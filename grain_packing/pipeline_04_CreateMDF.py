# if you try to directly import the MDF points into dream3d,
# it will superimpose them onto the existing mckenzie plot
# that dream3d uses as the default MDF.
# this script creates the MDF points by superimposing
# the desired MDF with the negative of the existing dream3d MDF
# plot so that when dream3d superimposes the mckenzie plot, 
# it is cancelled out and you're left with the true desired MDF.

import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

path_input_file  = "./pipeline_output/3-feature_attributes_ebsd.dream3d"
path_output_file = "./pipeline_output/4-mdf.txt"
path_input_data = "/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList"
# for some reason, this filter only acknowledges 13 points
bins = 20
mckenzie_weight = 5000
data_weight = 10000
axis = [0,1,0]
fontsize=20
subplots_adjustment=dict(left=0.20, right=0.95, top=0.75, bottom=0.15, hspace=0.25)

def import_data(path_input_file, path_input_data):
    
    with h5py.File(path_input_file,'r') as file:
        return file[path_input_data][...]
        
def get_data_points(data, bins):
    
    mckenzie_x_max = 62.4
    mckenzie_x_min =  2.4
    x_spacing      = (mckenzie_x_max-mckenzie_x_min)/(bins-1)
    xlim           = [mckenzie_x_min-x_spacing/2, mckenzie_x_max+x_spacing/2]

    weights = np.ones_like(data)/len(data)
    HIST_BINS = np.linspace(xlim[0], xlim[1], bins)

    fig, axis = plt.subplots(1,1)
    heights, edges, _ = plt.hist(x=data,bins=HIST_BINS,weights=weights,edgecolor="black",label="data")
    plt.close(fig)
    
    centers = 0.5*(edges[1:] + edges[:-1])
    
    return centers, heights
    
def get_mckenzie_points(data_x):
    
    mckenzie_approximation = \
        np.array([
            [ 2.4, 0.0020],
            [ 7.4, 0.0060],
            [12.4, 0.0170],
            [17.4, 0.0325],
            [22.4, 0.0510],
            [27.4, 0.0745],
            [32.4, 0.1080],
            [37.4, 0.1365],
            [42.4, 0.1750],
            [47.4, 0.1755],
            [52.4, 0.1330],
            [57.4, 0.0825],
            [62.4, 0.0085]
        ])
        
    mckenzie_y = np.zeros(data_x.shape)
    for i, data_x_i in enumerate(data_x):
    
        mckenzie_lower_x, mckenzie_lower_y = mckenzie_approximation[mckenzie_approximation[:,0] <= data_x_i, :][-1]
        mckenzie_upper_x, mckenzie_upper_y = mckenzie_approximation[mckenzie_approximation[:,0] >= data_x_i, :][ 0]
        
        mckenzie_y[i] = interpolate(data_x_i, mckenzie_lower_x, mckenzie_upper_x, mckenzie_lower_y, mckenzie_upper_y)
    
    return data_x, mckenzie_y
    
def interpolate(x, x_lower, x_upper, y_lower, y_upper):

    if y_upper == y_lower:
        return y_upper

    return (y_upper-y_lower)/(x_upper-x_lower)*(x-x_lower) + y_lower
    
def superimpose(data, weights):

    data_new = np.zeros(data[0].shape)
    for data_i, weight in zip(data, weights):
        data_new += data_i*weight

    return data_new
    
def export_data(xs, ys, axis, path_output_file):
    
    string = str(len(ys))+"\n"
    for x, y in zip(xs, ys):
        string += '%.10f %.10f %.10f %.10f %.10f\n' %(x, *axis, y)
        
    with open(path_output_file, 'w') as f:
        f.write(string)

def plot_image(
    x_desired,
        y_desired,
        x_mckenzie,
        y_mckenzie,
        x_superimposed,
        y_superimposed,
        path_output_file,
        subplots_adjustment=None
):

    fig, axis = plt.subplots(1,1)

    ### plot everything
    # desired
    plt.plot(x_desired, y_desired, marker="o", label="Desired")
    # mckenzie
    plt.plot(x_mckenzie, y_mckenzie, marker="o", label="McKenzie")
    # superimposed
    plt.plot(x_superimposed, y_superimposed, marker="o", label="Superimposed")

    # print everything
    axis.legend(                                        fontsize =fontsize)
    axis.set_xlabel ("Misorientation Angle ($^\circ$)", fontsize =fontsize)
    axis.set_ylabel ("Weighted PDF"                   , fontsize =fontsize)
    axis.tick_params(axis='both', which='major'       , labelsize=fontsize)

    # move the legend
    axis.get_legend().set_visible(False)
    handles = axis.get_legend().legend_handles
    labels  = [ text._text for text in axis.get_legend().texts]
    fig.legend(handles, labels)
    sns.move_legend(
        obj            = fig           ,
        loc            = "upper center",
        bbox_to_anchor = (0.5, 1.0)    ,
        ncol           = 2             ,
        title          = None          ,
        frameon        = False         ,
        fontsize       = 20
    )

    # set axis scaling
    if subplots_adjustment is None:
        fig.tight_layout()
    else:
        fig.subplots_adjust(**subplots_adjustment)

    # save
    path, ext = path_output_file.rsplit('.',1)
    path_output = f"{path}.png"
    os.makedirs(os.path.split(path_output)[0], exist_ok=True)
    fig.savefig(path_output)

if __name__ == "__main__":

    # import misorientation list from dream3d file
    data = import_data(path_input_file, path_input_data)

    # import the desired  MDF curve from the histogram of
    # the imported misorientation data
    x_desired, y_desired = get_data_points(data, bins)
    y_desired *= data_weight

    # create linear approximation of mckenzie curve
    # with the same number of points as the desired curve
    x_mckenzie, y_mckenzie = get_mckenzie_points(x_desired)
    y_mckenzie *= mckenzie_weight

    # create superimposed curve
    x_superimposed, y_superimposed = [x_desired, y_desired-y_mckenzie]

    # save superimposed curve to file
    export_data(x_superimposed, y_superimposed, axis, path_output_file)
    plot_image(
        x_desired,
        y_desired,
        x_mckenzie,
        y_mckenzie,
        x_superimposed,
        y_superimposed,
        path_output_file,
        subplots_adjustment
    )