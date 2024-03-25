# if you try to directly import the MDF points into dream3d,
# it will superimpose them onto the existing mckenzie plot
# that dream3d uses as the default MDF.
# this script creates creates the MDF points by superimposing
# the desired MDF with the negative of the existing dream3d MDF
# plot so that when dream3d superimposes the mckenzie plot, 
# you're left with the true desired MDF.

import h5py
import numpy as np
import matplotlib.pyplot as plt

path_input_file  = "./pipeline_output/5-feature_attributes_ebsd.dream3d"
path_output_file = "./pipeline_output/6-mdf.txt"
path_input_data = "/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList"
# for some reason, this filter only acknowledges 13 points
bins = 20
mckenzie_weight = 5000
data_weight = 10000
axis = [0,1,0]
plot = True

def import_data(path_input_file, path_input_data):
    
    with h5py.File(path_input_file,'r') as file:
        return file[path_input_data][...]
        
def get_data_points(data, bins, plot=False):
    
    mckenzie_x_max = 62.4
    mckenzie_x_min =  2.4
    x_spacing      = (mckenzie_x_max-mckenzie_x_min)/(bins-1)
    xlim           = [mckenzie_x_min-x_spacing/2, mckenzie_x_max+x_spacing/2]

    weights = np.ones_like(data)/len(data)
    HIST_BINS = np.linspace(xlim[0], xlim[1], bins)
    heights, edges, _ = plt.hist(x=data,bins=HIST_BINS,weights=weights,edgecolor="black",label="data")
    
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

data = import_data(path_input_file, path_input_data)

x_data, y_data = get_data_points(data, bins, plot=plot)
y_data *= data_weight
if plot:
    plt.plot(x_data, y_data, marker="o", label="data")
    plt.show(block=False)

x_mckenzie, y_mckenzie = get_mckenzie_points(x_data)
y_mckenzie *= mckenzie_weight
if plot:
    plt.plot(x_mckenzie, y_mckenzie, marker="o", label="mckenzie")
    plt.show(block=False)

x, y = [x_data, y_data-y_mckenzie]
if plot:
    plt.plot(x, y, marker="o", label="superimposed")
    plt.show(block=False)
    
export_data(x, y, axis, path_output_file)

if plot:
    plt.legend()
    input("Press enter to continue...")