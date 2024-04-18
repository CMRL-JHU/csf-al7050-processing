# UFGs are characterized by having recrystallized upon impact leading to a
# very low grain averaged kernel average misorientation (GAKAM) and grain size.
# The UFGs and CGs can be fairly effectively separated by finding a GAKAM
# threshold that divides the low valued UFGs from the higher valued CGs.
# Because the UFG size is approximately 1/10th the size of the average CG,
# a great deal of UFGs can appear in a very small volume and they can present
# on a GAKAM histogram as a very large spike near zero.
# This script attempts to find the thresheld by finding the minimum of the
# derivative of the GAKAM PDF.
# This problem is very data cleanliness dependent, so the choice of n_bins
# is of great importance. If it's determined that the derivative is just
# not doing a great job, threshold can be manually specified so that
# it can be found via trial and error through the gakam graphs.

import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

path_file_input  = './pipeline_output/2.2-segmented data preErodeDilate.dream3d'
path_file_output = './pipeline_output/3-gakam threshold.png'
path_gakam = '/DataContainers/ImageDataContainer/CellFeatureData/KernelAvgMisorientations'
n_bins = 100
threshold = None
xlabel = "GAKAM$^{\circ}$"
fontsize = 20
subplots_adjustment=dict(left=0.15, right=0.75, top=0.80, bottom=0.15, hspace=0.25)

# create a discrete pdf from the histogram of some data
def hist_to_points(data, n_bins=50):

    # create the figure
    fig, axis = plt.subplots(1,1)

    # plot the gakam histogram
    sns.histplot(
        data  = data,
        stat  = 'percent',
        bins  = n_bins,
        ax    = axis,
    )
    
    # gather the x and y locations from the histogram bins
    points = []
    for patch in axis.patches:
        points.append([patch.get_x(), patch.get_height()])
    points = np.asarray(points)

    # close the figure
    plt.close(fig)

    return points

# find derivative via finite differences
def find_derivative(data):

    x = data[:,0]
    y = data[:,1]

    dy = np.diff(y,1) / np.diff(x,1)
    dx = 0.5*(x[:-1]+x[1:])

    ddata = np.vstack((dx, dy)).T

    return ddata

# find the ufg/cg threshold via the derivative
def find_threshold(data, n_bins=50):

    points = hist_to_points(data, n_bins)
    points_diff = find_derivative(points)
    threshold = points_diff[points_diff[:,1] == min(points_diff[:,1])][0]

    return threshold[0]

def concat_legends(*axes):

    def get_handles(axis):
        return axis.get_legend().legend_handles
    def get_labels (axis):
        return [ text._text for text in axis.get_legend().texts]
    
    handles = []
    labels  = []
    for axis in axes:
        handles += get_handles(axis)
        labels  += get_labels (axis)

    return handles, labels

def plot_before(gakam, xlabel, threshold, path_output, n_bins=50, fontsize=20, subplots_adjustment=None):

    # create the figure
    fig, axis1 = plt.subplots(1,1)
    axis2 = axis1.twinx()

    # plot the gakam histogram
    sns.histplot(
        data  = gakam ,
        stat  = 'percent',
        bins  = n_bins,
        ax    = axis1,
        label = 'PDF (%)'
    )

    # set the ylimit to a little higher than the next highest bin after the threshold
    # otherwise the number of UFGs will dominate the graph to the point you can only
    # see the UFG bin
    x = []; y = []
    for patch in axis1.patches:
        x.append(patch.get_x     ())
        y.append(patch.get_height())
    y_limit = (0, 1.05*np.asarray(y)[np.asarray(x)>threshold].max())
    axis1.set_ylim  (y_limit)

    # plot the threshold as a vertical line
    axis1.axvline(x=threshold, color = 'r', label=f"T={np.around(threshold, 4)}",)

    # plot the derivative field
    points_diff = find_derivative(hist_to_points(gakam, n_bins))
    axis2.plot(points_diff[:,0], points_diff[:,1], label="$\dfrac{\partial Y}{\partial X}$", color='orange', linestyle='--')

    # continue the line colors from axis 1
    axis2._get_lines = axis1._get_lines
    axis2._get_patches_for_fill = axis1._get_patches_for_fill

    # label everything
    axis1.legend     (                                    fontsize =fontsize)
    axis1.set_xlabel (xlabel   ,                          fontsize =fontsize)
    axis1.set_ylabel ('PDF (%)',                          fontsize =fontsize)
    axis1.tick_params(axis='both', which='major',         labelsize=fontsize)
    axis2.legend     (                                    fontsize =fontsize)
    axis2.set_ylabel ("$\dfrac{\partial Y}{\partial X}$", fontsize =fontsize)
    axis2.tick_params(axis='both', which='major',         labelsize=fontsize)

    # create figure legend
    for axis in [axis1, axis2]:
        axis.get_legend().set_visible(False)
    fig.legend(*concat_legends(axis1, axis2))
    sns.move_legend(
        obj            = fig                     ,
        loc            = "upper center"          ,
        bbox_to_anchor = (0.5, 1.0)              ,
        ncol           = 3                       ,
        title          = None                    ,
        frameon        = False                   ,
        fontsize       = 20
    )

    # set axis scaling
    if subplots_adjustment is None:
        fig.tight_layout()
    else:
        fig.subplots_adjust(**subplots_adjustment)

    # save
    path, ext = path_output.rsplit('.',1)
    path_output = f"{path}_before.{ext}"
    os.makedirs(os.path.split(path_output)[0], exist_ok=True)
    fig.savefig(path_output)

    # close the figure
    plt.close(fig)

def plot_after(gakam, xlabel, threshold, path_output, n_bins=50, fontsize=20):

    # create the figure
    fig, axis = plt.subplots(1,1)

    # plot the gakam histogram
    sns.histplot(
        data  = gakam[gakam>threshold],
        stat  = 'percent',
        bins  = n_bins,
        ax    = axis,
        # label = 'PDF (%)'
    )

    # label everything
    # axis.legend     (                            fontsize =fontsize)
    axis.set_xlabel (xlabel   ,                  fontsize =fontsize)
    axis.set_ylabel ('PDF (%)',                  fontsize =fontsize)
    axis.tick_params(axis='both', which='major', labelsize=fontsize)

    # set axis scaling
    fig.tight_layout()

    # save
    path, ext = path_output.rsplit('.',1)
    path_output = f"{path}_after.{ext}"
    fig.savefig(path_output)

    # close the figure
    plt.close(fig)

if __name__ == "__main__":

    # read in gakam
    with h5py.File(path_file_input, 'r') as f:
        gakam = f[path_gakam][...]

    # if the threshold is set to none,
    # guess at the appropriate threshold
    if threshold is None:
        # by finding the minimum of the derivative
        threshold = find_threshold(gakam, n_bins)
        print(f"Approximating GAKAM threshold as: threshold > {threshold}")
    else:
        print(f"GAKAM manually set to: threshold > {threshold}")

    # plot histograms
    plot_before(gakam, xlabel, threshold, path_file_output, n_bins, fontsize, subplots_adjustment)
    plot_after (gakam, xlabel, threshold, path_file_output, n_bins, fontsize)

