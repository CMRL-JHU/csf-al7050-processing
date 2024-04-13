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

path_file_input  = './pipeline_output/2.2-segmented data preErodeDilate.dream3d'
path_file_output = './pipeline_output/3-gakam threshold.png'
path_gakam = '/DataContainers/ImageDataContainer/CellFeatureData/KernelAvgMisorientations'
n_bins = 100
threshold = None

# find derivative via finite differences
def find_derivative(data):

    x = data[:,0]
    y = data[:,1]

    dy = np.diff(y,1) / np.diff(x,1)
    dx = 0.5*(x[:-1]+x[1:])

    ddata = np.vstack((dx, dy)).T

    return ddata

# create a discrete pdf from the histogram of some data
def hist_to_points(data, n_bins=50):
    
    # plot the histogram and gather the x and y locations from it
    fig, axis = plot_hist(data, 'Before', n_bins)
    points = []
    for patch in axis.patches:
        points.append([patch.get_x(), patch.get_height()])
    points = np.asarray(points)
    plt.close(fig)

    return points

# find the ufg/cg threshold via the derivative
def find_threshold(data, n_bins=50):

    points = hist_to_points(data, n_bins)
    points_diff = find_derivative(points)
    threshold = points_diff[points_diff[:,1] == min(points_diff[:,1])][0]

    return threshold[0]

def plot_hist(data, title, n_bins=50):

    # plot gakam
    fig, axis = plt.subplots(1,1)
    sns.histplot(
        data  = data ,
        stat  = 'probability',
        bins  = n_bins,
        ax    = axis,
        label = 'GAKAM' 
    )

    # print labels
    axis.set_title(title)
    axis.legend(fontsize=20)
    axis.set_xlabel('GAKAM', fontsize=20)
    axis.set_ylabel('PDF'  , fontsize=20)

    return fig, axis

def plot_before(gakam, threshold, path_output, n_bins=50):

    fig, axis = plot_hist(gakam, 'Before', n_bins)

    # set the ylimit to a little higher than the next highest bin after the threshold
    # otherwise the number of UFGs will dominate the graph to the point you can only
    # see the UFG bin
    x = []; y = []
    for patch in axis.patches:
        x.append(patch.get_x     ())
        y.append(patch.get_height())
    y_limit = (0, 1.05*np.asarray(y)[np.asarray(x)>threshold].max())
    axis.set_ylim  (y_limit)

    # plot the threshold as a vertical line
    threshold_nudge = (gakam.max()-gakam.min())/(n_bins*2)
    axis.axvline(x = threshold+threshold_nudge, color = 'r', label=f"Threshold = {np.around(threshold, 4)}",)

    # save
    path, ext = path_output.rsplit('.',1)
    path_output = f"{path}_before.{ext}"
    fig.savefig(path_output)
    plt.close(fig)

def plot_after(gakam, threshold, path_output, n_bins=50):

    fig, axis = plot_hist(gakam[gakam>threshold], 'After', n_bins)

    # save
    path, ext = path_output.rsplit('.',1)
    path_output = f"{path}_after.{ext}"
    fig.savefig(path_output)
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
    plot_before(gakam, threshold, path_file_output, n_bins)
    plot_after (gakam, threshold, path_file_output, n_bins)

