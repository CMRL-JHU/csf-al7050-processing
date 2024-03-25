import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

path_file_input  = './pipeline_output/2.2-segmented data preErodeDilate.dream3d'
path_file_output = './pipeline_output/3-gakam threshold'
path_gakam = '/DataContainers/ImageDataContainer/CellFeatureData/KernelAvgMisorientations'
n_bins = 100
threshold = None

# read in gakam
with h5py.File(path_file_input, 'r') as f:
    gakam = f[path_gakam][...]

# if the threshold is set to none,
# guess at the appropriate threshold
# by finding the next most common value after zero
if threshold is None:
    threshold = stats.mode(gakam[gakam>0]).mode.item()
    print(f"Approximating threshold as: threshold > {threshold}")

# plot gakam
fig, axis = plt.subplots(1,1)
sns.histplot(
    data        = gakam ,
    stat        = 'probability',
    bins        = n_bins,
    ax          = axis
)

# find the ylimit
x = []; y = []
for patch in axis.patches:
    x.append(patch.get_x     ())
    y.append(patch.get_height())
y_limit = (0, 1.05*np.asarray(y)[np.asarray(x)>threshold].max())

# plot the threshold
threshold_nudge = (gakam.max()-gakam.min())/n_bins
axis.axvline(x = threshold+threshold_nudge, color = 'r', label=f"Threshold = {np.around(threshold, 4)}",)

axis.legend(fontsize=20)
axis.set_xlabel('GAKAM', fontsize=20)
axis.set_ylabel('PDF'  , fontsize=20)
axis.set_ylim  (y_limit)

fig.savefig(path_file_output)

