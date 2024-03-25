import h5py
import numpy as np
import os
from shutil import copyfile
import utils_dream3d

path_input = "pipeline_output/10-synthetic_grains.dream3d"
path_output = "pipeline_output/11-synthetic_grains_inserted_improved_mask.dream3d"
path_Geometry = "/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY"
path_CellData = "/DataContainers/ImageDataContainer/CellData"
path_featureids = path_CellData+"/"+"FeatureIds"
path_mask = path_CellData+"/"+"Mask"

# the volumetric fraction of the grain inside of the mask
# if more that this fraction is inside, mask over the whole grain 
threshold = 0.1

print(f"Creating new mask. Keeping grains where >= {threshold*100}% is inside mask.")
with h5py.File(path_input,"r") as file_input:
    featureids  = file_input[path_featureids][...]
    mask        = file_input[path_mask][...]

featureids_masked, mask_updated = utils_dream3d.remove_masked_grains(featureids,mask,threshold)

padding = "   "
print("Mask Morphology Changes:")
print(padding+"original mask volumetric fraction: ",np.count_nonzero(mask)/featureids.size)
print(padding+"updated mask volumetric fraction : ",np.count_nonzero(mask_updated)/featureids.size)

print("Copying reference data")
os.makedirs(path_output.rsplit('/',1)[0], exist_ok=True)
copyfile(path_input, path_output)

print("Exporting new mask")
# delete attribute array if it exists already
with h5py.File(path_output,"r+") as file_output:
    if "mask_updated" in file_output[path_CellData]:
        del file_output[path_CellData+"/"+"mask_updated"]
# write out new mask attribute array
utils_dream3d.insert_attribute_array(
    path_output           = path_output      ,
    path_group            = path_CellData    ,
    name                  = "mask_updated"   ,
    data                  = mask_updated     ,
    dtype                 = "DataArray<bool>",
    attribute_matrix_type = "Cell"
)

print("Creating *.xdmf file")
utils_dream3d.make_xdmf(path_output)
