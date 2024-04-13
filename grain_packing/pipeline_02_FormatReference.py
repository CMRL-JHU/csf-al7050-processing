# Delete any extraneous data from the reference data.
# We don't know what was done to the reference data before this point,
# but we do know what information should be contained in a
# newly imported EBSD image, so we can remove all but those data arrays.
# This way, whatever variables were created in the reference file
# before this point will not be able to interfere with later scripts.

import utils_dream3d
import h5py
import os
from shutil import copyfile

path_input            = "pipeline_input/6-yz_small_cleaned_grains.dream3d"
path_output           = "pipeline_output/2-reference.dream3d"
path_CellData         = "/DataContainers/ImageDataContainer/CellData"
path_CellFeatureData  = "/DataContainers/ImageDataContainer/CellFeatureData"

whitelist_CellData = [
    "BC",
    "BS",
    "Bands",
    "Error",
    "EulerAngles",
    "MAD",
    "Mask_Clear",
    "Phases",
    "X",
    "Y"
]

print("Copying reference data")
os.makedirs(path_output.rsplit('/',1)[0], exist_ok=True)
copyfile(path_input, path_output)

print("Removing unneeded data")
with h5py.File(path_output, "r+") as file_output:
    del file_output[path_CellFeatureData]
    for item in file_output[path_CellData]:
        if not item in whitelist_CellData:
            del file_output[path_CellData+"/"+item]
    
print("Creating *.xdmf file")
utils_dream3d.make_xdmf(path_output)
