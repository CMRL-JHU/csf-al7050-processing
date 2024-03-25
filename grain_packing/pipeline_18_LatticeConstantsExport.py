import os
from shutil import copyfile
import h5py
import numpy as np
import utils_dream3d

### User input
## Paths
# Files
path_file_input_latticeconstants  = "./pipeline_output/5-feature_attributes_ebsd.dream3d"
path_file_input_dream3d           = "./pipeline_output/17-synthetic_grains.dream3d"
path_file_output                  = "./pipeline_output/18-synthetic_grains.dream3d"
# HDF5
path_CellEnsembleData = "/DataContainers/ImageDataContainer/CellEnsembleData"
path_latticeconstants = path_CellEnsembleData+"/"+"LatticeConstants"
## Lattice constants to keep
lattice_constant_indecies = [0,1]

print(f"Copying reference data from {path_file_input_dream3d.rsplit('/',1)[-1]}")
os.makedirs(path_file_output.rsplit('/',1)[0], exist_ok=True)
copyfile(path_file_input_dream3d, path_file_output)

print(f"Importing lattice constants from {path_file_input_latticeconstants.rsplit('/',1)[-1]}")
with h5py.File(path_file_input_latticeconstants, "r") as file_input:
    data = file_input[path_latticeconstants][lattice_constant_indecies,:].astype(np.float32)

print(f"Exporting lattice constants from {path_file_output.rsplit('/',1)[-1]}")
utils_dream3d.insert_attribute_array(
        path_output           = path_file_output,
        path_group            = path_latticeconstants.rsplit('/',1)[0],
        name                  = path_latticeconstants.rsplit('/',1)[-1],
        data                  = data,
        dtype                 = "DataArray<float>",
        attribute_matrix_type = "CellEnsemble"
    )
    
utils_dream3d.make_xdmf(path_file_output)