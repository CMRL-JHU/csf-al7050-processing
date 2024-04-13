# This script copies the material lattice constants from the EBSD file
# to the SEVM. Some phases should have disappeared during the cleaning stages
# of the pipeline up to this point, so the script attempts to find the 
# lattice constants it should copy by finding the indices of matching
# phase names.
# Unfortunately the "StatsGenerator" filter in pipeline 7 has also renamed some
# of these, so it will never be able to match up for the example script.
# instead, type in the indices that match pretty closely.
# EX 1:                      EX 2:                      EX 3:                     
#   EBSD :                     EBSD :                     EBSD :                  
#      0) Invalid Phase           0) Invalid Phase           0) Invalid Phase     
#      1) Aluminium               1) Aluminium               1) Zinc              
#      2) Zinc                    2) Zinc                    2) Magnesium         
#      3) Magnesium               3) Magnesium               3) Mg Zn2            
#      4) Mg Zn2                  4) Mg Zn2                  
#   SEVM :                     SEVM :                     SEVM :                  
#      0) Unknown Phase Type      0) Aluminum                0) Unknown Phase Type
#      1) Aluminum                1) Unknown Phase Type      1) Aluminum          
#   Correct indices are: 0,1   Correct indices are: 1,0   Correct indices are: None



import os
from shutil import copyfile
import h5py
import numpy as np
import utils_dream3d

### User input
## Paths
# Files
path_file_input_ebsd  = "./pipeline_output/3-feature_attributes_ebsd.dream3d"
path_file_input_sevm  = "./pipeline_output/17-synthetic_grains.dream3d"
path_file_output      = "./pipeline_output/18-synthetic_grains.dream3d"
# HDF5
path_phasename_ebsd = "/DataContainers/ImageDataContainer/CellEnsembleData/MaterialName"
path_phasename_sevm = "/DataContainers/ImageDataContainer/CellEnsembleData/PhaseName"
path_latticeconstants_ebsd = "/DataContainers/ImageDataContainer/CellEnsembleData/LatticeConstants"

def find_matching_indices(list1, list2):

    map_value_index = { value: index for index, value in enumerate(list1) }
    matching = [
        (index, map_value_index[value])
        for index, value 
        in enumerate(list2)
        if value in map_value_index
    ]
    return matching

def import_latticeconstants(
        path_file_input_ebsd,
        path_file_input_sevm,
        path_phasename_ebsd ,
        path_phasename_sevm ,
        padding = "   "
    ):

    # find the indices of matching phase names
    with h5py.File(path_file_input_ebsd, "r") as file_ebsd, h5py.File(path_file_input_sevm, "r") as file_sevm:
        phasenames_ebsd = [name.decode() for name in file_ebsd[path_phasename_ebsd][...]]
        phasenames_sevm = [name.decode() for name in file_sevm[path_phasename_sevm][...]]
    matching_indices = find_matching_indices(phasenames_ebsd, phasenames_sevm)

    # print matching indices for verification
    newline = f"\n{padding*2}"
    print(f"{padding}EBSD :\n{padding*2}{newline.join([str(i)+') '+name for i, name in enumerate(phasenames_ebsd)])}")
    print(f"{padding}SEVM :\n{padding*2}{newline.join([str(i)+') '+name for i, name in enumerate(phasenames_sevm)])}")
    if len(matching_indices) == 0:
        print(f"{padding}Match:\n{padding*2}None")
    else:
        print(f"{padding}Match:\n{padding*2}{newline.join([str(i)+' '+phasenames_ebsd[i] for i, j in matching_indices])}")

    # error out
    if len(phasenames_sevm) > len(phasenames_ebsd):
        raise ValueError("More phoses appear in SEVM than in EBSD which is impossible")

    # if no matches found, get the indices manually
    if len(matching_indices) == 0:
        while True:
            indices = input(f"{padding}No matching indices found, please specify which ebsd phases to keep via comma separated indices: ")
            try:
                latticeconstants_indices = [int(index) for index in indices.split(',')]
                break
            except ValueError as e:
                print(repr(e))
    # otherwise use the matching phase names to find the lattice constants
    else:
        latticeconstants_indices = [i for i, j in matching_indices]

    # get the lattice constants for the sevm
    with h5py.File(path_file_input_ebsd, "r") as file_ebsd:
        latticeconstants_ebsd = file_ebsd[path_latticeconstants_ebsd][...]
    lattice_constants_sevm = latticeconstants_ebsd[latticeconstants_indices,:]

    return lattice_constants_sevm

# find path for lattice constants in the sevm
path_latticeconstants_sevm = f"{path_phasename_sevm.rsplit('/',1)[0]}/{path_latticeconstants_ebsd.rsplit('/',1)[-1]}"

print(f"Copying reference data from {path_file_input_sevm.rsplit('/',1)[-1]}")
os.makedirs(path_file_output.rsplit('/',1)[0], exist_ok=True)
copyfile(path_file_input_sevm, path_file_output)

print(f"Finding lattice constants")
lattice_constants_sevm = import_latticeconstants(
    path_file_input_ebsd,
    path_file_input_sevm,
    path_phasename_ebsd ,
    path_phasename_sevm
)

print(f"Exporting lattice constants from {path_file_output.rsplit('/',1)[-1]}")
utils_dream3d.insert_attribute_array(
        path_output           = path_file_output,
        path_group            = path_latticeconstants_sevm.rsplit('/',1)[0],
        name                  = path_latticeconstants_sevm.rsplit('/',1)[-1],
        data                  = lattice_constants_sevm,
        dtype                 = "DataArray<float>",
        attribute_matrix_type = "CellEnsemble"
    )
    
utils_dream3d.make_xdmf(path_file_output)