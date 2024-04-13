# writing to hdf5 files through matlab such that the data is still readable in dream3d
# is possible, but challenging.
# instead, we can import the text file containing the euler angles determined
# by the grain swapping algorithm and feed them into the dream3d file using python.

import numpy as np
import h5py
import os
from shutil import copyfile
import utils_dream3d

### paths
# files
path_hdf5_file_input  = "./pipeline_output/15-synthetic_grains.dream3d"
path_text_file_input  = "./pipeline_output/16-eulerangles.txt"
path_hdf5_file_output = "./pipeline_output/17-synthetic_grains.dream3d"
# groups
path_hdf5_celldata        = "/DataContainers/ImageDataContainer/CellData"
path_hdf5_cellfeaturedata = "/DataContainers/ImageDataContainer/CellFeatureData"
# datasets
path_hdf5_featureids          = path_hdf5_celldata       +"/"+"FeatureIds"
path_hdf5_feature_eulerangles = path_hdf5_cellfeaturedata+"/"+"AvgEulerAngles"

### vars
blacklist_name_dataarray = ["quats","eulerangles","ipf"]

#################################################################################################
############################################ EXECUTE ############################################
#################################################################################################

# load feature eulerangles from matlab generated text file
print(f"Loading feature\n   var : {path_hdf5_feature_eulerangles.rsplit('/',1)[-1]}\n   from: {path_text_file_input}")
feature_eulerangles = np.loadtxt(path_text_file_input, delimiter=",", usecols=range(3)).astype(np.float32)

# replace first row with zeros (if needed)
if not all(feature_eulerangles[0,:] == 0):
    feature_eulerangles[0,:] = [0,0,0]

# copy previous file containing crystallography
print(f"Copying reference data\n   from: {path_hdf5_file_input}\n   to  : {path_hdf5_file_output}")
os.makedirs(path_hdf5_file_output.rsplit('/',1)[0], exist_ok=True)
copyfile(path_hdf5_file_input, path_hdf5_file_output)

with h5py.File(path_hdf5_file_output, 'r+') as file_output:
    
    # error out if the shapes don't match.
    # otherwise there will be no way to reconcile with featureids
    if not file_output[path_hdf5_feature_eulerangles][...].shape == feature_eulerangles.shape:
        raise Exception(f" \
            Input shape of {feature_eulerangles.shape} \
            did not match existing shape of {file_output[path_hdf5_feature_eulerangles][...].shape} \
        ")

    # remove existing crystallography data
    print(f"Removing current crystallographic data arrays")
    for attributematrix in [path_hdf5_celldata, path_hdf5_cellfeaturedata]:
        for attributearray in file_output[attributematrix]:
            for item in blacklist_name_dataarray:
                if item.lower() in attributearray.lower():
                    print(f"   Deleting {attributematrix}/{attributearray}")
                    del file_output[f"{attributematrix}/{attributearray}"]
    
# export grain averaged eulerangles
print(f"Exporting data array\n   var : {path_hdf5_feature_eulerangles}\n   to  : {path_hdf5_file_output}")
utils_dream3d.insert_attribute_array(
    path_output           = path_hdf5_file_output,
    path_group            = path_hdf5_feature_eulerangles.rsplit("/",1)[ 0],
    name                  = path_hdf5_feature_eulerangles.rsplit("/",1)[-1],
    data                  = feature_eulerangles,
    dtype                 = "DataArray<float>",
    attribute_matrix_type = "CellFeature"
)

# export voxel eulerangles
print(f"Exporting data array\n   var : {path_hdf5_featureids.rsplit('/',1)[0]}/{path_hdf5_feature_eulerangles.rsplit('/',1)[-1]}\n   to  : {path_hdf5_file_output}")
with h5py.File(path_hdf5_file_output, 'r') as file_output:
    featureids = file_output[path_hdf5_featureids][...]
utils_dream3d.insert_attribute_array(
    path_output           = path_hdf5_file_output,
    path_group            = path_hdf5_featureids.rsplit('/',1)[0],
    name                  = path_hdf5_feature_eulerangles.rsplit("/",1)[-1],
    data                  = utils_dream3d.create_element_data_from_feature_data(featureids, feature_eulerangles),
    dtype                 = utils_dream3d.get_datatype_dream3d(feature_eulerangles.dtype),
    attribute_matrix_type = "Cell"
)

print('Creating xdmf file')
utils_dream3d.make_xdmf(path_hdf5_file_output)