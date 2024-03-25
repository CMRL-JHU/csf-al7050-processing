import numpy as np
import h5py
import os
import json
from shutil import copyfile
import utils_dream3d

### paths
# files
path_hdf5_file_input  = "./pipeline_output/15-synthetic_grains.dream3d"
path_text_file_input  = "./pipeline_output/16-eulerangles.txt"
path_hdf5_file_output = "./pipeline_output/17-synthetic_grains.dream3d"
path_json_file_output = "./15GrainSwapExport.json"
path_dream3d          = "../../DREAM3D"
# groups
path_hdf5_celldata        = "/DataContainers/ImageDataContainer/CellData"
path_hdf5_cellfeaturedata = "/DataContainers/ImageDataContainer/CellFeatureData"
# datasets
path_hdf5_featureids          = path_hdf5_celldata       +"/"+"FeatureIds"
path_hdf5_feature_eulerangles = path_hdf5_cellfeaturedata+"/"+"AvgEulerAngles"

# create a dream3d script that creates element eulerangles from feature eulerangles
def create_dream3d_script(path_hdf5_file_output, path_hdf5_featureids, path_hdf5_feature_eulerangles):

    _, name_datacontainer, name_celldata, name_featureids = path_hdf5_featureids.rsplit("/", 3)
    _, name_cellfeaturedata, name_eulerangles             = path_hdf5_feature_eulerangles.rsplit("/", 2)

    with h5py.File(path_hdf5_file_output, 'r') as file_output:

        data_json = \
        {
            "0": {
                "FilterVersion": "1.2.815",
                "Filter_Enabled": True,
                "Filter_Human_Label": "Read DREAM.3D Data File",
                "Filter_Name": "DataContainerReader",
                "Filter_Uuid": "{043cbde5-3878-5718-958f-ae75714df0df}",
                "InputFile": os.path.realpath(path_hdf5_file_output),
                "InputFileDataContainerArrayProxy": utils_dream3d.create_expected(file_output)["InputFileDataContainerArrayProxy"],
                "OverwriteExistingDataContainers": 0
            },
            "1": {
                "FilterVersion": "6.5.141",
                "Filter_Enabled": True,
                "Filter_Human_Label": "Convert Orientation Representation",
                "Filter_Name": "ConvertOrientations",
                "Filter_Uuid": "{e5629880-98c4-5656-82b8-c9fe2b9744de}",
                "InputOrientationArrayPath": {
                    "Attribute Matrix Name": name_cellfeaturedata,
                    "Data Array Name": name_eulerangles,
                    "Data Container Name": name_datacontainer
                },
                "InputType": 0,
                "OutputOrientationArrayName": "AvgQuats",
                "OutputType": 2
            },
            "2": {
                "CreatedArrayName": name_eulerangles,
                "FeatureIdsArrayPath": {
                    "Attribute Matrix Name": name_celldata,
                    "Data Array Name": name_featureids,
                    "Data Container Name": name_datacontainer
                },
                "FilterVersion": "1.2.815",
                "Filter_Enabled": True,
                "Filter_Human_Label": "Create Element Array from Feature Array",
                "Filter_Name": "CopyFeatureArrayToElementArray",
                "Filter_Uuid": "{99836b75-144b-5126-b261-b411133b5e8a}",
                "SelectedFeatureArrayPath": {
                    "Attribute Matrix Name": name_cellfeaturedata,
                    "Data Array Name": name_eulerangles,
                    "Data Container Name": name_datacontainer
                }
            },
            "3": {
                "CreatedArrayName": "AvgQuats",
                "FeatureIdsArrayPath": {
                    "Attribute Matrix Name": name_celldata,
                    "Data Array Name": name_featureids,
                    "Data Container Name": name_datacontainer
                },
                "FilterVersion": "1.2.815",
                "Filter_Enabled": True,
                "Filter_Human_Label": "Create Element Array from Feature Array",
                "Filter_Name": "CopyFeatureArrayToElementArray",
                "Filter_Uuid": "{99836b75-144b-5126-b261-b411133b5e8a}",
                "SelectedFeatureArrayPath": {
                    "Attribute Matrix Name": name_cellfeaturedata,
                    "Data Array Name": "AvgQuats",
                    "Data Container Name": name_datacontainer
                }
            },
            "4": {
                "FilterVersion": "1.2.815",
                "Filter_Enabled": True,
                "Filter_Human_Label": "Write DREAM.3D Data File",
                "Filter_Name": "DataContainerWriter",
                "Filter_Uuid": "{3fcd4c43-9d75-5b86-aad4-4441bc914f37}",
                "OutputFile": os.path.realpath(path_hdf5_file_output),
                "WriteTimeSeries": 0,
                "WriteXdmfFile": 1
            },
            "PipelineBuilder": {
                "Name": "15GrainSwapExport",
                "Number_Filters": 5,
                "Version": 6
            }
        }
        
        with open(path_json_file_output,"w") as f:
            json.dump(data_json, f, indent=4, separators=(",", ": "), sort_keys=True)
        
# unused here, but possibly useful for the future
# create element data from feature data
def create_element_data_from_feature_data():
    # Import element data
    featureids          = file_output[path_hdf5_featureids][...]
    element_eulerangles = file_output[path_hdf5_element_eulerangles][...]
    
    # find element geometry
    dims                = featureids.shape[:-1]
    
    # export element eulerangles
    for featureid in range(feature_eulerangles.shape[0]):
    
        mask = featureids == featureid
        mask = np.broadcast_to(mask, list(mask.shape[:-1])+[3])
        
        element_eulerangles_replacement = np.broadcast_to(feature_eulerangles[featureid], mask.shape)
        
        np.putmask(element_eulerangles, mask, element_eulerangles_replacement)
        
        print(featureid, " of ", feature_eulerangles.shape[0])

# load feature eulerangles from matlab generated text file
print(f"Loading feature\n   var : {path_hdf5_feature_eulerangles.rsplit('/',1)[-1]}\n   from: {path_text_file_input}")
feature_eulerangles = np.loadtxt(path_text_file_input, delimiter=",", usecols=range(3)).astype(np.float32)

# replace first row with zeros (if needed)
if not all(feature_eulerangles[0,:] == 0):
    feature_eulerangles[0,:] = [0,0,0]

# copy previous file containing crystallography
print(f"Copying reference data\n   from: {path_hdf5_file_input}\n   to  : {path_hdf5_file_output}")
path = path_hdf5_file_output.rsplit('/',1)[0]
if not os.path.exists(path):
    os.makedirs(path)
copyfile(path_hdf5_file_input, path_hdf5_file_output)

with h5py.File(path_hdf5_file_output, 'r+') as file_output:
    
    # error out if the shapes don't match.
    # otherwise there will be no way to reconcile with featureids
    if not file_output[path_hdf5_feature_eulerangles][...].shape == feature_eulerangles.shape:
        raise Exception(f" \
            Input shape of {feature_eulerangles.shape} \
            did not match existing shape of {file_output[path_hdf5_feature_eulerangles][...].shape} \
        ")

    # remove crystallography data
    for attributematrix in [path_hdf5_celldata, path_hdf5_cellfeaturedata]:
        for attributearray in [file_output[attributematrix+"/"+key].name for key in file_output[attributematrix]]:
            blacklist = ["quats","eulerangles","ipf"]
            for item in blacklist:
                if item in attributearray.lower():
                    print(f"Deleting {attributearray}")
                    del file_output[attributearray]
    
    # export feature eulerangles
    print(f"Exporting feature\n   var : {path_hdf5_feature_eulerangles.rsplit('/',1)[-1]}\n   to  : {path_hdf5_file_output}")
    # file_output[path_hdf5_feature_eulerangles][...] == feature_eulerangles
    utils_dream3d.insert_attribute_array(
        path_output           = path_hdf5_file_output,
        path_group            = path_hdf5_feature_eulerangles.rsplit("/",1)[ 0],
        name                  = path_hdf5_feature_eulerangles.rsplit("/",1)[-1],
        data                  = feature_eulerangles,
        dtype                 = "DataArray<float>",
        attribute_matrix_type = "CellFeature"
    )
    
# export element eulerangles
# create dream3d script to convert feature eulerangles to element eulerangles
print(f"Creating file: {path_json_file_output}")
create_dream3d_script(path_hdf5_file_output, path_hdf5_featureids, path_hdf5_feature_eulerangles)
# call dream.3d to create 3d microstructure
utils_dream3d.call_dream3d(path_dream3d, path_json_file_output)
# delete dream3d script
print(f"Removing file: {path_json_file_output}")
os.remove(path_json_file_output)
        
