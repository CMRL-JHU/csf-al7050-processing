import utils_dream3d
import h5py
from shutil import copyfile

path_input_statistics        = "pipeline_output/7-xy-statistics-StatsGenerator.dream3d"
path_input_boundaries        = "pipeline_output/8-subdomain.dream3d"
path_output                  = "pipeline_output/9-xy-prepack.dream3d"
path_CellData                = "/DataContainers/ImageDataContainer/CellData"
path_input_CellEnsembleData  = "/DataContainers/StatsGeneratorDataContainer/CellEnsembleData"
path_output_CellEnsembleData = "/DataContainers/ImageDataContainer/CellEnsembleData"

print("Copying boundary data")
copyfile(path_input_boundaries, path_output)

print("Removing uneeded data")
with h5py.File(path_output, "r+") as file_output:

    remove_list = [
        path_CellData+"/"+"Phases",
        path_output_CellEnsembleData
        ]
    
    for item in remove_list:
        path, name = item.rsplit("/",1)
        if name in file_output[path].keys():
            del file_output[path+"/"+name]

print("Copying statistics")
with h5py.File(path_input_statistics, "r") as file_input_statistics, h5py.File(path_output, "a") as file_output:
    file_output.copy(file_input_statistics[path_input_CellEnsembleData], path_output_CellEnsembleData+"Reference")
    
print("Creating *.xdmf file")
utils_dream3d.make_xdmf(path_output)
