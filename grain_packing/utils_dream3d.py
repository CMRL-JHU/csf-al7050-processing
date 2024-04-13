import h5py, json, os, sys, subprocess, platform, math
import numpy as np
import shlex

#make a spherical mask in a volume
def make_mask_sphere(dims, center=None, radius=None, origin=[0,0,0]):
    
    if center is None:
        center = [.5*dim for dim in dims]
    if radius is None:
        radius = .25*min(dims)

    coordinates = np.mgrid[[range(dim) for dim in dims]]
    r = np.sqrt(sum([(x-x0)**2 for x, x0 in zip(coordinates, center)]))
    r[r > radius] = 0
    r[r > 0     ] = 1
    
    dtype = "DataArray<bool>"
    r = r.reshape(dims+[1])
    r = r.astype(get_datatype(dtype))
    
    return r
    
#make a hyperellipsoidal mask in a volume
def make_mask_hyperellipsoid(dims, center=None, rotation_matrix=None, diameters=[1,1,1], exponent=2):
    
    if center is None:
        center = [.5*dim for dim in dims]
    if rotation_matrix is None:
        rotation_matrix = np.eye(3)
    
    x  = np.asarray(np.mgrid[[range(dim) for dim in dims]]) 
    x0 = np.broadcast_to( np.asarray(center).reshape(3,1,1,1), [3]+list(dims) )
    coordinates = np.einsum('ab,bijk', rotation_matrix, x-x0)
        
    r = np.sqrt( sum( [ (x/d)**exponent for x, d in zip(coordinates, diameters) ] ) )
    
    r[r > 1] = 0
    r[r > 0] = 1
    
    dtype = "DataArray<bool>"
    r = r.reshape(dims+[1])
    r = r.astype(get_datatype(dtype))
    
    return r

#make a spherical mask in a volume
def make_mask_superellipsoid(dims, epsilon1, epsilon2, center=None, rotation_matrix=None, radius=None, origin=[0,0,0]):
    
    if center is None:
        center = [.5*dim for dim in dims]
    if radius is None:
        radius = [.25*min(dims)]*3
    if rotation_matrix is None:
        rotation_matrix = np.eye(3)
    if epsilon1 is None:
        epsilon1 = 1.0
    if epsilon2 is None:
        epsilon2 = 1.0
    
    print(dims)
    
    #coordinates = np.asarray(np.mgrid[[range(dim) for dim in dims]])
    coordinates = np.mgrid[[range(dim) for dim in dims]] 
    print(coordinates)
    print(coordinates.shape)
    
    x0 = np.broadcast_to( np.asarray(center).reshape(3,1,1,1), [3]+list(dims) )
    print(x0)
    print(x0.shape)
    x = np.einsum('ab,bijk', rotation_matrix.as_matrix(), coordinates-x0)
    print(x)
    print(x.shape)
    print(radius)
#    exit()
#    r_a = np.sign(x[0])*( np.abs((x[0])/radius) )**(2.0/epsilon2)
#    r_b = np.sign(x[1])*( np.abs((x[1])/radius) )**(2.0/epsilon2)
#    r_c = np.sign(x[2])*( np.abs((x[2])/radius) )**(2.0/epsilon1)
    r_a = ( np.abs(x[2]/radius[2]) )**(2.0/epsilon2)
    r_b = ( np.abs(x[1]/radius[1]) )**(2.0/epsilon2)
    r_c = ( np.abs(x[0]/radius[0]) )**(2.0/epsilon1)
    r = (r_a +  r_b)**(epsilon2/epsilon1) + r_c
    
    r[r > 1] = 0
    r[r > 0     ] = 1
#    print("r ______________ ")
#    print(r)
#    print("r ______________ ")
    dtype = "DataArray<bool>"
    r = r.reshape(dims+[1])
    r = r.astype(get_datatype(dtype))
    
    return r
    
#make a cubic mask in a volume
def make_mask_cube(dims_outer, center=None, dims_inner=None):
    
    if center is None:
        center     = [int(0.5*dim) for dim in dims_outer]
    if dims_inner is None:
        dims_inner = [int(0.5*dim) for dim in dims_outer]
        
    indecies_inner = tuple(np.ogrid[[range(mid-dim//2, mid+dim//2) for dim, mid in zip(dims_inner, center)]])
    cube_outer = np.zeros(dims_outer)
    cube_outer[indecies_inner] = 1
    
    dtype = "DataArray<bool>"
    cube_outer = cube_outer.reshape(dims_outer+[1])
    cube_outer = cube_outer.astype(get_datatype(dtype))
    
    return cube_outer

#returns numpy class of datatype
def get_datatype(datatype):
    if datatype in ["float32", format_string("DataArray<float>"   ), "DataArray<float>"   ]:
        return np.float32
    if datatype in ["int32"  , format_string("DataArray<int32_t>" ), "DataArray<int32_t>" ]:
        return np.int32
    if datatype in ["int32"  , format_string("DataArray<uint32_t>"), "DataArray<uint32_t>"]:
        return np.uint32
    if datatype in ["uint8"  , format_string("DataArray<bool>"    ), "DataArray<bool>"    ]:
        return np.uint8
    if datatype in ["object" , format_string("StringDataArray"    ), "StringDataArray"    ]:
        return None
        
def get_datatype_dream3d(datatype):
    if datatype in ["float32", np.float32]:
        return "DataArray<float>"
    if datatype in ["int32"  , np.int32]:
        return "DataArray<int32_t>"
    if datatype in ["int32"  , np.uint32]:
        return "DataArray<uint32_t>"
    if datatype in ["uint8"  , np.uint8]:
        return "DataArray<bool>"
    if datatype in ["object"]:
        return "StringDataArray"
    
# Find all filter ids in a DREAM.3D *.JSON file with name
def get_filter_ids(data_json, name_filter):

    filter_ids = []
    for filter_id in data_json.keys():
    
        if "Filter_Name" in data_json[filter_id].keys() \
            and data_json[filter_id]["Filter_Enabled"]  \
            and data_json[filter_id]["Filter_Name"] == name_filter:
            
            filter_ids += [filter_id]
            
    return filter_ids

def create_element_data_from_feature_data(featureids, data_feature_array):
    
    # find element geometry
    dims       = list(featureids.shape[:-1])
    components = data_feature_array.shape[-1]

    # assign the correct euler angles to each voxel
    data_cell_array = np.zeros((np.prod(dims+[components]),))
    for i, featureid in enumerate(featureids.flatten()):
        data_cell_array[i*components:(i+1)*components] = data_feature_array[featureid]
    data_cell_array = data_cell_array.reshape(dims+[components]).astype(data_feature_array.dtype)

    # VERY SLOW
    # # find the voxels belonging to this grain
    # # and assign the correct eulerangle into those voxels
    # cell_eulerangles      = np.zeros(dims+[components])
    # for featureid in range(feature_eulerangles.shape[0]):
    #     mask = featureids == featureid
    #     cell_eulerangles[mask.squeeze(), :] = feature_eulerangles[featureid]

    return data_cell_array

#slicegan has no knowlegde of crystallography data, resolution, or naming schema, so it is imported from the .dream3d file
def import_resolution(ebsd_paths, name_planes, path_Geometry):

    size_voxels = []
    for ebsd_path in ebsd_paths:
            
        # Search through all data containers to find reference information
        with h5py.File(ebsd_path, 'r') as ebsd_file:

            # Each image is 2D, so the spacing in Z dimension is inconsequential
            size_voxels += [ebsd_file[path_Geometry+"/"+"SPACING"][:-1].tolist()]
        
    # If the input contains 3 planes, we are trying to reconstruct a 3D microstructure.
    # The side resolutions should be the same in order to reconstruct properly.
    # This function allows us to determine if the side resolutions are compatible
    # And gives us the common resolutions [x,y,z] if they exist
    if len(name_planes) == 3:
        size_voxels = find_common_global_value(name_planes, size_voxels)
        check_error_resolution(size_voxels)
    # If the input contains 1 plane, we're trying to reconstruct a 2D microstructure.
    # No compatibility check is required, and the resolution is given as [x,y,1]
    else:
        size_voxels = size_voxels+[1]

    return size_voxels
    
#filter a string with blacklist and whitelist
def filter_string(string, blacklist=None, whitelist=None):
    if not blacklist is None:
        for entry in blacklist:
            string = string.replace(entry, "")
    if not whitelist is None:
        string = "".join(list(filter(lambda char: char in set(whitelist), string)))
    return string  
        
#check if there is a common voxel size for each side. if there is not, then there will be no way to reconcile the dream3d data with the slicegan data at the end.
def check_error_resolution(size_voxels):
    if isinstance(size_voxels, str):
        raise ValueError(size_voxels)
    
#check if the crystallography matches exactly. if it does not, then there will be no way to reconcile the dream3d data with the slicegan data at the end.
def check_error_crystallography(ebsd_paths, path_CellEnsembleData, orientations_types):

    # If the reference data should contain multiple images crystallographic data
    if len(ebsd_paths) == 3 and is_crystallographic(orientations_types):

        # Open all the ebsd files
        ebsd_files = [h5py.File(ebsd_path, 'r') for ebsd_path in ebsd_paths]
        
        # Gather all the CellEnsembleData groups into a list
        CellEnsembleData = [file[path_CellEnsembleData] for file in ebsd_files]

        # Check if the groups are exactly equal
        result = groups_are_equal(CellEnsembleData)

        # Close all the ebsd files
        [file.close() for file in ebsd_files]
            
        if not result:
            raise RuntimeError("Crystallography data does not match between all 3 planes. No way to link crystallography")
            
#returns true if the value is a crystallographic data type
def is_crystallographic(orientations_types):
    if type(orientations_types) in [list, tuple]:
        result = []
        for orientations_type in orientations_types:
            result += [is_crystallographic(orientations_type)]
        return any(result)
    return orientations_types=="Quats" or orientations_types=="EulerAngles"
    
#EBSD scans contain two directional dimensions, but dream3d fills in the third to voxelize them
#dream3d seems to always keep the images oriented in the xy plane, so the z axis is filled with 1
#dream3d also seems to organize its directional dimensions backwards (z,y,x), so we remove dim 0
def remove_empty_z_container(array):
    return array.reshape(tuple(list(array.shape)[1:]))
    
#Replace the empty z container to put arrays back into dream3d
def replace_empty_z_container(array):
    return array.reshape(tuple([1]+list(array.shape)))
    
#rotate the 2-D array 90 degrees
#dream 3d keeps images oriented in a local xy plane, so to rotate, just swap the xy axes
#then flip either the x or y data to rotate counter-clockwise or clockwise respectively
def rotate_90_degrees(array, direction=1, n=1):
    for i in range(n):
        if direction == -1:
            array = array.swapaxes(0,1)[:,::-1,:] 
        else:
            array = array.swapaxes(0,1)[::-1,:,:]
    return array
    
#rotate all the CellData arrays within dream3d
def rotate_all(path_input ,n, path_CellData, path_Geometry):

    with h5py.File(path_input, 'r+') as file_input:

        # rotate n times
        for i in range(n):
            
            # gather geometry data
            dims        = file_input[path_Geometry+"/"+"DIMENSIONS"][:]
            size_voxels = file_input[path_Geometry+"/"+"SPACING"   ][:]

            # swap x and y geometry data
            dims[0], dims[1] = dims[1], dims[0]
            size_voxels[0], size_voxels[1] = size_voxels[1], size_voxels[0]

            # replace geometry data
            file_input[path_CellData].attrs["TupleDimensions"] = dims.astype(np.uint64)
            file_input[path_Geometry+"/"+"DIMENSIONS"][:] = dims.astype(np.uint64)
            file_input[path_Geometry+"/"+"SPACING"   ][:] = size_voxels

            for item in file_input[path_CellData]:

                # replace geometry data
                dims_string = ','.join(['='.join([name,str(value)]) for name,value in zip(['x','y','z'],dims)])
                file_input[path_CellData+"/"+item].attrs["TupleDimensions"      ] = dims.astype(np.uint64)
                file_input[path_CellData+"/"+item].attrs["Tuple Axis Dimensions"] = dims_string

                # rotate data
                data = file_input[path_CellData+"/"+item][:]
                data = remove_empty_z_container(data)
                data = rotate_90_degrees(data, n=n)
                data = replace_empty_z_container(data)
                data = data.astype(get_datatype(file_input[path_CellData+"/"+item].attrs["ObjectType"]))
                file_input[path_CellData+"/"+item][:] = data

# Create an empty DREAM.3D file with just the bare minimum to be recognized by DREAM.3D
def create_file_dream3d(path_output):

    # make sure the path exists
    path = os.path.split(path_output)[0]
    if not os.path.exists(path):
        os.makedirs(path)

    with h5py.File(path_output, 'w') as file_output:
    
        # create hdf file and required data structure
        file_output.attrs["DREAM3D Version"] = format_string("1.2.815.6bed39e95")
        file_output.attrs["FileVersion"]     = format_string("7.0")
        file_output.create_group("DataContainerBundles")
        Pipeline = file_output.create_group("Pipeline")
        Pipeline.attrs["Current Pipeline"] = format_string("Pipeline")
        Pipeline.attrs["Pipeline Version"] = np.int32([2])
        padding = 4*" "
        s  = 0*padding + "{"                                       + "\n"
        s += 1*padding + "\"PipelineBuilder\": " + "{"             + "\n"
        s += 2*padding + "\"Name\": "            + "\"Pipeline\"," + "\n"
        s += 2*padding + "\"Number_Filters\": "  + "0,"            + "\n"
        s += 2*padding + "\"Version\": "         + "6"             + "\n"
        s += 1*padding + "}"                                       + "\n"
        s += 0*padding + "}"
        Pipeline.create_dataset("Pipeline", data=s)

# Geometry should be a dictionary formatted as:
# {
#    "dims"       : [x,y,z]
#    "origin"     : [x,y,z] # normally [0,0,0]
#    "size_voxels": [x,y,z]
# }
def insert_geometry(
    path_output,
    geometry,
    path_Geometry="/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY"
    ):

    with h5py.File(path_output, 'a') as file_output:

        # export geometry data necessary for preserving dimensionality
        Geometry = file_output.create_group(path_Geometry)
        Geometry.attrs["GeometryName"]          = format_string("ImageGeometry")
        Geometry.attrs["GeometryType"]          = np.uint32([0])
        Geometry.attrs["GeometryTypeName"]      = format_string("ImageGeometry")
        Geometry.attrs["SpatialDimensionality"] = np.uint32([len(geometry["dims"])])
        Geometry.attrs["UnitDimensionality"]    = np.uint32([len(geometry["dims"])])
        Geometry.create_dataset("DIMENSIONS", data=np.asarray(geometry["dims"       ]).astype(np.int64  ))
        Geometry.create_dataset("ORIGIN"    , data=np.asarray(geometry["origin"     ]).astype(np.float32))
        Geometry.create_dataset("SPACING"   , data=np.asarray(geometry["size_voxels"]).astype(np.float32))
        
# Crystallography information comes from EBSD files of many types
# While you could create your own EBSD data from an EBSD file (.ang, .ctf, etc...)
# DREAM.3D handles this pretty well already, so it's likely best to import it with DREAM.3D
# This function can be used to copy already formatted crystallography information from a reference *.DREAM3D file
def copy_crystallography(
    path_dream3d_output,
    path_dream3d_input,
    path_CellEnsembleData_output="/DataContainers/ImageDataContainer/CellEnsembleData",
    path_CellEnsembleData_input="/DataContainers/ImageDataContainer/CellEnsembleData"
    ):
    
    with h5py.File(path_dream3d_input, 'r') as file_dream3d_input, h5py.File(path_dream3d_output, 'a') as file_dream3d_output:
    
        # export cell ensemble data necessary for crystallographic analysis
        file_dream3d_output.copy(file_dream3d_input[path_CellEnsembleData_input], path_CellEnsembleData_output)

# # Insert dream3d formatted attribute array with specified data into specified hdf5 group
# # Should should be a numpy array of shape: [z,y,x,components] or [x,components]
# # A missing component dimension can be inferred, but it's best to be explicit
# def insert_attribute_array(
#     path_output,
#     name,
#     data,
#     dtype="DataArray<float>",
#     path_CellData="/DataContainers/ImageDataContainer/CellData",
#     path_Geometry="/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY"
#     ):

#     ## Find dimension and component information
#     # Assume the data is 2 or 3 dimensional
#     if len(data.shape) > 2:
#         # Assume the data is missing the component dimension
#         if not len(data.shape) > 3:
#             data = data.reshape(np.append(data.shape,1).astype(tuple))
#     # Assume the data is 1 dimensional
#     else:
#         # Assume the data is missing the component dimension
#         if not len(data.shape) > 1:
#             data = data.reshape(np.append(data.shape,1).astype(tuple))
#     dims                  = list(data.shape[:-1][::-1])
#     components            = [data.shape[-1]]
#     tuple_axis_dimensions = ','.join(['='.join([name,str(value)]) for name,value in zip(['x','y','z'],dims)])
        
#     with h5py.File(path_output, 'a') as file_output:
        
#         # create the CellData group and its required attributes if it does not exist
#         group = file_output.require_group(path_CellData)
#         attributes = {
#             "AttributeMatrixType": np.uint32([len(file_output[path_Geometry+"/"+"DIMENSIONS"][...])]),
#             "TupleDimensions"    : np.uint64(file_output[path_Geometry+"/"+"DIMENSIONS"][...])
#         }
#         for key, val in zip(attributes.keys(), attributes.values()):
#             if not key in [i for i in file_output[path_CellData].attrs.keys()]:
#                 group.attrs[key] = val
        
#         # export attribute array
#         dataset = group.create_dataset(name, data=data)
#         dataset.attrs["ComponentDimensions"]   = np.uint64(components)
#         dataset.attrs["DataArrayVersion"]      = np.int32([2])
#         dataset.attrs["ObjectType"]            = format_string(dtype)
#         dataset.attrs["Tuple Axis Dimensions"] = format_string(tuple_axis_dimensions)
#         dataset.attrs["TupleDimensions"]       = np.uint64(dims)

def insert_attribute_array(
    path_output,
    path_group,
    name,
    data,
    dtype,
    attribute_matrix_type="Cell",
    linked_numneighbors_dataset=None,
    tuple_dimensions=None
    ):

    ## Find dimension and component information
    if len(data.shape) == 1:
        dims                        = np.asarray(data.shape).astype(np.uint64)
        components                  = np.asarray(1).reshape((1,)).astype(np.uint64)
    else:
        dims                  = np.asarray(data.shape[:-1][::-1]).astype(np.uint64)
        components            = np.asarray(data.shape[-1]).reshape((1,)).astype(np.uint64)
    if not linked_numneighbors_dataset is None:
        dims = tuple_dimensions
        linked_numneighbors_dataset = format_string(linked_numneighbors_dataset)
    else:
        tuple_axis_dimensions = format_string(','.join(['='.join([name,str(value)]) for name,value in zip(['x','y','z'],dims.tolist())]))

    ## Create attributes for the attribute matrix
    group_attributes = {
        "AttributeMatrixType": np.uint32(get_attribute_matrix_id(attribute_matrix_type)),
        "TupleDimensions"    : dims
    }

    ## Create attributes common to all datasets
    dataset_attributes = {
        "ComponentDimensions"  : components,
        "DataArrayVersion"     : np.asarray(2).reshape((1,)).astype(np.int32),
        "ObjectType"           : format_string(dtype),
        "TupleDimensions"      : dims
    }
    # Fill in attributes for special exceptions
    if not linked_numneighbors_dataset is None:
        dataset_attributes["Linked NumNeighbors Dataset"] = linked_numneighbors_dataset
    else:
        dataset_attributes["Tuple Axis Dimensions"      ] = tuple_axis_dimensions
    
    ## Write out to file
    with h5py.File(path_output, 'a') as file_output:

        # create the group and its required attributes if it does not exist
        group = file_output.require_group(path_group)
        for key, val in zip(group_attributes.keys(), group_attributes.values()):
            if not key in [i for i in group.attrs.keys()]:
                group.attrs[key] = val
        
        # export attribute array
        dataset = group.create_dataset(name, data=data)
        for key, val in zip(dataset_attributes.keys(), dataset_attributes.values()):
            if not key in [i for i in dataset.attrs.keys()]:
                dataset.attrs[key] = val

def get_attribute_matrix_id(attribute_matrix_type):
    map_attribute_matrix_type_id = {
        "Vertex" 	    :  0,
        "Edge" 	        :  1,
        "Face" 	        :  2,
        "Cell" 	        :  3,
        "VertexFeature" :  4,
        "EdgeFeature"   :  5,
        "FaceFeature"   :  6,
        "CellFeature"   :  7,
        "VertexEnsemble":  8,
        "EdgeEnsemble"  :  9,
        "FaceEnsemble"  : 10,
        "CellEnsemble"  : 11,
        "MetaData" 	    : 12,
        "Generic" 	    : 13,
        "Unknown" 	    :999
    }
    return map_attribute_matrix_type_id[attribute_matrix_type]
        
def make_xdmf(dream3d_path, padding="   ", n_padding=0):
    
    xdmf_path    = dream3d_path.rsplit(".",1)[0]+".xdmf"
    dream3d_name = os.path.split(dream3d_path)[-1]

    def insert_datacontainer(datacontainer, dream3d_name, padding, n_padding):
    
        def insert_geometry(datacontainer, padding, n_padding):
        
            dimensions = " ".join(str(i+1) for i in datacontainer["_SIMPL_GEOMETRY"+"/"+"DIMENSIONS"][...][::-1].astype(int))
            n_dims     = str(len(datacontainer["_SIMPL_GEOMETRY"+"/"+"DIMENSIONS"][...])) 
            origin     = " ".join(str(i+0) for i in datacontainer["_SIMPL_GEOMETRY"+"/"+"ORIGIN"    ][...][::-1].astype(int))
            spacing    = " ".join(str(i+0) for i in datacontainer["_SIMPL_GEOMETRY"+"/"+"SPACING"   ][...][::-1])
        
            string = ""
            string += (0+n_padding)*padding+f"<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"{dimensions}\">"+"\n"
            string += (0+n_padding)*padding+"</Topology>"+"\n"
            string += (0+n_padding)*padding+"<Geometry Type=\"ORIGIN_DXDYDZ\">"+"\n"
            string += (1+n_padding)*padding+"<!-- Origin  Z, Y, X -->"+"\n"
            string += (1+n_padding)*padding+f"<DataItem Format=\"XML\" Dimensions=\"{n_dims}\">"+"\n"
            string += (2+n_padding)*padding+f"{origin}"+"\n"
            string += (1+n_padding)*padding+"</DataItem>"+"\n"
            string += (1+n_padding)*padding+"<!-- DxDyDz (Spacing/Resolution) Z, Y, X -->"+"\n"
            string += (1+n_padding)*padding+f"<DataItem Format=\"XML\" Dimensions=\"{n_dims}\">"+"\n"
            string += (2+n_padding)*padding+f"{spacing}"+"\n"
            string += (1+n_padding)*padding+"</DataItem>"+"\n"
            string += (0+n_padding)*padding+"</Geometry>"
            
            return string
            
        def insert_attributearray(attributearray, dream3d_name, padding, n_padding):
        
            datatypes_xdmf    = ["Float", "Int", "UChar", "UChar", "UChar"]
            datatypes_dream3d = ["DataArray<float>","DataArray<int32_t>","DataArray<bool>","DataArray<uint8_t>","DataArray<int8_t>"]
            
            path          = attributearray.name
            name          = path.rsplit("/",1)[-1]
            shape         = " ".join(str(i+0) for i in attributearray[...].shape)
            attributetype = "Scalar" if attributearray[...].shape[-1] == 1 else "Vector"
            try:
                datatype  = datatypes_xdmf[datatypes_dream3d.index(attributearray.attrs["ObjectType"]   .decode("utf-8"))]
                precision = 4 if attributearray.attrs["ObjectType"]   .decode("utf-8") == "DataArray<float>" else 1
            except AttributeError: # hacky workaround for matlab string attributes
                datatype  = datatypes_xdmf[datatypes_dream3d.index(attributearray.attrs["ObjectType"][0].decode("utf-8"))]
                precision = 4 if attributearray.attrs["ObjectType"][0].decode("utf-8") == "DataArray<float>" else 1
            
            string = ""
            string += (0+n_padding)*padding+f"<Attribute Name=\"{name}\" AttributeType=\"{attributetype}\" Center=\"Cell\">"+"\n"
            string += (1+n_padding)*padding+f"<DataItem Format=\"HDF\" Dimensions=\"{shape}\" NumberType=\"{datatype}\" Precision=\"{precision}\" >"+"\n"
            string += (2+n_padding)*padding+f"{dream3d_name}:{path}"+"\n"
            string += (1+n_padding)*padding+"</DataItem>"+"\n"
            string += (0+n_padding)*padding+"</Attribute>"
            
            return string
            
        name_datacontainer = datacontainer.name.rsplit("/",1)[-1]
            
        string = ""
        string += (0+n_padding)*padding+f"<!-- *************** START OF {name_datacontainer} *************** -->"+"\n"
        string += (0+n_padding)*padding+f"<Grid Name=\"{name_datacontainer}\" GridType=\"Uniform\">"+"\n"
        
        string += insert_geometry(datacontainer, padding, 1+n_padding)+"\n"
        
        for attributematrix in datacontainer:
            if "TupleDimensions" in datacontainer[attributematrix].attrs:
                attributematrix_dims = " ".join(str(i) for i in datacontainer[attributematrix].attrs["TupleDimensions"].astype(int))
                geometry_dims        = " ".join(str(i) for i in datacontainer["_SIMPL_GEOMETRY"+"/"+"DIMENSIONS"][...].astype(int))
                if attributematrix_dims == geometry_dims:
                    name_attributematrix = datacontainer[attributematrix].name.rsplit('/',1)[-1]
                    string += (1+n_padding)*padding+f"<!-- *************** START OF {name_attributematrix} *************** -->"+"\n"
                    for attributearray in datacontainer[attributematrix]:
                        string += insert_attributearray(datacontainer[attributematrix+"/"+attributearray], dream3d_name, padding, 1+n_padding)+"\n"
                    string += (1+n_padding)*padding+f"<!-- *************** END OF {name_attributematrix} *************** -->"+"\n"

        string += (0+n_padding)*padding+"</Grid>"+"\n"
        string += (0+n_padding)*padding+f"<!-- *************** END OF {name_datacontainer} *************** -->"
    
        return string
    
    with h5py.File(dream3d_path, 'r') as dream3d_file, open(xdmf_path, 'w') as xdmf_file:
    
        string = ""
        string += (0+n_padding)*padding+"<?xml version=\"1.0\"?>"+"\n"
        string += (0+n_padding)*padding+"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]>"+"\n"
        string += (0+n_padding)*padding+"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">"+"\n"
        string += (1+n_padding)*padding+padding+"<Domain>"+"\n"
        
        for name_datacontainer in dream3d_file["DataContainers"]:
            datacontainer = dream3d_file[f"DataContainers/{name_datacontainer}"]
            if '_SIMPL_GEOMETRY' in datacontainer.keys() and len(datacontainer['_SIMPL_GEOMETRY'].keys()) > 0:
                string += insert_datacontainer(dream3d_file[f"DataContainers/{name_datacontainer}"], dream3d_name, padding, 2)+"\n"
        
        string += (1+n_padding)*padding+"</Domain>"+"\n"
        string += (0+n_padding)*padding+"</Xdmf>"
        
        xdmf_file.write(string)

#call dream3d file
def call_dream3d(dream3d_path, json_path):

    def quote(string):
        return "\""+string+"\""
    def path(string):
        if platform.system() == "Windows":
            return quote(os.path.realpath(string))
        else:
            return os.path.realpath(string)
    def subpath():
        if platform.system() == "Windows":
            return ""
        else:
            return "/bin"
    def extension():
        if platform.system() == "Windows":
            return ".exe"
        else:
            return ""
    
    cmd = path(dream3d_path+subpath()+"/PipelineRunner"+extension())+" -p "+path(json_path)
    return subprocess.call(shlex.split(cmd))

#replace the file paths in dream3d file because dream3d cannot do relative pathing
def replace_json_paths(json_path,input_path=None,output_path=None,image_path=None):
    
    # Read the DREAM.3D *.JSON file
    with open(json_path, 'r') as json_file:
        data = json.load(json_file)

    # Replace input paths
    if input_path is not None:

        for filter_type in ["DataContainerReader", "ReadCtfData"]:

            filter_ids = get_filter_ids(data, filter_type)
            if len(filter_ids) > 0:
                data[filter_ids[0]]["InputFile"] = input_path
                break
            
    # Replace *.DREAM3D output paths
    if output_path is not None:
        filter_ids = get_filter_ids(data, "DataContainerWriter")
        data[filter_ids[0]]["OutputFile"] = output_path
        
    # Replace *.PNG output paths
    if image_path is not None:
        filter_ids = get_filter_ids(data, "ITKImageWriter")[0]
        data[filter_ids[0]]["FileName"] = image_path
    
    # Write the modified DREAM.3D *.JSON file
    with open(json_path, 'w') as json_file:
        json.dump(data, json_file, indent=4, separators=(",", ": "), sort_keys=True)

# For reasons entirely beyond my comprehension, DREAM.3D does NOT check for attributes from *.DREAM3D files it reads
# Instead, it contains a massive dictionary of expected arrays from the first time the *.DREAM3D file was loaded
# If the previous *.DREAM3D file was changed such that new data exists, DREAM.3D will not read it because it's not in the dictionary
# This function:
#     Reads through the DREAM.3D *.JSON file to find input paths,
#     Scours those *.DREAM3D files for data arrays
#     Adds those data arrays back into the DREAM.3D *.JSON file
################## if there are issues with this function later, it may be because it also needs to scour CellEnsembleData and _SIMPL_GEOMETRY ################
def update_attribute_arrays_expected(path_json):
    
    # variables
    path_CellData = "/DataContainers/ImageDataContainer/CellData"

    with open(path_json,"r") as f:
        data_json = json.load(f)
    
    # find the data container reader ids
    filter_ids = get_filter_ids(data_json, "DataContainerReader")
    
    for filter_id in filter_ids:

        # Open the path to the *.DREAM3D file in the data container reader
        with h5py.File(data_json[filter_id]["InputFile"], 'r') as file_input:
        
            # Create a list of expected data arrays using the attribute arrays from the previous DREAM.3D files
            data_arrays = []
            for name, properties in zip(file_input[path_CellData].keys(), file_input[path_CellData].values()):
                data_arrays += [{
                    "Component Dimensions": file_input[path_CellData][name].attrs["ComponentDimensions"].tolist()       ,
                    "Flag"                : 2                                                                           , # no idea what flag is
                    "Name"                : name                                                                        ,
                    "Object Type"         : file_input[path_CellData][name].attrs["ObjectType"         ].decode('UTF-8'), # maybe should be "ascii"
                    "Path"                : path_CellData                                                               ,
                    "Tuple Dimensions"    : file_input[path_CellData][name].attrs["TupleDimensions"    ].tolist()       ,
                    "Version"             : file_input[path_CellData][name].attrs["DataArrayVersion"   ].tolist()[0]
                }]
                
        for i, data_container in enumerate(data_json[filter_id]["InputFileDataContainerArrayProxy"]["Data Containers"]):
            for j, attribute_matrix in enumerate(data_json[filter_id]["InputFileDataContainerArrayProxy"]["Data Containers"][i]["Attribute Matricies"]):
                if attribute_matrix["Name"] == "CellData":

                    # replace the expected attribute arrays
                    data_json[filter_id]["InputFileDataContainerArrayProxy"]["Data Containers"][i]["Attribute Matricies"][j]["Data Arrays"] = data_arrays
                    
                    # push inputs to dream3d pipeline
                    with open(path_json,"w") as f:
                        json.dump(data_json, f, indent=4, separators=(",", ": "), sort_keys=True)
                    
def update_expected(path_json):
    
    # pull dream3d pipeline
    with open(path_json,"r") as f:
            data_json = json.load(f)
    
    for filter_id in get_filter_ids(data_json, "DataContainerReader"):
    
        file_input = data_json[filter_id]["InputFile"]
        with h5py.File(file_input, 'r') as file_input:
            expected = create_expected(file_input)
        
        data_json[filter_id]["InputFileDataContainerArrayProxy"] = expected["InputFileDataContainerArrayProxy"]
        
    # push inputs to dream3d pipeline
    with open(path_json,"w") as f:
            json.dump(data_json, f, indent=4, separators=(",", ": "), sort_keys=True)
    
    
def create_expected(items, depth=0):
    
    if   depth == 0: # file
        expected = {
            "InputFileDataContainerArrayProxy": create_expected(items["DataContainers"], depth+1)
        }
    elif depth == 1: # root
        expected = {
            "Data Containers": [
                create_expected(items[item], depth+1) for item in [*items.keys()]
            ]
        }
    elif depth == 2: # data containers
        excluded = ["_SIMPL_GEOMETRY"]
        expected = {
            "Attribute Matricies": [
                create_expected(items[item], depth+1) for item in [*items.keys()] if not item in excluded
            ],
            "Flag"               : 2,
            "Name"               : items.name.rsplit('/',1)[1],
            "Type"               : 0
        }
    elif depth == 3: # attribute matricies
        expected = {
            "Data Arrays": [
                create_expected(items[item], depth+1) for item in [*items.keys()]
            ],
            "Flag"       : 2,
            "Name"       : items.name.rsplit('/',1)[1],
            "Type"       : int(items.attrs["AttributeMatrixType"])
        }
    elif depth == 4: # data sets
        expected = {
            "Component Dimensions": [int(dim) for dim in items.attrs["ComponentDimensions"]],
            "Flag"                : 2,
            "Name"                : items.name.rsplit('/',1)[1],
            "Object Type"         : items.attrs["ObjectType"].decode(),
            "Path"                : items.name.rsplit('/',1)[0],
            "Tuple Dimensions"    : [int(dim) for dim in items.attrs["TupleDimensions"    ]],
            "Version"             : int(items.attrs["DataArrayVersion"])
        }
            
    return expected

#you can either get null terminated (as opposed to null padded) strings with s.encode()
#or fixed length (as opposed to vlen or variable length) strings with numpy arrays
#but not both as shown here: https://forum.hdfgroup.org/t/nullpad-nullterm-strings/9107
#the docs explicitly state that null bytes will always cause h5py to throw an error
#dream3d requires non-arbitrary-length strings, so we choose to differ from dream3d
#created strings with null padding
def format_string(s):
    ascii_type = h5py.string_dtype('ascii', len(s)+1)
    b = np.array(s.encode("utf-8"), dtype=ascii_type)
    return b

#map values from local 2d orientations to global 3d orientation
def find_global_values(planes, local_values):
    global_values = []
    for plane in planes:
        global_value = []
        for direction in ["x","y","z"]:
            #this assumes that the orientation in the name matches the orientation in the value
            #if you're having trouble with errors, make sure your names match your actual orientation
            if direction in str(plane).lower():
                global_value.append(local_values[0][0])
                del local_values[0][0]
            else:
                global_value.append(0)
        del local_values[0]
        global_values.append(global_value)
    return np.array(global_values)

#find common value in global 3d orientation given values in local 2d orientations
def find_common_global_value(planes, local_values):
    global_values = find_global_values(planes, local_values)
    global_value = np.zeros((2,3))
    i = 0
    for direction in global_values.transpose():
        direction = direction[np.nonzero(direction)]
        if not all(direction == direction[0]):
            conflicting_planes = np.array(planes)[np.nonzero(direction)]
            direction_names = ["x","y","z"]
            error_message = "Global values do not match: Planes "+str(conflicting_planes).lower()+" disagree in "+direction_names[i]+" direction: "+str(direction)
            return error_message
        global_value[:,i] = direction
        i += 1
    return global_value[0,:]

#determine if multiple list-type objects are equal
def lists_are_equal(lists):
    logical_array = np.array([])
    #use chain logic: if x=y and y=z then x=z 
    for i in range(len(lists)-1):
        #convert list-type objects to numpy arrays and flatten to check equality. will return false if out of order
        logical_array = np.append(logical_array, np.array(lists[i]).flatten() == np.array(lists[i+1]).flatten())
    return all(logical_array)

#determine if multiple hdf5 groups are equal
def groups_are_equal(items):

    #check names
    if not lists_are_equal( [item.name for item in items] ):
        print("names unequal")
        return False

    #check attributes
    if not lists_are_equal( [{**item.attrs} for item in items] ):
        print("attributes unequal")
        return False
    
    # we've already checked the names and attributes, now to check
    # if it's a group or a dataset
    isdataset  = [ isinstance(item, h5py.Dataset) for item in items ]
    
    #check dataset contents
    if all(isdataset):
        if not lists_are_equal( [item[()] for item in items] ):
            print("dataset contents unequal")
            return False
        return True

    #check for mismatched types
    if any(isdataset):
        print("dataset/group types unequal")
        return False
    
    #must be a group, now check the subgroups
    keys = [ item.keys() for item in items ]

    #check group contents
    if not lists_are_equal(keys):
        print("group contents unequal")
        return False
    
    #check subgroups recursively
    for key in keys[0]:
        if not groups_are_equal( [item[key] for item in items] ):
            print('group contents unequal')
            return False
    return True
    
def find_mask_intersected_grains(featureids,mask):
    
    inside_mask  = np.unique(featureids[mask == 1])
    outside_mask = np.unique(featureids[mask == 0])

    return np.intersect1d(inside_mask, outside_mask)

def find_voxel_fraction_mask_intersected_grain(featureid,featureids,mask):
    
    n_voxels        = np.count_nonzero(featureids[featureids == featureid])
    n_voxels_masked = np.count_nonzero(featureids[np.logical_and(featureids == featureid, mask == 1)])

    return n_voxels_masked/n_voxels

def remove_masked_grains(featureids,mask,threshold):

    import time

    time_start = time.time()
    n_grains = len(np.unique(featureids))

    for i, featureid in enumerate(np.unique(featureids), 1):
        fraction = find_voxel_fraction_mask_intersected_grain(featureid,featureids,mask)

        # print ETA
        if i%25 == 0:
            time_elapsed = time.time()-time_start
            time_per_grain = time_elapsed/i
            ETA = time_per_grain*(n_grains-i)
            hrs = int(ETA // 3600)
            mins = int((ETA / 3600 % 1) * 60)
            print(f"Grain {i} of {n_grains}. Estimated time remaining: {hrs} hours, {mins} minutes")

        # delete the grain if too much is outside the mask
        if fraction < threshold:
            featureids[featureids == featureid] = 0

    mask = np.ones_like(featureids).astype(np.uint8)
    mask[featureids == 0] = 0

    return featureids, mask
    
# this attempts to find celldata attribe matricies by matching the attributematrix "TupleDimension"s
# with the _SIMPL_GEOMETRY dimensions
# returns a list of hdf5 paths to the celldata
# Ex: ["/DataContainers/VolumeDataContainer1/CellData","/DataContainers/VolumeDataContainer2/CellData"]
def find_celldata(path_input, path_datacontainer=None):

    paths_celldata = []

    is_celldata = lambda group:                                                  \
        "AttributeMatrixType" in group.attrs.keys()                              \
        and group.attrs["AttributeMatrixType"] == get_attribute_matrix_id("Cell")
    
    with h5py.File(path_input , 'r') as file_input:
        
        # recursively find celldata in every datacontainer if none is specified
        if path_datacontainer is None:
            for name_datacontainer in file_input["/DataContainers"]:
                path_datacontainer = f"/DataContainers/{name_datacontainer}"
                paths_celldata.extend(find_celldata(path_input, path_datacontainer))
                
        # find celldata
        else:
            datacontainer = file_input[path_datacontainer]
            for name_attributematrix in datacontainer:
                attributematrix = datacontainer[name_attributematrix]
                if is_celldata(attributematrix):
                    paths_celldata.append(attributematrix.name)
                        
    return paths_celldata
    
# this attempts to find "featureids" fields by searching in celldata attribute matricies
# for data arrays that contain signed 32 bit integers and whose names contain the substring "featureids"
# by default, it finds celldata using find_celldata, celldata can be specified as well
# celldata should be a list of pathnames.
# Ex: ["/DataContainers/VolumeDataContainer1/CellData","/DataContainers/VolumeDataContainer2/CellData"]
# returns a list of hdf5 paths to the featureids
# Ex: ["/DataContainers/VolumeDataContainer1/CellData/FeatureIds","/DataContainers/VolumeDataContainer2/CellData/FeatureIds"]
def find_featureids(path_input, paths_celldata=None, path_datacontainer=None):
    
    featureids = []
    
    if paths_celldata is None and path_datacontainer is None:
        paths_celldata = find_celldata(path_input)
    elif paths_celldata is None:
        paths_celldata = find_celldata(path_input, path_datacontainer=path_datacontainer)
    elif isinstance(paths_celldata, str):
        paths_celldata = [paths_celldata]
    elif not type(paths_celldata) in ["<class 'list'>", "<class 'tuple'>", "<class 'numpy.ndarray'>"]:
        raise('Celldata must be string or list-like')
    # regardless the value of datacontainer, if celldata is specified, it is used

    is_featureid = lambda dataarray:                                       \
        "ObjectType" in dataarray.attrs                                    \
        and dataarray.attrs["ObjectType"].decode() == "DataArray<int32_t>" \
        and "featureids" in dataarray.name.rsplit("/",1)[-1].lower()
    
    with h5py.File(path_input , 'r') as file_input:
    
        for path_attributematrix in paths_celldata:
            attributematrix = file_input[path_attributematrix]
            
            for name_dataarray in attributematrix:
                dataarray = attributematrix[name_dataarray]

                if is_featureid(dataarray):
                    featureids.append(dataarray.name)
                    
    return featureids
    
# this attempts to find cellfeaturedata linked to featureids fields by finding
# 1) the number of unique featureids in the featureids field
# 2) the length of the candidate cellfeaturedata "TupleDimensions"
# if these dimensions match, then the two are at least compatible if not intentionally linked
# by default, if finds featureids using find_featureids, but featureids can be specified as well
# featureids should be a list of hdf5 path names.
# Ex: ["/DataContainers/VolumeDataContainer1/CellData/FeatureIds","/DataContainers/VolumeDataContainer2/CellData/FeatureIds"]
# returns a dictionary mapping of each featureid with all the potential cellfeaturedata matches
# Ex: {"[...]/CellData/FeatureIds":["[...]/CellFeatureData1","[...]/CellFeatureData2"], [...]}
def find_map_featureids_cellfeaturedata(path_input, paths_featureids=None, paths_celldata=None, path_datacontainer=None):
    
    if paths_featureids is None and paths_celldata is None and path_datacontainer is None:
        paths_featureids = find_featureids(path_input)
    elif paths_featureids is None and not paths_celldata is None:
        paths_featureids = find_featureids(path_input, paths_celldata=paths_celldata)
    elif paths_featureids is None:
        paths_featureids = find_featureids(path_input, path_datacontainer=path_datacontainer)
    # regardless the value of datacontainer, if featureids is specified, it is used

    is_cellfeaturedata   = lambda group:                                                \
        "AttributeMatrixType" in group.attrs                                            \
        and group.attrs["AttributeMatrixType"] == get_attribute_matrix_id("CellFeature")
    has_equal_dimensions = lambda group, dims_reference:                      \
        "TupleDimensions" in group.attrs                                      \
         and all(group.attrs["TupleDimensions"].astype(int) == dims_reference)
        
    map_featureids_cellfeaturedata = {path_featureids:[] for path_featureids in paths_featureids}
    
    with h5py.File(path_input , 'r') as file_input:
    
        for path_featureids in map_featureids_cellfeaturedata.keys():
            
            # find the number of unique featureids
            featureids           = np.unique(file_input[path_featureids][...])
            dims_cellfeaturedata = np.array([(featureids.size+1, featureids.size)[0 in featureids]]).astype(int)
            
            # find the datacontainer housing the featureids
            path_datacontainer = path_featureids.rsplit("/",2)[0]
            datacontainer      = file_input[path_datacontainer]
            
            # find the any cellfeaturedata matricies that could be associated with the featureids
            for name_attributematrix in datacontainer:
                attributematrix = datacontainer[name_attributematrix]

                # associate the attribute matrix with the featureids data array if:
                # it is cellfeaturedata matrix and it has the same number of dimensions as the featureids array has unique features
                if is_cellfeaturedata(attributematrix) and has_equal_dimensions(attributematrix, dims_cellfeaturedata):
                    map_featureids_cellfeaturedata[path_featureids].append(attributematrix.name)
                    
    return map_featureids_cellfeaturedata

def crop(path_input, path_output, points):

    # accept only list-like objects for points
    types_accepted = ["<class 'list'>", "<class 'tuple'>", "<class 'numpy.ndarray'>"]
    if not str(type(points)) in types_accepted:
        raise TypeError(f"Points were given as {type(points)}. Accepted types are {', '.join(types_accepted)}")
    
    # cast points as a numpy array
    points = np.asarray(points)

    # catch selection errors
    if not points.shape == (2,3):
        raise ValueError(f"Points must be of shape [[x1,y1,z1],[x2,y2,z2]]")
    if any( points[1,:]-points[0,:] <  0 ):
        points_inverted = np.broadcast_to(points[1,:]-points[0,:] < 0, points.shape)
        points[points_inverted] = points[points_inverted].reshape(2,-1)[::-1,:].flatten()
        print(f"Axis was inverted in points. Swapped automatically")
    if any( points[1,:]-points[0,:] == 0):
        raise ValueError("Each axis in points must span a length greater than zero")
    
    def is_mutible(file_input, path_attributematrix):
        
        if '_SIMPL_GEOMETRY' in path_attributematrix:
            return True
        
        if file_input[path_attributematrix].attrs["AttributeMatrixType"] in [get_attribute_matrix_id(matrix_type) for matrix_type in ["CellFeature", "Cell"]]:
            return True
        
        return False
    
    def insert_celldata(file_input, path_output, path_celldata, points_bounded):

        map_featureids_cellfeaturedata = find_map_featureids_cellfeaturedata(path_input, paths_celldata=path_celldata)
        
        for name_dataarray in file_input[path_celldata]:
            path_dataarray = f"{path_celldata}/{name_dataarray}"

            # grab the selection from the array
            x, y, z = points_bounded.T
            data = file_input[path_dataarray][ z[0]:z[1], y[0]:y[1], x[0]:x[1], ... ]
            
            # find the array data type
            type_data = file_input[path_dataarray].attrs["ObjectType"].decode('UTF-8')
            
            # if it's a featureid, renumber and insert its associated cellfeaturedata
            if path_dataarray in [*map_featureids_cellfeaturedata.keys()]:

                ### re-order featureids
                # find the featureids that remain after the crop
                # sort required so that the values in cellfeaturedata matches
                # np.unique pre-sorts values, no not needed here explicitly
                featureids_old = np.unique(data)
                # the featureids must have contiguous unique ids 
                # such that they can be used to reference cellfeaturedata dataarrays
                # create the new values, but DON'T start at zero unless it already exists in features_old
                # this is because 0 has a special significance as an error value throughout dream.3d
                featureids_new = np.arange(0, featureids_old.size)+(1,0)[0 in data]
                # cycle through each old and new value as they're already sorted and of the same size
                # replace any value of data that matches the old value and replace it with the new value
                for featureid_old, featureid_new in zip(featureids_old, featureids_new): data[data==featureid_old]=featureid_new

                # insert cellfeaturedata
                insert_cellfeaturedata(file_input, path_output, path_dataarray, map_featureids_cellfeaturedata, featureids_old, featureids_new)
            
            # insert the array
            insert_attribute_array(
                path_output   ,
                path_celldata ,
                name_dataarray,
                data          ,
                type_data     ,
                "Cell"
            )

            # ### DEBUG ####
            # # test that the proposed and assigned featureids align
            # if path_celldata+"/"+name_dataarray in [key for key in map_featureids_cellfeaturedata.keys()]:
            #     with h5py.File(path_output, 'r') as f:
            #         features_actual = np.unique(f[f"{path_celldata}/{name_dataarray}"][...])
            #     if not all(featureids_new == features_actual):
            #         raise ValueError("feature ids do not align")

    def find_linked_datasets(file_input, path_attribute_matrix):

        map_linked_datasets = {}

        attribute_matrix = file_input[path_attribute_matrix]
        for name_dataarray in attribute_matrix:
            dataarray = attribute_matrix[name_dataarray]

            if "Linked NumNeighbors Dataset" in dataarray.attrs:
                map_ = dataarray.attrs["Linked NumNeighbors Dataset"].decode('UTF-8')
                if not map_ in map_linked_datasets: map_linked_datasets[map_]     = [name_dataarray]
                else                              : map_linked_datasets[map_].append(name_dataarray)

        return map_linked_datasets

    def find_neighborlist(file_input, path_cellfeaturedata, linked_datasets, name_numneighbors, path_featureids):

        linked_dataset_maps = [name_numneighbors]+[item for item in linked_datasets if not item == name_numneighbors]
        linked_dataset_vals = linked_datasets[name_numneighbors]

        names_neighborlist = []

        is_neighborlist = lambda dataset, featureids_unique: \
            dataset[...].dtype == np.int32 \
            and set(np.unique(dataset[...])).issubset(featureids_unique)
        
        featureids_unique = np.unique(file_input[path_featureids][...])

        for i, linked_dataset_map in enumerate(linked_dataset_maps, 1):

            linked_dataset_vals = linked_datasets[linked_dataset_map]

            for name_linked_dataset_val in linked_dataset_vals:
                path_linked_dataset_val = f"{path_cellfeaturedata}/{name_linked_dataset_val}"
                linked_dataset_val      = file_input[path_linked_dataset_val]

                if is_neighborlist(linked_dataset_val, featureids_unique): 
                    names_neighborlist.append(name_linked_dataset_val)

            if len(names_neighborlist) >= 1:
                break

            if len(linked_dataset_maps) > 1 and i < len(linked_dataset_maps):
                print(f'No neighborlists found in linked list. Checking next linked list')

        return(names_neighborlist)

    
    def insert_cellfeaturedata(file_input, path_output, path_featureids, map_featureids_cellfeaturedata, featureids_old, featureids_new):

        map_featureids_old_new = { old:new for old, new in zip(featureids_old, featureids_new) }
        
        # export CellFeatureData group(s)
        for path_cellfeaturedata in map_featureids_cellfeaturedata[path_featureids]:

            featureids_remaining = sorted((np.append([0],np.unique(featureids_old)),np.unique(featureids_old))[0 in featureids_old])

            # find the cellfeature arrays that are part of a linked list
            linked_datasets = find_linked_datasets(file_input, path_cellfeaturedata)
            linked_dataset_maps = [*linked_datasets.keys()]
            linked_dataset_vals = sum([*linked_datasets.values()],[])
            not_linked_datasets = list(set([*file_input[path_cellfeaturedata].keys()]).difference([*linked_dataset_maps, *linked_dataset_vals]))    
            
            # insert all the cellfeature arrays that are not part of a linked list
            for name_dataarray in not_linked_datasets:
                path_dataarray = f"{path_cellfeaturedata}/{name_dataarray}"
                dataarray      = file_input[path_dataarray]
            
                data      = dataarray[featureids_remaining]
                type_data = dataarray.attrs["ObjectType"].decode('UTF-8')
                
                insert_attribute_array(
                    path_output                ,
                    path_cellfeaturedata       ,
                    name_dataarray             ,
                    data                       ,
                    type_data                  ,
                    "CellFeature"
                )

            # restructure the linked lists
            for name_numneighbors in linked_datasets:
                path_numneighbors   = f"{path_cellfeaturedata}/{name_numneighbors}"
                numneighbors        = file_input[path_numneighbors]

                # find neighborlist array
                names_neighborlist = find_neighborlist(file_input, path_cellfeaturedata, linked_datasets, name_numneighbors, path_featureids)
                if len(names_neighborlist) == 0:
                    print(f"Skipping incomplete linked list. No possible neighborlists found for numneighbor: {name_numneighbors}")
                    continue
                if len(names_neighborlist) >  1:
                    print(f"Skipping ambiguous linked list. Multiple possible neighborlists found for numneighbor: {name_numneighbors}:{names_neighborlist}")
                    continue
                name_neighborlist = names_neighborlist[0]
                path_neighborlist = f"{path_cellfeaturedata}/{name_neighborlist}"
                neighborlist      = file_input[path_neighborlist]
                
                # not_neighborlists do not provide enough information to be processed on their own
                # so they have to be processed using the indecies of neighborlist
                names_not_neighborlist = list(set([*linked_datasets[name_numneighbors]]).difference(['NeighborList']))
                

                # initialize new arrays
                numneighbors_new = []
                neighborlist_new = []
                vals_not_neighborlist = [ [] for name in names_not_neighborlist ]

                # for each feature (grain), find the number of neighbors and what those neighbors are
                for index_numneighbors in range(0, numneighbors[...].size):
                    
                    # initialize the number of neighbors for this feature (grain) as zero
                    numneighbors_new_i = 0

                    # if this feature is not in the new domain, don't add it to numneighbors
                    if not index_numneighbors in featureids_remaining:
                        continue
                    
                    # find the index range for the current grain's neighbors in neighborlist
                    # by finding the sum of all previous grain's number of neighbors
                    start = sum(numneighbors[...].flatten()[:index_numneighbors])
                    stop  = start+numneighbors[...].flatten()[index_numneighbors]
                    for index_val in range(start, stop):
                        
                        # if this feature is not in the new domain, don't add it to neighborlist
                        # and don't increase the number of neighbors for the current feature (grain)
                        if not neighborlist[index_val] in featureids_remaining:
                            continue

                        numneighbors_new_i += 1
                        neighborlist_new.append( map_featureids_old_new[neighborlist[index_val]] )
                        for name, val in zip(names_not_neighborlist, vals_not_neighborlist):
                            val.append( file_input[f"{path_cellfeaturedata}/{name}"][index_val] )

                    numneighbors_new.append(numneighbors_new_i)

                # insert numneighbors
                data      = np.asarray(numneighbors_new).reshape((-1,1)).astype(numneighbors[...].dtype)
                type_data = numneighbors.attrs["ObjectType"].decode('UTF-8')
                insert_attribute_array(
                    path_output                ,
                    path_cellfeaturedata       ,
                    name_numneighbors          ,
                    data                       ,
                    type_data                  ,
                    "CellFeature"
                )

                # if this neighborlist was borrowed from another linked list, don't export it
                if name_neighborlist in linked_datasets[name_numneighbors]:

                    # insert neighborlist
                    data      = np.asarray(neighborlist_new).astype(neighborlist[...].dtype)
                    type_data = neighborlist.attrs["ObjectType"].decode('UTF-8')
                    linked_numneighbors_dataset = neighborlist.attrs["Linked NumNeighbors Dataset"].decode('UTF-8')
                    tuple_dimensions            = neighborlist.attrs["TupleDimensions"]
                    insert_attribute_array(
                        path_output                ,
                        path_cellfeaturedata       ,
                        name_neighborlist          ,
                        data                       ,
                        type_data                  ,
                        "CellFeature"              ,
                        linked_numneighbors_dataset,
                        tuple_dimensions
                    )

                # insert other linked lists
                for name, val in zip(names_not_neighborlist, vals_not_neighborlist):
                    data_old = file_input[f"{path_cellfeaturedata}/{name}"]
                    data      = np.asarray(val).astype(data_old[...].dtype)
                    type_data = data_old.attrs["ObjectType"].decode('UTF-8')
                    linked_numneighbors_dataset = data_old.attrs["Linked NumNeighbors Dataset"].decode('UTF-8')
                    tuple_dimensions            = data_old.attrs["TupleDimensions"]
                    insert_attribute_array(
                        path_output                ,
                        path_cellfeaturedata       ,
                        name                       ,
                        data                       ,
                        type_data                  ,
                        "CellFeature"              ,
                        linked_numneighbors_dataset,
                        tuple_dimensions
                    )

            # if the cellfeaturedata path was linked to any other featureids, unlink it
            for path_featureids_i, paths_cellfeaturedata_i in zip(map_featureids_cellfeaturedata.keys(), map_featureids_cellfeaturedata.values()):
                if not path_featureids == path_featureids_i and path_cellfeaturedata in paths_cellfeaturedata_i:
                    paths_cellfeaturedata_i.remove(path_cellfeaturedata)


        
    # create a new empty file
    create_file_dream3d(path_output)

    with h5py.File(path_input , 'r') as file_input:

        # copy select data from each datacontainer
        for name_datacontainer in file_input["/DataContainers"]:
            path_datacontainer = f"/DataContainers/{name_datacontainer}"
            
            # find the original geometry
            path_geometry  = f"{path_datacontainer}/_SIMPL_GEOMETRY"
            dims_domain    = file_input[f"{path_geometry}/DIMENSIONS"][...]
            origin_domain  = file_input[f"{path_geometry}/ORIGIN"    ][...]
            spacing_domain = file_input[f"{path_geometry}/SPACING"]

            # bound the selection to find the new geometry
            # origin refers to the physical coordinates, but we're referring to the 
            # indeces which should always start at zero
            # points_min     = np.max(np.vstack((origin_domain            , points[0,:])),axis=0)
            # points_max     = np.min(np.vstack((origin_domain+dims_domain, points[1,:])),axis=0)
            points_min     = np.max(np.vstack((np.zeros(dims_domain.shape), points[0,:])),axis=0)
            points_max     = np.min(np.vstack((dims_domain                , points[1,:])),axis=0)
            points_bounded = np.vstack((points_min, points_max)).astype(int)
            dims_subdomain = points_bounded[1,:]-points_bounded[0,:]
            geometry = {
                "dims"       : dims_subdomain,
                "origin"     : origin_domain ,
                "size_voxels": spacing_domain
            }
            
            # if any axis in the selection became zero length or inverted due to boundary restrictions
            # then the selection was not in the domain 
            if any( dims_subdomain <= 0 ):
                print(f"Selection: {points} is not in domain: {np.append(origin_domain.reshape(1,-1), dims_domain.reshape(1,-1), axis=0)} for datacontainer: {name_datacontainer}")
                continue
            
            # insert the geometry group
            insert_geometry(
                path_output,
                geometry,
                path_geometry
            )
            
            # insert the CellData group(s)
            # and any CellFeatureData groups associated with a FeatureId
            for path_celldata in find_celldata(path_input, path_datacontainer=path_datacontainer):
                insert_celldata(file_input, path_output, path_celldata, points_bounded)
            
            # export all other groups if they are not mutible types
            for name_attributematrix in file_input[path_datacontainer]:
                path_attributematrix = f"{path_datacontainer}/{name_attributematrix}"
                if not is_mutible(file_input, path_attributematrix):
                    with h5py.File(path_output, 'a') as file_output:
                        file_output.copy(file_input[path_attributematrix], path_attributematrix)
    
    # make the xdmf file to visualize in paraview
    make_xdmf(path_output)

    #Tested working

    #def unique2D_subarray(a):
    #    dtype1 = np.dtype((np.void, a.dtype.itemsize * np.prod(a.shape[1:])))
    #    b = np.ascontiguousarray(a.reshape(a.shape[0],-1)).view(dtype1)
    #    return a[np.unique(b, return_index=1)[1]]

    #with h5py.File(path_output , 'r') as file_output:

    #    path_celldata = "/DataContainers/ImageDataContainer/CellData"
    #    path_cellfeaturedata = "/DataContainers/ImageDataContainer/CellFeatureData"
    #    
    #    path_featureids = path_celldata+"/"+"FeatureIds"
    #    path_element_eulerangles = path_celldata+"/"+"AvgEulerAngles"
    #    path_feature_eulerangles = path_cellfeaturedata+"/"+"AvgEulerAngles"

    #    featureids = file_output[path_featureids][...]
    #    element_eulerangles = file_output[path_element_eulerangles][...]
    #    feature_eulerangles = file_output[path_feature_eulerangles][...]
    #    
    #    featureids_unique = np.unique(featureids)
    #    
    #    for i in range(featureids_unique.size):
    #    
    #        element_eulerangle_i = element_eulerangles[element_eulerangles.shape[-1]*[True]*featureids==i].reshape((-1,element_eulerangles.shape[-1]))
    #        element_eulerangle_i_unique = unique2D_subarray(element_eulerangle_i)
    #        
    #        feature_eulerangle_i = feature_eulerangles[i]
    #        
    #        if element_eulerangle_i_unique.shape[0] > 1:
    #            print("error! more than one eulerangle at {i}. maybe feature average not selected?")
    #            continue
    #        if not all(element_eulerangle_i_unique[0] == feature_eulerangle_i):
    #            print(f"failed! not equal at {i}")
    #            print("   "+"element: "+f"{element_eulerangle_i_unique[0]}")
    #            print("   "+"feature: "+f"{feature_eulerangle_i}")