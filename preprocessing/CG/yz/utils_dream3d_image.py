import h5py
import numpy as np
from PIL import Image
import cv2

import os, sys
new_path = os.path.dirname(__file__)
if new_path not in sys.path:
    sys.path.append(new_path)

# you can either get null terminated (as opposed to null padded) strings with s.encode()
# or fixed length (as opposed to vlen or variable length) strings with numpy arrays
# but not both as shown here: https://forum.hdfgroup.org/t/nullpad-nullterm-strings/9107
# the docs explicitly state that null bytes will always cause h5py to throw an error
# dream3d requires non-arbitrary-length strings, so we choose to differ from dream3d
# created strings with null padding
def format_string(s):
    ascii_type = h5py.string_dtype('ascii', len(s)+1)
    return np.array(s.encode("utf-8"), dtype=ascii_type)

def insert_mask(path_file, path_dataset, data):

    ## Find dimension and component information
    dims                  = np.asarray(data.shape[:-1][::-1]).astype(np.uint64)
    components            = np.asarray(data.shape[-1]).reshape((1,)).astype(np.uint64)
    tuple_axis_dimensions = format_string(','.join(['='.join([name,str(value)]) for name,value in zip(['x','y','z'],dims.tolist())]))

    ## Create attributes for the attribute matrix
    group_attributes = {
        "AttributeMatrixType": np.uint32(3),
        "TupleDimensions"    : dims
    }

    ## Create attributes common to all datasets
    dataset_attributes = {
        "ComponentDimensions"  : components,
        "DataArrayVersion"     : np.asarray(2).reshape((1,)).astype(np.int32),
        "ObjectType"           : format_string("DataArray<bool>"),
        "Tuple Axis Dimensions": tuple_axis_dimensions,
        "TupleDimensions"      : dims
    }

    ## Write out to file
    path_group, name_dataset = path_dataset.rsplit('/',1)
    with h5py.File(path_file, 'a') as file_output:

        # create the group and its required attributes if it does not exist
        group = file_output.require_group(path_group)
        for key, val in group_attributes.items():
            if not key in group.attrs:
                group.attrs[key] = val

        # export attribute array
        dataset = group.create_dataset(name_dataset, data=data)
        for key, val in dataset_attributes.items():
            if not key in dataset.attrs:
                dataset.attrs[key] = val

def mask_to_image(path_file_input, path_dataset_input, threshold=None):

    # import data from dream3d to numpy array
    with h5py.File(path_file_input, 'r') as file_input:
        dataset_input = file_input[path_dataset_input][...].squeeze()

    # threshold data and create boolean map
    if threshold is not None:
        dataset_input[dataset_input >= threshold] = 1
        dataset_input[dataset_input <  threshold] = 0

    # scale to image data range [0, 255]
    dataset_input = dataset_input/dataset_input.max()*255

    # convert numpy array to data
    image = Image.fromarray(dataset_input)
    image = image.convert('RGB')

    return image

def clean_mask(
    image_input           ,
    gauss_blur       =  10, # Gaussian blur level
    contrast_alpha   =   2, # Contrast multiplier
    contrast_beta    =  20, # Contrast shift
    binary_threshold = 200, # Binarization threshold (0-255)
    show_images      = False
):
    
    array_input = np.array(image_input)

    # input image data
    array_input = cv2.cvtColor(array_input, cv2.COLOR_RGB2BGR)
    if show_images: cv2.imshow('input',array_input); cv2.waitKey(0)

    # blur
    array_blurred = cv2.blur(array_input,(gauss_blur,gauss_blur))
    if show_images: cv2.imshow('blurred',array_blurred); cv2.waitKey(0)

    # contrast
    array_contrasted = np.clip(contrast_alpha*array_blurred.astype(np.float64)+contrast_beta, 0, 255).astype(np.uint8)
    if show_images: cv2.imshow('contrasted',array_contrasted); cv2.waitKey(0)

    # binarize
    array_binarized = (array_contrasted > binary_threshold).astype(np.uint8)*255
    if show_images: cv2.imshow('binarized',array_binarized); cv2.waitKey(0)

    return Image.fromarray(array_binarized)

def image_to_mask(image_input,path_file_output,path_dataset_output):

    array_input = np.array(image_input)[:,:,0]
    array_input = array_input.reshape(([1]+list(array_input.shape)+[1])).astype(np.uint8)

    #write in the dataset
    insert_mask(path_file_output, path_dataset_output, array_input)



