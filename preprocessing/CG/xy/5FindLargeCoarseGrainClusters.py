# Find clusters of large grains by creating a mask out of the indexed voxels,
# using a gaussian blur to smooth out the inter-particle grain boundaries,
# using a contrast filter further accuntuate the difference between the
# CG and UFG regions,
# and finally binarizing the image using a threshold to obtain the CG/UFG mask.
# This process is more art than science and will require a great deal of time
# tweaking the image cleaning variables.
# When the mask is created, it is then inserted into the new *.dream3d file.

import utils_dream3d_image
import shutil

# data import variables
path_file_input = "pipeline_output/4-Separated_domains.dream3d"
path_file_output = "pipeline_output/5-Separated_domains.dream3d"
path_dataset_input = "DataContainers/ImageDataContainer/CellData/Mask_Particle_Domain"
path_dataset_output = "DataContainers/ImageDataContainer/CellData/Mask_Particle_Domain_Cleaned"

# image creation variables
threshold        = 1     # binary threshold. (1 if val >= threshold else 0)
# image cleaning variables
gauss_blur       =  10   # Gaussian blur level           
contrast_alpha   =   2   # Contrast multiplier           
contrast_beta    = -20   # Contrast shift                
binary_threshold = 200   # Binarization threshold (0-255)
show_images      = True # Pop up images of each stage

# import data from dream3d to image
image = utils_dream3d_image.mask_to_image(
    path_file_input,
    path_dataset_input,
    threshold=threshold
)

# clean images
mask = utils_dream3d_image.clean_mask(
    image,
    gauss_blur=gauss_blur,
    contrast_alpha=contrast_alpha,
    contrast_beta=contrast_beta,
    binary_threshold=binary_threshold,
    show_images=show_images
)

# export data from image to dream3d
shutil.copy2(path_file_input, path_file_output)
utils_dream3d_image.image_to_mask(mask, path_file_output, path_dataset_output)
