import utils_dream3d_image
import shutil

# data import variables
path_file_input = "pipeline_output/4-Separated_domains.dream3d"
path_file_output = "pipeline_output/5-Separated_domains.dream3d"
path_dataset_input = "DataContainers/ImageDataContainer/CellData/Mask_Particle_Domain"
path_dataset_output = "DataContainers/ImageDataContainer/CellData/Mask_Particle_Domain_Cleaned"
threshold=None

# image cleaning variables
gauss_blur=10
contrast_alpha=2#1.75
contrast_beta=-20
binary_threshold=200
show_images=False

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
