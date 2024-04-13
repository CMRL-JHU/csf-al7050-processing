# sample windows of the 2d EBSD data for each orthoganal plane,
# and sample slices of each orthogonal direction for the 3d SEVMS,
# then use those samples create histograms of the UFG area fractions.
# we can then use these samples to compare the area fractions of our
# SEVMS with the area fractions of the EBSD data.

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as plt_norm
from matplotlib.cm import ScalarMappable as plt_sm
import seaborn as sns
import pandas as pd
import numpy as np
import h5py
from scipy.ndimage import affine_transform

import os, sys
new_path = os.path.dirname(__file__)
if new_path not in sys.path:
    sys.path.append(new_path)

from utils_dict import create_dict
import utils_python

path_config = {
    "project name": "13-AreaFractions",
    "input": {
        "EBSD"            :{
            "files":{
                "XY": {
                    "path"      : "pipeline_input/6.3-xy-cleaned_grains.dream3d",
                    "datasets"  : {
                        "featureids": "/DataContainers/ImageDataContainer/CellData/FeatureIds"
                    }
                },
                "ZX": {
                    "path"      : "pipeline_input/6.3-yz-cleaned_grains.dream3d",
                    "datasets"  : {
                        "featureids": "/DataContainers/ImageDataContainer/CellData/FeatureIds"
                    }
                },
                "YZ": {
                    "path"      : "pipeline_input/6.3-yz-cleaned_grains.dream3d",
                    "datasets"  : {
                        "featureids": "/DataContainers/ImageDataContainer/CellData/FeatureIds"
                    }
                }
            }
        },
        "GAN"        :{
            "file": {
                "path"      : "pipeline_input/2-slicegan_data.dream3d",
                "datasets"  : {
                    "featureids": "/DataContainers/ImageDataContainer/CellData/Mask"
                }
            }
        },
        "Improved"        :{
            "file": {
                "path"      : "pipeline_output/12-synthetic_grains_applied_improved_mask.dream3d",
                "datasets"  : {
                    "featureids": "/DataContainers/ImageDataContainer/CellData/mask_updated"
                }
            }
        }
    },
    "operations": {
        # "plot_slices"        : {
        #     "n_samples"  : 5,
        #     "path_output": "./pipeline_output/[self: project name]/[name source]/[name plane]@[name axis]=[coordinate axis].png"
        # },
        "find_area_fractions": {
            "n_samples"  : 300,
            "path_output": None
        }
    },
    "output": {
        "path_output" : "./pipeline_output/13-[type plot]_[attribute].eps",
        "map_names"   : {
            "find_area_fractions": "Area Fraction",
        }
    }
}

def Dream3dDataContainer(name, config):

    config = create_dict(config)
    data   = create_dict({})

    # 2d type
    if 'files' in config.get():
        for plane in config.get('files'):
            with h5py.File(config.get(f'files/{plane}/path'), 'r') as f:
                for dataset in config.get(f'files/{plane}/datasets'):
                    # reshape to [y,x,components]
                    array = f[config.get(f'files/{plane}/datasets/{dataset}')][...]
                    array = array.reshape(array.shape[1:])
                    data.set(
                        path      = f'{dataset}/{plane}',
                        set_value = array
                    )
        return Dream3dDataContainer2d(name, config, data)

    # 3d type
    elif 'file'  in config.get():
        with h5py.File(config.get(f'file/path'), 'r') as f:
            for dataset in config.get(f'file/datasets'):
                data.set(
                    path      = f'{dataset}',
                    set_value = f[config.get(f'file/datasets/{dataset}')][...]
                )
        return Dream3dDataContainer3d(name, config, data)

    else:
        raise ValueError(
            f"no recognized file structure to import in:\n{utils_python.pretty_file_dict(config.get())}"
        )

class Dream3dDataContainer2d():
    
    def __init__(self, name, config, data):

        self.name   = name
        self.config = config
        self.data   = data

        self.SliceSampler = SliceSampler2d

        ##################################
        print(f"{self.name = }")
        for plane in self.data.get('featureids'):
            featureids = self.data.get(f"featureids/{plane}")
            area_fraction = np.count_nonzero(featureids == 0)/featureids.size
            print(f"   {plane}: {area_fraction = }")
        ##################################

    def __str__(self):
        
        string =  "name:\n"
        string += f"   {self.name}"
        string += "type:\n"
        string += "   2d\n"
        string += "config:\n"
        string += utils_python.pretty_file_dict(self.config.get(), i=2)

        return string

    def find_area_fractions(self, n_samples=5, sample_shape=None, path_output=None):

        slice_sampler  = self.SliceSampler(data=self.data.get('featureids'), n=n_samples, sample_shape=sample_shape)
        area_fractions = {}

        for data_slice, name_plane, coordinates_plane in slice_sampler:

            # print_progress
            # print(f"{self.name}@find_area_fractions: {name_plane}={coordinates_plane}")

            # find the area fraction in this plane (zero values to all values)
            if not name_plane in area_fractions: area_fractions[name_plane] = []
            area_fractions[name_plane].append( np.count_nonzero(data_slice == 0)/data_slice.size )

        # save output
        if path_output is not None:

            for name_plane in area_fractions:

                # fill name placeholders
                replacements = {
                    '[name source]': self.name,
                    '[name plane]' : name_plane
                }
                path_output_finalized = path_output
                for placeholder, replacement in replacements.items():
                    path_output_finalized = path_output_finalized.replace(placeholder, replacement)

                # create directory tree and file
                os.makedirs(os.path.split(path_output_finalized)[0], exist_ok=True)
                with open(path_output_finalized, 'w') as f:
                    for area_fraction in area_fractions[name_plane]:
                        f.write(str(area_fraction))

        return area_fractions

class Dream3dDataContainer3d():

    def __init__(self, name, config, data):

        self.name   = name
        self.config = config
        self.data   = data

        self.SliceSampler = SliceSampler3d

    def __str__(self):

        string =  "name:\n"
        string += f"   {self.name}"
        string += "type:\n"
        string += "   3d\n"
        string += "config:\n"
        string += utils_python.pretty_file_dict(self.config.get(), i=2)

        return string

    def plot_slices(self, path_output=".", n_samples=5):

        slice_plotter = SlicePlotter     (data=self.data.get('featureids'))
        slice_sampler = self.SliceSampler(data=self.data.get('featureids'), n=n_samples)

        for data_slice, name_plane, coordinate_axis in slice_sampler:

            # print progress
            name_axis  = 'XYZ'.translate({ord(i): None for i in name_plane.upper()})
            # print(f"{self.name}@plot_slices: {name_plane}@{name_axis}={coordinate_axis}")

            # fill name placeholders
            replacements = {
                '[name source]'    : self.name           ,
                '[name plane]'     : name_plane          ,
                '[name axis]'      : name_axis           ,
                '[coordinate axis]': str(coordinate_axis)
            }
            path_output_finalized = path_output
            for placeholder, replacement in replacements.items():
                path_output_finalized = path_output_finalized.replace(placeholder, replacement)

            # plot image to file
            slice_plotter.plot_slice(data_slice, path_output_finalized)

    def find_area_fractions(self, n_samples=5, path_output=None):
        
        slice_sampler  = self.SliceSampler(data=self.data.get('featureids'), n=n_samples)
        area_fractions = {}

        for data_slice, name_plane, coordinate_axis in slice_sampler:

            # print progress
            name_axis = 'XYZ'.translate({ord(i): None for i in name_plane.upper()})
            # print(f"{self.name}@find_area_fractions: {name_plane}@{name_axis}={coordinate_axis}")

            # find the area fraction in this plane (zero values to all values)
            if not name_plane in area_fractions: area_fractions[name_plane] = []
            area_fractions[name_plane].append( np.count_nonzero(data_slice == 0)/data_slice.size )

        # save output
        if path_output is not None:

            for name_plane in area_fractions:

                # fill name placeholders
                replacements = {
                    '[name source]': self.name,
                    '[name plane]' : name_plane
                }
                path_output_finalized = path_output
                for placeholder, replacement in replacements.items():
                    path_output_finalized = path_output_finalized.replace(placeholder, replacement)

                # create directory tree and file
                os.makedirs(os.path.split(path_output_finalized)[0], exist_ok=True)
                with open(path_output_finalized, 'w') as f:
                    for area_fraction in area_fractions[name_plane]:
                        f.write(str(area_fraction))

        return area_fractions

class SlicePlotter():

    def __init__(self, data):

        if not len(data.shape) == 4: raise ValueError("data must be formatted as an array of [z,y,x,components]")

        self.n_colors     = np.unique(data.reshape(-1, data.shape[-1]), axis=0).shape[0]
        self.n_components = data.shape[3]
        
        self.init_vars_color(self.n_components)
            
    def init_vars_color(self, n_color_values):

        match n_color_values:
            case 1:
                self.type_color = 'scalar'
                self.cmap       = plt.get_cmap(name='jet', lut=self.n_colors)
                self.norm       = plt_norm(vmin=0, vmax=self.n_colors)
                self.sm         = plt_sm(cmap=self.cmap, norm=self.norm)
                self.plot_func  = self.plot_scalar
            case 3:
                self.type_color = 'color'
                raise ValueError('color data not implemented')
            case _:
                raise ValueError('unknown color type')

    def plot_slice(self, data, path_output):

        self.plot_func(data, path_output)
        
    def plot_scalar(self, data, path_output):

        # plot the slice
        fig, ax = plt.subplots()
        ax.imshow(self.cmap(data.reshape(data.shape[:-1])), origin='lower')
        fig.colorbar(self.sm, ax=ax)

        # save figure
        os.makedirs(os.path.split(path_output)[0], exist_ok=True)
        fig.savefig(path_output)

        # close the figure to free up memory
        plt.close(fig)

class SliceSampler3d():

    def __init__(self, data, n=5):

        self.data               = data
        self.n_samples_per_axis = n
        self.labels_plane       = ["XY", "ZX", "YZ"]

    def __iter__(self):

        for label_plane in self.labels_plane:
            for coordinate_axis in self.get_coordinates_axis():

                data_slice = self.data[coordinate_axis, ...]
                yield data_slice, label_plane, coordinate_axis

            # self.rotate_about_vector_by_theta(k= [1, 1, 1], theta = 120, center=False)
            self.rotate_about_111_by_120()
            
    def rotate_about_vector_by_theta(self, k, theta, center=True):

        # convert to radians
        theta = theta * np.pi/180
        # normalize
        k     = np.asarray(k)/(sum([i**2 for i in k])**(1/2))

        # create rodrigues rotation matrix 
        # with extra row and column for component dimension
        K = np.asarray([
            [    0, -k[2],  k[1]],
            [ k[2],     0, -k[0]],
            [-k[1],  k[0],     0],
        ])
        rotation         = np.eye(4)
        rotation[:3, :3] = np.eye(3) + np.sin(theta)*K + (1-np.cos(theta))*np.matmul(K,K)

        # rotation happens along 0,0,0
        # shift to centroid to 0,0,0 and unshift
        # to ratate about centroid
        shift = np.array(
            [
                [1, 0, 0, self.data.shape[2]/2],
                [0, 1, 0, self.data.shape[1]/2],
                [0, 0, 1, self.data.shape[0]/2],
                [0, 0, 0,                    1]
            ]
        )
        unshift = np.array(
            [
                [1, 0, 0, -self.data.shape[2]/2],
                [0, 1, 0, -self.data.shape[1]/2],
                [0, 0, 1, -self.data.shape[0]/2],
                [0, 0, 0,                     1]
            ]
        )

        if center:
            self.data = affine_transform(self.data, shift @ rotation @ unshift, offset=0)
        else:
            self.data = affine_transform(self.data, rotation, offset=0)

    def rotate_about_111_by_120(self):

        # dream3d volume data is stored as [z,y,x]
        # which corresponds to the following:
        #     ^ Y
        #     |__\ 
        #    /   / X
        # Z L
        # this function equates to rotating 120 degrees
        # about the [1,1,1] vector.
        # running this function results in a 3 state loop:
        # XY -> ZX:
        #        ^ Y    |      ^ X |     ^ X   
        #        |__\   | Y /__|   |     |__\  
        #       /   / X |   \ /    |    /   / Z
        #    Z L        |  Z L     | Y L       
        # ZX -> YZ
        #        ^ X    |      ^ Z |     ^ Z   
        #        |__\   | X /__|   |     |__\  
        #       /   / Z |   \ /    |    /   / Y
        #    Y L        |  Y L     | X L       
        # YZ -> XY
        #        ^ Z    |      ^ Y |     ^ Y   
        #        |__\   | Z /__|   |     |__\  
        #       /   / Y |   \ /    |    /   / X
        #    X L        |  X L     | Z L       

        rotations = {'axes':(2,1), 'k': 1}, {'axes':(0,2), 'k': 1}

        for rotation in rotations:
            self.data = np.rot90(self.data, **rotation)

    def get_coordinates_axis(self, exclude_edges=False):

        if exclude_edges:
            # use if also using rotate_about_111_by_120
            # rodrigues rotation uses floating point values
            # for coordinates which causes data loss at edges
            return np.linspace(1, self.data.shape[0]-2, self.n_samples_per_axis).astype(int)
        else:
            return np.linspace(0, self.data.shape[0]-1, self.n_samples_per_axis).astype(int)
        
class SliceSampler2d():

    def __init__(self, data, n=5, sample_shape=None):
        
        self.data                = data
        self.sample_shape        = sample_shape
        self.n_samples_per_plane = n

        # if self.sample_shape is None: self.sample_shape = [64, 64]
        if self.sample_shape is None: self.sample_shape = [128, 128]

    def __iter__(self):

        for label_plane, data in self.data.items():

            for n_sample in range(self.n_samples_per_plane):

                x, y       = self.get_random_coordinates(data)
                data_slice = data[y:y+self.sample_shape[0], x:x+self.sample_shape[1], :]

                yield data_slice, label_plane, [x, y]

    def get_random_coordinates(self, data):

        x = np.random.randint(0, data.shape[1] - self.sample_shape[1])
        y = np.random.randint(0, data.shape[0] - self.sample_shape[0])

        return x, y

def plot_histograms(data, path_output, map_names, n_bins=25, dpi=200, font_scale=2):

    # set seaborn figure properties
    # sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.grid' : False})
    # sns.set(font_scale=font_scale)
    # sns.set_palette("Pastel2")

    for attribute in data.attribute.unique():

        # find the row data
        data_i = data[data.attribute==attribute]
        row    = 0

        # create the figure and axes
        n_cols    = len(data_i.plane.unique())
        n_rows    = len(data_i.attribute.unique())
        figsize   = [6.4*n_cols,4.8*n_rows]
        fig, axes = plt.subplots( n_rows, n_cols, figsize=figsize, sharex='row', sharey='row' )
        axes      = np.asarray(axes).reshape((n_rows, n_cols))

        for col, plane in enumerate(data.plane.unique()):
        
            # find the column data
            data_j = data[data_i.plane==plane]
            axis   = axes[row, col]

            # find the bin sizes
            xmin = data_j.value.min()
            xmax = data_j.value.max()
            binwidth = (xmax - xmin) / n_bins
            axis.set_xlim([xmin, xmax])

            # do the histogram plot
            sns.histplot(
                data        = data_j   ,
                x           = "value"  ,
                stat        = 'percent',
                hue         = 'source' ,
                # bins        = n_bins   ,
                binwidth    = binwidth ,
                multiple    = 'dodge'  ,
                common_norm = False    ,
                shrink      = 0.9      ,
                ax          = axis
            )

            axis.grid(visible=True, linestyle=':', linewidth=1)
            # remap the feature attribute names in the x axis
            axis.set_xlabel (f"{map_names[attribute]} ({plane})", fontsize  = 20)
            axis.set_ylabel ('Fraction [%]'                     , fontsize  = 20)
            axis.tick_params(axis = 'both', which = 'major'     , labelsize = 20)

            # set only one legend on the figure
            axis.get_legend().set_visible(False)
            if col==0:
                # handles, labels = axis.get_legend_handles_labels()
                handles = axis.get_legend().legend_handles
                labels  = [ text._text for text in axis.get_legend().texts]
                fig.legend(handles, labels)
                sns.move_legend(
                    obj            = fig                     ,
                    loc            = "upper center"          ,
                    bbox_to_anchor = (0.5, 1.0)              ,
                    ncol           = len(data.source.unique()),
                    title          = None                    ,
                    frameon        = False                   ,
                    fontsize       = 20
                )
            
            if col > 0:
                axis.set_ylabel('', fontsize=0)

        fig.subplots_adjust(bottom=0.15, right=0.95, top=0.85, hspace=0.25)

        # fill name placeholders
        replacements = {
            '[attribute]': attribute,
            '[type plot]': 'histogram'
        }
        path_output_finalized = path_output
        for placeholder, replacement in replacements.items():
            path_output_finalized = path_output_finalized.replace(placeholder, replacement)

        # save the figure
        os.makedirs(os.path.split(path_output_finalized)[0], exist_ok=True)
        fig.savefig(path_output_finalized, dpi=dpi)
        plt.close(fig)

def plot_boxes(data, path_output, map_names, n_bins=25, dpi=200, font_scale=2):

    # assume only one attribute
    attribute = data.attribute.unique()[0]

    # create the figure and axes
    figsize   = [6.4, 4.8]
    fig, axis = plt.subplots( 1, 1, figsize=figsize, sharex='row', sharey='row' )

    # do the box plot
    sns.boxplot(data=data, x='plane', y='value', hue='source')

    # #########################################################
    # ### Sneak in the true (non sampled) 2D area fractions ###
    # #########################################################

    # ys_raw = [
    #     0.6228888888888889,
    #     0.6004888888888888,
    #     0.6004888888888888
    # ]

    # n_planes = len(data.plane.unique())
    # xmin, xmax = axis.get_xlim()
    # length_line = (xmax-xmin)/n_planes
    # xlims = [
    #     [start, stop]
    #     for start, stop 
    #     in zip( 
    #         np.linspace(xmin,xmax-length_line,n_planes), 
    #         np.linspace(xmin+length_line,xmax,n_planes) 
    #     )
    # ]

    # for plane in range(n_planes):
    #     axis.hlines(
    #         ys_raw[plane],
    #         *xlims[plane],
    #         colors='blue',
    #         linestyles='--',
    #         label='EBSD True'
    #     )

    
    # i = -(n_planes-1)*2
    # axis.legend()
    # handles = axis.get_legend().legend_handles
    # labels  = [ text._text for text in axis.get_legend().texts ]
    # axis.legend(handles[:i], labels[:i])

    # #########################################################
    # #########################################################
    # #########################################################

    # set only one legend on the figure
    axis.get_legend().set_visible(False)
    # handles, labels = axis.get_legend_handles_labels()
    handles = axis.get_legend().legend_handles
    labels  = [ text._text for text in axis.get_legend().texts ]
    fig.legend(handles, labels)
    sns.move_legend(
        obj            = fig                         ,
        loc            = "upper center"              ,
        bbox_to_anchor = (0.5, 1.0)                  ,
        # ncol           = len(data.source.unique())//2,
        ncol           = len(data.source.unique())   ,
        title          = None                        ,
        frameon        = False                       ,
        fontsize       = 20
        # fontsize       = 18
    )

    axis.grid(visible=True, linestyle=':', linewidth=1)
    axis.set_xlabel('')
    axis.set_ylabel(map_names[attribute], fontsize=20)
    axis.tick_params(axis = 'both', which = 'major'     , labelsize = 20)

    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.80, hspace=0.25)
    # fig.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.70, hspace=0.25)

    # fill name placeholders
    replacements = {
        '[attribute]': attribute,
        '[type plot]': 'box'
    }
    path_output_finalized = path_output
    for placeholder, replacement in replacements.items():
        path_output_finalized = path_output_finalized.replace(placeholder, replacement)

    # save the figure
    os.makedirs(os.path.split(path_output_finalized)[0], exist_ok=True)
    fig.savefig(path_output_finalized, dpi=dpi)
    plt.close(fig)

#############################################################
###################### BEGIN EXECUTION ######################
#############################################################

if __name__ == "__main__":

    # import configuration file
    config = create_dict(path_config)

    data = pd.DataFrame()

    # import the datasets from dream3d
    for name_input, data_input in config.get('/input').items():
        data_container = Dream3dDataContainer(name_input, config.get(f'input/{name_input}'))

        # perform the operations on the data
        for name_operation, data_operation in config.get('/operations').items():
            try:
                return_operation = getattr(data_container, name_operation)(**data_operation)
            except AttributeError as e:
                print(f"{name_input}: {e}")
                continue
            
            # compile the return values into a pandas dataframe
            if return_operation is not None:
                for name_plane, data_plane in return_operation.items():
                    data = pd.concat(
                        [
                            data,
                            pd.DataFrame({
                                'source'   : name_input,
                                'plane'    : name_plane,
                                'attribute': name_operation,
                                'value'    : data_plane
                            })
                        ],
                        ignore_index=True
                    )

    # plot the histograms
    for attribute in data.attribute.unique():
        plot_histograms(
            data[data.attribute == attribute],
            **config.get('output')
        )
        plot_boxes(
            data[data.attribute == attribute],
            **config.get('output')
        )