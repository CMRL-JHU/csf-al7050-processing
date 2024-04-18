# If the mask from SliceGAN is especially large, it may be beneficial
# to only sample a rectangular subdomain from somewhere in larger mask.
# This script does exactly that. If run as is, it samples the entire subdomain,
# however if run specifying an --x, --y, and/or --z, it samples only the
# portions of the domain within the specified bounds.
# EX:
#    python pipeline_08_ChooseSubdomain.py --x 0 30 --y 0 30 --z 0 30
#    samples only the portion of the original domain between voxels 0 and 20
#    along every axis
# It is highly recommended that a very small sample size such as the example
# is used in test cases


import h5py
import argparse
import utils_dream3d

path_input            = "pipeline_input/2-slicegan_data.dream3d"
path_output           = "pipeline_output/8-subdomain.dream3d"
path_Geometry         = "/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY"
path_CellData         = "/DataContainers/ImageDataContainer/CellData"
path_CellEnsembleData = "/DataContainers/ImageDataContainer/CellEnsembleData"

# gui sample
# x = [None , None ]
# y = [None , None ]
# z = [None , None ]
# entire domain
x = ['min', 'max']
y = ['min', 'max']
z = ['min', 'max']

parser = argparse.ArgumentParser()
parser.add_argument("--inp", dest="path_input" , default=path_input , type=str ,            help="Specify input file location")
parser.add_argument("--out", dest="path_output", default=path_output, type=str ,            help="Specify input file location")
parser.add_argument("--x"  , dest="x"          , default=x          , type=list, nargs="+", help="Specify bounds of subdomain cube")
parser.add_argument("--y"  , dest="y"          , default=y          , type=list, nargs="+", help="Specify bounds of subdomain cube")
parser.add_argument("--z"  , dest="z"          , default=z          , type=list, nargs="+", help="Specify bounds of subdomain cube")
args        = parser.parse_args()
path_input  = args.path_input
path_output = args.path_output
x           = args.x
y           = args.y
z           = args.z

# if dimensions not specified through command line,
# open GUI for editing
if x == y == z == [None,None]:

    import tkinter as tk
    
    # convert input strings to integers for x,y,z
    def submit(entries,window):
        
        global x; x = [int(i.strip()) for i in entries[0].split(',')]
        global y; y = [int(i.strip()) for i in entries[1].split(',')]
        global z; z = [int(i.strip()) for i in entries[2].split(',')]
        
        window.destroy()
    
    window = tk.Tk()
    
    label_instructions = tk.Label(window, text="Input Min,Max values for each dimension.\n   Ex: 0, 50")
    label_instructions.grid(row=0, column=1)
    
    label_x = tk.Label(window, text="X: ");label_x.grid(row=1, column=0)
    entry_x = tk.Entry(window)            ;entry_x.grid(row=1, column=1)
    label_y = tk.Label(window, text="Y: ");label_y.grid(row=2, column=0)
    entry_y = tk.Entry(window)            ;entry_y.grid(row=2, column=1)
    label_z = tk.Label(window, text="Z: ");label_z.grid(row=3, column=0)
    entry_z = tk.Entry(window)            ;entry_z.grid(row=3, column=1)
    
    button = tk.Button( window, text="Select", command = lambda: submit([entry_x.get(),entry_y.get(),entry_z.get()],window) )
    button.grid(row=4, column=1)
    
    window.mainloop()

# if we want the whole cube
elif x == y == z == ['min', 'max']:

    with h5py.File(path_input , 'r') as file_input:
        dims = file_input[path_Geometry+"/"+"DIMENSIONS"][...]

    x = [0, dims[0]]
    y = [0, dims[1]]
    z = [0, dims[2]]

# if dimensions were specified through command line,
# convert string to integer
else:

    x = [int("".join(val)) for val in x]
    y = [int("".join(val)) for val in y]
    z = [int("".join(val)) for val in z]

points = [
    [x[0], y[0], z[0]],
    [x[1], y[1], z[1]]
]

bounds = [x,y,z]
dims   = [bound[1]-bound[0] for bound in bounds]

with h5py.File(path_input , 'r') as file_input:
    
    geometry = {
        "dims"       : dims,
        "origin"     : file_input[path_Geometry+"/"+"ORIGIN"    ],
        "size_voxels": file_input[path_Geometry+"/"+"SPACING"   ]
    }
    
    utils_dream3d.create_file_dream3d(
        path_output
    )
    utils_dream3d.insert_geometry(
        path_output,
        geometry,
        path_Geometry
    )
    utils_dream3d.copy_crystallography(
        path_output,
        path_input,
        path_CellEnsembleData,
        path_CellEnsembleData
    )
    for name in file_input[path_CellData].keys():
        utils_dream3d.insert_attribute_array(
            path_output                                                                      ,
            path_CellData                                                                    ,
            name                                                                             ,
            data = file_input[path_CellData+"/"+name][...][z[0]:z[1],y[0]:y[1],x[0]:x[1],...],
            dtype = file_input[path_CellData+"/"+name].attrs["ObjectType"].decode('UTF-8')   ,
            attribute_matrix_type = "Cell"
        )
        
utils_dream3d.make_xdmf(path_output)
