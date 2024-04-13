# This script aids the user in curve fitting the microstructural characteristics
# such that the fitting parameters can be input into Dream3D's "Stats Generator"
# filter.
# To use, open this script with python, and a very small window should appear
# with a button prompting you to open a file. Use this navigation panel to open
# a *.dream3d file containing a "Statistics" matrix created by "Stats Generator".
# if you're following along with the example script, this file should be located
# under "./pipeline_output/5-statistics-StatsGenerator.dream3d".
# Once the file has opened, you will be able to see a colomn of variables you can
# cycle left and right through just below the open/close file button. The third
# such variable down is the selector that lets you choose which distribution you'd
# like to get fitting parameters for. The most noteworthy among the being:
#    "FeatureSize Distribution" for equivalent diameters
#    "FeatureSize Vs B Over A Distribution" for aspect ratios
# The selectors further below allow you to choose which data shows in the histogram
# plots to the bottom right. Because the data array variables have variable names
# and no signifiers denoting where they came from, this script has to present the
# user with all possible options for the correct data array based on the variable
# type (int/float) and the number of components (diameter = 1, aspect ratio = 2 ).
# Once the correct distribution and the correct data array have been selected,
# the user can either click the fit all button on the top right to allow scipy
# to attempt to automatically fit to the data, or if the results were unsatisfactory,
# the fitting parameters can be added manually and the graphs will be updated
# when a new value has been added to a field. It's also important to choose the
# correct number of bins, because that will drastically affect the fidelity
# of the distribution's pdf and therefor the fidelity of the fitting parameters.
# After the user is happy with the results, these fitting parameters will have to
# be added manually to the "Stats Generator" filter of pipeline component #7.

#for importing dream3d files, copying dictionaries by value instead of reference
import h5py, copy
#for translating dream3d/python data
import utils_dream3d
import json
#for creating best fit lines
import numpy as np

from scipy import stats

from scipy.stats import kstest, ks_2samp
from scipy.stats import beta
from scipy.stats import lognorm as lognormal
from scipy.stats import exponpow as power #### maybe supposed to be powerlaw? not implemented anyways
#gui
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#for running dream3d during recursive packing
import shlex, os, platform, time, subprocess

# To Do:
    #create warning handler for all the warnings that need output
    #calculate values for ALL statistical types
    #for some reason, within the dream3d file, power distributions
    #   are called "Beta Distribution"
    #   but with ["Average", "Minimum Value"]
    #   figure out a way to recognize and deal with them properly
    
    #recursive algorithm doesn't work.
    #   for some reason, changing the statistical distribution parameters
    #   and feeding them directly into dream3d will not result in a
    #   change in the attribute distribution. Although using the
    #   candidate distribution parameters from this program and typing
    #   them manually into dream.3d's statsgenerator DOES work.
    #   I haven't a clue why this may be or how to fix it, but the
    #   process that currently works is to use this program to find the
    #   next set of candidate distribution parameters, type them into
    #   statsgenerator manually, recreate the microstructure, then
    #   use this program to find the next set of candidate distribution
    #   parameters. rinse and repeat.
# NOTE! Abandoned:
    #because the recursive algorithm didn't work for aspectratios
    #   (what i was most interested in fixing), i never bothered to
    #   develop it further
    
class utils_packing():

    def __init__(self):
    
        locations = {
            "file": {
                "location"    : [0,0],
                "span"        : [1,1],
                "button width": 8
            },
            "stats": {
                "location"    : [1,0],
                "span"        : [1,1],
                "label width" : 30,
                "button width": 2
            },
            "bins": {
                "location"    : [2,0],
                "span"        : [1,1],
                "label width" : 30,
                "button width": 2
            },
            "features": {
                "location"    : [3,0],
                "span"        : [1,1],
                "label width" : 30,
                "button width": 2
            },
            "parameters": {
                "location"    : [0,1],
                "span"        : [1,1],
                "label width" : 20,
                "button width": 20
            },
            "figures": {
                "location"    : [1,1],
                "span"        : [3,1],
                "label width" : 20,
                "button width": 20
            },
            "recursive_fit": {
                "location"    : [4,0],
                "span"        : [1,1],
                "label width" : 20,
                "button width": 20,
                "button height": 13
            }
        }
        
        self.window            = tk.Tk()
        self.events            = events_handler            ()
        self.file              = file_handler              (self.events, dict({"widget":self.window}, **locations["file"         ]))
        self.stats             = stats_handler             (self.events, dict({"widget":self.window}, **locations["stats"        ]), self.file)
        self.features          = features_handler          (self.events, dict({"widget":self.window}, **locations["features"     ]), self.file, self.stats, {"self":"features","parent":"stats"})
        self.bins              = bins_handler              (self.events, dict({"widget":self.window}, **locations["bins"         ]), self.file, self.stats, self.features)
        self.distribution_data = distribution_data_handler (self.events, self.stats, self.features, self.bins)
        self.parameters        = parameters_handler        (self.events, dict({"widget":self.window}, **locations["parameters"   ]), self.file, self.stats, self.bins, self.distribution_data)
        self.figures           = figures_handler           (self.events, dict({"widget":self.window}, **locations["figures"      ]))
        self.plots             = plot_data                 (self.events, self.file, self.stats, self.features, self.bins, self.distribution_data, self.parameters, self.figures)
        # self.recursive_fit     = recursive_fit             (self.events, dict({"widget":self.window}, **locations["recursive_fit"]), self.file, self.stats, self.parameters, self.distribution_data)
        
        self.run()
        
    def close(self):
        
        self.events.event("quit")
        self.window.destroy()
        quit()
        
    def run(self):
        
        self.window.protocol("WM_DELETE_WINDOW", self.close)
        self.window.mainloop()

class canvas_handler():

    def __init__(self, window):
    
        self.window   = window["widget"  ]
        self.location = window["location"]
        self.span     = window["span"    ]
        self.elements = {}
        
        self.canvas   = tk.Canvas(self.window)
        location_row, location_column = self.location
        span_row    , span_column     = self.span
        self.canvas.grid(row=location_row, column=location_column, rowspan=span_row, columnspan=span_column)
        
    def add(self, name, element):
        
        self.elements[name] = element
        
    def display(self):
    
        # Place interface items
        for element in self.elements:
            
            # Check if the item is already in the canvas
            if self.elements[element]["widget"].winfo_ismapped():
                continue
            
            # Create all elements
            location_row, location_column = self.elements[element]["location"]
            span_row    , span_column     = self.elements[element]["span"    ]
            self.elements[element]["widget"].grid(row=location_row, column=location_column, rowspan=span_row, columnspan=span_column)
            
    def clear(self):
        
        for element in [*self.elements.keys()]:
            self.kill(element)
            
    def kill(self, element):
        
        if not element in [*self.elements.keys()]:
            return
        
        # figures stay in memory if not explicitly closed
        if "figure" in [*self.elements[element].keys()]:
            plt.close(self.elements[element]["figure"])
        
        # widgets stick around if not explicitly destroyed
        self.elements[element]["widget"].destroy()
        
        del self.elements[element]
            
    def destroy(self):
    
        self.clear()
        self.canvas.destroy()
        
class events_handler():

    def __init__(self):
    
        self.events = {}
        
    def add_event(self, event):
        
        if not event in [*self.events.keys()]:
            self.events[event] = []
        else:
            print(event+" already exists.")
        
    def add_action(self, event, action):
    
        if not event in [*self.events.keys()]:
            self.events[event] = []
        self.events[event].append(action)
        
    def event(self, event):
            
            for event in self.events[event]:
                event()
            
class file_handler():

    def __init__(self, events, window):
        
        self.events = events
        self.window = window
        self.canvas = None
        self.path   = None
        self.file   = None
        
        self.create_events()
        self.create_canvas()
        
    def create_events(self):
        
        self.events.add_event("file open" )
        self.events.add_event("file close")
        self.events.add_event("file export")
        
        self.events.add_action("quit", self.destroy)
        
    def create_canvas(self):
    
        self.canvas = canvas_handler(self.window)
    
        button =  {
            "text"    : tk.StringVar(),
            "location": [0,0],
            "span"    : [1,1]
        }
        
        button["text"].set("Open File")
        button["widget"] = \
            tk.Button(
                self.canvas.canvas,
                textvariable=button["text"],
                command=lambda: self.file_open(),
                width=self.window["button width"]
            )
        self.canvas.add("button", button)
        
        self.canvas.display()

    def file_open(self):

        # open file browser
        file_path = filedialog.askopenfilename(
            initialdir = ".",
            title = "Select a File",
            filetypes = (
                    ("DREAM.3D files", "*.dream3d*"),
                    ("all files","*.*")
            ),
        )
        
        if file_path:
        
            # update the file path and import the file data
            self.path = file_path
            
            # open the file
            self.file = h5py.File(self.path, 'r+')
            
            # swap action of explore button
            self.canvas.elements["button"]["text"].set("Close File")
            self.canvas.elements["button"]["widget"].config(command=lambda: self.file_close())
            
            # trigger event
            self.events.event("file open")
            
    def file_close(self):
    
        # close the file
        self.file.close()
    
        # swap action of explore button
        self.canvas.elements["button"]["text"].set("Select a File")
        self.canvas.elements["button"]["widget"].config(command=self.file_open)
        
        # trigger event
        self.events.event("file close")
        
    def export(self, stats):
        
        # list of valid parameters for each distribution type
        parameter_map = {
            "Log Normal Distribution": [
                "Average",
                "Standard Deviation"
            ],
            "Beta Distribution": [
                "Alpha",
                "Beta"
            ],
            "Power Distribution": [
                "Alpha",
                "Minimum Value"
            ]
        }
        
        # export all new parameters
        for statistic in [*stats.data["Modified"].keys()]:
            for phase     in [*stats.data["Modified"][statistic].keys()]:
                for property_ in [*stats.data["Modified"][statistic][phase    ].keys()]:
                    
                    # hdf5 can use either dictionary addressing or OS like paths
                    # build an OS like path here for easy addressing
                    path = "/".join([statistic, phase, property_])
                
                    # Change attribute type
                    self.file[path].attrs["Distribution Type"] = \
                        utils_dream3d.format_string(
                            stats.data["Modified"][statistic][phase][property_]["Distribution Type"]
                        )
                
                    for parameter_original, parameter_modified in \
                        zip(
                            [*stats.data["Original"][statistic][phase    ][property_].keys()],
                            [*stats.data["Modified"][statistic][phase    ][property_].keys()]
                        ):
                    
                        # extra parameters ("Distribution Type", "bins") were put into data for ease of use
                        # skip any parameters that don't match this map of valid parameters
                        if not parameter_modified in sum([*parameter_map.values()],[]):
                            continue
                        
                        # Delete old parameters
                        del self.file[path+"/"+parameter_original]
                        
                        # Gather the user entry data and format it for hdf5
                        data = stats.data["Modified"][statistic][phase][property_][parameter_modified]
                        data = data.reshape((data.size,1))
                        
                        dims                  = list(data.shape[:-1][::-1])
                        components            = [data.shape[-1]]
                        tuple_axis_dimensions = ','.join(['='.join([name,str(value)]) for name,value in zip(['x','y','z'],dims)])
                        
                        # Export attribute array
                        group = self.file[path]
                        dataset = group.create_dataset(parameter_modified, data=data)
                        dataset.attrs["ComponentDimensions"]   = np.uint64(components)
                        dataset.attrs["DataArrayVersion"]      = np.int32([2])
                        dataset.attrs["ObjectType"]            = utils_dream3d.format_string("DataArray<float>")
                        dataset.attrs["Tuple Axis Dimensions"] = utils_dream3d.format_string(tuple_axis_dimensions)
                        dataset.attrs["TupleDimensions"]       = np.uint64(dims)
                        
        self.events.event("file export")
        
    def destroy(self):
        
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("file close")
        
class stats_handler():

    def __init__(self, events, window, file):

        self.file   = file
        self.events = events
        self.window = window
        self.canvas = None

        self.data    = {
            "Original": {},
            "Modified": {}
        }
        self.current = {
            "path components":{
                "statistic": "",
                "phase"    : "",
                "property" : ""
            },
            "path":{
                "statistic": "",
                "phase"    : "",
                "property" : ""
            },
            "name list" : [],
            "value list": [],
            "required features": {}
        }
        
        # create initial distribution current indecies
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event ("stats cycle")
        self.events.add_event ("stats open" )
        self.events.add_event ("stats close")
    
        self.events.add_action("quit"       , self.destroy    )
        self.events.add_action("file close" , self.destroy    )
        self.events.add_action("file open"  , self.initialize )
        self.events.add_action("file export", self.import_data)
        
    def import_data(self):
    
        # initialize data container
        self.data = {
            "Original": {},
            "Modified": {}
        }
        
        file = self.file.file
        
        # create distribution data from hdf5 file
        for statistic in self.find_statistics(file):
            phases = {key:{} for key in file[statistic].keys() if isinstance(file[statistic+"/"+key], h5py.Group)}
            self.data["Original"][statistic] = phases
            for phase in phases:
                properties = {key:{} for key in file[statistic+"/"+phase].keys() if isinstance(file[statistic+"/"+phase+"/"+key], h5py.Group)}
                self.data["Original"][statistic][phase] = properties
                for property_ in properties:
                    parameters = {key:value[...].flatten() for key, value in zip(file[statistic+"/"+phase+"/"+property_].keys(),file[statistic+"/"+phase+"/"+property_].values())}
                    parameters["Distribution Type"] = file[statistic+"/"+phase+"/"+property_].attrs["Distribution Type"].decode()
                    self.data["Original"][statistic][phase][property_] = parameters
                    self.data["Original"][statistic][phase][property_]["bins"] = 50
                        
        self.data["Modified"] = copy.deepcopy(self.data["Original"])
        
    def find_statistics(self, items):
    
        # recursively dig through the groups in the hdf5 file to find
        # groups with the attribute "ObjectType" equal "Statistics"
        # these groups contain our statistics information
        hits = []
        for item in items:
            if isinstance(items[item], h5py.Dataset):
                continue
            if "ObjectType" in [*items[item].attrs.keys()]:
                if items[item].attrs["ObjectType"].decode() == "Statistics":
                    return items[item].name
            hit = self.find_statistics(items[item])
            hits += hit if type(hit) in [list,tuple] else [hit]
        if any(hits):
            return [hit for hit in hits if hit is not []]
        return hits
        
    def initialize(self):
    
        # import all valid statistical data from file
        self.import_data()
    
        # fail out gracefully if no stats
        if len([*self.data["Original"].keys()]) == 0:
            ##### TO DO #####
            
            # Output no statistics found error
            
            ##### TO DO #####
            return
        
        # set initial path to the first valid path found
        statistic = self.current["path components"]["statistic"] = [*self.data["Original"].keys()][0]
        phase     = self.current["path components"]["phase"    ] = [*self.data["Original"][statistic].keys()][0]
        property_ = self.current["path components"]["property" ] = [*self.data["Original"][statistic][phase].keys()][0]
        
        # update current["path", "name list", "value list"]
        self.update_current()
        
        # create the user interface module that allows for cycling through valid statistics paths
        self.create_canvas()
        
        # trigger event
        self.events.event("stats open")
        
    def create_canvas(self):
    
        self.canvas = canvas_handler(self.window)
    
        label =  {
            "text"    : tk.StringVar(),
            "location": [0, 0],
            "span"    : [1, 3]
        }
        label["text"].set("Select Distribution to Compare")
        label["widget"] = \
            tk.Label (
                self.canvas.canvas,
                textvariable=label["text"],
                font='Helvetica 12 bold',
                width=self.window["label width"]
            )
        self.canvas.add("label", label)
    
        for i, name in enumerate(self.current["name list"]):
    
            button = {
                "widget": tk.Button(
                    self.canvas.canvas,
                    text = "←",
                    command=lambda name_=name: self.cycle(name_, -1),
                    width = self.window["button width"] 
                ),
                "location": [1+i, 0],
                "span"    : [1  , 1]
            }
            self.canvas.add(f"button cycle negative {name}", button)
            
            label = {
                "text"    : tk.StringVar(),
                "location": [1+i, 1],
                "span"    : [1  , 1]
            }
            label["text"].set(self.current["path components"][name])
            label["widget"] = \
                tk.Label (
                    self.canvas.canvas,
                    textvariable=label["text"],
                    width = self.window["label width"]
                )
            self.canvas.add(f"label {name}", label)
            
            button = {
                "widget": tk.Button(
                    self.canvas.canvas,
                    text = "→",
                    command=lambda name_=name: self.cycle(name_, 1),
                    width = self.window["button width"] 
                ),
                "location": [1+i, 2],
                "span"    : [1  , 1]
            }
            self.canvas.add(f"button cycle positive {name}", button)
            
        self.canvas.display()
        
    def cycle(self, name, direction):
        
        # grab the list of values for each group
        statistic, phase, property_ = self.current["value list"]
        if name == "statistic":
            keys = [key for key in self.data["Original"].keys()]
        elif name == "phase":
            keys = [key for key in self.data["Original"][statistic].keys()]
        elif name == "property":
            keys = [key for key in self.data["Original"][statistic][phase].keys()]
        else:
            return
            
        # Rotate the array of possible values to current +/- desired direction,
        # and return zeroth value of the rotated array as the new current index
        current  = self.current["path components"][name]
        index    = keys.index(current)
        keys     = keys[index+direction:]+keys[:index+direction]
        name_new = keys[0]
        
        # commit new name
        self.current["path components"][name] = name_new
        
        # update the current hdf5 path
        self.update_current()
        self.update_labels ()
        
        # trigger cycle event
        self.events.event("stats cycle")
        
    def parameter_map(self):
    
        # list of valid parameters for each distribution type
        return {
            "Log Normal Distribution": [
                "Average",
                "Standard Deviation"
            ],
            "Beta Distribution": [
                "Alpha",
                "Beta"
            ],
            "Power Distribution": [
                "Alpha",
                "Minimum Value"
            ]
        }
        
    def get_required_features(self):
        
        statistic, phase, property_ = self.current["value list"]
        if   property_ == "FeatureSize Distribution":
            self.current["required features"] = {
                "Equivalent Diameter": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                }
            }
        elif property_ == "FeatureSize Vs B Over A Distributions":
            self.current["required features"] = {
                "Equivalent Diameter": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                },
                "Aspect Ratio (B/A)": {
                    "index"   :0,
                    "shape"   : [2,2],
                    "datatype": "DataArray<float>"
                }
            }
        elif property_ == "FeatureSize Vs C Over A Distributions":
            self.current["required features"] = {
                "Equivalent Diameter": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                },
                "Aspect Ratio (C/A)": {
                    "index"   : 1,
                    "shape"   : [2,2],
                    "datatype": "DataArray<float>"
                }
            }
        elif property_ == "FeatureSize Vs Neighbors Distributions":
            self.current["required features"] = {
                "Equivalent Diameter": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                }
                ##### TO DO ######
                
                # figure out whatever other information is needed here
                
                ##### TO DO ######
            }
        elif property_ == "FeatureSize Vs Omega3 Distributions":
            self.current["required features"] = {
                "Equivalent Diameter": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                },
                "Omega 3": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<float>"
                }
            }
        
    def update_current(self):
    
        # update paths
        names  = []
        values = []
        for name, value in zip(self.current["path components"].keys(), self.current["path components"].values()):
            names  += [name]
            values += [value]
            self.current["path"][name] = "/".join(values)
        self.current["name list" ] = names
        self.current["value list"] = values
        
        # update required features
        self.get_required_features()
        
    def update_labels(self):
    
        for name, value in zip(self.current["name list"], self.current["value list"]):
            self.canvas.elements[f"label {name}"]["text"].set(value)

    def destroy(self):
    
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("stats close")
    
class features_handler():

    def __init__(self, events, window, file, stats, names):
        
        self.events      = events
        self.window      = window
        self.canvas      = None
        self.file        = file
        self.stats       = stats
        # specifies the heirarchy of this class
        #     self.names = {'self':str(), 'parent':str()}
        self.names       = names
        
        self.possible = {}
        self.current  = {}
        
        # create initial distribution current indecies
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event (f"{self.names['self']} cycle")
        self.events.add_event (f"{self.names['self']} open" )
        self.events.add_event (f"{self.names['self']} close")
    
        self.events.add_action("quit"       , self.destroy     )
        self.events.add_action("file close" , self.destroy     )
        self.events.add_action(f"{self.names['parent']} close", self.destroy     )
        self.events.add_action(f"{self.names['parent']} open" , self.initialize  )
        self.events.add_action(f"{self.names['parent']} cycle", self.reinitialize)
        
        
    def initialize(self):
        
        # get features with the correct shape, number of components, and datatype
        # according to stats.current["required features"]
        self.possible = {}
        self.get_possible(self.file.file)
        
        # fail out gracefully if no features meet requirements
        if len([*self.possible.keys()]) == 0:
            ##### TO DO #####
            
            # Output no sizes found error
            
            ##### TO DO #####
            return
        
        # create initial pick for each required feature as first possible
        self.current = {}
        for feature in [*self.possible.keys()]:
        
            # create datastructure
            self.current[feature]  = {
                "path components": {
                    "datacontainer"  : None,
                    "attributematrix": None,
                    "dataset"        : None
                },
                "path": {
                    "datacontainer"  : None,
                    "attributematrix": None,
                    "dataset"        : None
                },
                "name list" : [],
                "value list": [],
                "value"     : None
            }
            
            # set initial path to the first valid path found
            datacontainer   = self.current[feature]["path components"]["datacontainer"  ] = [*self.possible[feature].keys()][0]
            attributematrix = self.current[feature]["path components"]["attributematrix"] = [*self.possible[feature][datacontainer].keys()][0]
            dataset         = self.current[feature]["path components"]["dataset"        ] = [*self.possible[feature][datacontainer][attributematrix].keys()][0]
            
        # update current["path", "name list", "value list", value]
        self.update_current()
        
        # create the user interface module that allows for cycling through valid statistics paths
        self.create_canvas()
        
        # trigger event
        self.events.event(f"{self.names['self']} open")
        
    def reinitialize(self):
        
        # get features with the correct shape, number of components, and datatype
        # according to stats.current["required features"]
        self.possible = {}
        self.get_possible(self.file.file)
        
        # fail out gracefully if no features meet requirements
        if len([*self.possible.keys()]) == 0:
            ##### TO DO #####
            
            # Output no sizes found error
            
            ##### TO DO #####
            return
        
        # create initial pick for each required feature as first possible
        self.current = {}
        for feature in [*self.possible.keys()]:
        
            # create datastructure
            self.current[feature]  = {
                "path components": {
                    "datacontainer"  : None,
                    "attributematrix": None,
                    "dataset"        : None
                },
                "path": {
                    "datacontainer"  : None,
                    "attributematrix": None,
                    "dataset"        : None
                },
                "name list": [],
                "value list": [],
                "value": None
            }
            
            # set initial path to the first valid path found
            datacontainer   = self.current[feature]["path components"]["datacontainer"  ] = [*self.possible[feature].keys()][0]
            attributematrix = self.current[feature]["path components"]["attributematrix"] = [*self.possible[feature][datacontainer].keys()][0]
            dataset         = self.current[feature]["path components"]["dataset"        ] = [*self.possible[feature][datacontainer][attributematrix].keys()][0]
            
        # update current["path", "name list", "value list", value]
        self.update_current()
        
        # update the user interface module that allows for cycling through valid statistics paths
        self.update_labels()
        
    def get_possible(self, items, depth=1, max_depth=4):
    
        # recursively dig through the groups in the hdf5 file to find
        # datasets at the right depth with the required attributes
        for item in items:
            item = items[item]
            
            if isinstance(item, h5py.Dataset):
                
                # if it's a dataset, not a group AND it's hit the correct depth
                if depth == max_depth:
                    
                    # check if all attribute requirements are met
                    for feature in [*self.stats.current["required features"].keys()]:
                        shape        = self.stats.current["required features"][feature]["shape"   ][0]
                        n_components = self.stats.current["required features"][feature]["shape"   ][1]
                        if n_components == "any":
                            n_components = item[:].shape[-1]
                        datatype     = self.stats.current["required features"][feature]["datatype"]   
                    
                        # if it has the correct shape, number of components, and data type
                        # we have a possible match! return it
                        if \
                            len(item[:].shape) == shape \
                            and item[:].shape[-1] == n_components \
                            and "ObjectType" in item.attrs.keys() \
                            and item.attrs["ObjectType"].decode("utf-8") == datatype:
                            
                            datacontainer, attributematrix, dataset = item.name[1:].split("/")[1:]
                            
                            if not feature in [*self.possible.keys()]:
                                self.possible[feature] = {}
                            if not datacontainer   in [*self.possible[feature].keys()]:
                                self.possible[feature][datacontainer  ] = {}
                            if not attributematrix in [*self.possible[feature][datacontainer  ]]:
                                self.possible[feature][datacontainer  ][attributematrix] = {}
                            if not dataset         in [*self.possible[feature][datacontainer  ][attributematrix]]:
                                self.possible[feature][datacontainer  ][attributematrix][dataset        ] = self.file.file[item.name][...]
                    
            elif depth < max_depth:
                self.get_possible(item, depth=depth+1)
                
    def create_canvas(self):
    
        self.canvas = canvas_handler(self.window)
    
        # the length of this is variable on purpose just in case a distribution
        # requires more than one feature. it must be cleared every update.
        self.update_labels()
        
    def cycle(self, feature, name, direction):
        
        # grab the list of values for each group
        datacontainer, attributematrix, dataset = self.current[feature]["value list"]
        if   name == "datacontainer":
            keys = [*self.possible[feature].keys()]
        elif name == "attributematrix":
            keys = [*self.possible[feature][datacontainer].keys()]
        elif name == "dataset":
            keys = [*self.possible[feature][datacontainer][attributematrix].keys()]
        else:
            return
            
        # Rotate the array of possible values to current +/- desired direction,
        # and return zeroth value of the rotated array as the new current index
        current  = self.current[feature]["path components"][name]
        index    = keys.index(current)
        keys     = keys[index+direction:]+keys[:index+direction]
        name_new = keys[0]
        
        # commit new name
        self.current[feature]["path components"][name] = name_new
        
        # set initial path to the first valid path found
        if name == "datacontainer":
            datacontainer   = name_new
            attributematrix = self.current[feature]["path components"]["attributematrix"] = [*self.possible[feature][datacontainer].keys()][0]
            dataset         = self.current[feature]["path components"]["dataset"        ] = [*self.possible[feature][datacontainer][attributematrix].keys()][0]
        if name == "attributematrix":
            attributematrix = name_new
            dataset         = self.current[feature]["path components"]["dataset"        ] = [*self.possible[feature][datacontainer][attributematrix].keys()][0]
        
        # update the current hdf5 path
        self.update_current()
        self.update_labels ()
        
        # save entries to self.data_new
        self.events.event(f"{self.names['self']} cycle")
        
    def update_current(self):
        
        # update path information
        for feature in [*self.current.keys()]:
            names  = []
            values = []
            for name, value in zip(self.current[feature]["path components"].keys(), self.current[feature]["path components"].values()):
                names  += [name]
                values += [value]
                self.current[feature]["path"][name] = "/".join(values)
            self.current[feature]["name list" ] = names
            self.current[feature]["value list"] = values
        
        # update size value
        for feature in [*self.current.keys()]:
            self.current[feature]["value"] = self.file.file["/DataContainers"+"/"+self.current[feature]["path"]["dataset"]][...]
        
    def update_labels(self):
        
        # clear out all canvas elements if they exist
        self.canvas.clear()
        
        # make (or remake) all canvas elements
        for i, feature in enumerate([*self.current.keys()]):
            
            ni = len(self.current[feature]["name list"])+1
        
            label =  {
                "text"    : tk.StringVar(),
                "location": [0+i*ni, 0],
                "span"    : [1     , 3]
            }
            label["text"].set(f"Select  {feature}")
            label["widget"] = \
                tk.Label (
                    self.canvas.canvas,
                    textvariable=label["text"],
                    font='Helvetica 12 bold',
                    width=self.window["label width"]
                )
            self.canvas.add(f"label {feature}", label)
        
            for j, name in enumerate(self.current[feature]["name list"]):
        
                button = {
                    "widget": tk.Button(
                        self.canvas.canvas,
                        text = "←",
                        command=lambda feature_=feature, name_=name: self.cycle(feature_, name_, -1),
                        width = self.window["button width"] 
                    ),
                    "location": [1+i*ni+j, 0],
                    "span"    : [1       , 1]
                }
                self.canvas.add(f"button {feature} cycle negative {name}", button)
                
                label = {
                    "text"    : tk.StringVar(),
                    "location": [1+i*ni+j, 1],
                    "span"    : [1       , 1]
                }
                label["text"].set(self.current[feature]["path components"][name])
                label["widget"] = \
                    tk.Label (
                        self.canvas.canvas,
                        textvariable=label["text"],
                        width = self.window["label width"]
                    )
                self.canvas.add(f"label {feature} {name}", label)
                
                button = {
                    "widget": tk.Button(
                        self.canvas.canvas,
                        text = "→",
                        command=lambda feature_=feature, name_=name: self.cycle(feature_, name_, 1),
                        width = self.window["button width"] 
                    ),
                    "location": [1+i*ni+j, 2],
                    "span"    : [1       , 1]
                }
                self.canvas.add(f"button {feature} cycle positive {name}", button)
            
        self.canvas.display()
            
    def destroy(self):
    
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event(f"{self.names['self']} close")
            
class bins_handler():

    def __init__(self, events, window, file, stats, features):
    
        self.events = events
        self.window = window
        self.canvas = None
        self.file   = file
        self.stats  = stats
        self.size   = features
        
        self.data = {
            "bounds": [],
            "size"  : []
        }
        self.current = {
            "bin"   : None,
            "bounds": None,
            "size"  : None
        }
        
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event ("bins cycle")
        self.events.add_event ("bins open" )
        self.events.add_event ("bins close")
    
        self.events.add_action("quit"          , self.destroy     )
        self.events.add_action("file close"    , self.destroy     )
        self.events.add_action("stats close"   , self.destroy     )
        self.events.add_action("features close", self.destroy     )
        self.events.add_action("features open" , self.initialize  )
        self.events.add_action("stats cycle"   , self.reinitialize)
        self.events.add_action("features cycle", self.reinitialize)
        
    def initialize(self):
        
        self.get_bins()
        
        self.update_current()
        
        self.create_canvas()
        
        self.events.event("bins open")
        
    def reinitialize(self):
    
        self.get_bins()
        
        self.update_current()
        self.update_labels()
        
    def get_bins(self):
        
        # Feature_Diameter_Info contains all the size bin information
        bin_size, range_max, range_min = self.file.file[self.stats.current["path"]["phase"]+"/"+"Feature_Diameter_Info"]
        
        # DREAM.3D picks an arbitrarily large number for bin size if it only wants one bin,
        # even if that number is larger than the maximum range. Normalize here.
        bin_size = min(bin_size, range_max)
        
        # find the bounds
        
        # get bin index
        # because size distribution always has just 1 bin, but everything else has multiple
        # we have to check what the number of values in a valid parameter of the current property
        statistic, phase, property_ = self.stats.current["value list"]
        for parameter in [*self.stats.data["Modified"][statistic][phase][property_].keys()]:
            if parameter in sum([*self.stats.parameter_map().values()],[]):
                n_parameter_values = len(self.stats.data["Modified"][statistic][phase][property_][parameter][...])
                break
        if n_parameter_values > 1:
            self.data["bounds"] = [[bin_min, bin_min+bin_size] for bin_min in np.arange(range_min, range_max, bin_size)]
        else:
            self.data["bounds"] = [[range_min, range_max]]
        
        # find the feature sizes that belong to each of those bounds
        # we kept the special zero values, so this captures every value of note in the size range for every bound
        self.data["size"] = []
        for bounds in self.data["bounds"]:
            logic = (self.size.current["Equivalent Diameter"]["value"] > bounds[0]) & (self.size.current["Equivalent Diameter"]["value"] <= bounds[1])
            self.data["size"] += [self.size.current["Equivalent Diameter"]["value"][logic]]
            
    def create_canvas(self):
    
        self.canvas = canvas_handler(self.window)
    
        label =  {
            "text"    : tk.StringVar(),
            "location": [0, 0],
            "span"    : [1, 1]
        }
        label["text"].set("Select Bin")
        label["widget"] = \
            tk.Label (
                self.canvas.canvas,
                textvariable=label["text"],
                font='Helvetica 12 bold',
                width=self.window["label width"]
            )
        self.canvas.add("label", label)
    
        combobox = {
            "value"   : tk.StringVar(),
            "location": [1, 0],
            "span"    : [1, 1]
        }
        combobox["value"].set(str(self.current["bounds"]))
        combobox["value"].trace("w", self.cycle_bin)
        combobox["widget"] = \
            ttk.Combobox(
                self.canvas.canvas,
                textvariable=combobox["value"]
            )
        combobox["widget"]["values"] = [str(bound) for bound in self.data["bounds"]]
        self.canvas.add("combobox", combobox)
            
        self.canvas.display()
        
    def cycle_bin(self, *args):
        
        # get new bound value from combobox box
        string = self.canvas.elements["combobox"]["value"].get()
        
        self.current["bounds"] = [float(val) for val in string.replace("[","").replace("]","").split(", ")]
        self.current["bin"] = self.data["bounds"].index(self.current["bounds"])
        
        # trigger event
        self.events.event("bins cycle")
        
    def update_current(self):
        
        self.current["bin"   ] = 0
        self.current["bounds"] = self.data["bounds"][0]
        self.current["size"  ] = self.data["size"  ][0]
        
    def update_labels(self):
    
        self.canvas.kill("combobox")
        
        combobox = {
            "value"   : tk.StringVar(),
            "location": [1, 0],
            "span"    : [1, 1]
        }
        combobox["value"].set(str(self.current["bounds"]))
        combobox["value"].trace("w", self.cycle_bin)
        combobox["widget"] = \
            ttk.Combobox(
                self.canvas.canvas,
                textvariable=combobox["value"]
            )
        combobox["widget"]["values"] = [str(bound) for bound in self.data["bounds"]]
        self.canvas.add("combobox", combobox)
        
        self.canvas.display()
        
    def destroy(self):
        
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("bins close")
        
        
class distribution_data_handler():

    def __init__(self, events, stats, features, bins):
        
        self.events        = events
        self.stats         = stats
        self.features      = features
        self.bins          = bins
        
        self.distribution  = None
        self.data = {
            "All"    : None,
            "Current": None
        }
        self.best_fits = {
            "All"    : None,
            "Current": None
        }
        
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event ("distributions open" )
        self.events.add_event ("distributions close")
    
        self.events.add_action("bins open"     , self.run    )
        self.events.add_action("stats cycle"   , self.run    )
        self.events.add_action("features cycle", self.run    )
        self.events.add_action("bins cycle"    , self.run    )
        
    def run(self):
        
        # Reset data and fits
        self.data = {}
        self.fits = {}
        
        # Get data for entire range and for the current bound
        for range_ in ["Current", "All"]:
            self.data[range_] = self.get_data(range_)
            
        # Find the best fit for the whole distribution
        self.find_overall_best_fit()
        
        # Trigger event to let dependents know the plots are ready
        self.events.event("distributions open")
        
    def get_data(self, bounds, data_source=None):
        
        if data_source is None:
            data_source = self.features.current
        
        if not "Equivalent Diameter" in [*data_source.keys()]:
            return
        
        if   bounds == "Current":
            bounds = self.bins.current["bounds"]
        elif bounds == "All":
            bounds = [self.bins.data["bounds"][0][0], self.bins.data["bounds"][-1][-1]]
        
        equivalent_diameters = data_source["Equivalent Diameter"]["value"].flatten()
        logic = (equivalent_diameters > bounds[0]) & (equivalent_diameters <= bounds[1])
        
        # Pull the prerequisite data and calculate the final form
        if   self.stats.current["path components"]["property"] == "FeatureSize Distribution":
            if "Equivalent Diameter" in [*data_source.keys()]:
                data = data_source["Equivalent Diameter"]["value"][...,0]
        elif self.stats.current["path components"]["property"] == "FeatureSize Vs B Over A Distributions":
            if "Aspect Ratio (B/A)"  in [*data_source.keys()]:
                data = data_source["Aspect Ratio (B/A)" ]["value"][...,0]
        elif self.stats.current["path components"]["property"] == "FeatureSize Vs C Over A Distributions":
            if "Aspect Ratio (C/A)"  in [*data_source.keys()]:
                data = data_source["Aspect Ratio (C/A)" ]["value"][...,1]
        elif self.stats.current["path components"]["property"] == "FeatureSize Vs Neighbors Distributions":
            ##### TO DO ######
            return
            ##### TO DO ######
        elif self.stats.current["path components"]["property"] == "FeatureSize Vs Omega3 Distributions":
            if "Omega 3"             in [*data_source.keys()]:
                data = data_source["Omega 3"            ]["value"][...,0]
        
        return data[logic]
        
    def fit_distribution(self, data, dist, name): 
        
        if data is None:
            return

        data = data[data>0]
        
        if len(data) == 0:
            return
            
        # find the best fit for the data using the given distribution
        if name == "Beta Distribution":
            try:
                fitted = [*dist.fit(data, floc=0, fscale=1)]
            except:
                fitted = [*dist.fit(data, floc=0.0)]
        else:
            fitted = [*dist.fit(data, floc=0.0)]
        
        # create the frozen random variate (contains pdf, cdf, and ppf)
        frv = dist(*fitted)
        
        # find Kolmogorov–Smirnov test statistic (higher means better fit)
        ks = kstest(data, dist.name, args=fitted, alternative='two-sided')[1]
        
        # find DREAM.3D parameters
        if   name == "Log Normal Distribution":
            parameters = {
                "Average"           : np.log(fitted[2]),
                "Standard Deviation": fitted[0]
            }
        elif name == "Beta Distribution":
            parameters = {
                "Alpha"             : fitted[0],
                "Beta"              : fitted[1]
            
            }
        elif name == "Power Distribution":
            parameters = {
                "Alpha"             : fitted[0],
                "Minimum Value"     : fitted[1]
            }
        
        distribution = {
            "frv"       : frv,
            "args"      : fitted,
            "parameters": parameters,
            "ks"        : ks
        }
        
        if not name is None:
            distribution["name"] = name
        
        return distribution
        
    def distribution_map(self):
        
        # map of DREAM.3D to scipy distribution names
        return {
            "Log Normal Distribution": "lognorm",
            "Beta Distribution"      : "beta"
            # "Power Distribution"     : "exponpow"
        }
        
    def find_overall_best_fit(self):
        
        # find the distribution that best fits the whole range
        best_fits = [ \
            self.fit_distribution(self.data["All"], getattr(stats, distribution), name) \
            for name, distribution in zip([*self.distribution_map().keys()], [*self.distribution_map().values()]) \
        ]
        ks_vals = [fit["ks"] for fit in best_fits if not fit is None]
        if len(ks_vals) == 0:
            return
        best_fit = best_fits[ks_vals.index(max(ks_vals))]
        self.distribution = best_fit
        
    def get_best_fits(self, distribution="Any"):
        
        return [self.get_best_fit(distribution=distribution, bounds=bounds) for bounds in self.bins.data["bounds"]]
        
    def get_best_fit(self, distribution="Current", bounds="Current"):
        
        if distribution == "Any":
            distribution = self.distribution["name"]
        elif distribution == "Current":
            statistic, phase, property_ = self.stats.current["value list"]
            distribution = self.stats.data["Modified"][statistic][phase][property_]["Distribution Type"]
        
        if   bounds == "Current":
            data = self.data["Current"]
        elif bounds == "All":
            data = self.get_data([self.bins.data["bounds"][0][0], self.bins.data["bounds"][-1][-1]])
        else:
            data = self.get_data(bounds)
            
        return self.fit_distribution(data, getattr(stats, self.distribution_map()[distribution]), distribution)
        
class parameters_handler():

    def __init__(self, events, window, file, stats, bins, distribution_data):
        
        self.events            = events
        self.window            = window
        self.canvas            = None
        self.file              = file
        self.stats             = stats
        self.bins              = bins
        self.distribution_data = distribution_data
        
        self.data    = {
            "Original":{},
            "Modified":{}
        }
        self.current = {
            "Original":{},
            "Modified":{}
        }
        
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event ("parameters cycle")
        self.events.add_event ("parameters open" )
        self.events.add_event ("parameters close")
    
        self.events.add_action("quit"            , self.destroy     )
        self.events.add_action("file close"      , self.destroy     )
        self.events.add_action("stats close"     , self.destroy     )
        self.events.add_action("bins close"      , self.destroy     )
        self.events.add_action("bins open"       , self.initialize  )
        self.events.add_action("file export"     , self.reinitialize)
        self.events.add_action("stats cycle"     , self.reinitialize)
        self.events.add_action("bins cycle"      , self.reinitialize)
        self.events.add_action("recursive_fit modify_stats done", self.reinitialize)
        
    def initialize(self):
    
        self.update_data()
        
        self.create_canvas()
        
        self.events.event("parameters open")
        
    def reinitialize(self):
        
        self.update_data()
        self.update_labels()
        
        self.events.event("parameters cycle")
        
    def update_data(self):
        
        statistic, phase, property_ = self.stats.current["value list"]
        bin_ = self.bins.current["bin"]
        
        # get data for all bins
        for version in self.current:
            self.data[version] = copy.deepcopy(self.stats.data[version][statistic][phase][property_])
            
        # get data for current bin
        self.current = copy.deepcopy(self.data)
        for version in self.current:
            for parameter in self.current[version]:
                # there should be a distribution type, a bins number, and two parameters in each property
                # we only want the parameters
                if parameter in ["Distribution Type", "bins"]:
                    continue
                # only apply current bin for parameters with multiple values
                if str(type(self.current[version][parameter])) in ["<class 'list'>", "<class 'tuple'>", "<class 'numpy.ndarray'>"]:
                    if len(self.current[version][parameter]) > 1:
                        value = self.current[version][parameter][bin_]
                    else:
                        value = self.current[version][parameter][0]
                else:
                    value = self.current[version][parameter]
                self.current[version][parameter] = value
        
    def create_canvas(self):
            
        self.canvas = canvas_handler(self.window)
        
        ni = len([*self.current.keys()])
        
        # distribution edit label
        label = {
            "text"    : tk.StringVar(),
            "location": [0,0],
            "span"    : [1,ni*2+1],
            "destroy" : False
        }
        label["text"].set("Original/Modified Best Fit Distributions")
        label["widget"] = \
            tk.Label (
                self.canvas.canvas,
                textvariable=label["text"],
                font='Helvetica 12 bold',
                width = 2*self.window["label width"]
            )
        self.canvas.add("label title ", label)
    
        for i, (version, data) in enumerate(zip([*self.current.keys()],[*self.current.values()])):
        
            # titles
            label = {
                "location": [1,0+i*2],
                "span"    : [1,2],
                "destroy" : False,
                "widget" : tk.Label (
                    self.canvas.canvas,
                    text=version,
                    width=self.window["label width"]
                )
            }
            self.canvas.add(f"label title {version}", label)
            
            # distribution title
            label = {
                "location": [2,0+i*2],
                "span"    : [1,1],
                "destroy" : False,
                "widget" : tk.Label (
                    self.canvas.canvas,
                    text="Distribution Type:",
                    width=self.window["label width"]
                )
            }
            self.canvas.add(f"label title {version} distribution", label)
            
            # bestfit title
            label = {
                "location": [1, 2+ni*2],
                "span"    : [1, 1    ],
                "destroy" : False,
                "widget" : tk.Label (
                    self.canvas.canvas,
                    text="Best Fit",
                    width=self.window["label width"]
                )
            }
            self.canvas.add(f"label title bestfit", label)
        
            # best fit current button
            button = {
                "widget": tk.Button(
                    self.canvas.canvas,
                    text = "Current Bin",
                    command=lambda: self.best_fit("Current"),
                    width=self.window["button width"]
                ),
                "location": [2, 2+ni*2],
                "span"    : [1, 1    ],
                "destroy" : False
            }
            self.canvas.add(f"button bestfit current", button)
            
            # best fit all button
            button = {
                "widget": tk.Button(
                    self.canvas.canvas,
                    text = "All Bins",
                    command=lambda: self.best_fit("All"),
                    width=self.window["button width"]
                ),
                "location": [3, 2+ni*2],
                "span"    : [1, 1    ],
                "destroy" : False
            }
            self.canvas.add(f"button bestfit all", button)
            
        self.update_labels()
        
    def update_labels(self):
    
        # this module is variable length by design.
        # destroy only the variable elements.
        for element in [*self.canvas.elements.keys()]:
            if self.canvas.elements[element]["destroy"]:
                self.canvas.elements[element]["widget"].destroy()
                del self.canvas.elements[element]
    
        # recreate the variable elements
        for i, (version, data) in enumerate(self.current.items()):
            
            # distribution type
            optionmenu = {
                "text"    : tk.StringVar(),
                "location": [2,1+i*2],
                "span"    : [1,1],
                "destroy" : True,
                "export"         : version=="Modified",
                "path components": {
                    "statistic": self.stats.current["path components"]["statistic"],
                    "phase"    : self.stats.current["path components"]["phase"    ],
                    "property ": self.stats.current["path components"]["property" ]
                }
            }
            optionmenu["text"].set(self.current[version]["Distribution Type"])
            optionmenu["text"].trace(
                "w", 
                lambda *args, optionmenu_=optionmenu: self.update_distribution(
                    {
                        "path" :optionmenu_["path components"],
                        "value":optionmenu_["text"].get()
                    }
                )
            )
            optionmenu["widget"] = \
                tk.OptionMenu (
                    self.canvas.canvas,
                    optionmenu["text"],
                    *["Log Normal Distribution", "Beta Distribution", "Power Distribution"]
                )
            if version == [*self.current.keys()][0]: # if original
                optionmenu["widget"].config(state="disabled")
            self.canvas.add(f"optionmenu {version}", optionmenu)
            
            # Distribution Parameters
            for j, (name, value) in enumerate(self.current[version].items(), 3):
            
                if not name == "Distribution Type":
                
                    # parameter label
                    label = {
                        "text"    : tk.StringVar(),
                        "location": [0+j,0+i*2],
                        "span"    : [1,1],
                        "destroy" : True
                    }
                    label["text"].set(name)
                    label["widget"] = \
                        tk.Label (
                            self.canvas.canvas,
                            textvariable=label["text"],
                            width=self.window["label width"]
                        )
                    self.canvas.add(f"label title {version} distribution name {name}", label)
                    
                    # parameter value
                    entry = {
                        "text"           : tk.StringVar(),
                        "location"       : [0+j,1+i*2],
                        "span"           : [1,1],
                        "destroy"        : True,
                        "export"         : version=="Modified",
                        "path"           : self.stats.current["path"]["property"]+"/"+name,
                        "path components": {
                            "statistic": self.stats.current["path components"]["statistic"],
                            "phase"    : self.stats.current["path components"]["phase"    ],
                            "property ": self.stats.current["path components"]["property" ],
                            "parameter": name,
                            "bin"      : self.bins .current["bin"]
                        }
                    }
                    entry["text"].set(value)
                    entry["widget"] = \
                        tk.Entry (
                            self.canvas.canvas,
                            textvariable=entry["text"],
                            width=self.window["label width"]
                        )
                    if version == [*self.current.keys()][0]:
                        #completely disable the original data entry boxes
                        entry["widget"].config(state="disabled")
                    else:
                        #update the modified data on the fly as it's typed
                        entry["text"].trace(
                            "w", 
                            lambda *args, entry_=entry: self.update_parameter(
                                {
                                    "path" :entry_["path components"],
                                    "value":entry_["text"].get()
                                }
                            )
                        )
                        #refresh everything down the pipeline that depends on the modified data
                        #once the enter button is hit in a modified data entry box
                        entry["widget"].bind('<Key-Return>', lambda *args:self.reinitialize())
                    self.canvas.add(f"entry title {version} distribution value {name}", entry)
        
        # submit button
        button = {
            "widget": tk.Button(
                self.canvas.canvas,
                text = "Export Distribution",
                command=self.export,
                width=2*self.window["button width"]
            ),
            "location": [1+j, 0],
            "span"    : [1  , 4],
            "destroy" : True
        }
        self.canvas.add(f"button export", button)
                    
        self.canvas.display()
        
    def update_distribution(self, *args):
        
        statistic, phase, property_ = [*args[0]["path" ].values()]
        value = args[0]["value"]
        
        # write out the new distribution type
        self.stats.data["Modified"][statistic][phase][property_]["Distribution Type"] = value
        
        # write out the new parameters for the new distribution type
        for parameter_current in [*self.stats.data["Modified"][statistic][phase][property_].keys()]:
        
            # if the current value is a distribution parameter value, 
            # replace it with the correct new distribution parameter value
            if parameter_current in sum([*self.stats.parameter_map().values()],[]):
        
                # find out if the parameter is a member of the first or second index,
                # then pop the current one out and replace it with the new one
                for index in [0,1]:
                    index_parameters = [val[index] for val in [*self.stats.parameter_map().values()]]
                    if parameter_current in index_parameters:
                        parameter_new = self.stats.parameter_map()[value][index]
                        self.stats.data["Modified"][statistic][phase][property_][parameter_new    ] = self.stats.data["Modified"][statistic][phase][property_].pop(parameter_current)
                    
            # if the current value is NOT a distribution parameter value, 
            # pop it out and re-append it just to keep the dict in the same order
            else:
                self.stats.data["Modified"][statistic][phase][property_][parameter_current] = self.stats.data["Modified"][statistic][phase][property_].pop(parameter_current)

        # change distribution parameter labels
        self.reinitialize()
        
        # trigger event
        self.events.event("parameters cycle")
        
    def update_parameter(self, *args):

        statistic, phase, property_, parameter = [*args[0]["path" ].values()][:-1]
        bin_ = [*args[0]["path" ].values()][ -1]
        value = args[0]["value"]

        # only apply current bin for parameters with multiple values
        if str(type(self.stats.data["Modified"][statistic][phase][property_][parameter])) in ["<class 'list'>", "<class 'tuple'>", "<class 'numpy.ndarray'>"]:
            if len(self.stats.data["Modified"][statistic][phase][property_][parameter]) > 1:
                type_ = type(self.stats.data["Modified"][statistic][phase][property_][parameter][bin_])
            else:
                type_ = type(self.stats.data["Modified"][statistic][phase][property_][parameter][0])
        else:
            type_ = type(self.stats.data["Modified"][statistic][phase][property_][parameter])

        # keep value type true to original data type
        # if the type is invalid, just don't save it
        if type_ is int:
            try:
                value = int(value)
            except:
                ##### TO DO #####
                
                # Warn value error
                
                ##### TO DO #####
                return
        else:
            try:
                value = np.float32(value)
            except:
                ##### TO DO #####
                
                # Warn value error
                
                ##### TO DO #####
                return
        
        # only apply current bin for parameters with multiple values
        if str(type(self.stats.data["Modified"][statistic][phase][property_][parameter])) in ["<class 'list'>", "<class 'tuple'>", "<class 'numpy.ndarray'>"]:
            if len(self.stats.data["Modified"][statistic][phase][property_][parameter]) > 1:
                self.stats.data["Modified"][statistic][phase][property_][parameter][bin_] = value
            else:
                self.stats.data["Modified"][statistic][phase][property_][parameter][0] = value
        else:
            self.stats.data["Modified"][statistic][phase][property_][parameter] = value
        
        # trigger event
        # self.events.event("parameters cycle")
        
    def best_fit(self, data_range):
        
        statistic, phase, property_ = self.stats.current["value list"]
        bin_ = self.bins.current["bin"]
            
        # if user selects "all bins" and this is a parameter with more than one value
        # iterate through every bin and best fit to the data
        if data_range == "All":
            
            best_fits = self.distribution_data.get_best_fits(distribution="Current")
            if best_fits is None or all(best_fit is None for best_fit in best_fits):
                return
            for bin_, best_fit in enumerate(best_fits):
                if best_fit is None:
                    continue
                for name, value in zip([*best_fit["parameters"].keys()], [*best_fit["parameters"].values()]):
                    if np.array(value) is None:
                        continue
                    self.stats.data["Modified"][statistic][phase][property_][name][bin_] = np.array(value)
                    
        else:
            
            bounds = "Current"
        
            best_fit = self.distribution_data.get_best_fit(distribution="Current", bounds=bounds)
            if best_fit is None:
                return
            for name, value in zip([*best_fit["parameters"].keys()], [*best_fit["parameters"].values()]):
                if np.array(value) is None:
                        continue
                self.stats.data["Modified"][statistic][phase][property_][name][bin_] = np.array(value)
                
        # change distribution parameter labels
        self.reinitialize()
                
        # trigger event
        self.events.event("parameters cycle")
        
    def export(self):
        
        self.file.export(self.stats)
                
    def destroy(self):
    
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("parameters close")
        
class recursive_fit():
    
    def __init__(self, events, window, file, stats, parameters, distribution_data):
    
        self.events            = events
        self.window            = window
        self.canvas            = None
        self.file              = file
        self.stats             = stats
        self.parameters        = parameters
        self.distribution_data = distribution_data
        
        self.paths = {
            "dream3d"      : None,
            "generator"    : {
                "file": {
                    "input" : "./utils_packing_recursive_3d_input.dream3d",
                    "json"  : "./utils_packing_recursive_3d.json",
                    "output": "./utils_packing_recursive_3d_output.dream3d"
                },
                "hdf5": {
                    "FeatureIds": "/DataContainers/SyntheticVolumeDataContainer/CellData/FeatureIds",
                    "Feature"   : "/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/AspectRatios",
                    "Geometry"  : "/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY"
                }
            } ,
            "discriminator": {
                "file": {
                    "input" : "./utils_packing_recursive_2d_input.dream3d",
                    "json"  : "./utils_packing_recursive_2d.json",
                    "output": "./utils_packing_recursive_2d_output.dream3d"
                },
                "hdf5"     : {
                    "FeatureIds": "/DataContainers/ImageDataContainer/CellData/FeatureIds",
                    "Feature"   : "/DataContainers/ImageDataContainer/CellFeatureData/AspectRatios",
                    "Geometry"  : "/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY"
                }
            }
        }
        
        # contains the data required to compare
        # the ebsd feature to the generated feature
        self.feature_data = {
            "Original": {
                "data": None,
                "fit" : None
            },
            "Modified": {
                "data": None,
                "fit" : None
            }
        }
        
        # Reuse features_handler to grab PhaseTyes and PhaseNames
        self.phases = None
        self.current = {
            "required features": {
                "Phase Names": {
                    "index"   : 0,
                    "shape"   : [1,'any'],
                    "datatype": "StringDataArray"
                },
                "Phase Types": {
                    "index"   : 0,
                    "shape"   : [2,1],
                    "datatype": "DataArray<uint32_t>"
                }
            }
        }
        
        #######################################################################################################################################################################
        
        # bounds for parameters
        self.bounds = {
            "mode": [],
            "var" : [],
            "skew": []
        }
        
        #######################################################################################################################################################################
        
        self.create_events()
        
    def create_events(self):
    
        self.events.add_event ("recursive_fit open"             )
        self.events.add_event ("recursive_fit close"            )
        self.events.add_event ("recursive_fit modify_stats done")
        self.events.add_event ("recursive_fit generate done"    )
        self.events.add_event ("recursive_fit discriminate done")
    
        self.events.add_action("quit"            , self.destroy     )
        self.events.add_action("file close"      , self.destroy     )
        self.events.add_action("parameters close", self.destroy     )
        self.events.add_action("parameters open" , self.initialize  )
        
        #######################################################################################################################################################################
        
        self.events.add_action("recursive_fit discriminate done", self.modify_stats)
        self.events.add_action("recursive_fit modify_stats done", self.generate    )
        
        #######################################################################################################################################################################
    
        self.events.add_action("recursive_fit generate done"    , self.discriminate)
        
    def initialize(self):
    
        self.create_canvas()
        
        # Reuse features_handler to grab PhaseTyes and PhaseNames
        # features_handler(self.events, self.window, self.file, self)
        
        self.events.event("recursive_fit open")
        
    def create_canvas(self):
            
        self.canvas = canvas_handler(self.window)
        
        label = {
            "text"    : "Recursive Fit",
            "location": [0, 0],
            "span"    : [1, 2]
        }
        label["widget"] = \
            tk.Label (
                self.canvas.canvas,
                text=label["text"],
                font='Helvetica 12 bold',
                width=self.window["label width"]
            )
        self.canvas.add("label title ", label)
        
        button = {
            "widget": tk.Button(
                self.canvas.canvas,
                text = "Begin",
                command=self.recursive_fit,
                width=self.window["button width"],
                height=self.window["button height"]
            ),
            "location": [1, 1],
            "span"    : [1, 1],
        }
        self.canvas.add(f"button export", button)
        
        # somethnig to grab the required PhaseNames and PhaseTypes
        locations = {
            "recursive_fit phases": {
                "location"     : [1,0],
                "span"         : [2,1],
                "label width"  : 20,
                "button width" : 2,
                "button height": 1
            }
        }
        self.phases = features_handler(self.events, dict({"widget":self.canvas.canvas}, **locations["recursive_fit phases"]), self.file, self, {"self":"recursive_fit phases", "parent":"recursive_fit"})
            
        self.canvas.display()
        
    def recursive_fit(self):
        
        # find dream.3d path
        self.browse_dream3d()
        
        # catch no path selected
        if self.paths["dream3d"] is None:
            return
        
        # initialize files to create first 3d microstructure
        self.create_input_3d()
        self.create_json_3d()
        
        # get fitment of original data
        self.feature_data["Original"]["data"] = self.distribution_data.data["All"]
        self.feature_data["Original"]["fit" ] = self.get_fit(self.feature_data["Original"]["data"])
        
        # run initial generator pass
        self.generate()
        
    #######################################################################################################################################################################
    
    def modify_stats(self):

        statistic, phase, property_ = self.stats.current["value list"]
        distribution_type = self.stats.data["Modified"][statistic][phase][property_]["Distribution Type"]

        if distribution_type == "Beta Distribution":
            
            a0 = self.feature_data["Original"]["fit"]["parameters"]["Alpha"]
            b0  = self.feature_data["Original"]["fit"]["parameters"]["Beta" ]
            mode0, var0, skew0 = self.get_shape_data(a0, b0, "Beta Distribution")
            # shape0 = np.array([mode0, var0, skew0])

            a1 = self.feature_data["Modified"]["fit"]["parameters"]["Alpha"]
            b1  = self.feature_data["Modified"]["fit"]["parameters"]["Beta" ]
            mode1, var1, skew1 = self.get_shape_data(a1, b1, "Beta Distribution")
            # shape1 = np.array([mode1, var1, skew1])
            
            if self.bounds["mode"] == []:
                
                dist_mode = mode1-mode0
                dist_var  = var1 -var0
                dist_skew = skew1-skew0
                
                c = 3
                self.bounds["mode"] = [mode0-max(0, c*dist_mode), mode0-min(0, c*dist_mode)]
                self.bounds["var" ] = [var0 -max(0, c*dist_var ), var0 -min(0, c*dist_var )]
                self.bounds["skew"] = [skew0-max(0, c*dist_skew), skew0-min(0, c*dist_skew)]
                
                # self.bounds["mode"] = [ 0.01, 0.99]
                # self.bounds["var" ] = [ 0.01, 0.99]
                # self.bounds["skew"] = [-0.99, 0.99]
            
            mode2 = sum(self.bounds["mode"])/2
            var2  = sum(self.bounds["var" ])/2
            skew2 = sum(self.bounds["skew"])/2
            
            if mode1 > mode0:
                self.bounds["mode"][1] = mode2
            else:
                self.bounds["mode"][0] = mode2
            if var1  > var0 :
                self.bounds["var" ][1] = var2
            else:
                self.bounds["var" ][0] = var2
            if skew1 > skew0:
                self.bounds["skew"][1] = skew2
            else:
                self.bounds["skew"][0] = skew2
            
            params = self.get_beta_parameters_given_mean(mode2, var2, skew2, a0, b0)
            a2, b2 = params[-1]
            
            # print(self.bounds)
            # print("mode  : ",mode0,", ",mode1)
            # print("var   : ",var0 ,", ",var1 )
            # print("skew  : ",skew0,", ",skew1)
            # print("params: ",a2   ,", ",b2   )
            
            # input()
            
            self.stats.data["Modified"][statistic][phase][property_]["Alpha"] = np.array([a2])
            self.stats.data["Modified"][statistic][phase][property_]["Beta" ] = np.array([b2])
            
        else:
            
            return
        
        self.events.event("recursive_fit modify_stats done")
        
        self.create_input_3d()
        
    ###############################################################################################################################################################
        
    def generate(self):
        
        path_json   = os.path.realpath(self.paths["generator"]["file"]["json"  ])
        path_output = os.path.realpath(self.paths["generator"]["file"]["output"])
        
        # call dream.3d to create 3d microstructure
        exit_code = utils_dream3d.call_dream3d(self.paths["dream3d"], path_json)
        
        # catch dream.3d pipeline failure
        if exit_code > 0:
        
            ###### TO DO ######
            # throw pipeline failure warning
            # make a tooltip for when this fails explaining failure prevention
            #     Ensure PhaseTypes is correctly chosen
            #     "Find Biased Features (BoundingBox)" filter always fails.
            #     Solution: Make a "Replace Value in Array" filter that changes values in "BiasedFeatures" from 1 to 0
            #     Make sure there are enough grains
            #     Make sure there are no small grains (<~9 voxels)
            ###### TO DO ######
            
            print(exit_code)
            return
        
        # # delete old datacontainer
        # with h5py.File(path_output, 'r+') as file_output:
            # del file_output["/".join(self.stats.current["path"]["statistic"].split("/",3)[:-1])]
            
        # # make xdmf file for paraview
        # utils_dream3d.make_xdmf(path_output)
        
        # flag next stage to start
        self.events.event("recursive_fit generate done")
        
    def discriminate(self):
        
        self.feature_data["Modified"] = {}
        
        # feature data from latest iteration
        # slice into 2d sections and gather cumulative feature distribution
        self.feature_data["Modified"]["data"] = self.get_data_2d()[:,0]
        
        ##### TEMP #####
        
        # plot the resultant 2d histogram
        plt.figure()
        
        # plot original data
        bins = np.linspace(self.feature_data["Original"]["data"].min(), self.feature_data["Original"]["data"].max(), 50)
        weights = np.ones_like(self.feature_data["Original"]["data"])/len(self.feature_data["Original"]["data"])
        plt.hist(self.feature_data["Original"]["data"], weights=weights, bins=bins, density=True, label="Original", alpha=0.6)
        
        # params = beta.fit(self.feature_data["Original"]["data"], floc=0)
        # frv = beta(*params)
        # x1 = np.linspace(frv.ppf(0.01), frv.ppf(0.99), 100)
        # y1 = frv.cdf(x1)
        # plt.plot(x1,y1, label="Original CDF")
        # print("original: ", params)
        
        # plot modified data
        bins = np.linspace(self.feature_data["Modified"]["data"].min(), self.feature_data["Modified"]["data"].max(), 50)
        weights = np.ones_like(self.feature_data["Modified"]["data"])/len(self.feature_data["Modified"]["data"])
        plt.hist(self.feature_data["Modified"]["data"], weights=weights, bins=bins, density=True, label="Modified", alpha=0.6)
        
        # params = beta.fit(self.feature_data["Modified"]["data"], floc=0)
        # frv = beta(*params)
        # x2 = np.linspace(frv.ppf(0.01), frv.ppf(0.99), 100)
        # y2 = frv.cdf(x2)
        # plt.plot(x2,y2, label="Modified CDF")
        # print("modified: ", params)
        
        plt.xlim([0,2])
        
        plt.legend()
        
        plt.show()
        
        ##### TEMP #####
        
        # # required for actually doing this correctly. for now, let's just get it working
        # required_features = self.stats.current["required features"]
        # distribution_modified = self.distribution_data.get_data(bounds, data_source=features_modified)
        
        # compare with ks test
        ks = ks_2samp(self.feature_data["Original"]["data"], self.feature_data["Modified"]["data"], alternative='two-sided')
        
        # find the feature fitment data for the latest iteration
        self.feature_data["Modified"]["fit"] = self.get_fit(self.feature_data["Modified"]["data"])
        print(ks)
        
        # flag next stage to start
        self.events.event("recursive_fit discriminate done")
            
    def get_shape_data(self,a,b,distribution_type):
        
        if distribution_type == "Beta Distribution":
            mode = (a-1)/(a+b-2)
            var  = a*b/(((a+b)**2)*(a+b+1))
            skew = 2*(b-a)*(a+b+1)**(1/2)/((a+b+2)*(a*b)**(1/2))
            
        return mode, var, skew
    
    def get_beta_parameters_given_mean(self, mode0, var0, skew0, alpha0, beta0, tolerance=.01, step=0.005):
    
        parameters = []
        for alpha1 in np.arange(0.01, alpha0, step):
            for beta1 in np.arange(0.01, beta0, step):
                mode1, var1, skew1 = self.get_shape_data(alpha1, beta1, "Beta Distribution")
                error_mode = abs(mode0 - mode1)
                error_var  = abs(var0  - var1 )
                error_skew = abs(skew0 - skew1)
                if \
                    error_mode < tolerance\
                    and error_var < tolerance\
                    and error_skew < tolerance:
                    
                    parameters.append([alpha1, beta1])
                    
        if parameters == []:
            print("No valid parameters found, relaxing tolerance")
            parameters = self.get_beta_parameters_given_mean(mode0, var0, skew0, alpha0, beta0, tolerance=tolerance*2, step=step)
        
        return parameters

    def get_fit(self, data):
        
        statistic, phase, property_ = self.stats.current["value list"]
        distribution_name_dream3d   = self.stats.data["Modified"][statistic][phase][property_]["Distribution Type"]
        distribution_name_scipy     = self.distribution_data.distribution_map()[distribution_name_dream3d]
        distribution_obj_scipy      = getattr(stats, distribution_name_scipy)
        
        fit = self.distribution_data.fit_distribution(data, distribution_obj_scipy, distribution_name_dream3d)
        
        return fit
    
    # def features_map(self):
        
        # # DREAM.3D <=> distribution_data
        # return {
            # "EquivalentDiameters": "Equivalent Diameter",
            # "AspectRatios"       : "Aspect Ratio (B/A)" ,
            # "AspectRatios"       : "Aspect Ratio (C/A)" ,
            # "Omega3s"            : "Omega 3"            
        # }
        
    def browse_dream3d(self):
        
        # open file browser
        path = filedialog.askdirectory(
            initialdir = ".",
            title = "Select DREAM.3D Directory"
        )
        
        if path:
        
            # update the file path and import the file data
            self.paths["dream3d"] = path
    
    def get_data_2d(self):
        
        # set up file paths
        path_output_3d = os.path.realpath(self.paths["generator"    ]["file"]["output"])
        path_input_2d  = os.path.realpath(self.paths["discriminator"]["file"]["input" ])
        path_json_2d   = os.path.realpath(self.paths["discriminator"]["file"]["json"  ])
        path_output_2d = os.path.realpath(self.paths["discriminator"]["file"]["output"])
        
        # set up dream.3d paths
        path_in_FeatureIds  = self.paths["generator"    ]["hdf5"]["FeatureIds"]
        path_in_Geometry    = self.paths["generator"    ]["hdf5"]["Geometry"  ]
        path_out_FeatureIds = self.paths["discriminator"]["hdf5"]["FeatureIds"]
        path_out_Feature    = self.paths["discriminator"]["hdf5"]["Feature"   ]
        path_out_Geometry   = self.paths["discriminator"]["hdf5"]["Geometry"  ]
        
        # only calculate feature for every "skip"th slice
        skips = 5
            
        def get_Geometry(path_input, path_Geometry):
        
            # gather 3d geometry
            with h5py.File(path_input, "r") as file_input:
                dims        = file_input[path_Geometry  ]["DIMENSIONS"][...].tolist()
                origin      = file_input[path_Geometry  ]["ORIGIN"    ][...].tolist()
                size_voxels = file_input[path_Geometry  ]["SPACING"   ][...].tolist()
                
            # turn the 3d geometry into a 2d geometry
            # must flip because dream3d needs 2d slices to be z-normal
            # and axes are [x,y,z]
            # flip      => [z,y,x]
            # ...[2]=1  => [z,y,1]
            geometry = {
                "dims"       : np.flip(dims       ),
                "origin"     : np.flip(origin     ),
                "size_voxels": np.flip(size_voxels)
            }
            geometry["dims"       ][2] = 1
            geometry["size_voxels"][2] = 1
            
            return geometry
            
        def get_FeatureIds(path_input, path_FeatureIds):
            
            with h5py.File(path_input, "r") as file_input:
                FeatureIds  = file_input[path_FeatureIds][...]
                
            return FeatureIds
            
        Geometry   = get_Geometry  (path_output_3d, path_in_Geometry  )
        FeatureIds = get_FeatureIds(path_output_3d, path_in_FeatureIds)
        
        # create the json script that finds the
        # feature data for the 2d slices
        self.create_json_2d()
        
        # find feature data for each slice of the 3d microstructure
        for i in range(FeatureIds.shape[0]):
            
            # skip some slices for faster runs
            if i%skips > 0:
                continue

            ### gather 2d slices of FeatureIds
            # dream3d needs 2d slices to be z-normal
            # and axes are [z,y,x,components]
            # [:,:,i,:] => [z,y,1,maybe_components]
            # swapaxes  => [1,y,z,maybe_components]
            # reshape   => [1,y,z,components]
            # this perfectly matches the [z,y,1] of geometry
            # because the two lists are opposites (flipped)
            slice_ = FeatureIds[:,:,i,:]
            slice_ = slice_.swapaxes(0,2)
            if not FeatureIds.shape[-1] > 1:
                slice_ = slice_.reshape(tuple(list(slice_.shape)+[1]))
            
            # push FeatureIds to .dream3d file and run ([scalar segment], [find shapes])
            self.create_input_2d(slice_, path_out_FeatureIds, Geometry, path_out_Geometry, path_input_2d)
            utils_dream3d.update_expected(path_json_2d)
            utils_dream3d.call_dream3d(self.paths["dream3d"], path_json_2d)
            utils_dream3d.make_xdmf(path_input_2d)
            
            # pull feature result
            with h5py.File(path_output_2d, 'r') as file_output_2d:
                if i == 0:
                    result = file_output_2d[path_out_Feature][...]
                else:
                    result = np.append(result, file_output_2d[path_out_Feature][...], axis=0)
                    
        return result
    
    def create_input_2d(self, FeatureIds, path_FeatureIds, Geometry, path_Geometry, path_output):
        
        # create a mask that excludes zero values
        mask = copy.deepcopy(FeatureIds)
        mask[mask>0] = 1
        
        # create the dream3d file that will act as
        # an input to the 2d slice feature finder
        utils_dream3d.create_file_dream3d(path_output)
        utils_dream3d.insert_geometry(path_output,Geometry,path_Geometry)
        utils_dream3d.insert_attribute_array(
            path_output           = path_output                     ,
            path_group            = path_FeatureIds.rsplit("/",1)[0],
            name                  = "FeatureIds"                    ,
            data                  = FeatureIds                      ,
            dtype                 = "DataArray<int32_t>"            ,
            attribute_matrix_type = "Cell"
        )
        utils_dream3d.insert_attribute_array(
            path_output           = path_output                     ,
            path_group            = path_FeatureIds.rsplit("/",1)[0],
            name                  = "Mask"                          ,
            data                  = mask                            ,
            dtype                 = "DataArray<bool>"                ,
            attribute_matrix_type = "Cell"
        )
    
    def create_input_3d(self):
        
        # setup file paths
        path_input  = os.path.realpath(self.paths["generator"]["file"]["input"])
        
        # setup hdf5 paths and gather variable names
        path_Statistics       = self.stats.current["path"]["statistic"]
        path_CellEnsembleData = path_Statistics.rsplit("/",1)[0]
        path_Geometry         = path_Statistics.rsplit("/",2)[0]+"/"+"_SIMPL_GEOMETRY"
        path_PhaseTypes       = "/DataContainers/" + "/".join(self.phases.current["Phase Types"]["value list"])
        path_PhaseNames       = "/DataContainers/" + "/".join(self.phases.current["Phase Names"]["value list"])

        print("PhaseTypes: ", path_PhaseTypes)
        print("PhaseNames: ", path_PhaseNames)
        
        # create the new hdf5 file
        utils_dream3d.create_file_dream3d(path_input)
        
        # copy over CellEnsembleData
        # can't use the utils_dream3d version because the reference file is still open by file_handler
        with h5py.File(path_input, 'a') as file_input:
            file_input.copy(self.file.file[path_CellEnsembleData], path_CellEnsembleData)
            # shouldn't be strictly required, but dream3d won't read the datacontainer without it
            file_input.copy(self.file.file[path_Geometry        ], path_Geometry        )

        # copy over the required variables for packing
        with h5py.File(path_input, 'r+') as file_input:
            for path_packing_req in [path_PhaseTypes, path_PhaseNames]:
                path, name = path_packing_req.rsplit('/',1)
                if not path == path_CellEnsembleData:
                    if name in file_input[path_CellEnsembleData]:
                        del file_input[path_CellEnsembleData+'/'+name]
                    file_input.copy(self.file.file[path_packing_req], path_CellEnsembleData+'/'+name)

        # change stats to modified stats
        with h5py.File(path_input, 'r+') as file_input:
            statistic, phase, property_ = self.stats.current["value list"]
            path = self.stats.current["path"]["property"]
            print(self.stats.data["Modified"][statistic][phase][property_])
            for name, value in self.stats.data["Modified"][statistic][phase][property_].items():
                if name in sum([*self.stats.parameter_map().values()],[]):
                    file_input[path+"/"+name][...] = value.reshape((value.size,1))
        
    def create_json_2d(self):
        
        # setup file paths
        path_input  = os.path.realpath(self.paths["discriminator"]["file"]["input" ])
        path_json   = os.path.realpath(self.paths["discriminator"]["file"]["json"  ])
        path_output = os.path.realpath(self.paths["discriminator"]["file"]["output"])
        
        data_json = \
            {
                "0": {
                    "FilterVersion": "1.2.815",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Read DREAM.3D Data File",
                    "Filter_Name": "DataContainerReader",
                    "Filter_Uuid": "{043cbde5-3878-5718-958f-ae75714df0df}",
                    "InputFile": path_input,
                    "InputFileDataContainerArrayProxy": {}
                },
                "1": {
                    "ActiveArrayName": "Active",
                    "CellFeatureAttributeMatrixName": "CellFeatureData",
                    "FeatureIdsArrayName": "FeatureIds_",
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Segment Features (Scalar)",
                    "Filter_Name": "ScalarSegmentFeatures",
                    "Filter_Uuid": "{2c5edebf-95d8-511f-b787-90ee2adf485c}",
                    "GoodVoxelsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "Mask",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "ScalarArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "ScalarTolerance": 1,
                    "UseGoodVoxels": 1
                },
                "2": {
                    "CentroidsArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "Centroids",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds_",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Centroids",
                    "Filter_Name": "FindFeatureCentroids",
                    "Filter_Uuid": "{6f8ca36f-2995-5bd3-8672-6b0b80d5b2ca}"
                },
                "3": {
                    "EquivalentDiametersArrayName": "EquivalentDiameters",
                    "FeatureAttributeMatrixName": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds_",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Sizes",
                    "Filter_Name": "FindSizes",
                    "Filter_Uuid": "{656f144c-a120-5c3b-bee5-06deab438588}",
                    "NumElementsArrayName": "NumElements",
                    "SaveElementSizes": 0,
                    "VolumesArrayName": "Volumes"
                },
                "4": {
                    "AspectRatiosArrayName": "AspectRatios",
                    "AxisEulerAnglesArrayName": "AxisEulerAngles",
                    "AxisLengthsArrayName": "AxisLengths",
                    "CellFeatureAttributeMatrixName": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "CentroidsArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "Centroids",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds_",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Shapes",
                    "Filter_Name": "FindShapes",
                    "Filter_Uuid": "{3b0ababf-9c8d-538d-96af-e40775c4f0ab}",
                    "Omega3sArrayName": "Omega3s",
                    "VolumesArrayName": "Volumes_"
                },
                "5": {
                    "FilterVersion": "1.2.815",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Write DREAM.3D Data File",
                    "Filter_Name": "DataContainerWriter",
                    "Filter_Uuid": "{3fcd4c43-9d75-5b86-aad4-4441bc914f37}",
                    "OutputFile": path_output,
                    "WriteTimeSeries": 0,
                    "WriteXdmfFile": 0
                },
                "PipelineBuilder": {
                    "Name": path_json.rsplit("/",1)[-1].rsplit(".",1)[0],
                    "Number_Filters": 6,
                    "Version": 6
                }
            }
        
        with open(path_json,"w") as f:
            json.dump(data_json, f, indent=4, separators=(",", ": "), sort_keys=True)
    
    def create_json_3d(self):
        
        # setup file paths
        path_input  = os.path.realpath(self.paths["generator"]["file"]["input" ])
        path_json   = os.path.realpath(self.paths["generator"]["file"]["json"  ])
        path_output = os.path.realpath(self.paths["generator"]["file"]["output"])
        
        ################ MAKE THESE INTO WIDGETS #################
        data_dimensions = [50, 50, 50]   ####### TO DO ###########
        data_resolution = [.2, .2, .2]   ####### TO DO ###########
        data_origin     = [ 0,  0,  0]   ####### TO DO ###########
        # these guys determine if you're using ellipsoids, hyper-ellipsoids, etc...
        # currently defaulting to ellipsoids
        data_shapetypes = self.phases.current["Phase Types"]["value"].flatten().tolist()
        data_shapetypes[1:] = np.zeros_like(data_shapetypes[1:]).tolist()
        ##########################################################
        
        with h5py.File(path_input, 'r') as file_input:
        
            data_json = \
            {
                "0": {
                    "FilterVersion": "1.2.815",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Read DREAM.3D Data File",
                    "Filter_Name": "DataContainerReader",
                    "Filter_Uuid": "{043cbde5-3878-5718-958f-ae75714df0df}",
                    "InputFile": path_input,
                    "InputFileDataContainerArrayProxy": utils_dream3d.create_expected(file_input)["InputFileDataContainerArrayProxy"],
                    "OverwriteExistingDataContainers": 0
                },
                "1": {
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Establish Shape Types",
                    "Filter_Name": "EstablishShapeTypes",
                    "Filter_Uuid": "{4edbbd35-a96b-5ff1-984a-153d733e2abb}",
                    "InputPhaseTypesArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.phases.current["Phase Types"]["path components"]["dataset"],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "ShapeTypeData": data_shapetypes,
                    "ShapeTypesArrayName": "ShapeTypes"
                },
                "2": {
                    "BoxDimensions": "X Range: 0 to 32 (Delta: 32)\nY Range: 0 to 32 (Delta: 32)\nZ Range: 0 to 32 (Delta: 32)",
                    "CellAttributeMatrixName": "CellData",
                    "DataContainerName": "SyntheticVolumeDataContainer",
                    "Dimensions": {
                        "x": data_dimensions[0],
                        "y": data_dimensions[1],
                        "z": data_dimensions[2]
                    },
                    "EnsembleAttributeMatrixName": "CellEnsembleData",
                    "EstimateNumberOfFeatures": 0,
                    "EstimatedPrimaryFeatures": "",
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Initialize Synthetic Volume",
                    "Filter_Name": "InitializeSyntheticVolume",
                    "Filter_Uuid": "{c2ae366b-251f-5dbd-9d70-d790376c0c0d}",
                    "InputPhaseTypesArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.phases.current["Phase Types"]["path components"]["dataset"],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "InputStatsArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.stats.current["path"]["statistic"].split("/")[4],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "Origin": {
                        "x": data_origin[0],
                        "y": data_origin[1],
                        "z": data_origin[2]
                    },
                    "Resolution": {
                        "x": data_resolution[0],
                        "y": data_resolution[1],
                        "z": data_resolution[2]
                    }
                },
                "3": {
                    "CellPhasesArrayName": "Phases",
                    "FeatureGeneration": 0,
                    "FeatureIdsArrayName": "FeatureIds",
                    "FeatureInputFile": path_input,
                    "FeaturePhasesArrayName": "Phases",
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Pack Primary Phases",
                    "Filter_Name": "PackPrimaryPhases",
                    "Filter_Uuid": "{84305312-0d10-50ca-b89a-fda17a353cc9}",
                    "InputPhaseNamesArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.phases.current["Phase Names"]["path components"]["dataset"],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "InputPhaseTypesArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.phases.current["Phase Types"]["path components"]["dataset"],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "InputShapeTypesArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": "ShapeTypes",
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "InputStatsArrayPath": {
                        "Attribute Matrix Name": self.stats.current["path"]["statistic"].split("/")[3],
                        "Data Array Name": self.stats.current["path"]["statistic"].split("/")[4],
                        "Data Container Name": self.stats.current["path"]["statistic"].split("/")[2]
                    },
                    "MaskArrayPath": {
                        "Attribute Matrix Name": "",
                        "Data Array Name": "",
                        "Data Container Name": ""
                    },
                    "NewAttributeMatrixPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "NumFeaturesArrayName": "NumFeatures",
                    "OutputCellAttributeMatrixPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "OutputCellEnsembleAttributeMatrixName": "CellEnsembleData",
                    "OutputCellFeatureAttributeMatrixName": "CellFeatureData",
                    "PeriodicBoundaries": 0,
                    "SaveGeometricDescriptions": 0,
                    "SelectedAttributeMatrixPath": {
                        "Attribute Matrix Name": "",
                        "Data Array Name": "",
                        "Data Container Name": ""
                    },
                    "UseMask": 0
                },
                "4": {
                    "CentroidsArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "Centroids",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Centroids",
                    "Filter_Name": "FindFeatureCentroids",
                    "Filter_Uuid": "{6f8ca36f-2995-5bd3-8672-6b0b80d5b2ca}"
                },
                "5": {
                    "AspectRatiosArrayName": "AspectRatios",
                    "AxisEulerAnglesArrayName": "AxisEulerAngles",
                    "AxisLengthsArrayName": "AxisLengths",
                    "CellFeatureAttributeMatrixName": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "CentroidsArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "Centroids",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Shapes",
                    "Filter_Name": "FindShapes",
                    "Filter_Uuid": "{3b0ababf-9c8d-538d-96af-e40775c4f0ab}",
                    "Omega3sArrayName": "Omega3s",
                    "VolumesArrayName": "Volumes"
                },
                "6": {
                    "EquivalentDiametersArrayName": "EquivalentDiameters",
                    "FeatureAttributeMatrixName": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Find Feature Sizes",
                    "Filter_Name": "FindSizes",
                    "Filter_Uuid": "{656f144c-a120-5c3b-bee5-06deab438588}",
                    "NumElementsArrayName": "NumElements",
                    "SaveElementSizes": 0,
                    "VolumesArrayName": "Volumes_"
                },
                "7": {
                    "ApplyToSinglePhase": 0,
                    "FeatureIdsArrayPath": {
                        "Attribute Matrix Name": "CellData",
                        "Data Array Name": "FeatureIds",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "FeaturePhasesArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "Phases",
                        "Data Container Name": "ImageDataContainer"
                    },
                    "FilterVersion": "6.5.141",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Minimum Size",
                    "Filter_Name": "MinSize",
                    "Filter_Uuid": "{53ac1638-8934-57b8-b8e5-4b91cdda23ec}",
                    "IgnoredDataArrayPaths": [
                    ],
                    "MinAllowedFeatureSize": 27,
                    "NumCellsArrayPath": {
                        "Attribute Matrix Name": "CellFeatureData",
                        "Data Array Name": "NumElements",
                        "Data Container Name": "SyntheticVolumeDataContainer"
                    },
                    "PhaseNumber": 0
                },
                "8": {
                    "FilterVersion": "1.2.815",
                    "Filter_Enabled": True,
                    "Filter_Human_Label": "Write DREAM.3D Data File",
                    "Filter_Name": "DataContainerWriter",
                    "Filter_Uuid": "{3fcd4c43-9d75-5b86-aad4-4441bc914f37}",
                    "OutputFile": path_output,
                    "WriteTimeSeries": 0,
                    "WriteXdmfFile": 0
                },
                "PipelineBuilder": {
                    "Name": path_json.rsplit("/",1)[-1].rsplit(".",1)[0],
                    "Number_Filters": 9,
                    "Version": 6
                }
            }

        with open(path_json,"w") as f:
            json.dump(data_json, f, indent=4, separators=(",", ": "), sort_keys=True)
        
    def destroy(self):
        
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("recursive_fit close")
        
class figures_handler():

    """
    # Example Data Structure
    data = {
        "[title1]":{
            "type": "hist",
            "data": {
                "bins": int(),
                "vals": np.array() # 1xn (flattened)
            }
        },
        "[title2]":{
            "type": "plot",
            "data": {
                "x" : np.array(), # 1xn (flattened)
                "y" : np.array()  # 1xn (flattened)
            }
        }
    }
    """

    def __init__(self, events, window):
        
        self.events = events
        self.window = window
        self.canvas = None
        self.figures   = {}
        
        self.create_events()
        
    def create_events(self):
        
        self.events.add_event ("figures open" )
        self.events.add_event ("figures close")
        
        self.events.add_action("quit"            , self.destroy     )
        self.events.add_action("file close"      , self.destroy     )
        self.events.add_action("parameters close", self.destroy     )
        self.events.add_action("parameters open" , self.initialize  )
        
    def initialize(self):
        
        self.create_canvas()
        
        self.events.event("figures open")
        
    def create_canvas(self):
        
        self.canvas = canvas_handler(self.window)
        
        self.canvas.display()
        
    def add_plots(self, name, plots):
    
        figure_size = (3,3)
        figure_dpi  = 100
        plot_alpha  = 0.6
        
        if name in [*self.canvas.elements.keys()]:
            location = self.canvas.elements[name]["location"][1]
            self.canvas.kill(name)
        
        else:
            location = len([*self.canvas.elements.keys()])
    
        # create the figure that will contain the plots
        figure  = plt.figure(figsize=figure_size, dpi = figure_dpi)
        axis    = figure.add_subplot()
        plt.title(name)
        
        # add the plot (as a child of above canvas) to the widgets list
        plot = {
            "figure"  : figure,
            "axis"    : axis,
            "location": [0, location],
            "span"    : [1,        1],
            "widget"  : FigureCanvasTkAgg(figure,self.canvas.canvas).get_tk_widget()
        }
        self.canvas.add(name, plot)
        
        # Plot
        for i, (title, type_, data) in \
            enumerate(
                zip(
                    [*plots.keys()],
                    [val["type"] for val in [*plots.values()]],
                    [val["data"] for val in [*plots.values()]]
                )
            ):
        
            if   type_ == "plot":
            
                # plot x,y data
                axis.plot(data["x"], data["y"], label=title, alpha=plot_alpha)
            
            elif type_ == "hist":
            
                if \
                    (
                        str(type(data["vals"])) in ["<class 'list'>", "<class 'tuple'>"] \
                        and len(data["vals"]) > 0 \
                    ) or \
                    (
                        str(type(data["vals"])) in ["<class 'numpy.ndarray'>"] \
                        and data["vals"].size > 1 \
                    ):
                    
                    # normalize histogram
                    HIST_BINS = np.linspace(data["vals"].min(), data["vals"].max(), data["bins"])
                    weights = np.ones_like(data["vals"])/len(data["vals"])
                    
                    # plot histogram data
                    axis.hist(x=data["vals"], bins=HIST_BINS, weights=weights, density=True, label=title, alpha=plot_alpha)
            
        axis.legend()
        
        self.canvas.display()

    def destroy(self):
    
        for fig in [*self.figures.values()]:
            plt.close(fig["figure"])
    
        if not self.canvas is None:
            self.canvas.destroy()
            
        self.events.event("figures close")
            
class plot_data():
    
    def __init__(self, events, file, stats, features, bins, distribution_data, parameters, figures):
        
        self.events            = events
        self.file              = file
        self.stats             = stats
        self.features          = features
        self.bins              = bins
        self.distribution_data = distribution_data
        self.parameters        = parameters
        self.figures           = figures
        
        self.create_events()
        
    def create_events(self):
        
        self.events.add_action("figures open"    , self.run  )
        self.events.add_action("stats cycle"     , self.run  )
        self.events.add_action("features cycle"  , self.run  )
        self.events.add_action("bins cycle"      , self.run  )
        self.events.add_action("parameters cycle", self.run  )
        
    def run(self):
        
        for range_ in ["All", "Current"]:

            data = {}
            
            # Plot the data as histograms
            data["Feature Data"] = {
                "type": "hist",
                "data": {
                    "bins": self.parameters.current["Modified"]["bins"],
                    "vals": self.distribution_data.data[range_]
                }
            }
            
            # Plot the data as a line
            for version in [*self.stats.data.keys()]:
                x, y = self.get_distribution(version, range_)
                data[version] = {
                    "type": "plot",
                    "data": {
                        "x" : x,
                        "y" : y
                    }
                }
            
            self.figures.add_plots(range_, data)
        
    def get_distribution(self, version, range_):
        
        statistic, phase, property_ = self.stats.current["value list"]
        properties = self.stats.data[version][statistic][phase][property_]
        bin_ = self.bins.current["bin"]
        
        if   properties["Distribution Type"] == "Beta Distribution":
            if   range_ == "Current":
                a, b = (properties["Alpha"][bin_],properties["Beta"][bin_])
            elif range_ == "All":
                a, b = (properties["Alpha"]      ,properties["Beta"]      )
            x = np.linspace(beta.ppf(0.01, a, b),beta.ppf(0.99, a, b), 100)
            y = beta.pdf(x, a, b)
        elif properties["Distribution Type"] == "Log Normal Distribution":
            if   range_ == "Current":
                u, s = (properties["Average"][bin_],properties["Standard Deviation"][bin_])
            elif range_ == "All":
                u, s = (properties["Average"]      ,properties["Standard Deviation"]      )
            x = np.linspace(lognormal.ppf(0.01, s),lognormal.ppf(0.99, s), 100)
            y = lognormal.pdf(x-u, s)
        elif properties["Distribution Type"] == "Power Distribution":
            ################## TO DO ###################
            
            pass
            
            ################## TO DO ###################
        else:
            ################## TO DO ###################
            
            # throw distribution type error
            pass
            
            ################## TO DO ###################
        
        return x, y
        

###################### run ######################

if __name__ == '__main__':
    utils_packing()
