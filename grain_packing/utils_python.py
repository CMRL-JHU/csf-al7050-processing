#pip install pickle
import pickle
import os
import numpy as np

def dump(obj):
    for attr in dir(obj):
        if hasattr( obj, attr ):
            print( "obj.%s = %s\n" % (attr, getattr(obj, attr)))
    
def pretty_file_list(data, i=1, padding="   "):
    
    # create the structure
    string = ""
    for val in data:
        if type(val) is dict:
            string +=                                                \
                i*padding+"{"+"\n"+                                  \
                pretty_file_dict(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"}"+","+"\n"
        elif type(val) is list:
            string +=                                                \
                i*padding+"["+"\n"+                                  \
                pretty_file_list(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"]"+","+"\n"
        elif type(val) is tuple:
            string +=                                                \
                i*padding+"("+"\n"+                                  \
                pretty_file_tuple(val, i=i+1, padding=padding)+"\n"+ \
                i*padding+")"+","+"\n"
        #elif type(val).__module__ == np.__name__:
        elif str(type(val)) in ["<class 'numpy.ndarray'>"]:
            val = val.tolist()
            string +=                                                \
                i*padding+"["+"\n"+                                  \
                pretty_file_list(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"]"+","+"\n"
        else:
            string += i*padding
            if type(val) is str:
                string += "\""+val+"\""+","
            else:
                string += str(val)+","
            string += "\n"
                
    # remove unnecessary delimiters (,\n)
    string = string[:-2]
    
    return string
    
def pretty_file_tuple(data, i=1, padding="   "):    
    return pretty_file_list(data, i=i+1, padding=padding)
    
def pretty_file_dict(data, i=1, padding="   "):

    # create the structure
    string = ""
    for key,val in zip(data.keys(),data.values()):
        if type(val) is dict:
            string +=                                                \
                i*padding+"\""+key+"\""+": "+"{"+"\n"+               \
                pretty_file_dict(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"}"+","+"\n"
        elif type(val) is list:
            string +=                                                \
                i*padding+"\""+key+"\""+": "+"["+"\n"+               \
                pretty_file_list(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"]"+","+"\n"
        elif type(val) is tuple:
            string +=                                                \
                i*padding+"\""+key+"\""+": "+"("+"\n"+               \
                pretty_file_tuple(val, i=i+1, padding=padding)+"\n"+ \
                i*padding+")"+","+"\n"
        #elif type(val).__module__ == np.__name__:
        elif str(type(val)) in ["<class 'numpy.ndarray'>"]:
            val = val.tolist()
            string +=                                                \
                i*padding+"\""+key+"\""+": "+"["+"\n"+               \
                pretty_file_list(val, i=i+1, padding=padding)+"\n"+  \
                i*padding+"]"+","+"\n"
        else:
            string += i*padding+"\""+key+"\""+": "
            if type(val) is str:
                string += "\""+val+"\""+","
            else:
                string += str(val)+","
            string += "\n"
                
    # remove unnecessary delimiters (,\n)
    string = string[:-2]
    
    return string
        
def data_load(path_file):

    if not os.path.exists(path_file):
        return
    
    with open(path_file, "rb") as file:
        return pickle.load(file)

def data_save(data, path_file):

    os.makedirs(os.path.split(path_file)[0], exist_ok=True)
    with open(path_file, "wb") as file:
        pickle.dump(data, file)
    
def data_copy(root_src, root_dst):
    import shutil, os
    
    if os.path.isdir(root_src):
        for src_dir, dirs, files in os.walk(root_src):
            dst_dir = src_dir.replace(root_src, root_dst, 1)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            for file_ in files:
                src_file = os.path.join(src_dir, file_)
                dst_file = os.path.join(dst_dir, file_)
                if os.path.exists(dst_file):
                    # in case of the src and dst are the same file
                    if os.path.samefile(src_file, dst_file):
                        continue
                    os.remove(dst_file)
                shutil.copy(src_file, dst_dir)
                
    elif os.path.isfile(root_src):
        if not os.path.exists(root_dst):
            os.makedirs(root_dst)
        shutil.copy(root_src, root_dst)
        
# find the shortest path segment that uniquely identifies
# each of n paths from one another
# ex:
#    paths = [
#       "1/2/3/a/b/c"
#       "1/2/4/b/c/a"
#       "1/2/3/c/a/b"
#    ]
#    result = [
#       "3/a"
#       "4"
#       "3/c"
#    ]
def find_shortest_unique_descriptor(paths, delimiter=os.sep, start=True):
    
    # if the solution is already unique, just return the solution
    if len(paths) == 1:
        return [paths[0].split(delimiter)[0]]
    
    # otherwise, recursively find the shortest unique solution
    solution = [None]*len(paths)
    segments = [ path.split(delimiter) for path in paths ]

    for idx_dir, dirs in enumerate( zip(*segments) ):

        dirs_unique = np.unique(dirs)

        # if the current directory level for each path has branching names
        if len(dirs_unique) > 1:

            # split the recursion into separate branches for each name
            for dir_unique in dirs_unique:
                
                # find all the path indecies in this branch
                indecies = [i for i, dir in enumerate(dirs) if dir==dir_unique]
                
                # reconstruct the path with only the the indicies corresponding to this name
                # and starting from the current directory level
                paths = [ delimiter.join(segments[idx_path][idx_dir:]) for idx_path in indecies ]
                
                # recursively find the solution for each branch
                for idx_path, descriptor in zip(indecies, find_shortest_unique_descriptor(paths, delimiter=delimiter, start=False)):
                    head = ( f"{delimiter.join(segments[idx_path][:idx_dir])}{delimiter}", "")[start]
                    solution[idx_path] = f"{head}{descriptor}"
            
            break
    
    return solution

def permute(list_nested):

    # if the list is empty, then we've come to the end
    # of the first list layer and it's time to start
    # going back up the chain
    # ex:
    # [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    # [[4, 5, 6], [7, 8, 9]]
    # [[7, 8, 9]]
    # []                                 <--- we are here
    if not list_nested:
        return

    list_permuted = []

    for x in list_nested[0]:

        vals = permute(list_nested[1:])

        # if we were just at the end of the chain,
        # we only have one set of values to return
        # ex:
        # [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        # [[4, 5, 6], [7, 8, 9]]
        # [[7, 8, 9]]                       <--- we are here
        # []
        if vals is None:
            list_permuted.append([x])

        # for every other part of the chain,
        # append the permutations of the rest
        # of the chain to the current set of values
        else:
            for val in vals:
                list_permuted.append([x, *val])

    return list_permuted

def is_numeric(x):
        
    # if it's a numpy object, unpack it to a standard python object
    if hasattr(x, 'item'): x = x.item()

    return isinstance(x, (int, float, complex)) and not isinstance(x, bool)

def is_list_like(obj):
    return hasattr(obj, '__len__') and not isinstance(obj, (str, dict))
    
def fill_screen(text="", char=1, pad=" || "):
    if   char == 1:
        char = "_.-'""`-._"
    elif char == 2:
        char = "~*:._.:*"

    # if full with is printed, a new blank line is added
    try:
        width_terminal = os.get_terminal_size()[0]-1
    except:
        # fall back to tty width
        width_terminal = 80-1
    
    if type(text) is str:
        text = [text]
        
    if type(text) in [list, tuple]:
    
        # fill each blurb with padding
        text = [pad+blurb+pad for blurb in text if len(blurb) > 0]
        
        # divide the screen width evenly between filler/text/filler/text/filler, etc...
        fill_length = (width_terminal-len("".join(text)))//(len(text)+1)
        
        # catch width lower than string size
        if fill_length < 0:
            return "".join([blurb for blurb in text])
        
        # create the filler string
        ## add all the full filler strings we can
        fill_string = (fill_length)//len(char)*char
        ## add whatever we can of the filler string to fill in the remaining space
        fill_string += char[0:(fill_length)%len(char)]
        
        # join filler strings and text blurbs together
        string = fill_string+"".join([blurb+fill_string for blurb in text])
        
        # add any missing filler characters from int division
        ## add all the full filler strings we can (probably never used)
        string += (width_terminal-len(string))//len(char)*char
        ## add whatever we can of the filler string to fill in the remaining space
        string += char[0:(width_terminal-len(string))%len(char)]
        
        return string