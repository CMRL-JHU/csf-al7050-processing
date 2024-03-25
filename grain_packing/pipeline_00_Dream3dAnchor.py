# Dream3D scripts require absolute paths in order to run,
# so if the scripts are moved or one of the containing directories is moved
# the scripts stop working.
# For large pipelines, this is a huge problem. Not only does the operator have
# to manually redeclare the i/o paths in each script, but certain elements of the
# scripts forget their state (e.g. threshold instructions, invert checkboxes, etc...)
# and information can be permanently lost. 
# This script attempts to fix that problem by going through each dream3d file in
# the current directory and anchoring every path string relative  to the current
# directory name or an explicitly given dir_name_old if the directory has been renamed
# and the old directory name is known.

import os, json, warnings, shlex

dir_name_old = "grain_packing"

def is_path_like(string):
    def is_dos_path_like(string):
        return len(list(filter(lambda x: x in string, ["\\\\","\\"]))) > 0 and string.startswith(":",1)
    def is_unix_path_like(string):
        return len(list(filter(lambda x: x in string, ["/"]))) > 0
    return  is_dos_path_like(string) or is_unix_path_like(string)
    
def replace_separators(string):
    separators = ["\\\\","\\","/"]
    for separator in separators:
        string = string.replace(separator,os.sep)
    return string
    
def replace_path_parent(file_path,dir_current,dir_name_old=None):


    dir_parent, dir_name = os.path.split(dir_current)
    file_path = replace_separators(file_path)

    # find the current directory name within the path string
    # then find the index to the first character after that
    # so we're 'anchoring' relative to the current directory name.
    if dir_name_old is None:

        index = file_path.rfind(os.sep+dir_name+os.sep)

        if index >= 0:
            index += len(dir_name)+2*len(os.sep)

    # we cannot anchor to the current directory name if the
    # directory name was changed.
    # so we're given the old directory name explicitly,
    # we can anchor to that instead.
    # additionally, we take whichever name appears later
    # in the string such that dir_name_old doesn't have
    # to be explicitly set to None in order to use the new path.
    else:

        index_new = file_path.rfind(os.sep+dir_name    +os.sep)
        index_old = file_path.rfind(os.sep+dir_name_old+os.sep)
        index     = max(index_new, index_old)
        dir_name_current  = [dir_name, dir_name_old][index_old > index_new]

        if index >= 0:
            index += len(dir_name_current)+2*len(os.sep)

    # if the index was -1, then we couldn't find either
    # the current directory name or the explicitly given
    # previous directory name in the path, so we cannot
    # anchor to the path name.
    if index < 0:
        warnings.warn(
            "\nCannot anchor path found outside of directory: "+"\n"+
            "   "+file_path
        )
        return file_path
    
    # rerurn the new path anchored relative to the current path
    file_path = dir_parent+os.sep+dir_name+os.sep+file_path[index:]
    return file_path
    
def find_io_paths(keys):
    whitelist = ["InputFile","InputPath","OutputFile","OutputPath","FileName"]
    blacklist = ["DataContainerArrayProxy"]
    return array_filter(keys, whitelist, blacklist)
    
def find_execute_paths(keys):
    whitelist = ["Arguments"]
    return array_filter(keys, whitelist)
    
def array_filter(keys, whitelist = None, blacklist = None):
    if not whitelist is None:
        keys = list(filter(lambda x: list(filter(lambda y: y in x, whitelist)), keys))
    if not blacklist is None:
        keys = list(filter(lambda x: list(filter(lambda y: y not in x, blacklist)), keys))
    return keys

def replace_json_paths(json_path,dir_current,dir_name_old=None):
    
    # read dream3d file
    with open(json_path, 'r') as json_file:
        data = json.load(json_file)

    # modify pathlike strings
    for filter_ in data:
        
        # modify strings behind a known file i/o key
        io_paths = find_io_paths(data[filter_].keys())
        if io_paths:
            for io_path in io_paths:
                data[filter_][io_path] = replace_path_parent(data[filter_][io_path],dir_current,dir_name_old)

        # modify strings behind a known executable key
        execute_paths = find_execute_paths(data[filter_].keys())
        if execute_paths:
            for execute_path in execute_paths:
                # not everything in an execute block will be a path, so we have to guess at what may be
                # ex:
                #    "/home/user/dream3d/dream3d.exe -p /home/user/do_stuff.json"
                #    ["/home/user/dream3d/dream3d.exe", "-p", "/home/user/do_stuff.json""]
                #    [True, False, True]
                components = shlex.split(data[filter_][execute_path],posix=False)
                components = [replace_path_parent(component,dir_current,dir_name_old) if is_path_like(component) else component for component in components]
                data[filter_][execute_path] = " ".join(components)

    # write dream3d file
    with open(json_path, 'w') as json_file:
        json.dump(data, json_file, indent=4, separators=(",", ": "), sort_keys=True)

if __name__ == '__main__':

    # find the current path
    dir_current, _ = os.path.split(os.path.realpath(__file__))

    # find all the json (dream3d) files in the current path
    json_paths = list(filter(lambda f: f.endswith(".json"), os.listdir(dir_current)))

    # anchor all path strings in the current dream3d file
    # relative to the current path name (or explicitly given dir_name_old)
    for json_path in json_paths:
        replace_json_paths(json_path,dir_current,dir_name_old)
