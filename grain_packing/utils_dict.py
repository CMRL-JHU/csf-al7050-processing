import regex as re
import json

import os, sys
new_path = os.path.dirname(__file__)
if new_path not in sys.path:
    sys.path.append(new_path)

import utils_python

class create_dict():
    
    data = {}
    
    def __init__(self, data, delimiter='/'):

        self.pattern_path_flag = r'\[self:.*?\]'
        self.pattern_path      = r'(?<=(\[self:)).*?(?=(\]))'
        self.delimiter         = delimiter
    
        # accept only dict or path string
        if not ( isinstance(data, str) or isinstance(data, dict) ):
            raise TypeError("User input must be a file path or a dictionary")
        
        # read in path file if type is string
        if isinstance(data, str):
            with open(data,"r") as f:
                self.data = json.load(f)
                
        # store if type is dict
        else: 
            self.data = data
        
        # fill in self-references (ex: [self: /some/path/in/this/file])
        self.data = self.fill_self_references_dict(self.data)
            
    def __iter__(self):
        
        return iter(self.data)
        
    def __str__(self):
        
        return utils_python.pretty_file_dict(self.data)
        
    # return the user input dictionary.
    #     path=None: the whole dictionary
    #     path=str : get the part of the dictionary under the string path
    def get(self, path='', optional=False):
        
        try:
            return self.reference_dict_by_string(
                data      = self.data,
                path      = path     ,
                operation = 'get'
            )
        except Exception as e:
            if optional:
                return None
            else:
                raise 
        
    def set(self, path, set_value):
        
        self.data = self.reference_dict_by_string(
            data      = self.data,
            path      = path     ,
            operation = 'set'    ,
            set_value = set_value
        )
        return self.data
    
    def delete(self, path=''):
        
        self.data = self.reference_dict_by_string(
            data      = self.data,
            path      = path     ,
            operation = 'delete'
        )
        return self.data
        
    def save(self, path_output, path_data=''):
        
        data = self.reference_dict_by_string(
            data      = self.data,
            path      = path_data,
            operation = 'get'
        )
        
        with open(path_output, 'w') as f:
            json.dump(data, f, indent=4, separators=(",", ": "), sort_keys=True)
            
    # rather than using[brackets][to][reference][a][dictionary],
    # reference the dictionary using('an/arbitrarily/delimited/string/path')
    def reference_dict_by_string(self, data, path, operation='get', set_value=None):
        
        # reference dictionary root path
        if path in ['', self.delimiter]:
            segments = []

        # reference a part of the dictionary
        else:
            segments = path.split(self.delimiter)
            # remove head delimiter
            if path.startswith(self.delimiter):
                segments = segments[1:]
            # check that the path exists
            if not segments[0] in data.keys() and not operation == 'set':
                raise ValueError(f"Path not found: {path}")

        # the entirety of the dictionary is being referenced
        if len(segments) == 0:
            if operation == 'get':
                return data
            if operation == 'set':
                data = set_value
                return data

        # end of path found, return result
        if len(segments) == 1:
            if operation == 'get':
                return data[segments[0]]
            if operation == 'set':
                data[segments[0]] = set_value
                return data
            if operation == 'delete':
                del data[segments[0]]
                return data

        # continue following path
        if len(segments) > 1:
            if operation == 'get':
                return self.reference_dict_by_string(
                    data      = data[segments[0]]           ,
                    path      = self.delimiter.join(segments[1:]),
                    operation = operation
                )
            if operation == 'set':
                if not segments[0] in data:
                    data[segments[0]] = {}
                data[segments[0]] = self.reference_dict_by_string(
                    data      = data[segments[0]]           ,
                    path      = self.delimiter.join(segments[1:]),
                    operation = operation                   ,
                    set_value = set_value
                )
                return data
            if operation == 'delete':
                data[segments[0]] = self.reference_dict_by_string(
                    data      = data[segments[0]]           ,
                    path      = self.delimiter.join(segments[1:]),
                    operation = operation                   ,
                    set_value = set_value
                )
                return data
    
    def fill_self_references_str(self, text):

        reference_flags = re.findall(self.pattern_path_flag, text)

        # if no self reference found, just return
        if not len(reference_flags) > 0:
            
            return text

        # if we don't have to embed the reference value into more text,
        # then the original value type can be anything
        if len(reference_flags) == 1 and reference_flags[0] == text:

            # find the reference path within the flag
            reference_path = re.search(self.pattern_path, reference_flags[0]).group().strip()

            # error out if the path is empty
            if len(reference_path) == 0:
                raise ValueError("Self reference path cannot be empty")
            
            val = self.reference_dict_by_string(self.data, reference_path, operation='get')

            if isinstance(val, str):
                val = self.fill_self_references_str(val)

            elif isinstance(val, dict):
                val = self.fill_self_references_dict(val)

            elif hasattr(val, '__len__'):
                val = self.fill_self_references_list(val)

            return val

        # if we do have to embed the reference value into more text...
        # replace all the self referencing paths with their corresponding values
        for reference_flag in reference_flags:

            # find the reference path within the flag
            reference_path = re.search(self.pattern_path, reference_flag).group().strip()
            
            # error out if the path is empty
            if len(reference_path) == 0:
                raise ValueError("Self reference path cannot be empty")

            # resolve the path
            reference_value = str(self.reference_dict_by_string(self.data, reference_path, operation='get'))
            
            # replace the path with the path value
            text = text.replace(reference_flag, reference_value)

        # check that the self referenced path did not contain another self reference
        text = self.fill_self_references_str(text)

        return text

    def fill_self_references_list(self, data):

        for i, val in enumerate(data):

            if isinstance(val, str):
                data[i] = self.fill_self_references_str(val)

            elif isinstance(val, dict):
                data[i] = self.fill_self_references_dict(val)

            elif hasattr(val, '__len__'):
                data[i] = self.fill_self_references_list(val)

        return data


    def fill_self_references_dict(self, data):

        for key, val in data.items():
            
            if isinstance(val, str):
                data[key] = self.fill_self_references_str(val)

            elif isinstance(val, dict):
                data[key] = self.fill_self_references_dict(val)

            elif hasattr(val, '__len__'):
                data[key] = self.fill_self_references_list(val)

        return data
