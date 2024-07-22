import numpy as np
import numpy
from collections import defaultdict
import sys, os

import constants as cst 


sys.tracebacklimit = 8



def folder_create(folder_name):
    """
    creates a folder

    Parameters
    ----------
    folder_name : TYPE: str
        DESCRIPTION: folder name.

    Returns
    -------
    None.

    """
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
        
        
        
class overload:
    """
    A class to allow function overloading
    
    """
    @staticmethod
    def determine_types(args, kwargs, class_name):
        """
        determines input types of a function

        """
        
        arg_types = [str(class_name.__name__) if isinstance(a, class_name) else type(a) for a in args]
        kwarg_types = sorted((k, class_name if isinstance(v, class_name) else type(v)) for k, v in kwargs.items())
        return tuple(arg_types),  tuple(kwarg_types)
                
                
    function_table = defaultdict(lambda: defaultdict(dict))
    
    
    @staticmethod
    def overload_class_module(args_types=(), kwargs_types=()):
        """
        decorator for a module

        """
        
        def wrap(func):
            class_name = func.__qualname__.split('.')[0]
            named_func = overload.function_table[class_name][func.__name__]
            named_func[args_types, kwargs_types] = func
            
            def call_function_by_signature(self, *args, **kwargs):
                sig = overload.determine_types(args, kwargs, type(self))
                if sig in named_func:
                    return named_func[sig](self, *args, **kwargs)
                else:
                    n_args = []
                    for types in named_func.keys():
                        n_args.append(len(types[0]))
                        if len(types[0]) == len(sig[0]):
                            error_msg = 'No matching function for provided types. Accepted types are:\n'
                            for types in named_func.keys():
                                error_msg += f'{types[0]}'
                            raise TypeError(error_msg)
                    n_args = np.asarray(n_args)
                    error_msg = f'function {func.__name__} accepts {n_args[0]+1}'
                    if len(n_args) > 1:
                        for n in n_args[1:]:
                            error_msg += f', {n+1}'
                    error_msg += f' arguments. {len(sig)} were provided.'
                    raise TypeError(error_msg)
                    
            return call_function_by_signature
        return wrap
    
    
    @staticmethod
    def overload_func(args_types=(), kwargs_types=()):
        """
        decorator for a function

        """
        
        def wrap(func):
            named_func = overload.function_table[func.__name__]
            named_func[args_types, kwargs_types] = func
            
            def call_function_by_signature(*args, **kwargs):
                sig = overload.determine_types(args, kwargs)
                if sig in named_func:
                    return named_func[sig](*args, **kwargs)
                else:
                    n_args = []
                    for types in named_func.keys():
                        n_args.append(len(types[0]))
                        if len(types[0]) == len(sig[0]):
                            error_msg = 'No matching function for provided types. Accepted types are:\n'
                            for types in named_func.keys():
                                error_msg += f'{types[0]}'
                            raise TypeError(error_msg)
                    n_args = set(np.asarray(n_args))
                    error_msg = f'function {func.__name__} accepts '
                    for n in n_args:
                        error_msg += f'{n} '
                    error_msg += f'arguments. {len(sig)} were provided.'
                    raise TypeError(error_msg)
                    
            return call_function_by_signature
        return wrap
        

        
class exceptions:
    """
    A class containing the essential checks
    
    ...
    
    Attributes
    ----------
    
    Methods
    -------
    type_check(Input, type_target)
        checks the type of an input
    isotope_check(isotope=str):
        checks if an isotope exists
    duplicate_check(array=list/array)
        checks for duplicates in an array
    nuclei_is_defined(nuclei_list=array, label=str)
        checks if a nuclei (label) is defined
    member_of_list_check(Input, Options=list/array)
        checks if an input is among a list of options
    """
    
    @staticmethod
    def type_check(Input, type_target):
        """
        checks the type of an object and raises a TypeError is False
    
        Parameters
        ----------
        Input : TYPE: should be the same as type_target
            DESCRIPTION: can be anything.
        type_target : TYPE: type or list of types
            DESCRIPTION: type to check against.
    
        Raises
        ------
        TypeError
            DESCRIPTION: raises an error if type(Input) doesn't match type_target.
    
        Returns
        -------
        bool
            DESCRIPTION: true if Input type matches target_type
    
        """
        
        if type_target == float:
            if not isinstance(Input, (int, float)):
                raise TypeError(f'{Input} should be of type {type_target}')
        else:
            if not isinstance(Input, type_target):
                raise TypeError(f'{Input} should be of type {type_target}')  
                
        return True
    
    
    @staticmethod
    def isotope_check(isotope):
        """
        checks if an input isotope is defined in PyRelax.
    
        Parameters
        ----------
        isotope : TYPE: str or list/array
            DESCRIPTION: isotope(s).
    
        Raises
        ------
        KeyError
            DESCRIPTION: raises an error is isotope cannot be found in the defined list.
        SyntaxError
            DESCRIPTION: raises an error is isotope is not a str or list/array of str.
    
        Returns
        -------
        None.
    
        """
        
        if isinstance(isotope, str):
            if isotope not in cst.spin_quantum_number.keys():
                raise KeyError(f'isotope {isotope} cannot be found in the list of possible isotopes.')
                
        elif isinstance(isotope, (list, numpy.ndarray)):
            for iso in isotope:
                exceptions.isotope_check(iso)
                
        else:
            raise SyntaxError('isotope must be a str or list/array of str')
                        
                           
    @staticmethod
    def duplicate_check(array):
        """
        checks for duplicate items in an array
    
        Parameters
        ----------
        array : TYPE: array or list
    
        Raises
        ------
        NameError
            DESCRIPTION: raises an error if it contains duplicate items.
    
        Returns
        -------
        None.
    
        """
        
        exceptions.type_check(array, (list, np.ndarray))
        if len(array) != len(set(array)):
            raise NameError(f'{array} contains duplicate items. \
                            This will create undesired effects.')
                            
                      
    @staticmethod
    def nuclei_is_defined(nuclei_list, label):
        """
        checks is a nuclei is defined in the spin system
    
        Parameters
        ----------
        nuclei_list : TYPE: list, array or dict
            DESCRIPTION: list of nuclei labels.
        label : TYPE: str
            DESCRIPTION: label to check against.
    
        Raises
        ------
        KeyError
            DESCRIPTION: raises an error if label in not in nuclei_list.
        TypeError
            DESCRIPTION: raises an error if nuclei_list in not a list, array or dict.
    
        Returns
        -------
        None.
    
        """
        
        exceptions.type_check(label, str)
        
        if isinstance(nuclei_list, (list, np.ndarray)):
            if label not in nuclei_list:
                raise KeyError(f'{label} is not part of the spin system')
        elif isinstance(nuclei_list, dict):
            if label not in nuclei_list.keys():
                raise KeyError(f'{label} is not part of the spin system')
        else:
            raise TypeError('A nuclei list has to be provided')
            
            
    @staticmethod
    def member_of_list_check(Input, Options):
        """
        checks if an input is found in a list/array
    
        Parameters
        ----------
        Input : TYPE: any
        Options : TYPE: list/array
    
        Raises
        ------
        NameError
            DESCRIPTION: raises an error if Input for function is not found in Options.
    
        Returns
        -------
        bool
            DESCRIPTION: true if Input is in Options
    
        """
        
        if Input not in Options:
            raise NameError(f'{Input} is not an option. \
                            Choose in {Options}.')
                            
        return True