"""
Utilities for Planar Multi-Body Dynamics

This module provides various utility functions to support 
the development of planar multi-body dynamic models. It includes tools 
for handling global variables, grouping class instances, and managing 
vector transformations.

Author: Giacomo Cangi
"""


import inspect
import numpy as np
from numpy.typing import *


def get_globals():
    """
    Retrieve the global variables of the calling module.

    This method uses the inspect module to access the stack and 
    fetch the global variables from the caller's context.

    Returns
    -------
    dict
        A dictionary of the caller's global variables.
    """
    for frame_info in inspect.stack():
        if frame_info.function == "<module>" and frame_info.frame.f_globals.get("__name__") == "__main__":
            return frame_info.frame.f_globals
        
    return inspect.stack()[1].frame.f_globals

def group_classes():
    """
    Looping through all the global variables this method returns a 
    dictionary, in which, all the defined instances are contained.

    Returns
    -------
    dict
        A dictionary where each key is a class name (string) and the corresponding 
        value is a list of instances of that class type.
    """
    grouped_instances = {}
    global_vars = get_globals()
    
    # loop through global vars
    for var_name, var_instance in global_vars.items():
        # check if the instance has a __class__ attribute
        if hasattr(var_instance, '__class__'):
            class_name = var_instance.__class__.__name__
            if class_name not in grouped_instances:
                grouped_instances[class_name] = []
            grouped_instances[class_name].append(var_instance)

    return grouped_instances

def as_column_property(name):
    """
    Creates a property that ensures the assigned value is always stored 
    as a column numpy array.

    This function converts the assigned value to a column numpy array of 
    shape (n, 1), regardless of whether the input is a list, a 1D numpy 
    array, or a 2D numpy array.

    Parameters
    ----------
    name (str)
        The name of the property to define. The actual value is stored in a 
        private attribute named with an underscore prefix.

    Returns
    -------
    property
        A property object that enforces column-vector format on assignment.
    """
    # Determine the private attribute name based on the input
    if name.startswith("_"):
        private_name = f"__{name.lstrip('_')}"
    else:
        private_name = f"__{name}"

    @property
    def prop(self):
        return getattr(self, private_name)

    @prop.setter
    def prop(self, value):
        # convert lists to numpy arrays
        if isinstance(value, list):
            value = np.array(value)

        # ensure the array is a column vector
        if isinstance(value, np.ndarray):
            if value.ndim == 1:
                value = value.reshape(-1, 1)
            elif value.ndim == 2 and value.shape[0] == 1:
                value = value.T
        else:
            raise ValueError("The value must be a list or a NumPy array. Received: {}".format(type(value)))

        setattr(self, private_name, value)

    return prop

def colvect(*args):
    """
    Create a column vector from the provided values.

    Args
    ----
    *args : float or int
        A variable number of arguments or a list representing the elements 
        of the column vector. If a list is provided, it will be used 
        to create the column vector.

    Returns
    -------
    numpy.ndarray
        A column vector with shape (n, 1), where n is the number of provided elements.

    Examples
    --------
    >>> colvect(1, 2, 3)
    array([[1],
           [2],
           [3]])

    >>> colvect([4, 5, 6])
    array([[4],
           [5],
           [6]])
    """
    # if a single list is passed, use it; otherwise, use *args
    if len(args) == 1 and isinstance(args[0], list):
        values = args[0]
    else:
        values = args
    return np.array(values).reshape(-1, 1)

def rowvect(*args):
    """
    Create a row vector from the provided values.

    Args
    ----
    *args : float or int
        A variable number of arguments or a list representing the elements 
        of the row vector. If a list is provided, it will be used 
        to create the row vector.

    Returns
    -------
    numpy.ndarray
        A row vector with shape (1, n), where n is the number of provided elements.

    Examples
    --------
    >>> rowvect(1, 2, 3)
    array([[1, 2, 3]])

    >>> rowvect([4, 5, 6])
    array([[4, 5, 6]])
    """
    # if a single list is passed, use it; otherwise, use *args
    if len(args) == 1 and isinstance(args[0], list):
        values = args[0]
    else:
        values = args
    return np.array(values).reshape(1, -1)

def validate_shape(vec: NDArray):
    """
    Check if the input vector is a column vector with 2 rows.
    
    Args
    ----
    vec (np.ndarray)
        The input vector to be checked.

    Raises
    ------
    ValueError: if the input is not a column vector with 2 rows.
    """
    if not isinstance(vec, np.ndarray):
        raise TypeError("Input must be a numpy ndarray.")
    
    if vec.shape != (2, 1):
        if vec.shape == (1, 2):
            raise ValueError("Input is a row vector with 2 elements, expected a column vector with shape (2, 1).")
        else:
            raise ValueError("Input must be a column vector with shape (2, 1).")