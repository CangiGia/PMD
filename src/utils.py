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
from scipy.interpolate import interp1d
from numba import njit, prange


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
        
def resample(x_source: NDArray,
                y_source: NDArray,
                x_target: NDArray,
                fill_value: bool = None, 
                kind: str = 'linear' # default interpolation method
            ) -> NDArray:
    
    """
    Resample a time series or signal matrix to a new target domain.

    This function interpolates the given `y_source` data (vector or matrix)
    originally sampled at `x_source`, to the new points in `x_target`.
    The result is a new matrix (or vector) with the same number of columns
    as `y_source` but resampled along `x_target`.

    Parameters
    ----------
    x_source : NDArray
        1D array of original time or sample points corresponding to `y_source`.
    y_source : NDArray
        2D array (shape: [n_samples, n_channels]) of data values to be downsampled.
    x_target : NDArray
        1D array of new time or sample points at which to evaluate the downsampled data.
    fill_value : bool, optional
        If True, allows extrapolation outside the bounds of `x_source`.
    kind : str, optional
        The type of interpolation to use. Default is 'linear'. Other options include 
        'nearest', 'zero', 'slinear', 'quadratic', 'cubic', etc. 
        See `scipy.interpolate.interp1d` for more options.

    Returns
    -------
    y_target : NDArray
        2D array of shape (len(x_target), n_channels), containing the downsampled values.

    Examples
    --------
    >>> import numpy as np

    >>> # original time vector and signal (2 channels)
    >>> x_source = np.linspace(0, 10, 100)
    >>> y_source = np.vstack([np.sin(x_source), np.cos(x_source)]).T  # shape (100, 2)

    >>> # new time vector with fewer points
    >>> x_target = np.linspace(0, 10, 20)

    >>> # downsample the signal with interpolation
    >>> y_target = downsampling(x_source, y_source, x_target)

    >>> y_target.shape()
    (20, 2)
    """

    PARALLEL_COL_THRESHOLD = 50  # minimum number of columns to enable parallel processing
    PARALLEL_ROW_THRESHOLD = 5000  # minimum number of rows to enable parallel processing
    
    #// ... some controls ...
    #// ... check for out-of-bounds values if extrapolation is not allowed ...
    if not fill_value:
        if x_target.min() < x_source.min() or x_target.max() > x_source.max():
            raise ValueError(
                " \t ... x_target contains values outside the range of x_source, "
                "but extrapolation is disabled! \n" \
                " Set fill_value = True to allow extrapolation ..."
            )
        
    #// ... check input ...
    if not isinstance(x_source, np.ndarray) or x_source.ndim != 1:
        raise ValueError("... `x_source` must be a 1D NumPy array ...")
    
    if not isinstance(x_target, np.ndarray) or x_target.ndim != 1:
        raise ValueError("... `x_target` must be a 1D NumPy array ...")

    if not isinstance(y_source, np.ndarray) or y_source.ndim != 2:
        raise ValueError("... `y_source` must be a 2D NumPy array (n_samples, n_channels) ...")

    if len(x_source) != y_source.shape[0]:
        raise ValueError("... The length of `x_source` must match the number of rows in `y_source` ...")

    if not np.all(np.diff(x_source) > 0):
        raise ValueError("... `x_source` must be monotonically increasing ...")

    #// ... check interpolation method ...
    valid_interpolations = {'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'}
    if kind not in valid_interpolations:
        raise ValueError(f"Invalid interpolation type `{kind}`. Choose from {valid_interpolations}.")

    #// ... pre-allocate the output array ...
    y_target = np.empty((len(x_target), y_source.shape[1]))  # more efficient than np.zeros()

    channels_number = y_source.shape[1]
    samples_number = x_source.shape[0]

    @njit(parallel=True)
    def parallel_resample(x_source, y_source, x_target):
        """
        Fast interpolation using NumPy's np.interp(), fully compatible with Numba.
        This function performs interpolation in parallel across multiple columns 
        to enance performance for large datasets.
        """
        n_channels = y_source.shape[1]
        y_target_parallel = np.empty((len(x_target), n_channels))

        for i in prange(n_channels):
            y_target_parallel[:, i] = np.interp(x_target, x_source, y_source[:, i])

        return y_target_parallel

    if channels_number >= PARALLEL_COL_THRESHOLD or samples_number >= PARALLEL_ROW_THRESHOLD:
        y_target = parallel_resample(x_source, y_source, x_target) # parallel processing
    else: # standard processing
        for i in range(y_source.shape[1]):
            if fill_value:
                f = interp1d(x_source, y_source[:, i], kind=kind, fill_value="extrapolate")
            else:
                f = interp1d(x_source, y_source[:, i], kind=kind)
            y_target[:, i] = f(x_target)

    return y_target