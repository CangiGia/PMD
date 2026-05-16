"""Utility helpers for planar multi-body dynamics.

Author: Giacomo Cangi
"""

import numpy as np
from numpy.typing import *


def as_column_property(name):
    """
    Creates a property that ensures the assigned value is always stored
    as a column numpy array.

    This function converts the assigned value to a column numpy array of
    shape (n, 1), regardless of whether the input is a list, a 1D numpy
    array, or a 2D numpy array.

    Parameters
    ----------
    name : str
        Name of the property to define. The value is stored in a private
        attribute with an underscore prefix.

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

    Parameters
    ----------
    *args : float, int, or list
        A variable number of scalars, or a single list whose elements
        become the column vector entries.

    Returns
    -------
    numpy.ndarray
        A column vector with shape (n, 1).

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
