import numpy as np
from numpy.typing import NDArray

def s_rot(vect: NDArray) -> NDArray:
    """
    Computes a 90-degree counterclockwise rotation of a 2D vector. 
    Given a 2D vector `vect` represented as a NumPy array, returns a new vector 
    that is rotated 90 degrees counterclockwise.

    Parameters
    ==========
    vect (ndarray)
        A NumPy array with two elements representing the vector to be rotated.

    Returns
    =======
    s_r (ndarray)
        A NumPy array containing the rotated vector.
    """
    
    s_r = np.array([-vect[1], vect[0]])
    return s_r

# function test
a = np.array([2, 3])
a_torated = s_rot(a)
eccomi = 1