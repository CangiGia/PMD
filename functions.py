import numpy as np
from numpy.typing import NDArray

def s_rot(vect: NDArray) -> NDArray:
    """
    Computes a 90-degree counterclockwise rotation of a 2D vector. 
    
    Given a 2D vector `vect` represented as a NumPy array, returns a new vector 
    that is rotated 90 degrees counterclockwise.

    Args
    ----
    vect (ndarray)
        A NumPy array with two elements representing the vector to be rotated.

    Returns
    -------
    s_r (ndarray)
        A NumPy array containing the rotated vector.
    """
    s_r = np.array([-vect[1], vect[0]])
    return s_r

def A_matrix(phi: float) -> NDArray:
    """
    Computes the 2D rotation matrix for a given angle.
    
    The rotation matrix is used to rotate a vector in 2D space 
    by an angle `phi` (in radians).

    Args
    ----
    phi (float)
        Angle of rotation in radians.

    Returns
    -------
    A (NDArray)
        A 2x2 rotation matrix.
    """
    cp = np.cos(phi)
    sp = np.sin(phi)
    A = np.array([[cp, -sp], [sp, cp]])
    return A

def r_Point(r: NDArray, s_P: NDArray) -> NDArray:
    """
    Computes the coordinates of point P relative to a global reference frame.
    
    This function calculates the global position of a point P that belongs to a body,
    using the fundamental relation of kinematics. The global position `r_P` is obtained 
    by adding the position vector `r` (which represents the position of the reference point of the body 
    in the global frame) to the local position vector `s_P` (which represents the position 
    of point P relative to the reference point of the body in the local frame).

    Args
    ----
    r (NDArray) 
        A vector representing the position of the body's reference point 
        in the global coordinate system.
    s_P (NDArray)
        A vector representing the position of point P relative to the body's 
        reference point, defined in the local coordinate system.

    Returns
    -------
    r_P (NDArray)
        A vector representing the global coordinates of point P in the global 
        reference frame.
    """    
    r_P = r + s_P
    return r_P

def r_Point_d(r_d: NDArray, s_P: NDArray, phi_d: float) -> NDArray:
    """
    Calculate the velocity of a point belonging to a rigid body using kinematic relations.

    Given a rigid body in motion, the velocity of a point P on the body can be expressed as a function
    of the velocity of a reference point (e.g., the origin or center of mass) of the body, the angular 
    velocity of the body, and the relative position of point P with respect to the reference point.

    Args
    ----
    r_d (NDArray)
        The velocity vector of the reference point (e.g., the center of mass) of the body.
    s_P (NDArray) 
        The position vector of point P relative to the reference point.
    phi_d (float) 
        The angular velocity of the body around the axis of rotation 
        (expressed in radians per second).

    Returns
    -------
    r_P_d (NDArray)
        The velocity vector of point P.
    """
    r_P_d = r_d + s_rot(s_P) * phi_d
    return r_P_d
    
def r_Point_dd(r_dd: NDArray, s_P: NDArray, phi_d: float, phi_dd: float) -> NDArray:
    """
    Calculate the acceleration of a point belonging to a rigid body using kinematic relations.

    Given a rigid body in motion, the acceleration of a point P on the body can be expressed as a 
    function of the acceleration of a reference point (e.g., the origin or center of mass) of the 
    body, the angular acceleration of the body, and the relative position of point P with respect 
    to the reference point.

    Args
    ----
    r_dd : NDArray
        The acceleration vector of the reference point (e.g., the center of mass) of the body.
    s_P : NDArray
        The position vector of point P relative to the reference point.
    phi_d : float
        The angular velocity of the body around the axis of rotation (expressed in radians per second).
    phi_dd : float
        The angular acceleration of the body around the axis of rotation (expressed in radians per second squared).

    Returns
    -------
    r_P_dd : NDArray
        The acceleration vector of point P.
    """
    r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    return r_P_dd

# function test
# s_rot
a = np.array([2, 3])
a_torated = s_rot(a)

# A_matrix
theta = np.pi/6
s_local = np.array([[2], [3]])
A = A_matrix(theta)
s = A@s_local

# example 3.1 
r = np.array([[2.5], [1.2]])
phi = np.radians(325)
s_A_local = np.array([[2.18], [0]])
s_B_local = np.array([[-1.8], [1.3]])
A = A_matrix(phi)
# a
s_A = A@s_A_local
s_B = A@s_B_local
# b
r_A = r_Point(r, s_A)
r_B = r_Point(r, s_B)
# c
s_BA = r_B - r_A

# example 3.2
r_d = np.array([[1],[-2]])
phi_d = 1

r_dd = np.array([[-0.4],[0.4]])
phi_dd = 4.65 

r_A_d = r_Point_d(r_d, s_A, phi_d)
r_A_dd = r_Point_dd(r_dd, s_A, phi_d, phi_dd)

ecchime = 1