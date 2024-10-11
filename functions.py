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

def r_Point(r: NDArray, s_P: NDArray) -> NDArray: #! replaced by "my_r_Point"
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

def my_r_Point(r: NDArray, s_P_local: NDArray, A: NDArray) -> NDArray:
    """
    Computes the coordinates of point P relative to a global reference frame.
    
    This function calculates the global position of a point P that belongs to a body,
    using the fundamental relation of kinematics. The global position `r_P` is obtained 
    by adding the position vector `r` (which represents the position of the local reference 
    frame with respect to the global reference frame) to the global position vector `s_P` 
    (which represents the position of point P relative to the local reference 
    frame, expressed in global coordinates).

    Args
    ----
    r (NDArray) 
        A vector representing the position of the body's reference point 
        in the global coordinate system.
    s_P_local (NDArray)
        A vector representing the position of point P expressed in the 
        local reference frame.
    A (NDArray)
        The rotation matrix that relates the local reference frame to 
        the global reference frame.

    Returns
    -------
    r_P (NDArray)
        A vector representing the global coordinates of point P in the global 
        reference frame.
    """ 
    s_P = A@s_P_local
    r_P = r + s_P
    return r_P

def r_Point_d(r_d: NDArray, s_P: NDArray, phi_d: float) -> NDArray: #! replaced by "my_r_Point_d"
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
    
def my_r_Point_d(r_d: NDArray, s_P_local: NDArray, A: NDArray, phi_d: float) -> NDArray:
    """
    Calculate the velocity of a point belonging to a rigid body using kinematic relations.

    Given a rigid body in motion, the velocity of a point P on the body can be expressed as a function
    of the velocity of a reference point (e.g., the origin of the local reference frame) of the body, 
    the angular velocity of the body, and the relative position of point P with respect to the reference point.

    Args
    ----
    r_d (NDArray)
        The velocity vector of the reference point (e.g., the center of mass) of the body.
    s_P_local (NDArray)
        A vector representing the position of point P expressed in the 
        local reference frame.
    A (NDArray)
        The rotation matrix that relates the local reference frame to 
        the global reference frame.
    phi_d (float) 
        The angular velocity of the body around the axis of rotation 
        (expressed in radians per second).

    Returns
    -------
    r_P_d (NDArray)
        The velocity vector of point P.
    """
    s_P = A@s_P_local
    r_P_d = r_d + (s_rot(s_P) * phi_d)
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

def pp_s(d, k, L0):
    """
    Calculate the force generated by a spring element in a planar multibody system.

    Args
    ----
    d (ndarray) 
        Displacement vector between two points in the system.
        It represents the relative position of the endpoints of the spring.
    k (float) 
        Stiffness coefficient of the spring, representing its resistance 
        to deformation (N/m).
    L0 (float) 
        Rest length of the spring when no force is applied (m).

    Returns
    -------
    f_s (ndarray) 
        Force vector applied by the spring element, oriented along the 
        direction of displacement and with a magnitude proportional 
        to the deformation from the rest length (N).
    """
    L = np.linalg.norm(d)
    u = d / L 
    f = k * (L - L0)
    f_s = f * u
    return f_s

def pp_sd(d: np.ndarray, d_d: np.ndarray, k: float, L0: float, dc: float) -> np.ndarray:
    """
    Calculate the force generated by a spring-damper element in a planar multibody system.

    Args
    ----
    d (ndarray) 
        Displacement vector between two points in the system.
        It represents the relative position of the endpoints of the element.
    d_d (ndarray) 
        Relative velocity vector between the two points, representing 
        the rate of change of the displacement (m/s).
    k (float) 
        Stiffness coefficient of the spring (N/m).
    L0 (float) 
        Rest length of the spring when no force is applied (m).
    dc (float) 
        Damping coefficient, representing the resistance to relative velocity 
        between the two points (Ns/m).

    Returns
    -------
    f_sd (ndarray) 
        Force vector applied by the spring-damper element, combining 
        both spring force and damping force (N).
    """
    L = np.linalg.norm(d)
    u = d / L 
    L_d = (d.T@d_d) / L
    f = (k * (L - L0)) + (dc * L_d)
    f_sd = f * u
    return f_sd

def pp_sda(d: NDArray, d_d: NDArray, k: float, L0: float, dc: float, fa: float) -> NDArray:
    """
    Calculate the force generated by a spring-damper-actuator element in a planar 
    multibody system.

    Args
    ----
    d (ndarray) 
        Displacement vector between two points in the system.
        It represents the relative position of the endpoints of the element.
    d_d (ndarray) 
        Relative velocity vector between the two points, representing 
        the rate of change of the displacement (m/s).
    k (float) 
        Stiffness coefficient of the spring (N/m).
    L0 (float) 
        Rest length of the spring when no force is applied (m).
    dc (float) 
        Damping coefficient, representing the resistance to relative velocity 
        between the two points (Ns/m).
    fa (float) 
        Actuator force, representing an externally applied force that modifies 
        the behavior of the spring-damper system (N).

    Returns
    -------
    f_sda (ndarray)
        Force vector applied by the spring-damper-actuator element, combining 
        spring force, damping force, and actuator force (N).
    """
    L = np.linalg.norm(d)
    u = d / L 
    L_d = ((d.T)@d_d) / L
    f = (k * (L - L0)) + (dc * L_d) + fa
    f_sda = f * u
    return f_sda