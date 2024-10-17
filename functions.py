import numpy as np
from numpy.typing import *


#* my func
import numpy as np


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

def unit_vector(u: NDArray) -> NDArray: 
    """
    Create a unit vector along the provided one. 
    
    Args
    ----
    u (NDArray)
        Vector along which create the unit vector. 

    Returns
    -------
    u_hat (NDArray)
        Unit vector required. 
    """    
    u_hat = u/np.linalg.norm(u)
    return u_hat 

#* Nikravesh
def s_rot(vect: NDArray) -> NDArray:
    """
    Computes a 90-degree counterclockwise rotation of a 2D vector. 
    
    Given a 2D vector `vect` represented as a NumPy array, returns a new
    vector that is rotated 90 degrees counterclockwise.

    Args
    ----
    vect (ndarray)
        A NumPy array with two elements representing the vector to be 
        rotated.

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

def r_Point(r: NDArray, s_P: NDArray) -> NDArray:  #! replaced by "my_r_Point"
    """
    Computes the coordinates of point `P` relative to a global reference 
    frame.
    
    This function calculates the global position of a point `P` that 
    belongs to a body, using the fundamental principles of kinematics. 
    The global position `r_P` is determined by adding the position 
    vector `r` (representing the location of the local reference 
    frame relative to the global reference frame) to the global 
    position vector `s_P` (representing the position of point `P` 
    relative to the body's reference point, defined in the local 
    coordinate system).

    Args
    ----
    r : NDArray
        A vector representing the position of the body's reference point 
        in the global coordinate system.
    s_P : NDArray
        A vector representing the position of point `P` relative to the 
        body's reference point, defined in the local coordinate system.

    Returns
    -------
    r_P : NDArray
        A vector representing the global coordinates of point `P` in the 
        global reference frame.
    """    
    r_P = r + s_P
    return r_P

def my_r_Point(r: NDArray, s_P_local: NDArray, A: NDArray, s_P: NDArray = None) -> NDArray:
    """
    Computes the coordinates of point `P` relative to a global reference 
    frame.
    
    This function calculates the global position of a point `P` that 
    belongs to a body, using the fundamental principles of kinematics. 
    If the position of point `P` relative to the local reference frame 
    is already known in global coordinates as the vector `s_P`, the 
    global position can be determined by simply adding `s_P` to the 
    position vector `r`.
    If the position of point `P` is known only relative to the local 
    reference frame (i.e., as `s_P_local`, and `s_P` is not directly known), 
    it must first be transformed into global coordinates using the rotation 
    matrix `A`, and then added to the vector `r`.
    
    Args
    ----
    r (NDArray) 
        A vector representing the position of the body's reference point 
        (e.g., the origin of the local reference frame) in the global 
        coordinate system. This defines where the local reference frame is 
        situated in global coordinates.
    s_P_local (NDArray)
        A vector representing the position of point `P` relative to the 
        body's reference point, expressed in the local reference frame.
        This defines the position of `P` within the local frame before 
        any transformation.
    A (NDArray)
        The rotation matrix that transforms vectors from the local reference 
        frame to the global reference frame. It is used to express `s_P_local` 
        in global coordinates.
    s_P (NDArray, optional)
        A vector representing the position of point `P`, already expressed 
        in the global reference frame. If this is provided, it is used directly 
        to compute the global position of `P` without the need for transformation.
        
    Returns
    -------
    r_P (NDArray)
        A vector representing the global coordinates of point `P` with 
        respect to the global reference frame.
    """ 

    if s_P is None: 
        s_P = A@s_P_local  # transform the local position s_P_local into global coordinates
        r_P = r + s_P
    else: 
        r_P = r + s_P
    return r_P

def r_Point_d(r_d: NDArray, s_P: NDArray, phi_d: float) -> NDArray:  #! replaced by "my_r_Point_d"
    """
    Calculate the velocity of a point belonging to a rigid body using 
    kinematic relations.

    Given a rigid body in motion, the velocity of a point `P` on the body
    can be expressed as the sum of the velocity of a reference point 
    (which could be the origin of the body’s local reference frame or 
    its center of mass) and the velocity due to the rotational motion 
    of the body around this reference point. The latter term is 
    calculated as the cross product of the body's angular velocity vector 
    and the relative position vector of point `P` with respect to the 
    reference point.

    This relationship is given by the equation:

        v_P = v_O + ω × r_{P/O}

    where:
    - `v_P`: velocity of point `P` in global coordinates
    - `v_O`: velocity of the reference point `O` in global coordinates
    - `ω`: angular velocity of the body (expressed in radians per second)
    - `r_{P/O}`: position vector of `P` relative to `O`, expressed in 
      the local reference frame and transformed to global coordinates.

    Args
    ----
    r_d (NDArray)
        The velocity vector of the reference point of the body, expressed 
        in global coordinates. This could be the velocity of the origin 
        of the body’s local reference frame or its center of mass.
    s_P (NDArray) 
        The position vector of point `P`, expressed in the global reference 
        frame. This represents the coordinates of point `P` in global 
        terms.
    phi_d (float) 
        The angular velocity of the body around the axis of rotation, 
        expressed in radians per second.

    Returns
    -------
    r_P_d (NDArray)
        The velocity vector of point `P`, expressed in global coordinates.
    """
    r_P_d = r_d + s_rot(s_P) * phi_d
    return r_P_d
    
def my_r_Point_d(r_d: NDArray, s_P_local: NDArray, A: NDArray, phi_d: float, s_P: NDArray = None) -> NDArray:
    """
    Calculate the velocity of a point belonging to a rigid body using 
    kinematic relations.

    Given a rigid body in motion, the velocity of a point `P` on the body
    can be expressed as the sum of the velocity of a reference point 
    (which could be the origin of the body’s local reference frame or 
    its center of mass) and the velocity due to the rotational motion 
    of the body around this reference point. The latter term is 
    calculated as the cross product of the body's angular velocity vector 
    and the relative position vector of point `P` with respect to the 
    reference point.

    This relationship is given by the equation:

        v_P = v_O + ω × r_{P/O}

    where:
    - `v_P`: velocity of point `P` in global coordinates
    - `v_O`: velocity of the reference point `O` in global coordinates
    - `ω`: angular velocity of the body (expressed in radians per second)
    - `r_{P/O}`: position vector of `P` relative to `O` in global coordinates.

    Args
    ----
    r_d (NDArray)
        The velocity vector of the reference point (e.g., the center of 
        mass or the origin of the local reference frame) of the body, 
        expressed in global coordinates.
    s_P_local (NDArray)
        A vector representing the position of point `P` expressed in the 
        local reference frame. This defines the position of `P` within 
        the local frame.
    A (NDArray)
        The rotation matrix that transforms vectors from the local reference 
        frame to the global reference frame. It is used to express `s_P_local` 
        in global coordinates.
    phi_d (float) 
        The angular velocity of the body around the axis of rotation, 
        expressed in radians per second.
    s_P (NDArray, optional)
        A vector representing the position of point `P` in global coordinates. 
        If this is provided, the function uses this value directly to 
        compute the velocity of `P`. If `None`, `s_P` will be calculated 
        using the rotation matrix `A` and the local position `s_P_local`.

    Returns
    -------
    r_P_d (NDArray)
        The velocity vector of point `P`, expressed in global coordinates.
    """
    if s_P is None: 
        s_P = A @ s_P_local
        r_P_d = r_d + (s_rot(s_P) * phi_d)
    else: 
        r_P_d = r_d + s_rot(s_P) * phi_d
    return r_P_d

def r_Point_dd(r_dd: NDArray, s_P_local: NDArray, A: NDArray, phi_d: float, phi_dd: float) -> NDArray: #! replaced by "my_r_Point_dd"
    """
    Calculate the acceleration of a point belonging to a rigid body using 
    kinematic relations.

    Given a rigid body in motion, the acceleration of a point `P` on the 
    body can be expressed as a function of the acceleration of a reference 
    point (e.g., the origin or center of mass) of the body, the angular 
    acceleration of the body, and the relative position of point `P` with 
    respect to the reference point.

    The acceleration is computed using the following equation:

        a_P = a_O + α × r_{P/O} + ω × (ω × r_{P/O})

    where:
    - `a_P`: acceleration of point `P` in global coordinates
    - `a_O`: acceleration of the reference point `O` in global coordinates
    - `α`: angular acceleration of the body
    - `r_{P/O}`: position vector of `P` relative to `O`, initially expressed 
      in the local reference frame and transformed to global coordinates.

    Args
    ----
    r_dd : NDArray
        The acceleration vector of the reference point (e.g., the center 
        of mass) of the body, expressed in global coordinates.
    s_P_local : NDArray
        The position vector of point `P` relative to the reference point,
        expressed in the local reference frame. This defines the position
        of `P` within the local frame.
    A : NDArray
        The rotation matrix that transforms vectors from the local 
        reference frame to the global reference frame. It is used to 
        express `s_P_local` in global coordinates.
    phi_d : float
        The angular velocity of the body around the axis of rotation, 
        expressed in radians per second.
    phi_dd : float
        The angular acceleration of the body around the axis of rotation, 
        expressed in radians per second squared.

    Returns
    -------
    r_P_dd : NDArray
        The acceleration vector of point `P`, expressed in global coordinates.
    """
    s_P = A @ s_P_local
    r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    return r_P_dd

def my_r_Point_dd(r_dd: NDArray, s_P_local: NDArray, A: NDArray, phi_d: float, phi_dd: float, s_P: NDArray = None) -> NDArray: 
    """
    Calculate the acceleration of a point belonging to a rigid body using 
    kinematic relations.

    Given a rigid body in motion, the acceleration of a point `P` on the 
    body can be expressed as a function of the acceleration of a reference 
    point (e.g., the origin of the local reference frame or center of mass) 
    of the body, the angular acceleration of the body, and the relative 
    position of point `P` with respect to the reference point.

    The acceleration is computed using the following equation:

        a_P = a_O + α × r_{P/O} + ω × (ω × r_{P/O})

    where:
    - `a_P`: acceleration of point `P` in global coordinates
    - `a_O`: acceleration of the reference point `O` in global coordinates
    - `α`: angular acceleration of the body
    - `r_{P/O}`: position vector of `P` relative to `O`, expressed in 
    global coordinates.

    Args
    ----
    r_dd (NDArray) 
        The acceleration vector of the reference point (e.g., the center 
        of mass or the local reference frame's origin) of the body, 
        expressed in global coordinates.
    s_P_local (NDArray)
        The position vector of point `P` relative to the reference point,
        expressed in the local reference frame. This defines the position
        of `P` within the local frame.
    A (NDArray)
        The rotation matrix that transforms vectors from the local 
        reference frame to the global reference frame. It is used to 
        express `s_P_local` in global coordinates.
    phi_d (float)
        The angular velocity of the body around the axis of rotation, 
        expressed in radians per second.
    phi_dd (float)
        The angular acceleration of the body around the axis of rotation, 
        expressed in radians per second squared.
    s_P (NDArray, optional)
        The position vector of point `P` relative to the reference point, 
        expressed in global coordinates. If not provided, it will be computed 
        from the local reference frame.

    Returns
    -------
    r_P_dd (NDArray)
        The acceleration vector of point `P`, expressed in global 
        coordinates.
    """
    if s_P is None:
        s_P = A @ s_P_local
        r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    else: 
        r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    return r_P_dd

def pp_s(d: NDArray, k: float, L0: float):
    """
    Calculate the force generated by a spring element in a planar 
    multibody system.

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
    u = d / L # unit vector
    f = k * (L - L0)
    f_s = f * u
    return f_s

def pp_sd(d: np.ndarray, d_d: np.ndarray, k: float, L0: float, dc: float) -> np.ndarray:
    """
    Calculate the force generated by a spring-damper element in a planar 
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
        Damping coefficient, representing the resistance to relative 
        velocity between the two points (Ns/m).

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
    Calculate the force generated by a spring-damper-actuator element in 
    a planar multibody system.

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
        Damping coefficient, representing the resistance to relative 
        velocity between the two points (Ns/m).
    fa (float) 
        Actuator force, representing an externally applied force that 
        modifies the behavior of the spring-damper system (N).

    Returns
    -------
    f_sda (ndarray)
        Force vector applied by the spring-damper-actuator element, 
        combining spring force, damping force, and actuator force (N).
    """
    L = np.linalg.norm(d)
    u = d / L 
    L_d = ((d.T)@d_d) / L
    f = (k * (L - L0)) + (dc * L_d) + fa
    f_sda = f * u
    return f_sda

def r_s(theta: float, k: float, theta0: float):
    """
    Calculate the torque generated by a torsional spring.

    The torque is generated by a torsional spring with stiffness `k` and 
    a rest angle `theta0`. It is proportional to the angular displacement 
    from the rest angle.

    Args
    ----
    theta (float) 
        Current angle (radians) of the torsional spring.
    k (float) 
        Stiffness of the torsional spring (Nm/rad).
    theta0 (float) 
        Rest angle of the torsional spring (radians).

    Returns
    -------
    T_s (float)
        The torque `T_s` generated by the torsional spring (Nm).
    """
    T_s = k * (theta - theta0)
    return T_s

def r_sd(theta: float, theta_d: float, k: float, theta0: float, dc: float):
    """
    Calculate the torque generated by a torsional spring and a torsional 
    damper.

    The total torque is generated by a torsional spring with stiffness `k` 
    and rest angle `theta0`, along with a torsional damper with damping 
    constant `dc` that resists angular velocity.

    Args
    ----
    theta (float) 
        Current angle (radians) of the system.
    theta_d (float) 
        Angular velocity (radians/second).
    k (float) 
        Stiffness of the torsional spring (Nm/rad).
    theta0 (float) 
        Rest angle of the torsional spring (radians).
    dc (float) 
        Damping coefficient of the torsional damper (Nm·s/rad).

    Returns
    -------
    T_sd (float) 
        The torque `T_sd` generated by the combined spring and damper (Nm).
    """
    T_sd = k * (theta - theta0) + dc * theta_d
    return T_sd

def r_sda(theta: float, theta_d: float, k: float, theta0: float, dc: float, Ta: float):
    """
    Calculate the torque generated by a torsional actuator system.

    The total torque is generated by a system consisting of a torsional 
    spring with stiffness `k` and rest angle `theta0`, a torsional damper 
    with damping constant `dc`, and an actuator that provides a constant 
    torque `Ta`.

    Args
    ----
    theta (float) 
        Current angle (radians) of the system.
    theta_d (float) 
        Angular velocity (radians/second).
    k (float) 
        Stiffness of the torsional spring (Nm/rad).
    theta0 (float) 
        Rest angle of the torsional spring (radians).
    dc (float) 
        Damping coefficient of the torsional damper (Nm·s/rad).
    Ta (float) 
        Constant torque provided by the actuator (Nm).

    Returns
    -------
    T_sda (float) 
        The torque `T_sda` generated by the combined spring, damper, and 
        actuator (Nm).
    """
    T_sda = k * (theta - theta0) + dc * theta_d + Ta
    return T_sda

def friction_A(mu_s: float, mu_d: float, v_s: float, p: float, k_t: float, v: float, fN: float) -> float:
    """
    Calculate the friction force based on the Anderson et al. model.

    This function computes the friction force using a model where the 
    viscous friction is not included. The formula takes into account the 
    transition from static to dynamic friction, with an exponential decay 
    controlled by the slip velocity `v` relative to a reference slip 
    velocity `v_s`.

    Args
    ----
    mu_s (float)
        The coefficient of static friction. Represents the friction 
        when the body is at rest.
    mu_d (float)
        The coefficient of dynamic friction. Represents the friction 
        when the body is in motion.
    v_s (float)
        The reference slip velocity where the transition between static 
        and dynamic friction occurs.
    p (float)
        The exponent parameter that controls the rate of exponential 
        decay in the friction model.
    k_t (float)
        A parameter that scales the hyperbolic tangent function, affecting 
        how sharply the friction changes with velocity.
    v (float)
        The relative slip velocity between the contacting surfaces.
    fN (float)
        The normal force, representing the force perpendicular to the 
        contact surface.

    Returns
    -------
    friction_force (float)
        The computed friction force.
        
    Note
    ----
    Typical values for the `v_s`, `p` and `k_t` are: `0.001` (m/s), `2` and 
    `10` respectively.
    """
    friction_force = fN * (mu_d + (mu_s - mu_d) * np.exp(-(abs(v) / v_s) ** p)) * np.tanh(k_t * v)
    return friction_force

def friction_B(mu_s: float, mu_d: float, mu_v: float, v_t: float, fnt: float, v: float, fN: float) -> float:
    """
    Calculate the friction force based on the Brown-McPhee model.

    This function computes the friction force using the Brown-McPhee model, 
    which includes the effects of viscous friction. The formula accounts for 
    both the transition from static to dynamic friction as well as the 
    influence of relative slip velocity and normal force.

    Args
    ----
    mu_s (float)
        The coefficient of static friction. Represents the friction 
        when the body is at rest.
    mu_d (float)
        The coefficient of dynamic friction. Represents the friction 
        when the body is in motion.
    mu_v (float)
        The coefficient of viscous friction, representing the friction 
        that is proportional to velocity.
    v_t (float)
        A reference velocity that scales the transition between static 
        and dynamic friction.
    fnt (float)
        A threshold normal force that affects the transition between 
        different friction regimes.
    v (float)
        The relative slip velocity between the contacting surfaces.
    fN (float)
        The normal force, representing the force perpendicular to the 
        contact surface.

    Returns
    -------
    friction_force (float)
        The computed friction force.
    """
    vr = v / v_t
    friction_force = (fN * (mu_d * np.tanh(4 * vr) + 
                    (mu_s - mu_d) * vr / (0.25 * vr**2 + 0.75)**2) +
                    mu_v * v * np.tanh(4 * fN / fnt))
    return friction_force