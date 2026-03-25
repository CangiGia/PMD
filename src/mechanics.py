"""
Mechanics Functions for Planar Multi-Body Dynamics

This module provides kinematic and force functions for planar multi-body 
dynamic systems. It includes rotation matrices, velocity/acceleration 
computations, spring-damper-actuator force models, torsional elements,
and friction models.

Author: Giacomo Cangi
"""

import numpy as np
from numpy.typing import *


def s_rot(vect: NDArray) -> NDArray:
    """
    Computes a 90-degree counterclockwise rotation of a 2D vector,
    preserving the input shape.

    Given a 2D vector `vect` represented as a NumPy array, returns a new
    vector that is rotated 90 degrees counterclockwise.

    Args:
        vect (ndarray): A NumPy array with two elements representing the vector to be
            rotated. Can be shape (2,), (2,1), or (1,2).

    Returns:
        ndarray: A NumPy array containing the rotated vector, same shape as input.
    """
    result = np.array([-vect[1], vect[0]])
    if vect.ndim == 2:
        return result.reshape(vect.shape)
    return result


def A_matrix(phi: float) -> NDArray:
    """
    Computes the 2D rotation matrix for a given angle.

    The rotation matrix is used to rotate a vector in 2D space 
    by an angle `phi` (in radians).

    Args:
        phi (float): Angle of rotation in radians.

    Returns:
        NDArray: A 2x2 rotation matrix.
    """
    cp = np.cos(phi)
    sp = np.sin(phi)
    A = np.array([[cp, -sp], [sp, cp]])
    return A


def unit_vector(u: NDArray) -> NDArray:
    """
    Create a unit vector along the provided one.

    Args:
        u (NDArray): Vector along which to create the unit vector.

    Returns:
        NDArray: Unit vector required.
    """
    u_hat = u / np.linalg.norm(u)
    return u_hat


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

    Args:
        r (NDArray): Position of the body's reference point (e.g., origin of the
            local reference frame) in global coordinates.
        s_P_local (NDArray): Position of point `P` relative to the body's reference
            point, expressed in the local reference frame.
        A (NDArray): Rotation matrix that transforms vectors from the local reference
            frame to the global reference frame.
        s_P (NDArray, optional): Position of point `P` already expressed in global
            coordinates. If provided, used directly without transformation.

    Returns:
        NDArray: Global coordinates of point `P`.
    """
    if s_P is None:
        s_P = A @ s_P_local
        r_P = r + s_P
    else:
        r_P = r + s_P
    return r_P


def my_r_Point_d(r_d: NDArray, s_P_local: NDArray, A: NDArray, phi_d: float, s_P: NDArray = None) -> NDArray:
    """
    Calculate the velocity of a point belonging to a rigid body using 
    kinematic relations.

    Given a rigid body in motion, the velocity of a point `P` on the body
    can be expressed as the sum of the velocity of a reference point 
    (which could be the origin of the body's local reference frame or 
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

    Args:
        r_d (NDArray): Velocity vector of the reference point (e.g., center of mass
            or origin of the local reference frame), in global coordinates.
        s_P_local (NDArray): Position of point `P` expressed in the local reference
            frame.
        A (NDArray): Rotation matrix that transforms vectors from the local reference
            frame to the global reference frame.
        phi_d (float): Angular velocity of the body around the axis of rotation (rad/s).
        s_P (NDArray, optional): Position of point `P` in global coordinates. If
            provided, used directly; otherwise computed from `A` and `s_P_local`.

    Returns:
        NDArray: Velocity vector of point `P`, expressed in global coordinates.
    """
    if s_P is None:
        s_P = A @ s_P_local
        r_P_d = r_d + (s_rot(s_P) * phi_d)
    else:
        r_P_d = r_d + s_rot(s_P) * phi_d
    return r_P_d


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

    Args:
        r_dd (NDArray): Acceleration vector of the reference point (e.g., center of
            mass or origin of the local reference frame), in global coordinates.
        s_P_local (NDArray): Position of point `P` relative to the reference point,
            expressed in the local reference frame.
        A (NDArray): Rotation matrix that transforms vectors from the local reference
            frame to the global reference frame.
        phi_d (float): Angular velocity of the body (rad/s).
        phi_dd (float): Angular acceleration of the body (rad/s²).
        s_P (NDArray, optional): Position of point `P` in global coordinates. If not
            provided, computed from `A` and `s_P_local`.

    Returns:
        NDArray: Acceleration vector of point `P`, expressed in global coordinates.
    """
    if s_P is None:
        s_P = A @ s_P_local
        r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    else:
        r_P_dd = r_dd + (s_rot(s_P) * phi_dd) - (s_P * (phi_d ** 2))
    return r_P_dd


def functData(Ci, Functs):
    """
    Compute and store function coefficients for analytical constraint 
    functions of type 'a', 'b', or 'c'.

    This function populates the `coeff` and `ncoeff` attributes of 
    the Function object at index `Ci` in the `Functs` list.

    Args:
        Ci (int): Index of the function in the Functs list.
        Functs (list): List of Function objects.
    """
    funct = Functs[Ci]
    funct_type = funct.type

    if funct_type == 'a':
        funct.ncoeff = 4
        funct.coeff[3] = 2 * funct.coeff[2]

    elif funct_type == 'b':
        funct.ncoeff = 9
        xe = funct.t_end - funct.t_start
        fe = funct.f_end - funct.f_start

        C = np.array([
            [xe**3, xe**4, xe**5],
            [3 * xe**2, 4 * xe**3, 5 * xe**4],
            [6 * xe, 12 * xe**2, 20 * xe**3]
        ])

        sol = np.linalg.solve(C, np.array([fe, 0, 0]))

        funct.coeff[0:3] = sol
        funct.coeff[3] = 3 * sol[0]
        funct.coeff[4] = 4 * sol[1]
        funct.coeff[5] = 5 * sol[2]
        funct.coeff[6] = 6 * sol[0]
        funct.coeff[7] = 12 * sol[1]
        funct.coeff[8] = 20 * sol[2]

    elif funct_type == 'c':
        funct.ncoeff = 9
        xe = funct.t_end - funct.t_start
        fpe = funct.dfdt_end

        C = np.array([
            [4 * xe**3, 5 * xe**4, 6 * xe**5],
            [12 * xe**2, 20 * xe**3, 30 * xe**4],
            [24 * xe, 60 * xe**2, 120 * xe**3]
        ])

        sol = np.linalg.solve(C, np.array([fpe, 0, 0]))

        funct.coeff[0:3] = sol
        funct.coeff[3] = 4 * sol[0]
        funct.coeff[4] = 5 * sol[1]
        funct.coeff[5] = 6 * sol[2]
        funct.coeff[6] = 12 * sol[0]
        funct.coeff[7] = 20 * sol[1]
        funct.coeff[8] = 30 * sol[2]

    else:
        raise ValueError(f"Unknown function type '{funct_type}'. Valid types: 'a', 'b', 'c'.")


def functEval(funct, t):
    """
    Evaluate a Function object at time ``t``.

    Returns ``(f, f_d, f_dd)`` — function value, first and second time
    derivatives — suitable for use in ``rel-rot`` / ``rel-tran`` joint
    constraint equations.

    Args:
        funct (Function): A Function object previously processed by ``functData``.
        t (float): Current simulation time.

    Returns:
        tuple[float, float, float]: A tuple ``(f, f_d, f_dd)`` where ``f`` is the
            function value at time ``t``, ``f_d`` is the first time derivative,
            and ``f_dd`` is the second time derivative.
    """
    ftype = funct.type
    c = funct.coeff

    if ftype == 'a':
        # Polynomial: f(t) = c[0] + c[1]*t + c[2]*t^2
        # After functData: c[3] = 2*c[2] is stored.
        f   = float(c[0] + c[1] * t + c[2] * t ** 2)
        f_d = float(c[1] + c[3] * t)   # c[3] = 2*c[2]
        f_dd = float(c[3])              # = 2*c[2] (constant)

    elif ftype in ('b', 'c'):
        t0 = float(funct.t_start)
        te = float(funct.t_end)
        fs = float(funct.f_start)

        if t < t0:
            f, f_d, f_dd = fs, 0.0, 0.0

        elif t >= te:
            tau = te - t0
            if ftype == 'b':
                # After t_end: motion has reached f_end, holds constant.
                f   = float(fs + c[0]*tau**3 + c[1]*tau**4 + c[2]*tau**5)
                f_d = 0.0
                f_dd = 0.0
            else:  # type 'c': after t_end, constant velocity = dfdt_end
                f_end_val = float(fs + c[0]*tau**4 + c[1]*tau**5 + c[2]*tau**6)
                fd_end    = float(c[3]*tau**3 + c[4]*tau**4 + c[5]*tau**5)
                f   = f_end_val + fd_end * (t - te)
                f_d = fd_end
                f_dd = 0.0

        else:
            tau = t - t0
            if ftype == 'b':
                f   = float(fs + c[0]*tau**3  + c[1]*tau**4  + c[2]*tau**5)
                f_d = float(c[3]*tau**2 + c[4]*tau**3 + c[5]*tau**4)
                f_dd = float(c[6]*tau   + c[7]*tau**2 + c[8]*tau**3)
            else:  # type 'c'
                f   = float(fs + c[0]*tau**4  + c[1]*tau**5  + c[2]*tau**6)
                f_d = float(c[3]*tau**3 + c[4]*tau**4 + c[5]*tau**5)
                f_dd = float(c[6]*tau**2 + c[7]*tau**3 + c[8]*tau**4)

    else:
        raise ValueError(f"Unknown function type '{ftype}'. Valid types: 'a', 'b', 'c'.")

    return f, f_d, f_dd


# --- Point-to-point force elements ---

def pp_s(d: NDArray, k: float, L0: float):
    """
    Calculate the force generated by a spring element in a planar 
    multibody system.

    Args:
        d (NDArray): Displacement vector between two points in the system,
            representing the relative position of the spring endpoints.
        k (float): Stiffness coefficient of the spring (N/m).
        L0 (float): Rest length of the spring when no force is applied (m).

    Returns:
        NDArray: Force vector applied by the spring element, oriented along the
            direction of displacement with magnitude proportional to deformation
            from rest length (N).
    """
    L = np.linalg.norm(d)
    u = d / L
    f = k * (L - L0)
    f_s = f * u
    return f_s


def pp_sd(d: np.ndarray, d_d: np.ndarray, k: float, L0: float, dc: float) -> np.ndarray:
    """
    Calculate the force generated by a spring-damper element in a planar 
    multibody system.

    Args:
        d (ndarray): Displacement vector between two points in the system.
        d_d (ndarray): Relative velocity vector between the two points (m/s).
        k (float): Stiffness coefficient of the spring (N/m).
        L0 (float): Rest length of the spring when no force is applied (m).
        dc (float): Damping coefficient (Ns/m).

    Returns:
        ndarray: Force vector applied by the spring-damper element (N).
    """
    L = np.linalg.norm(d)
    u = d / L
    L_d = (d.T @ d_d) / L
    f = (k * (L - L0)) + (dc * L_d)
    f_sd = f * u
    return f_sd


def pp_sda(d: NDArray, d_d: NDArray, k: float, L0: float, dc: float, fa: float) -> NDArray:
    """
    Calculate the force generated by a spring-damper-actuator element in 
    a planar multibody system.

    Args:
        d (NDArray): Displacement vector between two points in the system.
        d_d (NDArray): Relative velocity vector between the two points (m/s).
        k (float): Stiffness coefficient of the spring (N/m).
        L0 (float): Rest length of the spring when no force is applied (m).
        dc (float): Damping coefficient (Ns/m).
        fa (float): Actuator force (N).

    Returns:
        NDArray: Force vector applied by the spring-damper-actuator element (N).
    """
    L = np.linalg.norm(d)
    u = d / L
    L_d = ((d.T) @ d_d) / L
    f = (k * (L - L0)) + (dc * L_d) + fa
    f_sda = f * u
    return f_sda


# --- Torsional (rotational) force elements ---

def r_s(theta: float, k: float, theta0: float):
    """
    Calculate the torque generated by a torsional spring.

    Args:
        theta (float): Current angle of the torsional spring (rad).
        k (float): Stiffness of the torsional spring (Nm/rad).
        theta0 (float): Rest angle of the torsional spring (rad).

    Returns:
        float: Torque generated by the torsional spring (Nm).
    """
    T_s = k * (theta - theta0)
    return T_s


def r_sd(theta: float, theta_d: float, k: float, theta0: float, dc: float):
    """
    Calculate the torque generated by a torsional spring and damper.

    Args:
        theta (float): Current angle (rad).
        theta_d (float): Angular velocity (rad/s).
        k (float): Stiffness of the torsional spring (Nm/rad).
        theta0 (float): Rest angle (rad).
        dc (float): Damping coefficient (Nm·s/rad).

    Returns:
        float: Torque generated by the combined spring and damper (Nm).
    """
    T_sd = k * (theta - theta0) + dc * theta_d
    return T_sd


def r_sda(theta: float, theta_d: float, k: float, theta0: float, dc: float, Ta: float):
    """
    Calculate the torque generated by a torsional actuator system.

    Args:
        theta (float): Current angle (rad).
        theta_d (float): Angular velocity (rad/s).
        k (float): Stiffness of the torsional spring (Nm/rad).
        theta0 (float): Rest angle (rad).
        dc (float): Damping coefficient (Nm·s/rad).
        Ta (float): Constant torque provided by the actuator (Nm).

    Returns:
        float: Torque generated by the combined spring, damper, and actuator (Nm).
    """
    T_sda = k * (theta - theta0) + dc * theta_d + Ta
    return T_sda


# --- Friction models ---

def friction_A(mu_s: float, mu_d: float, v_s: float, p: float, k_t: float, v: float, fN: float) -> float:
    """
    Calculate the friction force based on the Anderson et al. model.

    This function computes the friction force using a model where the 
    viscous friction is not included. The formula takes into account the 
    transition from static to dynamic friction, with an exponential decay 
    controlled by the slip velocity `v` relative to a reference slip 
    velocity `v_s`.

    Args:
        mu_s (float): Coefficient of static friction.
        mu_d (float): Coefficient of dynamic friction.
        v_s (float): Reference slip velocity for friction transition (typical: 0.001 m/s).
        p (float): Exponent parameter controlling exponential decay rate (typical: 2).
        k_t (float): Parameter scaling the hyperbolic tangent function (typical: 10).
        v (float): Relative slip velocity between contacting surfaces.
        fN (float): Normal force perpendicular to the contact surface.

    Returns:
        float: The computed friction force.
    """
    friction_force = fN * (mu_d + (mu_s - mu_d) * np.exp(-(abs(v) / v_s) ** p)) * np.tanh(k_t * v)
    return friction_force


def friction_B(mu_s: float, mu_d: float, mu_v: float, v_t: float, fnt: float, v: float, fN: float) -> float:
    """
    Calculate the friction force based on the Brown-McPhee model.

    This function computes the friction force using the Brown-McPhee model, 
    which includes the effects of viscous friction.

    Args:
        mu_s (float): Coefficient of static friction.
        mu_d (float): Coefficient of dynamic friction.
        mu_v (float): Coefficient of viscous friction.
        v_t (float): Reference velocity for friction transition.
        fnt (float): Threshold normal force for friction regime transition.
        v (float): Relative slip velocity between contacting surfaces.
        fN (float): Normal force perpendicular to the contact surface.

    Returns:
        float: The computed friction force.
    """
    vr = v / v_t
    friction_force = (fN * (mu_d * np.tanh(4 * vr) +
                    (mu_s - mu_d) * vr / (0.25 * vr**2 + 0.75)**2) +
                    mu_v * v * np.tanh(4 * fN / fnt))
    return friction_force
