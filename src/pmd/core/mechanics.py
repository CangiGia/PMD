"""Mechanics utilities: rotations, analytical drivers, friction models.

Author: Giacomo Cangi
"""

import numpy as np
from numpy.typing import *


def rotate_90(vect: NDArray) -> NDArray:
    """Compute a 90-degree counterclockwise rotation of a 2D vector.

    Parameters
    ----------
    vect : NDArray
        A 2-element NumPy array. Can be shape (2,), (2, 1), or (1, 2).

    Returns
    -------
    NDArray
        The rotated vector, same shape as input.
    """
    result = np.array([-vect[1], vect[0]])
    if vect.ndim == 2:
        return result.reshape(vect.shape)
    return result


def rotation_matrix(phi: float) -> NDArray:
    """Compute the 2D rotation matrix for a given angle.

    Parameters
    ----------
    phi : float
        Angle of rotation in radians.

    Returns
    -------
    NDArray
        A 2×2 rotation matrix.
    """
    cp = np.cos(phi)
    sp = np.sin(phi)
    return np.array([[cp, -sp], [sp, cp]])



def functData(Ci, Functs):
    """
    Compute and store function coefficients for analytical constraint
    functions of type 'a', 'b', or 'c'.

    This function populates the `coeff` and `ncoeff` attributes of
    the Function object at index `Ci` in the `Functs` list.

    Parameters
    ----------
    Ci : int
        Index of the function in the Functs list.
    Functs : list of Function
        List of Function objects.
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

    Parameters
    ----------
    funct : Function
        A Function object previously processed by ``functData``.
    t : float
        Current simulation time.

    Returns
    -------
    tuple of float
        ``(f, f_d, f_dd)``: function value, first derivative, and second
        derivative at time ``t``.
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


# --- Friction models ---

def friction_A(mu_s: float, mu_d: float, v_s: float, p: float, k_t: float, v: float, fN: float) -> float:
    """
    Calculate the friction force based on the Anderson et al. model.

    This function computes the friction force using a model where the
    viscous friction is not included. The formula takes into account the
    transition from static to dynamic friction, with an exponential decay
    controlled by the slip velocity `v` relative to a reference slip
    velocity `v_s`.

    Parameters
    ----------
    mu_s : float
        Coefficient of static friction.
    mu_d : float
        Coefficient of dynamic friction.
    v_s : float
        Reference slip velocity for friction transition (typical: 0.001 m/s).
    p : float
        Exponent controlling exponential decay rate (typical: 2).
    k_t : float
        Parameter scaling the hyperbolic tangent function (typical: 10).
    v : float
        Relative slip velocity between contacting surfaces.
    fN : float
        Normal force perpendicular to the contact surface.

    Returns
    -------
    float
        The computed friction force.
    """
    friction_force = fN * (mu_d + (mu_s - mu_d) * np.exp(-(abs(v) / v_s) ** p)) * np.tanh(k_t * v)
    return friction_force
