"""Time-driven kinematic constraints (motions).

In MSC Adams terminology a *motion* is a constraint that prescribes a
time-history on a relative kinematic coordinate between two markers or
bodies.  Motions remove one DOF per equation (like joints) but depend
explicitly on time through a driving :class:`~pmd.core.forces.Function`.

Two concrete motions are provided:

* :class:`RotMotion` – prescribed relative rotation (angular driver).
* :class:`TranMotion` – prescribed relative translation distance (linear driver).

Author: Giacomo Cangi
"""

import numpy as np

from .model import Ground
from .mechanics import rotate_90, functEval
from .constraints import Joint


# ── Abstract base ─────────────────────────────────────────────────

class Motion(Joint):
    """Abstract base class for time-driven kinematic constraints.

    Inherits from :class:`~pmd.core.constraints.Joint` so that the solver
    treats motions and geometric joints uniformly when assembling the
    constraint Jacobian and right-hand side vectors.
    """
    pass


# ── Rotational Motion ─────────────────────────────────────────────

class RotMotion(Motion):
    """Prescribed relative rotation driven by an analytical function.

    Imposes the relative orientation angle between two bodies as a
    function of time, removing one DOF.  Equivalent to Adams
    ``MOTION/ROTATIONAL``.

    Parameters
    ----------
    iBody : Body, optional
        First body.  If omitted, ``iMarker.body`` is used.
    jBody : Body or Ground, optional
        Second body.  If omitted, ``jMarker.body`` is used.  Use
        ``Ground`` to drive a body relative to the inertial frame.
    iFunct : Function
        Driving function supplying the angle (rad) and its first two
        time derivatives at each instant.
    iMarker : Marker, optional
        Attachment marker on body *i*.
    jMarker : Marker, optional
        Attachment marker on body *j*.
    name : str, optional
        Human-readable identifier.

    Examples
    --------
    Constant angular velocity of 10 deg/s on *crank* relative to Ground::

        from pmd.core import RotMotion, Function

        fn = Function(type='a', coeff=[0.0, np.deg2rad(10), 0.0])
        drv = RotMotion(iBody=crank, jBody=Ground, iFunct=fn)
    """

    def __init__(self, *, iBody=None, jBody=None, iFunct=None,
                 iMarker=None, jMarker=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.iFunct = iFunct

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 1

    def compute_phi(self, model):
        fun, _, _ = functEval(self.iFunct, model.t)
        Bi, Bj = self.iBody, self.jBody
        if Bi is Ground:
            return -Bj.orientation - fun
        elif Bj is Ground:
            return Bi.orientation - fun
        return Bi.orientation - Bj.orientation - fun

    def compute_jacobian_i(self, model):
        return np.array([[0, 0, 1]])

    def compute_jacobian_j(self, model):
        return np.array([[0, 0, -1]])

    def compute_rhs_velocity(self, model):
        _, fun_d, _ = functEval(self.iFunct, model.t)
        return fun_d

    def compute_rhs_acceleration(self, model):
        _, _, fun_dd = functEval(self.iFunct, model.t)
        return fun_dd


# ── Translational Motion ──────────────────────────────────────────

class TranMotion(Motion):
    """Prescribed relative translation driven by an analytical function.

    Imposes the distance between two markers as a function of time,
    removing one DOF.  Equivalent to Adams ``MOTION/TRANSLATIONAL``.

    Parameters
    ----------
    iMarker : Marker
        Marker on body *i*.
    jMarker : Marker
        Marker on body *j*.
    iFunct : Function
        Driving function supplying the inter-marker distance L(t) (m)
        and its first two time derivatives.
    iBody : Body, optional
        Explicit body *i* (used when ``iMarker`` is not supplied).
    jBody : Body or Ground, optional
        Explicit body *j* (used when ``jMarker`` is not supplied).
    name : str, optional
        Human-readable identifier.

    Examples
    --------
    Prescribe a sinusoidal stroke of ±50 mm at 0.5 Hz::

        from pmd.core import TranMotion, Function
        import casadi as ca

        fn = Function(expr=lambda t: 0.05 * ca.sin(2 * ca.pi * 0.5 * t))
        act = TranMotion(iMarker=mk_barrel, jMarker=mk_rod, iFunct=fn)
    """

    def __init__(self, *, iMarker=None, jMarker=None, iFunct=None,
                 iBody=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.iFunct = iFunct

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 1

    def compute_phi(self, model):
        d = self.iMarker._rP - self.jMarker._rP
        fun, _, _ = functEval(self.iFunct, model.t)
        return (d.T @ d - fun ** 2) / 2

    def compute_jacobian_i(self, model):
        iPt, jPt = self.iMarker, self.jMarker
        d = iPt._rP - jPt._rP
        return np.block([d.T, (d.T @ iPt._sPr).reshape(1, 1)])

    def compute_jacobian_j(self, model):
        iPt, jPt = self.iMarker, self.jMarker
        d = iPt._rP - jPt._rP
        return np.block([-d.T, -(d.T @ jPt._sPr).reshape(1, 1)])

    def compute_rhs_velocity(self, model):
        fun, fun_d, _ = functEval(self.iFunct, model.t)
        return fun * fun_d

    def compute_rhs_acceleration(self, model):
        iPt, jPt = self.iMarker, self.jMarker
        Bi, Bj = self.iBody, self.jBody
        d = iPt._rP - jPt._rP
        dd = iPt._drP - jPt._drP
        fun, fun_d, fun_dd = functEval(self.iFunct, model.t)
        f = fun * fun_dd + fun_d ** 2
        if Bi is Ground:
            f = f + float(d.T @ rotate_90(jPt._dsP)) * Bj.angular_velocity
        elif Bj is Ground:
            f = (f - float(d.T @ rotate_90(iPt._dsP)) * Bi.angular_velocity
                 - float(dd.T @ dd))
        else:
            f = (f + float(d.T @ rotate_90(jPt._dsP)) * Bj.angular_velocity
                 - float(d.T @ rotate_90(iPt._dsP)) * Bi.angular_velocity
                 - float(dd.T @ dd))
        return f
