"""Force hierarchy and analytical driver functions.

This module contains all force elements and the :class:`Function` class
used to drive motions and actuators.

Force classes
-------------
* :class:`Force` – abstract base
* :class:`Weight` – gravitational load
* :class:`PtpForce` – point-to-point spring/damper/actuator
* :class:`RotSdaForce` – rotational spring/damper/actuator
* :class:`LocalForce` – constant force in body-local frame
* :class:`GlobalForce` – constant force in global frame
* :class:`Torque` – constant torque
* :class:`UserForce` – user-defined force via Python callback
* :class:`Actuator` – composite length- or force-controlled actuator

Driver functions
----------------
* :class:`Function` – polynomial or symbolic time-dependent driver

Data types
----------
* :class:`BodyState` – named-tuple snapshot of a body's kinematic state

Author: Giacomo Cangi
"""

import logging
from abc import ABC, abstractmethod
from collections import namedtuple

import numpy as np
from numpy.typing import *

from .utils import *
from .model import Base, Ground

# ─── BodyState ───────────────────────────────────────────────────────────────

BodyState = namedtuple(
    'BodyState',
    ['position', 'velocity', 'orientation', 'angular_velocity'],
)
"""Snapshot of a body's kinematic state passed to :class:`UserForce` callbacks.

Fields
------
position : ndarray, shape (2, 1)
    Global position vector ``[x, y]``.
velocity : ndarray, shape (2, 1)
    Global velocity vector ``[vx, vy]``.
orientation : float
    Body angle (rad).
angular_velocity : float
    Body angular velocity (rad/s).
"""

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
#  FORCE HIERARCHY
# ═══════════════════════════════════════════════════════════════════

class Force(ABC, Base):
    """Abstract base class for all force elements.

    Every concrete force must implement the :meth:`apply` method.

    Attributes
    ----------
    iMarker : Marker or None
        Marker on body *i*.
    jMarker : Marker or None
        Marker on body *j*.
    iBody : Body
        Body *i*.
    jBody : Body
        Body *j*.
    name : str or None
        Optional human-readable name for identification and legends.
    """

    def __init__(self, *, iMarker=None, jMarker=None,
                 iBody=None, jBody=None, name=None):
        """Initialize a Force.

        Parameters
        ----------
        iMarker : Marker, optional
            Marker on body *i*.
        jMarker : Marker, optional
            Marker on body *j*.
        iBody : Body, optional
            Explicit body *i* (used if *iMarker* is None).
        jBody : Body, optional
            Explicit body *j* (used if *jMarker* is None).
        name : str, optional
            Human-readable name.
        """
        super().__init__()
        self.iMarker = iMarker
        self.jMarker = jMarker
        self.iBody = iMarker.body if iMarker is not None else (iBody or Ground)
        self.jBody = jMarker.body if jMarker is not None else (jBody or Ground)
        self.name = name

    @abstractmethod
    def apply(self, bodies):
        """Apply this force to the relevant bodies.

        Parameters
        ----------
        bodies : list of Body
            List of Body objects in the model.
        """

    def __repr__(self):
        label = f"'{self.name}'" if self.name else ""
        body_info = f"iBody={self.iBody!r}, jBody={self.jBody!r}"
        if label:
            return f"{self.__class__.__name__}({label}, {body_info})"
        return f"{self.__class__.__name__}({body_info})"


# ── Weight ────────────────────────────────────────────────────────

class Weight(Force):
    """Gravitational force applied to all bodies.

    Attributes
    ----------
    gravity : float
        Gravitational acceleration magnitude (m/s²).
    gravity_direction : ndarray
        Unit direction of gravity, shape (2, 1).
    """

    DEFAULT_GRAVITY = 9.81
    DEFAULT_GRAVITY_VECTOR = colvect([0, -1])

    def __init__(self, gravity=None, gravity_direction=None, name=None):
        super().__init__(name=name)
        self.gravity = gravity if gravity is not None else self.DEFAULT_GRAVITY
        self.gravity_direction = (gravity_direction if gravity_direction is not None
                                  else self.DEFAULT_GRAVITY_VECTOR.copy())
        self._initialized = False

    def initialize_weights(self, bodies):
        """Pre-compute weight vectors for all bodies.

        Called once by ``PlanarMultibodyModel._initialize()``.

        Parameters
        ----------
        bodies : list of Body
            List of Body objects.
        """
        ug = self.gravity * self.gravity_direction
        for body in bodies:
            body._weight = body.mass * ug
        self._initialized = True

    def apply(self, bodies):
        for body in bodies:
            body._force += body._weight

    def __repr__(self):
        label = f"'{self.name}', " if self.name else ""
        return f"Weight({label}gravity={self.gravity})"


# ── Point-to-Point Force ─────────────────────────────────────────

class PtpForce(Force):
    """Point-to-point spring-damper-actuator force.

    Attributes
    ----------
    k : float
        Spring stiffness (N/m).
    L0 : float
        Undeformed spring length (m).
    dc : float
        Damping coefficient (N·s/m).
    f_a : float
        Constant actuator force (N).
    """

    def __init__(self, *, iMarker=None, jMarker=None, k=0, L0=0, dc=0,
                 f_a=0, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker, name=name)
        self.k = k
        self.L0 = L0
        self.dc = dc
        self.f_a = f_a

    def apply(self, bodies):
        iPt = self.iMarker
        jPt = self.jMarker
        Bi = self.iBody
        Bj = self.jBody
        d = iPt._rP - jPt._rP
        dd = iPt._drP - jPt._drP
        L = np.linalg.norm(d)
        dL = d.T @ dd / L
        delta = L - self.L0
        u = d / L
        f = self.k * delta + self.dc * dL + self.f_a
        fi = f * u
        if Bi is not Ground:
            Bi._force -= fi
            Bi._torque -= (iPt._sPr.T @ fi).item()
        if Bj is not Ground:
            Bj._force += fi
            Bj._torque += (jPt._sPr.T @ fi).item()

    def __repr__(self):
        label = f"'{self.name}', " if self.name else ""
        return f"PtpForce({label}k={self.k}, L0={self.L0}, dc={self.dc})"


# ── Rotational Spring-Damper-Actuator ─────────────────────────────

class RotSdaForce(Force):
    """Rotational spring-damper-actuator torque.

    Attributes
    ----------
    k : float
        Torsional stiffness (N·m/rad).
    theta0 : float
        Undeformed angle (rad).
    dc : float
        Damping coefficient (N·m·s/rad).
    T_a : float
        Constant actuator torque (N·m).
    """

    def __init__(self, *, iBody=None, jBody=None, k=0, theta0=0, dc=0,
                 T_a=0, iMarker=None, jMarker=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.k = k
        self.theta0 = theta0
        self.dc = dc
        self.T_a = T_a

    def apply(self, bodies):
        Bi = self.iBody
        Bj = self.jBody
        if Bi is Ground:
            theta = -Bj.orientation
            theta_d = -Bj.angular_velocity
            T = self.k * (theta - self.theta0) + self.dc * theta_d + self.T_a
            Bj._torque += T
        elif Bj is Ground:
            theta = Bi.orientation
            theta_d = Bi.angular_velocity
            T = self.k * (theta - self.theta0) + self.dc * theta_d + self.T_a
            Bi._torque -= T
        else:
            theta = Bi.orientation - Bj.orientation
            theta_d = Bi.angular_velocity - Bj.angular_velocity
            T = self.k * (theta - self.theta0) + self.dc * theta_d + self.T_a
            Bi._torque -= T
            Bj._torque += T


# ── Local Force ───────────────────────────────────────────────────

class LocalForce(Force):
    """Constant force in the body-local reference frame.

    Attributes
    ----------
    force_local : ndarray
        Force vector in local coordinates, shape (2, 1).
    """

    def __init__(self, *, iBody=None, force_local=None, iMarker=None,
                 jMarker=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.force_local = (force_local if force_local is not None
                            else colvect([0, 0]))

    # Ensure stored as column vector
    force_local = as_column_property("force_local")

    def apply(self, bodies):
        self.iBody._force += self.iBody._rotation_matrix @ self.force_local


# ── Global Force ──────────────────────────────────────────────────

class GlobalForce(Force):
    """Constant force in the global reference frame.

    Attributes
    ----------
    force_global : ndarray
        Force vector in global coordinates, shape (2, 1).
    """

    def __init__(self, *, iBody=None, force_global=None, iMarker=None,
                 jMarker=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.force_global = (force_global if force_global is not None
                             else colvect([0, 0]))

    # Ensure stored as column vector
    force_global = as_column_property("force_global")

    def apply(self, bodies):
        self.iBody._force += self.force_global


# ── Torque ────────────────────────────────────────────────────────

class Torque(Force):
    """Constant torque in the global reference frame.

    Attributes
    ----------
    torque_value : float
        Scalar torque value (N·m).
    """

    def __init__(self, *, iBody=None, torque_value=0, iMarker=None,
                 jMarker=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.torque_value = torque_value

    def apply(self, bodies):
        self.iBody._torque += self.torque_value


# ── User Force ────────────────────────────────────────────────────

class UserForce(Force):
    """User-defined force.

    Two mutually exclusive modes:

    **Callback mode** (``callback=...``)
        A Python callable invoked at every integration step with
        ``(t, iBody_state, jBody_state)`` and returning a dict::

            {'force': [Fx, Fy], 'torque': T}

        where ``iBody_state`` and ``jBody_state`` are
        :class:`~pmd.core.forces.BodyState` named-tuples exposing
        ``.position``, ``.velocity``, ``.orientation``,
        ``.angular_velocity``.  Use this mode when the force law
        cannot be expressed symbolically (e.g. look-up tables,
        external controllers, I/O).  Each integration step crosses
        the Python/C boundary — suitable for prototyping.

    **Symbolic mode** (``expression=...``) *(planned for Fase C)*
        A CasADi-compatible callable embedded in the symbolic
        equations of motion and JIT-compiled for maximum performance.

    Parameters
    ----------
    callback : callable, optional
        Callback ``(t, iBody_state, jBody_state) -> dict``.
    iMarker : Marker, optional
    jMarker : Marker, optional
    iBody : Body, optional
    jBody : Body, optional
    name : str, optional
    """

    def __init__(self, *, callback=None, iMarker=None, jMarker=None,
                 iBody=None, jBody=None, name=None, **kwargs):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
        self.callback = callback
        # Store any extra kwargs for user convenience (e.g. k, L0, dc)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def apply(self, bodies):
        if self.callback is None or not callable(self.callback):
            return

        t = getattr(self, '_t', 0.0)

        def _make_state(body):
            if body is Ground:
                return BodyState(
                    position=np.zeros((2, 1)),
                    velocity=np.zeros((2, 1)),
                    orientation=0.0,
                    angular_velocity=0.0,
                )
            return BodyState(
                position=body.position.copy(),
                velocity=body.velocity.copy(),
                orientation=float(body.orientation),
                angular_velocity=float(body.angular_velocity),
            )

        i_state = _make_state(self.iBody)
        j_state = _make_state(self.jBody)
        result = self.callback(t, i_state, j_state)
        if result is None:
            return

        iforce  = result.get('force',   None)
        itorque = result.get('torque',  None)
        jforce  = result.get('jforce',  None)
        jtorque = result.get('jtorque', None)

        if iforce is not None and self.iBody is not Ground:
            self.iBody._force  += np.asarray(iforce, dtype=float).reshape(2, 1)
        if itorque is not None and self.iBody is not Ground:
            self.iBody._torque += float(itorque)
        if jforce is not None and self.jBody is not Ground:
            self.jBody._force  += np.asarray(jforce, dtype=float).reshape(2, 1)
        if jtorque is not None and self.jBody is not Ground:
            self.jBody._torque += float(jtorque)

    def __repr__(self):
        label = f"'{self.name}', " if self.name else ""
        return f"UserForce({label}callback={self.callback!r})"


# ═══════════════════════════════════════════════════════════════════
#  FUNCTION CLASS
# ═══════════════════════════════════════════════════════════════════

class Function(Base):
    """Analytical driver function for motions and actuators.

    **Polynomial mode** (default)
        Specify ``type`` and ``coeff``; the function evaluates a
        polynomial up to degree 8 in time, with types ``'a'``, ``'b'``,
        ``'c'`` matching the legacy PMD conventions.

    **Symbolic mode**
        Pass ``expr=callable`` where the callable receives a CasADi SX
        symbol ``t`` and returns a scalar SX expression.  The solver
        embeds the expression directly in the compiled equations of
        motion.  First and second derivatives are computed via
        ``ca.jacobian`` automatically.  *(Implemented in Fase C.)*

    Parameters
    ----------
    type : {'a', 'b', 'c'}
        Polynomial function type.
    coeff : array_like, optional
        Polynomial coefficients ``[c0, c1, c2, ...]`` so that
        ``f(t) = c0 + c1*t + c2*t² + …``.
    expr : callable, optional
        Symbolic callable ``expr(t: ca.SX) -> ca.SX``.  When provided,
        ``type`` and ``coeff`` are ignored.

    Attributes
    ----------
    type : str
    t_start, f_start, t_end, f_end, dfdt_end : float
        Parameters for type ``'b'`` and ``'c'`` functions.
    ncoeff : int
        Number of active coefficients.
    coeff : numpy.ndarray
        Padded coefficient array (length 9).
    expr : callable or None
        Symbolic expression callable (None in polynomial mode).
    """

    def __init__(self, type='a', t_start=0, f_start=0, t_end=1, f_end=1,
                 dfdt_end=1, ncoeff=4, coeff=None, expr=None):
        super().__init__()
        self.type = type
        self.t_start = t_start
        self.f_start = f_start
        self.t_end = t_end
        self.f_end = f_end
        self.dfdt_end = dfdt_end
        self.ncoeff = ncoeff
        self.expr = expr
        if coeff is not None:
            c = np.asarray(coeff, dtype=float).flatten()
            padded = np.zeros(9)
            padded[:min(len(c), 9)] = c[:min(len(c), 9)]
            self.coeff = padded
        else:
            self.coeff = np.zeros(9)

    def __repr__(self):
        if self.expr is not None:
            return f"Function(expr={self.expr!r})"
        return f"Function(type='{self.type}')"


# ═══════════════════════════════════════════════════════════════════
#  ACTUATOR
# ═══════════════════════════════════════════════════════════════════

class Actuator(Force):
    """Composite length- or force-controlled actuator.

    Wraps a driver :class:`Function` law and expands into primitive
    model elements when the model is initialized.

    **Force control** (``control='force'``)
        Applies an axial force along the line connecting *iMarker* and
        *jMarker*.  Positive ``law(t)`` is a tensile (contracting) force.
        The actuator is evaluated numerically at each integration step.

    **Length control** (``control='length'``)
        Prescribes the distance between *iMarker* and *jMarker* via a
        :class:`~pmd.core.motion.TranMotion` joint.  During model
        initialization, :meth:`expand` appends the joint to the model.

    Parameters
    ----------
    iMarker : Marker
        Attachment marker on body *i*.
    jMarker : Marker
        Attachment marker on body *j*.
    control : {'force', 'length'}
        Actuator mode.
    law : Function or callable
        Driver law.  A plain callable is auto-wrapped as
        ``Function(expr=law)``.
    name : str, optional
    """

    def __init__(self, *, iMarker, jMarker, control='force', law, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker, name=name)
        if control not in ('force', 'length'):
            raise ValueError(
                f"Actuator control must be 'force' or 'length', got '{control}'.")
        self.control = control
        if callable(law) and not isinstance(law, Function):
            law = Function(expr=law)
        self.law = law

    def expand(self, model):
        """Register the internal joint for length-control mode.

        Called by :meth:`~pmd.core.solver.PlanarMultibodyModel._initialize`
        before joint initialization.  For force-control mode this is a no-op.

        Parameters
        ----------
        model : PlanarMultibodyModel
            The parent model.
        """
        if self.control == 'length':
            from .motion import TranMotion
            jt = TranMotion(
                iMarker=self.iMarker,
                jMarker=self.jMarker,
                iFunct=self.law,
                name=self.name,
            )
            model.Joints.append(jt)

    def apply(self, bodies):
        if self.control != 'force':
            return  # length mode → handled by the TranMotion joint

        from .mechanics import functEval
        t = getattr(self, '_t', 0.0)
        f_val, _, _ = functEval(self.law, t)

        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        L = float(np.linalg.norm(d))
        if L < 1e-12:
            return
        u = d / L
        fi = f_val * u  # tensile: pulls iBody toward jBody

        if self.iBody is not Ground:
            self.iBody._force  -= fi
            self.iBody._torque -= float((iPt._sPr.T @ fi).item())
        if self.jBody is not Ground:
            self.jBody._force  += fi
            self.jBody._torque += float((jPt._sPr.T @ fi).item())

    def __repr__(self):
        label = f"'{self.name}', " if self.name else ""
        return f"Actuator({label}control='{self.control}', law={self.law!r})"
