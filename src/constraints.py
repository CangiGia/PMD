"""Planar Multi-Body Dynamics Constraints and Forces.

This module provides polymorphic Joint and Force hierarchies using ABCs,
plus the Function class for analytical constraint drivers.

Author: Giacomo Cangi
"""

import logging
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import *

from .utils import *
from .model import Base, Ground
from .mechanics import s_rot, functEval

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
#  JOINT HIERARCHY
# ═══════════════════════════════════════════════════════════════════

class Joint(ABC, Base):
    """Abstract base class for all joint constraints.

    Every concrete joint must implement the abstract methods below.

    Attributes:
        iMarker: Marker on body i (or None).
        jMarker: Marker on body j (or None).
        iBody: Body i (derived from iMarker or explicit).
        jBody: Body j (derived from jMarker or explicit).
    """

    def __init__(self, *, iMarker=None, jMarker=None,
                 iBody=None, jBody=None):
        super().__init__()
        self.iMarker = iMarker
        self.jMarker = jMarker
        self.iBody = iMarker.body if iMarker is not None else (iBody or Ground)
        self.jBody = jMarker.body if jMarker is not None else (jBody or Ground)

        # Validation: markers must be on different bodies
        if iMarker is not None and jMarker is not None:
            if iMarker.body is jMarker.body:
                raise ValueError(
                    "iMarker and jMarker must be on different bodies")

        # Solver-internal bookkeeping (set by _initialize)
        self._nbody = 2
        self._mrows = 0
        self._rows = 0
        self._rowe = 0
        self._colis = 0
        self._colie = 0
        self._coljs = 0
        self._colje = 0
        self._lagrange = np.zeros([3, 1])

    # ── Abstract interface ────────────────────────────────────────

    @abstractmethod
    def initialize(self, model):
        """Compute derived quantities after model construction.

        Called once by PlanarMultibodyModel._initialize().

        Args:
            model: The PlanarMultibodyModel instance.
        """

    @abstractmethod
    def compute_phi(self, model):
        """Evaluate constraint equations.

        Args:
            model: The PlanarMultibodyModel instance.

        Returns:
            ndarray: Constraint residual, shape (mrows, 1).
        """

    @abstractmethod
    def compute_jacobian_i(self, model):
        """Compute the Jacobian block for body i.

        Args:
            model: The PlanarMultibodyModel instance.

        Returns:
            ndarray: Jacobian block for body i, shape (mrows, 3).
        """

    @abstractmethod
    def compute_jacobian_j(self, model):
        """Compute the Jacobian block for body j.

        Args:
            model: The PlanarMultibodyModel instance.

        Returns:
            ndarray: Jacobian block for body j, shape (mrows, 3).
        """

    @abstractmethod
    def compute_rhs_velocity(self, model):
        """Compute right-hand side of velocity constraint.

        Args:
            model: The PlanarMultibodyModel instance.

        Returns:
            ndarray or None: RHS velocity, shape (mrows, 1), or None if zero.
        """

    @abstractmethod
    def compute_rhs_acceleration(self, model):
        """Compute right-hand side of acceleration constraint (gamma).

        Args:
            model: The PlanarMultibodyModel instance.

        Returns:
            ndarray: RHS acceleration, shape (mrows, 1).
        """

    def fk_step(self, parent, child):
        """Forward kinematics: place child given parent state and q0.

        Default implementation does nothing (for joints without FK role).
        Override in subclasses that participate in assembly FK.

        Args:
            parent: The parent body.
            child: The child body.
        """

    def __repr__(self):
        return (f"{self.__class__.__name__}("
                f"iBody={self.iBody!r}, jBody={self.jBody!r})")


# ── Revolute Joint ────────────────────────────────────────────────

class RevJoint(Joint):
    """Revolute (pin) joint between two bodies.

    Attributes:
        fix: If 1, the relative rotation is also constrained.
        q0: Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, fix=0, q0=0,
                 iBody=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.fix = fix
        self.q0 = q0
        self._p0 = 0

    def initialize(self, model):
        self._mrows = 2
        self._nbody = 2
        self.iBody = self.iMarker.body
        self.jBody = self.jMarker.body
        if self.fix == 1:
            self._mrows = 3
            Bi = self.iBody
            Bj = self.jBody
            if Bi is Ground:
                self._p0 = -Bj.orientation
            elif Bj is Ground:
                self._p0 = Bi.orientation
            else:
                self._p0 = Bi.orientation - Bj.orientation

    def compute_phi(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        f = iPt._rP - jPt._rP
        if self.fix == 1:
            if self.iBody is Ground:
                f = np.append(f, (-self.jBody.orientation - self._p0))
            elif self.jBody is Ground:
                f = np.append(f, (self.iBody.orientation - self._p0))
            else:
                f = np.append(f, (self.iBody.orientation - self.jBody.orientation - self._p0))
        return f

    def compute_jacobian_i(self, model):
        iPt = self.iMarker
        Di = np.block([[np.eye(2), iPt._sPr.reshape(2, 1)]])
        if self.fix == 1:
            Di = np.vstack([Di, [0, 0, 1]])
        return Di

    def compute_jacobian_j(self, model):
        jPt = self.jMarker
        Dj = np.block([[-np.eye(2), -jPt._sPr.reshape(2, 1)]])
        if self.fix == 1:
            Dj = np.vstack([Dj, [0, 0, -1]])
        return Dj

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        Bi = self.iBody
        Bj = self.jBody
        if Bi is Ground:
            f = s_rot(jPt._dsP) * Bj.angular_velocity
        elif Bj is Ground:
            f = -s_rot(iPt._dsP) * Bi.angular_velocity
        else:
            f = (-s_rot(iPt._dsP) * Bi.angular_velocity
                 + s_rot(jPt._dsP) * Bj.angular_velocity)
        if self.fix == 1:
            f = np.vstack([f, [0]])
        return f

    def fk_step(self, parent, child):
        from .mechanics import A_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = A_matrix(parent.orientation)
        parent_sP = parent._rotation_matrix @ parent_marker.local_position.reshape(2, 1)
        r_joint = parent.position + parent_sP
        child.orientation = parent.orientation + q0
        child._rotation_matrix = A_matrix(child.orientation)
        child_sP = child._rotation_matrix @ child_marker.local_position.reshape(2, 1)
        child.position = r_joint - child_sP


# ── Translational Joint ──────────────────────────────────────────

class TranJoint(Joint):
    """Translational (prismatic) joint between two bodies.

    Attributes:
        fix: If 1, the relative translation is also constrained.
        q0: Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, fix=0, q0=0,
                 iBody=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.fix = fix
        self.q0 = q0
        self._p0 = 0

        # Validation: tran requires orientation on both markers
        if iMarker is not None and jMarker is not None:
            if not iMarker.has_orientation or not jMarker.has_orientation:
                raise ValueError(
                    "Translational joint requires both iMarker and jMarker "
                    "to have orientation (theta != None)")

    def initialize(self, model):
        self._mrows = 2
        self._nbody = 2
        self.iBody = self.iMarker.body
        self.jBody = self.jMarker.body
        if self.fix == 1:
            self._mrows = 3
            Bi = self.iBody
            Bj = self.jBody
            iPt = self.iMarker
            jPt = self.jMarker
            if Bi is Ground:
                self._p0 = np.linalg.norm(
                    iPt._rP - Bj.position - Bj._rotation_matrix @
                    jPt.local_position.reshape(2, 1))
            elif Bj is Ground:
                self._p0 = np.linalg.norm(
                    Bi.position + Bi._rotation_matrix @
                    iPt.local_position.reshape(2, 1) - jPt._rP)
            else:
                self._p0 = np.linalg.norm(
                    Bi.position + Bi._rotation_matrix @
                    iPt.local_position.reshape(2, 1) -
                    Bj.position - Bj._rotation_matrix @
                    jPt.local_position.reshape(2, 1))

    def compute_phi(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        ujr = self.jMarker._ur
        ui = self.iMarker._u
        d = iPt._rP - jPt._rP
        f = np.array([ujr.T @ d, ujr.T @ ui]).reshape(2, 1)
        if self.fix == 1:
            f = np.append(f, (ui.T @ d - self._p0) / 2).reshape(3, 1)
        return f

    def compute_jacobian_i(self, model):
        iPt = self.iMarker
        uj = self.jMarker._u
        ujr = self.jMarker._ur
        Di = np.block([
            [ujr.T, (uj.T @ iPt._sP).reshape(1, 1)],
            [np.array([0, 0, 1])]
        ])
        if self.fix == 1:
            Di = np.vstack([
                Di,
                [uj.T, (uj.T @ iPt._sPr).reshape(1)]
            ])
        return Di

    def compute_jacobian_j(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        uj = self.jMarker._u
        ujr = self.jMarker._ur
        d = iPt._rP - jPt._rP
        Dj = np.block([
            [-ujr.T, -(uj.T @ (jPt._sP + d)).reshape(1, 1)],
            [np.array([0, 0, -1])]
        ])
        if self.fix == 1:
            Dj = np.vstack([
                Dj,
                [-uj.T, -(uj.T @ jPt._sPr).reshape(1)]
            ])
        return Dj

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        Bi = self.iBody
        Bj = self.jBody
        iPt = self.iMarker
        jPt = self.jMarker
        ujd = self.jMarker._du
        ujdr = s_rot(ujd)

        if Bi is Ground:
            f2 = 0.0
        elif Bj is Ground:
            f2 = 0.0
        else:
            diffr = Bi.position - Bj.position
            dp_product = (ujd.T @ diffr).item() * Bi.angular_velocity
            diffdr = Bi.velocity - Bj.velocity
            f2 = dp_product - 2.0 * (ujdr.T @ diffdr).item()

        f = np.array([[f2], [0.0]])

        if self.fix == 1:
            d = iPt._rP - jPt._rP
            dd = iPt._drP - jPt._drP
            L = self._p0
            u = d / L
            du = dd / L
            f3 = -(du.T @ dd).item()

            if Bi is Ground:
                f3 += (u.T @ (s_rot(jPt._dsP) * Bj.angular_velocity)).item()
            elif Bj is Ground:
                f3 -= (u.T @ (s_rot(iPt._dsP) * Bi.angular_velocity)).item()
            else:
                term1 = iPt._dsP * Bi.angular_velocity
                term2 = jPt._dsP * Bj.angular_velocity
                f3 -= (u.T @ s_rot(term1 - term2)).item()

            f = np.vstack([f, [[f3]]])

        return f

    def fk_step(self, parent, child):
        from .mechanics import A_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = A_matrix(parent.orientation)
        child.orientation = parent.orientation
        child._rotation_matrix = A_matrix(child.orientation)
        if parent_marker.has_orientation:
            slide_dir = parent._rotation_matrix @ parent_marker._ulocal
        else:
            slide_dir = child._rotation_matrix @ child_marker._ulocal
        parent_sP = parent._rotation_matrix @ parent_marker.local_position.reshape(2, 1)
        r_joint = parent.position + parent_sP
        child_sP = child._rotation_matrix @ child_marker.local_position.reshape(2, 1)
        child.position = r_joint - child_sP + q0 * slide_dir


# ── Rev-Rev Joint ─────────────────────────────────────────────────

class RevRevJoint(Joint):
    """Constant-length link between two revolute joints.

    Attributes:
        L: Length of the link.
        q0: Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, L=0, q0=0,
                 iBody=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.L = L
        self.q0 = q0

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 2
        self.iBody = self.iMarker.body
        self.jBody = self.jMarker.body

    def compute_phi(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        L = self.L
        u = d / L
        return (u.T @ d - L) / 2

    def compute_jacobian_i(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        u = d / self.L
        return np.block([u.T, (u.T @ iPt._sPr).reshape(1, 1)])

    def compute_jacobian_j(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        u = d / self.L
        return np.block([-u.T, -(u.T @ jPt._sPr).reshape(1, 1)])

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        Bi = self.iBody
        Bj = self.jBody
        d = iPt._rP - jPt._rP
        dd = iPt._drP - jPt._drP
        L = self.L
        u = d / L
        ud = dd / L
        f = -ud.T @ dd
        if Bi is Ground:
            f = f + u.T @ s_rot(jPt._dsP) * Bj.angular_velocity
        elif Bj is Ground:
            f = f - u.T @ s_rot(iPt._dsP) * Bi.angular_velocity
        else:
            f = f - u.T @ s_rot(
                iPt._dsP * Bi.angular_velocity -
                jPt._dsP * Bj.angular_velocity)
        return f

    def fk_step(self, parent, child):
        from .mechanics import A_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = A_matrix(parent.orientation)
        parent_sP = parent._rotation_matrix @ parent_marker.local_position.reshape(2, 1)
        r_joint = parent.position + parent_sP
        child._rotation_matrix = A_matrix(child.orientation)
        child_sP = child._rotation_matrix @ child_marker.local_position.reshape(2, 1)
        L = self.L
        r_child_pivot = r_joint + L * np.array([[np.cos(q0)], [np.sin(q0)]])
        child.position = r_child_pivot - child_sP


# ── Rev-Tran Joint ────────────────────────────────────────────────

class RevTranJoint(Joint):
    """Revolute-translational composite joint.

    Attributes:
        L: Initial length parameter.
        q0: Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, L=0, q0=0,
                 iBody=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.L = L
        self.q0 = q0

        # Validation: rev-tran requires orientation on iMarker
        if iMarker is not None and not iMarker.has_orientation:
            raise ValueError(
                "RevTranJoint requires iMarker to have orientation "
                "(theta != None)")

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 2
        self.iBody = self.iMarker.body
        self.jBody = self.jMarker.body

    def compute_phi(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        uir = self.iMarker._ur
        d = iPt._rP - jPt._rP
        return uir.T @ d - self.L

    def compute_jacobian_i(self, model):
        iPt = self.iMarker
        ui = self.iMarker._u
        ui_r = self.iMarker._ur
        d = iPt._rP - self.jMarker._rP
        return np.block([ui_r.T, (ui.T @ (iPt._sP - d)).reshape(1, 1)])

    def compute_jacobian_j(self, model):
        jPt = self.jMarker
        ui = self.iMarker._u
        ui_r = self.iMarker._ur
        return np.block([-ui_r.T, -(ui.T @ jPt._sP).reshape(1, 1)])

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        Bi = self.iBody
        Bj = self.jBody
        ui = self.iMarker._u
        ui_d = self.iMarker._du
        d = iPt._rP - jPt._rP
        dd = iPt._drP - jPt._drP
        if Bi is Ground:
            f = ui.T @ jPt._dsP * Bj.angular_velocity
        elif Bj is Ground:
            f = (ui_d.T @ (d * Bi.angular_velocity + 2 * s_rot(dd))
                 - ui.T @ iPt._dsP * Bi.angular_velocity)
        else:
            f = (ui_d.T @ (d * Bi.angular_velocity + 2 * s_rot(dd))
                 - ui.T @ (iPt._dsP * Bi.angular_velocity
                           - jPt._dsP * Bj.angular_velocity))
        return f


# ── Rigid Joint ───────────────────────────────────────────────────

class RigidJoint(Joint):
    """Rigid (weld) joint between two bodies.

    Attributes:
        d0: Initial displacement vector.
    """

    def __init__(self, *, iBody=None, jBody=None, d0=None,
                 iMarker=None, jMarker=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.d0 = d0 if d0 is not None else []
        self._p0 = 0

    def initialize(self, model):
        self._mrows = 3
        self._nbody = 2
        Bi = self.iBody
        Bj = self.jBody
        if Bi is Ground:
            self.d0 = -Bj._rotation_matrix.T @ Bj.position
            self._p0 = -Bj.orientation
        elif Bj is Ground:
            self.d0 = Bi.position
            self._p0 = Bi.orientation
        else:
            self.d0 = Bj._rotation_matrix.T @ (Bi.position - Bj.position)
            self._p0 = Bi.orientation - Bj.orientation

    def compute_phi(self, model):
        Bi = self.iBody
        Bj = self.jBody
        if Bi is Ground:
            return np.vstack([
                -(Bj.position + Bj._rotation_matrix @ self.d0),
                -Bj.orientation - self._p0
            ])
        elif Bj is Ground:
            return np.vstack([
                Bi.position - self.d0,
                Bi.orientation - self._p0
            ])
        else:
            return np.vstack([
                Bi.position - (Bj.position + Bj._rotation_matrix @ self.d0),
                Bi.orientation - Bj.orientation - self._p0
            ])

    def compute_jacobian_i(self, model):
        return np.eye(3)

    def compute_jacobian_j(self, model):
        Bj = self.jBody
        if Bj is not Ground:
            return np.block([
                [-np.eye(2), -s_rot(Bj._rotation_matrix @ self.d0)],
                [np.array([0, 0, -1])]
            ])
        return np.zeros((3, 3))

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        Bj = self.jBody
        f = np.zeros(3)
        if Bj is not Ground:
            f = np.concatenate([
                -Bj._rotation_matrix @ self.d0 * Bj.angular_velocity ** 2,
                np.array([0])
            ])
        return f


# ── Disc Joint ────────────────────────────────────────────────────

class DiscJoint(Joint):
    """Rolling disc (wheel) constraint.

    Attributes:
        R: Disc radius.
        x0: Initial x-position.
    """

    def __init__(self, *, iBody=None, R=1, x0=0,
                 iMarker=None, jMarker=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.R = R
        self.x0 = x0
        self._p0 = 0

    def initialize(self, model):
        self._mrows = 2
        self._nbody = 1

    def compute_phi(self, model):
        Bi = self.iBody
        return np.vstack([
            Bi.position[1] - self.R,
            (Bi.position[0] - self.x0) + self.R * (Bi.orientation - self._p0)
        ])

    def compute_jacobian_i(self, model):
        return np.array([
            [0, 1, 0],
            [1, 0, self.R]
        ])

    def compute_jacobian_j(self, model):
        return np.zeros((2, 3))

    def compute_rhs_velocity(self, model):
        return None

    def compute_rhs_acceleration(self, model):
        return np.zeros(2)


# ── Relative Rotation Joint ──────────────────────────────────────

class RelRotJoint(Joint):
    """Relative rotation driven by an analytical function.

    Attributes:
        iFunct: Function object driving the constraint.
    """

    def __init__(self, *, iBody=None, jBody=None, iFunct=None,
                 iMarker=None, jMarker=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.iFunct = iFunct

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 1

    def compute_phi(self, model):
        fun, fun_d, fun_dd = functEval(self.iFunct, model.t)
        Bi = self.iBody
        Bj = self.jBody
        if Bi is Ground:
            return -Bj.orientation - fun
        elif Bj is Ground:
            return Bi.orientation - fun
        else:
            return Bi.orientation - Bj.orientation - fun

    def compute_jacobian_i(self, model):
        return np.array([[0, 0, 1]])

    def compute_jacobian_j(self, model):
        return np.array([[0, 0, -1]])

    def compute_rhs_velocity(self, model):
        fun, fun_d, _ = functEval(self.iFunct, model.t)
        return fun_d

    def compute_rhs_acceleration(self, model):
        fun, fun_d, fun_dd = functEval(self.iFunct, model.t)
        return fun_dd


# ── Relative Translation Joint ───────────────────────────────────

class RelTranJoint(Joint):
    """Relative translation driven by an analytical function.

    Attributes:
        iFunct: Function object driving the constraint.
    """

    def __init__(self, *, iMarker=None, jMarker=None, iFunct=None,
                 iBody=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.iFunct = iFunct

    def initialize(self, model):
        self._mrows = 1
        self._nbody = 1

    def compute_phi(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        fun, fun_d, fun_dd = functEval(self.iFunct, model.t)
        return (d.T @ d - fun ** 2) / 2

    def compute_jacobian_i(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        return np.block([d.T, (d.T @ iPt._sPr).reshape(1, 1)])

    def compute_jacobian_j(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        d = iPt._rP - jPt._rP
        return np.block([-d.T, -(d.T @ jPt._sPr).reshape(1, 1)])

    def compute_rhs_velocity(self, model):
        fun, fun_d, _ = functEval(self.iFunct, model.t)
        return fun * fun_d

    def compute_rhs_acceleration(self, model):
        iPt = self.iMarker
        jPt = self.jMarker
        Bi = self.iBody
        Bj = self.jBody
        d = iPt._rP - jPt._rP
        dd = iPt._drP - jPt._drP
        fun, fun_d, fun_dd = functEval(self.iFunct, model.t)
        f = fun * fun_dd + fun_d ** 2
        if Bi is Ground:
            f = f + d.T @ s_rot(jPt._dsP).T @ Bj.angular_velocity
        elif Bj is Ground:
            f = f - d.T @ s_rot(iPt._dsP).T @ Bi.angular_velocity - dd.T @ dd
        else:
            f = (f + d.T @ s_rot(jPt._dsP).T @ Bj.angular_velocity
                 - d.T @ s_rot(iPt._dsP).T @ Bi.angular_velocity - dd.T @ dd)
        return f


# ═══════════════════════════════════════════════════════════════════
#  FORCE HIERARCHY
# ═══════════════════════════════════════════════════════════════════

class Force(ABC, Base):
    """Abstract base class for all force elements.

    Every concrete force must implement the ``apply`` method.

    Attributes:
        iMarker: Marker on body i (or None).
        jMarker: Marker on body j (or None).
        iBody: Body i.
        jBody: Body j.
    """

    def __init__(self, *, iMarker=None, jMarker=None,
                 iBody=None, jBody=None):
        super().__init__()
        self.iMarker = iMarker
        self.jMarker = jMarker
        self.iBody = iMarker.body if iMarker is not None else (iBody or Ground)
        self.jBody = jMarker.body if jMarker is not None else (jBody or Ground)

    @abstractmethod
    def apply(self, bodies):
        """Apply this force to the relevant bodies.

        Args:
            bodies: List of Body objects in the model.
        """

    def __repr__(self):
        return f"{self.__class__.__name__}(iBody={self.iBody!r}, jBody={self.jBody!r})"


# ── Weight ────────────────────────────────────────────────────────

class Weight(Force):
    """Gravitational force applied to all bodies.

    Attributes:
        gravity: Gravitational acceleration magnitude (m/s^2).
        gravity_direction: Unit direction of gravity, shape (2, 1).
    """

    DEFAULT_GRAVITY = 9.81
    DEFAULT_GRAVITY_VECTOR = colvect([0, -1])

    def __init__(self, gravity=None, gravity_direction=None):
        super().__init__()
        self.gravity = gravity if gravity is not None else self.DEFAULT_GRAVITY
        self.gravity_direction = (gravity_direction if gravity_direction is not None
                                  else self.DEFAULT_GRAVITY_VECTOR.copy())
        self._initialized = False

    def initialize_weights(self, bodies):
        """Pre-compute weight vectors for all bodies.

        Called once by PlanarMultibodyModel._initialize().

        Args:
            bodies: List of Body objects.
        """
        ug = self.gravity * self.gravity_direction
        for body in bodies:
            body._weight = body.mass * ug
        self._initialized = True

    def apply(self, bodies):
        for body in bodies:
            body._force += body._weight

    def __repr__(self):
        return f"Weight(gravity={self.gravity})"


# ── Point-to-Point Force ─────────────────────────────────────────

class PtpForce(Force):
    """Point-to-point spring-damper-actuator force.

    Attributes:
        k: Spring stiffness (N/m).
        L0: Undeformed spring length (m).
        dc: Damping coefficient (N*s/m).
        f_a: Constant actuator force (N).
    """

    def __init__(self, *, iMarker=None, jMarker=None, k=0, L0=0, dc=0,
                 f_a=0):
        super().__init__(iMarker=iMarker, jMarker=jMarker)
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
        return f"PtpForce(k={self.k}, L0={self.L0}, dc={self.dc})"


# ── Rotational Spring-Damper-Actuator ─────────────────────────────

class RotSdaForce(Force):
    """Rotational spring-damper-actuator torque.

    Attributes:
        k: Torsional stiffness (N*m/rad).
        theta0: Undeformed angle (rad).
        dc: Damping coefficient (N*m*s/rad).
        T_a: Constant actuator torque (N*m).
    """

    def __init__(self, *, iBody=None, jBody=None, k=0, theta0=0, dc=0,
                 T_a=0, iMarker=None, jMarker=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
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

    Attributes:
        force_local: Force vector in local coordinates, shape (2, 1).
    """

    def __init__(self, *, iBody=None, force_local=None, iMarker=None,
                 jMarker=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.force_local = (force_local if force_local is not None
                            else colvect([0, 0]))

    # Ensure stored as column vector
    force_local = as_column_property("force_local")

    def apply(self, bodies):
        self.iBody._force += self.iBody._rotation_matrix @ self.force_local


# ── Global Force ──────────────────────────────────────────────────

class GlobalForce(Force):
    """Constant force in the global reference frame.

    Attributes:
        force_global: Force vector in global coordinates, shape (2, 1).
    """

    def __init__(self, *, iBody=None, force_global=None, iMarker=None,
                 jMarker=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.force_global = (force_global if force_global is not None
                             else colvect([0, 0]))

    # Ensure stored as column vector
    force_global = as_column_property("force_global")

    def apply(self, bodies):
        self.iBody._force += self.force_global


# ── Torque ────────────────────────────────────────────────────────

class Torque(Force):
    """Constant torque in the global reference frame.

    Attributes:
        torque_value: Scalar torque value (N*m).
    """

    def __init__(self, *, iBody=None, torque_value=0, iMarker=None,
                 jMarker=None, jBody=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.torque_value = torque_value

    def apply(self, bodies):
        self.iBody._torque += self.torque_value


# ── User Force ────────────────────────────────────────────────────

class UserForce(Force):
    """User-defined force via callback function.

    The callback is called with no arguments and should return a list of
    dicts, each with keys ``'body'``, ``'force'`` (2-element), and
    ``'torque'`` (scalar).

    Attributes:
        callback: Callable returning list of force contribution dicts.
    """

    def __init__(self, *, callback=None, iMarker=None, jMarker=None,
                 iBody=None, jBody=None, **kwargs):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody)
        self.callback = callback
        # Store any extra kwargs for user convenience (e.g. k, L0, dc)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def apply(self, bodies):
        if self.callback is not None and callable(self.callback):
            result = self.callback()
            if result is not None:
                for entry in result:
                    body = entry['body']
                    force = entry.get('force', None)
                    torque = entry.get('torque', None)
                    if force is not None:
                        f = np.asarray(force, dtype=float).reshape(2, 1)
                        body._force += f
                    if torque is not None:
                        body._torque += float(torque)

    def __repr__(self):
        return f"UserForce(callback={self.callback!r})"


# ═══════════════════════════════════════════════════════════════════
#  FUNCTION CLASS (unchanged)
# ═══════════════════════════════════════════════════════════════════

class Function(Base):
    """Analytical function for driving constraint equations.

    Attributes:
        type: Function type ``'a'``, ``'b'``, or ``'c'``.
        t_start: Start time for time-based functions.
        f_start: Start value for time-based functions.
        t_end: End time for time-based functions.
        f_end: End value for time-based functions.
        dfdt_end: Max first derivative value (type ``'c'``).
        ncoeff: Number of coefficients.
        coeff: Coefficient array.
    """

    def __init__(self, type='a', t_start=0, f_start=0, t_end=1, f_end=1,
                 dfdt_end=1, ncoeff=4, coeff=None):
        super().__init__()
        self.type = type
        self.t_start = t_start
        self.f_start = f_start
        self.t_end = t_end
        self.f_end = f_end
        self.dfdt_end = dfdt_end
        self.ncoeff = ncoeff
        if coeff is not None:
            c = np.asarray(coeff, dtype=float).flatten()
            padded = np.zeros(9)
            padded[:min(len(c), 9)] = c[:min(len(c), 9)]
            self.coeff = padded
        else:
            self.coeff = np.zeros(9)

    def __repr__(self):
        return f"Function(type='{self.type}')"
