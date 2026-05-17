"""Geometric joint constraints.

This module contains the joint hierarchy: abstract :class:`Joint` base and
all scleronomic (time-independent) concrete joints.  Time-driven kinematic
constraints (motions) live in :mod:`pmd.core.motion`; force elements live
in :mod:`pmd.core.forces`.

Author: Giacomo Cangi
"""

import logging
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import *

from .utils import *
from .model import Base, Ground
from .mechanics import rotate_90

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
#  JOINT HIERARCHY
# ═══════════════════════════════════════════════════════════════════

class Joint(ABC, Base):
    """Abstract base class for all joint constraints.

    Every concrete joint must implement the abstract methods below.

    Attributes
    ----------
    iMarker : Marker or None
        Marker on body i.
    jMarker : Marker or None
        Marker on body j.
    iBody : Body
        Body i (derived from iMarker or explicit).
    jBody : Body
        Body j (derived from jMarker or explicit).
    name : str or None
        Optional human-readable name for identification and legends.
    """

    def __init__(self, *, iMarker=None, jMarker=None,
                 iBody=None, jBody=None, name=None):
        """Initialize a Joint.

        Parameters
        ----------
        iMarker : Marker, optional
            Marker on body i.
        jMarker : Marker, optional
            Marker on body j.
        iBody : Body, optional
            Explicit body i (used if iMarker is None).
        jBody : Body, optional
            Explicit body j (used if jMarker is None).
        name : str, optional
            Human-readable name.
        """
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
        self._result_container = None
        self.name = name

    # ── Abstract interface ────────────────────────────────────────

    @abstractmethod
    def initialize(self, model):
        """Compute derived quantities after model construction.

        Called once by PlanarMultibodyModel._initialize().

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.
        """

    @abstractmethod
    def compute_phi(self, model):
        """Evaluate constraint equations.

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.

        Returns
        -------
        ndarray
            Constraint residual, shape (mrows, 1).
        """

    @abstractmethod
    def compute_jacobian_i(self, model):
        """Compute the Jacobian block for body i.

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.

        Returns
        -------
        ndarray
            Jacobian block for body i, shape (mrows, 3).
        """

    @abstractmethod
    def compute_jacobian_j(self, model):
        """Compute the Jacobian block for body j.

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.

        Returns
        -------
        ndarray
            Jacobian block for body j, shape (mrows, 3).
        """

    @abstractmethod
    def compute_rhs_velocity(self, model):
        """Compute right-hand side of velocity constraint.

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.

        Returns
        -------
        ndarray or None
            RHS velocity, shape (mrows, 1), or None if zero.
        """

    @abstractmethod
    def compute_rhs_acceleration(self, model):
        """Compute right-hand side of acceleration constraint (gamma).

        Parameters
        ----------
        model : PlanarMultibodyModel
            The model instance.

        Returns
        -------
        ndarray
            RHS acceleration, shape (mrows, 1).
        """

    def fk_step(self, parent, child):
        """Forward kinematics: place child given parent state and q0.

        Default implementation does nothing (for joints without FK role).
        Override in subclasses that participate in assembly FK.

        Parameters
        ----------
        parent : Body
            The parent body.
        child : Body
            The child body.
        """

    def get_reaction(self, index=None):
        """Get reaction force time history for a single constraint row.

        Parameters
        ----------
        index : int, optional
            0-based local row index within this joint's constraint rows.
            If None and the joint has exactly 1 row, returns that row.

        Returns
        -------
        NDArray
            1-D array of Lagrange multipliers over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        IndexError
            If index is out of range or ambiguous.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        reactions = self._result_container['reactions']
        if index is None:
            if reactions.shape[1] == 1:
                return reactions[:, 0]
            raise IndexError(
                f"Joint has {reactions.shape[1]} constraint rows; specify index.")
        return reactions[:, index]

    def get_reactions(self):
        """Get all reaction force time histories for this joint.

        Returns
        -------
        NDArray
            Array of shape (nSteps, mrows) with Lagrange multipliers.

        Raises
        ------
        RuntimeError
            If no results are available.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return self._result_container['reactions']

    def __repr__(self):
        label = f"'{self.name}'" if self.name else ""
        body_info = f"iBody={self.iBody!r}, jBody={self.jBody!r}"
        if label:
            return f"{self.__class__.__name__}({label}, {body_info})"
        return f"{self.__class__.__name__}({body_info})"


# ── Revolute Joint ────────────────────────────────────────────────

class RevJoint(Joint):
    """Revolute (pin) joint between two bodies.

    Attributes
    ----------
    fix : int
        If 1, the relative rotation is also constrained.
    q0 : float
        Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, fix=0, q0=0,
                 iBody=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
            f = rotate_90(jPt._dsP) * Bj.angular_velocity
        elif Bj is Ground:
            f = -rotate_90(iPt._dsP) * Bi.angular_velocity
        else:
            f = (-rotate_90(iPt._dsP) * Bi.angular_velocity
                 + rotate_90(jPt._dsP) * Bj.angular_velocity)
        if self.fix == 1:
            f = np.vstack([f, [0]])
        return f

    def fk_step(self, parent, child):
        from .mechanics import rotation_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = rotation_matrix(parent.orientation)
        parent_sP = parent._rotation_matrix @ parent_marker.local_position.reshape(2, 1)
        r_joint = parent.position + parent_sP
        child.orientation = parent.orientation + q0
        child._rotation_matrix = rotation_matrix(child.orientation)
        child_sP = child._rotation_matrix @ child_marker.local_position.reshape(2, 1)
        child.position = r_joint - child_sP


# ── Translational Joint ──────────────────────────────────────────

class TranJoint(Joint):
    """Translational (prismatic) joint between two bodies.

    Attributes
    ----------
    fix : int
        If 1, the relative translation is also constrained.
    q0 : float
        Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, fix=0, q0=0,
                 iBody=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
        ujdr = rotate_90(ujd)

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
                f3 += (u.T @ (rotate_90(jPt._dsP) * Bj.angular_velocity)).item()
            elif Bj is Ground:
                f3 -= (u.T @ (rotate_90(iPt._dsP) * Bi.angular_velocity)).item()
            else:
                term1 = iPt._dsP * Bi.angular_velocity
                term2 = jPt._dsP * Bj.angular_velocity
                f3 -= (u.T @ rotate_90(term1 - term2)).item()

            f = np.vstack([f, [[f3]]])

        return f

    def fk_step(self, parent, child):
        from .mechanics import rotation_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = rotation_matrix(parent.orientation)
        child.orientation = parent.orientation
        child._rotation_matrix = rotation_matrix(child.orientation)
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

    Attributes
    ----------
    L : float
        Length of the link.
    q0 : float
        Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, L=0, q0=0,
                 iBody=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
            f = f + u.T @ rotate_90(jPt._dsP) * Bj.angular_velocity
        elif Bj is Ground:
            f = f - u.T @ rotate_90(iPt._dsP) * Bi.angular_velocity
        else:
            f = f - u.T @ rotate_90(
                iPt._dsP * Bi.angular_velocity -
                jPt._dsP * Bj.angular_velocity)
        return f

    def fk_step(self, parent, child):
        from .mechanics import rotation_matrix
        q0 = self.q0
        im = self.iMarker
        jm = self.jMarker
        if id(im.body) == id(parent):
            parent_marker, child_marker = im, jm
        else:
            parent_marker, child_marker = jm, im
        parent._rotation_matrix = rotation_matrix(parent.orientation)
        parent_sP = parent._rotation_matrix @ parent_marker.local_position.reshape(2, 1)
        r_joint = parent.position + parent_sP
        child._rotation_matrix = rotation_matrix(child.orientation)
        child_sP = child._rotation_matrix @ child_marker.local_position.reshape(2, 1)
        L = self.L
        r_child_pivot = r_joint + L * np.array([[np.cos(q0)], [np.sin(q0)]])
        child.position = r_child_pivot - child_sP


# ── Rev-Tran Joint ────────────────────────────────────────────────

class RevTranJoint(Joint):
    """Revolute-translational composite joint.

    Attributes
    ----------
    L : float
        Initial length parameter.
    q0 : float
        Initial DOF value used by the assembler.
    """

    def __init__(self, *, iMarker=None, jMarker=None, L=0, q0=0,
                 iBody=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
            f = (ui_d.T @ (d * Bi.angular_velocity + 2 * rotate_90(dd))
                 - ui.T @ iPt._dsP * Bi.angular_velocity)
        else:
            f = (ui_d.T @ (d * Bi.angular_velocity + 2 * rotate_90(dd))
                 - ui.T @ (iPt._dsP * Bi.angular_velocity
                           - jPt._dsP * Bj.angular_velocity))
        return f


# ── Rigid Joint ───────────────────────────────────────────────────

class RigidJoint(Joint):
    """Rigid (weld) joint between two bodies.

    Attributes
    ----------
    d0 : ndarray
        Initial displacement vector.
    """

    def __init__(self, *, iBody=None, jBody=None, d0=None,
                 iMarker=None, jMarker=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
                [-np.eye(2), -rotate_90(Bj._rotation_matrix @ self.d0)],
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

    Attributes
    ----------
    R : float
        Disc radius.
    x0 : float
        Initial x-position.
    """

    def __init__(self, *, iBody=None, R=1, x0=0,
                 iMarker=None, jMarker=None, jBody=None, name=None):
        super().__init__(iMarker=iMarker, jMarker=jMarker,
                         iBody=iBody, jBody=jBody, name=name)
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
