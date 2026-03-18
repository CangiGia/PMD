"""
Planar Multi-Body Dynamics Model Components

This module provides the core model components: Ground (singleton),
Base (counter), Body (rigid body), and Marker (body-fixed reference frame).

Author: Giacomo Cangi
"""


import warnings
import numpy as np
from numpy.typing import *
from .utils import *
from .mechanics import s_rot


# ── Base class ────────────────────────────────────────────────────

class Base:
    """
    Base class for all multi-body simulation objects.

    This class provides a foundation for objects in the multi-body dynamic 
    simulation framework.
    """
    COUNT = 0  # class variable to track the number of instances

    def __new__(cls, *args, **kwargs):
        """
        Create a new instance of a subclass of `Base` and increment the 
        instance count.

        Parameters
        ----------
        args (tuple)
            Positional arguments passed to the constructor.
        kwargs (dict)
            Keyword arguments passed to the constructor.

        Returns
        -------
        instance
            A new instance of the class calling this method.
        """
        instance = super().__new__(cls)
        cls.COUNT += 1
        return instance

    @classmethod
    def get_count(cls):
        """
        Retrieve the current count of instances.

        Returns
        -------
        int
            The total count of instances derived from `Base`.
        """
        return cls.COUNT
    
    def get_type(self):
        """
        Get the name of the class of the instance.

        Returns
        -------
        str
            The class name of the instance.
        """
        return self.__class__.__name__


# ── Marker class ──────────────────────────────────────────────────

class Marker:
    def __init__(self, body, position, theta=None, name=None):
        pos = np.asarray(position, dtype=float).flatten()
        if pos.shape != (2,):
            raise ValueError(
                f"Marker position must be a 2-element array, got shape {pos.shape}")

        if theta is not None and not isinstance(theta, (int, float, np.floating)):
            raise TypeError(
                f"Marker theta must be a float or None, got {type(theta)}")

        self.body     = body
        self.position = pos        # flat (2,) array, local coords [xi, eta]
        self.theta    = theta
        self.name     = name

        # solver-internal: position kinematics
        self._sP   = np.zeros((2, 1))
        self._sPr  = np.zeros((2, 1))
        self._rP   = np.zeros((2, 1))
        self._dsP  = np.zeros((2, 1))
        self._drP  = np.zeros((2, 1))
        self._ddrP = np.zeros((2, 1))

        # solver-internal: orientation kinematics (only when theta is given)
        if theta is not None:
            self._ulocal = np.array([[np.cos(theta)], [np.sin(theta)]])
            self._u  = np.zeros((2, 1))
            self._ur = np.zeros((2, 1))
            self._du = np.zeros((2, 1))

    @property
    def has_orientation(self):
        return self.theta is not None

    @property
    def global_position(self):
        return self._rP

    @property
    def global_velocity(self):
        return self._drP


# ── Ground singleton ──────────────────────────────────────────────

class _GroundType:
    """Singleton representing the inertial ground/world frame.

    Ground is the immovable reference body with zero state.  Use the
    module-level ``Ground`` instance instead of instantiating this class.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            # auto-create an origin marker
            origin = Marker(body=cls._instance, position=[0, 0], theta=0.0, name='origin')
            pos_col = origin.position.reshape(2, 1)
            origin._sP = pos_col.copy()
            origin._sPr = s_rot(pos_col)
            origin._rP = pos_col.copy()
            origin._u = origin._ulocal.copy()
            origin._ur = s_rot(origin._ulocal)
            cls._markers = [origin]
            cls._instance.origin = origin
        return cls._instance

    # Fixed zero state — mirrors the Body interface that the solver reads.
    r   = np.zeros((2, 1))
    p   = 0.0
    dr  = np.zeros((2, 1))
    dp  = 0.0
    ddr = np.zeros((2, 1))
    ddp = 0.0
    _A  = np.eye(2)
    _bidx = 0          # internal body index: 0 ≡ ground

    _markers = []

    def add_marker(self, position, theta=None, name=None):
        marker = Marker(body=self, position=position, theta=theta, name=name)
        pos_col = marker.position.reshape(2, 1)
        marker._sP = pos_col.copy()
        marker._sPr = s_rot(pos_col)
        marker._rP = pos_col.copy()
        if marker.has_orientation:
            marker._u = marker._ulocal.copy()
            marker._ur = s_rot(marker._ulocal)
        self._markers.append(marker)
        return marker

    def __repr__(self):
        return "Ground"

    def __bool__(self):
        """Ground is falsy so ``if body:`` is False for Ground, True for Body."""
        return False


Ground = _GroundType()


# ── Body class ────────────────────────────────────────────────────

class Body(Base):
    """
    Create a `rigid body` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    m (float) 
        Mass of the body.
    J (float) 
        Moment of inertia of the body.
    r (NDArray) 
        Position vector (x, y) of the body.
    p (float) 
        Orientation angle (phi) of the body.
    dr (NDArray) 
        Time derivative of position (x, y).
    dp (float) 
        Time derivative of the orientation angle (phi).
    ddr (NDArray) 
        Second derivative of position (x, y).
    ddp (float) 
        Second time derivative of the orientation angle (phi).

    Notes
    -----
    _A (NDArray) 
        Rotational transformation matrix.
    _bidx (int)
        Internal body index (1-based), set by the solver during initialization.
    _irc (int) 
        Index of the first element of r in u.
    _irv (int) 
        Index of the first element of r_dot in u.
    _ira (int) 
        Index of the first element of r_dot2 in v.
    _invm (float) 
        Inverse of the mass.
    _invJ (float) 
        Inverse of the moment of inertia.
    _wgt (NDArray) 
        Weight of the body as a force vector.
    _f (NDArray) 
        Sum of forces acting on the body.
    _n (float) 
        Sum of moments acting on the body.
    _markers (list) 
        List of Marker objects associated with this body.
    """

    def __init__(self, m=1, J=1, r=None, p=None, dr=None, dp=0, ddr=None, ddp=0):
        super().__init__()  # call to the Base class constructor
        if m <= 0:
            raise ValueError(f"Body {self.COUNT}: mass must be positive, got {m}")
        if J < 0:
            raise ValueError(f"Body {self.COUNT}: moment of inertia cannot be negative, got {J}")
        self.m = m
        self.J = J
        self.r = r if r is not None else colvect(0, 0)
        self.p = p if p is not None else 0.0
        # Set the _given flags AFTER initial assignment to avoid
        # the setter marking defaults as user-provided.
        self._r_given = r is not None
        self._p_given = p is not None
        self.dr = dr if dr is not None else colvect(0, 0)
        self.dp = dp
        self.ddr = ddr if ddr is not None else colvect(0, 0)
        self.ddp = ddp
        self._A = np.eye(2)
        self._bidx = 0
        self._irc = 0
        self._irv = 0
        self._ira = 0
        self._invm = 1.0 / m
        self._invJ = 1.0 / J if J != 0 else float('inf')
        self._wgt = colvect(0, 0)
        self._f = colvect(0, 0)
        self._n = 0
        self._markers = []

    # new values will be automatically treated as column vectors, 
    # even if the user defines them as row vectors or lists!
    dr = as_column_property("dr")
    ddr = as_column_property("ddr")

    @property
    def r(self):
        return self.__r

    @r.setter
    def r(self, value):
        if isinstance(value, list):
            value = np.array(value)
        if isinstance(value, np.ndarray):
            if value.ndim == 1:
                value = value.reshape(-1, 1)
            elif value.ndim == 2 and value.shape[0] == 1:
                value = value.T
        self.__r = value
        # Mark as explicitly provided so auto-assembly won't overwrite it
        if hasattr(self, '_r_given'):
            self._r_given = True

    @property
    def p(self):
        return self.__p

    @p.setter
    def p(self, value):
        if isinstance(value, np.ndarray):
            value = float(value.flat[0])
        self.__p = value
        if hasattr(self, '_p_given'):
            self._p_given = True

    def add_marker(self, position, theta=None, name=None):
        marker = Marker(body=self, position=position, theta=theta, name=name)
        # V3: warn if duplicate position on same body
        for existing in self._markers:
            if np.linalg.norm(existing.position - marker.position) < 1e-10:
                warnings.warn(
                    f"Body already has a marker at position {marker.position}")
                break
        self._markers.append(marker)
        return marker

    def add_marker_at(self, reference_marker, offset=(0, 0), theta=None, name=None):
        marker = Marker(body=self, position=[0, 0], theta=theta, name=name)
        marker._deferred_ref = reference_marker
        marker._deferred_offset = offset
        self._markers.append(marker)
        return marker
