"""Planar Multi-Body Dynamics Model Components.

This module provides the core model components: Ground (singleton),
Base (counter), Body (rigid body), and Marker (body-fixed reference frame).

Author: Giacomo Cangi
"""

import logging
import warnings

import numpy as np
from numpy.typing import *

from .utils import *
from .mechanics import rotate_90

logger = logging.getLogger(__name__)


# ── Base class ────────────────────────────────────────────────────

class Base:
    """Base class for all multi-body simulation objects.

    Provides automatic instance counting for each subclass.

    Attributes:
        COUNT: Class variable tracking the number of instances created.
    """

    COUNT = 0

    def __new__(cls, *args, **kwargs):
        """Create a new instance and increment the instance count.

        Parameters
        ----------
        *args : tuple
            Positional arguments passed to the constructor.
        **kwargs : dict
            Keyword arguments passed to the constructor.

        Returns
        -------
        object
            A new instance of the class.
        """
        instance = super().__new__(cls)
        cls.COUNT += 1
        return instance

    @classmethod
    def get_count(cls):
        """Retrieve the current count of instances.

        Returns
        -------
        int
            The total count of instances.
        """
        return cls.COUNT

    def get_type(self):
        """Get the class name of this instance.

        Returns
        -------
        str
            The class name.
        """
        return self.__class__.__name__


# ── Marker class ──────────────────────────────────────────────────

class Marker:
    """Body-fixed reference frame (point + optional orientation).

    Attributes:
        body: The body this marker is attached to.
        local_position: Flat (2,) array of local coordinates [xi, eta].
        theta: Optional orientation angle in the local frame.
        name: Optional human-readable name.
    """

    def __init__(self, body, local_position, theta=None, name=None):
        """Initialize a Marker.

        Parameters
        ----------
        body : Body or _GroundType
            The body this marker is attached to.
        local_position : array_like
            2-element array of local coordinates.
        theta : float, optional
            Orientation angle (radians) in the local frame.
        name : str, optional
            Human-readable name.

        Raises
        ------
        ValueError
            If local_position is not a 2-element array.
        TypeError
            If theta is not a float or None.
        """
        pos = np.asarray(local_position, dtype=float).flatten()
        if pos.shape != (2,):
            raise ValueError(
                f"Marker position must be a 2-element array, got shape {pos.shape}")

        if theta is not None and not isinstance(theta, (int, float, np.floating)):
            raise TypeError(
                f"Marker theta must be a float or None, got {type(theta)}")

        self.body = body
        self.local_position = pos
        self.theta = theta
        self.name = name

        # solver-internal: position kinematics
        self._sP = np.zeros((2, 1))
        self._sPr = np.zeros((2, 1))
        self._rP = np.zeros((2, 1))
        self._dsP = np.zeros((2, 1))
        self._drP = np.zeros((2, 1))
        self._ddrP = np.zeros((2, 1))

        # solver-internal: orientation kinematics (only when theta is given)
        if theta is not None:
            self._ulocal = np.array([[np.cos(theta)], [np.sin(theta)]])
            self._u = np.zeros((2, 1))
            self._ur = np.zeros((2, 1))
            self._du = np.zeros((2, 1))

    @property
    def has_orientation(self):
        """bool: True if this marker has an orientation (theta is not None)."""
        return self.theta is not None

    @property
    def global_position(self):
        """ndarray: Global position of this marker, shape (2, 1)."""
        return self._rP

    @property
    def global_velocity(self):
        """ndarray: Global velocity of this marker, shape (2, 1)."""
        return self._drP

    def __repr__(self):
        label = f"'{self.name}'" if self.name else "unnamed"
        body_label = repr(self.body) if hasattr(self.body, '__repr__') else str(self.body)
        return f"Marker({label}, body={body_label}, local_position={self.local_position})"


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
            origin = Marker(body=cls._instance, local_position=[0, 0],
                            theta=0.0, name='origin')
            pos_col = origin.local_position.reshape(2, 1)
            origin._sP = pos_col.copy()
            origin._sPr = rotate_90(pos_col)
            origin._rP = pos_col.copy()
            origin._u = origin._ulocal.copy()
            origin._ur = rotate_90(origin._ulocal)
            cls._markers = [origin]
            cls._instance.origin = origin
        return cls._instance

    # Fixed zero state — mirrors the Body interface that the solver reads.
    position = np.zeros((2, 1))
    orientation = 0.0
    velocity = np.zeros((2, 1))
    angular_velocity = 0.0
    acceleration = np.zeros((2, 1))
    angular_acceleration = 0.0
    _rotation_matrix = np.eye(2)
    _body_index = 0

    _markers = []

    def add_marker(self, local_position, theta=None, name=None):
        """Create and attach a new Marker to Ground.

        Parameters
        ----------
        local_position : array_like
            2-element array of local coordinates.
        theta : float, optional
            Orientation angle (radians).
        name : str, optional
            Human-readable name.

        Returns
        -------
        Marker
            The newly created marker.
        """
        marker = Marker(body=self, local_position=local_position,
                        theta=theta, name=name)
        pos_col = marker.local_position.reshape(2, 1)
        marker._sP = pos_col.copy()
        marker._sPr = rotate_90(pos_col)
        marker._rP = pos_col.copy()
        if marker.has_orientation:
            marker._u = marker._ulocal.copy()
            marker._ur = rotate_90(marker._ulocal)
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
    """Rigid body in a planar multi-body dynamic simulation.

    Attributes:
        name: Optional human-readable name.
        mass: Mass of the body (kg).
        inertia: Moment of inertia of the body (kg*m^2).
        position: Position vector (x, y), shape (2, 1).
        orientation: Orientation angle phi (rad).
        velocity: Velocity vector (dx, dy), shape (2, 1).
        angular_velocity: Angular velocity dphi (rad/s).
        acceleration: Acceleration vector (ddx, ddy), shape (2, 1).
        angular_acceleration: Angular acceleration ddphi (rad/s^2).
        name: Optional human-readable name.
        shape: Optional shape descriptor for visualization.
    """

    def __init__(self, mass=1, inertia=1, position=None, orientation=None,
                 velocity=None, angular_velocity=0, acceleration=None,
                 angular_acceleration=0, name=None, shape=None):
        """Initialize a Body.

        Parameters
        ----------
        mass : float
            Mass (must be positive).
        inertia : float
            Moment of inertia (must be non-negative).
        position : array_like, optional
            Initial position [x, y]. Defaults to [0, 0].
        orientation : float, optional
            Initial orientation angle phi (rad). Defaults to 0.
        velocity : array_like, optional
            Initial velocity [dx, dy]. Defaults to [0, 0].
        angular_velocity : float
            Initial angular velocity dphi (rad/s).
        acceleration : array_like, optional
            Initial acceleration [ddx, ddy]. Defaults to [0, 0].
        angular_acceleration : float
            Initial angular acceleration ddphi (rad/s²).
        name : str, optional
            Human-readable name.
        shape : Rectangle or Circle or Polygon, optional
            Shape descriptor for visualization.

        Raises
        ------
        ValueError
            If mass is not positive or inertia is negative.
        """
        super().__init__()
        if mass <= 0:
            raise ValueError(
                f"Body {self.COUNT}: mass must be positive, got {mass}")
        if inertia < 0:
            raise ValueError(
                f"Body {self.COUNT}: moment of inertia cannot be negative, "
                f"got {inertia}")

        self.name = name
        self.shape = shape
        self.mass = mass
        self.inertia = inertia
        self.position = position if position is not None else colvect(0, 0)
        self.orientation = orientation if orientation is not None else 0.0
        # Set the _given flags AFTER initial assignment to avoid
        # the setter marking defaults as user-provided.
        self._position_given = position is not None
        self._orientation_given = orientation is not None
        self.velocity = velocity if velocity is not None else colvect(0, 0)
        self.angular_velocity = angular_velocity
        self.acceleration = (acceleration if acceleration is not None
                             else colvect(0, 0))
        self.angular_acceleration = angular_acceleration
        self._rotation_matrix = np.eye(2)
        self._body_index = 0
        self._index_position = 0
        self._index_velocity = 0
        self._index_acceleration = 0
        self._inv_mass = 1.0 / mass
        self._inv_inertia = 1.0 / inertia if inertia != 0 else float('inf')
        self._weight = colvect(0, 0)
        self._force = colvect(0, 0)
        self._torque = 0
        self._markers = []
        self._result_container = None

    # Ensure assigned values are stored as column vectors
    velocity = as_column_property("velocity")
    acceleration = as_column_property("acceleration")

    @property
    def position(self):
        """ndarray: Position vector (x, y), shape (2, 1)."""
        return self.__position

    @position.setter
    def position(self, value):
        if isinstance(value, list):
            value = np.array(value)
        if isinstance(value, np.ndarray):
            if value.ndim == 1:
                value = value.reshape(-1, 1)
            elif value.ndim == 2 and value.shape[0] == 1:
                value = value.T
        self.__position = value
        if hasattr(self, '_position_given'):
            self._position_given = True

    @property
    def orientation(self):
        """float: Orientation angle phi (rad)."""
        return self.__orientation

    @orientation.setter
    def orientation(self, value):
        if isinstance(value, np.ndarray):
            value = float(value.flat[0])
        self.__orientation = value
        if hasattr(self, '_orientation_given'):
            self._orientation_given = True

    def add_marker(self, local_position, theta=None, name=None):
        """Create and attach a new Marker to this body.

        Parameters
        ----------
        local_position : array_like
            2-element array of local coordinates.
        theta : float, optional
            Orientation angle (radians) in the local frame.
        name : str, optional
            Human-readable name.

        Returns
        -------
        Marker
            The newly created marker.
        """
        marker = Marker(body=self, local_position=local_position,
                        theta=theta, name=name)
        for existing in self._markers:
            if np.linalg.norm(existing.local_position - marker.local_position) < 1e-10:
                warnings.warn(
                    f"Body already has a marker at position "
                    f"{marker.local_position}")
                break
        self._markers.append(marker)
        return marker

    def add_marker_at(self, reference_marker, offset=(0, 0), theta=None,
                      name=None):
        """Create a marker whose position is deferred until assembly.

        The marker's local position will be computed from a reference marker
        on another body plus a local offset.

        Parameters
        ----------
        reference_marker : Marker
            A Marker on another body to use as reference.
        offset : array_like
            2-element offset in the reference body's local frame.
        theta : float, optional
            Orientation angle (radians).
        name : str, optional
            Human-readable name.

        Returns
        -------
        Marker
            The newly created (deferred) marker.
        """
        marker = Marker(body=self, local_position=[0, 0], theta=theta,
                        name=name)
        marker._deferred_ref = reference_marker
        marker._deferred_offset = offset
        self._markers.append(marker)
        return marker

    def __repr__(self):
        label = f"'{self.name}'" if self.name else f"#{self.COUNT}"
        return f"Body({label}, mass={self.mass}, inertia={self.inertia})"

    def get_position(self, component):
        """Get a single position component time history.

        Parameters
        ----------
        component : {'x', 'y', 'phi'}
            The position component to retrieve.

        Returns
        -------
        NDArray
            1-D array over time steps.

        Raises
        ------
        RuntimeError
            If no results are available (solve() not called).
        KeyError
            If component is invalid.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return self._result_container['positions'][component]

    def get_positions(self):
        """Get all position component time histories.

        Returns
        -------
        dict
            Dictionary with keys ``'x'``, ``'y'``, ``'phi'``, each
            mapping to a 1-D NDArray over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return dict(self._result_container['positions'])

    def get_velocity(self, component):
        """Get a single velocity component time history.

        Parameters
        ----------
        component : {'dx', 'dy', 'dphi'}
            The velocity component to retrieve.

        Returns
        -------
        NDArray
            1-D array over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        KeyError
            If component is invalid.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return self._result_container['velocities'][component]

    def get_velocities(self):
        """Get all velocity component time histories.

        Returns
        -------
        dict
            Dictionary with keys ``'dx'``, ``'dy'``, ``'dphi'``, each
            mapping to a 1-D NDArray over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return dict(self._result_container['velocities'])

    def get_acceleration(self, component):
        """Get a single acceleration component time history.

        Parameters
        ----------
        component : {'ddx', 'ddy', 'ddphi'}
            The acceleration component to retrieve.

        Returns
        -------
        NDArray
            1-D array over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        KeyError
            If component is invalid.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return self._result_container['accelerations'][component]

    def get_accelerations(self):
        """Get all acceleration component time histories.

        Returns
        -------
        dict
            Dictionary with keys ``'ddx'``, ``'ddy'``, ``'ddphi'``, each
            mapping to a 1-D NDArray over time steps.

        Raises
        ------
        RuntimeError
            If no results are available.
        """
        if self._result_container is None:
            raise RuntimeError("No results available. Run solve() first.")
        return dict(self._result_container['accelerations'])
