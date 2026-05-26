"""Core model components: Ground, Base, Body, Marker.

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

    Attributes
    ----------
    COUNT : int
        Class variable tracking the number of instances created.
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

    Attributes
    ----------
    body : Body or _GroundType
        The body this marker is attached to.
    local_position : ndarray
        Flat (2,) array of local coordinates [xi, eta].
    theta : float or None
        Optional orientation angle in the local frame.
    name : str or None
        Optional human-readable name.
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

    Mass properties can be specified in two ways:

    1. **Explicit** – provide ``mass`` and ``inertia`` directly.  ``shape``
       is optional and used only for visualisation.
    2. **Auto-derived** – provide ``shape`` and ``density``; mass and moment
       of inertia are computed from the geometry and ``thickness_z``.

    When ``density`` is supplied together with explicit ``mass``/``inertia``,
    the auto-derived values take precedence unless ``mass_override=True``.

    Attributes
    ----------
    name : str or None
        Optional human-readable name.
    mass : float
        Mass of the body (kg).
    inertia : float
        Moment of inertia of the body (kg·m²).
    mass_override : bool
        ``True`` when mass and inertia were supplied explicitly by the user
        (overriding auto-derivation from shape + density).
    position : ndarray
        Position vector (x, y), shape (2, 1).
    orientation : float
        Orientation angle phi (rad).
    velocity : ndarray
        Velocity vector (dx, dy), shape (2, 1).
    angular_velocity : float
        Angular velocity dphi (rad/s).
    acceleration : ndarray
        Acceleration vector (ddx, ddy), shape (2, 1).
    angular_acceleration : float
        Angular acceleration ddphi (rad/s²).
    shape : Rectangle or Circle or Plate or Link or None
        Optional shape descriptor for visualisation.
    """

    def __init__(self, *, shape=None, density=None, thickness_z=0.01,
                 mass=None, inertia=None, mass_override=False,
                 position=None, orientation=None,
                 velocity=None, angular_velocity=0,
                 acceleration=None, angular_acceleration=0,
                 name=None, color=None):
        """Initialize a Body.

        Parameters
        ----------
        shape : Rectangle, Circle, Plate, or Link, optional
            Geometry used for visualisation and, when ``density`` is given,
            for automatic mass-property computation.
        density : float, optional
            Material density (kg/m³).  Required when no explicit
            ``mass``/``inertia`` are given.
        thickness_z : float, optional
            Out-of-plane thickness (m) used in auto-derivation.
            Defaults to 0.01 m (10 mm).
        mass : float, optional
            Mass (kg, must be positive).  Required when ``density`` is not
            given, or when ``mass_override=True``.
        inertia : float, optional
            Moment of inertia (kg·m², must be non-negative).  Same rules as
            ``mass``.
        mass_override : bool, optional
            When ``True`` and ``density`` is also provided, the explicit
            ``mass``/``inertia`` values are kept and density is ignored for
            mass computation.  Defaults to ``False``.
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

        Raises
        ------
        ValueError
            If mass is not positive, inertia is negative, or if the
            combination of arguments is inconsistent.

        Examples
        --------
        Classic (explicit mass/inertia, shape for visualisation only)::

            b = Body(mass=1.5, inertia=0.02)
            b = Body(mass=1.5, inertia=0.02, shape=Rectangle(0.4, 0.1))

        Auto-derived from geometry::

            b = Body(shape=Rectangle(0.4, 0.1), density=7850)
            b = Body(shape=Plate(...), density=2700, thickness_z=0.005)

        Explicit override when density is present (GUI workflow)::

            b = Body(shape=Rectangle(0.4, 0.1), density=7850,
                     mass=2.0, inertia=0.05, mass_override=True)
        """
        super().__init__()

        # ── Resolve mass and inertia ───────────────────────────────
        user_gave_mass = (mass is not None) or (inertia is not None)

        if shape is None:
            # Classic mode: explicit mass and inertia required.
            if mass is None or inertia is None:
                raise ValueError(
                    f"Body {self.COUNT}: 'mass' and 'inertia' are required "
                    f"when 'shape' is not provided.")
        elif density is not None:
            # Auto-derive from shape + density (unless mass_override).
            if mass_override:
                if mass is None or inertia is None:
                    raise ValueError(
                        f"Body {self.COUNT}: mass_override=True requires "
                        f"explicit 'mass' and 'inertia'.")
            else:
                from .shapes import compute_mass_props
                auto_m, auto_J = compute_mass_props(
                    shape, density=density, thickness_z=thickness_z)
                if user_gave_mass:
                    warnings.warn(
                        f"Body {self.COUNT}: density given; auto-computed "
                        f"mass/inertia overwrite the user-supplied values. "
                        f"Pass mass_override=True to keep user values.",
                        UserWarning, stacklevel=2)
                mass, inertia = auto_m, auto_J
        else:
            # Shape present, no density → explicit mass and inertia required.
            if mass is None or inertia is None:
                raise ValueError(
                    f"Body {self.COUNT}: 'shape' without 'density' requires "
                    f"explicit 'mass' and 'inertia' (shape is used for "
                    f"visualisation only).  Provide 'density' to auto-compute.")

        # ── Validate resolved values ───────────────────────────────
        if mass <= 0:
            raise ValueError(
                f"Body {self.COUNT}: mass must be positive, got {mass}")
        if inertia < 0:
            raise ValueError(
                f"Body {self.COUNT}: moment of inertia cannot be negative, "
                f"got {inertia}")
        if shape is not None:
            from .shapes import _validate_shape
            _validate_shape(shape)

        self.name = name
        self.color = color
        self.shape = shape
        self.mass = mass
        self.inertia = inertia
        self.mass_override = mass_override
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

        # Auto-create a marker at the body's centre of mass.  Lives in
        # local coords [0, 0] regardless of how the body was built (direct
        # constructor or geometry factory) and is named ``<body>.cm`` so
        # it can be referenced unambiguously from joints, forces, etc.
        _cm_label = self.name if self.name else f"body_{self.COUNT}"
        self.cm = self.add_marker([0.0, 0.0], name=f"{_cm_label}.cm")

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

    # ── Geometry-driven factories ─────────────────────────────────
    @classmethod
    def as_link(cls, p1, p2, *, thickness, marker_theta=None, **body_kwargs):
        """Build a Link-shaped body from its two world-frame endpoints.

        Parameters
        ----------
        p1, p2 : array_like
            World coordinates ``(x, y)`` of the link endpoints.
        thickness : float
            Link visual/geometric thickness (m).
        marker_theta : float or None, optional
            Local orientation (rad) applied to both endpoint markers.
            Defaults to ``None`` (markers without orientation kinematics).
            Set to ``0.0`` to give the endpoint markers an x-axis aligned
            with the link itself (required by translational joints,
            actuators with ``control='length'``, etc.).
        **body_kwargs
            Forwarded to :class:`Body` (``mass``, ``inertia``, ``density``,
            ``name``, ...).  ``position``, ``orientation`` and ``shape``
            are computed automatically and must not be supplied.

        Returns
        -------
        body : Body
        mk_p1, mk_p2 : Marker
            Markers attached at ``p1`` and ``p2`` (local coords
            ``[-L/2, 0]`` and ``[+L/2, 0]``).
        """
        for k in ("position", "orientation", "shape"):
            if k in body_kwargs:
                raise TypeError(f"Body.from_link() does not accept {k!r}")
        from .shapes import Link
        p1 = np.asarray(p1, dtype=float).flatten()
        p2 = np.asarray(p2, dtype=float).flatten()
        d = p2 - p1
        length = float(np.linalg.norm(d))
        if length < 1e-12:
            raise ValueError("Link endpoints coincide.")
        body = cls(
            position=(p1 + p2) / 2.0,
            orientation=float(np.arctan2(d[1], d[0])),
            shape=Link(length=length, thickness=thickness),
            **body_kwargs,
        )
        mk_p1 = body.add_marker([-length / 2.0, 0.0], theta=marker_theta)
        mk_p2 = body.add_marker([ length / 2.0, 0.0], theta=marker_theta)
        return body, mk_p1, mk_p2

    @classmethod
    def as_plate(cls, v1, v2, v3, *, orientation=0.0, **body_kwargs):
        """Build a triangular Plate body from its three world-frame vertices.

        Parameters
        ----------
        v1, v2, v3 : array_like
            World coordinates of the triangle vertices in
            **counter-clockwise** order.
        orientation : float, optional
            Body orientation angle (rad).  Defaults to ``0`` so the local
            frame is axis-aligned with the world frame.
        **body_kwargs
            Forwarded to :class:`Body`.  ``position`` and ``shape`` are
            computed automatically and must not be supplied.

        Returns
        -------
        body : Body
        mk_v1, mk_v2, mk_v3 : Marker
            Markers attached at the three vertices.
        """
        for k in ("position", "shape"):
            if k in body_kwargs:
                raise TypeError(f"Body.from_plate() does not accept {k!r}")
        from .shapes import Plate
        shape, centroid = Plate.from_world_points(
            tuple(np.asarray(v1, dtype=float).flatten()),
            tuple(np.asarray(v2, dtype=float).flatten()),
            tuple(np.asarray(v3, dtype=float).flatten()),
        )
        body = cls(
            position=centroid,
            orientation=orientation,
            shape=shape,
            **body_kwargs,
        )
        mk_v1 = body.add_marker((np.asarray(v1, dtype=float).flatten() - centroid).tolist())
        mk_v2 = body.add_marker((np.asarray(v2, dtype=float).flatten() - centroid).tolist())
        mk_v3 = body.add_marker((np.asarray(v3, dtype=float).flatten() - centroid).tolist())
        return body, mk_v1, mk_v2, mk_v3

    @classmethod
    def as_rectangle(cls, c1, c2, c3, c4, **body_kwargs):
        """Build a Rectangle body from its four world-frame corners (CCW).

        The local frame is oriented so that the edge ``c1→c2`` lies along
        the local +x axis (``width``) and ``c2→c3`` along the local +y
        axis (``height``).

        Parameters
        ----------
        c1, c2, c3, c4 : array_like
            World coordinates of the four corners in
            **counter-clockwise** order.
        **body_kwargs
            Forwarded to :class:`Body`.  ``position``, ``orientation`` and
            ``shape`` are computed automatically and must not be supplied.

        Returns
        -------
        body : Body
        mk_c1, mk_c2, mk_c3, mk_c4 : Marker
            Markers attached at the four corners.
        """
        for k in ("position", "orientation", "shape"):
            if k in body_kwargs:
                raise TypeError(f"Body.from_rectangle() does not accept {k!r}")
        from .shapes import Rectangle
        corners = [np.asarray(c, dtype=float).flatten() for c in (c1, c2, c3, c4)]
        d12 = corners[1] - corners[0]
        d23 = corners[2] - corners[1]
        d34 = corners[3] - corners[2]
        d41 = corners[0] - corners[3]
        width  = float(np.linalg.norm(d12))
        height = float(np.linalg.norm(d23))
        if width < 1e-12 or height < 1e-12:
            raise ValueError("Rectangle has a degenerate (zero-length) side.")
        atol = 1e-9 * max(width, height)
        if (not np.isclose(np.linalg.norm(d34), width,  atol=atol) or
            not np.isclose(np.linalg.norm(d41), height, atol=atol) or
            not np.isclose(float(np.dot(d12, d23)), 0.0, atol=atol * max(width, height))):
            raise ValueError("Four corners do not form a rectangle.")
        centroid = sum(corners) / 4.0
        body = cls(
            position=centroid,
            orientation=float(np.arctan2(d12[1], d12[0])),
            shape=Rectangle(width=width, height=height),
            **body_kwargs,
        )
        mk_c1 = body.add_marker([-width / 2.0, -height / 2.0])
        mk_c2 = body.add_marker([ width / 2.0, -height / 2.0])
        mk_c3 = body.add_marker([ width / 2.0,  height / 2.0])
        mk_c4 = body.add_marker([-width / 2.0,  height / 2.0])
        return body, mk_c1, mk_c2, mk_c3, mk_c4

    @classmethod
    def as_circle(cls, center, radius, *, orientation=0.0, **body_kwargs):
        """Build a Circle body from its world-frame centre and radius.

        Parameters
        ----------
        center : array_like
            World coordinates ``(x, y)`` of the disc centre.
        radius : float
            Disc radius (m).
        orientation : float, optional
            Body orientation angle (rad).  A circle is rotationally
            symmetric, but the local frame is still meaningful for
            attaching markers later via :meth:`add_marker`.
        **body_kwargs
            Forwarded to :class:`Body`.  ``position`` and ``shape`` are
            computed automatically and must not be supplied.

        Returns
        -------
        body : Body
            No natural markers are returned (the disc is symmetric).
            Attach markers explicitly via :meth:`add_marker` if needed.
        """
        for k in ("position", "shape"):
            if k in body_kwargs:
                raise TypeError(f"Body.from_circle() does not accept {k!r}")
        from .shapes import Circle
        body = cls(
            position=np.asarray(center, dtype=float).flatten(),
            orientation=orientation,
            shape=Circle(radius=radius),
            **body_kwargs,
        )
        return body

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
