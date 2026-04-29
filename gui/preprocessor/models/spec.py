"""Dataclass specifications for the pre-processor model graph.

The GUI manipulates these light, serializable objects. A separate
*builder* module (to be implemented) will translate a ``ModelSpec``
into a runtime ``PMD.src.model.PlanarMultibodyModel`` at solve time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal
from uuid import uuid4


def _new_id(prefix: str) -> str:
    return f"{prefix}_{uuid4().hex[:8]}"


# ────────────────────────────────────────────────────────────────────
#  Geometry
# ────────────────────────────────────────────────────────────────────

ShapeKind = Literal["rectangle", "circle", "polygon", "link"]


@dataclass
class ShapeSpec:
    """Visual / inertial shape attached to a body.

    Attributes
    ----------
    kind : str
        ``"rectangle"``, ``"circle"``, ``"polygon"``, or ``"link"``
        (a simple rod between two endpoints).
    params : dict
        Shape-specific parameters. Examples:

        * rectangle: ``{"width": w, "height": h}``
        * circle:    ``{"radius": r}``
        * polygon:   ``{"vertices": [(xi, eta), ...]}`` (local coords)
        * link:      ``{"length": L, "thickness": t}``
    """

    kind: ShapeKind = "rectangle"
    params: dict = field(default_factory=dict)


# ────────────────────────────────────────────────────────────────────
#  Markers
# ────────────────────────────────────────────────────────────────────


@dataclass
class MarkerSpec:
    """Body-fixed reference frame.

    Attributes
    ----------
    id : str
        Unique identifier.
    name : str
        Human-readable label.
    body_id : str
        Owning body id (use ``"ground"`` for the ground frame).
    local_position : tuple of float
        ``(xi, eta)`` in the body's local frame.
    theta : float
        Orientation angle (rad) in the body's local frame.
    """

    id: str = field(default_factory=lambda: _new_id("mrk"))
    name: str = ""
    body_id: str = "ground"
    local_position: tuple[float, float] = (0.0, 0.0)
    theta: float = 0.0


# ────────────────────────────────────────────────────────────────────
#  Bodies
# ────────────────────────────────────────────────────────────────────


@dataclass
class BodySpec:
    """Rigid body specification.

    Attributes
    ----------
    id : str
        Unique identifier.
    name : str
        Human-readable label.
    mass : float
        Mass (kg).
    inertia : float
        Moment of inertia about the body centroid (kg·m²).
    position : tuple of float
        Initial CoM position ``(x, y)`` in the global frame.
    orientation : float
        Initial orientation angle phi (rad).
    velocity : tuple of float
        Initial linear velocity ``(dx, dy)``.
    angular_velocity : float
        Initial angular velocity (rad/s).
    shape : ShapeSpec or None
        Visual / inertial shape.
    """

    id: str = field(default_factory=lambda: _new_id("body"))
    name: str = ""
    mass: float = 1.0
    inertia: float = 1.0
    position: tuple[float, float] = (0.0, 0.0)
    orientation: float = 0.0
    velocity: tuple[float, float] = (0.0, 0.0)
    angular_velocity: float = 0.0
    shape: ShapeSpec | None = None
    # If True, geometry-defining fields (position, orientation, shape
    # dimensions) are considered frozen by the GUI: editors render them
    # read-only with a red highlight so the user can't accidentally
    # invalidate the markers tied to this body.
    locked: bool = False


# ────────────────────────────────────────────────────────────────────
#  Joints
# ────────────────────────────────────────────────────────────────────

JointKind = Literal[
    "RevJoint",
    "TranJoint",
    "RevRevJoint",
    "RevTranJoint",
    "RigidJoint",
    "DiscJoint",
    "RelRotJoint",
    "RelTranJoint",
]


@dataclass
class JointSpec:
    """Joint specification (Adams-style: marker-based).

    Attributes
    ----------
    id : str
        Unique identifier.
    name : str
        Human-readable label.
    kind : str
        Joint type (see :data:`JointKind`).
    i_marker_id : str
        Marker on body i.
    j_marker_id : str
        Marker on body j.
    params : dict
        Joint-specific parameters (e.g. ``{"fix": 1}`` for RevJoint,
        ``{"L": 0.3}`` for RevRevJoint, ``{"R": 0.05}`` for DiscJoint).
    """

    id: str = field(default_factory=lambda: _new_id("jnt"))
    name: str = ""
    kind: JointKind = "RevJoint"
    i_marker_id: str = ""
    j_marker_id: str = ""
    params: dict = field(default_factory=dict)


# ────────────────────────────────────────────────────────────────────
#  Forces
# ────────────────────────────────────────────────────────────────────

ForceKind = Literal[
    "Weight",
    "PtpForce",
    "RotSdaForce",
    "LocalForce",
    "GlobalForce",
    "Torque",
    "UserForce",
]


@dataclass
class ForceSpec:
    """Force / torque element specification.

    Attributes
    ----------
    id : str
        Unique identifier.
    name : str
        Human-readable label.
    kind : str
        Force type (see :data:`ForceKind`).
    i_marker_id : str
        Marker on body i (use ``""`` if not applicable, e.g. Weight).
    j_marker_id : str
        Marker on body j (use ``""`` if not applicable).
    params : dict
        Force-specific parameters (e.g. ``{"k": 1e3, "L0": 0.2,
        "dc": 10.0, "f_a": 0.0}`` for PtpForce).
    """

    id: str = field(default_factory=lambda: _new_id("frc"))
    name: str = ""
    kind: ForceKind = "Weight"
    i_marker_id: str = ""
    j_marker_id: str = ""
    params: dict = field(default_factory=dict)


# ────────────────────────────────────────────────────────────────────
#  Top-level container
# ────────────────────────────────────────────────────────────────────


@dataclass
class ModelSpec:
    """Whole-project specification edited by the pre-processor.

    Attributes
    ----------
    name : str
        Project name.
    bodies : list of BodySpec
        Moving rigid bodies (ground is implicit).
    markers : list of MarkerSpec
        All markers, on both ground and moving bodies.
    joints : list of JointSpec
        Constraint elements.
    forces : list of ForceSpec
        Force / torque elements.
    units : dict
        Unit-system overrides for display (e.g.
        ``{"length": "mm", "force": "N"}``).
    """

    name: str = "Untitled"
    bodies: list[BodySpec] = field(default_factory=list)
    markers: list[MarkerSpec] = field(default_factory=list)
    joints: list[JointSpec] = field(default_factory=list)
    forces: list[ForceSpec] = field(default_factory=list)
    units: dict = field(default_factory=dict)

    # ── Lookup helpers ─────────────────────────────────────────
    def body(self, body_id: str) -> BodySpec | None:
        return next((b for b in self.bodies if b.id == body_id), None)

    def marker(self, marker_id: str) -> MarkerSpec | None:
        return next((m for m in self.markers if m.id == marker_id), None)

    def joint(self, joint_id: str) -> JointSpec | None:
        return next((j for j in self.joints if j.id == joint_id), None)

    def force(self, force_id: str) -> ForceSpec | None:
        return next((f for f in self.forces if f.id == force_id), None)
