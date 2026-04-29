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
#  Materials database (used by the GUI to auto-derive mass props)
# ────────────────────────────────────────────────────────────────────

#: Common engineering materials → density in kg/m³.
MATERIALS: dict[str, float] = {
    "Steel":         7850.0,
    "Stainless":     7900.0,
    "Cast Iron":     7200.0,
    "Aluminum":      2700.0,
    "Titanium":      4500.0,
    "Brass":         8500.0,
    "Copper":        8960.0,
    "Plastic (PLA)": 1240.0,
    "Plastic (ABS)": 1050.0,
    "Wood (Oak)":     750.0,
    "Custom":           0.0,
}


def compute_mass_props(spec: "BodySpec") -> tuple[float, float] | None:
    """Return ``(mass, inertia)`` derived from a body's shape + density.

    The model is treated as a planar slab of out-of-plane thickness
    ``spec.thickness_z`` and uniform density ``spec.density``. Returns
    ``None`` if the body has no shape (so the caller can keep the
    user-typed values).
    """
    if spec.shape is None:
        return None
    p = spec.shape.params
    kind = spec.shape.kind
    rho = float(spec.density)
    tz  = float(spec.thickness_z)
    if kind == "rectangle":
        w = float(p.get("width", 0.0)); h = float(p.get("height", 0.0))
        area = w * h
        mass = rho * area * tz
        inertia = mass * (w * w + h * h) / 12.0
        return mass, inertia
    if kind == "link":
        L = float(p.get("length", 0.0)); t = float(p.get("thickness", 0.0))
        area = L * t
        mass = rho * area * tz
        inertia = mass * L * L / 12.0
        return mass, inertia
    if kind == "circle":
        r = float(p.get("radius", 0.0))
        import math as _m
        area = _m.pi * r * r
        mass = rho * area * tz
        inertia = 0.5 * mass * r * r
        return mass, inertia
    return None


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
    # Material / mass-properties strategy (GUI-only — the runtime model
    # only ever reads ``mass`` and ``inertia``). When ``mass_override``
    # is False the GUI re-derives mass + inertia from ``density`` and
    # the shape footprint times ``thickness_z`` (out-of-plane depth);
    # otherwise the user-typed mass/inertia values are kept verbatim.
    material: str = "Steel"
    density: float = 7850.0       # kg / m³
    thickness_z: float = 0.01     # m  (out-of-plane depth)
    mass_override: bool = False


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


# Reserved id of the implicit ground body.
GROUND_BODY_ID = "ground"


def _make_ground_body() -> "BodySpec":
    """Return the canonical implicit ground body.

    The ground body is read-only by convention: ``locked`` is ``True``
    so the GUI editors cannot edit its position / shape, and it carries
    a reserved id (``"ground"``) the runtime builder filters out before
    constructing the planar model.
    """
    return BodySpec(
        id=GROUND_BODY_ID,
        name="ground",
        mass=0.0,
        inertia=0.0,
        position=(0.0, 0.0),
        orientation=0.0,
        shape=None,
        locked=True,
        material="Custom",
        density=0.0,
        thickness_z=0.0,
        mass_override=True,
    )


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

    def __post_init__(self) -> None:
        # The ground body is *always* present in a model. It is created
        # implicitly so the user can attach markers to it without having
        # to define it explicitly. The runtime builder skips it (the
        # solver-side ``Ground`` singleton represents the same concept).
        if not any(b.id == GROUND_BODY_ID for b in self.bodies):
            self.bodies.insert(0, _make_ground_body())

    # ── Lookup helpers ─────────────────────────────────────────
    def body(self, body_id: str) -> BodySpec | None:
        return next((b for b in self.bodies if b.id == body_id), None)

    def marker(self, marker_id: str) -> MarkerSpec | None:
        return next((m for m in self.markers if m.id == marker_id), None)

    def joint(self, joint_id: str) -> JointSpec | None:
        return next((j for j in self.joints if j.id == joint_id), None)

    def force(self, force_id: str) -> ForceSpec | None:
        return next((f for f in self.forces if f.id == force_id), None)
