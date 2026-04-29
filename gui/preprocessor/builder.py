"""Compile a :class:`ModelSpec` into a runnable :class:`PlanarMultibodyModel`.

The GUI manipulates lightweight dataclass specs; this module is the
single place where they are translated into the runtime objects from
:mod:`PMD.src`. Keep the mapping in sync with the spec definitions.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple

from PMD.src import (
    Body,
    Ground,
    PlanarMultibodyModel,
    Weight,
    PtpForce,
    RevJoint,
    TranJoint,
    RevRevJoint,
    RevTranJoint,
    RigidJoint,
    DiscJoint,
    RelRotJoint,
    RelTranJoint,
    LocalForce,
    GlobalForce,
    Torque,
    RotSdaForce,
)
from .models import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
)


_JOINT_CLASSES = {
    "RevJoint":     RevJoint,
    "TranJoint":    TranJoint,
    "RevRevJoint":  RevRevJoint,
    "RevTranJoint": RevTranJoint,
    "RigidJoint":   RigidJoint,
    "DiscJoint":    DiscJoint,
    "RelRotJoint":  RelRotJoint,
    "RelTranJoint": RelTranJoint,
}

_FORCE_CLASSES = {
    "Weight":      Weight,
    "PtpForce":    PtpForce,
    "RotSdaForce": RotSdaForce,
    "LocalForce":  LocalForce,
    "GlobalForce": GlobalForce,
    "Torque":      Torque,
    # "UserForce" requires a callable → not buildable from the GUI alone.
}


class BuilderError(RuntimeError):
    """Raised when a spec cannot be compiled into a runtime model."""


# ──────────────────────────────────────────────────────────────────
def build_model(spec: ModelSpec) -> PlanarMultibodyModel:
    """Translate ``spec`` into a fully wired :class:`PlanarMultibodyModel`.

    Parameters
    ----------
    spec : ModelSpec
        Source-of-truth specification edited in the pre-processor.

    Returns
    -------
    PlanarMultibodyModel
        Ready to call ``.solve(...)``.
    """
    # 1) Bodies (and remember the runtime instance per spec id)
    body_by_id: Dict[str, Body] = {}
    for b in spec.bodies:
        body_by_id[b.id] = _make_body(b)

    # 2) Markers (Ground or owning body)
    marker_by_id: Dict[str, Any] = {}
    for m in spec.markers:
        marker_by_id[m.id] = _make_marker(m, body_by_id)

    # 3) Joints
    joints = [
        _make_joint(j, marker_by_id) for j in spec.joints
    ]

    # 4) Forces
    forces = [
        _make_force(f, marker_by_id) for f in spec.forces
    ]

    return PlanarMultibodyModel(
        bodies=list(body_by_id.values()),
        joints=joints,
        forces=forces,
    )


# ──────────────────────────────────────────────────────────────────
def _make_body(b: BodySpec) -> Body:
    return Body(
        mass=b.mass,
        inertia=b.inertia,
        position=list(b.position),
        orientation=b.orientation,
        name=b.name or b.id,
    )


def _make_marker(m: MarkerSpec, body_by_id: Dict[str, Body]):
    if m.body_id == "ground":
        owner = Ground
    else:
        owner = body_by_id.get(m.body_id)
        if owner is None:
            raise BuilderError(
                f"Marker {m.id!r} references unknown body {m.body_id!r}")
    # Marker.add_marker accepts an optional theta
    return owner.add_marker(list(m.local_position), theta=m.theta)


def _make_joint(j: JointSpec, marker_by_id: Dict[str, Any]):
    cls = _JOINT_CLASSES.get(j.kind)
    if cls is None:
        raise BuilderError(f"Unsupported joint kind: {j.kind!r}")
    mi = _resolve_marker(j.i_marker_id, marker_by_id, "i", j)
    mj = _resolve_marker(j.j_marker_id, marker_by_id, "j", j)
    kwargs: Dict[str, Any] = {"iMarker": mi, "jMarker": mj,
                              "name": j.name or j.id}
    # Pass through any user-supplied numerical params (k, R, L, fix, q0…)
    for k, v in j.params.items():
        kwargs.setdefault(k, v)
    return cls(**kwargs)


def _make_force(f: ForceSpec, marker_by_id: Dict[str, Any]):
    cls = _FORCE_CLASSES.get(f.kind)
    if cls is None:
        raise BuilderError(f"Unsupported force kind: {f.kind!r}")

    if f.kind == "Weight":
        # Weight is global; markers are ignored.
        return cls(name=f.name or f.id, **f.params)

    mi = _resolve_marker(f.i_marker_id, marker_by_id, "i", f) \
        if f.i_marker_id else None
    mj = _resolve_marker(f.j_marker_id, marker_by_id, "j", f) \
        if f.j_marker_id else None
    kwargs: Dict[str, Any] = {"name": f.name or f.id, **f.params}
    if mi is not None:
        kwargs.setdefault("iMarker", mi)
    if mj is not None:
        kwargs.setdefault("jMarker", mj)
    return cls(**kwargs)


def _resolve_marker(mid: str, marker_by_id: Dict[str, Any],
                    side: str, owner) -> Any:
    if not mid:
        raise BuilderError(
            f"{type(owner).__name__} {getattr(owner, 'id', '?')!r}"
            f" has no {side}-marker assigned")
    m = marker_by_id.get(mid)
    if m is None:
        raise BuilderError(
            f"{type(owner).__name__} {getattr(owner, 'id', '?')!r}"
            f" references unknown marker {mid!r}")
    return m
