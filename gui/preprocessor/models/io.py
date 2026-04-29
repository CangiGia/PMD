"""JSON serialization for :class:`ModelSpec`.

Tuples are encoded as JSON lists; on load, fields known to be 2-tuples
are coerced back. The on-disk format is a versioned dict::

    {"format": "pmdmodel", "version": 1, "model": {...}}
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .spec import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    ShapeSpec,
)

FORMAT_TAG = "pmdmodel"
FORMAT_VERSION = 1


# ──────────────────────────────────────────────────────────────────
# Save
# ──────────────────────────────────────────────────────────────────

def to_dict(spec: ModelSpec) -> dict[str, Any]:
    """Convert a :class:`ModelSpec` into a JSON-friendly dict."""
    return {
        "format": FORMAT_TAG,
        "version": FORMAT_VERSION,
        "model": asdict(spec),
    }


def save_model(spec: ModelSpec, path: str | Path) -> None:
    """Serialize ``spec`` to a JSON file at ``path``."""
    p = Path(path)
    p.write_text(json.dumps(to_dict(spec), indent=2), encoding="utf-8")


# ──────────────────────────────────────────────────────────────────
# Load
# ──────────────────────────────────────────────────────────────────

def _xy(v: Any, default: tuple[float, float] = (0.0, 0.0)) -> tuple[float, float]:
    if isinstance(v, (list, tuple)) and len(v) == 2:
        return (float(v[0]), float(v[1]))
    return default


def _shape_from_dict(d: dict | None) -> ShapeSpec | None:
    if not d:
        return None
    return ShapeSpec(kind=d.get("kind", "rectangle"),
                     params=dict(d.get("params") or {}))


def _body_from_dict(d: dict) -> BodySpec:
    return BodySpec(
        id=d.get("id") or "",
        name=d.get("name", ""),
        mass=float(d.get("mass", 1.0)),
        inertia=float(d.get("inertia", 1.0)),
        position=_xy(d.get("position")),
        orientation=float(d.get("orientation", 0.0)),
        velocity=_xy(d.get("velocity")),
        angular_velocity=float(d.get("angular_velocity", 0.0)),
        shape=_shape_from_dict(d.get("shape")),
        locked=bool(d.get("locked", False)),
    )


def _marker_from_dict(d: dict) -> MarkerSpec:
    return MarkerSpec(
        id=d.get("id") or "",
        name=d.get("name", ""),
        body_id=d.get("body_id", "ground"),
        local_position=_xy(d.get("local_position")),
        theta=float(d.get("theta", 0.0)),
    )


def _joint_from_dict(d: dict) -> JointSpec:
    return JointSpec(
        id=d.get("id") or "",
        name=d.get("name", ""),
        kind=d.get("kind", "RevJoint"),
        i_marker_id=d.get("i_marker_id", ""),
        j_marker_id=d.get("j_marker_id", ""),
        params=dict(d.get("params") or {}),
    )


def _force_from_dict(d: dict) -> ForceSpec:
    return ForceSpec(
        id=d.get("id") or "",
        name=d.get("name", ""),
        kind=d.get("kind", "Weight"),
        i_marker_id=d.get("i_marker_id", ""),
        j_marker_id=d.get("j_marker_id", ""),
        params=dict(d.get("params") or {}),
    )


def from_dict(payload: dict) -> ModelSpec:
    """Rebuild a :class:`ModelSpec` from a payload produced by :func:`to_dict`.

    ``payload`` may either be the wrapped envelope or the inner ``model``
    dict directly (for forward-compat).
    """
    if "model" in payload and "format" in payload:
        if payload.get("format") != FORMAT_TAG:
            raise ValueError(f"Unsupported file format: {payload.get('format')!r}")
        m = payload["model"]
    else:
        m = payload

    return ModelSpec(
        name=m.get("name", "Untitled"),
        bodies=[_body_from_dict(b) for b in m.get("bodies", [])],
        markers=[_marker_from_dict(x) for x in m.get("markers", [])],
        joints=[_joint_from_dict(x) for x in m.get("joints", [])],
        forces=[_force_from_dict(x) for x in m.get("forces", [])],
        units=dict(m.get("units") or {}),
    )


def load_model(path: str | Path) -> ModelSpec:
    """Load a :class:`ModelSpec` from a JSON file."""
    p = Path(path)
    payload = json.loads(p.read_text(encoding="utf-8"))
    return from_dict(payload)
