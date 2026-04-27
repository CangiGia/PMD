"""Serializable model specifications for the PMD pre-processor.

These dataclasses are the *source of truth* for the GUI: editing,
serialization, undo/redo all operate on them. They are compiled to a
runtime ``PlanarMultibodyModel`` only at solve time via :func:`build`.
"""

from .spec import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    ShapeSpec,
)

__all__ = [
    "BodySpec",
    "ForceSpec",
    "JointSpec",
    "MarkerSpec",
    "ModelSpec",
    "ShapeSpec",
]
