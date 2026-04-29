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
    MATERIALS,
    GROUND_BODY_ID,
    compute_mass_props,
)
from .io import load_model, save_model, to_dict, from_dict

__all__ = [
    "BodySpec",
    "ForceSpec",
    "JointSpec",
    "MarkerSpec",
    "ModelSpec",
    "ShapeSpec",
    "MATERIALS",
    "GROUND_BODY_ID",
    "compute_mass_props",
    "load_model",
    "save_model",
    "to_dict",
    "from_dict",
]
