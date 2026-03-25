"""Shape descriptors for body visualization.

These are pure data objects — they carry no rendering logic.
The GUI layer decides how to draw them.

Author: Giacomo Cangi
"""

import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True)
class Rectangle:
    """Rectangular body shape centered on the body's CoM.

    Args:
        width: Extent along the body's local x-axis (m).
        height: Extent along the body's local y-axis (m).
    """
    width: float
    height: float


@dataclass(frozen=True)
class Circle:
    """Circular body shape centered on the body's CoM.

    Args:
        radius: Circle radius (m).
    """
    radius: float


@dataclass(frozen=True)
class Polygon:
    """Arbitrary polygon shape in the body's local frame.

    Args:
        vertices: (N, 2) array of vertices in local frame, CCW order.
    """
    vertices: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, 'vertices', np.asarray(self.vertices, dtype=float))
