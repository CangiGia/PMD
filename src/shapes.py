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

    Attributes
    ----------
    width : float
        Extent along the body's local x-axis (m).
    height : float
        Extent along the body's local y-axis (m).
    """
    width: float
    height: float


@dataclass(frozen=True)
class Circle:
    """Circular body shape centered on the body's CoM.

    Attributes
    ----------
    radius : float
        Circle radius (m).
    """
    radius: float


@dataclass(frozen=True)
class Polygon:
    """Arbitrary polygon shape in the body's local frame.

    Attributes
    ----------
    vertices : numpy.ndarray
        Array of shape (N, 2) with vertices in local frame, CCW order.
    """
    vertices: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, 'vertices', np.asarray(self.vertices, dtype=float))
