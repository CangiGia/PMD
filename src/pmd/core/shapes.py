"""Shape descriptors for body visualization.

Author: Giacomo Cangi
"""

import math

import numpy as np
from dataclasses import dataclass

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
}


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
class Link:
    """Rod-shaped body (rectangle along the body's local x-axis).

    Typically used for connecting links, cranks, and rods. The shape is
    centered on the body's CoM with the long axis aligned with local x.

    Attributes
    ----------
    length : float
        Length along the body's local x-axis (m).
    thickness : float
        Cross-section height along the body's local y-axis (m).
    """
    length: float
    thickness: float


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


def _validate_shape(shape) -> None:
    """Raise ``ValueError`` if *shape* has non-physical parameters.

    Parameters
    ----------
    shape : Rectangle, Circle, Link, or Polygon
        Shape instance to check.

    Raises
    ------
    ValueError
        If any dimension is non-positive, or a Polygon is not CCW / has
        fewer than 3 vertices.
    TypeError
        If *shape* is not a recognised shape type.
    """
    if isinstance(shape, Rectangle):
        if shape.width <= 0 or shape.height <= 0:
            raise ValueError(
                f"Rectangle width and height must be positive, "
                f"got width={shape.width}, height={shape.height}")
    elif isinstance(shape, Circle):
        if shape.radius <= 0:
            raise ValueError(
                f"Circle radius must be positive, got {shape.radius}")
    elif isinstance(shape, Link):
        if shape.length <= 0 or shape.thickness <= 0:
            raise ValueError(
                f"Link length and thickness must be positive, "
                f"got length={shape.length}, thickness={shape.thickness}")
    elif isinstance(shape, Polygon):
        verts = shape.vertices
        if len(verts) < 3:
            raise ValueError(
                f"Polygon must have at least 3 vertices, got {len(verts)}")
        if np.any(np.isnan(verts)):
            raise ValueError("Polygon vertices contain NaN values")
        n = len(verts)
        signed_area = sum(
            verts[i, 0] * verts[(i + 1) % n, 1]
            - verts[(i + 1) % n, 0] * verts[i, 1]
            for i in range(n)
        )
        if signed_area <= 0:
            raise ValueError(
                "Polygon vertices must be in counter-clockwise order "
                "(signed area must be positive)")
    else:
        raise TypeError(f"Unsupported shape type: {type(shape).__name__}")


def compute_mass_props(shape, *, density: float,
                       thickness_z: float = 0.01) -> tuple[float, float]:
    """Return ``(mass, inertia)`` for a planar body with uniform density.

    The body is modelled as a planar slab of out-of-plane thickness
    *thickness_z* and uniform material *density*.

    Parameters
    ----------
    shape : Rectangle, Circle, Link, or Polygon
        Geometry of the body.
    density : float
        Material density (kg/m³).
    thickness_z : float
        Out-of-plane thickness (m). Defaults to 0.01 m (10 mm).

    Returns
    -------
    tuple[float, float]
        ``(mass, moment_of_inertia_about_z)`` in SI units.
    """
    rho = float(density)
    tz  = float(thickness_z)

    if isinstance(shape, Rectangle):
        w, h = shape.width, shape.height
        area = w * h
        mass = rho * area * tz
        inertia = mass * (w * w + h * h) / 12.0
        return mass, inertia

    if isinstance(shape, Link):
        L, t = shape.length, shape.thickness
        area = L * t
        mass = rho * area * tz
        inertia = mass * L * L / 12.0
        return mass, inertia

    if isinstance(shape, Circle):
        r = shape.radius
        area = math.pi * r * r
        mass = rho * area * tz
        inertia = 0.5 * mass * r * r
        return mass, inertia

    if isinstance(shape, Polygon):
        verts = shape.vertices
        n = len(verts)
        signed_area_2 = 0.0
        Jz = 0.0
        for i in range(n):
            xi, yi = verts[i]
            xj, yj = verts[(i + 1) % n]
            cross = xi * yj - xj * yi
            signed_area_2 += cross
            Jz += (xi*xi + xi*xj + xj*xj + yi*yi + yi*yj + yj*yj) * cross
        area = abs(signed_area_2) / 2.0
        mass = rho * area * tz
        inertia = mass * abs(Jz) / (12.0 * area) if area > 0 else 0.0
        return mass, inertia

    raise TypeError(f"Unsupported shape type: {type(shape).__name__}")
