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
class Plate:
    """Triangular plate shape defined by three vertices in the body local frame.

    The frame origin is always the centroid (CoM) of the triangle, which is
    required by the PMD solver (mass matrix is diagonal only when the body
    origin coincides with the CoM).  Vertices passed by the user are
    automatically re-centred: the stored :attr:`vertices` are shifted so that
    their mean is exactly ``(0, 0)``.

    Parameters
    ----------
    v1, v2, v3 : tuple of float
        Triangle vertices ``(x, y)`` in the body's local frame, listed in
        **counter-clockwise** order.

    Attributes
    ----------
    vertices : numpy.ndarray, shape (3, 2)
        Vertices re-centred on the centroid.  Used for rendering and inertia
        computation.
    input_vertices : numpy.ndarray, shape (3, 2)
        Original vertices exactly as supplied by the user.

    Examples
    --------
    >>> plate = Plate(v1=(0, 0), v2=(2, 0), v3=(1, 1.5))
    >>> plate.vertices.mean(axis=0)   # centroid at origin
    array([0., 0.])

    To create a Plate from three world-frame pin positions use the
    :meth:`from_world_points` factory, which also returns the centroid
    world coordinates so you can set ``body.position`` consistently:

    >>> plate, centroid = Plate.from_world_points((0,0), (2,0), (1,1.5))
    """

    v1: tuple[float, float]
    v2: tuple[float, float]
    v3: tuple[float, float]

    def __post_init__(self):
        raw = np.array([self.v1, self.v2, self.v3], dtype=float)
        if raw.shape != (3, 2):
            raise ValueError(
                "Each vertex must be a 2-element (x, y) sequence.")
        centroid = raw.mean(axis=0)
        object.__setattr__(self, 'input_vertices', raw.copy())
        object.__setattr__(self, 'vertices', raw - centroid)

    @classmethod
    def from_world_points(
        cls,
        p1: tuple[float, float],
        p2: tuple[float, float],
        p3: tuple[float, float],
    ) -> tuple["Plate", np.ndarray]:
        """Create a :class:`Plate` from three **world-frame** pin positions.

        Computes the centroid, expresses the vertices in the local frame
        (centred on the centroid), and returns both the :class:`Plate`
        instance and the centroid world coordinates so that you can set
        ``body.position`` correctly.

        Parameters
        ----------
        p1, p2, p3 : tuple of float
            World-frame coordinates ``(x, y)`` of the three vertices, CCW.

        Returns
        -------
        plate : Plate
        centroid : numpy.ndarray, shape (2,)
            World position of the centroid (set this as ``body.position``).

        Examples
        --------
        >>> plate, centroid = Plate.from_world_points((0, 0), (2, 0), (1, 1.5))
        >>> body = Body(shape=plate, density=7850.0)
        >>> body.position = centroid.reshape(2, 1)
        """
        pts = np.array([p1, p2, p3], dtype=float)
        centroid = pts.mean(axis=0)
        local = pts - centroid
        return cls(
            v1=tuple(local[0]),
            v2=tuple(local[1]),
            v3=tuple(local[2]),
        ), centroid


def _validate_shape(shape) -> None:
    """Raise ``ValueError`` if *shape* has non-physical parameters.

    Parameters
    ----------
    shape : Rectangle, Circle, Link, or Plate
        Shape instance to check.

    Raises
    ------
    ValueError
        If any dimension is non-positive, or a Plate is not CCW / has
        zero (or negative) area.
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
    elif isinstance(shape, Plate):
        raw = shape.input_vertices
        if np.any(np.isnan(raw)):
            raise ValueError("Plate vertices contain NaN values")
        (x0, y0), (x1, y1), (x2, y2) = raw
        signed_area_2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if signed_area_2 <= 0:
            raise ValueError(
                "Plate vertices must be in counter-clockwise order "
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
    shape : Rectangle, Circle, Link, or Plate
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

    if isinstance(shape, Plate):
        # vertices are already re-centred on the CoM (centroid)
        verts = shape.vertices   # shape (3, 2)
        signed_area_2 = 0.0
        Jz = 0.0
        for i in range(3):
            xi, yi = verts[i]
            xj, yj = verts[(i + 1) % 3]
            cross = xi * yj - xj * yi
            signed_area_2 += cross
            Jz += (xi*xi + xi*xj + xj*xj + yi*yi + yi*yj + yj*yj) * cross
        area = abs(signed_area_2) / 2.0
        mass = rho * area * tz
        inertia = mass * abs(Jz) / (12.0 * area) if area > 0 else 0.0
        return mass, inertia

    raise TypeError(f"Unsupported shape type: {type(shape).__name__}")
