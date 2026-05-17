"""Tests for pmd.core.shapes: geometry, mass properties, validation."""

import math

import numpy as np
import pytest

from pmd.core.shapes import (
    Circle,
    Link,
    Plate,
    Rectangle,
    _validate_shape,
    compute_mass_props,
)
from pmd.core.model import Body


# ---------------------------------------------------------------------------
# compute_mass_props — analytical checks
# ---------------------------------------------------------------------------

class TestComputeMassPropsRectangle:
    def test_mass(self):
        shape = Rectangle(width=0.4, height=0.3)
        mass, _ = compute_mass_props(shape, density=7850.0, thickness_z=0.01)
        expected_mass = 7850.0 * (0.4 * 0.3) * 0.01
        assert math.isclose(mass, expected_mass, rel_tol=1e-12)

    def test_inertia(self):
        w, h = 0.4, 0.3
        shape = Rectangle(width=w, height=h)
        mass, inertia = compute_mass_props(shape, density=7850.0, thickness_z=0.01)
        expected_J = mass * (w**2 + h**2) / 12.0
        assert math.isclose(inertia, expected_J, rel_tol=1e-12)

    def test_default_thickness(self):
        shape = Rectangle(width=1.0, height=1.0)
        mass, _ = compute_mass_props(shape, density=1000.0)
        assert math.isclose(mass, 1000.0 * 1.0 * 0.01, rel_tol=1e-12)


class TestComputeMassPropsLink:
    def test_mass(self):
        L, t = 0.40, 0.015
        shape = Link(length=L, thickness=t)
        mass, _ = compute_mass_props(shape, density=7850.0, thickness_z=0.015)
        expected_mass = 7850.0 * (L * t) * 0.015
        assert math.isclose(mass, expected_mass, rel_tol=1e-12)

    def test_inertia(self):
        L, t = 0.40, 0.015
        shape = Link(length=L, thickness=t)
        mass, inertia = compute_mass_props(shape, density=7850.0, thickness_z=0.015)
        expected_J = mass * L**2 / 12.0
        assert math.isclose(inertia, expected_J, rel_tol=1e-12)


class TestComputeMassPropsCircle:
    def test_mass(self):
        r = 0.05
        shape = Circle(radius=r)
        mass, _ = compute_mass_props(shape, density=2700.0, thickness_z=0.02)
        expected_mass = 2700.0 * math.pi * r**2 * 0.02
        assert math.isclose(mass, expected_mass, rel_tol=1e-12)

    def test_inertia(self):
        r = 0.05
        shape = Circle(radius=r)
        mass, inertia = compute_mass_props(shape, density=2700.0, thickness_z=0.02)
        assert math.isclose(inertia, 0.5 * mass * r**2, rel_tol=1e-12)


class TestComputeMassPropsPlate:
    """Right-isosceles triangle with legs = 1 (CCW), centred on centroid.

    Vertices before centering: (0,0), (1,0), (0,1).
    Centroid: (1/3, 1/3).
    Recentred:  v0=(-1/3,-1/3), v1=(2/3,-1/3), v2=(-1/3, 2/3).
    Area = 0.5.
    """

    def _plate(self):
        return Plate(v1=(0.0, 0.0), v2=(1.0, 0.0), v3=(0.0, 1.0))

    def test_mass(self):
        plate = self._plate()
        mass, _ = compute_mass_props(plate, density=1000.0, thickness_z=0.01)
        expected_mass = 1000.0 * 0.5 * 0.01
        assert math.isclose(mass, expected_mass, rel_tol=1e-12)

    def test_inertia_positive(self):
        plate = self._plate()
        _, inertia = compute_mass_props(plate, density=1000.0, thickness_z=0.01)
        assert inertia > 0.0

    def test_inertia_analytical(self):
        # Centroidal J_z for right-isosceles triangle, legs a=b=1:
        # J_z = m * (a^2 + b^2) / 18  (well-known formula for a right triangle)
        # Here a = b = 1, m = rho * 0.5 * tz.
        plate = self._plate()
        mass, inertia = compute_mass_props(plate, density=1000.0, thickness_z=0.01)
        expected_J = mass * (1.0**2 + 1.0**2) / 18.0
        assert math.isclose(inertia, expected_J, rel_tol=1e-10)

    def test_from_world_points_mass(self):
        plate, _ = Plate.from_world_points((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        mass, _ = compute_mass_props(plate, density=1000.0, thickness_z=0.01)
        assert math.isclose(mass, 1000.0 * 0.5 * 0.01, rel_tol=1e-12)

    def test_from_world_points_centroid(self):
        plate, centroid = Plate.from_world_points((0.0, 0.0), (3.0, 0.0),
                                                  (0.0, 3.0))
        assert math.isclose(centroid[0], 1.0, rel_tol=1e-12)
        assert math.isclose(centroid[1], 1.0, rel_tol=1e-12)
        # Stored vertices must be centred on origin
        assert math.isclose(plate.vertices.mean(axis=0)[0], 0.0, abs_tol=1e-12)
        assert math.isclose(plate.vertices.mean(axis=0)[1], 0.0, abs_tol=1e-12)

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError):
            compute_mass_props("not_a_shape", density=1000.0)


# ---------------------------------------------------------------------------
# Body(shape=..., density=...) — auto-derived mass properties
# ---------------------------------------------------------------------------

class TestBodyShapeConstructor:
    def test_rectangle_mass_inertia(self):
        shape = Rectangle(width=0.4, height=0.3)
        density = 7850.0
        tz = 0.01
        body = Body(shape=shape, density=density, thickness_z=tz)
        expected_mass = 7850.0 * 0.4 * 0.3 * 0.01
        expected_J = expected_mass * (0.4**2 + 0.3**2) / 12.0
        assert math.isclose(body.mass, expected_mass, rel_tol=1e-12)
        assert math.isclose(body.inertia, expected_J, rel_tol=1e-12)

    def test_shape_stored(self):
        shape = Link(length=0.1, thickness=0.015)
        body = Body(shape=shape, density=7850.0, thickness_z=0.01)
        assert body.shape is shape

    def test_name_forwarded(self):
        shape = Circle(radius=0.05)
        body = Body(shape=shape, density=2700.0, name="wheel")
        assert body.name == "wheel"

    def test_auto_matches_manual(self):
        shape = Rectangle(width=0.4, height=0.3)
        auto = Body(shape=shape, density=7850.0, thickness_z=0.01)
        m, J = compute_mass_props(shape, density=7850.0, thickness_z=0.01)
        manual = Body(mass=m, inertia=J, shape=shape)
        assert math.isclose(auto.mass, manual.mass, rel_tol=1e-12)
        assert math.isclose(auto.inertia, manual.inertia, rel_tol=1e-12)

    def test_mass_override_keeps_user_values(self):
        shape = Rectangle(width=0.4, height=0.3)
        body = Body(shape=shape, density=7850.0,
                    mass=99.0, inertia=9.9, mass_override=True)
        assert math.isclose(body.mass, 99.0)
        assert math.isclose(body.inertia, 9.9)
        assert body.mass_override is True

    def test_density_warns_when_mass_also_given(self):
        shape = Rectangle(width=0.4, height=0.3)
        with pytest.warns(UserWarning, match="mass_override"):
            Body(shape=shape, density=7850.0, mass=5.0, inertia=0.1)

    def test_shape_without_density_requires_mass_inertia(self):
        with pytest.raises(ValueError):
            Body(shape=Rectangle(width=0.4, height=0.3))

    def test_no_shape_requires_mass_inertia(self):
        with pytest.raises(ValueError):
            Body(mass=1.0)  # inertia missing


# ---------------------------------------------------------------------------
# _validate_shape
# ---------------------------------------------------------------------------

class TestValidateShape:
    def test_rectangle_ok(self):
        _validate_shape(Rectangle(width=1.0, height=0.5))  # no raise

    def test_rectangle_zero_width_raises(self):
        with pytest.raises(ValueError):
            _validate_shape(Rectangle(width=0.0, height=0.5))

    def test_rectangle_negative_height_raises(self):
        with pytest.raises(ValueError):
            _validate_shape(Rectangle(width=1.0, height=-0.1))

    def test_circle_ok(self):
        _validate_shape(Circle(radius=0.1))  # no raise

    def test_circle_zero_radius_raises(self):
        with pytest.raises(ValueError):
            _validate_shape(Circle(radius=0.0))

    def test_link_ok(self):
        _validate_shape(Link(length=0.4, thickness=0.01))  # no raise

    def test_link_zero_thickness_raises(self):
        with pytest.raises(ValueError):
            _validate_shape(Link(length=0.4, thickness=0.0))

    def test_plate_ccw_ok(self):
        _validate_shape(Plate(v1=(0.0, 0.0), v2=(1.0, 0.0),
                              v3=(0.0, 1.0)))  # no raise

    def test_plate_cw_raises(self):
        # CW: same triangle in reverse order
        with pytest.raises(ValueError, match="counter-clockwise"):
            _validate_shape(Plate(v1=(0.0, 0.0), v2=(0.0, 1.0),
                                  v3=(1.0, 0.0)))

    def test_plate_degenerate_collinear_raises(self):
        # Three collinear points → zero area → CW or zero signed area
        with pytest.raises(ValueError):
            _validate_shape(Plate(v1=(0.0, 0.0), v2=(1.0, 0.0),
                                  v3=(2.0, 0.0)))

    def test_plate_nan_raises(self):
        with pytest.raises(ValueError, match="NaN"):
            _validate_shape(Plate(v1=(0.0, 0.0), v2=(1.0, 0.0),
                                  v3=(float("nan"), 1.0)))

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError):
            _validate_shape("rectangle")


# ---------------------------------------------------------------------------
# Body.__init__ shape validation wiring
# ---------------------------------------------------------------------------

class TestBodyShapeValidation:
    def test_body_with_valid_shape_ok(self):
        Body(mass=1.0, inertia=0.1, shape=Rectangle(width=0.4, height=0.3))

    def test_body_with_cw_plate_raises(self):
        with pytest.raises(ValueError):
            Body(mass=1.0, inertia=0.1,
                 shape=Plate(v1=(0.0, 0.0), v2=(0.0, 1.0), v3=(1.0, 0.0)))

    def test_body_without_shape_ok(self):
        Body(mass=1.0, inertia=0.1)  # no raise
