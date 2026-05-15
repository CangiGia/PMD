"""Tests for pmd.core.shapes: geometry, mass properties, validation."""

import math

import numpy as np
import pytest

from pmd.core.shapes import (
    Circle,
    Link,
    Polygon,
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


class TestComputeMassPropsPolygon:
    def _unit_square(self):
        """1×1 square, CCW, centred at origin."""
        return Polygon(vertices=[[-0.5, -0.5], [0.5, -0.5],
                                  [0.5,  0.5], [-0.5,  0.5]])

    def test_mass_unit_square(self):
        shape = self._unit_square()
        mass, _ = compute_mass_props(shape, density=1000.0, thickness_z=0.01)
        assert math.isclose(mass, 1000.0 * 1.0 * 0.01, rel_tol=1e-10)

    def test_inertia_unit_square(self):
        # Centroidal J_z of 1×1 square = m*(1²+1²)/12
        shape = self._unit_square()
        mass, inertia = compute_mass_props(shape, density=1000.0, thickness_z=0.01)
        expected_J = mass * (1.0**2 + 1.0**2) / 12.0
        assert math.isclose(inertia, expected_J, rel_tol=1e-10)

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError):
            compute_mass_props("not_a_shape", density=1000.0)


# ---------------------------------------------------------------------------
# Body.from_shape — round-trip
# ---------------------------------------------------------------------------

class TestBodyFromShape:
    def test_rectangle_mass_inertia(self):
        shape = Rectangle(width=0.4, height=0.3)
        density = 7850.0
        tz = 0.01
        body = Body.from_shape(shape, density=density, thickness_z=tz)
        expected_mass = 7850.0 * 0.4 * 0.3 * 0.01
        expected_J = expected_mass * (0.4**2 + 0.3**2) / 12.0
        assert math.isclose(body.mass, expected_mass, rel_tol=1e-12)
        assert math.isclose(body.inertia, expected_J, rel_tol=1e-12)

    def test_shape_stored(self):
        shape = Link(length=0.1, thickness=0.015)
        body = Body.from_shape(shape, density=7850.0, thickness_z=0.01)
        assert body.shape is shape

    def test_name_forwarded(self):
        shape = Circle(radius=0.05)
        body = Body.from_shape(shape, density=2700.0, name="wheel")
        assert body.name == "wheel"

    def test_matches_manual_body(self):
        shape = Rectangle(width=0.4, height=0.3)
        auto = Body.from_shape(shape, density=7850.0, thickness_z=0.01)
        m, J = compute_mass_props(shape, density=7850.0, thickness_z=0.01)
        manual = Body(mass=m, inertia=J, shape=shape)
        assert math.isclose(auto.mass, manual.mass, rel_tol=1e-12)
        assert math.isclose(auto.inertia, manual.inertia, rel_tol=1e-12)


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

    def test_polygon_ccw_ok(self):
        verts = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]
        _validate_shape(Polygon(vertices=verts))  # no raise

    def test_polygon_cw_raises(self):
        # CW = same square in reverse
        verts = [[-0.5, -0.5], [-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]]
        with pytest.raises(ValueError, match="counter-clockwise"):
            _validate_shape(Polygon(vertices=verts))

    def test_polygon_too_few_vertices_raises(self):
        with pytest.raises(ValueError, match="at least 3"):
            _validate_shape(Polygon(vertices=[[0, 0], [1, 0]]))

    def test_polygon_nan_raises(self):
        verts = [[0, 0], [1, 0], [float("nan"), 1]]
        with pytest.raises(ValueError, match="NaN"):
            _validate_shape(Polygon(vertices=verts))

    def test_unsupported_type_raises(self):
        with pytest.raises(TypeError):
            _validate_shape("rectangle")


# ---------------------------------------------------------------------------
# Body.__init__ shape validation wiring
# ---------------------------------------------------------------------------

class TestBodyShapeValidation:
    def test_body_with_valid_shape_ok(self):
        Body(mass=1.0, inertia=0.1, shape=Rectangle(width=0.4, height=0.3))

    def test_body_with_cw_polygon_raises(self):
        cw_verts = [[-0.5, -0.5], [-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]]
        with pytest.raises(ValueError):
            Body(mass=1.0, inertia=0.1, shape=Polygon(vertices=cw_verts))

    def test_body_without_shape_ok(self):
        Body(mass=1.0, inertia=0.1)  # no raise
