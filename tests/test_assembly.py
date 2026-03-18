"""
Unit tests for the topological assembler (builder.py).

Verifies that _assemble() correctly computes body positions and orientations
from joint connectivity and q0 values.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:  python -m pytest PMD/tests/test_assembly.py -v
"""
import pytest
import numpy as np

from PMD.src.model import Body, Ground, _GroundType
from PMD.src.constraints import Joint
from PMD.src.builder import _assemble


@pytest.fixture(autouse=True)
def _reset():
    """Reset Ground markers and Body counter between tests."""
    yield
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]
    Body.COUNT = 0


class TestRevJointAssembly:
    """Test assembly with revolute joints."""

    def test_single_pendulum(self):
        """Single body attached to Ground by a rev joint at the origin."""
        B1 = Body(m=1, J=1)
        mk_G = Ground.add_marker([0.0, 0.0])
        mk_B = B1.add_marker([0.5, 0.0])  # pivot at left end of body
        j = Joint(type='rev', iMarker=mk_B, jMarker=mk_G, q0=-np.pi / 2)
        _assemble([B1], [j])
        # Body hangs straight down: r should place the CoM at (0, -0.5)
        # with p = -pi/2
        np.testing.assert_allclose(B1.p, -np.pi / 2, atol=1e-12)
        expected_r = np.array([[0.0], [-0.5]]) - B1._A @ np.array([[0.5], [0.0]]) + np.array([[0.0], [-0.5]])
        # Actually: joint at origin => B1.r = r_joint - A(p)*sP_child
        A = np.array([[np.cos(-np.pi/2), -np.sin(-np.pi/2)],
                       [np.sin(-np.pi/2),  np.cos(-np.pi/2)]])
        expected_r = np.array([[0.0], [0.0]]) - A @ np.array([[0.5], [0.0]])
        np.testing.assert_allclose(B1.r, expected_r, atol=1e-12)

    def test_double_pendulum(self):
        """Two bodies in series via rev joints."""
        B1 = Body(m=1, J=1)
        B2 = Body(m=1, J=1)
        mk_G = Ground.add_marker([0.0, 0.0])
        mk_B1_top = B1.add_marker([0.0, 0.5])
        mk_B1_bot = B1.add_marker([0.0, -0.5])
        mk_B2_top = B2.add_marker([0.0, 0.5])
        j0 = Joint(type='rev', iMarker=mk_B1_top, jMarker=mk_G, q0=0)
        j1 = Joint(type='rev', iMarker=mk_B2_top, jMarker=mk_B1_bot, q0=0)
        _assemble([B1, B2], [j0, j1])
        np.testing.assert_allclose(B1.p, 0.0, atol=1e-12)
        np.testing.assert_allclose(B2.p, 0.0, atol=1e-12)
        np.testing.assert_allclose(B1.r, np.array([[0.0], [-0.5]]), atol=1e-12)
        np.testing.assert_allclose(B2.r, np.array([[0.0], [-1.5]]), atol=1e-12)


class TestTranJointAssembly:
    """Test assembly with translational joints."""

    def test_slider(self):
        """Body sliding along horizontal axis from Ground."""
        B1 = Body(m=1, J=1)
        mk_G = Ground.add_marker([0.0, 0.0], theta=0.0)   # horizontal
        mk_B = B1.add_marker([0.0, 0.0], theta=0.0)
        j = Joint(type='tran', iMarker=mk_B, jMarker=mk_G, q0=2.0)
        _assemble([B1], [j])
        np.testing.assert_allclose(B1.p, 0.0, atol=1e-12)
        np.testing.assert_allclose(B1.r, np.array([[2.0], [0.0]]), atol=1e-12)


class TestRevTranJointAssembly:
    """Assembly involving rev-rev joints."""

    def test_coupler(self):
        """Single rev-rev coupler between Ground and a body."""
        B1 = Body(m=1, J=1, p=0.0)
        mk_G = Ground.add_marker([0.0, 0.0])
        mk_B = B1.add_marker([0.0, 0.0])
        j = Joint(type='rev-rev', iMarker=mk_B, jMarker=mk_G, L=1.0, q0=np.pi / 2)
        _assemble([B1], [j])
        # Child pivot at distance L along angle q0 from parent pivot
        expected_r = np.array([[0.0], [1.0]])
        np.testing.assert_allclose(B1.r, expected_r, atol=1e-12)
