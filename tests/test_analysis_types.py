"""
Tests for solve(analysis=...) parameter and new analysis types.

These tests verify:
  1. Backward compatibility: solve() without 'analysis' behaves identically
     to solve(analysis="dynamic").
  2. Input validation: unknown analysis type raises ValueError.
  3. DOF check: kinematic analysis on a model with DOF≠0 raises ValueError.
  4. Kinematic sanity: constraints satisfied at every step.
  5. Static sanity: equilibrium residual is below tolerance.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:  python -m pytest pmd/tests/test_analysis_types.py -v
"""

import numpy as np
import pytest

from pmd.core.model import Body, Ground, _GroundType
from pmd.core.constraints import RevJoint
from pmd.core.motion import RotMotion
from pmd.core.forces import Weight, Function
from pmd.core.solver import PlanarMultibodyModel


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset():
    """Reset global state between tests."""
    yield
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]
    Body.COUNT = 0


def _make_pendulum():
    """Simple single pendulum: 1 body, 1 RevJoint at the pivot, Weight.

    DOF = 3*1 - 2 = 1  (one rotational DOF).
    """
    B = Body(mass=1.0, inertia=0.1, position=[0.0, -0.5], orientation=0.0)
    mk_B = B.add_marker([0.0, 0.5])   # pivot at top of body
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    model = PlanarMultibodyModel(bodies=[B], joints=[j], forces=[w])
    return model, B


def _make_driven_crank():
    """Single body driven by a RotMotion — fully kinematic (DOF = 0).

    The crank rotates at constant angular velocity via Function type 'a'.
    DOF = 3*1 - 3 = 0.
    """
    B = Body(mass=1.0, inertia=0.1, position=[0.0, 0.0], orientation=0.0)
    mk_pivot = B.add_marker([0.0, 0.0])  # body origin = pivot
    mk_G = Ground.add_marker([0.0, 0.0])

    # Pin the body at origin (2 constraints)
    j_pin = RevJoint(iMarker=mk_pivot, jMarker=mk_G)

    # Drive rotation: phi(t) = 0 + 1*t  (omega = 1 rad/s)
    # Function type 'a': f(t) = c0 + c1*t + c2*t^2
    f = Function(type='a', coeff=[0.0, 1.0, 0.0])  # f = t
    j_rot = RotMotion(iMarker=mk_pivot, jMarker=mk_G, iFunct=f)

    model = PlanarMultibodyModel(
        bodies=[B], joints=[j_pin, j_rot], forces=[], functions=[f]
    )
    return model, B


# ---------------------------------------------------------------------------
# 1. Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """solve() with no 'analysis' arg must equal solve(analysis='dynamic')."""

    def test_default_equals_dynamic(self):
        model1, _ = _make_pendulum()
        T1, uT1 = model1.solve(t_final=0.5, dt=0.1)

        # Reset and rebuild identical model
        _GroundType._instance._markers = [_GroundType._instance.origin]
        Body.COUNT = 0
        model2, _ = _make_pendulum()
        T2, uT2 = model2.solve(analysis="dynamic", t_final=0.5, dt=0.1)

        np.testing.assert_array_equal(T1, T2, err_msg="T differs")
        np.testing.assert_array_equal(uT1, uT2, err_msg="uT differs")

    def test_analysis_param_case_insensitive(self):
        """'Dynamic', 'DYNAMIC', 'dynamic' are all accepted."""
        for variant in ("Dynamic", "DYNAMIC", "dynamic"):
            _GroundType._instance._markers = [_GroundType._instance.origin]
            Body.COUNT = 0
            model, _ = _make_pendulum()
            T, uT = model.solve(analysis=variant, t_final=0.1, dt=0.1)
            assert len(T) > 0


# ---------------------------------------------------------------------------
# 2. Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_unknown_analysis_raises(self):
        model, _ = _make_pendulum()
        with pytest.raises(ValueError, match="Unknown analysis type"):
            model.solve(analysis="magic", t_final=1.0, dt=0.1)

    def test_kinematic_on_underconstrained_model_raises(self):
        """DOF = 1 model must raise ValueError for kinematic analysis."""
        model, _ = _make_pendulum()  # DOF = 1
        with pytest.raises(ValueError, match="DOF = 0"):
            model.solve(analysis="kinematic", t_final=1.0, dt=0.1)


# ---------------------------------------------------------------------------
# 3. nDOF property
# ---------------------------------------------------------------------------

class TestNDOF:

    def test_pendulum_dof(self):
        model, _ = _make_pendulum()
        # 3*1 bodies - 2 constraints (RevJoint) = 1
        assert model.nDOF == 1

    def test_driven_crank_dof(self):
        model, _ = _make_driven_crank()
        # 3*1 - (2 + 1) = 0
        assert model.nDOF == 0


# ---------------------------------------------------------------------------
# 4. Kinematic analysis sanity
# ---------------------------------------------------------------------------

class TestKinematicAnalysis:

    def test_constraints_satisfied_at_all_steps(self):
        """max |Phi(q_k)| < 1e-9 at every output step."""
        model, B = _make_driven_crank()
        T, uT = model.solve(analysis="kinematic", t_final=1.0, dt=0.01,
                            ic_correct=True)
        assert len(T) > 1

        # Re-check constraints for each stored state
        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            model.t = float(T[k])          # <-- must set time before eval
            u_k = uT[k]
            model._u2bodies(u_k)
            model._update_position()
            phi = model._compute_constraints()
            max_phi = max(max_phi, float(np.linalg.norm(phi)))

        assert max_phi < 1e-9, (
            f"Kinematic constraints not satisfied: max |Phi| = {max_phi:.3e}"
        )

    def test_output_shape(self):
        model, _ = _make_driven_crank()
        T, uT = model.solve(analysis="kinematic", t_final=0.5, dt=0.05)
        nB3 = 3 * model.nB
        assert uT.shape == (len(T), 2 * nB3)

    def test_result_container_populated(self):
        """Body._result_container must have positions/velocities/accelerations."""
        model, B = _make_driven_crank()
        model.solve(analysis="kinematic", t_final=0.3, dt=0.05)
        rc = B._result_container
        for key in ("positions", "velocities", "accelerations"):
            assert key in rc, f"Missing key '{key}' in _result_container"


# ---------------------------------------------------------------------------
# 5. Static analysis sanity
# ---------------------------------------------------------------------------

class TestStaticAnalysis:

    def test_equilibrium_residual(self):
        """||Phi_q^T lam - h_a|| < 1e-8 and ||Phi|| < 1e-10 at equilibrium.

        A single pendulum under gravity.  The equilibrium is the body hanging
        straight down (body CoM directly below the pivot); constraints must be
        satisfied and the force balance must hold.
        """
        model, B = _make_pendulum()
        T, uT = model.solve(analysis="static", ic_correct=True)

        assert T.shape == (1,)
        assert float(T[0]) == 0.0

        # Unpack equilibrium state
        model.t = 0.0
        model._u2bodies(uT[0])
        model._update_position()

        # Constraints must be satisfied
        phi_eq = model._compute_constraints()
        assert float(np.linalg.norm(phi_eq)) < 1e-10, (
            f"Constraints violated at static equilibrium: "
            f"||Phi|| = {float(np.linalg.norm(phi_eq)):.3e}"
        )

        # Force balance: D^T lam = h_a
        D = model._compute_jacobian()
        h_a = model._compute_force().flatten()
        nConst = model.nC
        lam_vec = np.zeros(nConst)
        for joint in model.Joints:
            rs = joint._rows
            re = joint._rowe
            lam_vec[rs:re] = joint._result_container['reactions'].flatten()

        residual = float(np.linalg.norm(D.T @ lam_vec - h_a))
        assert residual < 1e-8, (
            f"Static equilibrium force residual too large: {residual:.3e}"
        )

    def test_static_output_shape(self):
        model, _ = _make_pendulum()
        T, uT = model.solve(analysis="static")
        assert T.shape == (1,)
        assert uT.shape[0] == 1

    def test_static_velocities_zero(self):
        """Velocity component of uT must be zero for static result."""
        model, _ = _make_pendulum()
        T, uT = model.solve(analysis="static")
        nB3 = 3 * model.nB
        velocities = uT[0, nB3:]
        np.testing.assert_array_equal(velocities, np.zeros(nB3))
