"""Tests for the CasADi DAE dynamic solver.

These tests verify:
  1. Pendulum dynamics: constraint satisfaction at every step.
  2. Constraint satisfaction: Phi(q_k) ≈ 0 at every output step.
  3. Energy consistency: a conservative pendulum preserves total energy.
  4. Output shape: T, uT conform to the expected format.
  5. Result containers: Body and Joint result containers are populated.
  6. Driven crank: a fully kinematic model works as a DAE too.
  7. Spring–damper force: PtpForce handled correctly.
  8. Torque force: constant Torque handled correctly.
  9. RotSdaForce: torsional spring–damper handled correctly.
 10. No constraints: pure ODE fallback (no joints) works.
 11. UserForce: callback-based force handled correctly.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:  python -m pytest pmd/tests/test_solver.py -v
"""

import numpy as np
import pytest

from pmd.core.model import Body, Ground, _GroundType
from pmd.core.constraints import (
    RevJoint, Weight, RelRotJoint, Function,
    PtpForce, RotSdaForce, Torque, UserForce,
)
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


# ---------------------------------------------------------------------------
# Model builders
# ---------------------------------------------------------------------------

def _make_pendulum():
    """Single pendulum: 1 body, 1 RevJoint, Weight.  DOF = 1."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, -0.5], orientation=0.0)
    mk_B = B.add_marker([0.0, 0.5])
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    model = PlanarMultibodyModel(bodies=[B], joints=[j], forces=[w])
    return model, B


def _make_driven_crank():
    """Crank driven by RelRotJoint: DOF = 0."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, 0.0], orientation=0.0)
    mk_pivot = B.add_marker([0.0, 0.0])
    mk_G = Ground.add_marker([0.0, 0.0])
    j_pin = RevJoint(iMarker=mk_pivot, jMarker=mk_G)
    f = Function(type='a', coeff=[0.0, 1.0, 0.0])
    j_rot = RelRotJoint(iMarker=mk_pivot, jMarker=mk_G, iFunct=f)
    model = PlanarMultibodyModel(
        bodies=[B], joints=[j_pin, j_rot], forces=[], functions=[f]
    )
    return model, B


def _make_pendulum_with_spring():
    """Pendulum with PtpForce spring–damper.  DOF = 1."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, -0.5], orientation=0.0)
    mk_B = B.add_marker([0.0, 0.5])
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    # Spring from body tip to a fixed point
    mk_tip = B.add_marker([0.0, -0.5])
    mk_anchor = Ground.add_marker([1.0, 0.0])
    spring = PtpForce(iMarker=mk_tip, jMarker=mk_anchor,
                       k=10.0, dc=0.5, L0=1.0, f_a=0.0)
    model = PlanarMultibodyModel(bodies=[B], joints=[j],
                                  forces=[w, spring])
    return model, B


def _make_free_body():
    """Body with no joints, only Weight.  nC = 0 → pure ODE."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, 0.0], orientation=0.0)
    w = Weight()
    model = PlanarMultibodyModel(bodies=[B], joints=[], forces=[w])
    return model, B


def _make_pendulum_with_torque():
    """Pendulum with constant Torque.  DOF = 1."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, -0.5], orientation=0.0)
    mk_B = B.add_marker([0.0, 0.5])
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    trq = Torque(iBody=B, torque_value=0.5)
    model = PlanarMultibodyModel(bodies=[B], joints=[j], forces=[w, trq])
    return model, B


def _make_pendulum_with_rotsda():
    """Pendulum with RotSdaForce torsional spring–damper.  DOF = 1."""
    B = Body(mass=1.0, inertia=0.1, position=[0.0, -0.5], orientation=0.0)
    mk_B = B.add_marker([0.0, 0.5])
    mk_G = Ground.add_marker([0.0, 0.0])
    j = RevJoint(iMarker=mk_B, jMarker=mk_G)
    w = Weight()
    sda = RotSdaForce(iBody=B, jBody=Ground, k=5.0, dc=0.2,
                       theta0=0.0, T_a=0.0)
    model = PlanarMultibodyModel(bodies=[B], joints=[j], forces=[w, sda])
    return model, B


# ---------------------------------------------------------------------------
# 2. Constraint satisfaction
# ---------------------------------------------------------------------------

class TestConstraintSatisfaction:

    def test_constraints_satisfied_pendulum(self):
        """max |Phi(q_k)| < 1e-7 at every output step."""
        model, _ = _make_pendulum()
        T, uT = model.solve(t_final=1.0, dt=0.01)

        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            model.t = float(T[k])
            model._u2bodies(uT[k])
            model._update_position()
            phi = model._compute_constraints()
            max_phi = max(max_phi, float(np.linalg.norm(phi)))

        assert max_phi < 1e-7, (
            f"Constraint violation too large: max |Phi| = {max_phi:.3e}"
        )

    def test_constraints_satisfied_driven_crank(self):
        """Driven crank with DOF=0 must satisfy constraints tightly."""
        model, _ = _make_driven_crank()
        T, uT = model.solve(t_final=1.0, dt=0.01)

        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            model.t = float(T[k])
            model._u2bodies(uT[k])
            model._update_position()
            phi = model._compute_constraints()
            max_phi = max(max_phi, float(np.linalg.norm(phi)))

        assert max_phi < 1e-7, (
            f"Driven crank constraint violation: max |Phi| = {max_phi:.3e}"
        )


# ---------------------------------------------------------------------------
# 4. Energy consistency (undamped pendulum)
# ---------------------------------------------------------------------------

class TestEnergyConsistency:

    def test_energy_conserved(self):
        """Total energy of an undamped pendulum must not drift by > 1%."""
        model, B = _make_pendulum()
        T, uT = model.solve(t_final=2.0, dt=0.01)

        nB3 = 3 * model.nB
        g = 9.81
        energies = []
        for k in range(len(T)):
            q = uT[k, :nB3]
            v = uT[k, nB3:]
            y = q[1]
            KE = 0.5 * B.mass * (v[0]**2 + v[1]**2) + 0.5 * B.inertia * v[2]**2
            PE = B.mass * g * y
            energies.append(KE + PE)

        energies = np.array(energies)
        E0 = energies[0]
        drift = np.max(np.abs(energies - E0)) / abs(E0) if E0 != 0 else 0.0
        assert drift < 0.01, f"Energy drift = {drift:.4f} (> 1%)"


# ---------------------------------------------------------------------------
# 5. Output shape
# ---------------------------------------------------------------------------

class TestOutputShape:

    def test_shape_pendulum(self):
        model, _ = _make_pendulum()
        T, uT = model.solve(t_final=0.5, dt=0.1)

        nB3 = 3 * model.nB
        assert uT.shape == (len(T), 2 * nB3)
        assert len(T) > 1

    def test_t_eval(self):
        """Provide explicit t_eval array."""
        model, _ = _make_pendulum()
        t_eval = np.array([0.0, 0.1, 0.25, 0.5])
        T, uT = model.solve(t_eval=t_eval)

        np.testing.assert_allclose(T, t_eval, atol=1e-12)
        assert uT.shape[0] == len(t_eval)


# ---------------------------------------------------------------------------
# 6. Result containers
# ---------------------------------------------------------------------------

class TestResultContainers:

    def test_body_result_container(self):
        model, B = _make_pendulum()
        model.solve(t_final=0.3, dt=0.1)

        rc = B._result_container
        for key in ("positions", "velocities", "accelerations"):
            assert key in rc, f"Missing key '{key}'"

    def test_joint_result_container(self):
        model, _ = _make_pendulum()
        model.solve(t_final=0.3, dt=0.1)

        for joint in model.Joints:
            assert 'reactions' in joint._result_container


# ---------------------------------------------------------------------------
# 7. Driven crank as DAE
# ---------------------------------------------------------------------------

class TestDrivenCrank:

    def test_driven_crank_phi_matches(self):
        """Angle should follow phi(t) = t for the driven crank."""
        model, B = _make_driven_crank()
        T, uT = model.solve(t_final=1.0, dt=0.1)

        nB3 = 3 * model.nB
        phi_vals = uT[:, 2]  # orientation of single body
        np.testing.assert_allclose(phi_vals, T, atol=1e-4,
                                   err_msg="Driven crank angle != t")


# ---------------------------------------------------------------------------
# 8. PtpForce spring–damper
# ---------------------------------------------------------------------------

class TestPtpForce:

    def test_spring_damper_runs(self):
        model, _ = _make_pendulum_with_spring()
        T, uT = model.solve(t_final=0.5, dt=0.05)
        assert len(T) > 1

    def test_spring_constraints_satisfied(self):
        model, _ = _make_pendulum_with_spring()
        T, uT = model.solve(t_final=0.5, dt=0.05)

        nB3 = 3 * model.nB
        max_phi = 0.0
        for k in range(len(T)):
            model.t = float(T[k])
            model._u2bodies(uT[k])
            model._update_position()
            phi = model._compute_constraints()
            max_phi = max(max_phi, float(np.linalg.norm(phi)))

        assert max_phi < 1e-7


# ---------------------------------------------------------------------------
# 9. Torque
# ---------------------------------------------------------------------------

class TestTorque:

    def test_torque_runs(self):
        model, _ = _make_pendulum_with_torque()
        T, uT = model.solve(t_final=0.5, dt=0.05)
        assert len(T) > 1


# ---------------------------------------------------------------------------
# 10. RotSdaForce
# ---------------------------------------------------------------------------

class TestRotSdaForce:

    def test_rotsda_runs(self):
        model, _ = _make_pendulum_with_rotsda()
        T, uT = model.solve(t_final=0.5, dt=0.05)
        assert len(T) > 1


# ---------------------------------------------------------------------------
# 11. No constraints (pure ODE fallback)
# ---------------------------------------------------------------------------

class TestNoConstraints:

    def test_free_body_gravity(self):
        """Free body under gravity: y(t) = -g/2 * t^2.  No joints."""
        model, B = _make_free_body()
        T, uT = model.solve(t_final=0.5, dt=0.05)

        g = 9.81
        y_expected = -0.5 * g * T**2
        y_actual = uT[:, 1]

        np.testing.assert_allclose(y_actual, y_expected, atol=1e-5,
                                   err_msg="Free-body gravity y(t) mismatch")


# ---------------------------------------------------------------------------
# 12. UserForce (callback-based force via numeric parameters)
# ---------------------------------------------------------------------------

def _make_quarter_car():
    """Quarter-car model (AA) with a UserForce wheel contact."""
    _GroundType._instance._markers = [_GroundType._instance.origin]
    Body.COUNT = 0

    B1 = Body(mass=2, inertia=0.5, position=[0.4398, 0.2512], orientation=-0.0367)
    B2 = Body(mass=30, inertia=2.5, position=[0.6817, 0.3498], orientation=0.0783)
    B3 = Body(mass=1, inertia=0.5, position=[0.4463, 0.4308], orientation=6.5222)

    q1 = B1.add_marker([-0.24, 0.0])
    a1 = B1.add_marker([0.18, 0.0])
    a2 = B2.add_marker([-0.07, -0.10])
    b2 = B2.add_marker([-0.10, 0.12])
    b3 = B3.add_marker([0.13, 0.0])
    o3 = B3.add_marker([-0.13, 0.0])
    o0 = Ground.add_marker([0.32, 0.40])
    q0 = Ground.add_marker([0.20, 0.26])
    e1 = B1.add_marker([0.0, 0.0])
    f0 = Ground.add_marker([0.38, 0.43])

    j1 = RevJoint(iMarker=q1, jMarker=q0)
    j2 = RevJoint(iMarker=a1, jMarker=a2)
    j3 = RevJoint(iMarker=b2, jMarker=b3)
    j4 = RevJoint(iMarker=o3, jMarker=o0)

    s1 = PtpForce(iMarker=e1, jMarker=f0, k=90000, L0=0.23, dc=1100)

    _k_wh, _L0_wh, _dc_wh = 50000, 0.35, 1000

    def wheel_contact():
        dely = B2.position[1] - _L0_wh
        if dely < 0:
            fy = (_k_wh * dely + _dc_wh * B2.velocity[1]).item()
            return [{'body': B2, 'force': [0, -fy], 'torque': 0}]
        return []

    s2 = UserForce(callback=wheel_contact)
    s3 = Weight()

    model = PlanarMultibodyModel(
        bodies=[B1, B2, B3],
        joints=[j1, j2, j3, j4],
        forces=[s1, s2, s3])
    return model, B2


class TestUserForce:

    def test_user_force_runs(self):
        """Quarter car with UserForce completes without error."""
        model, B2 = _make_quarter_car()
        T, uT = model.solve(t_final=1.0,
                             t_eval=np.linspace(0, 1, 101),
                             ic_correct=True)
        assert T.shape[0] == 101
        assert uT.shape == (101, 18)

    def test_user_force_constraints_satisfied(self):
        """Constraint violation stays small with UserForce."""
        model, _ = _make_quarter_car()
        T, uT = model.solve(t_final=0.5,
                             t_eval=np.linspace(0, 0.5, 51),
                             ic_correct=True)
        nB3 = 3 * model.nB
        for k in range(len(T)):
            q_k = uT[k, :nB3]
            for Bi, body in enumerate(model.Bodies):
                ir = 3 * Bi
                body.position = q_k[ir:ir + 2].reshape(2, 1)
                body.orientation = float(q_k[ir + 2])
            model.t = float(T[k])
            model._update_position()
            Phi_k = model._compute_constraints().flatten()
            assert np.max(np.abs(Phi_k)) < 1e-6, \
                f"Constraint violation {np.max(np.abs(Phi_k)):.2e} at step {k}"
