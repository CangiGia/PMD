"""
Tests for the CasADi + IDAS DAE solver (method="CASADI-DAE").

These tests verify:
  1. Import guard: solve() raises ImportError when CasADi is unavailable.
  2. Pendulum dynamics: CASADI-DAE agrees with LSODA within tolerance.
  3. Constraint satisfaction: Phi(q_k) ≈ 0 at every output step.
  4. Energy consistency: a conservative pendulum preserves total energy.
  5. Output shape: T, uT conform to the expected format.
  6. Result containers: Body and Joint result containers are populated.
  7. Driven crank: a fully kinematic model works as a DAE too.
  8. Spring–damper force: PtpForce handled correctly.
  9. Torque force: constant Torque handled correctly.
 10. RotSdaForce: torsional spring–damper handled correctly.
 11. No constraints: pure ODE fallback (no joints) works.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:  python -m pytest PMD/tests/test_dae_casadi.py -v
"""

import numpy as np
import pytest

_has_casadi = True
try:
    import casadi  # noqa: F401
except ImportError:
    _has_casadi = False

from PMD.src.model import Body, Ground, _GroundType
from PMD.src.constraints import (
    RevJoint, Weight, RelRotJoint, Function,
    PtpForce, RotSdaForce, Torque,
)
from PMD.src.solver import PlanarMultibodyModel

casadi_required = pytest.mark.skipif(
    not _has_casadi, reason="CasADi not installed"
)


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
# 1. Import guard
# ---------------------------------------------------------------------------

class TestImportGuard:

    def test_casadi_import_error(self, monkeypatch):
        """If casadi is not importable, solve must raise ImportError."""
        import builtins
        _original_import = builtins.__import__

        def _mock_import(name, *args, **kwargs):
            if name == 'casadi':
                raise ImportError("mocked")
            return _original_import(name, *args, **kwargs)

        model, _ = _make_pendulum()
        monkeypatch.setattr(builtins, '__import__', _mock_import)
        with pytest.raises(ImportError, match="CasADi is required"):
            model.solve(method="CASADI-DAE", t_final=0.1, dt=0.1)


# ---------------------------------------------------------------------------
# 2. Pendulum agreement with LSODA
# ---------------------------------------------------------------------------

class TestPendulumAgreement:

    @casadi_required
    def test_casadi_vs_lsoda(self):
        """CASADI-DAE should match LSODA within 1e-3 on a short run."""
        model1, _ = _make_pendulum()
        T1, uT1 = model1.solve(method="LSODA", t_final=0.5, dt=0.05)

        _GroundType._instance._markers = [_GroundType._instance.origin]
        Body.COUNT = 0

        model2, _ = _make_pendulum()
        T2, uT2 = model2.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)

        np.testing.assert_allclose(T1, T2, atol=1e-12)
        np.testing.assert_allclose(uT1, uT2, atol=1e-3,
                                   err_msg="CASADI-DAE vs LSODA mismatch")


# ---------------------------------------------------------------------------
# 3. Constraint satisfaction
# ---------------------------------------------------------------------------

class TestConstraintSatisfaction:

    @casadi_required
    def test_constraints_satisfied_pendulum(self):
        """max |Phi(q_k)| < 1e-7 at every output step."""
        model, _ = _make_pendulum()
        T, uT = model.solve(method="CASADI-DAE", t_final=1.0, dt=0.01)

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

    @casadi_required
    def test_constraints_satisfied_driven_crank(self):
        """Driven crank with DOF=0 must satisfy constraints tightly."""
        model, _ = _make_driven_crank()
        T, uT = model.solve(method="CASADI-DAE", t_final=1.0, dt=0.01)

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

    @casadi_required
    def test_energy_conserved(self):
        """Total energy of an undamped pendulum must not drift by > 1%."""
        model, B = _make_pendulum()
        T, uT = model.solve(method="CASADI-DAE", t_final=2.0, dt=0.01)

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

    @casadi_required
    def test_shape_pendulum(self):
        model, _ = _make_pendulum()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.1)

        nB3 = 3 * model.nB
        assert uT.shape == (len(T), 2 * nB3)
        assert len(T) > 1

    @casadi_required
    def test_t_eval(self):
        """Provide explicit t_eval array."""
        model, _ = _make_pendulum()
        t_eval = np.array([0.0, 0.1, 0.25, 0.5])
        T, uT = model.solve(method="CASADI-DAE", t_eval=t_eval)

        np.testing.assert_allclose(T, t_eval, atol=1e-12)
        assert uT.shape[0] == len(t_eval)


# ---------------------------------------------------------------------------
# 6. Result containers
# ---------------------------------------------------------------------------

class TestResultContainers:

    @casadi_required
    def test_body_result_container(self):
        model, B = _make_pendulum()
        model.solve(method="CASADI-DAE", t_final=0.3, dt=0.1)

        rc = B._result_container
        for key in ("positions", "velocities", "accelerations"):
            assert key in rc, f"Missing key '{key}'"

    @casadi_required
    def test_joint_result_container(self):
        model, _ = _make_pendulum()
        model.solve(method="CASADI-DAE", t_final=0.3, dt=0.1)

        for joint in model.Joints:
            assert 'reactions' in joint._result_container


# ---------------------------------------------------------------------------
# 7. Driven crank as DAE
# ---------------------------------------------------------------------------

class TestDrivenCrank:

    @casadi_required
    def test_driven_crank_phi_matches(self):
        """Angle should follow phi(t) = t for the driven crank."""
        model, B = _make_driven_crank()
        T, uT = model.solve(method="CASADI-DAE", t_final=1.0, dt=0.1)

        nB3 = 3 * model.nB
        phi_vals = uT[:, 2]  # orientation of single body
        np.testing.assert_allclose(phi_vals, T, atol=1e-4,
                                   err_msg="Driven crank angle != t")


# ---------------------------------------------------------------------------
# 8. PtpForce spring–damper
# ---------------------------------------------------------------------------

class TestPtpForce:

    @casadi_required
    def test_spring_damper_runs(self):
        model, _ = _make_pendulum_with_spring()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)
        assert len(T) > 1

    @casadi_required
    def test_spring_constraints_satisfied(self):
        model, _ = _make_pendulum_with_spring()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)

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

    @casadi_required
    def test_torque_runs(self):
        model, _ = _make_pendulum_with_torque()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)
        assert len(T) > 1

    @casadi_required
    def test_torque_vs_lsoda(self):
        """Torque model: CASADI-DAE agrees with LSODA."""
        model1, _ = _make_pendulum_with_torque()
        T1, uT1 = model1.solve(method="LSODA", t_final=0.5, dt=0.05)

        _GroundType._instance._markers = [_GroundType._instance.origin]
        Body.COUNT = 0

        model2, _ = _make_pendulum_with_torque()
        T2, uT2 = model2.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)

        np.testing.assert_allclose(uT1, uT2, atol=1e-3)


# ---------------------------------------------------------------------------
# 10. RotSdaForce
# ---------------------------------------------------------------------------

class TestRotSdaForce:

    @casadi_required
    def test_rotsda_runs(self):
        model, _ = _make_pendulum_with_rotsda()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)
        assert len(T) > 1

    @casadi_required
    def test_rotsda_vs_lsoda(self):
        """RotSdaForce model: CASADI-DAE agrees with LSODA."""
        model1, _ = _make_pendulum_with_rotsda()
        T1, uT1 = model1.solve(method="LSODA", t_final=0.5, dt=0.05)

        _GroundType._instance._markers = [_GroundType._instance.origin]
        Body.COUNT = 0

        model2, _ = _make_pendulum_with_rotsda()
        T2, uT2 = model2.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)

        np.testing.assert_allclose(uT1, uT2, atol=1e-3)


# ---------------------------------------------------------------------------
# 11. No constraints (pure ODE fallback)
# ---------------------------------------------------------------------------

class TestNoConstraints:

    @casadi_required
    def test_free_body_gravity(self):
        """Free body under gravity: y(t) = -g/2 * t^2.  No joints."""
        model, B = _make_free_body()
        T, uT = model.solve(method="CASADI-DAE", t_final=0.5, dt=0.05)

        g = 9.81
        y_expected = -0.5 * g * T**2
        y_actual = uT[:, 1]

        np.testing.assert_allclose(y_actual, y_expected, atol=1e-5,
                                   err_msg="Free-body gravity y(t) mismatch")
