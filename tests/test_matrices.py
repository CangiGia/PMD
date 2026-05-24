"""Static matrix regression tests for PMD validation models.

For each canonical validation model (MSD, DCP, FBL, CS) we verify:

  - ``M_matrix``: global (diagonal) mass matrix has exact expected values.
  - ``Phi``: constraint vector is zero at t=0 for models with consistent ICs.
  - ``Phi_q``: constraint Jacobian at t=0 matches a hardcoded reference to
    6 decimal places.

These tests act as numerical regression guards: any refactoring that
accidentally alters the mathematical structure of the constraint Jacobian
or the mass assembly logic will be detected immediately.

Pattern inspired by petrobras/ross — hardcoded reference arrays +
``assert_almost_equal`` — but at ``decimal=6`` (tighter than ROSS's 5).

Notes
-----
* MSD and DCP have open-chain topologies: Phi == 0 exactly at IC.
* CS has explicit initial positions with negligible (~1e-7) IC residuals
  due to floating-point representation; Phi is tested with atol=1e-5.
* FBL is a closed-loop mechanism — forward-kinematics assembly does NOT
  guarantee loop-closure; Phi is intentionally NOT tested for FBL.
  Only M_matrix and Phi_q are tested.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\PMD
Command:  python -m pytest tests/test_matrices.py -v
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_allclose

from pmd.core.constraints import RevJoint, TranJoint
from pmd.core.forces import Function, PtpForce, Weight
from pmd.core.model import Body, Ground, _GroundType
from pmd.core.motion import RotMotion
from pmd.core.solver import PlanarMultibodyModel


# ---------------------------------------------------------------------------
# Global-state reset (required between tests)
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset():
    """Reset singleton ground and body counter between tests."""
    yield
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]
    Body.COUNT = 0


# ---------------------------------------------------------------------------
# Model fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def msd():
    """Mass–Spring–Damper: 1 body, 1 TranJoint.  nB=1, nC=2, nDOF=1."""
    mass = Body(mass=1000.0, inertia=0.5,
                position=[0.0, 0.0], velocity=[0.0, 10.0], name="mass")
    mk_g_spring = Ground.add_marker([0.0, -2.0])
    mk_g_tran   = Ground.add_marker([0.0,  0.0], theta=np.pi / 2)
    mk_m_spring = mass.add_marker([0.0, -0.5])
    mk_m_tran   = mass.add_marker([0.0,  0.0], theta=np.pi / 2)
    f_sd = PtpForce(iMarker=mk_m_spring, jMarker=mk_g_spring,
                    k=1.0e6, L0=1.5, dc=6500.0)
    j_gs = TranJoint(iMarker=mk_m_tran, jMarker=mk_g_tran)
    model = PlanarMultibodyModel(bodies=[mass], joints=[j_gs], forces=[f_sd])
    model.t = 0.0
    return model


@pytest.fixture
def dcp():
    """Double Compound Pendulum: 2 bodies, 2 RevJoints.  nB=2, nC=4, nDOF=2."""
    _L1 = 250.0e-3
    _L2 = 400.0e-3
    _th1 = np.deg2rad(-12.9946167919)
    _th2 = np.deg2rad( 13.7697807428)

    link_1 = Body(
        mass=1.756261, inertia=0.011854,
        position=[_L1 / 2 * np.cos(_th1), _L1 / 2 * np.sin(_th1)],
        orientation=_th1, name="link_1",
    )
    link_2 = Body(
        mass=2.692381, inertia=0.042171,
        position=[
            _L1 * np.cos(_th1) + _L2 / 2 * np.cos(_th2),
            _L1 * np.sin(_th1) + _L2 / 2 * np.sin(_th2),
        ],
        orientation=_th2, name="link_2",
    )

    g_l1  = Ground.add_marker([0.0, 0.0])
    l1_g  = link_1.add_marker([-_L1 / 2, 0.0])
    l1_l2 = link_1.add_marker([ _L1 / 2, 0.0])
    l2_l1 = link_2.add_marker([-_L2 / 2, 0.0])

    j_gl1  = RevJoint(iMarker=g_l1,  jMarker=l1_g,  q0=_th1,        name="ground_link_1")
    j_l1l2 = RevJoint(iMarker=l1_l2, jMarker=l2_l1, q0=_th2 - _th1, name="link_1_link_2")
    fw = Weight()
    model = PlanarMultibodyModel(bodies=[link_1, link_2],
                                  joints=[j_gl1, j_l1l2], forces=[fw])
    model.t = 0.0
    return model


@pytest.fixture
def fbl():
    """Four-Bar Linkage: 3 bodies, 4 RevJoints + RotMotion.  nB=3, nC=9, nDOF=0.

    Closed-loop: auto-assembly does not guarantee loop-closure;
    Phi is NOT tested for this model.
    """
    _L1 = 40.0e-3
    _L2 = 120.0e-3
    _L3 = 80.0e-3
    _T1 = np.deg2rad( 40.000)
    _T2 = np.deg2rad( 20.298)
    _T3 = np.deg2rad(-122.605)
    _O4 = [0.100, 0.000]
    _OMEGA = np.deg2rad(-25.0)

    # No explicit positions — auto-assembly places bodies via joint q0 values
    link_1 = Body(mass=1.1528090607e-02, inertia=7.4696175803e-05, name="link_1")
    link_2 = Body(mass=3.1498650607e-02, inertia=0.2078332421e-05,  name="link_2")
    link_3 = Body(mass=2.1513370607e-02, inertia=0.141264709e-03,   name="link_3")

    g_l1  = Ground.add_marker([0.000, 0.0])
    g_l3  = Ground.add_marker(_O4)
    l1_g  = link_1.add_marker([-_L1 / 2, 0.0])
    l1_l2 = link_1.add_marker([ _L1 / 2, 0.0])
    l2_l1 = link_2.add_marker([-_L2 / 2, 0.0])
    l2_l3 = link_2.add_marker([ _L2 / 2, 0.0])
    l3_l2 = link_3.add_marker([ _L3 / 2, 0.0])
    l3_g  = link_3.add_marker([-_L3 / 2, 0.0])

    fn_mot = Function(type="a", coeff=[_T1, _OMEGA, 0.0])
    j_gl1  = RevJoint(iMarker=g_l1,  jMarker=l1_g,  q0=_T1,          name="ground_link1")
    j_l1l2 = RevJoint(iMarker=l1_l2, jMarker=l2_l1, q0=_T2 - _T1,    name="link1_link2")
    j_l2l3 = RevJoint(iMarker=l2_l3, jMarker=l3_l2, q0=_T3 - _T2,    name="link2_link3")
    j_l3g  = RevJoint(iMarker=l3_g,  jMarker=g_l3,  q0=np.pi - _T3,  name="link3_ground")
    j_mot  = RotMotion(iBody=link_1, jBody=Ground, iFunct=fn_mot,      name="motor")
    fw = Weight()
    model = PlanarMultibodyModel(
        bodies=[link_1, link_2, link_3],
        joints=[j_gl1, j_l1l2, j_l2l3, j_l3g, j_mot],
        forces=[fw], functions=[fn_mot],
    )
    model.t = 0.0
    return model


@pytest.fixture
def cs():
    """Crank–Slider: 3 bodies, 3 RevJoints + TranJoint + RotMotion.  nB=3, nC=9, nDOF=0."""
    crankshaft = Body(
        mass=0.8169313258, inertia=1.294122161962e-4,
        position=[0.03535533906, 0.03535533906],
        orientation=np.deg2rad(45.0), name="crankshaft",
    )
    rod = Body(
        mass=2.689171326, inertia=4.414522158e-4,
        position=[0.2675608549, 0.03535522707],
        orientation=np.deg2rad(-10.1821), name="rod",
    )
    slider = Body(
        mass=0.2614140191, inertia=4.182624305e-5,
        position=[0.4644110316, 0.0], name="slider",
    )

    mk_g_crank  = Ground.add_marker([0.0, 0.0])
    mk_g_slider = Ground.add_marker([0.4644110316, 0.0], theta=0.0)
    mk_c_ground = crankshaft.add_marker([-0.05, 0.0])
    mk_c_rod    = crankshaft.add_marker([ 0.05, 0.0])
    mk_r_crank  = rod.add_marker([-0.2, 0.0])
    mk_r_slider = rod.add_marker([ 0.2, 0.0])
    mk_s_rod    = slider.add_marker([0.0, 0.0])
    mk_s_ground = slider.add_marker([0.0, 0.0], theta=0.0)

    fn_mot = Function(type="a",
                      coeff=[np.deg2rad(55.1821), np.deg2rad(20.0), 0.0])
    j_gc  = RevJoint(iMarker=mk_g_crank,  jMarker=mk_c_ground, name="ground_crankshaft")
    j_cr  = RevJoint(iMarker=mk_c_rod,    jMarker=mk_r_crank,  name="crankshaft_rod")
    j_rs  = RevJoint(iMarker=mk_s_rod,    jMarker=mk_r_slider, name="rod_slider")
    j_gs  = TranJoint(iMarker=mk_s_ground, jMarker=mk_g_slider, name="ground_slider")
    j_mot = RotMotion(iBody=crankshaft, jBody=rod, iFunct=fn_mot, name="crank_rod_motor")
    fw = Weight(gravity=9.80665)
    model = PlanarMultibodyModel(
        bodies=[crankshaft, rod, slider],
        joints=[j_gc, j_cr, j_rs, j_gs, j_mot],
        forces=[fw], functions=[fn_mot],
    )
    model.t = 0.0
    return model


# ===========================================================================
# MSD — Mass–Spring–Damper
# ===========================================================================

class TestMSD:

    def test_dof(self, msd):
        assert msd.nB == 1
        assert msd.nC == 2
        assert msd.nDOF == 1

    def test_mass_matrix(self, msd):
        # fmt: off
        M_ref = np.array([
            [1000.000000,    0.000000,    0.000000],
            [   0.000000, 1000.000000,    0.000000],
            [   0.000000,    0.000000,    0.500000],
        ])
        # fmt: on
        assert_almost_equal(msd.M_matrix, M_ref, decimal=6)

    def test_phi_at_ic(self, msd):
        """Open-chain model: constraints satisfied exactly at IC."""
        phi = msd._compute_constraints().ravel()
        assert_allclose(phi, np.zeros(msd.nC), atol=1e-12)

    def test_phi_q_at_ic(self, msd):
        # fmt: off
        PhiQ_ref = np.array([
            [-1.000000,  0.000000,  0.000000],
            [ 0.000000,  0.000000,  1.000000],
        ])
        # fmt: on
        assert_almost_equal(msd._compute_jacobian(), PhiQ_ref, decimal=6)


# ===========================================================================
# DCP — Double Compound Pendulum
# ===========================================================================

class TestDCP:

    def test_dof(self, dcp):
        assert dcp.nB == 2
        assert dcp.nC == 4
        assert dcp.nDOF == 2

    def test_mass_matrix(self, dcp):
        # fmt: off
        M_ref = np.array([
            [1.756261, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
            [0.000000, 1.756261, 0.000000, 0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 0.011854, 0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 0.000000, 2.692381, 0.000000, 0.000000],
            [0.000000, 0.000000, 0.000000, 0.000000, 2.692381, 0.000000],
            [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.042171],
        ])
        # fmt: on
        assert_almost_equal(dcp.M_matrix, M_ref, decimal=6)

    def test_phi_at_ic(self, dcp):
        """Open-chain model: constraints satisfied exactly at IC."""
        phi = dcp._compute_constraints().ravel()
        assert_allclose(phi, np.zeros(dcp.nC), atol=1e-12)

    def test_phi_q_at_ic(self, dcp):
        # fmt: off
        PhiQ_ref = np.array([
            [-1.000000,  0.000000,  0.028107, 0.000000,  0.000000,  0.000000],
            [ 0.000000, -1.000000,  0.121799, 0.000000,  0.000000,  0.000000],
            [ 1.000000,  0.000000,  0.028107,-1.000000,  0.000000, -0.047604],
            [ 0.000000,  1.000000,  0.121799, 0.000000, -1.000000,  0.194252],
        ])
        # fmt: on
        assert_almost_equal(dcp._compute_jacobian(), PhiQ_ref, decimal=6)


# ===========================================================================
# FBL — Four-Bar Linkage  (closed-loop: Phi test skipped)
# ===========================================================================

class TestFBL:

    def test_dof(self, fbl):
        assert fbl.nB == 3
        assert fbl.nC == 9
        assert fbl.nDOF == 0

    def test_mass_matrix(self, fbl):
        # fmt: off
        M_ref = np.diag([
            0.011528090607,
            0.011528090607,
            0.000074696176,
            0.031498650607,
            0.031498650607,
            0.000002078332,
            0.021513370607,
            0.021513370607,
            0.000141264709,
        ])
        # fmt: on
        assert_almost_equal(fbl.M_matrix, M_ref, decimal=6)

    def test_phi_q_at_ic(self, fbl):
        # fmt: off
        PhiQ_ref = np.array([
            [-1.000000,  0.000000, -0.012856,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000],
            [ 0.000000, -1.000000,  0.015321,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000],
            [ 1.000000,  0.000000, -0.012856, -1.000000,  0.000000, -0.020601,  0.000000,  0.000000,  0.000000],
            [ 0.000000,  1.000000,  0.015321,  0.000000, -1.000000,  0.056352,  0.000000,  0.000000,  0.000000],
            [ 0.000000,  0.000000,  0.000000,  1.000000,  0.000000, -0.020601, -1.000000,  0.000000, -0.033665],
            [ 0.000000,  0.000000,  0.000000,  0.000000,  1.000000,  0.056352,  0.000000, -1.000000, -0.021603],
            [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  1.000000,  0.000000, -0.033665],
            [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  1.000000, -0.021603],
            [ 0.000000,  0.000000,  1.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000],
        ])
        # fmt: on
        assert_almost_equal(fbl._compute_jacobian(), PhiQ_ref, decimal=6)


# ===========================================================================
# CS — Crank–Slider
# ===========================================================================

class TestCS:

    def test_dof(self, cs):
        assert cs.nB == 3
        assert cs.nC == 9
        assert cs.nDOF == 0

    def test_mass_matrix(self, cs):
        # fmt: off
        M_ref = np.diag([
            0.816931325800,
            0.816931325800,
            0.000129412216,
            2.689171326000,
            2.689171326000,
            0.000441452216,
            0.261414019100,
            0.261414019100,
            0.000041826243,
        ])
        # fmt: on
        assert_almost_equal(cs.M_matrix, M_ref, decimal=6)

    def test_phi_at_ic(self, cs):
        """Explicit ICs for CS have negligible floating-point residuals (~1e-7)."""
        phi = cs._compute_constraints().ravel()
        assert_allclose(phi, np.zeros(cs.nC), atol=1e-5)

    def test_phi_q_at_ic(self, cs):
        # fmt: off
        PhiQ_ref = np.array([
            [-1.000000,  0.000000, -0.035355,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000],
            [ 0.000000, -1.000000,  0.035355,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000],
            [ 1.000000,  0.000000, -0.035355, -1.000000,  0.000000,  0.035355,  0.000000,  0.000000,  0.000000],
            [ 0.000000,  1.000000,  0.035355,  0.000000, -1.000000,  0.196850,  0.000000,  0.000000,  0.000000],
            [ 0.000000,  0.000000,  0.000000, -1.000000,  0.000000, -0.035355,  1.000000,  0.000000,  0.000000],
            [ 0.000000,  0.000000,  0.000000,  0.000000, -1.000000, -0.196850,  0.000000,  1.000000,  0.000000],
            [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  1.000000,  0.000000],
            [ 0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  1.000000],
            [ 0.000000,  0.000000,  1.000000,  0.000000,  0.000000, -1.000000,  0.000000,  0.000000,  0.000000],
        ])
        # fmt: on
        assert_almost_equal(cs._compute_jacobian(), PhiQ_ref, decimal=6)
