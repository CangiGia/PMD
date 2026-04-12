import time
import numpy as np
import scipy as sc
from PMD.src import *
from PMD.src.model import _GroundType
import matplotlib.pyplot as plt
from PMD.examples._plot_utils import plot_comparison
from PMD.gui.app import PostProcessor

#%% bodies
B1 = Body(mass=2, inertia=0.5, position=[0.4398, 0.2512], orientation=-0.0367)   # lower suspension arm
B2 = Body(mass=30, inertia=2.5, position=[0.6817, 0.3498], orientation=0.0783)  # wheel assembly
B3 = Body(mass=1, inertia=0.5, position=[0.4463, 0.4308], orientation=6.5222)   # upper suspension arm

#%% markers
q1 = B1.add_marker([-0.24, 0.0])       # Q point - lower suspension arm side
a1 = B1.add_marker([0.18, 0.0])        # A point - lower suspension arm side
a2 = B2.add_marker([-0.07, -0.10])     # A point - wheel assembly side
b2 = B2.add_marker([-0.10, 0.12])      # B point - wheel assembly side
b3 = B3.add_marker([0.13, 0.0])        # B point - upper suspension arm side
o3 = B3.add_marker([-0.13, 0.0])       # O point - upper suspension arm side
o0 = Ground.add_marker([0.32, 0.40])   # O point - ground
q0 = Ground.add_marker([0.20, 0.26])   # Q point - ground
e1 = B1.add_marker([0.0, 0.0])         # E point - Lower suspension arm side
f0 = Ground.add_marker([0.38, 0.43])   # F point - Ground

#%% joints
j1 = RevJoint(iMarker=q1, jMarker=q0)  # Revolute joint in Q
j2 = RevJoint(iMarker=a1, jMarker=a2)  # Revolute joint in A
j3 = RevJoint(iMarker=b2, jMarker=b3)  # Revolute joint in B
j4 = RevJoint(iMarker=o3, jMarker=o0)  # Revolute joint in O

#%% forces
s1 = PtpForce(iMarker=e1, jMarker=f0, k=90000, L0=0.23, dc=1100)  # spring-damper force

# custom wheel contact force
_k_wh = 50000
_L0_wh = 0.35
_dc_wh = 1000

def wheel_contact():
    """Custom force used to define the wheel contact condition."""
    dely = B2.position[1] - _L0_wh
    if dely < 0:
        fy = (_k_wh * dely + _dc_wh * B2.velocity[1]).item()
        return [{'body': B2, 'force': [0, -fy], 'torque': 0}]
    return []

s2 = UserForce(callback=wheel_contact)

s3 = Weight()  # gravity force

# ── Benchmark helpers ──────────────────────────────────────────────────────

def _make_model():
    """Rebuild the quarter-car from scratch (resets Ground singleton + Body counter)."""
    _GroundType._instance._markers = [_GroundType._instance.origin]
    Body.COUNT = 0
    _B1 = Body(mass=2,  inertia=0.5, position=[0.4398, 0.2512], orientation=-0.0367)
    _B2 = Body(mass=30, inertia=2.5, position=[0.6817, 0.3498], orientation=0.0783)
    _B3 = Body(mass=1,  inertia=0.5, position=[0.4463, 0.4308], orientation=6.5222)
    _q1 = _B1.add_marker([-0.24,  0.0]);  _a1 = _B1.add_marker([ 0.18,  0.0])
    _a2 = _B2.add_marker([-0.07, -0.10]); _b2 = _B2.add_marker([-0.10,  0.12])
    _b3 = _B3.add_marker([ 0.13,  0.0]); _o3 = _B3.add_marker([-0.13,  0.0])
    _o0 = Ground.add_marker([0.32, 0.40]); _q0 = Ground.add_marker([0.20, 0.26])
    _e1 = _B1.add_marker([0.0, 0.0]);     _f0 = Ground.add_marker([0.38, 0.43])
    _j1 = RevJoint(iMarker=_q1, jMarker=_q0)
    _j2 = RevJoint(iMarker=_a1, jMarker=_a2)
    _j3 = RevJoint(iMarker=_b2, jMarker=_b3)
    _j4 = RevJoint(iMarker=_o3, jMarker=_o0)
    _s1 = PtpForce(iMarker=_e1, jMarker=_f0, k=90000, L0=0.23, dc=1100)
    def _wc():
        dely = _B2.position[1] - 0.35
        if dely < 0:
            fy = (50000 * dely + 1000 * _B2.velocity[1]).item()
            return [{'body': _B2, 'force': [0, -fy], 'torque': 0}]
        return []
    _s2 = UserForce(callback=_wc)
    _s3 = Weight()
    return PlanarMultibodyModel(bodies=[_B1, _B2, _B3],
                                joints=[_j1, _j2, _j3, _j4],
                                forces=[_s1, _s2, _s3])


def _constraint_violation(model, T, uT):
    """Return array of max|Phi(q_k)| at each output time step."""
    nB3 = 3 * model.nB
    phi_max = np.zeros(len(T))
    for k in range(len(T)):
        q_k = uT[k, :nB3]
        for Bi, body in enumerate(model.Bodies):
            ir = 3 * Bi
            body.position    = q_k[ir:ir + 2].reshape(2, 1)
            body.orientation = float(q_k[ir + 2])
        model.t = float(T[k])
        model._update_position()
        phi_max[k] = float(np.max(np.abs(model._compute_constraints())))
    return phi_max


#%% solution
quarter_car = PlanarMultibodyModel(
    bodies=[B1, B2, B3],
    joints=[j1, j2, j3, j4],
    forces=[s1, s2, s3])
T, uT = quarter_car.solve(method='CASADI-DAE', t_final=10.0, t_eval=np.linspace(0, 10, 10001),
                          ic_correct=True)

# if __name__ == '__main__':
#     plot_comparison(T, uT, matlab_filename='AA.txt', model_title='AA')

post_proc = PostProcessor(model=quarter_car, T=T, uT=uT)
post_proc.show()


#%% ─────────────────────────────────────────────────────────────────────────
#   BENCHMARK: CasADi-DAE vs ODE solvers
#   run directly:  python _test_AA.py
#   (starts after the main PostProcessor window is closed)
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':

    BENCH_METHODS = [
        ('CASADI-DAE', dict(method='CASADI-DAE', ic_correct=True)),
        ('Radau',      dict(method='Radau',       ic_correct=True)),
        ('BDF',        dict(method='BDF',         ic_correct=True)),
        ('LSODA',      dict(method='LSODA',       ic_correct=True)),
    ]
    T_FINAL     = 5.0
    BENCH_TEVAL = np.linspace(0, T_FINAL, 501)

    print('\n' + '=' * 62)
    print('  BENCHMARK  —  quarter-car model  (_test_AA)')
    print('=' * 62)
    print(f'  t_final = {T_FINAL} s   |   output steps = {len(BENCH_TEVAL)}\n')

    bench = {}
    for label, kwargs in BENCH_METHODS:
        m = _make_model()
        t0 = time.perf_counter()
        _T, _uT = m.solve(t_eval=BENCH_TEVAL, **kwargs)
        elapsed = time.perf_counter() - t0
        phi = _constraint_violation(m, _T, _uT)
        bench[label] = {'T': _T, 'uT': _uT, 'time': elapsed, 'phi': phi}
        print(f'  {label:12s}  CPU: {elapsed:7.2f} s   max|Phi| = {phi.max():.2e}')

    print('=' * 62)

    # ── Figure ──────────────────────────────────────────────────────────
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    fig_b, (ax_t, ax_phi) = plt.subplots(1, 2, figsize=(13, 5))
    fig_b.suptitle(
        f'Quarter-car benchmark  (t_final = {T_FINAL} s,  {len(BENCH_TEVAL)} steps)',
        fontsize=13)

    # Left panel: wall-clock time
    labels = list(bench.keys())
    times  = [bench[k]['time'] for k in labels]
    bars   = ax_t.bar(labels, times, color=colors, edgecolor='k', width=0.5)
    ax_t.bar_label(bars, fmt='%.2f s', padding=4, fontsize=9)
    ax_t.set_ylabel('Wall-clock time [s]')
    ax_t.set_title('Simulation time')
    ax_t.set_ylim(0, max(times) * 1.35)

    # Right panel: constraint violation history
    for (label, res), c in zip(bench.items(), colors):
        ax_phi.semilogy(res['T'], np.maximum(res['phi'], 1e-18),
                        label=label, color=c, lw=1.5)
    ax_phi.set_xlabel('Time [s]')
    ax_phi.set_ylabel('max |Φ(q)| [-]')
    ax_phi.set_title('Constraint violation history')
    ax_phi.legend()
    ax_phi.grid(True, which='both', ls='--', alpha=0.4)

    plt.tight_layout()
    out_path = 'benchmark_AA.png'
    # plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\n  Plot saved → {out_path}')
    plt.show()
