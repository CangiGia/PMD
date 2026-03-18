"""
Planar Multi-Body Dynamics Simulation Solver

This module provides algorithms and numerical methods for solving 
planar multi-body dynamic models, including rigid bodies, constraints, 
and external forces. It is designed for academic research and engineering 
applications.

INDEX CONVENTION:
================
Body indices:  0 = ground (fixed), 1..nB = moving bodies
Point indices: 0-based (Python standard)
Joint indices: 0-based (Python standard)

Internal state vector u:
  u[0:nB3]      = positions  [r1x, r1y, p1, r2x, r2y, p2, ...]
  u[nB3:2*nB3]  = velocities [dr1x, dr1y, dp1, ...]

Internal index attributes (0-based):
  body._irc   = 3*Bi       (start index for position in u, inclusive)
  body._irv   = nB3 + 3*Bi (start index for velocity in u, inclusive)
  joint._rows = constraint row start (0-based, inclusive)
  joint._rowe = constraint row end (0-based, exclusive)
  joint._colis/_coljs = Jacobian column start (0-based, inclusive)
  joint._colie/_colje = Jacobian column end (0-based, exclusive)

All slicing uses standard Python [start:end) convention directly.

Author: Giacomo Cangi
"""


import os
import numpy as np
import scipy as sc
import numpy.linalg as lng
from .utils import *
from .mechanics import *
from .model import Ground
from scipy.integrate import solve_ivp
from tqdm import tqdm


class SolResult:
    """Simulation result container with lazy post-processing and plot methods.

    Supports both tuple unpacking (``T, uT = result``) and attribute access
    (``result.t``, ``result.y`` for a scipy-like interface).

    Attributes
    ----------
    t : ndarray, shape (n,)
        Time points.
    y : ndarray, shape (2*nB3, n)
        State matrix (positions + velocities), scipy ``solve_ivp`` convention.
    uT : ndarray, shape (n, 2*nB3)
        State matrix, legacy convention (rows = time steps).
    nB : int
        Number of moving bodies.
    """

    def __init__(self, t, uT, model=None, dense_sol=None):
        self.t = t
        self.y = uT.T   # scipy-like: shape (2*nB3, n_timesteps)
        self.uT = uT    # legacy: shape (n_timesteps, 2*nB3)
        self._model = model
        self._dense_sol = dense_sol
        self._accelerations = None
        self._reactions = None

        if model is not None:
            self.nB = model.nB
            self._nB3 = 3 * model.nB
            self._body_names = [f'Body {i+1}' for i in range(model.nB)]
        else:
            self.nB = 0
            self._nB3 = 0
            self._body_names = []

    def __iter__(self):
        return iter((self.t, self.uT))

    def __repr__(self):
        return f"SolResult(t: {self.t.shape}, y: {self.y.shape})"

    # ── Lazy post-processing properties ────────────────────────────

    @property
    def accelerations(self):
        """Generalized accelerations, shape (nSteps, nB3). Computed lazily on first access."""
        if self._accelerations is None:
            if self._model is None:
                raise RuntimeError("No model reference available for post-processing.")
            self._accelerations, self._reactions = self._model._post_process(self.t, self.uT)
        return self._accelerations

    @property
    def reactions(self):
        """Lagrange multipliers, shape (nSteps, nConstraints). Computed lazily on first access."""
        if self._reactions is None:
            if self._model is None:
                raise RuntimeError("No model reference available for post-processing.")
            self._accelerations, self._reactions = self._model._post_process(self.t, self.uT)
        return self._reactions

    # ── Data access helpers ────────────────────────────────────────

    def get_body_states(self, body_index):
        """Extract x, y, phi, dx, dy, dphi for a specific body.

        Parameters
        ----------
        body_index : int
            1-based body index.

        Returns
        -------
        dict
            Keys: 'x', 'y', 'phi', 'dx', 'dy', 'dphi', each ndarray shape (nSteps,).
        """
        i = body_index - 1  # 0-based
        nB3 = self._nB3
        return {
            'x':    self.uT[:, 3*i],
            'y':    self.uT[:, 3*i + 1],
            'phi':  self.uT[:, 3*i + 2],
            'dx':   self.uT[:, nB3 + 3*i],
            'dy':   self.uT[:, nB3 + 3*i + 1],
            'dphi': self.uT[:, nB3 + 3*i + 2],
        }

    def get_body_accelerations(self, body_index):
        """Extract ddx, ddy, ddphi for a specific body.

        Parameters
        ----------
        body_index : int
            1-based body index.

        Returns
        -------
        dict
            Keys: 'ddx', 'ddy', 'ddphi', each ndarray shape (nSteps,).
        """
        i = body_index - 1
        acc = self.accelerations  # trigger lazy compute
        return {
            'ddx':   acc[:, 3*i],
            'ddy':   acc[:, 3*i + 1],
            'ddphi': acc[:, 3*i + 2],
        }

    # ── Dense output resampling ────────────────────────────────────

    def resample(self, t_start=None, t_end=None, dt=None):
        """Create a new SolResult on a different time grid using dense output.

        Parameters
        ----------
        t_start : float, optional
            Start time (default: self.t[0]).
        t_end : float, optional
            End time (default: self.t[-1]).
        dt : float, optional
            Time step (default: same as original).

        Returns
        -------
        SolResult
            New result object on the resampled grid.

        Raises
        ------
        RuntimeError
            If dense_output was not enabled during solve().
        """
        if self._dense_sol is None:
            raise RuntimeError(
                "Dense output not available. Call solve() with dense_output=True "
                "to enable resampling."
            )
        if t_start is None:
            t_start = self.t[0]
        if t_end is None:
            t_end = self.t[-1]
        if dt is None:
            dt = self.t[1] - self.t[0]

        T_new = np.arange(t_start, t_end + dt / 2, dt)
        uT_new = self._dense_sol(T_new).T
        return SolResult(T_new, uT_new, self._model, dense_sol=self._dense_sol)

    # ── Plot methods ───────────────────────────────────────────────

    def plot_displacements(self, bodies=None, figsize=(12, 4)):
        """Plot x, y, φ for selected bodies.

        Parameters
        ----------
        bodies : list of int or None
            1-based body indices. None = all bodies.
        figsize : tuple
            Figure size per body row.
        """
        import matplotlib.pyplot as plt

        if bodies is None:
            bodies = list(range(1, self.nB + 1))

        fig, axes = plt.subplots(len(bodies), 3,
                                 figsize=(figsize[0], figsize[1] * len(bodies)),
                                 squeeze=False)

        labels = ['x [m]', 'y [m]', '\u03c6 [rad]']
        keys = ['x', 'y', 'phi']

        for row, bi in enumerate(bodies):
            states = self.get_body_states(bi)
            for col, (key, label) in enumerate(zip(keys, labels)):
                axes[row, col].plot(self.t, states[key], 'b-', linewidth=1.0)
                axes[row, col].set_xlabel('t [s]')
                axes[row, col].set_ylabel(label)
                axes[row, col].set_title(f'{self._body_names[bi-1]} \u2014 {key}')
                axes[row, col].grid(True, alpha=0.3)

        fig.suptitle('Displacements', fontsize=14, fontweight='bold')
        fig.tight_layout()
        plt.show()

    def plot_velocities(self, bodies=None, figsize=(12, 4)):
        """Plot dx, dy, dφ for selected bodies.

        Parameters
        ----------
        bodies : list of int or None
            1-based body indices. None = all bodies.
        figsize : tuple
            Figure size per body row.
        """
        import matplotlib.pyplot as plt

        if bodies is None:
            bodies = list(range(1, self.nB + 1))

        fig, axes = plt.subplots(len(bodies), 3,
                                 figsize=(figsize[0], figsize[1] * len(bodies)),
                                 squeeze=False)

        labels = ['\u1e8b [m/s]', '\u1e8f [m/s]', '\u03c6\u0307 [rad/s]']
        keys = ['dx', 'dy', 'dphi']

        for row, bi in enumerate(bodies):
            states = self.get_body_states(bi)
            for col, (key, label) in enumerate(zip(keys, labels)):
                axes[row, col].plot(self.t, states[key], 'b-', linewidth=1.0)
                axes[row, col].set_xlabel('t [s]')
                axes[row, col].set_ylabel(label)
                axes[row, col].set_title(f'{self._body_names[bi-1]} \u2014 {key}')
                axes[row, col].grid(True, alpha=0.3)

        fig.suptitle('Velocities', fontsize=14, fontweight='bold')
        fig.tight_layout()
        plt.show()

    def plot_accelerations(self, bodies=None, figsize=(12, 4)):
        """Plot ddx, ddy, ddφ for selected bodies. Triggers lazy computation.

        Parameters
        ----------
        bodies : list of int or None
            1-based body indices. None = all bodies.
        figsize : tuple
            Figure size per body row.
        """
        import matplotlib.pyplot as plt

        if bodies is None:
            bodies = list(range(1, self.nB + 1))

        fig, axes = plt.subplots(len(bodies), 3,
                                 figsize=(figsize[0], figsize[1] * len(bodies)),
                                 squeeze=False)

        labels = ['\u1e8d [m/s\u00b2]', '\u00ff [m/s\u00b2]', '\u03c6\u0308 [rad/s\u00b2]']
        keys = ['ddx', 'ddy', 'ddphi']

        for row, bi in enumerate(bodies):
            acc = self.get_body_accelerations(bi)
            for col, (key, label) in enumerate(zip(keys, labels)):
                axes[row, col].plot(self.t, acc[key], 'b-', linewidth=1.0)
                axes[row, col].set_xlabel('t [s]')
                axes[row, col].set_ylabel(label)
                axes[row, col].set_title(f'{self._body_names[bi-1]} \u2014 {key}')
                axes[row, col].grid(True, alpha=0.3)

        fig.suptitle('Accelerations', fontsize=14, fontweight='bold')
        fig.tight_layout()
        plt.show()

    def plot_reactions(self, joints=None, figsize=(12, 3)):
        """Plot Lagrange multipliers for selected constraint rows.

        Parameters
        ----------
        joints : list of int or None
            0-based constraint row indices. None = all.
        figsize : tuple
            Figure size per row.
        """
        import matplotlib.pyplot as plt

        reactions = self.reactions  # trigger lazy compute
        nConst = reactions.shape[1]

        if joints is None:
            indices = list(range(nConst))
        else:
            indices = [j for j in joints if j < nConst]

        nPlots = len(indices)
        fig, axes = plt.subplots(nPlots, 1,
                                 figsize=(figsize[0], figsize[1] * nPlots),
                                 squeeze=False)
        for row, k in enumerate(indices):
            axes[row, 0].plot(self.t, reactions[:, k], 'r-', linewidth=1.0)
            axes[row, 0].set_xlabel('t [s]')
            axes[row, 0].set_ylabel(f'\u03bb_{k+1} [N or N\u00b7m]')
            axes[row, 0].set_title(f'Constraint reaction {k+1}')
            axes[row, 0].grid(True, alpha=0.3)

        fig.suptitle('Constraint Reactions', fontsize=14, fontweight='bold')
        fig.tight_layout()
        plt.show()

    def plot_phase(self, body_index, dof='x', figsize=(6, 6)):
        """Phase portrait (q vs dq) for a specific DOF.

        Parameters
        ----------
        body_index : int
            1-based body index.
        dof : str
            'x', 'y', or 'phi'.
        figsize : tuple
            Figure size.
        """
        import matplotlib.pyplot as plt

        states = self.get_body_states(body_index)
        q_key = dof
        dq_key = 'd' + dof

        q_labels = {'x': 'x [m]', 'y': 'y [m]', 'phi': '\u03c6 [rad]'}
        dq_labels = {'x': '\u1e8b [m/s]', 'y': '\u1e8f [m/s]', 'phi': '\u03c6\u0307 [rad/s]'}

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.plot(states[q_key], states[dq_key], 'b-', linewidth=0.8)
        ax.set_xlabel(q_labels[dof])
        ax.set_ylabel(dq_labels[dof])
        ax.set_title(f'{self._body_names[body_index-1]} \u2014 Phase portrait ({dof})')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='datalim')
        fig.tight_layout()
        plt.show()

    def plot_energy(self, figsize=(10, 4)):
        """Plot total kinetic energy vs time.

        Parameters
        ----------
        figsize : tuple
            Figure size.
        """
        import matplotlib.pyplot as plt

        nB3 = self._nB3
        dq = self.uT[:, nB3:]  # velocities, shape (nSteps, nB3)

        # Build mass vector [m1, m1, J1, m2, m2, J2, ...]
        mass_vec = []
        for b in self._model.Bodies:
            mass_vec.extend([b.m, b.m, b.J])
        mass_vec = np.array(mass_vec)

        # KE = 0.5 * sum(m_i * dq_i^2)
        KE = 0.5 * np.sum(mass_vec * dq**2, axis=1)

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.plot(self.t, KE, 'g-', linewidth=1.0)
        ax.set_xlabel('t [s]')
        ax.set_ylabel('T [J]')
        ax.set_title('Total Kinetic Energy')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        plt.show()


class PlanarMultibodyModel:
    def __init__(self, bodies, joints=None, forces=None, functions=None,
                 verbose=False):
        self.verbose = verbose
        self.Bodies = list(bodies)
        self.Joints = list(joints) if joints else []
        self.Forces = list(forces) if forces else []
        self.Functs = list(functions) if functions else []

        # Auto-assembly trigger: only if a body has neither r nor p provided
        from .builder import _assemble
        needs_assembly = any(not b._r_given and not b._p_given for b in self.Bodies)
        if needs_assembly:
            _assemble(self.Bodies, self.Joints)

        # initialize the model for simulation automatically
        self._initialize()

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def nB(self):
        """Number of moving bodies."""
        return len(self.Bodies)

    @property
    def nC(self):
        """Total number of constraint equations."""
        return self.Joints[-1]._rowe if self.Joints else 0

    def _initialize(self):
        """
        Initialize the multi-body model considering the values defined 
        by the user.
        """
        # initialize variables
        nB = len(self.Bodies)
        nB3 = 3 * nB

        #// bodies — assign internal index and compute derived quantities
        for Bi, body in enumerate(self.Bodies):
            body._bidx = Bi + 1          # 1-based: ground=0, bodies=1..nB
            body._irc = 3 * Bi
            body._irv = nB3 + 3 * Bi
            body._invm = 1 / body.m
            body._invJ = 1 / body.J
            body._A = A_matrix(body.p)

        # mass (inertia) array and pre-computed diagonal matrix
        self.M_array = np.zeros(nB3)
        self.invM_array = np.zeros(nB3)
        for Bi, body in enumerate(self.Bodies):
            is_ = 3 * Bi
            ie_ = is_ + 3
            self.M_array[is_:ie_] = np.array([body.m, body.m, body.J])
            self.invM_array[is_:ie_] = np.array([body._invm, body._invm, body._invJ])
        self.M_matrix = np.diag(self.M_array)

        #// markers — Ground markers
        for marker in Ground._markers:
            pos_col = marker.position.reshape(2, 1)
            marker._sP  = pos_col.copy()
            marker._sPr = s_rot(pos_col)
            marker._rP  = pos_col.copy()
            if marker.has_orientation:
                marker._u  = marker._ulocal.copy()
                marker._ur = s_rot(marker._ulocal)

        #// markers — Body markers
        for body in self.Bodies:
            for marker in body._markers:
                pos_col = marker.position.reshape(2, 1)
                marker._sP  = body._A @ pos_col
                marker._sPr = s_rot(marker._sP)
                marker._rP  = body.r + marker._sP
                if marker.has_orientation:
                    marker._u  = body._A @ marker._ulocal
                    marker._ur = s_rot(marker._u)

        #// force elements
        for force in self.Forces:
            if force.type == 'weight':
                ug = force._gravity * force._wgt
                for body in self.Bodies:
                    body._wgt = body.m * ug
            elif force.type == 'ptp':
                # derive body references from marker attachments
                force.iBody = force.iMarker.body
                force.jBody = force.jMarker.body

        #// joints
        for joint in self.Joints:
            match joint.type:
                case 'rev':
                    joint._mrows = 2
                    joint._nbody = 2
                    joint.iBody = joint.iMarker.body
                    joint.jBody = joint.jMarker.body
                    if joint.fix == 1:
                        joint._mrows = 3
                        Bi = joint.iBody
                        Bj = joint.jBody
                        if Bi is Ground:
                            joint._p0 = -Bj.p
                        elif Bj is Ground:
                            joint._p0 = Bi.p
                        else:
                            joint._p0 = Bi.p - Bj.p

                case 'tran':
                    joint._mrows = 2
                    joint._nbody = 2
                    joint.iBody = joint.iMarker.body
                    joint.jBody = joint.jMarker.body
                    if joint.fix == 1:
                        joint._mrows = 3
                        Bi = joint.iBody
                        Bj = joint.jBody
                        iPt = joint.iMarker
                        jPt = joint.jMarker
                        if Bi is Ground:
                            joint._p0 = np.linalg.norm(iPt._rP -
                                                    Bj.r - Bj._A @
                                                    jPt.position.reshape(2,1))
                        elif Bj is Ground:
                            joint._p0 = np.linalg.norm(Bi.r +
                                                    Bi._A @
                                                    iPt.position.reshape(2,1) -
                                                    jPt._rP)
                        else:
                            joint._p0 = np.linalg.norm(Bi.r +
                                                    Bi._A @
                                                    iPt.position.reshape(2,1) -
                                                    Bj.r - Bj._A @
                                                    jPt.position.reshape(2,1))

                case 'rev-rev':
                    joint._mrows = 1
                    joint._nbody = 2
                    joint.iBody = joint.iMarker.body
                    joint.jBody = joint.jMarker.body

                case 'rev-tran':
                    joint._mrows = 1
                    joint._nbody = 2
                    joint.iBody = joint.iMarker.body
                    joint.jBody = joint.jMarker.body

                case 'rel-rot' | 'rel-tran':
                    joint._mrows = 1
                    joint._nbody = 1
                case 'disc':
                    joint._mrows = 2
                    joint._nbody = 1

                case 'rigid':
                    joint._mrows = 3
                    joint._nbody = 2
                    Bi = joint.iBody
                    Bj = joint.jBody
                    if Bi is Ground:
                        joint.d0 = -Bj._A.T @ Bj.r
                        joint._p0 = -Bj.p
                    elif Bj is Ground:
                        joint.d0 = Bi.r
                        joint._p0 = Bi.p
                    else:
                        joint.d0 = Bj._A.T @ (Bi.r - Bj.r)
                        joint._p0 = Bi.p - Bj.p

                case _:
                    raise ValueError("Joint type doesn't supported!")

        # Validation V7: check all joint bodies are in self.Bodies or are Ground
        body_ids = {id(b) for b in self.Bodies}
        body_ids.add(id(Ground))
        for joint in self.Joints:
            if id(joint.iBody) not in body_ids:
                raise ValueError(
                    f"Joint iBody {joint.iBody} is not in the model's bodies list")
            if id(joint.jBody) not in body_ids:
                raise ValueError(
                    f"Joint jBody {joint.jBody} is not in the model's bodies list")

        #// functions
        if self.Functs:
            nFc = len(self.Functs)
            for Ci in range(nFc):
                functData(Ci, self.Functs)

        # compute number of constraints and determine row/column pointers
        nConst = 0
        for joint in self.Joints:
            joint._rows = nConst
            joint._rowe = nConst + joint._mrows
            nConst = joint._rowe
            Bi_idx = joint.iBody._bidx
            if Bi_idx != 0:
                joint._colis = 3 * (Bi_idx - 1)
                joint._colie = 3 * Bi_idx
            Bj_idx = joint.jBody._bidx
            if Bj_idx != 0:
                joint._coljs = 3 * (Bj_idx - 1)
                joint._colje = 3 * Bj_idx

        if self.verbose:
            print("\n")
            print("\t... model has been created and initialized correctly ...")
            print("\n")
            
            print("\t... values after initializzation ...")
            print(f"-----")
            print(f"bodies ")
            print(f"-----")
            for i, body in enumerate(self.Bodies, start=1):
                print(f"\t... body: {i}")
                print(f"\t... mass: {body.m}")
                print(f"\t... moment of inertia: {body.J}")
                print(f"\t... position: {', '.join(map(str, body.r.flatten()))}")
                print(f"\t... orientation: {body.p}")
                print(f"\t... velocity: {', '.join(map(str, body.dr.flatten()))}")
                print(f"\t... angular velocity: {body.dp}")
                print(f"\t... acceleration: {', '.join(map(str, body.ddr.flatten()))}")
                print(f"\t... angular acceleration: {body.ddp}")
                print(f"\t... rotational matrix: {', '.join(map(str, body._A.flatten()))}")
                print(f"\t... inverse mass: {body._invm}")
                print(f"\t... inverse moment of inertia: {body._invJ}")
                print(f"\t... weight: {', '.join(map(str, body._wgt.flatten()))}")
                print(f"\t... force: {', '.join(map(str, body._f.flatten()))}")
                print(f"\t... torque: {body._n}")
                print(f"\t... markers: {body._markers}")
                print(f"\t... irc: {body._irc}")
                print(f"\t... irv: {body._irv}")
                print(f"\t... ira: {body._ira}")
                print("\n")
            
            print(f"-----")
            print(f"markers ")
            print(f"-----")
            for i, body in enumerate(self.Bodies, start=1):
                for j, marker in enumerate(body._markers):
                    print(f"\t... body {i}, marker {j}: {marker.name}")
                    print(f"\t... position (local): {marker.position}")
                    print(f"\t... theta: {marker.theta}")
                    print(f"\t... _rP: {', '.join(map(str, marker._rP.flatten()))}")
                    print(f"\t... _sP: {', '.join(map(str, marker._sP.flatten()))}")
                    print(f"\t... _sPr: {', '.join(map(str, marker._sPr.flatten()))}")
                    if marker.has_orientation:
                        print(f"\t... _u: {', '.join(map(str, marker._u.flatten()))}")
                        print(f"\t... _ur: {', '.join(map(str, marker._ur.flatten()))}")
                    print("\n")
            
            print(f"-----")
            print(f"forces ")
            print(f"-----")
            for i, force in enumerate(self.Forces, start=1):
                print(f"\t... force: {i}")
                print(f"\t... type: {force.type}")
                print(f"\t... head marker: {force.iMarker}")
                print(f"\t... tail marker: {force.jMarker}")
                print(f"\t... head body: {force.iBody}")
                print(f"\t... tail body: {force.jBody}")
                print(f"\t... spring stiffness: {force.k}")
                print(f"\t... undeformed spring length: {force.L0}")
                print(f"\t... undeformed torsional spring angle: {force.theta0}")
                print(f"\t... damping coefficient: {force.dc}")
                print(f"\t... constant actuator force: {force.f_a}")
                print(f"\t... constant actuator torque: {force.T_a}")
                print(f"\t... local force: {', '.join(map(str, force.flocal.flatten()))}")
                print(f"\t... global force: {', '.join(map(str, force.f.flatten()))}")
                print(f"\t... torque: {force.T}")
                print(f"\t... gravity: {force._gravity}")
                print(f"\t... weight: {', '.join(map(str, force._wgt.flatten()))}")
                print(f"\t... function index: {force._iFunct}")
                print("\n")
            
            print(f"-----")
            print(f"joints ")
            print(f"-----")
            for i, joint in enumerate(self.Joints, start=1):
                print(f"\t... joint: {i}")
                print(f"\t... type: {joint.type}")
                print(f"\t... body i: {joint.iBody}")
                print(f"\t... body j: {joint.jBody}")
                print(f"\t... iMarker: {joint.iMarker}")
                print(f"\t... jMarker: {joint.jMarker}")
                print(f"\t... function index: {joint.iFunct}")
                print(f"\t... length: {joint.L}")
                print(f"\t... radius: {joint.R}")
                print(f"\t... initial condition x: {joint.x0}")
                print(f"\t... initial condition d: {', '.join(map(str, joint.d0)) if hasattr(joint.d0, '__iter__') else joint.d0}")
                print(f"\t... fix: {joint.fix}")
                print(f"\t... q0: {joint.q0}")
                print(f"\t... initial condition phi: {joint._p0}")
                print(f"\t... number of bodies: {joint._nbody}")
                print(f"\t... number of rows: {joint._mrows}")
                print(f"\t... row start: {joint._rows}")
                print(f"\t... row end: {joint._rowe}")
                print(f"\t... column i start: {joint._colis}")
                print(f"\t... column i end: {joint._colie}")
                print(f"\t... column j start: {joint._coljs}")
                print(f"\t... column j end: {joint._colje}")
                print(f"\t... lagrange multipliers: {', '.join(map(str, joint._lagrange.flatten()))}")
                print("\n")

    def _update_position(self):
        """
        Update position entities.
        """
        for body in self.Bodies:
            body._A = A_matrix(body.p)

        for body in self.Bodies:
            for marker in body._markers:
                pos_col = marker.position.reshape(2, 1)
                marker._sP  = body._A @ pos_col
                marker._sPr = s_rot(marker._sP)
                marker._rP  = body.r + marker._sP
                if marker.has_orientation:
                    marker._u  = body._A @ marker._ulocal
                    marker._ur = s_rot(marker._u)

    def _update_velocity(self):
        """
        Compute sP_dot and rP_dot vectors and update velocity components.
        """
        for body in self.Bodies:
            for marker in body._markers:
                marker._dsP = marker._sPr * body.dp
                marker._drP = body.dr + marker._dsP
                if marker.has_orientation:
                    marker._du = marker._ur * body.dp
            
    def _compute_constraints(self):
        if not self.Joints:
            return np.zeros((0, 1))
        nConst = self.Joints[-1]._rowe
        phi = np.zeros([nConst, 1])

        for joint in self.Joints:
            match joint.type:
                case 'rev':
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    
                    # compute relative positions of the points
                    f = iPt._rP - jPt._rP

                    if joint.fix == 1:
                        if joint.iBody is Ground:
                            f = np.append(f, (-joint.jBody.p - joint._p0))
                        elif joint.jBody is Ground:
                            f = np.append(f, (joint.iBody.p - joint._p0))
                        else:
                            f = np.append(f, (joint.iBody.p - joint.jBody.p - joint._p0))

                case 'tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    ujr = joint.jMarker._ur
                    ui = joint.iMarker._u
                    d = iPt._rP - jPt._rP

                    # compute constraint equations
                    f = np.array([ujr.T @ d, ujr.T @ ui]).reshape(2,1)

                    # additional constraint if fixed
                    if joint.fix == 1:
                        f = np.append(f, (ui.T @ d - joint._p0) / 2).reshape(3,1)

                case 'rev-rev':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    d = iPt._rP - jPt._rP
                    L = joint.L
                    u = d/L
                    # compute constraint equations
                    f = (u.T @ d - L)/2

                case 'rev-tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    uir = joint.iMarker._ur
                    d = iPt._rP - jPt._rP

                    # compute constraint equations
                    f = (uir.T @ d - joint.L)

                case 'rigid':
                    Bi = joint.iBody
                    Bj = joint.jBody
                    
                    if Bi is Ground:
                        f = np.vstack([
                            -(Bj.r + Bj._A @ joint.d0),
                            -Bj.p - joint._p0
                        ])
                    elif Bj is Ground:
                        f = np.vstack([
                            Bi.r - joint.d0,
                            Bi.p - joint._p0
                        ])
                    else:
                        f = np.vstack([
                            Bi.r - (Bj.r + Bj._A @ joint.d0),
                            Bi.p - Bj.p - joint._p0
                        ])

                case 'disc':
                    Bi = joint.iBody
                    f = np.vstack([
                        Bi.r[1] - joint.R,
                        (Bi.r[0] - joint.x0) + joint.R * (Bi.p - joint._p0)
                    ])
                    
                case 'rel-rot':
                    fun, fun_d, fun_dd = functEval(joint.iFunct, self.t)
                    Bi = joint.iBody
                    Bj = joint.jBody

                    if Bi is Ground:
                        f = -Bj.p - fun
                    elif Bj is Ground:
                        f = Bi.p - fun
                    else:
                        f = Bi.p - Bj.p - fun

                case 'rel-tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    d = iPt._rP - jPt._rP
                    fun, fun_d, fun_dd = functEval(joint.iFunct, self.t)

                    f = (d.T @ d - fun**2)/2

                case _:
                    raise ValueError(f"Joint type '{joint.type}' is not supported.")

            rs = joint._rows
            re = joint._rowe
            phi[rs:re] = f
            
        return phi

    def _compute_jacobian(self):
        """
        Calculate the Jacobian matrix D for the system constraints.

        Returns
        -------
        D (NDArray)
            The Jacobian matrix of shape (nConst, nB3).
        """
        nConst = self.Joints[-1]._rowe
        nB3 = 3 * len(self.Bodies)
        D = np.zeros((nConst, nB3))
        
        for Ji, joint in enumerate(self.Joints):
            match joint.type:
                case 'rev':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    Di = np.block([
                        [np.eye(2), iPt._sPr.reshape(2, 1)]
                    ])
                    Dj = np.block([
                        [-np.eye(2), -jPt._sPr.reshape(2, 1)]
                    ])

                    if joint.fix == 1:
                        Di = np.vstack([
                            Di,
                            [0, 0, 1]
                        ])
                        Dj = np.vstack([
                            Dj,
                            [0, 0, -1]
                        ])

                case 'tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    uj = joint.jMarker._u
                    ujr = joint.jMarker._ur
                    d = iPt._rP - jPt._rP

                    Di = np.block([
                        [ujr.T, (uj.T @ iPt._sP).reshape(1, 1)],
                        [np.array([0, 0, 1])]
                    ])
                    Dj = np.block([
                        [-ujr.T, -(uj.T @ (jPt._sP + d)).reshape(1, 1)],
                        [np.array([0, 0, -1])]
                    ])

                    if joint.fix == 1:
                        Di = np.vstack([
                            Di,
                            [uj.T, (uj.T @ iPt._sPr).reshape(1)]
                        ])
                        Dj = np.vstack([
                            Dj,
                            [-uj.T, -(uj.T @ jPt._sPr).reshape(1)]
                        ])

                case 'rev-rev':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    d = iPt._rP - jPt._rP
                    L = joint.L
                    u = d/L

                    Di = np.block([
                        u.T, (u.T @ iPt._sPr).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -u.T, -(u.T @ jPt._sPr).reshape(1, 1)
                        ])
                    
                case 'rev-tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    ui = joint.iMarker._u
                    ui_r = joint.iMarker._ur
                    d = iPt._rP - jPt._rP

                    Di = np.block([
                        ui_r.T, (ui.T @ (iPt._sP - d)).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -ui_r.T, -(ui.T @ jPt._sP).reshape(1, 1)
                        ])

                case 'rigid':
                    Bj = joint.jBody

                    Di = np.eye(3)
                    if Bj is not Ground:
                        Dj = np.block([
                            [-np.eye(2), -s_rot(Bj._A @ joint.d0)],
                            [np.array([0, 0, -1])]
                        ])
                        
                case 'disc':
                    Di = np.array([
                        [0, 1, 0],
                        [1, 0, joint.R]
                    ])
                    
                case 'rel-rot':
                    Di = np.array([
                        [0, 0, 1]
                        ])
                    Dj = np.array([
                        [0, 0, -1]
                        ])

                case 'rel-tran':
                    iPt = joint.iMarker
                    jPt = joint.jMarker

                    d = iPt._rP - jPt._rP

                    Di = np.block([
                        d.T, (d.T @ iPt._sPr).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -d.T, -(d.T @ jPt._sPr).reshape(1, 1)
                        ])
                
                case _:
                    raise ValueError(f"Joint type '{joint.type}' is not supported.")

            # row indices for the current joint in the Jacobian matrix
            rs = joint._rows
            re = joint._rowe

            # column indices for body i 
            if joint.iBody is not Ground:
                cis = joint._colis
                cie = joint._colie
                D[rs:re, cis:cie] = Di

            # column indices for body j
            if joint.jBody is not Ground:
                cjs = joint._coljs
                cje = joint._colje
                D[rs:re, cjs:cje] = Dj

        return D

    def _rhs_velocity(self):
        """
        Calculate the right-hand side velocity vector for the system constraints.
        
        Returns
        -------
        rhsv : numpy.ndarray
            Right-hand side of velocity constraints.
        """
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        rhsv = np.zeros((nConst, 1))

        for joint in self.Joints:
            match joint.type:
                case 'rel-rot':
                    fun, fun_d, _ = functEval(joint.iFunct, self.t)
                    f = fun_d

                case 'rel-tran':
                    fun, fun_d, _ = functEval(joint.iFunct, self.t)
                    d = joint.iMarker._rP - joint.jMarker._rP
                    f = fun * fun_d

                case _:
                    continue

            rs = joint._rows
            re = joint._rowe
            rhsv[rs:re] = f

        return rhsv

    def _rhs_acceleration(self):
        """
        Compute the right-hand side of acceleration constraints.

        Returns
        -------
        numpy.ndarray
            Right-hand side of acceleration constraints (gamma).
        """
        nConst = self.Joints[-1]._rowe
        rhsa = np.zeros([nConst, 1])
        
        for Ji, joint in enumerate(self.Joints):
            joint_type = joint.type

            match joint_type:
                case "rev":
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    Bi = joint.iBody
                    Bj = joint.jBody

                    if Bi is Ground:
                        f = s_rot(jPt._dsP) * Bj.dp
                    elif Bj is Ground:
                        f = -s_rot(iPt._dsP) * Bi.dp
                    else:
                        f = (
                            -s_rot(iPt._dsP) * Bi.dp
                            + s_rot(jPt._dsP) * Bj.dp
                        )

                    if joint.fix == 1:
                        f = np.vstack([f, [0]])
                
                case "tran":
                    Bi = joint.iBody
                    Bj = joint.jBody
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    ujd = joint.jMarker._du
                    ujdr = s_rot(ujd)

                    if Bi is Ground:
                        f2 = 0.0
                    elif Bj is Ground:
                        f2 = 0.0
                    else:
                        diffr = Bi.r - Bj.r
                        dp_product = (ujd.T @ diffr).item() * Bi.dp
                        diffdr = Bi.dr - Bj.dr
                        f2 = dp_product - 2.0 * (ujdr.T @ diffdr).item()

                    f = np.array([[f2], [0.0]])

                    if joint.fix == 1:
                        d = iPt._rP - jPt._rP
                        dd = iPt._drP - jPt._drP
                        L = joint._p0 
                        u = d / L
                        du = dd / L
                        f3 = -(du.T @ dd).item()

                        if Bi is Ground:
                            f3 += (u.T @ (s_rot(jPt._dsP) * Bj.dp)).item()
                        elif Bj is Ground:
                            f3 -= (u.T @ (s_rot(iPt._dsP) * Bi.dp)).item()
                        else:
                            term1 = iPt._dsP * Bi.dp
                            term2 = jPt._dsP * Bj.dp
                            f3 -= (u.T @ s_rot(term1 - term2)).item()

                        f = np.vstack([f, [[f3]]])
                    
                case "rev-rev":
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    Bi = joint.iBody
                    Bj = joint.jBody
                    
                    d = iPt._rP - jPt._rP
                    dd = iPt._drP - jPt._drP
                    
                    L = joint.L
                    u = d/L
                    ud = dd/L
                    
                    f = -ud.T @ dd
                    
                    if Bi is Ground:
                        f = f + u.T @ s_rot(jPt._dsP) * Bj.dp
                    elif Bj is Ground:
                        f = f - u.T @ s_rot(iPt._dsP) * Bi.dp
                    else:
                        f = f - u.T @ s_rot(
                            iPt._dsP * Bi.dp - 
                            jPt._dsP * Bj.dp
                    )

                case "rev-tran":
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    Bi = joint.iBody
                    Bj = joint.jBody

                    ui = joint.iMarker._u
                    ui_d = joint.iMarker._du
                    d = iPt._rP - jPt._rP
                    dd = iPt._drP - jPt._drP

                    if Bi is Ground:
                        f = ui.T @ jPt._dsP * Bj.dp
                    elif Bj is Ground:
                        f = ui_d.T @ (d * Bi.dp + 2 * s_rot(dd)) - \
                            ui.T @ iPt._dsP * Bi.dp
                    else:
                        f = ui_d.T @ (d * Bi.dp + 2 * s_rot(dd)) - \
                            ui.T @ (iPt._dsP * Bi.dp - \
                                jPt._dsP * Bj.dp)
                
                case "rigid":
                    Bj = joint.jBody

                    f = np.zeros(3)
                    if Bj is not Ground:
                        f = np.concatenate([
                            -Bj._A @ joint.d0 * Bj.dp**2,
                            np.array([0])
                        ])
                
                case "disc":
                    f = np.zeros(2)
                    
                case "rel-rot":
                    fun, fun_d, fun_dd = functEval(joint.iFunct, self.t)
                    f = fun_dd

                case "rel-tran":
                    iPt = joint.iMarker
                    jPt = joint.jMarker
                    Bi = joint.iBody
                    Bj = joint.jBody

                    d = iPt._rP - jPt._rP
                    dd = iPt._drP - jPt._drP

                    fun, fun_d, fun_dd = functEval(joint.iFunct, self.t)

                    f = fun * fun_dd + fun_d**2

                    if Bi is Ground:
                        f = f + d.T @ s_rot(jPt._dsP).T @ Bj.dp
                    elif Bj is Ground:
                        f = f - d.T @ s_rot(iPt._dsP).T @ Bi.dp - dd.T @ dd
                    else:
                        f = f + d.T @ s_rot(jPt._dsP).T @ Bj.dp \
                            - d.T @ s_rot(iPt._dsP).T @ Bi.dp - dd.T @ dd
            
            rs = joint._rows
            re = rs + joint._mrows 
            rhsa[rs:re] = np.asarray(f).reshape(-1, 1)
        
        return rhsa

    def _bodies2u(self):
        """ 
        Pack coordinates and velocities into the u array.
        """
        nB = len(self.Bodies)
        u = np.zeros([(3 * nB * 2), 1])

        for Bi in range(nB):
            ir = self.Bodies[Bi]._irc
            ird = self.Bodies[Bi]._irv
            u[ir:ir+3] = np.block([[self.Bodies[Bi].r],[self.Bodies[Bi].p]])
            u[ird:ird+3] = np.block([[self.Bodies[Bi].dr], [self.Bodies[Bi].dp]])
        
        return u

    def _bodies2ud(self):
        """ 
        Pack velocities and accelerations into ud. 
        """
        nB6 = 6 * len(self.Bodies)
        ud = np.zeros([nB6, 1])

        for Bi, body in enumerate(self.Bodies):
            ir = body._irc
            ird = body._irv
            ud[ir:ir + 3] = np.vstack([body.dr, body.dp]).reshape(3, 1)
            ud[ird:ird + 3] = np.vstack([body.ddr, body.ddp]).reshape(3, 1)

        return ud
    
    def _u2bodies(self, u):
        """
        Unpack u into coordinate and velocity sub-arrays.
        """ 
        # check on "u" shape to avoid errors during the simulation
        if u.ndim != 2:
            u = u.reshape(-1, 1)
            
        nB = len(self.Bodies)
        for Bi in range(nB): 
            ir = self.Bodies[Bi]._irc
            ird = self.Bodies[Bi]._irv
            self.Bodies[Bi].r  = u[ir:ir+2]
            self.Bodies[Bi].p  = u[ir+2][0]
            self.Bodies[Bi].dr = u[ird:ird+2]
            self.Bodies[Bi].dp = u[ird+2][0]

    def _compute_force(self):
        """
        Compute and return the array of forces acting on the system at time t.
        """
        for body in self.Bodies:
            body._f = colvect([0.0, 0.0]) # initialize body force vectors
            body._n = 0.0                 # initialize body torque (moment) scalar

        # loop over all forces and apply them to the appropriate bodies
        for Fi, force in enumerate(self.Forces):
            match force.type:
                case 'weight':
                    for body in self.Bodies:
                        body._f += body._wgt

                case 'ptp':
                    iPt = force.iMarker
                    jPt = force.jMarker
                    Bi = force.iBody
                    Bj = force.jBody
                    d = iPt._rP - jPt._rP
                    dd = iPt._drP - jPt._drP
                    L = np.linalg.norm(d)
                    dL = d.T @ dd / L
                    delta = L - force.L0
                    u = d / L

                    f = force.k * delta + force.dc * dL + force.f_a
                    fi = f * u

                    if Bi is not Ground:
                        Bi._f -= fi
                        Bi._n -= (iPt._sPr.T @ fi).item()
                    
                    if Bj is not Ground:
                        Bj._f += fi
                        Bj._n += (jPt._sPr.T @ fi).item()

                case 'rot-sda':
                    Bi = force.iBody
                    Bj = force.jBody

                    if Bi is Ground:
                        theta = -Bj.p
                        theta_d = -Bj.dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        Bj._n += T
                    elif Bj is Ground:
                        theta = Bi.p
                        theta_d = Bi.dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        Bi._n -= T
                    else:
                        theta = Bi.p - Bj.p
                        theta_d = Bi.dp - Bj.dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        Bi._n -= T
                        Bj._n += T

                case 'flocal':
                    force.iBody._f += force.iBody._A @ force.flocal

                case 'f':
                    force.iBody._f += force.f

                case 'T':
                    force.iBody._n += force.T

                case 'user' if force.callback is not None and callable(force.callback):
                    force.callback()

                case _:
                    raise ValueError(f"Unsupported force type: '{force.type}'. Please check your input.")

        # Build force array
        nB3 = 3 * len(self.Bodies)
        g = np.zeros([nB3, 1])
        for Bi, body in enumerate(self.Bodies):
            ks = body._irc
            ke = ks + 3
            g[ks:ke] = np.vstack([body._f, body._n])

        return g

    def _ic_correct(self):
        """
        Corrects initial conditions on the body coordinates and velocities.
        """
        flag = False

        # position correction
        for _ in range(50):
            self._update_position()            # update position entities
            Phi = self._compute_constraints()  # evaluate constraints
            D = self._compute_jacobian()       # evaluate Jacobian
            ff = np.sqrt(Phi.T @ Phi)           # are the constraints violated?

            if ff < 1.0e-10:
                flag = True
                break

            # solve for corrections
            delta_c = -D.T @ np.linalg.solve(D @ D.T, Phi)

            # correct estimates
            nB = len(self.Bodies)
            for Bi in range(nB):
                ir = 3 * Bi
                self.Bodies[Bi].r = self.Bodies[Bi].r + delta_c[ir:ir + 2]
                self.Bodies[Bi].p = self.Bodies[Bi].p + delta_c[ir + 2][0]

        if not flag:
            raise ValueError("Convergence failed in Newton-Raphson!")

        # velocity correction
        nB = len(self.Bodies)
        Phi = np.zeros([3 * nB, 1])
        for Bi in range(nB):
            ir = 3 * Bi
            Phi[ir:ir + 2] = self.Bodies[Bi].dr
            Phi[ir + 2] = self.Bodies[Bi].dp

        rhsv = self._rhs_velocity()
        
        # solve for corrections
        delta_v = -D.T @ np.linalg.solve(D @ D.T, D @ Phi - rhsv)  

        # move corrected velocities to sub-arrays
        for Bi in range(nB):
            ir = 3 * Bi
            self.Bodies[Bi].dr = self.Bodies[Bi].dr + delta_v[ir:ir + 2]
            self.Bodies[Bi].dp = self.Bodies[Bi].dp + delta_v[ir + 2][0]

        coords = np.zeros((nB, 3))
        vels = np.zeros((nB, 3))
        for Bi in range(nB):
            coords[Bi, :] = np.hstack((self.Bodies[Bi].r.T, np.array(self.Bodies[Bi].p).reshape(-1, 1)))
            vels[Bi, :] = np.hstack((self.Bodies[Bi].dr.T, np.array(self.Bodies[Bi].dp).reshape(-1, 1)))

        if self.verbose:
            print("\t... initial conditions corrected ...")
            print(f"-----")
            print(f"bodies ")
            print(f"-----")
            for i, body in enumerate(self.Bodies, start=1):
                print(f"\t... body: {i}")
                print(f"\t... mass: {body.m}")
                print(f"\t... moment of inertia: {body.J}")
                print(f"\t... position: {', '.join(map(str, body.r.flatten()))} # [UPDATED]")
                print(f"\t... orientation: {body.p} # [UPDATED]")
                print(f"\t... velocity: {', '.join(map(str, body.dr.flatten()))} # [UPDATED]")
                print(f"\t... angular velocity: {body.dp} # [UPDATED]")
                print(f"\t... acceleration: {', '.join(map(str, body.ddr.flatten()))}")
                print(f"\t... angular acceleration: {body.ddp}")
                print(f"\t... rotational matrix: {', '.join(map(str, body._A.flatten()))}")
                print(f"\t... inverse mass: {body._invm}")
                print(f"\t... inverse moment of inertia: {body._invJ}")
                print(f"\t... weight: {', '.join(map(str, body._wgt.flatten()))}")
                print(f"\t... force: {', '.join(map(str, body._f.flatten()))}")
                print(f"\t... torque: {body._n}")
                print(f"\t... markers: {body._markers}")
                print(f"\t... irc: {body._irc}")
                print(f"\t... irv: {body._irv}")
                print(f"\t... ira: {body._ira}")
                print("\n")
            
            print(f"-----")
            print(f"markers ")
            print(f"-----")
            for i, body in enumerate(self.Bodies, start=1):
                for j, marker in enumerate(body._markers):
                    print(f"\t... body {i}, marker {j}: {marker.name}")
                    print(f"\t... position (local): {marker.position}")
                    print(f"\t... theta: {marker.theta}")
                    print(f"\t... _rP: {', '.join(map(str, marker._rP.flatten()))}")
                    print(f"\t... _sP: {', '.join(map(str, marker._sP.flatten()))}")
                    print(f"\t... _sPr: {', '.join(map(str, marker._sPr.flatten()))}")
                    if marker.has_orientation:
                        print(f"\t... _u: {', '.join(map(str, marker._u.flatten()))}")
                        print(f"\t... _ur: {', '.join(map(str, marker._ur.flatten()))}")
                    print("\n")
            
            print(f"-----")
            print(f"forces ")
            print(f"-----")
            for i, force in enumerate(self.Forces, start=1):
                print(f"\t... force: {i}")
                print(f"\t... type: {force.type}")
                print(f"\t... head marker: {force.iMarker}")
                print(f"\t... tail marker: {force.jMarker}")
                print(f"\t... head body: {force.iBody}")
                print(f"\t... tail body: {force.jBody}")
                print(f"\t... spring stiffness: {force.k}")
                print(f"\t... undeformed spring length: {force.L0}")
                print(f"\t... undeformed torsional spring angle: {force.theta0}")
                print(f"\t... damping coefficient: {force.dc}")
                print(f"\t... constant actuator force: {force.f_a}")
                print(f"\t... constant actuator torque: {force.T_a}")
                print(f"\t... local force: {', '.join(map(str, force.flocal.flatten()))}")
                print(f"\t... global force: {', '.join(map(str, force.f.flatten()))}")
                print(f"\t... torque: {force.T}")
                print(f"\t... gravity: {force._gravity}")
                print(f"\t... weight: {', '.join(map(str, force._wgt.flatten()))}")
                print(f"\t... function index: {force._iFunct}")
                print("\n")
            
            print(f"-----")
            print(f"joints ")
            print(f"-----")
            for i, joint in enumerate(self.Joints, start=1):
                print(f"\t... joint: {i}")
                print(f"\t... type: {joint.type}")
                print(f"\t... body i: {joint.iBody}")
                print(f"\t... body j: {joint.jBody}")
                print(f"\t... iMarker: {joint.iMarker}")
                print(f"\t... jMarker: {joint.jMarker}")
                print(f"\t... function index: {joint.iFunct}")
                print(f"\t... length: {joint.L}")
                print(f"\t... radius: {joint.R}")
                print(f"\t... initial condition x: {joint.x0}")
                print(f"\t... initial condition d: {', '.join(map(str, joint.d0)) if hasattr(joint.d0, '__iter__') else joint.d0}")
                print(f"\t... fix: {joint.fix}")
                print(f"\t... q0: {joint.q0}")
                print(f"\t... initial condition phi: {joint._p0}")
                print(f"\t... number of bodies: {joint._nbody}")
                print(f"\t... number of rows: {joint._mrows}")
                print(f"\t... row start: {joint._rows}")
                print(f"\t... row end: {joint._rowe}")
                print(f"\t... column i start: {joint._colis}")
                print(f"\t... column i end: {joint._colie}")
                print(f"\t... column j start: {joint._coljs}")
                print(f"\t... column j end: {joint._colje}")
                print(f"\t... lagrange multipliers: {', '.join(map(str, joint._lagrange.flatten()))}")
                print("\n")
        else:
            print("\n\t Corrected coordinates")
            print("\t", f"{'x':^12}{'y':^12}{'phi':^12}")
            for row in coords:
                print(f"\t {row[0]:^12.5f}{row[1]:^12.5f}{row[2]:^12.5f}")

            print("\n\t Corrected velocities")
            print("\t", f"{'x-dot':^12}{'y-dot':^12}{'phi-dot':^12}")
            for row in vels:
                print(f"\t {row[0]:^12.5f}{row[1]:^12.5f}{row[2]:^12.5f}")
            print("\n")

    def _analysis(self, t, u):
        """
        Solve the constrained equations of motion at time t with the standard
        Lagrange multiplier method.
        """        
        self._num += 1                  # increment the number of function evaluations
        self.t = t                      # store current time for force/constraint callbacks
        nB3 = 3 * len(self.Bodies)
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        self._u2bodies(u)               # unpack u into coordinate and velocity sub-arrays
        self._update_position()
        self._update_velocity()
        h_a = self._compute_force()     # array of applied forces

        if nConst == 0:
            ddc = self.invM_array.reshape(-1, 1) * h_a  # solve for accelerations
            Lambda = np.array([])  # no constraints, no multipliers
        else:
            D = self._compute_jacobian()
            rhsA = self._rhs_acceleration()  # right-hand side of acceleration constraints (gamma)

            # GGL regularization: replace zero block with (1/μ)·I
            if self._ggl_mu > 0:
                reg_block = (1.0 / self._ggl_mu) * np.eye(nConst)
            else:
                reg_block = np.zeros((nConst, nConst))

            DMD = np.block([
                [self.M_matrix, -D.T], 
                [D, reg_block]
            ])

            # Baumgarte stabilization: add -2α·dΦ - β²·Φ to RHS
            if self._baumgarte_alpha > 0 or self._baumgarte_beta > 0:
                Phi = np.asarray(self._compute_constraints()).flatten()
                dq = u[nB3:]  # velocity part of state vector
                dPhi = np.asarray(D @ dq).flatten()
                gamma_stab = (rhsA.flatten()
                              - 2.0 * self._baumgarte_alpha * dPhi
                              - self._baumgarte_beta**2 * Phi)
            else:
                gamma_stab = rhsA.flatten()

            rhs = np.concatenate([h_a.flatten(), gamma_stab])

            #* check on conditioned index of the coefficient matrix
            cond_number = np.linalg.cond(DMD)
            if cond_number > 1e12:
                print(f"Warning: DMD matrix is poorly conditioned with condition number {cond_number}")
    
            # solve the system of equations
            sol = np.linalg.solve(DMD, rhs).reshape(-1, 1)
            ddc = sol[:nB3]
            Lambda = sol[nB3:]

        # update accelerations for each body
        for Bi, body in enumerate(self.Bodies):
            ir = body._irc
            i2 = ir + 2
            i3 = i2
            body.ddr = ddc[ir:i2]
            body.ddp = ddc[i3][0]

        ud = self._bodies2ud()              # pack velocities and accelerations into ud
        return ud.flatten()

    def _taqaddum(self, t_initial, t_final, pbar):
        """
        Restituisce una funzione wrapper per _analysis con progresso ottimizzato
        """
        last_progress = 0
        
        def _wrapp_analysis(t, u):
            nonlocal last_progress
            progress = min(100, int(100 * (t - t_initial) / (t_final - t_initial)))
            if progress > last_progress:
                pbar.n = progress
                pbar.refresh()
                last_progress = progress
            return self._analysis(t, u)
        return _wrapp_analysis

    def solve(self, method="LSODA", t_final=None, dt=None, ic_correct=False,
              t_eval=None, t_span=None,
              baumgarte_alpha=5.0, baumgarte_beta=5.0, ggl_penalty=1e8):
        """Solve equations of motion.

        Supports both interactive (legacy) and programmatic (non-interactive)
        modes. Pass ``t_final`` directly to skip all ``input()`` prompts.

        Parameters
        ----------
        method : str, optional
            ODE solver method (default "LSODA"). See scipy.integrate.solve_ivp.
        t_final : float, optional
            Final simulation time. If None, prompts user interactively.
        dt : float, optional
            Output time step. Used only when ``t_eval`` is not given.
        ic_correct : bool, optional
            Whether to correct initial conditions before solving
            (default False). Only applies in non-interactive mode.
        t_eval : array-like, optional
            Explicit array of output time points. Overrides ``dt``.
        t_span : tuple, optional
            ``(t_start, t_end)`` shorthand; sets ``t_final = t_span[1]``.
        baumgarte_alpha : float, optional
            Baumgarte velocity-level stabilization gain (default 5.0).
            Controls damping of constraint velocity drift.
        baumgarte_beta : float, optional
            Baumgarte position-level stabilization gain (default 5.0).
            Controls correction of constraint position drift.
        ggl_penalty : float, optional
            GGL (Gear-Gupta-Leimkuhler) regularization parameter (default 1e8).
            Replaces the zero block with ``(1/ggl_penalty)*I`` in the augmented
            matrix, preventing singularity with redundant constraints.
            Set to 0 to disable regularization.

        Returns
        -------
        SolResult
            Object supporting tuple unpacking (``T, uT = sol``) and
            attribute access (``sol.t``, ``sol.y`` in scipy convention).
        """
        self.method = method

        # Store stabilization parameters for _analysis and _post_process
        self._baumgarte_alpha = baumgarte_alpha
        self._baumgarte_beta = baumgarte_beta
        self._ggl_mu = ggl_penalty

        # Handle t_span shorthand
        if t_span is not None and t_final is None:
            t_final = t_span[1]

        nConst = self.Joints[-1]._rowe if self.Joints else 0

        if t_final is None:
            # --- Interactive (legacy) mode ---
            print("\n")
            ans = input("\t... Do you want to correct the initial conditions? [(y)es/(n)o] ").lower()
        else:
            # --- Programmatic (non-interactive) mode ---
            ans = 'y' if ic_correct else 'n'

        if nConst != 0:
            self.t = 0.0  # needed by _compute_constraints / _rhs_velocity when ic_correct runs
            if ans == 'y':
                self._ic_correct()
            D = self._compute_jacobian()
            redund = np.linalg.matrix_rank(D)
            if redund < nConst:
                print("\n\t...Redundancy in the constraints")

        u = self._bodies2u()
        if self.verbose:
            header = "... initial u vector ..."
            print(f"\n\t{header}")
            header_width = len(header)
            formatted_u = [f"{float(element):.2f}" for element in u]
            for element in formatted_u:
                print(f"\t{element:^{header_width}}")
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            raise ValueError("\t ... check initial conditions, \"u\" vector contains NaN or Inf values.")

        t_initial = 0.0
        self._num = 0  # initialize the number of function evaluations

        if t_final is None:
            t_final = float(input("\n\t ...Final time = ? "))

        dense_sol = None
        if t_final == 0:
            self._analysis(0, u)
            T = np.array([0.0])
            uT = u.T
        else:
            if t_eval is not None:
                Tspan = np.asarray(t_eval, dtype=float)
            elif dt is not None:
                Tspan = np.arange(t_initial, t_final + dt * 0.5, dt)
            else:
                dt_input = float(input("\t ...Reporting time-step = ? "))
                Tspan = np.arange(t_initial, t_final, dt_input)

            u0 = u.flatten()
            options = {'rtol': 1e-6, 'atol': 1e-9,
                       'max_step': float(Tspan[1] - Tspan[0])}

            pbar = tqdm(total=100, desc="         ...Simulation progress",
                        bar_format="{l_bar}{bar}| [Elapsed time: {elapsed}, Remaining time: {remaining}]",
                        colour="green")

            _wrapp_analysis = self._taqaddum(t_initial, t_final, pbar)

            try:
                _sol = solve_ivp(_wrapp_analysis,
                                 [t_initial, t_final],
                                 u0,
                                 t_eval=Tspan,
                                 method=self.method,
                                 dense_output=True,
                                 **options)
            finally:  # ensure progress bar is closed even on error
                pbar.close()

            T = _sol.t
            uT = _sol.y.T
            dense_sol = _sol.sol

        print(f"\n ")
        print(f"\t ...Number of function evaluations: {self._num}")
        print(f"\t ...Simulation completed successfully!")
        print(f"\n ")
        return SolResult(T, uT, model=self, dense_sol=dense_sol)

    def _post_process(self, T, uT):
        """
        Recalculate accelerations and Lagrange multipliers on exact t_eval grid.

        Parameters
        ----------
        T : np.ndarray, shape (nSteps,)
            Time vector from solve_ivp (exact t_eval points).
        uT : np.ndarray, shape (nSteps, 2*nB3)
            State matrix from solve_ivp.

        Returns
        -------
        accelerations : np.ndarray, shape (nSteps, nB3)
            Generalized accelerations at each time step.
        reactions : np.ndarray, shape (nSteps, nConstraints)
            Lagrange multipliers at each time step.
        """
        nB3 = 3 * self.nB
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        nSteps = len(T)

        accelerations = np.zeros((nSteps, nB3))
        reactions = np.zeros((nSteps, nConst))

        for i in range(nSteps):
            t_i = T[i]
            u_i = uT[i]

            self.t = t_i
            self._u2bodies(u_i)
            self._update_position()
            self._update_velocity()
            h_a = self._compute_force()

            if nConst == 0:
                ddc = self.invM_array.reshape(-1, 1) * h_a
                Lambda = np.array([])
            else:
                D = self._compute_jacobian()
                rhsA = self._rhs_acceleration()

                # GGL regularization: replace zero block with (1/μ)·I
                if self._ggl_mu > 0:
                    reg_block = (1.0 / self._ggl_mu) * np.eye(nConst)
                else:
                    reg_block = np.zeros((nConst, nConst))

                DMD = np.block([
                    [self.M_matrix, -D.T],
                    [D, reg_block]
                ])

                # Baumgarte stabilization: add -2α·dΦ - β²·Φ to RHS
                if self._baumgarte_alpha > 0 or self._baumgarte_beta > 0:
                    Phi = np.asarray(self._compute_constraints()).flatten()
                    dq_i = u_i[nB3:]
                    dPhi = np.asarray(D @ dq_i).flatten()
                    gamma_stab = (rhsA.flatten()
                                  - 2.0 * self._baumgarte_alpha * dPhi
                                  - self._baumgarte_beta**2 * Phi)
                else:
                    gamma_stab = rhsA.flatten()

                rhs = np.concatenate([h_a.flatten(), gamma_stab])
                sol = np.linalg.solve(DMD, rhs)
                ddc = sol[:nB3]
                Lambda = sol[nB3:]

            accelerations[i] = ddc.flatten()
            reactions[i] = Lambda.flatten()

        return accelerations, reactions

    # ── Deprecated methods ─────────────────────────────────────────

    def get_reactions(self):
        """DEPRECATED: Use SolResult.reactions property instead."""
        import warnings
        warnings.warn(
            "get_reactions() is deprecated. Use the .reactions property "
            "on the SolResult object returned by solve().",
            DeprecationWarning,
            stacklevel=2
        )
        return None

    def get_accelerations(self):
        """DEPRECATED: Use SolResult.accelerations property instead."""
        import warnings
        warnings.warn(
            "get_accelerations() is deprecated. Use the .accelerations property "
            "on the SolResult object returned by solve().",
            DeprecationWarning,
            stacklevel=2
        )
        return None
