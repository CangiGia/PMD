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
import inspect
from .utils import *
from .mechanics import *
from .builder import *
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
        """Plot x, y, \u03c6 for selected bodies.

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
        """Plot dx, dy, d\u03c6 for selected bodies.

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
        """Plot ddx, ddy, dd\u03c6 for selected bodies. Triggers lazy computation.

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
    def __init__(self, verbose = False):
        grouped_calsses = group_classes()
        self.verbose = verbose
        self.Bodies = grouped_calsses.get("Body", [])
        self.Points = grouped_calsses.get("Point", [])
        self.uVectors = grouped_calsses.get("uVector", [])
        self.Forces = grouped_calsses.get("Force", [])
        self.Joints = grouped_calsses.get("Joint", [])
        self.Functs = grouped_calsses.get("Function", []) if "Function" in grouped_calsses else []
        
        # initialize the model for simulation automatically
        self.__initialize()

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

    def __initialize(self):
        """
        Initializi the multi-body model considering the values defined 
        by the user.
        """
        # initialize variables
        nB = len(self.Bodies)
        nB3 = 3 * nB
        nB6 = 6 * nB

        #// bodies
        for Bi in range(nB):
            body = self.Bodies[Bi]
            body._irc = 3 * Bi
            body._irv = nB3 + 3 * Bi
            body._invm = 1 / body.m
            body._invJ = 1 / body.J
            body._A = A_matrix(body.p)

        # mass (inertia) array and pre-computed diagonal matrix
        self.M_array = np.zeros(nB3)
        self.invM_array = np.zeros(nB3)
        for Bi in range(nB):
            is_ = 3 * Bi
            ie_ = is_ + 3
            self.M_array[is_:ie_] = np.array([self.Bodies[Bi].m, self.Bodies[Bi].m, self.Bodies[Bi].J])
            self.invM_array[is_:ie_] = np.array([self.Bodies[Bi]._invm, self.Bodies[Bi]._invm, self.Bodies[Bi]._invJ])
        self.M_matrix = np.diag(self.M_array)

        #// points
        nPtot = len(self.Points)
        for Pi in range(nPtot):
            point = self.Points[Pi]
            if point.Bindex == 0:
                point._sP = point.sPlocal
                point._sPr = s_rot(point._sP)
                point._rP = point._sP

            for Bi in range(nB):
                if point.Bindex == (Bi+1):
                    self.Bodies[Bi]._pts.append(Pi)  # append point index to the body's points

        #// unit vectors
        nU = len(self.uVectors)
        for Vi in range(nU):
            unit_vector = self.uVectors[Vi]
            if unit_vector.Bindex == 0:
                unit_vector._u = unit_vector.ulocal
                unit_vector._ur = s_rot(unit_vector._u)

        #// force elements
        nF = len(self.Forces)
        for Fi in range(nF):
            force = self.Forces[Fi]
            if force.type == 'weight':
                ug = force._gravity * force._wgt
                for Bi in range(nB):
                    self.Bodies[Bi]._wgt = self.Bodies[Bi].m * ug
            elif force.type == 'ptp':
                Pi = force.iPindex
                Pj = force.jPindex
                force.iBindex = self.Points[Pi].Bindex
                force.jBindex = self.Points[Pj].Bindex

        #// joints
        nJ = len(self.Joints)
        cfriction = 0

        for Ji in range(nJ):
            joint = self.Joints[Ji]
            joint_type = joint.type

            match joint_type:
                case 'rev':
                    joint._mrows = 2
                    joint._nbody = 2
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    joint.iBindex = self.Points[Pi].Bindex
                    joint.jBindex = self.Points[Pj].Bindex
                    if joint.fix == 1:
                        joint._mrows = 3
                        if joint.iBindex == 0:
                            joint._p0 = -self.Bodies[joint.jBindex].p
                        elif joint.jBindex == 0:
                            joint._p0 = self.Bodies[joint.iBindex].p
                        else:
                            joint._p0 = self.Bodies[joint.iBindex].p - self.Bodies[joint.jBindex].p

                case 'tran':
                    joint._mrows = 2
                    joint._nbody = 2
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    joint.iBindex = self.Points[Pi].Bindex
                    joint.jBindex = self.Points[Pj].Bindex
                    if joint.fix == 1:
                        joint._mrows = 3
                        if joint.iBindex == 0:
                            joint._p0 = np.linalg.norm(self.Points[Pi]._rP - 
                                                    self.Bodies[joint.jBindex].r - 
                                                    self.Bodies[joint.jBindex]._A @ 
                                                    self.Points[Pj].sPlocal)
                        elif joint.jBindex == 0:
                            joint._p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + 
                                                    self.Bodies[joint.iBindex]._A @ 
                                                    self.Points[Pi].sPlocal - 
                                                    self.Points[Pj]._rP)
                        else:
                            joint._p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + 
                                                    self.Bodies[joint.iBindex]._A @ 
                                                    self.Points[Pi].sPlocal - 
                                                    self.Bodies[joint.jBindex].r - 
                                                    self.Bodies[joint.jBindex]._A @ 
                                                    self.Points[Pj].sPlocal)

                case 'rev-rev':
                    joint._mrows = 1
                    joint._nbody = 2
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    joint.iBindex = self.Points[Pi].Bindex
                    joint.jBindex = self.Points[Pj].Bindex

                case 'rev-tran':
                    joint._mrows = 1
                    joint._nbody = 2
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    joint.iBindex = self.Points[Pi].Bindex
                    joint.jBindex = self.Points[Pj].Bindex

                case 'rel-rot' | 'rel-tran':
                    joint._mrows = 1
                    joint._nbody = 1
                case 'disc':
                    joint._mrows = 2
                    joint._nbody = 1

                case 'rigid':
                    joint._mrows = 3
                    joint._nbody = 2
                    Bi = joint.iBindex
                    Bj = joint.jBindex
                    if Bi == 0:
                        joint.d0 = -self.Bodies[Bj]._A.T @ self.Bodies[Bj].r
                        joint._p0 = -self.Bodies[Bj].p
                    elif Bj == 0:
                        joint.d0 = self.Bodies[Bi].r
                        joint._p0 = self.Bodies[Bi].p
                    else:
                        joint.d0 = self.Bodies[Bj]._A.T @ (self.Bodies[Bi].r - self.Bodies[Bj].r)
                        joint._p0 = self.Bodies[Bi].p - self.Bodies[Bj].p

                case _:
                    raise ValueError("Joint type doesn't supported!")

        #// functions
        if self.Functs:
            nFc = len(self.Functs)
            for Ci in range(nFc):
                functData(Ci, self.Functs)
        else:
            pass

        # compute number of constraints and determine row/column pointers
        nConst = 0
        for Ji in range(nJ):
            joint = self.Joints[Ji]
            joint._rows = nConst
            joint._rowe = nConst + joint._mrows
            nConst = joint._rowe
            Bi = joint.iBindex
            if Bi != 0:
                joint._colis = 3 * (Bi - 1)
                joint._colie = 3 * Bi
            Bj = joint.jBindex
            if Bj != 0:
                joint._coljs = 3 * (Bj - 1)
                joint._colje = 3 * Bj

        # // ---
        # // ... some check if required ...
        # // ---
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
                print(f"\t... points: {', '.join(map(str, body._pts))}")
                print(f"\t... irc: {body._irc}")
                print(f"\t... irv: {body._irv}")
                print(f"\t... ira: {body._ira}")
                print("\n")
            
            print(f"-----")
            print(f"points ")
            print(f"-----")
            for i, point in enumerate(self.Points):
                print(f"\t... point: {i}")
                print(f"\t... body index: {point.Bindex}")
                print(f"\t... local coordinates: {', '.join(map(str, point.sPlocal.flatten()))}")
                print(f"\t... global coordinates: {', '.join(map(str, point._rP.flatten()))}")
                print(f"\t... sP: {', '.join(map(str, point._sP.flatten()))}")
                print(f"\t... sPr: {', '.join(map(str, point._sPr.flatten()))}")
                print(f"\t... dsP: {', '.join(map(str, point._dsP.flatten()))}")
                print(f"\t... drP: {', '.join(map(str, point._drP.flatten()))}")
                print(f"\t... ddrP: {', '.join(map(str, point._ddrP.flatten()))}")
                print("\n")
            
            print(f"-----")
            print(f"vectors ")
            print(f"-----")
            for i, uvector in enumerate(self.uVectors, start=1):
                print(f"\t... uVector: {i}")
                print(f"\t... body index: {uvector.Bindex}")
                print(f"\t... local vector: {', '.join(map(str, uvector.ulocal.flatten()))}")
                print(f"\t... global vector: {', '.join(map(str, uvector._u.flatten()))}")
                print(f"\t... ur: {', '.join(map(str, uvector._ur.flatten()))}")
                print(f"\t... du: {', '.join(map(str, uvector._du.flatten()))}")
                print("\n")
            
            print(f"-----")
            print(f"forces ")
            print(f"-----")
            for i, force in enumerate(self.Forces, start=1):
                print(f"\t... force: {i}")
                print(f"\t... type: {force.type}")
                print(f"\t... head point index: {force.iPindex}")
                print(f"\t... tail point index: {force.jPindex}")
                print(f"\t... head body index: {force.iBindex}")
                print(f"\t... tail body index: {force.jBindex}")
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
                print(f"\t... body i index: {joint.iBindex}")
                print(f"\t... body j index: {joint.jBindex}")
                print(f"\t... point i index: {joint.iPindex}")
                print(f"\t... point j index: {joint.jPindex}")
                print(f"\t... unit vector i index: {joint.iUindex}")
                print(f"\t... unit vector j index: {joint.jUindex}")
                print(f"\t... function index: {joint.iFunct}")
                print(f"\t... length: {joint.L}")
                print(f"\t... radius: {joint.R}")
                print(f"\t... initial condition x: {joint.x0}")
                print(f"\t... initial condition d: {', '.join(map(str, joint.d0))}")
                print(f"\t... fix: {joint.fix}")
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

    def __update_position(self):
        """
        Update position entities.
        """
        # update rotation matrix of the body
        nB = len(self.Bodies)
        for Bi in range(nB):
            body = self.Bodies[Bi]
            body._A = A_matrix(body.p)

        # compute sP = A * sP_prime; rP = r + sP
        nP = len(self.Points)
        for Pi in range(nP):
            point = self.Points[Pi]
            Bi = point.Bindex
            if Bi != 0:
                body = self.Bodies[Bi-1]
                point._sP = body._A @ point.sPlocal
                point._sPr = s_rot(point._sP)
                point._rP = body.r + point._sP

        # compute u = _A * up
        nU = len(self.uVectors)
        for Vi in range(nU):
            unit_vector = self.uVectors[Vi]
            Bi = unit_vector.Bindex
            if Bi != 0:
                body = self.Bodies[Bi - 1]
                unit_vector._u = body._A @ unit_vector.ulocal
                unit_vector._ur = s_rot(unit_vector._u)

    def __update_velocity(self):
        """
        Compute sP_dot and rP_dot vectors and update velocity components.
        """
        for Pi, point in enumerate(self.Points):
            Bi = point.Bindex
            if Bi != 0:
                point._dsP = point._sPr * self.Bodies[Bi-1].dp
                point._drP = self.Bodies[Bi-1].dr + point._dsP

        for Vi, uvector in enumerate(self.uVectors):
            Bi = uvector.Bindex
            if Bi != 0:
                uvector._du = uvector._ur * self.Bodies[Bi-1].dp
            
    def __compute_constraints(self):
        if not self.Joints:
            return np.zeros((0, 1))
        nConst = self.Joints[-1]._rowe
        phi = np.zeros([nConst, 1])

        for joint in self.Joints:
            match joint.type:
                case 'rev':
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    iBindex = joint.iBindex
                    jBindex = joint.jBindex
                    
                    # compute relative positions of the points
                    rPi = self.Points[Pi]._rP
                    rPj = self.Points[Pj]._rP
                    f = rPi - rPj

                    if joint.fix == 1:
                        if iBindex == 0:  # body i is ground (fixed)
                            f = np.append(f, (-self.Bodies[jBindex].p - joint._p0))
                        elif jBindex == 0:  # body j is ground (fixed)
                            f = np.append(f, (self.Bodies[iBindex].p - joint._p0))
                        else:
                            f = np.append(f, (self.Bodies[iBindex].p - self.Bodies[jBindex].p - joint._p0))

                case 'tran':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    ujr = self.uVectors[joint.jUindex]._ur
                    ui = self.uVectors[joint.iUindex]._u
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP

                    # compute constraint equations
                    f = np.array([ujr.T @ d, ujr.T @ ui]).reshape(2,1)

                    # additional constraint if fixed
                    if joint.fix == 1:
                        f = np.append(f, (ui.T @ d - joint._p0) / 2).reshape(3,1)

                case 'rev-rev':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    L = joint.L
                    u = d/L
                    # compute constraint equations
                    f = (u.T @ d - L)/2

                case 'rev-tran':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    uir = self.uVectors[joint.iUindex]._ur
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP

                    # compute constraint equations
                    f = (uir.T @ d - joint.L)

                case 'rigid':
                    Bi = joint.iBindex
                    Bj = joint.jBindex
                    
                    if Bi == 0:
                        f = np.vstack([
                            -(self.Bodies[Bj].r + self.Bodies[Bj]._A @ joint.d0),
                            -self.Bodies[Bj].p - joint._p0
                        ])
                    elif Bj == 0:
                        f = np.vstack([
                            self.Bodies[Bi].r - joint.d0,
                            self.Bodies[Bi].p - joint._p0
                        ])
                    else:
                        f = np.vstack([
                            self.Bodies[Bi].r - (self.Bodies[Bj].r + self.Bodies[Bj]._A @ joint.d0),
                            self.Bodies[Bi].p - self.Bodies[Bj].p - joint._p0
                        ])

                case 'disc':
                    Bi = joint.iBindex
                    f = np.vstack([
                        self.Bodies[Bi - 1].r[1] - joint.R,
                        (self.Bodies[Bi - 1].r[0] - joint.x0) + joint.R * (self.Bodies[Bi - 1].p - joint._p0)
                    ])
                    
                case 'rel-rot':
                    fun, fun_d, fun_dd = functEval(self.Functs[joint.iFunct - 1], self.t)
                    Bi = joint.iBindex
                    Bj = joint.jBindex

                    if Bi == 0:
                        f = -self.Bodies[Bj - 1].p - fun
                    elif Bj == 0:
                        f = self.Bodies[Bi - 1].p - fun
                    else:
                        f = self.Bodies[Bi - 1].p - self.Bodies[Bj - 1].p - fun

                case 'rel-tran':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    fun, fun_d, fun_dd = functEval(self.Functs[joint.iFunct - 1], self.t)

                    f = (d.T @ d - fun**2)/2

                case _:
                    raise ValueError(f"Joint type '{joint.type}' is not supported.")

            rs = joint._rows
            re = joint._rowe
            phi[rs:re] = f
            
        return phi

    def __compute_jacobian(self):
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
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    Di = np.block([
                        [np.eye(2), self.Points[Pi]._sPr.reshape(2, 1)]
                    ])
                    Dj = np.block([
                        [-np.eye(2), -self.Points[Pj]._sPr.reshape(2, 1)]
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
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    uj = self.uVectors[joint.jUindex]._u
                    ujr = self.uVectors[joint.jUindex]._ur
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP

                    Di = np.block([
                        [ujr.T, (uj.T @ self.Points[Pi]._sP).reshape(1, 1)],
                        [np.array([0, 0, 1])]
                    ])
                    Dj = np.block([
                        [-ujr.T, -(uj.T @ (self.Points[Pj]._sP + d)).reshape(1, 1)],
                        [np.array([0, 0, -1])]
                    ])

                    if joint.fix == 1:
                        Di = np.vstack([
                            Di,
                            [uj.T, (uj.T @ self.Points[Pi]._sPr).reshape(1)]
                        ])
                        Dj = np.vstack([
                            Dj,
                            [-uj.T, -(uj.T @ self.Points[Pj]._sPr).reshape(1)]
                        ])

                case 'rev-rev':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    L = joint.L
                    u = d/L

                    Di = np.block([
                        u.T, (u.T @ self.Points[Pi]._sPr).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -u.T, -(u.T @ self.Points[Pj]._sPr).reshape(1, 1)
                        ])
                    
                case 'rev-tran':
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    ui = self.uVectors[joint.iUindex]._u
                    ui_r = self.uVectors[joint.iUindex]._ur
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP

                    Di = np.block([
                        ui_r.T, (ui.T @ (self.Points[Pi]._sP - d)).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -ui_r.T, -(ui.T @ self.Points[Pj]._sP).reshape(1, 1)
                        ])

                case 'rigid':
                    Bj = joint.jBindex

                    Di = np.eye(3)
                    if Bj != 0:
                        Dj = np.block([
                            [-np.eye(2), -s_rot(self.Bodies[Bj]._A @ joint.d0)],
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
                    Pi = joint.iPindex
                    Pj = joint.jPindex

                    d = self.Points[Pi]._rP - self.Points[Pj]._rP

                    Di = np.block([
                        d.T, (d.T @ self.Points[Pi]._sPr).reshape(1, 1)
                        ])
                    Dj = np.block([
                        -d.T, -(d.T @ self.Points[Pj]._sPr).reshape(1, 1)
                        ])
                
                case _:
                    raise ValueError(f"Joint type '{joint.type}' is not supported.")

            # row indices for the current joint in the Jacobian matrix
            rs = joint._rows
            re = joint._rowe

            # column indices for body i 
            if joint.iBindex != 0:
                cis = joint._colis
                cie = joint._colie
                D[rs:re, cis:cie] = Di

            # column indices for body j
            if joint.jBindex != 0:
                cjs = joint._coljs
                cje = joint._colje
                D[rs:re, cjs:cje] = Dj

        return D

    def __rhs_velocity(self):
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
                    fun, fun_d, _ = functEval(self.Functs[joint.iFunct - 1], self.t)
                    f = fun_d

                case 'rel-tran':
                    fun, fun_d, _ = functEval(self.Functs[joint.iFunct - 1], self.t)
                    d = self.Points[joint.iPindex]._rP - self.Points[joint.jPindex]._rP
                    f = fun * fun_d

                case _:
                    continue

            rs = joint._rows
            re = joint._rowe
            rhsv[rs:re] = f

        return rhsv

    def __rhs_acceleration(self): #// - Check required on rel-tran and rel-rot joints -
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
                    Pi, Pj = joint.iPindex, joint.jPindex
                    Bi, Bj = self.Points[Pi].Bindex, self.Points[Pj].Bindex

                    if Bi == 0:
                        f = s_rot(self.Points[Pj]._dsP) * self.Bodies[Bj - 1].dp
                    elif Bj == 0:
                        f = -s_rot(self.Points[Pi]._dsP) * self.Bodies[Bi - 1].dp
                    else:
                        f = (
                            -s_rot(self.Points[Pi]._dsP) * self.Bodies[Bi - 1].dp
                            + s_rot(self.Points[Pj]._dsP) * self.Bodies[Bj - 1].dp
                        )

                    if joint.fix == 1:
                        f = np.vstack([f, [0]])
                
                case "tran":
                    Bi, Bj = joint.iBindex, joint.jBindex
                    Pi, Pj = joint.iPindex, joint.jPindex
                    ujd = self.uVectors[joint.jUindex]._du
                    ujdr = s_rot(ujd)

                    if Bi == 0:
                        f2 = 0.0
                    elif Bj == 0:
                        f2 = 0.0
                    else:
                        diffr = self.Bodies[Bi-1].r - self.Bodies[Bj-1].r
                        dp_product = (ujd.T @ diffr).item() * self.Bodies[Bi-1].dp
                        diffdr = self.Bodies[Bi-1].dr - self.Bodies[Bj-1].dr
                        f2 = dp_product - 2.0 * (ujdr.T @ diffdr).item()

                    f = np.array([[f2], [0.0]])

                    if joint.fix == 1:
                        d = self.Points[Pi]._rP - self.Points[Pj]._rP
                        dd = self.Points[Pi]._drP - self.Points[Pj]._drP
                        L = joint._p0 
                        u = d / L
                        du = dd / L
                        f3 = -(du.T @ dd).item()

                        if Bi == 0:
                            f3 += (u.T @ (s_rot(self.Points[Pj]._dsP) * self.Bodies[Bj-1].dp)).item()
                        elif Bj == 0:
                            f3 -= (u.T @ (s_rot(self.Points[Pi]._dsP) * self.Bodies[Bi-1].dp)).item()
                        else:
                            term1 = self.Points[Pi]._dsP * self.Bodies[Bi-1].dp
                            term2 = self.Points[Pj]._dsP * self.Bodies[Bj-1].dp
                            f3 -= (u.T @ s_rot(term1 - term2)).item()

                        f = np.vstack([f, [[f3]]])
                    
                case "rev-rev":
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    Bi = joint.iBindex
                    Bj = joint.jBindex
                    
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    dd = self.Points[Pi]._drP - self.Points[Pj]._drP
                    
                    L = joint.L
                    u = d/L
                    ud = dd/L
                    
                    f = -ud.T @ dd
                    
                    if Bi == 0:
                        f = f + u.T @ s_rot(self.Points[Pj]._dsP) * self.Bodies[Bj-1].dp
                    elif Bj == 0:
                        f = f - u.T @ s_rot(self.Points[Pi]._dsP) * self.Bodies[Bi-1].dp
                    else:
                        f = f - u.T @ s_rot(
                            self.Points[Pi]._dsP * self.Bodies[Bi-1].dp - 
                            self.Points[Pj]._dsP * self.Bodies[Bj-1].dp
                    )

                case "rev-tran":
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    Bi = joint.iBindex
                    Bj = joint.jBindex

                    ui = self.uVectors[joint.iUindex]._u
                    ui_d = self.uVectors[joint.iUindex]._du
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    dd = self.Points[Pi]._drP - self.Points[Pj]._drP

                    if Bi == 0:
                        f = ui.T @ self.Points[Pj]._dsP * self.Bodies[Bj-1].dp
                    elif Bj == 0:
                        f = ui_d.T @ (d * self.Bodies[Bi-1].dp + 2 * s_rot(dd)) - \
                            ui.T @ self.Points[Pi]._dsP * self.Bodies[Bi-1].dp
                    else:
                        f = ui_d.T @ (d * self.Bodies[Bi-1].dp + 2 * s_rot(dd)) - \
                            ui.T @ (self.Points[Pi]._dsP * self.Bodies[Bi-1].dp - \
                                self.Points[Pj]._dsP * self.Bodies[Bj-1].dp)
                
                case "rigid":
                    Bj = joint.jBindex

                    f = np.zeros(3)
                    if Bj != 0:
                        f = np.concatenate([
                            -self.Bodies[Bj-1]._A @ joint.d0 * self.Bodies[Bj-1].dp**2,
                            np.array([0])
                        ])
                
                case "disc":
                    f = np.zeros(2)
                    
                case "rel-rot":
                    fun, fun_d, fun_dd = functEval(self.Functs[joint.iFunct - 1], self.t)
                    f = fun_dd

                case "rel-tran":
                    Pi = joint.iPindex
                    Pj = joint.jPindex
                    Bi = joint.iBindex
                    Bj = joint.jBindex

                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    dd = self.Points[Pi]._drP - self.Points[Pj]._drP

                    fun, fun_d, fun_dd = functEval(self.Functs[joint.iFunct - 1], self.t)

                    f = fun * fun_dd + fun_d**2

                    if Bi == 0:
                        f = f + d.T @ s_rot(self.Points[Pj]._dsP).T @ self.Bodies[Bj - 1].dp
                    elif Bj == 0:
                        f = f - d.T @ s_rot(self.Points[Pi]._dsP).T @ self.Bodies[Bi - 1].dp - dd.T @ dd
                    else:
                        f = f + d.T @ s_rot(self.Points[Pj]._dsP).T @ self.Bodies[Bj - 1].dp \
                            - d.T @ s_rot(self.Points[Pi]._dsP).T @ self.Bodies[Bi - 1].dp - dd.T @ dd
            
            rs = joint._rows
            re = rs + joint._mrows 
            rhsa[rs:re] = np.asarray(f).reshape(-1, 1)
        
        return rhsa

    def __bodies2u(self):
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

    def __bodies2ud(self):
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
    
    def __u2bodies(self, u):
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

    def __compute_force(self):
        """
        Compute and return the array of forces acting on the system at time t.
        """
        for body in self.Bodies:
            body._f = colvect([0.0, 0.0]) # initialize body force vectors
            body._n = 0.0                 # initialize body torque (moment) scalar

        #! the below code used to build the force array need to be optimized, 
        #! the class Force should keep inside all the method required to build
        #! the specific type of force -> code easier and more readable.
        # loop over all forces and apply them to the appropriate bodies
        for Fi, force in enumerate(self.Forces):
            match force.type:
                case 'weight':
                    for body in self.Bodies:
                        body._f += body._wgt

                case 'ptp':
                    Pi, Pj = force.iPindex, force.jPindex
                    Bi, Bj = force.iBindex, force.jBindex
                    d = self.Points[Pi]._rP - self.Points[Pj]._rP
                    dd = self.Points[Pi]._drP - self.Points[Pj]._drP
                    L = np.linalg.norm(d)
                    dL = d.T @ dd / L
                    delta = L - self.Forces[Fi].L0
                    u = d / L

                    f = self.Forces[Fi].k * delta + self.Forces[Fi].dc * dL + self.Forces[Fi].f_a
                    fi = f * u

                    if Bi != 0:
                        self.Bodies[Bi-1]._f -= fi
                        self.Bodies[Bi-1]._n -= (self.Points[Pi]._sPr.T @ fi).item()
                    
                    if Bj != 0:
                        self.Bodies[Bj-1]._f += fi
                        self.Bodies[Bj-1]._n += (self.Points[Pj]._sPr.T @ fi).item()

                case 'rot-sda':
                    Bi = force.iBindex
                    Bj = force.jBindex

                    if Bi == 0:
                        theta = -self.Bodies[Bj-1].p
                        theta_d = -self.Bodies[Bj-1].dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        self.Bodies[Bj-1]._n += T
                    elif Bj == 0:
                        theta = self.Bodies[Bi-1].p
                        theta_d = self.Bodies[Bi-1].dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        self.Bodies[Bi-1]._n -= T
                    else:
                        theta = self.Bodies[Bi-1].p - self.Bodies[Bj-1].p
                        theta_d = self.Bodies[Bi-1].dp - self.Bodies[Bj-1].dp
                        T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                        self.Bodies[Bi-1]._n -= T
                        self.Bodies[Bj-1]._n += T

                case 'flocal':
                    Bi = force.iBindex
                    self.Bodies[Bi-1]._f += self.Bodies[Bi-1]._A @ force.flocal

                case 'f':
                    Bi = force.iBindex
                    self.Bodies[Bi-1]._f += force.f

                case 'T':
                    Bi = force.iBindex
                    self.Bodies[Bi-1]._n += force.T

                case 'user' if force.callback is not None and callable(force.callback):
                    global_vars = get_globals()
                    params = list(inspect.signature(force.callback).parameters.keys())
                    args = [global_vars[name] for name in params if name in global_vars]
                    force.callback(*args)

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

    def __ic_correct(self):
        """
        Corrects initial conditions on the body coordinates and velocities.
        """
        flag = False

        # position correction
        for _ in range(50): #! 20 is an arbitrary value ... could be a parameter!
            self.__update_position()            # update position entities
            Phi = self.__compute_constraints()  # evaluate constraints
            D = self.__compute_jacobian()       # evaluate Jacobian
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
                self.Bodies[Bi].p = self.Bodies[Bi].p + delta_c[ir + 2][0] # [0] because I need to extract the single value

        if not flag:
            raise ValueError("Convergence failed in Newton-Raphson!")

        # velocity correction
        nB = len(self.Bodies)
        Phi = np.zeros([3 * nB, 1])
        for Bi in range(nB):
            ir = 3 * Bi
            Phi[ir:ir + 2] = self.Bodies[Bi].dr
            Phi[ir + 2] = self.Bodies[Bi].dp

        rhsv = self.__rhs_velocity()
        
        # solve for corrections
        delta_v = -D.T @ np.linalg.solve(D @ D.T, D @ Phi - rhsv)  

        # move corrected velocities to sub-arrays
        for Bi in range(nB):
            ir = 3 * Bi
            self.Bodies[Bi].dr = self.Bodies[Bi].dr + delta_v[ir:ir + 2]
            self.Bodies[Bi].dp = self.Bodies[Bi].dp + delta_v[ir + 2][0] # [0] because I need to extract the single value

        coords = np.zeros((nB, 3))
        vels = np.zeros((nB, 3))
        for Bi in range(nB):
            coords[Bi, :] = np.hstack((self.Bodies[Bi].r.T, np.array(self.Bodies[Bi].p).reshape(-1, 1)))
            vels[Bi, :] = np.hstack((self.Bodies[Bi].dr.T, np.array(self.Bodies[Bi].dp).reshape(-1, 1)))

        # // ---
        # // ... some check if required ...
        # // ---
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
                print(f"\t... points: {', '.join(map(str, body._pts))}")
                print(f"\t... irc: {body._irc}")
                print(f"\t... irv: {body._irv}")
                print(f"\t... ira: {body._ira}")
                print("\n")
            
            print(f"-----")
            print(f"points ")
            print(f"-----")
            for i, point in enumerate(self.Points):
                print(f"\t... point: {i}")
                print(f"\t... body index: {point.Bindex}")
                print(f"\t... local coordinates: {', '.join(map(str, point.sPlocal.flatten()))}")
                print(f"\t... global coordinates: {', '.join(map(str, point._rP.flatten()))}")
                print(f"\t... sP: {', '.join(map(str, point._sP.flatten()))}")
                print(f"\t... sPr: {', '.join(map(str, point._sPr.flatten()))}")
                print(f"\t... dsP: {', '.join(map(str, point._dsP.flatten()))}")
                print(f"\t... drP: {', '.join(map(str, point._drP.flatten()))}")
                print(f"\t... ddrP: {', '.join(map(str, point._ddrP.flatten()))}")
                print("\n")
            
            print(f"-----")
            print(f"vectors ")
            print(f"-----")
            for i, uvector in enumerate(self.uVectors, start=1):
                print(f"\t... uVector: {i}")
                print(f"\t... body index: {uvector.Bindex}")
                print(f"\t... local vector: {', '.join(map(str, uvector.ulocal.flatten()))}")
                print(f"\t... global vector: {', '.join(map(str, uvector._u.flatten()))}")
                print(f"\t... ur: {', '.join(map(str, uvector._ur.flatten()))}")
                print(f"\t... du: {', '.join(map(str, uvector._du.flatten()))}")
                print("\n")
            
            print(f"-----")
            print(f"forces ")
            print(f"-----")
            for i, force in enumerate(self.Forces, start=1):
                print(f"\t... force: {i}")
                print(f"\t... type: {force.type}")
                print(f"\t... head point index: {force.iPindex}")
                print(f"\t... tail point index: {force.jPindex}")
                print(f"\t... head body index: {force.iBindex}")
                print(f"\t... tail body index: {force.jBindex}")
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
                print(f"\t... body i index: {joint.iBindex}")
                print(f"\t... body j index: {joint.jBindex}")
                print(f"\t... point i index: {joint.iPindex}")
                print(f"\t... point j index: {joint.jPindex}")
                print(f"\t... unit vector i index: {joint.iUindex}")
                print(f"\t... unit vector j index: {joint.jUindex}")
                print(f"\t... function index: {joint.iFunct}")
                print(f"\t... length: {joint.L}")
                print(f"\t... radius: {joint.R}")
                print(f"\t... initial condition x: {joint.x0}")
                print(f"\t... initial condition d: {', '.join(map(str, joint.d0))}")
                print(f"\t... fix: {joint.fix}")
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

    def __analysis(self, t, u):
        """
        Solve the constrained equations of motion at time t with the standard
        Lagrange multiplier method.
        """        
        self.__num += 1                 # increment the number of function evaluations
        self.t = t                      # store current time for force/constraint callbacks
        nB3 = 3 * len(self.Bodies)
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        self.__u2bodies(u)              # unpack u into coordinate and velocity sub-arrays
        self.__update_position()
        self.__update_velocity()
        h_a = self.__compute_force()    # array of applied forces

        if nConst == 0:
            ddc = self.invM_array.reshape(-1, 1) * h_a  # solve for accelerations
            Lambda = np.array([])  # no constraints, no multipliers
        else:
            D = self.__compute_jacobian()
            rhsA = self.__rhs_acceleration()  # right-hand side of acceleration constraints (gamma)

            # construct the matrix system to solve
            DMD = np.block([
                [self.M_matrix, -D.T], 
                [D, np.zeros([nConst, nConst])]
            ])
            rhs = np.concatenate([h_a, rhsA])

            #* check on conditioned index of the coefficient matrix
            cond_number = np.linalg.cond(DMD)
            if cond_number > 1e12:
                print(f"Warning: DMD matrix is poorly conditioned with condition number {cond_number}")
    
            # solve the system of equations
            sol = np.linalg.solve(DMD, rhs)
            ddc = sol[:nB3]
            Lambda = sol[nB3:]

        # update accelerations for each body
        for Bi, body in enumerate(self.Bodies):
            ir = body._irc
            i2 = ir + 2
            i3 = i2
            body.ddr = ddc[ir:i2]
            body.ddp = ddc[i3][0]

        ud = self.__bodies2ud()             # pack velocities and accelerations into ud
        return ud.flatten()

    def __taqaddum(self, t_initial, t_final, pbar):
        """
        Restituisce una funzione wrapper per __analysis con progresso ottimizzato
        """
        last_progress = 0  # Ora è un intero invece di un dizionario
        
        def __wrapp_analysis(t, u):
            nonlocal last_progress
            progress = min(100, int(100 * (t - t_initial) / (t_final - t_initial)))
            if progress > last_progress:
                pbar.n = progress
                pbar.refresh()
                last_progress = progress
            return self.__analysis(t, u)
        return __wrapp_analysis

    def solve(self, method="LSODA", t_final=None, dt=None, ic_correct=False,
              t_eval=None, t_span=None):
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

        Returns
        -------
        SolResult
            Object supporting tuple unpacking (``T, uT = sol``) and
            attribute access (``sol.t``, ``sol.y`` in scipy convention).
        """
        self.method = method

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
            self.t = 0.0  # needed by __compute_constraints / __rhs_velocity when ic_correct runs
            if ans == 'y':
                self.__ic_correct()
            D = self.__compute_jacobian()
            redund = np.linalg.matrix_rank(D)
            if redund < nConst:
                print("\n\t...Redundancy in the constraints")

        u = self.__bodies2u()
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
        self.__num = 0  # initialize the number of function evaluations

        if t_final is None:
            t_final = float(input("\n\t ...Final time = ? "))

        dense_sol = None
        if t_final == 0:
            self.__analysis(0, u)
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

            __wrapp_analysis = self.__taqaddum(t_initial, t_final, pbar)

            try:
                _sol = solve_ivp(__wrapp_analysis,
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
        print(f"\t ...Number of function evaluations: {self.__num}")
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
            self.__u2bodies(u_i)
            self.__update_position()
            self.__update_velocity()
            h_a = self.__compute_force()

            if nConst == 0:
                ddc = self.invM_array.reshape(-1, 1) * h_a
                Lambda = np.array([])
            else:
                D = self.__compute_jacobian()
                rhsA = self.__rhs_acceleration()

                DMD = np.block([
                    [self.M_matrix, -D.T],
                    [D, np.zeros((nConst, nConst))]
                ])
                rhs = np.concatenate([h_a, rhsA])
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