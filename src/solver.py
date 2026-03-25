"""Planar Multi-Body Dynamics Simulation Solver.

This module provides algorithms and numerical methods for solving
planar multi-body dynamic models, including rigid bodies, constraints,
and external forces.

INDEX CONVENTION:
================
Body indices:  0 = ground (fixed), 1..nB = moving bodies
Point indices: 0-based (Python standard)
Joint indices: 0-based (Python standard)

Internal state vector u:
  u[0:nB3]      = positions  [r1x, r1y, p1, r2x, r2y, p2, ...]
  u[nB3:2*nB3]  = velocities [dr1x, dr1y, dp1, ...]

Internal index attributes (0-based):
  body._index_position    = 3*Bi       (start index for position in u)
  body._index_velocity    = nB3 + 3*Bi (start index for velocity in u)
  joint._rows             = constraint row start (0-based, inclusive)
  joint._rowe             = constraint row end (0-based, exclusive)
  joint._colis/_coljs     = Jacobian column start (0-based)
  joint._colie/_colje     = Jacobian column end (0-based, exclusive)

All slicing uses standard Python [start:end) convention directly.

Author: Giacomo Cangi
"""

import logging
import os

import numpy as np
import numpy.linalg as lng
import scipy as sc
from scipy.integrate import solve_ivp
from tqdm import tqdm

from .constraints import Weight
from .mechanics import *
from .model import Ground
from .utils import *

logger = logging.getLogger(__name__)



class PlanarMultibodyModel:
    """Planar multi-body dynamics model and solver.

    Attributes:
        Bodies: List of Body objects.
        Joints: List of Joint objects.
        Forces: List of Force objects.
        Functs: List of Function objects.
    """

    def __init__(self, bodies, joints=None, forces=None, functions=None):
        self.Bodies = list(bodies)
        self.Joints = list(joints) if joints else []
        self.Forces = list(forces) if forces else []
        self.Functs = list(functions) if functions else []

        # Auto-assembly trigger
        from .builder import _assemble
        needs_assembly = any(
            not b._position_given and not b._orientation_given
            for b in self.Bodies
        )
        if needs_assembly:
            _assemble(self.Bodies, self.Joints)

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
        """Initialize the multi-body model after construction."""
        nB = len(self.Bodies)
        nB3 = 3 * nB

        # Bodies — assign internal indices and compute derived quantities
        for Bi, body in enumerate(self.Bodies):
            body._body_index = Bi + 1
            body._index_position = 3 * Bi
            body._index_velocity = nB3 + 3 * Bi
            body._index_acceleration = nB3 + 3 * Bi  # same offset in accel vector
            body._inv_mass = 1 / body.mass
            body._inv_inertia = 1 / body.inertia
            body._rotation_matrix = rotation_matrix(body.orientation)

        # Mass (inertia) array and pre-computed diagonal matrix
        self.M_array = np.zeros(nB3)
        self.invM_array = np.zeros(nB3)
        for Bi, body in enumerate(self.Bodies):
            is_ = 3 * Bi
            ie_ = is_ + 3
            self.M_array[is_:ie_] = [body.mass, body.mass, body.inertia]
            self.invM_array[is_:ie_] = [body._inv_mass, body._inv_mass, body._inv_inertia]
        self.M_matrix = np.diag(self.M_array)

        # Ground markers
        for marker in Ground._markers:
            pos_col = marker.local_position.reshape(2, 1)
            marker._sP = pos_col.copy()
            marker._sPr = rotate_90(pos_col)
            marker._rP = pos_col.copy()
            if marker.has_orientation:
                marker._u = marker._ulocal.copy()
                marker._ur = rotate_90(marker._ulocal)

        # Body markers
        for body in self.Bodies:
            for marker in body._markers:
                pos_col = marker.local_position.reshape(2, 1)
                marker._sP = body._rotation_matrix @ pos_col
                marker._sPr = rotate_90(marker._sP)
                marker._rP = body.position + marker._sP
                if marker.has_orientation:
                    marker._u = body._rotation_matrix @ marker._ulocal
                    marker._ur = rotate_90(marker._u)

        # Force elements — initialize weights
        for force in self.Forces:
            if isinstance(force, Weight):
                force.initialize_weights(self.Bodies)

        # Joints — polymorphic initialization
        for joint in self.Joints:
            joint.initialize(self)

        # Validation: all joint bodies must be in self.Bodies or Ground
        body_ids = {id(b) for b in self.Bodies}
        body_ids.add(id(Ground))
        for joint in self.Joints:
            if id(joint.iBody) not in body_ids:
                raise ValueError(
                    f"Joint iBody {joint.iBody} is not in the model's bodies list")
            if id(joint.jBody) not in body_ids:
                raise ValueError(
                    f"Joint jBody {joint.jBody} is not in the model's bodies list")

        # Functions
        if self.Functs:
            for Ci in range(len(self.Functs)):
                functData(Ci, self.Functs)

        # Compute constraint row/column pointers
        nConst = 0
        for joint in self.Joints:
            joint._rows = nConst
            joint._rowe = nConst + joint._mrows
            nConst = joint._rowe
            Bi_idx = joint.iBody._body_index
            if Bi_idx != 0:
                joint._colis = 3 * (Bi_idx - 1)
                joint._colie = 3 * Bi_idx
            Bj_idx = joint.jBody._body_index
            if Bj_idx != 0:
                joint._coljs = 3 * (Bj_idx - 1)
                joint._colje = 3 * Bj_idx

        logger.info("Model initialized: %d bodies, %d joints, %d forces",
                     nB, len(self.Joints), len(self.Forces))

    def _update_position(self):
        """Update position-dependent kinematic quantities."""
        for body in self.Bodies:
            body._rotation_matrix = rotation_matrix(body.orientation)

        for body in self.Bodies:
            for marker in body._markers:
                pos_col = marker.local_position.reshape(2, 1)
                marker._sP = body._rotation_matrix @ pos_col
                marker._sPr = rotate_90(marker._sP)
                marker._rP = body.position + marker._sP
                if marker.has_orientation:
                    marker._u = body._rotation_matrix @ marker._ulocal
                    marker._ur = rotate_90(marker._u)

    def _update_velocity(self):
        """Update velocity-dependent kinematic quantities."""
        for body in self.Bodies:
            for marker in body._markers:
                marker._dsP = marker._sPr * body.angular_velocity
                marker._drP = body.velocity + marker._dsP
                if marker.has_orientation:
                    marker._du = marker._ur * body.angular_velocity

    def _compute_constraints(self):
        """Evaluate all constraint equations.

        Returns:
            ndarray of shape (nConst, 1), or (0, 1) if no joints.
        """
        if not self.Joints:
            return np.zeros((0, 1))
        nConst = self.Joints[-1]._rowe
        phi = np.zeros([nConst, 1])

        for joint in self.Joints:
            f = joint.compute_phi(self)
            rs = joint._rows
            re = joint._rowe
            phi[rs:re] = np.asarray(f).reshape(-1, 1)

        return phi

    def _compute_jacobian(self):
        """Compute the Jacobian matrix for all constraints.

        Returns:
            ndarray of shape (nConst, nB3).
        """
        nConst = self.Joints[-1]._rowe
        nB3 = 3 * len(self.Bodies)
        D = np.zeros((nConst, nB3))

        for joint in self.Joints:
            rs = joint._rows
            re = joint._rowe

            if joint.iBody is not Ground:
                Di = joint.compute_jacobian_i(self)
                cis = joint._colis
                cie = joint._colie
                D[rs:re, cis:cie] = Di

            if joint.jBody is not Ground:
                Dj = joint.compute_jacobian_j(self)
                cjs = joint._coljs
                cje = joint._colje
                D[rs:re, cjs:cje] = Dj

        return D

    def _rhs_velocity(self):
        """Compute the right-hand side of velocity constraints.

        Returns:
            ndarray of shape (nConst, 1).
        """
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        rhsv = np.zeros((nConst, 1))

        for joint in self.Joints:
            f = joint.compute_rhs_velocity(self)
            if f is not None:
                rs = joint._rows
                re = joint._rowe
                rhsv[rs:re] = np.asarray(f).reshape(-1, 1)

        return rhsv

    def _rhs_acceleration(self):
        """Compute the right-hand side of acceleration constraints (gamma).

        Returns:
            ndarray of shape (nConst, 1).
        """
        nConst = self.Joints[-1]._rowe
        rhsa = np.zeros([nConst, 1])

        for joint in self.Joints:
            f = joint.compute_rhs_acceleration(self)
            rs = joint._rows
            re = rs + joint._mrows
            rhsa[rs:re] = np.asarray(f).reshape(-1, 1)

        return rhsa

    def _bodies2u(self):
        """Pack coordinates and velocities into the state vector u."""
        nB = len(self.Bodies)
        u = np.zeros([(3 * nB * 2), 1])

        for Bi in range(nB):
            ir = self.Bodies[Bi]._index_position
            ird = self.Bodies[Bi]._index_velocity
            u[ir:ir+3] = np.block([[self.Bodies[Bi].position],
                                    [self.Bodies[Bi].orientation]])
            u[ird:ird+3] = np.block([[self.Bodies[Bi].velocity],
                                      [self.Bodies[Bi].angular_velocity]])

        return u

    def _bodies2ud(self):
        """Pack velocities and accelerations into ud."""
        nB6 = 6 * len(self.Bodies)
        ud = np.zeros([nB6, 1])

        for Bi, body in enumerate(self.Bodies):
            ir = body._index_position
            ird = body._index_velocity
            ud[ir:ir + 3] = np.vstack([body.velocity,
                                        body.angular_velocity]).reshape(3, 1)
            ud[ird:ird + 3] = np.vstack([body.acceleration,
                                          body.angular_acceleration]).reshape(3, 1)

        return ud

    def _u2bodies(self, u):
        """Unpack state vector u into body attributes."""
        if u.ndim != 2:
            u = u.reshape(-1, 1)

        nB = len(self.Bodies)
        for Bi in range(nB):
            ir = self.Bodies[Bi]._index_position
            ird = self.Bodies[Bi]._index_velocity
            self.Bodies[Bi].position = u[ir:ir+2]
            self.Bodies[Bi].orientation = u[ir+2][0]
            self.Bodies[Bi].velocity = u[ird:ird+2]
            self.Bodies[Bi].angular_velocity = u[ird+2][0]

    def _compute_force(self):
        """Compute the array of forces acting on the system.

        Returns:
            ndarray of shape (nB3, 1).
        """
        for body in self.Bodies:
            body._force = colvect([0.0, 0.0])
            body._torque = 0.0

        # Polymorphic force application
        for force in self.Forces:
            force.apply(self.Bodies)

        # Build force array
        nB3 = 3 * len(self.Bodies)
        g = np.zeros([nB3, 1])
        for Bi, body in enumerate(self.Bodies):
            ks = body._index_position
            ke = ks + 3
            g[ks:ke] = np.vstack([body._force, body._torque])

        return g

    def _ic_correct(self):
        """Correct initial conditions on body coordinates and velocities."""
        flag = False

        # Position correction
        for _ in range(50):
            self._update_position()
            Phi = self._compute_constraints()
            D = self._compute_jacobian()
            ff = np.sqrt(Phi.T @ Phi)

            if ff < 1.0e-10:
                flag = True
                break

            delta_c = -D.T @ np.linalg.solve(D @ D.T, Phi)

            nB = len(self.Bodies)
            for Bi in range(nB):
                ir = 3 * Bi
                self.Bodies[Bi].position = self.Bodies[Bi].position + delta_c[ir:ir + 2]
                self.Bodies[Bi].orientation = self.Bodies[Bi].orientation + delta_c[ir + 2][0]

        if not flag:
            raise ValueError("Convergence failed in Newton-Raphson!")

        # Velocity correction
        nB = len(self.Bodies)
        Phi = np.zeros([3 * nB, 1])
        for Bi in range(nB):
            ir = 3 * Bi
            Phi[ir:ir + 2] = self.Bodies[Bi].velocity
            Phi[ir + 2] = self.Bodies[Bi].angular_velocity

        rhsv = self._rhs_velocity()
        delta_v = -D.T @ np.linalg.solve(D @ D.T, D @ Phi - rhsv)

        for Bi in range(nB):
            ir = 3 * Bi
            self.Bodies[Bi].velocity = self.Bodies[Bi].velocity + delta_v[ir:ir + 2]
            self.Bodies[Bi].angular_velocity = (
                self.Bodies[Bi].angular_velocity + delta_v[ir + 2][0])

        coords = np.zeros((nB, 3))
        vels = np.zeros((nB, 3))
        for Bi in range(nB):
            coords[Bi, :] = np.hstack((
                self.Bodies[Bi].position.T,
                np.array(self.Bodies[Bi].orientation).reshape(-1, 1)))
            vels[Bi, :] = np.hstack((
                self.Bodies[Bi].velocity.T,
                np.array(self.Bodies[Bi].angular_velocity).reshape(-1, 1)))

        logger.info("Initial conditions corrected")
        logger.debug("Corrected coordinates:\n%s", coords)
        logger.debug("Corrected velocities:\n%s", vels)

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
        """Solve constrained equations of motion at time t.

        Args:
            t: Current time.
            u: State vector.

        Returns:
            Flattened derivative vector ud.
        """
        self._num += 1
        self.t = t
        nB3 = 3 * len(self.Bodies)
        nConst = self.Joints[-1]._rowe if self.Joints else 0
        self._u2bodies(u)
        self._update_position()
        self._update_velocity()
        h_a = self._compute_force()

        if nConst == 0:
            ddc = self.invM_array.reshape(-1, 1) * h_a
            Lambda = np.array([])
        else:
            D = self._compute_jacobian()
            rhsA = self._rhs_acceleration()

            if self._ggl_mu > 0:
                reg_block = (1.0 / self._ggl_mu) * np.eye(nConst)
            else:
                reg_block = np.zeros((nConst, nConst))

            DMD = np.block([
                [self.M_matrix, -D.T],
                [D, reg_block]
            ])

            if self._baumgarte_alpha > 0 or self._baumgarte_beta > 0:
                Phi = np.asarray(self._compute_constraints()).flatten()
                dq = u[nB3:]
                dPhi = np.asarray(D @ dq).flatten()
                gamma_stab = (rhsA.flatten()
                              - 2.0 * self._baumgarte_alpha * dPhi
                              - self._baumgarte_beta**2 * Phi)
            else:
                gamma_stab = rhsA.flatten()

            rhs = np.concatenate([h_a.flatten(), gamma_stab])

            cond_number = np.linalg.cond(DMD)
            if cond_number > 1e12:
                logger.warning("DMD matrix poorly conditioned: cond = %e", cond_number)

            sol = np.linalg.solve(DMD, rhs).reshape(-1, 1)
            ddc = sol[:nB3]
            Lambda = sol[nB3:]

        for Bi, body in enumerate(self.Bodies):
            ir = body._index_position
            i2 = ir + 2
            body.acceleration = ddc[ir:i2]
            body.angular_acceleration = ddc[i2][0]

        ud = self._bodies2ud()
        return ud.flatten()

    def _taqaddum(self, t_initial, t_final, pbar):
        """Return wrapped analysis function with progress tracking."""
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

    def _distribute_results(self, T, uT):
        """Push simulation results into Body and Joint _result_container dicts.

        Args:
            T: Time vector, shape (nSteps,).
            uT: State matrix, shape (nSteps, 2*nB3).
        """
        nB3 = 3 * self.nB
        accelerations, reactions = self._post_process(T, uT)
        for Bi, body in enumerate(self.Bodies):
            i = 3 * Bi
            body._result_container = {
                'positions': {
                    'x': uT[:, i],
                    'y': uT[:, i + 1],
                    'phi': uT[:, i + 2],
                },
                'velocities': {
                    'dx': uT[:, nB3 + i],
                    'dy': uT[:, nB3 + i + 1],
                    'dphi': uT[:, nB3 + i + 2],
                },
                'accelerations': {
                    'ddx': accelerations[:, i],
                    'ddy': accelerations[:, i + 1],
                    'ddphi': accelerations[:, i + 2],
                },
            }
        for joint in self.Joints:
            rs = joint._rows
            re = joint._rowe
            joint._result_container = {'reactions': reactions[:, rs:re]}

    def solve(self, method="LSODA", t_final=None, dt=None, ic_correct=False,
              t_eval=None, t_span=None,
              baumgarte_alpha=5.0, baumgarte_beta=5.0, ggl_penalty=1e8):
        """Solve equations of motion.

        Args:
            method: ODE solver method (default "LSODA").
            t_final: Final simulation time. If None, prompts interactively.
            dt: Output time step. Used only when t_eval is not given.
            ic_correct: Whether to correct initial conditions before solving.
            t_eval: Explicit array of output time points.
            t_span: (t_start, t_end) shorthand.
            baumgarte_alpha: Velocity-level stabilization gain (default 5.0).
            baumgarte_beta: Position-level stabilization gain (default 5.0).
            ggl_penalty: GGL regularization parameter (default 1e8).

        Returns:
            Tuple (T, uT): time vector (nSteps,) and state matrix (nSteps, 2*nB3).
            Results are also distributed into each Body and Joint via
            ``_result_container``; use ``body.get_position()`` etc. to retrieve them.
        """
        self.method = method
        self._baumgarte_alpha = baumgarte_alpha
        self._baumgarte_beta = baumgarte_beta
        self._ggl_mu = ggl_penalty

        if t_span is not None and t_final is None:
            t_final = t_span[1]

        nConst = self.Joints[-1]._rowe if self.Joints else 0

        if t_final is None:
            print("\n")
            ans = input("\t... Do you want to correct the initial conditions? [(y)es/(n)o] ").lower()
        else:
            ans = 'y' if ic_correct else 'n'

        if nConst != 0:
            self.t = 0.0
            if ans == 'y':
                self._ic_correct()
            D = self._compute_jacobian()
            redund = np.linalg.matrix_rank(D)
            if redund < nConst:
                logger.warning("Redundancy in the constraints")
                print("\n\t...Redundancy in the constraints")

        u = self._bodies2u()
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            raise ValueError("Check initial conditions: u vector contains NaN or Inf values.")

        t_initial = 0.0
        self._num = 0

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
            finally:
                pbar.close()

            T = _sol.t
            uT = _sol.y.T
            dense_sol = _sol.sol

        print(f"\n ")
        print(f"\t ...Number of function evaluations: {self._num}")
        print(f"\t ...Simulation completed successfully!")
        print(f"\n ")
        self._distribute_results(T, uT)
        return T, uT

    def _post_process(self, T, uT):
        """Recalculate accelerations and Lagrange multipliers on t_eval grid.

        Args:
            T: Time vector, shape (nSteps,).
            uT: State matrix, shape (nSteps, 2*nB3).

        Returns:
            Tuple of (accelerations, reactions) arrays.
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

                if self._ggl_mu > 0:
                    reg_block = (1.0 / self._ggl_mu) * np.eye(nConst)
                else:
                    reg_block = np.zeros((nConst, nConst))

                DMD = np.block([
                    [self.M_matrix, -D.T],
                    [D, reg_block]
                ])

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

