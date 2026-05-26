"""Planar multi-body dynamics simulation solver.

Author: Giacomo Cangi
"""

# ---------------------------------------------------------------------------
# INDEX CONVENTION (internal)
# ---------------------------------------------------------------------------
# Body indices:  0 = ground (fixed), 1..nB = moving bodies
# Point indices: 0-based (Python standard)
# Joint indices: 0-based (Python standard)
#
# Internal state vector u:
#   u[0:nB3]      = positions  [r1x, r1y, p1, r2x, r2y, p2, ...]
#   u[nB3:2*nB3]  = velocities [dr1x, dr1y, dp1, ...]
#
# Internal index attributes (0-based):
#   body._index_position    = 3*Bi       (start index for position in u)
#   body._index_velocity    = nB3 + 3*Bi (start index for velocity in u)
#   joint._rows             = constraint row start (0-based, inclusive)
#   joint._rowe             = constraint row end (0-based, exclusive)
#   joint._colis/_coljs     = Jacobian column start (0-based)
#   joint._colie/_colje     = Jacobian column end (0-based, exclusive)
#
# All slicing uses standard Python [start:end) convention directly.
# ---------------------------------------------------------------------------

import logging
import os

import casadi as ca
import numpy as np
import numpy.linalg as lng
from tqdm import tqdm

from .forces import Weight
from .mechanics import *
from .model import Ground
from .utils import *

logger = logging.getLogger(__name__)



class PlanarMultibodyModel:
    """Planar multi-body dynamics model and solver.

    Attributes
    ----------
    Bodies : list of Body
        List of Body objects.
    Joints : list of Joint
        List of Joint objects.
    Forces : list of Force
        List of Force objects.
    Functs : list of Function
        List of Function objects.
    """

    def __init__(self, bodies, joints=None, forces=None, functions=None, units=None):
        from .units import UnitSystem
        self.units = units if units is not None else UnitSystem()
        self.Bodies = list(bodies)
        self.Joints = list(joints) if joints else []
        self.Forces = list(forces) if forces else []
        self.Functs = list(functions) if functions else []

        # Auto-assembly trigger
        from .builder import _assemble, _resolve_deferred_markers
        needs_assembly = any(
            not b._position_given and not b._orientation_given
            for b in self.Bodies
        )
        if needs_assembly:
            _assemble(self.Bodies, self.Joints)
        else:
            _resolve_deferred_markers(self.Bodies)

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

    @property
    def nDOF(self):
        """Grübler degrees of freedom: 3*nB - nC."""
        return 3 * len(self.Bodies) - self.nC

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

        # Actuator expansion: length-mode actuators register a TranMotion joint
        from .forces import Actuator
        for force in list(self.Forces):
            if isinstance(force, Actuator):
                force.expand(self)

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

        Returns
        -------
        ndarray
            Array of shape (nConst, 1), or (0, 1) if no joints.
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

        Returns
        -------
        ndarray
            Array of shape (nConst, nB3).
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

        Returns
        -------
        ndarray
            Array of shape (nConst, 1).
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

        Returns
        -------
        ndarray
            Array of shape (nConst, 1).
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

        Returns
        -------
        ndarray
            Array of shape (nB3, 1).
        """
        for body in self.Bodies:
            body._force = colvect([0.0, 0.0])
            body._torque = 0.0

        # Polymorphic force application (inject current time for UserForce/Actuator)
        for force in self.Forces:
            force._t = self.t
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

            if ff < 1.0e-8:
                flag = True
                break

            delta_c = np.linalg.lstsq(D, -Phi, rcond=None)[0]

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
        delta_v = np.linalg.lstsq(D, rhsv - D @ Phi, rcond=None)[0]

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

    def solve(self, analysis="dynamic", t_final=None, dt=None,
              ic_correct=False, t_eval=None, t_span=None):
        """Solve the model equations.

        The type of analysis is selected via the ``analysis`` parameter:

        * ``"dynamic"``   — full time-integration via CasADi DAE collocation.
        * ``"kinematic"`` — position/velocity/acceleration analysis at each
          time step without inertia; requires DOF = 0.
        * ``"static"``    — static equilibrium search; returns a single-step
          result (``T = [0.0]``).

        Parameters
        ----------
        analysis : str
            Type of analysis: ``"dynamic"`` (default), ``"kinematic"``
            or ``"static"``. Case-insensitive.
        t_final : float, optional
            Final simulation time. If *None* and interactive, prompts the user.
            Not used for static analysis.
        dt : float, optional
            Output time step. Used only when *t_eval* is not given.
        ic_correct : bool
            Whether to project initial conditions onto the constraint manifold
            before solving.
        t_eval : array-like, optional
            Explicit array of output time points.
        t_span : tuple of float, optional
            ``(t_start, t_end)`` shorthand; sets *t_final* when given.

        Returns
        -------
        tuple
            ``(T, uT)``: time vector of shape ``(nSteps,)`` and state matrix of
            shape ``(nSteps, 2*nB3)``.  Results are also distributed into each
            ``Body`` and ``Joint`` via ``_result_container``.
        """
        analysis = analysis.lower()
        _valid = ("dynamic", "kinematic", "static")
        if analysis not in _valid:
            raise ValueError(
                f"Unknown analysis type '{analysis}'. "
                f"Valid options: {_valid}"
            )

        if analysis == "dynamic":
            return self._solve_dynamic(
                t_final=t_final, dt=dt,
                ic_correct=ic_correct,
                t_eval=t_eval, t_span=t_span,
            )
        elif analysis == "kinematic":
            return self._solve_kinematic(
                t_final=t_final, dt=dt,
                ic_correct=ic_correct,
                t_eval=t_eval, t_span=t_span,
            )
        else:  # "static"
            return self._solve_static(ic_correct=ic_correct)

    # ------------------------------------------------------------------
    # Private: kinematic solver
    # ------------------------------------------------------------------

    def _solve_kinematic(self, t_final=None, dt=None, ic_correct=False,
                         t_eval=None, t_span=None):
        """Position/velocity/acceleration analysis at each time step.

        Requires DOF = 0 (fully determined model).  At least one driven
        joint (``RotMotion`` or ``TranMotion`` with a ``Function``)
        must exist to supply kinematic input at each instant.

        Algorithm per step *k*:

        1. **Positions** – Newton-Raphson: :math:`\\Phi(q) = 0`
        2. **Velocities** – linear solve: :math:`\\Phi_q \\dot{q} = \\nu(t_k)`
        3. **Accelerations** – linear solve: :math:`\\Phi_q \\ddot{q} = \\gamma(t_k, q, \\dot{q})`
        4. **Reactions** – inverse dynamics: :math:`\\Phi_q^T \\lambda = M\\ddot{q} - h_a`
        """
        # ---- DOF check ----
        nDOF = self.nDOF
        if nDOF != 0:
            raise ValueError(
                f"Kinematic analysis requires DOF = 0, but this model has "
                f"DOF = {nDOF}. Add driven joints (RotMotion / TranMotion) "
                f"or remove free bodies."
            )

        nB = len(self.Bodies)
        nB3 = 3 * nB
        nConst = self.nC

        # ---- Build time vector ----
        if t_span is not None and t_final is None:
            t_final = t_span[1]
        if t_final is None:
            t_final = float(input("\n\t ...Final time = ? "))
        if t_eval is not None:
            T = np.asarray(t_eval, dtype=float)
        elif dt is not None:
            T = np.arange(0.0, t_final + dt * 0.5, dt)
        else:
            dt_input = float(input("\t ...Reporting time-step = ? "))
            T = np.arange(0.0, t_final + dt_input * 0.5, dt_input)

        nSteps = len(T)
        uT = np.zeros((nSteps, 2 * nB3))

        # ---- Correct initial conditions if requested ----
        self.t = float(T[0])
        if ic_correct:
            self._ic_correct()

        # ---- Storage for post-processed quantities ----
        accelerations = np.zeros((nSteps, nB3))
        reactions = np.zeros((nSteps, nConst))

        print(f"\n\t ...Kinematic analysis: {nSteps} steps")

        pbar = tqdm(total=nSteps - 1,
                    desc="         ...Simulation progress",
                    bar_format="{l_bar}{bar}| [{n_fmt}/{total_fmt} steps, "
                               "Elapsed: {elapsed}, Remaining: {remaining}]",
                    colour="green")

        for k, t_k in enumerate(T):
            self.t = float(t_k)

            # --- 1. Position: Newton-Raphson Phi(q) = 0 ---
            _pos_converged = False
            for _iter in range(100):
                self._update_position()
                Phi = self._compute_constraints()
                res = float(np.linalg.norm(Phi))
                if res < 1.0e-8:
                    _pos_converged = True
                    break
                D = self._compute_jacobian()
                # Minimum-norm correction: delta_q = -D^+ Phi  (SVD/lstsq)
                delta_q = np.linalg.lstsq(D, -Phi.flatten(), rcond=None)[0]
                for Bi in range(nB):
                    ir = 3 * Bi
                    self.Bodies[Bi].position = (
                        self.Bodies[Bi].position + delta_q[ir:ir + 2].reshape(2, 1))
                    self.Bodies[Bi].orientation = (
                        self.Bodies[Bi].orientation + float(delta_q[ir + 2]))

            if not _pos_converged:
                raise RuntimeError(
                    f"Kinematic position Newton-Raphson did not converge "
                    f"at t = {t_k:.6g} (step {k})."
                )

            # --- 2. Velocity: Phi_q * dq = nu(t) ---
            D = self._compute_jacobian()
            nu = self._rhs_velocity()  # shape (nConst, 1)
            dq = np.linalg.lstsq(D, nu.flatten(), rcond=None)[0]  # 1-D
            for Bi in range(nB):
                ir = 3 * Bi
                self.Bodies[Bi].velocity = dq[ir:ir + 2].reshape(2, 1)
                self.Bodies[Bi].angular_velocity = float(dq[ir + 2])

            # --- 3. Acceleration: Phi_q * ddq = gamma(t, q, dq) ---
            self._update_velocity()
            gamma = self._rhs_acceleration()  # shape (nConst, 1)
            ddq = np.linalg.lstsq(D, gamma.flatten(), rcond=None)[0]  # 1-D
            for Bi in range(nB):
                ir = 3 * Bi
                self.Bodies[Bi].acceleration = ddq[ir:ir + 2].reshape(2, 1)
                self.Bodies[Bi].angular_acceleration = float(ddq[ir + 2])

            # --- 4. Reactions: Phi_q^T lambda = M*ddq - h_a ---
            h_a = self._compute_force()
            M_ddq = (self.M_array * ddq).reshape(-1, 1)
            rhs_lam = (M_ddq - h_a).flatten()
            # Solve D^T lambda = rhs_lam  (least-squares; D is square here)
            lam = np.linalg.lstsq(D.T, rhs_lam, rcond=None)[0]

            # --- Pack state vector row ---
            u_k = np.zeros(2 * nB3)
            for Bi in range(nB):
                ir = 3 * Bi
                u_k[ir]     = self.Bodies[Bi].position.flat[0]
                u_k[ir + 1] = self.Bodies[Bi].position.flat[1]
                u_k[ir + 2] = float(self.Bodies[Bi].orientation)
                u_k[nB3 + ir]     = self.Bodies[Bi].velocity.flat[0]
                u_k[nB3 + ir + 1] = self.Bodies[Bi].velocity.flat[1]
                u_k[nB3 + ir + 2] = float(self.Bodies[Bi].angular_velocity)

            uT[k] = u_k
            accelerations[k] = ddq.flatten()
            reactions[k] = lam.flatten()
            if k < nSteps - 1:
                pbar.update(1)

        pbar.close()
        print(f"\t ...Kinematic analysis completed successfully!")
        print(f"\n ")
        self._distribute_results(T, uT, accelerations, reactions)
        return T, uT

    # ------------------------------------------------------------------
    # Private: static solver
    # ------------------------------------------------------------------

    def _solve_static(self, ic_correct=False):
        """Find static equilibrium: :math:`\\ddot{q} = 0`, :math:`\\dot{q} = 0`.

        Solves the nonlinear system

        .. math::

            \\begin{bmatrix} \\Phi(q) \\\\ \\Phi_q^T \\lambda - h_a(q) \\end{bmatrix} = 0

        via Newton-Raphson with a finite-difference Jacobian (no analytic
        second derivatives required).

        Returns a single-step result: ``T = [0.0]``,
        ``uT`` of shape ``(1, 2*nB3)`` with zero velocities.
        """
        nB = len(self.Bodies)
        nB3 = 3 * nB
        nConst = self.nC

        # ---- Zero velocities (static assumption) ----
        for body in self.Bodies:
            body.velocity = np.zeros((2, 1))
            body.angular_velocity = 0.0

        # ---- Correct initial positions if requested ----
        self.t = 0.0
        if ic_correct:
            self._ic_correct()

        # ---- Build initial guess vector z = [q; lambda] ----
        def _pack_q():
            q = np.zeros(nB3)
            for Bi, body in enumerate(self.Bodies):
                ir = 3 * Bi
                q[ir]     = body.position.flat[0]
                q[ir + 1] = body.position.flat[1]
                q[ir + 2] = float(body.orientation)
            return q

        def _unpack_q(q):
            for Bi in range(nB):
                ir = 3 * Bi
                self.Bodies[Bi].position = q[ir:ir + 2].reshape(2, 1)
                self.Bodies[Bi].orientation = float(q[ir + 2])

        def _residual(z):
            """Residual F(q, lambda) = [Phi(q); Phi_q^T lam - h_a(q)]."""
            q = z[:nB3]
            lam = z[nB3:].reshape(-1, 1)
            _unpack_q(q)
            self._update_position()
            Phi = self._compute_constraints().flatten()  # (nConst,)
            D = self._compute_jacobian()                 # (nConst, nB3)
            h_a = self._compute_force().flatten()        # (nB3,)
            eq_forces = (D.T @ lam).flatten() - h_a     # (nB3,)
            return np.concatenate([Phi, eq_forces])

        q0 = _pack_q()
        if nConst > 0:
            # Initial Lagrange multiplier guess via least squares
            self._update_position()
            D0 = self._compute_jacobian()
            h_a0 = self._compute_force().flatten()
            lam0 = np.linalg.lstsq(D0.T, h_a0, rcond=None)[0]
        else:
            lam0 = np.zeros(0)

        z0 = np.concatenate([q0, lam0])

        # ---- Newton-Raphson with finite-difference Jacobian ----
        _tol = 1.0e-10
        _max_iter = 100
        _eps = 1.0e-7  # FD perturbation
        n_z = len(z0)
        z = z0.copy()
        converged = False

        print(f"\n\t ...Static equilibrium search")
        for _iter in range(_max_iter):
            F = _residual(z)
            norm_F = float(np.linalg.norm(F))
            if norm_F < _tol:
                converged = True
                break

            # Finite-difference Jacobian
            J = np.zeros((n_z, n_z))
            for j in range(n_z):
                dz = np.zeros(n_z)
                dz[j] = _eps
                J[:, j] = (_residual(z + dz) - F) / _eps

            try:
                delta = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                delta = np.linalg.lstsq(J, -F, rcond=None)[0]

            z = z + delta

        if not converged:
            raise RuntimeError(
                f"Static equilibrium Newton-Raphson did not converge after "
                f"{_max_iter} iterations (final residual = {norm_F:.3e})."
            )

        # ---- Unpack final solution ----
        q_eq = z[:nB3]
        lam_eq = z[nB3:]
        _unpack_q(q_eq)
        for body in self.Bodies:
            body.velocity = np.zeros((2, 1))
            body.angular_velocity = 0.0
            body.acceleration = np.zeros((2, 1))
            body.angular_acceleration = 0.0

        # ---- Build output arrays ----
        T = np.array([0.0])
        u_eq = np.zeros(2 * nB3)
        for Bi in range(nB):
            ir = 3 * Bi
            u_eq[ir]     = self.Bodies[Bi].position.flat[0]
            u_eq[ir + 1] = self.Bodies[Bi].position.flat[1]
            u_eq[ir + 2] = float(self.Bodies[Bi].orientation)
            # velocities remain 0

        uT = u_eq.reshape(1, -1)

        # accelerations = 0 everywhere; reactions from lam_eq
        accelerations = np.zeros((1, nB3))
        reactions = lam_eq.reshape(1, -1) if nConst > 0 else np.zeros((1, 0))

        print(f"\t ...Static equilibrium found (iter = {_iter + 1})")
        print(f"\n ")
        self._distribute_results(T, uT, accelerations, reactions)
        return T, uT

    # ------------------------------------------------------------------
    # Private: distribute results for kinematic / static analyses
    # ------------------------------------------------------------------

    def _distribute_results(self, T, uT, accelerations, reactions):
        """Populate ``_result_container`` from pre-computed arrays.

        Used by ``_solve_kinematic`` and ``_solve_static`` which compute
        accelerations and reactions analytically (not via ``_post_process``).

        Parameters
        ----------
        T : ndarray
            Time vector, shape (nSteps,).
        uT : ndarray
            State matrix, shape (nSteps, 2*nB3).
        accelerations : ndarray
            Acceleration matrix, shape (nSteps, nB3).
        reactions : ndarray
            Reaction matrix, shape (nSteps, nConst).
        """
        nB3 = 3 * self.nB
        for Bi, body in enumerate(self.Bodies):
            i = 3 * Bi
            body._result_container = {
                'positions': {
                    'x':   uT[:, i],
                    'y':   uT[:, i + 1],
                    'phi': uT[:, i + 2],
                },
                'velocities': {
                    'dx':   uT[:, nB3 + i],
                    'dy':   uT[:, nB3 + i + 1],
                    'dphi': uT[:, nB3 + i + 2],
                },
                'accelerations': {
                    'ddx':   accelerations[:, i],
                    'ddy':   accelerations[:, i + 1],
                    'ddphi': accelerations[:, i + 2],
                },
            }
        for joint in self.Joints:
            rs = joint._rows
            re = joint._rowe
            joint._result_container = {'reactions': reactions[:, rs:re]}

    # ------------------------------------------------------------------
    # Private: CasADi DAE solver (Radau collocation)
    # ------------------------------------------------------------------

    def _build_casadi_phi(self, ca, q_sym, t_sym):
        """Build symbolic constraint vector Phi(q, t) using CasADi SX.

        Walks ``self.Joints`` and translates each constraint equation into
        a pure CasADi SX expression.  The Jacobian is **not** built here —
        CasADi computes it automatically via ``ca.jacobian``.

        Supported joint types: RevJoint, TranJoint, RevRevJoint,
        RevTranJoint, RigidJoint, DiscJoint, RotMotion, TranMotion.

        Parameters
        ----------
        ca : module
            The CasADi module.
        q_sym : ca.SX
            Symbolic position vector of shape (nB3, 1).
        t_sym : ca.SX
            Symbolic time variable.

        Returns
        -------
        ca.SX
            Column vector of shape (nC, 1).
        """
        from .constraints import (RevJoint, TranJoint, RevRevJoint,
                                  RevTranJoint, RigidJoint, DiscJoint)
        from .motion import RotMotion, TranMotion
        from .model import Ground

        phi_blocks = []

        for joint in self.Joints:
            iBody = joint.iBody
            jBody = joint.jBody

            # --- helper: extract symbolic body state from q_sym ---
            def _body_state(body):
                if body is Ground:
                    return (ca.SX.zeros(2, 1), ca.SX(0),
                            ca.SX.eye(2))  # pos, phi, R
                Bi = body._body_index - 1
                ir = 3 * Bi
                pos = q_sym[ir:ir + 2]
                phi = q_sym[ir + 2]
                R = ca.vertcat(
                    ca.horzcat(ca.cos(phi), -ca.sin(phi)),
                    ca.horzcat(ca.sin(phi),  ca.cos(phi)),
                )
                return pos, phi, R

            def _marker_rP(marker, body):
                """Global position of a marker on *body*."""
                pos, phi, R = _body_state(body)
                sP_local = marker.local_position.reshape(2, 1)
                sP = R @ sP_local
                return pos + sP

            def _marker_u(marker, body):
                """Global unit direction of a marker on *body*."""
                _, phi, R = _body_state(body)
                return R @ marker._ulocal

            # ---- Joint-specific symbolic Phi ----

            if isinstance(joint, RevJoint):
                rPi = _marker_rP(joint.iMarker, iBody)
                rPj = _marker_rP(joint.jMarker, jBody)
                block = rPi - rPj  # (2, 1)
                if joint.fix == 1:
                    _, phi_i, _ = _body_state(iBody)
                    _, phi_j, _ = _body_state(jBody)
                    if iBody is Ground:
                        extra = -phi_j - joint._p0
                    elif jBody is Ground:
                        extra = phi_i - joint._p0
                    else:
                        extra = phi_i - phi_j - joint._p0
                    block = ca.vertcat(block, extra)

            elif isinstance(joint, TranJoint):
                rPi = _marker_rP(joint.iMarker, iBody)
                rPj = _marker_rP(joint.jMarker, jBody)
                uj = _marker_u(joint.jMarker, jBody)
                ujr = ca.vertcat(-uj[1], uj[0])
                ui = _marker_u(joint.iMarker, iBody)
                d = rPi - rPj
                block = ca.vertcat(ujr.T @ d, ujr.T @ ui)
                if joint.fix == 1:
                    block = ca.vertcat(block, (ui.T @ d - joint._p0) / 2)

            elif isinstance(joint, RevRevJoint):
                rPi = _marker_rP(joint.iMarker, iBody)
                rPj = _marker_rP(joint.jMarker, jBody)
                d = rPi - rPj
                L = joint.L
                block = (d.T @ d - L**2) / (2 * L)

            elif isinstance(joint, RevTranJoint):
                rPi = _marker_rP(joint.iMarker, iBody)
                rPj = _marker_rP(joint.jMarker, jBody)
                ui = _marker_u(joint.iMarker, iBody)
                uir = ca.vertcat(-ui[1], ui[0])
                d = rPi - rPj
                block = uir.T @ d - joint.L

            elif isinstance(joint, RigidJoint):
                _, phi_i, _ = _body_state(iBody)
                _, phi_j, Rj = _body_state(jBody)
                pos_i, _, _ = _body_state(iBody)
                pos_j, _, _ = _body_state(jBody)
                d0 = joint.d0
                if not isinstance(d0, np.ndarray):
                    d0 = np.array(d0).reshape(-1, 1)
                if iBody is Ground:
                    block = ca.vertcat(
                        -(pos_j + Rj @ d0),
                        -phi_j - joint._p0,
                    )
                elif jBody is Ground:
                    block = ca.vertcat(
                        pos_i - d0,
                        phi_i - joint._p0,
                    )
                else:
                    block = ca.vertcat(
                        pos_i - (pos_j + Rj @ d0),
                        phi_i - phi_j - joint._p0,
                    )

            elif isinstance(joint, DiscJoint):
                pos_i, phi_i, _ = _body_state(iBody)
                block = ca.vertcat(
                    pos_i[1] - joint.R,
                    (pos_i[0] - joint.x0) + joint.R * (phi_i - joint._p0),
                )

            elif isinstance(joint, RotMotion):
                # Build symbolic function value from the Function object
                fun_sx = self._casadi_functEval(ca, joint.iFunct, t_sym)
                _, phi_i, _ = _body_state(iBody)
                _, phi_j, _ = _body_state(jBody)
                if iBody is Ground:
                    block = -phi_j - fun_sx
                elif jBody is Ground:
                    block = phi_i - fun_sx
                else:
                    block = phi_i - phi_j - fun_sx

            elif isinstance(joint, TranMotion):
                rPi = _marker_rP(joint.iMarker, iBody)
                rPj = _marker_rP(joint.jMarker, jBody)
                d = rPi - rPj
                fun_sx = self._casadi_functEval(ca, joint.iFunct, t_sym)
                block = (d.T @ d - fun_sx**2) / 2

            else:
                raise NotImplementedError(
                    f"CasADi symbolic Phi not implemented for "
                    f"{type(joint).__name__}")

            phi_blocks.append(block)

        if not phi_blocks:
            return ca.SX.zeros(0, 1)
        return ca.vertcat(*phi_blocks)

    @staticmethod
    def _casadi_functEval(ca, funct, t_sym):
        """Build CasADi SX expression for a Function value f(t).

        Supports types 'a', 'b', 'c' and the callable *expr* mode
        (``funct.expr(t_sym)``).  Mirrors ``mechanics.functEval`` but
        as a single symbolic SX graph using ``ca.if_else``.
        """
        # Symbolic/callable mode: let the callable build the SX graph
        if funct.expr is not None:
            return funct.expr(t_sym)

        ftype = funct.type
        c = funct.coeff

        if ftype == 'a':
            return c[0] + c[1] * t_sym + c[2] * t_sym**2

        elif ftype in ('b', 'c'):
            t0 = float(funct.t_start)
            te = float(funct.t_end)
            fs = float(funct.f_start)
            tau = t_sym - t0

            if ftype == 'b':
                f_mid  = fs + c[0]*tau**3 + c[1]*tau**4 + c[2]*tau**5
                tau_e  = te - t0
                f_end  = float(fs + c[0]*tau_e**3 + c[1]*tau_e**4 + c[2]*tau_e**5)
            else:  # 'c'
                f_mid    = fs + c[0]*tau**4 + c[1]*tau**5 + c[2]*tau**6
                tau_e    = te - t0
                f_end_v  = float(fs + c[0]*tau_e**4 + c[1]*tau_e**5 + c[2]*tau_e**6)
                fd_end_v = float(c[3]*tau_e**3 + c[4]*tau_e**4 + c[5]*tau_e**5)
                f_end    = f_end_v + fd_end_v * (t_sym - te)

            return ca.if_else(t_sym < t0, fs,
                              ca.if_else(t_sym >= te, f_end, f_mid))

        else:
            raise ValueError(f"Unknown function type '{ftype}'.")

    def _build_casadi_forces(self, ca, q_sym, v_sym, t_sym):
        """Build symbolic generalised force vector Q(q, v, t) using CasADi SX.

        Supported force types: Weight, PtpForce, RotSdaForce, LocalForce,
        GlobalForce, Torque.

        Parameters
        ----------
        ca : module
            The CasADi module.
        q_sym : ca.SX
            Symbolic position vector of shape (nB3, 1).
        v_sym : ca.SX
            Symbolic velocity vector of shape (nB3, 1).
        t_sym : ca.SX
            Symbolic time variable.

        Returns
        -------
        ca.SX
            Column vector of shape (nB3, 1).
        """
        from .forces import (Weight, PtpForce, RotSdaForce,
                              LocalForce, GlobalForce, Torque,
                              UserForce, Actuator)
        from .model import Ground

        nB3 = 3 * len(self.Bodies)
        Q = ca.SX.zeros(nB3, 1)
        has_user_forces = False

        def _body_sym(body, section='q'):
            if body is Ground:
                return None, None, None, None
            Bi = body._body_index - 1
            ir = 3 * Bi
            if section == 'q':
                pos = q_sym[ir:ir + 2]
                phi = q_sym[ir + 2]
            else:
                pos = v_sym[ir:ir + 2]
                phi = v_sym[ir + 2]
            R = ca.vertcat(
                ca.horzcat(ca.cos(q_sym[3*Bi + 2]), -ca.sin(q_sym[3*Bi + 2])),
                ca.horzcat(ca.sin(q_sym[3*Bi + 2]),  ca.cos(q_sym[3*Bi + 2])),
            )
            return pos, phi, R, ir

        def _marker_sPr_sym(marker, body):
            """Symbolic rotate_90(R @ sP_local) for torque computation."""
            _, phi_b, R, _ = _body_sym(body, 'q')
            if R is None:
                sP = ca.SX(marker.local_position.reshape(2, 1))
                return ca.vertcat(-sP[1], sP[0])
            sP = R @ marker.local_position.reshape(2, 1)
            return ca.vertcat(-sP[1], sP[0])

        def _marker_rP_sym(marker, body):
            if body is Ground:
                return ca.SX(marker.local_position.reshape(2, 1))
            _, _, R, ir = _body_sym(body, 'q')
            pos = q_sym[ir:ir + 2]
            sP = R @ marker.local_position.reshape(2, 1)
            return pos + sP

        def _marker_drP_sym(marker, body):
            """Symbolic velocity of marker global position."""
            if body is Ground:
                return ca.SX.zeros(2, 1)
            Bi = body._body_index - 1
            ir = 3 * Bi
            vel = v_sym[ir:ir + 2]
            omega = v_sym[ir + 2]
            _, _, R, _ = _body_sym(body, 'q')
            sP = R @ marker.local_position.reshape(2, 1)
            sPr = ca.vertcat(-sP[1], sP[0])
            dsP = sPr * omega
            return vel + dsP

        for force in self.Forces:
            if isinstance(force, Weight):
                g_val = force.gravity
                g_dir = force.gravity_direction.flatten()
                for body in self.Bodies:
                    Bi = body._body_index - 1
                    ir = 3 * Bi
                    w = body.mass * g_val * ca.SX(g_dir.reshape(2, 1))
                    Q[ir:ir + 2] += w

            elif isinstance(force, PtpForce):
                rPi = _marker_rP_sym(force.iMarker, force.iBody)
                rPj = _marker_rP_sym(force.jMarker, force.jBody)
                d = rPi - rPj
                L = ca.norm_2(d)
                u = d / L

                drPi = _marker_drP_sym(force.iMarker, force.iBody)
                drPj = _marker_drP_sym(force.jMarker, force.jBody)
                dd = drPi - drPj
                dL = (d.T @ dd) / L
                delta = L - force.L0
                f_mag = force.k * delta + force.dc * dL + force.f_a
                fi = f_mag * u

                if force.iBody is not Ground:
                    Bi = force.iBody._body_index - 1
                    ir = 3 * Bi
                    Q[ir:ir + 2] -= fi
                    sPr_i = _marker_sPr_sym(force.iMarker, force.iBody)
                    Q[ir + 2] -= (sPr_i.T @ fi)

                if force.jBody is not Ground:
                    Bj = force.jBody._body_index - 1
                    jr = 3 * Bj
                    Q[jr:jr + 2] += fi
                    sPr_j = _marker_sPr_sym(force.jMarker, force.jBody)
                    Q[jr + 2] += (sPr_j.T @ fi)

            elif isinstance(force, RotSdaForce):
                iB = force.iBody
                jB = force.jBody
                if iB is Ground:
                    Bj = jB._body_index - 1
                    jr = 3 * Bj
                    theta = -q_sym[jr + 2]
                    theta_d = -v_sym[jr + 2]
                    T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                    Q[jr + 2] += T
                elif jB is Ground:
                    Bi = iB._body_index - 1
                    ir = 3 * Bi
                    theta = q_sym[ir + 2]
                    theta_d = v_sym[ir + 2]
                    T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                    Q[ir + 2] -= T
                else:
                    Bi = iB._body_index - 1
                    Bj = jB._body_index - 1
                    ir = 3 * Bi
                    jr = 3 * Bj
                    theta = q_sym[ir + 2] - q_sym[jr + 2]
                    theta_d = v_sym[ir + 2] - v_sym[jr + 2]
                    T = force.k * (theta - force.theta0) + force.dc * theta_d + force.T_a
                    Q[ir + 2] -= T
                    Q[jr + 2] += T

            elif isinstance(force, LocalForce):
                Bi = force.iBody._body_index - 1
                ir = 3 * Bi
                _, _, R, _ = _body_sym(force.iBody, 'q')
                f_local = force.force_local
                if not isinstance(f_local, np.ndarray):
                    f_local = np.array(f_local).reshape(2, 1)
                Q[ir:ir + 2] += R @ f_local

            elif isinstance(force, GlobalForce):
                Bi = force.iBody._body_index - 1
                ir = 3 * Bi
                f_glob = force.force_global
                if not isinstance(f_glob, np.ndarray):
                    f_glob = np.array(f_glob).reshape(2, 1)
                Q[ir:ir + 2] += ca.SX(f_glob)

            elif isinstance(force, Torque):
                Bi = force.iBody._body_index - 1
                ir = 3 * Bi
                Q[ir + 2] += force.torque_value

            elif isinstance(force, UserForce):
                has_user_forces = True  # handled via numeric parameters

            elif isinstance(force, Actuator):
                if force.control == 'force':
                    has_user_forces = True  # evaluated numerically each step
                # 'length' mode is handled by the TranMotion joint, skip

            else:
                raise NotImplementedError(
                    f"CasADi symbolic Q not implemented for "
                    f"{type(force).__name__}")

        return Q, has_user_forces

    def _solve_dynamic(self, t_final=None, dt=None, ic_correct=False,
                       t_eval=None, t_span=None):
        """Time-integrate the equations of motion via CasADi DAE collocation.

        Builds the constraint and force equations as CasADi SX
        expressions, then integrates step-by-step with an implicit
        collocation scheme (Radau IIA) that handles the index-2 DAE
        arising from position-level constraints.

        The formulation uses differential states ``x = [q, v]`` and
        algebraic states ``z = [lam]`` in the semi-explicit form::

            q_dot = v
            M v_dot = Q(t, q, v) - Phi_q^T(q, t) lam
            0 = Phi(q, t)

        Parameters
        ----------
        t_final : float, optional
            Final simulation time.
        dt : float, optional
            Output time step (used when *t_eval* is absent).
        ic_correct : bool, optional
            Project initial conditions before integrating.
        t_eval : array-like, optional
            Explicit array of output time points.
        t_span : tuple of float, optional
            ``(t0, tf)`` shorthand; overrides *t_final*.

        Returns
        -------
        tuple
            ``(T, uT)``.
        """
        nB  = len(self.Bodies)
        nB3 = 3 * nB
        nC  = self.nC

        # ---- Build time grid ----
        if t_span is not None and t_final is None:
            t_final = t_span[1]
        if t_eval is not None:
            T = np.asarray(t_eval, dtype=float)
        elif t_final is not None and dt is not None:
            T = np.arange(0.0, t_final + dt * 0.5, dt)
        elif t_final is None:
            t_final = float(input("\n\t ...Final time = ? "))
            if dt is None:
                dt = float(input("\t ...Reporting time-step = ? "))
            T = np.arange(0.0, t_final + dt * 0.5, dt)
        else:
            dt_input = float(input("\t ...Reporting time-step = ? "))
            T = np.arange(0.0, t_final + dt_input * 0.5, dt_input)
        nSteps = len(T)

        # ---- Correct initial conditions if requested ----
        self.t = float(T[0])
        if ic_correct:
            self._ic_correct()

        # ---- Pack initial state ----
        q0 = np.zeros(nB3)
        v0 = np.zeros(nB3)
        for Bi, body in enumerate(self.Bodies):
            ir = 3 * Bi
            q0[ir]     = body.position.flat[0]
            q0[ir + 1] = body.position.flat[1]
            q0[ir + 2] = float(body.orientation)
            v0[ir]     = body.velocity.flat[0]
            v0[ir + 1] = body.velocity.flat[1]
            v0[ir + 2] = float(body.angular_velocity)

        self._update_position()
        self._update_velocity()
        Q0_num = self._compute_force().flatten()
        if nC > 0:
            D0 = self._compute_jacobian()
            # Velocity-level IC correction: project v0 so Phi_q*v + Phi_t = 0
            Phi0 = self._compute_constraints().flatten()
            eps_t = 1e-8
            self.t = float(T[0]) + eps_t
            self._update_position()  # re-eval for t+eps
            Phi0p = self._compute_constraints().flatten()
            self.t = float(T[0]) - eps_t
            self._update_position()
            Phi0m = self._compute_constraints().flatten()
            self.t = float(T[0])
            self._update_position()
            Phi_t = (Phi0p - Phi0m) / (2 * eps_t)
            rhs = -(D0 @ v0 + Phi_t)
            if np.linalg.norm(rhs) > 1e-12:
                v0 += np.linalg.lstsq(D0, rhs, rcond=None)[0]
                # Update body velocities with corrected v0
                for Bi, body in enumerate(self.Bodies):
                    ir = 3 * Bi
                    body.velocity         = v0[ir:ir + 2].reshape(2, 1)
                    body.angular_velocity = float(v0[ir + 2])
                self._update_velocity()
                Q0_num = self._compute_force().flatten()
            lam0 = np.linalg.lstsq(D0.T, Q0_num, rcond=None)[0]
        else:
            lam0 = np.zeros(0)

        # ---- Detect uniform time grid ----
        # A uniform grid lets us compile the CasADi integrator exactly ONCE.
        # We introduce "tau" as an explicit state variable that tracks
        # absolute time so the step map [q, v, tau] -> [q', v', tau'] is
        # identical for every step (only the initial tau changes).
        if nSteps > 1:
            _diffs = np.diff(T)
            _uniform = bool(np.allclose(_diffs, _diffs[0], rtol=1e-9, atol=0.0))
            _h = float(_diffs[0]) if _uniform else 0.0
        else:
            _uniform = False
            _h = 0.0

        # ---- Build CasADi symbolic DAE ----
        q_sym   = ca.SX.sym('q',   nB3)
        v_sym   = ca.SX.sym('v',   nB3)
        qd_sym  = ca.SX.sym('qd',  nB3)
        vd_sym  = ca.SX.sym('vd',  nB3)
        t_sym   = ca.SX.sym('t')

        if _uniform:
            # tau is the absolute time state; forces/constraints use it
            # instead of the CasADi integration-time symbol.
            tau_sym   = ca.SX.sym('tau')
            _time_arg = tau_sym
        else:
            tau_sym   = None
            _time_arg = t_sym

        Phi_sx = self._build_casadi_phi(ca, q_sym, _time_arg)       # (nC, 1)
        Q_sx, has_user_forces = self._build_casadi_forces(
            ca, q_sym, v_sym, _time_arg)                              # (nB3, 1)

        M_vec  = ca.SX(self.M_array.reshape(-1, 1))
        invM   = ca.SX(self.invM_array.reshape(-1, 1))

        # If any UserForce exists, add a numeric parameter vector
        if has_user_forces:
            p_user = ca.SX.sym('p_user', nB3)
            Q_sx = Q_sx + p_user
        else:
            p_user = ca.SX()

        if nC > 0:
            lam_sym = ca.SX.sym('lam', nC)
            D_sx = ca.jacobian(Phi_sx, q_sym)  # (nC, nB3) — automatic!

            if _uniform:
                # State: [q, v, tau].  ode: [v, M^{-1}(Q - D^T lam), 1]
                # No 't' key: absolute time tracked by tau, not by CasADi.
                ode = ca.vertcat(v_sym, invM * (Q_sx - D_sx.T @ lam_sym),
                                 ca.SX.ones(1))
                dae = {
                    'x':   ca.vertcat(q_sym, v_sym, tau_sym),
                    'z':   lam_sym,
                    'p':   p_user,
                    'ode': ode,
                    'alg': Phi_sx,
                }
            else:
                # Semi-explicit DAE:
                #   x = [q, v]   (differential)
                #   z = [lam]    (algebraic)
                #   ode: dx/dt = [v, M^{-1}(Q - D^T lam)]
                #   alg: 0 = Phi(q, t)
                ode = ca.vertcat(v_sym, invM * (Q_sx - D_sx.T @ lam_sym))
                alg = Phi_sx
                dae = {
                    'x':   ca.vertcat(q_sym, v_sym),
                    'z':   lam_sym,
                    'p':   p_user,
                    't':   t_sym,
                    'ode': ode,
                    'alg': alg,
                }
        else:
            if _uniform:
                # State: [q, v, tau], pure ODE
                ode = ca.vertcat(v_sym, invM * Q_sx, ca.SX.ones(1))
                dae = {
                    'x':   ca.vertcat(q_sym, v_sym, tau_sym),
                    'z':   ca.SX(),
                    'p':   p_user,
                    'ode': ode,
                    'alg': ca.SX(),
                }
            else:
                # No constraints — pure ODE
                ode = ca.vertcat(v_sym, invM * Q_sx)
                dae = {
                    'x':   ca.vertcat(q_sym, v_sym),
                    'z':   ca.SX(),
                    'p':   p_user,
                    't':   t_sym,
                    'ode': ode,
                    'alg': ca.SX(),
                }

        # ---- Collocation integrator options ----
        opts = {
            'rootfinder':                'kinsol',
            'collocation_scheme':        'radau',
            'number_of_finite_elements': 1,
        }

        # ---- Storage ----
        uT            = np.zeros((nSteps, 2 * nB3))
        accelerations = np.zeros((nSteps, nB3))
        reactions     = np.zeros((nSteps, nC))

        # Row 0 — initial state
        uT[0] = np.concatenate([q0, v0])
        reactions[0] = lam0
        if nC > 0:
            accelerations[0] = (Q0_num - D0.T @ lam0) / self.M_array
        else:
            accelerations[0] = Q0_num / self.M_array

        # ---- Build integrator once for uniform grids ----
        if _uniform:
            _integrator = ca.integrator('col_step', 'collocation', dae, 0.0, _h, opts)

        print(f"\n\t ...DAE analysis (CasADi-collocation): {nSteps} steps"
              f"{'  [one-shot integrator]' if _uniform else ''}")

        pbar = tqdm(total=nSteps - 1,
                    desc="         ...Simulation progress",
                    bar_format="{l_bar}{bar}| [{n_fmt}/{total_fmt} steps, "
                               "Elapsed: {elapsed}, Remaining: {remaining}]",
                    colour="green")

        # ---- Helper: evaluate UserForce/Actuator contributions numerically ----
        if has_user_forces:
            from .forces import UserForce, Actuator

            def _eval_user_forces():
                """Zero body forces, apply only UserForce/Actuator, collect vector."""
                for body in self.Bodies:
                    body._force = colvect([0.0, 0.0])
                    body._torque = 0.0
                for force in self.Forces:
                    if isinstance(force, (UserForce, Actuator)):
                        force._t = self.t
                        force.apply(self.Bodies)
                g = np.zeros(nB3)
                for Bi, body in enumerate(self.Bodies):
                    ir = 3 * Bi
                    g[ir]     = float(body._force.flat[0])
                    g[ir + 1] = float(body._force.flat[1])
                    g[ir + 2] = float(body._torque)
                return g

        # ---- Step-by-step integration ----
        x_curr = np.concatenate([q0, v0])
        z_curr = lam0.copy()

        for k in range(1, nSteps):
            t0_k = float(T[k - 1])
            tf_k = float(T[k])
            pbar.update(1)

            # Evaluate UserForce contributions at current state
            if has_user_forces:
                q_k = x_curr[:nB3]
                v_k = x_curr[nB3:]
                for Bi, body in enumerate(self.Bodies):
                    ir = 3 * Bi
                    body.position         = q_k[ir:ir + 2].reshape(2, 1)
                    body.orientation      = float(q_k[ir + 2])
                    body.velocity         = v_k[ir:ir + 2].reshape(2, 1)
                    body.angular_velocity = float(v_k[ir + 2])
                self.t = t0_k
                self._update_position()
                self._update_velocity()
                p_num = _eval_user_forces()
            else:
                p_num = np.zeros(0)

            if _uniform:
                # Augment state with current absolute time (tau initial cond.),
                # run the precompiled integrator, then strip tau from the result.
                x0_aug = np.concatenate([x_curr, [t0_k]])
                res    = _integrator(x0=x0_aug, z0=z_curr, p=p_num)
                xf     = np.array(res['xf']).flatten()
                x_curr = xf[:-1]  # drop tau state
            else:
                res    = ca.integrator(
                    'col_step', 'collocation', dae, t0_k, tf_k, opts
                )(x0=x_curr, z0=z_curr, p=p_num)
                x_curr = np.array(res['xf']).flatten()

            z_curr = np.array(res['zf']).flatten() if nC > 0 else np.zeros(0)

            uT[k] = x_curr
            reactions[k] = z_curr

            # Compute accelerations: vdot = M^{-1}(Q - D^T lam)
            q_k = x_curr[:nB3]
            v_k = x_curr[nB3:]
            for Bi, body in enumerate(self.Bodies):
                ir = 3 * Bi
                body.position         = q_k[ir:ir + 2].reshape(2, 1)
                body.orientation      = float(q_k[ir + 2])
                body.velocity         = v_k[ir:ir + 2].reshape(2, 1)
                body.angular_velocity = float(v_k[ir + 2])
            self.t = tf_k
            self._update_position()
            self._update_velocity()
            Q_k = self._compute_force().flatten()
            if nC > 0:
                D_k = self._compute_jacobian()
                accelerations[k] = (Q_k - D_k.T @ z_curr) / self.M_array
            else:
                accelerations[k] = Q_k / self.M_array

        # ---- Compute max constraint violation ----
        pbar.close()
        max_phi = 0.0
        if nC > 0:
            for k in range(nSteps):
                q_k = uT[k, :nB3]
                for Bi, body in enumerate(self.Bodies):
                    ir = 3 * Bi
                    body.position    = q_k[ir:ir + 2].reshape(2, 1)
                    body.orientation = float(q_k[ir + 2])
                self.t = float(T[k])
                self._update_position()
                Phi_k = self._compute_constraints().flatten()
                max_phi = max(max_phi, float(np.max(np.abs(Phi_k))))

        print(f"\t ...DAE analysis (collocation) completed successfully!")
        print(f"\t ...Max constraint violation: {max_phi:.3e}")
        print(f"\n ")

        self._distribute_results(T, uT, accelerations, reactions)
        return T, uT

