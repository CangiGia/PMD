"""
- Planar Multi-body Dynamics simulation solver - 
== 

This Python module provides the necessary algorithms and functions to solve 
planar multi-body dynamic models.

Author: - Giacomo Cangi, PhD student @ UniPG -
"""

import os
import numpy as np
import scipy as sc
import numpy.linalg as lng
from pmd_functions import *
from scipy.integrate import solve_ivp

# for debugging porpouse
import pdb 


class PlanarDynamicModel:
    def __init__(self):
        grouped_calsses = group_classes()
        self.Bodies = grouped_calsses.get("Body", [])
        self.Points = grouped_calsses.get("Point", [])
        self.uVectors = grouped_calsses.get("uVector", [])
        self.Forces = grouped_calsses.get("Force", [])
        self.Joints = grouped_calsses.get("Joint", [])
        self.Functs = grouped_calsses.get("Function", []) if "Function" in grouped_calsses else []

        # initialize the model for simulation automatically
        self.__initialize()

    def __initialize(self):
        """
        Initializi the multi-body model considering the values typed 
        by the user.
        """
        # initialize variables
        nB = len(self.Bodies)
        nB3 = 3 * nB
        nB6 = 6 * nB

        # bodies
        ##### ##### ##### ##### #####
        for Bi in range(nB):
            body = self.Bodies[Bi]
            body.irc = 3 * Bi + 1
            body.irv = nB3 + 3 * Bi + 1
            body.m_inv = 1 / body.m
            body.J_inv = 1 / body.J
            body.A = A_matrix(body.p)

        # mass (inertia) matrix as an array
        self.M_array = np.zeros((nB3, 1))
        self.M_inv_array = np.zeros((nB3, 1))
        for Bi in range(nB):
            is_ = 3 * Bi
            ie_ = is_ + 3
            self.M_array[is_:ie_] = np.array([[self.Bodies[Bi].m], [self.Bodies[Bi].m], [self.Bodies[Bi].J]])
            self.M_inv_array[is_:ie_] = np.array([[self.Bodies[Bi].m_inv], [self.Bodies[Bi].m_inv], [self.Bodies[Bi].J_inv]])
        ##### ##### ##### ##### #####  

        # points #! CHECK CON MATLAB - Probabile problema di indicizzazione
        ##### ##### ##### ##### #####
        nPtot = len(self.Points)
        for Pi in range(nPtot):
            point = self.Points[Pi]
            if point.Bindex == 0:
                point.sP = point.sPlocal
                point.sP_r = s_rot(point.sP)
                point.rP = point.sP

            for Bi in range(nB):
                if point.Bindex == (Bi+1): #! aggiunto il +1 per correggere l'assegnazione dei punti, check con MATLAB
                    self.Bodies[Bi].pts.append(Pi)  # append point index to the body's points
        ##### ##### ##### ##### #####

        # unit vectors
        ##### ##### ##### ##### #####
        nU = len(self.uVectors)
        for Vi in range(nU):
            unit_vector = self.uVectors[Vi]
            if unit_vector.Bindex == 0:
                unit_vector.u = unit_vector.ulocal
                unit_vector.u_r = s_rot(unit_vector.u)
        ##### ##### ##### ##### #####

        # force elements
        ##### ##### ##### ##### #####
        nF = len(self.Forces)
        for Fi in range(nF):
            force = self.Forces[Fi]
            if force.type == 'weight':
                ug = force.gravity * force.wgt
                for Bi in range(nB):
                    self.Bodies[Bi].wgt = self.Bodies[Bi].m * ug
            elif force.type == 'ptp':
                Pi = force.iPindex
                Pj = force.jPindex
                force.iBindex = self.Points[Pi].Bindex
                force.jBindex = self.Points[Pj].Bindex
        ##### ##### ##### ##### #####

        # joints
        ##### ##### ##### ##### #####
        nJ = len(self.Joints)
        cfriction = 0

        for Ji in range(nJ):
            joint = self.Joints[Ji]
            if joint.type == 'rev':
                joint.mrows = 2
                joint.nbody = 2
                Pi = joint.iPindex
                Pj = joint.jPindex
                joint.iBindex = self.Points[Pi].Bindex
                joint.jBindex = self.Points[Pj].Bindex
                if joint.fix == 1:
                    joint.mrows = 3
                    if joint.iBindex == 0:
                        joint.p0 = -self.Bodies[joint.jBindex].p
                    elif joint.jBindex == 0:
                        joint.p0 = self.Bodies[joint.iBindex].p
                    else:
                        joint.p0 = self.Bodies[joint.iBindex].p - self.Bodies[joint.jBindex].p

            elif joint.type == 'tran':
                joint.mrows = 2
                joint.nbody = 2
                Pi = joint.iPindex
                Pj = joint.jPindex
                joint.iBindex = self.Points[Pi].Bindex
                joint.jBindex = self.Points[Pj].Bindex
                if joint.fix == 1:
                    joint.mrows = 3
                    if joint.iBindex == 0:
                        joint.p0 = np.linalg.norm(self.Points[Pi].rP - 
                                                self.Bodies[joint.jBindex].r - 
                                                self.Bodies[joint.jBindex].A @ 
                                                self.Points[Pj].sPlocal)
                    elif joint.jBindex == 0:
                        joint.p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + 
                                                self.Bodies[joint.iBindex].A @ 
                                                self.Points[Pi].sPlocal - 
                                                self.Points[Pj].rP)
                    else:
                        joint.p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + 
                                                self.Bodies[joint.iBindex].A @ 
                                                self.Points[Pi].sPlocal - 
                                                self.Bodies[joint.jBindex].r - 
                                                self.Bodies[joint.jBindex].A @ 
                                                self.Points[Pj].sPlocal)

            elif joint.type == 'rev-rev':
                joint.mrows = 1
                joint.nbody = 2
                Pi = joint.iPindex
                Pj = joint.jPindex
                joint.iBindex = self.Points[Pi].Bindex
                joint.jBindex = self.Points[Pj].Bindex

            elif joint.type == 'rev-tran':
                joint.mrows = 1
                joint.nbody = 2
                Pi = joint.iPindex
                Pj = joint.jPindex
                joint.iBindex = self.Points[Pi].Bindex
                joint.jBindex = self.Points[Pj].Bindex

            elif joint.type in {'rel-rot', 'rel-tran'}:
                joint.mrows = 1
                joint.nbody = 1

            elif joint.type == 'disc':
                joint.mrows = 2
                joint.nbody = 1

            elif joint.type == 'rigid':
                joint.mrows = 3
                joint.nbody = 2
                Bi = joint.iBindex
                Bj = joint.jBindex
                if Bi == 0:
                    joint.d0 = -self.Bodies[Bj].A.T @ self.Bodies[Bj].r
                    joint.p0 = -self.Bodies[Bj].p
                elif Bj == 0:
                    joint.d0 = self.Bodies[Bi].r
                    joint.p0 = self.Bodies[Bi].p
                else:
                    joint.d0 = self.Bodies[Bj].A.T @ (self.Bodies[Bi].r - self.Bodies[Bj].r)
                    joint.p0 = self.Bodies[Bi].p - self.Bodies[Bj].p
            else:
                raise ValueError("Joint type doesn't supported!")
        ##### ##### ##### ##### #####

        # functions
        ##### ##### ##### ##### #####
        if self.Functs:
            nFc = len(self.Functs)
            for Ci in range(nFc):
                functData(Ci)
        else:
            pass
        ##### ##### ##### ##### #####

        # compute number of constraints and determine row/column pointers
        ##### ##### ##### ##### #####
        nConst = 0
        for Ji in range(nJ):
            joint = self.Joints[Ji]
            joint.rows = nConst + 1
            joint.rowe = nConst + joint.mrows
            nConst = joint.rowe
            Bi = joint.iBindex
            if Bi != 0:
                joint.colis = 3 * (Bi - 1) + 1
                joint.colie = 3 * Bi
            Bj = joint.jBindex
            if Bj != 0:
                joint.coljs = 3 * (Bj - 1) + 1
                joint.colje = 3 * Bj
        ##### ##### ##### ##### #####

    def __update_position(self):
        """
        Update position entities.
        """
        # update rotation matrix of the body
        nB = len(self.Bodies)
        for Bi in range(nB):
            body = self.Bodies[Bi]
            body.A = A_matrix(body.p)

        # compute sP = A * sP_prime; rP = r + sP
        nP = len(self.Points)
        for Pi in range(nP):
            point = self.Points[Pi]
            Bi = point.Bindex
            if Bi != 0:
                body = self.Bodies[Bi-1]
                point.sP = body.A @ point.sPlocal
                point.sP_r = s_rot(point.sP)
                point.rP = body.r + point.sP

        # compute u = A * up
        nU = len(self.uVectors)
        for Vi in range(nU):
            unit_vector = self.uVectors[Vi]
            Bi = unit_vector.Bindex
            if Bi != 0:
                body = self.Bodies[Bi - 1]
                unit_vector.u = body.A @ unit_vector.ulocal
                unit_vector.u_r = s_rot(unit_vector.u)

    def __update_velocity(self):
        """
        Compute sP_dot and rP_dot vectors and update velocity components.
        """
        for Pi, point in enumerate(self.Points):
            Bi = point.Bindex
            if Bi != 0:
                point.dsp = point.sP_r * self.Bodies[Bi-1].p_d
                point.drp = self.Bodies[Bi-1].r_d + point.dsp

        for Vi, uvector in enumerate(self.uVectors):
            Bi = uvector.Bindex
            if Bi != 0:
                uvector.u_d = uvector.u_r * self.Bodies[Bi-1].p_d
            
    def __compute_constraints(self): #! - To be completed for all joint type -
        nConst = self.Joints[-1].rowe
        phi = np.zeros([nConst, 1])

        for joint_data in self.Joints:
            joint_type = joint_data.type

            if joint_type == 'rev':
                Pi = joint_data.iPindex
                Pj = joint_data.jPindex
                iBindex = joint_data.iBindex
                jBindex = joint_data.jBindex
                
                # compute relative positions of the points
                rP_i = self.Points[Pi].rP
                rP_j = self.Points[Pj].rP
                f = rP_i - rP_j

                if joint_data.fix == 1:
                    if iBindex == 0:  # body i is ground (fixed)
                        f = np.append(f, (-self.Bodies[jBindex].p - joint_data.p0))
                    elif jBindex == 0:  # body j is ground (fixed)
                        f = np.append(f, (self.Bodies[iBindex].p - joint_data.p0))
                    else:
                        f = np.append(f, (self.Bodies[iBindex].p - self.Bodies[jBindex].p - joint_data.p0))

            elif joint_type == 'tran':
                Pi = joint_data.iPindex
                Pj = joint_data.jPindex

                uj_r = self.uVectors[joint_data.jUindex].u_r
                ui = self.uVectors[joint_data.iUindex].u
                d = self.Points[Pi].rP - self.Points[Pj].rP

                # compute constraint equations
                f = np.array([uj_r.T @ d, uj_r.T @ ui]).reshape(2,1)

                # additional constraint if fixed
                if joint_data.fix == 1:
                    f = np.append(f, (ui.T @ d - joint_data.p0) / 2).reshape(3,1)

            elif joint_type == 'rev-rev':
                pass
            elif joint_type == 'rev-tran':
                pass
            elif joint_type == 'rigid':
                pass
            elif joint_type == 'disc':
                pass
            elif joint_type == 'rel-rot':
                pass
            elif joint_type == 'rel-tran':
                pass
            else:
                raise ValueError(f"Joint type '{joint_type}' is not supported.")

            rs = joint_data.rows - 1
            re = joint_data.rowe
            phi[rs:re] = f
            
        return phi

    def __compute_jacobian(self): #! - To be completed for all joint type -
        """
        Calculate the Jacobian matrix D for the system constraints.

        Returns
        -------
        D (NDArray)
            The Jacobian matrix of shape (nConst, nB3).
        """
        nConst = self.Joints[-1].rowe
        nB3 = 3 * len(self.Bodies)
        D = np.zeros((nConst, nB3))
        
        for Ji in range(len(self.Joints)):
            joint = self.Joints[Ji]

            if joint.type == 'rev':
                Pi = joint.iPindex
                Pj = joint.jPindex

                Di = np.block([
                    [np.eye(2), self.Points[Pi].sP_r.reshape(2, 1)]
                ])
                Dj = np.block([
                    [-np.eye(2), -self.Points[Pj].sP_r.reshape(2, 1)]
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

            elif joint.type == 'tran':
                Pi = joint.iPindex
                Pj = joint.jPindex

                uj = self.uVectors[joint.jUindex].u
                uj_r = self.uVectors[joint.jUindex].u_r
                d = self.Points[Pi].rP - self.Points[Pj].rP

                Di = np.block([
                    [uj_r.T, np.dot(uj.T, self.Points[Pi].sP).reshape(1, 1)],
                    [np.array([0, 0, 1])]
                ])
                Dj = np.block([
                    [-uj_r.T, -np.dot(uj.T, (self.Points[Pj].sP + d)).reshape(1, 1)],
                    [np.array([0, 0, -1])]
                ])

                if joint.fix == 1:
                    Di = np.vstack([
                        Di,
                        [uj.T, np.dot(uj.T, self.Points[Pi].sP_r).reshape(1)]
                    ])
                    Dj = np.vstack([
                        Dj,
                        [-uj.T, -np.dot(uj.T, self.Points[Pj].sP_r).reshape(1)]
                    ])

            elif joint.type == 'rev-rev':
                pass
            elif joint.type == 'rev-tran':
                pass
            elif joint.type == 'rigid':
                pass
            elif joint.type == 'disc':
                pass
            elif joint.type == 'rel-rot':
                pass
            elif joint.type == 'rel-tran':
                pass
            else:
                raise ValueError(f"Joint type '{joint.type}' is not supported.")

            # row indices for the current joint in the Jacobian matrix
            rs = joint.rows - 1
            re = joint.rowe

            #! Possibile problema con gli indici, risolvere a monte, 
            #! quando attribuisco i valori a cis e cie !!!
            # column indices for body i 
            if joint.iBindex != 0:
                cis = joint.colis - 1
                cie = joint.colie
                D[rs:re, cis:cie] = Di

            # column indices for body j
            if joint.jBindex != 0:
                cjs = joint.coljs - 1
                cje = joint.colje
                D[rs:re, cjs:cje] = Dj

        return D

    #! ?? utilizzato solo in casi particolari di vincoli ??
    def __rhs_velocity(self):
        """
        Calculate the right-hand side velocity vector for the system 
        constraints.

        Returns
        -------
        rhsv : NDArray
            A column vector representing the right-hand side of the 
            velocity equations.
        """

        nConst = self.Joints[-1].rowe
        rhsv = np.zeros((nConst, 1))

        for Ji in range(len(self.Joints)):
            joint = self.Joints[Ji]

            if joint.type == 'rel-rot':
                f = self.V_rel_rot(joint)
                rhsv[joint.rows - 1:joint.rowe] = f
                
            elif joint.type == 'rel-tran':
                f = self.V_rel_tran(joint)
                rhsv[joint.rows - 1:joint.rowe] = f

        return rhsv

    def __rhs_acceleration(self):
        """
        Calculate the right-hand side acceleration vector for the system 
        constraints.

        Returns
        -------
        rhsa : NDArray
            A column vector representing the right-hand side of the 
            acceleration equations.
        """
        nConst = self.Joints[-1].rowe
        rhsa = np.zeros([nConst, 1])
        
        for Ji, joint in enumerate(self.Joints):
            if joint.type == "rev":
                Pi, Pj = joint.iPindex, joint.jPindex
                Bi, Bj = self.Points[Pi].Bindex, self.Points[Pj].Bindex

                if Bi == 0:
                    f = s_rot(self.Points[Pj].dsP) * self.Bodies[Bj - 1].p_d
                elif Bj == 0:
                    f = -s_rot(self.Points[Pi].dsP) * self.Bodies[Bi - 1].p_d
                else:
                    f = -s_rot(self.Points[Pi].dsP) * self.Bodies[Bi - 1].p_d 
                    + s_rot(self.Points[Pj].dsP) * self.Bodies[Bj - 1].p_d

                if joint.fix == 1:
                    f = np.vstack([f, [0]])
            
            elif joint.type == "tran":
                Bi, Bj = joint.iBindex, joint.jBindex
                Pi, Pj = joint.iPindex, joint.jPindex
                ujd = self.uVectors[joint.jUindex].u_d
                ujd_r = s_rot(ujd)

                if Bi == 0:
                    f2 = 0
                elif Bj == 0:
                    f2 = 0
                else:
                    r_diff = self.Bodies[Bi].r - self.Bodies[Bj].r
                    p_d_product = np.dot(ujd.T, r_diff) * self.Bodies[Bi].p_d
                    r_d_diff = self.Bodies[Bi].r_d - self.Bodies[Bj].r_d
                    f2 = p_d_product - 2 * np.dot(ujd_r.T, r_d_diff)

                f = np.array([[f2], [0]])

                if joint.fix == 1:
                    d = self.Points[Pi].rP - self.Points[Pj].rP
                    d_d = self.Points[Pi].rP_d - self.Points[Pj].rP_d
                    L = joint.p0
                    u = d / L
                    u_d = d_d / L
                    f3 = -np.dot(u_d.T, d_d)

                    if Bi == 0:
                        f3 += np.dot(u.T, self.s_rot(self.Points[Pj].sP_d) * self.Bodies[Bj].p_d)
                    elif Bj == 0:
                        f3 -= np.dot(u.T, self.s_rot(self.Points[Pi].sP_d) * self.Bodies[Bi].p_d)
                    else:
                        term1 = self.Points[Pi].sP_d * self.Bodies[Bi].p_d
                        term2 = self.Points[Pj].sP_d * self.Bodies[Bj].p_d
                        f3 -= np.dot(u.T, self.s_rot(term1 - term2))

                    f = np.vstack([f, [f3]])
                
            elif joint.type == "rev-rev": 
                pass
            elif joint.type == "rev-tran": 
                pass
            elif joint.type == "rigid": 
                pass
            elif joint.type == "disc": 
                pass
            elif joint.type == "rel-rot": 
                pass
            elif joint.type == "rel-tran": 
                pass
        
            #! -1 should be fixed during the constraint definition, to avoid
            #! this continuos managing od indexis for slicing operations
            rs = joint.rows - 1 
            re = rs + 2 
            rhsa[rs:re] = f
        
        return rhsa

    def __bodies2u(self): 
        """ 
        Pack coordinates and velocities into the u array.
        """
        nB = len(self.Bodies)
        u = np.zeros([(3 * nB * 2), 1])

        for Bi in range(nB):
            ir = self.Bodies[Bi].irc - 1
            ird = self.Bodies[Bi].irv - 1
            u[ir:ir+3] = np.block([[self.Bodies[Bi].r],[self.Bodies[Bi].p]])
            u[ird:ird+3] = np.block([[self.Bodies[Bi].r_d], [self.Bodies[Bi].p_d]])
        
        return u

    def __bodies2ud(self):
        """ 
        Pack velocities and accelerations into ud. 
        """
        nB6 = 6 * len(self.Bodies)
        ud = np.zeros([nB6, 1])

        for Bi, body in enumerate(self.Bodies):
            ir = body.irc - 1
            ird = body.irv - 1
            ud[ir:ir + 3] = np.vstack([body.r_d, body.p_d]).reshape(3, 1)
            ud[ird:ird + 3] = np.vstack([body.r_dd, body.p_dd]).reshape(3, 1)

        return ud
    
    def __u2bodies(self, u):
        """
        Unpack u into coordinate and velocity sub-arrays.
        """ 
        #* for debugging purpose
        pdb.set_trace()
        
        nB = len(self.Bodies)
        for Bi in range(nB): 
            ir = self.Bodies[Bi].irc - 1
            ird = self. Bodies[Bi].irv - 1
            self.Bodies[Bi].r   = u[ir:ir+2] #! attention to the index
            self.Bodies[Bi].p   = u[ir+2][0]
            self.Bodies[Bi].r_d = u[ird:ird+2]
            self.Bodies[Bi].p_d = u[ird+2][0]

    def __compute_force(self, t):
        """
        Compute and return the array of forces acting on the system at time t.
        """
        #! !!! THIS PART OF CODE NEED TO BE COMPELTED !!!
        #! is it necessary ??
        for body in self.Bodies:
            body.f = colvect([0.0, 0.0]) # initialize body force vectors
            body.n = 0.0                 # initialize torque (moment) scalar

        #! the below code used to build the force array need to be optimized, 
        #! the class Force should keep inside all the method required to build
        #! the specific type of force -> code easier and more readable.
        # loop over all forces and apply them to the appropriate bodies
        for Fi, force in enumerate(self.Forces):
            if force.type == 'weight':      # apply weight to each body
                for body in self.Bodies: 
                    body.f += body.wgt

            elif force.type == 'ptp':       # call method or function for point-to-point force
                Pi, Pj = force.iPindex, force.jPindex
                Bi, Bj = force.iBindex, force.jBindex
                d = self.Points[Pi].rP - self.Points[Pj].rP
                dd = self.Points[Pi].drP - self.Points[Pj].drP
                L = np.sqrt(d.T@d)
                dL = d.T@dd/L
                delta = L - self.Forces[Fi].L0
                u = d/L 

                f = self.Forces[Fi].k*delta + self.Forces[Fi].dc*dL + self.Forces[Fi].f_a
                fi = f*u

                if Bi != 0:
                    self.Bodies[Bi-1].f = self.Bodies[Bi-1].f - fi
                    self.Bodies[Bi-1].n = self.Bodies[Bi-1].n - ((self.Points[Pi].sP_r).T@fi).item() # extracting the single value
                
                if Bj != 0: 
                    self.Bodies[Bj-1].f = self.Bodies[Bj-1].f + fi
                    self.Bodies[Bj-1].n = self.Bodies[Bj-1].n + ((self.Points[Pj].sP_r).T@fi).item() # extracting the single value

            elif force.type == 'rot-sda':   # call method or function for rotational SDA force
                # self.SDA_rot()
                pass
            elif force.type == 'flocal':    # apply local force to the specified body
                # Bi = force.iBindex
                # self.Bodies[Bi].f += self.Bodies[Bi].A @ force.flocal
                pass
            elif force.type == 'f':         # apply a global force to the specified body
                # Bi = force.iBindex
                # self.Bodies[Bi].f += force.f
                pass
            elif force.type == 'T':         # apply a torque to the specified body
                # Bi = force.iBindex
                # self.Bodies[Bi].n += force.T
                pass
            elif force.type == 'user':      # call a user-defined force function
                # self.user_force()
                pass

        nB3 = 3 * len(self.Bodies)
        g = np.zeros([nB3, 1])
        for Bi, body in enumerate(self.Bodies):
            ks = body.irc - 1
            ke = ks + 3
            g[ks:ke] = np.vstack([body.f, body.n])

        return g

    def __ic_correct(self):
        """
        Corrects initial conditions on the body coordinates and velocities.
        """
        flag = False

        # position correction
        for _ in range(20): #! 20 is an arbitrary value ... could be a parameter!
            self.__update_position()            # update position entities
            Phi = self.__compute_constraints()  # evaluate constraints
            D = self.__compute_jacobian()       # evaluate Jacobian
            ff = np.sqrt(np.dot(Phi.T, Phi))    # are the constraints violated?

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
            Phi[ir:ir + 2] = self.Bodies[Bi].r_d
            Phi[ir + 2] = self.Bodies[Bi].p_d

        rhsv = self.__rhs_velocity()  
        
        # solve for corrections
        delta_v = -D.T @ np.linalg.solve(D @ D.T, D @ Phi - rhsv)  

        # move corrected velocities to sub-arrays
        for Bi in range(nB):
            ir = 3 * Bi
            self.Bodies[Bi].r_d = self.Bodies[Bi].r_d + delta_v[ir:ir + 2]
            self.Bodies[Bi].p_d = self.Bodies[Bi].p_d + delta_v[ir + 2][0] # [0] because I need to extract the single value

        coords = np.zeros((nB, 3))
        vels = np.zeros((nB, 3))
        for Bi in range(nB):
            coords[Bi, :] = np.hstack((self.Bodies[Bi].r.T, np.array(self.Bodies[Bi].p).reshape(-1, 1)))
            vels[Bi, :] = np.hstack((self.Bodies[Bi].r_d.T, np.array(self.Bodies[Bi].p_d).reshape(-1, 1)))

        # print("\nCorrected coordinates")
        # print(" x           y           phi")
        # print(coords)
        # print("Corrected velocities")
        # print(" x-dot       y-dot       phi-dot")
        # print(vels)
        # print()

        #! orded print of the correction
        print("\nCorrected coordinates")
        print(" x           y           phi")
        for row in coords:
            print(f"{row[0]:<12.5f}{row[1]:<12.5f}{row[2]:<12.5f}")
        print("\nCorrected velocities")
        print(" x-dot       y-dot       phi-dot")
        for row in vels:
            print(f"{row[0]:<12.5f}{row[1]:<12.5f}{row[2]:<12.5f}")

    def __analysis(self, t, u):
        """
        Solve the constrained equations of motion at time t with the standard
        Lagrange multiplier method.
        """
        nB3 = 3 * len(self.Bodies)
        nConst = self.Joints[-1].rowe
        self.__u2bodies(u)  # unpack u into coordinate and velocity sub-arrays
        self.__update_position()
        self.__update_velocity()
        h_a = self.__compute_force(t)  # array of applied forces

        if nConst == 0:
            c_dd = self.M_inv_array * h_a # solve for accelerations
        else:
            D = self.__compute_jacobian()
            rhsA = self.__rhs_acceleration()  # right-hand side of acceleration constraints (gamma)

            # construct the matrix system to solve
            #! devo estrarre i singoli valori dal vettore M_array perchè
            #! è un vettore di array, da modificare sopra sulla creazione 
            #! del vettore stesso !!!
            DMD = np.block([
                [np.diag([arr[0] for arr in self.M_array]), -D.T], 
                [D, np.zeros([nConst, nConst])]
            ])
            rhs = np.concatenate([h_a, rhsA])

            #* check on conditioned index of the coefficient matrix
            cond_number = np.linalg.cond(DMD)
            if cond_number > 1e12:
                print(f"Warning: DMD matrix is poorly conditioned with condition number {cond_number}")
    
            # solve the system of equations
            sol = np.linalg.solve(DMD, rhs)
            c_dd = sol[:nB3]
            Lambda = sol[nB3:]

        # update accelerations for each body
        for Bi, body in enumerate(self.Bodies):
            ir = body.irc - 1
            i2 = ir + 2
            i3 = i2
            body.r_dd = c_dd[ir:i2]
            body.p_dd = c_dd[i3]

        ud = self.__bodies2ud()             # pack velocities and accelerations into ud
        self.__num += 1                     # increment the number of function evaluations

        if self.__showtime == 1:
            if self.__t10 % 100 == 0:
                #* print to check - only for debuggin purpose
                print(f"Time: {t}\n")

                print("Positions:")
                for i, body in enumerate(self.Bodies):
                    print(f"Body {i+1}")
                    for row in np.atleast_2d(body.r):
                        for value in row:
                            print(f"    [{value:>10.6f}]") 
                
                print("\nVelocities:")
                for i, body in enumerate(self.Bodies):
                    print(f"Body {i+1}")
                    for row in np.atleast_2d(body.r_d):
                        for value in row:
                            print(f"    [{value:>10.6f}]")
                
                print("\nAccelerations:")
                for i, body in enumerate(self.Bodies):
                    print(f"Body {i+1}")
                    for row in np.atleast_2d(body.r_dd):
                        for value in row:
                            print(f"    [{value:>10.6f}]")
        return ud

    def solve(self):
        """
        Solve the EQMs of the planar multi-body system.
        """
        # initial conditions and Jacobian matrix definition
        nConst = self.Joints[-1].rowe
        nB = len(self.Bodies)
        nB6 = 6 * nB
        ans = input("Do you want to correct the initial conditions? [(y)es/(n)o] ").lower()

        if nConst != 0:
            if ans == 'y':
                self.__ic_correct()
            D = self.__compute_jacobian()
            redund = np.linalg.matrix_rank(D) # check the rank of D for redundancy
            if redund < nConst:
                print("Redundancy in the constraints")

        # pack coordinates and velocities ito u array
        u = self.__bodies2u()
        
        #* check on u array
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            raise ValueError("Initial conditions contain NaN or Inf values.")

        t_initial = 0
        t_final = float(input("Final time = ? "))

        # utils to check the convergence
        self.__showtime = 1
        self.__num = 0
        self.__t10 = 0
        
        if t_final == 0:
            self.__analysis(0, u)
            T = 0
            uT = u.T
        else:
            dt = float(input("Reporting time-step = ? "))
            Tspan = np.arange(t_initial, t_final + dt, dt)

            usol = u.flatten()
            options = {'rtol': 1e-6, 'atol': 1e-9}
            sol = solve_ivp(self.__analysis, [t_initial, t_final], usol, t_eval=Tspan, method='DOP853', **options)
            T = sol.t
            uT = sol.y.T

        num_evals = self._PlanarDynamicModel__num
        print(f"Number of function evaluations = {num_evals}")
        print(f"Simulation completed!")
        
        return T, uT