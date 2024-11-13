"""
This Python module provides the necessary algorithms and functions to solve 
planar multi-body dynamic models.

Author: Giacomo Cangi
"""

import numpy as np
import scipy as sc
import numpy.linalg as lng
from PMD.pmd_functions import *


class PmdDynamicModel:
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
        M_array = np.zeros((nB3, 1))
        M_inv_array = np.zeros((nB3, 1))
        for Bi in range(nB):
            is_ = 3 * Bi
            ie_ = is_ + 3
            M_array[is_:ie_] = np.array([[self.Bodies[Bi].m], [self.Bodies[Bi].m], [self.Bodies[Bi].J]])
            M_inv_array[is_:ie_] = np.array([[self.Bodies[Bi].m_inv], [self.Bodies[Bi].m_inv], [self.Bodies[Bi].J_inv]])
        ##### ##### ##### ##### ##### 

        # points #! CHECK CON MATLAB - Probabile problema di indicizzazione
        ##### ##### ##### ##### #####
        nPtot = len(self.Points)
        for Pi in range(nPtot):
            point = self.Points[Pi]
            if point.Bindex == 0:
                point.sP = point.sPlocal
                point.sP_rotated = s_rot(point.sP)
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
                point.sP_rotated = s_rot(point.sP)
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

    def __jacobian(self): #! - To be completed for all joint type -
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
                    [np.eye(2), self.Points[Pi].sP_rotated.reshape(2, 1)]
                ])
                Dj = np.block([
                    [-np.eye(2), -self.Points[Pj].sP_rotated.reshape(2, 1)]
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
                        [uj.T, np.dot(uj.T, self.Points[Pi].sP_rotated).reshape(1)]
                    ])
                    Dj = np.vstack([
                        Dj,
                        [-uj.T, -np.dot(uj.T, self.Points[Pj].sP_rotated).reshape(1)]
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

    def __rhs_velocity(self):
        """
        Calculate the right-hand side velocity vector for the system constraints.

        Returns
        -------
        rhs : NDArray
            A column vector representing the right-hand side of the velocity equations.
        """

        nConst = self.Joints[-1].rowe
        rhs = np.zeros((nConst, 1))

        for Ji in range(len(self.Joints)):
            joint = self.Joints[Ji]

            if joint.type == 'rel-rot':
                f = self.V_rel_rot(joint)
                rhs[joint.rows - 1:joint.rowe] = f
                
            elif joint.type == 'rel-tran':
                f = self.V_rel_tran(joint)
                rhs[joint.rows - 1:joint.rowe] = f

        return rhs

    def __ic_correct(self):
        """
        Corrects initial conditions on the body coordinates and velocities.
        """
        flag = False

        # position correction
        for _ in range(20): #! 20 is an arbitrary value ... could be a parameter!
            self.__update_position()            # update position entities
            Phi = self.__compute_constraints()  # evaluate constraints
            D = self.__jacobian()               # evaluate Jacobian
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

        #! orded print
        print("\nCorrected coordinates")
        print(" x           y           phi")
        for row in coords:
            print(f"{row[0]:<12.5f}{row[1]:<12.5f}{row[2]:<12.5f}")

        print("\nCorrected velocities")
        print(" x-dot       y-dot       phi-dot")
        for row in vels:
            print(f"{row[0]:<12.5f}{row[1]:<12.5f}{row[2]:<12.5f}")

    def solve(self):
        """
        Solve the EQMs of the planar multi-body system.
        """
        # initial conditions and Jacobian matrix definition
        nConst = self.Joints[-1].rowe
        ans = input("Do you want to correct the initial conditions? [(y)es/(n)o] ").lower()
        if nConst != 0:
            if ans == 'y':
                self.__ic_correct()
            D = self.__jacobian()
            redund = np.linalg.matrix_rank(D) # check the rank of D for redundancy
            if redund < nConst:
                print("Redundancy in the constraints")