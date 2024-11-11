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

    def __ic_correct(self):
        """
        Corrects initial conditions on the body coordinates and velocities.
        """
        flag = False

        for _ in range(20): #! 20 is an arbitrary value ... could be a parameter!
            self.__update_position()            # update position entities
            Phi = self.__constraints()         # evaluate constraints
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
                self.Bodies[Bi].r += delta_c[ir:ir + 2]
                self.Bodies[Bi].p += delta_c[ir + 2]

        if not flag:
            raise ValueError("Convergence failed in Newton-Raphson!")

        # velocity correction
        nB = len(self.Bodies)
        Phi = np.zeros([3 * nB, 1])
        for Bi in range(nB):
            ir = 3 * Bi
            Phi[ir:ir + 2, 0] = self.Bodies[Bi].r_d
            Phi[ir + 2, 0] = self.Bodies[Bi].p_d

        rhs = self.__rhs_velocity(0)  # Compute rhs
        delta_v = -D.T @ np.linalg.solve(D @ D.T, D @ Phi - rhs)  # Compute corrections

        # Move corrected velocities to sub-arrays
        for Bi in range(nB):
            ir = 3 * Bi
            self.Bodies[Bi].r_d += delta_v[ir:ir + 2]
            self.Bodies[Bi].p_d += delta_v[ir + 2]

        # Report corrected coordinates and velocities
        coords = np.zeros((nB, 3))
        vels = np.zeros((nB, 3))
        for Bi in range(nB):
            coords[Bi, :] = np.hstack((self.Bodies[Bi].r, self.Bodies[Bi].p))
            vels[Bi, :] = np.hstack((self.Bodies[Bi].r_d, self.Bodies[Bi].p_d))

        print("\nCorrected coordinates")
        print(" x           y           phi")
        print(coords)
        print("Corrected velocities")
        print(" x-dot       y-dot       phi-dot")
        print(vels)
        print()
        
    def solve(self):
        """
        
        """
        # initial conditions and Jacobian matrix definition
        nConst = self.Joints[-1].rowe
        if nConst != 0:
            ans = input("Do you want to correct the initial conditions? [(y)es/(n)o] ").lower()
            if ans == 'y':
                self.__ic_correct()

            D = self.__jacobian
            redund = np.linalg.matrix_rank(D) # check the rank of D for redundancy

            if redund < nConst:
                print("Redundancy in the constraints")