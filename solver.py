"""
This Python module provides the necessary algorithms and functions to solve 
planar multi-body dynamic models.

Author: Giacomo Cangi
"""

import numpy as np
import scipy as sc
from pmd_functions import *


class PmdDynamicModel:
    def __init__(self, Bodies, Points, Uvectors, Forces, Joints):
        self.Bodies = Bodies
        self.Points = Points
        self.Uvectors = Uvectors
        self.Forces = Forces
        self.Joints = Joints

        # initialize the model
        self.initialize()

    def initialize(self):
        # initialize variables
        nB = len(self.Bodies)
        nB3 = 3 * nB
        nB6 = 6 * nB

        for Bi in range(nB):
            body = self.Bodies[Bi]
            body.irc = 3 * Bi + 1
            body.irv = nB3 + 3 * Bi + 1
            body.m_inv = 1 / body.m
            body.J_inv = 1 / body.J
            body.A = A_matrix(body.p)

        # Mass (inertia) matrix as an array
        M_array = np.zeros((nB3, 1))
        M_inv_array = np.zeros((nB3, 1))
        for Bi in range(nB):
            is_ = 3 * Bi
            M_array[is_:is_ + 3] = np.array([[self.Bodies[Bi].m], [self.Bodies[Bi].m], [self.Bodies[Bi].J]])
            M_inv_array[is_:is_ + 3] = np.array([[self.Bodies[Bi].m_inv], [self.Bodies[Bi].m_inv], [self.Bodies[Bi].J_inv]])

        # Points
        nPtot = len(self.Points)
        for Pi in range(nPtot):
            point = self.Points[Pi]
            if point.Bindex == 0:
                point.sP = point.sPlocal
                point.sP_rotated = s_rot(point.sP)
                point.rP = point.sP

            for Bi in range(nB):
                if point.Bindex == Bi:
                    self.Bodies[Bi].pts.append(Pi)  # Append point index to the body's points

        # Unit vectors  
        nU = len(self.Uvectors)
        for Vi in range(nU):
            unit_vector = self.Uvectors[Vi]
            if unit_vector.Bindex == 0:
                unit_vector.u = unit_vector.ulocal
                unit_vector.u_r = s_rot(unit_vector.u)

        # Force elements
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

        # Joints
        nJ = len(self.Joints)
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
                        joint.p0 = np.linalg.norm(self.Points[Pi].rP - self.Bodies[joint.jBindex].r - self.Bodies[joint.jBindex].A @ self.Points[Pj].sPlocal)
                    elif joint.jBindex == 0:
                        joint.p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + self.Bodies[joint.iBindex].A @ self.Points[Pi].sPlocal - self.Points[Pj].rP)
                    else:
                        joint.p0 = np.linalg.norm(self.Bodies[joint.iBindex].r + self.Bodies[joint.iBindex].A @ self.Points[Pi].sPlocal - self.Bodies[joint.jBindex].r - self.Bodies[joint.jBindex].A @ self.Points[Pj].sPlocal)

# Usage example:
# model = DynamicModel(Bodies, Points, Points_anim, Uvectors, Forces, Joints)
