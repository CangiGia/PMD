# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 16:56:28 2023

@author: simone.lucertini
"""

"""
         !!!!!!!!!!!  PROBABLY IT IS BETTER TO INCLUDE THIS ONE IN INITIALIZE SO WE HAVE IT IN ONE PLACE?
"""


def Jacobian():
    import DAP_BC_Functions as Functions
    import numpy as np
    
    import initialize
    nB3 = initialize.nB3
    nConst = initialize.nConst
    Bodies = initialize.Bodies 
    Joints = initialize.Joints
    Points = initialize.Points
    Uvectors = initialize.Uvectors
    nJ = initialize.nJ
    #s_rot = initialize.s_rot
    
    D = np.zeros((nConst, nB3))
    for Ji in range(nJ):
        jointType = Joints[Ji]["type"]
        if jointType == "rev": # case of REVOLUTE JOINT [OK]
            
            # revolute_J from file: J_rev()
            # Jacobian sub-matrices for a revolute joint
            
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            
            Di = np.hstack((np.eye(2), np.array([Points[Pi-1]["sP_r"]]).T)) #<<<<<<< CONTROLLA GLI INDICI, FORSE VANNO BENE
            Dj = np.hstack((-np.eye(2), -np.array([Points[Pj-1]["sP_r"]]).T)) #<<<<<<< CONTROLLA GLI INDICI, FORSE VANNO BENE  
            
            if Joints[Ji]["fix"] == 1:
                Di = np.vstack((Di, [0, 0, 1]))
                Dj = np.vstack((Dj, [0, 0, -1]))
    
        elif jointType == "tran": # case of SLIDER JOINT [OK]
            
            # translational_J from file: J_tran()
            # Jacobian sub-matrices for a translational joint
            
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            i = Joints[Ji]["jUindex"]-1 #correct Python index for iUindex and jUindex TO BE VERIFIED!!!
            uj = np.array(Uvectors[i]["u"])   # was: uj = Uvectors[Joints[Ji]["jUindex"]]["u"]
            uj_r = np.array(Uvectors[i]["u_r"]) # was: uj_r = Uvectors[Joints[Ji]["jUindex"]]["u_r"]  
            del i
           
            d = np.array(Points[Pi-1]["rP"])- np.array(Points[Pj-1]["rP"]) 
            
            Di = np.vstack((np.array([[uj_r[0],  uj_r[1],  float(uj.T@np.array(Points[Pi-1]["sP"]))]])
                           , np.array([[0,  0,  +1]])))
            
            Dj = np.vstack((np.array([[-uj_r[0],  -uj_r[1],  -float(uj.T@(np.array(Points[Pj-1]["sP"])+d))]])
                           , np.array([[0,  0,  -1]])))                
            
    
            if Joints[Ji]["fix"] == 1:
                Di = np.vstack((Di, [uj.T, uj.T @ Points[Pi]["sP_r"]])) # !!!??? 
                Dj = np.vstack((Dj, [-uj.T, -uj.T @ Points[Pj]["sP_r"]])) # !!!??? 
    
        elif jointType == "rev-rev": # case of REVOLUTE-REVOLUTE (COMPOSITE) JOINT [!!!TO BE CHECKED!!!]
            
            # rev_rev_J from file J_rev_rev()
            # Jacobian sub-matrices for a revolute_revolute joint
            
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            L = Joints[Ji]["L"]
            u = d / L
            
            Di = np.array([u.T, u.T @ Points[Pi]["sP_r"]]) # !!!??? VEDI JOINT TRAN
            Dj = np.array([-u.T, -u.T @ Points[Pj]["sP_r"]]) # !!!??? VEDI JOINT TRAN
    
        elif jointType == "rev-tran": # case of REVOLUTE+SLIDER JOINT [!!!TO BE CHECKED!!!]
            
            # rev_tran_J from file: J_rev_tran()
            # Jacobian sub-matrices for a revolute-translational joint
            
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            ui = Uvectors[Joints[Ji]["iUindex"]]["u"]
            ui_r = Uvectors[Joints[Ji]["iUindex"]]["u_r"]
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            
            Di = np.array([[ui_r.T, ui.T @ (Points[Pi]["sP"] - d)]])
            Dj = np.array([[-ui_r.T, -ui.T @ Points[Pj]["sP"]]])
    
        elif jointType == "rigid": # case of RIGID JOINT [!!!TO BE CHECKED!!!]
            
            # rigid_J from file: J_rigid()
            # Jacobian sub-matrix for a Rigid joint
            
            Bj = Joints[Ji]["jBindex"]
            
            Di = np.eye(3)
            if Bj != 0:
                Dj = np.array([[-np.eye(2), -Functions.s_rot(Bodies[Bj]["A"] @ Joints[Ji]["d0"])],
                               [0, 0, -1]])
    
        elif jointType == "disc": # case of DISC JOINT [!!!TO BE CHECKED!!!]
            
            # disc_J, from file J_disc()
            # Jacobian sub-matrices for a disk joint between a body and the ground
            
            Di = np.array([[0, 1, 0],
                           [1, 0, Joints[Ji]["R"]]])
    
        elif jointType == "rel-rot": # case of RELATIVE ROTATION [!!!TO BE CHECKED!!!]
           
            # rel-rot_J from file :  J_rel_rot()
            # Jacobian sub-matrices for relative-rotational constraint 
            
            Di = np.array([0, 0, 1])
            Dj = np.array([0, 0, -1])
    
        elif jointType == "rel-tran": # case of RELATIVE TRANSLATION [!!!TO BE CHECKED!!!]
            
            # translational_J from file: J_rel_tran()
            # Jacobian sub-matrices for a translational joint
            
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            uj = Uvectors[Joints[Ji]["jUindex"]]["u"]
            uj_r = Uvectors[Joints[Ji]["jUindex"]]["u_r"]
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            
            Di = np.array([[uj_r.T, np.dot(uj.T, Points[Pi]["sP"])],
                           [0, 0, 1]])
            Dj = np.array([[-uj_r.T, -np.dot(uj.T, Points[Pj]["sP"] + d)],
                           [0, 0, -1]])
            
            if Joints[Ji]["fix"] == 1:
                Di = np.vstack((Di, [uj.T, np.dot(uj.T, Points[Pi]["sP_r"])]))
                Dj = np.vstack((Dj, [-uj.T, -np.dot(uj.T, Points[Pj]["sP_r"])]))
    
    # Buildind [D] Matrix!
        rs = Joints[Ji]["rows"]-1 #corrected to start from 0 
        re = Joints[Ji]["rowe"]
        
        if Joints[Ji]["iBindex"] != 0:
            cis = Joints[Ji]["colis"]-1 #corrected to start from 0 
            cie = Joints[Ji]["colie"]
            D[rs:re, cis:cie] = Di 
    
        if Joints[Ji]["jBindex"] != 0:
            cjs = Joints[Ji]["coljs"]-1 #corrected to start from 0 
            cje = Joints[Ji]["colje"]
            D[rs:re, cjs:cje] = Dj 
    
    return D
