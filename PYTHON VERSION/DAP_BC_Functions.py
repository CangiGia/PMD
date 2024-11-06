# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:43:24 2023

@author: Simone.Lucertini
"""

#%% FUNCTIONs FOR: ROTATIONS calcs

def Matrix_A(x): # This function computes the rotational transformation matrix A
    import numpy as np

    c = np.cos(x)
    s = np.sin(x)
    A = np.array(
                [[c, -s],
                 [s, c]])
    return A

def s_rot(s):
    import numpy as np
    # This function rotates a vector 90 degrees positively
    s_r = np.array([-s[1], s[0]])
    return s_r

#%%
# Bodies_to_u_d
# Pack velocities and accelerations into u_dot

""" WIP FROM HERE, to be checked!!!!!"""

def Bodies_to_u_d(nB6, nB, Bodies):
    import numpy as np
    u_d = np.zeros(nB6)
    for Bi in range(nB):
        ir = Bodies[Bi]["irc"]
        ird = Bodies[Bi]["irv"]
        u_d[ir:ir+3] = np.concatenate((Bodies[Bi]["r_d"], Bodies[Bi]["p_d"]))
        u_d[ird:ird+3] = np.concatenate((Bodies[Bi]["r_dd"], Bodies[Bi]["p_dd"]))

#%%
def RHSAcc(t, nJ, nConst, Bodies, Joints, Points, Uvectors):
    import numpy as np
    
    rhs = np.zeros(nConst)

    for Ji in range(nJ):
        joint_type = Joints[Ji]['type']

        if joint_type == "rev": #==================================================================
            pass
        elif joint_type == "tran":
            

            Bi = Joints[Ji]['iBindex']
            Bj = Joints[Ji]['jBindex']
            Pi = Joints[Ji]['iPindex']
            Pj = Joints[Ji]['jPindex']
            ujd = Uvectors[Joints[Ji]['jUindex']]['u_d']
            ujd_r = s_rot(ujd)
        
            if Bi == 0:
                f2 = 0
            elif Bj == 0:
                f2 = 0
            else:
                f2 = np.dot(ujd, (Bodies[Bi]['r'] - Bodies[Bj]['r'])) * Bodies[Bi]['p_d'] - 2 * np.dot(ujd_r, (Bodies[Bi]['r_d'] - Bodies[Bj]['r_d']))
            
            f = np.array((f2, 0))
        
            if Joints[Ji]['fix'] == 1:
                d = Points[Pi]['rP'] - Points[Pj]['rP']
                d_d = Points[Pi]['rP_d'] - Points[Pj]['rP_d']
                L = Joints[Ji]['p0']
                u = d / L
                u_d = d_d / L
                f3 = -np.dot(u_d, d_d)
        
                if Bi == 0:
                    f3 = f3 + np.dot(u, np.dot(s_rot(Points[Pj]['sP_d']), Bodies[Bj]['p_d']))
                elif Bj == 0:
                    f3 = f3 - np.dot(u, np.dot(s_rot(Points[Pi]['sP_d']), Bodies[Bi]['p_d']))
                else:
                    f3 = f3 - np.dot(u, (np.dot(s_rot(Points[Pi]['sP_d']), Bodies[Bi]['p_d']) - np.dot(Points[Pj]['sP_d'], Bodies[Bj]['p_d'])))
        
                f.extend([f3])


        elif joint_type == "rev-rev": #==================================================================
            
            # from file: A_rev_rev()

            Pi = Joints[Ji]['iPindex']
            Pj = Joints[Ji]['jPindex']
            Bi = Joints[Ji]['iBindex']
            Bj = Joints[Ji]['jBindex']
        
            d = Points[Pi]['rP'] - Points[Pj]['rP']
            d_d = Points[Pi]['rP_d'] - Points[Pj]['rP_d']
        
            L = Joints[Ji]['L']
            u = d / L
            u_d = d_d / L
        
            f = -np.dot(u_d, d_d)
        
            if Bi == 0:
                f = f + np.dot(u, s_rot(Points[Pj]["sP_d"])) * Bodies[Bj]["p_d"]
            elif Bj == 0:
                f = f - np.dot(u, s_rot(Points[Pi]["sP_d"])) * Bodies[Bi]["p_d"]
            else:
                f = f - np.dot(u, (s_rot(Points[Pi]["sP_d"]) * Bodies[Bi]["p_d"] - Points[Pj]["sP_d"] * Bodies[Bj]["p_d"]))

        elif joint_type == "rev-tran": #==================================================================
            #from file A_rev_tran()

            # r-h-s of acc. constraint for a revolute translational joint
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            Bi = Joints[Ji]["iBindex"]
            Bj = Joints[Ji]["jBindex"]
        
            ui = Uvectors[Joints[Ji]["iUindex"]]["u"]
            ui_d = Uvectors[Joints[Ji]["iUindex"]]["u_d"]
        
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            d_d = Points[Pi]["rP_d"] - Points[Pj]["rP_d"]
        
            if Bi == 0:
                f = np.dot(ui, Points[Pj]["sP_d"]) * Bodies[Bj]["p_d"]
            elif Bj == 0:
                f = np.dot(ui_d, (d * Bodies[Bi]["p_d"] + 2 * s_rot(d_d))) - np.dot(ui, Points[Pi]["sP_d"]) * Bodies[Bi]["p_d"]
            else:
                f = np.dot(ui_d, (d * Bodies[Bi]["p_d"] + 2 * s_rot(d_d))) - np.dot(ui, (Points[Pi]["sP_d"] * Bodies[Bi]["p_d"] - Points[Pj]["sP_d"] * Bodies[Bj]["p_d"]))

        elif joint_type == "rigid": #==================================================================
            #FROM A_rigid() 
            # r-h-s of acc. constraints for a Rigid joint
            Bj = Joints[Ji]["jBindex"]
            
            f = np.array([0, 0, 0])
            if Bj != 0:
                f = np.array([-Bodies[Bj]["A"] * Joints[Ji]["d0"] * Bodies[Bj]["p_d"]**2, 0])
                
        elif joint_type == "disc": #==================================================================
            # from A_disc()
            # The r-h-s of acc. constraints for a disk joint
            
            f = np.array([0, 0])    
            
        elif joint_type == "rel-rot": #==================================================================
            #f rom A_rel_rot()
            fun, fun_d, fun_dd = functs(Joints[Ji]["iFunct"], t)
            f = fun_dd
            
        elif joint_type == "rel-tran": #==================================================================
            # from A_rel_tran()
            # r-h-s of acc. constraint for relative-rotational constraint
            Pi = Joints[Ji]["iPindex"]
            Pj = Joints[Ji]["jPindex"]
            Bi = Joints[Ji]["iBindex"]
            Bj = Joints[Ji]["jBindex"]
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            d_d = Points[Pi]["rP_d"] - Points[Pj]["rP_d"]
            fun, fun_d, fun_dd = functs(Joints[Ji]["iFunct"], t)
            
            f = fun * fun_dd + fun_d**2
            if Bi == 0:
                f = f + np.dot(d, np.transpose(s_rot(Points[Pj]["sP_d"]))) * Bodies[Bj]["p_d"]
            elif Bj == 0:
                f = f - np.dot(d, np.transpose(s_rot(Points[Pi]["sP_d"]))) * Bodies[Bi]["p_d"] - np.dot(d_d, d_d)
            else:
                f = f + np.dot(d, np.transpose(s_rot(Points[Pj]["sP_d"]))) * Bodies[Bj]["p_d"] - np.dot(d, np.transpose(s_rot(Points[Pi]["sP_d"]))) * Bodies[Bi]["p_d"] - np.dot(d_d, d_d)            


        rs = Joints[Ji]["rows"]-1 # <---------- VERIFY INDEX
        re = Joints[Ji]["rowe"] # <---------- VERIFY INDEX
        rhs[rs:re] = f

    return rhs
