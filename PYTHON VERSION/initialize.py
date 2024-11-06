# initialize
"""
Created on Fri Sep 29 11:03:00 2023

@author: simone.lucertini
"""

import DAP_BC_Functions as Functions
import numpy as np
import inData
Bodies = inData.Bodies
Points = inData.Points
Forces = inData.Forces
Joints = inData.Joints
Functs = inData.Functs
Points_anim = inData.Points_anim
Uvectors = inData.Uvectors

bodycolor = ["r", "g", "b", "c", "m"]

num = 0 # number of function evaluations
t10 = 0
flags = np.zeros((10, 1))
pen_d0 = np.zeros((10, 1))
    
#%% BODIES
nB = len(Bodies) # retrieves the QTY of BODIES
nB3 = 3*nB
nB6 = 6*nB


for Bi in range(nB):
    Bodies[Bi]["irc"]    = 3*Bi+1 # index of the 1st element of r in u or u_dot     ??? SHOULD BE CORRECTED FOR PYTHON???
    Bodies[Bi]["irv"]    = nB3 + 3*Bi + 1 # index of the 1st element of r_dot in u or u_dot  ??? SHOULD BE CORRECTED FOR PYTHON???
    Bodies[Bi]["m_inv"]  = 1/Bodies[Bi]["m"] # mass inverse
    Bodies[Bi]["J_inv"]  = 1/Bodies[Bi]["J"] # inverse of moment of inertia
    Bodies[Bi]["A"]      = Functions.Matrix_A(Bodies[Bi]["p"]) #rotational transformation matrix
#         Bodies(Bi).color  = bodycolor{mod(Bi-1, 5) + 1}; DISATTIVATO ADLL"ORIGINE !!!?????????????????????

#%% MASS(INERTIA) MATRIX MATRIX
                                      
M_array = np.zeros((nB3, 1)) # Init of mass matrix [M]
M_inv_array = np.zeros((nB3, 1)) # Init ofinverse mass matrix

for Bi in range(nB):
    is_ = 3*Bi
    
    M_array[is_ + 0] = Bodies[Bi]["m"]
    M_array[is_ + 1] = Bodies[Bi]["m"]
    M_array[is_ + 2] = Bodies[Bi]["J"]
    M_inv_array [is_ + 0] = Bodies[Bi]["m_inv"]
    M_inv_array [is_ + 1] = Bodies[Bi]["m_inv"]
    M_inv_array [is_ + 2] = Bodies[Bi]["J_inv"]
    


#%% POINTS

nP = len(Points)
nPanim = len(Points_anim)
nPtot = nP + nPanim
Points = Points + Points_anim

for Pi in range(nPtot):
    if Points[Pi]["Bindex"] == 0: # if Bindex == 0 this is the Body "0" so the point is the "ground".
        Points[Pi]["sP"] = Points[Pi]["sPlocal"]
        Points[Pi]["sP_r"] = Functions.s_rot((Points[Pi]["sP"]))
        Points[Pi]["rP"] = Points[Pi]["sP"]
    
    for Bi in range(nB):
        if Points[Pi]["Bindex"] == Bi:
            #nun serve? len_ = len(Bodies[Bi]["pts"])       !!! questo non dovrebbe servire aavendo messo append dopo
            Bodies[Bi]["pts"].append(Pi)            #  !!!! CAPIRE SE VANNO BENE GLI INDICI DEI PUNTI O VANNO INCREMENTATI DI 1
            

#%% UNIT VECTORS
# Unit vectors
nU = len(Uvectors)
for Vi in range(nU):
    if Uvectors[Vi]["Bindex"] == 0:
        Uvectors[Vi]["u"] = Uvectors[Vi]["ulocal"]
        Uvectors[Vi]["u_r"] = Functions.s_rot(Uvectors[Vi]["u"])

#%% FORCE ELEMS
nF = len(Forces)
for Fi in range(nF):
    force_type = Forces[Fi]["type"]
    
    if force_type == "weight": # Force type = weight
        ug = [elem * Forces[Fi]["gravity"] for elem in Forces[Fi]["wgt"]]
        for Bi in range(nB):
            Bodies[Bi]["wgt"] = [elem * Bodies[Bi]["m"] for elem in ug] 
            
    elif force_type == "ptp":  # Force type = point-to-point (i.e. spring)
        Pi = Forces[Fi]["iPindex"]
        Pj = Forces[Fi]["jPindex"]
        
        Forces[Fi]["iBindex"] = Points[Pi-1]["Bindex"]
        Forces[Fi]["jBindex"] = Points[Pj-1]["Bindex"]


#%% JOINTS
# Joints
nJ = len(Joints)
cfriction = 0

for Ji in range(nJ):
    joint_type = Joints[Ji]["type"]
    
    if joint_type == "rev": #if REVOLUTE JOINT....
        Joints[Ji]["mrows"] = 2
        Joints[Ji]["nbody"] = 2
        Pi = Joints[Ji]["iPindex"]
        Pj = Joints[Ji]["jPindex"] 
        Bi = Points[Pi-1]["Bindex"]
        Bj = Points[Pj-1]["Bindex"]
        Joints[Ji]["iBindex"] = Bi
        Joints[Ji]["jBindex"] = Bj
        
        if Joints[Ji]["fix"] == 1: 
            Joints[Ji]["mrows"] = 3
            if Bi == 0:
                Joints[Ji]["p0"] = -Bodies[Bj]["p"]
            elif Bj == 0:
                Joints[Ji]["p0"] = Bodies[Bi]["p"]
            else:
                Joints[Ji]["p0"] = Bodies[Bi]["p"] - Bodies[Bj]["p"]
                
    elif joint_type == "tran": #if SLIDER JOINT....
        Joints[Ji]["mrows"] = 2
        Joints[Ji]["nbody"] = 2
        Pi = Joints[Ji]["iPindex"]
        Pj = Joints[Ji]["jPindex"]
        Bi = Points[Pi-1]["Bindex"]
        Joints[Ji]["iBindex"] = Bi
        Bj = Points[Pj-1]["Bindex"]
        Joints[Ji]["jBindex"] = Bj
        
        if Joints[Ji]["fix"] == 1:
            Joints[Ji]["mrows"] = 3
            if Bi == 0:
                Joints[Ji]["p0"] = np.linalg.norm(Points[Pi]["rP"] - Bodies[Bj]["r"] - Bodies[Bj]["A"] @ Points[Pj]["sPlocal"])
            elif Bj == 0:
                Joints[Ji]["p0"] = np.linalg.norm(Bodies[Bi]["r"] + Bodies[Bi]["A"] @ Points[Pi]["sPlocal"] - Points[Pj]["rP"])
            else:
                Joints[Ji]["p0"] = np.linalg.norm(Bodies[Bi]["r"] + Bodies[Bi]["A"] @ Points[Pi]["sPlocal"] - Bodies[Bj]["r"] - Bodies[Bj]["A"] @ Points[Pj]["sPlocal"])
    
    elif joint_type == "rev-rev": # if REVOLUTE_REVOLUTE ....
        Joints[Ji]["mrows"] = 1
        Joints[Ji]["nbody"] = 2
        Pi = Joints[Ji]["iPindex"]
        Pj = Joints[Ji]["jPindex"]
        Joints[Ji]["iBindex"] = Points[Pi-1]["Bindex"]
        Joints[Ji]["jBindex"] = Points[Pj-1]["Bindex"]
    
    elif joint_type == "rev-tran": # if REVOLUTE-SLIDER JOINT....
        Joints[Ji]["mrows"] = 1
        Joints[Ji]["nbody"] = 2
        Pi = Joints[Ji]["iPindex"]
        Pj = Joints[Ji]["jPindex"]
        Joints[Ji]["iBindex"] = Points[Pi-1]["Bindex"]
        Joints[Ji]["jBindex"] = Points[Pj-1]["Bindex"]
    
    elif joint_type in ["rel-rot", "rel-tran"]: # if relative DISPLACEMENT....
        Joints[Ji]["mrows"] = 1
        Joints[Ji]["nbody"] = 1
    
    elif joint_type == "disc": # if DISC....
        Joints[Ji]["mrows"] = 2
        Joints[Ji]["nbody"] = 1
    
    elif joint_type == "rigid": # if RIGID....
        Joints[Ji]["mrows"] = 3
        Joints[Ji]["nbody"] = 2
        Bi = Joints[Ji]["iBindex"]
        Bj = Joints[Ji]["jBindex"]
        
        if Bi == 0:
            Joints[Ji]["d0"] = -Bodies[Bj]["A"].T @ Bodies[Bj]["r"] # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Joints[Ji]["p0"] = -Bodies[Bj]["p"]
        elif Bj == 0:
            Joints[Ji]["d0"] = Bodies[Bi]["r"]
            Joints[Ji]["p0"] = Bodies[Bi]["p"]
        else:
            Joints[Ji]["d0"] = Bodies[Bj]["A"].T @ (Bodies[Bi]["r"] - Bodies[Bj]["r"])  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Joints[Ji]["p0"] = Bodies[Bi]["p"] - Bodies[Bj]["p"]
    
    else:
        print("ERROR: UNDEFINED JOINT TYPE!")

#%% FUNCTIONS
def functData(Ci): #mi sa che vanno convertti tutti in Dict es Functs[Ci].type -> Functs[Ci]["type"]
                                # !!! QUI PROBABILMENTE VA CORRETTA UN SACCO DI ROBA 
    if Functs[Ci].type == 'a':
        Functs[Ci].ncoeff = 4
        Functs[Ci].coeff[3] = 2 * Functs[Ci].coeff[2]
    elif Functs[Ci].type == 'b':
        Functs[Ci].ncoeff = 9
        xe = Functs[Ci].t_end - Functs[Ci].t_start
        fe = Functs[Ci].f_end - Functs[Ci].f_start
        C = [[xe**3, xe**4, xe**5],
             [3*xe**2, 4*xe**3, 5*xe**4],
             [6*xe, 12*xe**2, 20*xe**3]]
        sol = np.linalg.solve(C, [fe, 0, 0])
        Functs[Ci].coeff[0:3] = sol.transpose()
        Functs[Ci].coeff[3] = 3 * sol[0]
        Functs[Ci].coeff[4] = 4 * sol[1]
        Functs[Ci].coeff[5] = 5 * sol[2]
        Functs[Ci].coeff[6] = 6 * sol[0]
        Functs[Ci].coeff[7] = 12 * sol[1]
        Functs[Ci].coeff[8] = 20 * sol[2]
    elif Functs[Ci].type == 'c':
        Functs[Ci].ncoeff = 9
        xe = Functs[Ci].t_end - Functs[Ci].t_start
        fpe = Functs[Ci].dfdt_end
        C = [[4*xe**3, 5*xe**4, 6*xe**5],
             [12*xe**2, 20*xe**3, 30*xe**4],
             [24*xe, 60*xe**2, 120*xe**3]]
        sol = np.linalg.solve(C, [fpe, 0, 0])
        Functs[Ci].coeff[0:3] = sol.transpose()
        Functs[Ci].coeff[3] = 4 * sol[0]
        Functs[Ci].coeff[4] = 5 * sol[1]
        Functs[Ci].coeff[5] = 6 * sol[2]
        Functs[Ci].coeff[6] = 12 * sol[0]
        Functs[Ci].coeff[7] = 20 * sol[1]
        Functs[Ci].coeff[8] = 30 * sol[2]


nFc = len(Functs)
for Ci in range(nFc):
    functData(Ci)

# %% ------------------------------------------------------------
# Compute number of constraints and determine row/column pointers

nConst = 0
for Ji in range(nJ):
    Joints[Ji]["rows"] = nConst + 1
    Joints[Ji]["rowe"] = nConst + Joints[Ji]["mrows"]
    nConst = Joints[Ji]["rowe"]
    Bi = Joints[Ji]["iBindex"]
    if Bi != 0:
        Joints[Ji]["colis"] = 3*(Bi - 1) + 1
        Joints[Ji]["colie"] = 3*Bi
    Bj = Joints[Ji]["jBindex"]
    if Bj != 0:
        Joints[Ji]["coljs"] = 3*(Bj - 1) + 1
        Joints[Ji]["colje"] = 3*Bj

def init(): # Vars to be shared with parent scripts
    global nConst