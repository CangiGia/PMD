# =DAP
"""
Created on Fri Sep 29 16:37:20 2023

@author: simone.lucertini
"""
import DAP_BC_Functions as Functions
import numpy as np
from scipy.integrate import odeint #I import ODEint solver to integrate the equations of motion.

import initialize
nB3 = initialize.nB3
nB6 = initialize.nB6
nB = initialize.nB
nP = initialize.nP
nU = initialize.nU
nF = initialize.nF
num = initialize.num
Bodies = initialize.Bodies 
Joints = initialize.Joints
Points = initialize.Points
Forces = initialize.Forces
Uvectors = initialize.Uvectors
nJ = initialize.nJ
nConst = initialize.nConst
M_array = initialize.M_array
M_inv_array = initialize.M_inv_array


if nConst != 0:
    #ans = input("Do you want to correct the initial conditions? [(y)es/(n)o] ")
    ##### if ans == "y":   TEMPORARY DISABLED, RE-ENABLE PLZ!!!!!!
     ######   ic_correct()  # Correct initial conditions # XXXX SEE REFERENCED FILE--------------------
    import Jacobian                           # asassarebbe mejo portarlo qui dentro direttametne???  da decidere
    D = Jacobian.Jacobian() 
    
    redund = np.linalg.matrix_rank(D)  # Determine redundancy between constraints
    if redund < nConst:
        print("WARNING! Redundancy in the constraints")

# Pack coordinates and velocities into u array
#   (ri, pi); for i = 1:nB and (ri_dot, pi_dot); for i = 1:nB 
u = np.zeros(nB6)

#%% BODIES to "u"       # refer to file Bodies_to_u() that has been modev here

for Bi in range(nB):
    ir = Bodies[Bi]["irc"]-1 # <<<<<<<<<<<<<<<<<<<<<<<------------------to be double checked
    ird = Bodies[Bi]["irv"]-1 # <<<<<<<<<<<<<<<<<<<<<<<-----------------to be double checked
   
    u[ir:ir+2] = np.array((Bodies[Bi]["r"])) 
    u[ir+2] = Bodies[Bi]["p"]
    u[ird:ird+2] = np.array((Bodies[Bi]["r_d"])) 
    u[ird+2] = Bodies[Bi]["p_d"]
    
#%% ANALYSIS (from Analysis.m) ---------------- Anakysis.py moved here

# import analysis now included here
t_initial = 0 #set intial time to 0
t_final = 1 #float(input("Final time = ? "))  # ASK final time from user         #<<<<<<<<<<<<< ACTUALLY FORCED TO 1 TO SIMPLIFY DEBUG
showtime = 1




""" DEBUGGA DA QUI """



#%%  FUNCTIONs FOR INITIALIZATION OF BODY FORCE VECTORS (originally "force_array.m")

def Force_array(t):

    # initialize body force vectors
    for Bi in range(nB):
        Bodies[Bi]["f"] = np.array([0, 0])
        Bodies[Bi]["n"] = 0

    for Fi in range(nF): 
        force_type = Forces[Fi]["type"]

        if force_type == 'weight':
            for Bi in range(nB):
                Bodies[Bi]["f"] = Bodies[Bi]["f"] + Bodies[Bi]["wgt"]

        elif force_type == 'ptp': # Point-to-point spring-damper-actuator-------------------------------------------
            # PTP Spring-Damper ACTUATOR from file (SDA_ptp.m)-

            Pi = Forces[Fi]["iPindex"]
            Pj = Forces[Fi]["jPindex"]
            Bi = Forces[Fi]["iBindex"]
            Bj = Forces[Fi]["jBindex"]
            d = Points[Pi]["rP"] - Points[Pj]["rP"]
            d_dot = Points[Pi]["rP_d"] - Points[Pj]["rP_d"]
            L = np.sqrt(np.dot(d, d))
            L_dot = np.dot(d, d_dot) / L
            del_ = L - Forces[Fi]["L0"]
            u = d / L

            f = Forces[Fi]["k"] * del_ + Forces[Fi]["dc"] * L_dot + Forces[Fi]["f_a"]
            fi = f * u
            if Bi != 0:
                
                Bodies[Bi]["f"] = Bodies[Bi]["f"] - fi
                Bodies[Bi]["n"] = Bodies[Bi]["n"] - np.dot(Points[Pi]["sP_r"], fi)
                
            if Bj != 0:
                Bodies[Bj]["f"] = Bodies[Bj]["f"] + fi
                Bodies[Bj]["n"] = Bodies[Bj]["n"] + np.dot(Points[Pj]["sP_r"], fi)


        elif force_type == 'rot-sda': # Rotational spring-damper-actuator (FROM FILE sda_ROT.m)---------------------
                        
            Bi = Forces[Fi]["iBindex"]
            Bj = Forces[Fi]["jBindex"]
        
            if Bi == 0:
                theta = -Bodies[Bj]["p"]
                theta_d = -Bodies[Bj]["p_d"]
                T = Forces[Fi]["k"] * (theta - Forces[Fi]["theta0"]) + Forces[Fi]["dc"] * theta_d + Forces[Fi]["T_a"]
                Bodies[Bj]["n"] += T
            elif Bj == 0:
                theta = Bodies[Bi]["p"]
                theta_d = Bodies[Bi]["p_d"]
                T = Forces[Fi]["k"] * (theta - Forces[Fi]["theta0"]) + Forces[Fi]["dc"] * theta_d + Forces[Fi]["T_a"]
                Bodies[Bi]["n"] -= T
            else:
                theta = Bodies[Bi]["p"] - Bodies[Bj]["p"]
                theta_d = Bodies[Bi]["p_d"] - Bodies[Bj]["p_d"]
                T = Forces[Fi]["k"] * (theta - Forces[Fi]["theta0"]) + Forces[Fi]["dc"] * theta_d + Forces[Fi]["T_a"]
                Bodies[Bi]["n"] -= T
                Bodies[Bj]["n"] += T


        elif force_type == 'flocal': # ------------------------------------------------------------------------------
            Bi = Forces[Fi]["iBindex"]
            Bodies[Bi]["f"] += np.dot(Bodies[Bi]["A"], Forces[Fi]["flocal"])

        elif force_type == 'f': # -----------------------------------------------------------------------------------
            Bi = Forces[Fi]["iBindex"]
            Bodies[Bi]["f"] += Forces[Fi]["f"] 

        elif force_type == 'T': #------------------------------------------------------------------------------------
            Bi = Forces[Fi]["iBindex"]
            Bodies[Bi]["n"] += Forces[Fi]["T"]

        elif force_type == 'user': # ----------------------------- (FROM FILE user_force.m)

            def Friction_A(mu_s, mu_d, v_s, p, k_t, v, fy):
                # Anderson et al. friction model
                if abs(v) <= v_s:
                    ff = -mu_s * fy * (1 - (abs(v) / v_s) ** p)
                else:
                    ff = -mu_d * fy
                return ff
            
            mu_d = 0.15; mu_s = 0.2; mu_v = 0.0; v_s = 0.001; p = 2; k_t = 10000; fy = 9.81  # normal force
            v_conv = 0.1
            v = v_conv - Bodies[0]["r_d"][0]
            ff = Friction_A(mu_s, mu_d, v_s, p, k_t, v, fy)
            fx = ff + mu_v * v * fy
            fs = [fx, 0]
            Bodies[0]["f"] += fs



    g = np.zeros(nB3)
    for Bi in range(nB):
        ks = Bodies[Bi]["irc"]-1
        ke = ks + 3
        g[ks:ke] = np.concatenate( (Bodies[Bi]["f"], [Bodies[Bi]["n"]]) )

    return g


#%%    CORE Analysis function (passed to solver!)

def analysis(t, u): # function Analysis portata qui ma era su un file dedicato """
    # Solve the constrained equations of motion at time t with the standard Lagrange multiplier method
    
    
    #=== u to BODIES FROM "u_to_Bodies.m" finction ()-------------------------------------------------------------
    
    #    Unpack u into coordinate and velocity sub-arrays
    for Bi in range(nB): #              
        ir = Bodies[Bi]["irc"]-1
        ird = Bodies[Bi]["irv"]-1
        Bodies[Bi]["r"] = u[ir:ir+2]
        Bodies[Bi]["p"] = u[ir+2]
        Bodies[Bi]["r_d"] = u[ird:ird+2]
        Bodies[Bi]["p_d"] = u[ird+2]
    
    
    #=== update Positions FROM #def Update_Position():------------------------------------------------------------
    
    # Compute A's
    for Bi in range(nB): #                      
        #Bodies[Bi]["A"] = Matrix_A(Bodies[Bi]["p"])
        Bodies[Bi]["A"] =Functions.Matrix_A(Bodies[Bi]["p"])
        
    # Compute sP = A * sP_prime; rP = r + sP
    for Pi in range(nP):
        Bi = Points[Pi]["Bindex"]-1 #☺!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONFERMA SE IL -1 VA BENE!
        if Bi != 0:
            Points[Pi]["sP"] = Bodies[Bi]["A"] @ np.array((Points[Pi]["sPlocal"])) 
            Points[Pi]["sP_r"] = Functions.s_rot(Points[Pi]["sP"])
            Points[Pi]["rP"] = Bodies[Bi]["r"] + Points[Pi]["sP"]
            
    # Compute u = A * up
    for Vi in range(nU):
        Bi = Uvectors[Vi]["Bindex"]-1  #☺!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONFERMA SE IL -1 VA BENE!
        if Bi != 0:
            Uvectors[Vi]["u"] = Bodies[Bi]["A"] @ np.array((Uvectors[Vi]["ulocal"] ))
            Uvectors[Vi]["u_r"] = Functions.s_rot(Uvectors[Vi]["u"])

    
    
    # ===UPDATE VELOCITY, FROM def Update_Velocity():------------------------------------------------------------
    
    #   Compute sP_dot and rP_dot vectors
    for Pi in range(nP): #                   FROM Update_Velocity()
        Bi = Points[Pi]["Bindex"]-1  #☺!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONFERMA SE IL -1 VA BENE!
        if Bi != 0:
            Points[Pi]["sP_d"] = np.array((Points[Pi]["sP_r"] ))* np.array((Bodies[Bi]["p_d"]))
            Points[Pi]["rP_d"] = Bodies[Bi]["r_d"] + Points[Pi]["sP_d"]
            
    # Compute u_dot vectors
    for Vi in range(nU):
        Bi = Uvectors[Vi]["Bindex"]-1  #☺!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONFERMA SE IL -1 VA BENE!
        if Bi != 0:
            Uvectors[Vi]["u_d"] = Uvectors[Vi]["u_r"] * Bodies[Bi]["p_d"]

    
    
    h_a = Force_array(t)  # CALL, array of applied forces
    
    if nConst == 0:
        c_dd = M_inv_array * h_a  # solve for accelerations
    else:
        # Phi = Constraints(t)
        D = Jacobian.Jacobian()   # Jacobian matrix
        
        rhsA = Functions.RHSAcc(t, nJ, nConst, Bodies, Joints, Points, Uvectors) # r-h-s of acc constraints (gamma) #!!!!!!!!!! SEE FILE RHSAcc
        
        DMD = np.block([[np.diag(M_array), -D.T],
                        [D, np.zeros((nConst, nConst))]])
        rhs = np.concatenate((h_a, rhsA))
        sol = np.linalg.solve(DMD, rhs)
        c_dd = sol[:nB3]
        Lambda = sol[nB3:]
    
    for Bi in range(nB):
        ir = Bodies[Bi]["irc"]
        i2 = ir + 1
        i3 = i2 + 1
        Bodies[Bi]["r_dd"] = c_dd[ir:i2]
        Bodies[Bi]["p_dd"] = c_dd[i3]
    
    Functions.Bodies_to_u_d(nB6, nB, Bodies)  # pack velocities and accelerations into ud ========================== 
    
    
    num = num + 1  # number of function evaluations
    
    if showtime == 1:
        # Inform the user of the progress:
        # show the time once every 100 function evaluations
        if t10 % 100 == 0:
            print(t)
        t10 = t10 + 1


      
    #%%



PROVA = analysis(0, u) # ------------------------------------eliminare, mi serve solo per entrare nella FUNCT""" manualmente



#%% CALL THE SOLVER------
if t_final == 0: # one step function evaluation
    analysis(0, u)  # <<<<<<<<<<<<<<<<<<<<<<<-------------------------------------- !!! STILL WIP! SEE analysis.py
    T = np.array([0])
    uT = u.reshape(1, -1)
else:
    dt = 0.1 # float(input('Reporting time-step = ? ')) #<<<<<<<<<<<<< ACTUALLY FORCED TO 1 TO SIMPLIFY DEBUG
    Tspan = np.arange(t_initial, t_final+dt, dt)
    options = {'rtol': 1e-6, 'atol': 1e-9}
    T, uT = odeint(analysis, u, Tspan, options) # XXXXXXXXXXXXX FIX---------------------------------------------

num = len(T)
print("Number of function evaluations =", num)
print("\a")

#%% END
import winsound
frequency = 900; duration = 500; winsound.Beep(frequency, duration)




