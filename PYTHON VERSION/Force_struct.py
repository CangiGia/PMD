# Force_struct
"""
Created on Thu Sep 28 15:25:37 2023

@author: simone.lucertini
"""




# creating force DICT

Force = {
    "type"   : "ptp",  # element type: ptp, rot_sda, weight, fp, f, T
    "iPindex": 0    ,  # index of the head (arrow) point
    "jPindex": 0    ,  #index of the tail point
    "iBindex": 0    ,  # index of the head (arrow) body
    "jBindex": 0    ,  # index of the tail body
    "k"      : 0    ,  # spring stiffness
    "L0"     : 0    ,  # undeformed length
    "theta0" : 0    ,  #undeformed angle
    "dc"     : 0    ,  # damping coefficient
    "f_a"    : 0    ,  # constant actuator force
    "T_a"    : 0    ,  #constant actuator torque
    "gravity": 9.81 ,  # gravitational constant
    "wgt"    : [0,-1], # gravitational direction
    "flocal" : [0,0],  # constant force in local frame
    "f"      : [0,0],  # constant force in x-y frame
    "T"      : 0    ,  # constant torque in x-y frame
    "iFunct" : 0       # analytical function index  
}

