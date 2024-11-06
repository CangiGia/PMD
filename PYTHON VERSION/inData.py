# inData
"""
Created on Fri Sep 29 09:52:07 2023

@author: simone.lucertini
"""

# main.py


with open("inBodies.py") as f: #     inAnimate è stato messo qui dentro, verificare vada bene!!!!
    exec(f.read())
with open("inPoints.py") as f:
    exec(f.read())
with open("inUvectors.py") as f:
    exec(f.read())    
with open("inForces.py") as f:
    exec(f.read())
with open("inJoints.py") as f:
    exec(f.read())
with open("inFuncts.py") as f:
    exec(f.read())
# with open("inAnimate.py") as f: # !!! l'ho spostata in inBodies
#     exec(f.read())    

def init(): # Vars to be shared with parent scripts
    global Bodies, Points, Vectors, Forces, Joints, Functs, Points_anim, Uvectors