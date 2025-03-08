import numpy as np
import matplotlib.pyplot as plt

uT_python = np.loadtxt("_models/uT_python.txt")
T_python = np.loadtxt("_models/T_python.txt")

uT_matlab = np.loadtxt(r"C:\Users\Giacomo\OneDrive\Documents\03_Ing. Mec. PhD\03_Third Year\THESIS\Nikravesh - Planar Multibody Dynamics 2nd Edition\DAP-BC_11-2018\DAP_BC 11_2018\uT_matlab.txt")
T_matlab = np.loadtxt(r"C:\Users\Giacomo\OneDrive\Documents\03_Ing. Mec. PhD\03_Third Year\THESIS\Nikravesh - Planar Multibody Dynamics 2nd Edition\DAP-BC_11-2018\DAP_BC 11_2018\T_matlab.txt")

# ho 12 colonne nei file sopra
# le prime 6 sono le posizioni
# le seconde 6 sono le velocità 
# X1 Y1 Phi1 X2 Y2 Phi2 || dX1 dY1 dPhi1 dX2 dY2 dPhi2
plt.figure(figsize = (15, 4.5))
plt.plot(T_python, uT_python[:,11], label = "pitone")
plt.plot(T_matlab, uT_matlab[:,11], label = "matolabo")
plt.legend()
plt.show()
