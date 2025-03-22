import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt

current_directory = os.path.dirname(os.path.abspath(__file__))

uT_python = np.loadtxt(os.path.join(current_directory, "uT_python.txt"))
T_python = np.loadtxt(os.path.join(current_directory, "T_python.txt"))

uT_python_LSODA = np.loadtxt(os.path.join(current_directory, "uT_python_LSODA.txt"))
T_python_LSODA = np.loadtxt(os.path.join(current_directory, "T_python_LSODA.txt"))

uT_python_RK45 = np.loadtxt(os.path.join(current_directory, "uT_python_RK45.txt"))
T_python_RK45 = np.loadtxt(os.path.join(current_directory, "T_python_RK45.txt"))

uT_matlab = np.loadtxt(os.path.join(current_directory, "uT_matlab.txt"))
T_matlab = np.loadtxt(os.path.join(current_directory, "T_matlab.txt"))

# -----------------------------
# ho 12 colonne nei file sopra
# le prime 6 sono le posizioni
# le seconde 6 sono le velocità 
# X1 Y1 Phi1 X2 Y2 Phi2 || dX1 dY1 dPhi1 dX2 dY2 dPhi2
# -----------------------------

plt.figure(figsize = (15, 4.5))
# plt.plot(T_python, uT_python[:,5], label = "python")
plt.plot(T_matlab, uT_matlab[:,3], label = "matlab")
# plt.plot(T_python_LSODA, uT_python_LSODA[:,5], label = "python")
plt.plot(T_python_RK45, uT_python_RK45[:,3], label = "python")
plt.legend()
plt.tight_layout()
plt.show()