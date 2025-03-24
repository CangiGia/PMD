import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt

current_directory = os.path.dirname(os.path.abspath(__file__))

uT_python = np.loadtxt(os.path.join(current_directory, "uT_python.txt"))
T_python = np.loadtxt(os.path.join(current_directory, "T_python.txt"))

uT_matlab = np.loadtxt(os.path.join(current_directory, "uT_matlab.txt"))
T_matlab = np.loadtxt(os.path.join(current_directory, "T_matlab.txt"))

# -----------------------------
# ho 12 colonne nei file sopra
# le prime 6 sono le posizioni
# le seconde 6 sono le velocità 
# X1 Y1 Phi1 X2 Y2 Phi2 || dX1 dY1 dPhi1 dX2 dY2 dPhi2
# -----------------------------

fig, ax = plt.subplots(figsize=(15, 4.5))
ax.plot(T_python, uT_python[:, 3], label="python")
ax.plot(T_matlab, uT_matlab[:, 3], label="matlab")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Dispacement")
ax.set_title("Comparison between Python and Matlab")
ax.legend()
plt.show()