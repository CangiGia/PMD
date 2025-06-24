import os
import numpy as np
import matplotlib.pyplot as plt

# current_directory = os.path.dirname(os.path.abspath(__file__))

# time_ADAMS = np.loadtxt(os.path.join(current_directory, "time_ADAMS.txt"), skiprows=6)
# x1_ADAMS = np.loadtxt(os.path.join(current_directory, "x1_ADAMS.txt"), skiprows=6)
# x2_ADAMS = np.loadtxt(os.path.join(current_directory, "x2_ADAMS.txt"), skiprows=6)

# uT_python = np.loadtxt(os.path.join(current_directory, "uT_python.txt"))
# T_python = np.loadtxt(os.path.join(current_directory, "T_python.txt"))

# uT_matlab = np.loadtxt(os.path.join(current_directory, "uT_matlab.txt"))
# T_matlab = np.loadtxt(os.path.join(current_directory, "T_matlab.txt"))

# -----------------------------
# ho 12 colonne nei file sopra
# le prime 6 sono le posizioni
# le seconde 6 sono le velocità 
# X1 Y1 Phi1 X2 Y2 Phi2 || dX1 dY1 dPhi1 dX2 dY2 dPhi2
# -----------------------------

# fig, ax = plt.subplots(figsize=(15, 4.5))
# ax.plot(time_ADAMS, x1_ADAMS, label="adams")
# ax.plot(T_matlab, uT_matlab[:, 0], label="matlab")
# ax.plot(T_python, uT_python[:, 0], label="python")
# ax.set_xlabel("Time [s]")
# ax.set_ylabel("Dispacement")
# ax.set_title("Comparison between Adams, Python and Matlab")
# ax.legend()
# plt.show()

lam = np.loadtxt("AA_lambda_matlab.txt")
time = np.loadtxt("AA_time_matlab.txt").T

fig, ax = plt.subplots(figsize=(15, 4.5))
ax.plot(time, lam[:,0], label="Fx")
ax.plot(time, lam[:,1], label="Fy")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Force [N]")
ax.set_title("Reactioin force at Q revolute joint")
ax.legend()
plt.show()