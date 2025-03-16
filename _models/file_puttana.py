import numpy as np
import matplotlib.pyplot as plt

uT_python = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\uT_python.txt")
T_python = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\T_python.txt")

uT_python_LSODA = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\uT_python_LSODA.txt")
T_python_LSODA = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\T_python_LSODA.txt")

uT_python_RK45 = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\uT_python_RK45.txt")
T_python_RK45 = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\T_python_RK45.txt")

uT_matlab = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\uT_matlab.txt")
T_matlab = np.loadtxt(r"C:\Users\Giacomo\anaconda3\envs\GiacoEnv\PMD\\T_matlab.txt")

# -----------------------------
# ho 12 colonne nei file sopra
# le prime 6 sono le posizioni
# le seconde 6 sono le velocità 
# X1 Y1 Phi1 X2 Y2 Phi2 || dX1 dY1 dPhi1 dX2 dY2 dPhi2
# -----------------------------

plt.figure(figsize = (15, 4.5))
#plt.plot(T_python, uT_python[:,5], label = "python")
plt.plot(T_matlab, uT_matlab[:,4], label = "matlab")
#plt.plot(T_python_LSODA, uT_python_LSODA[:,5], label = "python")
plt.plot(T_python_RK45, uT_python_RK45[:,4], label = "python")
plt.legend()
plt.tight_layout()
plt.show()