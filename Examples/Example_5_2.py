import numpy as np
from scipy.optimize import fsolve
from PMD_functions import * 


#* fourbar function definition
def fourbar(x, theta1):
    u1 = np.array([np.cos(theta1), np.sin(theta1)])
    u2 = np.array([np.cos(x[0]), np.sin(x[0])])
    u3 = np.array([np.cos(x[1]), np.sin(x[1])])
    phi = L1 * u1 + L2 * u2 - L3 * u3 - np.array([L0, 0])
    return phi

#* data
L0, L1, L2, L3 = 0.2, 0.1, 0.3, 0.22 
theta1 = np.pi / 4
theta23_initial = np.array([0.5, 1.3])

#* solution
# positions
angles = fsolve(fourbar, theta23_initial, args=(theta1,), xtol=1e-9)
angles_deg = np.degrees(angles)
theta2 = angles[0]
theta3 = angles[1]
print(f"Theta2 (in radians): {theta2}")
print(f"Theta3 (in radians): {theta3}")
# print(f"Theta2 (in degrees): {angles_deg[0]}")
# print(f"Theta3 (in degrees): {angles_deg[1]}")

# velocities
theta1_d = 1.5
u1 = colvect(np.cos(theta1), np.sin(theta1))
u2 = colvect(np.cos(theta2), np.sin(theta2))
u3 = colvect(np.cos(theta3), np.sin(theta3))
C = np.column_stack((L2*u2, -L3*u3)) # jacobian
rhsv = -L1*u1*theta1_d # right-hand-side velocity
solution_v = np.linalg.solve(C,rhsv)
theta2_d = solution_v[0]
theta3_d = solution_v[1]
print(f"Theta2 (in radians/sec): {theta2_d}")
print(f"Theta3 (in radians/sec): {theta3_d}")

# accelerations 
theta1_dd = -1.2
rhsa = (-L1 * s_rot(u1) * theta1_dd +
        L1 * u1 * (theta1_d ** 2) +
        L2 * u2 * (theta2_d ** 2) -
        L3 * u3 * (theta3_d ** 2))
solution_a = np.linalg.solve(C,rhsa)
theta2_dd = solution_a[0]
theta3_dd = solution_a[1]
print(f"Theta2 (in radians/sec^2): {theta2_dd}")
print(f"Theta3 (in radians/sec^2): {theta3_dd}")

Lpa = 0.2
beta2 = np.pi / 6
upa = np.array([[np.cos(theta2 + beta2)], [np.sin(theta2 + beta2)]])
rP = L1 * u1 + Lpa * upa
rP_d = L1 * theta1_d * s_rot(u1) + Lpa * theta2_d * s_rot(upa)
rP_dd = (L1 * (theta1_dd * s_rot(u1) - theta1_d ** 2 * u1) + 
          Lpa * (theta2_dd * s_rot(upa) - theta2_d ** 2 * upa))

ecchime = 1