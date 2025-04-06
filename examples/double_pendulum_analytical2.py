import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Pendulum rod lengths (m) and masses (kg).
L1, L2 = 1, 1.5
m1, m2 = 1, 2
# The gravitational acceleration (m.s-2).
g = 9.81
# Check that the integration conserves total energy to within this absolute
# tolerance.
EDRIFT = 0.05

def deriv(t, y, L1, L2, m1, m2):
    """Return the derivatives of y = theta1, p_theta1, theta2, p_theta2.

    These are the generalized coordinates (here, angles) and generalized
    momenta for the two rigid rods.

    """

    theta1, p_theta1, theta2, p_theta2 = y

    M = m1 + 3*m2
    c, s = np.cos(theta1 - theta2), np.sin(theta1 - theta2)
    Lr = L1 / L2
    den = 4 * M - 9 * m2 * c**2

    theta1dot = 6 / L1**2 * (2*p_theta1 - 3 * Lr * c * p_theta2) / den
    theta2dot = 6 / m2 / L2**2 * (
                    (2 * p_theta2 * M - 3 * m2 / Lr * c * p_theta1) / den)
    term = m2 * L1 * L2 / 2 * theta1dot * theta2dot * s
    p_theta1dot = -term - (m1/2 + m2) * g * L1 * np.sin(theta1)
    p_theta2dot = term - m2/2 * g * L2 * np.sin(theta2)

    return theta1dot, p_theta1dot, theta2dot, p_theta2dot

def calc_H(y, L1, L2, m1, m2):
    """Calculate the Hamiltonian at y = theta1, p_theta1, theta2, p_theta2."""
    theta1, p_theta1, theta2, p_theta2 = y

    theta1dot, p_theta1dot, theta2dot, p_theta2dot = deriv(None, y, L1, L2,
                                                           m1, m2)
    # The Lagrangian
    c = np.cos(theta1 - theta2)
    L = ( m1 * (L1 * theta1dot)**2 / 6 + m2 * (L2 * theta2dot)**2 / 6
             + m2 / 2 * ((L1 * theta1dot)**2 + L1*L2*theta1dot*theta2dot * c)
             + g * L1 * np.cos(theta1) * (m1 / 2 + m2)
             + g * L2 * np.cos(theta2) * m2 / 2
        )

    # The Hamiltonian
    H = p_theta1 * theta1dot + p_theta2 * theta2dot - L
    return H

# Maximum time, time point spacings and the time grid (all in s).
tmax, dt = 10, 0.01
t = np.arange(0, tmax+dt, dt)

# Initial conditions:
# angles: theta1, theta2 and generalized momenta: p_theta1, p_theta2
theta1_0, theta2_0 = np.pi, np.pi/2
# theta1_0, theta2_0 = 0, 0
p_theta1_0, p_theta2_0 = 0, 0
y0 = np.array([theta1_0, p_theta1_0, theta2_0, p_theta2_0])
# Could call calc_H, but since the initial p_thetai are zero, and the
# Hamiltonian is conserved (since the Langrangian has no explicit time-
# dependence, H0 is just the potential energy of the initial configuration:
H0 = -g * (L1 * np.cos(theta1_0) * (m1 / 2 + m2) +
           L2 * np.cos(theta2_0) * m2 / 2)

# Do the numerical integration of the equations of motion.
y = solve_ivp(deriv, (0, tmax), y0, method='LSODA', dense_output=True,
              args=(L1, L2, m1, m2))

# Check that the Hamiltonian didn't drift too much.
H = calc_H(y.y, L1, L2, m1, m2)
if any(abs(H-H0) > EDRIFT):
    print('Maximum energy drift exceeded')

# Unpack dynamical variables as a function of time.
theta1, p_theta1, theta2, p_theta2 = y.sol(t)

# Convert to Cartesian coordinates of the two rods.
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# # Plot a trail of the m2 bob's position for the last trail_secs seconds.
# trail_secs = 1
# # This corresponds to max_trail time points.
# max_trail = int(trail_secs / dt)

# def make_plot(i):
#     """
#     Plot and save an image of the double pendulum configuration for time
#     point i.
#     """

#     # The pendulum rods (thick, black).
#     ax.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], lw=8, c='k')

#     # The trail will be divided into ns segments and plotted as a fading line.
#     ns = 20
#     s = max_trail // ns

#     for j in range(ns):
#         imin = i - (ns-j)*s
#         if imin < 0:
#             continue
#         imax = imin + s + 1
#         # The fading looks better if we square the fractional length along the
#         # trail.
#         alpha = (j/ns)**2
#         # Two trails, initiating at the centres of mass of the two rods.
#         ax.plot(x1[imin:imax]/2, y1[imin:imax]/2, c='b', solid_capstyle='butt',
#                 lw=2, alpha=alpha)
#         ax.plot((x1[imin:imax]+x2[imin:imax])/2,
#                 (y1[imin:imax]+y2[imin:imax])/2,
#                 c='r', solid_capstyle='butt', lw=2, alpha=alpha)

#     # Centre the image on the fixed anchor point, and ensure the axes are equal
#     ax.set_xlim(-L1-L2, L1+L2)
#     ax.set_ylim(-L1-L2, L1+L2)
#     ax.set_aspect('equal', adjustable='box')
#     plt.axis('off')
#     plt.savefig('frames/_img{:04d}.png'.format(i//di), dpi=72)
#     plt.cla()


# # Make an image every di time points, corresponding to a frame rate of fps
# # frames per second.
# # Frame rate, s-1
# fps = 10
# di = int(1/fps/dt)
# fig = plt.figure(figsize=(8.3333, 6.25), dpi=72)
# ax = fig.add_subplot(111)

# for i in range(0, t.size, di):
#     print(i // di, '/', t.size // di)
#     make_plot(i)

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot x1 and y1 in the first subplot
ax1.plot(t, x1, label='x1')
ax1.plot(t, y1, label='y1')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Coordinates (m)')
ax1.legend()
ax1.set_title('Coordinates of the First Pendulum')

# Plot x2 and y2 in the second subplot
ax2.plot(t, x2, label='x2')
ax2.plot(t, y2, label='y2')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Coordinates (m)')
ax2.legend()
ax2.set_title('Coordinates of the Second Pendulum')

plt.tight_layout()
plt.show()