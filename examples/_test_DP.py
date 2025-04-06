import numpy as np
import scipy as sc
import matplotlib as mpl
from PMD.src.builder import *
from PMD.src.functions import *
from PMD.src.solver import *
import matplotlib.pyplot as plt


#%% bodies
b1 = Body(m=1, J=0.0833, r=[0.0,0.5], p=0.0)
b2 = Body(m=2, J=0.375, r=[0.75,1.0], p=0.0)

#%% points
p0 = Point(Bindex=0, sPlocal=np.array([0.0, 0.0]))     #// revolute joint between b1 and b0 - ground side
p1 = Point(Bindex=1, sPlocal=np.array([0.0, -0.5]))    #// revolute joint between b1 and b0 - body 1 side
p2 = Point(Bindex=1, sPlocal=np.array([0.0, 0.5]))     #// revolute joint between b1 and b2 - body 1 side
p3 = Point(Bindex=2, sPlocal=np.array([-0.75, 0]))     #// revolute joint between b1 and b2 - body 2 side

#%% joints
j1 = Joint(type="rev", iPindex=0, jPindex=1) #// Revolute joint between b1 and b0
j2 = Joint(type="rev", iPindex=2, jPindex=3) #// Revolute joint between b1 and b2

#%% forces
s3 = Force()
s3 = Force(type="weight") #// only weight force, acting along -y axis

#%% double pendulum model creation
double_pendulum = PlanarDynamicModel()
time, solution = double_pendulum.solve(method='Radau')

# fig, ax = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

# ax[0].plot(time, solution[:, 0], label=r'$x_1$', color='b', linestyle='-', linewidth=2)
# ax[0].plot(time, solution[:, 1], label=r'$y_1$', color='r', linestyle='-', linewidth=2)
# ax[0].set_ylabel('Displacement (m)', fontsize=12)
# ax[0].legend()

# ax[1].plot(time, solution[:, 3], label=r'$x_2$', color='b', linestyle='-', linewidth=2)
# ax[1].plot(time, solution[:, 4], label=r'$y_2$', color='r', linestyle='-', linewidth=2)
# ax[1].set_ylabel('Displacement (m)', fontsize=12)
# ax[1].legend()

# # ax[2].plot(time, solution[:, 2], label=r'$\psi_1$', color='b', linestyle='-', linewidth=2)
# # ax[2].plot(time, solution[:, 5], label=r'$\psi_2$', color='r', linestyle='-', linewidth=2)
# # ax[2].set_xlabel('Time (s)', fontsize=12)
# # ax[2].set_ylabel('Displacement (rad)', fontsize=12)
# # ax[2].legend()

# fig.suptitle('System responses', fontsize=14)
# plt.show()

#%% analytical solution
# Pendulum rod lengths (m) and masses (kg).
L1, L2 = 1, 1.5
m1, m2 = 1, 2
g = 9.81
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

tmax, dt = 6, 0.01
t = np.arange(0, tmax+dt, dt)

# Initial conditions:
# angles: theta1, theta2 and generalized momenta: p_theta1, p_theta2
theta1_0, theta2_0 = np.pi, np.pi/2
p_theta1_0, p_theta2_0 = 0, 0
y0 = np.array([theta1_0, p_theta1_0, theta2_0, p_theta2_0])
H0 = -g * (L1 * np.cos(theta1_0) * (m1 / 2 + m2) +
           L2 * np.cos(theta2_0) * m2 / 2)

y = solve_ivp(deriv, (0, tmax), y0, method='Radau', dense_output=True,
              args=(L1, L2, m1, m2))

H = calc_H(y.y, L1, L2, m1, m2)
if any(abs(H-H0) > EDRIFT):
    print('Maximum energy drift exceeded')

theta1, p_theta1, theta2, p_theta2 = y.sol(t)
x1 = L1/2 * np.sin(theta1)
y1 = -L1/2 * np.cos(theta1)
x2 = (L1 * np.sin(theta1)) + L2/2 * np.sin(theta2)
y2 = (-L1 * np.cos(theta1)) - L2/2 * np.cos(theta2)

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
# fps = 12
# di = int(1/fps/dt)
# fig = plt.figure(figsize=(8.3333, 6.25), dpi=72)
# ax = fig.add_subplot(111)

# for i in range(0, t.size, di):
#     print(i // di, '/', t.size // di)
#     make_plot(i)


# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

# ax1.plot(t, x1, label='x1')
# ax1.plot(t, y1, label='y1')
# ax1.set_xlabel('Time (s)')
# ax1.set_ylabel('Coordinates (m)')
# ax1.legend()
# ax1.set_title('Coordinates of the First Pendulum')

# ax2.plot(t, x2, label='x2')
# ax2.plot(t, y2, label='y2')
# ax2.set_xlabel('Time (s)')
# ax2.set_ylabel('Coordinates (m)')
# ax2.legend()
# ax2.set_title('Coordinates of the Second Pendulum')

# plt.tight_layout()
# plt.show()

mpl.rcParams["font.family"] = "Comic Sans MS"
fig, axs = plt.subplots(2, 2, figsize=(15, 4.5), sharex=True)

axs[0, 0].plot(t, x1, color='greenyellow', linestyle='-', linewidth=7)
axs[0, 0].plot(time, solution[:, 0], color='k', linestyle='--', linewidth=2)
axs[0, 0].set_ylabel('Displacement (m)', fontsize=12)
axs[0, 0].set_title(r'CoG $x_1$ coordinate', fontsize=12)

axs[0, 1].plot(t, y1, color='greenyellow', linestyle='-', linewidth=7)
axs[0, 1].plot(time, solution[:, 1], color='k', linestyle='--', linewidth=2)
axs[0, 1].set_ylabel('Displacement (m)', fontsize=12)
axs[0, 1].set_title(r'CoG $y_1$ coordinate', fontsize=12)

axs[1, 0].plot(t, x2, color='greenyellow', linestyle='-', linewidth=7)
axs[1, 0].plot(time, solution[:, 3], color='k', linestyle='--', linewidth=2)
axs[1, 0].set_xlabel('Time (s)', fontsize=12)
axs[1, 0].set_ylabel('Displacement (m)', fontsize=12)
axs[1, 0].set_title(r'CoG $x_2$ coordinate', fontsize=12)

axs[1, 1].plot(t, y2, label='analytical model', color='greenyellow', linestyle='-', linewidth=7)
axs[1, 1].plot(time, solution[:, 4], label='pmd library', color='k', linestyle='--', linewidth=2)
axs[1, 1].set_xlabel('Time (s)', fontsize=12)
axs[1, 1].set_ylabel('Displacement (m)', fontsize=12)
axs[1, 1].set_title(r'CoG $y_2$ coordinate', fontsize=12)

fig.legend(loc='outside lower center', fontsize=12, ncol=2)
fig.suptitle('Responses comparison: analytical model VS pmd library', fontsize=14)
plt.tight_layout()
fig.savefig('pmdVSanalytical.png', dpi=900)
plt.show()
