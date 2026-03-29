import numpy as np
import scipy as sc
import matplotlib as mpl
from PMD.src import *
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from PMD.examples._plot_utils import plot_comparison


#%% bodies
b1 = Body(mass=1, inertia=0.0833, position=[0.0, 0.5], orientation=0.0)
b2 = Body(mass=2, inertia=0.375, position=[0.75, 1.0], orientation=0.0)

#%% markers
p0 = Ground.add_marker([0.0, 0.0])    # revolute joint between b1 and ground
p1 = b1.add_marker([0.0, -0.5])       # revolute joint between b1 and ground
p2 = b1.add_marker([0.0, 0.5])        # revolute joint between b1 and b2
p3 = b2.add_marker([-0.75, 0])        # revolute joint between b1 and b2

#%% joints
j1 = RevJoint(iMarker=p0, jMarker=p1)  # Revolute joint between b1 and ground
j2 = RevJoint(iMarker=p2, jMarker=p3)  # Revolute joint between b1 and b2

#%% forces
s3 = Weight()  # only weight force, acting along -y axis

#%% double pendulum model creation
double_pendulum = PlanarMultibodyModel(
    bodies=[b1, b2],
    joints=[j1, j2],
    forces=[s3])
T, uT = double_pendulum.solve(method='Radau', t_final=10.0, t_eval=np.linspace(0, 10, 10001))

if __name__ == '__main__':
    plot_comparison(T, uT, matlab_filename='DP.txt', model_title='DP')

# =============================================================================
# CALCULATION AND PLOT OF CONSTRAINT VIOLATIONS (PAPER STYLE)
# =============================================================================

# 0. Paper style settings (clear dark style)
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Comic Sans MS']
plt.rcParams['mathtext.fontset'] = 'cm' # Use TeX-style math font

# 1. Array initialization
n_steps = len(T)
phi_total_norm = np.zeros(n_steps)
phi_j1_norm = np.zeros(n_steps)
phi_j2_norm = np.zeros(n_steps)

# 2. Error calculation
for i in range(n_steps):
    double_pendulum._u2bodies(uT[i])
    double_pendulum._update_position()
    phi = double_pendulum._compute_constraints()
    
    phi_total_norm[i] = np.linalg.norm(phi)
    phi_j1_norm[i] = np.linalg.norm(phi[j1._rows : j1._rowe])
    phi_j2_norm[i] = np.linalg.norm(phi[j2._rows : j2._rowe])

# 3. Figure creation
fig, ax = plt.subplots(figsize=(10, 4.5)) # Aspect ratio matching your plots

# 4. Plotting lines
# - Total error: thick, light line in the background
ax.plot(T, phi_total_norm, label='Total Error ($||\mathbf{\Phi}||$)', 
        color='tab:blue', alpha=0.4, linewidth=5)

# - Joint 1: Green (matches body_1 color in the schematic)
ax.plot(T, phi_j1_norm, label='Revolute Joint 1 ($P_{10}-P_{11}$)', linewidth=1.8)

# - Joint 2: Blue (matches body_2 color in the schematic)
ax.plot(T, phi_j2_norm, label='Revolute Joint 2 ($P_{21}-P_{22}$)', linewidth=1.8)

# 5. Axes formatting
ax.set_yscale('log')

ax.set_xlabel('Time (s)', fontsize=14)
ax.set_ylabel('Constraint Violation (m)', fontsize=14)

ax.set_title(
    'Constraint Violation: PMD library model',
    fontsize=14,
)

# Grid (including minor ticks for logarithmic scale)
ax.grid(True, which='both', linestyle=':', color='gray', alpha=0.5)

# Tick label size (molto importante con Comic Sans)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=10)

# 6. Legend
ax.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.15),
    ncol=3,
    framealpha=1,
    edgecolor='lightgray',
    fontsize=12
)

# Adjust margins to make room for the bottom legend
plt.tight_layout()
plt.subplots_adjust(bottom=0.22) 

if __name__ == '__main__':
    # Uncomment if you also want to run your original comparison plot:
    # plot_comparison(T, uT, matlab_filename='DP.txt', model_title='DP')
    
    # Save the vector figure for LaTeX (highly recommended)
    plt.savefig('constraint_violation.pdf', format='pdf', bbox_inches='tight')
    plt.show()
