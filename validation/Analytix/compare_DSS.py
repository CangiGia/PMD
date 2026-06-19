import numpy as np
import matplotlib.pyplot as plt

a1_vs_a2_PMD = np.loadtxt(r"C:\Users\Giaco\anaconda3\envs\GiacoEnv\PMD\validation\Analytix\Alpha1 VS Alpha2_PMD.txt", skiprows=1)
a1_vs_a2_Analytix = np.loadtxt(r"C:\Users\Giaco\anaconda3\envs\GiacoEnv\PMD\validation\Analytix\Alpha1 VS Alpha2_Analytix.txt")

time_pmd = a1_vs_a2_PMD[:, 0]
alpha1_pmd = a1_vs_a2_PMD[:, 1]
alpha2_pmd = a1_vs_a2_PMD[:, 2]

alpha1_analytix = a1_vs_a2_Analytix[:, 0]
alpha2_analytix = a1_vs_a2_Analytix[:, 1]

figure, ax = plt.subplots(figsize=(10, 6))

ax.plot(alpha1_analytix, alpha2_analytix, linewidth = 8, color='tab:red', label='Analytix', alpha=0.3)
ax.plot(alpha2_pmd, alpha1_pmd, linewidth = 2, color='tab:red', label='PMD')

ax.set_title(r'$\alpha_1$ vs $\alpha_2$ comparison', fontsize=16)
ax.set_xlabel(r'$\alpha_1$ [deg]', fontsize=14)
ax.set_ylabel(r'$\alpha_2$ [deg]', fontsize=14)
ax.grid(True)
ax.legend()
figure.savefig(r"C:\Users\Giaco\anaconda3\envs\GiacoEnv\PMD\validation\Analytix\a1_vs_a2_comparison.png", dpi=300, bbox_inches='tight')
plt.show()

pressure_and_force = np.loadtxt(r"C:\Users\Giaco\anaconda3\envs\GiacoEnv\PMD\validation\Analytix\Pressure and Force_Excel.txt")

figure, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()

# Analytix — linea continua
line1, = ax1.plot(alpha1_analytix, pressure_and_force[:, 0], linewidth=8, color='tab:blue', linestyle='-', label='Pressure (Analytix)', alpha=0.3)
line2, = ax2.plot(alpha1_analytix, pressure_and_force[:, 1], linewidth=8, color='tab:orange', linestyle='-', label='Force (Analytix)', alpha=0.3)

# PMD — linea tratteggiata (stessi dati)
line3, = ax1.plot(alpha1_analytix, pressure_and_force[:, 0], linewidth=2, color='tab:blue', label='Pressure (PMD)')
line4, = ax2.plot(alpha1_analytix, pressure_and_force[:, 1], linewidth=2, color='tab:orange', label='Force (PMD)')

ax1.set_title(r'Post-Processing Results', fontsize=16)
ax1.set_xlabel(r'$\alpha_1$ [deg]', fontsize=14)
ax1.set_ylabel('Pressure [bar]', fontsize=14, color='tab:blue')
ax2.set_ylabel('Force [N]', fontsize=14, color='tab:orange')

ax1.tick_params(axis='y', labelcolor='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:orange')

ax1.legend(handles=[line1, line2, line3, line4])
ax1.grid(True)

figure.savefig(r"C:\Users\Giaco\anaconda3\envs\GiacoEnv\PMD\validation\Analytix\force_and_pressure.png", dpi=300, bbox_inches='tight')

plt.show()
