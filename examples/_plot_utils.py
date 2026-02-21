"""
Shared plotting utility for PMD example scripts.
Generates comparison plots between Python and MATLAB simulation results.
"""
import os
import numpy as np
import matplotlib.pyplot as plt


# Path to MATLAB reference results
MATLAB_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    '_references',
    'Nikravesh - Planar Multibody Dynamics 2nd Edition',
    'DAP_BC_12_2018',
    '_results'
)


def plot_comparison(T: np.ndarray, uT: np.ndarray, matlab_filename: str, model_title: str) -> None:
    """
    Plot Python vs MATLAB displacements for each body.

    Creates a figure with nB rows and 3 columns:
      - Column 1: x displacement for each body
      - Column 2: y displacement for each body
      - Column 3: phi (rotation) for each body

    Each subplot shows Python result (solid line) vs MATLAB result (dashed line).

    Args:
        T:               Time vector from Python simulation, shape (nSteps,)
        uT:              State matrix from Python simulation, shape (nSteps, 2*nB3)
        matlab_filename: Name of the MATLAB .txt file (e.g., 'AA.txt')
        model_title:     Title for the figure (e.g., 'AA')
    """
    # --- Derive number of bodies ---
    nB3 = uT.shape[1] // 2          # nB3 = 3 * nB
    nB  = nB3 // 3                   # number of moving bodies

    # --- Displacements only (first half of uT) ---
    disp = uT[:, :nB3]               # shape (nSteps, nB3)

    # --- Load MATLAB reference ---
    matlab_path = os.path.join(MATLAB_RESULTS_DIR, matlab_filename)
    if not os.path.exists(matlab_path):
        print(f"[plot_comparison] WARNING: MATLAB file not found: {matlab_path}")
        mat_data = None
    else:
        mat_data = np.loadtxt(matlab_path, delimiter=',')
        # mat_data columns: [t, x1, y1, p1, ..., xN, yN, pN, dx1, dy1, dp1, ...]
        T_mat    = mat_data[:, 0]
        disp_mat = mat_data[:, 1:nB3 + 1]   # displacement columns only

    # --- Build figure ---
    dof_labels = ['x [m]', 'y [m]', r'$\varphi$ [rad]']
    dof_keys   = ['x', 'y', 'φ']

    fig, axes = plt.subplots(
        nrows=nB, ncols=3,
        figsize=(14, 3.5 * nB),
        squeeze=False
    )
    fig.suptitle(f'PMD — {model_title}: Python vs MATLAB', fontsize=14, fontweight='bold')

    for bi in range(nB):              # body index (0-based)
        body_label = f'Body {bi + 1}'

        for dof in range(3):          # 0=x, 1=y, 2=phi
            ax      = axes[bi][dof]
            col_idx = 3 * bi + dof    # column in disp array

            # Python result
            ax.plot(T, disp[:, col_idx],
                    color='tab:blue', linewidth=1.5,
                    label='Python')

            # MATLAB result (if available)
            if mat_data is not None:
                ax.plot(T_mat, disp_mat[:, col_idx],
                        color='tab:orange', linewidth=1.5,
                        linestyle='--', label='MATLAB')

            # Labels
            ax.set_xlabel('t [s]', fontsize=9)
            ax.set_ylabel(dof_labels[dof], fontsize=9)
            ax.set_title(f'{body_label} — {dof_keys[dof]}', fontsize=10)
            ax.legend(fontsize=8)
            ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.show()
