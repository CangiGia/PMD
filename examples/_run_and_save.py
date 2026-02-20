"""
Runner script: executes a test model and saves T, uT as .npy files.
Usage: python -m PMD.examples._run_and_save <model_name>
Example: python -m PMD.examples._run_and_save _test_SP
Must be run from C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
"""
import sys
import os
import importlib
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend to avoid plt.show() blocking

model_name = sys.argv[1]  # e.g., "_test_SP"
out_dir_name = sys.argv[2] if len(sys.argv) > 2 else 'npy_ref'

# Import the module — this runs the model at module level
mod = importlib.import_module(f"PMD.examples.{model_name}")

T = mod.T
uT = mod.uT

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', out_dir_name)
os.makedirs(out_dir, exist_ok=True)

np.save(os.path.join(out_dir, f'{model_name}_T.npy'), T)
np.save(os.path.join(out_dir, f'{model_name}_uT.npy'), uT)
print(f"Saved {model_name}_T.npy and {model_name}_uT.npy to {out_dir_name}/")
