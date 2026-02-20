"""
Shared fixtures for PMD test suite.

All tests must be run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command: python -m pytest PMD/tests/ -v
"""
import pytest
import os
import sys
import subprocess
import numpy as np

# Path constants
PMD_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES_DIR = os.path.join(PMD_ROOT, 'examples')
NPY_REF_DIR = os.path.join(EXAMPLES_DIR, 'results', 'npy_ref')
WORKING_DIR = os.path.dirname(PMD_ROOT)  # C:\Users\Giaco\anaconda3\envs\GiacoEnv\

# List of all models to test (exclude _test_Rod which has known issues)
ALL_MODELS = [
    '_test_AA',
    '_test_Cart_A',
    '_test_Cart_B',
    '_test_Cart_C',
    '_test_Cart_D',
    '_test_CB',
    '_test_DP',
    '_test_MP_A',
    '_test_MP_B',
    '_test_MP_C',
    '_test_SP',
]


def run_model_subprocess(model_name: str) -> tuple:
    """
    Run a PMD model in an isolated subprocess and return (T, uT).

    Each model is run in a separate process to avoid global state
    contamination (Base.COUNT and subclass counters).
    Matplotlib is forced to 'Agg' backend to prevent plt.show() blocking.

    Args:
        model_name: Model name without .py extension (e.g., '_test_SP')

    Returns:
        Tuple of (T, uT) numpy arrays

    Raises:
        RuntimeError: If the subprocess fails
    """
    # Temp files for output
    tmp_dir = os.path.join(PMD_ROOT, 'tests', '_tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    t_path = os.path.join(tmp_dir, f'{model_name}_T.npy')
    ut_path = os.path.join(tmp_dir, f'{model_name}_uT.npy')

    # Script to execute in subprocess
    script = f"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import importlib

mod = importlib.import_module('PMD.examples.{model_name}')
np.save(r'{t_path}', mod.T)
np.save(r'{ut_path}', mod.uT)
"""

    result = subprocess.run(
        [sys.executable, '-c', script],
        cwd=WORKING_DIR,
        capture_output=True,
        text=True,
        timeout=300  # 5 minutes max per model
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Model {model_name} failed:\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )

    T = np.load(t_path)
    uT = np.load(ut_path)

    # Cleanup temp files
    os.remove(t_path)
    os.remove(ut_path)

    return T, uT


def load_reference(model_name: str) -> tuple:
    """
    Load golden reference data for a model.

    Args:
        model_name: Model name without .py extension (e.g., '_test_SP')

    Returns:
        Tuple of (T_ref, uT_ref) numpy arrays

    Raises:
        FileNotFoundError: If reference files don't exist
    """
    t_path = os.path.join(NPY_REF_DIR, f'{model_name}_T.npy')
    ut_path = os.path.join(NPY_REF_DIR, f'{model_name}_uT.npy')

    if not os.path.exists(t_path) or not os.path.exists(ut_path):
        raise FileNotFoundError(
            f"Reference files not found for {model_name}. "
            f"Expected:\n  {t_path}\n  {ut_path}"
        )

    return np.load(t_path), np.load(ut_path)
