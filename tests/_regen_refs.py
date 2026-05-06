"""Regenerate golden .npy reference files for all models."""
import os
import sys
import subprocess
import shutil

models = [
    '_test_AA', '_test_Cart_A', '_test_Cart_B', '_test_Cart_C', '_test_Cart_D',
    '_test_CB', '_test_DP', '_test_MP_A', '_test_MP_B', '_test_MP_C', '_test_SP',
]

PMD_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPY_REF = os.path.join(PMD_ROOT, 'examples', 'results', 'npy_ref')
TMP_DIR = os.path.join(PMD_ROOT, 'tests', '_tmp')
os.makedirs(TMP_DIR, exist_ok=True)
WORKING = os.path.dirname(PMD_ROOT)

for m in models:
    t_path = os.path.join(TMP_DIR, f'{m}_T.npy')
    ut_path = os.path.join(TMP_DIR, f'{m}_uT.npy')
    script = (
        f'import matplotlib; matplotlib.use("Agg"); '
        f'import numpy as np; import importlib; '
        f'mod = importlib.import_module("examples.{m}"); '
        f'np.save(r"{t_path}", mod.T); '
        f'np.save(r"{ut_path}", mod.uT)'
    )
    print(f'Running {m}...', end=' ', flush=True)
    result = subprocess.run(
        [sys.executable, '-c', script],
        cwd=WORKING, capture_output=True, text=True, timeout=600
    )
    if result.returncode != 0:
        print(f'FAILED!\n  STDERR: {result.stderr[:500]}')
    else:
        shutil.copy2(t_path, os.path.join(NPY_REF, f'{m}_T.npy'))
        shutil.copy2(ut_path, os.path.join(NPY_REF, f'{m}_uT.npy'))
        os.remove(t_path)
        os.remove(ut_path)
        print('OK')

print('All done!')
