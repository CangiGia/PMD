"""
Regression test: compares current model outputs against saved reference data.
Usage: python -m PMD.examples.test_regression
Must be run from C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
"""
import subprocess
import sys
import os
import shutil
import numpy as np

models = [
    '_test_AA', '_test_Cart_A', '_test_Cart_B', '_test_Cart_C', '_test_Cart_D',
    '_test_CB', '_test_DP', '_test_MP_A', '_test_MP_B', '_test_MP_C', '_test_SP'
]

ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'npy_ref')
tmp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'npy_tmp')
os.makedirs(tmp_dir, exist_ok=True)

cwd = r'C:\Users\Giaco\anaconda3\envs\GiacoEnv'
all_passed = True

for model in models:
    print(f"\n{'='*60}")
    print(f"Testing {model}...")
    print(f"{'='*60}")

    # Run model in a separate subprocess, save to tmp
    result = subprocess.run(
        [sys.executable, '-m', 'PMD.examples._run_and_save', model, 'npy_tmp'],
        cwd=cwd,
        capture_output=True, text=True,
        timeout=600
    )

    if result.returncode != 0:
        print(f"  FAILED to run:")
        print(result.stderr[-500:] if len(result.stderr) > 500 else result.stderr)
        all_passed = False
        continue

    # Load reference and new results
    T_ref = np.load(os.path.join(ref_dir, f'{model}_T.npy'))
    T_new = np.load(os.path.join(tmp_dir, f'{model}_T.npy'))
    uT_ref = np.load(os.path.join(ref_dir, f'{model}_uT.npy'))
    uT_new = np.load(os.path.join(tmp_dir, f'{model}_uT.npy'))

    # Compare bit-for-bit
    t_match = np.array_equal(T_ref, T_new)
    ut_match = np.array_equal(uT_ref, uT_new)

    if t_match and ut_match:
        print(f"  PASSED (bit-for-bit identical)")
    else:
        all_passed = False
        if not t_match:
            diff_t = np.max(np.abs(T_ref - T_new))
            print(f"  FAILED: T differs, max_diff = {diff_t:.2e}")
        if not ut_match:
            diff_ut = np.max(np.abs(uT_ref - uT_new))
            print(f"  FAILED: uT differs, max_diff = {diff_ut:.2e}")

print(f"\n{'='*60}")
if all_passed:
    print("ALL 11 MODELS PASSED (bit-for-bit identical)")
else:
    print("SOME MODELS FAILED")
print(f"{'='*60}")

# Cleanup tmp
shutil.rmtree(tmp_dir, ignore_errors=True)

sys.exit(0 if all_passed else 1)
