"""
Generate reference .npy files for all test models.
Usage: python -m PMD.examples.generate_references
Must be run from C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
"""
import subprocess
import sys
import os

models = [
    '_test_AA', '_test_Cart_A', '_test_Cart_B', '_test_Cart_C', '_test_Cart_D',
    '_test_CB', '_test_DP', '_test_MP_A', '_test_MP_B', '_test_MP_C', '_test_SP'
]

cwd = r'C:\Users\Giaco\anaconda3\envs\GiacoEnv'
failed = []

for model in models:
    print(f"\n{'='*60}")
    print(f"Running {model}...")
    print(f"{'='*60}")
    result = subprocess.run(
        [sys.executable, '-m', 'PMD.examples._run_and_save', model, 'npy_ref'],
        cwd=cwd,
        capture_output=True, text=True,
        timeout=600  # 10 min timeout per model
    )
    print(result.stdout)
    if result.returncode != 0:
        print(f"ERROR in {model}:")
        print(result.stderr)
        failed.append(model)

print(f"\n{'='*60}")
if failed:
    print(f"FAILED models: {', '.join(failed)}")
    sys.exit(1)
else:
    # Verify file count
    ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'npy_ref')
    npy_files = [f for f in os.listdir(ref_dir) if f.endswith('.npy')]
    print(f"All reference data generated successfully!")
    print(f"Files in npy_ref/: {len(npy_files)} (expected {len(models) * 2})")
    sys.exit(0)
