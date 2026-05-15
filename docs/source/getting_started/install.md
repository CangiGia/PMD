# Installation

## Requirements

- Python ≥ 3.10
- NumPy, CasADi, tqdm
- PyQt6 (optional — required for the GUI)

## Install from source

PMD is not yet published on PyPI.  Install directly from the GitHub repository:

```bash
git clone https://github.com/CangiGia/PMD.git
cd PMD
pip install -e ".[gui,dev]"
```

The `gui` extra pulls in PyQt6 and related packages.
The `dev` extra adds pytest, Sphinx, and other development tools.

## Verify the installation

```python
import pmd
print(pmd.__version__)   # should print the current version
```

Or run the test suite:

```bash
pytest tests/
```
