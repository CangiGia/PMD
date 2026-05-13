# Contributor's Guide

Thank you for your interest in improving PMD!

## Reporting bugs

Open an [Issue](https://github.com/CangiGia/PMD/issues) and include:

- A short description of the problem.
- A minimal reproducible example (model definition + `solve()` call).
- The full traceback if an exception is raised.
- Your Python and PMD version (`pip show pmd`).

## Suggesting features

Open an [Issue](https://github.com/CangiGia/PMD/issues) labelled
`enhancement` describing the feature and its motivation.

## Submitting a pull request

1. Fork the repository on GitHub.
2. Clone your fork and create a feature branch:

   ```bash
   git checkout -b feat/my-feature
   ```

3. Install in editable mode with all dev extras:

   ```bash
   pip install -e ".[gui,dev]"
   ```

4. Make your changes and add tests in `tests/`.

5. Verify the test suite passes:

   ```bash
   pytest tests/
   ```

6. Build the docs and check for warnings:

   ```bash
   cd docs
   python -m sphinx -b html source _build/html -W
   ```

7. Push your branch and open a pull request against `main`.

## Improving the documentation

Documentation source files live in `docs/source/`. They use
[MyST Markdown](https://myst-parser.readthedocs.io/) and are built with
[Sphinx](https://www.sphinx-doc.org/). Edit the relevant `.md` or `.rst`
file, rebuild locally, and submit a pull request.
