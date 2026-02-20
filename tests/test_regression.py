"""
Regression tests for PMD models.

Each model is run in a subprocess and its outputs (T, uT) are compared
bit-for-bit with golden reference .npy files.

Run from: C:\\Users\\Giaco\\anaconda3\\envs\\GiacoEnv\\
Command:  python -m pytest PMD/tests/test_regression.py -v
"""
import pytest
import numpy as np

from .conftest import ALL_MODELS, run_model_subprocess, load_reference


class TestModelRegression:
    """Regression tests: each model must reproduce its golden reference exactly."""

    @pytest.fixture()
    def setup_model(self, request):
        """
        Run the model and load its reference data.

        Stores results as instance attributes:
            self.T, self.uT       - current run
            self.T_ref, self.uT_ref - golden reference
        """
        model_name = request.param
        self.T, self.uT = run_model_subprocess(model_name)
        self.T_ref, self.uT_ref = load_reference(model_name)

    @pytest.mark.parametrize('setup_model', ALL_MODELS, indirect=True)
    def test_time_vector_identical(self, setup_model):
        """Time vector T must be bit-for-bit identical to reference."""
        np.testing.assert_array_equal(
            self.T, self.T_ref,
            err_msg="Time vector T differs from reference"
        )

    @pytest.mark.parametrize('setup_model', ALL_MODELS, indirect=True)
    def test_state_matrix_identical(self, setup_model):
        """State matrix uT must be bit-for-bit identical to reference."""
        np.testing.assert_array_equal(
            self.uT, self.uT_ref,
            err_msg="State matrix uT differs from reference"
        )
