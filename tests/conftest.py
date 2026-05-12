"""Shared fixtures for the pmd test-suite."""
import pytest


@pytest.fixture(autouse=True)
def _reset_ground_markers():
    """Reset ``Ground._markers`` to ``[origin]`` between tests.

    Prevents marker accumulation when multiple models are imported inside
    the same process (e.g. all ``test_numeric_snapshots`` cases).
    """
    from pmd.core.model import _GroundType
    yield
    gt = _GroundType._instance
    if gt is not None:
        _GroundType._markers = [gt.origin]
