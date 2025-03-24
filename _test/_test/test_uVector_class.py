import pytest
import numpy as np
from src.maker import uVector

def colvect(x, y):
    return np.array([[x], [y]])

# test creating an instance of uVector
def test_uvector_creation():
    uvector = uVector()
    assert uvector.Bindex == 0
    assert np.array_equal(uvector.ulocal, colvect(1, 0))
    assert np.array_equal(uvector._u, colvect(0, 0))
    assert np.array_equal(uvector._ur, colvect(0, 0))
    assert np.array_equal(uvector._du, colvect(0, 0))

# test with custom parameters
def test_uvector_custom_values():
    uvector = uVector(Bindex=3, ulocal=colvect(2, 3))
    assert uvector.Bindex == 3
    assert np.array_equal(uvector.ulocal, colvect(2, 3))

# test if ulocal and other vector properties are correctly converted to column vectors
def test_uvector_column_properties():
    uvector = uVector(ulocal=[4, 5])
    assert uvector.ulocal.shape == (2, 1)
    assert np.array_equal(uvector.ulocal, colvect(4, 5))

if __name__ == "__main__":
    pytest.main()