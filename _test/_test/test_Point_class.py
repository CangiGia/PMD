import pytest
import numpy as np
from PMD.maker import Point

def colvect(x, y):
    return np.array([[x], [y]])

# test creating an instance of Point
def test_point_creation():
    point = Point()
    assert point.Bindex == 0
    assert np.array_equal(point.sPlocal, colvect(0, 0))
    assert np.array_equal(point._sP, colvect(0, 0))
    assert np.array_equal(point._sPr, colvect(0, 0))
    assert np.array_equal(point._rP, colvect(0, 0))
    assert np.array_equal(point._dsP, colvect(0, 0))
    assert np.array_equal(point._drP, colvect(0, 0))
    assert np.array_equal(point._ddrP, colvect(0, 0))

# test with custom parameters
def test_point_custom_values():
    point = Point(Bindex=2, sPlocal=colvect(3, 4))
    assert point.Bindex == 2
    assert np.array_equal(point.sPlocal, colvect(3, 4))

# test if sPlocal and other vector properties are correctly converted to column vectors
def test_point_column_properties():
    point = Point(sPlocal=[5, 6])
    assert point.sPlocal.shape == (2, 1)
    assert np.array_equal(point.sPlocal, colvect(5, 6))

if __name__ == "__main__":
    pytest.main()