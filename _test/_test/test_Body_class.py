import pytest
import numpy as np
from src.maker import Base, Body

def colvect(x, y):
    return np.array([[x], [y]])

# test creating an instance of Body
def test_body_creation():
    body = Body()
    assert body.m == 1
    assert body.J == 1
    assert np.array_equal(body.r, colvect(0, 0))
    assert body.p == 0
    assert np.array_equal(body.dr, colvect(0, 0))
    assert body.dp == 0
    assert np.array_equal(body.ddr, colvect(0, 0))
    assert body.ddp == 0
    assert np.array_equal(body._A, np.eye(2))
    assert body._invm == 1.0
    assert body._invJ == 1.0
    assert np.array_equal(body._wgt, colvect(0, 0))
    assert np.array_equal(body._f, colvect(0, 0))
    assert body._n == 0
    assert body._pts == []

# test with custom parameters
def test_body_custom_values():
    body = Body(m=2, J=4, r=colvect(1, 2), p=0.5, dr=colvect(3, 4), dp=1, ddr=colvect(5, 6), ddp=2)
    assert body.m == 2
    assert body.J == 4
    assert np.array_equal(body.r, colvect(1, 2))
    assert body.p == 0.5
    assert np.array_equal(body.dr, colvect(3, 4))
    assert body.dp == 1
    assert np.array_equal(body.ddr, colvect(5, 6))
    assert body.ddp == 2
    assert body._invm == 0.5
    assert body._invJ == 0.25

# test for error if m is zero
def test_body_mass_zero():
    with pytest.raises(ValueError, match=r"Mass \(m\) cannot be zero"):
        Body(m=0)

# test if r, dr, and ddr properties are correctly converted to column vectors
def test_body_column_properties():
    body = Body(r=[1, 2], dr=[3, 4], ddr=[5, 6])
    assert body.r.shape == (2, 1)
    assert body.dr.shape == (2, 1)
    assert body.ddr.shape == (2, 1)
    assert np.array_equal(body.r, colvect(1, 2))
    assert np.array_equal(body.dr, colvect(3, 4))
    assert np.array_equal(body.ddr, colvect(5, 6))

if __name__ == "__main__":
    pytest.main()