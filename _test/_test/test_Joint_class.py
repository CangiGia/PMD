import pytest
import numpy as np
from src.maker import Joint

def test_joint_initialization_default():
    # test default initialization of the Joint class
    joint = Joint()
    assert joint.type == 'rev'
    assert joint.iBindex == 0
    assert joint.jBindex == 0
    assert joint.iPindex == 0
    assert joint.jPindex == 0
    assert joint.iUindex == 0
    assert joint.jUindex == 0
    assert joint.iFunct == 0
    assert joint.L == 0
    assert joint.R == 1
    assert joint.x0 == 0
    assert joint.d0 == []
    assert joint.fix == 0
    assert joint._p0 == 0
    assert joint._nbody == 2
    assert joint._mrows == 2
    assert joint._rows == 0
    assert joint._rowe == 0
    assert joint._colis == 0
    assert joint._colie == 0
    assert joint._coljs == 0
    assert joint._colje == 0
    assert np.array_equal(joint._lagrange, np.zeros([3, 1]))

def test_joint_initialization_custom_values():
    # test custom initialization of the Joint class
    joint = Joint(type='tran', iBindex=1, jBindex=2, iPindex=3, jPindex=4, iUindex=5, jUindex=6, 
                  iFunct=7, L=8, R=9, x0=10, d0=np.array([1, 2, 3]), fix=1)
    assert joint.type == 'tran'
    assert joint.iBindex == 1
    assert joint.jBindex == 2
    assert joint.iPindex == 3
    assert joint.jPindex == 4
    assert joint.iUindex == 5
    assert joint.jUindex == 6
    assert joint.iFunct == 7
    assert joint.L == 8
    assert joint.R == 9
    assert joint.x0 == 10
    assert np.array_equal(joint.d0, np.array([1, 2, 3]))
    assert joint.fix == 1
    assert joint._p0 == 0
    assert joint._nbody == 2
    assert joint._mrows == 2
    assert joint._rows == 0
    assert joint._rowe == 0
    assert joint._colis == 0
    assert joint._colie == 0
    assert joint._coljs == 0
    assert joint._colje == 0
    assert np.array_equal(joint._lagrange, np.zeros([3, 1]))

def test_joint_d0_default():
    # test default value of d0
    joint = Joint()
    assert joint.d0 == []

def test_joint_d0_custom():
    # test custom value of d0
    joint = Joint(d0=np.array([1, 2, 3]))
    assert np.array_equal(joint.d0, np.array([1, 2, 3]))

def test_joint_lagrange_default():
    # test if _lagrange is initialized as a column vector with default value
    joint = Joint()
    assert np.array_equal(joint._lagrange, np.zeros([3, 1]))

def test_joint_lagrange_custom():
    # test custom value of _lagrange
    joint = Joint()
    joint._lagrange = np.array([[5], [6], [7]])
    assert np.array_equal(joint._lagrange, np.array([[5], [6], [7]]))

if __name__ == "__main__":
    pytest.main()