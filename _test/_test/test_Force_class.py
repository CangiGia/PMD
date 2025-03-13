import pytest
import numpy as np
from PMD.maker import Force

def colvect(x, y):
    return np.array([[x], [y]])

def test_force_initialization_default():
    # test default initialization of the Force class
    force = Force()
    assert force.type == 'ptp'
    assert force.iPindex == 0
    assert force.jPindex == 0
    assert force.iBindex == 0
    assert force.jBindex == 0
    assert force.k == 0
    assert force.L0 == 0
    assert force.theta0 == 0
    assert force.dc == 0
    assert force.f_a == 0
    assert force.T_a == 0
    assert np.array_equal(force.flocal, colvect(0, 0))  # Expect column vector
    assert np.array_equal(force.f, colvect(0, 0))  # Expect column vector
    assert force.T == 0
    assert force._gravity == Force.DEFAULT_GRAVITY
    assert np.array_equal(force._wgt, Force.DEFAULT_GRAVITY_VECTOR)
    assert force._iFunct == 0

def test_force_initialization_custom_values():
    # test custom initialization of the Force class
    force = Force(type='f', iPindex=1, jPindex=2, iBindex=3, jBindex=4, k=10, L0=5, theta0=0.1, dc=0.5, f_a=15, T_a=20, flocal=np.array([1, 1]), f=np.array([2, 2]), T=25)
    assert force.type == 'f'
    assert force.iPindex == 1
    assert force.jPindex == 2
    assert force.iBindex == 3
    assert force.jBindex == 4
    assert force.k == 10
    assert force.L0 == 5
    assert force.theta0 == 0.1
    assert force.dc == 0.5
    assert force.f_a == 15
    assert force.T_a == 20
    assert np.array_equal(force.flocal, colvect(1, 1))  # Expect column vector
    assert np.array_equal(force.f, colvect(2, 2))  # Expect column vector
    assert force.T == 25
    assert force._gravity == Force.DEFAULT_GRAVITY
    assert np.array_equal(force._wgt, Force.DEFAULT_GRAVITY_VECTOR)
    assert force._iFunct == 0

def test_force_flocal_default():
    # test if flocal is initialized as a column vector with default value
    force = Force()
    assert np.array_equal(force.flocal, colvect(0, 0))  # Expect column vector

def test_force_flocal_custom():
    # test if flocal is initialized correctly with custom value
    force = Force(flocal=np.array([3, 4]))
    assert np.array_equal(force.flocal, colvect(3, 4))  # Expect column vector

def test_force_f_default():
    # test if f is initialized as a column vector with default value
    force = Force()
    assert np.array_equal(force.f, colvect(0, 0))  # Expect column vector

def test_force_f_custom():
    # test if f is initialized correctly with custom value
    force = Force(f=np.array([7, 8]))
    assert np.array_equal(force.f, colvect(7, 8))  # Expect column vector

def test_force_wgt_default():
    # test if _wgt is initialized with default value
    force = Force()
    assert np.array_equal(force._wgt, Force.DEFAULT_GRAVITY_VECTOR)

def test_force_gravity_constant():
    # test if _gravity is initialized with the default gravity value
    force = Force()
    assert force._gravity == Force.DEFAULT_GRAVITY

def test_force_column_property_wgt():
    # test the as_column_property behavior for wgt
    force = Force()
    force.wgt = np.array([1, 2])
    assert np.array_equal(force.wgt, colvect(1, 2))  # Expect column vector

def test_force_column_property_flocal():
    # test the as_column_property behavior for flocal
    force = Force()
    force.flocal = np.array([3, 4])
    assert np.array_equal(force.flocal, colvect(3, 4))  # Expect column vector

def test_force_column_property_f():
    # test the as_column_property behavior for f
    force = Force()
    force.f = np.array([5, 6])
    assert np.array_equal(force.f, colvect(5, 6))  # Expect column vector
    
if __name__ == "__main__":
    pytest.main()
