import pytest
import numpy as np
from PMD.maker import Function

def test_function_initialization_default():
    # test default initialization of the Function class
    func = Function()
    assert func.type == 'a'
    assert func.t_start == 0
    assert func.f_start == 0
    assert func.t_end == 1
    assert func.f_end == 1
    assert func.dfdt_end == 1
    assert func.ncoeff == 4
    assert func.coeff == []

def test_function_initialization_custom_values():
    # test custom initialization of the Function class
    coeff = np.array([1.0, 2.0, 3.0, 4.0])
    func = Function(type='b', t_start=2, f_start=3, t_end=4, f_end=5, dfdt_end=6, ncoeff=5, coeff=coeff)
    assert func.type == 'b'
    assert func.t_start == 2
    assert func.f_start == 3
    assert func.t_end == 4
    assert func.f_end == 5
    assert func.dfdt_end == 6
    assert func.ncoeff == 5
    assert np.array_equal(func.coeff, coeff)

def test_function_type_default():
    # test default value of type
    func = Function()
    assert func.type == 'a'

def test_function_type_custom():
    # test custom value of type
    func = Function(type='c')
    assert func.type == 'c'

def test_function_coeff_default():
    # test default value of coeff
    func = Function()
    assert func.coeff == []

def test_function_coeff_custom():
    # test custom value of coeff
    coeff = np.array([1.0, 2.0, 3.0, 4.0])
    func = Function(coeff=coeff)
    assert np.array_equal(func.coeff, coeff)

def test_function_ncoeff_default():
    # test default value of ncoeff
    func = Function()
    assert func.ncoeff == 4

def test_function_ncoeff_custom():
    # test custom value of ncoeff
    func = Function(ncoeff=6)
    assert func.ncoeff == 6

def test_function_dfdt_end_default():
    # test default value of dfdt_end
    func = Function()
    assert func.dfdt_end == 1

def test_function_dfdt_end_custom():
    # test custom value of dfdt_end
    func = Function(dfdt_end=10)
    assert func.dfdt_end == 10

if __name__ == "__main__":
    pytest.main()
