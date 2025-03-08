import unittest
import numpy as np
from PMD.maker import Base, Body, Point, uVector, Force, Joint, Function

class TestBase(unittest.TestCase):
    def test_instance_count(self):
        initial_count = Base.get_count()
        instance = Base()
        self.assertEqual(Base.get_count(), initial_count + 1)
        del instance
        self.assertEqual(Base.get_count(), initial_count + 1)

    def test_get_type(self):
        instance = Base()
        self.assertEqual(instance.get_type(), 'Base')

class TestBody(unittest.TestCase):
    def test_initialization(self):
        body = Body(m=5.0, J=4.0, r=np.array([[1.0], [0.2]]))
        self.assertEqual(body.m, 5.0)
        self.assertEqual(body.J, 4.0)
        np.testing.assert_array_equal(body.r, np.array([[1.0], [0.2]]))
        self.assertEqual(body.p, 0)
        np.testing.assert_array_equal(body.dr, np.array([[0.0], [0.0]]))
        self.assertEqual(body.dp, 0)
        np.testing.assert_array_equal(body.ddr, np.array([[0.0], [0.0]]))
        self.assertEqual(body.ddp, 0)
        np.testing.assert_array_equal(body._A, np.eye(2))
        self.assertEqual(body._irc, 0)
        self.assertEqual(body._irv, 0)
        self.assertEqual(body._ira, 0)
        self.assertEqual(body._invm, 1 / 5.0)
        self.assertEqual(body._invJ, 1 / 4.0)
        np.testing.assert_array_equal(body._wgt, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(body._f, np.array([[0.0], [0.0]]))
        self.assertEqual(body._n, 0)
        self.assertEqual(body._pts, [])

    def test_invalid_mass(self):
        with self.assertRaises(ValueError):
            Body(m=0)

class TestPoint(unittest.TestCase):
    def test_initialization(self):
        point = Point()
        self.assertEqual(point.Bindex, 0)
        np.testing.assert_array_equal(point.sPlocal, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._sP, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._sPr, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._rP, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._dsP, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._drP, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(point._ddrP, np.array([[0.0], [0.0]]))

class TestuVector(unittest.TestCase):
    def test_initialization(self):
        uvector = uVector()
        self.assertEqual(uvector.Bindex, 0)
        np.testing.assert_array_equal(uvector.ulocal, np.array([[1.0], [0.0]]))
        np.testing.assert_array_equal(uvector._u, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(uvector._ur, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(uvector._du, np.array([[0.0], [0.0]]))

class TestForce(unittest.TestCase):
    def test_initialization(self):
        force = Force()
        self.assertEqual(force.type, 'ptp')
        self.assertEqual(force.iPindex, 0)
        self.assertEqual(force.jPindex, 0)
        self.assertEqual(force.iBindex, 0)
        self.assertEqual(force.jBindex, 0)
        self.assertEqual(force.k, 0)
        self.assertEqual(force.L0, 0)
        self.assertEqual(force.theta0, 0)
        self.assertEqual(force.dc, 0)
        self.assertEqual(force.f_a, 0)
        self.assertEqual(force.T_a, 0)
        np.testing.assert_array_equal(force.flocal, np.array([[0.0], [0.0]]))
        np.testing.assert_array_equal(force.f, np.array([[0.0], [0.0]]))
        self.assertEqual(force.T, 0)
        self.assertEqual(force._gravity, Force.DEFAULT_GRAVITY)
        np.testing.assert_array_equal(force._wgt, Force.DEFAULT_GRAVITY_VECTOR)
        self.assertEqual(force._iFunct, 0)

class TestJoint(unittest.TestCase):
    def test_initialization(self):
        joint = Joint()
        self.assertEqual(joint.type, 'rev')
        self.assertEqual(joint.iBindex, 0)
        self.assertEqual(joint.jBindex, 0)
        self.assertEqual(joint.iPindex, 0)
        self.assertEqual(joint.jPindex, 0)
        self.assertEqual(joint.iUindex, 0)
        self.assertEqual(joint.jUindex, 0)
        self.assertEqual(joint.iFunct, 0)
        self.assertEqual(joint.L, 0)
        self.assertEqual(joint.R, 1)
        self.assertEqual(joint.x0, 0)
        self.assertEqual(joint.d0, [])
        self.assertEqual(joint.fix, 0)
        self.assertEqual(joint._p0, 0)
        self.assertEqual(joint._nbody, 2)
        self.assertEqual(joint._mrows, 2)
        self.assertEqual(joint._rows, 0)
        self.assertEqual(joint._rowe, 0)
        self.assertEqual(joint._colis, 0)
        self.assertEqual(joint._colie, 0)
        self.assertEqual(joint._coljs, 0)
        self.assertEqual(joint._colje, 0)
        np.testing.assert_array_equal(joint._lagrange, np.zeros([3, 1]))

class TestFunction(unittest.TestCase):
    def test_initialization(self):
        function = Function()
        self.assertEqual(function.type, 'a')
        self.assertEqual(function.t_start, 0)
        self.assertEqual(function.f_start, 0)
        self.assertEqual(function.t_end, 1)
        self.assertEqual(function.f_end, 1)
        self.assertEqual(function.dfdt_end, 1)
        self.assertEqual(function.ncoeff, 4)
        np.testing.assert_array_equal(function.coeff, np.array([]))

if __name__ == '__main__':
    unittest.main()