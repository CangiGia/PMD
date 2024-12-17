"""
- Planar Multi-body Dynamics simulation maker - 
== 

This Python module encompasses all the essential components for constructing a 
planar multi-body dynamic model.

Author: - Giacomo Cangi, PhD student @ UniPG -
"""

import numpy as np
from numpy.typing import *
from pmd_functions import *


class Base:
    """
    Base class for all multi-body simulation objects.

    This class provides a foundation for objects in the multi-body dynamic 
    simulation framework. It includes functionality for instance counting 
    and type retrieval, supporting introspection and organized tracking 
    of different object types.
    """
    COUNT = 0  # class variable to track the number of instances

    def __new__(cls, *args, **kwargs):
        """
        Create a new instance of a subclass of `Base` and increment the 
        instance count.

        Parameters
        ----------
        args (tuple)
            Positional arguments passed to the constructor.
        kwargs (dict)
            Keyword arguments passed to the constructor.

        Returns
        -------
        instance
            A new instance of the class calling this method.
        """
        instance = super().__new__(cls)
        cls.COUNT += 1
        return instance

    @classmethod
    def get_count(cls):
        """
        Retrieve the current count of instances.

        Returns
        -------
        int
            The total count of instances derived from `Base`.
        """
        return cls.COUNT
    
    def get_type(self):
        """
        Get the name of the class of the instance.

        Returns
        -------
        str
            The class name of the instance.
        """
        return self.__class__.__name__

class Body(Base):
    """
    Create a `rigid body` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    m (float) 
        Mass of the body.
    J (float) 
        Moment of inertia of the body.
    r (NDArray) 
        Position vector (x, y) of the body.
    p (float) 
        Orientation angle (phi) of the body.
    dr (NDArray) 
        Time derivative of position (x, y).
    dp (float) 
        Time derivative of the orientation angle (phi).
    ddr (NDArray) 
        Second derivative of position (x, y).
    ddp (float) 
        Second time derivative of the orientation angle (phi).

    Notes
    -----
    _A (NDArray) 
        Rotational transformation matrix.
    _irc (int) 
        Index of the first element of r in u.
    _irv (int) 
        Index of the first element of r_dot in u.
    _ira (int) 
        Index of the first element of r_dot2 in v.
    _invm (float) 
        Inverse of the mass.
    _invJ (float) 
        Inverse of the moment of inertia.
    _wgt (NDArray) 
        Weight of the body as a force vector.
    _f (NDArray) 
        Sum of forces acting on the body.
    _n (float) 
        Sum of moments acting on the body.
    _pts (list) 
        List of point indexes associated with this body.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.m = 1
        self.J = 1
        self.r = colvect(0, 0)
        self.p = 0
        self.dr = colvect(0, 0)    # default value
        self.dp = 0                # default value
        self.ddr = colvect(0, 0)   # default value
        self.ddp = 0               # default value
        self._A = np.eye(2)
        self._irc = 0
        self._irv = 0
        self._ira = 0
        self._invm = 1
        self._invJ = 1
        self._wgt = colvect(0, 0)
        self._f = colvect(0, 0)
        self._n = 0
        self._pts = []
    
    # new values will be automatically treated as column vectors, 
    # even if the user defines them as row vectors or lists!
    r = as_column_property("r")
    dr = as_column_property("dr")
    ddr = as_column_property("ddr")
    
class Point(Base): 
    """
    Create a body-fixed `point` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    Bindex (int)
        Body index to which the point belong.
    sPlocal (NDArray)
        Point coordinates with respect to the local reference frame.
    
    Notes
    -----
    _sP (NDArray)
        `x` and `y` components of the vector `s`.
    _sPr (NDArray)            
        `x` and `y` components of the vector `s` rotated.
    _rP (NDArray)
        `x` and `y` coordinates of the point. 
    _dsP (NDArray)
        First time derivative of the vector `s`. 
    _drP (NDArray)
        `x` and `y` first time derivative of the vector `r`. 
    _ddrP (NDArray)
        `x` and `y` second time derivative of the vector `r`.         
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.Bindex = 0
        self.sPlocal = colvect(0, 0)
        self._sP = colvect(0, 0)
        self._sPr = colvect(0, 0)
        self._rP = colvect(0, 0)
        self._dsP = colvect(0, 0)
        self._drP = colvect(0, 0)
        self._ddrP = colvect(0, 0)

    # new values will be automatically treated as column vectors, 
    # even if the user defines them as row vectors or lists!
    sPlocal = as_column_property("sPlocal")
    _sP = as_column_property("_sP")
    _sPr = as_column_property("_sPr")
    _rP = as_column_property("_rP")
    _dsP = as_column_property("_dsP")
    _drP = as_column_property("_drP")
    _ddrP = as_column_property("_ddrP")
    
class uVector(Base): 
    """
    Crate a body-fixed `unit vector` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    Bindex (int)
        Body index to which the vector belong.
    ulocal (NDArray)
        `xi` and `eta` components of the vector in the local reference frame.

    Notes
    -----
    _u (NDArray)
        `x` and `y` componenets of the vector in the global reference frame.
    _ur (NDArray)
        Components of the vector `u`, rotated. 
    _du (NDArray)
        First time derivative of the vector `u`.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.Bindex = 0
        self.ulocal = colvect(1, 0) # default value
        self._u = colvect(0, 0)
        self._ur = colvect(0, 0)
        self._du = colvect(0, 0)
    
    # new values will be automatically treated as column vectors, 
    # even if the user defines them as row vectors or lists!
    ulocal = as_column_property("ulocal")
    _u = as_column_property("_u")
    _ur = as_column_property("_ur")
    _du = as_column_property("_du")

class Force(Base):
    """
    Crate a force in a planar multi-body dynamic simulation.

    Attributes
    ----------
    type (str)
        Element type: ptp, rot_sda, weight, fp, f, T.
    iPindex (int)
        Index of the head (arrow) point.
    jPindex (int)
        Index of the tail point.
    iBindex (int)
        Index of the head (arrow) body.
    jBindex (int)
        Index of the tail body.
    k (float)
        Spring stiffness.
    L0 (float)
        Undeformed spring length.
    theta0 (float)
        Undeformed torsional spring angle.
    dc (float)
        Damping coefficient.
    f_a (float)
        Constant actuator force.
    T_a (float)
        Constant actuator torque.
    flocal (NDArray)
        Constant force in the local reference frame.

    Notes
    -----        
    _gravity (float)
        Gravitational constant.
    _wgt (NDArray)
        Gravitational direction (default `-y`).
    _f (NDArray)
        Constant force in the `x-y` reference frame.
    _T (float)
        Constant torque in the `x-y` refernce frame.
    _iFunct (int)
        Analytical function index.
    """

    # class-level constants
    DEFAULT_GRAVITY = 9.81
    DEFAULT_GRAVITY_VECTOR = colvect([0, -1])  # default gravitational force vector

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.type = 'ptp'
        self.iPindex = 0
        self.jPindex = 0
        self.iBindex = 0
        self.jBindex = 0
        self.k = 0
        self.L0 = 0
        self.theta0 = 0
        self.dc = 0
        self.f_a = 0
        self.T_a = 0
        self.flocal = colvect([0, 0])
        self._gravity = Force.DEFAULT_GRAVITY  
        self._wgt = Force.DEFAULT_GRAVITY_VECTOR 
        self._f = colvect([0, 0])
        self._T = 0
        self._iFunct = 0
    
    # new values will be automatically defined as column vector
    wgt = as_column_property("wgt")
    flocal = as_column_property("flocal")
    _f = as_column_property("f")

class Joint(Base):
    """
    Create a joint in a planar multi-body dynamic simulation.

    Attributes
    ----------
    type (str)
        Joint type: `rev`, `tran`, `rev-rev`, `rev-tran`, `rigid`, `disc`, 
        `rel-rot`, `rel-tran`.
    iBindex (int)
        Index of body `i`. 
    jBindex (int)
        Index of body `j`. 
    iPindex (int)
        Index of point `Pi`. 
    jPindex (int)
        Index of point `Pj`.
    iUindex (int)
        Index of unit vector `u_i`.
    jUindex (int)
        Index of unit vector `u_i`.
    iFunct (int)
        Analytical function index. 
    L (float)
        Constant length.
    R (float)
        Constant radius.
    x0 (float)
        Initial condition `x`, for disc.
    d0 (NDArray)
        Initial condition for `d` (rigid). 
    fix (int)
        Fix relative dof (if = 1, rev or tran). 
    _p0 (float)
        Initial condition `phi` for a disc (or rigid).
    _nbody (int)
        Number of moving bodies involved.
    _mrows (int)
        Number of rows (constraints).
    _rows (int)
        Row index start.
    _rowe (int)
        Row index end.
    _colis (int)
        Comlumn index for body i-start. 
    _colie (int)
        Column index for body j-start. 
    _coljs (int)
        Column index for body j-start.
    _colje (int)
        Column index for body i-start. 
    _lagrange (NDArray)
        Lagrange multipliers.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.type = 'rev'
        self.iBindex = 0
        self.jBindex = 0 
        self.iPindex = 0
        self.jPindex = 0
        self.iUindex = 0
        self.jUindex = 0
        self.iFunct = 0
        self.L = 0
        self.R = 1
        self.x0 = 0
        self.d0 = []
        self.fix = 0
        self._p0 = 0
        self._nbody = 2
        self._mrows = 2
        self._rows = 0
        self._rowe = 0
        self._colis = 0
        self._colie = 0
        self._coljs = 0
        self._colje = 0
        self._lagrange = np.zeros([3,1])
    
class Function(Base):
    """
    Create a user defined function in a planar multi-body dynamic simulation.
    
    Attributes
    ----------
    type (str)
        Function type `a`, `b` or `c`.
    t_start (float)
        Initial value for the time vector, i.e. the x axis of the function. 
        Required for funtction type `b` and `c`.
    f_start (float)
        Initial value for the function, i.e. initial value of the y axis. 
        Required for funtction type `b` and `c`.
    t_end (float)
        Final value for the time vector, i.e. the x axis of the function. 
        Required for funtction type `b` and `c`.
    f_end (float)
        Final value for the function, i.e. initial value of the y axis. 
        Required for funtction type `b` and `c`.
    dfdt_end (float) 
        Max value for the first time derivative function. Required for function 
        type `c`. 
    ncoeff (int) 
        Number of coefficients. 
    coeff (NDArray)
        Coefficient required only for function type `a`.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.type = 'a'
        self.t_start = 0
        self.f_start = 0 
        self.t_end = 1
        self.f_end = 1
        self.dfdt_end = 1
        self.ncoeff = 4
        self.coeff = []