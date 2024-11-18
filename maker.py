"""
This Python module encompasses all the essential components for constructing a 
planar multi-body dynamic model.

Author: Giacomo Cangi
"""

import numpy as np
from numpy.typing import *
from PMD.pmd_functions import *


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
    r_d (NDArray) 
        Time derivative of position (x, y).
    p_d (float) 
        Time derivative of the orientation angle (phi).
    A (NDArray) 
        Rotational transformation matrix.
    r_dd (NDArray) 
        Second derivative of position (x, y).
    p_dd (float) 
        Second time derivative of the orientation angle (phi).
    irc (int) 
        Index of the first element of r in u.
    irv (int) 
        Index of the first element of r_dot in u.
    ira (int) 
        Index of the first element of r_dot2 in v.
    m_inv (float) 
        Inverse of the mass.
    J_inv (float) 
        Inverse of the moment of inertia.
    wgt (NDArray) 
        Weight of the body as a force vector.
    f (NDArray) 
        Sum of forces acting on the body.
    n (float) 
        Sum of moments acting on the body.
    shape (str) 
        Shape of the body ('circle', 'rect', 'line').
    R (float) 
        Radius (for circular bodies).
    circ (list) 
        Points on the circumference of the body.
    W (float) 
        Width of the body (for rectangular bodies).
    H (float) 
        Height of the body (for rectangular bodies).
    color (str) 
        Color of the body for visualization.
    P4 (list) 
        List of 4 corners of the rectangle.
    pts (list) 
        List of point indexes associated with this body.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.m = 1
        self.J = 1
        self.r = colvect(0, 0)
        self.p = 0
        self.r_d = colvect(0, 0)
        self.p_d = 0
        self.A = np.eye(2)
        self.r_dd = colvect(0, 0)
        self.p_dd = 0
        self.irc = 0
        self.irv = 0
        self.ira = 0
        self.m_inv = 1
        self.J_inv = 1
        self.wgt = colvect(0, 0)
        self.f = colvect(0, 0)
        self.n = 0
        self.shape = ''
        self.R = 1
        self.circ = []
        self.W = 0
        self.H = 0
        self.color = 'k'
        self.P4 = []
        self.pts = []
    
    # new values will be automatically defined as column vector
    r = as_column_property("r")
    r_d = as_column_property("r_d")
    r_dd = as_column_property("r_dd")
    
class Point(Base): 
    """
    Create a body-fixed `point` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    Bindex (int)
        Body index to which the point belong.
    sPlocal (NDArray)
        Point coordinates with respect to the local reference frame.
    sP (NDArray)
        `x` and `y` components of the vector `s`.
    sP_r (NDArray)            
        `x` and `y` components of the vector `s` rotated.
    rP (NDArray)
        `x` and `y` coordinates of the point. 
    dsP (NDArray)
        First time derivative of the vector `s`. 
    drP (NDArray)
        `x` and `y` first time derivative of the vector `r`. 
    ddrP (NDArray)
        `x` and `y` second time derivative of the vector `r`.         
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.Bindex = 0
        self.sPlocal = colvect(0, 0)
        self.sP = colvect(0, 0)
        self.sP_r = colvect(0, 0)
        self.rP = colvect(0, 0)
        self.dsP = colvect(0, 0)
        self.drP = colvect(0, 0)
        self.ddrP = colvect(0, 0)

    # new values will be automatically defined as column vector
    sPlocal = as_column_property("sPlocal")
    sP = as_column_property("sP")
    sP_r = as_column_property("sP_r")
    rP = as_column_property("rP")
    dsP = as_column_property("dsP")
    drP = as_column_property("drP")
    ddrP = as_column_property("ddrP")
    
class uVector(Base): 
    """
    Crate a body-fixed `unit vector` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    Bindex (int)
        Body index to which the vector belong.
    ulocal (NDArray)
        `xi` and `eta` components of the vector in the local reference frame.
    u (NDArray)
        `x` and `y` componenets of the vector in the global reference frame.
    u_r (NDArray)
        Components of the vector `u`, rotated. 
    u_d (NDArray)
        First time derivative of the vector `u`.
    """

    def __init__(self):
        super().__init__()  # call to the Base class constructor
        self.Bindex = 0
        self.ulocal = colvect(1, 0)
        self.u = colvect(0, 0)
        self.u_r = colvect(0, 0)
        self.u_d = colvect(0, 0)
    
    # new values will be automatically defined as column vector
    ulocal = as_column_property("ulocal")
    u = as_column_property("u")
    u_r = as_column_property("u_r")
    u_d = as_column_property("u_d")

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
    gravity (float)
        Gravitational constant.
    wgt (NDArray)
        Gravitational direction (default `-y`)
    flocal (NDArray)
        Constant force in the local reference frame.
    f (NDArray)
        Constant force in the `x-y` reference frame.
    T (float)
        Constant torque in the `x-y` refernce frame.
    iFunct (int)
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
        self.gravity = Force.DEFAULT_GRAVITY  
        self.wgt = Force.DEFAULT_GRAVITY_VECTOR 
        self.flocal = colvect([0, 0])
        self.f = colvect([0, 0])
        self.T = 0
        self.iFunct = 0
    
    # new values will be automatically defined as column vector
    wgt = as_column_property("wgt")
    flocal = as_column_property("flocal")
    f = as_column_property("f")

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
    p0 (float)
        Initial condition `phi` for a disc (or rigid).
    d0 (NDArray)
        Initial condition for `d` (rigid). 
    fix (int)
        Fix relative dof (if = 1, rev or tran). 
    nbody (int)
        Number of moving bodies involved.
    mrows (int)
        Number of rows (constraints).
    rows (int)
        Row index start.
    rowe (int)
        Row index end.
    colis (int)
        Comlumn index for body i-start. 
    colie (int)
        Column index for body j-start. 
    coljs (int)
        Column index for body j-start.
    colje (int)
        Column index for body i-start. 
    lagrange (NDArray)
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
        self.p0 = 0
        self.d0 = []
        self.fix = 0
        self.nbody = 2
        self.mrows = 2
        self.rows = 0
        self.rowe = 0
        self.colis = 0
        self.colie = 0
        self.coljs = 0
        self.colje = 0
        self.lagrange = np.zeros([3,1])
    
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