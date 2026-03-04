"""
Planar Multi-Body Dynamics Model Builder

This module provides essential tools and algorithms for constructing 
and simulating planar multi-body dynamic systems. It includes support 
for defining rigid bodies, joints, constraints, and external forces.

Author: Giacomo Cangi
"""


import numpy as np
from numpy.typing import *
from .utils import *


# ── Ground singleton ──────────────────────────────────────────────

class _GroundType:
    """Singleton representing the inertial ground/world frame.

    Ground is the immovable reference body with zero state.  Use the
    module-level ``Ground`` instance instead of instantiating this class.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    # Fixed zero state — mirrors the Body interface that the solver reads.
    r   = np.zeros((2, 1))
    p   = 0.0
    dr  = np.zeros((2, 1))
    dp  = 0.0
    ddr = np.zeros((2, 1))
    ddp = 0.0
    _A  = np.eye(2)
    _bidx = 0          # internal body index: 0 ≡ ground

    def __repr__(self):
        return "Ground"

    def __bool__(self):
        """Ground is falsy so ``if body:`` is False for Ground, True for Body."""
        return False


Ground = _GroundType()


# ── Base class ────────────────────────────────────────────────────

class Base:
    """
    Base class for all multi-body simulation objects.

    This class provides a foundation for objects in the multi-body dynamic 
    simulation framework.
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
    _bidx (int)
        Internal body index (1-based), set by the solver during initialization.
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
        List of Point objects associated with this body.
    """

    def __init__(self, m=1, J=1, r=None, p=0, dr=None, dp=0, ddr=None, ddp=0):
        super().__init__()  # call to the Base class constructor
        if m <= 0:
            raise ValueError(f"Body {self.COUNT}: mass must be positive, got {m}")
        if J < 0:
            raise ValueError(f"Body {self.COUNT}: moment of inertia cannot be negative, got {J}")
        self.m = m
        self.J = J
        self.r = r if r is not None else colvect(0, 0)
        self.p = p
        self.dr = dr if dr is not None else colvect(0, 0)
        self.dp = dp
        self.ddr = ddr if ddr is not None else colvect(0, 0)
        self.ddp = ddp
        self._A = np.eye(2)
        self._bidx = 0
        self._irc = 0
        self._irv = 0
        self._ira = 0
        self._invm = 1.0 / m
        self._invJ = 1.0 / J if J != 0 else float('inf')
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
    body : Body or Ground
        The body to which this point is attached.
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

    def __init__(self, body=None, sPlocal=None):
        super().__init__()  # call to the Base class constructor
        self.body = body if body is not None else Ground
        self.sPlocal = sPlocal if sPlocal is not None else colvect(0, 0)
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
    Create a body-fixed `unit vector` in a planar multi-body dynamic simulation.

    Attributes
    ----------
    body : Body or Ground
        The body to which this vector is attached.
    ulocal (NDArray)
        `xi` and `eta` components of the vector in the local reference frame.

    Notes
    -----
    _u (NDArray)
        `x` and `y` components of the vector in the global reference frame.
    _ur (NDArray)
        Components of the vector `u`, rotated. 
    _du (NDArray)
        First time derivative of the vector `u`.
    """

    def __init__(self, body=None, ulocal=None):
        super().__init__()  # call to the Base class constructor
        self.body = body if body is not None else Ground
        self.ulocal = ulocal if ulocal is not None else colvect(1, 0) # default value
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
    Create a force in a planar multi-body dynamic simulation.

    Attributes
    ----------
    type : str
        Element type: ``'ptp'``, ``'rot-sda'``, ``'weight'``,
        ``'flocal'``, ``'f'``, ``'T'``, ``'user'``.
    iPoint : Point or None
        Head (arrow) point.
    jPoint : Point or None
        Tail point.
    iBody : Body, Ground, or None
        Head (arrow) body.
    jBody : Body, Ground, or None
        Tail body.
    k : float
        Spring stiffness.
    L0 : float
        Undeformed spring length.
    theta0 : float
        Undeformed torsional spring angle.
    dc : float
        Damping coefficient.
    f_a : float
        Constant actuator force.
    T_a : float
        Constant actuator torque.
    flocal : NDArray
        Constant force in the local reference frame.
    f : NDArray
       Constant force in the ``x-y`` reference frame.
    T : float
        Constant torque in the ``x-y`` reference frame.
    callback : callable or None
        User-defined force callback (no arguments, uses closures).

    Notes
    -----        
    _gravity : float
        Gravitational constant.
    _wgt : NDArray
        Gravitational direction (default ``-y``).
    _iFunct : int
        Analytical function index.
    """

    # class-level constants
    DEFAULT_GRAVITY = 9.81
    DEFAULT_GRAVITY_VECTOR = colvect([0, -1])  # default gravitational force vector

    VALID_TYPES = ['ptp', 'rot-sda', 'weight', 'flocal', 'f', 'T', 'user']

    def __init__(self, type='ptp', iPoint=None, jPoint=None,
                 iBody=None, jBody=None,
                 k=0, L0=0, theta0=0, dc=0, f_a=0, T_a=0,
                 flocal=None, f=None, T=0, callback=None):
        super().__init__()  # call to the Base class constructor
        if type not in Force.VALID_TYPES:
            raise ValueError(f"Force {self.COUNT}: unknown type '{type}', valid types: {Force.VALID_TYPES}")
        self.type = type
        self.iPoint = iPoint
        self.jPoint = jPoint
        self.iBody = iBody if iBody is not None else Ground
        self.jBody = jBody if jBody is not None else Ground
        self.k = k
        self.L0 = L0
        self.theta0 = theta0
        self.dc = dc
        self.f_a = f_a
        self.T_a = T_a
        self.flocal = flocal if flocal is not None else colvect([0, 0])
        self.f = f if f is not None else colvect([0, 0])
        self.T = T
        self.callback = callback
        self._gravity = Force.DEFAULT_GRAVITY  
        self._wgt = Force.DEFAULT_GRAVITY_VECTOR
        self._iFunct = 0

    # new values will be automatically defined as column vector
    wgt = as_column_property("wgt")
    flocal = as_column_property("flocal")
    f = as_column_property("f")

class Joint(Base):
    """
    Create a joint in a planar multi-body dynamic simulation.

    Attributes
    ----------
    type : str
        Joint type: ``'rev'``, ``'tran'``, ``'rev-rev'``, ``'rev-tran'``,
        ``'rigid'``, ``'disc'``, ``'rel-rot'``, ``'rel-tran'``.
    iBody : Body or Ground
        Body ``i``.
    jBody : Body or Ground
        Body ``j``.
    iPoint : Point or None
        Point ``Pi``.
    jPoint : Point or None
        Point ``Pj``.
    iUvec : uVector or None
        Unit vector ``ui``.
    jUvec : uVector or None
        Unit vector ``uj``.
    iFunct : Function or None
        Analytical function object.
    L : float
        Constant length.
    R : float
        Constant radius.
    x0 : float
        Initial condition ``x``, for disc.
    d0 : NDArray
        Initial condition for ``d`` (rigid). 
    fix : int
        Fix relative dof (if = 1, rev or tran). 

    Notes
    -----    
    _p0 : float
        Initial condition ``phi`` for a disc (or rigid).
    _nbody : int
        Number of moving bodies involved.
    _mrows : int
        Number of rows (constraints).
    _rows : int
        Row index start.
    _rowe : int
        Row index end.
    _colis : int
        Column index for body i-start. 
    _colie : int
        Column index for body i-end. 
    _coljs : int
        Column index for body j-start.
    _colje : int
        Column index for body j-end. 
    _lagrange : NDArray
        Lagrange multipliers.
    """

    VALID_TYPES = ['rev', 'tran', 'rev-rev', 'rev-tran', 'rigid', 'disc', 'rel-rot', 'rel-tran']

    def __init__(self, type='rev', iBody=None, jBody=None,
                 iPoint=None, jPoint=None,
                 iUvec=None, jUvec=None,
                 iFunct=None, L=0, R=1, x0=0, d0=None, fix=0):
        super().__init__()  # call to the Base class constructor
        if type not in Joint.VALID_TYPES:
            raise ValueError(f"Joint {self.COUNT}: unknown type '{type}', valid types: {Joint.VALID_TYPES}")
        self.type = type
        self.iBody = iBody if iBody is not None else Ground
        self.jBody = jBody if jBody is not None else Ground
        self.iPoint = iPoint
        self.jPoint = jPoint
        self.iUvec = iUvec
        self.jUvec = jUvec
        self.iFunct = iFunct
        self.L = L
        self.R = R
        self.x0 = x0
        self.d0 = d0 if d0 is not None else []
        self.fix = fix
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

    def __init__(self, type='a', t_start=0, f_start=0, t_end=1, f_end=1, dfdt_end=1, ncoeff=4, coeff=None):
        super().__init__()  # call to the Base class constructor
        self.type = type
        self.t_start = t_start
        self.f_start = f_start
        self.t_end = t_end
        self.f_end = f_end
        self.dfdt_end = dfdt_end
        self.ncoeff = ncoeff
        # Always allocate a 9-element numpy array so functData (types a, b, c)
        # can set any index [0..8] without IndexError.
        if coeff is not None:
            c = np.asarray(coeff, dtype=float).flatten()
            padded = np.zeros(9)
            padded[:min(len(c), 9)] = c[:min(len(c), 9)]
            self.coeff = padded
        else:
            self.coeff = np.zeros(9)