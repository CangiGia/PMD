"""
Planar Multi-Body Dynamics Constraints and Forces

This module provides Force, Joint, and Function classes for defining
constraints and external forces in planar multi-body dynamic systems.

Author: Giacomo Cangi
"""


import numpy as np
from numpy.typing import *
from .utils import *
from .model import Base, Ground


class Force(Base):
    """
    Create a force in a planar multi-body dynamic simulation.

    Attributes
    ----------
    type : str
        Element type: ``'ptp'``, ``'rot-sda'``, ``'weight'``,
        ``'flocal'``, ``'f'``, ``'T'``, ``'user'``.
    iMarker : Marker or None
        Head (arrow) marker.
    jMarker : Marker or None
        Tail marker.
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

    def __init__(self, type='ptp', iMarker=None, jMarker=None,
                 iBody=None, jBody=None,
                 k=0, L0=0, theta0=0, dc=0, f_a=0, T_a=0,
                 flocal=None, f=None, T=0, callback=None):
        super().__init__()  # call to the Base class constructor
        if type not in Force.VALID_TYPES:
            raise ValueError(f"Force {self.COUNT}: unknown type '{type}', valid types: {Force.VALID_TYPES}")
        self.type = type
        self.iMarker = iMarker
        self.jMarker = jMarker
        self.iBody = iMarker.body if iMarker is not None else (iBody or Ground)
        self.jBody = jMarker.body if jMarker is not None else (jBody or Ground)
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
    iMarker : Marker or None
        Marker on body i.
    jMarker : Marker or None
        Marker on body j.
    iBody : Body or Ground
        Body ``i``.
    jBody : Body or Ground
        Body ``j``.
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
    q0 : float
        Initial DOF value (used by assembler).

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

    def __init__(self, type='rev',
                 iMarker=None, jMarker=None,
                 iBody=None, jBody=None,
                 iFunct=None, L=0, R=1, x0=0, d0=None, fix=0,
                 q0=0):
        super().__init__()  # call to the Base class constructor
        if type not in Joint.VALID_TYPES:
            raise ValueError(f"Joint {self.COUNT}: unknown type '{type}', valid types: {Joint.VALID_TYPES}")

        self.type = type
        self.iMarker = iMarker
        self.jMarker = jMarker

        # For marker-based joints, derive iBody/jBody from markers
        if iMarker is not None:
            self.iBody = iMarker.body
        else:
            self.iBody = iBody if iBody is not None else Ground

        if jMarker is not None:
            self.jBody = jMarker.body
        else:
            self.jBody = jBody if jBody is not None else Ground

        # Validation V5: tran and rev-tran require orientation on both markers
        if type in ('tran',) and iMarker is not None and jMarker is not None:
            if not iMarker.has_orientation or not jMarker.has_orientation:
                raise ValueError(
                    f"Joint type '{type}' requires both iMarker and jMarker "
                    f"to have orientation (theta != None)")

        if type in ('rev-tran',) and iMarker is not None:
            if not iMarker.has_orientation:
                raise ValueError(
                    f"Joint type '{type}' requires iMarker "
                    f"to have orientation (theta != None)")

        # Validation V6: markers must be on different bodies
        if iMarker is not None and jMarker is not None:
            if iMarker.body is jMarker.body:
                raise ValueError(
                    f"iMarker and jMarker must be on different bodies")

        self.iFunct = iFunct
        self.L = L
        self.R = R
        self.x0 = x0
        self.d0 = d0 if d0 is not None else []
        self.fix = fix
        self.q0 = q0
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
