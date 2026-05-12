"""pmd - Planar Multi-Body Dynamics Library."""

import sys
import os

sys.path.append(os.path.dirname(__file__))

from .utils import *
from .mechanics import *
from .model import Base, Ground, _GroundType, Body, Marker
from .shapes import Rectangle, Circle, Polygon
from .constraints import (
    Joint,
    RevJoint,
    TranJoint,
    RevRevJoint,
    RevTranJoint,
    RigidJoint,
    DiscJoint,
    RelRotJoint,
    RelTranJoint,
    Force,
    Weight,
    PtpForce,
    RotSdaForce,
    LocalForce,
    GlobalForce,
    Torque,
    UserForce,
    Function,
)
from .builder import *
from .solver import PlanarMultibodyModel
from .units import UnitSystem

__all__ = [
    # Model primitives
    "Body",
    "Ground",
    "Marker",
    # Geometric shapes (for Body.shape)
    "Rectangle",
    "Circle",
    "Polygon",
    # Joints
    "Joint",
    "RevJoint",
    "TranJoint",
    "RevRevJoint",
    "RevTranJoint",
    "RigidJoint",
    "DiscJoint",
    "RelRotJoint",
    "RelTranJoint",
    # Forces
    "Force",
    "Weight",
    "PtpForce",
    "RotSdaForce",
    "LocalForce",
    "GlobalForce",
    "Torque",
    "UserForce",
    # Driver function
    "Function",
    # Solver
    "PlanarMultibodyModel",
    # Units
    "UnitSystem",
]