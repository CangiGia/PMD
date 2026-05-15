# Overview

PMD (*Planar Multi-Body Dynamics*) is an open-source Python library for
building and simulating planar rigid-body mechanical systems.

## Key concepts

| Concept | Description |
|---------|-------------|
| **Body** | A rigid body described by mass, inertia, and initial conditions. |
| **Marker** | A body-fixed reference point used to connect joints and forces. |
| **Joint** | A kinematic constraint between two bodies (or a body and ground). |
| **Force** | An applied load: gravity, spring-damper, torque, or user-defined. |
| **Function** | An analytical driver for prescribed motion constraints. |
| **PlanarMultibodyModel** | The model container and solver. |

## How it works

1. Create `Body` objects with mass, inertia, and initial conditions.
2. Add `Marker` objects to bodies to define attachment points.
3. Connect bodies with `Joint` constraints (revolute, translational, etc.).
4. Add `Force` elements (gravity, spring, etc.) and optionally `Function` drivers.
5. Pass everything to `PlanarMultibodyModel` and call `.solve()`.

## Supported analyses

- **Dynamic** — full time-integration of the equations of motion via
  CasADi DAE collocation (Radau IIA).
- **Kinematic** — position/velocity/acceleration analysis for fully-constrained
  systems (DOF = 0).
- **Static** — static equilibrium via Newton-Raphson iteration.

## Visualisation

Use the built-in `PostProcessor` to plot results, or the interactive
`PreProcessor` GUI to build models visually.
