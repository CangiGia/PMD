# Analyses

`PlanarMultibodyModel.solve()` supports three analysis types selected by the
`analysis=` parameter.

## Dynamic analysis (default)

Full time-integration of the equations of motion:

$$M\ddot{q} + \Phi_q^T \lambda = Q, \quad \Phi(q, t) = 0$$

```python
result = model.solve(analysis="dynamic", method="LSODA", t_final=10.0, dt=0.01)
```

The result dict contains:

| Key | Shape | Description |
|-----|-------|-------------|
| `T` | `(nSteps,)` | Time array |
| `U` | `(nSteps, 6·nB)` | State vector (positions + velocities) |

Individual body trajectories can be extracted from `U` using
`body._index_position` and `body._index_velocity`.

## Kinematic analysis

For fully-constrained systems (DOF = 0), the acceleration-level equations can
be solved without integrating:

```python
result = model.solve(analysis="kinematic", t_final=10.0, dt=0.001)
```

Positions and velocities are propagated step-by-step from the initial
conditions using constraint equations only.

## Static analysis

Finds the static equilibrium position by iterating Newton-Raphson until the
constraint residual is below tolerance:

```python
result = model.solve(analysis="static")
```

The result has `T = [0.0]` and a single-step `U`.

## Initial condition correction

Setting `ic_correct=True` projects the initial conditions onto the constraint
manifold before time integration — useful when the user-supplied ICs are
slightly off:

```python
result = model.solve(ic_correct=True, ...)
```
