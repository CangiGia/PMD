# Solver

`PlanarMultibodyModel` is the main simulation object.  Once constructed it
calls `_initialize()` automatically.

## Constructor

```python
from pmd.core import PlanarMultibodyModel

model = PlanarMultibodyModel(
    bodies=[body1, body2, ...],
    joints=[j1, j2, ...],
    forces=[f1, f2, ...],
    functions=[drv],          # optional
    units=UnitSystem(),       # optional
)
```

## Solving

```python
result = model.solve(
    analysis="dynamic",     # "dynamic" | "kinematic" | "static"
    method="LSODA",         # "LSODA" | "CASADI-DAE"
    t_final=10.0,
    dt=0.001,
    ic_correct=False,
    baumgarte_alpha=5.0,
    baumgarte_beta=5.0,
)
```

The returned `result` dict contains at minimum `T` (time array) and `U`
(state array, shape `(nSteps, 2·3·nB)`).

## Solver methods

| `method` | Description |
|----------|-------------|
| `"LSODA"` | SciPy `solve_ivp` with LSODA integrator; variable-order Adams/BDF. |
| `"CASADI-DAE"` | CasADi DAE solver (Radau collocation); symbolically-built system. Recommended for stiff problems. |

## Baumgarte stabilisation

For `"LSODA"`, constraint drift is suppressed by Baumgarte stabilisation:

$$\ddot{\Phi} + 2\alpha\dot{\Phi} + \beta^2\Phi = 0$$

Increase `baumgarte_alpha` and `baumgarte_beta` (default 5.0) for tighter
constraint enforcement at the cost of stiffness.

## Index convention (internal)

The following naming convention is used throughout the solver internals:

- **Body indices**: 0 = ground (fixed), 1..nB = moving bodies.
- **State vector** `u`:
  - `u[0 : 3·nB]` — positions `[r1x, r1y, φ1, r2x, r2y, φ2, …]`
  - `u[3·nB : 6·nB]` — velocities `[ṙ1x, ṙ1y, φ̇1, …]`
- **Internal attributes** (0-based):
  - `body._index_position = 3·Bi` — start index in position slice.
  - `body._index_velocity = 3·nB + 3·Bi` — start index in velocity slice.
  - `joint._rows / _rowe` — constraint row range `[rows, rowe)`.
  - `joint._colis / _colie`, `_coljs / _colje` — Jacobian column ranges.
