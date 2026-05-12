# Driver Functions

`Function` objects prescribe the motion of driver joints (`RelRotJoint`,
`RelTranJoint`).  They encode the time-history of a kinematic input as a
polynomial or piecewise-polynomial expression.

## Function types

| `type` | Description |
|--------|-------------|
| `'a'` | Polynomial: `f(t) = coeff[0] + coeff[1]*t + coeff[2]*t² + …` |
| `'b'` | Cubic spline between `(t_start, f_start)` and `(t_end, f_end)` |
| `'c'` | Modified-sine ramp with bounded acceleration |

## Type `'a'` — polynomial

```python
from pmd.core import Function

# Constant angular position of 55° (converted to radians)
import math
f = Function(type='a', coeff=[55.1821 * math.pi/180, 0.0])

# Linear ramp: 20 deg/s
f_ramp = Function(type='a', coeff=[0.0, 20.0 * math.pi/180])
```

The `coeff` array is always padded to 9 elements internally; unused
coefficients default to zero.

## Type `'b'` — cubic spline

```python
f = Function(type='b',
             t_start=0.0, f_start=0.0,
             t_end=2.0,   f_end=math.pi)
```

The function smoothly connects the two endpoints with zero velocity at both
ends.

## Using a Function with a driver joint

```python
from pmd.core import RelRotJoint

drv = RelRotJoint(iBody=crank, jBody=Ground, iFunct=f)
```

The solver evaluates `f(t)`, `f'(t)`, and `f''(t)` at each time step to
compute the constraint right-hand side.
