# Forces

Force elements apply loads to bodies without imposing position constraints.
They are passed to `PlanarMultibodyModel` via the `forces=` list.

## Weight — gravity

`Weight` applies gravitational force to **all** bodies in the model.

```python
from pmd.core import Weight

g = Weight(gravity=9.81)          # standard gravity, downward
```

Default direction is `[0, -1]` (negative Y). Override with `gravity_direction`.

## Point-to-Point Force — `PtpForce`

Linear spring-damper-actuator between two markers:

$$F = k \cdot (L - L_0) + c \cdot \dot{L} + F_a$$

```python
from pmd.core import PtpForce

spring = PtpForce(iMarker=m1, jMarker=m2, k=1e6, L0=1.5, dc=6500)
```

| Parameter | Description |
|-----------|-------------|
| `k` | Spring stiffness (N/m) |
| `L0` | Natural length (m) |
| `dc` | Damping coefficient (N·s/m) |
| `f_a` | Constant actuator force (N) |

## Rotational SDA — `RotSdaForce`

Torsional spring-damper-actuator between two bodies:

$$T = k \cdot (\theta - \theta_0) + c \cdot \dot\theta + T_a$$

```python
from pmd.core import RotSdaForce

torsion = RotSdaForce(iBody=arm, jBody=Ground, k=500, theta0=0.0, dc=10)
```

## Local Force — `LocalForce`

Constant force vector expressed in the body-local frame.

## Global Force — `GlobalForce`

Constant force vector expressed in the global (inertial) frame.

## Torque

Constant scalar torque applied to a body.

## User-Defined Force — `UserForce`

Arbitrary force via a Python callback function:

```python
def my_force():
    return [{'body': link, 'force': [0, -500], 'torque': 0}]

uf = UserForce(callback=my_force)
```
