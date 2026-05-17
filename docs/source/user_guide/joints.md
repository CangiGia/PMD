# Joints

Joints impose kinematic constraints between two bodies (or between a body and
`Ground`).  They are passed to `PlanarMultibodyModel` via the `joints=` list.

## Revolute Joint — `RevJoint`

Constrains two marker positions to coincide (2 DOF removed).  Optional `fix=1`
also locks relative rotation (3 DOF removed → rigid connection).

```python
from pmd.core import RevJoint

j = RevJoint(iMarker=m_crank, jMarker=m_rod)
```

| Parameter | Description |
|-----------|-------------|
| `iMarker`, `jMarker` | Markers on the two bodies |
| `fix` | `0` (default) = free rotation; `1` = fixed rotation |
| `q0` | Initial relative angle (rad) used by the assembler |

## Translational Joint — `TranJoint`

Constrains motion to one axis (2 DOF removed).  Requires `iMarker` to have
`theta` set (axis direction).

```python
j = TranJoint(iMarker=slider_top, jMarker=guide_top)
```

## Rev-Rev Joint — `RevRevJoint`

Constant-distance link between two revolute points (1 DOF removed).

```python
from pmd.core import RevRevJoint

j = RevRevJoint(iMarker=m1, jMarker=m2, L=0.12)
```

## Rev-Tran Joint — `RevTranJoint`

Revolute-translational composite joint (1 DOF removed).

## Rigid Joint — `RigidJoint`

Welds two bodies together (3 DOF removed).

## Disc Joint — `DiscJoint`

Rolling-without-slipping disc on a flat surface (2 DOF removed).

```python
from pmd.core import DiscJoint

j = DiscJoint(iBody=wheel, R=0.3, x0=0.0)
```

## Rotational Motion — `RotMotion`

Prescribed relative rotation between two bodies, driven by a `Function`.

```python
from pmd.core import RotMotion, Function

f   = Function(type='a', coeff=[0.0, 10.0 * 3.14159/180, 0.0])  # 10 deg/s ramp
drv = RotMotion(iBody=crank, jBody=Ground, iFunct=f)
```

## Translational Motion — `TranMotion`

Prescribed relative translation between two markers, driven by a `Function`.
