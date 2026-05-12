# Unit System

`UnitSystem` declares the physical units used when defining a PMD model.
The solver always works in SI internally; this object drives display conversions
and plot axis labels.

## Creating a unit system

```python
from pmd.core import UnitSystem

us = UnitSystem(length='mm', force='N', angle='deg')
```

Pass it to `PlanarMultibodyModel`:

```python
model = PlanarMultibodyModel(bodies=[...], joints=[...], units=us)
```

## Supported units

| Dimension | Choices | Default |
|-----------|---------|---------|
| `length` | `"m"`, `"mm"`, `"cm"`, `"in"`, `"ft"` | `"m"` |
| `force` | `"N"`, `"kN"`, `"lbf"` | `"N"` |
| `angle` | `"rad"`, `"deg"` | `"rad"` |

## Conversion helpers

```python
us = UnitSystem(length='mm', force='N', angle='deg')

factor = us.to_si('length')   # 0.001  (mm → m)
label  = us.latex_unit('length')  # r'\mathrm{mm}'
moment = us.moment_unit        # 'N·mm'
```
