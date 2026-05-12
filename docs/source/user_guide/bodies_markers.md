# Bodies and Markers

## Body

A `Body` represents a rigid body in the planar simulation.  It carries:

- **mass** and **inertia** — scalar physical properties.
- **position** `[x, y]` and **orientation** `phi` — initial configuration.
- **velocity** `[dx, dy]` and **angular_velocity** `dphi` — initial velocities.
- An optional **shape** (`Rectangle`, `Circle`, `Polygon`) for visualisation.

```python
from pmd.core import Body

link = Body(
    mass=2.689,
    inertia=1.5e-3,
    position=[0.05, 0.0],
    orientation=0.0,
)
```

Markers are attached via `Body.add_marker()`:

```python
m1 = link.add_marker([0.0, 0.0], name="left_pin")
m2 = link.add_marker([0.1, 0.0], name="right_pin")
```

## Ground

`Ground` is the singleton inertial reference frame.  It is always available at
module level and does not need to be created:

```python
from pmd.core import Ground

g1 = Ground.add_marker([0.0, 0.0], theta=0.0, name="origin_pin")
```

`Ground` is *falsy* — `if body` evaluates to `False` for `Ground` and `True`
for any `Body`, which simplifies conditional logic inside joints and forces.

## Marker

A `Marker` is a body-fixed reference frame (a point, and optionally an
orientation).  Joints and forces are attached between pairs of markers.

```python
from pmd.core import Marker

# Equivalent to Body.add_marker — rarely needed directly
m = Marker(body=link, local_position=[0.05, 0.0], theta=0.0, name="mid")
```

`theta` is required for joints that use an orientation axis (e.g. `TranJoint`).
