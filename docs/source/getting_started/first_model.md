# Your First Model — Mass-Spring-Damper

This walkthrough builds a single-DOF mass-spring-damper system: a 1 000 kg
body constrained to translate vertically, connected to ground by a spring
(k = 1 MN/m, L₀ = 1.5 m) and a damper (c = 6 500 N·s/m).

## Complete example

```python
from pmd.core import (
    Body, Ground, Marker,
    TranJoint, PtpForce,
    PlanarMultibodyModel,
)

# 1. Create the body
mass = Body(
    mass=1000, inertia=0.5,
    position=[0, 0],
    velocity=[0, 10],       # initial upward velocity 10 m/s
)

# 2. Add markers
g_bottom = Ground.add_marker([0, -1.5], name="ground_bottom")
g_top    = Ground.add_marker([0,  0.0], name="ground_top")
m_bottom = mass.add_marker([0, -0.25], name="mass_bottom")
m_top    = mass.add_marker([0,  0.25], name="mass_top")

# 3. Add joints and forces
joint  = TranJoint(iMarker=m_top,    jMarker=g_top)
spring = PtpForce(iMarker=m_bottom,  jMarker=g_bottom, k=1e6, L0=1.5, dc=6500)

# 4. Build and solve
model = PlanarMultibodyModel(
    bodies=[mass],
    joints=[joint],
    forces=[spring],
)
result = model.solve(t_final=3.0, dt=0.001)
```

## What happens step by step

1. `Body` — defines the moving mass with its inertia and initial state.
2. `Ground.add_marker` / `mass.add_marker` — attach named reference points.
3. `TranJoint` — constrains the mass to slide along the vertical axis.
4. `PtpForce` — applies a spring-damper force between two markers.
5. `PlanarMultibodyModel` — assembles the system and calls `_initialize()`.
6. `solve()` — integrates the equations of motion and returns a results dict.

## Visualising results

```python
from pmd.gui import PostProcessor
pp = PostProcessor(model)
pp.show()
```

See {doc}`../user_guide/gui` for more details on the post-processor.
