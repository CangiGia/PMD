# PMD — Planar Multi-Body Dynamics

**PMD** is an open-source Python library for modelling, simulating, and post-processing planar rigid-body mechanical systems.

---

## What PMD does

PMD lets you define a mechanical system as a set of rigid bodies connected by kinematic constraints and force elements, then solve its motion numerically and inspect the results interactively.

### Model definition
- **Rigid bodies** — each with mass, inertia, and body-fixed reference markers
- **Joints** — revolute, translational, rigid, disc-contact, relative-rotation/translation, and higher-pair composite joints (RevRev, RevTran)
- **Force elements** — gravity, point-to-point spring-damper, rotational spring-damper-actuator, local/global applied forces, torques, and fully user-defined force laws
- **Driven constraints** — analytical functions (sine, cosine, polynomial, …) to prescribe any coordinate over time
- **Unit system** — SI or mixed-unit models (m/mm/cm/in/ft, N/kN/lbf, rad/deg) with automatic conversion

### Analysis types
| Mode | What it computes |
|---|---|
| **Dynamic** | Full equations of motion — positions, velocities, accelerations, reaction forces over time |
| **Kinematic** | Position and velocity analysis only, for fully driven systems |

### Solvers
- **LSODA** (default) — adaptive-step ODE integration via SciPy
- **Radau / RK45 / DOP853** — any SciPy `solve_ivp` method
- **CASADI-DAE** — implicit DAE collocation via CasADi (optional dependency), suited for stiff and highly constrained systems
- Automatic initial-condition correction (Newton–Raphson)
- Optional Baumgarte stabilisation and GGL regularisation

### Interactive GUI (optional, requires PySide6)
- **PreProcessor** — drag-and-drop canvas to build the model graphically: create bodies, place markers, add joints and forces, set initial conditions, launch the solver
- **PostProcessor** — plot time histories of positions, velocities, accelerations and reaction forces across multiple simulation sessions; animated playback with zoom and export

---

## Installation

```bash
# core only
pip install -e .

# with GUI support
pip install -e ".[gui]"

# with GUI + dev tools
pip install -e ".[gui,dev]"
```

Python ≥ 3.10 required.

---

## Quick example

```python
import pmd
from pmd.core import Body, Marker, RevJoint, Weight, PlanarMultibodyModel, UnitSystem

us = UnitSystem(length="m", force="N")

link = Body(name="link", mass=1.0, inertia=0.01)
A    = Marker(body=link, local_position=[0.0, 0.0])
B    = Marker(body=link, local_position=[0.5, 0.0])

pin  = RevJoint(iMarker=pmd.core.Ground.origin, jMarker=A)
grav = Weight(g=9.81, unit_system=us)

model = PlanarMultibodyModel(Bodies=[link], Joints=[pin], Forces=[grav])
T, uT = model.solve(analysis="dynamic", t_final=2.0, dt=0.01)
```

---

## Project structure

```
src/pmd/
    core/   ← model, solver, constraints, mechanics, units
    gui/
        preprocessor/   ← interactive model builder
        postprocessor/  ← results viewer and animation
```

---

## License

[MIT](LICENSE)
