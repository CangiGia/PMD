# PMD — Planar Multi-Body Dynamics

**PMD** is an open-source Python library for modelling, simulating, and post-processing planar rigid-body mechanical systems.

Given a set of rigid bodies connected by kinematic constraints and force elements, PMD assembles and solves the constrained equations of motion numerically, then lets you inspect positions, velocities, accelerations, and reaction forces over time — both programmatically and through an interactive GUI.

---

## 🏗️ Model definition

A PMD model is built from three building blocks:

**Bodies and markers**
- Rigid bodies with mass, moment of inertia, and arbitrary body-fixed reference markers
- A singleton ground body acts as the inertial reference frame

**Joints (kinematic constraints)**
- Revolute, translational, rigid, disc-contact
- Higher-pair composites: RevRev, RevTran
- Relative-rotation and relative-translation constraints
- Analytical driver functions (sine, cosine, polynomial, …) to prescribe any coordinate over time

**Force elements**
- Gravity (Weight)
- Point-to-point spring-damper (PtpForce)
- Rotational spring-damper-actuator (RotSdaForce)
- Applied forces in local or global frame, torques
- Fully user-defined force laws (UserForce)

---

## 📐 Analysis and solvers

| Analysis | What is computed |
|---|---|
| **Dynamic** | Full equations of motion — positions, velocities, accelerations, constraint reactions |
| **Kinematic** | Position and velocity solution for fully driven systems |

Available integration methods:

- **LSODA** (default) — adaptive-step ODE via SciPy, handles stiff/non-stiff switching automatically
- **Radau / RK45 / DOP853** — any SciPy `solve_ivp` backend
- **CasADi-DAE** — implicit DAE collocation for stiff or highly constrained systems *(optional dependency)*

All analyses include automatic Newton–Raphson initial-condition correction. Optional Baumgarte stabilisation and GGL regularisation are available for long-duration simulations.

A configurable `UnitSystem` (m / mm / cm / in / ft · N / kN / lbf · rad / deg) handles unit conversion transparently.

---

## 🖥️ Interactive GUI

PMD ships two optional graphical interfaces (requires PySide6):

**PreProcessor** — an interactive canvas to build the model visually:
- Create bodies and place markers by click or snap
- Add joints and forces from a ribbon toolbar
- Set initial conditions and launch the solver directly from the UI
- Save and reload models to/from JSON

**PostProcessor** — a results viewer for one or more simulation sessions:
- Time-history plots of positions, velocities, accelerations, and constraint reactions
- Animated playback of the mechanism motion
- Configurable display units, zoom, and figure export

---

## ⚡ Quick example

```python
from pmd.core import Body, Marker, RevJoint, Weight, PlanarMultibodyModel, UnitSystem
from pmd.core import Ground

us   = UnitSystem(length="m", force="N")
link = Body(name="pendulum", mass=1.0, inertia=0.01)

A   = Marker(body=link, local_position=[0.0, 0.0])
pin = RevJoint(iMarker=Ground.origin, jMarker=A)
g   = Weight(g=9.81, unit_system=us)

model    = PlanarMultibodyModel(Bodies=[link], Joints=[pin], Forces=[g])
T, uT    = model.solve(analysis="dynamic", t_final=2.0, dt=0.01)
```

After `solve`, results are accessible on each `Body` and `Joint` object and can be passed directly to the PostProcessor GUI.

---

## ❓ Support and questions

For questions, ideas, or general discussion please open a [Discussion](../../discussions).

For bugs or feature requests please open an [Issue](../../issues) describing the problem and, where possible, a minimal reproducible example.

---

## 🤝 Contributing

PMD is under active development. Contributions, bug reports, and suggestions are welcome. Please open an issue or pull request — all feedback helps.

---

## 📜 License

[MIT](LICENSE)
