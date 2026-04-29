"""Smoke test: build a tiny pendulum spec and solve dynamically."""

from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np

from PMD.gui.preprocessor.builder import build_model
from PMD.gui.preprocessor.models import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    ShapeSpec,
)


def _make_pendulum() -> ModelSpec:
    spec = ModelSpec(name="pendulum")
    L = 0.4
    bar = BodySpec(
        id="bar", name="bar",
        mass=1.0,
        inertia=L * L / 12.0,
        position=(L / 2.0, 0.0),
        orientation=0.0,
        shape=ShapeSpec(kind="link",
                        params={"length": L, "thickness": 0.02}),
    )
    spec.bodies.append(bar)
    g = MarkerSpec(id="m_g", name="m_g", body_id="ground",
                   local_position=(0.0, 0.0))
    p = MarkerSpec(id="m_p", name="m_p", body_id="bar",
                   local_position=(-L / 2.0, 0.0))
    spec.markers.extend([g, p])
    spec.joints.append(JointSpec(id="j_r", name="hinge",
                                 kind="RevJoint",
                                 i_marker_id="m_g", j_marker_id="m_p"))
    spec.forces.append(ForceSpec(id="grav", name="g",
                                 kind="Weight",
                                 params={"gravity": 9.80665}))
    return spec


def main() -> int:
    spec = _make_pendulum()
    model = build_model(spec)
    print("model built: bodies =", len(model.Bodies),
          "joints =", len(model.Joints),
          "forces =", len(model.Forces))

    T, uT = model.solve(
        analysis="dynamic",
        method="CASADI-DAE",
        t_span=(0.0, 0.5),
        t_eval=np.linspace(0.0, 0.5, 200),
        ic_correct=True,
    )
    assert T.shape[0] == 200
    assert uT.shape[0] == 200
    print(f"solve OK: {len(T)} steps, last t = {T[-1]:.3f}")
    print("ALL OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
