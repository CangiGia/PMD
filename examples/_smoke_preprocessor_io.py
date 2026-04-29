"""Smoke test for save/load round-trip and post-load scene rebuild."""

from __future__ import annotations

import os
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from PMD.gui.preprocessor.app import PreProcessor
from PMD.gui.preprocessor.models import (
    BodySpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    ShapeSpec,
    load_model,
    save_model,
)


def _make_spec() -> ModelSpec:
    spec = ModelSpec(name="probe")
    b = BodySpec(id="body_a", name="link",
                 mass=1.5, inertia=0.02,
                 position=(0.3, 0.1), orientation=0.5,
                 shape=ShapeSpec(kind="link",
                                 params={"length": 0.4, "thickness": 0.02}))
    spec.bodies.append(b)
    m1 = MarkerSpec(id="mk_g", name="g",  body_id="ground",
                    local_position=(0.0, 0.0))
    m2 = MarkerSpec(id="mk_a", name="a",  body_id="body_a",
                    local_position=(0.2, 0.0), theta=0.1)
    spec.markers.extend([m1, m2])
    spec.joints.append(JointSpec(id="jt_a", name="rev",
                                 kind="RevJoint",
                                 i_marker_id="mk_a", j_marker_id="mk_g"))
    return spec


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)

    spec = _make_spec()

    # ---------------- round-trip ----------------
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "test.pmdmodel.json")
        save_model(spec, p)
        loaded = load_model(p)

    assert loaded.name == spec.name
    assert len(loaded.bodies) == 1 and loaded.bodies[0].id == "body_a"
    assert isinstance(loaded.bodies[0].position, tuple)
    assert loaded.bodies[0].shape.kind == "link"
    assert loaded.bodies[0].shape.params["length"] == 0.4
    assert len(loaded.markers) == 2
    assert loaded.markers[1].body_id == "body_a"
    assert isinstance(loaded.markers[1].local_position, tuple)
    assert len(loaded.joints) == 1 and loaded.joints[0].i_marker_id == "mk_a"
    print("round-trip OK")

    # ---------------- in-window reload ----------------
    from PMD.gui.preprocessor.main_window import PreProcessorWindow
    win = PreProcessorWindow(ModelSpec())

    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "test.pmdmodel.json")
        save_model(spec, p)
        win.load_project(p)

    assert len(win._body_items)   == 1
    assert len(win._marker_items) == 2
    assert len(win._joint_items)  == 1
    assert win._project_path is not None
    print("scene rebuild OK")

    # New project clears state
    win._on_new_project()
    assert win._body_items   == {}
    assert win._marker_items == {}
    assert win._joint_items  == {}
    assert win._project_path is None
    print("new-project clear OK")

    print("ALL OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
