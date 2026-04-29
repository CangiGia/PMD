"""Smoke test for the typed Inspector editors (offscreen)."""

from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from PMD.gui.preprocessor.models import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    ShapeSpec,
)
from PMD.gui.preprocessor.panels.inspector_panel import InspectorPanel


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)

    spec = ModelSpec(name="probe")
    body = BodySpec(
        id="body_1", name="link_a",
        mass=2.0, inertia=0.05,
        position=(0.1, 0.2), orientation=0.0,
        shape=ShapeSpec(kind="rectangle",
                        params={"width": 0.4, "height": 0.1}),
    )
    spec.bodies.append(body)
    marker = MarkerSpec(id="mk_1", name="m_a", body_id=body.id,
                        local_position=(0.0, 0.0), theta=0.0)
    spec.markers.append(marker)
    g_marker = MarkerSpec(id="mk_g", name="m_g", body_id="ground",
                          local_position=(0.0, 0.0), theta=0.0)
    spec.markers.append(g_marker)
    joint = JointSpec(id="jt_1", name="j_a", kind="RevJoint",
                      i_marker_id=marker.id, j_marker_id=g_marker.id)
    spec.joints.append(joint)
    force = ForceSpec(id="fr_1", name="grav",
                      kind="Weight", i_marker_id=marker.id,
                      params={"gx": 0.0, "gy": -9.81})
    spec.forces.append(force)

    insp = InspectorPanel()
    insp.set_spec(spec)

    # ── Body ──────────────────────────────────────────────
    edits: list[tuple[str, str]] = []
    insp.spec_changed.connect(lambda k, i: edits.append((k, i)))

    insp.show_item("body", body.id)
    ed = insp._editor
    assert ed is not None and ed.spec_kind == "body"
    # Enable manual override so mass / inertia become editable.
    ed.override.setChecked(True)
    ed.mass.setValue(3.5)
    ed.x.setValue(1.0)
    ed.phi_deg.setValue(90.0)
    ed._shape_widgets["width"].setValue(0.6)
    assert math.isclose(body.mass, 3.5)
    assert body.mass_override is True
    assert body.position[0] == 1.0
    assert math.isclose(body.orientation, math.pi / 2)
    assert body.shape.params["width"] == 0.6

    # Switch back to auto: mass must come from density × area × depth.
    ed.override.setChecked(False)
    ed.material.setCurrentText("Steel")
    expected_mass = 7850.0 * (0.6 * body.shape.params["height"]) * 0.01
    assert math.isclose(body.mass, expected_mass, rel_tol=1e-9)
    print("body edits OK", edits[-1])

    # ── Marker ────────────────────────────────────────────
    insp.show_item("marker", marker.id)
    ed = insp._editor
    ed.xi.setValue(0.05)
    ed.theta_deg.setValue(45.0)
    assert marker.local_position[0] == 0.05
    assert math.isclose(marker.theta, math.pi / 4)
    print("marker edits OK")

    # ── Joint ─────────────────────────────────────────────
    insp.show_item("joint", joint.id)
    ed = insp._editor
    ed.kind.setCurrentText("TranJoint")
    assert joint.kind == "TranJoint"
    print("joint edits OK")

    # ── Force ─────────────────────────────────────────────
    insp.show_item("force", force.id)
    ed = insp._editor
    print("force editor =", ed, "kind =", getattr(ed, 'spec_kind', None))
    assert ed is not None and ed.spec_kind == "force"
    print("force editor mounted OK")

    # ── Clear ─────────────────────────────────────────────
    insp.clear()
    assert insp._editor is None
    print("clear OK")

    print("ALL OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
