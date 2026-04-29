"""Smoke test for the new pre-processor capabilities:
delete, rename cascade, undo/redo, auto CM marker, material/density.

Run headless with::

    QT_QPA_PLATFORM=offscreen python examples/_smoke_preprocessor_actions.py
"""

from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent))

from PySide6.QtCore import QPointF, Qt
from PySide6.QtWidgets import QApplication

from PMD.gui.preprocessor.main_window import PreProcessorWindow
from PMD.gui.preprocessor.models import ModelSpec, MATERIALS


class FakeEv:
    def __init__(self, x, y, btn=Qt.LeftButton):
        self._p = QPointF(x, y)
        self._b = btn

    def scenePos(self): return self._p
    def button(self):   return self._b
    def accept(self):   pass


def main() -> int:
    app = QApplication.instance() or QApplication([])
    w = PreProcessorWindow(ModelSpec())
    w.show()

    # ── 1. Create a link → CM marker auto-created ──────────────
    w.set_active_tool("body_link")
    t = w._scene.active_tool
    t.mouse_press(FakeEv(0.0, 0.0))   # anchor A
    t.mouse_press(FakeEv(0.4, 0.0))   # commit
    body = w.spec.bodies[-1]
    cm = next(m for m in w.spec.markers if m.body_id == body.id
              and m.name.endswith(".cm"))
    assert cm.local_position == (0.0, 0.0)
    print("CM marker auto-created:", cm.name)

    # ── 2. Auto mass / inertia from steel density ──────────────
    L = body.shape.params["length"]; thk = body.shape.params["thickness"]
    expected = MATERIALS["Steel"] * L * thk * 0.01
    assert math.isclose(body.mass, expected, rel_tol=1e-9), (
        body.mass, expected)
    print(f"auto mass = {body.mass:.4f} kg  (steel)")

    # ── 3. Rename cascades to child markers ────────────────────
    w._inspector.show_item("body", body.id)
    ed = w._inspector._editor
    ed.name.setText("crank")
    ed._on_name()
    names = sorted(m.name for m in w.spec.markers if m.body_id == body.id)
    assert names == ["crank.A", "crank.B", "crank.cm"], names
    print("rename cascade OK:", names)

    # ── 4. Delete one marker via the tree pathway ──────────────
    mark_a = next(m for m in w.spec.markers if m.name == "crank.A")
    w._on_delete_requested("marker", mark_a.id)
    assert all(m.name != "crank.A" for m in w.spec.markers)
    print("marker delete OK")

    # ── 5. Undo restores the deleted marker ────────────────────
    w._undo()
    assert any(m.name == "crank.A" for m in w.spec.markers)
    print("undo OK")
    w._redo()
    assert all(m.name != "crank.A" for m in w.spec.markers)
    print("redo OK")

    # ── 6. Delete the body cascades to all its markers ─────────
    w._on_delete_requested("body", body.id)
    assert not w.spec.bodies
    assert not [m for m in w.spec.markers if m.body_id == body.id]
    print("body delete cascades to markers OK")

    print("ALL OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
