"""Functional smoke test for preprocessor tools (offscreen)."""
import os
os.environ["QT_QPA_PLATFORM"] = "offscreen"

from PySide6.QtCore import QPointF, Qt
from PySide6.QtWidgets import QApplication

from PMD.gui.style import apply_light_theme
from PMD.gui.preprocessor.main_window import PreProcessorWindow
from PMD.gui.preprocessor.models import ModelSpec


class FakeEv:
    def __init__(self, x, y, btn=Qt.LeftButton):
        self._p = QPointF(x, y)
        self._b = btn

    def scenePos(self): return self._p
    def button(self): return self._b
    def accept(self): pass


def main():
    app = QApplication([])
    apply_light_theme(app)
    w = PreProcessorWindow(ModelSpec())
    w.show()
    app.processEvents()

    w.set_active_tool("body_rect")
    t = w._scene.active_tool
    t.mouse_press(FakeEv(0.0, 0.0))   # anchor centre
    t.mouse_press(FakeEv(0.1, 0.0))   # second click → orientation
    print("after rect:  bodies =", len(w.spec.bodies))
    assert w.spec.bodies[-1].locked is True

    w.set_active_tool("body_link")
    t = w._scene.active_tool
    t.mouse_press(FakeEv(0.5, 0.5))   # anchor A
    t.mouse_press(FakeEv(0.8, 0.5))   # orientation
    print("after link:  bodies =", len(w.spec.bodies),
          "markers =", len(w.spec.markers))
    assert w.spec.bodies[-1].locked is True

    w.set_active_tool("marker")
    t = w._scene.active_tool
    t.mouse_press(FakeEv(-0.4, -0.2))
    print("after free marker: markers =", len(w.spec.markers))

    w.set_active_tool("joint_rev")
    t = w._scene.active_tool
    m_items = list(w._marker_items.values())
    print("marker items in scene:", len(m_items))
    if len(m_items) >= 2:
        a, b = m_items[0], m_items[1]
        t.mouse_press(FakeEv(a.scenePos().x(), a.scenePos().y()))
        t.mouse_press(FakeEv(b.scenePos().x(), b.scenePos().y()))
    print("after joint: joints =", len(w.spec.joints))

    # PtP spring between same two markers
    w.set_active_tool("force_ptp")
    t = w._scene.active_tool
    if len(m_items) >= 2:
        a, b = m_items[0], m_items[1]
        t.mouse_press(FakeEv(a.scenePos().x(), a.scenePos().y()))
        t.mouse_press(FakeEv(b.scenePos().x(), b.scenePos().y()))
    print("after ptp:   forces =", len(w.spec.forces))
    assert any(f.kind == "PtpForce" for f in w.spec.forces)

    # Gravity toggle ON
    w._on_action("force_grav")
    assert any(f.kind == "Weight" for f in w.spec.forces)
    # Gravity toggle OFF
    w._on_action("force_grav")
    assert not any(f.kind == "Weight" for f in w.spec.forces)
    print("gravity toggle OK")
    print("OK")


if __name__ == "__main__":
    main()
