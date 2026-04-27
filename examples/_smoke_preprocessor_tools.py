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
    t.mouse_press(FakeEv(0.0, 0.0))
    t.mouse_release(FakeEv(0.0, 0.0))
    print("after rect:  bodies =", len(w.spec.bodies))

    w.set_active_tool("body_link")
    t = w._scene.active_tool
    t.mouse_press(FakeEv(0.5, 0.5))
    t.mouse_press(FakeEv(0.8, 0.5))
    print("after link:  bodies =", len(w.spec.bodies),
          "markers =", len(w.spec.markers))

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
    print("OK")


if __name__ == "__main__":
    main()
