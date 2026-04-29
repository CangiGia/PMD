"""JointTool — Adams-style 2-click joint creation between markers."""

from __future__ import annotations

from typing import ClassVar

from PySide6.QtCore import Qt

from ..models import JointSpec
from ..widgets import MarkerItem
from .tool_base import Tool


class JointTool(Tool):
    """Click marker_i → click marker_j → joint of the configured kind."""

    name = "joint"
    cursor = Qt.PointingHandCursor

    KIND: ClassVar[str] = "RevJoint"

    def __init__(self, window):
        super().__init__(window)
        self._first: MarkerItem | None = None

    def activate(self):
        super().activate()
        self.window.statusBar().showMessage(
            f"{self.KIND}: pick first marker (Esc to cancel)…", 0)

    def deactivate(self):
        super().deactivate()
        if self._first is not None:
            self._first.set_highlighted(False)
        self._first = None
        self.window.statusBar().clearMessage()

    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton:
            return False
        marker = self._marker_under(event.scenePos())
        if marker is None:
            # Give the user feedback so they know the click was seen
            # but no marker was close enough.
            which = "first" if self._first is None else "second"
            self.window.statusBar().showMessage(
                f"{self.KIND}: no marker near the cursor \u2014 click on a "
                f"marker to pick the {which} one (Esc to cancel)\u2026",
                0)
            return True  # consume click but do nothing

        if self._first is None:
            self._first = marker
            marker.set_highlighted(True)
            self.window.statusBar().showMessage(
                f"{self.KIND}: pick second marker (Esc to cancel)…", 0)
            return True

        if marker is self._first:
            return True  # ignore double-click on same marker

        joint = JointSpec(
            name=f"jnt{len(self.spec.joints) + 1}",
            kind=self.KIND,  # type: ignore[arg-type]
            i_marker_id=self._first.spec.id,
            j_marker_id=marker.spec.id,
            params={},
        )
        self.spec.joints.append(joint)
        self.window.add_joint_visual(joint)
        self._first.set_highlighted(False)
        self._first = None
        self.window.statusBar().clearMessage()
        self._commit()
        # One-shot: revert to Select after the joint is created.
        self.window.set_active_tool("select")
        return True

    # ─────────────────────────────────────────────────────────
    # Pick radius around a marker (pixels). A click anywhere inside
    # this radius selects the closest marker — so the user doesn't
    # have to land exactly on the small origin dot.
    _PICK_PX: ClassVar[float] = 16.0

    def _marker_under(self, scene_pt) -> MarkerItem | None:
        # 1. Exact hit on the marker's pick shape — always wins.
        for it in self.scene.items(scene_pt):
            if isinstance(it, MarkerItem):
                return it
        # 2. Otherwise, fall back to nearest-marker within tolerance
        #    so the user doesn't need pixel-perfect aim.
        view = self.scene.views()[0] if self.scene.views() else None
        ppu = abs(view.transform().m11()) if view else 1.0
        tol2 = (self._PICK_PX / max(ppu, 1e-9)) ** 2
        px, py = float(scene_pt.x()), float(scene_pt.y())
        best: MarkerItem | None = None
        best_d2 = tol2
        for it in self.scene.items():
            if not isinstance(it, MarkerItem):
                continue
            sp = it.scenePos()
            dx = sp.x() - px
            dy = sp.y() - py
            d2 = dx * dx + dy * dy
            if d2 <= best_d2:
                best_d2 = d2
                best = it
        return best


class RevJointTool(JointTool):
    name = "joint_rev"
    KIND = "RevJoint"


class TranJointTool(JointTool):
    name = "joint_tran"
    KIND = "TranJoint"
