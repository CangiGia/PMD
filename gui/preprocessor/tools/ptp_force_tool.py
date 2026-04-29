"""PtpForceTool — point-to-point spring/damper between two markers."""

from __future__ import annotations

from PySide6.QtCore import Qt

from ..models import ForceSpec
from ..widgets import MarkerItem
from .tool_base import Tool


class PtpForceTool(Tool):
    """Click marker_i → click marker_j → PtpForce (k=0, L0=auto, dc=0)."""

    name = "force_ptp"
    cursor = Qt.PointingHandCursor

    def __init__(self, window):
        super().__init__(window)
        self._first: MarkerItem | None = None

    def deactivate(self):
        super().deactivate()
        if self._first is not None:
            self._first.set_highlighted(False)
        self._first = None

    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton:
            return False
        marker = self._marker_under(event.scenePos())
        if marker is None:
            return True

        if self._first is None:
            self._first = marker
            marker.set_highlighted(True)
            self.window.statusBar().showMessage(
                "PtpForce: pick second marker (Esc to cancel)\u2026", 0)
            return True

        if marker is self._first:
            return True

        # Use current world distance between the two markers as L0.
        a = self._first.scenePos()
        b = marker.scenePos()
        L0 = float(((a.x() - b.x()) ** 2 + (a.y() - b.y()) ** 2) ** 0.5)

        force = ForceSpec(
            name=f"ptp{len(self.spec.forces) + 1}",
            kind="PtpForce",
            i_marker_id=self._first.spec.id,
            j_marker_id=marker.spec.id,
            params={"k": 0.0, "L0": L0, "dc": 0.0, "f_a": 0.0},
        )
        self.spec.forces.append(force)
        self._first.set_highlighted(False)
        self._first = None
        self.window.statusBar().clearMessage()
        self._commit()
        return True

    # ─────────────────────────────────────────────────────────
    def _marker_under(self, scene_pt) -> MarkerItem | None:
        for it in self.scene.items(scene_pt):
            if isinstance(it, MarkerItem):
                return it
        return None
