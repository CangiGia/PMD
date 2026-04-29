"""BodyRectTool — click-drag to create a rectangular rigid body."""

from __future__ import annotations

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsRectItem

from ..models import BodySpec, ShapeSpec
from .tool_base import Tool


class BodyRectTool(Tool):
    """Click on the canvas to drop a default rectangular body.

    Behaviour
    ---------
    * Press → record start point (snapped).
    * Drag  → live preview rectangle.
    * Release → commit:
        - if the drag is degenerate (≈ click), use the default
          size 200 × 100 mm centred on the press point;
        - otherwise the rectangle's bounding box defines the body.
    """

    name = "body_rect"
    cursor = Qt.CrossCursor

    _DEFAULT_W = 0.20
    _DEFAULT_H = 0.10
    _MIN_SIZE  = 0.01  # below this the drag is treated as a click

    def __init__(self, window):
        super().__init__(window)
        self._press_pt = None
        self._preview: QGraphicsRectItem | None = None

    def deactivate(self):
        super().deactivate()
        self._clear_preview()
        self._press_pt = None

    # ─────────────────────────────────────────────────────────
    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton:
            return False
        self._press_pt = self.scene.snap(event.scenePos())
        self._make_preview(self._press_pt, self._press_pt)
        return True

    def mouse_move(self, event) -> bool:
        if self._press_pt is None or self._preview is None:
            return False
        end = self.scene.snap(event.scenePos())
        self._update_preview(self._press_pt, end)
        return True

    def mouse_release(self, event) -> bool:
        if event.button() != Qt.LeftButton or self._press_pt is None:
            return False
        end = self.scene.snap(event.scenePos())
        self._clear_preview()

        x0, y0 = self._press_pt
        x1, y1 = end
        w = abs(x1 - x0)
        h = abs(y1 - y0)
        if w < self._MIN_SIZE or h < self._MIN_SIZE:
            cx, cy = x0, y0
            w, h = self._DEFAULT_W, self._DEFAULT_H
        else:
            cx = 0.5 * (x0 + x1)
            cy = 0.5 * (y0 + y1)

        spec = BodySpec(
            name=f"body{len(self.spec.bodies) + 1}",
            mass=1.0, inertia=0.01,
            position=(cx, cy),
            shape=ShapeSpec(kind="rectangle",
                            params={"width": w, "height": h}),
        )
        self.spec.bodies.append(spec)
        self.window.add_body_item(spec)
        self._press_pt = None
        self._commit()
        # One-shot: revert to Select so the body cannot be accidentally
        # moved / re-edited by stray clicks right after creation.
        self.window.set_active_tool("select")
        return True

    # ─────────────────────────────────────────────────────────
    def _make_preview(self, p0, p1):
        rect = self._rect_from_points(p0, p1)
        self._preview = QGraphicsRectItem(rect)
        pen = QPen(QColor("#3f8cff"))
        pen.setCosmetic(True)
        pen.setWidthF(1.0)
        pen.setStyle(Qt.DashLine)
        self._preview.setPen(pen)
        self._preview.setBrush(QBrush(QColor(63, 140, 255, 60)))
        self._preview.setZValue(1000)
        self.scene.addItem(self._preview)

    def _update_preview(self, p0, p1):
        if self._preview is not None:
            self._preview.setRect(self._rect_from_points(p0, p1))

    def _clear_preview(self):
        if self._preview is not None:
            self.scene.removeItem(self._preview)
            self._preview = None

    @staticmethod
    def _rect_from_points(p0, p1) -> QRectF:
        x0, y0 = p0
        x1, y1 = p1
        return QRectF(min(x0, x1), min(y0, y1),
                      abs(x1 - x0), abs(y1 - y0))
