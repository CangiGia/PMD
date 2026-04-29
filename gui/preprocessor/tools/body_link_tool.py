"""BodyLinkTool — 2-click placement of a link body with end markers."""

from __future__ import annotations

import math

from PySide6.QtCore import QLineF, Qt
from PySide6.QtGui import QColor, QPen
from PySide6.QtWidgets import QGraphicsLineItem

from ..models import BodySpec, MarkerSpec, ShapeSpec
from .tool_base import Tool


class BodyLinkTool(Tool):
    """Two-click tool: pick endpoint A, then endpoint B.

    On commit, creates:

    * a ``BodySpec`` with shape ``"link"`` (length, thickness),
      positioned at the segment midpoint with orientation matching the
      segment direction;
    * two markers at the link endpoints (local coords ``±L/2, 0``).
    """

    name = "body_link"
    cursor = Qt.CrossCursor

    _THICKNESS = 0.02   # default link visual thickness (m)
    _MIN_LEN   = 0.02

    def __init__(self, window):
        super().__init__(window)
        self._a = None
        self._preview: QGraphicsLineItem | None = None

    def deactivate(self):
        super().deactivate()
        self._clear_preview()
        self._a = None

    # ─────────────────────────────────────────────────────────
    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton:
            return False
        p = self.scene.snap(event.scenePos())
        if self._a is None:
            self._a = p
            self._make_preview(p, p)
        else:
            self._commit_link(self._a, p)
            self._a = None
            self._clear_preview()
        return True

    def mouse_move(self, event) -> bool:
        if self._a is None or self._preview is None:
            return False
        b = self.scene.snap(event.scenePos())
        self._preview.setLine(QLineF(self._a[0], self._a[1], b[0], b[1]))
        return True

    # ─────────────────────────────────────────────────────────
    def _commit_link(self, a, b):
        ax, ay = a
        bx, by = b
        dx, dy = bx - ax, by - ay
        L = math.hypot(dx, dy)
        if L < self._MIN_LEN:
            return
        cx, cy = 0.5 * (ax + bx), 0.5 * (ay + by)
        phi = math.atan2(dy, dx)

        body = BodySpec(
            name=f"link{len(self.spec.bodies) + 1}",
            mass=1.0, inertia=L * L / 12.0,
            position=(cx, cy), orientation=phi,
            shape=ShapeSpec(kind="link",
                            params={"length": L, "thickness": self._THICKNESS}),
        )
        self.spec.bodies.append(body)

        m_a = MarkerSpec(name=f"{body.name}.A", body_id=body.id,
                         local_position=(-L / 2, 0.0), theta=0.0)
        m_b = MarkerSpec(name=f"{body.name}.B", body_id=body.id,
                         local_position=(+L / 2, 0.0), theta=0.0)
        self.spec.markers.extend([m_a, m_b])

        self.window.add_body_item(body)
        self.window.add_marker_item(m_a)
        self.window.add_marker_item(m_b)
        self._commit()
        # One-shot: revert to Select after a successful link creation.
        self.window.set_active_tool("select")

    # ─────────────────────────────────────────────────────────
    def _make_preview(self, p0, p1):
        self._preview = QGraphicsLineItem(QLineF(p0[0], p0[1], p1[0], p1[1]))
        pen = QPen(QColor("#3f8cff"))
        pen.setCosmetic(True)
        pen.setWidthF(2.0)
        pen.setStyle(Qt.DashLine)
        self._preview.setPen(pen)
        self._preview.setZValue(1000)
        self.scene.addItem(self._preview)

    def _clear_preview(self):
        if self._preview is not None:
            self.scene.removeItem(self._preview)
            self._preview = None
