"""BodyRectTool \u2014 inspector-driven 2-click placement of a rectangle.

Workflow mirrors :class:`BodyLinkTool`:

1. Pick the *Rect* tool. A draft :class:`BodySpec` is shown in the
   Inspector so the user can edit width / height / mass / inertia.
2. First click anchors the body's centre at the click point.
3. Mouse rotates the rectangle about its centre.
4. Second click commits the body, locks geometry and reverts to Select.
"""

from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsRectItem

from ..models import BodySpec, ShapeSpec, compute_mass_props
from .tool_base import Tool


class BodyRectTool(Tool):
    name = "body_rect"
    cursor = Qt.CrossCursor

    _DEF_W = 0.20
    _DEF_H = 0.10

    def __init__(self, window):
        super().__init__(window)
        self._draft: BodySpec | None = None
        self._phase: str = "anchor"
        self._a: tuple[float, float] | None = None
        self._preview: QGraphicsRectItem | None = None

    # ─────────────────────────────────────────────────────────
    def activate(self) -> None:
        super().activate()
        # Exclude the implicit ground body from the running counter.
        n = sum(1 for b in self.spec.bodies if b.id != "ground") + 1
        self._draft = BodySpec(
            name=f"body{n}",
            mass=1.0,
            inertia=(self._DEF_W ** 2 + self._DEF_H ** 2) / 12.0,
            position=(0.0, 0.0),
            orientation=0.0,
            shape=ShapeSpec(kind="rectangle",
                            params={"width": self._DEF_W,
                                    "height": self._DEF_H}),
        )
        self._phase = "anchor"
        self._a = None
        self.window._inspector.show_draft_body(self._draft)
        self.window.statusBar().showMessage(
            "Rect: set the dimensions in the Inspector, then click on the "
            "canvas to place the body's centre (Esc to cancel)\u2026", 0)

    def deactivate(self) -> None:
        super().deactivate()
        self._clear_preview()
        self._draft = None
        self._a = None
        self._phase = "anchor"
        self.window.statusBar().clearMessage()

    # ─────────────────────────────────────────────────────────
    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton or self._draft is None:
            return False
        p = self.scene.snap(event.scenePos())
        if self._phase == "anchor":
            self._a = p
            self._phase = "rotate"
            self._make_preview(p)
            self.window.statusBar().showMessage(
                "Rect: move the mouse to set the orientation, click to "
                "confirm (Esc to cancel)\u2026", 0)
            return True
        # commit
        self._commit_rect(self._a, p)
        self.window.set_active_tool("select")
        return True

    def mouse_move(self, event) -> bool:
        if (self._phase != "rotate" or self._a is None
                or self._preview is None):
            return False
        p = self.scene.snap(event.scenePos())
        ax, ay = self._a
        phi = math.atan2(p[1] - ay, p[0] - ax)
        self._preview.setRotation(math.degrees(phi))
        return True

    # ─────────────────────────────────────────────────────────
    def _commit_rect(self, a, mouse_b) -> None:
        assert self._draft is not None
        ax, ay = a
        phi = math.atan2(mouse_b[1] - ay, mouse_b[0] - ax)
        body = self._draft
        body.position = (ax, ay)
        body.orientation = phi
        # Re-derive inertia from the (possibly user-edited) sides.
        w = float(body.shape.params.get("width", self._DEF_W))
        h = float(body.shape.params.get("height", self._DEF_H))
        body.inertia = body.mass * (w * w + h * h) / 12.0
        body.locked = True
        if not body.mass_override:
            mp = compute_mass_props(body)
            if mp is not None:
                body.mass, body.inertia = mp
        self.spec.bodies.append(body)
        self.window.add_body_item(body)
        # Canonical centre-of-mass marker.
        self.window.auto_create_cm_marker(body)
        self._clear_preview()
        self._commit()

    # ─────────────────────────────────────────────────────────
    def _make_preview(self, anchor) -> None:
        w = float(self._draft.shape.params.get("width", self._DEF_W))
        h = float(self._draft.shape.params.get("height", self._DEF_H))
        rect = QRectF(-w / 2, -h / 2, w, h)
        self._preview = QGraphicsRectItem(rect)
        pen = QPen(QColor("#3f8cff"))
        pen.setCosmetic(True)
        pen.setWidthF(1.5)
        pen.setStyle(Qt.DashLine)
        self._preview.setPen(pen)
        self._preview.setBrush(QBrush(QColor(63, 140, 255, 60)))
        self._preview.setZValue(1000)
        self.scene.addItem(self._preview)
        ax, ay = anchor
        self._preview.setPos(ax, ay)
        self._preview.setRotation(0.0)

    def _clear_preview(self) -> None:
        if self._preview is not None and self._preview.scene() is not None:
            self.scene.removeItem(self._preview)
        self._preview = None
