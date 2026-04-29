"""BodyLinkTool \u2014 inspector-driven 2-click placement of a link body.

Workflow
--------
1. User selects the *Link* tool. A draft :class:`BodySpec` is built and
   exposed in the Inspector so the user can tune length / thickness /
   mass / initial velocity *before* anchoring it on the canvas.
2. First left-click on the canvas anchors marker A at the click point.
3. Mouse moves rotate the link around A, keeping the user-chosen length
   fixed.
4. Second left-click commits the body, places markers at the two ends,
   marks the body as ``locked`` and reverts to the *Select* tool.

Esc cancels the in-progress placement and discards the draft.
"""

from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsRectItem

from ..models import BodySpec, MarkerSpec, ShapeSpec, compute_mass_props
from .tool_base import Tool


class BodyLinkTool(Tool):
    name = "body_link"
    cursor = Qt.CrossCursor

    _DEF_LEN = 0.40
    _DEF_THK = 0.02

    def __init__(self, window):
        super().__init__(window)
        self._draft: BodySpec | None = None
        self._phase: str = "anchor"     # "anchor" | "rotate"
        self._a: tuple[float, float] | None = None
        self._preview: QGraphicsRectItem | None = None

    # ─────────────────────────────────────────────────────────
    def activate(self) -> None:
        super().activate()
        n = len(self.spec.bodies) + 1
        self._draft = BodySpec(
            name=f"link{n}",
            mass=1.0,
            inertia=self._DEF_LEN ** 2 / 12.0,
            position=(0.0, 0.0),
            orientation=0.0,
            shape=ShapeSpec(kind="link",
                            params={"length": self._DEF_LEN,
                                    "thickness": self._DEF_THK}),
        )
        self._phase = "anchor"
        self._a = None
        self.window._inspector.show_draft_body(self._draft)
        self.window.statusBar().showMessage(
            "Link: set the dimensions in the Inspector, then click on the "
            "canvas to place marker A (Esc to cancel)\u2026", 0)

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
                "Link: move the mouse to set the orientation, click to "
                "confirm (Esc to cancel)\u2026", 0)
            return True

        # phase == "rotate" → commit
        self._commit_link(self._a, p)
        self.window.set_active_tool("select")
        return True

    def mouse_move(self, event) -> bool:
        if (self._phase != "rotate" or self._a is None
                or self._preview is None or self._draft is None):
            return False
        p = self.scene.snap(event.scenePos())
        self._update_preview_to(p)
        return True

    # ─────────────────────────────────────────────────────────
    def _commit_link(self, a, mouse_b) -> None:
        assert self._draft is not None
        L = float(self._draft.shape.params.get("length", self._DEF_LEN))
        ax, ay = a
        phi = math.atan2(mouse_b[1] - ay, mouse_b[0] - ax)
        cx = ax + 0.5 * L * math.cos(phi)
        cy = ay + 0.5 * L * math.sin(phi)

        body = self._draft
        body.position = (cx, cy)
        body.orientation = phi
        body.inertia = body.mass * L * L / 12.0
        body.locked = True
        # Re-derive mass / inertia from material density once the
        # final shape parameters are settled.
        if not body.mass_override:
            mp = compute_mass_props(body)
            if mp is not None:
                body.mass, body.inertia = mp
        self.spec.bodies.append(body)

        m_a = MarkerSpec(name=f"{body.name}.A", body_id=body.id,
                         local_position=(-L / 2, 0.0), theta=0.0)
        m_b = MarkerSpec(name=f"{body.name}.B", body_id=body.id,
                         local_position=(+L / 2, 0.0), theta=0.0)
        self.spec.markers.extend([m_a, m_b])

        self.window.add_body_item(body)
        self.window.add_marker_item(m_a)
        self.window.add_marker_item(m_b)
        # Canonical centre-of-mass marker.
        self.window.auto_create_cm_marker(body)
        self._clear_preview()
        self._commit()

    # ─────────────────────────────────────────────────────────
    def _make_preview(self, anchor) -> None:
        L = float(self._draft.shape.params.get("length", self._DEF_LEN))
        thk = float(self._draft.shape.params.get("thickness", self._DEF_THK))
        rect = QRectF(-L / 2, -thk / 2, L, thk)
        self._preview = QGraphicsRectItem(rect)
        pen = QPen(QColor("#3f8cff"))
        pen.setCosmetic(True)
        pen.setWidthF(1.5)
        pen.setStyle(Qt.DashLine)
        self._preview.setPen(pen)
        self._preview.setBrush(QBrush(QColor(63, 140, 255, 70)))
        self._preview.setZValue(1000)
        self.scene.addItem(self._preview)
        ax, ay = anchor
        self._preview.setPos(ax + L / 2, ay)
        self._preview.setRotation(0.0)

    def _update_preview_to(self, mouse_pt) -> None:
        L = float(self._draft.shape.params.get("length", self._DEF_LEN))
        ax, ay = self._a
        phi = math.atan2(mouse_pt[1] - ay, mouse_pt[0] - ax)
        cx = ax + 0.5 * L * math.cos(phi)
        cy = ay + 0.5 * L * math.sin(phi)
        self._preview.setPos(cx, cy)
        self._preview.setRotation(math.degrees(phi))

    def _clear_preview(self) -> None:
        if self._preview is not None and self._preview.scene() is not None:
            self.scene.removeItem(self._preview)
        self._preview = None
