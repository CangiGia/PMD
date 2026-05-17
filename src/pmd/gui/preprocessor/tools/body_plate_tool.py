"""BodyPlateTool — 3-click placement of a triangular plate body.

Workflow
--------
1. User selects the *Plate* tool. A draft :class:`BodySpec` is built and
   shown in the Inspector.
2. First left-click → anchors vertex v1.
3. Second left-click → anchors vertex v2 (a line preview tracks the mouse
   between clicks 1 and 2).
4. Third left-click → anchors vertex v3; the triangle is committed
   automatically and the tool reverts to Select.

The centroid of the three clicked world-frame points is used as
``body.position`` (required by the PMD solver's diagonal mass matrix).
Local-frame vertices stored in ``spec.shape.params["vertices"]`` are
already centred on the centroid.

Esc at any point cancels the placement and discards the draft.
"""

from __future__ import annotations

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QBrush, QColor, QPen, QPolygonF
from PySide6.QtWidgets import QGraphicsPolygonItem

from ..models import BodySpec, ShapeSpec, compute_mass_props
from .tool_base import Tool


class BodyPlateTool(Tool):
    name = "body_plate"
    cursor = Qt.CrossCursor

    def __init__(self, window):
        super().__init__(window)
        self._draft: BodySpec | None = None
        self._phase: str = "v1"               # "v1" | "v2" | "v3"
        self._pts: list[tuple[float, float]] = []
        self._preview: QGraphicsPolygonItem | None = None

    # ─────────────────────────────────────────────────────────
    def activate(self) -> None:
        super().activate()
        n = sum(1 for b in self.spec.bodies if b.id != "ground") + 1
        self._draft = BodySpec(
            name=f"plate{n}",
            mass=1.0,
            inertia=0.01,
            position=(0.0, 0.0),
            orientation=0.0,
            shape=ShapeSpec(kind="plate", params={"vertices": []}),
        )
        self._phase = "v1"
        self._pts = []
        self.window._inspector.show_draft_body(self._draft)
        self.window.statusBar().showMessage(
            "Plate: click the first vertex (CCW order, Esc to cancel)\u2026", 0)

    def deactivate(self) -> None:
        super().deactivate()
        self._clear_preview()
        self._draft = None
        self._pts = []
        self._phase = "v1"
        self.window.statusBar().clearMessage()

    # ─────────────────────────────────────────────────────────
    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton or self._draft is None:
            return False
        p = self.scene.snap(event.scenePos())
        self._pts.append(p)

        if self._phase == "v1":
            self._phase = "v2"
            self.window.statusBar().showMessage(
                "Plate: click the second vertex (Esc to cancel)\u2026", 0)
            return True

        if self._phase == "v2":
            self._phase = "v3"
            self.window.statusBar().showMessage(
                "Plate: click the third vertex CCW to commit "
                "(Esc to cancel)\u2026", 0)
            return True

        # phase == "v3" → commit
        self._commit_plate()
        self.window.set_active_tool("select")
        return True

    def mouse_move(self, event) -> bool:
        if self._draft is None or self._phase == "v1":
            return False
        p = self.scene.snap(event.scenePos())
        self._update_preview(p)
        return True

    # ─────────────────────────────────────────────────────────
    def _update_preview(self, cursor: tuple[float, float]) -> None:
        self._clear_preview()
        pts = self._pts + [cursor]   # 2 or 3 points
        qpts = [QPointF(x, y) for x, y in pts]
        poly = QPolygonF(qpts)
        item = QGraphicsPolygonItem(poly)
        pen = QPen(QColor("#3f8cff"))
        pen.setCosmetic(True)
        pen.setWidthF(1.5)
        pen.setStyle(Qt.DashLine)
        item.setPen(pen)
        item.setBrush(
            QBrush(QColor(63, 140, 255, 60))
            if len(pts) == 3
            else QBrush(Qt.NoBrush)
        )
        item.setZValue(1000)
        self.scene.addItem(item)
        self._preview = item

    def _clear_preview(self) -> None:
        if self._preview is not None and self._preview.scene() is not None:
            self.scene.removeItem(self._preview)
        self._preview = None

    def _commit_plate(self) -> None:
        assert self._draft is not None
        assert len(self._pts) == 3
        (x0, y0), (x1, y1), (x2, y2) = self._pts
        cx = (x0 + x1 + x2) / 3.0
        cy = (y0 + y1 + y2) / 3.0
        local_verts = [
            (x0 - cx, y0 - cy),
            (x1 - cx, y1 - cy),
            (x2 - cx, y2 - cy),
        ]

        body = self._draft
        body.position = (cx, cy)
        body.orientation = 0.0
        body.shape.params["vertices"] = local_verts
        body.locked = True

        if not body.mass_override:
            mp = compute_mass_props(body)
            if mp is not None:
                body.mass, body.inertia = mp

        self.spec.bodies.append(body)
        self.window.add_body_item(body)
        self.window.auto_create_cm_marker(body)
        self._clear_preview()
        self._commit()
