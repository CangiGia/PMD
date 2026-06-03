"""BodyItem — graphical representation of a BodySpec on the canvas."""

from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import BodySpec


class RotationHandle(QGraphicsObject):
    """Blue dot child of BodyItem that lets the user rotate a body by dragging.

    Positioned at the top edge of the parent body (local +y direction).
    Created lazily by :meth:`BodyItem.show_rotation_handle`.
    """

    _R_SCENE = 0.012  # handle radius in scene units (metres)

    def __init__(self, parent_body: "BodyItem", on_done=None):
        super().__init__(parent=parent_body)
        self._body = parent_body
        self._on_done = on_done
        self._dragging = False
        gap = max(parent_body._h * 0.15, 0.020)
        self.setPos(0.0, parent_body._h / 2.0 + gap)
        self.setZValue(200)
        self.setCursor(Qt.ClosedHandCursor)
        self.setFlag(QGraphicsItem.ItemIsMovable, False)
        self.setFlag(QGraphicsItem.ItemIsSelectable, False)
        self.setAcceptedMouseButtons(Qt.LeftButton)

    # ── QGraphicsItem interface ────────────────────────────────

    def boundingRect(self) -> QRectF:
        r = self._R_SCENE
        return QRectF(-r, -r, 2.0 * r, 2.0 * r)

    def paint(self, painter, option, widget=None):
        r = self._R_SCENE
        painter.setBrush(QBrush(QColor("#3f8cff")))
        painter.setPen(QPen(QColor("white"), 0))
        painter.drawEllipse(QRectF(-r, -r, 2.0 * r, 2.0 * r))

    # ── Mouse handling ────────────────────────────────────────

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._dragging = True
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self._dragging:
            super().mouseMoveEvent(event)
            return
        center = self._body.scenePos()
        cursor = event.scenePos()
        dx = cursor.x() - center.x()
        dy = cursor.y() - center.y()
        phi = math.atan2(-dx, dy)
        self._body.setRotation(math.degrees(phi))
        # spec.orientation is updated automatically by BodyItem.itemChange
        event.accept()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton and self._dragging:
            self._dragging = False
            if self._on_done is not None:
                self._on_done(self._body)
            event.accept()
            return
        super().mouseReleaseEvent(event)


class BodyItem(QGraphicsObject):
    """Selectable / movable rectangular body on the editing canvas.

    The item draws itself in *local* coordinates centred on the body's
    centroid; its position in the scene is the body's CoM position.

    Parameters
    ----------
    spec : BodySpec
        The underlying body specification. The item edits this in place
        when the user drags it.
    """

    _FILL          = QColor("#cfe1ff")
    _FILL_SELECTED = QColor("#9ec5ff")
    _STROKE        = QColor("#3f8cff")
    _STROKE_W_PX   = 1.5

    def __init__(self, spec: BodySpec):
        super().__init__()
        self.spec = spec
        # Movability is controlled by the active tool. By default a body
        # is *not* draggable so the user can't accidentally displace it
        # while editing — the dedicated MoveTool re-enables this flag.
        self.setFlag(QGraphicsItem.ItemIsMovable, False)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
        self.setAcceptHoverEvents(True)
        self.setPos(spec.position[0], spec.position[1])
        self.setRotation(math.degrees(spec.orientation))

        kind = spec.shape.kind if spec.shape else "rectangle"
        params = spec.shape.params if spec.shape else {}
        if kind == "link":
            self._w = float(params.get("length",    0.30))
            self._h = float(params.get("thickness", 0.02))
        elif kind == "rectangle":
            self._w = float(params.get("width",  0.20))
            self._h = float(params.get("height", 0.10))
        elif kind == "plate":
            raw = params.get("vertices", [[-0.1, -0.1], [0.1, -0.1], [0.0, 0.1]])
            self._plate_verts = [(float(v[0]), float(v[1])) for v in raw]
            xs = [v[0] for v in self._plate_verts]
            ys = [v[1] for v in self._plate_verts]
            self._w = max(xs) - min(xs) if xs else 0.20
            self._h = max(ys) - min(ys) if ys else 0.10
        else:
            self._w = 0.20
            self._h = 0.10
        self._kind = kind
        self._rotation_handle: RotationHandle | None = None

    # ──────────────────────────────────────────────────────────
    # Rotation handle
    # ──────────────────────────────────────────────────────────

    def show_rotation_handle(self, on: bool, on_done=None) -> None:
        """Show or hide the rotation drag handle.

        The ground body never gets a handle. The handle is created
        lazily on first call with ``on=True`` and merely hidden
        afterwards — so the same object is reused across tool switches.
        """
        if on and self.spec.id != "ground":
            if self._rotation_handle is None:
                self._rotation_handle = RotationHandle(self, on_done=on_done)
            else:
                self._rotation_handle._on_done = on_done
            self._rotation_handle.setVisible(True)
        elif self._rotation_handle is not None:
            self._rotation_handle.setVisible(False)

    # ──────────────────────────────────────────────────────────
    # QGraphicsItem mandatory overrides
    # ──────────────────────────────────────────────────────────

    def boundingRect(self) -> QRectF:
        m = self._STROKE_W_PX  # leave room for the stroke
        return QRectF(-self._w / 2 - m,
                      -self._h / 2 - m,
                       self._w + 2 * m,
                       self._h + 2 * m)

    def shape(self):
        # Picking is constrained to the *actual* body rectangle (in
        # local coords). Without this the default shape() == bounding
        # rect, which carries the stroke-width margin in scene units
        # (i.e. ~1 m halo around every body) and makes everything
        # near a body select that body instead of what the user
        # actually clicked on.
        from PySide6.QtGui import QPainterPath
        path = QPainterPath()
        path.addRect(QRectF(-self._w / 2, -self._h / 2,
                            self._w, self._h))
        return path

    def paint(self, painter, option, widget=None):
        fill = self._FILL_SELECTED if self.isSelected() else self._FILL
        pen = QPen(self._STROKE)
        pen.setCosmetic(True)
        pen.setWidthF(self._STROKE_W_PX)

        painter.setBrush(QBrush(fill))
        painter.setPen(pen)

        if self._kind == "plate":
            from PySide6.QtGui import QPolygonF
            from PySide6.QtCore import QPointF
            poly = QPolygonF([QPointF(x, y) for x, y in self._plate_verts])
            painter.drawPolygon(poly)
            # CoM cross at origin (centroid = local (0,0))
            r = max(self._w, self._h) * 0.05
            painter.setPen(QPen(QColor("#1c2033"), 0))
            painter.drawLine(-r, 0.0, r, 0.0)
            painter.drawLine(0.0, -r, 0.0, r)
        else:
            rect = QRectF(-self._w / 2, -self._h / 2, self._w, self._h)
            painter.drawRect(rect)
            # CoM cross
            painter.setPen(QPen(QColor("#1c2033"), 0))
            painter.drawLine(-self._w * 0.05, 0, self._w * 0.05, 0)
            painter.drawLine(0, -self._h * 0.05, 0, self._h * 0.05)

    # ──────────────────────────────────────────────────────────
    # Sync spec on move
    # ──────────────────────────────────────────────────────────

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemPositionHasChanged:
            self.spec.position = (float(self.pos().x()), float(self.pos().y()))
        elif change == QGraphicsItem.ItemRotationHasChanged:
            self.spec.orientation = math.radians(float(self.rotation()))
        return super().itemChange(change, value)
