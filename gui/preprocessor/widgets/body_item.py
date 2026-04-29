"""BodyItem — graphical representation of a BodySpec on the canvas."""

from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import BodySpec


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
        else:
            self._w = 0.20
            self._h = 0.10
        self._kind = kind

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
        rect = QRectF(-self._w / 2, -self._h / 2, self._w, self._h)

        fill = self._FILL_SELECTED if self.isSelected() else self._FILL
        pen = QPen(self._STROKE)
        pen.setCosmetic(True)
        pen.setWidthF(self._STROKE_W_PX)

        painter.setBrush(QBrush(fill))
        painter.setPen(pen)
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
