"""BodyItem — graphical representation of a BodySpec on the canvas."""

from __future__ import annotations

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
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
        self.setAcceptHoverEvents(True)
        self.setPos(spec.position[0], spec.position[1])
        self.setRotation(0.0)  # rotation handled separately

        # Default rectangle dimensions (used when shape is None).
        if spec.shape and spec.shape.kind == "rectangle":
            self._w = float(spec.shape.params.get("width",  0.20))
            self._h = float(spec.shape.params.get("height", 0.10))
        else:
            self._w = 0.20
            self._h = 0.10

    # ──────────────────────────────────────────────────────────
    # QGraphicsItem mandatory overrides
    # ──────────────────────────────────────────────────────────

    def boundingRect(self) -> QRectF:
        m = self._STROKE_W_PX  # leave room for the stroke
        return QRectF(-self._w / 2 - m,
                      -self._h / 2 - m,
                       self._w + 2 * m,
                       self._h + 2 * m)

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
        return super().itemChange(change, value)
