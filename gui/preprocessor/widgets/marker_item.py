"""MarkerItem — body-attached reference frame visual."""

from __future__ import annotations

import math

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import MarkerSpec


class MarkerItem(QGraphicsObject):
    """Small triad indicating a body-fixed marker.

    Drawn as a filled circle plus a short red/green axis cross. The
    item is meant to be parented either to a :class:`BodyItem` (so it
    follows the body's pose) or directly to the scene for ground
    markers.
    """

    _RADIUS_PX  = 4.0
    _AXIS_LEN_PX = 14.0
    _FILL       = QColor("#ffd24d")
    _FILL_HI    = QColor("#ff8c00")
    _STROKE     = QColor("#1c2033")

    def __init__(self, spec: MarkerSpec, parent_body=None):
        super().__init__(parent_body)
        self.spec = spec
        self._highlight = False
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.ItemIgnoresTransformations, False)
        self.setZValue(500)
        self.setPos(*spec.local_position)
        self.setRotation(math.degrees(spec.theta))
        self.setAcceptHoverEvents(True)

    # ─────────────────────────────────────────────────────────
    def boundingRect(self) -> QRectF:
        # Use a generous box in scene units (we draw cosmetic px shapes).
        m = 0.04
        return QRectF(-m, -m, 2 * m, 2 * m)

    def paint(self, painter, option, widget=None):
        # Switch to device-pixel scale so the marker glyph stays at a
        # constant size on screen regardless of zoom.
        view = self.scene().views()[0] if self.scene() and self.scene().views() else None
        s = abs(view.transform().m11()) if view else 400.0
        r  = self._RADIUS_PX  / s
        ax = self._AXIS_LEN_PX / s

        fill = self._FILL_HI if self._highlight else self._FILL
        painter.setBrush(QBrush(fill))
        pen = QPen(self._STROKE)
        pen.setCosmetic(True)
        pen.setWidthF(1.0)
        painter.setPen(pen)
        painter.drawEllipse(QRectF(-r, -r, 2 * r, 2 * r))

        # X axis (red)
        pen_x = QPen(QColor("#d65a5a"))
        pen_x.setCosmetic(True)
        pen_x.setWidthF(1.5)
        painter.setPen(pen_x)
        painter.drawLine(0, 0, ax, 0)
        # Y axis (green)
        pen_y = QPen(QColor("#3aa35a"))
        pen_y.setCosmetic(True)
        pen_y.setWidthF(1.5)
        painter.setPen(pen_y)
        painter.drawLine(0, 0, 0, ax)

    # ─────────────────────────────────────────────────────────
    def set_highlighted(self, on: bool) -> None:
        self._highlight = bool(on)
        self.update()
