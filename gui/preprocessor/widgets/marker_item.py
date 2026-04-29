"""MarkerItem \u2014 body-attached reference frame visual."""

from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QFont, QPainterPath, QPen, QPolygonF
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import MarkerSpec


class MarkerItem(QGraphicsObject):
    """Reference-frame glyph for a body-fixed marker.

    Drawn as two arrows (red X, green Y) emanating from the marker
    origin, with a small picking dot at the origin and a label. Sizes
    are screen-constant via cosmetic pens and the view's current scale.
    """

    # Screen-space sizes (pixels)
    _AXIS_LEN_PX  = 22.0
    _ARROW_W_PX   = 6.0
    _ARROW_H_PX   = 8.0
    _ORIGIN_R_PX  = 2.5
    _PICK_R_PX    = 8.0
    _LABEL_DX_PX  = 4.0
    _LABEL_DY_PX  = -4.0

    # Colors
    _C_X      = QColor("#d65a5a")
    _C_Y      = QColor("#3aa35a")
    _C_X_HI   = QColor("#ff5050")
    _C_Y_HI   = QColor("#22c55e")
    _C_DOT    = QColor("#1c2033")
    _C_DOT_HI = QColor("#ff8c00")
    _C_LABEL  = QColor("#1c2033")

    def __init__(self, spec: MarkerSpec, parent_body=None):
        super().__init__(parent_body)
        self.spec = spec
        self._highlight = False
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setZValue(500)
        self.setPos(*spec.local_position)
        self.setRotation(math.degrees(spec.theta))
        self.setAcceptHoverEvents(True)

    # ─────────────────────────────────────────────────────────
    def _scale(self) -> float:
        view = self.scene().views()[0] if self.scene() and self.scene().views() else None
        return abs(view.transform().m11()) if view else 400.0

    def boundingRect(self) -> QRectF:
        s = self._scale()
        L = (self._AXIS_LEN_PX + self._ARROW_H_PX + 16.0) / s
        return QRectF(-L, -L, 2 * L, 2 * L)

    def shape(self):
        # Picking constrained to a small disc at the origin (axes are
        # visual only, can't be clicked through).
        s = self._scale()
        r = self._PICK_R_PX / s
        path = QPainterPath()
        path.addEllipse(QRectF(-r, -r, 2 * r, 2 * r))
        return path

    def paint(self, painter, option, widget=None):
        s = self._scale()
        L  = self._AXIS_LEN_PX / s
        aw = self._ARROW_W_PX / s
        ah = self._ARROW_H_PX / s
        r0 = self._ORIGIN_R_PX / s

        cx = self._C_X_HI if self._highlight else self._C_X
        cy = self._C_Y_HI if self._highlight else self._C_Y
        cd = self._C_DOT_HI if (self._highlight or self.isSelected()) else self._C_DOT

        # ── X axis (red) ──────────────────────────────────
        pen_x = QPen(cx); pen_x.setCosmetic(True); pen_x.setWidthF(1.8)
        painter.setPen(pen_x)
        painter.drawLine(QPointF(0.0, 0.0), QPointF(L - ah, 0.0))
        painter.setBrush(QBrush(cx))
        painter.drawPolygon(QPolygonF([
            QPointF(L,       0.0),
            QPointF(L - ah, +aw / 2.0),
            QPointF(L - ah, -aw / 2.0),
        ]))

        # ── Y axis (green) ────────────────────────────────
        pen_y = QPen(cy); pen_y.setCosmetic(True); pen_y.setWidthF(1.8)
        painter.setPen(pen_y)
        painter.drawLine(QPointF(0.0, 0.0), QPointF(0.0, L - ah))
        painter.setBrush(QBrush(cy))
        painter.drawPolygon(QPolygonF([
            QPointF(0.0,      L),
            QPointF(+aw / 2.0, L - ah),
            QPointF(-aw / 2.0, L - ah),
        ]))

        # ── Origin dot ────────────────────────────────────
        pen_d = QPen(self._C_DOT); pen_d.setCosmetic(True); pen_d.setWidthF(1.0)
        painter.setPen(pen_d)
        painter.setBrush(QBrush(cd))
        painter.drawEllipse(QRectF(-r0, -r0, 2 * r0, 2 * r0))

        # ── Selection / pick highlight ring ───────────────
        if self.isSelected() or self._highlight:
            ring = QPen(QColor("#ff8c00"))
            ring.setCosmetic(True); ring.setWidthF(1.2)
            ring.setStyle(Qt.DashLine)
            painter.setPen(ring)
            painter.setBrush(Qt.NoBrush)
            rs = (self._PICK_R_PX + 2.0) / s
            painter.drawEllipse(QRectF(-rs, -rs, 2 * rs, 2 * rs))

        # ── Label ─────────────────────────────────────────
        name = self.spec.name or self.spec.id
        if name:
            painter.save()
            painter.translate(self._LABEL_DX_PX / s, self._LABEL_DY_PX / s)
            # Cancel the world-units transform AND the canvas' y-flip
            # so the label is always horizontal and right-side up.
            painter.scale(1.0 / s, -1.0 / s)
            font = QFont()
            font.setPointSizeF(7.5)
            painter.setFont(font)
            pen_l = QPen(self._C_LABEL); pen_l.setCosmetic(True)
            painter.setPen(pen_l)
            painter.drawText(QPointF(0.0, 0.0), name)
            painter.restore()

    # ─────────────────────────────────────────────────────────
    def set_highlighted(self, on: bool) -> None:
        self._highlight = bool(on)
        self.update()
