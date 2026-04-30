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
    _AXIS_LEN_PX  = 44.0
    _ARROW_W_PX   = 11.0
    _ARROW_H_PX   = 14.0
    _ORIGIN_R_PX  = 4.0
    _PICK_R_PX    = 12.0
    _LABEL_DX_PX  = -2.0   # small horizontal nudge
    _LABEL_DY_PX  = 22.0   # below the origin (positive = screen-down)
    _LABEL_PT     = 11.0
    _AXIS_W_PX    = 2.2
    # Scene-pixel radius within which two markers are considered
    # co-located for the purpose of label stacking. Roughly the size
    # of the picking disc so a joint's two markers always group.
    _LABEL_GROUP_PX = 14.0

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
        pen_x = QPen(cx); pen_x.setCosmetic(True); pen_x.setWidthF(self._AXIS_W_PX)
        painter.setPen(pen_x)
        painter.drawLine(QPointF(0.0, 0.0), QPointF(L - ah, 0.0))
        painter.setBrush(QBrush(cx))
        painter.drawPolygon(QPolygonF([
            QPointF(L,       0.0),
            QPointF(L - ah, +aw / 2.0),
            QPointF(L - ah, -aw / 2.0),
        ]))

        # ── Y axis (green) ──────────────────────────────────────────────
        pen_y = QPen(cy); pen_y.setCosmetic(True); pen_y.setWidthF(self._AXIS_W_PX)
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
            # Total scene-rotation = parent body rotation + marker theta.
            # We *cancel* it so the label stays parallel to the global X
            # axis regardless of how the body is oriented. The y-flip of
            # the canvas view is also cancelled so the text reads upright.
            parent = self.parentItem()
            total_deg = float(self.rotation())
            if parent is not None:
                total_deg += float(parent.rotation())

            # Anti-overlap: when several markers share (almost) the same
            # scene origin (typical at a joint endpoint), stack their
            # labels vertically so they don't all draw on top of each
            # other. Each marker keeps its label close to its own dot.
            row = self._label_row_index()

            painter.save()
            painter.rotate(-total_deg)
            painter.scale(1.0 / s, -1.0 / s)        # painter local = pixels
            line_h = self._LABEL_PT * 1.35
            painter.translate(self._LABEL_DX_PX,
                              self._LABEL_DY_PX + row * line_h)
            font = QFont()
            font.setPointSizeF(self._LABEL_PT)
            painter.setFont(font)
            pen_l = QPen(self._C_LABEL); pen_l.setCosmetic(True)
            painter.setPen(pen_l)
            painter.drawText(QPointF(0.0, 0.0), name)
            painter.restore()

    # ─────────────────────────────────────────────────────────
    def _label_row_index(self) -> int:
        """Return this marker's row when stacking labels at the same spot.

        Two markers are considered co-located if their scene origins
        are within :data:`_LABEL_GROUP_PX` pixels. Within a group the
        rows are assigned by a stable sort (id, then python id) so the
        layout doesn't flicker between repaints.
        """
        scene = self.scene()
        if scene is None:
            return 0
        my_pos = self.scenePos()
        s = self._scale()
        tol = self._LABEL_GROUP_PX / max(s, 1e-9)
        peers: list[MarkerItem] = []
        for it in scene.items():
            if not isinstance(it, MarkerItem):
                continue
            if not (it.spec.name or it.spec.id):
                continue
            d = it.scenePos() - my_pos
            if abs(d.x()) <= tol and abs(d.y()) <= tol:
                peers.append(it)
        if len(peers) <= 1:
            return 0
        peers.sort(key=lambda m: (str(m.spec.id), id(m)))
        try:
            return peers.index(self)
        except ValueError:
            return 0

    # ─────────────────────────────────────────────────────────
    def set_highlighted(self, on: bool) -> None:
        self._highlight = bool(on)
        self.update()
