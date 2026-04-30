"""JointItem — graphical glyph for a joint between two markers."""

from __future__ import annotations

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import JointSpec


class JointItem(QGraphicsObject):
    """Generic joint glyph (placed at the i-marker's world position).

    The base class draws a *filled*, brightly coloured disc whose hue
    depends on the joint kind, with a thin dark outline so it stays
    visible on top of the marker glyphs that necessarily overlap it.
    """

    _RADIUS_PX = 9.0

    KIND_COLORS = {
        "RevJoint":     "#ff2bd6",   # bright fuchsia (hinge)
        "TranJoint":    "#7a5bf2",
        "RevRevJoint":  "#3a8d6a",
        "RevTranJoint": "#3a8d8d",
        "RigidJoint":   "#1c2033",
        "DiscJoint":    "#cf8a2a",
    }

    def __init__(self, spec: JointSpec):
        super().__init__()
        self.spec = spec
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        # Above markers (z=500) and bodies (default 0) so the joint
        # glyph is never visually buried under a marker.
        self.setZValue(600)

    # ─────────────────────────────────────────────────────────
    def boundingRect(self) -> QRectF:
        m = 0.04
        return QRectF(-m, -m, 2 * m, 2 * m)

    def paint(self, painter, option, widget=None):
        view = self.scene().views()[0] if self.scene() and self.scene().views() else None
        s = abs(view.transform().m11()) if view else 400.0
        r = self._RADIUS_PX / s

        color = QColor(self.KIND_COLORS.get(self.spec.kind, "#1c2033"))
        # Filled bright disc with a thin dark outline for contrast.
        painter.setBrush(QBrush(color))
        outline = QPen(QColor("#1c2033"))
        outline.setCosmetic(True)
        outline.setWidthF(1.2)
        painter.setPen(outline)
        painter.drawEllipse(QRectF(-r, -r, 2 * r, 2 * r))

        # Selection halo.
        if self.isSelected():
            halo = QPen(QColor("#ff8c00"))
            halo.setCosmetic(True)
            halo.setWidthF(1.4)
            halo.setStyle(Qt.DashLine)
            painter.setPen(halo)
            painter.setBrush(Qt.NoBrush)
            rr = (self._RADIUS_PX + 3.0) / s
            painter.drawEllipse(QRectF(-rr, -rr, 2 * rr, 2 * rr))
