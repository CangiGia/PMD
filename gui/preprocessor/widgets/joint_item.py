"""JointItem — graphical glyph for a joint between two markers."""

from __future__ import annotations

from PySide6.QtCore import QRectF
from PySide6.QtGui import QBrush, QColor, QPen
from PySide6.QtWidgets import QGraphicsItem, QGraphicsObject

from ..models import JointSpec


class JointItem(QGraphicsObject):
    """Generic joint glyph (placed at the i-marker's world position).

    Subclasses can refine the visual (a circle for revolute, two
    parallel bars for translational, etc.). The base class draws a
    ring whose colour depends on the joint kind.
    """

    _RADIUS_PX = 7.0

    KIND_COLORS = {
        "RevJoint":     "#5b8cf2",
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
        self.setZValue(450)

    # ─────────────────────────────────────────────────────────
    def boundingRect(self) -> QRectF:
        m = 0.04
        return QRectF(-m, -m, 2 * m, 2 * m)

    def paint(self, painter, option, widget=None):
        view = self.scene().views()[0] if self.scene() and self.scene().views() else None
        s = abs(view.transform().m11()) if view else 400.0
        r = self._RADIUS_PX / s

        color = QColor(self.KIND_COLORS.get(self.spec.kind, "#1c2033"))
        pen = QPen(color)
        pen.setCosmetic(True)
        pen.setWidthF(2.0)
        painter.setPen(pen)
        painter.setBrush(QBrush(QColor(255, 255, 255, 200)))
        painter.drawEllipse(QRectF(-r, -r, 2 * r, 2 * r))
