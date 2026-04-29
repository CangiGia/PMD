"""CanvasScene — QGraphicsScene with engineering grid background.

Coordinate convention
---------------------
* Scene units are *meters*; one Qt scene unit equals one meter.
* +x is to the right, +y is **upward** (engineering convention).
* :class:`CanvasView` flips the Y axis at render time so user coords
  match the scene coords visually.
"""

from __future__ import annotations

from PySide6.QtCore import QLineF, QRectF, Qt
from PySide6.QtGui import QColor, QPen
from PySide6.QtWidgets import QGraphicsScene


class CanvasScene(QGraphicsScene):
    """2D editing scene with adaptive grid background.

    Parameters
    ----------
    grid_step : float
        Grid spacing in scene units (meters). Default 0.05 m (50 mm).
    parent : QObject or None
        Optional parent.
    """

    # Visual tokens — kept in sync with the light theme palette.
    _GRID_MINOR  = QColor("#e4e7ec")
    _GRID_MAJOR  = QColor("#c8cdd6")
    _AXIS_X      = QColor("#d65a5a")  # red
    _AXIS_Y      = QColor("#3aa35a")  # green
    _BG          = QColor("#fafbfc")

    def __init__(self, grid_step: float = 0.05, parent=None):
        super().__init__(parent)
        self.grid_step = grid_step
        self.major_every = 10  # one major line every N minor steps

        self.active_tool = None  # set by PreProcessorWindow

        # Generous initial scene rect — expanded automatically by items.
        self.setSceneRect(-1.0, -1.0, 2.0, 2.0)
        self.setBackgroundBrush(self._BG)

    # ───────────────────────────────────────────────────────
    # Tool dispatch
    # ───────────────────────────────────────────────────────

    def mousePressEvent(self, event):
        if self.active_tool is not None and self.active_tool.mouse_press(event):
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.active_tool is not None and self.active_tool.mouse_move(event):
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.active_tool is not None and self.active_tool.mouse_release(event):
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        if self.active_tool is not None and self.active_tool.key_press(event):
            event.accept()
            return
        super().keyPressEvent(event)

    # ──────────────────────────────────────────────────────────
    # Background grid
    # ──────────────────────────────────────────────────────────

    def drawBackground(self, painter, rect: QRectF):
        super().drawBackground(painter, rect)

        # Adapt grid density to current zoom: skip minor grid when
        # too dense (less than 4 px between lines).
        view = self.views()[0] if self.views() else None
        if view is None:
            return
        px_per_unit = view.transform().m11()
        if px_per_unit <= 0:
            return

        step = self.grid_step
        # Scale step up if the minor grid would render too tight.
        while step * px_per_unit < 6:
            step *= self.major_every

        major_step = step * self.major_every

        left   = rect.left()
        right  = rect.right()
        top    = rect.top()
        bottom = rect.bottom()

        # Snap to step
        x0 = step * (int(left  / step) - 1)
        y0 = step * (int(top   / step) - 1)

        minor_pen = QPen(self._GRID_MINOR)
        minor_pen.setCosmetic(True)
        minor_pen.setWidth(0)

        major_pen = QPen(self._GRID_MAJOR)
        major_pen.setCosmetic(True)
        major_pen.setWidth(0)

        # Vertical lines
        minor_lines: list[QLineF] = []
        major_lines: list[QLineF] = []
        x = x0
        eps = step * 1e-3
        while x <= right:
            line = QLineF(x, top, x, bottom)
            if abs(x % major_step) < eps or abs(major_step - (x % major_step)) < eps:
                major_lines.append(line)
            else:
                minor_lines.append(line)
            x += step

        # Horizontal lines
        y = y0
        while y <= bottom:
            line = QLineF(left, y, right, y)
            if abs(y % major_step) < eps or abs(major_step - (y % major_step)) < eps:
                major_lines.append(line)
            else:
                minor_lines.append(line)
            y += step

        painter.setPen(minor_pen)
        painter.drawLines(minor_lines)
        painter.setPen(major_pen)
        painter.drawLines(major_lines)

        # ── World axes ─────────────────────────────────────────
        if left <= 0 <= right:
            axis_y = QPen(self._AXIS_Y)
            axis_y.setCosmetic(True)
            axis_y.setWidth(0)
            painter.setPen(axis_y)
            painter.drawLine(QLineF(0, top, 0, bottom))
        if top <= 0 <= bottom:
            axis_x = QPen(self._AXIS_X)
            axis_x.setCosmetic(True)
            axis_x.setWidth(0)
            painter.setPen(axis_x)
            painter.drawLine(QLineF(left, 0, right, 0))

    # ──────────────────────────────────────────────────────────
    # Snap helpers
    # ──────────────────────────────────────────────────────────

    # Pixel radius around the cursor inside which an existing marker
    # acts as a snap target (priority over the grid).
    _MARKER_SNAP_PX = 12.0

    def _px_per_unit(self) -> float:
        v = self.views()[0] if self.views() else None
        return abs(v.transform().m11()) if v else 1.0

    def _nearest_marker_pos(self, point) -> tuple[float, float] | None:
        """Return the world (x, y) of the marker closest to ``point``
        when within the snap radius, else ``None``."""
        # Local import to avoid a circular dependency at module load.
        from .marker_item import MarkerItem  # noqa: WPS433

        ppu = self._px_per_unit() or 1.0
        tol = self._MARKER_SNAP_PX / ppu
        tol2 = tol * tol
        best: tuple[float, float] | None = None
        best_d2 = tol2
        px, py = float(point.x()), float(point.y())
        for it in self.items():
            if not isinstance(it, MarkerItem):
                continue
            sp = it.scenePos()
            dx = sp.x() - px
            dy = sp.y() - py
            d2 = dx * dx + dy * dy
            if d2 <= best_d2:
                best_d2 = d2
                best = (sp.x(), sp.y())
        return best

    def snap(self, point) -> tuple[float, float]:
        """Snap a scene point to the nearest marker (within a small
        pixel tolerance) or, failing that, to the grid intersection."""
        m = self._nearest_marker_pos(point)
        if m is not None:
            return m
        step = self.grid_step
        return (round(point.x() / step) * step,
                round(point.y() / step) * step)
