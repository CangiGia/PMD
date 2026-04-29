"""CanvasView — interactive QGraphicsView with pan / zoom / Y-flip."""

from __future__ import annotations

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QPainter, QTransform
from PySide6.QtWidgets import QGraphicsView


class CanvasView(QGraphicsView):
    """QGraphicsView with engineering Y-axis (up) and mouse-wheel zoom.

    Signals
    -------
    cursor_moved(x, y) : float, float
        Emitted whenever the mouse moves over the scene; coordinates
        are expressed in scene units (meters).
    """

    cursor_moved = Signal(float, float)

    _ZOOM_IN_FACTOR  = 1.20
    _ZOOM_OUT_FACTOR = 1.0 / 1.20
    _MIN_SCALE       = 5.0       # px / m
    _MAX_SCALE       = 5_000.0   # px / m

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setRenderHint(QPainter.Antialiasing)
        self.setRenderHint(QPainter.SmoothPixmapTransform)
        self.setDragMode(QGraphicsView.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setMouseTracking(True)
        # Make sure the view can receive key events (e.g. Esc to cancel
        # the active tool). Without StrongFocus the scene never sees the
        # key press unless something inside it already has focus.
        self.setFocusPolicy(Qt.StrongFocus)

        # Engineering convention: +y points UP.
        # Default scale: 400 px / m (≈ 50 mm grid step shown ~20 px).
        self.setTransform(QTransform().scale(400, -400))

        self._panning = False
        self._pan_start = QPointF()

    # ──────────────────────────────────────────────────────────
    # Wheel zoom
    # ──────────────────────────────────────────────────────────

    def wheelEvent(self, event):
        factor = (self._ZOOM_IN_FACTOR
                  if event.angleDelta().y() > 0
                  else self._ZOOM_OUT_FACTOR)
        # Clamp scale
        current = abs(self.transform().m11())
        new = current * factor
        if new < self._MIN_SCALE or new > self._MAX_SCALE:
            return
        self.scale(factor, factor)
        self.scene().update()  # redraw adaptive grid

    # ──────────────────────────────────────────────────────────
    # Middle-mouse pan
    # ──────────────────────────────────────────────────────────

    def mousePressEvent(self, event):
        # Take keyboard focus so Esc (and other tool shortcuts) reach
        # the active tool via the scene's keyPressEvent dispatch.
        self.setFocus(Qt.MouseFocusReason)
        if event.button() == Qt.MiddleButton:
            self._panning = True
            self._pan_start = event.position()
            self.setCursor(Qt.ClosedHandCursor)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._panning:
            delta = event.position() - self._pan_start
            self._pan_start = event.position()
            h = self.horizontalScrollBar()
            v = self.verticalScrollBar()
            h.setValue(h.value() - int(delta.x()))
            v.setValue(v.value() - int(delta.y()))
            event.accept()
        else:
            scene_pt = self.mapToScene(event.position().toPoint())
            self.cursor_moved.emit(scene_pt.x(), scene_pt.y())
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MiddleButton and self._panning:
            self._panning = False
            self.unsetCursor()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    # ──────────────────────────────────────────────────────────
    # Convenience
    # ──────────────────────────────────────────────────────────

    def zoom_to_fit(self):
        """Fit all scene items in the view, with a 10 % margin."""
        bounds = self.scene().itemsBoundingRect()
        if bounds.isEmpty():
            bounds = self.scene().sceneRect()
        bounds = bounds.adjusted(-0.1 * bounds.width(),
                                  -0.1 * bounds.height(),
                                   0.1 * bounds.width(),
                                   0.1 * bounds.height())
        self.fitInView(bounds, Qt.KeepAspectRatio)
        self.scene().update()
