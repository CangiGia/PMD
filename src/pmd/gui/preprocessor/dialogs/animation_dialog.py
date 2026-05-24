"""AnimationDialog \u2014 replay a solved run on a preprocessor-style canvas.

Unlike the postprocessor's :class:`AnimationCanvas` (which falls back
to plain coloured circles when bodies have no shape attached at the
runtime level), this dialog rebuilds the *actual* preprocessor scene
\u2014 same body shapes (rect / link / circle), same marker glyphs, same
joint discs \u2014 and animates it by calling
``BodyItem.setPos / setRotation`` per frame from the solver's
``_result_container``. Markers (parented to bodies) and joints
(parented to their *i*-marker) follow automatically.
"""

from __future__ import annotations

import math
import time
from typing import Iterable

from PySide6.QtCore import Qt, QSize, QTimer, Signal
from PySide6.QtGui import QColor, QKeySequence, QPainter, QShortcut, QTransform
from PySide6.QtWidgets import (
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..models import ModelSpec
from ..widgets import BodyItem, JointItem, MarkerItem
from ... import icons as _icons


class _AnimationGraphicsView(QGraphicsView):
    """QGraphicsView with mouse-wheel zoom and middle-button pan.

    Left-button drag pans (``ScrollHandDrag``) and the wheel zooms
    around the cursor. Pan / zoom remain available *during* playback so
    the user can inspect the model while the animation runs.
    """

    _MIN_SCALE = 1.0      # px / m
    _MAX_SCALE = 1.0e5    # px / m
    _ZOOM_STEP = 1.15

    # Emitted on the first user pan or zoom; the parent canvas uses
    # this to disable the auto-refit that fires on showEvent /
    # resizeEvent (we don't want to wipe the user's framing).
    user_interacted = Signal()

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorViewCenter)
        self.setDragMode(QGraphicsView.ScrollHandDrag)

    def wheelEvent(self, event):  # noqa: D401 - Qt override
        delta = event.angleDelta().y()
        if delta == 0:
            return
        factor = self._ZOOM_STEP if delta > 0 else 1.0 / self._ZOOM_STEP
        cur = abs(self.transform().m11())
        new = cur * factor
        if new < self._MIN_SCALE or new > self._MAX_SCALE:
            return
        # Preserve the y-down flip while scaling.
        sign_y = -1.0 if self.transform().m22() < 0 else 1.0
        # Compose: scale around the cursor.
        # Easiest path: use scale() which respects the transformation
        # anchor we set in __init__.
        self.scale(factor, factor)
        # Re-apply the y-flip if it was lost (scale() preserves sign).
        _ = sign_y
        self.user_interacted.emit()
        event.accept()

    def mousePressEvent(self, event):  # noqa: D401 - Qt override
        # Any mouse press on the view is treated as a user interaction
        # (the ScrollHandDrag mode pans on left-button drag).
        self.user_interacted.emit()
        super().mousePressEvent(event)


class PreprocessorAnimationCanvas(QWidget):
    """Reusable widget that replays ``(T, uT)`` on a preprocessor scene.

    Used both inside :class:`AnimationDialog` (preprocessor "Animate"
    button) and inside the postprocessor's animation pane when the
    session was launched from the preprocessor.

    Parameters
    ----------
    spec : ModelSpec
        The preprocessor model spec used to rebuild body / marker /
        joint visuals (shapes, names, layout).
    model : PlanarMultibodyModel
        Solved runtime model. Each body's ``_result_container`` is read
        per-frame to drive the visuals.
    T : ndarray
        Time vector ``(n_steps,)``.
    parent : QWidget or None
        Optional parent widget.
    """

    # Emitted whenever the displayed frame changes; carries the
    # current simulation time (s). Listeners (e.g. the postprocessor's
    # PlotCanvas) can use it to draw a synchronised time cursor.
    time_changed = Signal(float)

    def __init__(self, spec: ModelSpec, model, T, parent=None):
        super().__init__(parent)

        self._spec    = spec
        self._model   = model
        self._T       = T
        self._n       = int(len(T))
        self._step    = 0
        self._playing = False

        # ---- map runtime bodies by name (== spec.name or spec.id) ----
        # The builder uses ``name=b.name or b.id`` so we can recover
        # the spec\u2192runtime mapping unambiguously.
        self._body_by_specid: dict[str, object] = {}
        for b_spec in spec.bodies:
            if b_spec.id == "ground":
                continue
            target = b_spec.name or b_spec.id
            for body in getattr(model, "Bodies", []):
                if getattr(body, "name", None) == target:
                    self._body_by_specid[b_spec.id] = body
                    break

        # ---- scene / view (mirrors the preprocessor canvas) ----
        self._scene = QGraphicsScene(self)
        self._scene.setBackgroundBrush(QColor("#fafafa"))
        self._view = _AnimationGraphicsView(self._scene, self)
        self._view.setRenderHints(
            QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        # Same y-up convention as the preprocessor canvas.
        self._view.setTransform(QTransform().scale(400, -400))

        self._body_items:   dict[str, BodyItem]   = {}
        self._marker_items: dict[str, MarkerItem] = {}
        self._joint_items:  dict[str, JointItem]  = {}
        self._build_scene()

        # Freeze the scene rect to the full *trajectory* envelope so
        # that bodies moving across the canvas during playback never
        # cause the view to re-fit / rescale itself (which produced
        # the constantly-shifting axes the user complained about).
        self._freeze_scene_rect()

        # Fit once after the scene is populated.
        QTimer.singleShot(0, self._fit_view)
        # Track whether the user has manually pan/zoomed since the last
        # auto-fit. While untouched, the canvas re-fits whenever it
        # actually receives a usable size (showEvent / resizeEvent),
        # which is what fixes the "tiny model in a huge pane" you see
        # when the postprocessor first toggles the animation pane on.
        self._user_zoomed = False
        self._view.user_interacted.connect(self._on_user_interacted)

        # ---- transport controls ----
        # Larger, modern icon buttons (24 px) with mdi6 glyphs.
        _BTN_SIZE = QSize(40, 40)
        _ICN_SIZE = QSize(24, 24)

        def _make_transport(icon_name: str, tip: str, slot) -> QToolButton:
            b = QToolButton()
            b.setIcon(_icons.icon(icon_name))
            b.setIconSize(_ICN_SIZE)
            b.setFixedSize(_BTN_SIZE)
            b.setAutoRaise(True)
            b.setCursor(Qt.PointingHandCursor)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            return b

        self._prev_btn = _make_transport(
            "mdi6.skip-previous", "Previous frame",
            lambda: self._goto(self._step - 1))
        self._play_btn = _make_transport(
            "mdi6.play", "Play / Pause (Space)", self._toggle_play)
        self._next_btn = _make_transport(
            "mdi6.skip-next", "Next frame",
            lambda: self._goto(self._step + 1))

        self._slider = QSlider(Qt.Horizontal)
        self._slider.setRange(0, max(0, self._n - 1))
        self._slider.valueChanged.connect(self._on_slider)

        # Playback speed multiplier (× real-time). Decoupled from the
        # display refresh rate: we always render at ~60 Hz and *skip*
        # frames as needed to stay locked to wall-clock time. This
        # makes "Speed" actually behave like a speed knob, even when
        # the simulation has thousands of steps over a few seconds
        # (where any sensible frame-per-tick scheme bottoms out at the
        # display refresh rate and looks frozen at high values).
        self._fps_spin = QDoubleSpinBox()
        self._fps_spin.setDecimals(2)
        self._fps_spin.setRange(0.05, 100.0)
        self._fps_spin.setSingleStep(0.25)
        self._fps_spin.setValue(1.0)
        self._fps_spin.setSuffix("\u00d7")
        self._fps_spin.setMinimumWidth(96)
        self._fps_spin.setToolTip(
            "Playback speed multiplier (\u00d7 real-time).\n"
            "1\u00d7 plays the simulation at its physical wall-clock rate.")
        self._fps_spin.valueChanged.connect(self._on_speed_changed)

        self._time_lbl = QLabel(self._fmt_time(0))

        self._fit_btn = QPushButton("Fit")
        self._fit_btn.clicked.connect(self._fit_view)

        # Marker visibility toggle.  Default ON.  Same UX as the
        # postprocessor's AnimationCanvas: "V" shortcut mirrors Adams.
        self._markers_btn = QToolButton()
        self._markers_btn.setIcon(_icons.icon("mdi6.eye-outline"))
        self._markers_btn.setIconSize(_ICN_SIZE)
        self._markers_btn.setFixedSize(_BTN_SIZE)
        self._markers_btn.setAutoRaise(True)
        self._markers_btn.setCursor(Qt.PointingHandCursor)
        self._markers_btn.setToolTip("Show / hide markers (V)")
        self._markers_btn.setCheckable(True)
        self._markers_btn.setChecked(True)
        self._markers_btn.toggled.connect(self._on_toggle_markers)

        self._export_btn = QPushButton("Export\u2026")
        self._export_btn.setToolTip(
            "Export the animation as a sequence of PNG frames")
        self._export_btn.clicked.connect(self._on_export)

        ctrl = QHBoxLayout()
        for w in (self._prev_btn, self._play_btn, self._next_btn):
            ctrl.addWidget(w)
        ctrl.addWidget(self._slider, stretch=1)
        ctrl.addWidget(self._time_lbl)
        ctrl.addSpacing(8)
        ctrl.addWidget(QLabel("Speed:"))
        ctrl.addWidget(self._fps_spin)
        ctrl.addSpacing(8)
        ctrl.addWidget(self._fit_btn)
        ctrl.addWidget(self._markers_btn)
        ctrl.addWidget(self._export_btn)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self._view, stretch=1)
        lay.addLayout(ctrl)

        # ---- playback timer ----
        self._timer = QTimer(self)
        self._timer.setTimerType(Qt.PreciseTimer)
        self._timer.timeout.connect(self._advance)
        # Wall-clock anchors (set on each play / speed change).
        self._wall_t0: float = 0.0
        self._sim_t0: float = 0.0
        self._update_timer_interval()

        # "V" toggles marker visibility (mirrors Adams / the
        # postprocessor's AnimationCanvas).  Window-scope so the
        # shortcut works regardless of which child currently has
        # keyboard focus.
        self._markers_shortcut = QShortcut(QKeySequence(Qt.Key_V), self)
        self._markers_shortcut.setContext(Qt.WindowShortcut)
        self._markers_shortcut.activated.connect(self._markers_btn.click)

        # First frame.
        self._apply_step(0)
    # ------------------------------------------------------------------
    # Scene build
    # ------------------------------------------------------------------
    def _build_scene(self) -> None:
        # Bodies, plus the same MarkerItem glyphs used on the
        # preprocessor's editing canvas so position *and* orientation
        # of each frame are legible during playback.  Visibility is
        # user-controlled (toolbar button or shortcut "V") -- toggle
        # off for a clean "bodies only" view.
        for b in self._spec.bodies:
            if b.id == "ground":
                continue
            item = BodyItem(b)
            # Make non-interactive in the playback scene.
            item.setFlag(QGraphicsItem.ItemIsMovable, False)
            item.setFlag(QGraphicsItem.ItemIsSelectable, False)
            self._scene.addItem(item)
            self._body_items[b.id] = item

        # Markers: parent to their owning BodyItem so they follow the
        # body automatically through setPos / setRotation.  Ground
        # markers (if any) are added at scene level.
        for m in self._spec.markers:
            parent_item = self._body_items.get(m.body_id) if m.body_id != "ground" else None
            mk_item = MarkerItem(m, parent_body=parent_item)
            mk_item.setFlag(QGraphicsItem.ItemIsSelectable, False)
            if parent_item is None:
                self._scene.addItem(mk_item)
            self._marker_items[m.id] = mk_item

    # ------------------------------------------------------------------
    # Auto-fit lifecycle
    # ------------------------------------------------------------------
    def _on_user_interacted(self) -> None:
        self._user_zoomed = True

    def showEvent(self, event):  # noqa: D401 - Qt override
        super().showEvent(event)
        # Re-fit on first show (and on subsequent shows while the user
        # hasn't manually framed the view). Without this, when the
        # postprocessor's animation pane is hidden at construction
        # time the initial QTimer.singleShot fit runs against a 0x0
        # viewport and the model ends up as a single pixel.
        if not self._user_zoomed:
            QTimer.singleShot(0, self._fit_view)

    def resizeEvent(self, event):  # noqa: D401 - Qt override
        super().resizeEvent(event)
        if not self._user_zoomed:
            # Defer to next event loop tick so the QGraphicsView has
            # already adopted its new viewport size before we fit.
            QTimer.singleShot(0, self._fit_view)

    def _fit_view(self) -> None:
        rect = self._scene.sceneRect()
        if rect.isEmpty():
            rect = self._scene.itemsBoundingRect()
        if rect.isEmpty():
            return
        self._view.fitInView(rect, Qt.KeepAspectRatio)
        # fitInView would reset our flipped Y; reapply the y-down flip
        # while preserving the computed scale.
        s = abs(self._view.transform().m11())
        self._view.setTransform(QTransform().scale(s, -s))
        self._view.centerOn(rect.center())

    def _freeze_scene_rect(self) -> None:
        """Compute the trajectory envelope and lock the scene rect to it.

        With a frozen scene rect the view never auto-rescales as the
        bodies move, eliminating the "axes constantly resizing" effect
        observed during playback.
        """
        # Start from a rough envelope of all bodies' positions across
        # the entire trajectory, then pad by the largest body's
        # bounding rect so geometry never clips the canvas edge.
        xs: list[float] = []
        ys: list[float] = []
        for spec_id, body in self._body_by_specid.items():
            try:
                rc = body._result_container
                xs.extend(float(v) for v in rc["positions"]["x"])
                ys.extend(float(v) for v in rc["positions"]["y"])
            except Exception:
                continue

        # Body extent (largest body half-diagonal, in metres).
        extent = 0.0
        for item in self._body_items.values():
            br = item.boundingRect()
            extent = max(extent,
                         math.hypot(br.width(), br.height()) * 0.5)

        if not xs or not ys:
            # Fall back to whatever the items currently span.
            r = self._scene.itemsBoundingRect()
            if r.isEmpty():
                r.setRect(-0.5, -0.5, 1.0, 1.0)
            pad = 0.10 * max(r.width(), r.height(), 0.1)
            r.adjust(-pad, -pad, pad, pad)
            self._scene.setSceneRect(r)
            return

        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)
        span = max(x_max - x_min, y_max - y_min, 0.1)
        pad = 0.15 * span + 1.5 * extent
        from PySide6.QtCore import QRectF
        rect = QRectF(x_min - pad, y_min - pad,
                      (x_max - x_min) + 2 * pad,
                      (y_max - y_min) + 2 * pad)
        self._scene.setSceneRect(rect)

    # ------------------------------------------------------------------
    # Transport
    # ------------------------------------------------------------------
    def _toggle_play(self) -> None:
        if self._playing:
            self._timer.stop()
            self._play_btn.setIcon(_icons.icon("mdi6.play"))
        else:
            if self._step >= self._n - 1:
                self._step = 0
            # Re-anchor wall clock so playback resumes from the
            # current step instead of jumping forward.
            self._anchor_playback()
            self._timer.start()
            self._play_btn.setIcon(_icons.icon("mdi6.pause"))
        self._playing = not self._playing

    def _anchor_playback(self) -> None:
        """Pin (wall_now, sim_now) so wall-clock advancing resumes here."""
        self._wall_t0 = time.monotonic()
        try:
            self._sim_t0 = float(self._T[self._step])
        except Exception:
            self._sim_t0 = 0.0

    def _on_speed_changed(self, _val: float) -> None:
        # Re-anchor so the speed change takes effect from "now" without
        # a discontinuous jump in displayed time.
        if self._playing:
            self._anchor_playback()

    def _update_timer_interval(self) -> None:
        # Fixed display refresh; "Speed" is decoupled (see _advance).
        # ~60 Hz is plenty for smooth perceived motion and keeps CPU
        # usage in check on slow machines.
        self._timer.setInterval(16)

    def _advance(self) -> None:
        # Wall-clock driven: figure out which step corresponds to the
        # elapsed wall time at the current playback speed, snap to it,
        # and stop at the end of the trajectory.
        speed = max(1e-6, float(self._fps_spin.value()))
        elapsed = time.monotonic() - self._wall_t0
        sim_t = self._sim_t0 + elapsed * speed
        # Map sim_t -> nearest step index (assumes T is monotonic).
        try:
            t_end = float(self._T[-1])
        except Exception:
            t_end = 0.0
        if sim_t >= t_end:
            self._goto(self._n - 1)
            self._timer.stop()
            self._playing = False
            self._play_btn.setIcon(_icons.icon("mdi6.play"))
            return
        # Linear search from current step (cheap; T is monotonic and we
        # never go backwards within _advance).
        s = self._step
        n = self._n
        while s + 1 < n and float(self._T[s + 1]) <= sim_t:
            s += 1
        if s != self._step:
            self._goto(s)

    def _goto(self, step: int) -> None:
        step = max(0, min(self._n - 1, int(step)))
        if step == self._step:
            return
        self._slider.blockSignals(True)
        self._slider.setValue(step)
        self._slider.blockSignals(False)
        self._apply_step(step)

    def _on_slider(self, v: int) -> None:
        self._apply_step(int(v))

    def set_step(self, step: int) -> None:
        """Public hook (used by the postprocessor's plot canvas)."""
        self._goto(step)

    def _on_toggle_markers(self, checked: bool) -> None:
        """Show / hide all marker glyphs (shortcut: V)."""
        for mk in self._marker_items.values():
            mk.setVisible(checked)
        self._markers_btn.setIcon(_icons.icon(
            "mdi6.eye-outline" if checked else "mdi6.eye-off-outline"))

    def _apply_step(self, step: int) -> None:
        self._step = step
        for spec_id, body in self._body_by_specid.items():
            item = self._body_items.get(spec_id)
            if item is None:
                continue
            try:
                rc = body._result_container
                x   = float(rc["positions"]["x"][step])
                y   = float(rc["positions"]["y"][step])
                phi = float(rc["positions"]["phi"][step])
            except Exception:
                continue
            item.setPos(x, y)
            item.setRotation(math.degrees(phi))
        self._time_lbl.setText(self._fmt_time(step))
        try:
            self.time_changed.emit(float(self._T[step]))
        except Exception:
            pass

    def _fmt_time(self, step: int) -> str:
        try:
            t = float(self._T[step])
        except Exception:
            t = 0.0
        return f"t = {t:.4f} s   ({step + 1}/{self._n})"

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    def _on_export(self) -> None:
        was_playing = self._playing
        if was_playing:
            self._toggle_play()

        d = QFileDialog.getExistingDirectory(
            self, "Choose a folder for the PNG frames")
        if not d:
            return

        from pathlib import Path
        out = Path(d)
        out.mkdir(parents=True, exist_ok=True)
        # Render every frame to a PNG.
        from PySide6.QtGui import QImage
        size = self._view.viewport().size()
        for i in range(self._n):
            self._apply_step(i)
            img = QImage(size, QImage.Format_ARGB32)
            img.fill(Qt.white)
            painter = QPainter(img)
            painter.setRenderHint(QPainter.Antialiasing, True)
            self._view.render(painter)
            painter.end()
            img.save(str(out / f"frame_{i:05d}.png"))
        # Reset to the first frame at the end.
        self._goto(0)

    # ------------------------------------------------------------------
    # Keyboard shortcuts
    # ------------------------------------------------------------------
    def keyPressEvent(self, event):
        k = event.key()
        if k == Qt.Key_Space:
            self._toggle_play(); event.accept(); return
        if k in (Qt.Key_Right, Qt.Key_Period):
            self._goto(self._step + 1); event.accept(); return
        if k in (Qt.Key_Left, Qt.Key_Comma):
            self._goto(self._step - 1); event.accept(); return
        if k == Qt.Key_Home:
            self._goto(0); event.accept(); return
        if k == Qt.Key_End:
            self._goto(self._n - 1); event.accept(); return
        super().keyPressEvent(event)


class AnimationDialog(QDialog):
    """Modal-less dialog wrapping :class:`PreprocessorAnimationCanvas`."""

    def __init__(self, spec: ModelSpec, model, T, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Animation \u2014 {spec.name}")
        self.resize(1000, 760)
        self._canvas = PreprocessorAnimationCanvas(spec, model, T, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self._canvas)


__all__ = ["AnimationDialog", "PreprocessorAnimationCanvas"]
