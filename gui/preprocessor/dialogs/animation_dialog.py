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
from typing import Iterable

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QPainter, QTransform
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
        self._view = QGraphicsView(self._scene, self)
        self._view.setRenderHints(
            QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        # Same y-up convention as the preprocessor canvas.
        self._view.setTransform(QTransform().scale(400, -400))
        self._view.setDragMode(QGraphicsView.ScrollHandDrag)

        self._body_items:   dict[str, BodyItem]   = {}
        self._marker_items: dict[str, MarkerItem] = {}
        self._joint_items:  dict[str, JointItem]  = {}
        self._build_scene()

        # Fit once after the scene is populated.
        QTimer.singleShot(0, self._fit_view)

        # ---- transport controls ----
        self._play_btn = QToolButton()
        self._play_btn.setText("\u25b6")     # play
        self._play_btn.setToolTip("Play / Pause (Space)")
        self._play_btn.clicked.connect(self._toggle_play)

        self._prev_btn = QToolButton()
        self._prev_btn.setText("\u23ee")
        self._prev_btn.setToolTip("Previous frame")
        self._prev_btn.clicked.connect(lambda: self._goto(self._step - 1))

        self._next_btn = QToolButton()
        self._next_btn.setText("\u23ed")
        self._next_btn.setToolTip("Next frame")
        self._next_btn.clicked.connect(lambda: self._goto(self._step + 1))

        self._slider = QSlider(Qt.Horizontal)
        self._slider.setRange(0, max(0, self._n - 1))
        self._slider.valueChanged.connect(self._on_slider)

        self._fps_spin = QDoubleSpinBox()
        self._fps_spin.setRange(1.0, 120.0)
        self._fps_spin.setValue(25.0)
        self._fps_spin.setSuffix(" fps")
        self._fps_spin.valueChanged.connect(self._update_timer_interval)

        self._time_lbl = QLabel(self._fmt_time(0))

        self._fit_btn = QPushButton("Fit")
        self._fit_btn.clicked.connect(self._fit_view)

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
        ctrl.addWidget(self._export_btn)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self._view, stretch=1)
        lay.addLayout(ctrl)

        # ---- playback timer ----
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._advance)
        self._update_timer_interval()

        # First frame.
        self._apply_step(0)
    # ------------------------------------------------------------------
    # Scene build
    # ------------------------------------------------------------------
    def _build_scene(self) -> None:
        # Bodies
        for b in self._spec.bodies:
            if b.id == "ground":
                continue
            item = BodyItem(b)
            # Make non-interactive in the playback scene.
            item.setFlag(QGraphicsItem.ItemIsMovable, False)
            item.setFlag(QGraphicsItem.ItemIsSelectable, False)
            self._scene.addItem(item)
            self._body_items[b.id] = item

        # Markers (parented to body so they follow it; ground markers
        # go straight into the scene).
        for m in self._spec.markers:
            parent = self._body_items.get(m.body_id)
            if parent is None:
                item = MarkerItem(m)
                item.setPos(*m.local_position)
                self._scene.addItem(item)
            else:
                item = MarkerItem(m, parent_body=parent)
            item.setFlag(QGraphicsItem.ItemIsSelectable, False)
            self._marker_items[m.id] = item

        # Joints \u2014 parent to the i-marker so the glyph follows the
        # owning body without per-frame work. The JointItem's local
        # origin already sits at the marker, so no extra setPos.
        for j in self._spec.joints:
            i_marker = self._marker_items.get(j.i_marker_id)
            if i_marker is None:
                continue
            item = JointItem(j)
            item.setParentItem(i_marker)
            item.setFlag(QGraphicsItem.ItemIsSelectable, False)
            self._joint_items[j.id] = item

    def _fit_view(self) -> None:
        rect = self._scene.itemsBoundingRect()
        if rect.isEmpty():
            return
        # Pad a bit so bodies don't sit on the edge.
        pad = 0.10 * max(rect.width(), rect.height(), 0.1)
        rect.adjust(-pad, -pad, pad, pad)
        self._view.fitInView(rect, Qt.KeepAspectRatio)
        # fitInView would reset our flipped Y; reapply the y-down flip
        # while preserving the computed scale.
        s = abs(self._view.transform().m11())
        self._view.setTransform(QTransform().scale(s, -s))
        self._view.centerOn(rect.center())

    # ------------------------------------------------------------------
    # Transport
    # ------------------------------------------------------------------
    def _toggle_play(self) -> None:
        if self._playing:
            self._timer.stop()
            self._play_btn.setText("\u25b6")
        else:
            if self._step >= self._n - 1:
                self._step = 0
            self._timer.start()
            self._play_btn.setText("\u23f8")
        self._playing = not self._playing

    def _update_timer_interval(self) -> None:
        fps = max(1.0, float(self._fps_spin.value()))
        self._timer.setInterval(int(round(1000.0 / fps)))

    def _advance(self) -> None:
        nxt = self._step + 1
        if nxt >= self._n:
            self._timer.stop()
            self._playing = False
            self._play_btn.setText("\u25b6")
            return
        self._goto(nxt)

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
