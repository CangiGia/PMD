"""AnimationCanvas — 2-D model animation with shape rendering."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.patches import Circle as MplCircle, FancyBboxPatch, Polygon as MplPolygon
from matplotlib.transforms import Affine2D

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QToolBar,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pmd.core.shapes import Rectangle, Circle, Plate, Link
from pmd.core.constraints import RevJoint, TranJoint, PtpForce
from pmd.core.mechanics import rotation_matrix
from ... import icons as _icons
from ...style import CANVAS_BG_DARK, CANVAS_BG_LIGHT, CANVAS_FG_DARK, CANVAS_FG_LIGHT

# Colour palette for bodies (tab10)
_BODY_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def _body_pos(body, step):
    """Return (x, y, phi) of *body* at *step*."""
    rc = body._result_container
    return (
        float(rc["positions"]["x"][step]),
        float(rc["positions"]["y"][step]),
        float(rc["positions"]["phi"][step]),
    )


def _marker_global(marker, step):
    """Return global (x, y) ndarray of *marker* at *step*."""
    body = marker.body
    if not body:                         # Ground
        return marker.local_position.ravel()[:2].copy()
    x, y, phi = _body_pos(body, step)
    return np.array([x, y]) + rotation_matrix(phi) @ marker.local_position.ravel()[:2]


class AnimationCanvas(QWidget):
    """2-D animation of the multi-body model at each time step.

    Parameters
    ----------
    sessions : list[Session]
        Solved simulation sessions.
    parent : QWidget or None
        Optional parent widget.
    """

    time_changed = Signal(int)

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        self._sessions = sessions
        self._T = sessions[0].T
        self._n_steps = len(self._T)
        self._step = 0
        self._playing = False

        # --- matplotlib widgets ---
        self._figure = Figure(tight_layout=True)
        self._ax = self._figure.add_subplot(111)
        self._ax.set_aspect("equal")
        self._canvas = FigureCanvasQTAgg(self._figure)

        # Backend nav toolbar (hidden) — used only for navigate state/history
        self._nav = NavigationToolbar2QT(self._canvas, self)
        self._nav.hide()

        # Custom visible toolbar (built after instance vars)
        self._toolbar = self._build_toolbar()

        # --- transport controls ---
        self._play_btn = QPushButton()
        self._play_btn.setIcon(_icons.icon("mdi6.play"))
        self._play_btn.setToolTip("Play / Pause")
        self._play_btn.setFixedWidth(36)
        self._play_btn.clicked.connect(self._on_play_pause)

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setRange(0, self._n_steps - 1)
        self._slider.setValue(0)
        self._slider.valueChanged.connect(self._on_slider_changed)

        self._time_lbl = QLabel(self._time_text(0))

        ctrl = QHBoxLayout()
        ctrl.addWidget(self._play_btn)
        ctrl.addWidget(self._slider, stretch=1)
        ctrl.addWidget(self._time_lbl)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)
        layout.addLayout(ctrl)

        # --- artist storage ---
        self._body_patches = []
        self._body_info = []       # list of (body, session_idx)
        self._marker_dots = []
        self._marker_info = []     # (marker, session_idx)
        self._joint_markers = []
        self._joint_info = []      # (joint, session_idx)
        self._force_lines = []
        self._force_info = []      # (force, session_idx)

        self._init_artists()

        # --- animation timer ---
        self._timer = self._canvas.new_timer(interval=30)
        self._timer.add_callback(self._advance_frame)

    # ------------------------------------------------------------------
    # Custom toolbar
    # ------------------------------------------------------------------

    def _build_toolbar(self) -> QToolBar:
        tb = QToolBar()
        tb.setMovable(False)
        self._btn_home = self._make_btn(tb, "Home",    "mdi6.home",                    self._on_home)
        self._btn_back = self._make_btn(tb, "Back",    "mdi6.arrow-left",              self._on_back)
        self._btn_fwd  = self._make_btn(tb, "Forward", "mdi6.arrow-right",             self._on_fwd)
        tb.addSeparator()
        self._btn_pan  = self._make_btn(tb, "Pan",     "mdi6.hand-back-right-outline", self._on_pan,  checkable=True)
        self._btn_zoom = self._make_btn(tb, "Zoom",    "mdi6.magnify",                 self._on_zoom, checkable=True)
        tb.addSeparator()
        self._btn_save = self._make_btn(tb, "Save",    "mdi6.content-save",            self._on_save)
        return tb

    @staticmethod
    def _make_btn(tb: QToolBar, text: str, icon_name: str, slot,
                  *, checkable: bool = False) -> QToolButton:
        btn = QToolButton()
        btn.setIcon(_icons.icon(icon_name))
        btn.setToolTip(text)
        btn.setCheckable(checkable)
        btn.clicked.connect(slot)
        tb.addWidget(btn)
        return btn

    def set_icon_theme(self, dark: bool) -> None:
        """Re-apply icon colours after a theme toggle."""
        self._btn_home.setIcon(_icons.icon("mdi6.home"))
        self._btn_back.setIcon(_icons.icon("mdi6.arrow-left"))
        self._btn_fwd.setIcon( _icons.icon("mdi6.arrow-right"))
        self._btn_pan.setIcon( _icons.icon("mdi6.hand-back-right-outline"))
        self._btn_zoom.setIcon(_icons.icon("mdi6.magnify"))
        self._btn_save.setIcon(_icons.icon("mdi6.content-save"))
        play_icon = "mdi6.pause" if self._playing else "mdi6.play"
        self._play_btn.setIcon(_icons.icon(play_icon))

    def set_dark(self, enabled: bool) -> None:
        """Switch the figure between dark and light appearance."""
        bg = CANVAS_BG_DARK  if enabled else CANVAS_BG_LIGHT
        fg = CANVAS_FG_DARK  if enabled else CANVAS_FG_LIGHT
        self._figure.set_facecolor(bg)
        self._ax.set_facecolor(bg)
        self._ax.tick_params(colors=fg, which="both")
        self._ax.xaxis.label.set_color(fg)
        self._ax.yaxis.label.set_color(fg)
        self._ax.title.set_color(fg)
        for spine in self._ax.spines.values():
            spine.set_edgecolor(fg)
        self._canvas.draw_idle()

    def _on_home(self):
        self._nav.home()
        self._canvas.draw_idle()

    def _on_back(self):
        self._nav.back()
        self._canvas.draw_idle()

    def _on_fwd(self):
        self._nav.forward()
        self._canvas.draw_idle()

    def _on_pan(self, checked: bool):
        if checked:
            self._btn_zoom.setChecked(False)
        self._nav.pan()

    def _on_zoom(self, checked: bool):
        if checked:
            self._btn_pan.setChecked(False)
        self._nav.zoom()

    def _on_save(self):
        self._nav.save_figure()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _time_text(self, step):
        return f"t = {self._T[step]:.4f} s"

    # ------------------------------------------------------------------
    # _init_artists  — draw everything at frame 0
    # ------------------------------------------------------------------

    def _init_artists(self):
        ax = self._ax
        color_idx = 0

        for si, session in enumerate(self._sessions):
            model = session.model

            # ---- Bodies ----
            for body in model.Bodies:
                if not body:
                    continue                 # skip Ground
                col = _BODY_COLORS[color_idx % len(_BODY_COLORS)]
                color_idx += 1

                x, y, phi = _body_pos(body, 0)
                shape = body.shape

                if isinstance(shape, Rectangle):
                    w, h = shape.width, shape.height
                    # FancyBboxPatch anchor is lower-left corner
                    patch = FancyBboxPatch(
                        (-w / 2, -h / 2), w, h,
                        boxstyle="round,pad=0.01",
                        facecolor=col, edgecolor=col,
                        alpha=0.30, linewidth=1.5,
                    )
                    patch.set_edgecolor(col)
                    patch.set_alpha(None)            # let face/edge alphas rule
                    patch.set_facecolor((*matplotlib.colors.to_rgb(col), 0.30))
                    patch.set_edgecolor((*matplotlib.colors.to_rgb(col), 1.0))
                    t = Affine2D().rotate(phi).translate(x, y) + ax.transData
                    patch.set_transform(t)
                    ax.add_patch(patch)
                    self._body_patches.append(patch)

                elif isinstance(shape, Circle):
                    patch = MplCircle(
                        (x, y), shape.radius,
                        facecolor=(*matplotlib.colors.to_rgb(col), 0.30),
                        edgecolor=(*matplotlib.colors.to_rgb(col), 1.0),
                        linewidth=1.5,
                    )
                    ax.add_patch(patch)
                    self._body_patches.append(patch)

                elif isinstance(shape, Plate):
                    verts = shape.vertices.copy()
                    patch = MplPolygon(
                        verts, closed=True,
                        facecolor=(*matplotlib.colors.to_rgb(col), 0.30),
                        edgecolor=(*matplotlib.colors.to_rgb(col), 1.0),
                        linewidth=1.5,
                    )
                    t = Affine2D().rotate(phi).translate(x, y) + ax.transData
                    patch.set_transform(t)
                    ax.add_patch(patch)
                    self._body_patches.append(patch)

                elif isinstance(shape, Link):
                    w, h = shape.length, shape.thickness
                    patch = FancyBboxPatch(
                        (-w / 2, -h / 2), w, h,
                        boxstyle="round,pad=0.002",
                        facecolor=(*matplotlib.colors.to_rgb(col), 0.30),
                        edgecolor=(*matplotlib.colors.to_rgb(col), 1.0),
                        linewidth=1.5,
                    )
                    t = Affine2D().rotate(phi).translate(x, y) + ax.transData
                    patch.set_transform(t)
                    ax.add_patch(patch)
                    self._body_patches.append(patch)

                else:
                    # No shape → small circle at CoM
                    patch = MplCircle(
                        (x, y), 0.05,
                        facecolor=(*matplotlib.colors.to_rgb(col), 0.30),
                        edgecolor=(*matplotlib.colors.to_rgb(col), 1.0),
                        linewidth=1.5,
                    )
                    ax.add_patch(patch)
                    self._body_patches.append(patch)

                self._body_info.append((body, si))

                # markers attached to this body
                for mk in body._markers:
                    gp = _marker_global(mk, 0)
                    dot, = ax.plot(gp[0], gp[1], "k.", markersize=3)
                    self._marker_dots.append(dot)
                    self._marker_info.append((mk, si))

            # ---- Joints ----
            for joint in model.Joints:
                mk = joint.iMarker or joint.jMarker
                if mk is None:
                    continue
                gp = _marker_global(mk, 0)

                if isinstance(joint, RevJoint):
                    dot, = ax.plot(gp[0], gp[1], "o", color="#333333",
                                   markersize=8, markerfacecolor="none",
                                   markeredgewidth=1.5)
                elif isinstance(joint, TranJoint):
                    dot, = ax.plot(gp[0], gp[1], "s", color="#333333",
                                   markersize=8, markerfacecolor="none",
                                   markeredgewidth=1.5)
                else:
                    dot, = ax.plot(gp[0], gp[1], ".", color="#333333",
                                   markersize=6)

                self._joint_markers.append(dot)
                self._joint_info.append((joint, si))

            # ---- Forces ----
            for force in model.Forces:
                if isinstance(force, PtpForce) and force.iMarker and force.jMarker:
                    p1 = _marker_global(force.iMarker, 0)
                    p2 = _marker_global(force.jMarker, 0)
                    line, = ax.plot(
                        [p1[0], p2[0]], [p1[1], p2[1]],
                        "--", color="#555555", linewidth=1.0,
                    )
                    self._force_lines.append(line)
                    self._force_info.append((force, si))

        # ---- axis limits from first + last frame ----
        self._auto_limits()
        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # auto limits
    # ------------------------------------------------------------------

    def _auto_limits(self):
        """Fit the axes to everything that will ever be drawn.

        The previous implementation only sampled body CM positions at
        the first/last frame, which produced a near-degenerate window
        for models whose CMs barely move while the link *endpoints*
        swing widely (typical pendulum case). We now include every
        marker's world position at every timestep — markers sit at
        joint locations and link ends, so their swept envelope is a
        tight, complete bound for the visible geometry.
        """
        xs: list[float] = []
        ys: list[float] = []

        # Body CMs across all frames (covers shapeless bodies).
        for body, _si in self._body_info:
            rc = body._result_container
            xs.extend([float(rc["positions"]["x"].min()),
                       float(rc["positions"]["x"].max())])
            ys.extend([float(rc["positions"]["y"].min()),
                       float(rc["positions"]["y"].max())])

        # Markers across all frames — captures link endpoints and any
        # off-CM features the user attached a frame to.
        for mk, _si in self._marker_info:
            body = mk.body
            lp = mk.local_position.ravel()[:2]
            if not body:
                xs.append(float(lp[0])); ys.append(float(lp[1]))
                continue
            rc = body._result_container
            bx = np.asarray(rc["positions"]["x"], dtype=float)
            by = np.asarray(rc["positions"]["y"], dtype=float)
            ph = np.asarray(rc["positions"]["phi"], dtype=float)
            cos_p = np.cos(ph); sin_p = np.sin(ph)
            mx_t = bx + cos_p * lp[0] - sin_p * lp[1]
            my_t = by + sin_p * lp[0] + cos_p * lp[1]
            xs.extend([float(mx_t.min()), float(mx_t.max())])
            ys.extend([float(my_t.min()), float(my_t.max())])

        if not xs:
            return
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)
        # 10% margin on the larger span; keep a sensible floor so a
        # static model doesn't collapse to a zero-size window.
        span = max((x_max - x_min), (y_max - y_min), 0.1)
        mx = span * 0.10
        self._ax.set_xlim(x_min - mx, x_max + mx)
        self._ax.set_ylim(y_min - mx, y_max + mx)

    # ------------------------------------------------------------------
    # _update_artists — move everything to *step*
    # ------------------------------------------------------------------

    def _update_artists(self, step):
        ax = self._ax

        # bodies
        for patch, (body, _si) in zip(self._body_patches, self._body_info):
            x, y, phi = _body_pos(body, step)
            shape = body.shape

            if isinstance(shape, (Rectangle, Polygon, Link)):
                t = Affine2D().rotate(phi).translate(x, y) + ax.transData
                patch.set_transform(t)
            elif isinstance(shape, Circle):
                patch.set_center((x, y))
            else:
                patch.set_center((x, y))

        # marker dots
        for dot, (mk, _si) in zip(self._marker_dots, self._marker_info):
            gp = _marker_global(mk, step)
            dot.set_data([gp[0]], [gp[1]])

        # joints
        for dot, (joint, _si) in zip(self._joint_markers, self._joint_info):
            mk = joint.iMarker or joint.jMarker
            if mk is not None:
                gp = _marker_global(mk, step)
                dot.set_data([gp[0]], [gp[1]])

        # forces
        for line, (force, _si) in zip(self._force_lines, self._force_info):
            p1 = _marker_global(force.iMarker, step)
            p2 = _marker_global(force.jMarker, step)
            line.set_data([p1[0], p2[0]], [p1[1], p2[1]])

        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def set_step(self, step):
        """Jump to the given time step and refresh all artists."""
        self._step = step
        self._slider.blockSignals(True)
        self._slider.setValue(step)
        self._slider.blockSignals(False)
        self._time_lbl.setText(self._time_text(step))
        self._update_artists(step)
        self.time_changed.emit(step)

    # ------------------------------------------------------------------
    # slots
    # ------------------------------------------------------------------

    def _on_slider_changed(self, value):
        self.set_step(value)

    def _on_play_pause(self):
        if self._playing:
            self._timer.stop()
            self._play_btn.setIcon(_icons.icon("mdi6.play"))
        else:
            self._timer.start()
            self._play_btn.setIcon(_icons.icon("mdi6.pause"))
        self._playing = not self._playing

    def _advance_frame(self):
        self.set_step((self._step + 1) % self._n_steps)
