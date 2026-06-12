"""AnimationCanvas — 2-D model animation with shape rendering."""

import time

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.patches import Circle as MplCircle, FancyArrowPatch, FancyBboxPatch, Polygon as MplPolygon
from matplotlib.transforms import Affine2D

from PySide6.QtCore import Qt, QEvent, QTimer, Signal
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDoubleSpinBox,
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
from pmd.core.constraints import RevJoint, TranJoint
from pmd.core.forces import PtpForce
from pmd.core.motion import Motion
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


def _marker_global_angle(marker, step) -> float:
    """Return the marker's global X-axis orientation (rad) at *step*.

    A marker without an explicit ``theta`` still tracks the parent body's
    rotation (it is body-fixed), so the triad rotates with the body.
    """
    body = marker.body
    phi_body = 0.0 if not body else _body_pos(body, step)[2]
    return phi_body + float(getattr(marker, "theta", None) or 0.0)


# ---------------------------------------------------------------
# Marker triad rendering constants (mirror MarkerItem in the
# preprocessor's widgets.marker_item to keep both views visually
# consistent: red X / green Y arrows, yellow origin dot, label
# offset slightly below-left of the origin so it never overlaps
# the dot or the axes).
# ---------------------------------------------------------------
_TRIAD_C_X      = "#d65a5a"   # red  -- X axis
_TRIAD_C_Y      = "#3aa35a"   # green -- Y axis
_TRIAD_C_DOT_FC = "#ffd400"   # yellow origin fill
_TRIAD_C_DOT_EC = "#1c2033"   # origin edge in light theme
_TRIAD_C_DOT_ED = "#f0f0f0"   # origin edge in dark theme
_TRIAD_DOT_SIZE = 3.5         # markersize (pt)
_TRIAD_LW       = 1.4         # arrow line width (pt)
_TRIAD_LABEL_PT = 9.5         # label font size (mirrors MarkerItem)
_TRIAD_LABEL_DX = -2.0        # label offset (points, +x = right)
_TRIAD_LABEL_DY = -6.0        # label offset (points, +y = up)


class AnimationCanvas(QWidget):
    """2-D animation of the multi-body model at each time step.

    Parameters
    ----------
    sessions : list[Session]
        Solved simulation sessions.
    parent : QWidget or None
        Optional parent widget.
    """

    # Emits the current simulation time (seconds) every time the
    # displayed frame changes.  The postprocessor's PlotCanvas connects
    # to it to draw a synchronised vertical cursor on all subplots.
    time_changed = Signal(float)
    # Emitted when the locked-cursor spin box changes (to force-update
    # the PlotCanvas cursor even when it is in locked state).
    cursor_time_changed = Signal(float)
    # Emitted from the locked-play dialog when the user chooses
    # "Unlock & Play" so that main_window can call
    # PlotCanvas.unlock_time_cursor().
    request_cursor_unlock = Signal()

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        self._sessions = sessions
        self._T = sessions[0].T
        self._n_steps = len(self._T)
        self._step = 0
        self._playing = False
        self._cursor_locked = False  # True when the plot cursor is locked
        self._joint_legend = None    # matplotlib Legend for the joint-type key

        # --- matplotlib widgets ---
        self._figure = Figure()
        # Fixed margins prevent tight_layout from recalculating on every
        # draw, which caused the axes box to jump/resize during pan.
        self._figure.subplots_adjust(left=0.06, right=0.97,
                                     top=0.96,  bottom=0.06)
        self._ax = self._figure.add_subplot(111)
        self._ax.set_aspect("equal")
        self._canvas = FigureCanvasQTAgg(self._figure)
        self._canvas.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._canvas.installEventFilter(self)
        self._canvas.mpl_connect("scroll_event", self._on_scroll)

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

        # Playback speed multiplier (× real-time).  Mirrors the
        # preprocessor's AnimationCanvas: rendering happens at a fixed
        # ~60 Hz and the wall-clock-driven _advance_frame() skips frames
        # as needed so that "Speed" actually behaves like a speed knob.
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

        # Display refresh rate — how many render calls per second.
        # Lowering this frees render budget per frame and eliminates
        # choppiness on complex scenes or slow GPUs.  Orthogonal to
        # Speed: Speed governs how much sim-time elapses per real
        # second; Display governs how often the canvas is redrawn.
        self._disp_fps_spin = QDoubleSpinBox()
        self._disp_fps_spin.setDecimals(0)
        self._disp_fps_spin.setRange(5.0, 144.0)
        self._disp_fps_spin.setSingleStep(5.0)
        self._disp_fps_spin.setValue(60.0)
        self._disp_fps_spin.setSuffix(" fps")
        self._disp_fps_spin.setMinimumWidth(80)
        self._disp_fps_spin.setToolTip(
            "Display refresh rate (frames per second).\n"
            "Higher values give smoother playback; lower values reduce\n"
            "CPU/GPU load on complex models.\n"
            "Does not affect playback speed — only how often the canvas redraws.")
        self._disp_fps_spin.valueChanged.connect(self._on_disp_fps_changed)

        # Locked-cursor controls (hidden when cursor is free)
        self._cursor_lock_lbl = QLabel("\U0001f512 t =")
        self._cursor_lock_lbl.setToolTip("Plot cursor locked — edit time to jump to a configuration")
        self._cursor_lock_lbl.setVisible(False)
        self._cursor_time_spin = QDoubleSpinBox()
        self._cursor_time_spin.setDecimals(4)
        _T0 = float(self._T[0])
        _T1 = float(self._T[-1])
        _dT = float(self._T[1] - self._T[0]) if len(self._T) > 1 else 0.001
        self._cursor_time_spin.setRange(_T0, _T1)
        self._cursor_time_spin.setSingleStep(_dT)
        self._cursor_time_spin.setMinimumWidth(100)
        self._cursor_time_spin.setToolTip(
            "Locked cursor time (s) — edit to jump to a specific configuration")
        self._cursor_time_spin.setVisible(False)
        self._cursor_time_spin.valueChanged.connect(self._on_locked_time_changed)

        ctrl = QHBoxLayout()
        ctrl.addWidget(self._play_btn)
        ctrl.addWidget(self._slider, stretch=1)
        ctrl.addWidget(self._time_lbl)
        ctrl.addSpacing(4)
        ctrl.addWidget(self._cursor_lock_lbl)
        ctrl.addWidget(self._cursor_time_spin)
        ctrl.addSpacing(8)
        ctrl.addWidget(QLabel("Speed:"))
        ctrl.addWidget(self._fps_spin)
        ctrl.addSpacing(4)
        ctrl.addWidget(QLabel("Display:"))
        ctrl.addWidget(self._disp_fps_spin)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)
        layout.addLayout(ctrl)

        # --- artist storage ---
        self._body_patches = []
        self._body_info = []           # (body, session_idx)
        self._body_labels = []         # Text at body CM (always visible, colour = body edge)
        self._marker_dots = []
        self._marker_arrows_x = []     # FancyArrowPatch, parallel to _marker_dots
        self._marker_arrows_y = []     # FancyArrowPatch, parallel to _marker_dots
        self._marker_labels = []       # text annotations, parallel to _marker_dots
        self._marker_info = []         # (marker, session_idx)
        self._joint_patches = []       # MplCircle (RevJoint) or FancyBboxPatch (TranJoint)
        self._joint_transforms = []    # Affine2D for TranJoint squares; None otherwise
        self._joint_rail1 = []         # Line2D upper rail (TranJoint); None otherwise
        self._joint_rail2 = []         # Line2D lower rail (TranJoint); None otherwise
        self._joint_info = []          # (joint, session_idx)
        self._joint_coord_labels = []  # Annotation showing global (x, y) per joint
        self._force_lines = []
        self._force_info = []          # (force, session_idx)

        # Triad length (data units) is fixed at init from the full
        # trajectory envelope so the visual size stays stable through
        # the animation and across sessions.  4% of the largest span
        # is small enough to never dominate the body footprint yet big
        # enough to convey orientation.
        self._triad_L = self._compute_triad_length()
        # Joint glyph radius (data-units): decoupled from the scene span so
        # it can be resized at runtime via View → Joint glyph size.
        self._joint_radius: float = self._triad_L * 0.10
        # Artists hidden while playback runs; restored when animation ends.
        self._play_hidden: list = []

        self._init_artists()

        # --- animation timer (wall-clock driven) ---
        # Display refresh at ~60 Hz; the wall-clock anchors are
        # (re-)set every time playback starts or the speed changes.
        self._timer = QTimer(self)
        self._timer.setTimerType(Qt.PreciseTimer)
        self._timer.setInterval(16)   # 60 fps default; user-adjustable via Display spin
        self._timer.timeout.connect(self._advance_frame)
        self._wall_t0 = 0.0
        self._sim_t0 = 0.0

        # "V" toggles marker visibility (mirrors Adams).  The shortcut
        # is wired to the toolbar button so its checked state stays in
        # sync and the slot logic lives in a single place.  Window-
        # scope means it fires whenever this widget's window is
        # active, regardless of which child currently has focus (the
        # matplotlib canvas does not take keyboard focus by default).
        self._markers_shortcut = QShortcut(QKeySequence(Qt.Key_V), self)
        self._markers_shortcut.setContext(Qt.WindowShortcut)
        self._markers_shortcut.activated.connect(self._btn_markers.click)

        # T-key: hold T to activate temporary pan mode while the key is held.
        # Implemented via Qt eventFilter (registered above on self._canvas) so
        # that the shortcut works without the canvas needing matplotlib keyboard
        # focus first.
        self._t_pan_active: bool = False
        self._pre_t_pan_state: bool = False

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
        # Marker visibility toggle — icon is a reference-frame triad so
        # the button's meaning is immediately recognisable. Default ON.
        self._btn_markers = self._make_btn(
            tb, "Show / hide markers (V)", "mdi6.axis-arrow",
            self._on_toggle_markers, checkable=True)
        self._btn_markers.setChecked(True)
        # Joint coordinates toggle — shows the global (x, y) of each
        # joint next to its symbol. Default OFF.
        self._btn_coords = self._make_btn(
            tb, "Show / hide joint coordinates", "mdi6.crosshairs-gps",
            self._on_toggle_coords, checkable=True)
        self._btn_coords.setChecked(False)
        tb.addSeparator()
        self._btn_save = self._make_btn(tb, "Save",    "mdi6.content-save",            self._on_save)
        tb.addSeparator()
        self._btn_export_video = self._make_btn(
            tb, "Export animation video", "mdi6.video-outline",
            self._on_export_video)
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
        self._btn_markers.setIcon(_icons.icon("mdi6.axis-arrow"))
        self._btn_coords.setIcon(_icons.icon("mdi6.crosshairs-gps"))
        self._btn_export_video.setIcon(_icons.icon("mdi6.video-outline"))
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
        # Marker glyphs: keep the yellow fill (pops on both themes)
        # but flip the ring edge + label colour so they stay legible.
        edge = _TRIAD_C_DOT_ED if enabled else _TRIAD_C_DOT_EC
        for dot in self._marker_dots:
            dot.set_markeredgecolor(edge)
        for label in self._marker_labels:
            label.set_color(fg)
        for label in self._joint_coord_labels:
            label.set_color(fg)
        if self._joint_legend is not None:
            self._joint_legend.get_frame().set_facecolor(bg)
            self._joint_legend.get_frame().set_edgecolor(
                "#555555" if enabled else "#cccccc")
            for text in self._joint_legend.get_texts():
                text.set_color(fg)
            title = self._joint_legend.get_title()
            if title:
                title.set_color(fg)
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

    def _on_scroll(self, event):
        """Zoom in/out with the scroll wheel, centred on the cursor."""
        if event.inaxes is not self._ax:
            return
        factor = 0.85 if event.button == "up" else (1.0 / 0.85)
        x0, y0 = event.xdata, event.ydata
        xlim = self._ax.get_xlim()
        ylim = self._ax.get_ylim()
        self._ax.set_xlim([x0 + (x - x0) * factor for x in xlim])
        self._ax.set_ylim([y0 + (y - y0) * factor for y in ylim])
        self._canvas.draw_idle()

    def _on_save(self):
        self._nav.save_figure()

    def _on_toggle_markers(self, checked: bool):
        """Show / hide marker glyphs and their labels (shortcut: V)."""
        for dot in self._marker_dots:
            dot.set_visible(checked)
        for arr in self._marker_arrows_x:
            arr.set_visible(checked)
        for arr in self._marker_arrows_y:
            arr.set_visible(checked)
        for label in self._marker_labels:
            label.set_visible(checked)
        self._canvas.draw_idle()

    def _on_toggle_coords(self, checked: bool):
        """Show / hide global (x, y) coordinates next to each joint."""
        for lbl in self._joint_coord_labels:
            lbl.set_visible(checked)
        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _time_text(self, step):
        return f"t = {self._T[step]:.4f} s"

    def _compute_triad_length(self) -> float:
        """Pick a data-unit length for marker triads.

        Uses the full trajectory envelope (all bodies, all steps, all
        sessions) so the value is independent of the current frame and
        of the matplotlib autoscale state.  Falls back to a sensible
        default if no per-step data is available.
        """
        xs: list[float] = []
        ys: list[float] = []
        for s in self._sessions:
            for body in getattr(s.model, "Bodies", []):
                rc = getattr(body, "_result_container", None)
                if rc is None:
                    continue
                try:
                    xs.extend(float(v) for v in rc["positions"]["x"])
                    ys.extend(float(v) for v in rc["positions"]["y"])
                except Exception:
                    continue
        if not xs or not ys:
            return 0.05
        span = max(max(xs) - min(xs), max(ys) - min(ys), 0.05)
        # 8% of the scene span: large enough that the marker axes are
        # immediately recognisable as a reference frame (rather than a
        # decorative tick) even on dense models.
        return 0.08 * span

    # ------------------------------------------------------------------
    # _init_artists  — draw everything at frame 0
    # ------------------------------------------------------------------

    def _init_artists(self):
        ax = self._ax
        color_idx = 0

        # Co-location bookkeeping: when several markers share (almost)
        # the same scene origin -- typical at a joint endpoint where 2+
        # bodies meet -- their labels would draw on top of one another
        # and become illegible.  Mirror the preprocessor's MarkerItem
        # behaviour and stack them vertically: the first label sits at
        # the default offset, each subsequent label is pushed one line
        # further down.  Tolerance is a fraction of the triad length so
        # it scales with the model.
        label_stack: list[list[float]] = []  # [x, y, count_so_far]
        tol = max(self._triad_L * 0.6, 1e-9)
        tol2 = tol * tol

        def _label_row(px: float, py: float) -> int:
            for slot in label_stack:
                dx = px - slot[0]; dy = py - slot[1]
                if dx * dx + dy * dy <= tol2:
                    row = int(slot[2])
                    slot[2] = row + 1
                    return row
            label_stack.append([px, py, 1.0])
            return 0

        line_h = _TRIAD_LABEL_PT * 1.3   # points per stacked label row

        for si, session in enumerate(self._sessions):
            model = session.model

            # ---- Bodies ----
            for body in model.Bodies:
                if not body:
                    continue                 # skip Ground
                if getattr(body, 'color', None) is not None:
                    col = body.color
                else:
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

                # Body name label — follows the CM, keeps the body's
                # edge colour, never rotated (always readable).  Uses
                # annotate so _auto_place_labels() can set the offset.
                _bname = getattr(body, "name", "") or ""
                _blbl = ax.annotate(
                    _bname,
                    xy=(x, y),
                    xytext=(12, 12),
                    textcoords="offset points",
                    ha="left", va="bottom",
                    fontsize=9.5, fontweight="bold",
                    color=col, zorder=11,
                    annotation_clip=False,
                    visible=bool(_bname),
                    arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.6,
                                    shrinkA=0, shrinkB=2),
                )
                _blbl.draggable(True)
                self._body_labels.append(_blbl)

                # markers attached to this body, drawn as a small
                # red-X / green-Y triad with a yellow origin dot and a
                # name label offset slightly below-left.  Matches the
                # preprocessor's MarkerItem so both views read the
                # same.  Triad length is chosen as a fraction of the
                # trajectory envelope (see _compute_triad_length).
                L_triad = self._triad_L
                for mk_idx, mk in enumerate(body._markers):
                    gp = _marker_global(mk, 0)
                    a = _marker_global_angle(mk, 0)
                    head_x = (gp[0] + L_triad * np.cos(a),
                              gp[1] + L_triad * np.sin(a))
                    head_y = (gp[0] - L_triad * np.sin(a),
                              gp[1] + L_triad * np.cos(a))
                    arrow_x = FancyArrowPatch(
                        (gp[0], gp[1]), head_x,
                        arrowstyle="-|>",
                        mutation_scale=8,
                        color=_TRIAD_C_X,
                        linewidth=_TRIAD_LW,
                        shrinkA=0, shrinkB=0,
                        zorder=5,
                    )
                    arrow_y = FancyArrowPatch(
                        (gp[0], gp[1]), head_y,
                        arrowstyle="-|>",
                        mutation_scale=8,
                        color=_TRIAD_C_Y,
                        linewidth=_TRIAD_LW,
                        shrinkA=0, shrinkB=0,
                        zorder=5,
                    )
                    ax.add_patch(arrow_x)
                    ax.add_patch(arrow_y)
                    dot, = ax.plot(
                        gp[0], gp[1],
                        marker="o",
                        linestyle="none",
                        markersize=_TRIAD_DOT_SIZE,
                        markerfacecolor=_TRIAD_C_DOT_FC,
                        markeredgecolor=_TRIAD_C_DOT_EC,
                        markeredgewidth=0.8,
                        zorder=6,
                    )
                    # Show a label only for markers that carry an explicit
                    # name.  Structural markers created by as_link/as_plate
                    # have name=None and appear as unlabelled triads.
                    mk_name = getattr(mk, "name", "") or ""
                    row = _label_row(gp[0], gp[1])
                    label = ax.annotate(
                        mk_name,
                        xy=(gp[0], gp[1]),
                        xytext=(_TRIAD_LABEL_DX,
                                _TRIAD_LABEL_DY - row * line_h),
                        textcoords="offset points",
                        fontsize=_TRIAD_LABEL_PT,
                        color="#1c2033",
                        ha="right",
                        va="top",
                        zorder=7,
                        annotation_clip=False,
                        visible=bool(mk_name),
                        arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.6,
                                        shrinkA=0, shrinkB=2),
                    )
                    label.draggable(True)
                    self._marker_dots.append(dot)
                    self._marker_arrows_x.append(arrow_x)
                    self._marker_arrows_y.append(arrow_y)
                    self._marker_labels.append(label)
                    self._marker_info.append((mk, si))

            # ---- Joints ----
            _r_joint = self._joint_radius          # RevJoint disc radius
            _sh      = self._joint_radius          # TranJoint square half-side
            _rail_hw = self._joint_radius * 2.5    # TranJoint rail half-width
            _gap_y   = self._joint_radius * 1.6    # perpendicular offset of rails
            for joint in model.Joints:
                if isinstance(joint, Motion):
                    continue
                mk = joint.iMarker or joint.jMarker
                if mk is None:
                    continue
                gp  = _marker_global(mk, 0)
                ang = _marker_global_angle(mk, 0)
                t_sq = None; r1_art = None; r2_art = None

                if isinstance(joint, RevJoint):
                    patch = MplCircle(
                        gp, radius=_r_joint,
                        facecolor="#ff2bd6", edgecolor="#1c2033",
                        linewidth=1.0, zorder=7,
                    )
                    ax.add_patch(patch)
                elif isinstance(joint, TranJoint):
                    t_sq = Affine2D().rotate(ang).translate(gp[0], gp[1])
                    patch = FancyBboxPatch(
                        (-_sh, -_sh), 2 * _sh, 2 * _sh,
                        boxstyle="square,pad=0",
                        facecolor="#00bcd4", edgecolor="#006978",
                        linewidth=1.2, zorder=7,
                    )
                    patch.set_transform(t_sq + ax.transData)
                    ax.add_patch(patch)
                    # Two rails parallel to the joint translation axis,
                    # offset ±gap_y in the local-Y direction.
                    c, s = np.cos(ang), np.sin(ang)
                    ex = np.array([c, s]); ey = np.array([-s, c])
                    p1 = gp - _rail_hw * ex + _gap_y * ey
                    p2 = gp + _rail_hw * ex + _gap_y * ey
                    p3 = gp - _rail_hw * ex - _gap_y * ey
                    p4 = gp + _rail_hw * ex - _gap_y * ey
                    r1_art, = ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                                      color="#006978", linewidth=1.4,
                                      solid_capstyle="butt", zorder=7)
                    r2_art, = ax.plot([p3[0], p4[0]], [p3[1], p4[1]],
                                      color="#006978", linewidth=1.4,
                                      solid_capstyle="butt", zorder=7)
                else:
                    patch = MplCircle(
                        gp, radius=_r_joint * 0.5,
                        facecolor="#888888", edgecolor="#333333",
                        linewidth=0.8, zorder=7,
                    )
                    ax.add_patch(patch)

                self._joint_patches.append(patch)
                self._joint_transforms.append(t_sq)
                self._joint_rail1.append(r1_art)
                self._joint_rail2.append(r2_art)
                self._joint_info.append((joint, si))

                # Coordinate annotation — hidden until toolbar button ON.
                _jname = getattr(joint, "name", "") or ""
                _clbl = ax.annotate(
                    f"{_jname}\n({gp[0]:.1f}, {gp[1]:.1f})",
                    xy=(gp[0], gp[1]),
                    xytext=(6, 6), textcoords="offset points",
                    fontsize=8.5, color="#1c2033",
                    ha="left", va="bottom",
                    zorder=9, annotation_clip=False,
                    visible=False,
                )
                self._joint_coord_labels.append(_clbl)

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
        # Disable autoscale BEFORE computing limits: FancyArrowPatch
        # patches are added in data-coordinates and can silently
        # expand ax.dataLim after an explicit set_xlim/set_ylim call
        # when autoscale is still active.
        self._ax.autoscale(False)
        self._auto_limits()
        self._auto_place_labels()
        self._build_joint_legend()
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
        # off-CM features the user attached a frame to.  We also add
        # the four cardinal triad tip positions (origin ± L in both
        # data axes) so the worst-case arrow direction never escapes
        # the padded envelope regardless of body orientation.
        L = self._triad_L
        for mk, _si in self._marker_info:
            body = mk.body
            lp = mk.local_position.ravel()[:2]
            if not body:
                px, py = float(lp[0]), float(lp[1])
                xs.extend([px - L, px + L])
                ys.extend([py - L, py + L])
                continue
            rc = body._result_container
            bx = np.asarray(rc["positions"]["x"], dtype=float)
            by = np.asarray(rc["positions"]["y"], dtype=float)
            ph = np.asarray(rc["positions"]["phi"], dtype=float)
            cos_p = np.cos(ph); sin_p = np.sin(ph)
            mx_t = bx + cos_p * lp[0] - sin_p * lp[1]
            my_t = by + sin_p * lp[0] + cos_p * lp[1]
            xs.extend([float(mx_t.min()) - L, float(mx_t.max()) + L])
            ys.extend([float(my_t.min()) - L, float(my_t.max()) + L])

        if not xs:
            return
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)
        # 15% breathing room on the larger span.  The triad extents
        # are already baked into xs/ys above so no extra _triad_L term
        # is needed here.
        span = max((x_max - x_min), (y_max - y_min), 0.1)
        pad = span * 0.15
        self._ax.set_xlim(x_min - pad, x_max + pad)
        self._ax.set_ylim(y_min - pad, y_max + pad)

    # ------------------------------------------------------------------
    # _build_joint_legend
    # ------------------------------------------------------------------

    def _build_joint_legend(self) -> None:
        """Create (or refresh) a joint-type key in the upper-right corner.

        Shows one sample symbol per joint TYPE present in the model,
        with a human-readable label.  The legend is recreated from
        scratch whenever the model is loaded so it stays in sync.
        """
        from matplotlib.lines import Line2D as _Line2D

        # Gather which joint types are present across all sessions.
        has_rev  = False
        has_tran = False
        for session in self._sessions:
            for joint in session.model.Joints:
                if isinstance(joint, RevJoint):
                    has_rev = True
                elif isinstance(joint, TranJoint):
                    has_tran = True

        handles = []
        if has_rev:
            handles.append(_Line2D(
                [0], [0], marker="o", linestyle="none",
                markerfacecolor="#ff2bd6", markeredgecolor="#1c2033",
                markeredgewidth=0.8, markersize=9,
                label="Revolute Joint",
            ))
        if has_tran:
            handles.append(_Line2D(
                [0], [0], marker="s", linestyle="none",
                markerfacecolor="#00bcd4", markeredgecolor="#006978",
                markeredgewidth=0.8, markersize=9,
                label="Translational Joint",
            ))

        if not handles:
            return

        if self._joint_legend is not None:
            try:
                self._joint_legend.remove()
            except Exception:
                pass

        self._joint_legend = self._ax.legend(
            handles=handles,
            loc="upper right",
            fontsize=8.5,
            framealpha=0.88,
            edgecolor="#cccccc",
            title="Joint types",
            title_fontsize=8.0,
            borderpad=0.6,
            handlelength=1.2,
            handletextpad=0.5,
        )
        self._joint_legend.set_zorder(12)

    # ------------------------------------------------------------------
    # _auto_place_labels
    # ------------------------------------------------------------------

    def _auto_place_labels(self) -> None:
        """Greedy label placement: adjust each annotation's pixel offset so
        estimated bounding boxes have minimal overlap at frame 0.

        Priority order: joint-coord labels → marker labels → body names.
        Joint labels are reserved even when invisible so that toggling them
        on later does not suddenly overlap other labels.
        """
        dpi = self._figure.get_dpi()
        pt_to_px = dpi / 72.0

        # Candidate offsets (dx_pt, dy_pt) in priority order per category.
        JOINT_CANDS = [
            (  6,   6), ( -6,   6), (  6, -10), ( -6, -10),
            ( 12,   0), (-12,   0), (  0,  10), (  0, -12),
            ( 18,   6), (-18,   6), ( 18, -10), (-18, -10),
        ]
        MK_CANDS = [
            (  5,  -4), ( -4, -14), (  5,   8), ( -4,   8),
            ( 14,  -4), (-18,  -4), ( 14,   8), (-18,   8),
            (  5, -18), ( -4, -18), ( 24,  -4), (-24,  -4),
        ]
        BODY_CANDS = [
            ( 44,  24), (-44,  24), ( 44, -24), (-44, -24),
            ( 62,   0), (-62,   0), (  0,  44), (  0, -44),
            ( 60,  38), (-60,  38), ( 60, -38), (-60, -38),
            ( 80,  28), (-80,  28), ( 80, -28), (-80, -28),
        ]

        def _est_px(ann) -> tuple[float, float]:
            """Estimate rendered (width, height) in display pixels."""
            txt = ann.get_text() or " "
            lines = txt.split("\n")
            n_c = max((len(l) for l in lines), default=5)
            fs  = ann.get_fontsize()
            return (n_c * fs * 0.58 * pt_to_px,
                    len(lines) * fs * 1.38 * pt_to_px)

        def _overlap(a, b) -> float:
            ox = max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
            oy = max(0.0, min(a[3], b[3]) - max(a[1], b[1]))
            return ox * oy

        placed: list[tuple[float, float, float, float]] = []

        def _place(ann, cands, anchor_data) -> None:
            try:
                adx, ady = self._ax.transData.transform(anchor_data)
            except Exception:
                return
            w, h = _est_px(ann)
            best_off = cands[0]
            best_ov  = float("inf")
            for (ox_pt, oy_pt) in cands:
                bx0 = adx + ox_pt * pt_to_px
                by0 = ady + oy_pt * pt_to_px
                box = (bx0, by0, bx0 + w, by0 + h)
                ov  = sum(_overlap(box, pb) for pb in placed)
                if ov < best_ov:
                    best_ov  = ov
                    best_off = (ox_pt, oy_pt)
                    if ov == 0.0:
                        break
            ox_pt, oy_pt = best_off
            ann.xyann = (ox_pt, oy_pt)
            bx0 = adx + ox_pt * pt_to_px
            by0 = ady + oy_pt * pt_to_px
            placed.append((bx0, by0, bx0 + w, by0 + h))

        # 1. Joint coord labels (always reserved, even when invisible)
        for clbl, (joint, _) in zip(self._joint_coord_labels, self._joint_info):
            mk = joint.iMarker or joint.jMarker
            if mk is not None:
                _place(clbl, JOINT_CANDS, _marker_global(mk, 0))

        # 2. Marker labels
        for lbl, (mk, _) in zip(self._marker_labels, self._marker_info):
            _place(lbl, MK_CANDS, _marker_global(mk, 0))

        # 3. Body names (largest offset to clear the body boundary)
        for lbl, (body, _) in zip(self._body_labels, self._body_info):
            if not lbl.get_visible():
                continue
            x, y, _ = _body_pos(body, 0)
            _place(lbl, BODY_CANDS, (x, y))

    # ------------------------------------------------------------------
    # _update_artists — move everything to *step*
    # ------------------------------------------------------------------

    def _update_artists(self, step):
        ax = self._ax

        # bodies
        for patch, (body, _si) in zip(self._body_patches, self._body_info):
            x, y, phi = _body_pos(body, step)
            shape = body.shape

            if isinstance(shape, (Rectangle, Plate, Link)):
                t = Affine2D().rotate(phi).translate(x, y) + ax.transData
                patch.set_transform(t)
            elif isinstance(shape, Circle):
                patch.set_center((x, y))
            else:
                patch.set_center((x, y))

        # body name labels — follow the CM (xy anchor; xyann fixed by placement)
        for lbl, (body, _si) in zip(self._body_labels, self._body_info):
            bx, by, _ = _body_pos(body, step)
            lbl.xy = (bx, by)

        # marker triads (dot + X arrow + Y arrow + label)
        L_triad = self._triad_L
        for i, (mk, _si) in enumerate(self._marker_info):
            gp = _marker_global(mk, step)
            a  = _marker_global_angle(mk, step)
            head_x = (gp[0] + L_triad * np.cos(a),
                      gp[1] + L_triad * np.sin(a))
            head_y = (gp[0] - L_triad * np.sin(a),
                      gp[1] + L_triad * np.cos(a))
            self._marker_arrows_x[i].set_positions((gp[0], gp[1]), head_x)
            self._marker_arrows_y[i].set_positions((gp[0], gp[1]), head_y)
            self._marker_dots[i].set_data([gp[0]], [gp[1]])
            # Only move the anchor (`xy`).  ``set_position`` would
            # overwrite ``xyann`` and -- because we use
            # ``textcoords="offset points"`` -- silently destroy the
            # row-stacking offsets computed at init time.
            self._marker_labels[i].xy = (gp[0], gp[1])

        # joints
        _rail_hw = self._joint_radius * 2.5
        _gap_y   = self._joint_radius * 1.6
        for patch, t_sq, r1, r2, clbl, (joint, _si) in zip(
            self._joint_patches, self._joint_transforms,
            self._joint_rail1, self._joint_rail2,
            self._joint_coord_labels, self._joint_info,
        ):
            mk = joint.iMarker or joint.jMarker
            if mk is None:
                continue
            gp  = _marker_global(mk, step)
            ang = _marker_global_angle(mk, step)
            if isinstance(joint, RevJoint):
                patch.set_center(gp)
            elif isinstance(joint, TranJoint) and t_sq is not None:
                t_sq.clear().rotate(ang).translate(gp[0], gp[1])
                c, s = np.cos(ang), np.sin(ang)
                ex = np.array([c, s]); ey = np.array([-s, c])
                p1 = gp - _rail_hw * ex + _gap_y * ey
                p2 = gp + _rail_hw * ex + _gap_y * ey
                p3 = gp - _rail_hw * ex - _gap_y * ey
                p4 = gp + _rail_hw * ex - _gap_y * ey
                r1.set_data([p1[0], p2[0]], [p1[1], p2[1]])
                r2.set_data([p3[0], p4[0]], [p3[1], p4[1]])
            else:
                patch.set_center(gp)
            # Update coordinate annotation text and anchor.
            _jname = getattr(joint, "name", "") or ""
            clbl.xy = (gp[0], gp[1])
            clbl.set_text(f"{_jname}\n({gp[0]:.1f}, {gp[1]:.1f})")

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
        # Emit the simulation time (seconds), not the step index, so
        # PlotCanvas.set_time_cursor() receives the value it expects.
        self.time_changed.emit(float(self._T[step]))

    # ------------------------------------------------------------------
    # slots
    # ------------------------------------------------------------------

    def _on_slider_changed(self, value):
        self.set_step(value)

    def _on_play_pause(self):
        if not self._playing and self._cursor_locked:
            # Cursor is locked — ask before starting playback.
            from PySide6.QtWidgets import QMessageBox
            dlg = QMessageBox(self)
            dlg.setWindowTitle("Cursor locked")
            dlg.setText(
                "The plot cursor is currently locked.\n"
                "Unlock it to resume animation playback.")
            dlg.setIcon(QMessageBox.Icon.Warning)
            btn_unlock = dlg.addButton(
                "Unlock & Play", QMessageBox.ButtonRole.AcceptRole)
            dlg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            dlg.exec()
            if dlg.clickedButton() is not btn_unlock:
                return
            # Unlock locally and ask PlotCanvas to unlock too.
            self.set_cursor_locked(False)
            self.request_cursor_unlock.emit()
        if self._playing:
            self._timer.stop()
            self._play_btn.setIcon(_icons.icon("mdi6.play"))
        else:
            if self._step >= self._n_steps - 1:
                self._step = 0
            # Pin wall-clock to the current sim time so playback resumes
            # smoothly from here instead of jumping forward.
            self._anchor_playback()
            self._timer.start()
            self._play_btn.setIcon(_icons.icon("mdi6.pause"))
            self._hide_for_play()
        self._playing = not self._playing

    def _anchor_playback(self):
        """Pin (wall_now, sim_now) so wall-clock advance resumes here."""
        self._wall_t0 = time.monotonic()
        try:
            self._sim_t0 = float(self._T[self._step])
        except Exception:
            self._sim_t0 = 0.0

    def set_cursor_locked(self, locked: bool, t: float | None = None) -> None:
        """Show / hide the locked-cursor UI and update internal state.

        Called by main_window when PlotCanvas.cursor_lock_changed fires,
        and internally from the "Unlock & Play" dialog branch.
        """
        self._cursor_locked = locked
        self._cursor_lock_lbl.setVisible(locked)
        self._cursor_time_spin.setVisible(locked)
        if locked and t is not None:
            self._cursor_time_spin.blockSignals(True)
            self._cursor_time_spin.setValue(t)
            self._cursor_time_spin.blockSignals(False)

    def _on_locked_time_changed(self, t: float) -> None:
        """Handle spin-box edits while the cursor is locked.

        Jumps the animation to the step nearest to *t* and emits
        cursor_time_changed so PlotCanvas can force-update its cursor.
        """
        if not self._cursor_locked:
            return
        T = self._T
        step = int(np.argmin(np.abs(T - t)))
        step = max(0, min(step, self._n_steps - 1))
        self.set_step(step)
        self.cursor_time_changed.emit(float(T[step]))

    def _on_speed_changed(self, _val: float):
        # Re-anchor so the speed change takes effect from "now" without
        # a discontinuous jump in displayed time.
        if self._playing:
            self._anchor_playback()

    def _on_disp_fps_changed(self, fps: float) -> None:
        """Update the render timer interval when the Display spin changes."""
        self._timer.setInterval(max(1, int(1000.0 / max(fps, 1.0))))

    def _advance_frame(self):
        # Wall-clock driven: figure out which step corresponds to the
        # elapsed wall time at the current playback speed, snap to it,
        # and stop at the end of the trajectory.
        speed = max(1e-6, float(self._fps_spin.value()))
        elapsed = time.monotonic() - self._wall_t0
        sim_t = self._sim_t0 + elapsed * speed
        try:
            t_end = float(self._T[-1])
        except Exception:
            t_end = 0.0
        if sim_t >= t_end:
            self.set_step(self._n_steps - 1)
            self._timer.stop()
            self._playing = False
            self._play_btn.setIcon(_icons.icon("mdi6.play"))
            self._restore_play_hidden()
            return
        # Linear forward scan (T is monotonic; we never go backwards here).
        s = self._step
        n = self._n_steps
        while s + 1 < n and float(self._T[s + 1]) <= sim_t:
            s += 1
        if s != self._step:
            self.set_step(s)

    # ------------------------------------------------------------------
    # Video export
    # ------------------------------------------------------------------

    def _on_export_video(self) -> None:
        """Open the export dialog and run the offscreen render pipeline."""
        from ..dialogs import ExportVideoDialog
        from PySide6.QtWidgets import QMessageBox
        dlg = ExportVideoDialog(
            current_speed=float(self._fps_spin.value()),
            parent=self,
        )
        if dlg.exec() != ExportVideoDialog.DialogCode.Accepted:
            return
        try:
            self._export_video(
                output_path=dlg.output_path,
                video_fps=dlg.video_fps,
                speed=dlg.speed,
                layout=dlg.layout,
                dpi=dlg.dpi,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        QMessageBox.information(
            self, "Export complete",
            f"Video saved to:\n{dlg.output_path}"
        )

    def _export_video(
        self,
        output_path: str,
        video_fps: int = 30,
        speed: float = 1.0,
        layout: str = "anim",
        dpi: int = 100,
    ) -> None:
        """Render an offscreen video of the animation (and optionally the plots).

        Parameters
        ----------
        output_path : str
            Destination file path (.mp4 recommended).
        video_fps : int
            Output video frame rate.
        speed : float
            Simulation time per real second in the video.
        layout : str
            ``"anim"`` — animation canvas only.
            ``"combo"`` — animation + plot canvas side by side.
        dpi : int
            Rendering resolution for both canvases.
        """
        import shutil
        import subprocess
        import numpy as np
        from PySide6.QtWidgets import QProgressDialog, QApplication
        from PySide6.QtCore import Qt

        if not shutil.which("ffmpeg"):
            raise RuntimeError(
                "ffmpeg not found in PATH.\n"
                "Install ffmpeg and ensure it is accessible from the command line."
            )

        # -- compute which steps to render ------------------------------------
        T = self._T
        n_steps = self._n_steps
        if n_steps < 2:
            raise RuntimeError("Need at least 2 time steps to export a video.")
        dt_sim = float(T[-1] - T[0])
        dt_per_frame = speed / video_fps      # sim-time per video frame
        # Build a list of step indices, one per video frame
        frame_steps: list[int] = []
        t_cursor = float(T[0])
        while t_cursor <= float(T[-1]) + 1e-12:
            idx = int(np.argmin(np.abs(T - t_cursor)))
            frame_steps.append(idx)
            t_cursor += dt_per_frame
        if not frame_steps:
            raise RuntimeError("No frames to export.")

        # -- canvas geometry --------------------------------------------------
        # Pixel sizes are computed from the figure's logical size × export DPI.
        # We NEVER call set_dpi() or canvas.draw() on the live Qt figures:
        # those operations resize the embedded widget and redraw on screen,
        # causing the visible layout to shift during export.  Instead every
        # frame is rendered fully offscreen via savefig(format='raw') which
        # invokes the Agg renderer in memory without touching the Qt canvas.
        #
        # Two sets of dimensions are kept separate:
        #   *_raw  — what matplotlib's Agg actually renders: round(inches × dpi)
        #            → used for reshaping the savefig buffer (must match exactly)
        #   W / H  — ffmpeg frame size, rounded UP to nearest even pixel as
        #            required by yuv420p; any leftover edge pixels stay black.
        # Mixing the two caused the "cannot reshape array" crash.

        import io as _io

        def _even(n: int) -> int:
            """Round up to nearest even pixel (required by yuv420p)."""
            return n + n % 2

        fig_a   = self._figure
        w_a_raw = int(fig_a.get_figwidth()  * dpi)
        h_a_raw = int(fig_a.get_figheight() * dpi)

        plot_canvas_ref = getattr(self, "_plot_canvas_ref", None)
        use_combo = (layout == "combo") and (plot_canvas_ref is not None)
        if use_combo:
            fig_p   = plot_canvas_ref._figure
            w_p_raw = int(fig_p.get_figwidth()  * dpi)
            h_p_raw = int(fig_p.get_figheight() * dpi)
            # Mirror the GUI splitter order: PlotCanvas LEFT, AnimCanvas RIGHT.
            # Pad the TOTAL frame to even for yuv420p; individual bufs use raw dims.
            W = _even(w_p_raw + w_a_raw)
            H = _even(max(h_p_raw, h_a_raw))
        else:
            W = _even(w_a_raw)
            H = _even(h_a_raw)

        # -- stash current state ----------------------------------------------
        saved_step = self._step
        was_playing = self._playing
        if was_playing:
            self._timer.stop()
            self._playing = False
            self._play_btn.setIcon(_icons.icon("mdi6.play"))

        self._hide_for_play()

        # -- launch ffmpeg pipe -----------------------------------------------
        cmd = [
            "ffmpeg", "-y",
            "-f", "rawvideo", "-vcodec", "rawvideo",
            "-s", f"{W}x{H}",
            "-pix_fmt", "rgb24",
            "-r", str(video_fps),
            "-i", "pipe:0",
            "-vcodec", "libx264",
            "-pix_fmt", "yuv420p",
            output_path,
        ]
        pipe = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)

        # -- progress dialog --------------------------------------------------
        n_frames = len(frame_steps)
        prog = QProgressDialog(
            "Exporting video…", "Cancel", 0, n_frames, self
        )
        prog.setWindowTitle("Export animation video")
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)
        prog.setValue(0)

        cancelled = False
        try:
            for fi, step in enumerate(frame_steps):
                if prog.wasCanceled():
                    cancelled = True
                    break

                # Advance animation state; set_step emits time_changed which
                # updates the PlotCanvas time cursor (used in combo layout).
                self.set_step(step)

                # Render offscreen via savefig — pure Agg, no Qt canvas touched.
                # Reshape MUST use the *_raw dimensions (exact Agg output size);
                # the frame array is larger (W/H, even-padded) and stays black
                # in any padding pixels.
                _buf_a = _io.BytesIO()
                fig_a.savefig(_buf_a, format="raw", dpi=dpi, bbox_inches=None)
                buf_a = np.frombuffer(_buf_a.getvalue(), dtype=np.uint8).reshape(h_a_raw, w_a_raw, 4)[:, :, :3]

                if use_combo:
                    _buf_p = _io.BytesIO()
                    fig_p.savefig(_buf_p, format="raw", dpi=dpi, bbox_inches=None)
                    buf_p = np.frombuffer(_buf_p.getvalue(), dtype=np.uint8).reshape(h_p_raw, w_p_raw, 4)[:, :, :3]
                    # PlotCanvas LEFT, AnimCanvas RIGHT — mirrors GUI splitter order.
                    frame = np.zeros((H, W, 3), dtype=np.uint8)
                    frame[:h_p_raw, :w_p_raw]              = buf_p
                    frame[:h_a_raw, w_p_raw:w_p_raw + w_a_raw] = buf_a
                else:
                    frame = np.zeros((H, W, 3), dtype=np.uint8)
                    frame[:h_a_raw, :w_a_raw] = buf_a

                pipe.stdin.write(frame.tobytes())

                if fi % 5 == 0:
                    prog.setValue(fi)
                    QApplication.processEvents()

        finally:
            pipe.stdin.close()
            pipe.wait()
            prog.close()
            self._restore_play_hidden()
            self.set_step(saved_step)

        if cancelled:
            import os
            try:
                os.remove(output_path)
            except OSError:
                pass
            raise RuntimeError("Export cancelled by user.")

    # ------------------------------------------------------------------
    # Hide / restore artists during playback (Feature A)
    # ------------------------------------------------------------------

    def _hide_for_play(self) -> None:
        """Hide labels and marker triads while playback is running."""
        self._play_hidden.clear()
        for art in (
            self._body_labels
            + self._marker_arrows_x
            + self._marker_arrows_y
            + self._marker_dots
            + self._marker_labels
            + self._joint_coord_labels
        ):
            if art.get_visible():
                art.set_visible(False)
                self._play_hidden.append(art)
        self._canvas.draw_idle()

    def _restore_play_hidden(self) -> None:
        """Restore artists that were hidden by _hide_for_play."""
        for art in self._play_hidden:
            art.set_visible(True)
        self._play_hidden.clear()
        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # T-key temporary pan mode (Feature B) — Qt eventFilter
    # ------------------------------------------------------------------

    def eventFilter(self, obj, event) -> bool:
        """Intercept Qt key events on the matplotlib canvas.

        Holding **T** activates pan mode; releasing T restores the
        previous pan state.  Using Qt's event filter rather than
        ``mpl_connect("key_press_event")`` avoids the requirement that
        the matplotlib canvas already has Qt keyboard focus.
        """
        if obj is self._canvas:
            etype = event.type()
            if etype == QEvent.Type.KeyPress:
                if event.key() == Qt.Key.Key_T and not event.isAutoRepeat():
                    if not self._t_pan_active:
                        self._t_pan_active = True
                        self._pre_t_pan_state = self._btn_pan.isChecked()
                        self._btn_pan.setChecked(True)
                        self._on_pan(True)
                    return True
            elif etype == QEvent.Type.KeyRelease:
                if event.key() == Qt.Key.Key_T and not event.isAutoRepeat():
                    if self._t_pan_active:
                        self._t_pan_active = False
                        if not self._pre_t_pan_state:
                            self._btn_pan.setChecked(False)
                            self._on_pan(False)
                    return True
        return super().eventFilter(obj, event)

    # ------------------------------------------------------------------
    # Joint glyph resize (Feature E)
    # ------------------------------------------------------------------

    def _rebuild_joint_glyphs(self, new_radius: float) -> None:
        """Remove existing joint patches and recreate them with *new_radius*.

        Called by main_window when the user changes the joint glyph size
        via View → Joint glyph size.
        """
        self._joint_radius = new_radius
        ax = self._ax
        step = self._step

        # Remove old patches and rails from the axes.
        for p in self._joint_patches:
            try:
                p.remove()
            except Exception:
                pass
        for r in self._joint_rail1:
            if r is not None:
                try:
                    r.remove()
                except Exception:
                    pass
        for r in self._joint_rail2:
            if r is not None:
                try:
                    r.remove()
                except Exception:
                    pass

        self._joint_patches.clear()
        self._joint_transforms.clear()
        self._joint_rail1.clear()
        self._joint_rail2.clear()

        _r = new_radius
        _sh = _r
        _rail_hw = _r * 2.5
        _gap_y = _r * 1.6

        for joint, _si in self._joint_info:
            mk = joint.iMarker or joint.jMarker
            gp = _marker_global(mk, step)
            ang = _marker_global_angle(mk, step)
            t_sq = None
            r1_art = None
            r2_art = None

            if isinstance(joint, RevJoint):
                patch = MplCircle(
                    gp, radius=_r,
                    facecolor="#ff2bd6", edgecolor="#1c2033",
                    linewidth=1.0, zorder=7,
                )
                ax.add_patch(patch)
            elif isinstance(joint, TranJoint):
                t_sq = Affine2D().rotate(ang).translate(gp[0], gp[1])
                patch = FancyBboxPatch(
                    (-_sh, -_sh), 2 * _sh, 2 * _sh,
                    boxstyle="square,pad=0",
                    facecolor="#00bcd4", edgecolor="#006978",
                    linewidth=1.2, zorder=7,
                )
                patch.set_transform(t_sq + ax.transData)
                ax.add_patch(patch)
                c, s = np.cos(ang), np.sin(ang)
                ex = np.array([c, s])
                ey = np.array([-s, c])
                p1 = gp - _rail_hw * ex + _gap_y * ey
                p2 = gp + _rail_hw * ex + _gap_y * ey
                p3 = gp - _rail_hw * ex - _gap_y * ey
                p4 = gp + _rail_hw * ex - _gap_y * ey
                r1_art, = ax.plot(
                    [p1[0], p2[0]], [p1[1], p2[1]],
                    color="#006978", linewidth=1.4,
                    solid_capstyle="butt", zorder=7,
                )
                r2_art, = ax.plot(
                    [p3[0], p4[0]], [p3[1], p4[1]],
                    color="#006978", linewidth=1.4,
                    solid_capstyle="butt", zorder=7,
                )
            else:
                patch = MplCircle(
                    gp, radius=_r * 0.5,
                    facecolor="#888888", edgecolor="#333333",
                    linewidth=0.8, zorder=7,
                )
                ax.add_patch(patch)

            self._joint_patches.append(patch)
            self._joint_transforms.append(t_sq)
            self._joint_rail1.append(r1_art)
            self._joint_rail2.append(r2_art)

        self._canvas.draw_idle()
