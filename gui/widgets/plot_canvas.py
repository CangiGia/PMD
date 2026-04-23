"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from PySide6.QtCore import Signal

from ._zoom_inset import ZoomInset
from PySide6.QtWidgets import QToolBar, QToolButton, QVBoxLayout, QWidget

from ..models import CurveItem
from .. import icons as _icons
from ..style import CANVAS_BG_DARK, CANVAS_BG_LIGHT, CANVAS_FG_DARK, CANVAS_FG_LIGHT


class PlotCanvas(QWidget):
    """Embeds a matplotlib Figure + NavigationToolbar.

    Public API
    ----------
    update_plot(curves: list[CurveItem])
        Clear the axes and redraw only the given curves.
    """

    step_requested = Signal(int)  # emitted on click: nearest time-step index

    def __init__(self, parent=None):
        super().__init__(parent)

        self._figure = Figure(tight_layout=True)
        self._canvas = FigureCanvasQTAgg(self._figure)

        # Backend nav toolbar (hidden) — used only for navigate state/history
        self._nav = NavigationToolbar2QT(self._canvas, self)
        self._nav.hide()

        self._curves: list = []
        self._dark = False
        ax0 = self._figure.add_subplot(111)
        self._axes: list = [ax0]
        self._canvas.mpl_connect("button_press_event", self._on_click)

        # Zoom-inset state
        self._zoom_insets: list[ZoomInset] = []
        self._pending_selectors: list = []  # temporary selectors during "Add Zoom" mode

        # Data-cursor state
        self._cursor_artists: dict = {}     # ax -> {'vline','hline','text'}
        self._cid_cursor_motion: int | None = None

        # Custom visible toolbar (built after instance vars)
        self._toolbar = self._build_toolbar()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

    def _build_toolbar(self) -> QToolBar:
        tb = QToolBar()
        tb.setMovable(False)
        self._btn_home = self._make_btn(tb, "Home",    "mdi6.home",                    self._on_home)
        self._btn_back = self._make_btn(tb, "Back",    "mdi6.arrow-left",              self._on_back)
        self._btn_fwd  = self._make_btn(tb, "Forward", "mdi6.arrow-right",             self._on_fwd)
        tb.addSeparator()
        self._btn_pan      = self._make_btn(tb, "Pan",      "mdi6.hand-back-right-outline", self._on_pan,      checkable=True)
        self._btn_zoom     = self._make_btn(tb, "Zoom",     "mdi6.magnify",                 self._on_zoom,     checkable=True)
        self._btn_add_zoom = self._make_btn(tb, "Add Zoom", "mdi6.selection",               self._on_add_zoom, checkable=True)
        tb.addSeparator()
        self._btn_cursor   = self._make_btn(tb, "Cursor",   "mdi6.crosshairs-gps",          self._on_cursor,   checkable=True)
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
        self._btn_home.setIcon(    _icons.icon("mdi6.home"))
        self._btn_back.setIcon(    _icons.icon("mdi6.arrow-left"))
        self._btn_fwd.setIcon(     _icons.icon("mdi6.arrow-right"))
        self._btn_pan.setIcon(     _icons.icon("mdi6.hand-back-right-outline"))
        self._btn_zoom.setIcon(    _icons.icon("mdi6.magnify"))
        self._btn_add_zoom.setIcon(_icons.icon("mdi6.selection"))
        self._btn_cursor.setIcon(  _icons.icon("mdi6.crosshairs-gps"))

    # -- Toolbar actions ------------------------------------------------

    def _on_home(self):
        self._nav.home()
        self._canvas.draw_idle()

    def _on_back(self):
        self._nav.back()
        self._canvas.draw_idle()

    def _on_fwd(self):
        self._nav.forward()
        self._canvas.draw_idle()

    def _sync_nav_mode(self, target: str | None) -> None:
        """Force the matplotlib NavigationToolbar mode to ``target``.

        ``target`` is one of ``'pan'``, ``'zoom'`` or ``None`` (idle).
        ``NavigationToolbar2.pan()`` and ``zoom()`` are toggles, so we
        inspect ``nav.mode`` and only call them when needed. This avoids
        leaving matplotlib in pan/zoom mode when switching to the custom
        Add Zoom mode.
        """
        current = str(self._nav.mode)  # '', 'pan/zoom', 'zoom rect'
        if target == "pan":
            if current == "zoom rect":
                self._nav.zoom()           # turn off zoom rect
            if str(self._nav.mode) != "pan/zoom":
                self._nav.pan()            # turn on pan
        elif target == "zoom":
            if current == "pan/zoom":
                self._nav.pan()            # turn off pan
            if str(self._nav.mode) != "zoom rect":
                self._nav.zoom()           # turn on zoom rect
        else:  # idle
            if current == "pan/zoom":
                self._nav.pan()
            elif current == "zoom rect":
                self._nav.zoom()

    def _on_pan(self, checked: bool):
        if checked:
            self._btn_zoom.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._deactivate_pending_selectors()
            self._deactivate_cursor()
            self._sync_nav_mode("pan")
        else:
            self._sync_nav_mode(None)

    def _on_zoom(self, checked: bool):
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._deactivate_pending_selectors()
            self._deactivate_cursor()
            self._sync_nav_mode("zoom")
        else:
            self._sync_nav_mode(None)

    def _on_add_zoom(self, checked: bool):
        """Toggle the zoom-inset selection mode."""
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._sync_nav_mode(None)       # ensure matplotlib pan/zoom is OFF
            self._deactivate_cursor()
            self._activate_pending_selectors()
        else:
            self._deactivate_pending_selectors()

    def _on_cursor(self, checked: bool):
        """Toggle the data-cursor mode."""
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_zoom.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._sync_nav_mode(None)
            self._deactivate_pending_selectors()
            self._activate_cursor()
        else:
            self._deactivate_cursor()

    # -- Zoom-inset helpers ------------------------------------------------

    def _activate_pending_selectors(self) -> None:
        """Create one temporary RectangleSelector per subplot axis.

        Existing zoom-inset selectors are paused so they do not steal events
        while a new region is being drawn. Each pending selector also vetoes
        clicks that fall inside any existing inset.
        """
        self._deactivate_pending_selectors()

        # Pause all existing inset selectors
        for _z in self._zoom_insets:
            _z.pause_selector()

        for idx, ax in enumerate(self._axes):
            sel = RectangleSelector(
                ax,
                lambda ec, er, i=idx: self._on_temp_select(i, ec, er),
                useblit=False,
                button=[1],
                interactive=False,
                minspanx=5, minspany=5,
                spancoords="pixels",
                props=dict(
                    edgecolor="steelblue",
                    facecolor=(0.27, 0.51, 0.71, 0.06),
                    linestyle="--", linewidth=1,
                ),
            )
            ZoomInset.install_selector_veto(sel)
            self._pending_selectors.append(sel)

    def _deactivate_pending_selectors(self) -> None:
        for sel in self._pending_selectors:
            try:
                sel.set_active(False)
                sel.clear()
            except Exception:
                pass
        self._pending_selectors.clear()
        # Resume all existing inset selectors
        for _z in self._zoom_insets:
            _z.resume_selector()

    def _on_temp_select(self, ax_idx: int, eclick, erelease) -> None:
        """Called when the user finishes drawing the initial selection region."""
        if eclick.xdata is None or erelease.xdata is None:
            return
        x0, x1 = sorted([eclick.xdata, erelease.xdata])
        y0, y1 = sorted([eclick.ydata, erelease.ydata])
        if x1 - x0 < 1e-10 or y1 - y0 < 1e-10:
            return

        self._deactivate_pending_selectors()
        self._btn_add_zoom.setChecked(False)

        ax_parent = self._axes[ax_idx]

        # Identify the curves that belong to this subplot
        groups: dict[str, list] = {}
        for c in self._curves:
            groups.setdefault(c.unit or "Value", []).append(c)
        group_curves = list(groups.values())[ax_idx] if ax_idx < len(groups) else []

        def _remove_cb(z: ZoomInset) -> None:
            if z in self._zoom_insets:
                self._zoom_insets.remove(z)

        zoom = ZoomInset(
            parent_ax=ax_parent,
            curves=group_curves,
            x0=x0, x1=x1,
            y0=y0, y1=y1,
            dark=self._dark,
            on_remove=_remove_cb,
        )
        self._zoom_insets.append(zoom)
        self._canvas.draw_idle()

    def _on_save(self):
        # Kept as a no-op for backwards compatibility (button was removed).
        self._nav.save_figure()

    # -- Data cursor ---------------------------------------------------

    def _activate_cursor(self) -> None:
        """Install per-axes crosshair artists and connect the motion handler."""
        self._deactivate_cursor()
        fg = CANVAS_FG_DARK if self._dark else CANVAS_FG_LIGHT
        for ax in self._axes:
            vline = ax.axvline(x=0, color=fg, linestyle=":", linewidth=0.8, visible=False)
            hline = ax.axhline(y=0, color=fg, linestyle=":", linewidth=0.8, visible=False)
            text = ax.annotate(
                "", xy=(0, 0), xycoords="data",
                xytext=(8, 8), textcoords="offset points",
                fontsize=8, color=fg, visible=False, zorder=12,
                bbox=dict(boxstyle="round,pad=0.25",
                          facecolor=(CANVAS_BG_DARK if self._dark else CANVAS_BG_LIGHT),
                          edgecolor=fg, linewidth=0.6, alpha=0.85),
            )
            self._cursor_artists[ax] = {"vline": vline, "hline": hline, "text": text}
        self._cid_cursor_motion = self._canvas.mpl_connect(
            "motion_notify_event", self._on_cursor_motion
        )
        self._canvas.draw_idle()

    def _deactivate_cursor(self) -> None:
        """Remove cursor artists and disconnect the motion handler."""
        if self._cid_cursor_motion is not None:
            try:
                self._canvas.mpl_disconnect(self._cid_cursor_motion)
            except Exception:
                pass
            self._cid_cursor_motion = None
        for artists in self._cursor_artists.values():
            for a in artists.values():
                try:
                    a.remove()
                except Exception:
                    pass
        self._cursor_artists.clear()
        self._canvas.draw_idle()

    def _on_cursor_motion(self, event) -> None:
        ax = event.inaxes
        if ax is None or ax not in self._cursor_artists or event.xdata is None:
            # Hide all cursors when the pointer is outside any tracked axes
            changed = False
            for artists in self._cursor_artists.values():
                for a in artists.values():
                    if a.get_visible():
                        a.set_visible(False)
                        changed = True
            if changed:
                self._canvas.draw_idle()
            return

        # Find the curves that belong to this axes
        ax_idx = self._axes.index(ax)
        groups: dict[str, list] = {}
        for c in self._curves:
            groups.setdefault(c.unit or "Value", []).append(c)
        group_curves = list(groups.values())[ax_idx] if ax_idx < len(groups) else []
        if not group_curves:
            return

        # Snap to the closest curve in y at the mouse x
        x_target = event.xdata
        y_target = event.ydata
        best_x = best_y = None
        best_dy = float("inf")
        for c in group_curves:
            T = c.T
            i = int(np.searchsorted(T, x_target))
            if i >= len(T):
                i = len(T) - 1
            elif i > 0 and abs(T[i - 1] - x_target) < abs(T[i] - x_target):
                i -= 1
            cx = float(T[i])
            cy = float(c.data[i])
            dy = abs(cy - y_target) if y_target is not None else 0.0
            if dy < best_dy:
                best_dy = dy
                best_x, best_y = cx, cy

        if best_x is None:
            return

        # Hide cursors on other axes
        for other_ax, artists in self._cursor_artists.items():
            if other_ax is ax:
                continue
            for a in artists.values():
                if a.get_visible():
                    a.set_visible(False)

        artists = self._cursor_artists[ax]
        artists["vline"].set_xdata([best_x])
        artists["hline"].set_ydata([best_y])
        artists["text"].xy = (best_x, best_y)
        artists["text"].set_text(f"x={best_x:.4g}\ny={best_y:.4g}")
        for a in artists.values():
            a.set_visible(True)
        self._canvas.draw_idle()

    def update_plot(self, curves: list[CurveItem]):
        """Clear figure, create one subplot per unit group, plot curves."""
        self._curves = curves
        self._deactivate_pending_selectors()
        self._deactivate_cursor()
        self._btn_add_zoom.setChecked(False)
        self._btn_cursor.setChecked(False)
        # Remove each inset cleanly (disconnects events, clears selector)
        for _z in list(self._zoom_insets):
            _z.remove()
        self._zoom_insets.clear()
        self._figure.clear()

        if not curves:
            ax0 = self._figure.add_subplot(111)
            self._axes = [ax0]
            self.set_dark(self._dark)
            self._canvas.draw_idle()
            return

        # Group curves by unit (= shared Y-label); preserves insertion order
        groups: dict[str, list] = {}
        for c in curves:
            groups.setdefault(c.unit or "Value", []).append(c)

        n = len(groups)
        self._axes = []
        for idx, (ylabel, group) in enumerate(groups.items()):
            ax = self._figure.add_subplot(n, 1, idx + 1)
            self._axes.append(ax)
            for c in group:
                ax.plot(c.T, c.data, color=c.color, label=c.label)
            ax.set_ylabel(ylabel)
            ax.legend(fontsize="small", loc="best")
            ax.grid(True, linestyle="--", alpha=0.5)
            if idx == n - 1:
                ax.set_xlabel("Time [s]")

        self.set_dark(self._dark)
        self._canvas.draw_idle()

    def set_dark(self, enabled: bool):
        """Switch the existing Figure between dark and light appearance."""
        self._dark = enabled
        bg = CANVAS_BG_DARK  if enabled else CANVAS_BG_LIGHT
        fg = CANVAS_FG_DARK  if enabled else CANVAS_FG_LIGHT
        self._figure.set_facecolor(bg)
        for ax in self._axes:
            ax.set_facecolor(bg)
            ax.tick_params(colors=fg, which="both")
            ax.xaxis.label.set_color(fg)
            ax.yaxis.label.set_color(fg)
            ax.title.set_color(fg)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg)
            # Re-skin legend so text/frame match the theme
            leg = ax.get_legend()
            if leg is not None:
                for txt in leg.get_texts():
                    txt.set_color(fg)
                frame = leg.get_frame()
                frame.set_facecolor(bg)
                frame.set_edgecolor(fg)
        # Re-skin existing cursor artists, if any
        for artists in self._cursor_artists.values():
            artists["vline"].set_color(fg)
            artists["hline"].set_color(fg)
            artists["text"].set_color(fg)
            box = artists["text"].get_bbox_patch()
            if box is not None:
                box.set_facecolor(bg)
                box.set_edgecolor(fg)
        for _z in self._zoom_insets:
            _z.apply_theme(enabled)
        self._canvas.draw_idle()

    def _on_click(self, event):
        """Convert a matplotlib click to the nearest time-step index and emit it."""
        if (self._btn_pan.isChecked() or self._btn_zoom.isChecked()
                or self._btn_add_zoom.isChecked() or self._btn_cursor.isChecked()):
            return
        if event.button != 1:
            return
        if not self._curves or event.inaxes not in self._axes or event.xdata is None:
            return
        t_click = event.xdata
        T = self._curves[0].T
        step = int(np.argmin(np.abs(T - t_click)))
        self.step_requested.emit(step)
