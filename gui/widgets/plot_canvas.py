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
        self._vlines: list = [ax0.axvline(x=0, color="gray", linestyle="--",
                                           linewidth=0.8, visible=False)]
        self._canvas.mpl_connect("button_press_event", self._on_click)

        # Zoom-inset state
        self._zoom_insets: list[ZoomInset] = []
        self._pending_selectors: list = []  # temporary selectors during "Add Zoom" mode

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
        self._btn_add_zoom = self._make_btn(tb, "Add Zoom", "mdi6.magnify-plus-outline",    self._on_add_zoom, checkable=True)
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
        self._btn_home.setIcon(    _icons.icon("mdi6.home"))
        self._btn_back.setIcon(    _icons.icon("mdi6.arrow-left"))
        self._btn_fwd.setIcon(     _icons.icon("mdi6.arrow-right"))
        self._btn_pan.setIcon(     _icons.icon("mdi6.hand-back-right-outline"))
        self._btn_zoom.setIcon(    _icons.icon("mdi6.magnify"))
        self._btn_add_zoom.setIcon(_icons.icon("mdi6.magnify-plus-outline"))
        self._btn_save.setIcon(    _icons.icon("mdi6.content-save"))

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

    def _on_pan(self, checked: bool):
        if checked:
            self._btn_zoom.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._deactivate_pending_selectors()
        self._nav.pan()

    def _on_zoom(self, checked: bool):
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._deactivate_pending_selectors()
        self._nav.zoom()

    def _on_add_zoom(self, checked: bool):
        """Toggle the zoom-inset selection mode."""
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_zoom.setChecked(False)
            self._activate_pending_selectors()
        else:
            self._deactivate_pending_selectors()

    # -- Zoom-inset helpers ------------------------------------------------

    def _activate_pending_selectors(self) -> None:
        """Create one temporary RectangleSelector per subplot axis."""
        self._deactivate_pending_selectors()
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
            self._pending_selectors.append(sel)

    def _deactivate_pending_selectors(self) -> None:
        for sel in self._pending_selectors:
            try:
                sel.set_active(False)
                sel.clear()
            except Exception:
                pass
        self._pending_selectors.clear()

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
        self._nav.save_figure()

    def update_plot(self, curves: list[CurveItem]):
        """Clear figure, create one subplot per unit group, plot curves."""
        self._curves = curves
        self._deactivate_pending_selectors()
        self._btn_add_zoom.setChecked(False)
        # Remove each inset cleanly (disconnects events, clears selector)
        for _z in list(self._zoom_insets):
            _z.remove()
        self._zoom_insets.clear()
        self._figure.clear()
        self._vlines = []

        if not curves:
            ax0 = self._figure.add_subplot(111)
            self._axes = [ax0]
            self._vlines = [ax0.axvline(x=0, color="gray", linestyle="--",
                                         linewidth=0.8, visible=False)]
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
            vline = ax.axvline(x=0, color="gray", linestyle="--",
                               linewidth=0.8, visible=False)
            self._vlines.append(vline)
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
        for _z in self._zoom_insets:
            _z.apply_theme(enabled)
        self._canvas.draw_idle()

    def _on_click(self, event):
        """Convert a matplotlib click to the nearest time-step index and emit it."""
        if self._btn_pan.isChecked() or self._btn_zoom.isChecked() or self._btn_add_zoom.isChecked():
            return
        if event.button != 1:
            return
        if not self._curves or event.inaxes not in self._axes or event.xdata is None:
            return
        t_click = event.xdata
        T = self._curves[0].T
        step = int(np.argmin(np.abs(T - t_click)))
        for vline in self._vlines:
            vline.set_xdata([t_click])
            vline.set_visible(True)
        self._canvas.draw_idle()
        self.step_requested.emit(step)
