"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QToolBar, QToolButton, QVBoxLayout, QWidget

from ..models import CurveItem
from .. import icons as _icons


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
        self._nav.pan()

    def _on_zoom(self, checked: bool):
        if checked:
            self._btn_pan.setChecked(False)
        self._nav.zoom()

    def _on_save(self):
        self._nav.save_figure()

    def update_plot(self, curves: list[CurveItem]):
        """Clear figure, create one subplot per unit group, plot curves."""
        self._curves = curves
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
        bg = "#2b2b2b" if enabled else "white"
        fg = "#cccccc" if enabled else "black"
        self._figure.set_facecolor(bg)
        for ax in self._axes:
            ax.set_facecolor(bg)
            ax.tick_params(colors=fg, which="both")
            ax.xaxis.label.set_color(fg)
            ax.yaxis.label.set_color(fg)
            ax.title.set_color(fg)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg)
        self._canvas.draw_idle()

    def _on_click(self, event):
        """Convert a matplotlib click to the nearest time-step index and emit it."""
        if self._btn_pan.isChecked() or self._btn_zoom.isChecked():
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
