"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..models import CurveItem


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
        self._toolbar = NavigationToolbar2QT(self._canvas, self)

        self._curves: list = []
        self._dark = False
        ax0 = self._figure.add_subplot(111)
        self._axes: list = [ax0]
        self._vlines: list = [ax0.axvline(x=0, color="gray", linestyle="--",
                                           linewidth=0.8, visible=False)]
        self._canvas.mpl_connect("button_press_event", self._on_click)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

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
