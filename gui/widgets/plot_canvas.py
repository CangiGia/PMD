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
        self._ax = self._figure.add_subplot(111)
        self._canvas = FigureCanvasQTAgg(self._figure)
        self._toolbar = NavigationToolbar2QT(self._canvas, self)

        self._curves = []
        self._dark = False
        self._vline = self._ax.axvline(x=0, color="gray", linestyle="--",
                                        linewidth=0.8, visible=False)
        self._canvas.mpl_connect("button_press_event", self._on_click)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

    def update_plot(self, curves: list[CurveItem]):
        """Clear axes, plot each curve, refresh canvas."""
        self._curves = curves
        ax = self._ax
        ax.clear()
        self._vline = self._ax.axvline(x=0, color="gray", linestyle="--",
                                        linewidth=0.8, visible=False)

        for curve in curves:
            ax.plot(curve.T, curve.data, color=curve.color, label=curve.label)

        if curves:
            ax.set_xlabel("Time [s]")
            # Y-label: use shared unit if all curves agree, else leave empty
            units = {c.unit for c in curves if c.unit}
            ax.set_ylabel(next(iter(units)) if len(units) == 1 else "")
            ax.legend(fontsize="small", loc="best")
            ax.grid(True, linestyle="--", alpha=0.5)

        # apply dark/light style after ax.clear() wipes all colors
        self.set_dark(self._dark)

        self._canvas.draw_idle()

    def set_dark(self, enabled: bool):
        """Switch the existing Figure between dark and light appearance."""
        self._dark = enabled
        bg = "#2b2b2b" if enabled else "white"
        fg = "#cccccc" if enabled else "black"
        grid_c = "#444444" if enabled else "#cccccc"
        self._figure.set_facecolor(bg)
        ax = self._ax
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
        if event.inaxes is not self._ax:
            return
        if not self._curves:
            return
        t_click = event.xdata
        T = self._curves[0].T
        step = int(np.argmin(np.abs(T - t_click)))
        self._vline.set_xdata([t_click])
        self._vline.set_visible(True)
        self._canvas.draw_idle()
        self.step_requested.emit(step)
