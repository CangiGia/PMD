"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtWidgets import QVBoxLayout, QWidget

from .curve_item import CurveItem


class PlotCanvas(QWidget):
    """Embeds a matplotlib Figure + NavigationToolbar.

    Public API
    ----------
    update_plot(curves: list[CurveItem])
        Clear the axes and redraw only the given curves.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self._figure = Figure(tight_layout=True)
        self._ax = self._figure.add_subplot(111)
        self._canvas = FigureCanvasQTAgg(self._figure)
        self._toolbar = NavigationToolbar2QT(self._canvas, self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

    def update_plot(self, curves: list[CurveItem]):
        """Clear axes, plot each curve, refresh canvas."""
        ax = self._ax
        ax.clear()

        for curve in curves:
            ax.plot(curve.T, curve.data, color=curve.color, label=curve.label)

        if curves:
            ax.set_xlabel("Time [s]")
            ax.legend(fontsize="small", loc="best")
            ax.grid(True, linestyle="--", alpha=0.5)

        self._canvas.draw_idle()
