"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QToolBar, QVBoxLayout, QWidget

from ..models import CurveItem

# Matplotlib figure background colours (match QSS theme tokens)
_MPL_BG_LIGHT = "#f0f0f0"
_MPL_BG_DARK  = "#2a2a2a"
_MPL_FG_LIGHT = "#1a1a1a"
_MPL_FG_DARK  = "#e0e0e0"


class PlotCanvas(QWidget):
    """Embeds a matplotlib Figure + custom QToolBar.

    Public API
    ----------
    update_plot(curves: list[CurveItem])
        Clear the axes and redraw only the given curves.
    set_dark(enabled: bool)
        Switch between dark and light appearance.
    """

    step_requested = Signal(int)  # emitted on click: nearest time-step index

    def __init__(self, parent=None):
        super().__init__(parent)

        self._figure = Figure(tight_layout=True)
        self._canvas = FigureCanvasQTAgg(self._figure)

        # Hidden backend toolbar — manages pan/zoom mode state internally
        self._nav = NavigationToolbar2QT(self._canvas, self)
        self._nav.setVisible(False)

        self._curves: list = []
        self._dark = False
        ax0 = self._figure.add_subplot(111)
        self._axes: list = [ax0]
        self._vlines: list = [ax0.axvline(x=0, color="gray", linestyle="--",
                                           linewidth=0.8, visible=False)]
        self._canvas.mpl_connect("button_press_event", self._on_click)

        self._toolbar = self._build_toolbar()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

    # ------------------------------------------------------------------
    # Toolbar construction
    # ------------------------------------------------------------------

    def _build_toolbar(self) -> QToolBar:
        tb = QToolBar()
        tb.setMovable(False)
        tb.setFloatable(False)

        self._act_home = tb.addAction("⌂  Home")
        self._act_home.setToolTip("Reset to original view")
        self._act_home.triggered.connect(self._on_home)

        self._act_back = tb.addAction("←  Back")
        self._act_back.setToolTip("Back to previous view")
        self._act_back.triggered.connect(self._nav.back)

        self._act_fwd = tb.addAction("→  Forward")
        self._act_fwd.setToolTip("Forward to next view")
        self._act_fwd.triggered.connect(self._nav.forward)

        tb.addSeparator()

        self._act_pan = tb.addAction("✥  Pan")
        self._act_pan.setToolTip("Pan axes (left drag)")
        self._act_pan.setCheckable(True)
        self._act_pan.triggered.connect(self._on_pan)

        self._act_zoom = tb.addAction("⌕  Zoom")
        self._act_zoom.setToolTip("Zoom to rectangle")
        self._act_zoom.setCheckable(True)
        self._act_zoom.triggered.connect(self._on_zoom)

        tb.addSeparator()

        self._act_save = tb.addAction("⬇  Save")
        self._act_save.setToolTip("Save figure to file")
        self._act_save.triggered.connect(self._nav.save_figure)

        return tb

    # ------------------------------------------------------------------
    # Toolbar slots
    # ------------------------------------------------------------------

    def _on_home(self):
        self._nav.home()
        self._sync_mode_buttons()

    def _on_pan(self):
        self._nav.pan()
        self._sync_mode_buttons()

    def _on_zoom(self):
        self._nav.zoom()
        self._sync_mode_buttons()

    def _sync_mode_buttons(self):
        """Keep Pan/Zoom checked state in sync with matplotlib's internal mode."""
        mode = str(self._nav.mode).lower()
        self._act_pan.setChecked("pan" in mode)
        self._act_zoom.setChecked("zoom" in mode)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def update_plot(self, curves: list[CurveItem]):
        """Clear figure, create one subplot per unit group, plot curves."""
        self._curves = curves
        self._figure.clear()
        self._vlines = []
        self._sync_mode_buttons()

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
        bg = _MPL_BG_DARK  if enabled else _MPL_BG_LIGHT
        fg = _MPL_FG_DARK  if enabled else _MPL_FG_LIGHT
        grid_color = "#444444" if enabled else "#cccccc"

        self._figure.set_facecolor(bg)
        for ax in self._axes:
            ax.set_facecolor(bg)
            ax.tick_params(colors=fg, which="both")
            ax.xaxis.label.set_color(fg)
            ax.yaxis.label.set_color(fg)
            ax.title.set_color(fg)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg)
            # Update existing grid lines colour
            ax.grid(True, linestyle="--", alpha=0.5, color=grid_color)
            # Legend text colour
            legend = ax.get_legend()
            if legend:
                for text in legend.get_texts():
                    text.set_color(fg)
                legend.get_frame().set_facecolor(bg)
                legend.get_frame().set_edgecolor(fg)

        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # Click → seek
    # ------------------------------------------------------------------

    def _on_click(self, event):
        """Convert a matplotlib click to the nearest time-step index and emit it."""
        if not self._curves or event.inaxes not in self._axes or event.xdata is None:
            return
        # Only seek on plain left-click (not during pan/zoom)
        mode = str(self._nav.mode).lower()
        if "pan" in mode or "zoom" in mode:
            return
        t_click = event.xdata
        T = self._curves[0].T
        step = int(np.argmin(np.abs(T - t_click)))
        for vline in self._vlines:
            vline.set_xdata([t_click])
            vline.set_visible(True)
        self._canvas.draw_idle()
        self.step_requested.emit(step)

