"""MainWindow — primary QMainWindow shell for the PMD PostProcessor."""

import csv
import logging

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from .models import build_curves
from .panels import FilterPanel, ResultSetPanel, SimulationPanel
from .widgets import AnimationCanvas, PlotCanvas

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """Top-level window: menu bar, left panel | central area, status bar.

    Parameters
    ----------
    sessions : list[Session]
        One or more loaded simulation sessions to display.
    parent : QWidget or None
        Optional parent widget.
    """

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        self._sessions = sessions
        self._selection = []

        self.setWindowTitle("PMD PostProcessor")
        self.resize(1200, 700)

        self._build_menu_bar()
        self._build_central_area()
        self._build_status_bar()

    # ------------------------------------------------------------------
    # Menu bar
    # ------------------------------------------------------------------

    def _build_menu_bar(self):
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("&File")
        file_menu.addAction("&Close", self.close)

        view_menu = menu_bar.addMenu("&View")
        self._anim_action = QAction("Animation Pane", self, checkable=True, checked=False)
        self._anim_action.toggled.connect(self._on_toggle_animation)
        view_menu.addAction(self._anim_action)
        view_menu.addSeparator()
        self._dark_action = QAction("Dark Theme", self, checkable=True, checked=False)
        self._dark_action.toggled.connect(self._on_toggle_theme)
        view_menu.addAction(self._dark_action)

        export_menu = menu_bar.addMenu("E&xport")
        export_menu.addAction("Save Plot as Image…", self._on_export_plot)
        export_menu.addAction("Export Curves to CSV…", self._on_export_csv)

    # ------------------------------------------------------------------
    # Central area (splitter: SimulationPanel | FilterPanel + plot)
    # ------------------------------------------------------------------

    def _build_central_area(self):
        splitter = QSplitter(Qt.Orientation.Horizontal, self)

        # Left panel — SimulationPanel
        self._sim_panel = SimulationPanel(self._sessions)
        self._sim_panel.setMinimumWidth(200)
        self._sim_panel.setMaximumWidth(350)
        self._sim_panel.selection_changed.connect(self._on_selection_changed)

        # Right side — FilterPanel on top, vertical splitter (plot | animation), ResultSetPanel bottom
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self._filter_panel = FilterPanel()
        self._filter_panel.add_curves_requested.connect(self._on_add_curves)
        right_layout.addWidget(self._filter_panel)

        viz_splitter = QSplitter(Qt.Orientation.Vertical)
        self._plot_canvas = PlotCanvas()
        viz_splitter.addWidget(self._plot_canvas)
        self._anim_canvas = AnimationCanvas(self._sessions)
        viz_splitter.addWidget(self._anim_canvas)
        self._anim_canvas.setVisible(False)
        self._plot_canvas.step_requested.connect(self._anim_canvas.set_step)
        viz_splitter.setStretchFactor(0, 1)
        viz_splitter.setStretchFactor(1, 1)
        right_layout.addWidget(viz_splitter, stretch=1)

        self._result_set_panel = ResultSetPanel()
        self._result_set_panel.setMaximumHeight(200)
        self._result_set_panel.curves_changed.connect(self._on_curves_changed)
        right_layout.addWidget(self._result_set_panel)

        splitter.addWidget(self._sim_panel)
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.setCentralWidget(splitter)

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------

    def _build_status_bar(self):
        n_bodies = sum(len(s.model.Bodies) for s in self._sessions)
        n_joints = sum(len(s.model.Joints) for s in self._sessions)
        n_sessions = len(self._sessions)

        parts = [f"Sessions: {n_sessions}",
                 f"Bodies: {n_bodies}",
                 f"Joints: {n_joints}"]

        if n_sessions == 1:
            T = self._sessions[0].T
            n_steps = len(T)
            t_range = f"{T[0]:.4g} – {T[-1]:.4g} s"
            parts += [f"Steps: {n_steps}", f"T: {t_range}"]

        self._status_base = "  |  ".join(parts)
        self.statusBar().showMessage(self._status_base)

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_selection_changed(self, selection):
        """Receives list of checked descriptor dicts from SimulationPanel."""
        self._selection = selection
        self._filter_panel.update_from_selection(selection)
        msg = self._status_base + f"  |  Selected: {len(selection)}"
        self.statusBar().showMessage(msg)
        logger.debug("Selection changed: %s", [d["label"] for d in selection])

    def _on_add_curves(self, category, component, selection):
        """Build CurveItems and add them to the ResultSetPanel."""
        curves = build_curves(category, component, selection)
        if curves:
            self._result_set_panel.add_curves(curves)
        logger.debug(
            "Add curves: category=%s, component=%s, added=%d",
            category, component, len(curves),
        )

    def _on_curves_changed(self):
        """Refresh the plot with currently visible curves."""
        visible = self._result_set_panel.visible_curves()
        self._plot_canvas.update_plot(visible)
        logger.debug("Curves changed: %d visible", len(visible))

    def _on_toggle_animation(self, checked: bool):
        """Show/hide the AnimationCanvas from the View menu."""
        self._anim_canvas.setVisible(checked)

    def _on_toggle_theme(self, checked: bool):
        """Switch the PlotCanvas between dark and light appearance."""
        self._plot_canvas.set_dark(checked)
        visible = self._result_set_panel.visible_curves()
        self._plot_canvas.update_plot(visible)

    def _on_export_plot(self):
        """Save the current plot figure to PNG/SVG/PDF."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "",
            "PNG Image (*.png);;SVG Vector (*.svg);;PDF Document (*.pdf)"
        )
        if path:
            self._plot_canvas._figure.savefig(path, bbox_inches="tight", dpi=150)

    def _on_export_csv(self):
        """Export all visible curves to a single CSV file."""
        curves = self._result_set_panel.visible_curves()
        if not curves:
            QMessageBox.information(self, "Export CSV", "No visible curves to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Curves", "", "CSV file (*.csv)"
        )
        if not path:
            return
        T_ref = curves[0].T
        header = ["time_s"] + [c.label for c in curves]
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for i, t in enumerate(T_ref):
                row = [t] + [float(c.data[i]) if i < len(c.data) else "" for c in curves]
                writer.writerow(row)
