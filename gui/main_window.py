"""MainWindow — primary QMainWindow shell for the PMD PostProcessor."""

import logging

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QLabel,
    QMainWindow,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from .filter_panel import FilterPanel
from .simulation_panel import SimulationPanel

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """Top-level window: menu bar, left panel | central area, status bar.

    Parameters
    ----------
    session : Session
        The loaded simulation session to display.
    parent : QWidget or None
        Optional parent widget.
    """

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self._session = session
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

    # ------------------------------------------------------------------
    # Central area (splitter: SimulationPanel | FilterPanel + plot)
    # ------------------------------------------------------------------

    def _build_central_area(self):
        splitter = QSplitter(Qt.Orientation.Horizontal, self)

        # Left panel — SimulationPanel
        self._sim_panel = SimulationPanel(self._session)
        self._sim_panel.setMinimumWidth(200)
        self._sim_panel.setMaximumWidth(350)
        self._sim_panel.selection_changed.connect(self._on_selection_changed)

        # Right side — FilterPanel on top, plot placeholder below
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self._filter_panel = FilterPanel()
        self._filter_panel.add_curves_requested.connect(self._on_add_curves)
        right_layout.addWidget(self._filter_panel)

        self._plot_placeholder = QLabel("Plot area (placeholder)")
        self._plot_placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        right_layout.addWidget(self._plot_placeholder, stretch=1)

        splitter.addWidget(self._sim_panel)
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.setCentralWidget(splitter)

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------

    def _build_status_bar(self):
        model = self._session.model
        T = self._session.T
        n_bodies = len(model.Bodies)
        n_joints = len(model.Joints)
        n_steps = len(T)
        t_range = f"{T[0]:.4g} – {T[-1]:.4g} s"

        self._status_base = (
            f"Bodies: {n_bodies}  |  Joints: {n_joints}  |  "
            f"Steps: {n_steps}  |  T: {t_range}"
        )
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
        """Stub — will forward to PlotCanvas / ResultSetPanel in S4/S5."""
        logger.debug(
            "Add curves: category=%s, component=%s, items=%s",
            category, component, [d["label"] for d in selection],
        )
