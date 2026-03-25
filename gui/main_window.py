"""MainWindow — primary QMainWindow shell for the PMD PostProcessor."""

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QLabel,
    QMainWindow,
    QSplitter,
    QVBoxLayout,
    QWidget,
)


class MainWindow(QMainWindow):
    """Top-level window: menu bar, left panel | central area, status bar.

    Parameters
    ----------
    model : PlanarMultibodyModel
        The solved model whose results are displayed.
    T : ndarray
        Time vector, shape (nSteps,).
    parent : QWidget or None
        Optional parent widget.
    """

    def __init__(self, model, T, parent=None):
        super().__init__(parent)
        self._model = model
        self._T = T

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
    # Central area (splitter: left panel | central placeholder)
    # ------------------------------------------------------------------

    def _build_central_area(self):
        splitter = QSplitter(Qt.Orientation.Horizontal, self)

        # Left panel — placeholder
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.addWidget(QLabel("Result browser (placeholder)"))
        left_layout.addStretch()
        left_panel.setMinimumWidth(200)
        left_panel.setMaximumWidth(350)

        # Central area — placeholder
        central_widget = QWidget()
        central_layout = QVBoxLayout(central_widget)
        central_layout.addWidget(QLabel("Plot area (placeholder)"))
        central_layout.addStretch()

        splitter.addWidget(left_panel)
        splitter.addWidget(central_widget)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.setCentralWidget(splitter)

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------

    def _build_status_bar(self):
        n_bodies = len(self._model.Bodies)
        n_joints = len(self._model.Joints)
        n_steps = len(self._T)
        t_range = f"{self._T[0]:.4g} – {self._T[-1]:.4g} s"

        self.statusBar().showMessage(
            f"Bodies: {n_bodies}  |  Joints: {n_joints}  |  "
            f"Steps: {n_steps}  |  T: {t_range}"
        )
