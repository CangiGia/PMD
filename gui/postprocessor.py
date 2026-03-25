"""PostProcessor — entry-point façade for the PMD GUI."""

import sys

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow


class PostProcessor:
    """Launch a post-processing GUI for a solved PMD model.

    Parameters
    ----------
    model : PlanarMultibodyModel
        A solved (or unsolved) model instance.
    T : ndarray
        Time vector, shape (nSteps,).
    uT : ndarray
        State matrix, shape (nSteps, 2*nB3).
    """

    def __init__(self, model, T, uT):
        self._model = model
        self._T = T
        self._uT = uT
        # Ensure results are distributed into bodies/joints
        model._distribute_results(T, uT)

    def show(self):
        """Create (or reuse) a QApplication and display the main window."""
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        self._window = MainWindow(self._model, self._T)
        self._window.show()
        app.exec()
