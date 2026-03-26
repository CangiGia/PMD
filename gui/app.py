"""PostProcessor — entry-point façade for the PMD GUI."""

import sys

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow
from .models import Session


class PostProcessor:
    """Launch a post-processing GUI for a solved PMD model.

    Parameters
    ----------
    model : PlanarMultibodyModel
        A solved model instance.
    T : ndarray
        Time vector, shape (nSteps,).
    uT : ndarray
        State matrix, shape (nSteps, 2*nB3).
    name : str, optional
        Label for this simulation session shown in the UI.
    """

    def __init__(self, model, T, uT, name=None):
        model._distribute_results(T, uT)
        self._session = Session(model, T, uT, name=name)

    def show(self):
        """Create (or reuse) a QApplication and display the main window."""
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        self._window = MainWindow(self._session)
        self._window.show()
        app.exec()
