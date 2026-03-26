"""PostProcessor — entry-point façade for the PMD GUI."""

import sys

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow
from .models import Session


class PostProcessor:
    """Launch a post-processing GUI for one or more solved PMD models.

    Parameters
    ----------
    model : PlanarMultibodyModel, optional
        A solved model instance.  When provided, the first session is
        created immediately.
    T : ndarray, optional
        Time vector, shape (nSteps,).
    uT : ndarray, optional
        State matrix, shape (nSteps, 2*nB3).
    name : str, optional
        Label for this simulation session shown in the UI.
    """

    def __init__(self, model=None, T=None, uT=None, name=None):
        self._sessions: list[Session] = []
        if model is not None:
            self.add(model, T, uT, name=name)

    def add(self, model, T, uT, name=None):
        """Add a solved simulation to the post-processor.

        Returns *self* so calls can be chained.
        """
        model._distribute_results(T, uT)
        self._sessions.append(Session(model, T, uT, name=name))
        return self

    def show(self):
        """Create (or reuse) a QApplication and display the main window."""
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        self._window = MainWindow(self._sessions)
        self._window.show()
        app.exec()
