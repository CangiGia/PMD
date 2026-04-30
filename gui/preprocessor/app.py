"""PreProcessor — entry-point façade for the PMD model-builder GUI."""

from __future__ import annotations

import sys

from PySide6.QtWidgets import QApplication

from ..style import apply_light_theme
from .main_window import PreProcessorWindow
from .models import ModelSpec


class PreProcessor:
    """Launch the PMD pre-processor (interactive model builder).

    Parameters
    ----------
    spec : ModelSpec, optional
        Existing model specification to load. If ``None``, an empty
        project is created.
    """

    def __init__(self, spec: ModelSpec | None = None):
        self._spec = spec if spec is not None else ModelSpec()

    def show(self):
        """Create (or reuse) a QApplication and display the main window."""
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        apply_light_theme(app)
        self._window = PreProcessorWindow(self._spec)
        self._window.show()
        app.exec()
