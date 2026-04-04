"""FilterPanel — Add Curves action bar."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QPushButton, QWidget


class FilterPanel(QWidget):
    """Thin bar with only an *Add Curves* button.

    The curves to add are determined by what is checked in NavigationPanel.

    Signals
    -------
    add_curves_requested()
        Emitted when the user clicks *Add Curves*.
    """

    add_curves_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        self._add_btn = QPushButton("Add Curves")
        self._add_btn.clicked.connect(self.add_curves_requested)
        layout.addWidget(self._add_btn)
        layout.addStretch()
