"""UnitsToolbar — GUI toolbar for switching display unit system."""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QWidget

from pmd.core.units import UnitSystem


class UnitsToolbar(QWidget):
    """Horizontal toolbar that lets the user pick display units.

    Emits :attr:`units_changed` whenever any combo-box selection changes.

    Layout::

        Display as:  Length [m ▾]  Angle [rad ▾]  Force [N ▾]
    """

    units_changed = Signal(object)  # payload: UnitSystem

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("units_toolbar")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(8)

        layout.addWidget(QLabel("Display as:"))

        # --- Length ---
        layout.addWidget(QLabel("Length"))
        self._length_cb = QComboBox()
        self._length_cb.addItems(["m", "mm", "cm", "in", "ft"])
        layout.addWidget(self._length_cb)

        # --- Angle ---
        layout.addWidget(QLabel("Angle"))
        self._angle_cb = QComboBox()
        self._angle_cb.addItems(["rad", "deg"])
        layout.addWidget(self._angle_cb)

        # --- Force ---
        layout.addWidget(QLabel("Force"))
        self._force_cb = QComboBox()
        self._force_cb.addItems(["N", "kN", "lbf"])
        layout.addWidget(self._force_cb)

        layout.addStretch()

        self._length_cb.currentTextChanged.connect(self._emit)
        self._angle_cb.currentTextChanged.connect(self._emit)
        self._force_cb.currentTextChanged.connect(self._emit)

    # ------------------------------------------------------------------
    def _emit(self, _text: str = "") -> None:
        self.units_changed.emit(self.current_units())

    def current_units(self) -> UnitSystem:
        """Return a :class:`UnitSystem` reflecting the current selections."""
        return UnitSystem(
            length=self._length_cb.currentText(),
            angle=self._angle_cb.currentText(),
            force=self._force_cb.currentText(),
        )
