"""ResultSetPanel — list of active curves with visibility toggles."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QPixmap, QIcon
from PySide6.QtWidgets import (
    QHBoxLayout,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..models import CurveItem


def _color_icon(hex_color: str, size: int = 12) -> QIcon:
    """Create a small square icon filled with *hex_color*."""
    pm = QPixmap(size, size)
    pm.fill(QColor(hex_color))
    return QIcon(pm)


class ResultSetPanel(QWidget):
    """Panel listing all active curves with checkboxes and remove buttons.

    Signals
    -------
    curves_changed()
        Emitted after any add / remove / clear / visibility toggle.
    """

    curves_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._curves: list[CurveItem] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._curve_list = QListWidget()
        self._curve_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self._curve_list.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self._curve_list)

        btn_layout = QHBoxLayout()
        btn_layout.setContentsMargins(4, 2, 4, 2)

        self._remove_btn = QPushButton("Remove selected")
        self._remove_btn.setEnabled(False)
        self._remove_btn.clicked.connect(self.remove_selected)
        btn_layout.addWidget(self._remove_btn)

        self._clear_btn = QPushButton("Clear all")
        self._clear_btn.setEnabled(False)
        self._clear_btn.clicked.connect(self.clear)
        btn_layout.addWidget(self._clear_btn)

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        self._curve_list.itemSelectionChanged.connect(self._update_remove_btn)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add_curves(self, curves: list[CurveItem]):
        """Append *curves* to the list and create matching QListWidgetItems."""
        self._curve_list.blockSignals(True)
        for curve in curves:
            item = QListWidgetItem(_color_icon(curve.color), curve.label)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(
                Qt.CheckState.Checked if curve.visible else Qt.CheckState.Unchecked
            )
            self._curve_list.addItem(item)
            self._curves.append(curve)
        self._curve_list.blockSignals(False)
        self._clear_btn.setEnabled(bool(self._curves))
        self.curves_changed.emit()

    def remove_selected(self):
        """Remove curves whose QListWidgetItems are *selected* (highlighted)."""
        selected_rows = sorted(
            [self._curve_list.row(it) for it in self._curve_list.selectedItems()],
            reverse=True,
        )
        for row in selected_rows:
            self._curve_list.takeItem(row)
            del self._curves[row]
        self._clear_btn.setEnabled(bool(self._curves))
        self._update_remove_btn()
        self.curves_changed.emit()

    def clear(self):
        """Remove all curves."""
        self._curve_list.clear()
        self._curves.clear()
        self._clear_btn.setEnabled(False)
        self._update_remove_btn()
        self.curves_changed.emit()

    def visible_curves(self) -> list[CurveItem]:
        """Return only curves whose checkbox is checked."""
        result = []
        for i, curve in enumerate(self._curves):
            item = self._curve_list.item(i)
            if item and item.checkState() == Qt.CheckState.Checked:
                result.append(curve)
        return result

    # ------------------------------------------------------------------
    # Internal slots
    # ------------------------------------------------------------------

    def _on_item_changed(self, item: QListWidgetItem):
        row = self._curve_list.row(item)
        if 0 <= row < len(self._curves):
            checked = item.checkState() == Qt.CheckState.Checked
            self._curves[row].visible = checked
            self.curves_changed.emit()

    def _update_remove_btn(self):
        self._remove_btn.setEnabled(bool(self._curve_list.selectedItems()))
