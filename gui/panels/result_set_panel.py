"""ResultSetPanel — list of active curves with visibility toggles."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QPixmap, QIcon
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..models import CurveItem
from .. import icons as _icons


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
        # Each entry stores the original build_curves arguments for a batch
        # so that rebuild_all() can replay them with new display units.
        self._requests: list[dict] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QLabel("Active Curves")
        header.setObjectName("panel_header")
        layout.addWidget(header)

        self._curve_list = QListWidget()
        self._curve_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self._curve_list.itemChanged.connect(self._on_item_changed)
        layout.addWidget(self._curve_list)

        btn_layout = QHBoxLayout()
        btn_layout.setContentsMargins(4, 2, 4, 2)

        self._remove_btn = QPushButton("Remove selected")
        self._remove_btn.setIcon(_icons.icon("mdi6.minus"))
        self._remove_btn.setEnabled(False)
        self._remove_btn.clicked.connect(self.remove_selected)
        btn_layout.addWidget(self._remove_btn)

        self._clear_btn = QPushButton("Clear all")
        self._clear_btn.setIcon(_icons.icon("mdi6.trash-can-outline"))
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
        """Append *curves* to the list without storing a rebuild request.

        Use :meth:`add_curves_with_request` when curves should be rebuildable
        on unit-system changes.
        """
        self._append_curves(curves)
        self.curves_changed.emit()

    def add_curves_with_request(
        self,
        curves: list[CurveItem],
        category: str,
        component: str,
        selection: list[dict],
    ):
        """Append *curves* and store the original request for later rebuilds.

        Parameters
        ----------
        curves : list[CurveItem]
            Curves already built from the request (with current display units).
        category, component, selection :
            The arguments originally passed to :func:`build_curves`.
        """
        self._requests.append({
            "category": category,
            "component": component,
            "selection": selection,
        })
        self._append_curves(curves)
        self.curves_changed.emit()

    def rebuild_all(self, display_units, build_curves_fn) -> None:
        """Rebuild all stored requests with a new display unit system.

        Existing curves are replaced by freshly scaled counterparts.
        Colors and visibility states are preserved by matching curve labels.

        Parameters
        ----------
        display_units : UnitSystem
            The new display unit system.
        build_curves_fn : callable
            The ``build_curves`` factory, signature:
            ``(category, component, selection, display_units) -> list[CurveItem]``.
        """
        if not self._requests:
            return

        # Snapshot current label → (color, visible) for state preservation
        label_state: dict[str, tuple[str, bool]] = {}
        for i, curve in enumerate(self._curves):
            item = self._curve_list.item(i)
            vis = (item.checkState() == Qt.CheckState.Checked) if item else curve.visible
            label_state[curve.label] = (curve.color, vis)

        # Clear display without touching _requests
        self._curve_list.blockSignals(True)
        self._curve_list.clear()
        self._curves.clear()
        self._curve_list.blockSignals(False)

        # Replay each stored request
        for req in self._requests:
            new_curves = build_curves_fn(
                req["category"], req["component"], req["selection"], display_units
            )
            # Keep only curves that were present before (i.e. not user-removed)
            kept = [c for c in new_curves if c.label in label_state]
            for c in kept:
                c.color, c.visible = label_state[c.label]
            self._append_curves(kept)

        self._clear_btn.setEnabled(bool(self._curves))
        self._update_remove_btn()
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
        """Remove all curves and stored requests."""
        self._curve_list.clear()
        self._curves.clear()
        self._requests.clear()
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
    # Internal helpers
    # ------------------------------------------------------------------

    def _append_curves(self, curves: list[CurveItem]) -> None:
        """Add *curves* to the list widget and internal list (no signal)."""
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

    def _on_item_changed(self, item: QListWidgetItem):
        row = self._curve_list.row(item)
        if 0 <= row < len(self._curves):
            checked = item.checkState() == Qt.CheckState.Checked
            self._curves[row].visible = checked
            self.curves_changed.emit()

    def _update_remove_btn(self):
        self._remove_btn.setEnabled(bool(self._curve_list.selectedItems()))
