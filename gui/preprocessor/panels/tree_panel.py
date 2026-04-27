"""TreePanel — hierarchical browser of the model spec."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem, QWidget, QVBoxLayout

from ..models import ModelSpec


class TreePanel(QWidget):
    """Read-only tree view with categories: Bodies / Markers / Joints / Forces.

    Signals
    -------
    item_selected(kind, id) : str, str
        Emitted when the user clicks an item. ``kind`` is one of
        ``"body" | "marker" | "joint" | "force"``; ``id`` is the spec's
        unique id.
    """

    item_selected = Signal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(self.tree)

        self._spec: ModelSpec | None = None

    # ──────────────────────────────────────────────────────────
    # Refresh
    # ──────────────────────────────────────────────────────────

    def set_spec(self, spec: ModelSpec) -> None:
        self._spec = spec
        self.refresh()

    def refresh(self) -> None:
        if self._spec is None:
            return
        self.tree.clear()

        cats = [
            ("Bodies",  "body",   self._spec.bodies),
            ("Markers", "marker", self._spec.markers),
            ("Joints",  "joint",  self._spec.joints),
            ("Forces",  "force",  self._spec.forces),
        ]
        for label, kind, items in cats:
            root = QTreeWidgetItem(self.tree, [f"{label} ({len(items)})"])
            root.setExpanded(True)
            for it in items:
                node = QTreeWidgetItem(root, [it.name or it.id])
                node.setData(0, Qt.UserRole, (kind, it.id))

    # ──────────────────────────────────────────────────────────
    # Slots
    # ──────────────────────────────────────────────────────────

    def _on_item_clicked(self, item: QTreeWidgetItem, _column: int):
        data = item.data(0, Qt.UserRole)
        if data is not None:
            kind, _id = data
            self.item_selected.emit(kind, _id)
