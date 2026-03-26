"""SimulationPanel — left-sidebar result browser for a PMD Session."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)


class SimulationPanel(QWidget):
    """Tree-based browser of Bodies, Joints, and Forces for a Session.

    Emits ``selection_changed`` whenever a checkbox is toggled.
    Each item in the payload list is a dict with keys:
        kind   : "body" | "joint" | "force"
        index  : 0-based index into the model list
        label  : display label
        object : reference to the actual model object
    """

    selection_changed = Signal(object)

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self._session = session
        self._leaves: list[QTreeWidgetItem] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._tree = QTreeWidget()
        self._tree.setHeaderLabel(session.name)
        self._tree.setAnimated(True)
        layout.addWidget(self._tree)

        self._populate()
        self._tree.itemChanged.connect(self._on_item_changed)

    # ------------------------------------------------------------------
    # Tree construction
    # ------------------------------------------------------------------

    def _populate(self):
        self._tree.blockSignals(True)
        model = self._session.model

        root_bodies = self._make_root("Bodies")
        for i, body in enumerate(model.Bodies, start=1):
            label = body.name or f"Body_{i}"
            desc = {"kind": "body", "index": i - 1, "label": label, "object": body}
            self._add_leaf(root_bodies, label, desc)

        root_joints = self._make_root("Joints")
        for i, joint in enumerate(model.Joints, start=1):
            label = joint.name or f"{type(joint).__name__}_{i}"
            desc = {"kind": "joint", "index": i - 1, "label": label, "object": joint}
            self._add_leaf(root_joints, label, desc)

        root_forces = self._make_root("Forces")
        for i, force in enumerate(model.Forces, start=1):
            label = force.name or f"{type(force).__name__}_{i}"
            desc = {"kind": "force", "index": i - 1, "label": label, "object": force}
            self._add_leaf(root_forces, label, desc)

        self._tree.expandAll()
        self._tree.blockSignals(False)

    def _make_root(self, text):
        item = QTreeWidgetItem(self._tree, [text])
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        font = item.font(0)
        font.setBold(True)
        item.setFont(0, font)
        return item

    def _add_leaf(self, parent, label, descriptor):
        item = QTreeWidgetItem(parent, [label])
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(0, Qt.CheckState.Unchecked)
        item.setData(0, Qt.ItemDataRole.UserRole, descriptor)
        self._leaves.append(item)

    # ------------------------------------------------------------------
    # Signal handler
    # ------------------------------------------------------------------

    def _on_item_changed(self, item, column):
        if column != 0:
            return
        if not (item.flags() & Qt.ItemFlag.ItemIsUserCheckable):
            return
        self.selection_changed.emit(self.current_selection())

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def current_selection(self):
        """Return list of descriptor dicts for all currently checked items."""
        return [
            leaf.data(0, Qt.ItemDataRole.UserRole)
            for leaf in self._leaves
            if leaf.checkState(0) == Qt.CheckState.Checked
        ]
