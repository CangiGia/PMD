"""SimulationPanel — left-sidebar result browser for PMD Sessions."""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QLabel,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .. import icons as _icons

_ROOT_ICONS = {
    "Bodies": "mdi6.cube-outline",
    "Joints": "mdi6.link-variant",
}
_LEAF_ICONS = {
    "body":  "mdi6.cube-scan",
    "joint": "mdi6.vector-link",
    "force": "mdi6.lightning-bolt",
}


class SimulationPanel(QWidget):
    """Tree-based browser of Bodies, Joints, and Forces across Sessions.

    Supports one or more sessions.  When a single session is loaded the
    tree reads ``Session_name / Bodies / …``.  With multiple sessions each
    top-level node is a session.

    Emits ``selection_changed`` whenever a checkbox is toggled.
    Each item in the payload list is a dict with keys:
        kind    : "body" | "joint" | "force"
        index   : 0-based index into the model list
        label   : display label
        object  : reference to the actual model object
        session : the Session this item belongs to
    """

    selection_changed = Signal(object)

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        if not isinstance(sessions, list):
            sessions = [sessions]
        self._sessions = sessions
        self._leaves: list[QTreeWidgetItem] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QLabel("Simulations")
        header.setObjectName("panel_header")
        layout.addWidget(header)

        self._tree = QTreeWidget()
        self._tree.setIndentation(16)
        self._tree.header().hide()
        self._tree.setAnimated(True)
        layout.addWidget(self._tree)

        self._populate()
        self._tree.itemChanged.connect(self._on_item_changed)

    # ------------------------------------------------------------------
    # Tree construction
    # ------------------------------------------------------------------

    def _populate(self):
        self._tree.blockSignals(True)

        for session in self._sessions:
            model = session.model
            session_root = self._make_root(session.name)

            root_bodies = self._make_root("Bodies", parent=session_root)
            for i, body in enumerate(model.Bodies, start=1):
                label = body.name or f"Body_{i}"
                desc = {"kind": "body", "index": i - 1, "label": label,
                        "object": body, "session": session}
                self._add_leaf(root_bodies, label, desc)

            root_joints = self._make_root("Joints", parent=session_root)
            for i, joint in enumerate(model.Joints, start=1):
                label = joint.name or f"{type(joint).__name__}_{i}"
                desc = {"kind": "joint", "index": i - 1, "label": label,
                        "object": joint, "session": session}
                self._add_leaf(root_joints, label, desc)

        self._tree.expandAll()
        self._tree.blockSignals(False)

    def _make_root(self, text, parent=None):
        target = parent if parent is not None else self._tree
        item = QTreeWidgetItem(target, [text])
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        font = item.font(0)
        font.setBold(True)
        item.setFont(0, font)
        icon_name = _ROOT_ICONS.get(text, "mdi6.layers-outline")
        item.setIcon(0, _icons.icon(icon_name))
        return item

    def _add_leaf(self, parent, label, descriptor):
        item = QTreeWidgetItem(parent, [label])
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(0, Qt.CheckState.Unchecked)
        item.setData(0, Qt.ItemDataRole.UserRole, descriptor)
        kind = descriptor.get("kind", "")
        item.setIcon(0, _icons.icon(_LEAF_ICONS.get(kind, "mdi6.circle-outline")))
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

    def clear_selection(self):
        """Uncheck all leaf items without emitting selection_changed."""
        self._tree.blockSignals(True)
        for leaf in self._leaves:
            leaf.setCheckState(0, Qt.CheckState.Unchecked)
        self._tree.blockSignals(False)
        self.selection_changed.emit([])
