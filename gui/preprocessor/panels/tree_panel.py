"""TreePanel — Adams-style left-side model browser.

Layout
------

::

    ┌──────────────────────┐
    │  [.MODEL_1   ▼]      │  ← model selector
    ├──────────────────────┤
    │ ╭ Browse | Groups | Filters ╮
    │ │  Bodies (n)              │
    │ │  Connectors (n)          │
    │ │  Forces (n)              │
    │ │  Markers (n)             │
    │ ╰──────────────────────────╯
    ├──────────────────────┤
    │ Search: [          ] │
    └──────────────────────┘
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..models import ModelSpec


class TreePanel(QWidget):
    """Hierarchical browser of the model spec with Browse/Groups/Filters tabs.

    Signals
    -------
    item_selected(kind, id) : str, str
        Kind ∈ {body, marker, joint, force}; id is the spec id.
    item_delete_requested(kind, id) : str, str
        Right-click → Delete or Del key on a tree item.
    """

    item_selected         = Signal(str, str)
    item_delete_requested = Signal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # ── Model selector ────────────────────────────────────
        self._model_combo = QComboBox()
        self._model_combo.addItem(".MODEL_1")
        outer.addWidget(self._model_combo, 0)

        # ── Tabs ──────────────────────────────────────────────
        self._tabs = QTabWidget()
        self._tabs.setDocumentMode(True)
        outer.addWidget(self._tabs, 1)

        self._browse_tree = QTreeWidget()
        self._browse_tree.setHeaderHidden(True)
        self._browse_tree.itemClicked.connect(self._on_item_clicked)
        self._browse_tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self._browse_tree.customContextMenuRequested.connect(
            self._on_context_menu)
        # Del key on the tree triggers a delete request for the
        # currently focused spec node.
        self._del_action = QAction("Delete", self._browse_tree)
        self._del_action.setShortcut(QKeySequence.Delete)
        self._del_action.setShortcutContext(Qt.WidgetWithChildrenShortcut)
        self._del_action.triggered.connect(self._on_delete_shortcut)
        self._browse_tree.addAction(self._del_action)
        self._tabs.addTab(self._browse_tree, "Browse")

        groups_page = QWidget()
        gl = QVBoxLayout(groups_page)
        gl.addWidget(QLabel("<i>No groups defined.</i>"))
        gl.addStretch(1)
        self._tabs.addTab(groups_page, "Groups")

        filt_page = QWidget()
        fl = QVBoxLayout(filt_page)
        fl.addWidget(QLabel("<i>Filters: (not implemented)</i>"))
        fl.addStretch(1)
        self._tabs.addTab(filt_page, "Filters")

        # ── Search ────────────────────────────────────────────
        search_row = QHBoxLayout()
        search_row.setContentsMargins(0, 0, 0, 0)
        search_row.addWidget(QLabel("Search:"))
        self._search = QLineEdit()
        self._search.setPlaceholderText("Filter tree…")
        self._search.textChanged.connect(self._apply_search_filter)
        search_row.addWidget(self._search, 1)
        outer.addLayout(search_row, 0)

        self._spec: ModelSpec | None = None

    # ──────────────────────────────────────────────────────────
    def set_spec(self, spec: ModelSpec) -> None:
        self._spec = spec
        self.refresh()

    def refresh(self) -> None:
        if self._spec is None:
            return
        self._browse_tree.clear()

        # ── Bodies (with their child markers nested below) ──
        bodies_root = QTreeWidgetItem(
            self._browse_tree,
            [f"Bodies ({len(self._spec.bodies)})"]
        )
        bodies_root.setExpanded(True)
        for b in self._spec.bodies:
            node = QTreeWidgetItem(bodies_root, [b.name or b.id])
            node.setData(0, Qt.UserRole, ("body", b.id))
            # Nested markers belonging to this body.
            for m in self._spec.markers:
                if m.body_id == b.id:
                    leaf = QTreeWidgetItem(node, [m.name or m.id])
                    leaf.setData(0, Qt.UserRole, ("marker", m.id))
            node.setExpanded(False)

        # ── Connectors / forces ──────────────────────────────
        for label, kind, items in (
            ("Connectors", "joint", self._spec.joints),
            ("Forces",     "force", self._spec.forces),
        ):
            root = QTreeWidgetItem(self._browse_tree,
                                   [f"{label} ({len(items)})"])
            root.setExpanded(True)
            for it in items:
                node = QTreeWidgetItem(root, [it.name or it.id])
                node.setData(0, Qt.UserRole, (kind, it.id))

        self._apply_search_filter(self._search.text())

    # ──────────────────────────────────────────────────────────
    def _on_item_clicked(self, item: QTreeWidgetItem, _column: int):
        data = item.data(0, Qt.UserRole)
        if data is not None:
            kind, _id = data
            self.item_selected.emit(kind, _id)

    def _apply_search_filter(self, text: str) -> None:
        text = text.strip().lower()

        def filter_item(node: QTreeWidgetItem) -> bool:
            """Return True if ``node`` (or any descendant) matches."""
            self_match = (not text) or (text in node.text(0).lower())
            any_child = False
            for k in range(node.childCount()):
                if filter_item(node.child(k)):
                    any_child = True
            visible = bool(self_match or any_child)
            node.setHidden(not visible)
            return visible

        for i in range(self._browse_tree.topLevelItemCount()):
            filter_item(self._browse_tree.topLevelItem(i))

    # ──────────────────────────────────────────────────────────
    # Context menu / Delete shortcut
    # ──────────────────────────────────────────────────────────

    def _on_context_menu(self, pos):
        item = self._browse_tree.itemAt(pos)
        if item is None:
            return
        data = item.data(0, Qt.UserRole)
        if data is None:
            return
        kind, _id = data
        menu = QMenu(self._browse_tree)
        act_select = menu.addAction("Edit in Inspector")
        act_delete = menu.addAction("Delete")
        chosen = menu.exec(self._browse_tree.viewport().mapToGlobal(pos))
        if chosen is act_delete:
            self.item_delete_requested.emit(kind, _id)
        elif chosen is act_select:
            self.item_selected.emit(kind, _id)

    def _on_delete_shortcut(self):
        item = self._browse_tree.currentItem()
        if item is None:
            return
        data = item.data(0, Qt.UserRole)
        if data is None:
            return
        kind, _id = data
        self.item_delete_requested.emit(kind, _id)
