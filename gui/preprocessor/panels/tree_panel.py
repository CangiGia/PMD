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
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
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
    """

    item_selected = Signal(str, str)

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

        cats = [
            ("Bodies",     "body",   self._spec.bodies),
            ("Markers",    "marker", self._spec.markers),
            ("Connectors", "joint",  self._spec.joints),
            ("Forces",     "force",  self._spec.forces),
        ]
        for label, kind, items in cats:
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
        for i in range(self._browse_tree.topLevelItemCount()):
            root = self._browse_tree.topLevelItem(i)
            any_visible = False
            for j in range(root.childCount()):
                child = root.child(j)
                visible = (not text) or (text in child.text(0).lower())
                child.setHidden(not visible)
                any_visible = any_visible or visible
            root.setHidden(bool(text) and not any_visible)
