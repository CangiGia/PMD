"""NavigationPanel — vertical sidebar navigation for PMD PostProcessor."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QLabel,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..models import reaction_labels

# Body sub-categories: display label → (category_key, [component_keys])
_BODY_CATEGORIES: dict[str, tuple[str, list[str]]] = {
    "Positions":     ("positions",     ["x", "y", "phi"]),
    "Velocities":    ("velocities",    ["dx", "dy", "dphi"]),
    "Accelerations": ("accelerations", ["ddx", "ddy", "ddphi"]),
}

_COMP_LABEL: dict[str, str] = {
    "x": "x",   "y": "y",   "phi":   "\u03c6",
    "dx": "dx", "dy": "dy", "dphi":  "d\u03c6",
    "ddx": "ddx", "ddy": "ddy", "ddphi": "dd\u03c6",
}

_TYPE = Qt.ItemDataRole.UserRole
_DATA = Qt.ItemDataRole.UserRole + 1


class NavigationPanel(QWidget):
    """Vertical navigation sidebar: File / Simulation sections + Settings footer.

    Signals
    -------
    action_triggered(str)
        ``"close"`` | ``"export_plot"`` | ``"export_csv"``
    dark_theme_toggled(bool)
    anim_pane_toggled(bool)
    """

    action_triggered   = Signal(str)
    dark_theme_toggled = Signal(bool)
    anim_pane_toggled  = Signal(bool)

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        if not isinstance(sessions, list):
            sessions = [sessions]
        self._sessions = sessions
        self._leaves: list[QTreeWidgetItem] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # ── Navigation tree (fills available space) ──────────────────
        self._tree = QTreeWidget()
        self._tree.header().hide()
        self._tree.setIndentation(14)
        self._tree.setAnimated(True)
        layout.addWidget(self._tree, stretch=1)

        self._populate()
        self._tree.itemClicked.connect(self._on_item_clicked)

        # ── Settings footer (pinned to bottom) ───────────────────────
        footer = QWidget()
        footer.setObjectName("settings_footer")
        fl = QVBoxLayout(footer)
        fl.setContentsMargins(8, 6, 8, 10)
        fl.setSpacing(6)

        settings_hdr = QLabel("Settings")
        settings_hdr.setObjectName("panel_header")
        fl.addWidget(settings_hdr)

        self._dark_check = QCheckBox("Dark Theme")
        self._dark_check.toggled.connect(self.dark_theme_toggled)
        fl.addWidget(self._dark_check)

        self._anim_check = QCheckBox("Animation Pane")
        self._anim_check.toggled.connect(self.anim_pane_toggled)
        fl.addWidget(self._anim_check)

        layout.addWidget(footer)

    # ------------------------------------------------------------------
    # Tree construction
    # ------------------------------------------------------------------

    def _populate(self):
        self._tree.blockSignals(True)
        self._leaves.clear()

        # ── File section ─────────────────────────────────────────────
        file_root = self._add_section("File")
        self._add_action(file_root, "Close",                     "close")
        self._add_action(file_root, "Save Plot as Image\u2026",  "export_plot")
        self._add_action(file_root, "Export Curves to CSV\u2026","export_csv")
        file_root.setExpanded(True)

        # ── Simulation sections ──────────────────────────────────────
        for session in self._sessions:
            model = session.model
            sim_root = self._add_section(session.name or "Simulation")

            bodies_root = self._add_node(sim_root, "Bodies")
            for i, body in enumerate(model.Bodies, start=1):
                b_label = body.name or f"Body_{i}"
                b_desc = {
                    "kind": "body", "index": i - 1,
                    "label": b_label, "object": body, "session": session,
                }
                b_node = self._add_node(bodies_root, b_label)
                for cat_label, (cat_key, comps) in _BODY_CATEGORIES.items():
                    cat_node = self._add_node(b_node, cat_label)
                    for comp in comps:
                        self._add_leaf(
                            cat_node,
                            _COMP_LABEL.get(comp, comp),
                            {**b_desc, "category": cat_key, "component": comp},
                        )

            joints_root = self._add_node(sim_root, "Joints")
            for i, joint in enumerate(model.Joints, start=1):
                j_label = joint.name or f"{type(joint).__name__}_{i}"
                j_desc = {
                    "kind": "joint", "index": i - 1,
                    "label": j_label, "object": joint, "session": session,
                }
                j_node = self._add_node(joints_root, j_label)
                r_labels = reaction_labels(joint)
                if r_labels:
                    react_node = self._add_node(j_node, "Reactions")
                    for idx, r_lbl in enumerate(r_labels):
                        self._add_leaf(
                            react_node, r_lbl,
                            {**j_desc, "category": "reactions",
                             "component": str(idx)},
                        )

            sim_root.setExpanded(True)
            bodies_root.setExpanded(True)
            joints_root.setExpanded(True)

        self._tree.blockSignals(False)

    # -- helpers --------------------------------------------------------

    def _add_section(self, text: str) -> QTreeWidgetItem:
        item = QTreeWidgetItem(self._tree, [text])
        item.setData(0, _TYPE, "section")
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        f = item.font(0)
        f.setBold(True)
        item.setFont(0, f)
        return item

    def _add_node(self, parent: QTreeWidgetItem, text: str) -> QTreeWidgetItem:
        item = QTreeWidgetItem(parent, [text])
        item.setData(0, _TYPE, "node")
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        return item

    def _add_action(self, parent: QTreeWidgetItem, text: str, action_id: str):
        item = QTreeWidgetItem(parent, [text])
        item.setData(0, _TYPE, "action")
        item.setData(0, _DATA, action_id)
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)

    def _add_leaf(self, parent: QTreeWidgetItem, text: str, data: dict):
        item = QTreeWidgetItem(parent, [text])
        item.setData(0, _TYPE, "leaf")
        item.setData(0, _DATA, data)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(0, Qt.CheckState.Unchecked)
        self._leaves.append(item)

    # ------------------------------------------------------------------
    # Signal handler
    # ------------------------------------------------------------------

    def _on_item_clicked(self, item: QTreeWidgetItem, _column: int):
        if item.data(0, _TYPE) == "action":
            self.action_triggered.emit(item.data(0, _DATA))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_checked_items(self) -> list[dict]:
        """Return data dicts (with 'category' and 'component') for all checked
        leaf nodes."""
        return [
            leaf.data(0, _DATA)
            for leaf in self._leaves
            if leaf.checkState(0) == Qt.CheckState.Checked
        ]

    def clear_checks(self):
        """Uncheck all leaf items without emitting signals."""
        self._tree.blockSignals(True)
        for leaf in self._leaves:
            leaf.setCheckState(0, Qt.CheckState.Unchecked)
        self._tree.blockSignals(False)


# Backward-compatible alias
SimulationPanel = NavigationPanel
