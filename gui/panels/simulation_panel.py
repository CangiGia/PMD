"""SimulationPanel — left-sidebar result browser for PMD Sessions."""

from PySide6.QtCore import Qt, QPoint, QSize, Signal
from PySide6.QtWidgets import (
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .. import icons as _icons



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
    theme_toggle_requested = Signal(bool)

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        if not isinstance(sessions, list):
            sessions = [sessions]
        self._sessions = sessions
        self._leaves: list[QTreeWidgetItem] = []
        self._icon_items: list[tuple] = []  # (item, icon_name, dim)
        self._is_dark = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QLabel("Simulations")
        header.setObjectName("panel_header")
        layout.addWidget(header)

        # Vertical splitter: tree (resizable) | footer nav buttons
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.setObjectName("nav_splitter")
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(5)

        self._tree = QTreeWidget()
        self._tree.setIndentation(16)
        self._tree.setIconSize(QSize(16, 16))
        self._tree.header().hide()
        self._tree.setAnimated(True)
        splitter.addWidget(self._tree)

        footer = QWidget()
        footer_layout = QVBoxLayout(footer)
        footer_layout.setContentsMargins(0, 0, 0, 0)
        footer_layout.setSpacing(0)
        self._btn_settings    = self._make_nav_btn("Settings",    "mdi6.cog-outline")
        self._btn_information = self._make_nav_btn("Information", "mdi6.information-outline")
        self._btn_help        = self._make_nav_btn("Help",        "mdi6.help-circle-outline")
        for btn in (self._btn_settings, self._btn_information, self._btn_help):
            footer_layout.addWidget(btn)
        splitter.addWidget(footer)

        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        layout.addWidget(splitter, stretch=1)

        self._btn_settings.clicked.connect(self._on_settings)
        self._btn_information.clicked.connect(self._on_information)
        self._btn_help.clicked.connect(self._on_help)

        self._populate()
        self._tree.itemChanged.connect(self._on_item_changed)

    # ------------------------------------------------------------------
    # Nav footer helpers
    # ------------------------------------------------------------------

    def _make_nav_btn(self, text: str, icon_name: str) -> QPushButton:
        btn = QPushButton(f"  {text}")
        btn.setProperty("nav", "true")
        btn.setIcon(_icons.icon(icon_name, dim=True))
        btn.setIconSize(QSize(24, 24))
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._icon_items.append((btn, icon_name, True))
        return btn

    def set_dark(self, enabled: bool) -> None:
        self._is_dark = enabled

    def _on_settings(self):
        menu = QMenu(self)
        dark_action = menu.addAction("Dark Theme")
        dark_action.setCheckable(True)
        dark_action.setChecked(self._is_dark)
        dark_action.triggered.connect(self.theme_toggle_requested.emit)
        pos = self._btn_settings.mapToGlobal(
            QPoint(self._btn_settings.width(), 0)
        )
        menu.exec(pos)

    def _on_information(self):
        QMessageBox.information(self, "Information", "Under development")

    def _on_help(self):
        QMessageBox.information(self, "Help", "Under development")

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

    def clear_selection(self):
        """Uncheck all leaf items without emitting selection_changed."""
        self._tree.blockSignals(True)
        for leaf in self._leaves:
            leaf.setCheckState(0, Qt.CheckState.Unchecked)
        self._tree.blockSignals(False)
        self.selection_changed.emit([])

    def refresh_icons(self) -> None:
        """Re-apply icons from the current theme (call after set_dark)."""
        for item, icon_name, dim in self._icon_items:
            ic = _icons.icon(icon_name, dim=dim)
            if isinstance(item, QTreeWidgetItem):
                item.setIcon(0, ic)
            else:
                item.setIcon(ic)
