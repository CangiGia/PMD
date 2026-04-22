"""SimulationPanel — left-sidebar result browser for PMD Sessions."""

from __future__ import annotations

from PySide6.QtCore import Qt, QPoint, QSize, Signal
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QMenu,
    QMessageBox,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from .. import icons as _icons


class SimulationPanel(QWidget):
    """Sidebar with File / View / Simulations nav buttons and a footer.

    Emits ``selection_changed`` whenever a leaf checkbox is toggled.
    Each item in the payload is a dict with keys:
        kind    : "body" | "joint"
        index   : 0-based index into the model list
        label   : display label
        object  : reference to the actual model object
        session : the Session this item belongs to
    """

    selection_changed          = Signal(object)
    theme_toggle_requested     = Signal(bool)
    export_requested           = Signal(str)   # "plot" | "csv" | "txt"
    close_requested            = Signal()
    animation_toggle_requested = Signal(bool)

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        if not isinstance(sessions, list):
            sessions = [sessions]
        self._sessions = sessions
        self._leaf_actions: list[tuple[QAction, dict]] = []
        self._icon_items: list[tuple[QPushButton, str, bool]] = []
        self._is_dark = False
        self._is_anim = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # ── Top nav (File, View, Simulations) ───────────────────────────
        top = QWidget()
        top_layout = QVBoxLayout(top)
        top_layout.setContentsMargins(0, 0, 0, 0)
        top_layout.setSpacing(0)

        self._btn_file        = self._make_nav_btn("File",        "mdi6.home-outline")
        self._btn_view        = self._make_nav_btn("View",        "mdi6.binoculars")
        self._btn_simulations = self._make_nav_btn("Simulations", "mdi6.sine-wave")
        for btn in (self._btn_file, self._btn_view, self._btn_simulations):
            top_layout.addWidget(btn)
        top_layout.addStretch(1)

        # ── Footer (Settings, Information, Help) ────────────────────────
        footer = QWidget()
        footer_layout = QVBoxLayout(footer)
        footer_layout.setContentsMargins(0, 0, 0, 0)
        footer_layout.setSpacing(0)

        self._btn_settings    = self._make_nav_btn("Settings",    "mdi6.cog-outline")
        self._btn_information = self._make_nav_btn("Information", "mdi6.information-outline")
        self._btn_help        = self._make_nav_btn("Help",        "mdi6.help-circle-outline")
        for btn in (self._btn_settings, self._btn_information, self._btn_help):
            footer_layout.addWidget(btn)

        # ── Splitter: top (resizable) | footer (fixed) ──────────────────
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.setObjectName("nav_splitter")
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(5)
        splitter.addWidget(top)
        splitter.addWidget(footer)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        layout.addWidget(splitter, stretch=1)

        self._btn_file.clicked.connect(self._on_file)
        self._btn_view.clicked.connect(self._on_view)
        self._btn_simulations.clicked.connect(self._on_simulations)
        self._btn_settings.clicked.connect(self._on_settings)
        self._btn_information.clicked.connect(self._on_information)
        self._btn_help.clicked.connect(self._on_help)

        self._build_leaf_actions()

    # ------------------------------------------------------------------
    # Button factory
    # ------------------------------------------------------------------

    def _make_nav_btn(self, text: str, icon_name: str | None = None) -> QPushButton:
        btn = QPushButton(f"  {text}")
        btn.setProperty("nav", "true")
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        if icon_name:
            btn.setIcon(_icons.icon(icon_name, dim=True))
            btn.setIconSize(QSize(32, 32))
            self._icon_items.append((btn, icon_name, True))
        return btn

    # ------------------------------------------------------------------
    # Leaf QActions (checkable, persistent across menu opens)
    # ------------------------------------------------------------------

    def _build_leaf_actions(self):
        self._leaf_actions.clear()
        for session in self._sessions:
            model = session.model
            for i, body in enumerate(model.Bodies, start=1):
                label = body.name or f"Body_{i}"
                desc = {"kind": "body", "index": i - 1, "label": label,
                        "object": body, "session": session}
                action = QAction(label)
                action.setCheckable(True)
                action.toggled.connect(self._on_leaf_toggled)
                self._leaf_actions.append((action, desc))

            for i, joint in enumerate(model.Joints, start=1):
                label = joint.name or f"{type(joint).__name__}_{i}"
                desc = {"kind": "joint", "index": i - 1, "label": label,
                        "object": joint, "session": session}
                action = QAction(label)
                action.setCheckable(True)
                action.toggled.connect(self._on_leaf_toggled)
                self._leaf_actions.append((action, desc))

    # ------------------------------------------------------------------
    # Theme / state
    # ------------------------------------------------------------------

    def set_dark(self, enabled: bool) -> None:
        self._is_dark = enabled

    def set_animation(self, enabled: bool) -> None:
        self._is_anim = enabled

    # ------------------------------------------------------------------
    # Popup position helper
    # ------------------------------------------------------------------

    def _popup_below(self, btn: QPushButton) -> QPoint:
        """Global position at the bottom-left corner of *btn*."""
        return btn.mapToGlobal(QPoint(0, btn.height()))

    # ------------------------------------------------------------------
    # Top nav handlers
    # ------------------------------------------------------------------

    def _on_file(self):
        menu = QMenu(self)
        menu.addAction("Export Plot...").triggered.connect(
            lambda: self.export_requested.emit("plot")
        )
        menu.addAction("Export CSV...").triggered.connect(
            lambda: self.export_requested.emit("csv")
        )
        menu.addAction("Export TXT...").triggered.connect(
            lambda: self.export_requested.emit("txt")
        )
        menu.addSeparator()
        menu.addAction("Close").triggered.connect(self.close_requested.emit)
        menu.exec(self._popup_below(self._btn_file))

    def _on_view(self):
        menu = QMenu(self)
        anim_action = menu.addAction("Animation Pane")
        anim_action.setCheckable(True)
        anim_action.setChecked(self._is_anim)
        anim_action.toggled.connect(self._on_anim_toggled)
        menu.exec(self._popup_below(self._btn_view))

    def _on_anim_toggled(self, checked: bool):
        self._is_anim = checked
        self.animation_toggle_requested.emit(checked)

    def _on_simulations(self):
        menu = QMenu(self)
        for session in self._sessions:
            session_menu = menu.addMenu(session.name)
            bodies_menu = session_menu.addMenu("Bodies")
            for action, desc in self._leaf_actions:
                if desc["session"] is session and desc["kind"] == "body":
                    bodies_menu.addAction(action)
            joints_menu = session_menu.addMenu("Joints")
            for action, desc in self._leaf_actions:
                if desc["session"] is session and desc["kind"] == "joint":
                    joints_menu.addAction(action)
        menu.exec(self._popup_below(self._btn_simulations))

    # ------------------------------------------------------------------
    # Footer handlers
    # ------------------------------------------------------------------

    def _on_settings(self):
        menu = QMenu(self)
        dark_action = menu.addAction("Dark Theme")
        dark_action.setCheckable(True)
        dark_action.setChecked(self._is_dark)
        dark_action.toggled.connect(self.theme_toggle_requested.emit)
        pos = self._btn_settings.mapToGlobal(
            QPoint(self._btn_settings.width(), 0)
        )
        menu.exec(pos)

    def _on_information(self):
        QMessageBox.information(self, "Information", "Under development")

    def _on_help(self):
        QMessageBox.information(self, "Help", "Under development")

    # ------------------------------------------------------------------
    # Leaf toggle
    # ------------------------------------------------------------------

    def _on_leaf_toggled(self):
        self.selection_changed.emit(self.current_selection())

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def current_selection(self):
        """Return list of descriptor dicts for all currently checked items."""
        return [desc for action, desc in self._leaf_actions if action.isChecked()]

    def clear_selection(self):
        """Uncheck all leaf items without emitting selection_changed."""
        for action, _ in self._leaf_actions:
            action.blockSignals(True)
            action.setChecked(False)
            action.blockSignals(False)
        self.selection_changed.emit([])

    def refresh_icons(self) -> None:
        """Re-apply icons from the current theme (call after set_dark)."""
        for btn, icon_name, dim in self._icon_items:
            btn.setIcon(_icons.icon(icon_name, dim=dim))
