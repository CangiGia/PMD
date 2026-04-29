"""RibbonBar — Adams-style top ribbon with tabs and grouped tool buttons.

The ribbon publishes two signals:

* ``tool_changed(name)`` — a *modal* tool was activated (mutually
  exclusive with all the other tool buttons; ``"select"`` deactivates
  every tool).
* ``action_triggered(name)`` — a *one-shot* command was invoked
  (e.g. running a simulation, opening the post-processor).
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QSize, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


# ──────────────────────────────────────────────────────────────────
# Local stylesheet — ensures tab + tool-button text are readable
# regardless of what the global theme defines for QWidget.
# ──────────────────────────────────────────────────────────────────
_RIBBON_QSS = """
QTabWidget::pane {
    border: 1px solid #8a9bb4;
    border-top: none;
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #f4f6fa, stop:1 #dfe4ec);
}
QTabBar {
    background: transparent;
    qproperty-drawBase: 0;
}
QTabBar::tab {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #fbfcfe, stop:1 #cdd6e3);
    color: #1c2033;
    padding: 5px 18px 4px 18px;
    margin: 2px 1px 0 1px;
    border: 1px solid #8a9bb4;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    min-height: 18px;
    font-size: 9pt;
}
QTabBar::tab:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #ffffff, stop:1 #dde6f2);
    color: #1a5cbf;
}
QTabBar::tab:selected {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #ffffff, stop:1 #eef2f8);
    color: #1a5cbf;
    font-weight: 600;
    margin-top: 0px;
    padding-top: 7px;
    border-color: #5a6f8f;
}
QTabBar::tab:!selected {
    margin-top: 3px;
}
QToolButton {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #ffffff, stop:0.5 #eef2f8, stop:1 #c9d3e2);
    color: #1c2033;
    border: 1px solid #8a9bb4;
    border-radius: 4px;
    padding: 3px 8px;
    font-size: 9pt;
}
QToolButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #ffffff, stop:0.5 #e1ebf9, stop:1 #b6c8e3);
    border-color: #3f8cff;
    color: #1a5cbf;
}
QToolButton:pressed, QToolButton:checked {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
        stop:0 #b6c8e3, stop:0.5 #cddaee, stop:1 #e6ecf6);
    border: 1px solid #2f6fcf;
    color: #1a5cbf;
}
QWidget#RibbonGroup {
    border-right: 1px solid #b6bfd0;
    background: transparent;
}
"""


# ──────────────────────────────────────────────────────────────────
# Ribbon group (one labelled box inside a tab)
# ──────────────────────────────────────────────────────────────────

class RibbonGroup(QWidget):
    """A vertically-stacked group of buttons + a centred title label.

    Visually mimics MS-Office / Adams ribbon groups: buttons on top,
    group name in small caption below, thin separator on the right.
    """

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setObjectName("RibbonGroup")

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 2, 4, 1)
        outer.setSpacing(1)

        self._buttons_host = QWidget()
        self._row = QHBoxLayout(self._buttons_host)
        self._row.setContentsMargins(0, 0, 0, 0)
        self._row.setSpacing(1)
        outer.addWidget(self._buttons_host, 1)

        self._title = QLabel(title)
        self._title.setAlignment(Qt.AlignCenter)
        f = self._title.font()
        f.setPointSize(max(7, f.pointSize() - 2))
        self._title.setFont(f)
        self._title.setStyleSheet("color: #6b7280;")
        outer.addWidget(self._title, 0)

    # ----------------------------------------------------------
    def add_button(self, label: str, tooltip: str = "",
                   large: bool = True) -> QToolButton:
        """Add a (text-only for now) tool button to this group."""
        btn = QToolButton()
        btn.setText(label)
        btn.setToolTip(tooltip or label)
        btn.setAutoRaise(False)
        btn.setToolButtonStyle(
            Qt.ToolButtonTextUnderIcon if large else Qt.ToolButtonTextOnly
        )
        if large:
            btn.setMinimumSize(QSize(56, 46))
            btn.setMaximumHeight(50)
            btn.setIconSize(QSize(20, 20))
        else:
            btn.setMinimumHeight(22)
        btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self._row.addWidget(btn)
        return btn


# ──────────────────────────────────────────────────────────────────
# Ribbon bar (the QTabWidget at the top of the main window)
# ──────────────────────────────────────────────────────────────────

class RibbonBar(QTabWidget):
    """Top-level ribbon. Holds one ``QWidget`` page per category tab."""

    tool_changed     = Signal(str)
    action_triggered = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDocumentMode(True)
        self.setMovable(False)
        self.setTabPosition(QTabWidget.North)
        self.setMinimumHeight(86)
        self.setMaximumHeight(96)
        self.setStyleSheet(_RIBBON_QSS)

        self._tool_group = QButtonGroup(self)
        self._tool_group.setExclusive(True)
        self._tool_buttons: dict[str, QToolButton] = {}

        # Hidden silent button to allow "no tool active" state
        self._silent = QToolButton()
        self._silent.setCheckable(True)
        self._silent.setVisible(False)
        self._tool_group.addButton(self._silent)

        self._build_all_tabs()

    # ----------------------------------------------------------
    # Public API
    # ----------------------------------------------------------

    def set_tool(self, key: str) -> None:
        """Programmatically activate a tool. ``"select"`` disables all."""
        if key == "select":
            self._silent.setChecked(True)
            return
        btn = self._tool_buttons.get(key)
        if btn is not None:
            btn.setChecked(True)

    # ----------------------------------------------------------
    # Internal builders
    # ----------------------------------------------------------

    def _new_tab(self, title: str) -> QHBoxLayout:
        page = QWidget()
        lay = QHBoxLayout(page)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(2)
        self.addTab(page, title)
        return lay

    def _add_tool(self, group: RibbonGroup, key: str, label: str,
                  tooltip: str = "") -> QToolButton:
        btn = group.add_button(label, tooltip)
        btn.setCheckable(True)
        self._tool_group.addButton(btn)
        self._tool_buttons[key] = btn
        btn.toggled.connect(
            lambda checked, k=key: checked and self.tool_changed.emit(k)
        )
        return btn

    def _add_action(self, group: RibbonGroup, key: str, label: str,
                    tooltip: str = "") -> QToolButton:
        btn = group.add_button(label, tooltip)
        btn.setCheckable(False)
        btn.clicked.connect(lambda _=False, k=key: self.action_triggered.emit(k))
        return btn

    # ----------------------------------------------------------
    def _build_all_tabs(self) -> None:
        self._build_bodies_tab()
        self._build_connectors_tab()
        self._build_motions_tab()
        self._build_forces_tab()
        self._build_elements_tab()
        self._build_simulation_tab()
        self._build_results_tab()

    # ── Bodies ────────────────────────────────────────────────
    def _build_bodies_tab(self):
        lay = self._new_tab("Bodies")
        g = RibbonGroup("Solids"); lay.addWidget(g)
        self._add_tool(g, "body_rect",   "Rect",    "Rectangular rigid body")
        self._add_tool(g, "body_link",   "Link",    "Two-point link / rod")
        self._add_tool(g, "body_circle", "Circle",  "Circular rigid body")
        self._add_tool(g, "body_poly",   "Polygon", "Polygonal rigid body")

        g2 = RibbonGroup("Edit"); lay.addWidget(g2)
        self._add_tool(g2, "move", "Move", "Drag bodies to translate (orientation: Inspector only)")
        lay.addStretch(1)

    # ── Connectors (joints) ───────────────────────────────────
    def _build_connectors_tab(self):
        lay = self._new_tab("Connectors")
        g = RibbonGroup("Joints"); lay.addWidget(g)
        self._add_tool(g, "joint_rev",    "Revolute",      "Revolute joint")
        self._add_tool(g, "joint_tran",   "Translational", "Translational joint")
        self._add_tool(g, "joint_revrev", "Rev-Rev",       "Two-link couple")
        self._add_tool(g, "joint_revtran","Rev-Tran",      "Rev-Tran couple")
        self._add_tool(g, "joint_rigid",  "Rigid",         "Rigid (fixed) joint")
        self._add_tool(g, "joint_disc",   "Disc",          "Disc on plane")
        lay.addStretch(1)

    # ── Motions ───────────────────────────────────────────────
    def _build_motions_tab(self):
        lay = self._new_tab("Motions")
        g = RibbonGroup("Drivers"); lay.addWidget(g)
        self._add_tool(g, "motion_rev",  "Rotational",    "Prescribed rotation")
        self._add_tool(g, "motion_tran", "Translational", "Prescribed translation")
        lay.addStretch(1)

    # ── Forces ────────────────────────────────────────────────
    def _build_forces_tab(self):
        lay = self._new_tab("Forces")
        g1 = RibbonGroup("Global"); lay.addWidget(g1)
        self._add_action(g1, "force_grav", "Gravity", "Toggle global weight (–9.81 m/s² in Y)")

        g2 = RibbonGroup("Applied"); lay.addWidget(g2)
        self._add_tool(g2, "force_ptp",    "PtP Spring",  "Point-to-point spring/damper")
        self._add_tool(g2, "force_rotsda", "Rot SDA",     "Rotational spring/damper/actuator")
        self._add_tool(g2, "force_local",  "Local Force", "Force in body frame")
        self._add_tool(g2, "force_global", "Global Force","Force in world frame")
        self._add_tool(g2, "torque",       "Torque",      "Pure torque on a body")
        lay.addStretch(1)

    # ── Elements (markers, construction) ──────────────────────
    def _build_elements_tab(self):
        lay = self._new_tab("Elements")
        g = RibbonGroup("Construction"); lay.addWidget(g)
        self._add_tool(g, "marker", "Marker", "Place a marker on a body")
        lay.addStretch(1)

    # ── Simulation ────────────────────────────────────────────
    def _build_simulation_tab(self):
        lay = self._new_tab("Simulation")
        g = RibbonGroup("Run"); lay.addWidget(g)
        self._add_action(g, "sim_dynamic",   "Dynamic",   "Run dynamic analysis")
        self._add_action(g, "sim_kinematic", "Kinematic", "Run kinematic analysis")
        self._add_action(g, "sim_static",    "Static",    "Run static equilibrium")
        g2 = RibbonGroup("Settings"); lay.addWidget(g2)
        self._add_action(g2, "sim_settings", "Solver…", "Solver settings")
        lay.addStretch(1)

    # ── Results ───────────────────────────────────────────────
    def _build_results_tab(self):
        lay = self._new_tab("Results")
        g = RibbonGroup("Post-Processor"); lay.addWidget(g)
        self._add_action(g, "post_open", "Open", "Open the post-processor")
        lay.addStretch(1)
