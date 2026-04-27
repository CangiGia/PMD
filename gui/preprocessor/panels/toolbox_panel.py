"""ToolboxPanel — vertical tool selector (Select / Body / Marker / Joint / …)."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QButtonGroup, QPushButton, QVBoxLayout, QWidget


class ToolboxPanel(QWidget):
    """Vertical bar with mutually-exclusive tool toggles.

    Signals
    -------
    tool_changed(name) : str
        Emitted when the active tool changes. ``name`` is one of the
        keys registered in :data:`TOOLS` (or ``"select"`` for the
        default arrow tool).
    """

    tool_changed = Signal(str)

    # (key, label, tooltip) — extend as new tools come online.
    TOOLS = [
        ("select",     "Select",     "Selection tool (Esc)"),
        ("body_rect",  "Rect Body",  "Add a rectangular rigid body"),
        ("body_link",  "Link Body",  "Add a 2-point link (rod)"),
        ("body_poly",  "Polygon",    "Add a polygonal rigid body"),
        ("marker",     "Marker",     "Place a marker on a body"),
        ("joint_rev",  "Rev Joint",  "Revolute joint between two markers"),
        ("joint_tran", "Tran Joint", "Translational joint between two markers"),
        ("force_grav", "Gravity",    "Toggle gravity (Weight)"),
        ("force_ptp",  "Spring",     "Point-to-point spring/damper/actuator"),
    ]

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        self._group = QButtonGroup(self)
        self._group.setExclusive(True)
        self._buttons: dict[str, QPushButton] = {}

        for key, label, tip in self.TOOLS:
            btn = QPushButton(label)
            btn.setCheckable(True)
            btn.setToolTip(tip)
            btn.setMinimumHeight(28)
            self._group.addButton(btn)
            self._buttons[key] = btn
            btn.toggled.connect(
                lambda checked, k=key: checked and self.tool_changed.emit(k)
            )
            layout.addWidget(btn)

        layout.addStretch(1)

        # Default tool
        self._buttons["select"].setChecked(True)

    def set_tool(self, key: str) -> None:
        """Programmatically activate a tool (no-op if unknown)."""
        btn = self._buttons.get(key)
        if btn is not None:
            btn.setChecked(True)
