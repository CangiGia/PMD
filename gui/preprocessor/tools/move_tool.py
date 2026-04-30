"""MoveTool \u2014 drag bodies on the canvas; orientation only via Inspector."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QGraphicsItem

from ..widgets import BodyItem
from .tool_base import Tool


class MoveTool(Tool):
    """While active, body items become draggable; everything else stays put.

    Translation is the only DOF exposed: orientation must be edited from
    the Inspector (numeric field). On deactivation, all bodies are
    locked again so accidental drags cannot move them.
    """

    name = "move"
    cursor = Qt.OpenHandCursor

    # ── lifecycle ─────────────────────────────────────────────
    def activate(self) -> None:
        super().activate()
        for item in self.window._body_items.values():
            item.setFlag(QGraphicsItem.ItemIsMovable, True)

    def deactivate(self) -> None:
        super().deactivate()
        for item in self.window._body_items.values():
            item.setFlag(QGraphicsItem.ItemIsMovable, False)

    # ── events ────────────────────────────────────────────────
    def mouse_press(self, event) -> bool:
        # Let the scene handle the press so the body's normal drag kicks in.
        # We don't return True (consumed) so the default behaviour runs.
        return False

    def mouse_release(self, event) -> bool:
        # After a drag the body's itemChange() has already updated the
        # spec. Refresh the inspector / tree so coords reflect the new
        # position.
        body_item = None
        for it in self.scene.selectedItems():
            if isinstance(it, BodyItem):
                body_item = it
                break
        if body_item is not None:
            self.window._inspector.show_item("body", body_item.spec.id)
            self._commit()
        return False
