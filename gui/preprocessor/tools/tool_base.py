"""Tool — abstract base for interactive editing tools.

A *tool* is a stateful object that hijacks scene mouse events while it
is active. The default selection behaviour is implemented by the
:class:`SelectTool` (which simply forwards to ``QGraphicsScene``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from PySide6.QtCore import QObject, Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QGraphicsSceneMouseEvent

if TYPE_CHECKING:
    from ..main_window import PreProcessorWindow
    from ..widgets.canvas_scene import CanvasScene


class Tool(QObject):
    """Base class for canvas tools.

    Subclasses override the ``mouse_*`` and ``key_press`` hooks. Each
    hook should return ``True`` to indicate the event was consumed
    (suppressing the default scene behaviour).

    Parameters
    ----------
    window : PreProcessorWindow
        Owning main window — exposes ``spec``, ``_scene``, ``_view``,
        and refresh helpers.
    """

    name: str = "tool"
    cursor: Qt.CursorShape = Qt.CrossCursor

    def __init__(self, window: "PreProcessorWindow"):
        super().__init__(window)
        self.window = window
        self.spec   = window.spec
        self.scene: "CanvasScene" = window._scene
        self.view  = window._view

    # ── lifecycle ─────────────────────────────────────────────
    def activate(self) -> None:
        self.view.viewport().setCursor(self.cursor)

    def deactivate(self) -> None:
        self.view.viewport().unsetCursor()

    # ── event hooks (override in subclasses) ──────────────────
    def mouse_press(self, event: QGraphicsSceneMouseEvent) -> bool:
        return False

    def mouse_move(self, event: QGraphicsSceneMouseEvent) -> bool:
        return False

    def mouse_release(self, event: QGraphicsSceneMouseEvent) -> bool:
        return False

    def key_press(self, event: QKeyEvent) -> bool:
        # Esc: cancel any in-progress operation and revert to Select.
        if event.key() == Qt.Key_Escape:
            self.window.set_active_tool("select")
            return True
        return False

    # ── helpers ───────────────────────────────────────────────
    def _commit(self) -> None:
        """Refresh tree/inspector/status after spec mutation."""
        self.window._tree.refresh()
        self.window._update_count_label()
