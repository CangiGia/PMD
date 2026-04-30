"""SelectTool — the default arrow tool (selection + drag-move)."""

from __future__ import annotations

from PySide6.QtCore import Qt

from .tool_base import Tool


class SelectTool(Tool):
    """No-op tool: defers everything to the underlying QGraphicsScene."""

    name = "select"
    cursor = Qt.ArrowCursor

    # All handlers return False → scene processes events normally.
