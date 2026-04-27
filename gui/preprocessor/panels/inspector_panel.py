"""InspectorPanel — property editor for the currently selected item.

Stub implementation: shows id, kind and a JSON-like dump of the spec.
A future revision will provide proper typed editors per spec class.
"""

from __future__ import annotations

from dataclasses import asdict

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout, QWidget

from ..models import ModelSpec


class InspectorPanel(QWidget):
    """Read-only property dump for the selected spec (placeholder)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        self._title = QLabel("<i>No selection</i>")
        self._title.setTextFormat(Qt.RichText)
        layout.addWidget(self._title)

        self._dump = QPlainTextEdit()
        self._dump.setReadOnly(True)
        layout.addWidget(self._dump, 1)

        self._spec: ModelSpec | None = None

    def set_spec(self, spec: ModelSpec) -> None:
        self._spec = spec

    def show_item(self, kind: str, item_id: str) -> None:
        if self._spec is None:
            return
        getter = {
            "body":   self._spec.body,
            "marker": self._spec.marker,
            "joint":  self._spec.joint,
            "force":  self._spec.force,
        }.get(kind)
        if getter is None:
            return
        item = getter(item_id)
        if item is None:
            self.clear()
            return
        self._title.setText(
            f"<b>{kind.capitalize()}</b> &nbsp; "
            f"<span style='color:#6b7280'>{item.id}</span>"
        )
        self._dump.setPlainText(_pretty(asdict(item)))

    def clear(self) -> None:
        self._title.setText("<i>No selection</i>")
        self._dump.clear()


def _pretty(d: dict, indent: int = 0) -> str:
    pad = "  " * indent
    lines = []
    for k, v in d.items():
        if isinstance(v, dict):
            lines.append(f"{pad}{k}:")
            lines.append(_pretty(v, indent + 1))
        else:
            lines.append(f"{pad}{k}: {v!r}")
    return "\n".join(lines)
