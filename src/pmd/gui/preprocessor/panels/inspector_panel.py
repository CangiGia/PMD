"""InspectorPanel — typed property editor for the selected spec."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFrame,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ..models import ModelSpec
from .editors import BodyEditor, EditorBase, ForceEditor, JointEditor, MarkerEditor


class InspectorPanel(QWidget):
    """Swap-in typed editor depending on the selected spec kind.

    Signals
    -------
    spec_changed(kind, id) : str, str
        Forwarded from the active editor; the main window listens to
        this and refreshes scene + tree.
    """

    spec_changed = Signal(str, str)
    body_renamed = Signal(str, str, str)   # forwarded from BodyEditor

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 8, 8, 8)
        outer.setSpacing(6)

        # Header
        self._title = QLabel("<i>No selection</i>")
        self._title.setTextFormat(Qt.RichText)
        outer.addWidget(self._title)

        # Separator
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        outer.addWidget(line)

        # Scrollable host for the active editor
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.NoFrame)
        outer.addWidget(self._scroll, 1)

        self._spec: ModelSpec | None = None
        self._editor: EditorBase | None = None

    # ──────────────────────────────────────────────────────────
    def set_spec(self, spec: ModelSpec) -> None:
        self._spec = spec
        self.clear()

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

        editor = self._make_editor(kind, item)
        self._set_editor(editor)

    def clear(self) -> None:
        self._title.setText("<i>No selection</i>")
        self._set_editor(None)

    def show_draft_body(self, spec) -> None:
        """Edit a not-yet-committed body during canvas placement.

        The position/orientation fields are disabled (they get set by
        the next mouse clicks) but the user can freely tune mass,
        inertia, and the shape's principal dimensions before anchoring.
        """
        self._title.setText(
            "<b>Body</b> &nbsp; "
            "<span style='color:#a7350c;'>(draft — pick A, then orientation)</span>"
        )
        editor = BodyEditor(spec)
        editor.set_position_field_enabled(False)
        self._set_editor(editor)

    # ──────────────────────────────────────────────────────────
    def _make_editor(self, kind: str, item) -> EditorBase | None:
        assert self._spec is not None
        if kind == "body":
            return BodyEditor(item)
        if kind == "marker":
            return MarkerEditor(item, self._spec)
        if kind == "joint":
            return JointEditor(item, self._spec)
        if kind == "force":
            return ForceEditor(item, self._spec)
        return None

    def _set_editor(self, editor: EditorBase | None) -> None:
        # Detach old editor
        if self._editor is not None:
            try:
                self._editor.spec_changed.disconnect(self.spec_changed)
            except (RuntimeError, TypeError):
                pass
            try:
                if isinstance(self._editor, BodyEditor):
                    self._editor.body_renamed.disconnect(self.body_renamed)
            except (RuntimeError, TypeError):
                pass
        self._editor = editor

        if editor is None:
            self._scroll.setWidget(QWidget())
            return

        editor.spec_changed.connect(self.spec_changed)
        if isinstance(editor, BodyEditor):
            editor.body_renamed.connect(self.body_renamed)
        # Scroll area takes ownership of the widget.
        self._scroll.setWidget(editor)
