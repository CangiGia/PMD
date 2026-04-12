"""FilterPanel — category / component dropdown bar for result selection."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QWidget,
)

from ..models import reaction_labels
from .. import icons as _icons

# Map from display name → (result‑container key, list of component keys)
_BODY_CATEGORIES = {
    "Positions": ("positions", ["x", "y", "phi"]),
    "Velocities": ("velocities", ["dx", "dy", "dphi"]),
    "Accelerations": ("accelerations", ["ddx", "ddy", "ddphi"]),
}

_JOINT_CATEGORIES = {
    "Reactions": ("reactions", None),  # components built dynamically
}


class FilterPanel(QWidget):
    """Horizontal bar with category / component combos and an *Add curves* button.

    Signals
    -------
    add_curves_requested(category_key, component_key, selection)
        Emitted when the user clicks *Add curves*.
        *category_key*: ``"positions"`` | ``"velocities"`` | ``"accelerations"`` | ``"reactions"``
        *component_key*: ``"x"`` … ``"ddphi"`` for bodies, ``"0"`` … for reactions
        *selection*: list of descriptor dicts from SimulationPanel
    """

    add_curves_requested = Signal(str, str, list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._selection: list[dict] = []

        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 2, 4, 2)

        layout.addWidget(QLabel("Category:"))
        self._category_combo = QComboBox()
        self._category_combo.setMinimumWidth(130)
        layout.addWidget(self._category_combo)

        layout.addWidget(QLabel("Component:"))
        self._component_combo = QComboBox()
        self._component_combo.setMinimumWidth(100)
        layout.addWidget(self._component_combo)

        self._add_btn = QPushButton("Add curves")
        self._add_btn.setIcon(_icons.icon("mdi6.plus"))
        layout.addWidget(self._add_btn)

        layout.addStretch()

        # Internal connections
        self._category_combo.currentTextChanged.connect(self._on_category_changed)
        self._add_btn.clicked.connect(self._on_add_clicked)

        # Start disabled
        self._set_enabled(False)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def update_from_selection(self, selection: list[dict]):
        """Rebuild combos based on the current SimulationPanel selection."""
        self._selection = selection

        kinds = {d["kind"] for d in selection}
        has_body = "body" in kinds
        has_joint = "joint" in kinds

        # Determine available categories
        categories: list[str] = []
        if has_body:
            categories.extend(_BODY_CATEGORIES.keys())
        if has_joint:
            categories.extend(k for k in _JOINT_CATEGORIES.keys()
                              if k not in categories)

        self._category_combo.blockSignals(True)
        self._category_combo.clear()
        self._category_combo.addItems(categories)
        self._category_combo.blockSignals(False)

        if categories:
            self._set_enabled(True)
            self._on_category_changed(self._category_combo.currentText())
        else:
            self._set_enabled(False)
            self._component_combo.clear()

    # ------------------------------------------------------------------
    # Internal slots
    # ------------------------------------------------------------------

    def _on_category_changed(self, text: str):
        self._component_combo.blockSignals(True)
        self._component_combo.clear()

        if text in _BODY_CATEGORIES:
            _, components = _BODY_CATEGORIES[text]
            self._component_combo.addItems(components)
            self._component_combo.setEnabled(True)

        elif text in _JOINT_CATEGORIES:
            labels = self._get_reaction_labels()
            if labels:
                for i, lbl in enumerate(labels):
                    self._component_combo.addItem(lbl, userData=str(i))
                self._component_combo.setEnabled(True)
            else:
                self._component_combo.setEnabled(False)
        else:
            self._component_combo.setEnabled(False)

        self._component_combo.blockSignals(False)
        self._update_add_btn()

    def _on_add_clicked(self):
        cat_text = self._category_combo.currentText()
        comp_text = self._component_combo.currentText()

        # Resolve category key
        if cat_text in _BODY_CATEGORIES:
            cat_key = _BODY_CATEGORIES[cat_text][0]
            comp_key = comp_text
        elif cat_text in _JOINT_CATEGORIES:
            cat_key = _JOINT_CATEGORIES[cat_text][0]
            comp_key = self._component_combo.currentData() or comp_text
        else:
            return

        self.add_curves_requested.emit(cat_key, comp_key, list(self._selection))

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_reaction_labels(self) -> list[str]:
        """Return reaction labels from the first selected joint (or generic fallback)."""
        for d in self._selection:
            if d["kind"] == "joint":
                lbls = reaction_labels(d["object"])
                if lbls:
                    return lbls
        # fallback: use column count
        max_cols = self._get_max_reaction_cols()
        return [f"\u03bb_{i}" for i in range(max_cols)]

    def _get_max_reaction_cols(self) -> int:
        """Return the max number of reaction columns among selected joints."""
        max_cols = 0
        for d in self._selection:
            if d["kind"] == "joint":
                obj = d["object"]
                rc = obj._result_container
                if rc is not None:
                    cols = rc["reactions"].shape[1]
                    if cols > max_cols:
                        max_cols = cols
        return max_cols

    def _update_add_btn(self):
        enabled = (
            bool(self._selection)
            and self._category_combo.currentText() != ""
            and self._component_combo.isEnabled()
            and self._component_combo.currentText() != ""
        )
        self._add_btn.setEnabled(enabled)

    def _set_enabled(self, enabled: bool):
        self._category_combo.setEnabled(enabled)
        self._component_combo.setEnabled(enabled)
        self._add_btn.setEnabled(enabled)
