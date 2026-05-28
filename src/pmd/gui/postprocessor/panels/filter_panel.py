"""FilterPanel â€” category / component list selectors for result selection."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..models import reaction_labels
from ... import icons as _icons

# Map from display name â†’ (resultâ€‘container key, list of component keys)
_BODY_CATEGORIES = {
    "Positions": ("positions", ["x", "y", "phi"]),
    "Velocities": ("velocities", ["dx", "dy", "dphi"]),
    "Accelerations": ("accelerations", ["ddx", "ddy", "ddphi"]),
}

_JOINT_CATEGORIES = {
    "Reactions": ("reactions", None),  # components built dynamically
}


class FilterPanel(QWidget):
    """Horizontal bar with category / component list widgets and an *Add curves* button.

    Signals
    -------
    add_curves_requested(category_keys, component_keys, selection)
        Emitted when the user clicks *Add curves*.
        *category_keys*: list of selected category keys (e.g. ``["positions"]``)
        *component_keys*: list of selected component keys (e.g. ``["x","y"]``)
        *selection*: list of descriptor dicts from SimulationPanel
    """

    add_curves_requested = Signal(list, list, list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._selection: list[dict] = []
        # Maps display name â†’ (cat_key, comp_key) for each component item
        self._comp_key_map: dict[str, str] = {}
        self._cat_key_map: dict[str, str] = {}

        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 6)
        layout.setSpacing(8)

        # Category list
        cat_col = QVBoxLayout()
        cat_col.setSpacing(2)
        cat_col.addWidget(QLabel("Category:"))
        self._category_list = QListWidget()
        self._category_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection)
        self._category_list.setMaximumHeight(90)
        self._category_list.setMinimumWidth(130)
        cat_col.addWidget(self._category_list)
        layout.addLayout(cat_col)

        # Component list
        comp_col = QVBoxLayout()
        comp_col.setSpacing(2)
        comp_col.addWidget(QLabel("Component:"))
        self._component_list = QListWidget()
        self._component_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection)
        self._component_list.setMaximumHeight(90)
        self._component_list.setMinimumWidth(110)
        comp_col.addWidget(self._component_list)
        layout.addLayout(comp_col)

        self._add_btn = QPushButton("Add curves")
        self._add_btn.setProperty("primary", True)
        self._add_btn.setIcon(_icons.icon("mdi6.plus"))
        layout.addWidget(self._add_btn)

        layout.addStretch()

        # Internal connections
        self._category_list.itemSelectionChanged.connect(
            self._on_category_selection_changed)
        self._component_list.itemSelectionChanged.connect(self._update_add_btn)
        self._add_btn.clicked.connect(self._on_add_clicked)

        # Start disabled
        self._set_enabled(False)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def update_from_selection(self, selection: list[dict]):
        """Rebuild lists based on the current SimulationPanel selection."""
        self._selection = selection

        kinds = {d["kind"] for d in selection}
        has_body = "body" in kinds
        has_joint = "joint" in kinds

        categories: list[str] = []
        if has_body:
            categories.extend(_BODY_CATEGORIES.keys())
        if has_joint:
            categories.extend(k for k in _JOINT_CATEGORIES.keys()
                              if k not in categories)

        self._category_list.blockSignals(True)
        self._category_list.clear()
        self._cat_key_map.clear()
        for name in categories:
            if name in _BODY_CATEGORIES:
                self._cat_key_map[name] = _BODY_CATEGORIES[name][0]
            elif name in _JOINT_CATEGORIES:
                self._cat_key_map[name] = _JOINT_CATEGORIES[name][0]
            self._category_list.addItem(name)
        self._category_list.blockSignals(False)

        if categories:
            self._set_enabled(True)
            # Pre-select first category and rebuild components
            self._category_list.setCurrentRow(0)
            self._on_category_selection_changed()
        else:
            self._set_enabled(False)
            self._component_list.clear()

    # ------------------------------------------------------------------
    # Internal slots
    # ------------------------------------------------------------------

    def _on_category_selection_changed(self):
        """Rebuild the component list when the category selection changes."""
        selected_cats = [
            item.text() for item in self._category_list.selectedItems()
        ]
        # Use the first selected category to populate components.
        first_cat = selected_cats[0] if selected_cats else ""

        self._component_list.blockSignals(True)
        self._component_list.clear()
        self._comp_key_map.clear()

        if first_cat in _BODY_CATEGORIES:
            _, components = _BODY_CATEGORIES[first_cat]
            for comp in components:
                self._component_list.addItem(comp)
                self._comp_key_map[comp] = comp
            self._component_list.setEnabled(True)

        elif first_cat in _JOINT_CATEGORIES:
            labels = self._get_reaction_labels()
            for i, lbl in enumerate(labels):
                self._component_list.addItem(lbl)
                self._comp_key_map[lbl] = str(i)
            self._component_list.setEnabled(bool(labels))

        else:
            self._component_list.setEnabled(False)

        self._component_list.blockSignals(False)
        self._update_add_btn()

    def _on_add_clicked(self):
        selected_cat_names = [
            item.text() for item in self._category_list.selectedItems()
        ]
        selected_comp_names = [
            item.text() for item in self._component_list.selectedItems()
        ]

        cat_keys = [self._cat_key_map.get(n, n) for n in selected_cat_names]
        comp_keys = [self._comp_key_map.get(n, n) for n in selected_comp_names]

        if not cat_keys or not comp_keys:
            return

        self.add_curves_requested.emit(cat_keys, comp_keys, list(self._selection))

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
            and bool(self._category_list.selectedItems())
            and bool(self._component_list.selectedItems())
            and self._component_list.isEnabled()
        )
        self._add_btn.setEnabled(enabled)

    def _set_enabled(self, enabled: bool):
        self._category_list.setEnabled(enabled)
        self._component_list.setEnabled(enabled)
        self._add_btn.setEnabled(enabled)

    def refresh_icons(self) -> None:
        """Re-apply icons from the current theme (call after set_dark)."""
        self._add_btn.setIcon(_icons.icon("mdi6.plus"))
