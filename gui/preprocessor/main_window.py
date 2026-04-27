"""PreProcessorWindow — main shell of the model-builder GUI.

Layout (left → right):

    ┌──────────┬──────────────┬──────────────────┬─────────────┐
    │ Toolbox  │ Tree browser │  Canvas (centre) │  Inspector  │
    │ (dock L) │ (dock L)     │                  │  (dock R)   │
    └──────────┴──────────────┴──────────────────┴─────────────┘

All side panes are :class:`QDockWidget`-based so the user can
re-arrange / detach / hide them as desired.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QDockWidget,
    QLabel,
    QMainWindow,
    QStatusBar,
)

from .models import BodySpec, ModelSpec, ShapeSpec
from .panels import InspectorPanel, ToolboxPanel, TreePanel
from .widgets import BodyItem, CanvasScene, CanvasView


class PreProcessorWindow(QMainWindow):
    """Top-level QMainWindow for the PMD pre-processor.

    Parameters
    ----------
    spec : ModelSpec
        The project specification edited by this window.
    parent : QWidget or None
        Optional parent widget.
    """

    def __init__(self, spec: ModelSpec, parent=None):
        super().__init__(parent)
        self.spec = spec
        self._active_tool = "select"

        self.setWindowTitle("PMD Pre-Processor — Model Builder")
        self.resize(1500, 850)

        self._build_central_canvas()
        self._build_docks()
        self._build_menu_bar()
        self._build_status_bar()

        # Wire up
        self._toolbox.tool_changed.connect(self._on_tool_changed)
        self._tree.item_selected.connect(self._inspector.show_item)
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)

    # ──────────────────────────────────────────────────────────
    # Build helpers
    # ──────────────────────────────────────────────────────────

    def _build_central_canvas(self):
        self._scene = CanvasScene(grid_step=0.05)
        self._view = CanvasView(self._scene)
        self._view.cursor_moved.connect(self._on_cursor_moved)
        self._view.scene().selectionChanged.connect(self._on_selection_changed)
        self.setCentralWidget(self._view)

    def _build_docks(self):
        # Toolbox (left)
        self._toolbox = ToolboxPanel()
        dock_tools = QDockWidget("Tools", self)
        dock_tools.setWidget(self._toolbox)
        dock_tools.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_tools)

        # Tree browser (left, beside toolbox)
        self._tree = TreePanel()
        dock_tree = QDockWidget("Model Tree", self)
        dock_tree.setWidget(self._tree)
        dock_tree.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_tree)
        self.splitDockWidget(dock_tools, dock_tree, Qt.Horizontal)

        # Inspector (right)
        self._inspector = InspectorPanel()
        dock_inspect = QDockWidget("Inspector", self)
        dock_inspect.setWidget(self._inspector)
        dock_inspect.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_inspect)

        # Default sizes
        self.resizeDocks([dock_tools, dock_tree, dock_inspect],
                         [120, 220, 280], Qt.Horizontal)

    def _build_menu_bar(self):
        mbar = self.menuBar()

        # File
        file_menu = mbar.addMenu("&File")
        act_new = QAction("&New Project", self)
        act_new.triggered.connect(self._on_new_project)
        file_menu.addAction(act_new)
        file_menu.addSeparator()
        act_quit = QAction("&Quit", self)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        # View
        view_menu = mbar.addMenu("&View")
        act_fit = QAction("Zoom to &Fit", self)
        act_fit.setShortcut("F")
        act_fit.triggered.connect(self._view.zoom_to_fit)
        view_menu.addAction(act_fit)

        # Simulate (placeholder)
        sim_menu = mbar.addMenu("&Simulate")
        act_run = QAction("&Run…", self)
        act_run.setEnabled(False)  # not implemented yet
        sim_menu.addAction(act_run)

    def _build_status_bar(self):
        sb = QStatusBar(self)
        self.setStatusBar(sb)

        self._lbl_coords = QLabel("x=0.000  y=0.000 m")
        self._lbl_tool   = QLabel("tool: select")
        self._lbl_count  = QLabel("bodies=0  joints=0  forces=0")

        sb.addWidget(self._lbl_coords, 1)
        sb.addPermanentWidget(self._lbl_tool)
        sb.addPermanentWidget(self._lbl_count)

    # ──────────────────────────────────────────────────────────
    # Slots
    # ──────────────────────────────────────────────────────────

    def _on_cursor_moved(self, x: float, y: float):
        self._lbl_coords.setText(f"x={x:+.3f}  y={y:+.3f} m")

    def _on_tool_changed(self, name: str):
        self._active_tool = name
        self._lbl_tool.setText(f"tool: {name}")
        # Click-to-add proof of concept: rectangle body
        if name == "body_rect":
            # Place a default body at origin for now (real tool comes later).
            self._add_demo_body()
            self._toolbox.set_tool("select")

    def _on_selection_changed(self):
        items = [it for it in self._scene.selectedItems()
                 if isinstance(it, BodyItem)]
        if items:
            self._inspector.show_item("body", items[0].spec.id)
        # Note: tree refresh on every change is fine for small models;
        # optimise later if needed.

    def _on_new_project(self):
        self._scene.clear()
        self.spec = ModelSpec()
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)
        self._inspector.clear()
        self._update_count_label()

    # ──────────────────────────────────────────────────────────
    # Demo / placeholder
    # ──────────────────────────────────────────────────────────

    def _add_demo_body(self):
        """Create a default rectangular body at the origin (placeholder)."""
        spec = BodySpec(
            name=f"body{len(self.spec.bodies) + 1}",
            mass=1.0,
            inertia=0.01,
            position=(0.0, 0.0),
            shape=ShapeSpec(kind="rectangle",
                            params={"width": 0.20, "height": 0.10}),
        )
        self.spec.bodies.append(spec)
        item = BodyItem(spec)
        self._scene.addItem(item)
        self._tree.refresh()
        self._update_count_label()

    def _update_count_label(self):
        self._lbl_count.setText(
            f"bodies={len(self.spec.bodies)}  "
            f"joints={len(self.spec.joints)}  "
            f"forces={len(self.spec.forces)}"
        )
