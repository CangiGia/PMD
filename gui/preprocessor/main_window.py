"""PreProcessorWindow — Adams-style shell of the model-builder GUI.

Layout
------

::

    ┌─────────────────────────────────────────────────────────────┐
    │  Menu bar                                                   │
    ├─────────────────────────────────────────────────────────────┤
    │  Ribbon  [Bodies | Connectors | Motions | Forces | …]       │
    ├──────────────┬──────────────────────────────┬───────────────┤
    │  Tree        │                              │               │
    │  browser     │      Canvas (QGraphicsView)  │  Inspector    │
    │  (dock L)    │                              │  (dock R)     │
    │              │                              │               │
    ├──────────────┴──────────────────────────────┴───────────────┤
    │  Status bar:  x=…  y=…  | tool: …    | counts               │
    └─────────────────────────────────────────────────────────────┘
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QDockWidget,
    QLabel,
    QMainWindow,
    QStatusBar,
    QToolBar,
)

from .models import BodySpec, ModelSpec, ShapeSpec
from .panels import InspectorPanel, RibbonBar, TreePanel
from .widgets import BodyItem, CanvasScene, CanvasView


class PreProcessorWindow(QMainWindow):
    """Top-level QMainWindow for the PMD pre-processor."""

    def __init__(self, spec: ModelSpec, parent=None):
        super().__init__(parent)
        self.spec = spec
        self._active_tool = "select"

        self.setWindowTitle("PMD Pre-Processor — Model Builder")
        self.resize(1500, 900)

        self._build_central_canvas()
        self._build_ribbon()
        self._build_docks()
        self._build_menu_bar()
        self._build_status_bar()

        # Wire up
        self._ribbon.tool_changed.connect(self._on_tool_changed)
        self._ribbon.action_triggered.connect(self._on_action)
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

    def _build_ribbon(self):
        """Place the ribbon as a non-movable top toolbar."""
        self._ribbon = RibbonBar(self)
        ribbon_bar = QToolBar("Ribbon", self)
        ribbon_bar.setObjectName("RibbonToolBar")
        ribbon_bar.setMovable(False)
        ribbon_bar.setFloatable(False)
        ribbon_bar.addWidget(self._ribbon)
        self.addToolBar(Qt.TopToolBarArea, ribbon_bar)

    def _build_docks(self):
        # Tree browser (left)
        self._tree = TreePanel()
        dock_tree = QDockWidget("Model Browser", self)
        dock_tree.setObjectName("DockModelBrowser")
        dock_tree.setWidget(self._tree)
        dock_tree.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_tree)

        # Inspector (right)
        self._inspector = InspectorPanel()
        dock_inspect = QDockWidget("Inspector", self)
        dock_inspect.setObjectName("DockInspector")
        dock_inspect.setWidget(self._inspector)
        dock_inspect.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_inspect)

        self.resizeDocks([dock_tree, dock_inspect], [260, 280], Qt.Horizontal)

    def _build_menu_bar(self):
        mbar = self.menuBar()

        file_menu = mbar.addMenu("&File")
        act_new = QAction("&New Project", self)
        act_new.triggered.connect(self._on_new_project)
        file_menu.addAction(act_new)
        file_menu.addSeparator()
        act_quit = QAction("&Quit", self)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        edit_menu = mbar.addMenu("&Edit")
        act_undo = QAction("Undo", self); act_undo.setEnabled(False)
        act_redo = QAction("Redo", self); act_redo.setEnabled(False)
        edit_menu.addAction(act_undo)
        edit_menu.addAction(act_redo)

        view_menu = mbar.addMenu("&View")
        act_fit = QAction("Zoom to &Fit", self)
        act_fit.setShortcut("F")
        act_fit.triggered.connect(self._view.zoom_to_fit)
        view_menu.addAction(act_fit)

        mbar.addMenu("&Settings")
        mbar.addMenu("&Tools")

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
        # TEMP: rectangular body proof-of-concept
        if name == "body_rect":
            self._add_demo_body()
            self._ribbon.set_tool("select")

    def _on_action(self, name: str):
        # Real handlers will be wired in subsequent steps.
        self.statusBar().showMessage(f"Action: {name}", 2000)

    def _on_selection_changed(self):
        items = [it for it in self._scene.selectedItems()
                 if isinstance(it, BodyItem)]
        if items:
            self._inspector.show_item("body", items[0].spec.id)

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
