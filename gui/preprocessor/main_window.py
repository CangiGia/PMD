"""PreProcessorWindow — Adams-style shell of the model-builder GUI."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QStatusBar,
    QToolBar,
)

from .models import (
    BodySpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    load_model,
    save_model,
)
from .panels import InspectorPanel, RibbonBar, TreePanel
from .tools  import make_tool, Tool
from .widgets import (
    BodyItem,
    CanvasScene,
    CanvasView,
    JointItem,
    MarkerItem,
)


class PreProcessorWindow(QMainWindow):
    """Top-level QMainWindow for the PMD pre-processor."""

    def __init__(self, spec: ModelSpec, parent=None):
        super().__init__(parent)
        self.spec = spec
        self._project_path: str | None = None

        # spec.id → graphics item registries
        self._body_items:   dict[str, BodyItem]   = {}
        self._marker_items: dict[str, MarkerItem] = {}
        self._joint_items:  dict[str, JointItem]  = {}

        self.setWindowTitle("PMD Pre-Processor — Model Builder")
        self.resize(1500, 900)

        self._build_central_canvas()
        self._build_ribbon()
        self._build_docks()
        self._build_menu_bar()
        self._build_status_bar()

        # Wire up
        self._ribbon.tool_changed.connect(self.set_active_tool)
        self._ribbon.action_triggered.connect(self._on_action)
        self._tree.item_selected.connect(self._inspector.show_item)
        self._inspector.spec_changed.connect(self._on_spec_edited)
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)

        # Default tool
        self._active_tool: Tool | None = None
        self.set_active_tool("select")

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
        self._ribbon = RibbonBar(self)
        ribbon_bar = QToolBar("Ribbon", self)
        ribbon_bar.setObjectName("RibbonToolBar")
        ribbon_bar.setMovable(False)
        ribbon_bar.setFloatable(False)
        ribbon_bar.addWidget(self._ribbon)
        self.addToolBar(Qt.TopToolBarArea, ribbon_bar)

    def _build_docks(self):
        self._tree = TreePanel()
        dock_tree = QDockWidget("Model Browser", self)
        dock_tree.setObjectName("DockModelBrowser")
        dock_tree.setWidget(self._tree)
        dock_tree.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_tree)

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
        act_new.setShortcut(QKeySequence.New)
        act_new.triggered.connect(self._on_new_project)
        file_menu.addAction(act_new)

        act_open = QAction("&Open…", self)
        act_open.setShortcut(QKeySequence.Open)
        act_open.triggered.connect(self._on_open)
        file_menu.addAction(act_open)

        act_save = QAction("&Save", self)
        act_save.setShortcut(QKeySequence.Save)
        act_save.triggered.connect(self._on_save)
        file_menu.addAction(act_save)

        act_save_as = QAction("Save &As…", self)
        act_save_as.setShortcut(QKeySequence.SaveAs)
        act_save_as.triggered.connect(self._on_save_as)
        file_menu.addAction(act_save_as)

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
        self._lbl_count  = QLabel("bodies=0  markers=0  joints=0  forces=0")

        sb.addWidget(self._lbl_coords, 1)
        sb.addPermanentWidget(self._lbl_tool)
        sb.addPermanentWidget(self._lbl_count)

    # ──────────────────────────────────────────────────────────
    # Tool dispatch
    # ──────────────────────────────────────────────────────────

    def set_active_tool(self, name: str) -> None:
        """Activate the tool registered under ``name`` (and update UI)."""
        if self._active_tool is not None:
            self._active_tool.deactivate()
        self._active_tool = make_tool(name, self)
        self._scene.active_tool = self._active_tool
        self._active_tool.activate()
        self._lbl_tool.setText(f"tool: {name}")
        self._ribbon.set_tool(name)

    # ──────────────────────────────────────────────────────────
    # Item registry helpers (used by tools)
    # ──────────────────────────────────────────────────────────

    def add_body_item(self, spec: BodySpec) -> BodyItem:
        item = BodyItem(spec)
        self._scene.addItem(item)
        self._body_items[spec.id] = item
        return item

    def add_marker_item(self, spec: MarkerSpec) -> MarkerItem:
        parent_body = self._body_items.get(spec.body_id)
        if parent_body is None:
            # ground (or stale) → place directly in scene
            item = MarkerItem(spec)
            item.setPos(*spec.local_position)
            self._scene.addItem(item)
        else:
            item = MarkerItem(spec, parent_body=parent_body)
        self._marker_items[spec.id] = item
        return item

    def add_joint_visual(self, spec: JointSpec) -> JointItem | None:
        i_marker = self._marker_items.get(spec.i_marker_id)
        if i_marker is None:
            return None
        item = JointItem(spec)
        # place at i-marker world position
        wp = i_marker.scenePos()
        item.setPos(wp.x(), wp.y())
        self._scene.addItem(item)
        self._joint_items[spec.id] = item
        return item

    # ──────────────────────────────────────────────────────────
    # Slots
    # ──────────────────────────────────────────────────────────

    def _on_cursor_moved(self, x: float, y: float):
        self._lbl_coords.setText(f"x={x:+.3f}  y={y:+.3f} m")

    def _on_action(self, name: str):
        # Real handlers will be wired in subsequent steps.
        self.statusBar().showMessage(f"Action: {name}", 2000)

    def _on_spec_edited(self, kind: str, item_id: str) -> None:
        """Re-sync graphics + tree after the inspector mutates a spec."""
        if kind == "body":
            self._refresh_body_visual(item_id)
        elif kind == "marker":
            self._refresh_marker_visual(item_id)
        elif kind == "joint":
            self._refresh_joint_visual(item_id)
        # forces have no visual yet
        self._tree.refresh()

    def _refresh_body_visual(self, body_id: str) -> None:
        import math
        body = self.spec.body(body_id)
        item = self._body_items.get(body_id)
        if body is None or item is None:
            return
        item.setPos(body.position[0], body.position[1])
        item.setRotation(math.degrees(body.orientation))
        # Shape size may have changed:
        params = body.shape.params if body.shape else {}
        if body.shape and body.shape.kind == "link":
            item._w = float(params.get("length",    item._w))
            item._h = float(params.get("thickness", item._h))
        elif body.shape and body.shape.kind == "rectangle":
            item._w = float(params.get("width",  item._w))
            item._h = float(params.get("height", item._h))
        item.prepareGeometryChange()
        item.update()

    def _refresh_marker_visual(self, marker_id: str) -> None:
        import math
        spec = self.spec.marker(marker_id)
        item = self._marker_items.get(marker_id)
        if spec is None or item is None:
            return
        # Re-parent if owning body changed
        new_parent = self._body_items.get(spec.body_id)
        if item.parentItem() is not new_parent:
            self._scene.removeItem(item)
            new_item = self.add_marker_item(spec)
            # Replace in registry (already done by add_marker_item)
            return
        item.setPos(*spec.local_position)
        item.setRotation(math.degrees(spec.theta))
        item.update()

    def _refresh_joint_visual(self, joint_id: str) -> None:
        spec = self.spec.joint(joint_id)
        item = self._joint_items.get(joint_id)
        if spec is None or item is None:
            return
        i_marker = self._marker_items.get(spec.i_marker_id)
        if i_marker is not None:
            wp = i_marker.scenePos()
            item.setPos(wp.x(), wp.y())
        item.update()

    def _on_selection_changed(self):
        items = self._scene.selectedItems()
        if not items:
            return
        top = items[0]
        if isinstance(top, BodyItem):
            self._inspector.show_item("body", top.spec.id)
        elif isinstance(top, MarkerItem):
            self._inspector.show_item("marker", top.spec.id)
        elif isinstance(top, JointItem):
            self._inspector.show_item("joint", top.spec.id)

    def _on_new_project(self):
        self._reset_scene()
        self.spec = ModelSpec()
        self._project_path = None
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)
        self._inspector.clear()
        self._update_count_label()
        self._update_title()

    # -- Project I/O ------------------------------------------------

    PROJECT_FILTER = "PMD Model (*.pmdmodel.json *.json);;All files (*.*)"

    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Project", "", self.PROJECT_FILTER)
        if not path:
            return
        try:
            self.load_project(path)
        except Exception as exc:  # surface to user
            QMessageBox.critical(self, "Open failed", str(exc))

    def _on_save(self):
        if not self._project_path:
            return self._on_save_as()
        self.save_project(self._project_path)

    def _on_save_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Project As", "untitled.pmdmodel.json",
            self.PROJECT_FILTER)
        if not path:
            return
        self.save_project(path)

    def save_project(self, path: str) -> None:
        save_model(self.spec, path)
        self._project_path = path
        self._update_title()
        self.statusBar().showMessage(f"Saved: {path}", 3000)

    def load_project(self, path: str) -> None:
        spec = load_model(path)
        self._reset_scene()
        self.spec = spec
        self._project_path = path
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)
        self._inspector.clear()

        # Rebuild scene from spec (bodies first, then markers, then joints).
        for b in self.spec.bodies:
            self.add_body_item(b)
        for m in self.spec.markers:
            self.add_marker_item(m)
        for j in self.spec.joints:
            self.add_joint_visual(j)

        self._update_count_label()
        self._update_title()
        self._view.zoom_to_fit()
        self.statusBar().showMessage(f"Loaded: {path}", 3000)

    def _reset_scene(self) -> None:
        self._scene.clear()
        self._body_items.clear()
        self._marker_items.clear()
        self._joint_items.clear()

    def _update_title(self) -> None:
        suffix = self._project_path or "<unsaved>"
        self.setWindowTitle(
            f"PMD Pre-Processor \u2014 {self.spec.name} [{suffix}]")

    # ──────────────────────────────────────────────────────────
    def _update_count_label(self):
        self._lbl_count.setText(
            f"bodies={len(self.spec.bodies)}  "
            f"markers={len(self.spec.markers)}  "
            f"joints={len(self.spec.joints)}  "
            f"forces={len(self.spec.forces)}"
        )
