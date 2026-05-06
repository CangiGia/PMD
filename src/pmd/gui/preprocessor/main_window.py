"""PreProcessorWindow — Adams-style shell of the model-builder GUI."""

from __future__ import annotations

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QStatusBar,
    QToolBar,
    QWidget,
)

from .models import (
    BodySpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    load_model,
    save_model,
    to_dict,
    from_dict,
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


# Maximum number of states retained in the undo/redo history.
_HISTORY_CAP = 100


# ----------------------------------------------------------------------
# Solver worker thread
# ----------------------------------------------------------------------
class _SolveWorker(QThread):
    """Run ``model.solve(...)`` off the GUI thread.

    The solver normally renders a tqdm progress bar to stdout. While
    the worker is running we monkey-patch ``pmd.core.solver.tqdm`` with
    a thin shim that emits a Qt ``progress(int)`` signal each time the
    solver updates its bar, so the GUI can show a real
    :class:`QProgressDialog` instead of the (invisible) console bar.
    """

    progress     = Signal(int)              # 0..100
    finished_ok  = Signal(object, object, object)  # (model, T, uT)
    failed       = Signal(str)

    def __init__(self, model, cfg, parent=None):
        super().__init__(parent)
        self._model   = model
        self._cfg     = cfg
        self._cancel  = False

    def request_cancel(self) -> None:
        self._cancel = True

    def run(self) -> None:                 # noqa: D401 (Qt override)
        import numpy as np
        try:
            from pmd.core import solver as _solver_mod
        except Exception as exc:           # pragma: no cover
            self.failed.emit(f"Could not import pmd solver: {exc}")
            return

        worker = self
        original_tqdm = _solver_mod.tqdm

        class _GuiTqdm:
            def __init__(self, total=100, **_kw):
                self.total = max(int(total or 1), 1)
                self._n = 0
                worker.progress.emit(0)

            @property
            def n(self):
                return self._n

            @n.setter
            def n(self, v):
                self._n = v
                pct = int(100 * v / self.total)
                worker.progress.emit(min(100, max(0, pct)))

            def refresh(self):
                pass

            def update(self, k=1):
                self.n = self._n + k

            def close(self):
                worker.progress.emit(100)

            def __enter__(self):
                return self

            def __exit__(self, *_a):
                self.close()

        _solver_mod.tqdm = _GuiTqdm
        try:
            cfg = self._cfg
            t_eval = np.linspace(cfg.t_start, cfg.t_end, cfg.n_steps)
            T, uT = self._model.solve(
                analysis=cfg.analysis,
                method=cfg.method,
                t_span=(cfg.t_start, cfg.t_end),
                t_eval=t_eval,
                ic_correct=cfg.ic_correct,
                baumgarte_alpha=cfg.baumgarte_alpha,
                baumgarte_beta=cfg.baumgarte_beta,
            )
            self.finished_ok.emit(self._model, T, uT)
        except Exception:
            import traceback
            self.failed.emit(traceback.format_exc())
        finally:
            _solver_mod.tqdm = original_tqdm


class PreProcessorWindow(QMainWindow):
    """Top-level QMainWindow for the pmd pre-processor."""

    def __init__(self, spec: ModelSpec, parent=None):
        super().__init__(parent)
        self.spec = spec
        self._project_path: str | None = None
        self._solver_settings = None  # lazy SolverSettings, last used
        self._last_run = None         # (model, T, uT) of the last solve
        self._post_window = None      # keep ref so it isn't GC'd

        # spec.id → graphics item registries
        self._body_items:   dict[str, BodyItem]   = {}
        self._marker_items: dict[str, MarkerItem] = {}
        self._joint_items:  dict[str, JointItem]  = {}
        # 'V' shortcut toggles this; new markers respect the current value.
        self._markers_visible: bool = True

        self.setWindowTitle("pmd Pre-Processor — Model Builder")
        self.resize(1500, 900)

        self._build_central_canvas()
        self._build_ribbon()
        self._build_docks()
        self._build_menu_bar()
        self._build_status_bar()

        # Wire up
        self._ribbon.tool_changed.connect(self.set_active_tool)
        self._ribbon.action_triggered.connect(self._on_action)
        self._tree.item_selected.connect(self._on_tree_item_selected)
        self._tree.item_selected.connect(self._inspector.show_item)
        self._tree.item_delete_requested.connect(self._on_delete_requested)
        self._tree.model_renamed.connect(self._on_model_renamed)
        self._inspector.spec_changed.connect(self._on_spec_edited)
        self._inspector.body_renamed.connect(self._on_body_renamed)
        self._view.delete_pressed.connect(self._on_canvas_delete_shortcut)
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)

        # Undo/redo history (list of JSON-serialisable snapshots).
        self._history: list[dict] = []
        self._history_idx: int = -1
        self._history_record(initial=True)

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
        self._act_undo = QAction("Undo", self)
        self._act_undo.setShortcut(QKeySequence.Undo)
        self._act_undo.triggered.connect(self._undo)
        edit_menu.addAction(self._act_undo)

        self._act_redo = QAction("Redo", self)
        self._act_redo.setShortcuts([QKeySequence.Redo,
                                     QKeySequence("Ctrl+Y")])
        self._act_redo.triggered.connect(self._redo)
        edit_menu.addAction(self._act_redo)

        edit_menu.addSeparator()
        self._act_delete = QAction("Delete selection", self)
        self._act_delete.setShortcut(QKeySequence.Delete)
        self._act_delete.setShortcutContext(Qt.ApplicationShortcut)
        self._act_delete.triggered.connect(self._on_canvas_delete_shortcut)
        edit_menu.addAction(self._act_delete)

        # Global Esc → cancel current tool / revert to Select. Using
        # ApplicationShortcut so it fires regardless of which child
        # widget currently has focus.
        self._act_cancel = QAction("Cancel tool", self)
        self._act_cancel.setShortcut(QKeySequence(Qt.Key_Escape))
        self._act_cancel.setShortcutContext(Qt.ApplicationShortcut)
        self._act_cancel.triggered.connect(
            lambda: self.set_active_tool("select"))
        self.addAction(self._act_cancel)

        view_menu = mbar.addMenu("&View")
        act_fit = QAction("Zoom to &Fit", self)
        act_fit.setShortcut("F")
        act_fit.triggered.connect(self._view.zoom_to_fit)
        view_menu.addAction(act_fit)

        # V \u2192 toggle visibility of all markers (and, transitively,
        # their labels and the joint-name labels). The marker glyphs
        # are visually noisy on dense models; pressing V gives a clean
        # body-only view for screenshots / inspection. Application-
        # scoped so it fires regardless of focus.
        self._act_toggle_markers = QAction("Toggle markers", self)
        self._act_toggle_markers.setShortcut(QKeySequence(Qt.Key_V))
        self._act_toggle_markers.setShortcutContext(Qt.ApplicationShortcut)
        self._act_toggle_markers.triggered.connect(self._toggle_markers_visible)
        self.addAction(self._act_toggle_markers)
        view_menu.addAction(self._act_toggle_markers)

        mbar.addMenu("&Settings")
        mbar.addMenu("&Tools")

    def _build_status_bar(self):
        sb = QStatusBar(self)
        self.setStatusBar(sb)
        # Slightly larger font so prompts ("pick first marker…") are
        # easy to read while the user's eye is on the canvas.
        sb.setStyleSheet(
            "QStatusBar { font-size: 12pt; }"
            "QStatusBar QLabel { font-size: 11pt; }"
        )

        self._lbl_coords = QLabel("x=0.000  y=0.000 m")
        self._lbl_tool   = QLabel("tool: select")
        self._lbl_count  = QLabel("bodies=0  markers=0  joints=0  forces=0")

        # Inline solver-progress group: a description label + a thin
        # progress bar + a small Cancel button. Lives in the status
        # bar (left side, just after the coordinate readout) so the
        # solver no longer steals focus with a modal popup. Hidden
        # while idle.
        self._sim_status_widget = QWidget()
        _sim_lay = QHBoxLayout(self._sim_status_widget)
        _sim_lay.setContentsMargins(0, 0, 0, 0)
        _sim_lay.setSpacing(8)
        self._lbl_sim = QLabel("")
        self._lbl_sim.setStyleSheet("QLabel { font-size: 11pt; }")
        self._sim_progress = QProgressBar()
        self._sim_progress.setRange(0, 100)
        self._sim_progress.setValue(0)
        self._sim_progress.setTextVisible(True)
        self._sim_progress.setFixedHeight(14)
        self._sim_progress.setFixedWidth(220)
        self._sim_cancel = QPushButton("Cancel")
        self._sim_cancel.setFixedHeight(20)
        self._sim_cancel.setFlat(True)
        _sim_lay.addWidget(self._lbl_sim)
        _sim_lay.addWidget(self._sim_progress)
        _sim_lay.addWidget(self._sim_cancel)
        self._sim_status_widget.hide()

        sb.addWidget(self._lbl_coords, 1)
        sb.addWidget(self._sim_status_widget)
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
        # Give keyboard focus to the canvas so Esc (and other shortcuts)
        # reach the active tool without requiring a click first.
        self._view.setFocus(Qt.OtherFocusReason)

    # ──────────────────────────────────────────────────────────
    # Item registry helpers (used by tools)
    # ──────────────────────────────────────────────────────────

    def add_body_item(self, spec: BodySpec) -> BodyItem | None:
        # The implicit ground body has no visual representation on the
        # canvas — it lives only in the tree as an attachment target.
        if spec.id == "ground":
            return None
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
        if not self._markers_visible:
            item.setVisible(False)
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
        if name == "sim_run":
            self._run_simulation()
        elif name == "sim_animate":
            self._open_animation_window()
        elif name == "sim_settings":
            # Legacy key (kept for backwards compatibility).
            self._open_solver_dialog()
        elif name == "post_open":
            self._open_post_processor()
        elif name == "force_grav":
            self._toggle_gravity()
        else:
            self.statusBar().showMessage(f"Action: {name}", 2000)

    def _toggle_gravity(self):
        """Add (or remove) a global Weight ForceSpec."""
        from .models import ForceSpec
        existing = next(
            (f for f in self.spec.forces if f.kind == "Weight"), None)
        if existing is not None:
            self.spec.forces.remove(existing)
            self.statusBar().showMessage("Gravity: OFF", 2000)
        else:
            self.spec.forces.append(ForceSpec(
                name="gravity", kind="Weight",
                params={"gravity": 9.80665},
            ))
            self.statusBar().showMessage("Gravity: ON  (9.80665 m/s², -Y)", 2000)
        self._tree.refresh()
        self._update_count_label()

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
        self._history_record()

    def _on_model_renamed(self, new_name: str) -> None:
        """User renamed the model from the tree's name editor.

        Updates :attr:`spec.name` (which feeds the window title and the
        default file name on Save As) and snapshots the change so it
        can be undone like any other edit.
        """
        if not new_name or new_name == self.spec.name:
            return
        self.spec.name = new_name
        self._update_title()
        self._history_record()

    def _on_body_renamed(self, body_id: str, old_name: str, new_name: str):
        """Rename markers prefixed with the body's old name to track it.

        Markers whose name starts with ``"<old_name>."`` (e.g. the
        auto-CM marker ``"<old_name>.cm"`` or end markers
        ``"<old_name>.A"`` / ``".B"``) are updated so they stay in sync
        with the body label shown in the tree.
        """
        if not old_name or new_name == old_name:
            return
        prefix_old = f"{old_name}."
        prefix_new = f"{new_name}."
        for m in self.spec.markers:
            if m.body_id != body_id:
                continue
            if m.name == old_name:
                m.name = new_name
            elif m.name.startswith(prefix_old):
                m.name = prefix_new + m.name[len(prefix_old):]
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
            # Nothing is selected anymore → reset the Inspector so the
            # last-selected entity's editor doesn't linger.
            self._inspector.clear()
            return
        top = items[0]
        if isinstance(top, BodyItem):
            self._inspector.show_item("body", top.spec.id)
        elif isinstance(top, MarkerItem):
            self._inspector.show_item("marker", top.spec.id)
        elif isinstance(top, JointItem):
            self._inspector.show_item("joint", top.spec.id)

    def _on_tree_item_selected(self, kind: str, item_id: str) -> None:
        """Bridge tree clicks into the active tool when relevant.

        Currently used so that, while the JointTool is active, the
        user can pick the two markers either by clicking on the canvas
        or by simply clicking them in the Model Browser tree.

        In addition, mirror the tree selection on the canvas so the
        clicked body / marker / joint visually highlights itself
        (BodyItem fill change, MarkerItem dashed ring + bright dot,
        JointItem orange halo). The scene's :pysignal:`selectionChanged`
        is temporarily blocked so we don't bounce back into the
        inspector / tree.
        """
        # Resolve the canvas item, if any, for this tree row.
        item = None
        if kind == "body":
            item = self._body_items.get(item_id)
        elif kind == "marker":
            item = self._marker_items.get(item_id)
        elif kind == "joint":
            item = self._joint_items.get(item_id)

        if item is not None:
            scene = self._view.scene()
            scene.blockSignals(True)
            try:
                # Single-selection semantics: clear previous, select
                # the just-clicked item, then ensure it's repainted
                # (selection halos are computed in paint()).
                for it in scene.selectedItems():
                    if it is not item:
                        it.setSelected(False)
                item.setSelected(True)
                item.update()
            finally:
                scene.blockSignals(False)

        if kind != "marker" or self._active_tool is None:
            return
        picker = getattr(self._active_tool, "pick_marker_by_id", None)
        if callable(picker):
            picker(item_id)

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

    PROJECT_FILTER = "pmd Model (*.pmdmodel.json *.json);;All files (*.*)"

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
        # Seed the dialog with the model's current name so the saved
        # file matches what the user sees in the tree's name editor.
        stem = (self.spec.name or "untitled").strip().replace(" ", "_") \
               or "untitled"
        default = f"{stem}.pmdmodel.json"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Project As", default,
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
            f"pmd Pre-Processor \u2014 {self.spec.name} [{suffix}]")

    # ──────────────────────────────────────────────────────────
    def _update_count_label(self):
        self._lbl_count.setText(
            f"bodies={len(self.spec.bodies)}  "
            f"markers={len(self.spec.markers)}  "
            f"joints={len(self.spec.joints)}  "
            f"forces={len(self.spec.forces)}"
        )

    # ──────────────────────────────────────────────────────────
    # Auto-marker / body creation hooks (called by tools)
    # ──────────────────────────────────────────────────────────

    def auto_create_cm_marker(self, body: BodySpec) -> MarkerSpec:
        """Create the canonical centre-of-mass marker for ``body``.

        Tools that commit a new body should call this so every body
        always has a ``"<name>.cm"`` marker at its origin.
        """
        cm = MarkerSpec(name=f"{body.name}.cm", body_id=body.id,
                        local_position=(0.0, 0.0), theta=0.0)
        self.spec.markers.append(cm)
        self.add_marker_item(cm)
        return cm

    # ──────────────────────────────────────────────────────────
    # Deletion
    # ──────────────────────────────────────────────────────────

    def _on_delete_requested(self, kind: str, item_id: str) -> None:
        self._delete_spec(kind, item_id)

    def _on_canvas_delete_shortcut(self) -> None:
        """Delete every currently selected scene item."""
        items = self._scene.selectedItems()
        if not items:
            return
        targets: list[tuple[str, str]] = []
        for it in items:
            if isinstance(it, BodyItem):
                targets.append(("body", it.spec.id))
            elif isinstance(it, MarkerItem):
                targets.append(("marker", it.spec.id))
            elif isinstance(it, JointItem):
                targets.append(("joint", it.spec.id))
        for kind, _id in targets:
            self._delete_spec(kind, _id, _record=False)
        self._history_record()

    def _delete_spec(self, kind: str, item_id: str,
                     *, _record: bool = True) -> None:
        """Remove a spec (and its dependents) from the model + scene."""
        if kind == "body":
            # The implicit ground body cannot be deleted.
            if item_id == "ground":
                self.statusBar().showMessage(
                    "The ground body is implicit and cannot be removed.",
                    4000)
                return
            body = self.spec.body(item_id)
            if body is None:
                return
            # Cascade: remove markers belonging to the body, then any
            # joints / forces referencing those markers.
            mids = {m.id for m in self.spec.markers if m.body_id == item_id}
            self.spec.markers = [m for m in self.spec.markers
                                 if m.body_id != item_id]
            self.spec.joints = [j for j in self.spec.joints
                                if j.i_marker_id not in mids
                                and j.j_marker_id not in mids]
            self.spec.forces = [f for f in self.spec.forces
                                if f.i_marker_id not in mids
                                and f.j_marker_id not in mids]
            self.spec.bodies = [b for b in self.spec.bodies
                                if b.id != item_id]
            # Scene cleanup
            for mid in mids:
                self._scene_remove_marker(mid)
            self._scene_remove_body(item_id)
            self._scene_remove_joints_orphans()
        elif kind == "marker":
            marker = self.spec.marker(item_id)
            if marker is None:
                return
            self.spec.markers = [m for m in self.spec.markers
                                 if m.id != item_id]
            self.spec.joints = [j for j in self.spec.joints
                                if j.i_marker_id != item_id
                                and j.j_marker_id != item_id]
            self.spec.forces = [f for f in self.spec.forces
                                if f.i_marker_id != item_id
                                and f.j_marker_id != item_id]
            self._scene_remove_marker(item_id)
            self._scene_remove_joints_orphans()
        elif kind == "joint":
            self.spec.joints = [j for j in self.spec.joints
                                if j.id != item_id]
            self._scene_remove_joint(item_id)
        elif kind == "force":
            self.spec.forces = [f for f in self.spec.forces
                                if f.id != item_id]
        else:
            return

        self._inspector.clear()
        self._tree.refresh()
        self._update_count_label()
        if _record:
            self._history_record()

    # ── Scene-removal helpers ─────────────────────────────────
    def _scene_remove_body(self, body_id: str) -> None:
        item = self._body_items.pop(body_id, None)
        if item is not None and item.scene() is not None:
            self._scene.removeItem(item)

    def _scene_remove_marker(self, marker_id: str) -> None:
        item = self._marker_items.pop(marker_id, None)
        if item is not None and item.scene() is not None:
            self._scene.removeItem(item)

    def _scene_remove_joint(self, joint_id: str) -> None:
        item = self._joint_items.pop(joint_id, None)
        if item is not None and item.scene() is not None:
            self._scene.removeItem(item)

    def _scene_remove_joints_orphans(self) -> None:
        live = {j.id for j in self.spec.joints}
        for jid in [j for j in self._joint_items if j not in live]:
            self._scene_remove_joint(jid)

    # ──────────────────────────────────────────────────────────
    # Undo / Redo
    # ──────────────────────────────────────────────────────────

    def _history_record(self, *, initial: bool = False) -> None:
        """Snapshot the current spec to the undo history."""
        snap = to_dict(self.spec)
        if (not initial and self._history
                and self._history_idx >= 0
                and self._history[self._history_idx] == snap):
            return  # no-op edit
        # Drop redo tail.
        del self._history[self._history_idx + 1:]
        self._history.append(snap)
        # Cap history to avoid unbounded memory growth.
        if len(self._history) > _HISTORY_CAP:
            self._history = self._history[-_HISTORY_CAP:]
        self._history_idx = len(self._history) - 1

    def _undo(self) -> None:
        if self._history_idx <= 0:
            return
        self._history_idx -= 1
        self._restore_snapshot(self._history[self._history_idx])

    def _redo(self) -> None:
        if self._history_idx >= len(self._history) - 1:
            return
        self._history_idx += 1
        self._restore_snapshot(self._history[self._history_idx])

    def _restore_snapshot(self, snap: dict) -> None:
        spec = from_dict(snap)
        self._reset_scene()
        self.spec = spec
        self._tree.set_spec(self.spec)
        self._inspector.set_spec(self.spec)
        for b in self.spec.bodies:
            self.add_body_item(b)
        for m in self.spec.markers:
            self.add_marker_item(m)
        for j in self.spec.joints:
            self.add_joint_visual(j)
        self._update_count_label()
        self.statusBar().showMessage(
            f"History: {self._history_idx + 1}/{len(self._history)}", 2000)

    # ──────────────────────────────────────────────────────────
    # Simulation bridge
    # ──────────────────────────────────────────────────────────

    def _toggle_markers_visible(self):
        # Defensive: callable before items exist.
        self._markers_visible = not getattr(self, "_markers_visible", True)
        for item in self._marker_items.values():
            item.setVisible(self._markers_visible)
        self.statusBar().showMessage(
            "Markers shown" if self._markers_visible else "Markers hidden",
            1500)

    def _open_solver_dialog(self):
        from .dialogs import SolverDialog
        dlg = SolverDialog(self, initial=self._solver_settings)
        if dlg.exec() == dlg.Accepted:
            self._solver_settings = dlg.settings()
            self.statusBar().showMessage("Solver settings updated", 2000)

    def _run_simulation(self, *, analysis: str | None = None):
        from .builder import build_model, BuilderError
        from .dialogs import SolverDialog, SolverSettings

        if not self.spec.bodies or all(b.id == "ground"
                                       for b in self.spec.bodies):
            QMessageBox.warning(self, "Run", "Model has no bodies.")
            return

        # Show dialog seeded with the previous settings (or defaults);
        # the dialog itself lets the user pick the analysis kind.
        seed = self._solver_settings or SolverSettings()
        if analysis:
            seed.analysis = analysis
        dlg = SolverDialog(self, initial=seed)
        from PySide6.QtWidgets import QDialog
        if dlg.exec() != QDialog.Accepted:
            return
        cfg = dlg.settings()
        self._solver_settings = cfg

        try:
            model = build_model(self.spec)
        except BuilderError as exc:
            QMessageBox.critical(self, "Build failed", str(exc))
            return

        # Run the solver in a worker thread; show an inline progress
        # bar in the status bar (next to the simulation-type label).
        # The solver's stdout-based tqdm bar is intercepted (by
        # monkey-patching ``solver.tqdm``) so the user sees the GUI
        # bar instead of a console bar that's invisible while the GUI
        # is in the foreground.
        worker = _SolveWorker(model, cfg)

        # Show inline progress group in the status bar.
        self._lbl_sim.setText(f"Solving ({cfg.analysis}, {cfg.method})\u2026")
        self._sim_progress.setValue(0)
        self._sim_status_widget.show()
        try:
            self._sim_cancel.clicked.disconnect()
        except (TypeError, RuntimeError):
            pass
        self._sim_cancel.clicked.connect(worker.request_cancel)
        worker.progress.connect(self._sim_progress.setValue)

        def _hide_sim_status() -> None:
            self._sim_status_widget.hide()
            self._lbl_sim.setText("")
            self._sim_progress.setValue(0)

        def _on_done(model_, T, uT):
            self._sim_progress.setValue(100)
            _hide_sim_status()
            self._last_run = (model_, T, uT)
            self.statusBar().showMessage(
                f"Solve OK: {len(T)} steps. Click Animate to replay, "
                "or Results \u2192 Open for plots.", 6000)

        def _on_failed(msg: str):
            _hide_sim_status()
            QMessageBox.critical(self, "Solve failed", msg)
            self.statusBar().clearMessage()

        worker.finished_ok.connect(_on_done)
        worker.failed.connect(_on_failed)
        # Keep a reference so the QThread isn't GC'd mid-run.
        self._sim_worker = worker
        worker.start()

    def _open_animation_window(self):
        if self._last_run is None:
            QMessageBox.information(
                self, "Animate",
                "No simulation results yet \u2014 click Run Simulation first.")
            return
        from .dialogs import AnimationDialog

        model, T, _uT = self._last_run
        try:
            dlg = AnimationDialog(self.spec, model, T, parent=self)
        except Exception:
            import traceback
            QMessageBox.critical(
                self, "Animate failed",
                f"Could not build the animation:\n\n{traceback.format_exc()}")
            return
        dlg.show()
        # Keep a reference so the dialog isn't GC'd when the function
        # returns; storing on self also lets a second click reuse it.
        self._anim_dialog = dlg

    def _open_post_processor(self):
        if self._last_run is None:
            QMessageBox.information(
                self, "Post-processor",
                "No simulation results yet \u2014 run an analysis first.")
            return
        model, T, uT = self._last_run
        # Build the post-processor's main window in our running event loop;
        # avoid PostProcessor.show() (which would call app.exec() again).
        from ..postprocessor import MainWindow as _PostMain, Session
        model._distribute_results(T, uT)
        sessions = [Session(model, T, uT, name=self.spec.name,
                            preprocessor_spec=self.spec)]
        self._post_window = _PostMain(sessions)
        self._post_window.show()
