"""MainWindow — primary QMainWindow shell for the pmd PostProcessor."""

import csv
import logging

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QMainWindow,
    QMessageBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from pmd.core.units import UnitSystem
from .models import build_curves
from .panels import FilterPanel, ResultSetPanel, SimulationPanel, UnitsToolbar
from ..style import apply_dark_theme, apply_light_theme
from .. import icons as _icons
from .widgets import AnimationCanvas, PlotCanvas

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """Top-level window: menu bar, left panel | central area, status bar.

    Parameters
    ----------
    sessions : list[Session]
        One or more loaded simulation sessions to display.
    parent : QWidget or None
        Optional parent widget.
    """

    def __init__(self, sessions, parent=None):
        super().__init__(parent)
        self._sessions = sessions
        self._selection = []
        self._display_units = UnitSystem()

        self.setWindowTitle("pmd PostProcessor")
        self.resize(1500, 720)

        self._build_central_area()
        self._build_status_bar()

    # Central area (splitter: NavigationPanel | FilterPanel + plot)

    def _build_central_area(self):
        splitter = QSplitter(Qt.Orientation.Horizontal, self)

        # Left panel — SimulationPanel
        self._sim_panel = SimulationPanel(self._sessions)
        self._sim_panel.setObjectName("sim_sidebar")
        self._sim_panel.setMinimumWidth(200)
        self._sim_panel.setMaximumWidth(350)
        self._sim_panel.selection_changed.connect(self._on_selection_changed)
        self._sim_panel.theme_toggle_requested.connect(self._on_toggle_theme)
        self._sim_panel.export_requested.connect(self._on_export_requested)
        self._sim_panel.close_requested.connect(self.close)
        self._sim_panel.animation_toggle_requested.connect(self._on_toggle_animation)
        self._sim_panel.math_toggle_requested.connect(self._on_toggle_math)
        self._sim_panel.joint_glyph_requested.connect(self._on_joint_glyph_size)

        # Right side — FilterPanel on top, vertical splitter (plot | animation), ResultSetPanel bottom
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(0)

        self._units_toolbar = UnitsToolbar()
        self._units_toolbar.units_changed.connect(self._on_units_changed)
        right_layout.addWidget(self._units_toolbar)

        self._filter_panel = FilterPanel()
        self._filter_panel.setObjectName("filter_card")
        self._filter_panel.add_curves_requested.connect(self._on_add_curves)
        right_layout.addWidget(self._filter_panel)

        # The animation pane sits *next to* the plots (horizontal split)
        # rather than below them: vertical stacking pushed the plots
        # into a tiny strip on common laptop screens.
        viz_splitter = QSplitter(Qt.Orientation.Horizontal)
        self._plot_canvas = PlotCanvas()
        self._plot_canvas.curve_computed.connect(self._on_curve_computed)
        viz_splitter.addWidget(self._plot_canvas)

        # Math bar (Feature F) — lives here so it spans the full right area
        # (plot canvas + animation canvas), keeping both vertically aligned.
        # Inserted AFTER self._plot_canvas is created.
        # Hidden until the user enables it via View → Math bar.
        self._math_bar_widget = self._plot_canvas._math_bar
        right_layout.addWidget(self._math_bar_widget)
        self._math_bar_widget.setVisible(False)
        # When all sessions were launched from the preprocessor, use the
        # preprocessor-style animation canvas (real body shapes, marker
        # axes, fuchsia revolute discs) instead of the plain coloured
        # circles fallback. Falls back to the legacy AnimationCanvas
        # otherwise (e.g. sessions loaded from a .pkl).
        specs = [getattr(s, "preprocessor_spec", None) for s in self._sessions]
        if specs and all(sp is not None for sp in specs):
            from ..preprocessor.dialogs import PreprocessorAnimationCanvas
            s0 = self._sessions[0]
            self._anim_canvas = PreprocessorAnimationCanvas(
                specs[0], s0.model, s0.T)
        else:
            self._anim_canvas = AnimationCanvas(self._sessions)
        viz_splitter.addWidget(self._anim_canvas)
        self._anim_canvas.setVisible(False)
        # Give AnimationCanvas a reference to PlotCanvas so that the video
        # export backend can draw synchronised plot frames (combo layout).
        if hasattr(self._anim_canvas, "_export_video"):
            self._anim_canvas._plot_canvas_ref = self._plot_canvas
        if hasattr(self._anim_canvas, "set_step"):
            self._plot_canvas.step_requested.connect(self._anim_canvas.set_step)
        # When the preprocessor-style canvas is in use, mirror its
        # current frame as a synchronised vertical cursor on every plot.
        if hasattr(self._anim_canvas, "time_changed"):
            self._anim_canvas.time_changed.connect(
                self._plot_canvas.set_time_cursor)
        # Cursor lock: PlotCanvas ↔ AnimationCanvas (both directions)
        if (hasattr(self._plot_canvas, "cursor_lock_changed")
                and hasattr(self._anim_canvas, "set_cursor_locked")):
            self._plot_canvas.cursor_lock_changed.connect(
                self._on_plot_cursor_lock_changed)
        if (hasattr(self._anim_canvas, "cursor_time_changed")
                and hasattr(self._plot_canvas, "force_time_cursor")):
            self._anim_canvas.cursor_time_changed.connect(
                self._plot_canvas.force_time_cursor)
        if (hasattr(self._anim_canvas, "request_cursor_unlock")
                and hasattr(self._plot_canvas, "unlock_time_cursor")):
            self._anim_canvas.request_cursor_unlock.connect(
                self._plot_canvas.unlock_time_cursor)
        viz_splitter.setStretchFactor(0, 2)
        viz_splitter.setStretchFactor(1, 3)
        self._viz_splitter = viz_splitter
        right_layout.addWidget(viz_splitter, stretch=1)

        self._result_set_panel = ResultSetPanel()
        self._result_set_panel.setMaximumHeight(200)
        self._result_set_panel.curves_changed.connect(self._on_curves_changed)
        right_layout.addWidget(self._result_set_panel)

        splitter.addWidget(self._sim_panel)
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.setCentralWidget(splitter)

    # Status bar

    def _build_status_bar(self):
        n_bodies = sum(len(s.model.Bodies) for s in self._sessions)
        n_joints = sum(len(s.model.Joints) for s in self._sessions)
        n_sessions = len(self._sessions)

        parts = [f"Sessions: {n_sessions}",
                 f"Bodies: {n_bodies}",
                 f"Joints: {n_joints}"]

        if n_sessions == 1:
            T = self._sessions[0].T
            n_steps = len(T)
            t_range = f"{T[0]:.4g} – {T[-1]:.4g} s"
            parts += [f"Steps: {n_steps}", f"T: {t_range}"]

        self._status_base = "  |  ".join(parts)
        self.statusBar().showMessage(self._status_base)

    # Slots

    def _on_plot_cursor_lock_changed(self, locked: bool) -> None:
        """Propagate PlotCanvas cursor lock state to AnimationCanvas."""
        if hasattr(self._anim_canvas, "set_cursor_locked"):
            t = None
            if locked and hasattr(self._plot_canvas, "time_cursor_t"):
                t = self._plot_canvas.time_cursor_t
            self._anim_canvas.set_cursor_locked(locked, t)

    def _on_selection_changed(self, selection):
        """Receives list of checked descriptor dicts from SimulationPanel."""
        self._selection = selection
        self._filter_panel.update_from_selection(selection)
        msg = self._status_base + f"  |  Selected: {len(selection)}"
        self.statusBar().showMessage(msg)
        logger.debug("Selection changed: %s", [d["label"] for d in selection])

    def _on_add_curves(self, cat_keys, comp_keys, selection):
        """Build CurveItems for each (cat, comp) pair and add them to the panel."""
        added = 0
        for cat_key in cat_keys:
            for comp_key in comp_keys:
                curves = build_curves(cat_key, comp_key, selection,
                                      self._display_units)
                if curves:
                    self._result_set_panel.add_curves_with_request(
                        curves, cat_key, comp_key, selection)
                    added += len(curves)
        if added:
            self._sim_panel.clear_selection()
        logger.debug(
            "Add curves: cats=%s, comps=%s, added=%d",
            cat_keys, comp_keys, added,
        )

    def _on_units_changed(self, unit_system: UnitSystem):
        """Rebuild all curves whenever the user changes the display unit system."""
        self._display_units = unit_system
        self._result_set_panel.rebuild_all(unit_system, build_curves)
        # rebuild_all emits curves_changed which triggers _on_curves_changed automatically

    def _on_curves_changed(self):
        """Refresh the plot with currently visible curves."""
        visible = self._result_set_panel.visible_curves()
        self._plot_canvas.update_plot(visible)
        logger.debug("Curves changed: %d visible", len(visible))

    def _on_toggle_animation(self, checked: bool):
        """Show/hide the AnimationCanvas."""
        self._anim_canvas.setVisible(checked)
        if checked:
            # Give the animation pane ~60% of the horizontal space and
            # the plot panel ~40%.  Do this every time the pane is
            # shown so the layout resets even if the user had manually
            # dragged the splitter to a different position.
            total = self._viz_splitter.width()
            if total > 0:
                self._viz_splitter.setSizes([
                    int(total * 2 / 5),
                    int(total * 3 / 5),
                ])
        # Show the time cursor on the plots only while the animation
        # pane is visible -- the two are conceptually one feature.
        if hasattr(self._plot_canvas, "set_time_cursor_visible"):
            self._plot_canvas.set_time_cursor_visible(checked)
        if checked and hasattr(self._anim_canvas, "_T") \
                and hasattr(self._anim_canvas, "_step") \
                and hasattr(self._plot_canvas, "set_time_cursor"):
            try:
                self._plot_canvas.set_time_cursor(
                    float(self._anim_canvas._T[self._anim_canvas._step]))
            except Exception:
                pass
        self._sim_panel.set_animation(checked)

    def _on_toggle_math(self, checked: bool) -> None:
        """Show/hide the math bar (View → Math bar)."""
        self._math_bar_widget.setVisible(checked)

    def _on_export_requested(self, kind: str):
        """Dispatch export actions triggered from the sidebar File menu."""
        if kind == "plot":
            self._on_export_plot()
        elif kind == "csv":
            self._on_export_csv()
        elif kind == "txt":
            self._on_export_txt()
        elif kind == "all":
            self._on_export_all()

    def _on_toggle_theme(self, checked: bool):
        """Switch both Qt widgets and Matplotlib between dark and light theme."""
        app = QApplication.instance()
        if checked:
            apply_dark_theme(app)
        else:
            apply_light_theme(app)
        _icons.set_dark(checked)
        self._sim_panel.set_dark(checked)
        # Refresh panel icons
        self._sim_panel.refresh_icons()
        self._result_set_panel.refresh_icons()
        self._filter_panel.refresh_icons()
        # Refresh toolbar icons on both canvases
        self._plot_canvas.set_icon_theme(checked)
        self._anim_canvas.set_icon_theme(checked)
        # Redraw canvases — set_dark preserves curves AND any active zoom insets
        self._plot_canvas.set_dark(checked)
        self._anim_canvas.set_dark(checked)

    def _on_export_plot(self):
        """Save the current plot figure to PNG/SVG/PDF."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "",
            "PNG Image (*.png);;SVG Vector (*.svg);;PDF Document (*.pdf)"
        )
        if path:
            self._plot_canvas._figure.savefig(path, bbox_inches="tight", dpi=150)

    def _on_export_csv(self):
        """Export all visible curves to a single CSV file."""
        curves = self._result_set_panel.visible_curves()
        if not curves:
            QMessageBox.information(self, "Export CSV", "No visible curves to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Curves", "", "CSV file (*.csv)"
        )
        if not path:
            return
        T_ref = curves[0].T
        header = ["time_s"] + [c.label for c in curves]
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for i, t in enumerate(T_ref):
                row = [t] + [float(c.data[i]) if i < len(c.data) else "" for c in curves]
                writer.writerow(row)

    def _on_export_txt(self):
        """Export all visible curves to a tab-separated text file."""
        curves = self._result_set_panel.visible_curves()
        if not curves:
            QMessageBox.information(self, "Export TXT", "No visible curves to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Curves", "", "Text file (*.txt)"
        )
        if not path:
            return
        T_ref = curves[0].T
        cols = ["time_s"] + [c.label for c in curves]
        with open(path, "w", encoding="utf-8") as f:
            f.write("\t".join(cols) + "\n")
            for i, t in enumerate(T_ref):
                row = [str(t)] + [
                    str(float(c.data[i])) if i < len(c.data) else ""
                    for c in curves
                ]
                f.write("\t".join(row) + "\n")

    def _on_export_all(self):
        """Export all loaded curves (visible and hidden) to a tab-separated text file."""
        curves = self._result_set_panel._curves
        if not curves:
            QMessageBox.information(self, "Export", "No curves loaded to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export All Curves", "", "Text file (*.txt)"
        )
        if not path:
            return
        T_ref = curves[0].T
        cols = ["time_s"] + [c.label for c in curves]
        with open(path, "w", encoding="utf-8") as f:
            f.write("\t".join(cols) + "\n")
            for i, t in enumerate(T_ref):
                row = [str(t)] + [
                    str(float(c.data[i])) if i < len(c.data) else ""
                    for c in curves
                ]
                f.write("\t".join(row) + "\n")

    def _on_curve_computed(self, curve) -> None:
        """Add a math-computed curve (from PlotCanvas) to the result set panel."""
        self._result_set_panel.add_curves([curve])

    def _on_joint_glyph_size(self) -> None:
        """Open a dialog to set the joint glyph radius in the animation canvas."""
        if not hasattr(self._anim_canvas, "_rebuild_joint_glyphs"):
            QMessageBox.information(
                self, "Joint glyph size",
                "Joint glyph resize is not available for the current animation canvas.")
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Joint glyph size")
        form = QFormLayout(dlg)

        spin = QDoubleSpinBox()
        spin.setDecimals(4)
        spin.setRange(1e-6, 1e6)
        spin.setSingleStep(0.01)
        current = getattr(self._anim_canvas, "_joint_radius",
                          getattr(self._anim_canvas, "_triad_L", 0.1) * 0.10)
        spin.setValue(current)
        form.addRow("Glyph size (model units):", spin)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        form.addRow(buttons)

        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._anim_canvas._rebuild_joint_glyphs(spin.value())
