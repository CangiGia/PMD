"""PlotCanvas — matplotlib Figure embedded in a PySide6 widget."""

import numpy as np

import matplotlib
matplotlib.use("qtagg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
from PySide6.QtCore import Signal, Qt

from ._zoom_inset import ZoomInset
from PySide6.QtWidgets import QToolBar, QToolButton, QVBoxLayout, QWidget

from ..models import CurveItem
from ... import icons as _icons
from ...style import CANVAS_BG_DARK, CANVAS_BG_LIGHT, CANVAS_FG_DARK, CANVAS_FG_LIGHT


class PlotCanvas(QWidget):
    """Embeds a matplotlib Figure + NavigationToolbar.

    Public API
    ----------
    update_plot(curves: list[CurveItem])
        Clear the axes and redraw only the given curves.
    """

    step_requested = Signal(int)  # emitted on click: nearest time-step index

    def __init__(self, parent=None):
        super().__init__(parent)

        self._figure = Figure(tight_layout=True)
        self._canvas = FigureCanvasQTAgg(self._figure)

        # Backend nav toolbar (hidden) — used only for navigate state/history
        self._nav = NavigationToolbar2QT(self._canvas, self)
        self._nav.hide()

        self._curves: list = []
        self._dark = False
        ax0 = self._figure.add_subplot(111)
        self._axes: list = [ax0]
        self._canvas.mpl_connect("button_press_event", self._on_click)
        # Allow the canvas to receive keyboard events (used to delete
        # the currently selected data tip with the Delete/Backspace key).
        self._canvas.setFocusPolicy(Qt.StrongFocus)
        self._canvas.mpl_connect("key_press_event", self._on_key_press)

        # Zoom-inset state
        self._zoom_insets: list[ZoomInset] = []
        self._pending_selectors: list = []  # temporary selectors during "Add Zoom" mode

        # Data-cursor state
        self._cursor_artists: dict = {}     # ax -> {'vline','hline','text'}
        self._cid_cursor_motion: int | None = None
        self._cid_cursor_click: int | None = None
        self._pinned_tips: list = []        # list of dicts (see _pin_datatip)
        self._selected_pin: dict | None = None
        self._grid_on: bool = True

        # Animation time cursor (vertical line per axis, in sync with
        # the postprocessor's animation pane). The cursor is only ever
        # *shown* while the animation pane is visible.
        self._time_cursor_t: float | None = None
        self._time_cursor_visible: bool = False
        self._time_cursor_lines: dict = {}  # ax -> Line2D

        # Custom visible toolbar (built after instance vars)
        self._toolbar = self._build_toolbar()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, stretch=1)

    def _build_toolbar(self) -> QToolBar:
        tb = QToolBar()
        tb.setMovable(False)
        self._btn_home = self._make_btn(tb, "Home",    "mdi6.home",                    self._on_home)
        self._btn_back = self._make_btn(tb, "Back",    "mdi6.arrow-left",              self._on_back)
        self._btn_fwd  = self._make_btn(tb, "Forward", "mdi6.arrow-right",             self._on_fwd)
        tb.addSeparator()
        self._btn_pan      = self._make_btn(tb, "Pan",      "mdi6.hand-back-right-outline", self._on_pan,      checkable=True)
        self._btn_zoom     = self._make_btn(tb, "Zoom",     "mdi6.magnify",                 self._on_zoom,     checkable=True)
        self._btn_add_zoom = self._make_btn(tb, "Add Zoom", "mdi6.selection",               self._on_add_zoom, checkable=True)
        tb.addSeparator()
        self._btn_subplots = self._make_btn(tb, "Configure subplots",                 "mdi6.tune-vertical",          self._on_subplots)
        self._btn_customize = self._make_btn(tb, "Edit axis, curve and image parameters", "mdi6.format-list-bulleted",  self._on_customize)
        self._btn_grid     = self._make_btn(tb, "Toggle grid",                          "mdi6.grid",                   self._on_grid,     checkable=True)
        self._btn_grid.setChecked(True)
        tb.addSeparator()
        self._btn_cursor   = self._make_btn(tb, "Data cursor: click a curve to pin a static tip. Click a pin to select it, then press Delete to remove. Right-click a pin to remove directly.", "mdi6.crosshairs-gps", self._on_cursor, checkable=True)
        self._btn_clear_tips = self._make_btn(tb, "Clear all data tips",                "mdi6.eraser",                 self._on_clear_tips)
        return tb

    @staticmethod
    def _make_btn(tb: QToolBar, text: str, icon_name: str, slot,
                  *, checkable: bool = False) -> QToolButton:
        btn = QToolButton()
        btn.setIcon(_icons.icon(icon_name))
        btn.setToolTip(text)
        btn.setCheckable(checkable)
        btn.clicked.connect(slot)
        tb.addWidget(btn)
        return btn

    def set_icon_theme(self, dark: bool) -> None:
        """Re-apply icon colours after a theme toggle."""
        self._btn_home.setIcon(    _icons.icon("mdi6.home"))
        self._btn_back.setIcon(    _icons.icon("mdi6.arrow-left"))
        self._btn_fwd.setIcon(     _icons.icon("mdi6.arrow-right"))
        self._btn_pan.setIcon(     _icons.icon("mdi6.hand-back-right-outline"))
        self._btn_zoom.setIcon(    _icons.icon("mdi6.magnify"))
        self._btn_add_zoom.setIcon(_icons.icon("mdi6.selection"))
        self._btn_subplots.setIcon( _icons.icon("mdi6.tune-vertical"))
        self._btn_customize.setIcon(_icons.icon("mdi6.format-list-bulleted"))
        self._btn_grid.setIcon(    _icons.icon("mdi6.grid"))
        self._btn_cursor.setIcon(  _icons.icon("mdi6.crosshairs-gps"))
        self._btn_clear_tips.setIcon(_icons.icon("mdi6.eraser"))

    # -- Toolbar actions ------------------------------------------------

    def _on_home(self):
        self._nav.home()
        self._canvas.draw_idle()

    def _on_back(self):
        self._nav.back()
        self._canvas.draw_idle()

    def _on_fwd(self):
        self._nav.forward()
        self._canvas.draw_idle()

    def _sync_nav_mode(self, target: str | None) -> None:
        """Force the matplotlib NavigationToolbar mode to ``target``.

        ``target`` is one of ``'pan'``, ``'zoom'`` or ``None`` (idle).
        ``NavigationToolbar2.pan()`` and ``zoom()`` are toggles, so we
        inspect ``nav.mode`` and only call them when needed. This avoids
        leaving matplotlib in pan/zoom mode when switching to the custom
        Add Zoom mode.
        """
        current = str(self._nav.mode)  # '', 'pan/zoom', 'zoom rect'
        if target == "pan":
            if current == "zoom rect":
                self._nav.zoom()           # turn off zoom rect
            if str(self._nav.mode) != "pan/zoom":
                self._nav.pan()            # turn on pan
        elif target == "zoom":
            if current == "pan/zoom":
                self._nav.pan()            # turn off pan
            if str(self._nav.mode) != "zoom rect":
                self._nav.zoom()           # turn on zoom rect
        else:  # idle
            if current == "pan/zoom":
                self._nav.pan()
            elif current == "zoom rect":
                self._nav.zoom()

    def _on_pan(self, checked: bool):
        if checked:
            self._btn_zoom.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._deactivate_pending_selectors()
            self._deactivate_cursor()
            self._sync_nav_mode("pan")
        else:
            self._sync_nav_mode(None)

    def _on_zoom(self, checked: bool):
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._deactivate_pending_selectors()
            self._deactivate_cursor()
            self._sync_nav_mode("zoom")
        else:
            self._sync_nav_mode(None)

    def _on_add_zoom(self, checked: bool):
        """Toggle the zoom-inset selection mode."""
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_zoom.setChecked(False)
            self._btn_cursor.setChecked(False)
            self._sync_nav_mode(None)       # ensure matplotlib pan/zoom is OFF
            self._deactivate_cursor()
            self._activate_pending_selectors()
        else:
            self._deactivate_pending_selectors()

    def _on_cursor(self, checked: bool):
        """Toggle the data-cursor mode."""
        if checked:
            self._btn_pan.setChecked(False)
            self._btn_zoom.setChecked(False)
            self._btn_add_zoom.setChecked(False)
            self._sync_nav_mode(None)
            self._deactivate_pending_selectors()
            self._activate_cursor()
        else:
            self._deactivate_cursor()

    # -- Zoom-inset helpers ------------------------------------------------

    def _activate_pending_selectors(self) -> None:
        """Create one temporary RectangleSelector per subplot axis.

        Existing zoom-inset selectors are paused so they do not steal events
        while a new region is being drawn. Each pending selector also vetoes
        clicks that fall inside any existing inset.
        """
        self._deactivate_pending_selectors()

        # Pause all existing inset selectors
        for _z in self._zoom_insets:
            _z.pause_selector()

        for idx, ax in enumerate(self._axes):
            sel = RectangleSelector(
                ax,
                lambda ec, er, i=idx: self._on_temp_select(i, ec, er),
                useblit=False,
                button=[1],
                interactive=False,
                minspanx=5, minspany=5,
                spancoords="pixels",
                props=dict(
                    edgecolor="steelblue",
                    facecolor=(0.27, 0.51, 0.71, 0.06),
                    linestyle="--", linewidth=1,
                ),
            )
            ZoomInset.install_selector_veto(sel)
            self._pending_selectors.append(sel)

    def _deactivate_pending_selectors(self) -> None:
        for sel in self._pending_selectors:
            try:
                sel.set_active(False)
                sel.clear()
            except Exception:
                pass
        self._pending_selectors.clear()
        # Resume all existing inset selectors
        for _z in self._zoom_insets:
            _z.resume_selector()

    def _on_temp_select(self, ax_idx: int, eclick, erelease) -> None:
        """Called when the user finishes drawing the initial selection region."""
        if eclick.xdata is None or erelease.xdata is None:
            return
        x0, x1 = sorted([eclick.xdata, erelease.xdata])
        y0, y1 = sorted([eclick.ydata, erelease.ydata])
        if x1 - x0 < 1e-10 or y1 - y0 < 1e-10:
            return

        self._deactivate_pending_selectors()
        self._btn_add_zoom.setChecked(False)

        ax_parent = self._axes[ax_idx]

        # Identify the curves that belong to this subplot
        groups: dict[str, list] = {}
        for c in self._curves:
            groups.setdefault(c.unit or "Value", []).append(c)
        group_curves = list(groups.values())[ax_idx] if ax_idx < len(groups) else []

        def _remove_cb(z: ZoomInset) -> None:
            if z in self._zoom_insets:
                self._zoom_insets.remove(z)

        zoom = ZoomInset(
            parent_ax=ax_parent,
            curves=group_curves,
            x0=x0, x1=x1,
            y0=y0, y1=y1,
            dark=self._dark,
            on_remove=_remove_cb,
        )
        self._zoom_insets.append(zoom)
        self._canvas.draw_idle()

    def _on_save(self):
        # Kept as a no-op for backwards compatibility (button was removed).
        self._nav.save_figure()

    def _on_subplots(self):
        """Open the matplotlib 'Configure subplots' dialog (margins, spacing)."""
        self._nav.configure_subplots()

    def _on_customize(self):
        """Open the matplotlib 'Edit axis, curve and image parameters' dialog.

        Lets the user toggle the grid, change line colour/style/width,
        rename axes/legend entries, set scales, etc.
        """
        self._nav.edit_parameters()

    def _on_grid(self, checked: bool):
        """Toggle the dashed grid on every subplot.

        Note: matplotlib's ``Axes.grid()`` enables the grid as soon as any
        line property is supplied, regardless of the boolean argument
        (UserWarning: 'First parameter to grid() is false, but line
        properties are supplied. The grid will be enabled.'). So we apply
        the line properties only when turning the grid ON, and call
        grid(False) without extras when turning it OFF.
        """
        self._grid_on = checked
        for ax in self._axes:
            if checked:
                ax.grid(True, linestyle="--", alpha=0.5)
            else:
                ax.grid(False)
        self._canvas.draw_idle()

    def _on_clear_tips(self):
        """Remove all pinned data tips."""
        self._clear_pinned_tips()
        self._canvas.draw_idle()

    # -- Data cursor ---------------------------------------------------

    def _activate_cursor(self) -> None:
        """Install per-axes vertical cursor + multi-curve readout artists.

        For each axes we create:
          * a vertical guide line that snaps to the mouse x;
          * one annotation per curve, colored like the curve, listing the
            (x, y) value at the cursor x;
          * a header annotation showing the cursor x.
        All annotations are anchored INSIDE the axes (blended transform:
        x in data, y in axes fraction) so they never push the figure
        layout outwards.
        """
        from matplotlib.transforms import blended_transform_factory

        self._deactivate_cursor()
        fg = CANVAS_FG_DARK if self._dark else CANVAS_FG_LIGHT
        bg = CANVAS_BG_DARK if self._dark else CANVAS_BG_LIGHT

        # Group curves once, indexed by axes
        groups: dict[str, list] = {}
        for c in self._curves:
            groups.setdefault(c.unit or "Value", []).append(c)
        groups_list = list(groups.values())

        for ax_idx, ax in enumerate(self._axes):
            curves = groups_list[ax_idx] if ax_idx < len(groups_list) else []

            vline = ax.axvline(
                x=0, color=fg, linestyle=":", linewidth=1.2,
                visible=False, zorder=11,
            )

            blended = blended_transform_factory(ax.transData, ax.transAxes)
            lines: list = []

            # Header: shows the cursor x
            header = ax.annotate(
                "", xy=(0, 0.98), xycoords=blended,
                xytext=(10, -10), textcoords="offset points",
                ha="left", va="top",
                fontsize=9, color=fg, visible=False, zorder=12,
                annotation_clip=True,
                bbox=dict(boxstyle="round,pad=0.30",
                          facecolor=bg, edgecolor=fg,
                          linewidth=0.7, alpha=0.92),
            )
            lines.append(header)

            # One colored line per curve
            for i, c in enumerate(curves):
                ann = ax.annotate(
                    "", xy=(0, 0.98), xycoords=blended,
                    xytext=(10, -10 - 16 * (i + 1)), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=9, color=c.color, visible=False, zorder=12,
                    annotation_clip=True,
                    bbox=dict(boxstyle="round,pad=0.30",
                              facecolor=bg, edgecolor=c.color,
                              linewidth=0.7, alpha=0.92),
                )
                lines.append(ann)

            self._cursor_artists[ax] = {
                "vline": vline,
                "lines": lines,      # [header, curve_0, curve_1, ...]
                "curves": curves,
            }

        self._cid_cursor_motion = self._canvas.mpl_connect(
            "motion_notify_event", self._on_cursor_motion
        )
        self._cid_cursor_click = self._canvas.mpl_connect(
            "button_press_event", self._on_cursor_click
        )
        self._canvas.draw_idle()

    def _deactivate_cursor(self) -> None:
        """Remove cursor artists and disconnect the motion handler."""
        if self._cid_cursor_motion is not None:
            try:
                self._canvas.mpl_disconnect(self._cid_cursor_motion)
            except Exception:
                pass
            self._cid_cursor_motion = None
        if self._cid_cursor_click is not None:
            try:
                self._canvas.mpl_disconnect(self._cid_cursor_click)
            except Exception:
                pass
            self._cid_cursor_click = None
        for entry in self._cursor_artists.values():
            try:
                entry["vline"].remove()
            except Exception:
                pass
            for a in entry.get("lines", []):
                try:
                    a.remove()
                except Exception:
                    pass
        self._cursor_artists.clear()
        self._canvas.draw_idle()

    def _hide_all_cursors(self) -> None:
        for entry in self._cursor_artists.values():
            entry["vline"].set_visible(False)
            for a in entry.get("lines", []):
                a.set_visible(False)

    def _on_cursor_motion(self, event) -> None:
        ax = event.inaxes
        if ax is None or ax not in self._cursor_artists or event.xdata is None:
            self._hide_all_cursors()
            self._canvas.draw_idle()
            return

        entry = self._cursor_artists[ax]
        curves = entry["curves"]
        if not curves:
            return

        x_target = event.xdata
        # Sample each curve at the column closest to x_target (shared T grid is
        # not assumed: each curve is sampled independently).
        samples: list[tuple[float, float]] = []
        for c in curves:
            T = c.T
            i = int(np.searchsorted(T, x_target))
            if i >= len(T):
                i = len(T) - 1
            elif i > 0 and abs(T[i - 1] - x_target) < abs(T[i] - x_target):
                i -= 1
            samples.append((float(T[i]), float(c.data[i])))

        # Use the first curve's x as the anchor (typically all share the time grid)
        anchor_x = samples[0][0]

        # Hide cursors on the other axes
        for other_ax, other in self._cursor_artists.items():
            if other_ax is ax:
                continue
            other["vline"].set_visible(False)
            for a in other.get("lines", []):
                a.set_visible(False)

        # Update the vline
        entry["vline"].set_xdata([anchor_x])
        entry["vline"].set_visible(True)

        # Decide left/right placement based on cursor pixel position to keep
        # the readout inside the axes.
        ax_bbox = ax.get_window_extent()
        place_left = (event.x - ax_bbox.x0) > (ax_bbox.width * 0.6)
        if place_left:
            ha = "right"
            ox = -10
        else:
            ha = "left"
            ox = 10

        # Header line (x value)
        header = entry["lines"][0]
        header.xy = (anchor_x, 0.98)
        header.set_text(f"x = {anchor_x:.4g}")
        header.set_ha(ha)
        header.xyann = (ox, -10)
        header.set_visible(True)

        # Curve lines (●  label  y = value)
        for i, (c, (cx, cy)) in enumerate(zip(curves, samples)):
            ann = entry["lines"][i + 1]
            ann.xy = (anchor_x, 0.98)
            ann.set_text(f"\u25CF  {c.label}   y = {cy:.4g}")
            ann.set_ha(ha)
            ann.xyann = (ox, -10 - 16 * (i + 1))
            ann.set_visible(True)

        self._canvas.draw_idle()

    # -- Data tip pinning ----------------------------------------------

    def _hit_pin(self, event) -> dict | None:
        """Return the pin whose annotation bbox contains the event, if any."""
        if event.x is None or event.y is None:
            return None
        for pin in self._pinned_tips:
            if pin["ax"] is not event.inaxes:
                continue
            for ann in pin["anns"]:
                bbox = ann.get_window_extent()
                if bbox.contains(event.x, event.y):
                    return pin
        return None

    def _on_cursor_click(self, event) -> None:
        """Handle clicks while Cursor mode is active.

        * Left click on an existing pin -> select it (highlighted).
        * Left click elsewhere on a tracked axes -> pin a new static tip.
        * Right click on/near a pin -> remove it.
        """
        ax = event.inaxes
        if ax is None or ax not in self._cursor_artists or event.xdata is None:
            self._select_pin(None)
            return

        if event.button == 3:
            self._remove_pin_near(event)
            return
        if event.button != 1:
            return

        # Click on an existing pin -> select it instead of creating a new one
        hit = self._hit_pin(event)
        if hit is not None:
            self._select_pin(hit)
            self._canvas.setFocus()        # ensure key events reach us
            return

        entry = self._cursor_artists[ax]
        curves = entry["curves"]
        if not curves:
            return

        x_target = event.xdata
        samples: list[tuple[float, float]] = []
        for c in curves:
            T = c.T
            i = int(np.searchsorted(T, x_target))
            if i >= len(T):
                i = len(T) - 1
            elif i > 0 and abs(T[i - 1] - x_target) < abs(T[i] - x_target):
                i -= 1
            samples.append((float(T[i]), float(c.data[i])))

        anchor_x = samples[0][0]
        self._pin_datatip(ax, anchor_x, curves, samples, event)
        self._select_pin(None)         # new pin starts unselected
        self._canvas.setFocus()
        self._relayout_pins(ax)
        self._canvas.draw_idle()

    def _pin_datatip(self, ax, anchor_x: float, curves: list,
                     samples: list, event) -> None:
        """Create a static data tip (vline + colored readout) at ``anchor_x``."""
        from matplotlib.transforms import blended_transform_factory
        fg = CANVAS_FG_DARK if self._dark else CANVAS_FG_LIGHT
        bg = CANVAS_BG_DARK if self._dark else CANVAS_BG_LIGHT
        blended = blended_transform_factory(ax.transData, ax.transAxes)

        # Pin vline
        vline = ax.axvline(
            x=anchor_x, color=fg, linestyle="-",
            linewidth=0.9, alpha=0.55, zorder=10,
        )

        ax_bbox = ax.get_window_extent()
        place_left = (event.x - ax_bbox.x0) > (ax_bbox.width * 0.6)
        ha, ox = ("right", -10) if place_left else ("left", 10)

        # Anchor near the bottom of the axes (so it does not overlap the live
        # readout pinned at the top).
        y_anchor = 0.05
        anns: list = []

        header = ax.annotate(
            f"x = {anchor_x:.4g}",
            xy=(anchor_x, y_anchor), xycoords=blended,
            xytext=(ox, 10 + 16 * len(curves)), textcoords="offset points",
            ha=ha, va="bottom",
            fontsize=9, color=fg, zorder=13, annotation_clip=True,
            bbox=dict(boxstyle="round,pad=0.30",
                      facecolor=bg, edgecolor=fg,
                      linewidth=0.9, alpha=0.95),
        )
        anns.append(header)

        for i, (c, (cx, cy)) in enumerate(zip(curves, samples)):
            ann = ax.annotate(
                f"\u25CF  {c.label}   y = {cy:.4g}",
                xy=(anchor_x, y_anchor), xycoords=blended,
                xytext=(ox, 10 + 16 * (len(curves) - 1 - i)),
                textcoords="offset points",
                ha=ha, va="bottom",
                fontsize=9, color=c.color, zorder=13, annotation_clip=True,
                bbox=dict(boxstyle="round,pad=0.30",
                          facecolor=bg, edgecolor=c.color,
                          linewidth=0.9, alpha=0.95),
            )
            anns.append(ann)

        self._pinned_tips.append({
            "ax": ax,
            "x": anchor_x,
            "vline": vline,
            "anns": anns,
            "y_anchor": y_anchor,
            "base_oy": 10,            # base vertical offset (px) of header
            "line_h": 16,             # vertical step (px) between rows
            "n_rows": len(curves) + 1,
            "slot": 0,                # vertical slot index, set by _relayout_pins
            "selected": False,
        })

    def _remove_pin_near(self, event) -> None:
        """Remove the pin whose annotation bbox contains the click, if any."""
        if event.x is None or event.y is None:
            return
        for pin in list(self._pinned_tips):
            if pin["ax"] is not event.inaxes:
                continue
            for ann in pin["anns"]:
                bbox = ann.get_window_extent()
                if bbox.contains(event.x, event.y):
                    self._remove_pin(pin)
                    self._canvas.draw_idle()
                    return
        # Fallback: snap to the closest pin within ~12 px of the vline
        best = None
        best_d = 12.0
        for pin in self._pinned_tips:
            if pin["ax"] is not event.inaxes:
                continue
            try:
                px, _ = pin["ax"].transData.transform((pin["x"], 0))
            except Exception:
                continue
            d = abs(px - event.x)
            if d < best_d:
                best_d = d
                best = pin
        if best is not None:
            self._remove_pin(best)
            self._canvas.draw_idle()

    def _remove_pin(self, pin: dict) -> None:
        if self._selected_pin is pin:
            self._selected_pin = None
        ax = pin["ax"]
        try:
            pin["vline"].remove()
        except Exception:
            pass
        for ann in pin["anns"]:
            try:
                ann.remove()
            except Exception:
                pass
        if pin in self._pinned_tips:
            self._pinned_tips.remove(pin)
        self._relayout_pins(ax)

    def _clear_pinned_tips(self) -> None:
        for pin in list(self._pinned_tips):
            self._remove_pin(pin)

    # -- Selection / keyboard -----------------------------------------

    def _select_pin(self, pin: dict | None) -> None:
        """Highlight ``pin`` (None to clear selection)."""
        if self._selected_pin is pin:
            return
        # De-highlight previously selected
        prev = self._selected_pin
        if prev is not None and prev in self._pinned_tips:
            for ann in prev["anns"]:
                box = ann.get_bbox_patch()
                if box is not None:
                    box.set_linewidth(0.9)
            prev["vline"].set_alpha(0.55)
            prev["vline"].set_linewidth(0.9)
            prev["selected"] = False
        # Highlight the new selection
        if pin is not None:
            for ann in pin["anns"]:
                box = ann.get_bbox_patch()
                if box is not None:
                    box.set_linewidth(2.0)
            pin["vline"].set_alpha(0.9)
            pin["vline"].set_linewidth(1.4)
            pin["selected"] = True
        self._selected_pin = pin
        self._canvas.draw_idle()

    def _on_key_press(self, event) -> None:
        """Delete / Backspace removes the currently selected pin."""
        if event.key in ("delete", "backspace"):
            if self._selected_pin is not None:
                self._remove_pin(self._selected_pin)
                self._selected_pin = None
                self._canvas.draw_idle()
        elif event.key == "escape":
            self._select_pin(None)

    # -- Pin layout (overlap avoidance) -------------------------------

    def _relayout_pins(self, ax) -> None:
        """Re-stack pins on ``ax`` so their annotation boxes never overlap.

        Each pin's annotations are positioned at
            xytext = (ox, base_oy + slot * (n_rows * line_h + gap) + row * line_h)
        For each axes, sort pins by x and assign the lowest free slot whose
        horizontal extent does not collide with already-placed pins.
        """
        pins = [p for p in self._pinned_tips if p["ax"] is ax]
        if not pins:
            return

        # We need rendered extents to know horizontal overlap. Force a draw.
        self._canvas.draw()

        # Sort by x position (left to right, time grows that way)
        pins.sort(key=lambda p: p["x"])

        # Group pins by their (vline) x coord and decide slots
        placed: list[tuple[float, float, int]] = []  # (x_left_px, x_right_px, slot)

        for pin in pins:
            # Compute combined horizontal extent of this pin's annotations
            xs0, xs1 = float("inf"), float("-inf")
            for ann in pin["anns"]:
                bb = ann.get_window_extent()
                xs0 = min(xs0, bb.x0)
                xs1 = max(xs1, bb.x1)

            # Pick the lowest slot index that does not collide horizontally
            slot = 0
            while True:
                collides = any(
                    not (xs1 < other_x0 or xs0 > other_x1) and slot == s
                    for other_x0, other_x1, s in placed
                )
                if not collides:
                    break
                slot += 1
            placed.append((xs0, xs1, slot))

            if pin["slot"] != slot:
                pin["slot"] = slot
                self._apply_pin_slot(pin)

    def _apply_pin_slot(self, pin: dict) -> None:
        """Update annotation xytext offsets to match pin['slot']."""
        base = pin["base_oy"]
        line_h = pin["line_h"]
        n_rows = pin["n_rows"]
        gap = 6
        slot_offset = pin["slot"] * (n_rows * line_h + gap)
        # anns are [header, curve_0, curve_1, ...]
        # header sits at top of the stack -> highest offset
        # curves go below in order
        anns = pin["anns"]
        # Header offset = base + slot_offset + (n_rows-1)*line_h
        # Curve i offset = base + slot_offset + (n_rows-2-i)*line_h
        for i, ann in enumerate(anns):
            # current ox is preserved; rebuild oy
            ox = ann.xyann[0]
            if i == 0:
                oy = base + slot_offset + (n_rows - 1) * line_h
            else:
                oy = base + slot_offset + (n_rows - 1 - i) * line_h
            ann.xyann = (ox, oy)

    def update_plot(self, curves: list[CurveItem]):
        """Clear figure, create one subplot per unit group, plot curves."""
        self._curves = curves
        self._deactivate_pending_selectors()
        self._deactivate_cursor()
        self._clear_pinned_tips()
        self._btn_add_zoom.setChecked(False)
        self._btn_cursor.setChecked(False)
        # Remove each inset cleanly (disconnects events, clears selector)
        for _z in list(self._zoom_insets):
            _z.remove()
        self._zoom_insets.clear()
        # Drop stale time-cursor artists; they will be recreated below
        # against the freshly cleared axes.
        self._time_cursor_lines.clear()
        self._figure.clear()

        if not curves:
            ax0 = self._figure.add_subplot(111)
            self._axes = [ax0]
            self.set_dark(self._dark)
            self._canvas.draw_idle()
            return

        # Group curves by unit (= shared Y-label); preserves insertion order
        groups: dict[str, list] = {}
        for c in curves:
            groups.setdefault(c.unit or "Value", []).append(c)

        n = len(groups)
        self._axes = []
        for idx, (ylabel, group) in enumerate(groups.items()):
            ax = self._figure.add_subplot(n, 1, idx + 1)
            self._axes.append(ax)
            for c in group:
                ax.plot(c.T, c.data, color=c.color, label=c.label)
            ax.set_ylabel(ylabel)
            ax.legend(fontsize="small", loc="best")
            if self._grid_on:
                ax.grid(True, linestyle="--", alpha=0.5)
            else:
                ax.grid(False)
            if idx == n - 1:
                ax.set_xlabel("Time [s]")

        # Recreate the animation time-cursor (if any) against the new
        # set of axes; preserves position and visibility across redraws.
        self._refresh_time_cursor()

        self.set_dark(self._dark)
        self._canvas.draw_idle()

    # ------------------------------------------------------------------
    # Animation time-cursor
    # ------------------------------------------------------------------
    def set_time_cursor(self, t: float | None):
        """Move the synced time cursor to absolute time ``t`` (seconds).

        Pass ``None`` to clear the cursor entirely. The cursor only
        renders while it is also marked visible
        (:meth:`set_time_cursor_visible`).
        """
        self._time_cursor_t = None if t is None else float(t)
        self._refresh_time_cursor()
        self._canvas.draw_idle()

    def set_time_cursor_visible(self, on: bool):
        """Show / hide the synced time cursor without dropping its position."""
        self._time_cursor_visible = bool(on)
        for line in self._time_cursor_lines.values():
            line.set_visible(self._time_cursor_visible
                             and self._time_cursor_t is not None)
        self._canvas.draw_idle()

    def _refresh_time_cursor(self):
        """Ensure exactly one cursor Line2D per axes, at the stored time."""
        # Drop entries whose ax was removed (e.g. by ``update_plot``).
        for ax in list(self._time_cursor_lines.keys()):
            if ax not in self._axes:
                self._time_cursor_lines.pop(ax, None)

        t = self._time_cursor_t
        show = self._time_cursor_visible and (t is not None)
        fg = CANVAS_FG_DARK if self._dark else CANVAS_FG_LIGHT

        for ax in self._axes:
            line = self._time_cursor_lines.get(ax)
            if line is None:
                # Use a distinctive, theme-friendly orange so it stands
                # apart from the curve colors and from the hover/pin
                # cursors (which are foreground-coloured).
                line = ax.axvline(
                    x=0.0, color="#ff8c00", linestyle="--",
                    linewidth=1.2, alpha=0.95, zorder=15, visible=False,
                )
                self._time_cursor_lines[ax] = line
            if t is not None:
                line.set_xdata([t, t])
            line.set_visible(show)

    def set_dark(self, enabled: bool):
        """Switch the existing Figure between dark and light appearance."""
        self._dark = enabled
        bg = CANVAS_BG_DARK  if enabled else CANVAS_BG_LIGHT
        fg = CANVAS_FG_DARK  if enabled else CANVAS_FG_LIGHT
        self._figure.set_facecolor(bg)
        for ax in self._axes:
            ax.set_facecolor(bg)
            ax.tick_params(colors=fg, which="both")
            ax.xaxis.label.set_color(fg)
            ax.yaxis.label.set_color(fg)
            ax.title.set_color(fg)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg)
            # Re-skin legend so text/frame match the theme
            leg = ax.get_legend()
            if leg is not None:
                for txt in leg.get_texts():
                    txt.set_color(fg)
                frame = leg.get_frame()
                frame.set_facecolor(bg)
                frame.set_edgecolor(fg)
        # Re-skin existing cursor artists, if any
        for entry in self._cursor_artists.values():
            entry["vline"].set_color(fg)
            for i, a in enumerate(entry.get("lines", [])):
                box = a.get_bbox_patch()
                if box is not None:
                    box.set_facecolor(bg)
                if i == 0:
                    a.set_color(fg)
                    if box is not None:
                        box.set_edgecolor(fg)
                # Curve lines keep their per-curve color and edge
        # Re-skin pinned data tips: header gets fg, curve pills keep curve color
        for pin in self._pinned_tips:
            pin["vline"].set_color(fg)
            for i, ann in enumerate(pin["anns"]):
                box = ann.get_bbox_patch()
                if box is not None:
                    box.set_facecolor(bg)
                if i == 0:
                    ann.set_color(fg)
                    if box is not None:
                        box.set_edgecolor(fg)
        for _z in self._zoom_insets:
            _z.apply_theme(enabled)
        self._canvas.draw_idle()

    def _on_click(self, event):
        """Convert a matplotlib click to the nearest time-step index and emit it."""
        if (self._btn_pan.isChecked() or self._btn_zoom.isChecked()
                or self._btn_add_zoom.isChecked() or self._btn_cursor.isChecked()):
            return
        if event.button != 1:
            return
        if not self._curves or event.inaxes not in self._axes or event.xdata is None:
            return
        t_click = event.xdata
        T = self._curves[0].T
        step = int(np.argmin(np.abs(T - t_click)))
        self.step_requested.emit(step)
