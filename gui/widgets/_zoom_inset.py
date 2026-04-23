"""_zoom_inset.py — Interactive zoom inset for PlotCanvas.

Each ZoomInset owns:
  - A RectangleSelector (interactive=True, drag_from_anywhere=True) on the
    parent axes — users can move or resize the selected region with handles.
  - An inset Axes with a fully opaque background (no parent curves bleeding
    through) and visible tick labels and axis scales (publication style).
  - Four ConnectionPatch connectors, one per corner of the selection rectangle,
    that update live as the selector or the inset moves.
  - A ✕ close label (top-right of inset) — click to remove.
  - Full-body drag: clicking and dragging anywhere on the inset body (outside
    the ✕ and the resize handle) repositions the inset panel.
  - A resize handle (bottom-right corner, small coloured square) — drag to
    resize the inset panel.
"""

from __future__ import annotations

from typing import Callable

from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import RectangleSelector

from ..style import CANVAS_BG_DARK, CANVAS_BG_LIGHT, CANVAS_FG_DARK, CANVAS_FG_LIGHT

_CONN_LIGHT = "steelblue"
_CONN_DARK  = "#7aa3cc"


class ZoomInset:
    """An interactive zoom-inset overlay on a single parent Axes."""

    _W0    = 0.38   # default width  in parent-axes fraction
    _H0    = 0.32   # default height in parent-axes fraction
    _MIN_W = 0.12
    _MIN_H = 0.10

    # Class-level registry of all live insets. Used to veto mouse clicks that
    # fall inside any existing inset on RectangleSelectors created by other
    # ZoomInsets (or by the PlotCanvas "Add Zoom" mode).
    _registry: list["ZoomInset"] = []

    # ------------------------------------------------------------------
    # Selector veto helper
    # ------------------------------------------------------------------

    @classmethod
    def install_selector_veto(cls, selector) -> None:
        """Patch a selector's ignore() to reject clicks inside any live inset.

        Used both by ZoomInset itself (on its persistent selector) and by
        PlotCanvas (on the temporary "Add Zoom" selectors).
        """
        original_ignore = selector.ignore

        def _patched(event):
            if original_ignore(event):
                return True
            if event.x is None or event.y is None:
                return False
            for z in cls._registry:
                try:
                    bbox = z._inset_ax.get_window_extent()
                except Exception:
                    continue
                if bbox.contains(event.x, event.y):
                    return True
            return False

        selector.ignore = _patched

    def __init__(
        self,
        parent_ax,
        curves: list,
        x0: float, x1: float,
        y0: float, y1: float,
        dark: bool = False,
        on_remove: Callable | None = None,
    ):
        self._parent_ax = parent_ax
        self._fig       = parent_ax.get_figure()
        self._curves    = curves
        self._dark      = dark
        self._on_remove = on_remove

        # Current zoom region in parent data coordinates
        self._x0, self._x1 = x0, x1
        self._y0, self._y1 = y0, y1

        # Inset position as [l, b, w, h] in parent-axes fraction
        self._pos: list[float] = self._initial_pos(x0, x1)

        # Interaction state
        self._btn_down      = False
        self._drag_active   = False
        self._resize_active = False
        self._start_px      = (0.0, 0.0)
        self._start_pos: list[float] = list(self._pos)

        # Build all matplotlib artists
        self._build_inset()
        self._build_selector(x0, x1, y0, y1)
        self._connectors = self._make_connectors()
        self._build_overlays()

        # Register this inset (must happen before connecting events)
        ZoomInset._registry.append(self)

        # Connect canvas events
        canvas = self._fig.canvas
        self._cids = [
            canvas.mpl_connect("button_press_event",   self._on_press),
            canvas.mpl_connect("motion_notify_event",  self._on_motion),
            canvas.mpl_connect("button_release_event", self._on_release),
        ]

        self.apply_theme(dark)

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    def _initial_pos(self, x0: float, x1: float) -> list[float]:
        """Place the inset in the horizontal half opposite to the selection."""
        xlim  = self._parent_ax.get_xlim()
        xspan = (xlim[1] - xlim[0]) or 1.0
        xmid  = ((x0 + x1) / 2.0 - xlim[0]) / xspan
        left  = 0.57 if xmid < 0.5 else 0.03
        return [left, 0.55, self._W0, self._H0]

    def _build_inset(self) -> None:
        self._inset_ax = self._parent_ax.inset_axes(self._pos)
        # Render the inset above all parent axes content.
        self._inset_ax.set_zorder(5)
        # Ensure a fully opaque background so parent curves do not bleed through.
        self._inset_ax.patch.set_alpha(1.0)
        for c in self._curves:
            self._inset_ax.plot(c.T, c.data, color=c.color, linewidth=0.9)
        self._inset_ax.set_xlim(self._x0, self._x1)
        self._inset_ax.set_ylim(self._y0, self._y1)
        self._inset_ax.tick_params(labelsize=7, length=2, width=0.6)
        self._inset_ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        self._inset_ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        for spine in self._inset_ax.spines.values():
            spine.set_linewidth(0.8)

    def _build_selector(self, x0, x1, y0, y1) -> None:
        self._selector = RectangleSelector(
            self._parent_ax,
            self._on_selector_release,
            useblit=False,
            button=[1],
            interactive=True,
            drag_from_anywhere=True,
            ignore_event_outside=True,   # clicks outside the rect are ignored,
                                          # so the rectangle never disappears
            minspanx=2, minspany=2,
            spancoords="pixels",
            props=dict(
                edgecolor="steelblue",
                facecolor="none",   # no fill — dashed border is enough
                linestyle="--",
                linewidth=1.0,
            ),
            handle_props=dict(
                markersize=5,
                markerfacecolor="white",
                markeredgecolor="steelblue",
                markeredgewidth=1.2,
            ),
        )
        self._selector.extents = (x0, x1, y0, y1)

        # Veto clicks that fall inside any live inset (including this one):
        # without this, the selector — which has drag_from_anywhere=True —
        # would steal events from our drag/resize/close handlers.
        ZoomInset.install_selector_veto(self._selector)

    def _make_connectors(self) -> list:
        """Four ConnectionPatch lines — one per corner of the selection region.

        Connectors are added as children of the parent axes (low zorder) so
        the inset axes (high zorder) renders on top of them and visually
        clips the lines at the inset border.
        """
        color = _CONN_DARK if self._dark else _CONN_LIGHT
        corners = [
            ((self._x0, self._y0), (0.0, 0.0)),
            ((self._x1, self._y0), (1.0, 0.0)),
            ((self._x0, self._y1), (0.0, 1.0)),
            ((self._x1, self._y1), (1.0, 1.0)),
        ]
        result = []
        for pa_pt, ins_pt in corners:
            cp = ConnectionPatch(
                xyA=pa_pt,  coordsA="data",          axesA=self._parent_ax,
                xyB=ins_pt, coordsB="axes fraction",  axesB=self._inset_ax,
                color=color, linestyle="--", linewidth=0.8, alpha=0.7,
                zorder=2, clip_on=False,
            )
            self._parent_ax.add_artist(cp)
            result.append(cp)
        return result

    def _build_overlays(self) -> None:
        # ✕ annotation — visible hit target in _on_press, no picker needed
        self._close_btn = self._inset_ax.annotate(
            "✕",
            xy=(1.0, 1.0), xycoords="axes fraction",
            xytext=(-3, -3), textcoords="offset points",
            ha="right", va="top",
            fontsize=8, fontweight="bold",
            color="#888888", zorder=10,
        )
        # Small coloured square at bottom-right — visual cue for the resize zone
        self._resize_hint = Rectangle(
            (0.86, 0.0), 0.14, 0.10,
            transform=self._inset_ax.transAxes,
            facecolor="steelblue", edgecolor="none",
            alpha=0.35, zorder=9,
        )
        self._inset_ax.add_patch(self._resize_hint)

    # ------------------------------------------------------------------
    # Selector callbacks
    # ------------------------------------------------------------------

    def _on_selector_release(self, eclick, erelease) -> None:
        self._sync_from_selector()

    def _sync_from_selector(self) -> None:
        """Read current selector extents and refresh the inset and connectors."""
        try:
            x0, x1, y0, y1 = self._selector.extents
        except Exception:
            return
        if x1 <= x0 or y1 <= y0:
            return
        self._x0, self._x1 = x0, x1
        self._y0, self._y1 = y0, y1
        self._inset_ax.set_xlim(x0, x1)
        self._inset_ax.set_ylim(y0, y1)
        self._update_connectors()
        self._fig.canvas.draw_idle()

    def _update_connectors(self) -> None:
        """Move connector start-points to the current selector corners."""
        corners = [
            (self._x0, self._y0), (self._x1, self._y0),
            (self._x0, self._y1), (self._x1, self._y1),
        ]
        for cp, pt in zip(self._connectors, corners):
            cp.xy1 = pt  # xy1 is read at draw time by ConnectionPatch

    # ------------------------------------------------------------------
    # Mouse event handlers
    # ------------------------------------------------------------------

    def _inset_frac(self, x_px: float, y_px: float):
        """Convert display-space pixel position to inset-axes fraction.

        Returns (xf, yf) or None on error.
        """
        try:
            return self._inset_ax.transAxes.inverted().transform((x_px, y_px))
        except Exception:
            return None

    def _on_press(self, event) -> None:
        if event.button == 1:
            self._btn_down = True  # used for live selector sync in _on_motion

        if event.button != 1 or event.x is None:
            return

        frac = self._inset_frac(event.x, event.y)
        if frac is None:
            return
        xf, yf = frac

        if not (0.0 <= xf <= 1.0 and 0.0 <= yf <= 1.0):
            return  # click is outside the inset panel

        # ✕ close button — top-right corner
        if xf > 0.80 and yf > 0.84:
            self.remove()
            return

        # Resize handle — bottom-right corner (matches the blue square)
        if xf > 0.85 and yf < 0.12:
            self._resize_active = True
            self._start_px  = (event.x, event.y)
            self._start_pos = list(self._pos)
            return

        # Any other point inside the inset body → drag the inset panel
        self._drag_active = True
        self._start_px  = (event.x, event.y)
        self._start_pos = list(self._pos)

    def _on_motion(self, event) -> None:
        if event.x is None or event.y is None:
            return

        ax_bbox = self._parent_ax.get_window_extent()
        if ax_bbox.width == 0 or ax_bbox.height == 0:
            return

        # --- Dragging the inset panel ---
        if self._drag_active:
            dx = (event.x - self._start_px[0]) / ax_bbox.width
            dy = (event.y - self._start_px[1]) / ax_bbox.height
            l0, b0, w0, h0 = self._start_pos
            new_l = max(0.0, min(1.0 - w0, l0 + dx))
            new_b = max(0.0, min(1.0 - h0, b0 + dy))
            self._apply_pos([new_l, new_b, w0, h0])
            return

        # --- Resizing the inset panel ---
        if self._resize_active:
            dx = (event.x - self._start_px[0]) / ax_bbox.width
            dy = (event.y - self._start_px[1]) / ax_bbox.height
            l0, b0, w0, h0 = self._start_pos
            new_w = max(self._MIN_W, min(1.0 - l0, w0 + dx))
            new_h = max(self._MIN_H, min(1.0 - b0, h0 - dy))  # drag down → taller
            self._apply_pos([l0, b0, new_w, new_h])
            return

        # --- Live sync while the selector is being moved on the parent axes ---
        if self._btn_down and event.inaxes is self._parent_ax:
            self._sync_from_selector()

    def _on_release(self, event) -> None:
        if event.button == 1:
            self._btn_down      = False
            self._drag_active   = False
            self._resize_active = False

    def _apply_pos(self, pos: list[float]) -> None:
        """Apply a new [l, b, w, h] position to the inset, updating the locator."""
        self._pos = pos
        loc = self._inset_ax.get_axes_locator()
        if loc is not None:
            # _TransformedBoundsLocator stores bounds as loc._bounds and reads
            # them at draw time, so a direct assignment is sufficient.
            loc._bounds = pos
        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Theme
    # ------------------------------------------------------------------

    def apply_theme(self, dark: bool) -> None:
        self._dark = dark
        bg    = CANVAS_BG_DARK  if dark else CANVAS_BG_LIGHT
        fg    = CANVAS_FG_DARK  if dark else CANVAS_FG_LIGHT
        c_col = _CONN_DARK if dark else _CONN_LIGHT

        self._inset_ax.set_facecolor(bg)
        self._inset_ax.patch.set_alpha(1.0)  # re-assert opaque after theme change
        self._inset_ax.tick_params(colors=fg, which="both")
        for spine in self._inset_ax.spines.values():
            spine.set_edgecolor(fg)
        self._close_btn.set_color("#aaaaaa" if dark else "#555555")
        for cp in self._connectors:
            cp.set_color(c_col)

    # ------------------------------------------------------------------
    # Pause / resume — used by PlotCanvas during the "Add Zoom" mode
    # ------------------------------------------------------------------

    def pause_selector(self) -> None:
        """Disable the persistent selector (e.g. while a new one is being drawn)."""
        try:
            self._selector.set_active(False)
            self._selector.set_visible(False)
        except Exception:
            pass

    def resume_selector(self) -> None:
        """Re-enable the persistent selector and restore its visible rectangle."""
        try:
            self._selector.set_active(True)
            self._selector.set_visible(True)
            # extents are preserved across set_active toggles
        except Exception:
            pass
        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Removal
    # ------------------------------------------------------------------

    def remove(self) -> None:
        """Disconnect events, remove all artists, notify the parent canvas."""
        # Unregister first so other selectors stop vetoing clicks at our bbox
        if self in ZoomInset._registry:
            ZoomInset._registry.remove(self)

        canvas = self._fig.canvas
        for cid in self._cids:
            try:
                canvas.mpl_disconnect(cid)
            except Exception:
                pass

        for cp in self._connectors:
            try:
                cp.remove()
            except Exception:
                pass

        try:
            self._selector.set_active(False)
            self._selector.clear()
        except Exception:
            pass

        try:
            self._inset_ax.remove()
        except Exception:
            pass

        if self._on_remove is not None:
            try:
                self._on_remove(self)
            except Exception:
                pass

        canvas.draw_idle()
