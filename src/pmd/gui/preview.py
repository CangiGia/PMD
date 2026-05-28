"""Pre-simulation model snapshot viewer.

Renders a static matplotlib figure of a :class:`PlanarMultibodyModel` at
its current state (typically right after construction, before any call to
:meth:`~pmd.core.solver.PlanarMultibodyModel.solve`).

This module is intentionally lightweight:

* No Qt import — uses whatever matplotlib backend the caller has active,
  so it works identically in scripts, Jupyter notebooks, and headless
  environments (``show=False``).
* Read-only — reads ``body.position``, ``body.orientation``, ``body.shape``,
  ``body.color``, and ``body._markers``; **no model state is mutated**.
* Optional — nothing in the solver path imports this module, so it has
  zero runtime cost when unused.

The visual style (body fills, joint glyphs, marker triads) mirrors the
post-processor :class:`~pmd.gui.postprocessor.widgets.animation_canvas.AnimationCanvas`
so both views look consistent.

Examples
--------
Basic blocking window::

    from pmd.gui import preview_model
    preview_model(model)

Hide marker triads::

    preview_model(model, show_markers=False)

Non-blocking, save to file::

    fig = preview_model(model, show=False)
    fig.savefig("initial_config.png", dpi=150)

Embed in an existing axes::

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    preview_model(model, ax=ax, show=False)
    plt.tight_layout()
    plt.show()

Author: Giacomo Cangi
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import (
    Circle as MplCircle,
)
from matplotlib.patches import (
    FancyArrowPatch,
    FancyBboxPatch,
)
from matplotlib.patches import (
    Polygon as MplPolygon,
)
from matplotlib.transforms import Affine2D

from pmd.core.constraints import RevJoint, TranJoint
from pmd.core.mechanics import rotation_matrix
from pmd.core.motion import Motion
from pmd.core.shapes import Circle, Link, Plate, Rectangle

# ── Colour constants (mirror animation_canvas) ────────────────────────────

# Same tab10 palette used by the post-processor so pre- and post-simulation
# views are visually consistent for bodies that have no explicit colour.
_BODY_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

_TRIAD_C_X = "#d65a5a"  # red  — X axis
_TRIAD_C_Y = "#3aa35a"  # green — Y axis
_TRIAD_C_DOT_FC = "#ffd400"  # yellow origin fill
_TRIAD_C_DOT_EC = "#1c2033"  # origin edge


# ── Internal geometry helpers ─────────────────────────────────────────────


def _body_state(body) -> tuple[float, float, float]:
    """Return (x, y, phi) of a body from its current position/orientation."""
    pos = np.asarray(body.position, dtype=float).ravel()
    return float(pos[0]), float(pos[1]), float(body.orientation)


def _marker_global(marker) -> np.ndarray:
    """Return the global (x, y) position of a marker."""
    body = marker.body
    if not body:  # Ground marker
        return marker.local_position.ravel()[:2].copy()
    x, y, phi = _body_state(body)
    return np.array([x, y]) + rotation_matrix(phi) @ marker.local_position.ravel()[:2]


def _marker_global_angle(marker) -> float:
    """Return the global orientation (rad) of a marker."""
    body = marker.body
    phi_body = 0.0 if not body else _body_state(body)[2]
    return phi_body + float(getattr(marker, "theta", None) or 0.0)


def _scene_span(model) -> float:
    """Return the largest x- or y-extent of body positions."""
    if not model.Bodies:
        return 1.0
    pts = np.array([_body_state(b)[:2] for b in model.Bodies])
    span_x = float(pts[:, 0].max() - pts[:, 0].min())
    span_y = float(pts[:, 1].max() - pts[:, 1].min())
    return max(span_x, span_y, 1.0)


# ── Internal drawing helpers ──────────────────────────────────────────────


def _draw_body(ax, body, col, triad_L, show_markers):
    """Draw one body patch plus optional label and marker triads."""
    x, y, phi = _body_state(body)
    shape = body.shape
    rgb = matplotlib.colors.to_rgb(col)
    face = (*rgb, 0.30)
    edge = (*rgb, 1.00)

    if isinstance(shape, Rectangle):
        patch = FancyBboxPatch(
            (-shape.width / 2, -shape.height / 2),
            shape.width,
            shape.height,
            boxstyle="round,pad=0.01",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.5,
        )
        patch.set_transform(Affine2D().rotate(phi).translate(x, y) + ax.transData)
        ax.add_patch(patch)

    elif isinstance(shape, Circle):
        ax.add_patch(
            MplCircle(
                (x, y),
                shape.radius,
                facecolor=face,
                edgecolor=edge,
                linewidth=1.5,
            )
        )

    elif isinstance(shape, Plate):
        patch = MplPolygon(
            shape.vertices.copy(),
            closed=True,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.5,
        )
        patch.set_transform(Affine2D().rotate(phi).translate(x, y) + ax.transData)
        ax.add_patch(patch)

    elif isinstance(shape, Link):
        patch = FancyBboxPatch(
            (-shape.length / 2, -shape.thickness / 2),
            shape.length,
            shape.thickness,
            boxstyle="round,pad=0.002",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.5,
        )
        patch.set_transform(Affine2D().rotate(phi).translate(x, y) + ax.transData)
        ax.add_patch(patch)

    else:
        # No recognised shape — small disc at the CoM
        ax.add_patch(
            MplCircle(
                (x, y),
                triad_L * 0.5,
                facecolor=face,
                edgecolor=edge,
                linewidth=1.5,
            )
        )

    if show_markers:
        for mk in body._markers:
            _draw_marker(ax, mk, triad_L)


def _draw_marker(ax, mk, triad_L: float) -> None:
    """Draw a single marker as a red-X / green-Y triad with a yellow dot."""
    gp = _marker_global(mk)
    a = _marker_global_angle(mk)
    head_x = (gp[0] + triad_L * np.cos(a), gp[1] + triad_L * np.sin(a))
    head_y = (gp[0] - triad_L * np.sin(a), gp[1] + triad_L * np.cos(a))
    ax.add_patch(
        FancyArrowPatch(
            (gp[0], gp[1]),
            head_x,
            arrowstyle="-|>",
            mutation_scale=8,
            color=_TRIAD_C_X,
            linewidth=1.4,
            shrinkA=0,
            shrinkB=0,
            zorder=5,
        )
    )
    ax.add_patch(
        FancyArrowPatch(
            (gp[0], gp[1]),
            head_y,
            arrowstyle="-|>",
            mutation_scale=8,
            color=_TRIAD_C_Y,
            linewidth=1.4,
            shrinkA=0,
            shrinkB=0,
            zorder=5,
        )
    )
    ax.plot(
        gp[0],
        gp[1],
        marker="o",
        linestyle="none",
        markersize=3.5,
        markerfacecolor=_TRIAD_C_DOT_FC,
        markeredgecolor=_TRIAD_C_DOT_EC,
        markeredgewidth=0.8,
        zorder=6,
    )
    mk_name = getattr(mk, "name", "") or ""
    if mk_name:
        ax.annotate(
            mk_name,
            xy=(gp[0], gp[1]),
            xytext=(-2.0, -6.0),
            textcoords="offset points",
            fontsize=9.0,
            color="#1c2033",
            ha="right",
            va="top",
            zorder=7,
            annotation_clip=False,
        )


def _draw_joint(ax, joint, r_joint, sh, rail_hw, gap_y) -> None:
    """Draw one joint glyph (RevJoint disc or TranJoint square + rails)."""
    if isinstance(joint, Motion):
        return
    mk = joint.iMarker or joint.jMarker
    if mk is None:
        return
    gp = _marker_global(mk)
    ang = _marker_global_angle(mk)

    if isinstance(joint, RevJoint):
        ax.add_patch(
            MplCircle(
                gp,
                radius=r_joint,
                facecolor="#ff2bd6",
                edgecolor="#1c2033",
                linewidth=1.0,
                zorder=7,
            )
        )

    elif isinstance(joint, TranJoint):
        patch = FancyBboxPatch(
            (-sh, -sh),
            2 * sh,
            2 * sh,
            boxstyle="square,pad=0",
            facecolor="#00bcd4",
            edgecolor="#006978",
            linewidth=1.2,
            zorder=7,
        )
        patch.set_transform(Affine2D().rotate(ang).translate(gp[0], gp[1]) + ax.transData)
        ax.add_patch(patch)
        c, s = np.cos(ang), np.sin(ang)
        ex, ey = np.array([c, s]), np.array([-s, c])
        for sign in (+1, -1):
            p1 = gp - rail_hw * ex + sign * gap_y * ey
            p2 = gp + rail_hw * ex + sign * gap_y * ey
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                color="#006978",
                linewidth=1.4,
                solid_capstyle="butt",
                zorder=7,
            )

    else:
        ax.add_patch(
            MplCircle(
                gp,
                radius=r_joint * 0.5,
                facecolor="#888888",
                edgecolor="#333333",
                linewidth=0.8,
                zorder=7,
            )
        )


# ── Public API ────────────────────────────────────────────────────────────


def preview_model(
    model,
    *,
    show_markers: bool = True,
    show_joints: bool = True,
    show_body_labels: bool = True,
    show_joint_legend: bool = True,
    ax=None,
    show: bool = True,
    block: bool = True,
):
    """Render a static snapshot of the model at its current state.

    Parameters
    ----------
    model : PlanarMultibodyModel
        A constructed model.  Called before or after :meth:`solve`; the
        function reads positions/orientations as they currently are.
    show_markers : bool, optional
        Draw the red-X / green-Y triad at every marker. Default ``True``.
    show_joints : bool, optional
        Draw a magenta disc (revolute) or cyan square (translational) at
        every joint. Default ``True``.
    show_body_labels : bool, optional
        Annotate each body with its name. Default ``True``.
    show_joint_legend : bool, optional
        Add a legend describing the joint glyphs. Default ``True``.
    ax : matplotlib.axes.Axes, optional
        Target axes to draw on.  When ``None`` (default) a new figure is
        created with ``figsize=(11, 6)``.
    show : bool, optional
        Call :func:`matplotlib.pyplot.show` before returning.
        Default ``True``.
    block : bool, optional
        Forwarded to ``plt.show(block=...)``.  Set to ``False`` for a
        non-blocking display (useful in interactive sessions).
        Default ``True``.

    Returns
    -------
    matplotlib.figure.Figure
        The figure containing the snapshot.  Useful for saving::

            fig = preview_model(model, show=False)
            fig.savefig("snap.png", dpi=150)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 6))
    else:
        fig = ax.figure

    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_title(
        f"Model preview — initial state ({len(model.Bodies)} bodies, {len(model.Joints)} joints)",
        fontsize=10,
    )

    span = _scene_span(model)
    triad_L = 0.04 * span  # marker arm length in data units
    r_joint = triad_L * 0.18
    sh = triad_L * 0.18
    rail_hw = triad_L * 0.45
    gap_y = sh * 1.6

    color_idx = 0

    # ── Bodies ────────────────────────────────────────────────────
    for body in model.Bodies:
        if not body:
            continue
        if getattr(body, "color", None) is not None:
            col = body.color
        else:
            col = _BODY_COLORS[color_idx % len(_BODY_COLORS)]
            color_idx += 1

        _draw_body(ax, body, col, triad_L, show_markers)

        if show_body_labels and getattr(body, "name", None):
            x, y, _ = _body_state(body)
            ax.annotate(
                body.name,
                xy=(x, y),
                xytext=(12, 12),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=9.5,
                fontweight="bold",
                color=col,
                zorder=11,
                annotation_clip=False,
            )

    # Ground markers
    if show_markers:
        from pmd.core.model import Ground

        for mk in Ground._markers:
            _draw_marker(ax, mk, triad_L)

    # ── Joints ────────────────────────────────────────────────────
    if show_joints:
        for joint in model.Joints:
            _draw_joint(ax, joint, r_joint, sh, rail_hw, gap_y)

        if show_joint_legend:
            handles = [
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="#ff2bd6",
                    markeredgecolor="#1c2033",
                    markersize=10,
                    label="Revolute Joint",
                ),
                Line2D(
                    [],
                    [],
                    marker="s",
                    linestyle="none",
                    markerfacecolor="#00bcd4",
                    markeredgecolor="#006978",
                    markersize=10,
                    label="Translational Joint",
                ),
            ]
            ax.legend(handles=handles, title="Joint types", loc="upper right", fontsize=9)

    ax.relim()
    ax.autoscale_view()
    ax.margins(0.08)

    if show:
        plt.show(block=block)
    return fig
