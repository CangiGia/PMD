"""CurveItem — lightweight data container for a single plottable curve."""

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

# 10-colour palette (tab10-inspired hex values)
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
_color_index = 0


def _next_color() -> str:
    global _color_index
    color = _PALETTE[_color_index % len(_PALETTE)]
    _color_index += 1
    return color


@dataclass
class CurveItem:
    """A single time-series curve ready for plotting.

    Attributes
    ----------
    label : str
        Human-readable label (e.g. ``"Body_1 / x"``).
    T : NDArray
        Time vector, shape (nSteps,).
    data : NDArray
        Values vector, shape (nSteps,).
    color : str
        CSS / hex colour string.
    visible : bool
        Whether the curve should be drawn.
    """

    label: str
    T: NDArray
    data: NDArray
    color: str = field(default_factory=_next_color)
    visible: bool = True


def build_curves(category: str, component: str,
                 selection: list[dict]) -> list[CurveItem]:
    """Build CurveItem instances from a FilterPanel request.

    Parameters
    ----------
    category : str
        ``"positions"`` | ``"velocities"`` | ``"accelerations"`` | ``"reactions"``
    component : str
        Component key (``"x"`` … ``"ddphi"`` for bodies, ``"0"`` … for reactions).
    selection : list[dict]
        Descriptor dicts from SimulationPanel (keys: kind, index, label,
        object, session).

    Returns
    -------
    list[CurveItem]
        One curve per compatible item in *selection*.
    """
    curves: list[CurveItem] = []

    # Detect multi-session to prefix labels
    sessions = {id(d["session"]) for d in selection}
    multi = len(sessions) > 1

    for desc in selection:
        kind = desc["kind"]
        obj = desc["object"]
        lbl = desc["label"]
        T = desc["session"].T

        if multi:
            lbl = f"{desc['session'].name} / {lbl}"

        if kind == "body" and category in ("positions", "velocities", "accelerations"):
            rc = obj._result_container
            if rc is None:
                continue
            data = rc[category][component]
            curves.append(CurveItem(
                label=f"{lbl} / {component}",
                T=T,
                data=data,
            ))

        elif kind == "joint" and category == "reactions":
            rc = obj._result_container
            if rc is None:
                continue
            col_idx = int(component)
            reactions = rc["reactions"]
            if col_idx >= reactions.shape[1]:
                continue
            data = reactions[:, col_idx]
            curves.append(CurveItem(
                label=f"{lbl} / \u03bb_{component}",
                T=T,
                data=data,
            ))
        # else: skip (force without data, or incompatible kind/category)

    return curves
