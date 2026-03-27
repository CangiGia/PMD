"""CurveItem — lightweight data container for a single plottable curve."""

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from PMD.src.constraints import RevJoint, TranJoint, RevRevJoint

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


_YLABEL_MAP: dict[tuple[str, str], str] = {
    ("positions",     "x"):      "Position [m]",
    ("positions",     "y"):      "Position [m]",
    ("positions",     "phi"):    "Orientation [rad]",
    ("velocities",    "dx"):     "Velocity [m/s]",
    ("velocities",    "dy"):     "Velocity [m/s]",
    ("velocities",    "dphi"):   "Angular velocity [rad/s]",
    ("accelerations", "ddx"):    "Acceleration [m/s\u00b2]",
    ("accelerations", "ddy"):    "Acceleration [m/s\u00b2]",
    ("accelerations", "ddphi"):  "Angular acceleration [rad/s\u00b2]",
}

_REACTION_YLABEL = {
    "Fx":      "Reaction [N]",
    "Fy":      "Reaction [N]",
    "Mz":      "Reaction [N\u00b7m]",
    "F_perp":  "Reaction [N]",
    "M":       "Reaction [N\u00b7m]",
    "F_slide": "Reaction [N]",
    "F_link":  "Reaction [N]",
}


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
    unit: str = ""


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
                unit=_YLABEL_MAP.get((category, component), ""),
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
            rxn_lbl = reaction_labels(obj)[col_idx]
            curves.append(CurveItem(
                label=f"{lbl} / {rxn_lbl}",
                T=T,
                data=data,
                unit=_REACTION_YLABEL.get(rxn_lbl, "Reaction"),
            ))
        # else: skip (force without data, or incompatible kind/category)

    return curves


def reaction_labels(joint) -> list[str]:
    """Return a human-readable label for each reaction column of *joint*.

    Returned labels use SI notation and reflect the physical meaning of
    each Lagrange multiplier for the given joint type.
    """
    if isinstance(joint, RevJoint):
        labels = ["Fx", "Fy"]
        if getattr(joint, "fix", 0) == 1:
            labels.append("Mz")
        return labels
    if isinstance(joint, TranJoint):
        labels = ["F_perp", "M"]
        if getattr(joint, "fix", 0) == 1:
            labels.append("F_slide")
        return labels
    if isinstance(joint, RevRevJoint):
        return ["F_link"]
    # generic fallback for any other joint type
    rc = joint._result_container
    if rc is not None:
        return [f"\u03bb_{i}" for i in range(rc["reactions"].shape[1])]
    return []
