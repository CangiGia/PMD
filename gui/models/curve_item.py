"""CurveItem — lightweight data container for a single plottable curve."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from PMD.src.constraints import RevJoint, TranJoint, RevRevJoint
from PMD.src.units import (
    UnitSystem,
    conversion_factor,
    ylabel_for,
    reaction_ylabel_for,
)

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


# Mapping from (category, component) to physical dimension used for
# unit conversion.  "velocity" and "acceleration" scale like length/time^n;
# angular rates scale like angle/time^n.
_DIMENSION: dict[tuple[str, str], str] = {
    ("positions",     "x"):      "length",
    ("positions",     "y"):      "length",
    ("positions",     "phi"):    "angle",
    ("velocities",    "dx"):     "velocity",
    ("velocities",    "dy"):     "velocity",
    ("velocities",    "dphi"):   "angle",
    ("accelerations", "ddx"):    "acceleration",
    ("accelerations", "ddy"):    "acceleration",
    ("accelerations", "ddphi"):  "angle",
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
    unit : str
        Y-axis group label (LaTeX-ready string).  Curves that share the same
        ``unit`` string are drawn in the same subplot.
    """

    label: str
    T: NDArray
    data: NDArray
    color: str = field(default_factory=_next_color)
    visible: bool = True
    unit: str = ""


def build_curves(
    category: str,
    component: str,
    selection: list[dict],
    display_units: UnitSystem | None = None,
) -> list[CurveItem]:
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
    display_units : UnitSystem, optional
        Unit system to use when displaying data.  If *None*, the model's own
        unit system is used (i.e. no conversion).

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
        session = desc["session"]
        T = session.T

        # Retrieve the unit system the model data is stored in
        model_us: UnitSystem = getattr(session, "units", UnitSystem())
        disp_us: UnitSystem = display_units if display_units is not None else model_us

        if multi:
            lbl = f"{session.name} / {lbl}"

        if kind == "body" and category in ("positions", "velocities", "accelerations"):
            rc = obj._result_container
            if rc is None:
                continue
            raw_data = rc[category][component]
            dim = _DIMENSION.get((category, component), "length")
            factor = conversion_factor(model_us, disp_us, dim)
            unit_label = ylabel_for(category, component, disp_us)
            curves.append(CurveItem(
                label=f"{lbl} / {component}",
                T=T,
                data=raw_data * factor,
                unit=unit_label,
            ))

        elif kind == "joint" and category == "reactions":
            rc = obj._result_container
            if rc is None:
                continue
            col_idx = int(component)
            reactions = rc["reactions"]
            if col_idx >= reactions.shape[1]:
                continue
            raw_data = reactions[:, col_idx]
            rxn_lbl = reaction_labels(obj)[col_idx]
            unit_label, dim = reaction_ylabel_for(rxn_lbl, disp_us)
            factor = conversion_factor(model_us, disp_us, dim)
            curves.append(CurveItem(
                label=f"{lbl} / {rxn_lbl}",
                T=T,
                data=raw_data * factor,
                unit=unit_label,
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
