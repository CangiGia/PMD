"""UnitSystem — declaration of physical units used in a PMD model.

The solver always works with the raw numbers supplied by the user.
``UnitSystem`` is a pure metadata object that records which physical
units those numbers represent so the GUI can apply the correct
conversion factors when displaying results.

Usage
-----
>>> from PMD.src.units import UnitSystem
>>> us = UnitSystem(length="mm", force="N")   # model built in mm / N
>>> us.factor("length")                        # 1e-3  (mm → m)
>>> us.factor("moment")                        # 1e-3  (N·mm → N·m)
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Registry: every supported unit → its SI equivalent value
# ---------------------------------------------------------------------------

#: Conversion factors to SI base units.
_TO_SI: dict[str, float] = {
    # length → metre
    "m":   1.0,
    "mm":  1e-3,
    "cm":  1e-2,
    "in":  0.0254,
    "ft":  0.3048,
    # force → newton
    "N":   1.0,
    "kN":  1e3,
    "lbf": 4.44822,
    # angle → radian  (input only; display can be "deg")
    "rad": 1.0,
    "deg": 3.141592653589793 / 180.0,
}

#: Human-readable LaTeX symbol for each unit, used in plot axis labels.
_LATEX: dict[str, str] = {
    "m":   r"\mathrm{m}",
    "mm":  r"\mathrm{mm}",
    "cm":  r"\mathrm{cm}",
    "in":  r"\mathrm{in}",
    "ft":  r"\mathrm{ft}",
    "N":   r"\mathrm{N}",
    "kN":  r"\mathrm{kN}",
    "lbf": r"\mathrm{lbf}",
    "rad": r"\mathrm{rad}",
    "deg": r"°",
}

# Derived moment symbols: force·length, built on demand.
def _moment_latex(force: str, length: str) -> str:
    return rf"{_LATEX[force]}{{\cdot}}{_LATEX[length]}"


# ---------------------------------------------------------------------------
# UnitSystem dataclass
# ---------------------------------------------------------------------------

@dataclass
class UnitSystem:
    """Declares the physical units used when defining a PMD model.

    Parameters
    ----------
    length : str
        Unit for positions / displacements.
        Choices: ``"m"`` (SI default), ``"mm"``, ``"cm"``, ``"in"``, ``"ft"``.
    force : str
        Unit for forces.
        Choices: ``"N"`` (SI default), ``"kN"``, ``"lbf"``.
    angle : str
        Unit for angles.
        Choices: ``"rad"`` (SI default), ``"deg"``.
        Note: the solver always uses radians internally; this field is
        informational and drives the display angle conversion.
    """

    length: str = "m"
    force:  str = "N"
    angle:  str = "rad"

    def __post_init__(self):
        valid = {
            "length": {"m", "mm", "cm", "in", "ft"},
            "force":  {"N", "kN", "lbf"},
            "angle":  {"rad", "deg"},
        }
        for dim, choices in valid.items():
            val = getattr(self, dim)
            if val not in choices:
                raise ValueError(
                    f"UnitSystem: invalid {dim} unit '{val}'. "
                    f"Valid choices: {sorted(choices)}"
                )

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def moment_unit(self) -> str:
        """String representation of the moment unit (force × length)."""
        return f"{self.force}·{self.length}"

    # ------------------------------------------------------------------
    # Conversion helpers
    # ------------------------------------------------------------------

    def to_si(self, dimension: str) -> float:
        """Return the factor to convert *this* system's unit to SI.

        Parameters
        ----------
        dimension : ``"length"`` | ``"force"`` | ``"angle"`` | ``"moment"``
            The physical dimension to convert.

        Returns
        -------
        float
            Multiplicative factor such that ``value_SI = value * factor``.
        """
        if dimension == "moment":
            return _TO_SI[self.force] * _TO_SI[self.length]
        unit = getattr(self, dimension)
        return _TO_SI[unit]

    def latex_unit(self, dimension: str) -> str:
        """Return the LaTeX string for the unit of *dimension*."""
        if dimension == "moment":
            return _moment_latex(self.force, self.length)
        unit = getattr(self, dimension)
        return _LATEX.get(unit, unit)

    def __repr__(self) -> str:
        return (
            f"UnitSystem(length='{self.length}', force='{self.force}', "
            f"angle='{self.angle}')"
        )


# ---------------------------------------------------------------------------
# Module-level helpers used by build_curves
# ---------------------------------------------------------------------------

def conversion_factor(model_us: "UnitSystem", display_us: "UnitSystem",
                      dimension: str) -> float:
    """Scalar to multiply raw model data to obtain display values.

    ``display_value = raw_value * conversion_factor(model, display, dim)``

    Parameters
    ----------
    model_us : UnitSystem
        Units in which the model data are stored.
    display_us : UnitSystem
        Units in which the GUI should display the data.
    dimension : ``"length"`` | ``"velocity"`` | ``"acceleration"`` \\
                | ``"force"`` | ``"moment"`` | ``"angle"``
        The physical dimension to convert.

    Returns
    -------
    float
        Multiplicative factor to convert model values to display values.
    """
    if dimension == "velocity":
        # [length_model/s]  →  [length_display/s]
        return model_us.to_si("length") / display_us.to_si("length")
    if dimension == "acceleration":
        return model_us.to_si("length") / display_us.to_si("length")
    if dimension in ("length", "force", "angle", "moment"):
        return model_us.to_si(dimension) / display_us.to_si(dimension)
    return 1.0


def ylabel_for(category: str, component: str, display_us: "UnitSystem") -> str:
    """Return a LaTeX-ready Y-axis label string for the given category/component.

    Examples
    --------
    >>> ylabel_for("positions", "x", UnitSystem(length="mm"))
    'Position $[\\\\mathrm{mm}]$'
    """
    ul = display_us.latex_unit("length")
    uf = display_us.latex_unit("force")
    ua = display_us.latex_unit("angle")
    uv = rf"{display_us.latex_unit('length')}/\mathrm{{s}}"
    uacc = rf"{display_us.latex_unit('length')}/\mathrm{{s}}^2"
    uav = rf"{display_us.latex_unit('angle')}/\mathrm{{s}}"
    uaacc = rf"{display_us.latex_unit('angle')}/\mathrm{{s}}^2"

    _MAP = {
        ("positions",     "x"):     (f"Position $[{ul}]$",            "length",   "Position"),
        ("positions",     "y"):     (f"Position $[{ul}]$",            "length",   "Position"),
        ("positions",     "phi"):   (f"Orientation $[{ua}]$",         "angle",    "Orientation"),
        ("velocities",    "dx"):    (f"Velocity $[{uv}]$",            "velocity", "Velocity"),
        ("velocities",    "dy"):    (f"Velocity $[{uv}]$",            "velocity", "Velocity"),
        ("velocities",    "dphi"):  (f"Angular velocity $[{uav}]$",   "angle",    "AngVel"),
        ("accelerations", "ddx"):   (f"Acceleration $[{uacc}]$",      "acceleration", "Acceleration"),
        ("accelerations", "ddy"):   (f"Acceleration $[{uacc}]$",      "acceleration", "Acceleration"),
        ("accelerations", "ddphi"): (f"Angular acc. $[{uaacc}]$",     "angle",    "AngAcc"),
    }
    entry = _MAP.get((category, component))
    if entry is None:
        return ""
    return entry[0]


def reaction_ylabel_for(rxn_label: str, display_us: "UnitSystem") -> tuple[str, str]:
    """Return (y_axis_label, dimension) for a given reaction physical label.

    The *dimension* string is used by ``conversion_factor``.
    """
    uf = display_us.latex_unit("force")
    um = display_us.latex_unit("moment")
    _force_labels = {"Fx", "Fy", "F_perp", "F_slide", "F_link"}
    _moment_labels = {"Mz", "M"}
    if rxn_label in _force_labels:
        return (f"Reaction $[{uf}]$", "force")
    if rxn_label in _moment_labels:
        return (f"Reaction $[{um}]$", "moment")
    return ("Reaction", "force")
