"""Typed editors for the Inspector panel.

Each editor exposes a Qt form (label + widget pairs) bound to one
dataclass instance from :mod:`pmd.gui.preprocessor.models`. Edits are
written back to the spec live and emit :sig:`spec_changed`.
"""

from __future__ import annotations

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QWidget,
)

from ..models import (
    BodySpec,
    ForceSpec,
    JointSpec,
    MarkerSpec,
    ModelSpec,
    MATERIALS,
    compute_mass_props,
)


from PySide6.QtWidgets import QCheckBox, QFrame


# QSS applied to widgets that are locked after creation: light red fill
# + dark red text + red border so the user clearly sees the field is
# read-only and cannot be changed without breaking the model geometry.
_LOCKED_QSS = (
    "QDoubleSpinBox[locked=\"true\"], QLineEdit[locked=\"true\"], "
    "QComboBox[locked=\"true\"] {"
    "  background: rgba(214, 90, 90, 0.12);"
    "  color: #8b1c1c;"
    "  border: 1px solid rgba(184, 50, 50, 0.55);"
    "  border-radius: 3px;"
    "}"
    "QDoubleSpinBox[autoval=\"true\"], QLineEdit[autoval=\"true\"] {"
    "  background: rgba(63, 140, 255, 0.10);"
    "  color: #1d3f73;"
    "  border: 1px solid rgba(63, 140, 255, 0.55);"
    "  border-radius: 3px;"
    "}"
    # Make sure the check indicator is always visible regardless
    # of the active palette / global stylesheet.
    "QCheckBox::indicator {"
    "  width: 14px; height: 14px;"
    "  border: 1px solid #5a6f8f;"
    "  border-radius: 3px;"
    "  background: #ffffff;"
    "}"
    "QCheckBox::indicator:hover { border-color: #3f8cff; }"
    "QCheckBox::indicator:checked {"
    "  background: #3f8cff;"
    "  border: 1px solid #1a5cbf;"
    "  image: none;"
    "}"
    "QCheckBox::indicator:disabled {"
    "  background: #e6e6ea; border-color: #b8bccb;"
    "}"
)


def _set_locked(widget, on: bool) -> None:
    """Tag a widget as 'locked' (red) and disable it accordingly."""
    widget.setProperty("locked", bool(on))
    widget.setProperty("autoval", False)
    widget.setReadOnly(on) if hasattr(widget, "setReadOnly") else None
    widget.setEnabled(not on)
    widget.style().unpolish(widget)
    widget.style().polish(widget)


def _set_autoval(widget, on: bool) -> None:
    """Tag a widget as auto-computed (blue) and disable user input."""
    widget.setProperty("autoval", bool(on))
    widget.setProperty("locked", False)
    if hasattr(widget, "setReadOnly"):
        widget.setReadOnly(on)
    widget.setEnabled(not on)
    widget.style().unpolish(widget)
    widget.style().polish(widget)


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

def _spinbox(value: float, *, decimals: int = 4, step: float = 0.01,
             vmin: float = -1e6, vmax: float = 1e6,
             suffix: str = "") -> QDoubleSpinBox:
    sb = QDoubleSpinBox()
    sb.setDecimals(decimals)
    sb.setRange(vmin, vmax)
    sb.setSingleStep(step)
    sb.setValue(value)
    if suffix:
        sb.setSuffix(f"  {suffix}")
    sb.setKeyboardTracking(False)
    return sb


# ──────────────────────────────────────────────────────────────────
# Base
# ──────────────────────────────────────────────────────────────────

class EditorBase(QWidget):
    """Common skeleton: form layout + spec_changed signal."""

    spec_changed = Signal(str, str)  # (kind, spec.id)

    def __init__(self, spec, kind: str, parent=None):
        super().__init__(parent)
        self.spec = spec
        self.spec_kind = kind
        self.form = QFormLayout(self)
        self.form.setContentsMargins(0, 0, 0, 0)
        self.form.setHorizontalSpacing(8)
        self.form.setVerticalSpacing(4)
        self._build()

    def _build(self) -> None:
        raise NotImplementedError

    def _emit(self) -> None:
        self.spec_changed.emit(self.spec_kind, self.spec.id)


# ──────────────────────────────────────────────────────────────────
# Body
# ──────────────────────────────────────────────────────────────────

class BodyEditor(EditorBase):
    """Editor for :class:`BodySpec` — name / inertia / state / shape."""

    # Notifies the main window that this body's name changed so it can
    # cascade-rename its child markers ("<old>.cm" → "<new>.cm" etc.).
    body_renamed = Signal(str, str, str)   # (body_id, old_name, new_name)

    def __init__(self, spec: BodySpec, parent=None):
        super().__init__(spec, "body", parent)

    def _build(self) -> None:
        s: BodySpec = self.spec
        self._old_name = s.name

        self.name = QLineEdit(s.name)
        self.name.editingFinished.connect(self._on_name)
        self.form.addRow("Name", self.name)

        # ── Material / density block ───────────────────────────
        self.material = QComboBox()
        for mat in MATERIALS:
            self.material.addItem(mat, mat)
        if s.material in MATERIALS:
            self.material.setCurrentText(s.material)
        else:
            self.material.setCurrentText("Custom")
        self.material.currentTextChanged.connect(self._on_material)
        self.form.addRow("Material", self.material)

        self.density = _spinbox(s.density, decimals=2, step=10.0,
                                vmin=0.0, vmax=2.5e4, suffix="kg/m³")
        self.density.valueChanged.connect(self._on_density)
        self.form.addRow("Density", self.density)

        # NOTE: ``thickness_z`` is shown lower down inside the Shape
        # block — it's a shape parameter (not part of the rigid-body
        # state) used only to derive mass / inertia.
        self.thickness_z = _spinbox(s.thickness_z, decimals=4, step=0.001,
                                    vmin=1e-6, vmax=10.0, suffix="m")
        self.thickness_z.valueChanged.connect(self._on_thickness_z)

        # Mass / inertia (auto by default, blue; or override → editable)
        self.mass = _spinbox(s.mass, decimals=4, step=0.1, vmin=1e-9,
                             vmax=1e9, suffix="kg")
        self.mass.valueChanged.connect(self._on_mass)
        self.form.addRow("Mass", self.mass)

        self.inertia = _spinbox(s.inertia, decimals=6, step=0.001, vmin=1e-12,
                                vmax=1e9, suffix="kg·m²")
        self.inertia.valueChanged.connect(self._on_inertia)
        self.form.addRow("Inertia", self.inertia)

        self.override = QCheckBox("Override mass / inertia manually")
        self.override.setChecked(bool(s.mass_override))
        self.override.toggled.connect(self._on_override_toggled)
        self.form.addRow("", self.override)

        # ── Position ───────────────────────────────────────────
        self.x = _spinbox(s.position[0], suffix="m")
        self.y = _spinbox(s.position[1], suffix="m")
        self.x.valueChanged.connect(self._on_pos)
        self.y.valueChanged.connect(self._on_pos)
        self.form.addRow("x", self.x)
        self.form.addRow("y", self.y)

        # Orientation (deg in UI, rad in spec)
        self.phi_deg = _spinbox(math.degrees(s.orientation),
                                decimals=3, step=1.0,
                                vmin=-3600, vmax=3600, suffix="°")
        self.phi_deg.valueChanged.connect(self._on_phi)
        self.form.addRow("Orientation", self.phi_deg)

        # ── Velocity ───────────────────────────────────────────
        self.vx = _spinbox(s.velocity[0], suffix="m/s")
        self.vy = _spinbox(s.velocity[1], suffix="m/s")
        self.omega = _spinbox(s.angular_velocity, decimals=4, step=0.1,
                              suffix="rad/s")
        self.vx.valueChanged.connect(self._on_vel)
        self.vy.valueChanged.connect(self._on_vel)
        self.omega.valueChanged.connect(self._on_omega)
        self.form.addRow("vx", self.vx)
        self.form.addRow("vy", self.vy)
        self.form.addRow("ω", self.omega)

        # ── Shape params ───────────────────────────────────────
        if s.shape is not None:
            self.form.addRow(QLabel(
                f"<b>Shape</b> &nbsp; <span style='color:#6b7280'>"
                f"{s.shape.kind}</span>"))
            self._shape_widgets: dict[str, QDoubleSpinBox] = {}
            for key, val in s.shape.params.items():
                w = _spinbox(float(val), suffix="m")
                w.valueChanged.connect(
                    lambda v, k=key: self._on_shape_param(k, v)
                )
                self._shape_widgets[key] = w
                self.form.addRow(key, w)
            # Depth (out-of-plane). Belongs to the shape conceptually
            # — it's only consumed by the mass-property calculator.
            self.form.addRow("depth (z)", self.thickness_z)
            note = QLabel(
                "Note: the z dimension is used only to compute mass "
                "properties; it has no other effect on the simulation."
            )
            note.setWordWrap(True)
            f = note.font(); f.setPointSize(max(7, f.pointSize() - 2))
            note.setFont(f)
            note.setStyleSheet("color: #6b7280;")
            self.form.addRow("", note)
        else:
            self._shape_widgets = {}
            # No shape: still expose depth so the user can tweak the
            # mass-property thickness directly.
            self.form.addRow("depth (z)", self.thickness_z)

        # Apply visual states last so all widgets exist.
        self.setStyleSheet(_LOCKED_QSS)
        self._refresh_mass_visuals()
        if getattr(s, "locked", False):
            self.set_locked(True)
        # Ground is fully read-only: lock every field including name,
        # material and mass / inertia.
        if s.id == "ground":
            self.name.setReadOnly(True)
            self.name.setEnabled(False)
            for w in (self.material, self.density, self.thickness_z,
                      self.mass, self.inertia, self.vx, self.vy,
                      self.omega, self.override):
                w.setEnabled(False)

    # ──────────────────────────────────────────────────────────
    def set_locked(self, on: bool) -> None:
        """Freeze geometry-defining fields once the body is committed."""
        for w in (self.x, self.y, self.phi_deg):
            _set_locked(w, on)
        for w in self._shape_widgets.values():
            _set_locked(w, on)

    def set_position_field_enabled(self, on: bool) -> None:
        """Enable/disable just the (x, y, phi) trio — used for drafts."""
        for w in (self.x, self.y, self.phi_deg):
            _set_locked(w, not on)

    # ----------------------------------------------------------
    def _on_name(self):
        new = self.name.text()
        old = self._old_name
        if new == old:
            return
        self.spec.name = new
        self._old_name = new
        self.body_renamed.emit(self.spec.id, old, new)
        self._emit()

    def _on_mass(self, v):
        if not self.spec.mass_override:
            return    # ignored — auto-computed
        self.spec.mass = float(v)
        self._emit()

    def _on_inertia(self, v):
        if not self.spec.mass_override:
            return
        self.spec.inertia = float(v)
        self._emit()

    def _on_pos(self, _=0):
        self.spec.position = (self.x.value(), self.y.value()); self._emit()

    def _on_phi(self, v):
        self.spec.orientation = math.radians(v); self._emit()

    def _on_vel(self, _=0):
        self.spec.velocity = (self.vx.value(), self.vy.value()); self._emit()

    def _on_omega(self, v):
        self.spec.angular_velocity = float(v); self._emit()

    def _on_shape_param(self, key: str, value: float):
        self.spec.shape.params[key] = float(value)
        self._recompute_mass_props()
        self._emit()

    # ── Material / density handlers ───────────────────────────
    def _on_material(self, name: str):
        self.spec.material = name
        if name in MATERIALS and name != "Custom":
            rho = MATERIALS[name]
            self.density.blockSignals(True)
            self.density.setValue(rho)
            self.density.blockSignals(False)
            self.spec.density = rho
        self._recompute_mass_props()
        self._emit()

    def _on_density(self, v: float):
        self.spec.density = float(v)
        # Free-form density => Custom material
        if self.material.currentText() != "Custom":
            self.material.blockSignals(True)
            self.material.setCurrentText("Custom")
            self.material.blockSignals(False)
            self.spec.material = "Custom"
        self._recompute_mass_props()
        self._emit()

    def _on_thickness_z(self, v: float):
        self.spec.thickness_z = float(v)
        self._recompute_mass_props()
        self._emit()

    def _on_override_toggled(self, on: bool):
        self.spec.mass_override = bool(on)
        if not on:
            self._recompute_mass_props()
        self._refresh_mass_visuals()
        self._emit()

    # ── Mass-props auto helpers ───────────────────────────────
    def _recompute_mass_props(self) -> None:
        if self.spec.mass_override:
            return
        out = compute_mass_props(self.spec)
        if out is None:
            return
        m, j = out
        self.spec.mass, self.spec.inertia = m, j
        self.mass.blockSignals(True);    self.mass.setValue(m);    self.mass.blockSignals(False)
        self.inertia.blockSignals(True); self.inertia.setValue(j); self.inertia.blockSignals(False)

    def _refresh_mass_visuals(self) -> None:
        auto = not self.spec.mass_override
        _set_autoval(self.mass, auto)
        _set_autoval(self.inertia, auto)


# ──────────────────────────────────────────────────────────────────
# Marker
# ──────────────────────────────────────────────────────────────────

class MarkerEditor(EditorBase):
    def __init__(self, spec: MarkerSpec, model: ModelSpec, parent=None):
        self._model = model
        super().__init__(spec, "marker", parent)

    def _build(self) -> None:
        s: MarkerSpec = self.spec

        self.name = QLineEdit(s.name)
        self.name.editingFinished.connect(self._on_name)
        self.form.addRow("Name", self.name)

        # Owning body (combo of all bodies + ground). Ground is
        # already part of ``self._model.bodies`` (auto-injected by
        # ``ModelSpec``), so we just iterate to avoid a duplicate.
        self.body = QComboBox()
        for b in self._model.bodies:
            self.body.addItem(b.name or b.id, b.id)
        idx = self.body.findData(s.body_id)
        if idx >= 0:
            self.body.setCurrentIndex(idx)
        self.body.currentIndexChanged.connect(self._on_body)
        self.form.addRow("Body", self.body)

        # Coordinate labels switch between (x, y) when the marker
        # belongs to ground (defined in the global frame) and
        # (ξ, η) otherwise (defined in the body's local frame).
        self.xi  = _spinbox(s.local_position[0])
        self.eta = _spinbox(s.local_position[1])
        self.xi.valueChanged.connect(self._on_pos)
        self.eta.valueChanged.connect(self._on_pos)
        self._lbl_xi  = QLabel()
        self._lbl_eta = QLabel()
        self.form.addRow(self._lbl_xi,  self.xi)
        self.form.addRow(self._lbl_eta, self.eta)
        self._refresh_pos_labels()

        self.theta_deg = _spinbox(math.degrees(s.theta),
                                  decimals=3, step=1.0,
                                  vmin=-3600, vmax=3600, suffix="°")
        self.theta_deg.valueChanged.connect(self._on_theta)
        self.form.addRow("θ", self.theta_deg)

    def _on_name(self):    self.spec.name = self.name.text();          self._emit()
    def _on_body(self, _):
        self.spec.body_id = self.body.currentData()
        self._refresh_pos_labels()
        self._emit()
    def _on_pos(self, _=0):
        self.spec.local_position = (self.xi.value(), self.eta.value())
        self._emit()
    def _on_theta(self, v): self.spec.theta = math.radians(v);          self._emit()

    def _refresh_pos_labels(self) -> None:
        """Re-label coordinate fields based on the owning body."""
        is_ground = self.spec.body_id == "ground"
        if is_ground:
            self._lbl_xi.setText("x")
            self._lbl_eta.setText("y")
            self.xi.setSuffix("  m")
            self.eta.setSuffix("  m")
        else:
            self._lbl_xi.setText("ξ (local x)")
            self._lbl_eta.setText("η (local y)")
            self.xi.setSuffix("  m (ξ)")
            self.eta.setSuffix("  m (η)")


# ──────────────────────────────────────────────────────────────────
# Joint
# ──────────────────────────────────────────────────────────────────

_JOINT_KINDS = [
    "RevJoint", "TranJoint", "RevRevJoint", "RevTranJoint",
    "RigidJoint", "DiscJoint", "RelRotJoint", "RelTranJoint",
]


class JointEditor(EditorBase):
    def __init__(self, spec: JointSpec, model: ModelSpec, parent=None):
        self._model = model
        super().__init__(spec, "joint", parent)

    def _build(self) -> None:
        s: JointSpec = self.spec

        self.name = QLineEdit(s.name)
        self.name.editingFinished.connect(self._on_name)
        self.form.addRow("Name", self.name)

        self.kind = QComboBox()
        for k in _JOINT_KINDS:
            self.kind.addItem(k, k)
        self.kind.setCurrentText(s.kind)
        self.kind.currentTextChanged.connect(self._on_kind)
        self.form.addRow("Kind", self.kind)

        # Marker selectors
        self.mi = self._marker_combo(s.i_marker_id)
        self.mj = self._marker_combo(s.j_marker_id)
        self.mi.currentIndexChanged.connect(self._on_mi)
        self.mj.currentIndexChanged.connect(self._on_mj)
        self.form.addRow("Marker i", self.mi)
        self.form.addRow("Marker j", self.mj)

        # Free-form params (key → spinbox)
        if s.params:
            self.form.addRow(QLabel("<b>Params</b>"))
            self._param_widgets: dict[str, QDoubleSpinBox] = {}
            for k, v in s.params.items():
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    continue
                w = _spinbox(fv, decimals=6, step=0.01)
                w.valueChanged.connect(lambda val, key=k: self._on_param(key, val))
                self._param_widgets[k] = w
                self.form.addRow(k, w)

    def _marker_combo(self, current_id: str) -> QComboBox:
        cb = QComboBox()
        cb.addItem("(none)", "")
        for m in self._model.markers:
            label = m.name or m.id
            owner = self._model.body(m.body_id)
            if owner is not None:
                label = f"{label}  [{owner.name or owner.id}]"
            cb.addItem(label, m.id)
        idx = cb.findData(current_id)
        if idx >= 0:
            cb.setCurrentIndex(idx)
        return cb

    def _on_name(self):     self.spec.name = self.name.text();           self._emit()
    def _on_kind(self, t):  self.spec.kind = t;                           self._emit()
    def _on_mi(self, _):    self.spec.i_marker_id = self.mi.currentData(); self._emit()
    def _on_mj(self, _):    self.spec.j_marker_id = self.mj.currentData(); self._emit()
    def _on_param(self, k: str, v: float):
        self.spec.params[k] = float(v)
        self._emit()


# ──────────────────────────────────────────────────────────────────
# Force
# ──────────────────────────────────────────────────────────────────

_FORCE_KINDS = [
    "Weight", "PtpForce", "RotSdaForce",
    "LocalForce", "GlobalForce", "Torque", "UserForce",
]


class ForceEditor(EditorBase):
    def __init__(self, spec: ForceSpec, model: ModelSpec, parent=None):
        self._model = model
        super().__init__(spec, "force", parent)

    def _build(self) -> None:
        s: ForceSpec = self.spec

        self.name = QLineEdit(s.name)
        self.name.editingFinished.connect(self._on_name)
        self.form.addRow("Name", self.name)

        self.kind = QComboBox()
        for k in _FORCE_KINDS:
            self.kind.addItem(k, k)
        self.kind.setCurrentText(s.kind)
        self.kind.currentTextChanged.connect(self._on_kind)
        self.form.addRow("Kind", self.kind)

        self.mi = self._marker_combo(s.i_marker_id)
        self.mj = self._marker_combo(s.j_marker_id)
        self.mi.currentIndexChanged.connect(self._on_mi)
        self.mj.currentIndexChanged.connect(self._on_mj)
        self.form.addRow("Marker i", self.mi)
        self.form.addRow("Marker j", self.mj)

        if s.params:
            self.form.addRow(QLabel("<b>Params</b>"))
            for k, v in s.params.items():
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    continue
                w = _spinbox(fv, decimals=6, step=0.01)
                w.valueChanged.connect(lambda val, key=k: self._on_param(key, val))
                self.form.addRow(k, w)

    def _marker_combo(self, current_id: str) -> QComboBox:
        cb = QComboBox()
        cb.addItem("(none)", "")
        for m in self._model.markers:
            cb.addItem(m.name or m.id, m.id)
        idx = cb.findData(current_id)
        if idx >= 0:
            cb.setCurrentIndex(idx)
        return cb

    def _on_name(self):     self.spec.name = self.name.text();             self._emit()
    def _on_kind(self, t):  self.spec.kind = t;                             self._emit()
    def _on_mi(self, _):    self.spec.i_marker_id = self.mi.currentData(); self._emit()
    def _on_mj(self, _):    self.spec.j_marker_id = self.mj.currentData(); self._emit()
    def _on_param(self, k: str, v: float):
        self.spec.params[k] = float(v)
        self._emit()
