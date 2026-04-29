"""Typed editors for the Inspector panel.

Each editor exposes a Qt form (label + widget pairs) bound to one
dataclass instance from :mod:`PMD.gui.preprocessor.models`. Edits are
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
)


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

    def __init__(self, spec: BodySpec, parent=None):
        super().__init__(spec, "body", parent)

    def _build(self) -> None:
        s: BodySpec = self.spec

        self.name = QLineEdit(s.name)
        self.name.editingFinished.connect(self._on_name)
        self.form.addRow("Name", self.name)

        self.mass = _spinbox(s.mass, decimals=4, step=0.1, vmin=1e-6,
                             suffix="kg")
        self.mass.valueChanged.connect(self._on_mass)
        self.form.addRow("Mass", self.mass)

        self.inertia = _spinbox(s.inertia, decimals=6, step=0.001, vmin=1e-9,
                                suffix="kg·m²")
        self.inertia.valueChanged.connect(self._on_inertia)
        self.form.addRow("Inertia", self.inertia)

        # Position
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

        # Velocity
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

        # Shape params (depend on kind)
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

    # ----------------------------------------------------------
    def _on_name(self):     self.spec.name = self.name.text();           self._emit()
    def _on_mass(self, v):  self.spec.mass = float(v);                   self._emit()
    def _on_inertia(self,v):self.spec.inertia = float(v);                self._emit()
    def _on_pos(self, _=0): self.spec.position = (self.x.value(), self.y.value()); self._emit()
    def _on_phi(self, v):   self.spec.orientation = math.radians(v);    self._emit()
    def _on_vel(self, _=0): self.spec.velocity = (self.vx.value(), self.vy.value()); self._emit()
    def _on_omega(self,v):  self.spec.angular_velocity = float(v);       self._emit()

    def _on_shape_param(self, key: str, value: float):
        self.spec.shape.params[key] = float(value)
        self._emit()


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

        # Owning body (combo of all bodies + ground)
        self.body = QComboBox()
        self.body.addItem("ground", "ground")
        for b in self._model.bodies:
            self.body.addItem(b.name or b.id, b.id)
        idx = self.body.findData(s.body_id)
        if idx >= 0:
            self.body.setCurrentIndex(idx)
        self.body.currentIndexChanged.connect(self._on_body)
        self.form.addRow("Body", self.body)

        self.xi  = _spinbox(s.local_position[0], suffix="m (ξ)")
        self.eta = _spinbox(s.local_position[1], suffix="m (η)")
        self.xi.valueChanged.connect(self._on_pos)
        self.eta.valueChanged.connect(self._on_pos)
        self.form.addRow("ξ (local x)", self.xi)
        self.form.addRow("η (local y)", self.eta)

        self.theta_deg = _spinbox(math.degrees(s.theta),
                                  decimals=3, step=1.0,
                                  vmin=-3600, vmax=3600, suffix="°")
        self.theta_deg.valueChanged.connect(self._on_theta)
        self.form.addRow("θ", self.theta_deg)

    def _on_name(self):    self.spec.name = self.name.text();          self._emit()
    def _on_body(self, _): self.spec.body_id = self.body.currentData(); self._emit()
    def _on_pos(self, _=0):
        self.spec.local_position = (self.xi.value(), self.eta.value())
        self._emit()
    def _on_theta(self, v): self.spec.theta = math.radians(v);          self._emit()


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
