"""Modal dialog to configure and launch a simulation from the GUI."""

from __future__ import annotations

from dataclasses import dataclass

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


@dataclass
class SolverSettings:
    analysis: str = "dynamic"      # dynamic | kinematic | static
    t_start: float = 0.0
    t_end:   float = 2.0
    n_steps: int   = 2000
    ic_correct: bool = True


class SolverDialog(QDialog):
    """Collect :class:`SolverSettings` from the user."""

    ANALYSES = ["dynamic", "kinematic", "static"]

    def __init__(self, parent: QWidget | None = None,
                 initial: SolverSettings | None = None):
        super().__init__(parent)
        self.setWindowTitle("Run Simulation")
        self.setModal(True)
        self._settings = initial or SolverSettings()

        outer = QVBoxLayout(self)
        form = QFormLayout()
        outer.addLayout(form)

        self.cmb_analysis = QComboBox()
        self.cmb_analysis.addItems(self.ANALYSES)
        self.cmb_analysis.setCurrentText(self._settings.analysis)
        form.addRow("Analysis", self.cmb_analysis)

        self.sb_t0 = self._spin(self._settings.t_start, suffix="s")
        self.sb_tf = self._spin(self._settings.t_end,   suffix="s")
        form.addRow("t start", self.sb_t0)
        form.addRow("t end",   self.sb_tf)

        self.sb_n = QSpinBox()
        self.sb_n.setRange(2, 1_000_000)
        self.sb_n.setValue(self._settings.n_steps)
        form.addRow("Output steps", self.sb_n)

        # IC correction toggle
        from PySide6.QtWidgets import QCheckBox
        self.chk_ic = QCheckBox("Project ICs onto constraint manifold")
        self.chk_ic.setChecked(self._settings.ic_correct)
        outer.addWidget(self.chk_ic)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        outer.addWidget(bb)

    # ------------------------------------------------------------
    @staticmethod
    def _spin(value: float, *, decimals: int = 4,
              vmin: float = -1e6, vmax: float = 1e6,
              suffix: str = "") -> QDoubleSpinBox:
        sb = QDoubleSpinBox()
        sb.setDecimals(decimals)
        sb.setRange(vmin, vmax)
        sb.setValue(value)
        if suffix:
            sb.setSuffix(f"  {suffix}")
        sb.setKeyboardTracking(False)
        return sb

    # ------------------------------------------------------------
    def settings(self) -> SolverSettings:
        return SolverSettings(
            analysis=self.cmb_analysis.currentText(),
            t_start=self.sb_t0.value(),
            t_end=self.sb_tf.value(),
            n_steps=self.sb_n.value(),
            ic_correct=self.chk_ic.isChecked(),
        )
