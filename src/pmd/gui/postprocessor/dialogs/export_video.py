"""ExportVideoDialog — parameters dialog for video export."""

from __future__ import annotations

import shutil

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class ExportVideoDialog(QDialog):
    """Modal dialog that collects video-export parameters.

    Parameters
    ----------
    current_speed : float
        The current playback speed multiplier shown in the AnimationCanvas,
        used as the default for the video export speed field.
    parent : QWidget or None
    """

    def __init__(self, current_speed: float = 1.0, parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle("Export animation video")
        self.setMinimumWidth(480)

        # ── ffmpeg check ──────────────────────────────────────────────────────
        self._ffmpeg_ok = shutil.which("ffmpeg") is not None

        # ── form ─────────────────────────────────────────────────────────────
        form = QFormLayout()
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)

        # Output path
        path_row = QWidget()
        path_hl  = QHBoxLayout(path_row)
        path_hl.setContentsMargins(0, 0, 0, 0)
        self._path_edit = QLineEdit("animation.mp4")
        self._path_edit.setPlaceholderText("output file path…")
        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(80)
        browse_btn.clicked.connect(self._browse)
        path_hl.addWidget(self._path_edit, stretch=1)
        path_hl.addWidget(browse_btn)
        form.addRow("Output file", path_row)

        # Video FPS
        self._fps_spin = QSpinBox()
        self._fps_spin.setRange(10, 60)
        self._fps_spin.setValue(30)
        self._fps_spin.setSuffix(" fps")
        self._fps_spin.setToolTip(
            "Frame rate of the output video.\n"
            "30 fps is standard for presentations; 60 fps for smoother visuals.")
        form.addRow("Video FPS", self._fps_spin)

        # Playback speed multiplier
        self._speed_spin = QDoubleSpinBox()
        self._speed_spin.setDecimals(2)
        self._speed_spin.setRange(0.05, 100.0)
        self._speed_spin.setSingleStep(0.25)
        self._speed_spin.setValue(current_speed)
        self._speed_spin.setSuffix("×")
        self._speed_spin.setToolTip(
            "How many seconds of simulation time appear in one real second of video.\n"
            "1× = real time.  10× = ten times faster than real time.")
        form.addRow("Speed", self._speed_spin)

        # Layout
        self._layout_combo = QComboBox()
        self._layout_combo.addItem("Animation only",      userData="anim")
        self._layout_combo.addItem("Animation + Plots (side by side)", userData="combo")
        self._layout_combo.setToolTip(
            "Animation only — exports only the 2-D scene.\n"
            "Animation + Plots — places the scene and the plot canvas side by side;\n"
            "  the plot cursor scrolls in sync with the animation.")
        form.addRow("Layout", self._layout_combo)

        # DPI
        self._dpi_spin = QSpinBox()
        self._dpi_spin.setRange(50, 300)
        self._dpi_spin.setValue(100)
        self._dpi_spin.setSuffix(" dpi")
        self._dpi_spin.setToolTip(
            "Rendering resolution (dots per inch).\n"
            "Higher values give sharper output but increase file size and render time.")
        form.addRow("Resolution", self._dpi_spin)

        # ── buttons ───────────────────────────────────────────────────────────
        bbox = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        bbox.accepted.connect(self._on_accept)
        bbox.rejected.connect(self.reject)
        self._ok_btn = bbox.button(QDialogButtonBox.StandardButton.Ok)
        self._ok_btn.setText("Export")

        # ── warning if ffmpeg missing ─────────────────────────────────────────
        self._warn_lbl = QLabel(
            "⚠  ffmpeg not found in PATH.  Install ffmpeg and make sure it is "
            "accessible from the command line before exporting."
        )
        self._warn_lbl.setWordWrap(True)
        self._warn_lbl.setStyleSheet("color: #b85c00; font-style: italic;")
        self._warn_lbl.setVisible(not self._ffmpeg_ok)
        if not self._ffmpeg_ok:
            self._ok_btn.setEnabled(False)

        # ── layout ────────────────────────────────────────────────────────────
        vl = QVBoxLayout(self)
        vl.addLayout(form)
        vl.addSpacing(4)
        vl.addWidget(self._warn_lbl)
        vl.addSpacing(8)
        vl.addWidget(bbox)

    # ------------------------------------------------------------------
    # slots
    # ------------------------------------------------------------------

    def _browse(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save video as", self._path_edit.text(),
            "MP4 video (*.mp4);;AVI video (*.avi);;All files (*)"
        )
        if path:
            self._path_edit.setText(path)

    def _on_accept(self) -> None:
        if not self._path_edit.text().strip():
            QMessageBox.warning(self, "Export video", "Please specify an output file path.")
            return
        self.accept()

    # ------------------------------------------------------------------
    # Result accessors (call after exec() returns Accepted)
    # ------------------------------------------------------------------

    @property
    def output_path(self) -> str:
        return self._path_edit.text().strip()

    @property
    def video_fps(self) -> int:
        return self._fps_spin.value()

    @property
    def speed(self) -> float:
        return self._speed_spin.value()

    @property
    def layout(self) -> str:
        """``"anim"`` or ``"combo"``."""
        return self._layout_combo.currentData()

    @property
    def dpi(self) -> int:
        return self._dpi_spin.value()
