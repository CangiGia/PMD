"""GUI dialogs (modal) for the PMD pre-processor."""

from .solver_dialog import SolverDialog, SolverSettings
from .animation_dialog import AnimationDialog, PreprocessorAnimationCanvas

__all__ = [
    "SolverDialog",
    "SolverSettings",
    "AnimationDialog",
    "PreprocessorAnimationCanvas",
]
