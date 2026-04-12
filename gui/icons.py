"""icons.py — theme-aware qtawesome icon factory for the PMD GUI."""

from __future__ import annotations

import qtawesome as qta
from PySide6.QtGui import QIcon

from .style import ICON_COLOR_LIGHT, ICON_COLOR_DARK, ICON_COLOR_DISABLED

_dark: bool = False


def set_dark(dark: bool) -> None:
    """Update the global dark-mode flag (call before requesting icons)."""
    global _dark
    _dark = dark


def icon(name: str) -> QIcon:
    """Return a theme-coloured QIcon for the given *mdi6.* icon name."""
    color = ICON_COLOR_DARK if _dark else ICON_COLOR_LIGHT
    return qta.icon(name, color=color, color_disabled=ICON_COLOR_DISABLED)
