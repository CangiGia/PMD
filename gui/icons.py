"""icons.py — theme-aware qtawesome icon factory for the PMD GUI."""

from __future__ import annotations

import qtawesome as qta
from PySide6.QtGui import QIcon

from .style import (
    ICON_COLOR_LIGHT, ICON_COLOR_DARK,
    ICON_COLOR_LIGHT_DIM, ICON_COLOR_DARK_DIM,
    ICON_COLOR_DISABLED, ICON_COLOR_DARK_DISABLED,
    ICON_COLOR_ACCENT,
)

_dark: bool = False
_cache: dict[tuple, QIcon] = {}


def set_dark(dark: bool) -> None:
    """Switch dark mode and invalidate the icon cache."""
    global _dark
    _dark = dark
    clear_cache()


def clear_cache() -> None:
    """Discard all cached QIcon instances."""
    _cache.clear()


def icon(name: str, *, dim: bool = False) -> QIcon:
    """Return a theme-aware, cached QIcon for the given mdi6 icon name.

    Parameters
    ----------
    name:
        An mdi6 icon name, e.g. ``'mdi6.play-outline'``.
    dim:
        Use the secondary (dimmed) colour. Pass ``True`` for tree category
        headers, status-bar glyphs and other non-interactive decorations.
    """
    key = (name, _dark, dim)
    cached = _cache.get(key)
    if cached is not None:
        return cached

    if _dark:
        color    = ICON_COLOR_DARK_DIM    if dim else ICON_COLOR_DARK
        disabled = ICON_COLOR_DARK_DISABLED
    else:
        color    = ICON_COLOR_LIGHT_DIM   if dim else ICON_COLOR_LIGHT
        disabled = ICON_COLOR_DISABLED

    ic = qta.icon(
        name,
        color=color,
        color_active=ICON_COLOR_ACCENT,
        color_selected=ICON_COLOR_ACCENT,
        color_disabled=disabled,
    )
    _cache[key] = ic
    return ic
