"""style.py — QSS themes for the PMD PostProcessor GUI.

Usage
-----
    from PMD.gui.style import apply_light_theme
    apply_light_theme(app)          # call once after QApplication is created
"""

from __future__ import annotations

from PySide6.QtWidgets import QApplication

# ---------------------------------------------------------------------------
# Light theme colour tokens
# ---------------------------------------------------------------------------
_ACCENT       = "#0060c0"
_ACCENT_HOVER = "#0052a3"
_ACCENT_PRESS = "#004080"
_BG           = "#ebebeb"
_SURFACE      = "#f6f6f6"   # card surfaces — lighter than bg
_SURFACE_HI   = "#eeeeee"   # hover state
_BORDER       = "#c8c8c8"   # input borders only
_BORDER_FOCUS = "#0060c0"
_TEXT         = "#1a1a1a"
_TEXT_DIM     = "#5a5a5a"
_TEXT_DISABLED= "#a8a8a8"
_RADIUS       = "6px"

_LIGHT_QSS = f"""
/* ── Window & base ─────────────────────────────────── */
QMainWindow, QWidget {{
    background-color: {_BG};
    color: {_TEXT};
    font-family: "Segoe UI", sans-serif;
    font-size: 10pt;
}}

/* ── Menu bar ───────────────────────────────────────── */
QMenuBar {{
    background: {_SURFACE};
    border-bottom: 1px solid {_BORDER};
    padding: 2px 6px;
}}
QMenuBar::item {{
    padding: 4px 10px;
    border-radius: 4px;
}}
QMenuBar::item:selected {{
    background: rgba(0, 96, 192, 0.1);
    color: {_ACCENT};
}}
QMenu {{
    background: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    padding: 4px 0;
}}
QMenu::item {{
    padding: 5px 28px 5px 16px;
}}
QMenu::item:selected {{
    background: rgba(0, 96, 192, 0.1);
    color: {_ACCENT};
}}
QMenu::separator {{
    height: 1px;
    background: {_BORDER};
    margin: 4px 8px;
}}

/* ── Buttons — ghost default, primary opt-in ─────── */
QPushButton {{
    background: transparent;
    color: {_TEXT};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    padding: 5px 14px;
    min-height: 26px;
}}
QPushButton:hover {{
    background: {_SURFACE_HI};
    border-color: #a0a0a0;
}}
QPushButton:pressed {{
    background: #e0e0e0;
}}
QPushButton:disabled {{
    color: rgba(26, 26, 26, 0.35);
    border-color: #dcdcdc;
}}
QPushButton[primary="true"] {{
    background: {_ACCENT};
    border-color: transparent;
    color: #ffffff;
}}
QPushButton[primary="true"]:hover {{
    background: {_ACCENT_HOVER};
}}
QPushButton[primary="true"]:pressed {{
    background: {_ACCENT_PRESS};
}}
QPushButton[primary="true"]:disabled {{
    background: #c0d4e8;
    color: rgba(255, 255, 255, 0.5);
}}

/* ── ComboBox — flat ─────────────────────────────── */
QComboBox {{
    background: {_SURFACE};
    border: 1px solid transparent;
    border-radius: {_RADIUS};
    padding: 4px 8px;
    min-height: 24px;
    color: {_TEXT};
}}
QComboBox:hover {{
    background: {_SURFACE_HI};
}}
QComboBox:focus {{
    border-color: {_BORDER_FOCUS};
}}
QComboBox::drop-down {{
    border: none;
    width: 20px;
}}
QComboBox::down-arrow {{
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid {_TEXT_DIM};
    width: 0;
    height: 0;
    margin-right: 6px;
}}
QComboBox QAbstractItemView {{
    background: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    selection-background-color: rgba(0, 96, 192, 0.1);
    selection-color: {_ACCENT};
    color: {_TEXT};
    outline: none;
}}

/* ── Labels ─────────────────────────────────────────── */
QLabel {{
    color: {_TEXT};
    background: transparent;
}}

/* ── Tree widget ────────────────────────────────────── */
QTreeWidget {{
    background: {_SURFACE};
    border: none;
    outline: none;
    color: {_TEXT};
}}
QTreeWidget::item {{
    padding: 4px 8px;
    border-radius: 3px;
    min-height: 22px;
}}
QTreeWidget::item:selected {{
    background: rgba(0, 96, 192, 0.1);
    color: {_TEXT};
    border-left: 2px solid {_ACCENT};
    padding-left: 6px;
}}
QTreeWidget::item:hover:!selected {{
    background: rgba(0, 0, 0, 0.04);
}}
QHeaderView::section {{
    background: {_BG};
    border: none;
    border-bottom: 1px solid {_BORDER};
    padding: 4px 8px;
    font-weight: 600;
    color: {_TEXT_DIM};
}}

/* ── List widget ────────────────────────────────────── */
QListWidget {{
    background: {_SURFACE};
    border: none;
    outline: none;
    color: {_TEXT};
}}
QListWidget::item {{
    padding: 4px 8px;
    border-radius: 3px;
    min-height: 22px;
}}
QListWidget::item:selected {{
    background: rgba(0, 96, 192, 0.1);
    color: {_TEXT};
    border-left: 2px solid {_ACCENT};
    padding-left: 6px;
}}
QListWidget::item:hover:!selected {{
    background: rgba(0, 0, 0, 0.04);
}}

/* ── Splitter — invisible, 4 px hit area, accent on hover ── */
QSplitter::handle {{
    background: transparent;
}}
QSplitter::handle:horizontal {{
    width: 4px;
}}
QSplitter::handle:vertical {{
    height: 4px;
}}
QSplitter::handle:hover {{
    background: rgba(0, 96, 192, 0.4);
}}

/* ── Toolbar ─────────────────────────────────────────── */
QToolBar {{
    background: {_BG};
    border: none;
    spacing: 2px;
    padding: 2px 4px;
}}
QToolBar QToolButton {{
    background: transparent;
    border: 1px solid transparent;
    border-radius: {_RADIUS};
    padding: 4px 6px;
    color: {_TEXT};
    min-width: 28px;
    min-height: 28px;
}}
QToolBar QToolButton:hover {{
    background: rgba(0, 96, 192, 0.08);
    border-color: {_BORDER};
}}
QToolBar QToolButton:checked {{
    background: rgba(0, 96, 192, 0.14);
    border-color: {_ACCENT};
}}
QToolBar QToolButton:pressed {{
    background: rgba(0, 96, 192, 0.2);
}}
QToolBar::separator {{
    width: 1px;
    background: {_BORDER};
    margin: 4px 3px;
}}

/* ── Status bar ─────────────────────────────────────── */
QStatusBar {{
    background: {_SURFACE};
    border-top: 1px solid {_BORDER};
    color: {_TEXT_DIM};
    font-size: 8pt;
    padding: 2px 8px;
}}

/* ── Scroll bars — overlay style ─────────────────── */
QScrollBar:vertical {{
    background: transparent;
    width: 8px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {_BORDER};
    border-radius: 4px;
    min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{
    background: #a0a0a0;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
QScrollBar:horizontal {{
    background: transparent;
    height: 8px;
    margin: 0;
}}
QScrollBar::handle:horizontal {{
    background: {_BORDER};
    border-radius: 4px;
    min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{
    background: #a0a0a0;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}

/* ── Checkboxes ──────────────────────────────────────── */
QTreeWidget::indicator,
QListWidget::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_BORDER};
    border-radius: 3px;
    background: white;
}}
QTreeWidget::indicator:hover,
QListWidget::indicator:hover {{
    border-color: {_ACCENT};
}}
QTreeWidget::indicator:checked,
QListWidget::indicator:checked {{
    background: {_ACCENT};
    border-color: {_ACCENT};
}}

/* ── Panel section headers ───────────────────────────── */
QLabel#panel_header {{
    color: {_TEXT_DIM};
    font-size: 8pt;
    font-weight: 600;
    padding: 8px 12px 6px 12px;
    background: transparent;
    border-bottom: 1px solid #d8d8d8;
}}

/* ── Sidebar ─────────────────────────────────────────── */
#sim_sidebar {{
    background: #e0e0e0;
    border-right: 1px solid #d0d0d0;
}}

/* ── Filter card ─────────────────────────────────────── */
#filter_card {{
    background: {_SURFACE};
    border-radius: {_RADIUS};
    margin: 2px 4px;
}}

/* ── Units toolbar ───────────────────────────────────── */
#units_toolbar {{
    background: {_BG};
    border-bottom: 1px solid #d8d8d8;
}}
"""

# ---------------------------------------------------------------------------
# Dark theme colour tokens — 3-level surface hierarchy
# ---------------------------------------------------------------------------
_D_ACCENT       = "#4a9eff"
_D_ACCENT_HOVER = "#6ab5ff"
_D_ACCENT_PRESS = "#2c80e0"
_D_BG_DEEP    = "#161616"   # sidebar, status bar
_D_BG         = "#1e1e1e"   # main window background
_D_SURFACE    = "#252525"   # raised surfaces (cards, lists)
_D_SURFACE_HI = "#2d2d2d"   # hover on raised surfaces
_D_DIVIDER    = "#2a2a2a"   # subtle separators
_D_BORDER     = "#3a3a3a"   # input borders
_D_BORDER_FOCUS = "#4a9eff"
_D_TEXT         = "#e0e0e0"
_D_TEXT_DIM     = "#888888"
_D_TEXT_DISABLED= "rgba(224, 224, 224, 0.35)"

_DARK_QSS = f"""
/* ── Window & base ──────────────────────────────────── */
QMainWindow, QWidget {{
    background-color: {_D_BG};
    color: {_D_TEXT};
    font-family: "Segoe UI", sans-serif;
    font-size: 10pt;
}}

/* ── Menu bar ─────────────────────────────────────── */
QMenuBar {{
    background: {_D_SURFACE};
    border-bottom: 1px solid {_D_DIVIDER};
    padding: 2px 6px;
}}
QMenuBar::item {{
    padding: 4px 10px;
    border-radius: 4px;
}}
QMenuBar::item:selected {{
    background: rgba(74, 158, 255, 0.15);
    color: {_D_ACCENT};
}}
QMenu {{
    background: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: 6px;
    padding: 4px 0;
}}
QMenu::item {{
    padding: 5px 28px 5px 16px;
}}
QMenu::item:selected {{
    background: rgba(74, 158, 255, 0.15);
    color: {_D_ACCENT};
}}
QMenu::separator {{
    height: 1px;
    background: {_D_DIVIDER};
    margin: 4px 8px;
}}

/* ── Buttons — ghost default, primary opt-in ─────── */
QPushButton {{
    background: transparent;
    color: {_D_TEXT};
    border: 1px solid {_D_BORDER};
    border-radius: 6px;
    padding: 5px 14px;
    min-height: 26px;
}}
QPushButton:hover {{
    background: {_D_SURFACE_HI};
    border-color: #555555;
}}
QPushButton:pressed {{
    background: #222222;
}}
QPushButton:disabled {{
    color: {_D_TEXT_DISABLED};
    border-color: {_D_DIVIDER};
}}
QPushButton[primary="true"] {{
    background: {_D_ACCENT};
    border-color: transparent;
    color: #ffffff;
}}
QPushButton[primary="true"]:hover {{
    background: {_D_ACCENT_HOVER};
}}
QPushButton[primary="true"]:pressed {{
    background: {_D_ACCENT_PRESS};
}}
QPushButton[primary="true"]:disabled {{
    background: #1e3a5a;
    color: rgba(255, 255, 255, 0.35);
}}

/* ── ComboBox — flat ─────────────────────────────── */
QComboBox {{
    background: {_D_SURFACE_HI};
    border: 1px solid transparent;
    border-radius: 6px;
    padding: 4px 8px;
    min-height: 24px;
    color: {_D_TEXT};
}}
QComboBox:hover {{
    background: #333333;
}}
QComboBox:focus {{
    border-color: {_D_BORDER_FOCUS};
    background: {_D_SURFACE};
}}
QComboBox::drop-down {{
    border: none;
    width: 20px;
}}
QComboBox::down-arrow {{
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid {_D_TEXT_DIM};
    width: 0;
    height: 0;
    margin-right: 6px;
}}
QComboBox QAbstractItemView {{
    background: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: 6px;
    selection-background-color: rgba(74, 158, 255, 0.15);
    selection-color: {_D_ACCENT};
    color: {_D_TEXT};
    outline: none;
}}

/* ── Labels ────────────────────────────────────────── */
QLabel {{
    color: {_D_TEXT};
    background: transparent;
}}

/* ── Tree widget ────────────────────────────────────── */
QTreeWidget {{
    background: {_D_BG};
    border: none;
    outline: none;
    color: {_D_TEXT};
}}
QTreeWidget::item {{
    padding: 4px 8px;
    border-radius: 3px;
    min-height: 22px;
}}
QTreeWidget::item:selected {{
    background: rgba(74, 158, 255, 0.12);
    color: {_D_TEXT};
    border-left: 2px solid {_D_ACCENT};
    padding-left: 6px;
}}
QTreeWidget::item:hover:!selected {{
    background: rgba(255, 255, 255, 0.04);
}}
QHeaderView::section {{
    background: {_D_BG};
    border: none;
    border-bottom: 1px solid {_D_DIVIDER};
    padding: 4px 8px;
    font-weight: 600;
    color: {_D_TEXT_DIM};
}}

/* ── List widget ─────────────────────────────────────── */
QListWidget {{
    background: {_D_SURFACE};
    border: none;
    outline: none;
    color: {_D_TEXT};
}}
QListWidget::item {{
    padding: 4px 8px;
    border-radius: 3px;
    min-height: 22px;
}}
QListWidget::item:selected {{
    background: rgba(74, 158, 255, 0.12);
    color: {_D_TEXT};
    border-left: 2px solid {_D_ACCENT};
    padding-left: 6px;
}}
QListWidget::item:hover:!selected {{
    background: rgba(255, 255, 255, 0.04);
}}

/* ── Splitter — invisible, 4 px hit area, accent on hover ── */
QSplitter::handle {{
    background: transparent;
}}
QSplitter::handle:horizontal {{
    width: 4px;
}}
QSplitter::handle:vertical {{
    height: 4px;
}}
QSplitter::handle:hover {{
    background: rgba(74, 158, 255, 0.4);
}}

/* ── Toolbar ─────────────────────────────────────────── */
QToolBar {{
    background: {_D_BG};
    border: none;
    spacing: 2px;
    padding: 2px 4px;
}}
QToolBar QToolButton {{
    background: transparent;
    border: 1px solid transparent;
    border-radius: 6px;
    padding: 4px 6px;
    color: {_D_TEXT};
    min-width: 28px;
    min-height: 28px;
}}
QToolBar QToolButton:hover {{
    background: rgba(74, 158, 255, 0.1);
    border-color: {_D_BORDER};
}}
QToolBar QToolButton:checked {{
    background: rgba(74, 158, 255, 0.18);
    border-color: {_D_ACCENT};
}}
QToolBar QToolButton:pressed {{
    background: rgba(74, 158, 255, 0.25);
}}
QToolBar::separator {{
    width: 1px;
    background: {_D_DIVIDER};
    margin: 4px 3px;
}}

/* ── Status bar ───────────────────────────────────────── */
QStatusBar {{
    background: {_D_BG_DEEP};
    border-top: 1px solid {_D_DIVIDER};
    color: {_D_TEXT_DIM};
    font-size: 8pt;
    padding: 2px 8px;
}}

/* ── Scroll bars — overlay style ──────────────────── */
QScrollBar:vertical {{
    background: transparent;
    width: 8px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {_D_BORDER};
    border-radius: 4px;
    min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{
    background: #555555;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
QScrollBar:horizontal {{
    background: transparent;
    height: 8px;
    margin: 0;
}}
QScrollBar::handle:horizontal {{
    background: {_D_BORDER};
    border-radius: 4px;
    min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{
    background: #555555;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}

/* ── Checkboxes ──────────────────────────────────────── */
QTreeWidget::indicator,
QListWidget::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_D_BORDER};
    border-radius: 3px;
    background: {_D_SURFACE};
}}
QTreeWidget::indicator:hover,
QListWidget::indicator:hover {{
    border-color: {_D_ACCENT};
}}
QTreeWidget::indicator:checked,
QListWidget::indicator:checked {{
    background: {_D_ACCENT};
    border-color: {_D_ACCENT};
}}

/* ── Panel section headers ───────────────────────────── */
QLabel#panel_header {{
    color: {_D_TEXT_DIM};
    font-size: 8pt;
    font-weight: 600;
    padding: 8px 12px 6px 12px;
    background: transparent;
    border-bottom: 1px solid {_D_DIVIDER};
}}

/* ── Sidebar — deep background, contrast-only separation ── */
#sim_sidebar {{
    background: {_D_BG_DEEP};
}}

/* ── Filter card — raised surface, no border ─────────── */
#filter_card {{
    background: {_D_SURFACE};
    border-radius: 6px;
    margin: 2px 4px;
}}

/* ── Units toolbar ───────────────────────────────────── */
#units_toolbar {{
    background: {_D_BG};
    border-bottom: 1px solid {_D_DIVIDER};
}}
"""

# ---------------------------------------------------------------------------
# Public icon-colour constants  (consumed by gui/icons.py)
# ---------------------------------------------------------------------------
ICON_COLOR_LIGHT          = _TEXT          # "#1a1a1a"  — default on light
ICON_COLOR_DARK           = _D_TEXT        # "#e0e0e0"  — default on dark
ICON_COLOR_LIGHT_DIM      = _TEXT_DIM      # "#5a5a5a"  — secondary on light
ICON_COLOR_DARK_DIM       = _D_TEXT_DIM    # "#888888"  — secondary on dark
ICON_COLOR_DISABLED       = _TEXT_DISABLED # "#a8a8a8"  — light disabled
ICON_COLOR_DARK_DISABLED  = "#5a5a5a"      # ~30% of #e0e0e0 — dark disabled
ICON_COLOR_ACCENT         = _D_ACCENT      # "#4a9eff"  — active/checked (both themes)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def apply_light_theme(app: QApplication) -> None:
    """Apply the flat light theme QSS to *app*."""
    app.setStyleSheet(_LIGHT_QSS)


def apply_dark_theme(app: QApplication) -> None:
    """Apply the flat dark theme QSS to *app*."""
    app.setStyleSheet(_DARK_QSS)
