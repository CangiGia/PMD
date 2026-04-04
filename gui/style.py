"""style.py — QSS themes for the PMD PostProcessor GUI.

Usage
-----
    from PMD.gui.style import apply_light_theme
    apply_light_theme(app)          # call once after QApplication is created
"""

from __future__ import annotations

from PySide6.QtWidgets import QApplication

# ---------------------------------------------------------------------------
# Colour tokens
# ---------------------------------------------------------------------------
_ACCENT       = "#0078d4"   # Fluent blue
_ACCENT_HOVER = "#106ebe"
_ACCENT_PRESS = "#005a9e"
_BG           = "#e8e8e8"   # window background (warm mid-grey, not pure white)
_SURFACE      = "#f0f0f0"   # widget surfaces (list, tree, input)
_BORDER       = "#c4c4c4"
_BORDER_FOCUS = "#0078d4"
_TEXT         = "#1a1a1a"
_TEXT_DIM     = "#606060"
_TEXT_DISABLED= "#a0a0a0"
_RADIUS       = "4px"
_SEL_BG       = "#cce4f7"   # selection background (light blue)
_SEL_TEXT     = "#003d6b"

_LIGHT_QSS = f"""
/* ── Window & base ─────────────────────────────────── */
QMainWindow, QWidget {{
    background-color: {_BG};
    color: {_TEXT};
    font-family: "Segoe UI", sans-serif;
    font-size: 9pt;
}}

/* ── Menu bar ───────────────────────────────────────── */
QMenuBar {{
    background-color: {_SURFACE};
    border-bottom: 1px solid {_BORDER};
    padding: 2px 4px;
}}
QMenuBar::item {{
    padding: 4px 10px;
    border-radius: {_RADIUS};
}}
QMenuBar::item:selected {{
    background-color: {_SEL_BG};
    color: {_SEL_TEXT};
}}
QMenu {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    padding: 4px 0px;
}}
QMenu::item {{
    padding: 5px 28px 5px 16px;
    border-radius: {_RADIUS};
}}
QMenu::item:selected {{
    background-color: {_SEL_BG};
    color: {_SEL_TEXT};
}}
QMenu::separator {{
    height: 1px;
    background: {_BORDER};
    margin: 4px 8px;
}}

/* ── Buttons ────────────────────────────────────────── */
QPushButton {{
    background-color: {_ACCENT};
    color: #ffffff;
    border: none;
    border-radius: {_RADIUS};
    padding: 6px 16px;
    min-height: 28px;
}}
QPushButton:hover {{
    background-color: {_ACCENT_HOVER};
}}
QPushButton:pressed {{
    background-color: {_ACCENT_PRESS};
}}
QPushButton:disabled {{
    background-color: #c8c8c8;
    color: {_TEXT_DISABLED};
}}

/* ── ComboBox ───────────────────────────────────────── */
QComboBox {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    padding: 3px 8px;
    min-height: 22px;
    selection-background-color: {_SEL_BG};
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
    width: 0px;
    height: 0px;
    margin-right: 6px;
}}
QComboBox QAbstractItemView {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    selection-background-color: {_SEL_BG};
    selection-color: {_SEL_TEXT};
    outline: none;
}}

/* ── Labels ─────────────────────────────────────────── */
QLabel {{
    color: {_TEXT};
    background: transparent;
}}

/* ── Tree widget (generic) ──────────────────────────── */
QTreeWidget {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    outline: none;
}}
QTreeWidget::item {{
    padding: 3px 4px;
    border-radius: 3px;
}}
QTreeWidget::item:selected {{
    background-color: {_SEL_BG};
    color: {_SEL_TEXT};
}}
QTreeWidget::item:hover {{
    background-color: #e5f1fb;
}}

/* ── Navigation sidebar tree ────────────────────────── */
#nav_tree {{
    background-color: transparent;
    border: none;
    outline: none;
}}
#nav_tree::item {{
    min-height: 32px;
    padding: 2px 6px;
    border-radius: 4px;
}}
#nav_tree::item:selected {{
    background-color: {_ACCENT};
    color: #ffffff;
}}
#nav_tree::item:hover:!selected {{
    background-color: rgba(0, 0, 0, 0.07);
}}
#nav_tree::branch {{
    width: 0px;
    height: 0px;
    margin: 0px;
    padding: 0px;
    border-image: none;
    image: none;
    background: transparent;
}}
QHeaderView::section {{
    background-color: {_BG};
    border: none;
    border-bottom: 1px solid {_BORDER};
    padding: 4px 6px;
    font-weight: bold;
    color: {_TEXT_DIM};
}}

/* ── List widget ────────────────────────────────────── */
QListWidget {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    outline: none;
}}
QListWidget::item {{
    padding: 4px 6px;
    border-radius: 3px;
}}
QListWidget::item:selected {{
    background-color: {_SEL_BG};
    color: {_SEL_TEXT};
}}
QListWidget::item:hover {{
    background-color: #e5f1fb;
}}

/* ── Splitter ───────────────────────────────────────── */
QSplitter::handle {{
    background-color: {_BORDER};
}}
QSplitter::handle:horizontal {{
    width: 4px;
}}
QSplitter::handle:vertical {{
    height: 4px;
}}
QSplitter::handle:hover {{
    background-color: {_ACCENT};
}}

/* ── Matplotlib NavigationToolbar ───────────────────── */
QToolBar {{
    background-color: {_BG};
    border: none;
    border-bottom: 1px solid {_BORDER};
    spacing: 2px;
    padding: 2px 4px;
}}
QToolBar QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: {_RADIUS};
    padding: 3px 5px;
    color: {_TEXT};
}}
QToolBar QToolButton:hover {{
    background-color: {_SEL_BG};
    border-color: {_BORDER};
}}
QToolBar QToolButton:checked {{
    background-color: {_SEL_BG};
    border-color: {_ACCENT};
}}
QToolBar QToolButton:pressed {{
    background-color: #b8d6ee;
}}
QToolBar::separator {{
    width: 1px;
    background: {_BORDER};
    margin: 4px 3px;
}}

/* ── Status bar ─────────────────────────────────────── */
QStatusBar {{
    background-color: {_SURFACE};
    border-top: 1px solid {_BORDER};
    color: {_TEXT_DIM};
    font-size: 8pt;
    padding: 2px 6px;
}}

/* ── Scroll bars ─────────────────────────────────────── */
QScrollBar:vertical {{
    background: {_BG};
    width: 10px;
    border-radius: 5px;
}}
QScrollBar::handle:vertical {{
    background: {_BORDER};
    border-radius: 5px;
    min-height: 20px;
}}
QScrollBar::handle:vertical:hover {{
    background: #b0b0b0;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0px;
}}
QScrollBar:horizontal {{
    background: {_BG};
    height: 10px;
    border-radius: 5px;
}}
QScrollBar::handle:horizontal {{
    background: {_BORDER};
    border-radius: 5px;
    min-width: 20px;
}}
QScrollBar::handle:horizontal:hover {{
    background: #b0b0b0;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0px;
}}

/* ── Checkboxes ──────────────────────────────────────── */
QTreeWidget::indicator,
QListWidget::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_BORDER};
    border-radius: 3px;
    background-color: white;
}}
QTreeWidget::indicator:hover,
QListWidget::indicator:hover {{
    border-color: {_ACCENT};
}}
QTreeWidget::indicator:checked,
QListWidget::indicator:checked {{
    background-color: {_ACCENT};
    border-color: {_ACCENT};
}}

/* ── Panel section headers ───────────────────────────── */
QLabel#panel_header {{
    color: {_TEXT};
    font-size: 10pt;
    font-weight: bold;
    padding: 6px 8px 4px 8px;
    background-color: transparent;
    border-bottom: 1px solid {_BORDER};
}}

/* ── Sidebar ─────────────────────────────────────────── */
#sim_sidebar {{
    background-color: #d2d2d2;
    border-right: 1px solid {_BORDER};
}}

/* ── Filter card ─────────────────────────────────────── */
#filter_card {{
    background-color: {_SURFACE};
    border: 1px solid {_BORDER};
    border-radius: {_RADIUS};
    padding: 0px 2px;
}}

/* ── QCheckBox ──────────────────────────────────────── */
QCheckBox {{
    color: {_TEXT};
    spacing: 6px;
}}
QCheckBox::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_BORDER};
    border-radius: 3px;
    background-color: white;
}}
QCheckBox::indicator:hover  {{ border-color: {_ACCENT}; }}
QCheckBox::indicator:checked {{
    background-color: {_ACCENT};
    border-color: {_ACCENT};
}}

/* ── Settings footer ────────────────────────────────── */
#settings_footer {{
    background-color: #c8c8c8;
    border-top: 1px solid {_BORDER};
}}
"""

# ---------------------------------------------------------------------------
# Dark theme tokens
# ---------------------------------------------------------------------------
_D_ACCENT       = "#4ca3e0"
_D_ACCENT_HOVER = "#6ab8f0"
_D_ACCENT_PRESS = "#2c85c4"
_D_BG           = "#1e1e1e"
_D_SURFACE      = "#2a2a2a"
_D_BORDER       = "#3d3d3d"
_D_BORDER_FOCUS = "#4ca3e0"
_D_TEXT         = "#e0e0e0"
_D_TEXT_DIM     = "#888888"
_D_TEXT_DISABLED= "#555555"
_D_SEL_BG       = "#1e4a72"
_D_SEL_TEXT     = "#d0eaff"

_DARK_QSS = f"""
/* ── Window & base ──────────────────────────────────── */
QMainWindow, QWidget {{
    background-color: {_D_BG};
    color: {_D_TEXT};
    font-family: "Segoe UI", sans-serif;
    font-size: 9pt;
}}

/* ── Menu bar ─────────────────────────────────────── */
QMenuBar {{
    background-color: {_D_SURFACE};
    border-bottom: 1px solid {_D_BORDER};
    padding: 2px 4px;
}}
QMenuBar::item {{
    padding: 4px 10px;
    border-radius: {_RADIUS};
}}
QMenuBar::item:selected {{
    background-color: {_D_SEL_BG};
    color: {_D_SEL_TEXT};
}}
QMenu {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    padding: 4px 0px;
}}
QMenu::item {{
    padding: 5px 28px 5px 16px;
    border-radius: {_RADIUS};
}}
QMenu::item:selected {{
    background-color: {_D_SEL_BG};
    color: {_D_SEL_TEXT};
}}
QMenu::separator {{
    height: 1px;
    background: {_D_BORDER};
    margin: 4px 8px;
}}

/* ── Buttons ────────────────────────────────────────── */
QPushButton {{
    background-color: {_D_ACCENT};
    color: #ffffff;
    border: none;
    border-radius: {_RADIUS};
    padding: 6px 16px;
    min-height: 28px;
}}
QPushButton:hover {{
    background-color: {_D_ACCENT_HOVER};
}}
QPushButton:pressed {{
    background-color: {_D_ACCENT_PRESS};
}}
QPushButton:disabled {{
    background-color: #3a3a3a;
    color: {_D_TEXT_DISABLED};
}}

/* ── ComboBox ───────────────────────────────────────── */
QComboBox {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    padding: 3px 8px;
    min-height: 22px;
    color: {_D_TEXT};
    selection-background-color: {_D_SEL_BG};
}}
QComboBox:focus {{
    border-color: {_D_BORDER_FOCUS};
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
    width: 0px;
    height: 0px;
    margin-right: 6px;
}}
QComboBox QAbstractItemView {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    selection-background-color: {_D_SEL_BG};
    selection-color: {_D_SEL_TEXT};
    color: {_D_TEXT};
    outline: none;
}}

/* ── Labels ────────────────────────────────────────── */
QLabel {{
    color: {_D_TEXT};
    background: transparent;
}}

/* ── Tree widget (generic) ──────────────────────────── */
QTreeWidget {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    color: {_D_TEXT};
    outline: none;
}}
QTreeWidget::item {{
    padding: 3px 4px;
    border-radius: 3px;
}}
QTreeWidget::item:selected {{
    background-color: {_D_SEL_BG};
    color: {_D_SEL_TEXT};
}}
QTreeWidget::item:hover {{
    background-color: #2e3e50;
}}

/* ── Navigation sidebar tree ────────────────────────── */
#nav_tree {{
    background-color: transparent;
    border: none;
    outline: none;
    color: {_D_TEXT};
}}
#nav_tree::item {{
    min-height: 32px;
    padding: 2px 6px;
    border-radius: 4px;
}}
#nav_tree::item:selected {{
    background-color: {_D_ACCENT};
    color: #ffffff;
}}
#nav_tree::item:hover:!selected {{
    background-color: rgba(255, 255, 255, 0.07);
}}
#nav_tree::branch {{
    width: 0px;
    height: 0px;
    margin: 0px;
    padding: 0px;
    border-image: none;
    image: none;
    background: transparent;
}}
QHeaderView::section {{
    background-color: {_D_BG};
    border: none;
    border-bottom: 1px solid {_D_BORDER};
    padding: 4px 6px;
    font-weight: bold;
    color: {_D_TEXT_DIM};
}}

/* ── List widget ─────────────────────────────────────── */
QListWidget {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    color: {_D_TEXT};
    outline: none;
}}
QListWidget::item {{
    padding: 4px 6px;
    border-radius: 3px;
}}
QListWidget::item:selected {{
    background-color: {_D_SEL_BG};
    color: {_D_SEL_TEXT};
}}
QListWidget::item:hover {{
    background-color: #2e3e50;
}}

/* ── Splitter ────────────────────────────────────────── */
QSplitter::handle {{
    background-color: {_D_BORDER};
}}
QSplitter::handle:horizontal {{
    width: 4px;
}}
QSplitter::handle:vertical {{
    height: 4px;
}}
QSplitter::handle:hover {{
    background-color: {_D_ACCENT};
}}
/* ── Matplotlib NavigationToolbar (dark) ────────────── */
QToolBar {{
    background-color: {_D_BG};
    border: none;
    border-bottom: 1px solid {_D_BORDER};
    spacing: 2px;
    padding: 2px 4px;
}}
QToolBar QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: {_RADIUS};
    padding: 3px 5px;
    color: {_D_TEXT};
}}
QToolBar QToolButton:hover {{
    background-color: {_D_SEL_BG};
    border-color: {_D_BORDER};
}}
QToolBar QToolButton:checked {{
    background-color: {_D_SEL_BG};
    border-color: {_D_ACCENT};
}}
QToolBar QToolButton:pressed {{
    background-color: #1a3a5c;
}}
QToolBar::separator {{
    width: 1px;
    background: {_D_BORDER};
    margin: 4px 3px;
}}
/* ── Status bar ───────────────────────────────────────── */
QStatusBar {{
    background-color: {_D_SURFACE};
    border-top: 1px solid {_D_BORDER};
    color: {_D_TEXT_DIM};
    font-size: 8pt;
    padding: 2px 6px;
}}

/* ── Scroll bars ───────────────────────────────────────── */
QScrollBar:vertical {{
    background: {_D_BG};
    width: 10px;
    border-radius: 5px;
}}
QScrollBar::handle:vertical {{
    background: {_D_BORDER};
    border-radius: 5px;
    min-height: 20px;
}}
QScrollBar::handle:vertical:hover {{
    background: #606060;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0px;
}}
QScrollBar:horizontal {{
    background: {_D_BG};
    height: 10px;
    border-radius: 5px;
}}
QScrollBar::handle:horizontal {{
    background: {_D_BORDER};
    border-radius: 5px;
    min-width: 20px;
}}
QScrollBar::handle:horizontal:hover {{
    background: #606060;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0px;
}}

/* ── Checkboxes ──────────────────────────────────────── */
QTreeWidget::indicator,
QListWidget::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_D_BORDER};
    border-radius: 3px;
    background-color: {_D_SURFACE};
}}
QTreeWidget::indicator:hover,
QListWidget::indicator:hover {{
    border-color: {_D_ACCENT};
}}
QTreeWidget::indicator:checked,
QListWidget::indicator:checked {{
    background-color: {_D_ACCENT};
    border-color: {_D_ACCENT};
}}

/* ── Panel section headers ───────────────────────────── */
QLabel#panel_header {{
    color: {_D_TEXT};
    font-size: 10pt;
    font-weight: bold;
    padding: 6px 8px 4px 8px;
    background-color: transparent;
    border-bottom: 1px solid {_D_BORDER};
}}

/* ── Sidebar ─────────────────────────────────────────── */
#sim_sidebar {{
    background-color: #161616;
    border-right: 1px solid {_D_BORDER};
}}

/* ── Filter card ─────────────────────────────────────── */
#filter_card {{
    background-color: {_D_SURFACE};
    border: 1px solid {_D_BORDER};
    border-radius: {_RADIUS};
    padding: 0px 2px;
}}

/* ── QCheckBox ──────────────────────────────────────── */
QCheckBox {{
    color: {_D_TEXT};
    spacing: 6px;
}}
QCheckBox::indicator {{
    width: 14px;
    height: 14px;
    border: 1.5px solid {_D_BORDER};
    border-radius: 3px;
    background-color: {_D_SURFACE};
}}
QCheckBox::indicator:hover  {{ border-color: {_D_ACCENT}; }}
QCheckBox::indicator:checked {{
    background-color: {_D_ACCENT};
    border-color: {_D_ACCENT};
}}

/* ── Settings footer ────────────────────────────────── */
#settings_footer {{
    background-color: #252525;
    border-top: 1px solid {_D_BORDER};
}}
"""

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def apply_light_theme(app: QApplication) -> None:
    """Apply the flat light theme QSS to *app*."""
    app.setStyleSheet(_LIGHT_QSS)


def apply_dark_theme(app: QApplication) -> None:
    """Apply the flat dark theme QSS to *app*."""
    app.setStyleSheet(_DARK_QSS)
