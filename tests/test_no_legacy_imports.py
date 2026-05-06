"""Garde-fou: ensures no legacy ``PMD.*`` imports survive the src-layout migration."""
from __future__ import annotations
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
LEGACY = re.compile(r"^\s*(?:from|import)\s+PMD\b", re.MULTILINE)
SKIP = {"_references", ".git", "__pycache__", ".pytest_cache", ".venv"}


def test_no_legacy_pmd_imports() -> None:
    bad: list[str] = []
    for p in ROOT.rglob("*.py"):
        if any(s in p.parts for s in SKIP):
            continue
        if LEGACY.search(p.read_text(encoding="utf-8", errors="ignore")):
            bad.append(str(p.relative_to(ROOT)))
    assert not bad, "Legacy PMD.* imports found:\n  " + "\n  ".join(bad)
