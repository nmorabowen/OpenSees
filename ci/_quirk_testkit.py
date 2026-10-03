"""Shared helpers for the quirk-lint self-tests (ci/test_quirk_<slug>.py, WP-162).

Pure functions and constants only -- no module-level state a test could mutate.
This file should rarely change: a new rule gets its own test file and its own
helpers there, so two PRs adding rules never edit the same file.
"""
from pathlib import Path

STAMP = "// LADRUNO-HEADER-START\n// LADRUNO-HEADER-END\n"


def tree(root: Path, files: dict) -> Path:
    """Write {relative path: text} under root; return root."""
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return root


def rel(root: Path):
    """The `rel` callback the check functions take: path -> root-relative posix string."""
    return lambda p: p.resolve().relative_to(root.resolve()).as_posix()
