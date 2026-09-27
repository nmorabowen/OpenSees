"""Shared helpers for the WP-120 code-health survey (read-only, dependency-free).

Reuses the C++ cleaner and function finder of the quirk lint
(`ci/check_quirk_patterns.py`: `clean()` blanks comments, string/char literals,
`#if 0` blocks and every preprocessor line, keeping line numbers; `functions()`
finds definitions by brace matching) instead of writing a second one.
"""
from __future__ import annotations

import importlib.util
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("WP120_ROOT", REPO))   # survey a git-archive tree (acceptance runs)
SUFFIXES = (".cpp", ".h", ".hpp", ".cc", ".cxx")
STAMP = "LADRUNO-HEADER-START"

_spec = importlib.util.spec_from_file_location("cqp", REPO / "ci" / "check_quirk_patterns.py")
cqp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cqp)
clean = cqp.clean
functions = cqp.functions


def rel(p: Path) -> str:
    return p.resolve().relative_to(ROOT).as_posix()


def read(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def all_sources():
    for p in sorted((ROOT / "SRC").rglob("*")):
        if p.is_file() and p.suffix in SUFFIXES:
            yield p


def split_fork_vanilla():
    """(fork, vanilla): fork = files carrying the LADRUNO-HEADER-START stamp."""
    fork, van = [], []
    for p in all_sources():
        (fork if STAMP in read(p) else van).append(p)
    return fork, van


TOKEN = re.compile(r"[A-Za-z_]\w*|\d[\w.]*(?:[eE][+-]?\d+)?|::|->|\+\+|--|<<=|>>=|<<|>>|[-+*/%&|^!=<>]=|&&|\|\||\S")
KEYWORDS = set("""alignas alignof and asm auto bool break case catch char class const constexpr
const_cast continue decltype default delete do double dynamic_cast else enum explicit
extern false float for friend goto if inline int long mutable namespace new noexcept
nullptr operator private protected public register reinterpret_cast return short signed
sizeof static static_assert static_cast struct switch template this throw true try typedef
typeid typename union unsigned using virtual void volatile while override final size_t
std Vector Matrix ID""".split())


def tokens(cleaned_lines):
    """[(token, line_index)] of cleaned source (comments/strings/preprocessor gone)."""
    out = []
    for li, ln in enumerate(cleaned_lines):
        for m in TOKEN.finditer(ln):
            out.append((m.group(0), li))
    return out


def code_lines(cleaned_lines):
    return sum(1 for ln in cleaned_lines if ln.strip())
