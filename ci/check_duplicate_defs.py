#!/usr/bin/env python3
"""Duplicate-definition gate (WP-162). Dependency-free.

  D1 redef      A `def` / `async def` / `class` whose name is defined AGAIN as a
                direct statement of the same scope (module, class or function
                body). Python keeps the last one silently, so every caller --
                including the tests written for the first -- runs the second.

The incident: WP-153 and WP-158 both took quirk rule "L9", and both added a
`def _l9(...)` helper to ci/test_check_quirk_patterns.py. `git merge` of the two
(the #899 merge-up of #901) completed with NO conflict and two `def _l9`; every
WP-153 test then called WP-158's helper.

Why not just `ruff check --select F811`? F811 ("redefinition of unused name")
does NOT fire on that merge: the first `_l9` is referenced by the test bodies
between the two definitions, so pyflakes counts it as used. Reproduced on the
real three-way merge of e2560dd2b and 5300da720 (ruff 0.16.10: "All checks
passed"). The CI step runs both: F811 catches what it catches (a duplicated
test name, a re-import), this catches the redefinition itself.

Not flagged: definitions inside `if`/`try`/`with` blocks (platform or
import-fallback alternatives); `@overload`; property `@x.setter/.getter/.deleter`;
`@x.register` (singledispatch). Waive one redefinition on its `def` line or the
line above, with a reason of at least 12 characters:
    # ladruno-lint: redef-ok <reason>

Usage:
    python ci/check_duplicate_defs.py [PATH ...]     # default: ci tests Ladruno_scripts
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ("ci", "tests", "Ladruno_scripts")
WAIVER = re.compile(r"#\s*ladruno-lint:\s*redef-ok\b(.*)$")
MIN_REASON = 12
DEFS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def _exempt(node) -> bool:
    for d in getattr(node, "decorator_list", []):
        target = d.func if isinstance(d, ast.Call) else d
        if isinstance(target, ast.Name) and target.id == "overload":
            return True
        if isinstance(target, ast.Attribute) and target.attr in ("setter", "getter", "deleter",
                                                                 "register", "overload"):
            return True
    return False


def _scopes(tree):
    """Every statement list that is a scope body: the module, each class, each function."""
    yield tree.body
    for node in ast.walk(tree):
        if isinstance(node, DEFS):
            yield node.body


def find(path: Path, root: Path = ROOT) -> list[str]:
    try:
        src = path.read_text(encoding="utf-8", errors="replace")
        tree = ast.parse(src, filename=str(path))
    except SyntaxError as e:
        return [f"D1 {_rel(path, root)}:{e.lineno}: cannot parse ({e.msg})"]
    raw = src.splitlines()
    out = []
    for body in _scopes(tree):
        first: dict[str, int] = {}
        for node in body:
            if not isinstance(node, DEFS) or _exempt(node):
                continue
            if node.name not in first:
                first[node.name] = node.lineno
                continue
            line = node.lineno - len(getattr(node, "decorator_list", []))
            waived = None
            for j in (node.lineno - 1, line - 2):
                if 0 <= j < len(raw):
                    m = WAIVER.search(raw[j])
                    if m:
                        waived = m.group(1).strip()
            if waived is not None and len(waived) >= MIN_REASON:
                continue
            why = "redef-ok waiver reason too short" if waived is not None else (
                f"'{node.name}' is defined again (first at line {first[node.name]}); Python keeps this "
                "one silently, so every caller of the first runs this. Rename one (a merge of two PRs "
                "that took the same name does exactly this), or waive with '# ladruno-lint: redef-ok <reason>'")
            out.append(f"D1 {_rel(path, root)}:{node.lineno}: {why}")
    return out


def _rel(p: Path, root: Path) -> str:
    try:
        return p.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return p.as_posix()


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    targets = [Path(a) for a in args] or [ROOT / d for d in DEFAULT]
    files = []
    for t in targets:
        files += [t] if t.is_file() else sorted(t.rglob("*.py"))
    findings = [f for p in files for f in find(p)]
    for f in findings:
        print(f)
    print(f"check_duplicate_defs: {len(findings)} finding(s) in {len(files)} file(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
