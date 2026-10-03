#!/usr/bin/env python3
"""Ledger-fragment gate (WP-161). Dependency-free.

The build-control ledgers are one fragment per entry under
`Ladruno_implementation/ledger/<kind>/` (see ci/ledger.py), so no two PRs write
the same file. This gate keeps the fragments well-formed:

  F1 schema     every fragment parses: front matter between `---` lines, one
                `key: value` per line (a JSON value or a plain bare string),
                known keys only, `wp` and `title` present; non-migrated
                fragments (no `legacy_seq`) also need `date: YYYY-MM-DD`, and an
                implementations row needs `status`. The body is non-empty and
                has the shape of its kind: one table row (implementations,
                section table), `- ` bullets (section history), a `### `
                heading first (quirks; migrated entries may be `## `), table
                rows only, at least 3 cells each (vanilla; 2 for a few migrated rows).
  F2 name       the file name is `<wp>-<slug>.md` for the fragment's own `wp`
                (`WP-<n>`, `ADR-<n>`, `PR-<n>` or `LEGACY`), slug lowercase
                `a-z0-9-`.
  F3 unique     no two fragments of a kind share a name (case-insensitively:
                Windows checkouts), and no two share a body or (quirks) a
                heading -- a duplicated entry is the union-merge failure this
                layout replaces. Duplicates that are BOTH migrated legacy are
                reported as warnings: they were already in the old ledger.
  F4 classTag   `class_tags` of an implementations fragment name symbols in
                SRC/classTags.h (parsed by ci/check_classtags.py): an unknown
                symbol is an error unless the entry says RESERVED or PLANNED
                (a warning for migrated entries), `SYMBOL=value` must match the
                header, and two current fragments may not claim one symbol.
  F5 banner     a `banner` value is a line of Ladruno_scripts/banner_features.txt.
  F6 template   each kind has `_template.md` with its placeholders, and
                `python ci/ledger.py build` renders all three ledgers.
  F7 stub       a committed `LEDGER_*.md` that is a stub is exactly the stub
                ci/ledger.py writes. A stale branch that union-merges rows into
                a stub, or an agent who edits it, fails here: write a fragment.

Usage:
    python ci/check_ledger_fragments.py            # exit 1 on any error
    python ci/check_ledger_fragments.py --root DIR
"""
from __future__ import annotations

import argparse
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_classtags  # noqa: E402
import ledger as L  # noqa: E402

ALLOWED = set(L.KEYS)
SECTIONS = {"implementations": {"table", "history"}, "vanilla": {"main", "upstreamable"}}


def _err(errors, kind, name, msg):
    errors.append(f"{kind}/{name}.md: {msg}")


def check_schema(f: L.Fragment, errors, warns):
    kind, name, meta, body = f.kind, f.name, f.meta, f.body
    for k in meta:
        if k not in ALLOWED:
            _err(errors, kind, name, f"F1 unknown front-matter key {k!r} (allowed: {', '.join(L.KEYS)})")
    for k in ("wp", "title"):
        if not isinstance(meta.get(k), str) or not meta.get(k).strip():
            _err(errors, kind, name, f"F1 front matter needs a non-empty {k!r}")
    if not body.strip():
        _err(errors, kind, name, "F1 empty body")
        return
    if not f.legacy:
        if not isinstance(meta.get("date"), str) or not L.DATE_RE.match(meta.get("date", "")):
            _err(errors, kind, name, "F1 a new fragment needs `date: YYYY-MM-DD`")
    elif not (isinstance(meta["legacy_seq"], int) or
              (isinstance(meta["legacy_seq"], list) and all(isinstance(x, int) for x in meta["legacy_seq"]))):
        _err(errors, kind, name, "F1 legacy_seq must be an integer (or a list of them, vanilla)")
    for k in ("files", "class_tags"):
        if k in meta and not (isinstance(meta[k], list) and all(isinstance(x, str) for x in meta[k])):
            _err(errors, kind, name, f"F1 {k} must be a JSON list of strings")
    lines = body.split("\n")
    if kind == "implementations":
        sec = meta.get("section", "table")
        if sec not in SECTIONS[kind]:
            _err(errors, kind, name, f"F1 section must be one of {sorted(SECTIONS[kind])}")
        elif sec == "table":
            if len(lines) != 1 or not L.is_row(lines[0]):
                _err(errors, kind, name, "F1 an implementations fragment is ONE table row `| ... |`")
            if not f.legacy and not meta.get("status"):
                _err(errors, kind, name, "F1 an implementations row needs `status`")
        elif not lines[0].startswith("- ") or any(ln.strip() and not ln.startswith(("- ", " ", "\t")) for ln in lines):
            _err(errors, kind, name, "F1 a history fragment is `- ` bullet(s)")
    elif kind == "quirks":
        first = lines[0]
        ok = first.startswith("### ") if not f.legacy else bool(L.HEADING_RE.match(first))
        if not ok:
            _err(errors, kind, name, "F1 a quirk starts with its `### ` heading")
        fence = None
        for ln in lines[1:]:
            fm = L.FENCE_RE.match(ln)
            if fm:
                fence = None if fence == fm.group(1) else (fence or fm.group(1))
            elif fence is None and L.HEADING_RE.match(ln) and not f.legacy:
                _err(errors, kind, name, f"F1 a second `##`/`###` heading inside one quirk ({ln[:60]!r}); "
                                         "use `####` or split it into two fragments")
                break
    else:
        if meta.get("table", "main") not in SECTIONS[kind]:
            _err(errors, kind, name, f"F1 table must be one of {sorted(SECTIONS[kind])}")
        rows = [ln for ln in lines if ln.strip()]
        need = 2 if f.legacy else 3     # a few migrated rows never had a PR cell
        bad = [ln for ln in rows if not L.is_row(ln) or len(L.cells(ln)) < need]
        if bad:
            _err(errors, kind, name, f"F1 a vanilla fragment holds table rows only (`| file | why | PR |`): {bad[0][:60]!r}")
        if "files" in meta and isinstance(meta["files"], list) and len(meta["files"]) != len(rows):
            _err(errors, kind, name, f"F1 files lists {len(meta['files'])} file(s) for {len(rows)} row(s)")


def check_name(f: L.Fragment, errors):
    m = L.NAME_RE.match(f.name)
    wp = f.meta.get("wp", "")
    if not isinstance(wp, str) or not L.KEY_RE.match(wp):
        _err(errors, f.kind, f.name, f"F2 wp {wp!r} is not WP-<n>, ADR-<n>, PR-<n> or LEGACY")
    elif not m:
        _err(errors, f.kind, f.name, f"F2 file name must be `{wp}-<slug>.md` (slug: a-z 0-9 -)")
    elif m.group("key") != wp:
        _err(errors, f.kind, f.name, f"F2 file name key {m.group('key')!r} does not match wp {wp!r}")


def check_unique(kind, frags, errors, warns):
    by_name = defaultdict(list)
    by_body = defaultdict(list)
    by_head = defaultdict(list)
    for f in frags:
        by_name[f.name.lower()].append(f)
        by_body[f.body.strip()].append(f)
        if kind == "quirks":
            by_head[f.body.split("\n", 1)[0].strip()].append(f)
    for name, fs in by_name.items():
        if len(fs) > 1:
            errors.append(f"{kind}: F3 {len(fs)} fragments named {name}.md (case-insensitive)")
    for label, table in (("body", by_body), ("heading", by_head)):
        for _, fs in table.items():
            if len(fs) < 2:
                continue
            names = ", ".join(sorted(f.name for f in fs))
            msg = f"{kind}: F3 duplicate {label} in {names}"
            (warns if all(f.legacy for f in fs) else errors).append(
                msg + (" (both migrated from the old ledger)" if all(f.legacy for f in fs)
                       else " -- edit the older fragment instead of adding a copy"))


def check_classtags_of(root, frags, errors, warns):
    header = check_classtags.parse(root / "SRC" / "classTags.h")
    claimed = defaultdict(list)
    for f in frags:
        for tag in f.meta.get("class_tags", []) if isinstance(f.meta.get("class_tags"), list) else []:
            sym, _, val = tag.partition("=")
            sym = sym.strip()
            text = (f.body + " " + str(f.meta.get("status", ""))).upper()
            if not re.fullmatch(r"[A-Z][A-Z0-9]*_TAG_\w+", sym):
                _err(errors, f.kind, f.name, f"F4 class_tags entry {tag!r} is not a SYMBOL or SYMBOL=value")
                continue
            if not f.legacy:
                claimed[sym].append(f.name)
            if sym not in header:
                if "RESERVED" in text or "PLANNED" in text:
                    continue
                (warns if f.legacy else errors).append(
                    f"implementations/{f.name}.md: F4 class tag {sym} is not in SRC/classTags.h "
                    "(mark the entry RESERVED/PLANNED, or fix the symbol)")
            elif val.strip() and int(val) != header[sym][0]:
                _err(errors, f.kind, f.name, f"F4 {sym}={val.strip()} but SRC/classTags.h has {header[sym][0]}")
    for sym, names in claimed.items():
        if len(names) > 1:
            errors.append(f"implementations: F4 class tag {sym} claimed by {', '.join(sorted(names))}")


def check_banner(root, frags, errors):
    p = root / "Ladruno_scripts" / "banner_features.txt"
    lines = {ln.strip() for ln in p.read_text(encoding="utf-8").splitlines()} if p.exists() else set()
    for f in frags:
        b = f.meta.get("banner")
        if b is not None and (not isinstance(b, str) or b.strip() not in lines):
            _err(errors, f.kind, f.name, f"F5 banner {b!r} is not a line of Ladruno_scripts/banner_features.txt")


def check_templates_and_build(root, errors):
    need = {"implementations": ["<!-- ledger:rows table -->", "<!-- ledger:history -->"],
            "quirks": ["<!-- ledger:entries -->"],
            "vanilla": ["<!-- ledger:rows main -->", "<!-- ledger:rows upstreamable -->"]}
    for kind, marks in need.items():
        try:
            tpl = L.load_template(root, kind)
        except L.LedgerError as e:
            errors.append(f"{kind}: F6 {e}")
            continue
        for mk in marks:
            if mk not in tpl.split("\n"):
                errors.append(f"{kind}: F6 _template.md lost its placeholder line {mk}")
    if any("F6" in e for e in errors):
        return
    with tempfile.TemporaryDirectory() as tmp:
        try:
            L.build(root, Path(tmp))
        except (L.LedgerError, OSError, ValueError) as e:
            errors.append(f"F6 `python ci/ledger.py build` fails: {e}")


def check_stubs(root, errors):
    for kind, spec in L.KINDS.items():
        p = root / L.IMPL_DIR / spec["ledger"]
        if not p.exists():
            errors.append(f"F7 {L.IMPL_DIR}/{spec['ledger']} is missing (keep the stub: wiki-links point at it)")
            continue
        text = p.read_text(encoding="utf-8").replace("\r\n", "\n")
        if L.STUB_MARK not in text:
            continue                    # not split yet: the old single-file ledger
        src = L.stub_source(text)
        if src is None or text != L.render_stub(kind, src):
            errors.append(f"F7 {L.IMPL_DIR}/{spec['ledger']} is generated -- do not edit it; write a fragment in "
                          f"{L.FRAG_DIR}/{kind}/ (restore the stub: `git checkout origin/ladruno -- {L.IMPL_DIR}/{spec['ledger']}`)")


def run(root: Path) -> tuple[list[str], list[str], int]:
    errors: list[str] = []
    warns: list[str] = []
    total = 0
    split_done = all((L.frag_dir(root, k) / L.TEMPLATE).exists() for k in L.KINDS)
    for kind in L.KINDS:
        frags = []
        d = L.frag_dir(root, kind)
        for p in sorted(d.glob("*.md")) if d.is_dir() else []:
            if p.name.startswith("_"):
                continue
            try:
                meta, body = L.parse_fragment(p.read_text(encoding="utf-8"))
            except (L.LedgerError, UnicodeDecodeError) as e:
                errors.append(f"{kind}/{p.name}: F1 {e}")
                continue
            frags.append(L.Fragment(kind, p.stem, meta, body, p))
        total += len(frags)
        for f in frags:
            check_schema(f, errors, warns)
            check_name(f, errors)
        check_unique(kind, frags, errors, warns)
        if kind == "implementations":
            check_classtags_of(root, frags, errors, warns)
        check_banner(root, frags, errors)
    if split_done:
        check_templates_and_build(root, errors)
    check_stubs(root, errors)
    return errors, warns, total


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Ledger-fragment gate (WP-161).")
    ap.add_argument("--root", type=Path, default=L.ROOT)
    ap.add_argument("--quiet", action="store_true", help="do not list warnings")
    args = ap.parse_args(argv)
    errors, warns, total = run(args.root.resolve())
    if not args.quiet:
        for w in warns:
            print("WARN  " + w)
    for e in errors:
        print("ERROR " + e)
    status = "FAIL" if errors else "OK"
    print(f"check_ledger_fragments: {status} ({total} fragments, {len(errors)} errors, {len(warns)} warns)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
