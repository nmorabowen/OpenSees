"""Self-test for ci/ledger.py and ci/check_ledger_fragments.py (WP-161).

The round-trip proof: build(split(ledger)) has the same entries and the same
template as the ledger, on synthetic ledgers that carry every shape the real
ones have (union-merge artefacts included) AND on the real ledgers -- the
working-tree files before the split, the split source the stubs record after
it. Then the gate's rules, `migrate`, and the two gates retargeted to the
fragments (quirk lint L3, viewer gate V1).
Run: pytest -q ci/test_ledger.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_ledger_fragments as gate  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402
import check_viewer_ledger as cvl  # noqa: E402
import ledger as L  # noqa: E402

FM = "---\ntitle: t\nproject: Ladruno\n---\n\n"

IMPL = FM + """# Ledger — implementations

Intro.

## Conventions

- **One row per feature**.

## Ledger

| Feature | Kind | Class tag | Files | Status | PR(s) |
|---|---|---|---|---|---|
| **WP-151 — a re-seat fix** (a `x \\| y` escaped pipe) | material | — | `a.cpp` | shipped | [#893](https://github.com/o/r/pull/893) |
| **LadrunoBrick** — brick | Element | `ELE_TAG_LadrunoBrick` 33002 | `b.cpp` | shipped | #65 |

| **ADR-97 P1 — closest point** | material | — | `c.cpp` | shipped | #800 |

## Documentation / ADR PRs + shipped-feature build history

Prose about the bullets.

- [#4](https://github.com/o/r/pull/4) — ADR: explicit integrators
- **Stage 1 SHIPPED:** a long bullet
  with an indented continuation line
| **WP-133 — an orphan row past the bullets** | material | — | `d.cpp` | shipped | #866 |

| **WP-133 — a second orphan after a blank** | material | — | `e.cpp` | shipped | #866 |
"""

QUIRKS = FM + """# Ledger — quirks

## Conventions

- One section per quirk.

## Quirks

### First quirk (WP-155, 2026-09-29)
- **Bites:** something.

```bash
## not a heading: inside a fence
### nor this
```

## Quirk: a level-two entry (ADR 46 P1)

Body.

### A sub-entry under it, with no id

Text.
---
"""

VANILLA = FM + """# Ledger — vanilla

## Conventions

- One row per (file, PR).

## Ledger

| Vanilla file | Why touched | PR |
|---|---|---|
| `SRC/a.cpp` | `// Ladruno WP-158` fix | [#901](https://github.com/o/r/pull/901) |
| `SRC/b.cpp` | other | [#754](https://github.com/o/r/pull/754) |

| `SRC/a.cpp` | `// Ladruno WP-158` second edit | [#901](https://github.com/o/r/pull/901) |
| `SRC/c.cpp` | a row with no PR cell |

> [!note] Upstreamable bugfixes
> Track those below.

| Vanilla file | Upstreamable fix | Fork PR |
|---|---|---|
| `SRC/d.cpp` | fix | [#7](https://github.com/o/r/pull/7) |
| `SRC/a.cpp` | appended at the end of the file | 793 |
"""

SYNTH = {"implementations": IMPL, "quirks": QUIRKS, "vanilla": VANILLA}


# ---------------------------------------------------------------- round trip
@pytest.mark.parametrize("kind", sorted(SYNTH))
def test_roundtrip_synthetic(kind):
    assert L.roundtrip(kind, SYNTH[kind]) == []


def test_synthetic_shapes_are_parsed_as_intended():
    _, impl = L.parse_ledger("implementations", IMPL)
    assert [e.section for e in impl] == ["table"] * 3 + ["history"] * 2 + ["table"] * 2
    assert "indented continuation" in impl[4].text
    _, q = L.parse_ledger("quirks", QUIRKS)
    assert [e.text.split("\n")[0] for e in q] == [
        "### First quirk (WP-155, 2026-09-29)", "## Quirk: a level-two entry (ADR 46 P1)",
        "### A sub-entry under it, with no id"]
    assert "## not a heading" in q[0].text
    tpl, v = L.parse_ledger("vanilla", VANILLA)
    assert [e.section for e in v] == ["main"] * 4 + ["upstreamable"] * 2
    assert "<!-- ledger:rows main -->" in tpl and "<!-- ledger:rows upstreamable -->" in tpl


def test_split_keys_names_and_vanilla_grouping():
    _, v = L.parse_ledger("vanilla", VANILLA)
    frags = L.split_entries("vanilla", v, set())
    by_key = {f.meta["wp"]: f for f in frags}
    wp158 = by_key["WP-158"]                 # one PR's rows, keyed by the WP the rows name
    assert wp158.body.count("\n") == 1 and wp158.meta["legacy_seq"] == [1, 3]
    assert wp158.meta["files"] == ["`SRC/a.cpp`", "`SRC/a.cpp`"]
    assert "PR-793" in by_key and by_key["PR-793"].meta["table"] == "upstreamable"
    for f in frags:
        assert L.NAME_RE.match(f.name) and f.name.startswith(f.meta["wp"] + "-"), f.name
    _, q = L.parse_ledger("quirks", QUIRKS)
    keys = [f.meta["wp"] for f in L.split_entries("quirks", q, set())]
    assert keys == ["WP-155", "ADR-46", "LEGACY"]


def test_build_groups_vanilla_rows_by_file_and_puts_new_entries_first(tmp_path):
    root = _split_tree(tmp_path)
    _frag(root, "vanilla", "WP-200-x", {"wp": "WP-200", "title": "x", "date": "2026-10-02"},
          "| `SRC/b.cpp` | `// Ladruno WP-200` | #999 |")
    _frag(root, "quirks", "WP-200-new-quirk", {"wp": "WP-200", "title": "New", "date": "2026-10-02"},
          "### New quirk (WP-200)\n- **Bites:** x.")
    _frag(root, "implementations", "WP-200-new-row",
          {"wp": "WP-200", "title": "New", "date": "2026-10-02", "status": "draft"},
          "| **WP-200 — new** | tool | — | `x.py` | draft | — |")
    out = L.build(root, tmp_path / "out")
    van = out["vanilla"].read_text(encoding="utf-8")
    rows = [ln for ln in van.split("\n") if ln.startswith("| `SRC/")]
    # a.cpp rows together, then b.cpp with its new row, then c.cpp
    assert [r.split("|")[1].strip() for r in rows[:4]] == ["`SRC/a.cpp`", "`SRC/a.cpp`", "`SRC/b.cpp`", "`SRC/b.cpp`"]
    assert "WP-200" in rows[3]
    q = out["quirks"].read_text(encoding="utf-8")
    assert q.index("### New quirk") < q.index("### First quirk")
    impl = out["implementations"].read_text(encoding="utf-8")
    assert impl.index("WP-200 — new") < impl.index("WP-151")
    assert L.GENERATED_RE.match(van.split("\n")[4])


@pytest.mark.parametrize("kind", sorted(L.KINDS))
def test_roundtrip_real_ledgers(kind):
    """THE proof on the fork's own ledgers: the working-tree file before the split,
    the split source the stub records after it."""
    try:
        text, src = L.ledger_source(L.ROOT, kind)
    except (L.LedgerError, OSError) as e:
        pytest.skip(f"real {kind} ledger not available here: {e}")
    assert L.roundtrip(kind, text) == [], f"round trip of {kind} from {src}"


# ---------------------------------------------------------------- the committed tree
def test_committed_fragments_pass_the_gate():
    errors, _, _ = gate.run(L.ROOT)
    assert errors == []


# ---------------------------------------------------------------- gate rules
def _split_tree(tmp_path):
    root = tmp_path / "repo"
    (root / L.IMPL_DIR).mkdir(parents=True)
    for kind, text in SYNTH.items():
        (root / L.IMPL_DIR / L.KINDS[kind]["ledger"]).write_text(text, encoding="utf-8")
    (root / "SRC").mkdir()
    (root / "SRC" / "classTags.h").write_text("#define ELE_TAG_LadrunoBrick 33002 // Ladruno\n", encoding="utf-8")
    assert L.cmd_split(root, None) == 0
    return root


def _frag(root, kind, name, meta, body):
    p = L.frag_dir(root, kind) / f"{name}.md"
    p.write_text(L.dump_fragment(meta, body), encoding="utf-8")
    return p


def _errors(root):
    return gate.run(root)[0]


def test_split_tree_is_clean_and_writes_stubs(tmp_path):
    root = _split_tree(tmp_path)
    assert _errors(root) == []
    stub = (root / L.IMPL_DIR / "LEDGER_quirks.md").read_text(encoding="utf-8")
    assert L.STUB_MARK in stub and len(stub.splitlines()) < 20
    with pytest.raises(L.LedgerError, match="already a stub"):
        L.cmd_split(root, None)


def test_split_keeps_a_wps_own_fragment_and_regenerates_legacy(tmp_path):
    root = _split_tree(tmp_path)
    own = _frag(root, "quirks", "WP-161-own", {"wp": "WP-161", "title": "own", "date": "2026-10-02"}, "### own")
    for kind, text in SYNTH.items():         # restore the full ledgers, as at merge time
        (root / L.IMPL_DIR / L.KINDS[kind]["ledger"]).write_text(text, encoding="utf-8")
    assert L.cmd_split(root, None) == 0
    assert own.exists() and _errors(root) == []


@pytest.mark.parametrize("name,meta,body,needle", [
    ("WP-200-x", {"wp": "WP-200", "title": "x"}, "### x", "needs `date"),
    ("WP-200-x", {"wp": "WP-201", "title": "x", "date": "2026-10-02"}, "### x", "does not match wp"),
    ("WP-200-Bad_Name", {"wp": "WP-200", "title": "x", "date": "2026-10-02"}, "### x", "F2"),
    ("WP-200-x", {"wp": "WP-200", "title": "x", "date": "2026-10-02", "colour": "red"}, "### x", "unknown"),
    ("WP-200-x", {"wp": "WP-200", "title": "x", "date": "2026-10-02"}, "no heading", "heading"),
    ("WP-200-x", {"wp": "WP-200", "title": "x", "date": "2026-10-02"}, "### a\n\n### b", "second"),
    ("WP-200-x", {"wp": "WP-200", "title": "x", "date": "2026-10-02", "banner": "not a line"}, "### x", "F5"),
])
def test_gate_flags_bad_quirk_fragments(tmp_path, name, meta, body, needle):
    root = _split_tree(tmp_path)
    _frag(root, "quirks", name, meta, body)
    errs = _errors(root)
    assert len(errs) == 1 and needle in errs[0], errs


def test_gate_flags_a_duplicated_entry(tmp_path):
    root = _split_tree(tmp_path)
    first = sorted(L.frag_dir(root, "quirks").glob("WP-155-*.md"))[0]
    meta, body = L.parse_fragment(first.read_text(encoding="utf-8"))
    _frag(root, "quirks", "WP-200-copy", {"wp": "WP-200", "title": "copy", "date": "2026-10-02"}, body)
    errs = _errors(root)
    assert errs and all("F3 duplicate" in e for e in errs), errs


def test_gate_flags_an_edited_stub(tmp_path):
    # the shape a stale branch produces when it union-merges a row into the stub
    root = _split_tree(tmp_path)
    p = root / L.IMPL_DIR / "LEDGER_vanilla_files.md"
    p.write_text(p.read_text(encoding="utf-8") + "| `SRC/x.cpp` | late row | #1 |\n", encoding="utf-8")
    errs = _errors(root)
    assert len(errs) == 1 and "F7" in errs[0], errs


def test_gate_checks_class_tags_against_classtags_h(tmp_path):
    root = _split_tree(tmp_path)
    row = "| **WP-200 — e** | Element | x | `e.cpp` | shipped | #1 |"
    base = {"wp": "WP-200", "title": "e", "date": "2026-10-02", "status": "shipped"}
    _frag(root, "implementations", "WP-200-e", {**base, "class_tags": ["ELE_TAG_LadrunoBrick=33002"]}, row)
    assert _errors(root) == []
    _frag(root, "implementations", "WP-200-e", {**base, "class_tags": ["ELE_TAG_LadrunoBrick=33003"]}, row)
    assert any("F4" in e and "33003" in e for e in _errors(root))
    _frag(root, "implementations", "WP-200-e", {**base, "class_tags": ["ELE_TAG_Nope"]}, row)
    assert any("F4" in e and "ELE_TAG_Nope" in e for e in _errors(root))
    _frag(root, "implementations", "WP-200-e",
          {**base, "status": "PLANNED", "class_tags": ["ELE_TAG_Nope"]}, row.replace("shipped", "PLANNED"))
    assert _errors(root) == []
    _frag(root, "implementations", "WP-201-f", {**base, "wp": "WP-201", "class_tags": ["ELE_TAG_LadrunoBrick"]},
          row.replace("WP-200", "WP-201"))
    _frag(root, "implementations", "WP-200-e", {**base, "class_tags": ["ELE_TAG_LadrunoBrick"]}, row)
    assert any("claimed by" in e for e in _errors(root))


def test_front_matter_round_trips_and_is_strict():
    meta = {"wp": "WP-161", "title": 'a "quoted": title', "date": "2026-10-02", "files": ["`a`", "b"],
            "legacy_seq": [3, 4]}
    assert L.parse_fragment(L.dump_fragment(meta, "body")) == (meta, "body")
    with pytest.raises(L.LedgerError):
        L.parse_fragment("---\ntitle: a: b\n---\nx")
    with pytest.raises(L.LedgerError):
        L.parse_fragment("no front matter")


# ---------------------------------------------------------------- migrate
def test_migrate_turns_a_branchs_ledger_edits_into_fragments(tmp_path):
    root = _split_tree(tmp_path)
    head = {
        "implementations": IMPL.replace(
            "|---|---|---|---|---|---|\n",
            "|---|---|---|---|---|---|\n| **WP-300 — new feature** | tool | — | `n.py` | draft | — |\n"
        ).replace("| shipped | #65 |", "| shipped, fixed | #65 |"),
        "quirks": QUIRKS + "\n### A new trap (WP-300)\n- **Bites:** y.\n",
        "vanilla": VANILLA.replace(
            "| `SRC/c.cpp` | a row with no PR cell |\n",
            "| `SRC/c.cpp` | a row with no PR cell |\n| `SRC/z.cpp` | `// Ladruno WP-300` | #300 |\n"),
    }
    made = {}
    for kind in L.KINDS:
        frags = L.load_fragments(root, kind)
        new, edits, notes = L.plan_migration(kind, SYNTH[kind], head[kind], frags, "WP-300", "2026-10-02",
                                             {f.name.lower() for f in frags})
        assert notes == [], notes
        made[kind] = (new, edits)
    new, edits = made["implementations"]
    assert [f.meta["wp"] for f in new] == ["WP-300"] and new[0].meta["status"] == "draft"
    assert len(edits) == 1 and "shipped, fixed" in edits[0][1] and edits[0][0].legacy
    new, edits = made["quirks"]
    assert len(new) == 1 and new[0].body.startswith("### A new trap") and not edits
    new, edits = made["vanilla"]
    assert len(new) == 1 and new[0].body == "| `SRC/z.cpp` | `// Ladruno WP-300` | #300 |"
    assert new[0].meta["pr"] == "#300"


# ---------------------------------------------------------------- retargeted gates
def test_l3_pointers_resolve_against_fragments(tmp_path):
    root = _split_tree(tmp_path)
    g = root / ".claude" / "skills" / "g" / "SKILL.md"
    g.parent.mkdir(parents=True)
    g.write_text('Quirks: "A sub-entry under it", "a heading nobody wrote".\n', encoding="utf-8")
    out = cq.check_pointers(root, lambda p: p.as_posix())
    assert len(out) == 1 and "a heading nobody wrote" in out[0]


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
                           "-c", "commit.gpgsign=false", "-c", "core.autocrlf=false", *args],
                          stdin=subprocess.DEVNULL, capture_output=True, text=True, check=True).stdout.strip()


def test_viewer_gate_accepts_a_fragment_line(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git not on PATH")
    _git(tmp_path, "init", "-q")
    (tmp_path / "Ladruno_tools" / "profiler_viewer").mkdir(parents=True)
    (tmp_path / "Ladruno_tools" / "profiler_viewer" / "a.py").write_text("x\n")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "base")
    base = _git(tmp_path, "rev-parse", "HEAD")
    (tmp_path / "Ladruno_tools" / "profiler_viewer" / "b.py").write_text("y\n")
    frag = tmp_path / cvl.FRAGMENTS / "WP-300-profiler.md"
    frag.parent.mkdir(parents=True)
    frag.write_text("---\nwp: WP-300\n---\n| `Ladruno_tools/profiler_viewer` gains b.py | tool |\n")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "tool + fragment")
    changes, added, removed = cvl.collect(tmp_path, base, "HEAD")
    assert cvl.find_violations(changes, added, removed) == []
