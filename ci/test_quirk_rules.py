"""Self-test for the quirk-lint rule registry (WP-162): slugs, aliases, the CLI,
and the ci/README.md rule table generated from RULES.
Run: pytest -q ci/test_quirk_rules.py
"""
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402

README = Path(__file__).resolve().parent / "README.md"
START, END = "<!-- quirk-rules:start (generated: python ci/check_quirk_patterns.py --rules-table) -->", \
    "<!-- quirk-rules:end -->"


def test_every_check_function_is_registered_once():
    checks = {n for n in dir(cq) if n.startswith("check_") and callable(getattr(cq, n))}
    assert checks - {"check_stale_waivers"} == {r.check.__name__ for r in cq._RULE_LIST}
    assert len(cq.RULES) == len(cq._RULE_LIST) == len(cq.ALIASES)


def test_slugs_are_kebab_case_and_aliases_are_the_old_numbers():
    for r in cq._RULE_LIST:
        assert re.fullmatch(r"[a-z][a-z0-9]*(?:-[a-z0-9]+)*", r.slug), r.slug
        assert re.fullmatch(r"L\d+", r.alias), r.alias


def test_only_takes_slugs_and_deprecated_aliases():
    assert cq.resolve_only("rayleigh,wipe") == (["rayleigh", "wipe"], [])
    slugs, notes = cq.resolve_only("L9, dead-decl ,L1")
    assert slugs == ["dead-decl", "rayleigh"] and len(notes) == 2 and "deprecated" in notes[0]
    with pytest.raises(ValueError, match="unknown rule 'L11'"):
        cq.resolve_only("L11")


DEAD_DECL = ("void M::f(const Vector& s)\n{\n    Vector r(6);\n    if (p > small)\n"
             "        Vector r = s / p;\n}\n")


def test_cli_labels_findings_by_slug(tmp_path, capsys):
    root = kit.tree(tmp_path, {"SRC/material/M.cpp": DEAD_DECL})
    assert cq.main(["--root", str(root), "--only", "dead-decl"]) == 1
    out = capsys.readouterr().out.splitlines()
    assert out[0].startswith("dead-decl SRC/material/M.cpp:4:"), out
    assert out[-1].endswith("[dead-decl]")


def test_cli_alias_still_works_and_warns(tmp_path, capsys):
    root = kit.tree(tmp_path, {"SRC/material/M.cpp": DEAD_DECL})
    assert cq.main(["--root", str(root), "--only", "L9"]) == 1
    cap = capsys.readouterr()
    assert "deprecated alias" in cap.err and cap.out.startswith("dead-decl ")


def test_cli_unknown_rule_is_exit_2(tmp_path, capsys):
    assert cq.main(["--root", str(tmp_path), "--only", "nope"]) == 2


def test_readme_rule_table_is_generated_from_RULES():
    text = README.read_text(encoding="utf-8").replace("\r\n", "\n")
    assert START in text and END in text, "ci/README.md lost its generated rule-table markers"
    table = text.split(START, 1)[1].split(END, 1)[0].strip("\n")
    assert table == cq.rules_table(), \
        "ci/README.md rule table is stale: paste the output of `python ci/check_quirk_patterns.py --rules-table`"
