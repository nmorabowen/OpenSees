"""Self-test for quirk-lint rule `unknown-token` (WP-167) in ci/check_quirk_patterns.py.

The incident: OPS_LadrunoRCConcrete ended its option ladder with
`// unknown tokens are ignored (forward-compat)` and silently dropped apeGmsh's
-crackedNu / -betaC. One file per rule (WP-162).
Run: pytest -q ci/test_quirk_unknown_token.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402

PATH = "SRC/material/nD/M.cpp"

# the pre-WP-167 LadrunoRCConcrete shape
SILENT = kit.STAMP + """void* OPS_M(void)
{
  int tag = 1; double Kc = 0.0; bool b = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if      (strcmp(opt, "-Kc") == 0)   { int nd = 1; OPS_GetDoubleInput(&nd, &Kc); }
    else if (strcmp(opt, "-beta") == 0) b = true;
    // unknown tokens are ignored (forward-compat)
  }
  return 0;
}
"""
FIXED = SILENT.replace("    // unknown tokens are ignored (forward-compat)\n",
                       '    else { opserr << "unknown option"; return 0; }\n')


def _run(tmp_path, files):
    root = kit.tree(tmp_path, files)
    used = set()
    return cq.check_unknown_token(root, kit.rel(root), used) + cq.check_stale_waivers(root, kit.rel(root), used)


def test_flags_the_ladrunoRCConcrete_incident_twice(tmp_path):
    out = _run(tmp_path, {PATH: SILENT})
    assert len(out) == 2, out
    # kit.STAMP is two lines: the ladder starts on line 8, the comment is line 10
    assert any("comment" in f and ":10:" in f for f in out), out
    assert any("no final `else`" in f and ":8:" in f for f in out), out


def test_passes_the_fail_closed_ladder(tmp_path):
    assert _run(tmp_path, {PATH: FIXED}) == []


def test_ladder_without_else_flagged_even_without_the_comment(tmp_path):
    src = SILENT.replace("    // unknown tokens are ignored (forward-compat)\n", "")
    out = _run(tmp_path, {PATH: src})
    assert len(out) == 1 and "no final `else`" in out[0], out


def test_ladder_followed_by_a_report_is_not_silent(tmp_path):
    # the chain is not the loop's last statement: something handles the fall-through
    src = SILENT.replace("    // unknown tokens are ignored (forward-compat)\n",
                         '    opserr << "unknown"; return 0;\n')
    assert _run(tmp_path, {PATH: src}) == []


def test_comment_variants(tmp_path):
    for i, c in enumerate(["// unknown options are ignored (forward-compatible)",
                           "/* unrecognised flags are silently skipped */",
                           "// ignore unknown arguments"]):
        out = _run(tmp_path / str(i), {PATH: FIXED.replace("  return 0;\n}", "  " + c + "\n  return 0;\n}")})
        assert len(out) == 1 and "comment" in out[0], (c, out)


def test_message_strings_and_vanilla_are_out_of_scope(tmp_path):
    # a WARNING string is a warn policy, not a comment; an unstamped (vanilla) file is skipped
    warn = FIXED.replace('opserr << "unknown option"', 'opserr << "unknown option ignored"')
    assert _run(tmp_path / "a", {PATH: warn}) == []
    assert _run(tmp_path / "b", {PATH: SILENT.replace(kit.STAMP, "")}) == []


def test_waiver_and_stale_waiver(tmp_path):
    waived = SILENT.replace("    // unknown tokens are ignored (forward-compat)\n", "").replace(
        "    if      (strcmp", "    // ladruno-lint: unknown-ok the tail is re-parsed by the section\n    if      (strcmp")
    assert _run(tmp_path / "a", {PATH: waived}) == []
    short = waived.replace("the tail is re-parsed by the section", "ok")
    out = _run(tmp_path / "b", {PATH: short})
    assert len(out) == 1 and "too short" in out[0], out
    stale = FIXED.replace("    if      (strcmp",
                          "    // ladruno-lint: unknown-ok the tail is re-parsed by the section\n    if      (strcmp")
    out = _run(tmp_path / "c", {PATH: stale})
    assert len(out) == 1 and "stale unknown-ok" in out[0], out


def test_legacy_ratchet_only_shrinks(tmp_path, monkeypatch):
    monkeypatch.setattr(cq, "UNKNOWN_TOKEN_LEGACY", {PATH: "recorded"})
    assert _run(tmp_path / "a", {PATH: SILENT}) == []           # recorded offender: tolerated
    out = _run(tmp_path / "b", {PATH: FIXED})                   # converted: the entry must go
    assert len(out) == 1 and "delete the entry" in out[0], out


def test_real_tree_is_clean():
    root = Path(__file__).resolve().parent.parent
    assert cq.main(["--root", str(root), "--only", "unknown-token"]) == 0
