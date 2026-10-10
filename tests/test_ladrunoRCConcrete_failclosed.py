"""Fail-closed unknown-token policy for the LadrunoRCConcrete family (WP-167).

Before WP-167 both `OPS_LadrunoRCConcrete` and `OPS_LadrunoRCFiniteStrain` ended
their option ladder with "unknown tokens are ignored (forward-compat)". apeGmsh
#1184 emitted `-crackedNu 0.2` / `-betaC 170` for this material; the build had
neither, so both were DROPPED without a word and the model ran with the default
physics. A misspelled flag (`-tensStif vc`) went the same way.

The policy now (SRC/material/LadrunoOptSpec.h): every accepted option is
declared in `kLadrunoRCOptions` (name, value count, `since`); any other token
makes the command FAIL and the message names the material, its tag, the token
and the build stamp.

The oracle is the declared grammar itself:
  * every option the table declares still parses, for both materials (so the
    table and the ladder agree, and no accepted option changed meaning -- the
    behavioural suites test_ladrunoRCConcrete_*.py pin the meanings);
  * the table in the header equals the list below (a new option must be added
    to both, with a parse case);
  * `-crackedNu 0.2`, `-betaC 170` (refused by maintainer decision 2026-10-05
    until fork PR #877 declares them) and a misspelling FAIL, name the token,
    and create no material -- under openseespy AND the classic Tcl exe.

Revert-proof: on unmodified `ladruno` the refusal tests fail (the tokens are
accepted silently).
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
HEADER = REPO / "SRC" / "material" / "LadrunoOptSpec.h"
FAMILY = ("LadrunoRCConcrete", "LadrunoRCFiniteStrain")

E, NU = 30000.0, 0.2
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0, 30.0, 5.0]
CD = [0.0, 0.0, 0.25, 1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0, 0.5]
TD = [0.0, 0.0, 0.9]
BACKBONE = ["-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Te", *TE, "-Ts", *TS, "-Td", *TD]

# Every declared option -> the values a valid command passes after it. The six
# backbone lists are in BACKBONE (every command needs -Ce/-Cs/-Te/-Ts).
ACCEPTED = {
    "-Ce": None, "-Cs": None, "-Cd": None, "-Te": None, "-Ts": None, "-Td": None,
    "-Kc": [2.0 / 3.0], "-betaFloor": [0.1], "-rho": [2.4e-9],
    "-beta": [], "-lublinerReduced": [], "-secant": [], "-numericalTangent": [],
    "-interlock": [], "-agg": [16.0], "-crackStrain": [1.0e-4], "-crackSpacing": [100.0],
    "-lch": [50.0], "-betaSrMin": [0.01],
    "-cyclic": [],
    "-xcrack": [], "-degKappa": [0.5], "-degSlipRef": [0.01], "-degMin": [0.1],
    "-implex": [], "-implexAlpha": [1.0], "-implexControl": [0.05, 0.01],
    "-shearRetention": ["const"], "-shearRetFactor": [0.4],
    "-tensStiff": ["vc"], "-tensStiffC": [500.0], "-tensStiffAlpha": [1.0],
    "-autoRegularization": [50.0],
}

REFUSED = {
    "crackedNu": ["-crackedNu", 0.2],     # apeGmsh #1184, open fork PR #877
    "betaC": ["-betaC", 170.0],           # apeGmsh #1184, open fork PR #877
    "misspelled": ["-tensStif", "vc"],    # one letter short of -tensStiff
    "bare-word": ["bogus"],
}


def _declared_in_header():
    text = HEADER.read_text(encoding="utf-8")
    body = text.split("kLadrunoRCOptions[] = {", 1)[1].split("};", 1)[0]
    body = re.sub(r"//[^\n]*", "", body)
    return re.findall(r'\{"(-[A-Za-z]+)",\s*(LIST|\d+),\s*"[^"]+"\}', body)


def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


def _cmd(mat, tag, extra):
    return [mat, tag, E, NU, *BACKBONE, *extra]


def test_header_table_is_the_accepted_list():
    declared = _declared_in_header()
    names = [n for n, _ in declared]
    assert len(names) == len(set(names)), f"duplicate option in kLadrunoRCOptions: {names}"
    assert set(names) == set(ACCEPTED), (
        "kLadrunoRCOptions and this test disagree: add the option to both, with a parse case. "
        f"only in header: {sorted(set(names) - set(ACCEPTED))}; "
        f"only in test: {sorted(set(ACCEPTED) - set(names))}")
    for name, nargs in declared:
        vals = ACCEPTED[name]
        if vals is None:
            assert nargs == "LIST", name
        else:
            assert int(nargs) == len(vals), f"{name}: header says {nargs} value(s), test passes {vals}"
    assert "-crackedNu" not in names and "-betaC" not in names, \
        "-crackedNu/-betaC are refused until fork PR #877 declares them with its own `since`"


@pytest.mark.parametrize("mat", FAMILY)
@pytest.mark.parametrize("opt", sorted(o for o, v in ACCEPTED.items() if v is not None))
def test_every_declared_option_still_parses(mat, opt):
    _fresh()
    ops.nDMaterial(*_cmd(mat, 11, [opt, *ACCEPTED[opt]]))   # raises if refused


@pytest.mark.parametrize("mat", FAMILY)
def test_all_declared_options_together_parse(mat):
    _fresh()
    extra = [x for o, v in ACCEPTED.items() if v is not None for x in (o, *v)]
    ops.nDMaterial(*_cmd(mat, 12, extra))


@pytest.mark.parametrize("mat", FAMILY)
@pytest.mark.parametrize("case", sorted(REFUSED))
def test_undeclared_token_fails_and_names_it(mat, case, capfd):
    _fresh()
    tag = 7
    tok = REFUSED[case][0]
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.nDMaterial(*_cmd(mat, tag, ["-beta", *REFUSED[case]]))
    err = capfd.readouterr().err
    assert f"{mat} {tag}: unknown option '{tok}'" in err, err
    if hasattr(ops, "ladrunoBuild"):
        assert f"(build {ops.ladrunoBuild()})" in err, err
    # nothing was registered under the tag: the same tag takes a valid material
    ops.nDMaterial(*_cmd(mat, tag, ["-beta"]))


@pytest.mark.parametrize("mat", FAMILY)
@pytest.mark.parametrize("opt", ["-Kc", "-betaFloor", "-rho", "-implexControl", "-shearRetention"])
def test_declared_option_missing_its_values_fails(mat, opt, capfd):
    _fresh()
    vals = ACCEPTED[opt][:-1]               # one value short, at the end of the line
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.nDMaterial(*_cmd(mat, 9, [opt, *vals]))
    ops.nDMaterial(*_cmd(mat, 9, ["-beta"]))   # nothing registered under the tag


@pytest.mark.parametrize("mat", FAMILY)
@pytest.mark.parametrize("opt", ["-Kc", "-betaFloor", "-rho"])
def test_non_numeric_scalar_value_fails(mat, opt, capfd):
    # before WP-167 these three reads were unchecked: `-Kc abc` kept the default
    # (and `-rho abc` silently meant mass 0). Now the command fails and says why.
    _fresh()
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.nDMaterial(*_cmd(mat, 8, [opt, "abc", "-beta"]))
    assert f"{opt} needs a value" in capfd.readouterr().err
    ops.nDMaterial(*_cmd(mat, 8, ["-beta"]))


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()


@pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_refuses_the_same_tokens(tmp_path):
    bb = " ".join(str(x) for x in BACKBONE)
    lines = ["model basic -ndm 3 -ndf 3"]
    for mat in FAMILY:
        lines.append(f'puts "OK-{mat} [catch {{nDMaterial {mat} 1{len(lines)} {E} {NU} {bb} -beta}}]"')
        for case, toks in sorted(REFUSED.items()):
            t = " ".join(str(x) for x in toks)
            lines.append(f'puts "REF-{mat}-{case} [catch {{nDMaterial {mat} 2{len(lines)} {E} {NU} {bb} {t}}}]"')
    deck = tmp_path / "failclosed.tcl"
    deck.write_text("\n".join(lines) + "\n", encoding="ascii")
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True, timeout=300,
                          env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"), cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    for mat in FAMILY:
        assert f"OK-{mat} 0" in out, out          # a valid command still succeeds
        for case, toks in REFUSED.items():
            assert f"REF-{mat}-{case} 1" in out, out
            assert f"unknown option '{toks[0]}'" in out, out
