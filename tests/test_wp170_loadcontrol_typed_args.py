"""WP-170 -- `integrator LadrunoLoadControl` reads typed openseespy args.

openseespy passes TYPED args. `OPS_GetString()` answers the placeholder
"Invalid String Input!" for an int/float, so a parser that peeks a slot which
may be a number must use `OPS_GetStringFromAll` (LEDGER quirk "openseespy
parsers must peek a maybe-numeric arg with OPS_GetStringFromAll").

`OPS_LadrunoLoadControl` peeked the optional `numIter minLambda maxLambda`
triple with `OPS_GetString` and tested `nxt[0] != '-'`. That was correct only
by accident -- the placeholder's first letter is 'I', not '-' -- so the triple
was read, but a stray number among the flags was reported as
`unknown option 'Invalid String Input!'` instead of by its value. WP-170 moves
both reads to `OPS_GetStringFromAll` and decides "is this the triple" by
whether the token parses as a number, the same rule in both interpreters.

HOW THE TRIPLE IS OBSERVED. With lambda=0.05, numIter=4, min=0.01, max=0.1 on
a linear model (1 iteration per step), LoadControl grows the increment to
min(0.05*4/1, 0.1) = 0.1 after step 1, so two steps reach t=0.15. If the
triple is dropped, min=max=lambda and two steps reach t=0.10. The stock
`LoadControl` row pins the arithmetic, so the oracle is not this file's own
guess.

Measured before the fix (2026-10-02 build): every triple case below already
passed; only `test_stray_number_is_named_by_its_value` failed.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

LAM, NUM_ITER, LMIN, LMAX = 0.05, 4, 0.01, 0.1
T_TRIPLE = 0.15      # increment grew to LMAX after step 1
T_FIXED = 0.10       # increment stayed at LAM
TRIPLE = (LAM, NUM_ITER, LMIN, LMAX)


def _two_steps(*integrator_args):
    """1-DOF truss with a driven end, linear. Returns (oks, time)."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 30, 0)
    ops.algorithm("Newton")
    ops.integrator(*integrator_args)
    ops.analysis("Static")
    oks = [ops.analyze(1) for _ in range(2)]
    return oks, ops.getTime()


def test_premise_stock_loadcontrol_reads_the_triple():
    """Pins T_TRIPLE on the stock integrator; without it the rest is vacuous."""
    oks, t = _two_steps("LoadControl", *TRIPLE)
    assert oks == [0, 0]
    assert t == pytest.approx(T_TRIPLE, rel=1e-12)


# (flags, expected tangentPredictor, expected extrapolate)
FLAGS = [
    pytest.param((), 0.0, 0.0, id="no-flag"),
    pytest.param(("-tangentPredictor",), 1.0, 0.0, id="tangentPredictor"),
    pytest.param(("-extrapolate", 0.5), 0.0, 0.5, id="extrapolate-float"),
    pytest.param(("-extrapolate", 1), 0.0, 1.0, id="extrapolate-int"),
]


@pytest.mark.parametrize("flags,tp,ex", FLAGS)
def test_numeric_triple_plus_flag(flags, tp, ex):
    """The typed triple is read AND the flag after it lands."""
    oks, t = _two_steps("LadrunoLoadControl", *TRIPLE, *flags)
    assert oks == [0, 0]
    assert t == pytest.approx(T_TRIPLE, rel=1e-12), "numIter/min/max triple was not read"
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(tp)
    assert ops.ladrunoLoadControl("extrapolate") == pytest.approx(ex)


@pytest.mark.parametrize("flags,tp,ex", FLAGS[1:])
def test_flag_without_triple_is_not_taken_for_the_triple(flags, tp, ex):
    """`lambda -flag ...` must not swallow the flag (or its value) as numIter."""
    oks, t = _two_steps("LadrunoLoadControl", LAM, *flags)
    assert oks == [0, 0]
    assert t == pytest.approx(T_FIXED, rel=1e-12)
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(tp)
    assert ops.ladrunoLoadControl("extrapolate") == pytest.approx(ex)


def test_stray_number_is_named_by_its_value(capfd):
    """An unknown numeric token is reported by its value, not the placeholder.

    This is the case that failed before WP-170: the warning read
    `unknown option 'Invalid String Input!' ignored`.
    """
    capfd.readouterr()
    oks, t = _two_steps("LadrunoLoadControl", LAM, "-tangentPredictor", 7)
    err = capfd.readouterr().err
    assert oks == [0, 0]
    assert t == pytest.approx(T_FIXED, rel=1e-12)
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(1.0)
    assert "unknown option '7'" in err, err
    assert "Invalid String Input" not in err, err


# --- classic Tcl: the new "is it a number" peek must read the same there ------
# Only the step time and the parser's warnings are observable: the
# `ladrunoLoadControl` query command is not registered in classic Tcl
# (LadrunoCommandTable.h: LADRUNO_NONE in the Tcl slot).

REPO = Path(__file__).resolve().parents[1]


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()

TCL_DECK = """
proc run {args} {
    wipe
    model basic -ndm 1 -ndf 1
    node 1 0.0; node 2 1.0; fix 1 1
    uniaxialMaterial Elastic 1 1.0
    element Truss 1 1 2 1.0 1
    timeSeries Linear 1
    pattern Plain 1 1 { sp 2 1 1.0 }
    constraints Transformation; numberer Plain; system BandGeneral
    test NormDispIncr 1e-10 30 0; algorithm Newton
    eval integrator LadrunoLoadControl $args
    analysis Static
    analyze 1; analyze 1
    puts "RESULT [getTime]"
}
run 0.05 4 0.01 0.1
run 0.05 4 0.01 0.1 -tangentPredictor
run 0.05 4 0.01 0.1 -extrapolate 0.5
run 0.05 -tangentPredictor
run 0.05 -extrapolate 0.5
run 0.05 -tangentPredictor 7
puts "DECK-DONE"
"""


@pytest.mark.skipif(TCL_EXE is None,
                    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
def test_classic_tcl_parity(tmp_path):
    deck = tmp_path / "wp170.tcl"
    deck.write_text(TCL_DECK)
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see tests/test_wp103_getstringfromall_tcl.py (WinError 6).
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=120, env=env, cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert "DECK-DONE" in out, out
    times = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("RESULT ")]
    expected = [T_TRIPLE, T_TRIPLE, T_TRIPLE, T_FIXED, T_FIXED, T_FIXED]
    assert times == pytest.approx(expected, rel=1e-12), out
    # a flag mistaken for the triple would trip the int/double readers
    assert "failed to read" not in out, out
    assert out.count("unknown option") == 1 and "unknown option '7'" in out, out
