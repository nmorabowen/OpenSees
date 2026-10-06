"""WP-170 -- fork integrator parsers read typed openseespy args as text safely.

openseespy passes TYPED args. `OPS_GetString()` answers the placeholder
"Invalid String Input!" for an int/float, so a parser that reads a slot which
may be a number must use `OPS_GetStringFromAll` (LEDGER quirk "openseespy
parsers must peek a maybe-numeric arg with OPS_GetStringFromAll").

Four fork parsers peeked an optional value with `OPS_GetString` and decided
"value or next flag" by `peek[0] == '-'`:

  * `LadrunoLoadControl`   -- the `numIter minLambda maxLambda` triple
  * `CentralDifferenceSMS`, `CentralDifferenceSMSConsistent`,
    `ExplicitBathe{SMS,SMSConsistent,LNVDSMS,LNVDSMSConsistent}`
                           -- the downgraded `-recompute N`
  * `ExplicitBathe -lnvd [alpha]`

Under openseespy that worked only by accident (the placeholder starts with
'I'); under classic Tcl a NEGATIVE value starts with '-' and was taken for a
flag. The flag loops also printed the placeholder for a stray number. WP-170
reads every token with `OPS_GetStringFromAll` and classifies by "does it parse
as a number".

HOW THE LoadControl TRIPLE IS OBSERVED. With lambda=0.05, numIter=4,
min=0.01, max=0.1 on a linear model (1 iteration per step), LoadControl grows
the increment to min(0.05*4/1, 0.1) = 0.1 after step 1, so two steps reach
t=0.15. If the triple is dropped, min=max=lambda and two steps reach t=0.10.
The stock `LoadControl` row pins the arithmetic.

Measured on the pre-fix build (2026-10-02): every "already worked" case passes;
the stray-number / non-numeric-value cases and the Tcl negative-value cases
fail. Each test says which group it is in.
"""
import os
import subprocess
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

PLACEHOLDER = "Invalid String Input"


# =========================================================================
# LadrunoLoadControl
# =========================================================================
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
def test_loadcontrol_numeric_triple_plus_flag(flags, tp, ex):
    """[already worked] The typed triple is read AND the flag after it lands."""
    oks, t = _two_steps("LadrunoLoadControl", *TRIPLE, *flags)
    assert oks == [0, 0]
    assert t == pytest.approx(T_TRIPLE, rel=1e-12), "numIter/min/max triple was not read"
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(tp)
    assert ops.ladrunoLoadControl("extrapolate") == pytest.approx(ex)


@pytest.mark.parametrize("flags,tp,ex", FLAGS[1:])
def test_loadcontrol_flag_without_triple(flags, tp, ex):
    """[already worked] `lambda -flag ...` must not swallow the flag as numIter."""
    oks, t = _two_steps("LadrunoLoadControl", LAM, *flags)
    assert oks == [0, 0]
    assert t == pytest.approx(T_FIXED, rel=1e-12)
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(tp)
    assert ops.ladrunoLoadControl("extrapolate") == pytest.approx(ex)


def test_loadcontrol_stray_number_named_by_value(capfd):
    """[failed pre-fix] warned `unknown option 'Invalid String Input!'`."""
    capfd.readouterr()
    oks, t = _two_steps("LadrunoLoadControl", LAM, "-tangentPredictor", 7)
    err = capfd.readouterr().err
    assert oks == [0, 0]
    assert t == pytest.approx(T_FIXED, rel=1e-12)
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(1.0)
    assert "unknown option '7'" in err, err
    assert PLACEHOLDER not in err, err


# =========================================================================
# explicit integrators: -recompute N (SMS downgrade), -lnvd [alpha], stray tokens
# =========================================================================
SMS_PREFIX = {
    "cd-sms": ("CentralDifferenceSMS", 0.1),
    "cd-sms-consistent": ("CentralDifferenceSMSConsistent", 0.1),
    "eb-sms": ("ExplicitBatheSMS", 0.54, 0.1),
    "eb-sms-consistent": ("ExplicitBatheSMSConsistent", 0.54, 0.1),
    "eb-lnvd-sms": ("ExplicitBatheLNVDSMS", 0.54, 0.6, 0.1),
}
ALL_PREFIX = dict(SMS_PREFIX, **{
    "eb": ("ExplicitBathe", 0.54),
    "eb-lnvd-alias": ("ExplicitBatheLNVD", 0.54, 0.6),
})


def _integrator_stderr(args, capfd):
    """Issue `integrator *args` on an empty model. Returns (accepted, stderr).

    `accepted` is only meaningful as True-means-no-exception: a parser that
    returns null does NOT raise under openseespy (see the -lnvd tests).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    capfd.readouterr()
    try:
        ops.integrator(*args)
        accepted = True
    except Exception:          # a null parser return surfaces as an OpenSees error
        accepted = False
    return accepted, capfd.readouterr().err


@pytest.mark.parametrize("n", [5, 5.0, -1], ids=["int", "float", "negative"])
@pytest.mark.parametrize("key", list(SMS_PREFIX))
def test_recompute_numeric_value_consumed(key, n, capfd):
    """[already worked under openseespy] any typed number after -recompute is its N."""
    ok, err = _integrator_stderr(SMS_PREFIX[key] + ("-recompute", n), capfd)
    assert ok, err
    assert "downgraded to REPORT-ONLY" in err, err
    assert "unknown option" not in err, err


@pytest.mark.parametrize("key", list(SMS_PREFIX))
def test_recompute_non_numeric_word_is_not_swallowed(key, capfd):
    """[failed pre-fix] `-recompute abc`: the word was silently eaten as N."""
    ok, err = _integrator_stderr(SMS_PREFIX[key] + ("-recompute", "abc"), capfd)
    assert ok, err
    assert "unknown option" in err and "abc" in err, err


@pytest.mark.parametrize("key", list(ALL_PREFIX))
def test_stray_number_named_by_value(key, capfd):
    """[failed pre-fix] the unknown-option warning printed the placeholder."""
    ok, err = _integrator_stderr(ALL_PREFIX[key] + ("-verbose", 7), capfd)
    assert ok, err
    assert "unknown option" in err and "7" in err, err
    assert PLACEHOLDER not in err, err


@pytest.mark.parametrize("key", ["cd-sms", "eb-sms", "eb", "eb-lnvd-alias"])
def test_numeric_lump_named_by_value(key, capfd):
    """[failed pre-fix] `-lump 3` warned `unknown -lump Invalid String Input!`."""
    ok, err = _integrator_stderr(ALL_PREFIX[key] + ("-lump", 3), capfd)
    assert ok, err
    assert "unknown -lump 3" in err, err
    assert PLACEHOLDER not in err, err


@pytest.mark.parametrize("alpha", [0.6, 0, 0.0], ids=["float", "int-zero", "float-zero"])
def test_lnvd_typed_alpha_accepted(alpha, capfd):
    """[already worked] a typed alpha after -lnvd is read, not taken for a flag."""
    ok, err = _integrator_stderr(("ExplicitBathe", 0.54, "-lnvd", alpha, "-verbose"), capfd)
    assert ok, err
    assert "unknown option" not in err, err


def test_lnvd_without_alpha_before_flag(capfd):
    """[already worked] `-lnvd -verbose`: -verbose is the next flag, alpha defaults."""
    ok, err = _integrator_stderr(("ExplicitBathe", 0.54, "-lnvd", "-verbose"), capfd)
    assert ok, err
    assert "unknown option" not in err, err


def test_lnvd_negative_alpha_is_refused(capfd):
    """[already worked] A negative alpha is a value; the [0,1) check refuses it.

    Asserted on stderr, not on an exception: `OPS_Integrator` returns 0 even
    when the factory returns null, so openseespy does not raise here.
    """
    _, err = _integrator_stderr(("ExplicitBathe", 0.54, "-lnvd", -0.5), capfd)
    assert "alpha must be in [0,1)" in err, err


def test_lnvd_non_numeric_value_still_fatal(capfd):
    """[already worked] `-lnvd abc` stays a refusal, not an 'unknown option'."""
    _, err = _integrator_stderr(("ExplicitBathe", 0.54, "-lnvd", "abc"), capfd)
    assert "unknown option" not in err, err
    assert "need an alpha value" in err, err


# =========================================================================
# classic Tcl: same classification there. Only step time and warnings are
# observable for LoadControl -- `ladrunoLoadControl` is not registered in
# classic Tcl (LadrunoCommandTable.h: LADRUNO_NONE in the Tcl slot).
# =========================================================================
REPO = Path(__file__).resolve().parents[1]


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(
    TCL_EXE is None,
    reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")


def _run_tcl(tmp_path, deck_text):
    deck = tmp_path / "wp170.tcl"
    deck.write_text(deck_text)
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1")
    # stdin=DEVNULL: see tests/test_wp103_getstringfromall_tcl.py (WinError 6).
    proc = subprocess.run([TCL_EXE, str(deck)], capture_output=True, text=True,
                          timeout=120, env=env, cwd=str(tmp_path),
                          stdin=subprocess.DEVNULL)
    out = proc.stdout + proc.stderr
    assert "DECK-DONE" in out, out
    return out


LOADCONTROL_DECK = """
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


@needs_exe
def test_classic_tcl_loadcontrol(tmp_path):
    """[already worked] Tcl was all-strings; pins the new peek reads the same."""
    out = _run_tcl(tmp_path, LOADCONTROL_DECK)
    times = [float(ln.split()[1]) for ln in out.splitlines() if ln.startswith("RESULT ")]
    expected = [T_TRIPLE, T_TRIPLE, T_TRIPLE, T_FIXED, T_FIXED, T_FIXED]
    assert times == pytest.approx(expected, rel=1e-12), out
    # a flag mistaken for the triple would trip the int/double readers
    assert "failed to read" not in out, out
    assert out.count("unknown option") == 1 and "unknown option '7'" in out, out


EXPLICIT_DECK = """
proc try {label args} {
    wipe
    model basic -ndm 1 -ndf 1
    set rc [catch {eval integrator $args} msg]
    puts "CASE $label rc=$rc"
}
try recompute-neg   CentralDifferenceSMS 0.1 -recompute -1
try recompute-flag  CentralDifferenceSMS 0.1 -recompute -verbose
try lnvd-neg        ExplicitBathe 0.54 -lnvd -0.5
try lnvd-flag       ExplicitBathe 0.54 -lnvd -verbose
puts "DECK-DONE"
"""


@needs_exe
def test_classic_tcl_negative_values(tmp_path):
    """[failed pre-fix] Tcl took `-1` / `-0.5` for a flag.

    `-recompute -1` fell into `unknown option -1`; `-lnvd -0.5` was un-read,
    warned as an unknown option, and ran with the DEFAULT alpha 0.8 instead of
    refusing the out-of-range value.
    """
    out = _run_tcl(tmp_path, EXPLICIT_DECK)
    assert "unknown option -1" not in out, out
    assert "unknown option -0.5" not in out, out
    assert "alpha must be in [0,1)" in out, out
    # the flag-after-flag forms still un-read the flag
    assert "unknown option -verbose" not in out, out
