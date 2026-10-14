"""WP-193: cost-aware refactorization trigger for factorization reuse
(`algorithm KrylovNewton -reuseFactor m -reuseCost on`).

ADR: Ladruno_implementation/193_ladruno_reuse_cost_trigger_adr.md. Reuses the WP-188
frame deck (tests/wp188_factor_reuse/frame_gm.tcl) and its runner.

What is pinned:
  * `-reuseCost off` is the WP-188 path byte for byte (and no cost line is printed);
    with LADRUNO_WP193_BASE_DIST pointing at a dist/ built from the base commit, the
    WP-188 outputs (option absent, -reuseFactor 30) are md5-equal to it as well;
  * the trigger follows the measured cost ratio on the same model, with the regimes
    made by test-only hooks (the frame alone has t_F ~ t_I):
      - iterations expensive (LADRUNO_WP193_TEST_ITER_DELAY_MS = 1: t_F/t_I ~ 0.06,
        the Paper 1 1D-cell regime): it refactors at most steps, like full Newton;
      - factorization expensive (LADRUNO_WP193_TEST_FACTOR_DELAY_MS = 20: t_F/t_I
        ~ 200, the 4D regime): it keeps the factor, like WP-188 with a large m;
    and its wall time is within 1.2x (+50 ms) of the better of the two fixed
    policies run by the same driver with the same hooks;
  * results stay within the test-tolerance scale of full Newton;
  * the parsers fail closed (classic Tcl and openseespy);
  * OpenSeesMP on 2 ranks: every rank prints identical statistics, cost line included
    (the timings are max-reduced before they enter the model).
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

from test_wp188_factor_reuse import (DIST, MP_EXE, MPIEXEC, SER_EXE, _edp_diff, _hist_diff,
                                     _run as _run188)

pytestmark = [pytest.mark.zone_a]

needs_serial = pytest.mark.skipif(not SER_EXE.exists(), reason="dist OpenSees not built")
needs_mp = pytest.mark.skipif(not (MP_EXE.exists() and MPIEXEC.exists()),
                              reason="dist OpenSeesMP / mpiexec not built")

COST = re.compile(
    r"reuse cost model \(WP-193\) ([\w ]+?)(?: \[rank (\d+)\])?: t_factor ([-\d.e+]+) s, t_iter ([-\d.e+]+) s, "
    r"ratio ([-\d.e+]+), n_fresh ([-\d.e+]+), n_kept1 ([-\d.e+]+), slope ([-\d.e+]+), starts kept by cost (\d+), "
    r"starts formed by cost (\d+), in-step cost refreshes (\d+), probes (\d+)")
CKEYS = ("tF", "tI", "ratio", "nF", "nK1", "slope", "kept", "formed", "instep", "probes")

DELAY = "LADRUNO_WP193_TEST_FACTOR_DELAY_MS"
ITER_DELAY = "LADRUNO_WP193_TEST_ITER_DELAY_MS"


def _cost_lines(text):
    out = []
    for m in COST.finditer(text):
        d = {"when": m.group(1), "rank": None if m.group(2) is None else int(m.group(2))}
        for k, v in zip(CKEYS, m.groups()[2:]):
            d[k] = int(v) if k in ("kept", "formed", "instep", "probes") else float(v)
        out.append(d)
    return out


def _run(tmp, name, np=1, delay_ms=None, iter_ms=None, cost=None, **opts):
    extra = {}
    if delay_ms is not None:
        extra[DELAY] = str(delay_ms)
    if iter_ms is not None:
        extra[ITER_DELAY] = str(iter_ms)
    if cost is not None:
        opts["reusextra"] = cost
    r = _run188(tmp, name, np=np, extra_env=extra or None, **opts)
    r["cost"] = []
    for k in range(np):
        f = Path(tmp) / (name + ".algo%d" % k)
        c = _cost_lines(f.read_text(errors="replace")) if f.exists() else []
        r["cost"].append(c[-1] if c else None)
    wf = Path(tmp) / (name + ".wall")
    r["wall"] = float(wf.read_text().split()[0]) if wf.exists() else None
    return r


# frame, the first 300 steps of the record (ramp + strong motion), tol 1e-6, -maxDim 20
BASE = dict(tol=1e-6, maxdim=20, nstep=300)


@needs_serial
def test_cost_off_is_the_wp188_path(tmp_path):
    a = _run(tmp_path, "r30", reuse=30, **BASE)
    b = _run(tmp_path, "r30off", reuse=30, cost="-reuseCost off", **BASE)
    assert a["edp"]["fail"] == 0
    assert a["md5"] == b["md5"]
    assert a["algo"][0] == b["algo"][0]
    assert a["cost"][0] is None and b["cost"][0] is None   # no cost line when off
    base = os.environ.get("LADRUNO_WP193_BASE_DIST")
    if base:
        for opts in ({}, {"reuse": 30}, {"reuse": 30, "model": "elastic"}):
            o = dict(BASE, **opts)
            x = _run188(tmp_path, "b%d" % len(opts), dist=base, **o)
            y = _run188(tmp_path, "n%d" % len(opts), **o)
            assert x["md5"] == y["md5"], f"output differs from the base build for {opts}"
            assert x["algo"][0] == y["algo"][0]


@needs_serial
def test_linear_model_keeps_like_wp188(tmp_path):
    # n_K = n_F on a linear model: keeping pays, the factor lives (nearly) to m.
    # The fitted slope of n_K(age) is timing noise around 0, so a few cost
    # refactorizations may happen; the answer is full Newton's to round-off anyway.
    full = _run(tmp_path, "full", model="elastic", tol=1e-10, nstep=400)
    cost = _run(tmp_path, "cost", model="elastic", tol=1e-10, nstep=400, reuse=50, cost="-reuseCost on")
    assert _hist_diff(cost, full) < 1e-12
    s, c = cost["algo"][0], cost["cost"][0]
    assert s["pct"] > 90.0 and s["retries"] == 0, (s, c)
    assert c["instep"] == 0


# The two fixed policies, run by the same driver with the same hooks:
#   always refactor: the cost rule with an infinite safety factor (every step forms
#                    once the warm-up is over; no probes) -- full Newton in effect;
#   always keep:     WP-188 -reuseFactor 1000, no cost model.
REFACTOR = "-reuseCost on -reuseCostSafety 1e6 -reuseCostProbe 0"


@needs_serial
@pytest.mark.parametrize("regime", ["iterations_expensive", "factorization_expensive"])
def test_trigger_follows_the_cost_ratio(tmp_path, regime):
    """The same frame, two cost regimes made with the test-only hooks:
      iterations_expensive    every iteration burns 1 ms (an expensive state
                              determination; the frame's factorization is ~0.06 ms):
                              t_F/t_I ~ 0.05, the 1D-cell regime -> refactor;
      factorization_expensive every factorization burns 20 ms (a large model;
                              an iteration is ~0.07 ms): t_F/t_I ~ 300, the 4D
                              regime -> keep.
    The trigger must choose the matching policy and cost no more wall time than
    the better fixed policy (a 1.2x + 50 ms margin for wall-clock noise)."""
    hook = {"iter_ms": 1} if regime == "iterations_expensive" else {"delay_ms": 20}
    full = _run(tmp_path, "full", **BASE)
    ref = _run(tmp_path, "refactor", reuse=1000, cost=REFACTOR, **hook, **BASE)
    keep = _run(tmp_path, "keep", reuse=1000, **hook, **BASE)
    cost = _run(tmp_path, "cost", reuse=1000, cost="-reuseCost on", **hook, **BASE)
    for r in (ref, keep, cost):
        assert r["edp"]["fail"] == 0 and r["edp"]["steps"] == BASE["nstep"]
        # results: the test-tolerance scale of full Newton (as WP-188)
        assert max(_edp_diff(r, full), _hist_diff(r, full)) < 1e-4
    s, c = cost["algo"][0], cost["cost"][0]
    walls = {"refactor": ref["wall"], "keep": keep["wall"], "cost": cost["wall"]}
    print("WP-193", regime, walls, s, c)
    if regime == "iterations_expensive":
        assert c["ratio"] < 0.5
        assert ref["algo"][0]["pct"] < 5.0                  # the reference really refactors
        assert keep["wall"] > 1.4 * ref["wall"], walls      # and here that is the better policy
        assert s["pct"] < 40.0, (s, c)                      # the trigger mostly refactors
        assert c["formed"] + c["instep"] > 0.5 * BASE["nstep"]
        best = ref["wall"]
    else:
        assert c["ratio"] > 50.0
        assert ref["wall"] > 5.0 * keep["wall"], walls      # keeping is the better policy
        assert s["pct"] > 90.0, (s, c)                      # the trigger keeps
        assert s["fact"] <= 10
        best = keep["wall"]
    assert cost["wall"] <= 1.2 * best + 0.05, walls


@needs_serial
def test_no_hook_frame_converges_within_tolerance(tmp_path):
    """The frame as is (t_F ~ t_I): whatever the trigger decides, the results stay at
    the test-tolerance scale of full Newton and the cost line is complete."""
    full = _run(tmp_path, "full", **BASE)
    cost = _run(tmp_path, "cost", reuse=1000, cost="-reuseCost on", **BASE)
    assert cost["edp"]["fail"] == 0 and cost["edp"]["steps"] == BASE["nstep"]
    assert max(_edp_diff(cost, full), _hist_diff(cost, full)) < 1e-4
    c = cost["cost"][0]
    assert c["tF"] > 0 and c["tI"] > 0 and c["nF"] >= 1 and c["nK1"] >= 1
    assert c["kept"] + c["formed"] + c["probes"] > 0


_REFUSE_DECK = r"""
model BasicBuilder -ndm 1 -ndf 1
node 1 0.0
node 2 1.0
fix 1 1
uniaxialMaterial Elastic 1 100.0
element truss 1 1 2 1.0 1
test NormDispIncr 1e-8 10
set r [catch {algorithm KrylovNewton %s} msg]
puts "RESULT $r"
"""


@needs_serial
@pytest.mark.parametrize("args,refused", [
    ("-reuseFactor 5 -reuseCost on", 0),
    ("-reuseFactor 5 -reuseCost off", 0),
    ("-reuseFactor 5 -reuseCost on -reuseCostSafety 1.5 -reuseCostWarmup 2 -reuseCostProbe 0", 0),
    ("-reuseCost on", 1),
    ("-reuseFactor 5 -reuseCost maybe", 1),
    ("-reuseFactor 5 -reuseCost", 1),
    ("-reuseFactor 5 -reuseCostSafety 2", 1),
    ("-reuseFactor 5 -reuseCost off -reuseCostWarmup 2", 1),
    ("-reuseFactor 5 -reuseCost on -reuseCostSafety 0", 1),
    ("-reuseFactor 5 -reuseCost on -reuseCostWarmup -1", 1),
    ("-reuseFactor 5 -reuseCost on -reuseCostProbe x", 1),
])
def test_classic_parser_fails_closed(tmp_path, args, refused):
    deck = tmp_path / "refuse.tcl"
    deck.write_text(_REFUSE_DECK % args)
    env = dict(os.environ)
    tcl = DIST / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    p = subprocess.run([str(SER_EXE), str(deck)], stdin=subprocess.DEVNULL, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=60, env=env)
    m = re.search(r"RESULT (\d)", p.stdout + p.stderr)
    assert m, p.stdout + p.stderr
    assert int(m.group(1)) == refused, p.stdout + p.stderr


def test_openseespy_parser_and_linear_cost_reuse():
    """The interpreter path (OPS_KrylovNewton: openseespy / openseesmp), typed args."""
    from _testbed import ops

    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.node(2, 1.0)
        ops.fix(1, 1)
        ops.mass(2, 1.0)
        ops.uniaxialMaterial("Elastic", 1, 100.0)
        ops.element("truss", 1, 1, 2, 1.0, 1)
        ops.timeSeries("Sine", 1, 0.0, 10.0, 0.5)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0)
        ops.constraints("Plain"); ops.numberer("Plain"); ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-12, 20)

    def drive(*algo):
        build()
        ops.algorithm("KrylovNewton", *algo)
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        out = []
        for _ in range(50):
            assert ops.analyze(1, 0.01) == 0
            out.append(ops.nodeDisp(2, 1))
        return out

    full = drive()
    for args in (("-reuseFactor", 10, "-reuseCost", "on"),
                 ("-reuseFactor", 10, "-reuseCost", 1, "-reuseCostSafety", 2, "-reuseCostWarmup", 2),
                 ("-reuseFactor", 10, "-reuseCost", "off")):
        reuse = drive(*args)
        assert max(abs(a - b) for a, b in zip(full, reuse)) < 1e-14, args
    for bad in (("-reuseCost", "on"), ("-reuseFactor", 3, "-reuseCost", "yes"),
                ("-reuseFactor", 3, "-reuseCostSafety", 2.0),
                ("-reuseFactor", 3, "-reuseCost", "on", "-reuseCostSafety", -1.0)):
        build()
        with pytest.raises(Exception):
            ops.algorithm("KrylovNewton", *bad)
    ops.wipe()


@needs_mp
def test_mp_ranks_agree_on_the_cost_model(tmp_path):
    full = _run(tmp_path, "full", np=2, **BASE)
    cost = _run(tmp_path, "cost", np=2, reuse=1000, cost="-reuseCost on", delay_ms=20, **BASE)
    assert full["edp"]["fail"] == 0 and cost["edp"]["fail"] == 0
    a0, a1 = cost["algo"]
    c0, c1 = cost["cost"]
    assert a0 is not None and a1 is not None and c0 is not None and c1 is not None
    strip = lambda d: {k: v for k, v in d.items() if k not in ("rank", "when")}
    assert strip(a0) == strip(a1), (a0, a1)
    # the timings are max-reduced before the model sees them: identical on every rank
    assert strip(c0) == strip(c1), (c0, c1)
    assert a0["pct"] > 80.0
    assert _hist_diff(cost, full) < 1e-4
