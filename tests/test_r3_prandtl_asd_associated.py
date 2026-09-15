"""R3 associated control, run on the OTHER Drucker-Prager -- the ASD twin.

SLOW TIER (`@pytest.mark.slow`, opt-in via `--runslow` / `LADRUNO_RUN_SLOW=1`).
MEASURED WALL TIME: see `_WALL_TIME` below.

WHY THIS FILE EXISTS
--------------------
`tests/test_r3_prandtl_collapse_gate.py` pins the Prandtl-Reissner collapse load
on the vanilla UW `nDMaterial DruckerPrager`, with an ASSOCIATED (psi = phi) leg
as its falsification control.  ADR-95 established that the fork's OTHER
Drucker-Prager, `ASDPlasticMaterial3D` with `DruckerPrager_YF`/`DruckerPrager_PF`,
is the SAME CONE (note 95 Sec. 3b: the two agree to the printed digit on the
non-associated leg, 1.0850 both, after #815 fixed the gradient and wp/94f fixed
the zero-dilatancy apex classification).

Two implementations of one cone must answer alike under BOTH flow rules or one
of them is wrong.  Under associated flow they did not: wp/94f UNIONED the yield
function's EUCLIDEAN apex region with the exact ELASTIC-METRIC one, and the union
keeps the WIDER region, which is the exact one only at zero dilatancy.  Under
associated flow the exact slope `K*etabar/G` is ~10x the Euclidean `eta`, so the
union kept a region ~10x too wide and silently apex-projected trial states whose
correct return is to the cone flank -- pinning them at `sigma_apex` with no
deviator and, under `tangent_type Continuum`, a ZERO tangent.  The fix (ADR-94
addendum, F8) makes the elastic-metric test REPLACE the Euclidean one.

WHAT IS ASSERTED
----------------
Both associated legs are run in this session on the same mesh, so the comparison
is not against a hard-coded number from another box:

1. the ASD associated leg is a CAPACITY under the gate's own three clauses;
2. its ratio matches the UW associated leg measured in the SAME session to the
   gate's own band half-width (+/- 3 %) -- one cone, two implementations;
3. it still satisfies the gate's falsification-control premise: the associated
   answer sits at least `ASSOC_MIN_SEPARATION` above the top of the
   NON-associated band, so the two flow rules answer differently;
4. the two associated LOAD-SETTLEMENT PATHS agree at matched settlement, not
   merely at the peak -- the pre-fix failure was a path failure;
5. the NON-associated ASD leg is unchanged and lands in the h0 = 1.0 band -- the
   regression guard for wp/94f, whose result this fix must not move.

The CHEAP variant of the same defect is `tests/test_f8_asd_dp_associated_apex.py`
(zone_a, a few seconds, four one-element subprocesses): it pins the misclassified
associated wedge against the closed-form cone return, and it FAILS on `9c2f964ea`
and passes here.  Run that first; this file is the deck-level confirmation, not
the primary gate.

MEASURED WALL TIME: **20.0 min** (1196.5 s) for the module, 5 tests, one process,
on the reference box (Windows 11, Intel oneAPI build, `system Pardiso`).  Slow
tier, opt-in, and it must never enter the default battery.

MEASURED (build `3324485f74f7ace1a7459c012e5c45167a541368`, h0 = 1.0, exact
`q_u = 138.907 kPa`):

    leg                     q_num    ratio    tail %       mode  ds/floor  CAP  s_end/B  wall s
    h1.0_assoc             268.75   1.9348     0.120     TARGET    4400.0  yes   0.1500    47.5
    h1.0_assoc_asd         268.38   1.9321     0.163     BUDGET    1250.0  yes   0.1367  1089.0
    h1.0_nonassoc_asd      150.71   1.0850     0.001     TARGET    4400.0  yes   0.1500    59.9

The two associated legs agree to 0.075 % worst over the whole common range, not
merely at the peak.  `BUDGET` is the gate's own "capacity WITH A NAMED
ALLOWANCE": the load had been flat for the last tenth of the run with the step
still 1250x the floor.

Pre-fix, on `9c2f964ea`, the ASD associated leg read **1.6307, mode BUDGET,
plateau NO (tail 8.270 %), free-advance NO (terminal step 50x the floor),
s/B 0.01685, 898 failed attempts, 894 s** -- and **zero refusals in the entire
log**, which is why the failure had to be found by classification rather than by
reading diagnostics.
"""
import os
import sys
import time

import pytest

_DIST = os.environ.get("LADRUNO_DIST_BIN")
if _DIST and os.path.isdir(_DIST):
    sys.path.insert(0, _DIST)

import test_r3_prandtl_collapse_gate as G  # noqa: E402

from _testbed import ops  # noqa: E402

pytestmark = [pytest.mark.zone_a, pytest.mark.slow, pytest.mark.t2a]

_WALL_TIME = "20.0 min (1196.5 s), 5 tests, one process -- see the module docstring"

# h0 = 1.0 is the cheapest rung of the gate's own refinement sequence and is
# where the associated wall was reported.  The gate's `CONTROL_H0` (0.5) rung is
# four times the cost and measures the same statement.
H0 = 1.0


@pytest.fixture(scope="module")
def assoc_pair(tmp_path_factory):
    out = tmp_path_factory.mktemp("r3_asd_assoc")
    build = ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "unknown"
    t0 = time.time()
    print(f"\n=== R3 associated control, UW vs ASD Drucker-Prager ==========")
    print(f"    engine build (ladrunoBuild)  : {build}")
    print(f"    openseespy module            : {os.path.abspath(ops.__file__)}")
    legs = {
        ("UW", True): G._run_leg(H0, True, str(out), material="UW"),
        ("ASD", True): G._run_leg(H0, True, str(out), material="ASD"),
        ("ASD", False): G._run_leg(H0, False, str(out), material="ASD"),
    }
    hdr = (f"{'leg':>22} {'q_num':>9} {'ratio':>8} {'tail %':>9} {'mode':>10} "
           f"{'ds/floor':>9} {'CAP':>5} {'s_end/B':>8} {'wall s':>8}")
    print(f"\n    exact q_u = {next(iter(legs.values()))['q_exact']:.3f} kPa")
    print(f"{hdr}\n{'-' * len(hdr)}")
    for r in legs.values():
        print(f"{r['tag']:>22} {r['qmax']:9.2f} {r['ratio']:8.4f} "
              f"{r['tail_pct']:9.3f} {r['mode']:>10} {r['headroom']:9.1f} "
              f"{'yes' if r['capacity'] else 'NO':>5} {r['s_end_over_B']:8.4f} "
              f"{r['wall_s']:8.1f}")
    print("\n    termination stories:")
    for r in legs.values():
        print(f"      {G._why(r)}")
    print(f"\n    TOTAL WALL TIME: {time.time() - t0:.1f} s")
    print("=" * 62)
    return dict(build=build, legs=legs)


def test_asd_associated_is_a_capacity(assoc_pair):
    """Clauses 1-3 of the gate's capacity rule, on the ASD material.  Pre-fix
    this leg stalled while STILL HARDENING, with no refusal of any kind -- the
    misclassified Gauss points reported success and a zero tangent."""
    r = assoc_pair["legs"][("ASD", True)]
    assert r["capacity"], (
        "the ASD associated leg is not a capacity (plateau="
        f"{r['plateau']}, free-advance={r['free']}, mode={r['mode']}). "
        "Pre-fix signature: no refusals at all, Gauss points pinned at "
        "sigma_apex with q = 0, step size collapsing while the load still "
        "climbs. " + G._why(r))


def test_asd_and_uw_associated_agree(assoc_pair):
    """One cone, two implementations, same mesh, same session.  The tolerance is
    the gate's own band half-width, not a fitted number."""
    a = assoc_pair["legs"][("ASD", True)]
    u = assoc_pair["legs"][("UW", True)]
    assert u["capacity"], (
        "the UW associated REFERENCE leg is not a capacity, so this comparison "
        "is inconclusive: " + G._why(u))
    rel = abs(a["ratio"] - u["ratio"]) / u["ratio"]
    assert rel <= G._BAND_HALFWIDTH, (
        f"the two Drucker-Pragers disagree under ASSOCIATED flow: ASD "
        f"{a['ratio']:.4f} vs UW {u['ratio']:.4f} ({100*rel:.2f} % apart, "
        f"allowed {100*G._BAND_HALFWIDTH:.0f} %). They agree to the printed "
        f"digit at psi = 0 (ADR-95 Sec. 3b), so a disagreement here is an "
        f"integrator defect, not physics.\n  ASD: {G._why(a)}\n  UW:  {G._why(u)}")


def test_asd_associated_is_separated_from_the_nonassociated_band(assoc_pair):
    """The gate's falsification premise, restated on the ASD material: the two
    flow rules must answer DIFFERENTLY, or the non-associated agreement carries
    no information."""
    r = assoc_pair["legs"][("ASD", True)]
    band_top = G._MEASURED[H0] * (1.0 + G._BAND_HALFWIDTH)
    assert r["ratio"] >= band_top + G.ASSOC_MIN_SEPARATION, (
        f"ASD associated landed at {r['ratio']:.4f}, within "
        f"{G.ASSOC_MIN_SEPARATION:.2f} of the non-associated band top "
        f"{band_top:.4f}. " + G._why(r))


def test_asd_and_uw_associated_curves_agree_at_matched_settlement(assoc_pair):
    """Stronger than the ratio: the two implementations must produce the SAME
    load-settlement PATH, not merely the same peak.  Both legs are read from
    their own CSVs and compared at matched s/B on a common grid.

    A peak-only comparison can be satisfied by two curves that disagree
    everywhere and happen to top out together; the pre-fix failure was a path
    failure (the ASD leg diverged into a wall at s/B 0.017 with the load still
    climbing), so the path is what this pins."""
    import csv

    def curve(r):
        with open(r["csv"], newline="") as fh:
            rows = list(csv.reader(fh))[1:]
        return [(float(x[1]), float(x[2])) for x in rows]

    a = curve(assoc_pair["legs"][("ASD", True)])
    u = curve(assoc_pair["legs"][("UW", True)])
    s_hi = min(a[-1][0], u[-1][0])
    grid = [s_hi * f for f in (0.1, 0.2, 0.4, 0.6, 0.8, 1.0)]

    def q_at(c, s):
        return min(c, key=lambda p: abs(p[0] - s))[1]

    worst, worst_s = 0.0, 0.0
    for s in grid:
        qa, qu = q_at(a, s), q_at(u, s)
        rel = abs(qa - qu) / abs(qu)
        if rel > worst:
            worst, worst_s = rel, s
    assert worst <= 0.01, (
        f"the two Drucker-Pragers follow DIFFERENT associated paths: worst "
        f"relative q difference {100*worst:.3f} % at s/B = {worst_s:.4f}, over "
        f"a common range up to s/B = {s_hi:.4f}. Measured 0.05 % on the fixed "
        f"build; pre-fix the ASD leg did not reach this range at all.")


def test_asd_nonassociated_leg_is_unchanged(assoc_pair):
    """wp/94f's own result is the thing this fix must not move: at etabar = 0
    the elastic-metric region strictly CONTAINS the Euclidean one, so replacing
    the union by the exact test is a no-op there -- by construction, and this
    checks it on the deck."""
    r = assoc_pair["legs"][("ASD", False)]
    lo, hi = G.BANDS[H0]
    assert r["capacity"], "the ASD psi = 0 leg stopped being a capacity: " + G._why(r)
    assert lo <= r["ratio"] <= hi, (
        f"the ASD psi = 0 leg moved out of its band: {r['ratio']:.4f} not in "
        f"[{lo:.4f}, {hi:.4f}] (ADR-95 / wp/94f measured 1.0850). " + G._why(r))
