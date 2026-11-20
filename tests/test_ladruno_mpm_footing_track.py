"""WP-244 (fix #1112): the footing box's half-cell lead, row P1 of wp-briefs "WP-244".

P1, the rigid-translation patch (the 232l T13 lesson: a rule that changes which grid nodes are
prescribed under motion is shown on a free translation before any footing number). A stress-free
elastic block is carried down by a tracked gridvel box spanning its full width, at dS = h/10 per
step, through 2.5 cells. The exact answer is a rigid translation: zero stress, zero condition
force, every particle displaced by s. Both leads (the box's lower bound at 0 and at -h/2) are run;
they differ in which rows are prescribed once s >= h/2, which this isolates from any soil.

Also checked: before s = h/2 the two leads prescribe the same node rows, so the footing records
that end before s = h/2 stay bit-identical (S3 r1 / r2, L4 r1, L3 r1 / r2).
"""
import json
from itertools import pairwise
from pathlib import Path

import pytest
from _testbed import ops

pytestmark = pytest.mark.zone_a

HOST = 244
H = 0.25
E, NU = 1.0e5, 0.3
DS = H / 10.0                   # settlement per step (dt = 1)
NSTEP = 25                      # s = 2.5 h
TOL_SIG = 1e-8 * E * DS / H
WIDTH = 2.0


def _l(v):
    return list(v) if isinstance(v, (list, tuple)) else [v]


def _build(lead, footing=False):
    """footing=False: P1's free block, the box spanning its width. footing=True: the same block
    on a fixed base (rollers at x = 0 and 2, base fixed), a footing of half-width 0.5 m."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    if footing:
        ops.mpm("grid", HOST, "-origin", 0.0, -2.0, "-h", H, "-cells", int(WIDTH / H),
                int(3.0 / H), "-basis", "gimp")
    else:
        ops.mpm("grid", HOST, "-origin", -1.0, -4.0, "-h", H, "-cells", int(4.0 / H),
                int(5.0 / H), "-basis", "gimp")
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 0.0, -2.0, WIDTH, 0.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    if footing:
        ops.mpm("fix", "xmin", "x")
        ops.mpm("fix", "xmax", "x")
        ops.mpm("fix", "ymin", "xy")
    hi_x = 0.5 if footing else WIDTH + H
    ops.mpm("gridvel", "footing", "box", -H, -lead * H, hi_x, 1.0e9, "v", 0.0, -DS,
            "start", 0.0, "track", "y")
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    if footing:                     # forces ~ E dS = 2.5e3 kN: a displacement test
        ops.test("NormDispIncr", 1e-12, 25)
    else:
        ops.test("NormUnbalance", 1e-8, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _run(lead):
    _build(lead)
    x0 = _l(ops.mpm("particleData", "x"))
    worst = {"sig": 0.0, "force": 0.0, "disp": 0.0}
    for k in range(1, NSTEP + 1):
        ops.mpm("analyze", 1, 1.0)
        s = DS * k
        sig = _l(ops.mpm("particleData", "sig"))
        x = _l(ops.mpm("particleData", "x"))
        f = _l(ops.mpm("condForce", "footing"))
        worst["sig"] = max(worst["sig"], max(abs(v) for v in sig))
        worst["force"] = max(worst["force"], abs(f[1]))
        worst["disp"] = max(worst["disp"],
                            max(abs(x[2 * p] - x0[2 * p]) for p in range(len(x) // 2)),
                            max(abs((x[2 * p + 1] - x0[2 * p + 1]) + s) for p in range(len(x) // 2)))
    return worst


LEAD0_MISS = (
    "pre-registered P1 row, lead 0 (the old box, #1112): the block lags exactly one step "
    "(displacement error = dS) from the step that starts at s = 0.9 h. There the only prescribed "
    "row (y = 0) overlaps the top particles by 0.1 h, its nodal mass falls below the mode-A mass "
    "cutoff, so it leaves the active set and the footing prescribes nothing for that step. "
    "With `mpm massCutoff 1e-14` the lag is gone. The lead keeps >= h/2 of soil on the lowest "
    "prescribed row's support.")


@pytest.mark.parametrize("lead", [pytest.param(0.0, marks=pytest.mark.xfail(
                                      strict=True, raises=AssertionError, reason=LEAD0_MISS)),
                                  0.5])
def test_p1_rigid_translation(lead):
    w = _run(lead)
    print(f"P1 lead {lead}: max|sig| {w['sig']:.3e} (tol {TOL_SIG:.3e}), max|condForce_y| "
          f"{w['force']:.3e} (tol {TOL_SIG * WIDTH:.3e}), max displacement error {w['disp']:.3e} m")
    assert w["sig"] <= TOL_SIG, w
    assert w["force"] <= TOL_SIG * WIDTH, w
    assert w["disp"] <= 1e-9 * H, w


def _footing_state(lead, nstep):
    _build(lead, footing=True)
    for _ in range(nstep):
        ops.mpm("analyze", 1, 1.0)
    return (_l(ops.mpm("particleData", "x")) + _l(ops.mpm("particleData", "sig"))
            + _l(ops.mpm("condForce", "footing")))


def test_lead_is_bit_identical_before_half_a_cell():
    """An elastic footing on a fixed base: before s = h/2 (4 steps, s = 0.4 h) the two leads
    prescribe the same rows and agree bit for bit; after it (6 steps, s = 0.6 h) the lead has
    added the row at -h and they differ, so the comparison can fail."""
    assert _footing_state(0.0, 4) == _footing_state(0.5, 4)
    assert _footing_state(0.0, 6) != _footing_state(0.5, 6)


# --------------------------------------------------------------------------- #
#  F1, the artifact window (wp-briefs "WP-244"): HS S3 at r2 / r3 to s = 2 h   #
#  The records are written by the WP-244 runs (`s3_run`, s_end = 2 h) and read  #
#  here; lead 0 is the #1112 reproduction, reported.                           #
# --------------------------------------------------------------------------- #
WP244 = Path(__file__).resolve().parent / "data" / "mpm" / "wp244"


def _worst_drop(curve, s_min=0.01):
    """The largest one-step drop of q after s_min, as a fraction of the running maximum."""
    run, worst = 0.0, 0.0
    for (_, q0), (s1, q1) in pairwise(curve):
        run = max(run, q0)
        if s1 > s_min and run > 0.0:
            worst = max(worst, (q0 - q1) / run)
    return worst


@pytest.mark.parametrize("r", [2, 3])
def test_f1_artifact_window(r):
    h = 0.5 / 2 ** (r - 1)
    new = json.loads((WP244 / f"f1_hs_r{r}_lead0.5.json").read_text(encoding="utf-8"))
    old = json.loads((WP244 / f"f1_hs_r{r}_lead0.0.json").read_text(encoding="utf-8"))
    w_new, w_old = _worst_drop(new["curve"]), _worst_drop(old["curve"])
    print(f"F1 r{r}: lead 0.5 reaches s/h {new['curve'][-1][0] / h:.2f}, worst drop "
          f"{100 * w_new:.1f} %; lead 0 (reported) reaches s/h {old['curve'][-1][0] / h:.2f}, "
          f"worst drop {100 * w_old:.1f} %")
    assert new["reached"] and abs(new["curve"][-1][0] - 2.0 * h) < 1e-9, new["stopped"]
    assert w_new <= 0.05, w_new
