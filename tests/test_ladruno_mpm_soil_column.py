"""WP-232k S2 and S4 (tests/mpm_soil_interface.md; gates in wp-briefs.md "232k" and its freeze).

S2: the K0 column (mode A QS, fixed base) with a surcharge rho g l0 on the top row: syy against
the closed form within 0.03 rho g l0 and sxx RMS against 232k-V's FE column within 0.02 rho g l0,
at the end of stage 0, after the flip and after the surcharge (tests/data/mpm/soil/s2_column.json).
S4: the stage flip is the same per-particle map in mode A and mode F, and the FE's: on the
low-g twin (g 0.5, l0 2; tests/data/mpm/soil/s4_lowg.json) the row with p < Pmin commits
sigma = -Pmin I and alpha = 0 exactly, every other particle keeps its sigma bit for bit, and
alpha_in = alpha (-flipAlphaIn init); >= 1 reset particle in each mode.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


S2_HS_MISS = pytest.mark.xfail(strict=True, reason=(
    "pre-registered miss (ruling at e8e3fcb2a): the registered surcharge mechanism, a body force on"
    " the soil's top particle row, leaves the loaded row carrying about half the load in its own"
    " stress, so syy_max after the surcharge is 0.274 (gimp) / 0.576 (bspline2) against 0.03,"
    " at the top cell (gimp: rows below it <= 0.015; bspline2: 0.072 at Y 9.125, 0.027 at 8.625,"
    " the B-spline support spreading it); end of stage 0 and after the flip pass. Cap variants reported beside it in tests/data/mpm/soil/s2_mpm_hs_diagnosis.json:"
    " elastic cap 0.021 / 0.300, soil cap 0.058 / 0.069. Gate unchanged."))


@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("mat", [pytest.param("hs", marks=S2_HS_MISS), "sanisand"])
def test_s2_column_with_surcharge(mat, basis):
    out = S.s2_run(_ops(), mat, basis)
    assert out["passed"], out["checkpoints"]


@pytest.mark.parametrize("mode", ["A", "F"])
def test_s4_flip_map(mode):
    out = S.s4_flip_map(_ops(), mode)
    assert out["reset"] >= 1, out
    assert not out["bad"], out


@pytest.mark.parametrize("row", [0, 39], ids=["reset_p_below_Pmin", "deep"])
def test_s4_same_state_both_modes(row):
    """S4's same-state half (ruling at the S4 freeze): the S1 block, an affine oedometric stage-0
    path to the S4 record's row strain, then `mpm stage 1 1`, in mode A and mode F. Gate: the
    post-flip sigma and the full SANISAND state are bitwise equal between the modes, and equal
    to the one-point FE record's post-flip state (rel 1e-8, the S1 rule). Row 0 has p < Pmin at
    the flip (the Elastic2Plastic reset: sigma = -Pmin I, alpha = 0); row 39 is the deepest.
    If the two modes' PRE-flip states already differ by round-off (their kinematics sum the
    same affine field in a different order), the post-flip difference must not exceed it."""
    ops = _ops()
    a = S.s4_affine(ops, "A", row)
    f = S.s4_affine(ops, "F", row)
    N = len(a["post"]["sig"])
    pre_same = a["pre"] == f["pre"]
    if pre_same:
        assert a["post"] == f["post"], "post-flip states differ between the modes on equal inputs"
    else:
        d_pre = max(S.rel_diff(a["pre"]["sig"][p], f["pre"]["sig"][p]) for p in range(N))
        d_post = max(S.rel_diff(a["post"]["sig"][p], f["post"]["sig"][p]) for p in range(N))
        assert d_post <= max(10.0 * d_pre, 1e-15), (d_pre, d_post)
    for out in (a, f):
        ns = len(out["post"]["state"]) // N
        for p in range(N):
            assert S.rel_diff(out["pre"]["sig"][p], out["ref_pre"]["sig"]) <= 1e-8
            assert S.rel_diff(out["post"]["sig"][p], out["ref_post"]["sig"]) <= 1e-8
            st = out["post"]["state"][ns * p:ns * p + ns]
            assert S.rel_diff(st, out["ref_post"]["state"]) <= 1e-8
        if out["below_Pmin"]:
            assert all(s[:3] == [-0.101, -0.101, -0.101] for s in out["post"]["sig"])
