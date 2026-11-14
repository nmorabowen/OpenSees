"""WP-232f: the fine-resolution P1 gates of MPM mode A (tests/mpm_implicit_cases.txt).

These are the pre-registered gate rows at their finest resolution (cantilever QS and transient
r2 / r3, the mild-cone footing r2 / r3, the plate with a hole): 10^4 - 5 * 10^4 particles over
up to 400 implicit steps, too long for Zone-A's 40-minute job. They run where the warrant is
assembled: set LADRUNO_MPM_GATES=1 (and dist/bin/opensees.pyd on PYTHONPATH). Without it every
test here is skipped with that reason — they are never silently passed.

An expected value still `TBD` (232f-V filling it) is skipped loudly, never passed.
Written from tests/mpm_implicit_interface.md and the cases file only.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_implicit_harness as H  # noqa: E402

pytestmark = [pytest.mark.zone_a,
              pytest.mark.skipif(os.environ.get("LADRUNO_MPM_GATES") != "1",
                                 reason="fine-resolution MPM gates: set LADRUNO_MPM_GATES=1 "
                                        "(too long for Zone-A)")]

GATES = ["cantilever_qs_r2", "cantilever_qs_r3", "cantilever_dyn_r2", "cantilever_dyn_r3",
         "plate_hole_dw050"]

_CACHE: dict = {}


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    return ops


def run(name):
    if name not in _CACHE:
        _CACHE[name] = H.run_python(H.case_by_name(name), _ops())
    return _CACHE[name]


@pytest.mark.parametrize("name", GATES)
def test_mode_a_gate(name):
    c = H.case_by_name(name)
    r = run(name)
    res = H.evaluate(c, r, runs=run)
    for lab, v, det in res:
        print(f"{v.upper():6s} {lab}: {det}")
    # the resolution trend is reported by test_cantilever_error_trend (a pre-registered miss)
    failed = [f"{lab}: {det}" for lab, v, det in res if v == "fail" and "relerrle" not in lab]
    assert not failed, "\n".join(failed)
    tbd = [lab for lab, v, _ in res if v == "tbd"]
    if tbd and all(v in ("tbd", "report") for _, v, _ in res):
        pytest.skip("expected values TBD (232f-V): " + "; ".join(tbd))


@pytest.mark.xfail(strict=True, reason=(
    "pre-registered miss (gate row 'large-deflection cantilever, QS': error non-increasing "
    "across resolutions). Measured 0.10 / 0.23 / 0.39 % at r1 / r2 / r3 (renormalize, c_m 5e-2): "
    "every resolution sits inside the 5 % band, but the MPM tip converges to about 0.4 % beyond "
    "the FE Richardson value, so the error grows as the coarse-grid cancellation fades. "
    "Model-form diagnosis: uGIMP particle quadrature at 2x2 per cell and the renormalised "
    "free-surface supports against a hypo-H8 slab -- a different discrete limit, not a "
    "divergence. The criterion is not loosened."))
def test_cantilever_error_trend():
    for name in ("cantilever_qs_r2", "cantilever_qs_r3"):
        c = H.case_by_name(name)
        res = H.evaluate(c, run(name), runs=run)
        trend = [(lab, v, det) for lab, v, det in res if "relerrle" in lab]
        assert trend and all(v == "pass" for _, v, _ in trend), trend


FOOTING_MISS = (
    "pre-registered miss, recorded 2026-11-14 (WP-232f footing diagnosis, PR #1083): the "
    "mild-cone footing does not reach the band [0.90, 1.10] x 138.907 kPa at r3. Two findings. "
    "(a) PLATEAU OVERSHOOT: r1 q_max 188.6 (1.36x), in the locked-FE range (232f-V: H8 without "
    "B-bar 1.577 / 1.290 at 4 / 8 cells across B, with B-bar 1.018 / 0.941); the measurement is "
    "exact (condForce = sigma * width, ratio 1.000000, test_condforce_is_the_reaction) and a "
    "lower elastic nu does not move it (nu 0.3 / 0.2: 190.3 / 190.5), so the diagnosis is "
    "plastic-incompressibility locking under isochoric flow (psi = 0) with 2x2 particles per "
    "cell -- the ADR's phase 4 (anti-locking), not phase 1. (b) NEWTON DIVERGENCE at r2 / r3 "
    "(r2 stops near t = 30, r3 at step 42): iterates invert particles and carry centres out of "
    "the grid; NewtonLineSearch (t = 28) and KrylovNewton (t = 29) do not help. The element "
    "tangent is not the cause: ||K - K_FD|| / ||K|| (FD converged over eps 1e-8..1e-10) is "
    "3e-6 elastic undeformed, 1e-3 at 0.2 % strain, mean 3.7 % on the plastic footing state; the "
    "worst particles (rel ~1) carry |K| ~ 75 against ~1e5 elastic: DruckerPrager at its apex "
    "(weightless soil, sigma_y 0.2 kPa near the free surface) returns its regularised 1e-3 Ce "
    "tangent where the response is perfectly plastic, so the nodes they support are nearly "
    "singular. The band is not touched.")


@pytest.mark.xfail(strict=True, reason=FOOTING_MISS)
@pytest.mark.parametrize("name", ["footing_mildcone_r2", "footing_mildcone_r3"])
def test_footing_gate(name):
    c = H.case_by_name(name)
    r = run(name)
    res = H.evaluate(c, r, runs=run)
    for lab, v, det in res:
        print(f"{v.upper():6s} {lab}: {det}")
    assert r["stopped"] is None, f"{name}: stopped at step {r['stopped']}"
    got = H.observable(c, r, [e for e in c["expects"] if e["obs"] == "cond_pressure_max"][0])[0]
    assert 0.90 * 138.907 <= got <= 1.10 * 138.907, got
