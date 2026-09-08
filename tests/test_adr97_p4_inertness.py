"""ADR-97 P1 (wp/97b) — gate 4: `Backward_Euler` is BYTE-IDENTICAL, and the
closest-point map agrees with it exactly where the theory says it must.

ADR-97 **D1**: every result the fork has published on `ASDPlasticMaterial3D`
(Cerro Lindo M-series, the ADR-84 battery, the finite-strain D3 oracle) was
obtained on `Backward_Euler`.  `Closest_Point` is a new member function and a
new switch case; it touches no line the old integrator reads.  This file makes
that promise checkable instead of asserted.

Two halves:

1. **Byte identity.**  23 decks / 282 committed-stress rows dumped from the
   PRE-CHANGE binary (`3622d6214`) and stored in
   `Ladruno_implementation/adr97_oracle/baselines/`.  Each deck is re-run in a
   FRESH SUBPROCESS -- the per-tag `INT_OPT_*`/`DBL_OPT_*` option maps are
   `std::map<int,...>` keyed by material tag and shared across every instance of
   that tag *by design*, and the pytest heap leaks state across tests in one
   process (ADR-94 `capfd` / `WinError 6`).  Comparison is `==` on the doubles,
   not `allclose`.

2. **Agreement and divergence.**  `Closest_Point` and `Backward_Euler` are the
   SAME map exactly when the flow direction does not rotate over the step
   (proportional loading), because then
   ``sum_k dl_k m(sigma^k) == dl m(sigma_{n+1})``.  They must therefore agree to
   Newton tolerance on non-rotating-normal perfectly plastic decks, and differ
   measurably on Armstrong-Frederick decks (where `Backward_Euler` integrates
   the recovery term explicitly inside an otherwise implicit loop).  A gate that
   only checked agreement would pass on an integrator that silently WAS
   `Backward_Euler`.

Zone-A: the fast half is ~15 s; the full 23-deck sweep is `@pytest.mark.slow`
(measured 4 min 10 s on the reference box) and runs under `--runslow`.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from _testbed import ops

import test_adr97_p1_smooth as S  # noqa: E402

pytestmark = [pytest.mark.zone_a]

HERE = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(HERE, os.pardir, "Ladruno_implementation",
                        "adr97_oracle", "baselines")
BASELINE = os.path.join(BASE_DIR, "be_secant_baseline_3622d6214.json")
DUMPER = os.path.join(BASE_DIR, "dump_hist.py")

#: a representative slice of the baseline: the DEFAULT integrator/tangent pair
#: on both a plastic and an elastic leg, the two refusal-carrying tet decks, and
#: one deck per remaining yield function family.
FAST_DECKS = [
    "cube/vm/BE/Secant/plastic",
    "cube/vm/BE/Secant/elastic",
    "tet/vm/BE/default",
    "tet/dp/BE/default",
    "tet/mc/BE/strict1",
]


@pytest.fixture(scope="module")
def baseline():
    if not os.path.exists(BASELINE):
        pytest.skip("ADR-97 gate-4 baseline not present: %s" % BASELINE)
    with open(BASELINE) as f:
        return json.load(f)


def _run_child(deck, tmp_path):
    """Run ONE deck in a fresh interpreter and return its recorded history."""
    out = os.path.join(str(tmp_path), "hist.json")
    env = dict(os.environ)
    env.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    p = subprocess.run([sys.executable, DUMPER, "--one", deck, out],
                       cwd=HERE, env=env, capture_output=True, text=True)
    if not os.path.exists(out):
        pytest.fail("child process produced no history for %r:\n%s"
                    % (deck, (p.stdout + p.stderr)[-2000:]))
    with open(out) as f:
        return json.load(f)


def _assert_identical(deck, ref, got):
    assert "stress" in ref, "baseline deck %r never built: %r" % (deck, ref)
    assert "stress" in got, "current binary failed deck %r: %r" % (deck, got)
    assert got["codes"] == ref["codes"], (
        "deck %r: analyze() return codes changed %s -> %s"
        % (deck, ref["codes"], got["codes"]))
    assert len(got["stress"]) == len(ref["stress"]), (
        "deck %r: %d committed steps, was %d"
        % (deck, len(got["stress"]), len(ref["stress"])))
    nbit, worst = 0, 0.0
    for a, b in zip(ref["stress"], got["stress"]):
        for x, y in zip(a, b):
            if x != y:
                nbit += 1
                worst = max(worst, abs(x - y))
    assert nbit == 0, (
        "deck %r: %d of %d committed stress components changed (worst |d| = "
        "%.3e). ADR-97 D1 says Backward_Euler is byte-identical -- either a"
        " shared code path was edited, or the baseline needs a deliberate,"
        " documented regeneration (never to make this gate green)."
        % (deck, nbit, len(ref["stress"]) * 6, worst))


@pytest.mark.parametrize("deck", FAST_DECKS)
def test_gate4_backward_euler_is_byte_identical_fast(baseline, tmp_path, deck):
    """The representative slice, on every push."""
    _assert_identical(deck, baseline[deck], _run_child(deck, tmp_path))


@pytest.mark.slow
def test_gate4_backward_euler_is_byte_identical_full(baseline, tmp_path):
    """All 23 decks / 282 committed-stress rows.  Measured 4 min 10 s (23 fresh
    interpreters, each constructing the material catalogue)."""
    changed = []
    for deck, ref in sorted(baseline.items()):
        got = _run_child(deck, tmp_path)
        try:
            _assert_identical(deck, ref, got)
        except AssertionError as exc:
            changed.append(str(exc))
    assert not changed, "\n".join(changed)


# ===========================================================================
# agreement / divergence
# ===========================================================================
NON_ROTATING = {
    "vm-triaxial": (S.VM_PATHS["triaxial"], "vm"),
    "vm-simple-shear": (S.VM_PATHS["simple-shear"], "vm"),
    "dp-compress": (S.DP_COMPRESS, "dp"),
}


@pytest.mark.parametrize("case", list(NON_ROTATING))
def test_gate4_cp_and_be_agree_on_non_rotating_perfectly_plastic_decks(case):
    """Where the flow direction does NOT rotate over the step, the cutting
    plane's fixed point ``sigma_tr - sum_k dl_k E m(sigma^k)`` and the closest
    point's ``sigma_tr - dl E m(sigma_{n+1})`` are the same point.  The decks are
    named non-rotating on purpose: on a rotating-normal path the two maps differ
    by the cutting-plane path error, which is gate 3's pinned contrast and not a
    gate-4 failure."""
    path, fam = NON_ROTATING[case]
    if fam == "vm":
        cp = S.drive(lambda t: S.mat_vm(t), path, nstep=10)
        be = S.drive(lambda t: S.mat_vm(t, method="Backward_Euler",
                                        tangent="Secant"), path, nstep=10)
    else:
        if not S._constructible(lambda t: S.mat_dp(t)):
            pytest.skip("DruckerPrager Closest_Point unavailable")
        cp = S.drive(lambda t: S.mat_dp(t), path, nstep=10)
        be = S.drive(lambda t: S.mat_dp(t, method="Backward_Euler",
                                        tangent="Secant"), path, nstep=10)
    assert all(c == 0 for c in cp["codes"] + be["codes"])
    d = S._rel(cp["sigma"][-1], be["sigma"][-1])
    print("gate 4 agreement %-16s CP vs BE rel difference = %.3e" % (case, d))
    assert d < 1e-8, (cp["sigma"][-1], be["sigma"][-1])


def test_gate4_cp_and_be_differ_measurably_on_armstrong_frederick():
    """The other half of the gate: an integrator that silently WAS
    `Backward_Euler` would pass every agreement test above.  Armstrong-Frederick
    is where the two maps must part company -- `Backward_Euler` evaluates AF's
    recovery term ``-c_r ||dev m||_eq alpha`` at the TRIAL alpha inside an
    otherwise implicit loop, so its answer depends on how many corrector steps
    the loop happened to take (P0 measured a spread of 9.15 on a stress of ~46
    over four equally valid iterate paths)."""
    if not S._constructible(lambda t: S.mat_vm_af(t)):
        pytest.skip("VonMises + ArmstrongFrederick specialization unavailable")
    path = S.VM_PATHS["rotating-normal"]
    cp = S.drive(lambda t: S.mat_vm_af(t), path, nstep=10)
    be = S.drive(lambda t: S.mat_vm_af(t, method="Backward_Euler",
                                       tangent="Secant"), path, nstep=10)
    assert all(c == 0 for c in cp["codes"] + be["codes"])
    d = S._rel(cp["sigma"][-1], be["sigma"][-1])
    print("gate 4 divergence  AF rotating-normal  CP vs BE rel difference "
          "= %.3e" % d)
    assert d > 1e-4, (
        "Closest_Point and Backward_Euler agree on an Armstrong-Frederick "
        "rotating-normal path -- the new integrator may be dispatching to the "
        "old one")
