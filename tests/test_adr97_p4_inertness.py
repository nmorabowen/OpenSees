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

Cross-platform note (Zone-A run 34174010743, the first Linux run of this file):
the baseline was dumped on Windows/MSVC.  GCC/Linux reproduces every deck to
5.8e-09 worst absolute (VM at Newton tolerance) and 2e-20 (DP, FP noise) -- a
compiler/libm difference, not a code path.  A byte-identity gate is only
meaningful on the platform that produced the baseline (ADR-94 lesson: 1e-9 pins
failed twice on Linux), so `==` is enforced on win32 and a 1e-6 RELATIVE bound
(the fork's cross-platform float-pin floor) everywhere else.  Regenerate the
baseline on Windows only, deliberately, and say so in the commit.

ADR-97 wp/97e (P4) baseline regeneration, 2026-09-07: re-pointing
`Numerical_Algorithmic_FirstOrder/SecondOrder` at the actual committed map
(see `ASDPlasticMaterial3D::numerical_tangent_of_committed_map()`) changed
`cube/vm/BE/Numerical_Algorithmic_FirstOrder/plastic` by up to 5.428e-09
absolute in 40 of its 60 committed-stress components -- EXPECTED, not a
regression, and the ONLY entry of the 23 that moved (every other deck,
including all 4 other `cube/vm/BE/*` tangent types, re-verified bit-identical
against the same binary). Root cause: `cube/*` decks are LOAD-CONTROLLED with
free DOFs (`test_adr94_hlist_numerics._cube_build`), not fully prescribed --
the tangent this option returns feeds the outer Newton's stiffness matrix at
EVERY outer iteration (`Backward_Euler` calls `ComputeTangentStiffness()` at
the tail of every `setTrialStrainIncr()`, not just the converged one), so a
tangent that is now materially different (FD of `Backward_Euler` itself,
instead of FD of the unrelated `compute_local_stress()` map) changes the
outer Newton's convergence PATH and therefore its converged point within the
`NormDispIncr` tolerance ball (`1e-12` here) -- the same order of magnitude
as the cross-platform compiler-noise floor measured above. This is exactly
what the WP intends (a materially better tangent) and does not touch
`Backward_Euler`'s own residual: the paired `.../elastic` leg of the SAME
deck is untouched (0 of 60 components changed) because an elastic trial's
map is exactly linear, so old and new FD agree exactly there, with no
path-dependence to expose. `Backward_Euler`'s own source is confirmed
untouched by `grep -n compute_local_stress` inside its body (ADR-97 D1).

Zone-A, 5.9 s for the whole file: a fresh-interpreter deck costs ~0.2 s, so
BOTH the representative slice and the full 23-deck sweep run on every push.
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
    # `stdin=subprocess.DEVNULL` is LOAD-BEARING on Windows: with the inherited
    # stdin, Popen tries to DuplicateHandle whatever pytest's fd-level capture
    # left there and raises `OSError: [WinError 6] The handle is invalid`
    # intermittently -- the same trap `_testbed/subprocess_run.py` documents,
    # and it looks exactly like the child crashing.  `errors="replace"` for the
    # same reason: an `opserr` line with a non-ASCII byte otherwise raises
    # UnicodeDecodeError inside the capture on a cp1252 console.
    p = subprocess.run([sys.executable, DUMPER, "--one", deck, out],
                       cwd=HERE, env=env, capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, encoding="utf-8",
                       errors="replace")
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
    nbit, worst, scale = 0, 0.0, 0.0
    for a, b in zip(ref["stress"], got["stress"]):
        for x, y in zip(a, b):
            scale = max(scale, abs(x))
            if x != y:
                nbit += 1
                worst = max(worst, abs(x - y))
    if sys.platform != "win32":
        # Not the baseline's platform: MSVC vs GCC/libm differ at 1e-9 absolute
        # (measured, see the module docstring).  Enforce the fork's
        # cross-platform floor instead of bit equality.
        assert worst <= 1e-6 * max(scale, 1.0), (
            "deck %r: committed stress differs from the Windows baseline by "
            "%.3e (scale %.3e) on %s -- beyond the 1e-6 cross-platform floor, "
            "so this is a code-path change, not compiler noise."
            % (deck, worst, scale, sys.platform))
        return
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


def test_gate4_backward_euler_is_byte_identical_full(baseline, tmp_path):
    """All 23 decks / 282 committed-stress rows, one fresh interpreter each.

    Measured 5.9 s for the whole file -- a child costs ~0.2 s, so there is no
    reason to hide this behind `--runslow`.  (It was written as `@slow` on the
    assumption that 23 interpreter starts would be minutes; measuring is
    cheaper than assuming.)"""
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
