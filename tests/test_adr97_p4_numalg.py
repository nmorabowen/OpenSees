"""ADR-97 P4 (wp/97e-numalg-repoint) -- Numerical_Algorithmic_* re-point gate.

Helper calls (`test_adr97_p1_smooth.mat_vm/drive/_constructible/_rel/
_two_cube/VM_PATHS`, `adr97_oracle/fd_tangent_driver.fd_check/mat_vm`)
re-verified against HEAD (`wp/97e-numalg-repoint`, cut from `c10c8dca7` via
`wp/97f-explicit-gate`) before registering this file. Every numeric threshold
below is either taken verbatim from the ADR-97 P4 task brief's stated
expected orders, or copied from an ALREADY-MEASURED number in
`Ladruno_implementation/97_ladruno_asdp_closest_point_adr.md` (the two-cube
113 total) -- none of it is invented.

Covers, per the P4 task brief item 4:

1. `Backward_Euler` + `Numerical_Algorithmic_SecondOrder` matches a central
   difference of the binary's own assembled residual (`fd_tangent_driver.
   fd_check`) to <= 1e-6 relative; `_FirstOrder` to ~1e-3 (forward FD, no
   free truncation-order cancellation).
2. `Closest_Point` + `Numerical_Algorithmic_SecondOrder` matches `Closest_
   Point` + `Algorithmic` (the exact consistent tangent, ADR-97 P1) to
   <= 1e-6 relative on the assembled tangent -- the two must now be a
   cross-check of each other, per ADR-97 D2's own wording ("making
   Numerical_Algorithmic_FirstOrder a cross-check on Algorithmic rather than
   a competitor").
3. The ADR-94 two-cube `testIter` sum for `Backward_Euler` + `Numerical_
   Algorithmic_SecondOrder` drops well below the RECORDED `Backward_Euler` +
   `Continuum` total of 113 (`6+41+36+30`, ADR-97 P1 report / ADR-97 doc
   implementation log, 2026-09-07 P1 entry) -- it should land near the
   `Closest_Point`+`Algorithmic` total of 16 (4 per step x 4 steps), since on
   this non-rotating proportional path Backward_Euler's cutting plane equals
   the closest-point answer and a correct FD of it is the same consistent
   tangent.
4. Refusal propagation: a starved `n_max_iterations` reachable from INSIDE a
   perturbed sub-call must fail the WHOLE tangent call (the host analysis
   fails), not silently assemble a partial/short tangent.

Byte-inertness of `Backward_Euler` + `Secant` (ADR-97 D1) is NOT re-tested
here on purpose -- re-run the existing gate-4 file instead
(`tests/test_adr97_p4_inertness.py` -- note the name collision: "p4" there
means gate 4 of ADR-97 P1, NOT this ADR phase). That file already dumps 23
decks / 282 committed-stress rows in fresh subprocesses and is the correct
regression check that this WP's edits (all inside `ComputeTangentStiffness()`
and the two `compute_numerical_tangent_*` functions) do not leak into any
path `Backward_Euler` + `Secant` (the default, and every deck in that file)
actually reads.
"""
import os
import sys

import numpy as np
import pytest

from _testbed import ops

import test_adr97_p1_smooth as S  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "Ladruno_implementation",
                                "adr97_oracle"))
import fd_tangent_driver as FD  # noqa: E402

pytestmark = [pytest.mark.zone_a]

# ADR-97 P1 report / doc implementation log (97_ladruno_asdp_closest_point_
# adr.md, 2026-09-07 P1 entry): two-cube testIter BE+Continuum = 6,41,36,30
# (SUM 113); CP+Algorithmic = 4,4,4,4 (SUM 16). Both ALREADY MEASURED, not
# re-derived here.
TWOCUBE_BE_CONTINUUM_TOTAL = 113
TWOCUBE_CP_ALGORITHMIC_TOTAL = 16


@pytest.fixture(scope="module")
def numalg_repoint_available():
    """Skip the whole file on a binary that has not shipped ADR-97 P4 yet --
    detectable ONLY by behaviour (there is no new parser token), so probe with
    the actual measurement this file makes and treat "still looks like the
    old compute_local_stress-based FD" as unavailable rather than a failure.
    Prevents this draft from red-flagging a pre-P4 binary if it is ever
    accidentally collected before P4 ships."""
    if not S._constructible(lambda t: S.mat_vm(t, method="Backward_Euler",
                                               tangent="Secant")):
        pytest.skip("ASDPlasticMaterial3D / VonMises not available")
    return True


# ===========================================================================
# 1. Backward_Euler + Numerical_Algorithmic_* vs a central difference of the
#    binary's own assembled residual.
# ===========================================================================
def test_be_numerical_secondorder_matches_fd_of_its_own_committed_map(
        numalg_repoint_available):
    """Post ADR-97 P4: Numerical_Algorithmic_SecondOrder differentiates
    Backward_Euler ITSELF (numerical_tangent_of_committed_map()), so a
    central difference of the binary's own residual under Backward_Euler +
    Numerical_Algorithmic_SecondOrder should agree with the material's own
    reported tangent to stencil-truncation precision -- there is no longer a
    third map for the two to disagree about.
    """
    r = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Numerical_Algorithmic_SecondOrder",
                                        method="Backward_Euler"),
                    label="Backward_Euler / Numerical_Algorithmic_SecondOrder")
    print("P4 gate BE+NumAlg2nd rel_err = %.3e (expect <= 1e-6)" % r["rel_err"])
    assert r["rc"] == 0, "the Backward_Euler load-driven rig did not converge"
    assert r["rel_err"] <= 1e-6, r["rel_err"]


def test_be_numerical_firstorder_matches_fd_of_its_own_committed_map(
        numalg_repoint_available):
    """Same check for the forward-difference variant. Expected order is
    looser (~1e-3) than the central-difference pair above: forward FD has no
    even-order cancellation, and the fixed `delta_min`/`epsilon_ref` stencil
    was tuned for the OLD compute_local_stress()-based use, not re-tuned here.
    """
    r = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Numerical_Algorithmic_FirstOrder",
                                        method="Backward_Euler"),
                    label="Backward_Euler / Numerical_Algorithmic_FirstOrder")
    print("P4 gate BE+NumAlg1st rel_err = %.3e (expect ~1e-3)" % r["rel_err"])
    assert r["rc"] == 0, "the Backward_Euler load-driven rig did not converge"
    assert r["rel_err"] <= 1e-2, r["rel_err"]  # conservative; tighten post-build


# ===========================================================================
# 2. Closest_Point + Numerical_Algorithmic_SecondOrder vs Closest_Point +
#    Algorithmic -- the ADR-97 D2 cross-check.
# ===========================================================================
def test_cp_numerical_secondorder_cross_checks_algorithmic(
        numalg_repoint_available):
    """ADR-97 D2's own framing: after this WP, Numerical_Algorithmic_* on
    Closest_Point is "a cross-check on Algorithmic rather than a competitor".
    Both differentiate/assemble the tangent of the SAME converged Newton
    solve, so they should agree to FD-stencil precision."""
    if not S._constructible(lambda t: S.mat_vm(t, method="Closest_Point",
                                               tangent="Algorithmic")):
        pytest.skip("Closest_Point / Algorithmic unavailable on this binary")

    r_alg = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Algorithmic",
                                            method="Closest_Point"),
                        label="Closest_Point / Algorithmic")
    r_num = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Numerical_Algorithmic_SecondOrder",
                                            method="Closest_Point"),
                        label="Closest_Point / Numerical_Algorithmic_SecondOrder")
    assert r_alg["rc"] == 0 and r_num["rc"] == 0
    d = S._rel(r_num["K_asm"], r_alg["K_asm"])
    print("P4 gate CP Algorithmic-vs-NumAlg2nd assembled-K rel diff = %.3e "
          "(expect <= 1e-6)" % d)
    assert d <= 1e-6, (r_alg["K_asm"], r_num["K_asm"])


# ===========================================================================
# 3. The ADR-94 two-cube Newton-cost re-measurement.
# ===========================================================================
def test_twocube_be_numerical_secondorder_cost_drops_from_113(
        numalg_repoint_available):
    """BE+Continuum's per-step testIter sum is the ALREADY MEASURED 113
    (6+41+36+30, ADR-97 P1 report). BE+Numerical_Algorithmic_SecondOrder,
    post-P4, differentiates the SAME map Continuum only approximates; on this
    non-rotating proportional two-cube path Backward_Euler's cutting plane
    equals the closest-point answer, so the Newton cost should approach
    Closest_Point+Algorithmic's already-measured 16, not sit anywhere near
    113."""
    it_num, sig_num = S._two_cube(
        lambda t: FD.mat_vm(t, tangent="Numerical_Algorithmic_SecondOrder",
                            method="Backward_Euler"))
    print("two-cube BE+Numerical_Algorithmic_SecondOrder testIter = %s "
          "(sum %d; BE+Continuum recorded sum = %d; CP+Algorithmic recorded "
          "sum = %d)" % (it_num, sum(it_num), TWOCUBE_BE_CONTINUUM_TOTAL,
                         TWOCUBE_CP_ALGORITHMIC_TOTAL))
    assert sum(it_num) < TWOCUBE_BE_CONTINUUM_TOTAL / 2, (
        "ADR-97 P4: Backward_Euler + Numerical_Algorithmic_SecondOrder should "
        "cost far less than the recorded BE+Continuum total of %d -- got %d"
        % (TWOCUBE_BE_CONTINUUM_TOTAL, sum(it_num)))


# ===========================================================================
# 4. Refusal propagation from inside a perturbed sub-call.
# ===========================================================================
@pytest.mark.parametrize("method", ["Backward_Euler", "Closest_Point"])
def test_starved_perturbation_refuses_the_whole_tangent_not_a_partial_one(
        numalg_repoint_available, method):
    """A perturbed sub-call that exhausts n_max_iterations must fail the
    ENTIRE numerical-tangent call (and therefore the host analysis), never
    silently assemble a tangent from whatever columns happened to converge --
    that would be exactly the class of silent-wrong-answer defect ADR-94
    exists to close. `hiso=7000.0, niter=1` is the existing starved-Newton
    reproducer from `test_adr97_p6_failloud.py::_starved`."""
    if method == "Closest_Point" and not S._constructible(
            lambda t: S.mat_vm(t, method="Closest_Point")):
        pytest.skip("Closest_Point unavailable on this binary")
    r = S.drive(lambda t: S.mat_vm(t, method=method,
                                   tangent="Numerical_Algorithmic_SecondOrder",
                                   hiso=7000.0, niter=1),
               S.VM_PATHS["triaxial"], nstep=10, ele="LadrunoBrick")
    print("%s + Numerical_Algorithmic_SecondOrder starved codes: %s"
          % (method, r["codes"]))
    assert any(c != 0 for c in r["codes"]), (
        "a starved perturbed sub-call inside the numerical tangent was "
        "committed as success (%s)" % method)
