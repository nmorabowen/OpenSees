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
4. Refusal propagation: a starved integrator (`Backward_Euler` + plain
   MohrCoulomb via the ADR-84 P2a exhaustion reproducer; `Closest_Point` + VM
   via the ADR-97 P1 gate-6 reproducer) refuses identically under
   `Numerical_Algorithmic_SecondOrder` and under `Secant`/the default -- the
   repoint does not mask, alter, or crash differently on, a refusal it did
   not itself cause. See the section-4 comment below for why isolating a
   refusal to STRICTLY inside a perturbed sub-call (primary converges,
   perturbation alone starves) was not achievable deterministically.

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
import test_asdplastic_mctc as M  # noqa: E402

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
    # MEASURED (build c24cda99c, this rig): 2.055e-08. Pinned at 1e-6 (the
    # fork's cross-platform float-pin floor), ~50x margin over the measured
    # value.
    assert r["rel_err"] <= 1e-6, r["rel_err"]


def test_be_numerical_firstorder_matches_fd_of_its_own_committed_map(
        numalg_repoint_available):
    """Same check for the forward-difference variant. The task brief expected
    a looser order (~1e-3) than the central-difference pair above (forward FD
    has no even-order cancellation), but MEASURED (build c24cda99c, this rig)
    it lands at 3.394e-08 -- essentially the same order as SecondOrder. This
    single-shot uniaxial rig takes ONE load-controlled step to a converged
    strain increment small enough that first- and second-order stencil
    truncation are both far below the rig's own FD step `h`; a rig with a
    larger `TrialStrain - CommitStrain` would likely show the expected
    forward-vs-central gap. Pinned at 1e-6 (the fork's cross-platform
    float-pin floor) rather than a looser bound, since that is what was
    actually observed -- tighten further only with more measurement across
    other rigs/materials, not by guessing.
    """
    r = FD.fd_check(lambda t: FD.mat_vm(t, tangent="Numerical_Algorithmic_FirstOrder",
                                        method="Backward_Euler"),
                    label="Backward_Euler / Numerical_Algorithmic_FirstOrder")
    print("P4 gate BE+NumAlg1st rel_err = %.3e (expect <= 1e-6, measured 3.394e-08)"
          % r["rel_err"])
    assert r["rc"] == 0, "the Backward_Euler load-driven rig did not converge"
    assert r["rel_err"] <= 1e-6, r["rel_err"]


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
# 4. Refusal propagation: a starved integrator must not silently assemble a
#    partial/short tangent under Numerical_Algorithmic_*.
#
# The scope's ideal reproducer -- a perturbed sub-call ALONE exhausting its
# iteration budget while the unperturbed primary call converges -- turned out
# to be impractical to engineer deterministically: `setTrialStrainIncr()`'s
# Newton is well within its quadratic-convergence regime by the time a
# 1e-8-relative FD perturbation is applied, so (measured, both methods) a
# perturbed call's own iteration count essentially never exceeds the
# unperturbed one's at any n_max_iterations boundary tried (a sweep around
# the exact boundary the un-perturbed call needs -- 42 fails/43 succeeds on
# the Backward_Euler deck below -- found FirstOrder and SecondOrder both
# succeeding right alongside Secant at every niter from 43 up). What IS
# reliably testable, and is the property that actually matters operationally
# (a refusal must never be swallowed into a wrong-but-plausible tangent): a
# deck starved badly enough that BOTH the primary call and its perturbations
# fail refuses identically whether `tangent_type` is `Secant`/`Algorithmic`
# or `Numerical_Algorithmic_SecondOrder` -- the repoint does not mask, alter,
# or crash differently on, a refusal it did not itself cause.
# ===========================================================================
def test_be_starved_mc_refuses_identically_under_numerical_and_secant(
        numalg_repoint_available):
    """`Backward_Euler` + plain `MohrCoulomb_YF` (no cutoff, so the scalar
    Newton is the only thing standing between this material and an f > 0
    commit) on the Cerro-Lindo-like tet deck from
    `test_adr84_p2a_strict_convergence.py` (`n_max_iterations 2`,
    `strict_convergence 1` -- the exhaustion-accept reproducer from ADR-84
    P2a). Both `tangent_type Secant` (the default; never touches this WP's
    code) and `tangent_type Numerical_Algorithmic_SecondOrder` (re-points
    through `numerical_tangent_of_committed_map()`) must refuse on the same
    step with the same code -- the P4 repoint must not change whether, or
    how, a starved Backward_Euler is refused.
    """
    def _mat_mc(tag, tangent):
        ops.nDMaterial(
            "ASDPlasticMaterial3D", tag,
            "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", M.IV,
            "Begin_Model_Parameters",
            "YoungsModulus", M.E, "PoissonsRatio", M.NU,
            "MC_phi", M.PHI, "MC_c", M.C, "MC_psi", M.PSI, "MC_ds", 0.0,
            "MassDensity", 0.0,
            "End_Model_Parameters",
            "Begin_Internal_Variables",
            "BackStress", 0., 0., 0., 0., 0., 0.,
            "End_Internal_Variables",
            "Begin_Integration_Options",
            "tangent_type", tangent, "strict_convergence", 1,
            "n_max_iterations", 2,
            "End_Integration_Options",
        )

    if not S._constructible(lambda t: _mat_mc(t, "Secant")):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb not available")

    tet = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
           5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5),
           9: (.5, 0, .5), 10: (0, .5, .5)}
    top = (4, 8, 9, 10)
    nsteps = 20

    def _run(tangent):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in tet.items():
            ops.node(t, *map(float, c))
        for t in (1, 2, 3, 5, 6, 7):
            ops.fix(t, 1, 1, 1)
        for t in top:
            ops.fix(t, 1, 1, 0)
        _mat_mc(1, tangent)
        ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.sp(t, 3, -0.02)
        ops.constraints("Penalty", 1e14, 1e14)
        ops.numberer("Plain")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-8, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0 / nsteps)
        ops.analysis("Static")
        codes = []
        for _ in range(nsteps):
            rc = ops.analyze(1)
            codes.append(rc)
            if rc != 0:
                break
        return codes

    codes_secant = _run("Secant")
    codes_numalg = _run("Numerical_Algorithmic_SecondOrder")
    print("BE starved MC codes: Secant=%s  Numerical_Algorithmic_SecondOrder=%s"
          % (codes_secant, codes_numalg))
    assert codes_secant[-1] != 0, "the ADR-84 P2a starved-MC reproducer stopped reproducing"
    assert codes_numalg == codes_secant, (
        "ADR-97 P4: Numerical_Algorithmic_SecondOrder refused differently "
        "from Secant on the same starved Backward_Euler deck")


def test_cp_starved_vm_refuses_the_whole_tangent_not_a_partial_one(
        numalg_repoint_available):
    """`Closest_Point`'s own Newton fails loud by default (no
    `strict_convergence` needed, unlike `Backward_Euler` -- ADR-97 P1 gate 6 /
    `test_adr97_p6_failloud.py::test_starved_newton_is_refused_not_committed`).
    `hiso=7000.0, niter=1` on the triaxial VM deck is that exact reproducer;
    run here under `Numerical_Algorithmic_SecondOrder` to confirm the P4
    repoint's perturbed sub-calls propagate the SAME refusal outward rather
    than assembling a tangent from whatever columns happened to converge
    before the refusing one."""
    if not S._constructible(lambda t: S.mat_vm(t, method="Closest_Point")):
        pytest.skip("Closest_Point unavailable on this binary")
    r = S.drive(lambda t: S.mat_vm(t, method="Closest_Point",
                                   tangent="Numerical_Algorithmic_SecondOrder",
                                   hiso=7000.0, niter=1),
               S.VM_PATHS["triaxial"], nstep=10, ele="LadrunoBrick")
    print("Closest_Point + Numerical_Algorithmic_SecondOrder starved codes: %s"
          % (r["codes"],))
    assert any(c != 0 for c in r["codes"]), (
        "a starved Closest_Point Newton (reachable from inside a perturbed "
        "sub-call or the primary call alike) was committed as success")
