"""ADR-84 P2(a) — the ``strict_convergence`` integration option, which fixes the
two UPSTREAM silent-accept defects in ``ASDPlasticMaterial3D::Backward_Euler``.

THE DEFECTS (both in the DEFAULT integrator, both affecting every ASDP material)
--------------------------------------------------------------------------------
1. **Exhaustion-accept.**  The scalar-Newton consistency loop is
   ``for (int iter = 0; iter < max_iter; ++iter) { ... if (|Phi| < tol) break; ... }``
   and falls out of ``max_iter`` with **no convergence check**, dropping straight
   through to ``ComputeTangentStiffness(); return 0;``.  A stalled Gauss point
   reports SUCCESS, and the global Newton then converges on a residual assembled
   from an inadmissible stress with nothing in the output saying so.  This is the
   persistence mechanism behind the Cerro Lindo ADR-0005 M3 finding (20 Gauss
   points at ``f/(2c cos phi) = +0.0299``).
2. **f-decreasing elastic early-exit.**  ``(yf_val_start - yf_val_end > tol_yf)``
   accepts a step as elastic on the sole grounds that f went DOWN, with no
   requirement that the end state be admissible.  From an admissible commit it
   cannot manufacture a violation; from an already-inadmissible one — exactly
   what defect 1 produces — it carries ``f > 0`` forward instead of correcting
   it.  The two compound.

THE FIX IS OPT-IN.  ``strict_convergence`` (int, default 0) is parsed in the
``Begin_Integration_Options`` block.  Default 0 is byte-identical to upstream
**by design**: turning it on changes convergence behaviour for every existing
VonMises / DruckerPrager / MohrCoulomb deck, so it cannot be the default in a
fork other people's models run on.

REPRODUCER — a STARVED ITERATION BUDGET, not the apex path
----------------------------------------------------------
The obvious candidate, plain ``MohrCoulomb_YF`` on a near-apex hydrostatic
tension path, does **not** reproduce on this build: measured 30/30 steps, worst
committed ``f_MC = 2.8e-14``, flag-on bit-identical.  (That path is still worth
keeping as ``test_asdplastic_mctc::test_plain_mc_apex_reference_run``; it simply
is not where the exhaustion-accept bites.)  What reproduces it deterministically
is starving ``n_max_iterations`` on a genuinely plastic leg — which is not a
contrivance but the defect stated precisely: *exhaustion is accepted*, and the
cheapest way to exhaust is a small budget.  Measured on the tet host below with
``n_max_iterations 2``: flag off completes 20/20 steps with **all 20 committed
states inadmissible** and worst ``f_MC = 77.6 kPa``; flag on refuses.

TWO MEASURED CAVEATS, both pinned by tests in this file
--------------------------------------------------------
* **``stdBrick`` swallows the refusal.**  ``Brick::update()`` writes
  ``success = materialPointers[i]->setTrialStrain(strain);`` and then
  unconditionally ``return 0;`` — ``success`` is assigned and never read.  So on
  a ``stdBrick`` host the material returns -1, prints its ``opserr`` line, and
  the analysis still reports success.  ``TenNodeTetrahedron`` accumulates and
  returns the sum (the fork already fixed that one — TIMs report item 8), so it
  is the host used for the contract tests here.  This is a THIRD silent-accept
  defect, in a vanilla element, and is deliberately NOT fixed in this PR:
  ``return success`` there is an unconditional behaviour change for every
  ``stdBrick`` + every material, which is exactly the blast radius the opt-in
  flag exists to avoid.  ``test_stdbrick_swallows_the_refusal`` pins it so it is
  not rediscovered.
* **``strict_convergence`` gates every global-Newton TRIAL, not just the commit.**
  A step's first global iterate can be far from the solution, and the local
  return map on that large trial increment may legitimately exhaust; flag-off
  lets it through and the step still converges.  So strict mode needs
  ``f_absolute_tol`` **scaled to the model's stress magnitude** — the option
  already exists for this.  Measured on the tet host (max |sigma| ~ 2.8e5 kPa)
  with ``n_max_iterations 100``: at ``f_absolute_tol`` 1e-6 / 1e-4 / 1e-2 strict
  mode refuses step 1 while flag-off runs clean, and at 1e-1 (a ~3.6e-7 RELATIVE
  demand) the two are bit-identical.  ``test_flag_on_is_inert_with_a_scaled_tolerance``
  pins that, and it is the operating guidance of record.

Parameters are the Cerro Lindo-like set of ``test_asdplastic_mctc`` (kPa), whose
oracles (``f_mc``, ``principals``) are imported rather than duplicated.
"""
import math

import numpy as np
import pytest

from _testbed import ops

import test_asdplastic_mctc as M

pytestmark = [pytest.mark.zone_a]

# Starved budget that makes exhaustion deterministic on a plastic leg.
NITER_STARVED = 2
# A budget the return map is comfortable with.
NITER_AMPLE = 100
# f_absolute_tol scaled to the tet rig's stress magnitude (see the docstring).
FTOL_SCALED = 1.0e-1

# VonMises control model (perfect plasticity), stresses O(sigma_y) = O(30 kPa),
# where the DEFAULT f_absolute_tol = 1e-6 is already a ~3e-8 relative demand.
E_VM, NU_VM, SY_VM = 70000.0, 0.3, 30.0
IV_VM = ("BackStress(TensorLinearHardeningFunction):"
         "YieldStress(ScalarLinearHardeningFunction):")


# ---------------------------------------------------------------------------
# materials — the integration-options block is what this file is about
# ---------------------------------------------------------------------------
def _int_opts(strict=None, niter=None, ftol=None):
    """``Begin_Integration_Options`` block, or NOTHING when every field is None
    (which is how every pre-P2a deck is spelled — the case flag-off must match).
    """
    out = []
    if strict is not None:
        out += ["strict_convergence", int(strict)]
    if niter is not None:
        out += ["n_max_iterations", int(niter)]
    if ftol is not None:
        out += ["f_absolute_tol", float(ftol)]
    if not out:
        return []
    return ["Begin_Integration_Options"] + out + ["End_Integration_Options"]


def mat_mc(tag, strict=None, niter=None, ftol=None):
    """Plain Mohr-Coulomb — no cutoff, so no ``special_return`` hook: the scalar
    Newton is the only thing standing between this material and an f > 0 commit.
    """
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
        *_int_opts(strict, niter, ftol),
    )


def mat_vm(tag, strict=None, niter=None, ftol=None):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV_VM,
        "Begin_Model_Parameters",
        "YoungsModulus", E_VM, "PoissonsRatio", NU_VM,
        "ScalarLinearHardeningParameter", 0.0,
        "TensorLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "YieldStress", SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *_int_opts(strict, niter, ftol),
    )


def f_vm(s, sy=SY_VM):
    """``||dev(sigma)|| - sqrt(2/3)*sigma_y`` — the arithmetic of
    ``VonMises_YF::YF`` with zero backstress (both hardening parameters are 0).
    """
    S = M._tensor(s)
    dev = S - (np.trace(S) / 3.0) * np.eye(3)
    return math.sqrt(float(np.sum(dev * dev))) - math.sqrt(2.0 / 3.0) * sy


# ---------------------------------------------------------------------------
# host A: TenNodeTetrahedron — PROPAGATES the material's return code, so the
# refusal is observable from Python.  Compression of a unit tet, base fixed,
# apex nodes driven down: a plainly plastic leg (max |sigma| ~ 2.8e5 kPa).
# ---------------------------------------------------------------------------
_TET = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
        5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5),
        9: (.5, 0, .5), 10: (0, .5, .5)}
_TET_TOP = (4, 8, 9, 10)
TET_NSTEPS = 20
TET_UTOP = -0.02


def _tet_build(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _TET.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    for t in _TET_TOP:
        ops.fix(t, 1, 1, 0)
    mat_fn(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in _TET_TOP:
        ops.sp(t, 3, TET_UTOP)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / TET_NSTEPS)
    ops.analysis("Static")


def _drive(build_fn, nsteps):
    """Step WITHOUT asserting convergence — a refused step is the point.

    Returns ``(codes, hist)``: every analyze code up to and including the first
    failure, and the committed GP-1 stress after each SUCCESSFUL step.  A failed
    step contributes a code and no history row: its state is never committed,
    which is exactly the property under test.
    """
    build_fn()
    codes, hist = [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")
        hist.append(list(ops.eleResponse(1, "stresses"))[0:6])
    return codes, (np.array(hist) if hist else np.zeros((0, 6)))


def _tet_run(strict=None, niter=None, ftol=None, mat=mat_mc):
    return _drive(lambda: _tet_build(lambda t: mat(t, strict, niter, ftol)),
                  TET_NSTEPS)


def _audit(hist, tol):
    """max f_MC over a committed history (nan on an empty history)."""
    if not len(hist):
        return float("nan")
    return float(np.max([M.f_mc(s) for s in hist]))


def _cond_tol(hist):
    """Conditioning-aware admissibility tolerance — f is O(|sigma|), so a fixed
    absolute bound is a meaningless relative demand on a 1e5 kPa state
    (test_asdplastic_mctc's rule, reused)."""
    smax = float(np.max(np.abs(hist))) if len(hist) else 0.0
    return 1.0e-6 * max(2.0 * M.C * math.cos(math.radians(M.PHI)), smax)


# ---------------------------------------------------------------------------
# availability guards
# ---------------------------------------------------------------------------
def _tet_constructible(mat_fn):
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in _TET.items():
            ops.node(t, *map(float, c))
        mat_fn(1)
        ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def mc_available():
    if not _tet_constructible(lambda t: mat_mc(t)):
        pytest.skip("ASDPlasticMaterial3D / MohrCoulomb_YF / "
                    "TenNodeTetrahedron not available in this build")


@pytest.fixture(scope="module")
def vm_available():
    if not _tet_constructible(lambda t: mat_vm(t)):
        pytest.skip("ASDPlasticMaterial3D / VonMises_YF not available in this "
                    "build")


# ===========================================================================
# 1. flag OFF is inert — the hard gate
# ===========================================================================
@pytest.mark.t0m
def test_flag_off_is_inert(mc_available):
    """Omitting the option entirely and passing ``strict_convergence 0`` must
    produce the SAME run, bit for bit — on the very path where the flag DOES
    change behaviour when set to 1 (test 3), so this is not a vacuous
    comparison of two clean runs.

    This is what makes the option safe to add to a shared fork: a deck that does
    not mention it and a deck that mentions it as 0 are the same deck.  (That
    flag-off also equals *upstream* is carried by the unchanged ASDP batteries,
    which know nothing about this option.)
    """
    codes_a, hist_a = _tet_run(strict=None, niter=NITER_STARVED)
    codes_b, hist_b = _tet_run(strict=0, niter=NITER_STARVED)

    assert codes_a == codes_b, (
        f"parsing 'strict_convergence 0' changed the analyze codes: {codes_a} "
        f"(option absent) vs {codes_b} (option = 0)")
    assert hist_a.shape == hist_b.shape, (
        f"parsing 'strict_convergence 0' changed the number of committed "
        f"steps: {hist_a.shape[0]} vs {hist_b.shape[0]}")
    if not np.array_equal(hist_a, hist_b):
        d = np.abs(hist_a - hist_b)
        i, j = np.unravel_index(int(np.argmax(d)), d.shape)
        pytest.fail(
            f"'strict_convergence 0' is NOT byte-identical to the option being "
            f"absent: max|diff| = {d.max():.6e} at step {i} comp {j} "
            f"({hist_a[i, j]!r} vs {hist_b[i, j]!r}). Flag-off must be inert.")


# ===========================================================================
# 2. flag OFF commits inadmissible states — the defect, measured
# ===========================================================================
@pytest.mark.t0m
def test_flag_off_commits_inadmissible_states(mc_available):
    """The upstream behaviour this option exists to make refusable, pinned.

    With a starved budget the return map cannot converge, and flag-off reports
    SUCCESS anyway: measured 20/20 steps "converged", every committed state
    outside the yield surface, worst ``f_MC`` ~ 77.6 kPa against an
    admissibility tolerance of ~2.8e-1.

    If this test ever fails because the run came back CLEAN, the defect was
    fixed by default somewhere — update ADR-84 §6b and delete this test rather
    than "repairing" it.  It asserts a bug on purpose; that is what makes
    test 3 non-vacuous.
    """
    codes, hist = _tet_run(strict=None, niter=NITER_STARVED)

    assert codes == [0] * TET_NSTEPS, (
        f"flag-off was expected to report SUCCESS on every step (that is the "
        f"defect); got codes = {codes}")

    tol = _cond_tol(hist)
    fmc = np.array([M.f_mc(s) for s in hist])
    n_bad = int(np.sum(fmc > tol))
    assert n_bad > 0, (
        f"flag-off committed NO inadmissible state (max f_MC = {fmc.max():.6e} "
        f"<= tol {tol:.3e}) — the starved-budget reproducer no longer "
        f"reproduces, so test 3 would be vacuous. Re-tune "
        f"n_max_iterations/the drive before trusting anything in this file.")

    print(f"\n[flag OFF, n_max_iterations={NITER_STARVED}] "
          f"steps completed = {len(hist)}/{TET_NSTEPS} (all reported success), "
          f"inadmissible committed states = {n_bad}/{len(hist)}, "
          f"worst committed f_MC = {fmc.max():.6e} vs tol {tol:.3e} "
          f"(= {fmc.max() / tol:.1f}x). This is the upstream silent accept.")


# ===========================================================================
# 3. flag ON refuses instead of committing — THE contract
# ===========================================================================
@pytest.mark.t0m
def test_flag_on_refuses_instead_of_committing(mc_available):
    """With ``strict_convergence 1`` the same starved run must either converge
    every step or FAIL a step — and either way NO COMMITTED STATE may carry
    ``f > tol``.

    This is the Cerro Lindo M3 acceptance ("no committed state may carry f > 0
    anywhere in the model") reduced to one material point and enforced BY the
    material instead of by a post-hoc audit.
    """
    codes_off, hist_off = _tet_run(strict=None, niter=NITER_STARVED)
    codes_on, hist_on = _tet_run(strict=1, niter=NITER_STARVED)

    # -- the contract: audit EVERY committed state in numpy
    tol = _cond_tol(hist_on) if len(hist_on) else _cond_tol(hist_off)
    if len(hist_on):
        fmc = np.array([M.f_mc(s) for s in hist_on])
        k = int(np.argmax(fmc))
        assert fmc[k] <= tol, (
            f"strict_convergence=1 STILL committed an inadmissible state: max "
            f"f_MC = {fmc[k]:.6e} > tol {tol:.3e} at committed step {k} "
            f"(0-based) of {len(hist_on)}; sigma = {hist_on[k]}, principals = "
            f"{M.principals(hist_on[k])}. analyze codes = {codes_on}")

    # -- it must have REFUSED: flag-off ran all 20 steps dirty (test 2), so a
    #    flag-on run that also reports 20 clean successes would mean the gate
    #    never fired.
    refused = any(c != 0 for c in codes_on)
    assert refused, (
        f"strict_convergence=1 did NOT refuse a step on a path where flag-off "
        f"commits {int(np.sum([M.f_mc(s) > tol for s in hist_off]))} "
        f"inadmissible states: codes = {codes_on}. The gate did not fire.")

    assert codes_on != codes_off, (
        f"strict_convergence=1 produced the same analyze codes as flag-off "
        f"({codes_on}) — the flag is a no-op on this path")

    print(f"\n[flag ON, n_max_iterations={NITER_STARVED}] "
          f"analyze codes = {codes_on} (refused at step {len(codes_on)}), "
          f"committed steps = {len(hist_on)} vs {len(hist_off)} flag-off, "
          f"worst committed f_MC = {_audit(hist_on, tol):.6e} (tol {tol:.3e}). "
          f"The non-converged state is refused, not committed.")


# ===========================================================================
# 4. flag ON must not perturb a path that already converges
# ===========================================================================
@pytest.mark.t0m
def test_flag_on_does_not_perturb_a_converging_path(vm_available):
    """A well-behaved VonMises leg must be BIT-IDENTICAL flag-on vs flag-off.

    The flag may only ever ADD a refusal on a path that was not converging; if
    it moves a converged answer by even one ulp it is not a convergence gate,
    it is a second constitutive model.  VonMises here carries stresses of order
    ``sigma_y`` = 30 kPa, where the default ``f_absolute_tol`` = 1e-6 is already
    a ~3e-8 relative demand — i.e. the tolerance is comfortably reachable, which
    is the regime the flag is meant to be silent in.
    """
    codes_off, hist_off = _tet_run(strict=0, mat=mat_vm)
    codes_on, hist_on = _tet_run(strict=1, mat=mat_vm)

    assert codes_off == [0] * TET_NSTEPS, (
        f"the VonMises control leg does not converge even with the flag OFF "
        f"(codes {codes_off}) — it cannot serve as the no-perturbation control")

    # non-vacuity: it must actually be PLASTIC
    fvm = np.array([f_vm(s) for s in hist_off])
    assert float(np.max(fvm)) > -1.0e-3 * SY_VM, (
        f"the VonMises control leg never reached yield (max f_VM = "
        f"{np.max(fvm):.6e}, sigma_y = {SY_VM}) — a purely elastic path cannot "
        f"show that a convergence gate is inert")

    assert codes_on == codes_off, (
        f"strict_convergence=1 changed the analyze codes on a CONVERGING "
        f"VonMises path: {codes_off} -> {codes_on}")
    assert hist_on.shape == hist_off.shape, (
        f"strict_convergence=1 changed the number of committed steps on a "
        f"converging VonMises path: {hist_off.shape[0]} -> {hist_on.shape[0]}")
    if not np.array_equal(hist_on, hist_off):
        d = np.abs(hist_on - hist_off)
        i, j = np.unravel_index(int(np.argmax(d)), d.shape)
        scale = max(1.0, float(np.max(np.abs(hist_off))))
        pytest.fail(
            f"strict_convergence=1 PERTURBED a converging VonMises path: "
            f"max|diff| = {d.max():.6e} at step {i} comp {j} "
            f"({hist_on[i, j]!r} vs {hist_off[i, j]!r}), ~"
            f"{d.max() / (scale * 2.220446049250313e-16):.1f} ulp of the "
            f"stress scale. The flag must only ever REFUSE a non-converged "
            f"step, never change a converged one.")


# ===========================================================================
# 5. the operating guidance: strict mode needs a SCALED f_absolute_tol
# ===========================================================================
@pytest.mark.t0m
def test_flag_on_is_inert_with_a_scaled_tolerance(mc_available):
    """``strict_convergence`` gates every global-Newton TRIAL evaluation, not
    only the committed state — so on a high-stress model the default
    ``f_absolute_tol`` = 1e-6 (an absolute bound in stress units) is an
    unreachable relative demand and strict mode refuses step 1 of a run that is
    perfectly healthy with the flag off.

    Measured on this rig (max |sigma| ~ 2.8e5 kPa, ``n_max_iterations`` 100):
    ``f_absolute_tol`` 1e-6 / 1e-4 / 1e-2 all refuse at step 1; 1e-1 — a ~3.6e-7
    RELATIVE demand — is bit-identical to flag-off.  That is the operating
    guidance of record: **scale f_absolute_tol to your stress magnitude before
    turning strict_convergence on.**  This test pins both halves.
    """
    # (a) the trap: default tolerance, ample budget, healthy flag-off run
    codes_off, hist_off = _tet_run(strict=0, niter=NITER_AMPLE)
    codes_on, _hist_on = _tet_run(strict=1, niter=NITER_AMPLE)
    assert codes_off == [0] * TET_NSTEPS, (
        f"the ample-budget control run does not converge flag-off: {codes_off}")
    assert any(c != 0 for c in codes_on), (
        f"the documented over-refusal at the DEFAULT f_absolute_tol no longer "
        f"happens (codes {codes_on}) — the caveat in this module's docstring "
        f"and in ADR-84 P2a is stale and must be re-measured")

    # (b) the remedy: a tolerance scaled to the stress magnitude is inert
    codes_off_s, hist_off_s = _tet_run(strict=0, niter=NITER_AMPLE,
                                       ftol=FTOL_SCALED)
    codes_on_s, hist_on_s = _tet_run(strict=1, niter=NITER_AMPLE,
                                     ftol=FTOL_SCALED)

    assert codes_on_s == codes_off_s == [0] * TET_NSTEPS, (
        f"with f_absolute_tol = {FTOL_SCALED} the run should complete both "
        f"ways: flag-off {codes_off_s}, flag-on {codes_on_s}")
    assert hist_on_s.shape == hist_off_s.shape and np.array_equal(hist_on_s,
                                                                 hist_off_s), (
        f"with a stress-scaled f_absolute_tol = {FTOL_SCALED}, "
        f"strict_convergence=1 must be bit-identical to flag-off; it is not")

    smax = float(np.max(np.abs(hist_off_s)))
    print(f"\n[scaled tolerance] max |sigma| = {smax:.4e} kPa; "
          f"f_absolute_tol = {FTOL_SCALED:.1e} is a {FTOL_SCALED / smax:.2e} "
          f"RELATIVE demand -> flag-on bit-identical to flag-off over "
          f"{len(hist_on_s)} steps. At the default 1e-6 "
          f"({1e-6 / smax:.2e} relative) flag-on refuses at step 1.")


# ===========================================================================
# 6. the third silent accept: stdBrick swallows the material's refusal
# ===========================================================================
@pytest.mark.t0m
def test_stdbrick_swallows_the_refusal(mc_available):
    """``Brick::update()`` writes ``success = ...->setTrialStrain(strain);`` and
    then unconditionally ``return 0;`` — the code is assigned and never read.

    So on a ``stdBrick`` host the material prints its refusal and returns -1 and
    the analysis still reports success.  This is a THIRD silent-accept defect,
    in a vanilla ELEMENT, found while measuring P2a.  It is deliberately not
    fixed here: ``return success`` is an unconditional behaviour change for
    every stdBrick + every material, which is precisely the blast radius the
    opt-in flag exists to avoid.  (``TenNodeTetrahedron`` already accumulates
    and returns the sum — fork fix, TIMs report item 8; that is why the contract
    tests above use the tet.)

    Pinned here so it is not rediscovered, and so that FIXING Brick.cpp shows up
    as a loud, informative failure of this test rather than as a mystery
    elsewhere.
    """
    ax = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -5.0e-3])
    lat = ([0.0, 0.2, 1.0], [0.0, -1.0e-4, -1.0e-4])
    nsteps = 50

    codes, hist = _drive(
        lambda: M._build(lambda t: mat_mc(t, 1, NITER_STARVED),
                         ax, lat, 1.0, 1.0 / nsteps),
        nsteps)

    assert codes == [0] * nsteps, (
        f"stdBrick now PROPAGATES the material's refusal (codes {codes}). If "
        f"Brick::update() was fixed to 'return success', that is good news — "
        f"delete this test, note it in ADR-84 P2a and LEDGER_quirks, and move "
        f"the contract tests back onto stdBrick.")

    tol = _cond_tol(hist)
    fmc = np.array([M.f_mc(s) for s in hist])
    n_bad = int(np.sum(fmc > tol))
    assert n_bad > 0, (
        f"expected stdBrick to commit inadmissible states despite the material "
        f"refusing every one of them (max f_MC = {fmc.max():.6e} <= tol "
        f"{tol:.3e}) — re-measure before trusting this pin")

    print(f"\n[stdBrick swallow] strict_convergence=1, "
          f"n_max_iterations={NITER_STARVED}: the material refused (see the "
          f"opserr lines) yet analyze() returned 0 on all {nsteps} steps and "
          f"{n_bad} committed states carry f_MC > tol (worst "
          f"{fmc.max():.6e} vs {tol:.3e}). Brick::update() drops the code.")
