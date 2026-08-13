"""ADR-84 P3 — ``MohrCoulombTensionCutoff`` in a CONFINED Mohr-Coulomb shear
state on the tension-cutoff face and the MC∩TC corner: the regime the composite
was requested for, and the one the shipped P0 verification never visited.

PROVENANCE
----------
Cerro Lindo SSI, rung **M5** (13 Aug 2026).  Their M3/M4 model runs to λ = 1
with plain MC; substituting ``MohrCoulombTensionCutoff_YF/PF`` and changing
nothing else stalls stage S1 (pure continuum) at λ ≈ 0.45-0.50 with a ~9e-5 m
displacement increment against a ~500 kN residual — on all three settings they
tried, including ``tangent_type = Continuum``.

The P0 battery (``test_asdplastic_mctc``) exercises the composite in **uniaxial
tension on a single stdBrick**, and its material-point driver prescribes every
DOF — zero free equations, so the global Newton never solves anything and the
tangent the material hands the assembler is unobservable.  Both gaps are why a
tangent defect could ship.  This module closes them.

WHAT WAS ACTUALLY WRONG (measured; ADR §9)
------------------------------------------
``special_return`` hardcoded ``secant_blend`` — ``(E + D)/2`` — on **every**
rung.  Two consequences, both fixed in P3:

1. **It silently overrode ``tangent_type``.**  Every Gauss point resolved by the
   hook got the blend no matter what the deck asked for.  The consumer's
   ``tangent_type = Continuum`` was inert on exactly the Gauss points that
   needed it — which is why their three settings all behaved the same.
2. **The blend fabricates stiffness that is not there.**  At the corner the true
   algorithmic tangent of the hook's own return is rank-deficient
   (``min|eig| ≈ 9e-4``, cond 6.8e9): the corner is a rigid attractor, the
   stress does not move for a 20x range of strain increment.  The blend reports
   ``min|eig| = 2.0e6 kPa`` in that same direction — 87% error in norm, and a
   *well-conditioned* matrix (cond 3.2) where the truth is nearly singular.
   For contrast, plain MC driven to a corner-adjacent state (on the isotropic
   in-situ deconfinement path of the numpy replica, not this driver) correctly
   hands out a rank-deficient tangent, cond ~1e16.  A global Newton told "push
   here and stress rises by 2e6·Δε", whose stress then does not move at all,
   produces precisely the reported signature: a vanishing increment against a
   residual that never falls.

Note the defect is **not corner-specific**: the Rankine FACE rung is blended too,
and on this path the face steps show the same 7-vs-2 iteration penalty.  The M5
model has far more face points than corner points.

P3 moves the tangent-operator policy to the integrator, where the generic path
already applies it.  ``Secant`` stays the default and is byte-identical;
``Continuum`` now actually delivers the raw active-set (Koiter) tangent.

DRIVER
------
Unit-cube ``stdBrick``, 1/8-symmetry restraints.  x-faces driven in COMPRESSION
and y-faces in TENSION by independent ``Path`` series; **the z-faces are left
FREE**.  That last part is the whole point: it leaves 4 real unknowns (the
z-DOFs at z = 1) whose stiffness comes from the material tangent, so tangent
quality is visible in ``testIter()``.  σ_zz = 0 throughout.

The path walks the Rankine FACE up to the MC surface and then sits on the
MC∩TC CORNER.  With the Cerro Lindo host rock (c = 1014 kPa, φ = 45.95°,
ψ = 11.49° — strongly non-associated — T = 24.7 kPa) the corner is at

    s1 = -4862.2 kPa,   s2 = 0,   s3 = T = +24.7 kPa

i.e. **4.86 MPa of confinement**, inside the σ₃ ≈ 1-6 MPa band the M5 report
asked for.  Measured on this build: steps 0-29 on the FACE, step 30 enters the
CORNER, steps 31-44 sit on it as a perfect attractor (stress frozen to the last
digit, so those steps converge in 1 iteration for every tangent type and carry
no information — the discriminating steps are the 31 that precede them).

THE GATE
--------
``test_tangent_type_is_not_ignored`` is the regression.  Before P3 the hook
assigned ``Stiffness = stiff_sr`` unconditionally and ``stiff_sr`` was always
secant-blended, so on this path — where every plastic step is hook-resolved —
``Secant``, ``Continuum`` and ``Elastic`` produced **identical** iteration
histories.  Identical histories are now a failure.

Measured on the fixed build: Secant 7 iterations/step, Continuum 2, Elastic 9.

Tolerances are conditioning-aware: |σ| ~ 5e3 kPa here, so the admissibility
audit uses the ``test_asdplastic_mctc`` rule (1e-6 scaled by
``max(2c·cosφ, max|σ|)``) rather than a fixed absolute.
"""
import math

import numpy as np
import pytest

from _testbed import ops

import test_asdplastic_mctc as M

pytestmark = [pytest.mark.zone_a]

# ---------------------------------------------------------------------------
# Cerro Lindo HOST ROCK (kPa) — the M5 parameters, not the softer EDZ set the
# P0 battery uses.  These are what put the corner at ~4.9 MPa of confinement.
# ---------------------------------------------------------------------------
E, NU = 2.0e6, 0.3
C, PHI, PSI = 1014.0, 45.95, 11.49          # psi << phi: strongly non-associated
T_ROCK = 24.7

IV = "BackStress(NullHardeningTensorFunction):"

# Closed form of the corner reached by this driver (s2 = 0 because z is free):
# the s1 solving f_MC(s1, 0, T) = 0.
S1_CORNER = -4862.2

NSTEPS = 45
EPS_XX_END = -3.6e-3        # compression: drives the confinement
EPS_YY_END = +6.0e-3        # tension: drives s3 onto the cutoff, then the corner


# ---------------------------------------------------------------------------
# material — one knob: tangent_type
# ---------------------------------------------------------------------------
def mat_mctc(tag, tangent_type=None, strict=None):
    opts = []
    if tangent_type is not None:
        opts += ["tangent_type", tangent_type]
    if strict is not None:
        opts += ["strict_convergence", int(strict)]
    block = (["Begin_Integration_Options"] + opts + ["End_Integration_Options"]
             if opts else [])
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF",
        "LinearIsotropic3D_EL", IV,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", NU,
        "MC_phi", PHI, "MC_c", C, "MC_psi", PSI, "MC_ds", 0.0,
        "TC_min_stress", T_ROCK, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *block,
    )


def mat_mc(tag, tangent_type=None):
    """Plain MC on the identical driver — the run the M5 model completes."""
    opts = ["Begin_Integration_Options", "tangent_type", tangent_type,
            "End_Integration_Options"] if tangent_type else []
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "MohrCoulomb_YF", "MohrCoulomb_PF", "LinearIsotropic3D_EL", IV,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", NU,
        "MC_phi", PHI, "MC_c", C, "MC_psi", PSI, "MC_ds", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "BackStress", 0., 0., 0., 0., 0., 0.,
        "End_Internal_Variables",
        *opts,
    )


# ---------------------------------------------------------------------------
# driver: x prescribed (compression), y prescribed (tension), z FREE
# ---------------------------------------------------------------------------
def _build(mat_fn, **mat_kw):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in M._CUBE.items():
        ops.node(t, *map(float, c))
    for t, m in M._FIX.items():
        ops.fix(t, *m)
    mat_fn(1, **mat_kw)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    M._path(1, [0.0, 1.0], [0.0, EPS_XX_END], 1.0)
    ops.pattern("Plain", 1, 1)
    for n in M._XFACE:
        ops.sp(n, 1, 1.0)

    M._path(2, [0.0, 1.0], [0.0, EPS_YY_END], 1.0)
    ops.pattern("Plain", 2, 2)
    for n in M._YFACE:
        ops.sp(n, 2, 1.0)

    # z-faces deliberately NOT prescribed => sigma_zz = 0, and the four z-DOFs
    # at z = 1 are genuine unknowns.  Without them testIter() is identically 1
    # and the tangent is unobservable — which is how P0 shipped.
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")


def _drive(mat_fn, nsteps=NSTEPS, **mat_kw):
    """Run the path WITHOUT asserting convergence (a refused step is data).

    Returns ``(codes, sig, iters)``: the analyze code of every attempted step,
    the committed GP-1 stress of every SUCCESSFUL step, and the Newton iteration
    count each successful step needed.
    """
    _build(mat_fn, **mat_kw)
    codes, sig, iters = [], [], []
    for _ in range(nsteps):
        rc = ops.analyze(1)
        codes.append(rc)
        if rc != 0:
            break
        ops.eleResponse(1, "forces")
        sig.append(list(ops.eleResponse(1, "stresses"))[0:6])
        iters.append(ops.testIter())
    return codes, np.array(sig), np.array(iters, dtype=int)


def _f_mc(s):
    return M.f_mc(s, phi_deg=PHI, c=C)


def _f_tc(s):
    return M.f_tc(s, T_ROCK)


def _scale(sig):
    return max(2.0 * C * math.cos(math.radians(PHI)), float(np.max(np.abs(sig))))


def _corner_mask(sig):
    """Steps on the MC∩TC corner: cutoff active AND MC surface active."""
    sc = _scale(sig)
    on_cut = np.array([abs(_f_tc(s)) <= 1e-6 * sc for s in sig])
    on_mc = np.array([abs(_f_mc(s)) <= 1e-6 * sc for s in sig])
    return on_cut & on_mc


def _assert_admissible(sig, label):
    tol = 1e-6 * _scale(sig)
    worst = max(max(_f_mc(s), _f_tc(s)) for s in sig)
    assert worst <= tol, (
        f"{label}: committed an INADMISSIBLE state, max f = {worst:.6g} > tol "
        f"{tol:.6g} (|sigma|max = {np.max(np.abs(sig)):.6g} kPa)")


@pytest.fixture(scope="module")
def mctc_available():
    if not M._constructible(M.mat_mctc):
        pytest.skip("ASDPlasticMaterial3D MohrCoulombTensionCutoff_YF/_PF not "
                    "available in this build")


@pytest.fixture(scope="module")
def runs(mctc_available):
    """The three tangent policies on the identical confined path."""
    out = {}
    for tt in ("Secant", "Continuum", "Elastic"):
        out[tt] = _drive(mat_mctc, tangent_type=tt)
    return out


# ===========================================================================
# 0. the path really does reach a confined MC-shear corner — coverage, measured
# ===========================================================================
@pytest.mark.t0m
def test_path_reaches_a_confined_mc_shear_corner(runs):
    """The state the M5 report says was never tested: strongly confined,
    non-associated, ON the MC-cutoff corner.  Everything else here is only
    meaningful if this holds, so it is asserted rather than assumed."""
    codes, sig, _ = runs["Continuum"]
    assert all(c == 0 for c in codes), (
        f"the confined-corner path did not run to completion: codes {codes}")
    assert len(sig) == NSTEPS

    # the cutoff is active on essentially the whole path (Rankine FACE rung)
    on_cut = np.array([abs(_f_tc(s)) <= 1e-6 * _scale(sig) for s in sig])
    assert on_cut.sum() >= 30, (
        f"only {on_cut.sum()} steps had the tension cutoff active - the "
        f"reproducer has drifted off the cutoff entirely")

    mask = _corner_mask(sig)
    n_corner = int(mask.sum())
    assert n_corner >= 5, (
        f"path spent only {n_corner} steps on the MC-TC corner - the "
        f"reproducer no longer covers the reported regime")

    s_corner = sig[mask]
    s1 = np.array([M.principals(s)[0] for s in s_corner])
    s3 = np.array([M.principals(s)[2] for s in s_corner])

    assert np.allclose(s3, T_ROCK, atol=1e-6 * abs(S1_CORNER)), (
        f"corner steps are not on the cutoff: s3 = {s3}")
    assert np.allclose(s1, S1_CORNER, rtol=5e-3), (
        f"corner s1 = {s1} does not match the closed form {S1_CORNER}")
    # the confinement band the M5 report asked for
    assert 1000.0 <= np.abs(s1).min() and np.abs(s1).max() <= 6000.0, (
        f"confinement {np.abs(s1).min()/1000:.2f}-{np.abs(s1).max()/1000:.2f} "
        f"MPa is outside the 1-6 MPa band the report asked for")

    _assert_admissible(sig, "confined-corner path")
    print(f"\n[coverage] {on_cut.sum()}/{NSTEPS} steps with the cutoff active, "
          f"{n_corner} on the MC-TC corner at {np.abs(s1).mean()/1000:.2f} MPa "
          f"confinement (s1 = {s1.mean():.1f} kPa)")


# ===========================================================================
# 1. THE REGRESSION — tangent_type must not be silently overridden
# ===========================================================================
@pytest.mark.t0m
def test_tangent_type_is_not_ignored(runs):
    """Before P3 the hook hardcoded ``secant_blend`` on every rung and the
    integrator assigned it unconditionally, so on this path — where every
    plastic step is hook-resolved — all three tangent policies produced
    IDENTICAL Newton histories.  The knob the M5 report reached for was inert
    on exactly the Gauss points that needed it.

    The three policies must now differ.  The committed STRESSES must not: the
    tangent decides how you get to the solution, never what it is.
    """
    codes_s, sig_s, it_s = runs["Secant"]
    codes_c, sig_c, it_c = runs["Continuum"]
    codes_e, sig_e, it_e = runs["Elastic"]

    for tag, codes in (("Secant", codes_s), ("Continuum", codes_c),
                       ("Elastic", codes_e)):
        assert all(c == 0 for c in codes), f"{tag} run failed: {codes}"

    assert not np.array_equal(it_s, it_c), (
        "tangent_type is being IGNORED: 'Secant' and 'Continuum' needed "
        f"identical Newton iterations on every step ({it_s.tolist()}). This is "
        "the ADR-84 P0 defect — special_return hardcoded the secant blend and "
        "the integrator assigned it regardless of the configured operator.")
    assert not np.array_equal(it_e, it_c), (
        "tangent_type is being IGNORED: 'Elastic' and 'Continuum' needed "
        f"identical Newton iterations on every step ({it_e.tolist()})")

    sc = _scale(sig_s)
    assert np.allclose(sig_s, sig_c, atol=1e-6 * sc), (
        "Secant and Continuum converged to DIFFERENT stress states; the "
        "tangent must affect the path to the solution, not the solution")
    assert np.allclose(sig_s, sig_e, atol=1e-6 * sc), (
        "Elastic and Secant converged to DIFFERENT stress states")

    for tag, s in (("Secant", sig_s), ("Continuum", sig_c), ("Elastic", sig_e)):
        _assert_admissible(s, tag)

    print(f"\n[tangent knob] total Newton iterations: "
          f"Secant {it_s.sum()}, Continuum {it_c.sum()}, Elastic {it_e.sum()}")


# ===========================================================================
# 2. the consequence that matters — the raw tangent is the one that converges
# ===========================================================================
@pytest.mark.t0m
def test_raw_tangent_converges_faster_than_the_blend(runs):
    """On the cutoff face and at the corner the true tangent is rank-deficient.
    ``Continuum`` reports that; ``Secant`` blends half the elastic stiffness
    back into the degenerate direction, which is what turns a global Newton
    into the M5 stall.

    Measured on the fixed build: 2 iterations/step vs 7.  The gate demands only
    a factor of 2, so it pins the mechanism without pinning the exact build.
    """
    _, sig_s, it_s = runs["Secant"]
    _, _, it_c = runs["Continuum"]

    # Steps where the corner attractor has locked the stress converge in one
    # iteration for every operator and carry no information — exclude them.
    active = it_s > 1
    assert active.sum() >= 10, (
        f"only {active.sum()} steps did real Newton work; the driver no longer "
        f"discriminates between tangent operators")

    tot_s = int(it_s[active].sum())
    tot_c = int(it_c[active].sum())
    assert tot_c * 2 <= tot_s, (
        f"the raw (Continuum) tangent did not beat the secant blend by the "
        f"expected margin on the {active.sum()} working steps: Continuum "
        f"{tot_c} iterations vs Secant {tot_s}. The active-set tangent is "
        f"rank-deficient in truth; a well-conditioned blend of it cannot be "
        f"the better Newton operator there.")
    print(f"\n[convergence] over {active.sum()} working steps: Continuum "
          f"{tot_c} iterations vs Secant {tot_s} "
          f"({tot_s/max(tot_c,1):.1f}x fewer)")


# ===========================================================================
# 3. the M5 headline symptom must not reproduce
# ===========================================================================
@pytest.mark.t0m
def test_composite_completes_the_confined_path(mctc_available):
    """The M5 headline is that the composite STALLS partway through a run.  On
    this driver it must complete every step.

    Plain MC is reported alongside as context, NOT asserted on: it fails this
    particular driver around step 19 for a reason that has nothing to do with
    the cutoff.  With no cutoff the MC envelope runs into the tensile quadrant
    to c*cot(phi) = 980.9 kPa, so the tensile leg walks plain MC into its
    hydrostatic apex, where the upstream scalar Newton is known to break down
    (ADR-84 section 2 -- the defect this material exists to avoid).  The
    composite completing a path plain MC cannot is the cutoff working as
    intended, so pinning plain MC's failure here would pin upstream behaviour
    we do not own.
    """
    codes_tc, sig_tc, it_tc = _drive(mat_mctc, tangent_type="Continuum")
    assert all(c == 0 for c in codes_tc), (
        f"the composite stalled on the confined-corner path - the M5 headline "
        f"symptom, reproduced: codes {codes_tc}")
    assert len(sig_tc) == NSTEPS
    _assert_admissible(sig_tc, "composite")

    codes_mc, sig_mc, _ = _drive(mat_mc, tangent_type="Continuum")
    n_mc = sum(1 for c in codes_mc if c == 0)
    tail = "(also completes)" if n_mc == NSTEPS else "(walks into its own MC apex)"
    print(f"\n[context] composite {len(sig_tc)}/{NSTEPS} steps; "
          f"plain MC (no cutoff) {n_mc}/{NSTEPS} {tail}")


# ===========================================================================
# 4. the conservative vertex projection is refusable (M5 request #3)
# ===========================================================================
@pytest.mark.t0m
def test_strict_mode_is_inert_on_an_exactly_resolved_path(mctc_available):
    """M5 asked for "a documented, loud refusal when the integrator cannot land
    on the corner, rather than a Newton stall that reads as a modelling
    problem".

    ``special_return``'s terminal rung replaces the whole deviatoric state with
    ``T_eff·δ`` when no exact cutoff feature validated.  It preserves the
    "never commit f > 0" guarantee, but it is a large unphysical stress drop.
    P3 reports that rung as ``SR_QUALITY_FALLBACK`` and the integrator refuses
    it under ``strict_convergence`` instead of committing it silently.

    Every step of THIS path lands on an exact feature (FACE or CORNER), so the
    fallback must never fire and strict mode must be completely inert — that is
    what is pinned here.  A firing fallback would be a different fixture.
    """
    codes_off, sig_off, it_off = _drive(mat_mctc, tangent_type="Continuum",
                                        strict=0)
    codes_on, sig_on, it_on = _drive(mat_mctc, tangent_type="Continuum",
                                     strict=1)

    assert codes_off == codes_on, (
        f"strict_convergence changed the analyze codes on a path whose every "
        f"step lands on an EXACT cutoff feature: {codes_off} vs {codes_on}")
    assert len(sig_off) == len(sig_on) == NSTEPS
    assert np.array_equal(sig_off, sig_on), (
        "strict_convergence perturbed a converging confined-corner path; it "
        "must only ever refuse, never change a committed state")
    assert np.array_equal(it_off, it_on), (
        "strict_convergence changed the Newton history on a path it should not "
        "touch at all")
