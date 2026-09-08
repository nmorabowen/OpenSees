"""ADR-97 P5 (`wp/97f-explicit-gate`, M9) -- StiffSoilShear/StiffSoilCap step-1
finite-ness and hardening-path admissibility, plus the ADR-97 D3 refusal for
`Closest_Point` on the whole StiffSoil family.

Root cause (`Ladruno_implementation/reviews/adr97_p5_report.md`,
`Ladruno_implementation/_adr94_components.md` finding 3): TWO independent
bugs in `StiffSoilShear_YF`/`_PF`, both NaN-producing:

* Bug A -- `StiffSoilShear_YF.h`'s ``qf`` used an unguarded ``cot(phi)``,
  which is ``Inf`` at ``phi == 0`` (a normal, legitimate cohesive-only /
  undrained-clay choice). ``Inf * 2*sin(phi)`` (``== Inf * 0`` at ``phi==0``)
  is the IEEE-754 indeterminate NaN. Fixed by an algebraic rewrite that
  removes ``cot`` entirely: ``qf = 2*(c*cos(phi) + sigma3*sin(phi)) /
  (1 - sin(phi))`` -- identical to the old formula for any ``phi != 0``
  (matched to ~1e-14 relative) and gives the correct Tresca limit
  ``qf -> 2c`` as ``phi -> 0`` instead of NaN.
* Bug B -- `StiffSoilShear_PF.h`'s numerically-differentiated flow direction
  normalizes by its own norm with no zero-guard. At an EXACTLY hydrostatic
  trial stress (the natural first step of any isotropic consolidation leg,
  and `InitialP0`'s own hydrostatic seed) `computeMobilizedDilatancy()`
  returns `psi_m == 0` exactly, which makes the central-difference numerator
  the exact zero vector on all six Voigt axes -- `norm == 0.0` and
  `vv_out /= norm` is the indeterminate 0/0. This IS the ADR-94 R3a "6/194
  non-finite cloud points" finding: those 6 points are exactly the
  hydrostatic-axis points, not "random general (non-diagonal) points" as
  originally logged. Fixed by a zero-guard (``setZero()`` fallback, the same
  precedent the DruckerPrager NaN*0 fix used).

`StiffSoil_EL` supplies no `ELASTICITY_STRESS_DERIVATIVE` (ADR-97 D6) -- NOT
required here, because `Closest_Point` refuses the whole StiffSoil family at
parse time regardless (ADR-97 D3); `Backward_Euler` never evaluates that
block. StiffSoilCap has no NaN defect of its own (neither bug reaches it);
it gets a smoke test only.

Zone-A.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]

MC_C = 10.0
MC_PSI = 5.0
SS_E50_REF = 20000.0
SS_EUR_REF = 60000.0
SS_RF = 0.9
SS_M = 0.5
SS_PREF = 100.0
NU = 0.2

RTOL = 1.0e-6  # testbed policy floor


# ===========================================================================
# deck builders
# ===========================================================================
def mat_stiffsoil_shear(tag, phi=30.0, method="Backward_Euler",
                         tangent="Secant", initial_p0=0.0, eps0=0.0,
                         niter=100):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "StiffSoilShear_YF", "StiffSoilShear_PF", "StiffSoil_EL",
        "EpsQpShear(StiffSoilShearHardening):",
        "Begin_Model_Parameters",
        "PoissonsRatio", NU,
        "MC_phi", phi, "MC_c", MC_C, "MC_psi", MC_PSI, "MC_ds", 1e-4,
        "SS_E50_ref", SS_E50_REF, "SS_Eur_ref", SS_EUR_REF,
        "SS_Rf", SS_RF, "SS_m", SS_M, "SS_pref", SS_PREF,
        "MassDensity", 0.0, "InitialP0", float(initial_p0),
        "End_Model_Parameters",
        "Begin_Internal_Variables", "EpsQpShear", float(eps0),
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", niter,
        "End_Integration_Options",
    )


def mat_stiffsoil_cap(tag, method="Backward_Euler", tangent="Secant",
                       initial_p0=0.0, phi=30.0, pc0=100.0, niter=100):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "StiffSoilCap_YF", "StiffSoilCap_PF", "StiffSoil_EL",
        "CapPressure(StiffSoilCapHardening):",
        "Begin_Model_Parameters",
        "PoissonsRatio", NU,
        "MC_phi", phi, "MC_c", MC_C, "MC_ds", 1e-4,
        "SS_Eur_ref", SS_EUR_REF, "SS_pref", SS_PREF, "SS_m", SS_M,
        "SS_alpha", 0.5, "SS_beta", 1.0,
        "MassDensity", 0.0, "InitialP0", float(initial_p0),
        "End_Model_Parameters",
        "Begin_Internal_Variables", "CapPressure", float(pc0),
        "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", niter,
        "End_Integration_Options",
    )


def _u_of(eps, x, y, z):
    e11, e22, e33, g12, g23, g13 = eps
    return (e11 * x + g12 * y + g13 * z, e22 * y + g23 * z, e33 * z)


def drive(mat_fn, legs, nstep=10, ele="LadrunoBrick", tol=1e-10, maxiter=60):
    """One LadrunoBrick unit cube, system UmfPack, prescribed homogeneous
    strain via sp + Path series (same idiom as tests/test_adr97_p1_smooth.py
    ::drive)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    mat_fn(1)
    ops.element(ele, 1, *range(1, 9), 1)

    U = []
    for eps in legs:
        d = {}
        for k, (x, y, z) in enumerate(NODES):
            u = _u_of(eps, float(x), float(y), float(z))
            for dd in (1, 2, 3):
                d[(k + 1, dd)] = u[dd - 1]
        U.append(d)

    nl = len(legs)
    owner = {}
    for key in U[0]:
        vals = [0.0] + [U[L][key] for L in range(nl)]
        changed = [L for L in range(nl) if abs(vals[L + 1] - vals[L]) > 1e-18]
        assert len(changed) <= 1, "path not component-disjoint at %s" % (key,)
        owner[key] = changed[0] if changed else 0

    for L in range(nl):
        if nl == 1:
            ops.timeSeries("Linear", 100 + L)
        else:
            times = [float(t) for t in range(nl + 2)]
            values = [0.0] * (L + 1) + [1.0] * (nl - L + 1)
            ops.timeSeries("Path", 100 + L, "-time", *times, "-values", *values)
        ops.pattern("Plain", 100 + L, 100 + L)
        for key, o in owner.items():
            if o == L:
                ops.sp(key[0], key[1], float(U[L][key]))

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxiter, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.analysis("Static")

    out = {"codes": [], "sigma": []}
    for _ in range(nl * nstep):
        rc = ops.analyze(1)
        out["codes"].append(rc)
        if rc != 0:
            break
        out["sigma"].append(np.array(list(ops.eleResponse(1, "stresses"))[0:6]))
    out["sigma"] = np.array(out["sigma"]) if out["sigma"] else np.zeros((0, 6))
    return out


def _constructible(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    try:
        mat_fn(1)
        ops.element("LadrunoBrick", 1, *range(1, 9), 1)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


@pytest.fixture(scope="module")
def stiffsoil_available():
    if not _constructible(lambda t: mat_stiffsoil_shear(t)):
        pytest.skip("ASDPlasticMaterial3D / StiffSoilShear / LadrunoBrick "
                    "not available in this build")


# ===========================================================================
# qf oracle -- SAME algebraic rewrite as the C++ fix (no cot(phi)).
# ===========================================================================
def _qa(sigma3, phi_deg=30.0, c=MC_C, rf=SS_RF):
    phi = math.radians(phi_deg)
    sin_phi, cos_phi = math.sin(phi), math.cos(phi)
    qf = 2.0 * (c * cos_phi + sigma3 * sin_phi) / (1.0 - sin_phi)
    return qf / rf


def _principal_q(sigma):
    """q = sigma1 - sigma3 (geotechnical/compression-positive convention,
    i.e. on the NEGATED tension-positive Voigt stress OpenSees stores)."""
    sig_geo = -np.asarray(sigma, dtype=float)
    mat = np.array([[sig_geo[0], sig_geo[3], sig_geo[5]],
                     [sig_geo[3], sig_geo[1], sig_geo[4]],
                     [sig_geo[5], sig_geo[4], sig_geo[2]]])
    w = np.sort(np.linalg.eigvalsh(mat))
    sigma3, sigma2, sigma1 = w[0], w[1], w[2]
    return sigma1 - sigma3, sigma3


# ===========================================================================
# Bug A regression: MC_phi == 0.0 must not NaN at step 1.
# ===========================================================================
def test_bug_a_phi_zero_gives_finite_step1(stiffsoil_available):
    out = drive(lambda t: mat_stiffsoil_shear(t, phi=0.0, initial_p0=-50.0),
                legs=[(0.0, 0.0, -0.001, 0.0, 0.0, 0.0)], nstep=5)
    assert out["codes"] and out["codes"][0] == 0, out["codes"]
    assert np.all(np.isfinite(out["sigma"])), out["sigma"]


# ===========================================================================
# Bug B regression: a purely isotropic (hydrostatic) leg must not NaN.
#
# MEASURED (this build): the fix removes every "NaN!" print and every
# committed value is finite, but the deeper isotropic path still cuts off
# with a clean refusal (rc=-3) a couple of steps in -- a SEPARATE, legitimate
# Newton/step-size limit unrelated to Bug B (see p5/run_iso.log from the
# pre-investigation probe, which shows the identical [0, 0, -3] shape once
# the zero-guard is applied). A clean refusal is not the defect; a NaN
# COMMIT is (ADR-94 B2's strong invariant: no committed state is ever
# non-finite, whatever the step codes are).
# ===========================================================================
def test_bug_b_isotropic_leg_gives_finite_step1(stiffsoil_available):
    out = drive(lambda t: mat_stiffsoil_shear(t, phi=30.0, initial_p0=0.0),
                legs=[(-0.001, -0.001, -0.001, 0.0, 0.0, 0.0)], nstep=5)
    assert np.all(np.isfinite(out["sigma"])), out["sigma"]
    assert out["codes"] and out["codes"][0] == 0, (
        "the very first step -- the exact hydrostatic trial that used to "
        "hit the Bug B 0/0 -- must still succeed: %r" % out["codes"])
    bad = [c for c in out["codes"] if c not in (0, -3)]
    assert not bad, (
        "unexpected analyze() codes on the isotropic leg: %r (0 = step "
        "taken, -3 = the global Newton gave up after a clean material "
        "refusal)" % out["codes"])


# ===========================================================================
# Non-trivial hardening path: isotropic consolidation, then triaxial
# compression path, run in enough increments to accumulate a non-trivial
# EpsQpShear well past first yield.  Invariant checked: q never exceeds the
# hyperbolic law's own asymptote q_a = qf/Rf, at every committed step -- true
# regardless of the accumulated EpsQpShear internal variable (the hardening
# only slides q up TOWARD q_a as eps_qp_shear -> infinity; it can never push
# q past it). This is the strongest hardening-state-independent invariant of
# the corrected qf formula.
#
# (A single Path leg, not two: ``drive()``'s multi-leg "owner" bookkeeping
# requires each dof's target to change in at most one leg-to-leg step, i.e.
# legs each ENGAGE a previously-idle dof rather than re-target one that is
# already moving -- a continuously-deepening e33 across two legs is not
# component-disjoint. A single proportional ramp from zero to a markedly
# non-isotropic target already mixes isotropic and deviatoric loading and
# drives the deck well past first yield.)
# ===========================================================================
def test_hardening_path_never_exceeds_the_hyperbolic_asymptote(
        stiffsoil_available):
    leg = (-0.002, -0.002, -0.02, 0.0, 0.0, 0.0)   # triaxial compression
    out = drive(lambda t: mat_stiffsoil_shear(t, phi=30.0, initial_p0=-50.0),
                legs=[leg], nstep=16)
    assert np.all(np.isfinite(out["sigma"]))
    assert len(out["sigma"]) == 16, out["codes"]

    for sig in out["sigma"]:
        q, sigma3 = _principal_q(sig)
        if sigma3 <= 1e-10:
            continue   # outside the hyperbolic law's own admissible range
        qa = _qa(sigma3)
        tol = RTOL * max(qa, 1.0)
        assert q <= qa + tol, (
            f"committed q={q:.6f} exceeded the hyperbolic asymptote "
            f"q_a={qa:.6f} at sigma3={sigma3:.6f} -- the corrected qf "
            f"formula (Bug A fix) or the hyperbolic law itself has "
            f"regressed")


# ===========================================================================
# StiffSoilCap: smoke test only (ADR-94 H11 -- zero prior Zone-A coverage).
# Neither Bug A nor Bug B is reachable from StiffSoilCap's own YF/PF, so this
# is a no-regression pin, not a bug reproducer.
# ===========================================================================
@pytest.mark.parametrize("phi", [0.0, 30.0])
@pytest.mark.parametrize("pc0", [0.0, 100.0])
def test_stiffsoilcap_smoke_finite(phi, pc0):
    if not _constructible(lambda t: mat_stiffsoil_cap(t, phi=phi, pc0=pc0)):
        pytest.skip("ASDPlasticMaterial3D / StiffSoilCap / LadrunoBrick not "
                    "available in this build")
    out = drive(lambda t: mat_stiffsoil_cap(t, initial_p0=-50.0, phi=phi,
                                            pc0=pc0),
                legs=[(-0.001, -0.001, -0.001, 0.0, 0.0, 0.0)], nstep=5)
    assert np.all(np.isfinite(out["sigma"])), (phi, pc0, out["sigma"])


# ===========================================================================
# ADR-97 D3: Closest_Point refuses the WHOLE StiffSoil family, checked BOTH
# WAYS -- the same deck (same YF/PF/EL/parameters) must build under
# Backward_Euler and be refused under Closest_Point.
# ===========================================================================
def test_closest_point_refuses_stiffsoilshear_both_ways(stiffsoil_available):
    assert _constructible(
        lambda t: mat_stiffsoil_shear(t, method="Backward_Euler")), (
        "control StiffSoilShear deck does not build under Backward_Euler -- "
        "fix the deck before reading anything into the refusal below")
    assert not _constructible(
        lambda t: mat_stiffsoil_shear(t, method="Closest_Point",
                                      tangent="Secant")), (
        "integration_method Closest_Point was accepted for StiffSoilShear -- "
        "ADR-97 D3 says the whole StiffSoil family stays refused (P5 adds "
        "only its Backward_Euler NaN fixes, not a closest-point return)")


def test_closest_point_refuses_stiffsoilcap_both_ways():
    if not _constructible(lambda t: mat_stiffsoil_cap(t)):
        pytest.skip("ASDPlasticMaterial3D / StiffSoilCap / LadrunoBrick not "
                    "available in this build")
    assert _constructible(
        lambda t: mat_stiffsoil_cap(t, method="Backward_Euler")), (
        "control StiffSoilCap deck does not build under Backward_Euler -- "
        "fix the deck before reading anything into the refusal below")
    assert not _constructible(
        lambda t: mat_stiffsoil_cap(t, method="Closest_Point",
                                    tangent="Secant")), (
        "integration_method Closest_Point was accepted for StiffSoilCap -- "
        "ADR-97 D3 says the whole StiffSoil family stays refused")
