"""LadrunoRebarBuckling (rebar-buckling WRAPPER over any UniaxialMaterial,
Dhakal-Maekawa, MAT_TAG 33001) -- Zone-A material battery.

The wrapper applies, in compression past a slenderness-dependent onset,
    sigma_buckled = r(e, lambda) * sigma_bare        (r in [resid, 1])
leaving the wrapped bar untouched. Tension and pre-onset compression pass
through; lsr<=0 is the exact identity gate.

Probes use the direct material harness (testUniaxialMaterial / setStrain /
getStress / getTangent). setStrain COMMITS each step, so:
  * path-dependent checks (B0/B1) drive a full stepwise path and compare
    step-by-step to the bare material;
  * formula / tangent checks (B2/B3/B5) probe each point as an INDEPENDENT
    one-step return from a fresh (zero-committed) material -- for a monotonic
    compression from virgin (e_cross = 0) that one-step map IS the backbone,
    which is exactly what the C++ committed-clone probe reconstructs.

Oracle (B2): the DM formula ported to Python below, fed the bare material's own
one-step stress sigma_bare(eps) and its onset backbone fStarL = sigma_bare(eStar)
-- i.e. the verbatim ReinforcingSteel::Buckled_stress_Dhakal modification.
See Ladruno_implementation/14_ladruno_rebar_buckling_adr.md.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# bar / wrapper parameters (steel, consistent N-mm units)
_FY = 420.0
_E = 200000.0
_B = 0.01          # Steel02 strain-hardening ratio


def _bar(tag, Fy=_FY, E=_E, b=_B):
    ops.uniaxialMaterial("Steel02", tag, Fy, E, b, 18.0, 0.925, 0.15)


def _wrap(tag, bartag, lsr, alpha=1.0, Fy=_FY, E=_E):
    ops.uniaxialMaterial("LadrunoRebarBuckling", tag, bartag,
                         "-lsr", lsr, "-alpha", alpha, "-fy", Fy, "-E", E)


# ==========================================================================
#  drivers
# ==========================================================================
def _drive(mat_fn, strains):
    """Stepwise (committing) path. Returns [(stress, tangent), ...]."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    mat_fn(1)
    ops.testUniaxialMaterial(1)
    out = []
    for e in strains:
        ops.setStrain(float(e))
        out.append((ops.getStress(), ops.getTangent()))
    return out


def _oneshot(mat_fn, eps, want="stress"):
    """Independent one-step return from a fresh material at strain eps."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    mat_fn(1)
    ops.testUniaxialMaterial(1)
    ops.setStrain(float(eps))
    return ops.getStress() if want == "stress" else ops.getTangent()


def _bar_fn(tag):
    _bar(tag)


def _wrap_fn(lsr, alpha=1.0):
    def f(tag):
        _bar(tag + 100)            # bar on a distinct tag
        _wrap(tag, tag + 100, lsr, alpha=alpha)
    return f


# ==========================================================================
#  Python port of ReinforcingSteel::Buckled_stress_Dhakal (monotonic, e_cross=0)
# ==========================================================================
def _eStar(lsr, Fy=_FY, E=_E):
    eY = Fy / E
    kfac = math.sqrt(eY * 2000.0)
    mag = max(7.0, 55.0 - 2.3 * kfac * lsr)
    return -mag * eY


def _dm_ref(eps, sBare, fStarL, lsr, alpha=1.0, Fy=_FY, E=_E):
    eY = Fy / E
    e = eps                                   # e_cross = 0 (monotonic from virgin)
    if e >= -eY:
        return sBare
    kfac = math.sqrt(eY * 2000.0)
    eStar = _eStar(lsr, Fy, E)
    fStar = fStarL * alpha * (1.1 - 0.016 * kfac * lsr)
    if fStar > -0.2 * Fy:
        fStar = -0.2 * Fy
    if e >= eStar:
        ratio = 1.0 - (1.0 - fStar / fStarL) * (e + eY) / (eStar + eY)
        return sBare * ratio
    num = fStar - 0.02 * E * (e - eStar)
    ave = sBare * num / fStarL
    if ave > -0.2 * Fy:
        ave = -0.2 * Fy
    return ave


# ==========================================================================
#  B0 -- lsr=0 pass-through is bit-for-bit the bare material (identity gate)
# ==========================================================================
@pytest.mark.t0m
def test_B0_identity_gate():
    eY = _FY / _E
    # full push-pull incl. deep compression and a reversal
    strains = [i * eY for i in
               [0, 1, 2, 4, 6, 3, 0, -1, -3, -6, -10, -15, -8, 0, 5]]
    bare = _drive(_bar_fn, strains)
    wrap = _drive(_wrap_fn(0.0), strains)
    for (sb, tb), (sw, tw) in zip(bare, wrap):
        assert sw == pytest.approx(sb, rel=1e-12, abs=1e-9)
        assert tw == pytest.approx(tb, rel=1e-12, abs=1e-9)


# ==========================================================================
#  B1 -- with lsr>0, tension and pre-onset (|e|<eY) compression pass through
# ==========================================================================
@pytest.mark.t0m
def test_B1_tension_and_preonset_passthrough():
    eY = _FY / _E
    # tension-only monotonic: always pass-through
    tens = [i * eY for i in [0, 1, 2, 4, 8, 12]]
    assert _drive(_bar_fn, tens) == pytest.approx(_drive(_wrap_fn(8.0), tens))

    # compression shallower than yield (from virgin): e = eps > -eY => pass-through
    comp = [i * eY for i in [0, -0.25, -0.5, -0.75, -0.99]]
    bare = _drive(_bar_fn, comp)
    wrap = _drive(_wrap_fn(8.0), comp)
    for (sb, _), (sw, _t) in zip(bare, wrap):
        assert sw == pytest.approx(sb, rel=1e-10, abs=1e-9)


# ==========================================================================
#  B2 -- buckled compression matches the ported DM formula, several lambda
# ==========================================================================
@pytest.mark.t1
@pytest.mark.parametrize("lsr", [6.0, 8.0, 12.0])
def test_B2_dm_formula(lsr):
    eY = _FY / _E
    fStarL = _oneshot(_bar_fn, _eStar(lsr))            # bare backbone at onset
    # sample the intermediate AND deep branches
    for k in [-1.5, -3, -6, -10, -15, -20, -25]:
        eps = k * eY
        sBare = _oneshot(_bar_fn, eps)
        sBuck = _oneshot(_wrap_fn(lsr), eps)
        ref = _dm_ref(eps, sBare, fStarL, lsr)
        assert sBuck == pytest.approx(ref, rel=1e-6, abs=1e-4), (
            f"lambda={lsr} eps={eps}: wrapper {sBuck} != DM ref {ref}")
        # a genuine knock-down (not pass-through) once past onset
        if eps < _eStar(lsr):
            assert abs(sBuck) < abs(sBare)


# ==========================================================================
#  B2b -- the multi-step committing path agrees with the formula too
#         (exercises commitState + the fStarL cache + the anchor)
# ==========================================================================
@pytest.mark.t1
def test_B2b_multistep_matches_formula():
    eY = _FY / _E
    lsr = 8.0
    fStarL = _oneshot(_bar_fn, _eStar(lsr))
    strains = [(-0.25 * i) * eY for i in range(0, 81)]   # 0 -> -20 eY, 80 steps
    wrap = _drive(_wrap_fn(lsr), strains)
    bare = _drive(_bar_fn, strains)
    for (sw, _), (sb, _b), eps in zip(wrap, bare, strains):
        ref = _dm_ref(eps, sb, fStarL, lsr)
        assert sw == pytest.approx(ref, rel=1e-6, abs=1e-3), (
            f"multistep eps={eps}: wrapper {sw} != DM ref {ref}")


# ==========================================================================
#  B3 -- larger lambda => earlier (shallower) onset + lower residual capacity
# ==========================================================================
@pytest.mark.t1
def test_B3_slenderness_trend():
    eY = _FY / _E
    # eStar (onset) is shallower for larger lambda
    assert abs(_eStar(12.0)) < abs(_eStar(8.0)) < abs(_eStar(6.0))

    # at a moderate compressive strain, more slender => more knocked down
    eps = -10.0 * eY
    s6 = _oneshot(_wrap_fn(6.0), eps)
    s8 = _oneshot(_wrap_fn(8.0), eps)
    s12 = _oneshot(_wrap_fn(12.0), eps)
    sbare = _oneshot(_bar_fn, eps)
    assert abs(s12) < abs(s8) < abs(s6) < abs(sbare)

    # the deep residual never falls less compressive than the -0.2 fy floor,
    # and at very large strain it clamps exactly onto it
    s12_mid = _oneshot(_wrap_fn(12.0), -28.0 * eY)
    assert s12_mid <= -0.2 * _FY + 1e-9              # at or below the floor
    s12_deep = _oneshot(_wrap_fn(12.0), -60.0 * eY)
    assert s12_deep == pytest.approx(-0.2 * _FY, rel=0.0, abs=1e-6)


# ==========================================================================
#  B5 -- analytic consistent tangent vs central finite difference
#        (each point an independent one-step return from a fresh material)
# ==========================================================================
@pytest.mark.t0m
def test_B5_consistent_tangent_fd():
    eY = _FY / _E
    lsr = 8.0
    d = 1.0e-7
    # -8 eY: intermediate branch;  -20 eY: deep (unfloored) branch
    for eps0 in (-8.0 * eY, -20.0 * eY):
        k = _oneshot(_wrap_fn(lsr), eps0, "tangent")
        sp = _oneshot(_wrap_fn(lsr), eps0 + d, "stress")
        sm = _oneshot(_wrap_fn(lsr), eps0 - d, "stress")
        kfd = (sp - sm) / (2.0 * d)
        assert k == pytest.approx(kfd, rel=1e-4, abs=1e-2), (
            f"consistent tangent {k} != FD {kfd} at eps={eps0}")


# ==========================================================================
#  Gomes–Appleton  (-model ga)
# ==========================================================================
_GA_C = 9.42477796076938       # 3*pi (upstream constant)


def _wrap_ga(lsr, reduction=0.0, ff=0.5):
    """Wrapper with -model ga. GA needs E (e_cross) but not fy."""
    def f(tag):
        _bar(tag + 100)
        ops.uniaxialMaterial("LadrunoRebarBuckling", tag, tag + 100,
                             "-model", "ga", "-lsr", lsr,
                             "-reduction", reduction, "-fsufrac", ff, "-E", _E)
    return f


def _ga_ref(eps, sBare, lsr, reduction=0.0, ff=0.5):
    """Python port of ReinforcingSteel::Buckled_stress_Gomes for the
    monotonic-from-virgin case (e_cross=0, fsup=0)."""
    e_cross, fsup = 0.0, 0.0
    if eps >= e_cross:
        return sBare
    fs_buck = math.sqrt(32.0 / (e_cross - eps)) / (_GA_C * lsr)
    beta_loc, gama_loc, Dft = 1.0, 0.1, 0.25
    sd = abs(fs_buck - 1.0)
    if sd <= Dft:
        beta_loc = 1.0 - gama_loc * (Dft - sd) / Dft
    m = min(1.0, fs_buck)
    factor = m * beta_loc * (1.0 - reduction) + reduction
    return fsup * ff - (factor + ff) * (fsup * ff - sBare) / (1.0 + ff)


# GA0 -- -lsr 0 identity gate (full push-pull == bare)
@pytest.mark.t0m
def test_GA0_identity_gate():
    eY = _FY / _E
    strains = [i * eY for i in [0, 2, 5, 2, 0, -3, -8, -15, -6, 0, 4]]
    bare = _drive(_bar_fn, strains)
    wrap = _drive(_wrap_ga(0.0), strains)
    for (sb, tb), (sw, tw) in zip(bare, wrap):
        assert sw == pytest.approx(sb, rel=1e-12, abs=1e-9)
        assert tw == pytest.approx(tb, rel=1e-12, abs=1e-9)


# GA1 -- tension pass-through (GA gate: eps >= e_cross == 0 from virgin)
@pytest.mark.t0m
def test_GA1_tension_passthrough():
    eY = _FY / _E
    tens = [i * eY for i in [0, 1, 3, 6, 10]]
    assert _drive(_bar_fn, tens) == pytest.approx(_drive(_wrap_ga(8.0), tens))


# GA2 -- buckled compression matches the ported Gomes formula
@pytest.mark.t1
@pytest.mark.parametrize("lsr,red", [(6.0, 0.0), (8.0, 0.0), (12.0, 0.0), (8.0, 0.3)])
def test_GA2_gomes_formula(lsr, red):
    eY = _FY / _E
    for k in [-1, -3, -6, -10, -15, -22]:
        eps = k * eY
        sBare = _oneshot(_bar_fn, eps)
        sBuck = _oneshot(_wrap_ga(lsr, reduction=red), eps)
        ref = _ga_ref(eps, sBare, lsr, reduction=red)
        assert sBuck == pytest.approx(ref, rel=1e-6, abs=1e-4), (
            f"λ={lsr} r={red} eps={eps}: wrapper {sBuck} != GA ref {ref}")


# GA3 -- reduction limits: r=1 => no buckling; r=0 => genuine knock-down;
#        more slender => more knock-down
@pytest.mark.t1
def test_GA3_reduction_and_slenderness():
    eY = _FY / _E
    eps = -10.0 * eY
    sbare = _oneshot(_bar_fn, eps)
    # r = 1 is the GA "no buckling" value -> pass-through
    assert _oneshot(_wrap_ga(8.0, reduction=1.0), eps) == pytest.approx(sbare, rel=1e-6)
    # r = 0 (full) knocks the stress down, more so for a more slender bar
    s6 = _oneshot(_wrap_ga(6.0), eps)
    s8 = _oneshot(_wrap_ga(8.0), eps)
    s12 = _oneshot(_wrap_ga(12.0), eps)
    assert abs(s12) < abs(s8) < abs(s6) < abs(sbare)


# GA5 -- consistent-tangent FD (smooth region, away from the fs_buck≈1 kink)
@pytest.mark.t0m
def test_GA5_consistent_tangent_fd():
    eY = _FY / _E
    lsr = 8.0
    d = 1.0e-7
    for eps0 in (-10.0 * eY, -25.0 * eY):
        k = _oneshot(_wrap_ga(lsr), eps0, "tangent")
        sp = _oneshot(_wrap_ga(lsr), eps0 + d, "stress")
        sm = _oneshot(_wrap_ga(lsr), eps0 - d, "stress")
        kfd = (sp - sm) / (2.0 * d)
        assert k == pytest.approx(kfd, rel=1e-3, abs=1e-2), (
            f"GA consistent tangent {k} != FD {kfd} at eps={eps0}")
