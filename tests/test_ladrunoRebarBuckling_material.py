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
from _testbed.roundtrip import database_roundtrip

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


# ==========================================================================
#  B6 -- composition  Fatigue ∘ RebarBuckling ∘ bar  builds, runs, and the
#        buckling reduction propagates through the outer Fatigue layer; the
#        Fatigue rupture still triggers on the buckled response.
# ==========================================================================
def _fatigue_buck(lsr):
    """Fatigue(tag) wraps RebarBuckling(tag+50) wraps Steel02(tag+100). The
    Fatigue -min/-max bounds give a DETERMINISTIC rupture for the test."""
    def f(tag):
        _bar(tag + 100)
        ops.uniaxialMaterial("LadrunoRebarBuckling", tag + 50, tag + 100,
                             "-lsr", lsr, "-fy", _FY, "-E", _E)
        ops.uniaxialMaterial("Fatigue", tag, tag + 50, "-min", -0.01, "-max", 0.01)
    return f


@pytest.mark.t1
def test_B6_composition_fatigue_buckling():
    # ... -0.008 is buckled but within the Fatigue -min; -0.02 trips Fatigue.
    strains = [0.0, -0.002, -0.004, -0.006, -0.008, -0.012, -0.02]
    runA = _drive(_fatigue_buck(8.0), strains)    # with buckling
    runB = _drive(_fatigue_buck(0.0), strains)    # identity gate (no buckling)

    # buckling shows through the Fatigue wrapper at the deep-but-unfailed point
    sA8, sB8 = runA[4][0], runB[4][0]             # eps = -0.008
    assert abs(sA8) < abs(sB8), (
        f"buckling did not propagate through Fatigue: {sA8} vs {sB8}")

    # Fatigue rupture still triggers on the buckled response (stress collapses)
    assert abs(runA[-1][0]) < 1.0, f"Fatigue did not rupture the buckled mat: {runA[-1][0]}"
    assert abs(runB[-1][0]) < 1.0


# ==========================================================================
#  B7 -- sendSelf/recvSelf + broker round-trip (nested material), via the
#        FE_Datastore. Covers the serialization defect class (D4) incl. GA's
#        extra (gaReduction, gaFsuFrac) fields. Skips if the build lacks a DB.
# ==========================================================================
def _b7_build(mat_fn):
    def f():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        mat_fn(1)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.element("Truss", 1, 1, 2, 1.0, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0, 0.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-12, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 2, 1, -0.01)   # one step into buckling
        ops.analysis("Static")
        ops.analyze(1)
    return f


@pytest.mark.t1
def test_B7_database_roundtrip_dm():
    database_roundtrip(_b7_build(_wrap_fn(8.0)), probe_nodes=[2], ndf=2)


@pytest.mark.t1
def test_B7_database_roundtrip_ga():
    database_roundtrip(_b7_build(_wrap_ga(8.0)), probe_nodes=[2], ndf=2)


# ==========================================================================
#  v2 CYCLIC re-straightening (B4 battery) -- residual-gap memory.
#
#  Model (h-convention): sigma = sigma_bare + g, with the RAISE g >= 0 (buckling
#  reduces the compressive magnitude => raises sigma toward 0). v1 monotonic is
#  exactly g == g_law = sigma_v1 - sigma_bare. On a tension reload the raise is
#  HELD on the compression side then decayed to 0 over L_rs (smoothstep) on the
#  tension side, rejoining bare. The RESTRAIGHTEN->BUCKLING re-entry is
#  C0-continuous (residual-gap carry).  See ADR 14 §9.
# ==========================================================================
def _wrap_rs(lsr, mode="lambda", c=1.0, alpha=1.0, Fy=_FY, E=_E):
    """DM wrapper with an explicit -restraighten span mode."""
    def f(tag):
        _bar(tag + 100)
        args = ["LadrunoRebarBuckling", tag, tag + 100,
                "-lsr", lsr, "-alpha", alpha, "-fy", Fy, "-E", E]
        if mode == "lambda":
            args += ["-restraighten", "lambda"]
        else:
            args += ["-restraighten", "c", c]
        ops.uniaxialMaterial(*args)
    return f


def _wrap_ga_rs(lsr, mode="lambda", c=1.0, reduction=0.0, ff=0.5, E=_E):
    def f(tag):
        _bar(tag + 100)
        args = ["LadrunoRebarBuckling", tag, tag + 100, "-model", "ga",
                "-lsr", lsr, "-reduction", reduction, "-fsufrac", ff, "-E", E]
        if mode == "lambda":
            args += ["-restraighten", "lambda"]
        else:
            args += ["-restraighten", "c", c]
        ops.uniaxialMaterial(*args)
    return f


def _Lrs_lambda(lsr, Fy=_FY, E=_E):
    eY = Fy / E
    f = 0.5 + 0.10 * max(0.0, lsr - 5.0)
    f = min(6.0, max(0.5, f))
    return f * eY


# B4a -- reduce-to-v1: pure monotonic compression reproduces the v1 DM/GA law
#        (the new state machine, with no reversal, is the v1 law verbatim).
@pytest.mark.t1
def test_B4a_reduce_to_v1():
    eY = _FY / _E
    lsr = 8.0
    fStarL = _oneshot(_bar_fn, _eStar(lsr))
    strains = [(-0.25 * i) * eY for i in range(0, 121)]   # 0 -> -30 eY
    wrap = _drive(_wrap_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    for (sw, _), (sb, _b), eps in zip(wrap, bare, strains):
        ref = _dm_ref(eps, sb, fStarL, lsr)
        assert sw == pytest.approx(ref, rel=1e-6, abs=1e-3), (
            f"B4a monotonic != v1 DM at eps={eps}: {sw} vs {ref}")
    # GA leg
    wrapg = _drive(_wrap_ga_rs(lsr, "lambda"), strains)
    bareg = _drive(_bar_fn, strains)
    for (sw, _), (sb, _b), eps in zip(wrapg, bareg, strains):
        ref = _ga_ref(eps, sb, lsr)
        assert sw == pytest.approx(ref, rel=1e-6, abs=1e-3), (
            f"B4a monotonic != v1 GA at eps={eps}: {sw} vs {ref}")


# B4b -- compression -> full tension reload rejoins the bare curve (PASS).
@pytest.mark.t1
def test_B4b_full_reload_rejoins_bare():
    eY = _FY / _E
    lsr = 8.0
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # 0 -> -25 eY (buckled)
    up = [(-25.0 + 0.5 * i) * eY for i in range(1, 57)]     # -25 eY -> +3 eY
    strains = down + up
    wrap = _drive(_wrap_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    # buckled at the reversal (raise above bare)
    sw_rev, sb_rev = wrap[len(down) - 1][0], bare[len(down) - 1][0]
    assert sw_rev > sb_rev + 1.0, f"not buckled at reversal: {sw_rev} vs {sb_rev}"
    # past e_cross_rs(=0) + L_rs the wrapper has rejoined bare
    assert wrap[-1][0] == pytest.approx(bare[-1][0], rel=1e-6, abs=1e-2)
    assert wrap[-1][1] == pytest.approx(bare[-1][1], rel=1e-4, abs=1e-1)


# B4c -- Phase 1 (compression-side reload) holds the raise parallel to bare:
#        the raise g = sigma_wrap - sigma_bare stays ~constant; tangent ~ bare.
@pytest.mark.t1
def test_B4c_phase1_holds_raise():
    eY = _FY / _E
    lsr = 8.0
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # -25 eY
    up = [(-25.0 + 0.5 * i) * eY for i in range(1, 46)]     # -25 -> -2.5 eY (Phase 1, e_cross_rs=0)
    strains = down + up
    wrap = _drive(_wrap_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    n = len(down)
    g_rev = wrap[n - 1][0] - bare[n - 1][0]
    assert g_rev > 1.0
    for i in range(n, len(strains)):
        g = wrap[i][0] - bare[i][0]
        assert g == pytest.approx(g_rev, rel=0.05, abs=1.0), (
            f"Phase-1 raise drifted at eps={strains[i]}: {g} vs {g_rev}")


# B4d -- Phase 2 bracket: strictly between bare and the held-raise line.
@pytest.mark.t1
def test_B4d_phase2_bracket():
    eY = _FY / _E
    lsr = 8.0
    L = _Lrs_lambda(lsr)
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # -25 eY
    up = []
    e = -25.0
    while e < 1.5:                                          # fine reload through Phase 2
        e += 0.05
        up.append(e * eY)
    strains = down + up
    wrap = _drive(_wrap_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    g_rev = wrap[len(down) - 1][0] - bare[len(down) - 1][0]
    found = 0
    for (sw, _), (sb, _b), eps in zip(wrap, bare, strains):
        if 1e-9 < eps < L - 1e-9:                          # strictly inside Phase 2 (e_cross_rs=0)
            assert sw > sb + 1e-6, f"not above bare at {eps}: {sw} vs {sb}"
            assert sw < sb + g_rev + 1e-3, f"not below held-raise at {eps}"
            found += 1
    assert found > 0, "no Phase-2 samples found"


# B4f -- partial reload then re-compress: the raise g is C0-continuous across
#        the RESTRAIGHTEN->BUCKLING re-entry (no stress jump), then re-buckles.
@pytest.mark.t1
def test_B4f_recompress_c0_continuous():
    eY = _FY / _E
    lsr = 8.0
    L = _Lrs_lambda(lsr)
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # -25 eY
    # reload partway INTO Phase 2 (to ~0.5*L past the crossing), then re-compress
    mid = 0.5 * L / eY
    up = [(-25.0 + 0.5 * i) * eY for i in range(1, int((25.0 + mid) / 0.5) + 1)]
    redown = [up[-1] + (-0.1 * i) * eY for i in range(1, 121)]   # re-compress to deep
    strains = down + up + redown
    wrap = _drive(_wrap_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    g = [w[0] - b[0] for w, b in zip(wrap, bare)]
    k_re = len(down) + len(up)                              # first re-compress step index
    g_rev = g[len(down) - 1]
    # C0 of the raise across the re-entry: no jump
    for i in range(k_re, k_re + 3):
        assert abs(g[i] - g[i - 1]) < 0.10 * g_rev + 1.0, (
            f"raise jumped at re-entry step {i}: {g[i-1]} -> {g[i]}")
    # genuinely re-buckles deeper (raise grows back above bare)
    assert g[-1] > 1.0, f"did not re-buckle on re-compression: g={g[-1]}"
    # nothing pathological
    for sw, _ in wrap:
        assert abs(sw) < 1e4 and sw == sw, "NaN / blow-up in re-compress cycle"


# B4-GA -- the law-agnostic RESTRAIGHTEN: GA buckle -> reload rejoins bare too.
@pytest.mark.t1
def test_B4_ga_cyclic_rejoins_bare():
    eY = _FY / _E
    lsr = 8.0
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # -25 eY (GA buckles past e_cross=0)
    up = [(-25.0 + 0.5 * i) * eY for i in range(1, 57)]     # -> +3 eY
    strains = down + up
    wrap = _drive(_wrap_ga_rs(lsr, "lambda"), strains)
    bare = _drive(_bar_fn, strains)
    assert wrap[len(down) - 1][0] > bare[len(down) - 1][0] + 1.0   # buckled raise
    assert wrap[-1][0] == pytest.approx(bare[-1][0], rel=1e-6, abs=1e-2)  # rejoined


# B4i -- L_rs floor: a tiny c clamps the span to eY (NOT c*|e_rev-e_cross_rs|),
#        so the recovery is finite/bounded and rejoins ~one yield strain past the
#        crossing, not instantly. Guards against D0/L_rs blow-up.
@pytest.mark.t1
def test_B4i_lrs_floor():
    eY = _FY / _E
    lsr = 8.0
    # buckle, then reload with fine steps; c=0.01 -> raw span 0.25eY would be tiny,
    # but the eY floor makes L_rs = eY (rejoin ~+1 eY past e_cross_rs = 0).
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # 0 -> -25 eY
    up = [(-25.0 + 0.25 * i) * eY for i in range(1, 117)]   # -> +4 eY (0.25 eY steps)
    strains = down + up
    wrap = _drive(_wrap_rs(lsr, mode="c", c=0.01), strains)
    bare = _drive(_bar_fn, strains)
    for (sw, kt), eps in zip(wrap, strains):
        assert sw == sw and abs(sw) < 1e4, f"blow-up at {eps}: {sw}"
        assert abs(kt) < 1e3 * _E, f"tangent blow-up at {eps}: {kt}"

    def at(target):
        j = min(range(len(strains)), key=lambda k: abs(strains[k] - target * eY))
        return wrap[j][0], bare[j][0]

    # floor binds at eY: still raised at +0.5 eY (q<1), rejoined by +2 eY (q>1).
    sw_half, sb_half = at(0.5)
    assert sw_half > sb_half + 1e-3, "floor too small: already rejoined at +0.5 eY"
    sw_two, sb_two = at(2.0)
    assert sw_two == pytest.approx(sb_two, rel=1e-6, abs=1e-1), "did not rejoin by +2 eY"


# B4h -- v2 serialization round-trip from a committed RESTRAIGHTEN state.
def _b7_build_cyclic(mat_fn, n_into_restraighten=8):
    """Compress into buckling then partial-reload so the committed branch is
    RESTRAIGHTEN, then hand off to the datastore round-trip."""
    def f():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        mat_fn(1)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.element("Truss", 1, 1, 2, 1.0, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0, 0.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-10, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 2, 1, -0.001)   # compress (buckle)
        ops.analysis("Static")
        ops.analyze(30)
        ops.integrator("DisplacementControl", 2, 1, +0.001)   # partial reload -> RESTRAIGHTEN
        ops.analyze(n_into_restraighten)
    return f


@pytest.mark.t1
def test_B4h_serialization_roundtrip_restraighten():
    database_roundtrip(_b7_build_cyclic(_wrap_rs(8.0, "lambda")),
                       probe_nodes=[2], ndf=2)


# B4e -- consistent-tangent FD on the RESTRAIGHTEN Phase-2 smoothstep branch.
#        The direct harness commits every setStrain, so the FD is done by
#        PATH-REPLAY: drive the identical committed path, then take the final
#        single step to e0, e0+d, e0-d from the SAME committed state (the
#        material is path-independent, so the committed pre-state is identical).
@pytest.mark.t1
def test_B4e_restraighten_tangent_fd():
    eY = _FY / _E
    lsr = 8.0
    d = 1.0e-7
    L = _Lrs_lambda(lsr)
    e0 = 0.5 * L                                            # mid Phase-2 (q=0.5), e_cross_rs=0
    down = [(-0.5 * i) * eY for i in range(0, 51)]          # 0 -> -25 eY (buckled)
    reload_pts = []
    e = -25.0 * eY
    while e < e0 - 1e-12:
        e = min(e + 0.25 * eY, e0)
        reload_pts.append(e)
    base = down + reload_pts                                # ends exactly at e0
    pre = base[:-1]
    s0 = _drive(_wrap_rs(lsr, "lambda"), base)[-1]
    sp = _drive(_wrap_rs(lsr, "lambda"), pre + [e0 + d])[-1][0]
    sm = _drive(_wrap_rs(lsr, "lambda"), pre + [e0 - d])[-1][0]
    k = s0[1]
    kfd = (sp - sm) / (2.0 * d)
    assert k == pytest.approx(kfd, rel=1e-4, abs=1e-2), (
        f"RESTRAIGHTEN Phase-2 tangent {k} != FD {kfd} at e0={e0}")
    # sanity: this point is genuinely inside Phase 2 (raised above bare, q in (0,1))
    sb0 = _drive(_bar_fn, base)[-1][0]
    assert s0[0] > sb0 + 1e-6


# B4g -- structural Newton across the reversal + the q=0/q=1 seams: a real Truss
#        under DisplacementControl must converge through buckle -> reload.
def _b4g_run(mat_fn, comp_steps, comp_dU, reload_steps, reload_dU):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    mat_fn(1)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    fails = 0
    ops.integrator("DisplacementControl", 2, 1, comp_dU)
    for _ in range(comp_steps):
        if ops.analyze(1) != 0:
            fails += 1
    ops.integrator("DisplacementControl", 2, 1, reload_dU)
    for _ in range(reload_steps):
        if ops.analyze(1) != 0:
            fails += 1
    return fails


@pytest.mark.t1
def test_B4g_structural_newton_seams():
    eY = _FY / _E
    # compress past onset into buckling, then reload through the crossing (0) and
    # past e_cross_rs + L_rs (~0.8 eY) so both Phase-2 seams are crossed.
    f_dm = _b4g_run(_wrap_rs(8.0, "lambda"),
                    comp_steps=60, comp_dU=-0.001,      # -> ~ -28.6 eY
                    reload_steps=80, reload_dU=+0.0005)  # -> ~ +0.4 eY ... +
    assert f_dm == 0, f"DM: Newton failed to converge on {f_dm} step(s)"
    f_ga = _b4g_run(_wrap_ga_rs(8.0, "lambda"),
                    comp_steps=60, comp_dU=-0.001,
                    reload_steps=80, reload_dU=+0.0005)
    assert f_ga == 0, f"GA: Newton failed to converge on {f_ga} step(s)"
