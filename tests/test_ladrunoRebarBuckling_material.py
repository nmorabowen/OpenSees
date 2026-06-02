"""LadrunoRebarBuckling (Dhakal-Maekawa rebar-buckling WRAPPER, MAT_TAG 33001)
— Zone-A material battery.

The wrapper applies a compression-side buckled-average reduction on top of any
bare-bar UniaxialMaterial: sigma_buckled = r(eps, lsr) * sigma_bare. Tension and
elastic/pre-onset compression pass through unchanged.

Driver: a unit Truss (L=1, A=1) along x, displacement-controlled into
compression. strain = nodeDisp(2,1); stress = -nodeReaction(1,1).

Oracles (Ladruno_implementation/14_ladruno_rebar_buckling_adr.md):
  * B0 identity:   -lsr 0 ≡ the wrapped bare material (exact)
  * B1/B2 buckling: past the compression onset the buckled |stress| drops below
    the bare |stress|; elastic compression still matches bare
  * B3 residual+trend: deep compression floors near 0.2*fy; more slender (higher
    lsr) ⇒ lower compressive capacity at a fixed strain
  * B5 tangent:    analytic consistent tangent vs one-step-from-zero central FD
  * B6 composition: Fatigue ∘ Buckling ∘ LadrunoUniaxialJ2 builds and runs
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_FY, _E, _B = 420.0, 200000.0, 0.01     # typical rebar (MPa); Steel02 bare bar


def _steel02(tag):
    ops.uniaxialMaterial("Steel02", tag, _FY, _E, _B, 18.0, 0.925, 0.15)


# --------------------------------------------------------------------------
# compression truss driver
# --------------------------------------------------------------------------
def _build(mat_fn):
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
    ops.analysis("Static")


def _run(mat_fn, eps_target, nsteps):
    _build(mat_fn)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"analyze failed heading to eps={eps_target}"
        ops.reactions()
        out.append((ops.nodeDisp(2, 1), -ops.nodeReaction(1, 1)))
    return out


# --------------------------------------------------------------------------
# B0: -lsr 0 pass-through ≡ bare material (exact)
# --------------------------------------------------------------------------
@pytest.mark.t0m
def test_identity_no_buckling():
    def bare(t):
        _steel02(t)

    def wrapped(t):
        _steel02(101)
        ops.uniaxialMaterial("LadrunoRebarBuckling", t, 101, "-lsr", 0.0)

    b = _run(bare, -0.05, 100)
    w = _run(wrapped, -0.05, 100)
    for (eb, sb), (ew, sw) in zip(b, w):
        assert sw == pytest.approx(sb, rel=1e-9, abs=1e-9), f"identity broke at eps={eb}"


# --------------------------------------------------------------------------
# B1/B2: buckling reduces compression past onset; elastic compression matches bare
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_buckling_reduces_compression():
    def bare(t):
        _steel02(t)

    def wrapped(t):
        _steel02(101)
        ops.uniaxialMaterial("LadrunoRebarBuckling", t, 101, "-lsr", 10.0, "-fy", _FY)

    b = _run(bare, -0.04, 80)
    w = _run(wrapped, -0.04, 80)

    # elastic compression (first step, |eps|<eyp=fy/E≈0.0021): identical to bare
    assert w[0][1] == pytest.approx(b[0][1], rel=1e-6), "elastic compression must match bare"
    # deep compression: buckled |stress| is strictly below the bare |stress|
    assert abs(w[-1][1]) < abs(b[-1][1]) - 1.0, (
        f"buckling should reduce compression: |w|={abs(w[-1][1]):.1f} vs |b|={abs(b[-1][1]):.1f}")
    # and the bare bar is genuinely yielded/hardening (guard against vacuous compare)
    assert abs(b[-1][1]) > _FY


# --------------------------------------------------------------------------
# B3: residual floor (~0.2 fy) + slenderness trend (higher lsr ⇒ lower capacity)
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_residual_and_slenderness_trend():
    cap = {}
    for lsr in (6.0, 10.0, 15.0):
        def wrapped(t, lsr=lsr):
            _steel02(101)
            ops.uniaxialMaterial("LadrunoRebarBuckling", t, 101, "-lsr", lsr, "-fy", _FY)
        cap[lsr] = abs(_run(wrapped, -0.06, 120)[-1][1])

    # more slender ⇒ lower compressive capacity at the same deep strain
    assert cap[15.0] < cap[10.0] < cap[6.0], f"slenderness trend broken: {cap}"
    # the most slender floors near the 0.2*fy residual
    assert cap[15.0] == pytest.approx(0.2 * _FY, rel=0.3), (
        f"residual {cap[15.0]:.1f} not near 0.2 fy = {0.2*_FY:.1f}")


# --------------------------------------------------------------------------
# B5: analytic consistent tangent vs one-step-from-zero central FD
# (setStrain commits, so probe each point as an independent one-step return;
#  see the J2 V6 rationale + LEDGER_quirks)
# --------------------------------------------------------------------------
def _one_step(eps, want):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    _steel02(101)
    ops.uniaxialMaterial("LadrunoRebarBuckling", 1, 101, "-lsr", 10.0, "-fy", _FY)
    ops.testUniaxialMaterial(1)
    ops.setStrain(float(eps))
    return ops.getStress() if want == "stress" else ops.getTangent()


@pytest.mark.t0m
def test_consistent_tangent_fd():
    d = 1.0e-7
    # eyp≈0.0021, eStar(lsr=10)≈-0.0166: -0.008 is in the smooth intermediate branch
    for eps0 in (-0.008, -0.012):
        k = _one_step(eps0, "tangent")
        sp = _one_step(eps0 + d, "stress")
        sm = _one_step(eps0 - d, "stress")
        kfd = (sp - sm) / (2.0 * d)
        assert k == pytest.approx(kfd, rel=2e-3, abs=5.0), (
            f"buckling consistent tangent {k} != FD {kfd} at eps={eps0}")


# --------------------------------------------------------------------------
# B6: composition Fatigue ∘ Buckling ∘ LadrunoUniaxialJ2 builds and runs
# --------------------------------------------------------------------------
@pytest.mark.t1
def test_composition_fatigue_buckling_j2():
    def stack(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", 101, _E,
                             "-iso", "voce", _FY, 0.0, 0.0, _E * _B, "-kin", 1, 4000.0, 30.0)
        ops.uniaxialMaterial("LadrunoRebarBuckling", 102, 101, "-lsr", 10.0, "-fy", _FY)
        ops.uniaxialMaterial("Fatigue", t, 102)

    # drive into compression — the 3-layer stack must run without error and stay finite
    res = _run(stack, -0.03, 60)
    assert all(abs(s) < 1e9 for _, s in res), "composition produced non-finite stress"
    assert abs(res[-1][1]) < _FY * 2.0, "composition stress unreasonable"
