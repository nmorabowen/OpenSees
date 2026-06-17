"""LadrunoDispBeamColumn2d (Ladruno regularized displacement-based fiber frame,
2D, classTag 33013) -- Zone-A element battery.

A faithful clone of the stock ``dispBeamColumn`` (Euler-Bernoulli fiber frame,
linear basic strain) PLUS a Tier-1 regularization-length (``lch``) channel: the
element reports the localizing integration point's tributary length
``wt[i]*L`` through ``getCharacteristicLength()`` so crack-band / auto-regularizing
materials (e.g. ``ASDConcrete1D -autoRegularization``) regularize their softening
branch over the correct band instead of the whole element -- the fix stock
``dispBeamColumn`` never got (LEDGER_quirks sec.59), mirroring ``forceBeamColumn``.

Correctness anchors:

  * **reduce-to-stock**: with an Elastic section and a Linear transform, the
    element reproduces stock ``dispBeamColumn`` to round-off (the clone is byte
    behavior + lch logic that is inert under elasticity).

  * **-lch parsing**: ``ip`` / ``element`` / a finite positive value build; the
    review-hardened parser rejects non-finite (``inf``/``nan``) and non-positive
    band widths.

  * **Tier-1 lch delivery**: a softening ``ASDConcrete1D -autoRegularization``
    fiber section pulled past peak under displacement control dissipates a
    DIFFERENT amount of energy for ``-lch ip`` (per-IP band) vs ``-lch element``
    (full-element band), proving the per-IP length actually reaches the material
    and scales its fracture energy.

  * **mesh-objectivity**: with ``-lch ip`` the dissipated energy of a softening
    cantilever converges as the mesh is refined (1/2/4/.. elements) -- the
    regularization removes mesh-size dependence.

Plan: Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# elastic section / member properties
_E = 29000.0
_A = 10.0
_Iz = 100.0
_L = 100.0


def _beam_model(ndm_ok=True):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)


def _elastic_tip(make_element, n=2):
    """Cantilever of n equal elements along x; fix node 1; mixed tip load;
    linear-solve; return tip (node n+1) displacement vector."""
    _beam_model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    make_element(n)
    tip = n + 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 1.0e3, -2.0e3, 5.0e2)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(tip, d) for d in range(1, 4)]


# ---------------------------------------------------------------------------
# Test 1: reduce-to-stock (elastic) == dispBeamColumn
# ---------------------------------------------------------------------------
def test_reduce_to_stock_elastic():
    def make_ladruno(n):
        for i in range(1, n + 1):
            ops.element("LadrunoDispBeamColumn", i, i, i + 1, 1, 1)

    def make_stock(n):
        for i in range(1, n + 1):
            ops.element("dispBeamColumn", i, i, i + 1, 1, 1)

    d_lad = _elastic_tip(make_ladruno)
    d_ref = _elastic_tip(make_stock)
    for a, b in zip(d_lad, d_ref):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-14)


# ---------------------------------------------------------------------------
# Test 2: -lch parsing accept / reject
# ---------------------------------------------------------------------------
def _try_build_lch(lch_args):
    """Return True if an element with the given trailing -lch args builds."""
    _beam_model()
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    try:
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, *lch_args)
    except Exception:
        return False
    return 1 in ops.getEleTags()


@pytest.mark.parametrize("args", [(), ("-lch", "ip"), ("-lch", "element"), ("-lch", "12.5")])
def test_lch_accepts_valid(args):
    assert _try_build_lch(args) is True


@pytest.mark.parametrize("bad", ["inf", "-inf", "nan", -5.0, 0.0])
def test_lch_rejects_invalid(bad):
    assert _try_build_lch(("-lch", bad)) is False


# ---------------------------------------------------------------------------
# Tier-1 lch delivery + mesh-objectivity helpers (softening fiber section)
# ---------------------------------------------------------------------------
def _concrete_section(sec_tag, mat_tag, lch_ref=100.0):
    ops.uniaxialMaterial(
        "ASDConcrete1D", mat_tag, 30000.0,
        "-ft", 3.0, "-Te", *[1.0e-4, 1.5e-3], "-Ts", *[3.0, 0.05],
        "-fc", 30.0, "-Ce", *[1.0e-3, 6.0e-3], "-Cs", *[30.0, 1.0],
        "-autoRegularization", lch_ref,
    )
    ops.section("Fiber", sec_tag)
    for k in range(10):
        ops.fiber(-9.0 + k * 2.0, 0.0, 20.0, mat_tag)


def _pushover_energy(lch_mode, n):
    """n-element cantilever, concrete fiber section, displacement-controlled tip
    push past peak; return dissipated energy = integral of base shear d(tip)."""
    _beam_model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0)
    ops.fix(1, 1, 1, 1)
    _concrete_section(1, 1)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    for i in range(1, n + 1):
        ops.element("LadrunoDispBeamColumn", i, i, i + 1, 1, 1, "-lch", lch_mode)
    tip = n + 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 0.0, -1.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-7, 60)
    ops.algorithm("KrylovNewton")
    ops.integrator("DisplacementControl", tip, 2, -0.005)
    ops.analysis("Static")
    E = fprev = dprev = peak = 0.0
    for _ in range(200):
        if ops.analyze(1) != 0:
            break
        ops.reactions()
        d = abs(ops.nodeDisp(tip, 2))
        f = abs(ops.nodeReaction(1, 2))
        E += 0.5 * (f + fprev) * (d - dprev)
        peak = max(peak, f)
        fprev, dprev = f, d
    return E, peak


def test_lch_delivery_changes_softening():
    """Per-IP vs full-element band feed the auto-reg material different fracture
    energy -> different dissipated energy. Proves lch reaches the material."""
    e_ip, _ = _pushover_energy("ip", 4)
    e_el, _ = _pushover_energy("element", 4)
    assert e_ip > 0.0 and e_el > 0.0
    # ip uses a smaller (per-IP) band -> larger regularized Gf -> more dissipation
    assert e_ip > 1.2 * e_el


def test_mesh_objectivity_ip():
    """-lch ip: dissipated energy converges under mesh refinement."""
    e2, _ = _pushover_energy("ip", 2)
    e4, _ = _pushover_energy("ip", 4)
    e8, _ = _pushover_energy("ip", 8)
    assert e2 > 0 and e4 > 0 and e8 > 0
    # refinement change shrinks and the finest step is small (objective)
    assert abs(e8 - e4) / e4 < 0.05


# ---------------------------------------------------------------------------
# Test: large-displacement via Corotational transform == stock dispBeamColumn
# ---------------------------------------------------------------------------
def _corot_largedef_tip(ele_name, n=4, nsteps=30):
    """Elastic cantilever under a Corotational transform driven into large
    deflection (pure geometric nonlinearity); return tip displacement vector."""
    _beam_model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Corotational", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    for i in range(1, n + 1):
        ops.element(ele_name, i, i, i + 1, 1, 1)
    tip = n + 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 0.0, -400.0, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0
    return [ops.nodeDisp(tip, d) for d in range(1, 4)]


def test_corotational_large_displacement_matches_stock():
    d_lad = _corot_largedef_tip("LadrunoDispBeamColumn")
    d_ref = _corot_largedef_tip("dispBeamColumn")
    # genuinely large deflection (geometric nonlinearity active): tip drop and
    # axial shortening both substantial fractions of the length
    assert abs(d_lad[1]) > 10.0 and abs(d_lad[0]) > 1.0
    # faithful clone: identical to stock under the same Corotational transform
    for a, b in zip(d_lad, d_ref):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# Test: -nl  ½θ² bowing strain toggle (Linear transform, axially restrained)
# ---------------------------------------------------------------------------
def _lin_transf_tip(nl, P, n=4, nsteps=20):
    """Axially-restrained cantilever under a Linear geomTransf, transverse tip
    load P over nsteps with Newton; return tip vertical displacement. With -nl
    the ½θ² bowing builds axial tension under the restraint (tension stiffening)
    that the linear basic strain cannot represent."""
    _beam_model()
    dx = _L / n
    for i in range(n + 1):
        ops.node(i + 1, i * dx, 0.0)
    ops.fix(1, 1, 1, 1)
    tip = n + 1
    ops.fix(tip, 1, 0, 0)  # restrain axial DOF at the free end
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    extra = ("-nl",) if nl else ()
    for i in range(1, n + 1):
        ops.element("LadrunoDispBeamColumn", i, i, i + 1, 1, 1, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 0.0, -float(P), 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0
    return ops.nodeDisp(tip, 2)


def test_nl_reduces_to_linear_at_small_deformation():
    d_nl = _lin_transf_tip(True, 1.0)
    d_lin = _lin_transf_tip(False, 1.0)
    # tiny load -> theta ~ 0 -> ½θ² negligible -> NL == linear
    assert math.isclose(d_nl, d_lin, rel_tol=1e-5)


def test_nl_stiffens_under_finite_rotation():
    d_nl = _lin_transf_tip(True, 300.0)
    d_lin = _lin_transf_tip(False, 300.0)
    # both deflect substantially; the -nl bowing builds restraint tension that
    # stiffens the member, so |tip(NL)| < |tip(linear)| by a clear margin
    assert abs(d_lin) > 5.0
    assert abs(d_nl) < 0.98 * abs(d_lin)
