"""LadrunoShellModifier (WP-91) -- ETABS-style stiffness-modifier decorator for any
order-8 plate section, Zone-A section/element-level battery.

    section LadrunoShellModifier $tag $innerSecTag
        <-f11 v> <-f22 v> <-f12 v> <-m11 v> <-m22 v> <-m12 v>
        <-v13 v> <-v23 v> <-mass v>

All nine flags optional, default 1.0 (there is no weight modifier). Index mapping onto
the OpenSees plate resultant order (0..7 = FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ) is
f11,f22,f12,m11,m22,m12,v13,v23. With s_i = sqrt(m_i), S = diag(s):

    eps_inner = S @ eps_outer                  (inner sees scaled strain)
    sig_outer = S @ sig_inner                  (outer stress = scaled inner stress)
    D_outer   = S @ D_inner @ S                 i.e. D'(i,j) = sqrt(m_i*m_j) * D(i,j)
    rho_outer = massMod * rho_inner
    getSectionDeformation() returns the OUTER strain (unchanged pass-through).

ElasticMembranePlateSection stores its plate (bending/shear) block NEGATIVE
(tangent(3,3) = -D, etc. -- see SRC/material/section/ElasticMembranePlateSection.cpp).
The congruence transform preserves that sign (S is diagonal-positive, congruence
cannot flip an eigen-sign), so every closed-form D used below carries the SAME
negative bending/shear block the vanilla section does -- do not "fix" the sign.

Gates in this file (see Ladruno_implementation ledgers for the WP-91 feature row):
  G1  identity            all nine modifiers = 1 <=> bit-identical to the bare inner
                          section, on ShellMITC4 AND ASDShellQ4.
  G2  congruence          one non-trivial modifier set, direct numpy verification of
                          sig_outer = S @ D_inner @ S @ eps_outer against the FE's own
                          strain read-back (eleResponse 'strains') -- see the note in
                          test_g2_congruence_all_nine_modifiers below for why this file
                          reads the actual FE strain rather than hand-deriving a
                          "uniform patch" nodal field (a deliberate, stronger adaptation
                          of the literal gate wording).
  G3  analytic structural  a uniform block-scalar exactly halves/doubles the reduced
                          stiffness of an isolated membrane- or bending-active DOF set,
                          so tip deflection scales by EXACTLY 1/r on a cantilever strip.
  G4  SPD                 f11=f22=0.01, f12=1.0, nu=0.25 keeps the membrane congruence
                          block SPD; a numpy-only companion shows the naive
                          diagonal-only scaling of the same block would NOT be.
  G5  Ep_mod equivalence  m11=m22=m12=v13=v23=r (f=1) on an Ep_mod=1 inner section must
                          reproduce ElasticMembranePlateSection's native Ep_mod=r exactly
                          (confirmed against the .cpp source below -- see the docstring
                          of test_g5_ep_mod_equivalence for the derivation).
  G7  refusals            non-order-8 inner, missing inner tag, negative modifier all
                          raise; a modifier of exactly 0.0 is accepted and produces a
                          singular (not erroring) mode.
  G8  nonlinear passthrough  wrapping a LayeredShell(LadrunoRCConcrete) section with all
                          modifiers = 1 must reproduce the bare section's nonlinear
                          push-to-cracking response exactly.

G6 (mass) lives in test_ladrunoShellModifier_structural.py.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


@pytest.fixture(scope="module", autouse=True)
def _build_stamp():
    """Fail loud on a stale/openseespy-wheel binary rather than silently passing --
    LadrunoShellModifier is a fork-only feature (see LEDGER_quirks / CLAUDE.md)."""
    if not hasattr(ops, "ladrunoBuild"):
        pytest.skip("openseespy wheel fallback in use -- LadrunoShellModifier needs the fork build")
    print(f"\n[ladrunoBuild] {ops.ladrunoBuild()}")


# ---------------------------------------------------------------------------
# shared constants / helpers
# ---------------------------------------------------------------------------
E0, NU0, H0, RHO0 = 200.0, 0.25, 0.2, 2.5
FIVE6 = 5.0 / 6.0

_QUAD = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0)}


def _mod_args(**mods):
    args = []
    for k, v in mods.items():
        args += [f"-{k}", float(v)]
    return args


def _wrap(tag, inner_tag, **mods):
    ops.section("LadrunoShellModifier", tag, inner_tag, *_mod_args(**mods))


def _bare_section(tag, E=E0, nu=NU0, h=H0, rho=RHO0, ep_mod=1.0):
    # ElasticMembranePlateSection ctor order: tag, E, nu, h, <rho>, <Ep_mod>
    ops.section("ElasticMembranePlateSection", tag, E, nu, h, rho, ep_mod)


def _closed_form_D(E, nu, h, ep_mod=1.0):
    """8x8 D from SRC/material/section/ElasticMembranePlateSection.cpp::getSectionTangent,
    INCLUDING its negative bending/shear sign convention. Row/col order matches the
    section's own getType(): FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ."""
    Em = E
    Ep = E * ep_mod
    M = Em / (1.0 - nu * nu) * h
    Gm = 0.5 * Em / (1.0 + nu) * h            # membrane (in-plane) shear stiffness
    Gt = Gm * (FIVE6 * (Ep / Em))             # transverse-shear stiffness (uses Ep)
    D = Ep * h ** 3 / 12.0 / (1.0 - nu * nu)  # bending modulus magnitude
    K = np.zeros((8, 8))
    K[0, 0] = M
    K[1, 1] = M
    K[0, 1] = K[1, 0] = nu * M
    K[2, 2] = Gm
    K[3, 3] = -D
    K[4, 4] = -D
    K[3, 4] = K[4, 3] = -nu * D
    K[5, 5] = -0.5 * D * (1.0 - nu)
    K[6, 6] = Gt
    K[7, 7] = Gt
    return K


def _build_quad(ele, sec_tag, sec_builder, nodes=_QUAD):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in nodes.items():
        ops.node(t, x, y, z)
    sec_builder(sec_tag)
    ops.element(ele, 1, *nodes.keys(), sec_tag)


def _static_setup():
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


# ===========================================================================
# G1 -- identity gate
# ===========================================================================
# arbitrary but fixed nodal loads (forces + moments, all six DOFs) at the three
# unrestrained corner nodes -- exercises membrane, bending, transverse shear AND
# the drilling DOF simultaneously.
_G1_LOADS = {
    2: (1.3, -0.7, 0.4, 0.2, -0.5, 0.1),
    3: (-0.6, 0.9, -0.3, 0.15, 0.25, -0.2),
    4: (0.5, 0.2, 0.6, -0.1, 0.3, 0.05),
}


def _g1_solve():
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy, fz, mx, my, mz) in _G1_LOADS.items():
        ops.load(n, fx, fy, fz, mx, my, mz)
    _static_setup()
    assert ops.analyze(1) == 0, "G1 identity solve failed"
    K = np.array(ops.printA("-ret"))
    F = np.array(ops.eleForce(1))
    return K, F


def _run_identity(ele):
    _build_quad(ele, 1, lambda t: _bare_section(t))
    K_bare, F_bare = _g1_solve()

    def _wrapped_builder(t):
        _bare_section(100)
        _wrap(t, 100, f11=1.0, f22=1.0, f12=1.0, m11=1.0, m22=1.0, m12=1.0,
              v13=1.0, v23=1.0, mass=1.0)

    _build_quad(ele, 2, _wrapped_builder)
    K_wrap, F_wrap = _g1_solve()

    # "ideally bit-identical" -- sqrt(1.0) == 1.0 exactly in IEEE754 and 1.0*x == x
    # exactly, so a real bit-identical implementation should pass at rtol=0; we assert
    # at 1e-13 to absorb any incidental fp reordering inside the wrapper's own math.
    np.testing.assert_allclose(K_wrap, K_bare, rtol=1e-13, atol=1e-13,
                                err_msg=f"{ele}: wrapped(all=1) K != bare K")
    np.testing.assert_allclose(F_wrap, F_bare, rtol=1e-13, atol=1e-13,
                                err_msg=f"{ele}: wrapped(all=1) F != bare F")


def test_g1_identity_shellmitc4():
    """All nine modifiers = 1 must be bit-identical to the bare inner section,
    both the assembled free-DOF tangent (printA, single-element patch -- ShellMITC4/
    ASDShellQ4 do not implement eleResponse('stiff') the way the Ladruno solid
    elements do, so this is the closest available proxy for 'element stiffness
    matrix') and the resisting force vector (eleForce)."""
    _run_identity("ShellMITC4")


def test_g1_identity_asdshellq4():
    _run_identity("ASDShellQ4")


# ===========================================================================
# G2 -- congruence gate
# ===========================================================================
_MODS_G2 = dict(f11=0.35, f22=0.70, f12=0.25, m11=0.50, m22=0.90, m12=0.15,
                v13=0.80, v23=0.60, mass=0.75)
_G2_ORDER = ("f11", "f22", "f12", "m11", "m22", "m12", "v13", "v23")

# a deliberately generic (non-uniform, non-"patch") nodal displacement field --
# see the docstring below for why this file does not attempt to hand-derive an
# exact constant-curvature/constant-shear MITC4 patch field.
_G2_DISP = {
    2: (0.0040, -0.0020, 0.0030, 0.0010, -0.0015, 0.0008),
    3: (-0.0030, 0.0025, -0.0018, -0.0012, 0.0009, -0.0006),
    4: (0.0020, 0.0015, 0.0021, 0.0007, -0.0011, 0.0013),
}


def test_g2_congruence_all_nine_modifiers():
    """Direct verification of sig_outer = S @ D_inner @ S @ eps_outer (equivalently
    D_outer = S @ D_inner @ S) for a modifier set where all nine values differ.

    ADAPTATION FROM THE LITERAL GATE WORDING: rather than hand-deriving a nodal
    displacement field that forces a "known uniform" curvature+shear state through
    ShellMITC4's MITC tying (which requires reverse-engineering its B-matrix sign
    convention -- SRC/element/shell/ShellMITC4.cpp computeBdrill/Bmembrane/Bbend --
    to get right), this test imposes an arbitrary, non-uniform nodal displacement
    field and reads the ACTUAL section strain the element computed at Gauss point 0
    via eleResponse(1, "strains"). The congruence identity
        sig_outer = S @ D_inner @ S @ eps_outer
    is a pointwise material law (D_inner is constant/elastic) that must hold for
    ANY eps, uniform or not -- so this is a strictly more general check than the
    prompt's literal "uniform state" construction, without the patch-test-fidelity
    risk. D_inner itself comes from the closed-form _closed_form_D (independently
    derived from the .cpp source, see its docstring), not from the code under test.
    """
    E, nu, h, rho = 190.0, 0.22, 0.18, 1.7
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _QUAD.items():
        ops.node(t, x, y, z)
    _bare_section(100, E=E, nu=nu, h=h, rho=rho, ep_mod=1.0)
    _wrap(1, 100, **_MODS_G2)
    ops.element("ShellMITC4", 1, *_QUAD.keys(), 1)

    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, vals in _G2_DISP.items():
        for d, v in enumerate(vals, start=1):
            ops.sp(n, d, v)
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "G2 imposed-displacement solve failed"

    ops.eleResponse(1, "forces")  # prime lazy state (see test_ladrunoRCConcrete_shell.py)
    eps = np.array(ops.eleResponse(1, "strains"))[:8]
    sig_actual = np.array(ops.eleResponse(1, "stresses"))[:8]

    D_inner = _closed_form_D(E, nu, h, ep_mod=1.0)
    s = np.sqrt([_MODS_G2[k] for k in _G2_ORDER])
    S = np.diag(s)
    D_outer = S @ D_inner @ S
    sig_expected = D_outer @ eps

    np.testing.assert_allclose(
        sig_actual, sig_expected, rtol=1.0e-7, atol=1.0e-7,
        err_msg=f"congruence mismatch: eps={eps}\nactual={sig_actual}\nexpected={sig_expected}",
    )

    # bonus (best-effort, not load-bearing): if LadrunoShellModifier does not override
    # setResponse it inherits SectionForceDeformation's "stiffness" keyword, which
    # returns the flattened 8x8 getSectionTangent(). If that keyword is not wired up
    # this block is skipped rather than failing the gate above -- the frozen contract
    # does not specify setResponse behavior, so this is treated as a nice-to-have.
    try:
        raw = ops.eleResponse(1, "material", 1, "stiffness")
    except Exception:
        raw = None
    if raw is not None and len(raw) == 64:
        D_fe = np.array(raw).reshape(8, 8)
        np.testing.assert_allclose(D_fe, D_outer, rtol=1.0e-8, atol=1.0e-8,
                                    err_msg="eleResponse('material',1,'stiffness') tangent "
                                            "!= closed-form S@D@S")


# ===========================================================================
# G3 -- analytic structural gate (exact block-scalar scaling)
# ===========================================================================
_NX = 4  # cantilever strip: NX elements along x, width 1 along y


def _strip_nodes(nx=_NX):
    return {(i, j): i * 2 + j + 1 for i in range(nx + 1) for j in (0, 1)}


def _build_strip(ele, sec_tag, sec_builder, nx=_NX):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nodes = _strip_nodes(nx)
    for (i, j), tag in nodes.items():
        ops.node(tag, float(i), float(j), 0.0)
    sec_builder(sec_tag)
    for i in range(nx):
        n1, n2, n3, n4 = nodes[(i, 0)], nodes[(i + 1, 0)], nodes[(i + 1, 1)], nodes[(i, 1)]
        ops.element(ele, i + 1, n1, n2, n3, n4, sec_tag)
    return nodes


def _membrane_tip_disp(sec_builder, nx=_NX):
    """Restrain uz/rotx/roty/drilling EVERYWHERE (isolate the membrane block
    entirely) and load the tip with in-plane force -- the reduced system only ever
    sees K_membrane, so scaling the whole 3x3 membrane D block by r scales this
    displacement by EXACTLY 1/r."""
    nodes = _build_strip("ShellMITC4", 1, sec_builder, nx)
    for (i, j), tag in nodes.items():
        if i == 0:
            ops.fix(tag, 1, 1, 1, 1, 1, 1)
        else:
            ops.fix(tag, 0, 0, 1, 1, 1, 1)
    tip = [nodes[(nx, 0)], nodes[(nx, 1)]]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in tip:
        ops.load(t, 10.0, 4.0, 0.0, 0.0, 0.0, 0.0)
    _static_setup()
    assert ops.analyze(1) == 0, "G3 membrane solve failed"
    return ops.nodeDisp(tip[0], 1), ops.nodeDisp(tip[0], 2)


def _bending_tip_disp(sec_builder, nx=_NX):
    """Restrain ux/uy/drilling EVERYWHERE (isolate bending+shear entirely) and load
    the tip out-of-plane -- scaling the ENTIRE out-of-plane (bending+shear) 5x5 D
    block by r scales this displacement by EXACTLY 1/r, regardless of the
    bending/shear split of the response."""
    nodes = _build_strip("ShellMITC4", 1, sec_builder, nx)
    for (i, j), tag in nodes.items():
        if i == 0:
            ops.fix(tag, 1, 1, 1, 1, 1, 1)
        else:
            ops.fix(tag, 1, 1, 0, 0, 0, 1)
    tip = [nodes[(nx, 0)], nodes[(nx, 1)]]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in tip:
        ops.load(t, 0.0, 0.0, 5.0, 0.3, 0.2, 0.0)
    _static_setup()
    assert ops.analyze(1) == 0, "G3 bending solve failed"
    return ops.nodeDisp(tip[0], 3)


def test_g3_membrane_modifier_exactly_doubles_inplane_tip():
    ux_bare, uy_bare = _membrane_tip_disp(lambda t: _bare_section(t))

    def _mod(t):
        _bare_section(t + 100)
        _wrap(t, t + 100, f11=0.5, f22=0.5, f12=0.5, m11=1.0, m22=1.0, m12=1.0,
              v13=1.0, v23=1.0, mass=1.0)

    ux_mod, uy_mod = _membrane_tip_disp(_mod)
    assert ux_bare != 0.0 and uy_bare != 0.0
    assert ux_mod / ux_bare == pytest.approx(2.0, rel=1.0e-9)
    assert uy_mod / uy_bare == pytest.approx(2.0, rel=1.0e-9)


def test_g3_bending_modifier_exactly_doubles_outofplane_tip():
    uz_bare = _bending_tip_disp(lambda t: _bare_section(t))

    def _mod(t):
        _bare_section(t + 100)
        _wrap(t, t + 100, f11=1.0, f22=1.0, f12=1.0, m11=0.5, m22=0.5, m12=0.5,
              v13=0.5, v23=0.5, mass=1.0)

    uz_mod = _bending_tip_disp(_mod)
    assert uz_bare != 0.0
    assert uz_mod / uz_bare == pytest.approx(2.0, rel=1.0e-9)


# ===========================================================================
# G4 -- SPD gate
# ===========================================================================
def _membrane_patch_K(sec_builder, ele="ShellMITC4"):
    """Single unit-quad element, node 1 pinned, nodes 2-4 restrained to their
    membrane translations only (uz/rotx/roty/drilling fixed everywhere) -- a 6-DOF
    reduced system driven purely by the section's 3x3 membrane block."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _QUAD.items():
        ops.node(t, x, y, z)
    sec_builder(1)
    ops.element(ele, 1, *_QUAD.keys(), 1)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    for t in (2, 3, 4):
        ops.fix(t, 0, 0, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (2, 3, 4):
        ops.load(t, 0.01, 0.007, 0.0, 0.0, 0.0, 0.0)
    _static_setup()
    try:
        ops.analyze(1)          # may legitimately fail to SOLVE for the G7 zero-modifier
    except Exception:           # case below -- the tangent is still formed beforehand
        pass
    K = np.array(ops.printA("-ret"))
    n = int(round(len(K) ** 0.5))
    return K.reshape(n, n)


def test_g4_membrane_congruence_is_spd():
    def build(tag):
        _bare_section(tag + 100, nu=0.25)
        _wrap(tag, tag + 100, f11=0.01, f22=0.01, f12=1.0,
              m11=1.0, m22=1.0, m12=1.0, v13=1.0, v23=1.0, mass=1.0)

    K = _membrane_patch_K(build)
    vals = np.linalg.eigvalsh(K)
    assert np.all(vals > 0), f"congruence membrane stiffness not SPD: eigenvalues={vals}"


def test_g4_companion_naive_diagonal_scaling_would_be_indefinite():
    """DESIGN-RATIONALE check, not a gate on the C++ code: pure numpy, demonstrating
    WHY the congruence form D'(i,j) = sqrt(m_i*m_j)*D(i,j) was chosen over a naive
    "scale the diagonal only, leave off-diagonal terms as-is" implementation. At
    f11=f22=0.01, nu=0.25 the naive Fxx/Fyy 2x2 submatrix [[0.01M, nu*M],[nu*M, 0.01M]]
    has eigenvalues 0.01M*(1 +/- nu/0.01) = M*(0.26, -0.24) -- indefinite -- while the
    congruence form [[0.01M, 0.01*nu*M],[0.01*nu*M, 0.01M]] = 0.01*[[M,nu*M],[nu*M,M]]
    stays SPD (it is just the original SPD block uniformly rescaled)."""
    nu, M = 0.25, 1.0  # M is an arbitrary positive scale (E*h/(1-nu^2)); sign only matters
    naive = np.array([[0.01 * M, nu * M], [nu * M, 0.01 * M]])
    congruent = np.array([[0.01 * M, 0.01 * nu * M], [0.01 * nu * M, 0.01 * M]])
    assert np.linalg.eigvalsh(naive).min() < 0.0, "companion check invalid: naive stayed SPD"
    assert np.linalg.eigvalsh(congruent).min() > 0.0


# ===========================================================================
# G5 -- upstream Ep_mod equivalence gate
# ===========================================================================
_G5_DISP = {
    2: (0.0021, -0.0011, 0.0031, 0.0009, -0.0007, 0.0004),
    3: (-0.0017, 0.0026, -0.0022, -0.0006, 0.0013, -0.0009),
    4: (0.0012, 0.0008, 0.0019, 0.0005, -0.0010, 0.0011),
}


def _g5_solve_generic(sec_builder):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _QUAD.items():
        ops.node(t, x, y, z)
    sec_builder(1)
    ops.element("ShellMITC4", 1, *_QUAD.keys(), 1)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, vals in _G5_DISP.items():
        for d, v in enumerate(vals, start=1):
            ops.sp(n, d, v)
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "G5 imposed-displacement solve failed"
    K = np.array(ops.printA("-ret"))
    F = np.array(ops.eleForce(1))
    return K, F


def test_g5_ep_mod_equivalence():
    """Derivation (see ElasticMembranePlateSection::getSectionTangent): with an
    Ep_mod=1 inner section (Ep=Em), wrapping with m11=m22=m12=v13=v23=r (f=1)
    scales D(3,3),D(4,4) by r, D(3,4) by sqrt(r*r)=r, D(5,5) by r, and D(6,6)/D(7,7)
    (which are proportional to Ep) by r as well -- termwise IDENTICAL to setting
    Ep_mod=r directly on the vanilla section (Ep=r*Em produces exactly the same five
    terms, since each is linear in Ep). The membrane block is untouched in both
    constructions. No discrepancy from the frozen contract was found; this test
    pins the equivalence to solver precision via full element K + eleForce identity
    (same technique as G1), not just a hand-picked stress state."""
    r = 0.4

    def direct(tag):
        _bare_section(tag, ep_mod=r)

    def wrapped(tag):
        _bare_section(tag + 100, ep_mod=1.0)
        _wrap(tag, tag + 100, m11=r, m22=r, m12=r, v13=r, v23=r)

    K_direct, F_direct = _g5_solve_generic(direct)
    K_wrap, F_wrap = _g5_solve_generic(wrapped)
    np.testing.assert_allclose(K_wrap, K_direct, rtol=1.0e-11, atol=1.0e-11,
                                err_msg="Ep_mod=r K != wrapper-equivalent K")
    np.testing.assert_allclose(F_wrap, F_direct, rtol=1.0e-11, atol=1.0e-11,
                                err_msg="Ep_mod=r F != wrapper-equivalent F")


# ===========================================================================
# G7 -- refusal gates
# ===========================================================================
def test_g7_refuses_non_order8_inner_section():
    ops.wipe()
    # ndm 2 / ndf 3: `section Elastic $tag $E $A $Iz` is the 2d (order-2) form.
    # In ndm 3 that same command wants E A Iz Iy G J and errors on 3 args, which
    # would abort the test in SETUP rather than exercising the R1 refusal.
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.section("Elastic", 1, 1000.0, 1.0, 1.0)  # order-2 beam section, not order 8
    with pytest.raises(Exception):
        ops.section("LadrunoShellModifier", 2, 1, "-f11", 1.0)


def test_g7_refuses_missing_inner_section_tag():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    with pytest.raises(Exception):
        ops.section("LadrunoShellModifier", 1, 9999)


def test_g7_refuses_negative_modifier():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    _bare_section(1)
    with pytest.raises(Exception):
        ops.section("LadrunoShellModifier", 2, 1, "-f11", -0.5)


def test_g7_zero_modifier_accepted_and_singular_not_error():
    """A modifier of exactly 0.0 is ETABS-legal and must NOT raise; it must instead
    zero the corresponding row/col of D_outer, leaving the patch rank-deficient.

    The contract (ADR 91 §6) says 'singular in that response mode' -- it does NOT
    promise a rank drop of exactly one. Killing f11 zeroes the whole FXX row AND
    column, including the Poisson coupling D(0,1), so the section resists only
    eps_yy and gamma_xy; the space of patch deformations with eps_yy = gamma_xy = 0
    but eps_xx != 0 is more than one-dimensional for this 6-free-DOF quad (measured:
    2). So assert what the contract actually claims -- the wrapped patch is singular
    and the unmodified control patch is not."""
    def build_zero(tag):
        _bare_section(tag + 100)
        _wrap(tag, tag + 100, f11=0.0)   # no exception expected here

    def build_control(tag):
        _bare_section(tag + 100)
        _wrap(tag, tag + 100, f11=1.0)

    K = _membrane_patch_K(build_zero)
    vals = np.sort(np.linalg.eigvalsh(K))
    scale = vals[-1]
    n_zero = int(np.sum(np.abs(vals) < 1.0e-8 * scale))
    assert n_zero >= 1, f"expected f11=0.0 to make the patch singular: {vals}"

    K_ctl = np.linalg.eigvalsh(_membrane_patch_K(build_control))
    scale_ctl = np.max(K_ctl)
    n_zero_ctl = int(np.sum(np.abs(K_ctl) < 1.0e-8 * scale_ctl))
    assert n_zero_ctl == 0, (
        f"control patch (f11=1.0) must be non-singular, found {n_zero_ctl} "
        f"zero modes: {K_ctl}")


# ===========================================================================
# G8 -- nonlinear passthrough smoke (LayeredShell / LadrunoRCConcrete)
# ===========================================================================
# same backbone as tests/test_ladrunoRCConcrete_shell.py (already-verified material)
_RC_E, _RC_NU, _RC_KC = 30000.0, 0.2, 2.0 / 3.0
_RC_CE = [0.0, 0.0007, 0.0020, 0.0100]
_RC_CS = [0.0, 24.0, 30.0, 5.0]
_RC_CD = [0.0, 0.0, 0.25, 1.0 - 5.0 / 45.0]
_RC_TE = [0.0, 0.0001, 0.0010]
_RC_TS = [0.0, 3.0, 0.5]
_RC_TD = [0.0, 0.0, 1.0 - 0.5 / 5.0]
_RC_H = 0.1
_RC_NLAYERS = 4


def _rc_mat(tag):
    ops.nDMaterial("LadrunoRCConcrete", tag, _RC_E, _RC_NU,
                   "-Ce", *_RC_CE, "-Cs", *_RC_CS, "-Cd", *_RC_CD,
                   "-Te", *_RC_TE, "-Ts", *_RC_TS, "-Td", *_RC_TD, "-Kc", _RC_KC)


def test_g8_layeredshell_nonlinear_passthrough():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    coordsA = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0)}
    coordsB = {11: (0.0, 0.0, 0.0), 12: (1.0, 0.0, 0.0), 13: (1.0, 1.0, 0.0), 14: (0.0, 1.0, 0.0)}
    for group in (coordsA, coordsB):
        for t, (x, y, z) in group.items():
            ops.node(t, x, y, z)

    _rc_mat(1)
    layers = []
    for _ in range(_RC_NLAYERS):
        layers += [1, _RC_H / _RC_NLAYERS]
    ops.section("LayeredShell", 10, _RC_NLAYERS, *layers)          # bare (elementA)
    ops.section("LayeredShell", 20, _RC_NLAYERS, *layers)          # inner (elementB)
    ops.section("LadrunoShellModifier", 21, 20,
                "-f11", 1.0, "-f22", 1.0, "-f12", 1.0,
                "-m11", 1.0, "-m22", 1.0, "-m12", 1.0,
                "-v13", 1.0, "-v23", 1.0, "-mass", 1.0)
    ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 10)
    ops.element("ASDShellQ4", 2, 11, 12, 13, 14, 21)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # same per-step increment (1e-5) as the already-verified 60-step push in
    # test_ladrunoRCConcrete_shell.py, just fewer total steps (this is a passthrough
    # smoke test, not a re-verification of the material's softening physics)
    exx_max, nsteps = 3.5e-4, 35
    for group in (coordsA, coordsB):
        for t, (x, y, z) in group.items():
            ops.sp(t, 1, exx_max * x)
            ops.sp(t, 2, 0.0)
            ops.sp(t, 3, 0.0)
            ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-7, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")

    for step in range(nsteps):
        assert ops.analyze(1) == 0, f"G8 push step {step} failed"
        ops.eleResponse(1, "forces"); ops.eleResponse(2, "forces")  # prime lazy state
        nxx_bare = ops.eleResponse(1, "stresses")[0]
        nxx_wrap = ops.eleResponse(2, "stresses")[0]
        assert nxx_wrap == pytest.approx(nxx_bare, rel=1.0e-10, abs=1.0e-10), (
            f"step {step}: wrapped(all=1) LayeredShell Nxx {nxx_wrap} != bare {nxx_bare}"
        )
