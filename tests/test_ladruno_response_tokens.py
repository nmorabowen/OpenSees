"""Ladruno recorder-token consistency — `SRC/element/LadrunoResponseTokens.h`.

Upstream every element hand-rolls its own ``strcmp`` chain in ``setResponse``,
so the accepted spelling of the SAME quantity drifts between elements (vanilla
``Brick`` takes ``stresses`` but not ``stress``; ``FourNodeQuad`` takes both).
The fork inherited that drift into every new element, and added its own.

This battery locks in the three things the shared table buys:

  1. **Alias folding.** Every registered spelling of a quantity resolves to the
     same branch AND returns bit-identical values (``stress`` == ``stresses``,
     ``charLength`` == ``characteristicLength`` == ``lch``, ``kt`` == ``k`` ==
     ``penalty``, ...). A token that silently returns an empty response is the
     failure this exists to prevent — an empty recorder file looks like "the
     quantity is zero", not like "you spelled it wrong".

  2. **The family contract.** Every fork continuum element answers
     ``force / stress / strain / stiff / stiffInitial / charLength`` — several
     of them answered only a subset before (``LadrunoSolidShell`` had no force
     response at all; ``LadrunoUP`` had no force, strain or charLength;
     ``LadrunoCSTPair`` had no strain; the plane family had no stiffness).

  3. **The base-class chain.** No fork element used to delegate to
     ``Element::setResponse``, so ``inertialForce`` / ``dampingForce`` /
     ``dynamicForce`` were unavailable everywhere. The delegation must come
     AFTER ``output.endTag()`` because the base opens its own ``ElementOutput``
     tag (see LEDGER_quirks).

``eleResponse`` returns ``[]`` for an unresolved token, which is exactly the
signal we assert on.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# hex20 mid-edge pairs in brcshl order (0-based corners) — copied from
# tests/hex20_reference.py:EDGES rather than imported, so this battery does not
# drag in that module's sympy dependency.
_HEX20_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7)]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _resp(tag, *args):
    r = ops.eleResponse(tag, *args)
    return list(r) if r else []


def _assert_resolves(tag, token, minlen=1):
    v = _resp(tag, token)
    assert len(v) >= minlen, f"token {token!r} did not resolve on element {tag}"
    return v


def _assert_same(tag, *tokens):
    """All spellings resolve AND return the same values.

    Compared to a tight tolerance rather than bit-exactly: some elements'
    getters are not idempotent (``LadrunoSolidShell::getResistingForce``
    re-runs the EAS inner Newton, so two consecutive calls differ in the last
    ULP). Anything looser than a shared branch would be off by orders of
    magnitude, not by 1e-15.
    """
    ref = _assert_resolves(tag, tokens[0])
    for t in tokens[1:]:
        got = _assert_resolves(tag, t)
        assert len(got) == len(ref), (
            f"alias width mismatch on element {tag}: {t!r} gave {len(got)}, "
            f"{tokens[0]!r} gave {len(ref)}"
        )
        assert got == pytest.approx(ref, rel=1e-9, abs=1e-12), (
            f"alias mismatch on element {tag}: {t!r} != {tokens[0]!r}\n"
            f"  {t}: {got[:6]}\n  {tokens[0]}: {ref[:6]}"
        )
    return ref


def _static_solve():
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# The family contract every fork continuum element must honour.
_CONTINUUM = ("force", "stress", "strain", "stiff", "stiffInitial", "charLength")

# What the Element base class provides once an element chains to it.
_BASE_TOKENS = ("globalForce", "dampingForce", "dynamicForce", "inertialForce")


def _check_continuum_contract(tag, ndof):
    for token in _CONTINUUM:
        _assert_resolves(tag, token)
    # the singular/plural fold, both directions
    _assert_same(tag, "stress", "stresses")
    _assert_same(tag, "strain", "strains")
    _assert_same(tag, "force", "forces")
    _assert_same(tag, "charLength", "characteristicLength", "lch")
    _assert_same(tag, "stiff", "stiffness", "tangent")
    _assert_same(tag, "stiffInitial", "initialStiffness")
    assert len(_resp(tag, "force")) == ndof
    assert len(_resp(tag, "stiff")) == ndof * ndof


def _check_base_chain(tag, ndof):
    for token in _BASE_TOKENS:
        v = _assert_resolves(tag, token)
        assert len(v) == ndof, f"{token} on element {tag}: {len(v)} != {ndof}"


# ---------------------------------------------------------------------------
# 3D solids — LadrunoBrick / LadrunoBrick20 / LadrunoSolidShell
# ---------------------------------------------------------------------------
_HEX = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]


def _build_hex8(ele_name, *extra):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.node(i, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)
    ops.element(ele_name, 1, *range(1, 9), 1, *extra)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 1.0, 0.5, -2.0)
    _static_solve()


def test_brick_tokens():
    _build_hex8("LadrunoBrick")
    _check_continuum_contract(1, 24)
    _check_base_chain(1, 24)
    # brick-only quantities keep their own alias sets
    _assert_same(1, "hourglassEnergy", "hgEnergy", "hourglassDissipation")
    assert len(_resp(1, "stress")) == 48          # 8 GP x 6
    assert len(_resp(1, "stress3D6")) == 6
    assert _resp(1, "material", 1, "stress")      # per-GP forwarding still works


def test_brick20_tokens():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = {i + 1: _HEX[i] for i in range(8)}
    for k, (a, b) in enumerate(_HEX20_EDGES):     # 0-based corner pairs
        ca, cb = _HEX[a], _HEX[b]
        coords[9 + k] = tuple(0.5 * (ca[d] + cb[d]) for d in range(3))
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)
    ops.element("LadrunoBrick20", 1, *range(1, 21), 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 1.0, 0.5, -2.0)
    _static_solve()

    _check_continuum_contract(1, 60)
    _check_base_chain(1, 60)


def test_solidshell_tokens():
    """SolidShell had NO force and NO charLength response before the sweep."""
    _build_hex8("LadrunoSolidShell")
    _check_continuum_contract(1, 24)
    _check_base_chain(1, 24)
    _assert_same(1, "alpha", "eas")


# ---------------------------------------------------------------------------
# plane family — LadrunoQuad / LadrunoCST / LadrunoLST / LadrunoCSTPair
# ---------------------------------------------------------------------------
def _plane_common():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)


def _plane_load(top, fx=1.0, fy=-2.0):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, fx, fy)
    _static_solve()


def test_quad_tokens():
    _plane_common()
    for tag, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(tag, float(x), float(y))
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain",
                "-thick", 1.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    _plane_load((3, 4))

    _check_continuum_contract(1, 8)
    _check_base_chain(1, 8)
    _assert_same(1, "stressPlaneStrain", "stressesPlaneStrain")


def test_cst_tokens():
    _plane_common()
    for tag, (x, y) in enumerate([(0, 0), (1, 0), (0, 1)], start=1):
        ops.node(tag, float(x), float(y))
    ops.element("LadrunoCST", 1, 1, 2, 3, 1, "-type", "PlaneStrain",
                "-thick", 1.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    _plane_load((3,))

    _check_continuum_contract(1, 6)
    _check_base_chain(1, 6)


def test_lst_tokens():
    _plane_common()
    corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    mids = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
    for tag, (x, y) in enumerate(corners + mids, start=1):
        ops.node(tag, x, y)
    ops.element("LadrunoLST", 1, *range(1, 7), 1, "-type", "PlaneStrain",
                "-thick", 1.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(4, 1, 1)
    _plane_load((3, 5, 6))

    _check_continuum_contract(1, 12)
    _check_base_chain(1, 12)


def test_cstpair_tokens():
    """CSTPair had stress but NO strain response before the sweep."""
    _plane_common()
    ops.nDMaterial("LogStrain2D", 2, 1)     # finite-strain only, by contract
    for tag, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(tag, float(x), float(y))
    ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 2, "-type", "PlaneStrain",
                "-thick", 1.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    _plane_load((3, 4))

    _check_continuum_contract(1, 8)
    _check_base_chain(1, 8)
    _assert_same(1, "Jbar", "jbar")
    assert len(_resp(1, "strain")) == 6           # 2 triangles x 3 components


# ---------------------------------------------------------------------------
# LadrunoUP — had no force / strain / stiff / charLength at all
# ---------------------------------------------------------------------------
def test_up_tokens():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for tag, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 2.2e6, "-poro", 0.4,
                "-rhoF", 1.0, "-perm", 1.0e-4, 1.0e-4)
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 1, 1, 0)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 0, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, -2.0, 0.0)
    ops.load(4, 1.0, -2.0, 0.0)
    _static_solve()

    _check_continuum_contract(1, 12)
    _check_base_chain(1, 12)
    _assert_same(1, "stressTotal", "stressesTotal")
    _assert_same(1, "flux", "darcyFlux")
    assert _resp(1, "porePressure")


# ---------------------------------------------------------------------------
# LadrunoIMKBeam — realigned on the DispBeamColumn vocabulary
# ---------------------------------------------------------------------------
# INCLINED on purpose: on a beam along +x the local and global frames coincide,
# so a "localForce" that secretly returns the global vector would pass. 3-4-5.
_IMK2D_L = 5.0


def _imk2d():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 3.0, 4.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.element("LadrunoIMKBeam", 1, 1, 2, 0.01, 2.0e8, 1.0e-4, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, -20.0, 5.0)
    _static_solve()


def test_imkbeam2d_tokens():
    _imk2d()
    # beam-family spellings, all newly accepted
    _assert_same(1, "force", "forces", "globalForce", "globalForces")
    _assert_same(1, "basicForce", "basicForces")
    _assert_same(1, "chordRotation", "chordDeformation", "basicDeformation",
                 "deformation", "deformations")
    _assert_same(1, "plasticRotation", "plasticDeformation", "hingeRotation",
                 "hingeRotations")
    _assert_same(1, "hingeMoment", "hingeMoments")
    _check_base_chain(1, 6)

    # localForce is now the ELEMENT-FRAME force, not a third spelling of the
    # global one. On a beam along +x the two frames coincide in magnitude but
    # the local vector is the 6 basic-derived end forces, and equilibrium must
    # hold: N_i = -N_j, and the two end moments balance the shear couple.
    loc = _assert_resolves(1, "localForce")
    glo = _resp(1, "force")
    assert len(loc) == 6 and len(glo) == 6
    assert loc != pytest.approx(glo, rel=1e-6), \
        "localForce is still returning the global vector"
    q = _resp(1, "basicForce")
    assert abs(loc[0] + q[0]) < 1e-9 and abs(loc[3] - q[0]) < 1e-9
    assert abs(loc[2] - q[1]) < 1e-9 and abs(loc[5] - q[2]) < 1e-9
    assert abs(loc[1] + loc[4]) < 1e-9            # shear pair
    assert abs(loc[1] * _IMK2D_L - (q[1] + q[2])) < 1e-6
    # local and global must still describe the same total force
    assert abs(sum(loc[i] ** 2 for i in (0, 1)) ** 0.5
               - sum(glo[i] ** 2 for i in (0, 1)) ** 0.5) < 1e-6


def test_imkbeam3d_tokens():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 3.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("LadrunoIMKBeam", 1, 1, 2, 0.01, 2.0e8, 8.0e7, 2.0e-4,
                1.0e-4, 1.0e-4, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, -20.0, 7.0, 3.0, 2.0, 1.0)
    _static_solve()

    _assert_same(1, "force", "globalForce", "globalForces")
    _assert_same(1, "chordRotation", "basicDeformation", "deformation")
    _assert_same(1, "plasticRotation", "hingeRotation")
    _check_base_chain(1, 12)

    loc = _assert_resolves(1, "localForce")
    q = _resp(1, "basicForce")
    assert len(loc) == 12
    assert abs(loc[0] + q[0]) < 1e-9 and abs(loc[6] - q[0]) < 1e-9   # axial
    assert abs(loc[3] + q[5]) < 1e-9 and abs(loc[9] - q[5]) < 1e-9   # torsion
    assert abs(loc[5] - q[1]) < 1e-9 and abs(loc[11] - q[2]) < 1e-9  # Mz ends
    assert abs(loc[4] - q[3]) < 1e-9 and abs(loc[10] - q[4]) < 1e-9  # My ends


# ---------------------------------------------------------------------------
# tie / coupling family
# ---------------------------------------------------------------------------
def _tie_model(build):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(_HEX, start=1):
        ops.node(i, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3, 2.0)
    ops.element("LadrunoBrick", 100, *range(1, 9), 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    build()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 1.0, 0.5, -2.0)
    _static_solve()


def test_embedded_node_tokens():
    def build():
        ops.node(50, 0.5, 0.5, 0.5)
        ops.element("LadrunoEmbeddedNode", 1, 50, "-host", 100,
                    "-xi", 0.0, 0.0, 0.0, "-k", 1.0e6)
    _tie_model(build)

    _assert_same(1, "penalty", "kt", "k")
    _assert_same(1, "lambda", "augLambda")
    _assert_same(1, "dtCritical", "dtcr")
    _assert_same(1, "massPenalty", "mpenalty")
    for t in ("force", "gap", "penaltyEnergy", "constraintViolation"):
        _assert_resolves(1, t)

    # localForce used to be a second spelling of "force" (global components),
    # which left the D9 local-frame branch permanently dead. It now resolves to
    # the local-frame one.
    _assert_same(1, "localForce", "frameForce", "interfaceForce")
    _assert_same(1, "frameGap", "localGap")


def test_embedded_rebar_tokens():
    def build():
        ops.node(50, 0.5, 0.5, 0.5)
        ops.element("LadrunoEmbeddedRebar", 1, 50, "-host", 100,
                    "-xi", 0.0, 0.0, 0.0, "-dir", 1.0, 0.0, 0.0,
                    "-perfect", 1.0e6, "-kt", 1.0e6)
    _tie_model(build)

    _assert_same(1, "penalty", "kt", "k")     # "k" was NOT accepted before
    _assert_same(1, "lambda", "augLambda")
    _assert_same(1, "dtCritical", "dtcr")
    _assert_same(1, "massPenalty", "mpenalty")
    _assert_same(1, "barAxis", "dir")
    _assert_same(1, "bondEnergy", "bondDissipation")
    for t in ("force", "gap", "slip", "bondStress", "axialForce",
              "penaltyEnergy", "constraintViolation"):
        _assert_resolves(1, t)


def test_kinematic_coupling_tokens():
    def build():
        ops.node(50, 2.0, 0.5, 0.5)
        ops.element("LadrunoKinematicCoupling", 1, 50, 1, 5, "-k", 1.0e7)
    _tie_model(build)

    _assert_same(1, "penalty", "kt", "k")     # "penalty" was NOT accepted before
    _assert_same(1, "lambda", "augLambda")
    _assert_same(1, "dtCritical", "dtcr")
    _assert_same(1, "force", "couplingForce")
    _assert_same(1, "tiedDOFs", "nGap")
    # newly added, to match LadrunoEmbeddedNode/Rebar
    for t in ("penaltyEnergy", "constraintViolation"):
        _assert_resolves(1, t)


def test_distributing_coupling_tokens():
    def build():
        # Two constraints on this element, both load-bearing:
        #  * the reference node MUST carry rotations — the element refuses ndf 3
        #    at setDomain, and a refused element is a hard exit(-1) in FE_Element;
        #  * the independent set must be NON-COLLINEAR. With just nodes 5 and 6
        #    the reference rotation about that line is unconstrained (the element
        #    says so: "degenerate independent set"), the tangent is singular, and
        #    whether the solve reports it is platform-dependent — Linux LAPACK
        #    fails the step, the Windows build did not. 5-6-7 spans a plane.
        ops.node(50, 2.0, 0.5, 0.5, "-ndf", 6)
        ops.element("LadrunoDistributingCoupling", 1, 50, 3, 5, 6, 7,
                    "-k", 1.0e7)
    _tie_model(build)

    _assert_same(1, "penalty", "kt", "k")     # "penalty" was NOT accepted before
    _assert_same(1, "lambda", "augLambda")
    _assert_same(1, "dtCritical", "dtcr")
    _assert_same(1, "massPenalty", "mpenalty")
    _assert_same(1, "force", "couplingForce")
    _assert_same(1, "rotationAxes", "nKept")
    for t in ("penaltyEnergy", "constraintViolation"):
        _assert_resolves(1, t)


# ---------------------------------------------------------------------------
# unregistered tokens must still behave exactly like upstream
# ---------------------------------------------------------------------------
def test_unknown_token_is_still_empty():
    """The fallback is an exact strcmp — a bogus token must NOT resolve to
    something merely because it shares a prefix with a registered one."""
    _build_hex8("LadrunoBrick")
    for bogus in ("stresss", "STRESS", "charLen", "forcess", "notAToken"):
        assert _resp(1, bogus) == [], f"{bogus!r} unexpectedly resolved"
