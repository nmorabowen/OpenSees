"""WP-124: the shared Element-contract helpers of the eight continuum element shells
(SRC/element/LadrunoElementShell.h), pinned on every element that uses them.

Each test has an ABSOLUTE oracle and was run against a one-line mutation of the helper it
guards (Ladruno_implementation/wp124_shells/mutation_rows.py); the table is in
Ladruno_implementation/124_continuum_shell_helpers.md.

Ki cache (cacheKi / dropKi)
  * on a linear-elastic model K == K0, so ModifiedNewton -initial must converge in <= 2
    iterations and land on the Newton answer: a zero, scaled or stale Ki cannot;
    (LadrunoCSTPair: its Ki is the UNRELIEVED per-triangle B'D0B while its tangent is the
    F-bar-patch operator, so Ki != K even at F = I by design -- there the check is only
    "converges to the Newton answer", which a zero Ki still fails);
  * Ki is formed ONCE (the vanilla convention): after halving E through a parameter,
    ModifiedNewton -initial must iterate against the frozen K0 -- a helper that forgot to
    cache would re-form Ki at the new E and converge at once;
  * C6: a restore into a LIVE element (database checkpoint) drops Ki (dropKi in recvSelf).
"""
import pytest

import opensees as ops

E, NU, RHO = 1000.0, 0.3, 2.0

QUAD = [(0.0, 0.0), (2.0, 0.0), (2.2, 1.1), (-0.1, 1.0)]
TRI = [(0.0, 0.0), (2.0, 0.0), (0.3, 1.2)]
TRI_MIDS = [(0, 1), (1, 2), (2, 0)]
HEX = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
       (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.1, 1.05, 1.1), (0.0, 1.0, 1.0)]
HEX20_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
               (0, 4), (1, 5), (2, 6), (3, 7)]
TET = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
TET_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]


_SUPPORTS = [True]          # False: the builders skip their fix() calls (rigid-body probes)


def _mids(verts, edges):
    return [tuple(0.5 * (verts[a][d] + verts[b][d]) for d in range(len(verts[0])))
            for a, b in edges]


def _plane(coords, fixes, finite=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, (x, y) in enumerate(coords):
        ops.node(i + 1, x, y)
    for n, fx in (fixes if _SUPPORTS[0] else []):
        ops.fix(n, *fx)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    if finite:
        ops.nDMaterial("LogStrain2D", 2, 1)
        return 2
    return 1


def _solid(coords, fixed):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, xyz in enumerate(coords):
        ops.node(i + 1, *xyz)
    for n in (fixed if _SUPPORTS[0] else []):
        ops.fix(n, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    return 1


# name -> builder returning (loaded node, load vector)
def _quad():
    m = _plane(QUAD, [(1, (1, 1)), (2, (0, 1))])
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, m, "-type", "PlaneStrain", "-thick", 0.8)
    return 3, (1.0, 0.3)


def _cst():
    m = _plane(TRI, [(1, (1, 1)), (2, (0, 1))])
    ops.element("LadrunoCST", 1, 1, 2, 3, m, "-type", "PlaneStrain", "-thick", 0.8)
    return 3, (1.0, 0.3)


def _lst():
    m = _plane(TRI + _mids(TRI, TRI_MIDS), [(1, (1, 1)), (2, (0, 1)), (4, (0, 1))])
    ops.element("LadrunoLST", 1, *range(1, 7), m, "-type", "PlaneStrain", "-thick", 0.8)
    return 3, (1.0, 0.3)


def _cstpair():
    m = _plane(QUAD, [(1, (1, 1)), (2, (0, 1))], finite=True)
    ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, m, "-thick", 0.8)
    return 3, (1.0e-3, 3.0e-4)          # finite-only element: stay in the linear range


def _brick():
    m = _solid(HEX, [1, 2, 3, 4])
    ops.element("LadrunoBrick", 1, *range(1, 9), m)
    return 7, (1.0, 0.3, 0.0)


def _brick20():
    reg = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
           (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
    m = _solid(reg + _mids(reg, HEX20_EDGES), [1, 2, 3, 4, 9, 10, 11, 12])
    ops.element("LadrunoBrick20", 1, *range(1, 21), m)
    return 7, (1.0, 0.3, 0.0)


def _tri6():
    m = _plane(TRI + _mids(TRI, TRI_MIDS), [(1, (1, 1)), (2, (0, 1)), (4, (0, 1))])
    ops.element("BezierTri6", 1, *range(1, 7), 0.8, "PlaneStrain", m, "-rho", RHO)
    return 3, (1.0, 0.3)


def _tet10():
    m = _solid(TET + _mids(TET, TET_EDGES), [1, 2, 3, 5, 6, 7])
    ops.element("BezierTet10", 1, *range(1, 11), m, "-rho", RHO)
    return 4, (1.0, 0.3, 0.0)


ELEMENTS = {"LadrunoQuad": _quad, "LadrunoCST": _cst, "LadrunoLST": _lst,
            "LadrunoCSTPair": _cstpair, "LadrunoBrick": _brick, "LadrunoBrick20": _brick20,
            "BezierTri6": _tri6, "BezierTet10": _tet10}


def _static(node, load, algo, tol=1.0e-12, iters=2, lf=1.0):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(node, *load)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm(*algo)
    ops.integrator("LoadControl", lf)
    ops.analysis("Static", "-noWarnings")


def _all_disp(n):
    return [x for i in range(1, n + 1) for x in ops.nodeDisp(i)]


@pytest.mark.parametrize("name", list(ELEMENTS))
def test_initial_stiffness_is_K0(name):
    """Elastic: ModifiedNewton -initial converges like Newton (Ki == K0 == K)."""
    node, load = ELEMENTS[name]()
    nn = len(ops.getNodeTags())
    _static(node, load, ("Newton",), iters=10)
    assert ops.analyze(1) == 0
    ref = _all_disp(nn)

    ELEMENTS[name]()
    iters = 200 if name == "LadrunoCSTPair" else 3
    _static(node, load, ("ModifiedNewton", "-initial"), iters=iters)
    assert ops.analyze(1) == 0, f"{name}: ModifiedNewton -initial did not converge in {iters}"
    got = _all_disp(nn)
    scale = max(abs(x) for x in ref)
    assert scale > 0.0
    rel = 1.0e-5 if name == "LadrunoCSTPair" else 1.0e-9   # pair: linear convergence
    assert max(abs(a - b) for a, b in zip(got, ref)) <= rel * scale


# not CSTPair: Ki != K by design there, so iteration counts cannot see a re-formed Ki
FREEZE = [n for n in ELEMENTS if n != "LadrunoCSTPair"]


@pytest.mark.parametrize("name", FREEZE)
def test_initial_stiffness_is_formed_once(name):
    """Ki is cached at first use: halving E afterwards leaves ModifiedNewton -initial
    iterating against the OLD K0 (error ratio ~1/2 per iteration)."""
    node, load = ELEMENTS[name]()
    ops.parameter(1, "element", 1, "E")
    _static(node, load, ("ModifiedNewton", "-initial"), tol=1.0e-12, iters=200, lf=0.5)
    assert ops.analyze(1) == 0                  # forms and caches Ki at E
    ops.updateParameter(1, 0.5 * E)
    assert ops.analyze(1) == 0
    n_iter = ops.testIter()
    assert n_iter > 8, f"{name}: Ki was re-formed after the E update ({n_iter} iterations)"


# C6: restore into a LIVE element (same domain, database restore of a checkpoint).
# Not CSTPair (Ki != K by design, so iteration counts cannot see a stale Ki).
# LadrunoBrick20 also needs C14: before it, the live restore left its geometry cache empty
# (recvSelf cleared it; no setDomain follows) and the system was singular.
LIVE = [n for n in ELEMENTS if n != "LadrunoCSTPair"]


@pytest.mark.parametrize("name", LIVE)
def test_live_restore_drops_Ki(name, tmp_path):
    """Save at E, move the live elements to E/2 and form Ki there, restore the checkpoint
    (live recvSelf back to E). Ki must follow: ModifiedNewton -initial converges at once.
    A kept Ki (E/2) gives the iteration factor 1 - K/K0 = -1 and never converges."""
    node, load = ELEMENTS[name]()
    ops.parameter(1, "element", 1, "E")
    ops.database("File", str(tmp_path / "db"))
    ops.save(1)
    ops.updateParameter(1, 0.5 * E)
    # form Ki at E/2 WITHOUT adding a load pattern: an object the checkpoint does not hold
    # makes the restore itself fail
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 5, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0
    assert ops.restore(1) in (0, None)
    ops.wipeAnalysis()
    ops.timeSeries("Linear", 7)
    ops.pattern("Plain", 7, 7)
    ops.load(node, *load)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    assert ops.analyze(1) == 0, f"{name}: stale Ki after a live restore"
    assert ops.testIter() <= 3


@pytest.mark.parametrize("name", list(ELEMENTS))
def test_live_restore_is_usable(name, tmp_path):
    """C14 oracle, independent of Ki: after a checkpoint restore into the LIVE domain, Newton
    must land on the answer of a freshly built model at the restored E (LadrunoBrick20 was
    singular: its geometry cache stayed empty because no setDomain follows a live restore)."""
    node, load = ELEMENTS[name]()
    _static(node, load, ("Newton",), iters=20)
    assert ops.analyze(1) == 0
    fresh = _all_disp(len(ops.getNodeTags()))

    ELEMENTS[name]()
    ops.database("File", str(tmp_path / "db"))
    ops.save(1)
    assert ops.restore(1) in (0, None)
    _static(node, load, ("Newton",), iters=20)
    assert ops.analyze(1) == 0, f"{name}: analysis fails after a live restore"
    got = _all_disp(len(ops.getNodeTags()))
    scale = max(abs(x) for x in fresh)
    assert max(abs(a - b) for a, b in zip(got, fresh)) <= 1.0e-9 * scale


# ---- ground-motion inertia (addGroundInertia) -------------------------------------------
# Rigid-body probe (ladruno-new-element guide): an UNSUPPORTED element under a constant
# UniformExcitation a_g moves rigidly, so every node's RELATIVE acceleration is exactly
# -a_g. A sign flip gives +a_g, a skipped load gives 0, and a full (consistent) mass reduced
# to its diagonal breaks the nodal balance. Horizontal and vertical ground motion.
def _tri6_cmass():
    m = _plane(TRI + _mids(TRI, TRI_MIDS), [(1, (1, 1)), (2, (0, 1)), (4, (0, 1))])
    ops.element("BezierTri6", 1, *range(1, 7), 0.8, "PlaneStrain", m, "-rho", RHO, "-cMass")
    return 3, (1.0, 0.3)


def _tet10_cmass():
    m = _solid(TET + _mids(TET, TET_EDGES), [1, 2, 3, 5, 6, 7])
    ops.element("BezierTet10", 1, *range(1, 11), m, "-rho", RHO, "-cMass")
    return 4, (1.0, 0.3, 0.0)


# Bezier defaults to LUMPED mass; the -cMass variants put a genuinely full M through the
# helper's full branch (Brick/Brick20 are consistent by default)
GROUND = dict(ELEMENTS, **{"BezierTri6-cMass": _tri6_cmass, "BezierTet10-cMass": _tet10_cmass})


@pytest.mark.parametrize("direction", [1, 2])
@pytest.mark.parametrize("name", list(GROUND))
def test_ground_inertia_rigid_body(name, direction):
    _SUPPORTS[0] = False
    try:
        GROUND[name]()
    finally:
        _SUPPORTS[0] = True
    ndf = len(ops.nodeDisp(ops.getNodeTags()[0]))
    ag = 1.5
    ops.timeSeries("Constant", 3, "-factor", ag)
    ops.pattern("UniformExcitation", 3, direction, "-accel", 3)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-14, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    for _ in range(3):
        assert ops.analyze(1, 0.01) == 0
    for n in ops.getNodeTags():
        acc = ops.nodeAccel(n)
        for d in range(ndf):
            want = -ag if d == direction - 1 else 0.0
            assert abs(acc[d] - want) <= 1.0e-9 * ag,                 f"{name} node {n} dof {d + 1}: relative accel {acc[d]!r}, want {want!r}"


# ---- parameter forwarding (forwardToMaterials / forwardToMaterialPoint) -------------------
# Linear elastic: halving E everywhere doubles every displacement EXACTLY (ratio 2), whether
# through the forall broadcast ("E") or one "material k E" parameter per Gauss point.
N_POINTS = {"LadrunoQuad": 4, "LadrunoCST": 1, "LadrunoLST": 3, "LadrunoCSTPair": 2,
            "LadrunoBrick": 8, "LadrunoBrick20": 27, "BezierTri6": 3, "BezierTet10": 4}
PARAM = list(ELEMENTS)   # C3 fixed: LadrunoCSTPair forwards parameters too


def _solve_disp(name, params=()):
    node, load = ELEMENTS[name]()
    for tag, argv in enumerate(params, start=1):
        ops.parameter(tag, "element", 1, *argv)
        ops.updateParameter(tag, 0.5 * E)
    _static(node, load, ("Newton",), iters=20)
    assert ops.analyze(1) == 0
    return _all_disp(len(ops.getNodeTags()))


def _ratio(ref, got):
    k = max(range(len(ref)), key=lambda i: abs(ref[i]))
    return got[k] / ref[k]


@pytest.mark.parametrize("name", PARAM)
def test_forall_parameter_reaches_every_material(name):
    ref = _solve_disp(name)
    r = _ratio(ref, _solve_disp(name, [("E",)]))
    tol = 1.0e-5 if name == "LadrunoCSTPair" else 1.0e-9   # finite-only element
    assert abs(r - 2.0) <= tol, f"{name}: forall E/2 gives displacement ratio {r!r}, want 2"


@pytest.mark.parametrize("name", PARAM)
def test_material_point_parameters_cover_every_point(name):
    ref = _solve_disp(name)
    n = N_POINTS[name]
    r = _ratio(ref, _solve_disp(name, [("material", str(k), "E") for k in range(1, n + 1)]))
    tol = 1.0e-5 if name == "LadrunoCSTPair" else 1.0e-9
    assert abs(r - 2.0) <= tol, f"{name}: 'material k E' for k=1..{n} gives ratio {r!r}, want 2"


# ---- response finalise (finishResponse + Element::getResponse fallback) -----------------
RESP = list(ELEMENTS)   # C2 fixed: Bezier chains to Element::setResponse too


def _hex_volume(X):
    """Exact volume of a trilinear hex (detJ is at most quadratic per axis: 2-pt Gauss)."""
    g = 1.0 / 3.0 ** 0.5
    s = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
         (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
    vol = 0.0
    for a in (-g, g):
        for b in (-g, g):
            for c in (-g, g):
                J = [[0.0] * 3 for _ in range(3)]
                for (si, ti, ui), x in zip(s, X):
                    d = (si * (1 + ti * b) * (1 + ui * c) / 8, ti * (1 + si * a) * (1 + ui * c) / 8,
                         ui * (1 + si * a) * (1 + ti * b) / 8)
                    for r in range(3):
                        for q in range(3):
                            J[r][q] += x[r] * d[q]
                vol += (J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                        - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                        + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]))
    return vol


def _shoelace(P):
    return 0.5 * abs(sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1]
                         for i in range(len(P))))


VOLUME = {"LadrunoQuad": _shoelace(QUAD) * 0.8, "LadrunoCST": _shoelace(TRI) * 0.8,
          "LadrunoLST": _shoelace(TRI) * 0.8, "LadrunoCSTPair": _shoelace(QUAD) * 0.8,
          "BezierTri6": _shoelace(TRI) * 0.8, "LadrunoBrick": _hex_volume(HEX),
          "LadrunoBrick20": 1.0, "BezierTet10": 1.0 / 6.0}


@pytest.mark.parametrize("name", RESP)
def test_base_response_vocabulary(name):
    """globalForce == force; on an unsupported element translating at v0 with alphaM only,
    the x-components of dampingForce sum to alphaM * rho * V * v0 (row sums of M: exact for
    lumped and consistent mass alike)."""
    _SUPPORTS[0] = False
    try:
        ELEMENTS[name]()
    finally:
        _SUPPORTS[0] = True
    tags = ops.getNodeTags()
    ndf = len(ops.nodeDisp(tags[0]))
    for i, n in enumerate(tags):                   # a non-trivial deformation for 'force'
        ops.setNodeDisp(n, 1, 1.0e-3 * (i % 3), "-commit")
    force = ops.eleResponse(1, "force")
    glob = ops.eleResponse(1, "globalForce")
    assert glob and len(glob) == len(force), f"{name}: globalForce records nothing"
    assert glob == force
    for n in tags:
        ops.setNodeDisp(n, 1, 0.0, "-commit")
    aM, v0 = 0.3, 2.0
    ops.rayleigh(aM, 0.0, 0.0, 0.0)
    for n in tags:
        ops.setNodeVel(n, 1, v0, "-commit")
    damp = ops.eleResponse(1, "dampingForce")
    assert damp, f"{name}: dampingForce records nothing"
    fx = sum(damp[i] for i in range(0, len(damp), ndf))
    want = aM * RHO * VOLUME[name] * v0
    assert abs(fx - want) <= 1.0e-9 * want, f"{name}: sum dampingForce_x {fx!r}, want {want!r}"

# ---- gap fixes C4, C5, C10 --------------------------------------------------------------
def test_cstpair_stress_plane_strain():
    """C4: the fourth plane element exposes the 4-component plane-strain stress per triangle
    ([sxx, syy, sxy, szz] x 2), in-plane parts identical to 'stress'."""
    node, load = _cstpair()
    _static(node, load, ("Newton",), iters=20)
    assert ops.analyze(1) == 0
    s = ops.eleResponse(1, "stress")
    s4 = ops.eleResponse(1, "stressPlaneStrain")
    assert len(s) == 6 and len(s4) == 8, (len(s), len(s4))
    for t in range(2):
        assert s4[4 * t:4 * t + 3] == s[3 * t:3 * t + 3]


@pytest.mark.parametrize("k", [2, 3, 4])
def test_quad_ssp_material_point_maps_to_the_live_slot(k):
    """C5: under -formulation ssp the only live material is slot 0 (setResponse and
    LadrunoBrick already map every k there), so 'material k E' must act like 'material 1 E'."""
    def solve(param):
        _plane(QUAD, [(1, (1, 1)), (2, (0, 1))])
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-formulation", "ssp", "-type", "PlaneStrain")
        if param:
            ops.parameter(1, "element", 1, *param)
            ops.updateParameter(1, 0.5 * E)
        _static(3, (1.0, 0.3), ("Newton",), iters=20)
        assert ops.analyze(1) == 0
        return ops.nodeDisp(3, 1)
    ref = solve(None)
    one = solve(("material", "1", "E")) / ref
    kth = solve(("material", str(k), "E")) / ref
    assert one != 1.0
    assert kth == one, f"'material {k} E' ratio {kth!r} != 'material 1 E' ratio {one!r}"


C10 = pytest.mark.xfail(strict=True, reason="C10: 'materialState' is swallowed by the 'material' branch")
STATE = [n if n.startswith("Bezier") else pytest.param(n, marks=C10)
         for n in ELEMENTS if n != "LadrunoCSTPair"]


@pytest.mark.parametrize("name", STATE)
def test_material_state_reaches_the_materials(name, capfd, monkeypatch):
    """C10: the UW staged-analysis switch ('parameter ... materialState', e.g. DruckerPrager
    elastic -> plastic) must reach the GP materials. The Ladruno elements' 'material'-prefix
    test swallowed it (argc < 3 -> -1); Bezier and SixNodeTri exclude it."""
    orig = ops.nDMaterial

    def dp(*a):
        if a[0] == "ElasticIsotropic":
            return orig("DruckerPrager", a[1], 833.3, 384.6, 1.0e3,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, RHO)
        return orig(*a)
    monkeypatch.setattr(ops, "nDMaterial", dp)
    ELEMENTS[name]()
    capfd.readouterr()
    ops.parameter(1, "element", 1, "materialState")
    err = capfd.readouterr().err
    assert "no objects were able to identify" not in err, f"{name}: materialState unclaimed"
