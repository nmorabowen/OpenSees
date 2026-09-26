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
    cache would re-form Ki at the new E and converge at once.
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


def _mids(verts, edges):
    return [tuple(0.5 * (verts[a][d] + verts[b][d]) for d in range(len(verts[0])))
            for a, b in edges]


def _plane(coords, fixes, finite=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, (x, y) in enumerate(coords):
        ops.node(i + 1, x, y)
    for n, fx in fixes:
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
    for n in fixed:
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
