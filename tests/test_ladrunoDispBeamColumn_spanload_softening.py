"""LadrunoDispBeamColumn -spanLoad: the bubble solve on cracked and softening sections.

A cracked section with one bar layer yielded leaves Kbb nearly rank one, and softening
fibers make it indefinite. The old inner Newton diverged there, stopped at its cap, and
returned success with an unbalanced bubble (residual up to 1e6 N against a test of 1e-3):
the element force was then not a function of the deformation, the global Newton cycled,
and the next re-solve from the new committed beta moved the force after a converged step.
The 3-story, 2-bay RC frame of the validation suite (04a) stopped at 80/141 and 21/141
steps at 4 and 8 elements per member; the same frame with -noSpanLoad completes.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_BAYS = (0.0, 6000.0, 12000.0)
_FLOORS = (0.0, 4000.0, 7500.0, 11000.0)
_FC, _EC = 35.0, 27805.57
_ROOF = 10
_GJ = 1.0e13


def _section(tag, b, h, n_bar_layer, as_layer, ndm):
    cov = 40.0
    extra = ("-GJ", _GJ) if ndm == 3 else ()
    ops.section("Fiber", tag, *extra)
    ops.patch("rect", 1, 10, 1, -(0.5 * h - cov), -0.5 * b + cov, 0.5 * h - cov, 0.5 * b - cov)
    ops.patch("rect", 2, 2, 1, -0.5 * h, -0.5 * b, -(0.5 * h - cov), 0.5 * b)
    ops.patch("rect", 2, 2, 1, 0.5 * h - cov, -0.5 * b, 0.5 * h, 0.5 * b)
    y = 0.5 * h - cov
    for s in (-y, y):
        ops.layer("straight", 3, n_bar_layer, as_layer / n_bar_layer, s, -0.5 * b + cov, s, 0.5 * b - cov)


def _frame(ndm, n_sub, span_flag):
    """Gravity-loaded frame, ready for the lateral push. span_flag: -spanLoad / -noSpanLoad."""
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", 3 if ndm == 2 else 6)

    def node(tag, x, y):
        if ndm == 2:
            ops.node(tag, x, y)
        else:
            ops.node(tag, x, y, 0.0)
            ops.fix(tag, 0, 0, 1, 1, 1, 0)

    for f, y in enumerate(_FLOORS):
        for b, x in enumerate(_BAYS):
            node(3 * f + b + 1, x, y)
    for n in (1, 2, 3):
        ops.fix(n, *((1, 1, 1) if ndm == 2 else (1, 1, 0, 0, 0, 1)))

    vec = () if ndm == 2 else (0.0, 0.0, 1.0)
    ops.geomTransf("PDelta", 1, *vec)
    ops.geomTransf("Linear", 2, *vec)

    eps0 = -2.0 * _FC / _EC
    ops.uniaxialMaterial("Concrete02", 1, -1.2 * _FC, 1.2 * eps0, -0.5 * _FC, -0.025, 0.10, 3.0, 1500.0)
    ops.uniaxialMaterial("Concrete02", 2, -_FC, eps0, -0.35 * _FC, -0.018, 0.10, 2.5, 1500.0)
    ops.uniaxialMaterial("Steel02", 3, 420.0, 200000.0, 0.012, 18.0, 0.925, 0.15)
    _section(1, 450.0, 450.0, 4, 1600.0, ndm)
    _section(2, 350.0, 550.0, 3, 1080.0, ndm)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.beamIntegration("Lobatto", 2, 2, 5)

    members = []
    for f in range(3):
        for b in range(3):
            members.append((3 * f + b + 1, 3 * (f + 1) + b + 1, 1, 1))
    for f in range(1, 4):
        for b in range(2):
            members.append((3 * f + b + 1, 3 * f + b + 2, 2, 2))

    beams = []
    next_node, next_ele = 1000, 1
    for ni, nj, trf, integ in members:
        ci, cj = ops.nodeCoord(ni), ops.nodeCoord(nj)
        chain = [ni]
        for k in range(1, n_sub):
            t = k / n_sub
            node(next_node, ci[0] + t * (cj[0] - ci[0]), ci[1] + t * (cj[1] - ci[1]))
            chain.append(next_node)
            next_node += 1
        chain.append(nj)
        for na, nb in zip(chain[:-1], chain[1:]):
            ops.element("LadrunoDispBeamColumn", next_ele, na, nb, trf, integ, span_flag, "-lch", "ip")
            if integ == 2:
                beams.append(next_ele)
            next_ele += 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in beams:
        ops.eleLoad("-ele", tag, "-type", "-beamUniform", -25.0, *((0.0,) if ndm == 2 else (0.0, 0.0)))
    for f in range(1, 4):
        for b, p in enumerate((-100.0e3, -200.0e3, -100.0e3)):
            load = (0.0, p, 0.0) if ndm == 2 else (0.0, p, 0.0, 0.0, 0.0, 0.0)
            ops.load(3 * f + b + 1, *load)

    ops.system("BandGeneral" if ndm == 2 else "UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.2)
    ops.analysis("Static")
    assert ops.analyze(5) == 0
    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for f, w in enumerate((1.0, 2.0, 3.0), start=1):
        for b in range(3):
            load = (w, 0.0, 0.0) if ndm == 2 else (w, 0.0, 0.0, 0.0, 0.0, 0.0)
            ops.load(3 * f + b + 1, *load)
    return list(ops.getEleTags())


def _step(du):
    """Newton, KrylovNewton, then 10 sub-steps: the ladder of the validation frame."""
    ops.integrator("DisplacementControl", _ROOF, 1, du)
    for alg in (("Newton",), ("KrylovNewton",)):
        ops.algorithm(*alg)
        if ops.analyze(1) == 0:
            ops.algorithm("Newton")
            return 0
    ops.integrator("DisplacementControl", _ROOF, 1, du / 10.0)
    for _ in range(10):
        for alg in (("Newton",), ("KrylovNewton",)):
            ops.algorithm(*alg)
            rc = ops.analyze(1)
            if rc == 0:
                break
        if rc != 0:
            return rc
    ops.algorithm("Newton")
    return 0


def _base_shear():
    ops.reactions()
    return -sum(ops.nodeReaction(n, 1) for n in (1, 2, 3)) / 1.0e3


def _push(targets):
    """Returns (steps completed, base shears)."""
    u, shear = 0.0, []
    for k, tgt in enumerate(targets):
        if _step(tgt - u) != 0:
            return k, shear
        u = ops.nodeDisp(_ROOF, 1)
        shear.append(_base_shear())
    return len(targets), shear


_MONO = [385.0 * (k + 1) / 140 for k in range(140)]


@pytest.mark.parametrize("ndm,n_sub", [(2, 4), (2, 8), (3, 4), (3, 8)])
def test_spanload_push_completes(ndm, n_sub):
    _frame(ndm, n_sub, "-spanLoad")
    done, _ = _push(_MONO)
    assert done == len(_MONO)


def test_spanload_matches_nospanload_on_a_fine_mesh():
    # Eight elements per member already resolve the span load: the flag must not move the push.
    _frame(2, 8, "-spanLoad")
    done_a, va = _push(_MONO)
    _frame(2, 8, "-noSpanLoad")
    done_b, vb = _push(_MONO)
    assert done_a == done_b == len(_MONO)
    assert max(abs(v) for v in va) == pytest.approx(max(abs(v) for v in vb), rel=0.01)
    assert va[-1] == pytest.approx(vb[-1], rel=0.02)


def test_spanload_cyclic_completes():
    _frame(2, 8, "-spanLoad")
    pts = [0.0]
    for p in (38.5, -38.5, 82.5, -82.5, 165.0, -165.0, 275.0, -275.0):
        start = pts[-1]
        pts += [start + (p - start) * (k + 1) / 25 for k in range(25)]
    done, _ = _push(pts[1:])
    assert done == len(pts) - 1


def test_force_does_not_move_after_a_converged_step():
    # Forming reactions re-applies the loads and re-solves the bubble from the new committed
    # beta. With a solved bubble that finds the same state; the old code moved the force.
    tags = _frame(2, 4, "-spanLoad")
    u = 0.0
    worst = 0.0
    for tgt in _MONO[:80]:
        assert _step(tgt - u) == 0
        u = ops.nodeDisp(_ROOF, 1)
        before = {t: list(ops.eleForce(t)) for t in tags}
        ops.reactions()
        for t in tags:
            for a, b in zip(before[t], ops.eleForce(t)):
                worst = max(worst, abs(a - b) / max(abs(a), 1.0e6))
    assert worst < 1.0e-6
