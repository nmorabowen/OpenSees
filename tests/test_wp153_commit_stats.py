"""WP-153 -- LadrunoSANISAND `commitStats` (response id 33102): the COMMITTED-path
alpha_in census.

`sasStats` counts the alpha_in re-seats of every integration CALL: every Newton
iterate of an implicit step, and both update passes of every
LadrunoDynamicRelaxation step. So its re-seat count is a property of the
solver, not of the path. The TIMs explicit campaign (Step 0) needs the path:
does the DR march COMMIT more re-seats than the implicit march over the same
settlement (guide 06, sec. 6 criterion 3)?

`commitStats` = [commits, commits whose alpha_in differs from the last committed
alpha_in, summed ||alpha_in - alpha_in_n|| over them].

  C-1  ORACLE: on a cyclic plane-strain column (LadrunoQuad, SAS-ME, a Newton
       ladder with step cuts) the count equals the number of commits at which
       the committed alpha_in read through `state` changed, per Gauss point,
       exactly; and the norm sum matches too.
  C-2  it is a PATH count: the sasStats re-seat count is >= it, and strictly
       greater here (the Newton iterates re-seat too).
  C-3  failed attempts do not count (a cut step is reverted, not committed).
  C-4  revertToStart zeroes it.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

TYR = [125.0, 0.05, 0.635, 1.25, 0.712, 0.019, 0.934, 0.7, 100.0, 0.01, 7.05,
       0.968, 1.1, 0.704, 3.5, 4.0, 600.0, 1.62]
R1 = ["-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5]
P0 = 50.0
GPS = [(e, g) for e in (1, 2) for g in (1, 2, 3, 4)]
_converged_steps = [0]


def _build():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for n, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2)], 1):
        ops.node(n, float(x), float(y))
    ops.nDMaterial("LadrunoSANISAND", 1, *TYR, 129, 0, 1, 1.0e-7, 1.0e-4,
                   "-Pmin", 0.0101, "-maxSubsteps", 2000, *R1)
    for e, nn in ((1, (1, 2, 4, 3)), (2, (3, 4, 6, 5))):
        ops.element("LadrunoQuad", e, *nn, 1, "-formulation", "bbar",
                    "-type", "PlaneStrain", "-thick", 1.0)
    ops.fix(1, 1, 1); ops.fix(2, 0, 1); ops.fix(3, 1, 0); ops.fix(5, 1, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, f in ((2, 0.5), (4, 1.0), (6, 0.5)):
        ops.load(n, -P0 * f, 0.0)
    for n in (5, 6):
        ops.load(n, 0.0, -P0 * 0.5)
    ops.constraints("Transformation"); ops.numberer("Plain")
    ops.system("FullGeneral"); ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton"); ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    assert ops.analyze(10) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)
    # vertical strain cycles on the top (nodes 5 6), lateral load held
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    uy = ops.nodeDisp(6, 2)
    for n in (5, 6):
        ops.sp(n, 2, 1.0)
    ops.setTime(uy)


def _alpha_in():
    return [tuple(ops.eleResponse(e, "material", g, "state")[18:24]) for e, g in GPS]


def _commit_stats():
    return [list(ops.eleResponse(e, "material", g, "commitStats")) for e, g in GPS]


def _sas_reseats():
    return [ops.eleResponse(e, "material", g, "sasStats")[11] for e, g in GPS]


def _march():
    """Compress to -0.02, unload to -0.005, reload to -0.025, unload: the
    reversals re-seat alpha_in. A failed Newton step is cut (and reverted)."""
    base = _commit_stats()
    a_prev = _alpha_in()
    changed = [0] * len(GPS)
    dnorm = [0.0] * len(GPS)
    cuts = 0
    _converged_steps[0] = 0
    u = ops.getTime()
    for target in (u - 0.02, u - 0.005, u - 0.025, u - 0.01):
        while abs(ops.getTime() - target) > 1e-12:
            ds = max(-2.0e-4, min(2.0e-4, target - ops.getTime()))
            while True:
                ops.integrator("LoadControl", ds)
                ops.test("NormDispIncr", 1e-10, 30, 0)
                if ops.analyze(1) == 0:
                    _converged_steps[0] += 1
                    break
                cuts += 1
                ds *= 0.5
                assert abs(ds) > 1e-9, "step cut to nothing"
            a = _alpha_in()
            for k in range(len(GPS)):
                d = math.sqrt(sum((x - y) ** 2 for x, y in zip(a[k], a_prev[k])))
                if d > 0.0:
                    changed[k] += 1
                    dnorm[k] += d
            a_prev = a
    return base, changed, dnorm, cuts


def test_C1_C2_C3_commit_stats_is_the_committed_path_count():
    _build()
    base, changed, dnorm, cuts = _march()
    cs = _commit_stats()
    for k, (e, g) in enumerate(GPS):
        n_re = cs[k][1] - base[k][1]
        assert n_re == changed[k], (e, g, n_re, changed[k])                       # C-1
        assert cs[k][2] - base[k][2] == pytest.approx(dnorm[k], rel=1e-12, abs=1e-15)
    tot_path = sum(changed)
    assert tot_path > 0, "the cycle re-seated nothing: the test is vacuous"
    assert sum(_sas_reseats()) > tot_path                                          # C-2
    # C-3: a cut attempt is reverted, never committed, so every GP's commit count
    # since the march began is the number of CONVERGED march steps, cuts or not
    assert len({c[0] for c in cs}) == 1
    assert cs[0][0] - base[0][0] == _converged_steps[0]


def test_C4_revertToStart_zeroes_it():
    _build()
    _march()
    assert sum(c[1] for c in _commit_stats()) > 0
    ops.reset()
    assert all(c == [0.0, 0.0, 0.0] for c in _commit_stats())
