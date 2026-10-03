"""WP-153 -- LadrunoDynamicRelaxation::revertToLastStep: a failed step retried
must march EXACTLY as if it had never failed.

Found by the TIMs explicit/DR footing campaign (Step 0, 2026-09-29): DR did not
override IncrementalIntegrator::revertToLastStep(), so it inherited a no-op. On a
failed step DirectIntegrationAnalysis reverts the Domain and calls it, but DR's
private leap-frog state (Ut, Vhalf, Aprev -- and M*, which -recompute and the
KE-peak auto-refresh rebuild INSIDE a step) stayed where the failed step left
it, so the retry marched from a state that was never committed. The fix
snapshots the march state at every successful commit (and at domainChanged) and
restores it.

The failure is injected through DR's own one-solve guard: under `algorithm
Newton` with an unreachable tolerance the second iteration calls update() a
second time, which DR refuses. By then newStep() AND one update() have moved
every private vector -- the stale case. The algorithm is then set back to
`Linear` (DirectIntegrationAnalysis::setAlgorithm touches the algorithm only,
never the integrator) and the march resumes.

  R-1  failed-then-retried == uninterrupted, BIT FOR BIT, on a plastic brick
       block under kinetic (Cundall) damping, with -recompute landing on the
       failed step, and under viscous damping.
  R-2  the failure is REAL: the failed attempt returns < 0 and leaves the
       Domain at the last commit (time and displacement).
  R-3  two consecutive failures, then the retry: still identical.

Non-vacuity was measured against the pre-fix build (8ebde5cbd, Esmeralda,
2026-09-30): R-1 fails there on every parametrization (see the PR body).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_E, _NU, _SIG0 = 2.0e5, 0.3, 250.0
_K, _G = _E / (3.0 * (1 - 2 * _NU)), _E / (2.0 * (1 + _NU))
_RATE = 1.0e-5                 # m per DR step on the top face: 4e-3 at step 400 (3x yield)
_N = 400                       # DR steps in the march
_FAIL_AT = 150                 # the failed attempt is DR step 150


def _build():
    """A 2x2x2 LadrunoBrick unit cube of LadrunoJ2, base fixed, the top face
    driven down in z by a Linear-series SP: u_z = -t (the DR clock)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, _K, _G, "-iso", "voce", _SIG0, 0.0, 0.0, 200.0)
    n, h = 3, 0.5
    tag = {}
    t = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                ops.node(t, i * h, j * h, k * h)
                tag[(i, j, k)] = t
                t += 1
    e = 1
    for k in range(n - 1):
        for j in range(n - 1):
            for i in range(n - 1):
                c = [tag[(i, j, k)], tag[(i + 1, j, k)], tag[(i + 1, j + 1, k)],
                     tag[(i, j + 1, k)], tag[(i, j, k + 1)], tag[(i + 1, j, k + 1)],
                     tag[(i + 1, j + 1, k + 1)], tag[(i, j + 1, k + 1)]]
                ops.element("LadrunoBrick", e, *c, 1)
                e += 1
    for (i, j, k), nd in tag.items():
        if k == 0:
            ops.fix(nd, 1, 1, 1)
    top = [nd for (i, j, k), nd in tag.items() if k == n - 1]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nd in top:
        ops.sp(nd, 3, -1.0)
    return sorted(tag.values())


def _dr(options):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation", "-dt", _RATE, *options)
    ops.analysis("Transient")


def _state(nodes):
    return [ops.getTime()] + [u for nd in nodes for u in ops.nodeDisp(nd)]


def _fail_once():
    """One DR step that fails AFTER newStep() and one update() moved the march."""
    ops.test("NormUnbalance", 1.0e-300, 2, 0)
    ops.algorithm("Newton")
    rc = ops.analyze(1, _RATE)
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    return rc


_OPTS = {
    "kinetic": ["-damping", "kinetic"],
    "kinetic_recompute": ["-damping", "kinetic", "-recompute", 50],
    "viscous": ["-damping", "viscous", 1.0],
}


def _uninterrupted(opts):
    nodes = _build()
    _dr(opts)
    assert ops.analyze(_N, _RATE) == 0
    return _state(nodes)


@pytest.mark.parametrize("name", sorted(_OPTS))
def test_R1_failed_then_retried_equals_uninterrupted(name):
    opts = _OPTS[name]
    ref = _uninterrupted(opts)
    nodes = _build()
    _dr(opts)
    assert ops.analyze(_FAIL_AT - 1, _RATE) == 0
    assert _fail_once() < 0
    assert ops.analyze(_N - (_FAIL_AT - 1), _RATE) == 0
    got = _state(nodes)
    assert got == ref, (name, max(abs(a - b) for a, b in zip(got, ref)))


def test_R2_the_failure_is_real_and_the_domain_sits_at_the_last_commit():
    nodes = _build()
    _dr(_OPTS["kinetic"])
    assert ops.analyze(_FAIL_AT - 1, _RATE) == 0
    before = _state(nodes)
    assert before[0] > 0.0
    assert _fail_once() < 0
    assert _state(nodes) == before


def test_R3_two_failures_in_a_row_then_the_retry():
    ref = _uninterrupted(_OPTS["kinetic_recompute"])
    nodes = _build()
    _dr(_OPTS["kinetic_recompute"])
    assert ops.analyze(_FAIL_AT - 1, _RATE) == 0
    assert _fail_once() < 0
    assert _fail_once() < 0
    assert ops.analyze(_N - (_FAIL_AT - 1), _RATE) == 0
    assert _state(nodes) == ref
