"""WP-240 -- LadrunoDynamicRelaxation updates the elements ONCE per step.

Ported from jaabell/OpenSees c8b0f3bc4 (José Abell). `update()` used to call
`updateDomain()` after `setResponse()`. The displacements it pushed were those
of `newStep()`, whose `updateDomain()` had already brought the elements to
them, so every element that integrates its materials in `update()` did each
material point twice. The second pass is gone.

  S-1  WITNESS. A truss with a rate-dependent material records the velocity
       its last `update()` saw. After a DR step that is the half-step velocity
       of `newStep()`, (u_{n+1} - u_n)/dt, and not the full-step velocity that
       `update()` pushes to the node afterwards. Fails on a build with the
       second element pass (it then sees the full-step velocity).
  S-2  SAME ANSWER. Rate-independent models give the numbers the two-pass
       build gave (references recorded with the build of ladruno f9864cfec,
       before the change): a plastic LadrunoBrick block driven by prescribed
       `sp` displacements under the Transformation handler (kinetic and
       viscous damping), and a plastic LadrunoQuad panel under nodal load.
  S-3  The same block with a failed step in the middle (WP-153
       `revertToLastStep`) lands on the same reference.

What changes, by design: an element whose `update()` reads the nodal velocity
commits the state of the half-step velocity, the one its force in the residual
was computed with (S-1). Quirk: ledger/quirks/WP-240-*.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# ---------------------------------------------------------------------------
# S-1  one element update per step, seen through a rate-dependent material
# ---------------------------------------------------------------------------
_E1, _C1, _P1, _DT1 = 1000.0, 5.0, 10.0, 1.0


def _viscous_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, _E1)
    ops.uniaxialMaterial("Viscous", 2, _C1, 1.0)
    ops.uniaxialMaterial("Parallel", 3, 1, 2)
    ops.element("Truss", 1, 1, 2, 1.0, 3)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, _P1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e30, 1, 0)
    ops.algorithm("Linear")
    # viscous damping: the velocities are never zeroed, so every step has a
    # half-step and a full-step velocity that differ
    ops.integrator("LadrunoDynamicRelaxation", "-damping", "viscous", 1.0)
    ops.analysis("Transient")


def test_S1_the_element_sees_the_half_step_velocity_once():
    _viscous_truss()
    u_prev = 0.0
    distinct = 0
    for step in range(1, 21):
        assert ops.analyze(1, _DT1) == 0
        u = ops.nodeDisp(2, 1)
        force = ops.eleResponse(1, "axialForce")[0]
        v_seen = (force - _E1 * u) / _C1          # the rate of the last update()
        v_half = (u - u_prev) / _DT1               # newStep(): u += dt v_{n+1/2}
        v_full = ops.nodeVel(2, 1)                 # pushed by update() afterwards
        scale = max(abs(v_half), abs(v_full), 1.0e-12)
        assert abs(v_seen - v_half) <= 1.0e-9 * scale, (step, v_seen, v_half, v_full)
        if abs(v_full - v_half) > 1.0e-3 * scale:
            distinct += 1
        u_prev = u
    # non-vacuity: the two velocities really differ, so the assertion above
    # tells one element pass from two
    assert distinct >= 10, distinct


# ---------------------------------------------------------------------------
# S-2 / S-3  rate-independent models: the answer of the two-pass build
# ---------------------------------------------------------------------------
_E, _NU, _SIG0 = 2.0e5, 0.3, 250.0
_K, _G = _E / (3.0 * (1 - 2 * _NU)), _E / (2.0 * (1 + _NU))
_RATE = 1.0e-5                 # m per DR step on the driven face
_N = 400
_FAIL_AT = 150


def _brick_block():
    """2x2x2 LadrunoBrick unit cube of LadrunoJ2, base fixed, top face driven
    down by a Linear-series `sp` (u_z = -t): 4e-3 at step 400, 3x yield."""
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
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for (i, j, k), nd in tag.items():
        if k == n - 1:
            ops.sp(nd, 3, -1.0)
    return tag


def _quad_panel():
    """4x4 LadrunoQuad plane-strain panel of LadrunoJ2, left edge fixed, the
    right edge pulled by nodal loads past yield."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("LadrunoJ2", 1, _K, _G, "-iso", "voce", _SIG0, 0.0, 0.0, 200.0)
    n, h = 5, 0.25
    tag = {}
    t = 1
    for j in range(n):
        for i in range(n):
            ops.node(t, i * h, j * h)
            tag[(i, j)] = t
            t += 1
    e = 1
    for j in range(n - 1):
        for i in range(n - 1):
            ops.element("LadrunoQuad", e, tag[(i, j)], tag[(i + 1, j)],
                        tag[(i + 1, j + 1)], tag[(i, j + 1)], 1,
                        "-type", "PlaneStrain", "-thick", 1.0)
            e += 1
    for j in range(n):
        ops.fix(tag[(0, j)], 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    total = 1.2 * _SIG0 * 1.0                 # 1.2 x the uniaxial yield force
    for j in range(n):
        w = 0.5 if j in (0, n - 1) else 1.0
        ops.load(tag[(n - 1, j)], total * w / (n - 1), 0.0)
    return tag


def _dr(handler, dt, options):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation", "-dt", dt, *options)
    ops.analysis("Transient")


def _brick_observables(tag):
    mid = tag[(1, 1, 1)]
    side = tag[(2, 1, 1)]
    ops.reactions()
    base = sum(ops.nodeReaction(nd, 3) for (i, j, k), nd in tag.items() if k == 0)
    return {
        "time": ops.getTime(),
        "mid_uz": ops.nodeDisp(mid, 3),
        "side_ux": ops.nodeDisp(side, 1),
        "base_rz": base,
        "residual": ops.ladrunoDR("residualNorm"),
    }


def _quad_observables(tag):
    return {
        "tip_ux": ops.nodeDisp(tag[(4, 2)], 1),
        "corner_uy": ops.nodeDisp(tag[(4, 4)], 2),
        "mid_ux": ops.nodeDisp(tag[(2, 2)], 1),
        "kinetic": ops.ladrunoDR("kineticEnergy"),
    }


def _run_brick(options, fail=False):
    tag = _brick_block()
    _dr("Transformation", _RATE, options)
    if not fail:
        assert ops.analyze(_N, _RATE) == 0
    else:
        assert ops.analyze(_FAIL_AT - 1, _RATE) == 0
        # WP-153's injection: a second update() in one step, which DR refuses
        ops.test("NormUnbalance", 1.0e-300, 2, 0)
        ops.algorithm("Newton")
        assert ops.analyze(1, _RATE) < 0
        ops.test("NormUnbalance", 1.0e-12, 1, 0)
        ops.algorithm("Linear")
        assert ops.analyze(_N - (_FAIL_AT - 1), _RATE) == 0
    return _brick_observables(tag)


def _run_quad(options):
    tag = _quad_panel()
    _dr("Plain", 1.0, options)
    assert ops.analyze(600, 1.0) == 0
    return _quad_observables(tag)


_KINETIC = ["-damping", "kinetic"]
_VISCOUS = ["-damping", "viscous", 1.0]

# Recorded with the two-pass build (ladruno f9864cfec, MSVC, 2026-11-19) by
# Ladruno_files/testbed/perf/wp240/single_update_reference.py.
_REFERENCE = {
    "brick_kinetic": {
        "time": 0.004000000000000005,
        "mid_uz": -0.0010768666030345013,
        "side_ux": 0.0009265060081457676,
        "base_rz": 292.1743406134454,
        "residual": 1.6296019983613979,
    },
    "brick_viscous": {
        "time": 0.004000000000000005,
        "mid_uz": -0.0010327960820380355,
        "side_ux": 0.0008799590892033896,
        "base_rz": 281.47677048461264,
        "residual": 12.518332461872081,
    },
    "quad_kinetic": {
        "tip_ux": 0.0027694429580592517,
        "corner_uy": -0.001805291745365371,
        "mid_ux": 0.0008213126918831219,
        "kinetic": 6.886430052033285e-26,
    },
    "quad_viscous": {
        "tip_ux": 0.0012713945342665612,
        "corner_uy": -0.00036142123352239235,
        "mid_ux": 0.0005618720838294374,
        "kinetic": 6.167704483284549e-06,
    },
}

_RTOL, _ATOL = 1.0e-8, 1.0e-12


def _assert_matches(name, got):
    ref = _REFERENCE[name]
    assert ref is not None, "reference for %s was never recorded" % name
    assert sorted(got) == sorted(ref), (sorted(got), sorted(ref))
    for key, value in ref.items():
        assert abs(got[key] - value) <= _RTOL * abs(value) + _ATOL, (
            name, key, got[key], value)


@pytest.mark.parametrize("name,options", [("brick_kinetic", _KINETIC),
                                          ("brick_viscous", _VISCOUS)])
def test_S2_brick_with_prescribed_sp_matches_the_two_pass_build(name, options):
    got = _run_brick(options)
    # the model is a probe: the top face moved, and the block is past yield
    assert abs(got["time"] - _N * _RATE) < 1.0e-12
    assert got["mid_uz"] < -1.0e-3
    _assert_matches(name, got)


@pytest.mark.parametrize("name,options", [("quad_kinetic", _KINETIC),
                                          ("quad_viscous", _VISCOUS)])
def test_S2_quad_under_load_matches_the_two_pass_build(name, options):
    got = _run_quad(options)
    assert got["tip_ux"] > 1.0e-4
    _assert_matches(name, got)


def test_S3_a_failed_step_in_the_middle_lands_on_the_same_reference():
    got = _run_brick(_KINETIC, fail=True)
    _assert_matches("brick_kinetic", got)
