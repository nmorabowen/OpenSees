"""LadrunoRebarBuckling -- fiber-section / structural integration (Zone-A).

The Zone-A material battery (test_ladrunoRebarBuckling_material.py) drives the
WRAPPER directly. These tests put it where it is actually used: as a rebar fiber
inside a fiber section, under a cyclic structural analysis through a real
beam/zero-length element and Newton solve.

  * B4j -- zeroLengthSection moment-curvature: cyclic curvature drives the extreme
    bar fibers past the buckling onset; the wrapped (lsr>0) section loses moment
    capacity on the buckled half-cycles vs the identity (lsr=0) section, while the
    tension-dominated small-curvature response is unchanged. Both solve cleanly.
  * B4k -- forceBeamColumn RC cantilever: an axially-loaded RC column (Concrete02
    core/cover + wrapped rebar layers) under a growing cyclic tip drift -- the
    integration robustness check (Newton converges every step) plus a soft
    capacity comparison vs the identity-gate column.

See Ladruno_implementation/14_ladruno_rebar_buckling_adr.md (§9).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_FY = 420.0
_E = 200000.0
_B = 0.01
_EY = _FY / _E


def _steel(tag):
    ops.uniaxialMaterial("Steel02", tag, _FY, _E, _B, 18.0, 0.925, 0.15)


# ==========================================================================
#  B4j -- zeroLengthSection moment-curvature, wrapped vs identity rebar
# ==========================================================================
def _bar_mat(tag, steel_tag, lsr):
    """Rebar fiber material: RebarBuckling(lsr) wrapping Steel02 (lsr=0 = identity)."""
    _steel(steel_tag)
    ops.uniaxialMaterial("LadrunoRebarBuckling", tag, steel_tag,
                         "-lsr", lsr, "-fy", _FY, "-E", _E,
                         "-restraighten", "lambda")


def _section_cyclic(lsr, kappa_path, ybar=0.23, Abar=500.0, nsub=25):
    """Drive a 2-bar fiber section (zeroLengthSection) through a curvature path
    under DisplacementControl on the rotation DOF. Returns [(kappa, M)] or None
    if Newton ever fails to converge."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 0)                       # axial + shear fixed => pure curvature on dof 3

    _bar_mat(1, 2, lsr)                       # bar material tag 1 (wraps Steel02 tag 2)
    ops.section("Fiber", 1)
    ops.fiber(+ybar, 0.0, Abar, 1)           # top bar
    ops.fiber(-ybar, 0.0, Abar, 1)           # bottom bar
    ops.element("zeroLengthSection", 1, 1, 2, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 1.0)               # reference moment on dof 3
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 3, kappa_path[0] / nsub)
    ops.analysis("Static")

    out = []
    cur = 0.0
    for target in kappa_path:
        d = (target - cur) / nsub
        ops.integrator("DisplacementControl", 2, 3, d)
        for _ in range(nsub):
            if ops.analyze(1) != 0:
                return None
        cur = target
        ops.reactions()
        out.append((target, ops.nodeReaction(1, 3)))   # base moment reaction
    return out


@pytest.mark.t1
def test_B4j_fiber_section_cyclic_buckling():
    lsr = 20.0                               # eStar floors at ~ -7 eY
    ybar = 0.23
    # DM reduction begins at the YIELD strain -eY (intermediate branch); the deep
    # branch starts at eStar (~ -7 eY). So pre-onset == extreme bar strain < eY.
    k_small = 0.8 * _EY / ybar               # extreme-bar strain ~ 0.8 eY (elastic, pass-through)
    k_big = 12.0 * _EY / ybar                # extreme-bar strain ~ 12 eY (deep buckling)
    path = [k_small, -k_small, k_big, -k_big, k_big, -k_big]

    wrap = _section_cyclic(lsr, path)
    iden = _section_cyclic(0.0, path)        # identity gate (no buckling)
    assert wrap is not None, "wrapped section: Newton failed to converge"
    assert iden is not None, "identity section: Newton failed to converge"

    # pre-onset peak (index 0, +k_small): bars not buckled => moments match
    assert abs(wrap[0][1]) == pytest.approx(abs(iden[0][1]), rel=1e-3), (
        "pre-onset moment should equal the identity section")

    # buckled peaks (the +/-k_big peaks, indices 2..5): the compressed bar sheds
    # capacity => the wrapped section moment magnitude is clearly below identity.
    for idx in (2, 3, 4, 5):
        mw, mi = abs(wrap[idx][1]), abs(iden[idx][1])
        assert mw < 0.85 * mi, (
            f"peak {idx}: buckling did not degrade the section "
            f"(|M_wrap|={mw:.1f} not < 0.85*|M_id|={mi:.1f})")

    # the wrapped section dissipates LESS per cycle than the un-buckled section
    # (the compressed bar sheds load): loop area over the two big cycles is smaller.
    def _area(series):
        a = 0.0
        for (k0, m0), (k1, m1) in zip(series, series[1:]):
            a += 0.5 * (m0 + m1) * (k1 - k0)
        return abs(a)
    assert _area(wrap[1:]) < _area(iden[1:]), "buckled loop should enclose less area"


# ==========================================================================
#  B4k -- forceBeamColumn RC cantilever cyclic (integration robustness)
# ==========================================================================
def _rc_section(sec_tag, bar_mat_tag):
    """Simple RC rectangular fiber section: Concrete02 core/cover + 2 wrapped
    rebar layers (top/bottom). b x h = 400 x 500 mm."""
    b, h, cov = 400.0, 500.0, 40.0
    fc, ec, fcu, ecu = -30.0, -0.002, -6.0, -0.01     # Concrete02 (N/mm^2)
    ops.uniaxialMaterial("Concrete02", 50, fc, ec, fcu, ecu, 0.1, 3.0, 1500.0)
    yh = h / 2.0
    yc = yh - cov
    ops.section("Fiber", sec_tag)
    ops.patch("rect", 50, 1, 12, -yh, -b / 2.0, yh, b / 2.0)        # concrete
    ops.layer("straight", bar_mat_tag, 3, 500.0, +yc, -b / 2.0 + cov, +yc, b / 2.0 - cov)
    ops.layer("straight", bar_mat_tag, 3, 500.0, -yc, -b / 2.0 + cov, -yc, b / 2.0 - cov)


def _advance(node, dof, d_total, base_steps):
    """Displacement-control advance with adaptive step-halving (the standard way
    to push a cyclic RC section). Returns True if the full increment is reached."""
    done = 0.0
    dU = d_total / base_steps
    floor = abs(d_total) / base_steps / 64.0
    guard = 0
    while abs(d_total) - abs(done) > 1e-9:
        guard += 1
        if guard > 100000:
            return False
        step = dU
        if abs(done + step) > abs(d_total):          # don't overshoot the target
            step = d_total - done
        ops.integrator("DisplacementControl", node, dof, step)
        if ops.analyze(1) == 0:
            done += step
            dU = d_total / base_steps                 # restore the nominal step after a success
        else:
            dU *= 0.5
            if abs(dU) < floor:
                return False
    return True


def _rc_cantilever_cyclic(lsr, drifts, H=3000.0, axial=-600.0e3, nsub=30):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, H)
    ops.fix(1, 1, 1, 1)

    _bar_mat(1, 2, lsr)
    _rc_section(1, 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.geomTransf("Linear", 1)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)

    # gravity / axial, applied gradually
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, axial, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        return None, False

    # cyclic lateral tip drift (DisplacementControl on node 2 dof 1, adaptive)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 1.0, 0.0, 0.0)
    peaks = []
    cur = 0.0
    ok = True
    for dr in drifts:
        target = dr * H
        if not _advance(2, 1, target - cur, nsub):
            ok = False
            break
        cur = target
        ops.reactions()
        peaks.append((dr, ops.nodeReaction(1, 1)))   # base shear
    return peaks, ok


@pytest.mark.t1
def test_B4k_rc_cantilever_cyclic_converges():
    # growing symmetric drift cycles up to ~3% (adaptive sub-stepping)
    drifts = [0.005, -0.005, 0.015, -0.015, 0.03, -0.03]
    wrap, ok_w = _rc_cantilever_cyclic(20.0, drifts)
    assert wrap is not None, "wrapped RC column: gravity diverged"
    assert ok_w and len(wrap) == len(drifts), (
        "wrapped RC column: cyclic analysis did not complete all drift cycles")

    iden, ok_i = _rc_cantilever_cyclic(0.0, drifts)
    assert iden is not None and ok_i and len(iden) == len(drifts)
    # soft capacity check: at the largest drift the buckled column carries no MORE
    # base shear than the identity column (bar buckling can only reduce capacity).
    assert abs(wrap[-1][1]) <= abs(iden[-1][1]) * 1.02 + 1.0
