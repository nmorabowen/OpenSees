"""ADR-57 E3 — edge-edge Coulomb/Tresca FRICTION (the 4th consumer of LadrunoFrictionKernel).

The EDGE_EDGE mode now adds a tangential traction on top of the E2 normal penalty: the slip is the
tangential relative DISPLACEMENT of the two closest points (the C3.1 trap — at the closest point the
relative POSITION is purely normal), run through the SHIPPED return map with the per-pair committed
slip on the EdgeEdgeState. -edgeMu/-edgeKt/-edgeCohesion/-edgeTauMax form the unified cone
min(μN+c, τmax); the symmetric tangent is solver-safe (the Coulomb Csl behind -edgeConsistentTan).
Mechanics validated oracle-first in contact_prototypes/proto_e3_friction.py.

Two perpendicular bars: a slave bar whose bottom edge (∥x) crosses a master bar's top edge (∥y),
n=ẑ. The bottom edge is held at a fixed penetration (z fixed ⇒ a CONSTANT normal pressure N=εN|gN|)
and RAKES in y across the master edge — so the friction acts on a clean 1-DOF tangential mode.

Gates:
  - EXPLICIT raking bar: Coulomb friction DECELERATES the rake (opposes motion), vs μ=0 which coasts;
  - IMPLICIT stick: a sub-cone tangential load sticks ⇒ static Newton converges, slip bounded;
  - the symmetric tangent is SOLVER-SAFE (ProfileSPD == FullGeneral);
  - μ=0 ⇒ byte-identical to the frictionless E2 path (the caller short-circuits the friction block).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e5      # edge-edge normal penalty εN  ⇒ N = εN·|ZC|
KT = 1.0e4       # tangential (stick) penalty
ZC = -1.0e-4     # held penetration ⇒ constant normal pressure N = EPS·|ZC| = 10
N_PRESS = EPS * abs(ZC)


def _build_bars():
    """Master (fixed) bar in the y-z plane, top edge m1-m2 ∥ y at z=0. Slave bar in the x-z plane,
    bottom edge s1-s2 ∥ x at z=ZC. n = ê_s×ê_m = ẑ."""
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.0, ZC + 0.3), (8, -1.0, 0.0, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)               # anchor the slave top edge


def _edge_contact(mu, kt=KT, friction=True, cohesion=0.0):
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    args = [1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if friction:
        if mu > 0.0:
            args += ["-edgeMu", mu]
        if cohesion > 0.0:
            args += ["-edgeCohesion", cohesion]
        args += ["-edgeKt", kt]
    ops.contact(*args)


def _add_y_springs(k=100.0):
    """Soft y-springs on the bottom nodes 5,6. A point-like edge pair restrains only the COMBINED
    tangential gap mode (weights (1−s),s shared over the two nodes), so the differential y-mode is a
    zero-energy mode under the IMPLICIT static tangent — regularize it (the friction stick stiffness
    still dominates the response). Under EXPLICIT the nodal mass regularizes it, so the rake needs no
    springs. (Same point-like trap E2 handled for the normal mode; inherent to penalty contact.)"""
    ops.uniaxialMaterial("Elastic", 1, k)
    for gnd, sn, x in [(405, 5, -1.0), (406, 6, 1.0)]:
        ops.node(gnd, x, 0.0, ZC)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 800 + sn, gnd, sn, "-mat", 1, "-dir", 2)


def _rake_velocity_history(mu, nsteps=60, dt=1.0e-3, v0=0.5):
    """EXPLICIT rake: bottom edge held at z=ZC (x,z fixed; constant N), y free with mass + a +y
    velocity. Returns the y-velocity history of node 5. μ>0 ⇒ decelerates; μ=0 ⇒ coasts at v0."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (5, 6):
        ops.fix(t, 1, 0, 1)               # x,z fixed (hold penetration); y free (the rake DOF)
        ops.mass(t, 1.0, 1.0, 1.0)
        ops.setNodeVel(t, 2, v0, "-commit")
    _edge_contact(mu, friction=(mu > 0.0))
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    vy = [ops.nodeVel(5, 2)]
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        vy.append(ops.nodeVel(5, 2))
    return vy


def test_e3_explicit_rake_coulomb_opposes_motion():
    """THE HEADLINE: a raking bottom edge with Coulomb friction DECELERATES (friction opposes the
    +y motion), while μ=0 coasts at v0. The gap between them is the assembled friction force."""
    v0 = 0.5
    vy_fric = _rake_velocity_history(mu=0.5, v0=v0)
    vy_free = _rake_velocity_history(mu=0.0, v0=v0)
    # μ=0: no friction ⇒ coasts at v0 (x,z fixed ⇒ no other y-force)
    assert abs(vy_free[-1] - v0) < 1.0e-6, f"μ=0 should coast at {v0}, got {vy_free[-1]}"
    # μ>0: friction opposes ⇒ velocity strictly decreases and ends well below the coast
    assert vy_fric[-1] < 0.8 * v0, f"Coulomb did not decelerate the rake ({vy_fric[-1]} vs {v0})"
    assert all(vy_fric[i + 1] <= vy_fric[i] + 1e-9 for i in range(len(vy_fric) - 1)), \
        "friction velocity not monotonically decreasing (sign error?)"
    assert vy_fric[-1] >= 0.0, "friction over-shot to negative velocity (energy injection?)"


def _stick_static(system="FullGeneral"):
    """IMPLICIT static: bottom edge held at z=ZC (x,z fixed) + soft y-springs (regularize the point-
    pair differential mode), a sub-cone +y load ⇒ friction STICKS ⇒ Newton converges at a small
    elastic slip. cap = μN + c = 0.3·10 + 20 = 23 > load 10 ⇒ stick. Returns the y-disp of node 5."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (5, 6):
        ops.fix(t, 1, 0, 1)               # x,z fixed; y free (held by friction stick + the y-spring)
    _add_y_springs()
    _edge_contact(0.3, cohesion=20.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 5.0, 0.0)                    # total +y load 10 < cap 23 ⇒ stick
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "implicit friction-stick Newton did not converge (tangent?)"
    return ops.nodeDisp(5, 2)


def test_e3_implicit_stick_converges():
    """A sub-cone tangential load STICKS ⇒ static Newton converges at a bounded elastic slip — the
    friction-stick tangent (the SYMMETRIC kt·P_t) is what makes the solve converge."""
    y = _stick_static("FullGeneral")
    # stick ⇒ elastic slip is small & bounded (kt-dominated), NOT a runaway slip across the cone.
    assert 0.0 < y < 0.05, f"stick slip not small/bounded ({y}) — did it slip the cone?"


def test_e3_symmetric_tangent_solver_safe():
    """The default friction tangent is SYMMETRIC ⇒ a symmetric solver (ProfileSPD) gives the SAME
    answer as the general solver (FullGeneral). (The non-symmetric Csl is behind -edgeConsistentTan.)"""
    y_spd = _stick_static("ProfileSPD")
    y_gen = _stick_static("FullGeneral")
    assert abs(y_spd - y_gen) < 1.0e-10, f"symmetric tangent not solver-safe ({y_spd} vs {y_gen})"


def test_e3_mu0_byte_identical():
    """μ=0 ⇒ the friction block short-circuits ⇒ BITWISE identical to the frictionless E2 path
    (no -edgeMu/-edgeKt flags). An explicit rake history compared bitwise."""
    vy_mu0 = _rake_velocity_history(mu=0.0)            # -edgeMu omitted (friction=False)

    def _run_explicit_mu0_flag():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _build_bars()
        for t in (5, 6):
            ops.fix(t, 1, 0, 1)
            ops.mass(t, 1.0, 1.0, 1.0)
            ops.setNodeVel(t, 2, 0.5, "-commit")
        # -edgeMu 0 explicitly: must take the SAME short-circuit path as omitting it
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
                    "-edgeMu", 0.0, "-edgeKt", KT, "-outward", 0.0, 0.0, 1.0)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
        ops.analysis("Transient")
        vy = [ops.nodeVel(5, 2)]
        for _ in range(60):
            assert ops.analyze(1, 1.0e-3) == 0
            vy.append(ops.nodeVel(5, 2))
        return vy
    assert vy_mu0 == _run_explicit_mu0_flag()
