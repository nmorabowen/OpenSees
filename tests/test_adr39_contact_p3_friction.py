"""ADR-39 P3 — Coulomb friction on the shipped frictionless NTS penalty contact
(P2b/P2.5). v1 = the EXPLICIT ship: friction FORCE only (the tangent is mass-only
under CentralDifferenceLadruno, so no friction tangent is assembled — P3.5). The
return map is validated oracle-first in contact_prototypes/proto_p3_friction.py (6/6);
here it runs in compiled OpenSees through the REAL assembly (CDL + addB), which is
the only way to catch the design-gate SIGN BLOCKER (friction must OPPOSE the motion;
a missing negation injects energy → a=g(sinθ+μcosθ)).

All tests are explicit (CDL) to match the v1 ship leg. Models are element-free (a mass
node + fixed master corners) — they also exercise the B1 element-free teardown fix.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
KT = 1.0e5


# --------------------------------------------------------------------------- utils
def _horiz_master(L=20.0, base=101):
    """4 FIXED corner nodes of a flat quad in the z=0 plane (normal +z)."""
    h = L / 2.0
    tags = []
    for i, (x, y) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        t = base + i
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags


def _tilt_master(theta, L=1.0, base=101):
    """4 FIXED corners of a quad rotated about the x-axis by theta. Outward normal
    n = (0, -sinθ, cosθ) (toward the +z half-space where the slave rests)."""
    c, s = math.cos(theta), math.sin(theta)
    h = L / 2.0
    tags = []
    for i, (x, y0) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        t = base + i
        ops.node(t, x, y0 * c, y0 * s)     # rotate (x,y0,0) about x by theta
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags


def _run_cdl(nsteps, dt, sample=None):
    """Step CDL nsteps; if sample (a callable) is given, record sample() each step."""
    rec = []
    for k in range(nsteps):
        assert ops.analyze(1, dt) == 0, f"step {k} failed"
        if sample is not None:
            rec.append(sample(k))
    return rec


# --------------------------------------------------------------------- the tests
def test_p3_incline_sliding_acceleration():
    """THE sign gate (design test #3): a block on a θ-incline under gravity slides at
    a = g(sinθ − μcosθ) — measured through real CDL+addB (NOT the oracle's hand
    subtraction). A missing friction negation would give a = g(sinθ + μcosθ). The
    frictionless leg (μ=0) slides at g·sinθ (covers the regression direction too)."""
    g, m, theta = 9.81, 2.0, math.radians(30.0)
    c, s = math.cos(theta), math.sin(theta)
    n = (0.0, -s, c)                               # outward normal
    dhat = (0.0, -c, -s)                           # unit down-slope direction
    W = m * g
    N = W * c                                      # held normal force
    dt, k1, k2 = 1.0e-3, 80, 180

    def slope_accel(mu):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        mt = _tilt_master(theta, L=1.0)
        dn = N / KN                                # normal-equilibrium pre-penetration
        ops.node(1, 0.0, dn * s, -dn * c)          # start at -dn*n (normal balanced)
        ops.mass(1, m, m, m)
        ops.contactSurface(10, "-master", 4, *mt)
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, KT, mu, "-outward", *n)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, 0.0, -W)                  # gravity
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
        ops.analysis("Transient")

        def vslope(_k):
            return (ops.nodeVel(1, 1)*dhat[0] + ops.nodeVel(1, 2)*dhat[1]
                    + ops.nodeVel(1, 3)*dhat[2])
        rec = _run_cdl(k2 + 1, dt, sample=vslope)
        return (rec[k2] - rec[k1]) / ((k2 - k1) * dt)

    a_fric = slope_accel(0.3)
    a_free = slope_accel(0.0)
    a_th_fric = g * (s - 0.3 * c)                  # 2.356
    a_th_free = g * s                              # 4.905
    assert abs(a_fric - a_th_fric) / a_th_fric < 0.03, f"a_fric={a_fric} vs {a_th_fric}"
    assert abs(a_free - a_th_free) / a_th_free < 0.03, f"a_free={a_free} vs {a_th_free}"


def _driven_block(mu, Q, P=1000.0, m=1.0, nsteps=150, k1=40, k2=140, L=20.0):
    """Horizontal plane (normal +z), slave pre-penetrated to the normal equilibrium
    (z=-P/kn ⇒ N=P, no vertical ringing), driven by a constant +x force Q and held
    down by -P. Returns (a_x, last vx, last x, max|x|) measured through real CDL."""
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _horiz_master(L=L)
    z0 = -P / KN
    ops.node(1, 0.0, 0.0, z0)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    vx, xmax = [], 0.0
    for k in range(nsteps):
        assert ops.analyze(1, dt) == 0
        vx.append(ops.nodeVel(1, 1))
        xmax = max(xmax, abs(ops.nodeDisp(1, 1)))
    a_x = (vx[k2] - vx[k1]) / ((k2 - k1) * dt)
    return a_x, vx[-1], ops.nodeDisp(1, 1), xmax


def test_p3_slip_caps_at_muN():
    """Slip caps the tangential force at μN (design tests #1-2): a driven block slides
    at a=(Q−μN)/m, so the friction force f=Q−m·a == μN — INDEPENDENT of Q. Two drives
    above the cone give the same μN (the cap does not scale with the drive)."""
    P, m, mu = 1000.0, 1.0, 0.3
    cap = mu * P                                   # 300
    for Q, a_th in [(600.0, 300.0), (900.0, 600.0)]:
        a_x, _, _, _ = _driven_block(mu, Q, P=P, m=m)
        assert abs(a_x - a_th) / a_th < 0.02, f"Q={Q}: a={a_x} vs {a_th}"
        f = Q - m * a_x                            # measured friction force
        assert abs(f - cap) / cap < 0.02, f"Q={Q}: friction {f} != μN {cap}"


def test_p3_stick_below_cone():
    """Below the cone (Q < μN, and 2Q < μN so the undamped overshoot stays sub-cone)
    the block is HELD by friction — it elastically oscillates about Q/kt and never
    slides away (design test #1 stick). max|x| stays within the undamped overshoot
    2Q/kt; a BROKEN no-friction impl free-accelerates to x≈O(1) (caught), and a
    LEAKING stick (spurious plastic slip) would exceed the tight bound below."""
    P, m, mu, Q = 1000.0, 1.0, 0.3, 100.0          # μN=300, 2Q=200 < 300 ⇒ stick
    _, _, _, xmax = _driven_block(mu, Q, P=P, m=m, nsteps=300, k1=100, k2=200)
    assert xmax < 1.2 * (2.0 * Q / KT), f"stick should hold tightly, max|x|={xmax}"


def test_p3_energy_dissipation():
    """Friction dissipates energy (design test #4). Two independent checks: (a) SIGN —
    the block's exit speed with friction is LESS than frictionless (a sign-flipped
    friction would speed it up); (b) ENERGY BALANCE from the SOLVED state variables —
    drive work Q·x equals kinetic energy ½m·v² plus frictional dissipation μN·x
    (work-energy theorem, measured from x and v directly, not from the asserted a)."""
    P, m, mu, Q = 1000.0, 1.0, 0.3, 600.0
    cap = mu * P
    _, v_fric, x_fric, _ = _driven_block(mu, Q, P=P, m=m)
    _, v_free, _, _ = _driven_block(0.0, Q, P=P, m=m)
    assert v_fric < v_free, f"friction must slow the block: {v_fric} !< {v_free}"
    W_diss = Q * x_fric - 0.5 * m * v_fric * v_fric     # drive work − KE = dissipation
    assert abs(W_diss - cap * x_fric) / (cap * x_fric) < 0.05, \
        f"energy balance: dissipation {W_diss} != μN·x {cap * x_fric}"


def _static_normal_pen(mu, P=1.0e3):
    """Static: slave just-penetrated under a fixed horizontal quad, normal load -P,
    tangential DOFs fixed (free z only). Returns the normal penetration. The friction
    branch IS entered for mu>0 (slot lookup + return map) but returns zero traction
    (no tangential motion) — it must not perturb the normal result. (A statically-
    ACTIVE friction leg is infeasible in v1: with no friction tangent a free tangential
    DOF is singular — which is exactly why the ship is explicit; active friction is
    covered by the explicit slide/stick/incline tests.)"""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _horiz_master(L=2.0)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)                            # free z only
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return -(-1.0e-8 + ops.nodeDisp(1, 3))


def test_p3_frictionless_regression():
    """μ never perturbs the NORMAL penetration (design test #5). μ=0 takes the short-
    circuit path (byte-identical to the shipped P2b); μ=0.5 enters the friction branch
    (return map runs, returns zero stick traction). Both give exactly P/kn, identical
    to each other — friction does not leak into the normal solve."""
    P = 1.0e3
    expect = P / KN
    pen_mu0 = _static_normal_pen(0.0, P=P)
    pen_mu5 = _static_normal_pen(0.5, P=P)
    for pen in (pen_mu0, pen_mu5):
        assert pen > 0 and abs(pen - expect) / expect < 1e-6, f"pen={pen} expect={expect}"
    assert abs(pen_mu0 - pen_mu5) / expect < 1e-9, f"μ perturbed the normal: {pen_mu0} vs {pen_mu5}"


def test_p3_commit_fires():
    """The friction state commits at the Domain::commit choke point: after N explicit
    steps of a frictional slide, ladrunoContactInfo reports N commits (design test #6
    — the hook that promotes gpT_trial → gpT each step actually fires)."""
    nsteps = 50
    P, m, mu, Q = 1000.0, 1.0, 0.3, 600.0
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _horiz_master(L=20.0)
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
    info = ops.ladrunoContactInfo()                # [numContacts, numCommits, numReverts]
    assert info[0] == 1 and info[1] == nsteps, f"contactInfo={info}, expected commits={nsteps}"


def test_p3_wipe_reanalyze():
    """Engine lifetime + friction-state GC: running a frictional slide, wiping, and
    re-running in the SAME process gives the SAME acceleration with no crash (the
    wipe deletes the engine; the adapter re-fetches it lazily; the GC prunes stale
    slots). Element-free model ⇒ also re-exercises the B1 teardown fix."""
    P, m, mu, Q = 1000.0, 1.0, 0.3, 600.0
    a1, _, _, _ = _driven_block(mu, Q, P=P, m=m)
    a2, _, _, _ = _driven_block(mu, Q, P=P, m=m)
    assert abs(a1 - 300.0) / 300.0 < 0.02, f"run1 a={a1}"
    assert abs(a2 - 300.0) / 300.0 < 0.02, f"run2 a={a2}"
    assert abs(a1 - a2) < 1.0e-6, f"reanalysis drifted: {a1} vs {a2}"


def test_p3_late_engagement_gt0():
    """ENGAGEMENT-referenced slip (design BLOCKER/MAJOR-1, gT0). A slave that drifts
    tangentially WHILE SEPARATED, then penetrates, must feel ZERO friction at the
    engagement step — the pre-contact drift is NOT a slip. (A global-reference slip
    would inject a spurious stick traction kt·Δ at first contact → a large tangential
    acceleration spike.) The slave flies in +x and falls; at the first penetrating
    step its x-acceleration must be ≈0 (no tangential load, gTeff=0 at capture)."""
    mu, m, gap0 = 0.3, 1.0, 0.01
    V, Vz, dt = 0.5, 2.0, 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _horiz_master(L=20.0)
    ops.node(1, 0.0, 0.0, gap0)                    # starts ABOVE the plane (separated)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.setNodeVel(1, 1, V, "-commit")             # drift +x while separated
    ops.setNodeVel(1, 3, -Vz, "-commit")           # fall toward the plane
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    ax_engage, x_engage, bug_scale = None, None, None
    for _ in range(60):
        assert ops.analyze(1, dt) == 0
        z = gap0 + ops.nodeDisp(1, 3)
        if z < 0.0 and ax_engage is None:          # first penetrating step
            ax_engage = ops.nodeAccel(1, 1)
            x_engage = ops.nodeDisp(1, 1)
            N = KN * (-z)
            bug_scale = KT * x_engage / m           # spurious a_x a global-ref slip gives
            break
    assert ax_engage is not None, "slave never engaged"
    assert x_engage > 1.0e-4, f"slave should have drifted before engaging, x={x_engage}"
    # correct: ~0; the global-reference bug would give ~ -KT·Δ/m (orders larger)
    assert abs(ax_engage) < 0.05 * abs(bug_scale), \
        f"spurious engagement friction: a_x={ax_engage} vs bug-scale {bug_scale}"


def test_p3_revert_path():
    """Friction state reverts on an aborted step (design MAJOR-7). Explicit CDL has no
    auto-retry, so an adaptive-step wrapper must drive Domain::revertToLastCommit() —
    here a deliberately non-converging implicit step triggers it. The revert hook must
    fire (numReverts increments), dropping the rejected trial slip."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _horiz_master(L=20.0)
    ops.node(1, 0.0, 0.0, -1.0e-3)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-30, 1, 0)          # impossible ⇒ the step aborts
    ops.algorithm("Newton")
    ops.analysis("Transient")
    rc = ops.analyze(1, 1.0e-3)
    info = ops.ladrunoContactInfo()                # [numContacts, numCommits, numReverts]
    assert rc != 0, "step should have failed to converge"
    assert info[2] >= 1, f"revert hook must fire on an aborted step, contactInfo={info}"


def _cube(base_node, z0, L, E, nu, ele_tag, mat_tag, rho=0.0):
    """8-node cube (base_node..+7), bottom at z0, edge L. Returns (bottom4, top4)."""
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        t = base_node + i
        ops.node(t, float(x), float(y), z0 + z)
        tags.append(t)
    ops.nDMaterial("ElasticIsotropic", mat_tag, E, nu, rho)
    ops.element("LadrunoBrick", ele_tag, *tags, mat_tag)
    return tags[:4], tags[4:]


def test_p3_deformable_master_drag():
    """The MASTER friction reaction block (resid_master_i += −N_i·tFric) — never
    exercised when the master is fixed. A slave slides under friction on the top face
    of a DEFORMABLE brick (fixed bottom); the friction must DRAG the brick top in the
    slide direction (+x), shearing it. A sign error in the master block would drag it
    the wrong way (or not at all). Explicit CDL + a little stiffness-Rayleigh damping
    to settle the shear; assert the brick top drifted +x and the slave slid forward
    relative to it."""
    E, nu, rho, L = 2.0e4, 0.0, 1.0, 1.0
    mu, P, Q, m = 0.3, 1.0e3, 6.0e2, 1.0
    dt = 2.0e-4
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    bot, top = _cube(101, 0.0, L, E, nu, ele_tag=1, mat_tag=1, rho=rho)
    for t in bot:
        ops.fix(t, 1, 1, 1)
    ops.node(1, L / 2, L / 2, L - P / KN)          # slave on top face, pre-penetrated
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)                         # Q>μP ⇒ slides
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("ProfileSPD")
    ops.integrator("CentralDifferenceLadruno")
    ops.rayleigh(0.0, 0.0, 1.0e-4, 0.0)            # stiffness damping settles the shear
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(1500):
        assert ops.analyze(1, dt) == 0
    top_x = sum(ops.nodeDisp(t, 1) for t in top) / 4.0
    slave_x = ops.nodeDisp(1, 1)
    assert top_x > 1.0e-5, f"friction must drag the master top +x, got {top_x}"
    assert slave_x > top_x, f"slave should slide forward of the master, {slave_x} !> {top_x}"
