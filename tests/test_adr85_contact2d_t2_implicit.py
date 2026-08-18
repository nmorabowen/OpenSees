"""ADR-85 T2 -- the implicit-transient 2D friction lane (deferral D3).

WHAT THIS FILE IS
    T1b's panel deferrals (_adr85_t2_design.md Sec... via the ADR's
    Implementation-log) assign D3 to T2: "the `-initial` (addKiToTang) 2D arm
    gains its first test in this PR (ModifiedNewton `-initial` patch rerun);
    the full implicit-transient 2D lane (addCtoTang under HHT/GN) is T2's
    gate." A StaticIntegrator never calls `addCtoTang` at all (contact's own
    comment in LadrunoContactFE.h: "contact has no damping -> no-op" -- true
    for the FRICTIONLESS/undamped case, but ADR-85 How/6 and the SEGMENT-mode
    code both show the muc*n(x)n `-visc` damping TANGENT is assembled THERE,
    LadrunoContactFE.cpp ~997), so addCtoTang's 2D arm is genuinely UNTESTED
    by every static test in this battery, friction or not -- only a
    TRANSIENT implicit integrator (Newmark / HHT) reaches it.

    "A stride-3 index on a 2-DOF node reads the WRONG equation" (the task
    brief's own words) is the sharp failure mode this file has to catch: 2D
    nodes have ndf==2, so any leftover `for (d=0; d<3; d++)` or hardcoded
    `(0,1,2)` DOF loop copied from the 3D SEGMENT branch reads one equation
    past a 2-DOF node's own block and corrupts an ADJACENT node's tangent/
    residual entry. That is not "slightly wrong" -- it either crashes
    (out-of-bounds), or silently assembles into the WRONG unknown, producing
    an answer with the right ORDER OF MAGNITUDE but the wrong VALUE (a
    coupling that should not exist). Every test below therefore checks
    against a closed form tight enough that "off by a misrouted DOF" fails
    loudly, not a vague "did it converge" check.

    Written against not-yet-built C++; every test here is expected to fail on
    the `-mu > 0` FATAL (or, for the -visc+Newmark test, to fail for a
    DIFFERENT reason -- see that test's docstring) until T2 lands.

COVERAGE MAP -> ADR-85 D3
    ModifiedNewton -initial (addKiToTang), 2D friction STICK
            test_2d_modifiednewton_initial_friction_stick
    Newmark implicit dynamic slide (addKtToTang under a transient
    integrator -- the "GN" alpha-family half of "HHT/GN")
            test_2d_newmark_implicit_friction_slide
    HHT-alpha implicit dynamic slide, PLUS a live `-visc` damping tangent
    (the only test in the whole 2D battery that puts a REAL (nonzero)
    contribution through addCtoTang, per the note above)
            test_2d_hht_visc_friction_addctotang
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
KT = 1.0e5


def _flat_master():
    ops.node(101, -10.0, 0.0)
    ops.node(102, 10.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)


# ---------------------------------------------------------------------------
# 1. ModifiedNewton -initial (addKiToTang) -- D3
# ---------------------------------------------------------------------------
def test_2d_modifiednewton_initial_friction_stick():
    """D3, the `-initial` half. Direct 2D port of the 3D precedent
    (test_adr39_contact_p35_implicit.py::test_p35_static_stick_converges,
    solved under `-initial` instead of full Newton) and of the T1b analogue
    that pinned the 2D addKiToTang arm frictionless
    (test_adr85_contact2d_t1b_nts.py::test_2d_patch_initial_tangent).

    A free 2-DOF slave (no anchor: friction stick provides the ONLY
    tangential stiffness), driven below the cone, under ModifiedNewton
    `-initial`. `K_initial == K_c` for a flat/fixed master (no geometric
    term, ADR-85 How/5), so `addKiToTang` should mirror `addKtToTang` exactly
    and converge to the IDENTICAL stick closed form x=Q/kt as full Newton --
    with MORE iterations allowed (modified Newton is linearly convergent
    where exact Newton is not, for a problem that is exactly linear once the
    stick/slip branch is fixed) but not a different answer.

    SHARPNESS AGAINST THE STRIDE-3 HAZARD. `addKiToTang` populating the WRONG
    equation (a stride-3 leftover) has two distinguishable failure shapes,
    both caught here: (a) the 2x2 local block lands on the wrong global DOF
    pair entirely, leaving the TRUE x/y DOFs under-stiffened -> the
    ModifiedNewton iteration either fails to converge within the raised
    iteration budget or drifts off the linear closed form by an O(1)
    fraction, not a rounding-level amount; (b) it silently drops the
    tangential (kt) contribution and keeps only kn -> the iteration matrix is
    RANK-DEFICIENT in the tangential direction and modified Newton cannot
    converge to a nonzero x at all (residual stalls). Neither shape passes
    the 1e-6 closed-form check below.
    """
    P, Q, mu = 1.0e3, 200.0, 0.5             # mu*P=500 > Q=200 -> stick
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(1, 0.0, -P / KN)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-10, 200, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.analysis("Static")
    rc = ops.analyze(1)
    assert rc == 0, "2D friction stick under ModifiedNewton -initial must converge"
    x = ops.nodeDisp(1, 1)
    x_ref = Q / KT
    assert abs(x - x_ref) / x_ref < 1.0e-6, (
        f"-initial stick x={x!r} != Q/kt={x_ref!r} (full-Newton closed form) "
        "-- addKiToTang's 2D arm assembled the wrong stiffness")
    z = ops.nodeDisp(1, 2)
    assert abs(z) < 1.0e-9, (
        f"the normal DOF moved under -initial ({z!r}) -- a stride-3 leftover "
        "coupling the tangential update into the normal equation")


# ---------------------------------------------------------------------------
# 2. Newmark implicit dynamic slide (addKtToTang under a transient
#    integrator) -- D3
# ---------------------------------------------------------------------------
def test_2d_newmark_implicit_friction_slide():
    """D3, direct 2D port of test_adr39_contact_p35_implicit.py::
    test_p35_newmark_dynamic_friction: a mass driven above the cone under
    implicit Newmark (gamma=0.5, beta=0.25 -- no numerical damping), slides
    at a=(Q-mu*N)/m, exactly as under explicit CDL
    (test_adr85_contact2d_t2_explicit.py::test_2d_explicit_friction_energy_balance
    uses the same closed form), but reached through addKtToTang under a
    TRANSIENT integrator rather than addMtoTang-only (explicit) or a single
    static LoadControl step (the stick tests above) -- the third of the
    three tangent-assembly shapes the ADR names (static / transient-implicit
    / explicit), and the one a StaticIntegrator-only battery cannot touch.

    Tolerance 5%: same numerical-scheme argument as the explicit energy-
    balance test and its 3D Newmark precedent (test_p35_newmark_dynamic_friction
    used the identical 5% bound for the identical measurement -- a
    finite-difference velocity slope over a sliding window, now with Newmark's
    own O(dt^2) local truncation error in place of central difference's).
    """
    P, Q, mu, m = 1.0e3, 6.0e2, 0.3, 1.0      # mu*P=300 < Q=600 -> slides at a=300
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(1, 0.0, -P / KN)
    ops.mass(1, m, m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Transient")
    vx = []
    for _ in range(120):
        assert ops.analyze(1, dt) == 0, "2D Newmark frictional step must converge"
        vx.append(ops.nodeVel(1, 1))
    a = (vx[110] - vx[40]) / (70 * dt)
    a_th = (Q - mu * P) / m
    assert abs(a - a_th) / a_th < 0.05, f"2D implicit slide a={a!r} != {a_th!r}"


# ---------------------------------------------------------------------------
# 3. HHT-alpha + a LIVE -visc damping tangent -- D3, addCtoTang for real
# ---------------------------------------------------------------------------
def test_2d_hht_visc_friction_addctotang():
    """D3, the sharpest addCtoTang test in this battery. Contact contributes
    NOTHING to addCtoTang when `-visc` is absent (the header's own "no
    damping -> no-op" comment) -- every OTHER transient test in this file
    would pass even if addCtoTang's 2D arm were entirely missing, because a
    no-op has no wrong answer to expose. This test is the one place a
    genuinely WRONG (or stride-3-corrupted) addCtoTang 2D implementation has
    something to corrupt: `-visc mu_c` assembles a real muc*(n_x)n tangent
    contribution there (ADR-41 D2 / the SEGMENT-mode comment near
    LadrunoContactFE.cpp:997), and HHT-alpha (unlike Newmark
    gamma=0.5/beta=0.25) forms its effective tangent as a nontrivial
    combination of K, C, and M (c1*K + c2*C + c3*M with c2 != 0), so a
    corrupted or missing 2D C-block changes the CONVERGED answer, not just
    the iteration count.

    RIG: the T1b `-visc` impact rig
    (test_adr85_contact2d_t1b_explicit.py::_impact), reused under HHT-alpha
    instead of CentralDifferenceLadruno, WITH friction (mu>0, oblique
    approach so the tangential channel is live) stacked on top of the
    viscous normal dashpot -- both new-to-T2 (friction) and T1b-shipped
    (visc) code share one assembled operator here.

    CLOSED FORM: the classic damped-spring-mass restitution used by the T1b
    -visc test, e = exp(-pi*zeta/sqrt(1-zeta^2)), zeta = mu_c/(2*sqrt(kn*m)).
    Bracket identical to the T1b precedent's own reasoning (a factor of 2
    either side of the closed form: HHT's numerical dissipation and the
    coarser dt relative to the contact period both add a few percent of
    integration error on top of the physical damping, same order as central
    difference's -- the bracket exists to separate "damped roughly right" from
    "not damped" or "damping the wrong DOF", not to pin the number tightly).
    """
    M_, KN_ = 1.0, 1.0e6
    zeta = 0.3
    muc = 2.0 * zeta * math.sqrt(KN_ * M_)
    mu, gap0, vy0, vx0 = 0.2, 0.01, -2.0, 0.3
    dt = 0.5 * 2.0 * math.sqrt(M_ / KN_)

    def _impact(with_visc):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        _flat_master()
        ops.node(1, 0.0, gap0)
        ops.mass(1, M_, M_)
        ops.contactSurface(20, "-slave", 1)
        extra = ("-visc", muc) if with_visc else ()
        ops.contact(1, 10, 20, KN_, KT, mu, "-outward", 0.0, 1.0, *extra)
        ops.setNodeVel(1, 1, vx0, "-commit")
        ops.setNodeVel(1, 2, vy0, "-commit")
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("HHT", 0.9)
        ops.test("NormDispIncr", 1.0e-10, 40, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        max_pen, v_out = 0.0, 0.0
        for _ in range(1200):
            assert ops.analyze(1, dt) == 0, "2D HHT visc+friction impact diverged"
            y = gap0 + ops.nodeDisp(1, 2)
            vy = ops.nodeVel(1, 2)
            if y < 0.0:
                max_pen = max(max_pen, -y)
            if y > 0.5 * gap0 and vy > 0.0:
                v_out = vy
        return v_out, max_pen

    v_free, _ = _impact(False)
    v_damped, _ = _impact(True)
    e_free = v_free / abs(vy0)
    e_damped = v_damped / abs(vy0)
    e_ref = math.exp(-math.pi * zeta / math.sqrt(1.0 - zeta * zeta))
    assert e_free > 0.85, (
        f"undamped HHT control not elastic enough (e={e_free:.4f}) -- HHT's own "
        "numerical dissipation (alpha=0.9) is already visible; retune alpha "
        "if this fails on a real build")
    assert 0.5 * e_ref < e_damped < 2.0 * e_ref, (
        f"HHT+visc restitution {e_damped:.4f} outside "
        f"[{0.5 * e_ref:.4f}, {2.0 * e_ref:.4f}] around the closed form "
        f"{e_ref:.4f} -- addCtoTang's 2D arm (with friction ALSO live) is not "
        "damping the normal channel correctly")
    assert e_damped < 0.8 * e_free, (
        f"-visc did not measurably dissipate under HHT+friction: damped "
        f"{e_damped:.4f} vs free {e_free:.4f}")
