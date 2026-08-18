"""ADR-85 T1b -- `-visc` and `-soft` on the 2D SEGMENT (NTS) lane.

WHY THIS FILE EXISTS
    The splash banner claims both options for the 2D lane and, until this file,
    nothing tested either of them there.  That gap is not symmetric with the
    others in this battery, because the two options sit on OPPOSITE sides of the
    ADR-85 How/6 free-vs-ported line and only one of them was ever going to fail
    loudly:

      `-visc` on SEGMENT is a PORT.  The dashpot loops are stride-3 hardcoded
        (LadrunoContactFE.cpp:700-703 and the addCtoTang twins), so in 2D an
        unported loop reads past the end of a 2-component array -- or, if the
        stride is merely wrong, applies the damper to the wrong component and
        produces a plausible-but-wrong dissipation.  (The RIGID_PLANE dashpot
        is free, loops `d < ndm`, and is already gated by the T0 file's
        test_2d_rigidplane_visc.  That test passing says nothing about this one.)

      `-soft` on SEGMENT rides `softKn`, which the ADR lists as free-but-verify.
        Free code that nobody exercises is how ADR-76's `-initial` finding
        happened, so it is verified rather than assumed.

    Both tests are re-pointed versions of shipped idioms, so their convergence
    and their constants are already known-good:
      -visc  <- tests/test_adr85_contact2d_t0_rigidplane.py::test_2d_rigidplane_visc
               (rebound comparison) on the 3D NTS impact rig of
               tests/test_adr39_contact_p2b.py::test_p2b_explicit_impact_restitution
      -soft  <- tests/test_adr39_contact_p4_soft.py::_segment_loaded /
               ::test_soft_segment_path_stabilized, transcribed to 2D

THE MASTER OVERHANGS, DELIBERATELY
    The master segment runs from x = -5 to +5 under a slave at x = 0, i.e. the
    slave sits at xi = 0.5 with a 5-length margin to either end.  That is the
    lesson the slice twin paid for: a slave at a master boundary makes the test
    measure END POLICY (the terminal acceptance window) instead of the thing
    under test.  Nothing here is about ends, so nothing here goes near one.

TOLERANCES ARE RELATIVE (ADR-85 How/1)
    Rebound speeds are compared to the impact speed and to each other,
    penetrations to the deck's impact penalty depth v*sqrt(m/kn) and to the
    soft-equilibrium P/k_soft, never to a bare absolute.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---- deck constants (every tolerance below is expressed against one of these)
M = 1.0             # slave mass
KN = 1.0e6          # accuracy-driven penalty for the impact rig
MX = 5.0            # master half-span: the slave rides at xi = 0.5, far from both ends
V0 = -2.0           # impact speed
GAP0 = 0.01         # drop height


def _master(kn_dt_safe=True):
    """Two FIXED master nodes spanning [-MX, MX] at y = 0.  Fixed and massless,
    exactly as the shipped 3D SEGMENT rigs do it: an infinite-mass master makes
    m_eff in the SOFT sizing equal the slave mass, which is what keeps the
    k_soft prediction below a closed form."""
    ops.node(101, -MX, 0.0)
    ops.node(102, MX, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)


def _explicit():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")


# --------------------------------------------------------------------------
# 1. -visc on the 2D SEGMENT lane -- the D2.1 dashpot PORT
# --------------------------------------------------------------------------
def _impact(muc, kn=KN, gap0=GAP0, v0=V0, nsteps=1200):
    """Slave mass dropped onto the fixed master segment, central difference.

    Returns the history summary the gate reads.  `-outward 0 1` is given
    explicitly: the slave starts above the master so the centroid vote would
    resolve, but pinning it removes one variable from a test about damping.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _master()
    ops.node(1, 0.0, gap0)
    ops.mass(1, M, M)
    ops.fix(1, 1, 0)                      # free y only (the pair is rank-1)
    ops.contactSurface(20, "-slave", 1)
    extra = ("-visc", muc) if muc > 0.0 else ()
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 1.0, *extra)
    ops.setNodeVel(1, 2, v0, "-commit")
    _explicit()
    dt = 0.5 * 2.0 * math.sqrt(M / kn)    # half the contact CFL, the shipped idiom
    max_pen, v_out, max_absv, blew = 0.0, 0.0, 0.0, False
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        y = gap0 + ops.nodeDisp(1, 2)
        vy = ops.nodeVel(1, 2)
        max_absv = max(max_absv, abs(vy))
        if y < 0.0:
            max_pen = max(max_pen, -y)
        if y > 0.5 * gap0 and vy > 0.0:
            v_out = vy
    return dict(blew=blew, max_pen=max_pen, v_out=v_out, max_absv=max_absv)


def test_2d_segment_visc_damps_rebound():
    """`-visc mu_c` on a 2D NTS pair dissipates, and dissipates by the right
    amount.

    THE PORT UNDER TEST.  The SEGMENT dashpot reads the normal approach rate as
    gdot = B.v over a stride-3 B (LadrunoContactFE.cpp:700-703).  In 2D that
    stride is wrong by construction, so the failure modes are (a) reading past
    a 2-component slave velocity, or (b) picking up the wrong component and
    damping something that is not the gap rate.  Mode (b) is the dangerous one:
    it still converges, still looks damped-ish, and is wrong.  So the gate is
    not "some energy was lost" but a QUANTITATIVE bracket.

    The pair is a linear spring-dashpot during contact, so the rebound ratio is
    the classic single-mode result

        e = exp(-pi*zeta/sqrt(1 - zeta^2)),   zeta = mu_c/(2*sqrt(kn*M))

    At zeta = 0.3 that is 0.372.  The bracket below is deliberately wide (a
    factor of 2 either side of the closed form) because the contact is
    discretised at omega*dt = 1 and the slave arrives with the master fixed, so
    the measured restitution carries a few percent of integration error -- but
    it is narrow enough to exclude both "no damping" (e = 1) and "damping the
    wrong component / double-counting" (e far below the bracket).

    Three further claims keep the bracket honest: the undamped control must be
    nearly elastic, the damper must never AMPLIFY, and neither leg may diverge.
    """
    zeta = 0.3
    muc = 2.0 * zeta * math.sqrt(KN * M)
    free = _impact(0.0)
    damped = _impact(muc)

    assert not free["blew"] and not damped["blew"], "a 2D -visc impact diverged"
    e_free = free["v_out"] / abs(V0)
    e_damped = damped["v_out"] / abs(V0)
    e_ref = math.exp(-math.pi * zeta / math.sqrt(1.0 - zeta * zeta))
    print(f"\n[2D SEGMENT -visc] zeta={zeta} mu_c={muc:.4e}  e_free={e_free:.4f} "
          f"e_damped={e_damped:.4f} (closed form {e_ref:.4f})")

    assert e_free > 0.9, (
        f"the UNDAMPED control is not elastic (e={e_free:.4f}) -- the rig is not "
        "measuring what this test thinks it is")
    assert 0.5 * e_ref < e_damped < 2.0 * e_ref, (
        f"damped restitution {e_damped:.4f} outside the bracket "
        f"[{0.5 * e_ref:.4f}, {2.0 * e_ref:.4f}] around the spring-dashpot "
        f"closed form {e_ref:.4f} -- the SEGMENT dashpot is damping the wrong "
        "quantity (the stride-3 loop was not ported to ndm)")
    assert e_damped < 0.7 * e_free, (
        f"-visc did not dissipate: e {e_damped:.4f} vs undamped {e_free:.4f}")
    assert damped["max_absv"] <= free["max_absv"] * 1.05, (
        f"-visc AMPLIFIED the response: max|v| {damped['max_absv']:.3e} vs "
        f"{free['max_absv']:.3e} -- a dashpot removes energy, never adds it")


# --------------------------------------------------------------------------
# 2. -soft on the 2D SEGMENT lane -- the B1 Courant-stable penalty
# --------------------------------------------------------------------------
KN_STIFF = 1.0e8    # omega*dt = sqrt(kn/M)*dt = 10 >> 2 at DT_STRUCT
DT_STRUCT = 1.0e-3  # a structural-scale step, stable for everything but contact
RUNAWAY = 1.0       # |u_y| beyond this (>> any physical penetration) == diverged


def _held(kn, soft, dt=DT_STRUCT, P=1.0e3, v=-0.5, nsteps=2000):
    """The p4_soft `_segment_loaded` rig in 2D: the slave starts ON the master
    under a SUSTAINED downward load (so contact stays engaged and an unstable
    contact mode compounds) and is kicked downward.  Returns (diverged, max|uy|).

    `-outward 0 1` is REQUIRED here, not optional: the slave sits exactly on the
    master, so the interface-level centroid datum is degenerate and the vote
    would refuse the deck (ADR-85 How/2).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _master()
    ops.node(1, 0.0, 0.0)
    ops.mass(1, M, M)
    ops.fix(1, 1, 0)
    ops.contactSurface(20, "-slave", 1)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 1.0, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    ops.setNodeVel(1, 2, v, "-commit")
    _explicit()
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs                     # NaN / Inf abort
        maxabs = max(maxabs, abs(ops.nodeDisp(1, 2)))
        if maxabs > RUNAWAY:
            return True, maxabs
    return False, maxabs


def test_2d_segment_soft_stabilizes_explicit():
    """`-soft` on a 2D NTS pair sizes the penalty from mass and dt, so a kn that
    DIVERGES at the structural step stays bounded.

    ADR-85 How/6 lists the normal SOFT sizing as free-but-verify: `softKn`'s
    inputs (gapModeInvMass / ladrunoNodeMass / ladrunoInvMassProj) are all
    size-guarded for 2-DOF nodes, and under the `ndf == ndm` guard slot 2 is
    structurally zero.  "Free" is a prediction; this is the measurement.  The
    argument k_soft = SOFSCL*4*m_eff/dt^2 is itself dimension-independent
    (a Courant bound), so if 2D needed anything the failure would be in the
    plumbing, not the formula.

    dt is IDENTICAL across the two legs -- only SOFT differs -- which is what
    makes the comparison a statement about the penalty and not about the step.
    The bound on the soft response is a closed form, not a fudge: with a fixed
    (infinite-mass) master m_eff = M, so

        k_soft = SOFSCL*4*M/dt^2 = 0.1*4*1/1e-6 = 4.0e5
        equilibrium penetration = P/k_soft = 2.5e-3

    and the response is that plus the kick amplitude v/omega, both well inside
    the 0.05 bound.  The 3D sibling (test_adr39_contact_p4_soft::
    test_soft_segment_path_stabilized) asserts the same three things on the same
    numbers.
    """
    div_stiff, _ = _held(KN_STIFF, 0.0)
    div_soft, maxabs = _held(KN_STIFF, 0.1)
    k_soft = 0.1 * 4.0 * M / (DT_STRUCT * DT_STRUCT)
    print(f"\n[2D SEGMENT -soft] k_soft={k_soft:.4e}, equilibrium "
          f"{1.0e3 / k_soft:.4e}, measured max|uy|={maxabs:.4e}")

    assert div_stiff, (
        "a stiff penalty must DIVERGE at the structural dt (omega*dt = 10 >> 2) "
        "-- if it does not, this deck no longer gates anything")
    assert not div_soft, (
        f"-soft must stabilize the 2D SEGMENT lane at the same dt "
        f"(max|uy|={maxabs:.3e})")
    assert maxabs < 0.05, (
        f"soft response {maxabs:.3e} is not bounded near the soft equilibrium "
        f"{1.0e3 / k_soft:.3e} -- k_soft was sized from something other than "
        "the slave mass and the step")
