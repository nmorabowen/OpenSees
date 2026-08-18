"""ADR-85 T2 -- explicit friction: energy balance, SOFT pounding, and the
`-visc` statics-inertness regression (G-T2 b / c / d).

WHAT THIS FILE IS
    T1b's `test_adr85_contact2d_t1b_explicit.py` gated `-visc` and `-soft` on
    the 2D SEGMENT lane FRICTIONLESS.  T2 is the first phase where friction is
    live on that same explicit lane, and three things become checkable for
    the first time: (b) the force-only explicit friction dissipates the RIGHT
    amount of energy; (c) `-soft`'s TANGENTIAL sizing (`softKt`) is provably
    DEAD CODE until a 2D pair can slip at all (LadrunoContactFE.cpp:773-783,
    "Unreachable until T2 wires 2D friction" -- _adr85_t2_design.md Sec 2)
    and this file is what makes it live; (d) `-visc`'s statics-inertness must
    still hold once a friction TANGENT is also being assembled under the same
    StaticIntegrator call (T1b only ever proved `-visc` inert alongside a
    frictionless tangent).

    All three tests are written against the not-yet-built C++ scalar friction
    path (see test_adr85_contact2d_t2_friction.py's module docstring for the
    pre-T2 decision) and are EXPECTED TO FAIL right now on the `-mu > 0`
    FATAL (LadrunoContactHandler.cpp:995-999) or its `-soft` pre-scan
    ordering quirk (deferral D6, see (c) below) -- never silently.

COVERAGE MAP -> ADR-85 G-T2 / deferrals
    (b) explicit sliding energy balance: dissipation == mu*N*slip
            test_2d_explicit_friction_energy_balance
    (c) SOFT pounding stable at dt = 0.9*dt_cr, SOFSCL default, softKt LIVE
            test_2d_soft_pounding_stable_with_friction
    (d) -visc statics-inert BYTE-IDENTITY, WITH a live friction tangent
            test_2d_visc_static_inert_with_friction_tangent
    D6 (mentioned, not separately gated here): the `-soft` pre-scan
        refusal-ordering cleanup is a MESSAGE-ATTRIBUTION issue (a 2D `-soft`
        deck is refused by the `-soft` subsystem before the ADR-85 named
        refusal, loud but wrongly attributed) -- it does not change whether
        test_2d_soft_pounding_stable_with_friction can pass once T2 lands, so
        it is not separately gated; if T2 leaves D6 unfixed, this test's
        FATAL will simply name the wrong subsystem, which is a message-text
        detail this test does not assert on (no REFUSALS substring check
        here -- it asserts on numeric stability, not on why a pre-friction
        build refuses to run at all).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# ---------------------------------------------------------------------------
# 1. explicit sliding energy balance -- G-T2(b)
# ---------------------------------------------------------------------------
KN = 1.0e6
KT = 1.0e5


def _driven_block_2d(mu, Q, P=1.0e3, m=1.0, nsteps=150, dt=1.0e-3):
    """The 2D analogue of 3D's test_adr39_contact_p3_friction.py::_driven_block:
    a flat 2-node master (normal +y), a mass node pre-penetrated to the normal
    equilibrium (y0 = -P/kn, so N == P with no vertical ringing), driven by a
    constant +x force above the cone. Fully free (no `fix`) -- the contact
    itself supplies the only normal restraint, exactly like the 3D rig."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    L = 20.0
    ops.node(101, -L / 2.0, 0.0)
    ops.node(102, L / 2.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.node(1, 0.0, -P / KN)
    ops.mass(1, m, m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "2D explicit driven-block step failed"
    return ops.nodeDisp(1, 1), ops.nodeVel(1, 1)


def test_2d_explicit_friction_energy_balance():
    """G-T2(b): dissipation == mu*N*slip, measured from the work-energy
    theorem on the SOLVED state (x, v), not from the return map's own
    internal bookkeeping -- an independent check on the force-only explicit
    friction (CentralDifferenceLadruno never forms a friction tangent, so this
    is purely the FORCE path, LadrunoContactFE.cpp's SEGMENT-mode addB).

    W_drive (Q*x) = KE (0.5*m*v^2) + dissipation, so
        dissipation_measured = Q*x - 0.5*m*v^2
    compared against mu*N*x (the cap force times the total slid distance --
    valid because Q > mu*N throughout, so the block is ALWAYS in the slip
    branch from the first step, never re-sticking, so x IS the accumulated
    slip distance, not merely the total displacement of a partly-elastic
    motion).

    TOLERANCE 5%, same as the 3D precedent this ports
    (test_adr39_contact_p3_friction.py::test_p3_energy_dissipation, which
    established this bound for the identical numerical scheme).  The
    justification carries over unchanged: dt=1e-3 central difference has
    O(dt^2) local truncation error in position but the FRICTION FORCE
    direction is decided from a half-step-lagged velocity/displacement
    estimate (force-only, no tangent), so the sliding-then-settling dynamics
    of a driven block accumulate a few percent of integration error over
    nsteps=150 -- the SAME numerical mechanism in 2D as in 3D (same scalar
    return map's force branch, same explicit scheme), so the same bound
    applies without re-deriving it.
    """
    P, m, mu, Q = 1.0e3, 1.0, 0.3, 6.0e2       # mu*P=300 < Q=600 -> continuous slip
    cap = mu * P
    x_fric, v_fric = _driven_block_2d(mu, Q, P=P, m=m)
    _x_free, v_free = _driven_block_2d(0.0, Q, P=P, m=m)
    assert v_fric < v_free, (
        f"friction must slow the 2D block: v_fric={v_fric!r} !< v_free={v_free!r} "
        "(a sign-flipped friction force would speed it up -- the ADR-39 P3 "
        "sign gate, replayed in 2D)")
    W_diss = Q * x_fric - 0.5 * m * v_fric * v_fric
    ref = cap * x_fric
    assert abs(W_diss - ref) / ref < 0.05, (
        f"energy balance: measured dissipation {W_diss:.6e} != mu*N*slip "
        f"{ref:.6e} (x={x_fric:.6e}, cap={cap:.6e})")


# ---------------------------------------------------------------------------
# 2. SOFT pounding with friction (softKt goes LIVE) -- G-T2(c)
# ---------------------------------------------------------------------------
M = 1.0
DT_CR = 1.0e-3       # the "structural" reference dt (ADR-85 How/6's k_soft is
                     # SIZED FROM whatever dt it is fed, ω·dt = 2√SOFSCL, so
                     # -soft's OWN stability is dt-invariant by construction --
                     # the number that matters here is a STRUCTURAL-scale dt
                     # chosen independent of the contact penalty, matching the
                     # T1b -soft precedent's own DT_STRUCT
                     # (test_adr85_contact2d_t1b_explicit.py), not a contact CFL)
DT = 0.9 * DT_CR
KN_RAW = 1.0e8       # an UNSIZED penalty representative of a naive -kn choice
                     # at this problem's scale; sqrt(KN_RAW/M)*DT ~ 30 >> 2, so
                     # it MUST diverge at DT without SOFT sizing -- the same
                     # idiom as the T1b -soft test's KN_STIFF
RUNAWAY = 5.0        # |disp| beyond this (>> any physical impact penetration
                     # or slip travel here) means the run blew up


def _pounding(soft, mu, kn=KN_RAW, kt=KN_RAW, dt=DT, nsteps=1500,
              gap0=0.02, vy0=-3.0, vx0=1.5):
    """An OBLIQUE 2D impact: a free mass approaches a fixed master line with
    BOTH a normal (vy0) and a tangential (vx0) approach velocity -- a
    pounding event is rarely purely normal, and the tangential component is
    what puts softKt on the hot path (a purely normal impact would only ever
    exercise softKn, which T1b already gated frictionless).  `-soft` at
    SOFSCL default (bare flag, no numeric arg -> 0.10,
    OpenSeesOutputCommands.cpp:632-650).  Returns (stable, maxabs) over both
    DOFs."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -5.0, 0.0)
    ops.node(102, 5.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.node(1, 0.0, gap0)
    ops.mass(1, M, M)
    ops.contactSurface(20, "-slave", 1)
    extra = ("-soft",) if soft else ()
    ops.contact(1, 10, 20, kn, kt, mu, "-outward", 0.0, 1.0, *extra)
    ops.setNodeVel(1, 1, vx0, "-commit")
    ops.setNodeVel(1, 2, vy0, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return False, maxabs
        maxabs = max(maxabs, abs(ops.nodeDisp(1, 1)), abs(ops.nodeDisp(1, 2)))
        if maxabs > RUNAWAY:
            return False, maxabs
    return True, maxabs


def test_2d_soft_pounding_stable_with_friction():
    """G-T2(c): an oblique 2D pounding impact stays STABLE at dt = 0.9*DT_CR
    under `-soft` (SOFSCL default) with friction ON, exercising softKt's 2D
    branch for the first time (it is provably dead code until a 2D pair can
    take the SLIP branch at all -- LadrunoContactFE.cpp:773-783's own comment,
    quoted in the module docstring).

    THE CONTRAST, mirroring the T1b -soft precedent
    (test_adr85_contact2d_t1b_explicit.py::test_2d_segment_soft_stabilizes_explicit):
    an UNSIZED raw penalty (kn=kt=KN_RAW=1e8) at the SAME dt, SAME oblique
    impact, SAME friction MUST diverge (omega*dt = sqrt(KN_RAW/M)*DT = 9.0,
    4.5x the central-difference stability limit of 2) -- that divergence
    is the proof DT_CR is genuinely "too aggressive for a raw penalty", not an
    accident of the rig, and the `-soft` leg's stability at the IDENTICAL dt
    with the IDENTICAL kn/kt inputs is therefore attributable to the SOFT
    sizing (both softKn on the normal approach AND softKt on the tangential
    slip -- the impact has both components by construction), not to a
    conveniently gentle dt choice.

    No closed-form bound on the bounded response magnitude is asserted here
    (unlike T1b's `-soft` test, which had a clean SOFSCL-derived equilibrium
    penetration under a SUSTAINED load): a genuine two-component IMPACT's
    peak penetration/slip depends on the impact geometry in a way a sustained
    load's equilibrium does not, so the gate is BOUNDEDNESS (stays inside
    RUNAWAY = 5.0, four orders above the mm-scale gap0/velocities here) rather
    than a numeric target -- the same "stable vs blows up" contrast the T1b
    precedent itself leads with before adding its closed-form check.
    """
    mu = 0.4
    raw_stable, raw_max = _pounding(soft=False, mu=mu)
    soft_stable, soft_max = _pounding(soft=True, mu=mu)
    assert not raw_stable, (
        f"the raw (unsized) penalty must DIVERGE at dt=0.9*DT_CR under an "
        f"oblique impact (max|disp|={raw_max:.3e}) -- if it does not, DT_CR is "
        "not actually aggressive relative to KN_RAW and this deck gates nothing")
    assert soft_stable, (
        f"-soft (SOFSCL default) with friction ON must stay STABLE at the "
        f"identical dt (max|disp|={soft_max:.3e}) -- softKt's 2D branch was "
        "either not wired or not Courant-sized")
    assert soft_max < RUNAWAY, f"soft response {soft_max:.3e} exceeded the runaway bound"


# ---------------------------------------------------------------------------
# 3. -visc statics-inertness WITH a live friction tangent -- G-T2(d)
# ---------------------------------------------------------------------------
# Mirrors tests/test_adr41_viscous_d2.py::test_d2_static_byte_identical (the
# shipped 3D D2 regression contract), extended to a configuration where a
# friction STICK tangent is ALSO being assembled in the same static solve --
# T1b only ever proved -visc inert alongside a frictionless (kn*B^T*B-only)
# tangent; this is the first static deck where the two new/newish T2 code
# paths (visc's addCtoTang no-op path is unaffected by transient-vs-static,
# but the ASSEMBLED tangent it would-be contribute to now also carries the
# friction stick block) run side by side.
VISC_KN = 1.0e6
VISC_KT = 1.0e5


def _visc_static_state(muc):
    """A 2D NTS pair in STICK (mu comfortably above any trial tangential
    force) under a StaticIntegrator, with `-visc muc` declared. Velocity is
    identically zero throughout a static solve (nodal velocities are never
    advanced by a StaticIntegrator), so the viscous term p_visc = muc*gdot
    must be IDENTICALLY zero regardless of muc -- this is a structural
    argument (v=0 makes the product zero, not a numerical coincidence), so
    the comparison is asserted at STRICT byte-identity (`==`), the same
    strength as the shipped ADR-41 D2 precedent, not a tolerance."""
    P, Q, mu = 1.0e3, 50.0, 0.5              # mu*P=500 >> Q=50 -> stick, well inside the cone
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -2.0, 0.0)
    ops.node(102, 2.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.node(1, 0.0, -P / VISC_KN)
    ops.contactSurface(20, "-slave", 1)
    extra = ("-visc", muc) if muc > 0.0 else ()
    ops.contact(1, 10, 20, VISC_KN, VISC_KT, mu, "-outward", 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "visc+friction stick static solve did not converge"
    return ops.nodeDisp(1, 1), ops.nodeDisp(1, 2)


def test_2d_visc_static_inert_with_friction_tangent():
    """G-T2(d): under a StaticIntegrator, `-visc mu_c` is IDENTICALLY inert
    (bit-for-bit) on a 2D SEGMENT pair that is SIMULTANEOUSLY sticking under
    friction -- the previously-untested combination of two T2-adjacent code
    paths sharing one assembled tangent. `ladrunoContactForce` is NOT used
    here (it is the NTS-lane NORMAL-force query only, per ADR-85 How/7 and
    the contact_dump.py harness convention); the comparison is on the raw
    nodal displacement, which is what a genuine viscous perturbation (however
    small) would move first.
    """
    ux0, uy0 = _visc_static_state(0.0)
    ux1, uy1 = _visc_static_state(5.0e3)
    assert ux0 == ux1 and uy0 == uy1, (
        f"-visc perturbed a static (stick+friction) 2D solve: "
        f"({ux0!r},{uy0!r}) vs ({ux1!r},{uy1!r}) with mu_c=5e3 -- velocity is "
        "identically zero under statics, so this must be bit-identical")
