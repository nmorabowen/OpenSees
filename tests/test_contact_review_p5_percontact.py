"""Contact review fix P5 (2026-07) — per-contact NormalField sign re-key + hygiene batch.

THE HEADLINE (MED, five-agent review): the ADR-63 `-smoothNormal` nodal-normal field was
keyed by MASTER SURFACE tag, but the frozen global outward SIGN it carries is a per-CONTACT
datum (parity with the faceted per-contact orientDir). A SECOND contact sharing the master —
the two-sided baffle: slaves pressing on BOTH faces of one declared plate — inherited the
FIRST contact's frozen sign, and because `signCaptured` skipped the whole vote, even its own
explicit `-outward` was IGNORED. The wrong-signed field reads the second contact's slaves as
"already penetrating" from the wrong side and EJECTS them THROUGH the plate (silent
pass-through). Post-fix the field store is keyed by CONTACT tag: each contact votes and
freezes its own sign.

THE HYGIENE BATCH gated here too:
  - one-time warnings moved from process-lifetime function-local statics to per-(contact,
    topic) engine latches (a second contact with the same defect now warns too; a wipe/new
    model re-warns);
  - `-visc` is DISABLED (with a warning) under a StaticIntegrator — trial velocities are
    stale STATE there (a static stage after a transient one keeps the last velocities), so
    the dashpot injected an unphysical, tangent-less force;
  - mortar AUTO orientation (scen − mcen) degenerating on COINCIDENT conforming facets now
    warns (the sign was silently left to facet-winding luck; ties exempt — sign-independent).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e5


# ===================================================== (1) two-sided baffle, shared master
def _run_two_sided_baffle():
    """One flat quad master (a baffle plate at z=0) shared by TWO -smoothNormal contacts:
    contact 1 = slave A ABOVE (auto sign — the per-slave robust vote gives +z), contact 2 =
    slave B BELOW with an explicit -outward (0,0,-1). Both slaves are pressed INTO the plate
    (A down, B up), x,y fixed (only z free — the sign-isolation discipline). Returns the
    extreme z-excursions (min zA, max zB): a correct per-contact sign HOLDS each slave at its
    own face (|z| ~ P/kn); the pre-fix shared-master sign EJECTS slave B through the plate."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x, y in [(1, -1.0, -1.0), (2, 1.0, -1.0), (3, 1.0, 1.0), (4, -1.0, 1.0)]:
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
    zA0, zB0 = 1.0e-3, -1.0e-3
    ops.node(100, 0.3, 0.0, zA0); ops.mass(100, 1.0, 1.0, 1.0); ops.fix(100, 1, 1, 0)
    ops.node(200, -0.3, 0.0, zB0); ops.mass(200, 1.0, 1.0, 1.0); ops.fix(200, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(20, "-slave", 100)
    ops.contactSurface(30, "-slave", 200)
    # contact 1 (processed FIRST — it captures the shared-master field pre-fix)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-smoothNormal")
    # contact 2: same master, slave below, its own -outward toward ITS side (down)
    ops.contact(2, 10, 30, KN, 0.0, 0.0, "-outward", 0.0, 0.0, -1.0, "-smoothNormal")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(100, 0.0, 0.0, -20.0)   # A pressed DOWN onto the top face
    ops.load(200, 0.0, 0.0, 20.0)    # B pressed UP onto the bottom face
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    zA_min, zB_max = zA0, zB0
    for _ in range(3000):
        if ops.analyze(1, dt) != 0:
            break
        zA = zA0 + ops.nodeDisp(100, 3)
        zB = zB0 + ops.nodeDisp(200, 3)
        if not (math.isfinite(zA) and math.isfinite(zB)):
            break
        zA_min = min(zA_min, zA)
        zB_max = max(zB_max, zB)
    return zA_min, zB_max


def test_p5_two_sided_baffle_per_contact_sign():
    """Each contact sharing the master gets ITS OWN frozen sign: slave A is held at the top
    face AND slave B at the bottom face. Pre-fix, contact 2 inherited contact 1's +z field
    (its -outward ignored), read slave B as penetrating from the wrong side, and ejected it
    UP through the plate (zB → large positive)."""
    zA_min, zB_max = _run_two_sided_baffle()
    assert zA_min > -0.05, f"slave A passed through the baffle from above (min z = {zA_min})"
    assert zB_max < 0.05, (
        f"slave B was ejected through the baffle (max z = {zB_max}) — the second contact "
        "inherited the first contact's frozen sign (shared-master keying)")


# ============================================ (2) warning latches are per-CONTACT (engine)
def test_p5_warning_latch_is_per_contact(capfd):
    """TWO mortar contacts, each with friction but no -epsT: EACH must emit the
    'friction requested without -epsT' warning. Pre-fix the function-local static latched
    after the FIRST contact anywhere in the process, silencing the second (and every later
    model in the same interpreter)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nid = 1
    for c, x0 in [(1, 0.0), (2, 10.0)]:      # two independent facet pairs, far apart
        mt, st = [], []
        for x, y in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
            ops.node(nid, x0 + x, y, 0.0); ops.fix(nid, 1, 1, 1); mt.append(nid); nid += 1
        for x, y in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
            ops.node(nid, x0 + x, y, 0.5); ops.mass(nid, 1.0, 1.0, 1.0); st.append(nid); nid += 1
        ops.contactSurface(10 * c, "-master", 4, *mt)
        ops.contactSurface(10 * c + 1, "-slave-segments", 4, *st)
        ops.contact(c, 10 * c, 10 * c + 1, "-mortar", "-epsN", KN, "-mu", 0.3,
                    "-outward", 0.0, 0.0, 1.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-4) == 0
    err = capfd.readouterr().err
    assert "mortar contact 1: friction requested without -epsT" in err, \
        f"contact 1 -epsT warning missing (got: {err[-400:]!r})"
    assert "mortar contact 2: friction requested without -epsT" in err, (
        "contact 2 -epsT warning missing — the one-time latch is process-wide, not "
        f"per-contact (got: {err[-400:]!r})")


# ======================================================= (3) -visc under a StaticIntegrator
def _static_plane_press(muc, stale_vz):
    """The shipped D2 static rig (node ON the plane, x,y fixed, pressed down) + a STALE
    committed velocity (as a static stage following a transient one would have).
    Returns the converged z-displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 0)
    if stale_vz != 0.0:
        ops.setNodeVel(1, 3, stale_vz, "-commit")
    ops.contactSurface(20, "-slave", 1)
    args = [1, 20, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN]
    if muc > 0.0:
        args += ["-visc", muc]
    ops.contactPlane(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "static plane press did not converge"
    return ops.nodeDisp(1, 3)


def test_p5_visc_inert_under_static_integrator_with_stale_velocity(capfd):
    """THE REGRESSION: with a non-zero COMMITTED velocity (stale state from a prior transient
    stage), a static solve with -visc must equal the no-visc solve — under statics the trial
    velocity is not a rate, and the damping tangent is never assembled. Pre-fix the dashpot
    read the stale velocity and injected a constant unphysical force, shifting the static
    equilibrium. A one-time warning must surface the disablement."""
    z_off = _static_plane_press(muc=0.0, stale_vz=-5.0)
    z_on = _static_plane_press(muc=50.0, stale_vz=-5.0)
    err = capfd.readouterr().err
    assert z_on == z_off, (
        f"-visc perturbed a STATIC solve via a stale committed velocity: {z_on!r} vs {z_off!r}")
    assert "-visc" in err and "static integrator" in err, \
        f"the -visc static disablement did not warn (got: {err[-300:]!r})"


# ========================================= (4) coincident-facet mortar AUTO orientation warn
def _coincident_pair(outward, tie=False):
    """A conforming COINCIDENT mortar facet pair (slave facet node-coincident with the master
    facet — the zero-gap tie-style start where scen − mcen degenerates). One explicit step."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt, st = [], []
    nid = 1
    for x, y in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
        ops.node(nid, x, y, 0.0); ops.fix(nid, 1, 1, 1); mt.append(nid); nid += 1
    for x, y in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
        ops.node(nid, x, y, 0.0); ops.mass(nid, 1.0, 1.0, 1.0); st.append(nid); nid += 1
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(11, "-slave-segments", 4, *st)
    args = [1, 10, 11, "-mortar", "-epsN", KN]
    if tie:
        args.append("-tie")
    if outward:
        args += ["-outward", 0.0, 0.0, 1.0]
    ops.contact(*args)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-5) == 0
    return


def test_p5_coincident_mortar_auto_orientation_warns(capfd):
    """Coincident conforming facets WITHOUT -outward: the auto orientation scen − mcen is
    degenerate and the sign falls to facet-winding luck — the handler must say so (pre-fix:
    silence; the PR-3 gate follow-up)."""
    _coincident_pair(outward=False)
    err = capfd.readouterr().err
    assert "auto orientation is DEGENERATE" in err, \
        f"coincident-facet auto-orientation warning missing (got: {err[-300:]!r})"


def test_p5_staggered_zero_gap_mortar_auto_orientation_warns(capfd):
    """The gate lens-B extension: a zero-gap STAGGERED (non-conforming, laterally offset)
    pair — the canonical mortar interface — has orientDir TANGENTIAL to the interface, so
    the facet-normal flip test is equally on roundoff. The warning predicate tests the
    NORMAL COMPONENT n̂_m·orientDir, not |orientDir|, so this must warn too."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt, st = [], []
    nid = 1
    for x, y in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
        ops.node(nid, x, y, 0.0); ops.fix(nid, 1, 1, 1); mt.append(nid); nid += 1
    for x, y in [(0.4, 0.0), (1.4, 0.0), (1.4, 1.0), (0.4, 1.0)]:   # same plane, offset +0.4x
        ops.node(nid, x, y, 0.0); ops.mass(nid, 1.0, 1.0, 1.0); st.append(nid); nid += 1
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(11, "-slave-segments", 4, *st)
    ops.contact(1, 10, 11, "-mortar", "-epsN", KN)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-5) == 0
    err = capfd.readouterr().err
    assert "auto orientation is DEGENERATE" in err, \
        f"staggered zero-gap auto-orientation warning missing (got: {err[-300:]!r})"


def test_p5_coincident_mortar_outward_and_tie_stay_silent(capfd):
    """The warning is targeted: an explicit -outward (the documented escape) and a -tie pair
    (sign-independent full 3-vec bond — the COMMON coincident case) must NOT warn."""
    _coincident_pair(outward=True)
    err1 = capfd.readouterr().err
    _coincident_pair(outward=False, tie=True)
    err2 = capfd.readouterr().err
    assert "auto orientation is DEGENERATE" not in err1, "warned despite an explicit -outward"
    assert "auto orientation is DEGENERATE" not in err2, "warned on a -tie (sign-independent) pair"
