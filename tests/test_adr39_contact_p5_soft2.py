"""ADR-39 P5 / capstone B2 — SOFT=2 segment-based explicit penalty.

The shipped B1 SOFT=1 penalty is NODE-to-segment (NTS): it checks slave NODES against master
SEGMENTS. That MISSES contact in three regimes — (1) two facets overlapping (penetrating) in an
interior strip where NO slave node projects in-bounds of the master, (2) a slave node sliding just
off a master segment edge (the projection goes out-of-bounds ⇒ the pair goes inert), (3) a slave
node at a T-intersection with no unique owning master segment. B2 adds the SOFT=2 SEGMENT-BASED
penalty: `contact ... -mortar -soft <SOFSCL>` detects slave-facet ∩ master-facet overlap (the
shipped clip→Gauss mortar kernel) and distributes a Courant-stable penalty traction
(k_soft = SOFSCL·4·m_eff/dt², m_eff the segment-based gap-mode mass) over the clipped overlap — so
explicit contact survives corners/edges/recontact at the STRUCTURAL timestep (the progressive-
collapse / element-removal enabler). Derivation + FD/energy/stability checks: the numpy oracle
contact_prototypes/proto_b2_soft2_segment.py.

Gates (capstone B2 row):
  (1) CORNER/EDGE ROBUSTNESS — a partial-overlap case the NTS node-to-segment lane gets WRONG
      (no slave node in-bounds ⇒ pass-through) is CAUGHT by SOFT=2 (the slave is restrained);
  (2) EXPLICIT STABILITY + ENERGY BALANCE — a stiff segment penalty that DIVERGES at the structural
      dt runs STABLY with -soft (ω·dt = 2√SOFSCL); a free impact rebounds energy-balanced;
  (3) SOFT=2 EXPLICIT-ONLY — under an implicit integrator -soft is inert (the configured epsN is
      used ⇒ byte-identical to a plain -mortar penalty); -soft-ABSENT mortar is unchanged;
  (5) command-surface refusals (-mortar -soft is frictionless MVP: ⊥ -tie/-mu/-cohesion/-tauMax/
      -visc; needs a base penalty; SOFT=1 vs SOFT=2 selected by the formulation).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

M = 1.0
DT = 1.0e-3          # the "structural" step; ω(kn=1e8)·dt ≫ 2 ⇒ a stiff penalty is unstable here
KN_STIFF = 1.0e8     # a stiff (accuracy-driven) penalty that throttles the explicit dt
RUNAWAY = 1.0


def _fixed_master():
    """A FIXED unit master quad [0,1]² at z=0 (nodes 1..4), outward +z."""
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


# ===================================================================================
# GATE 1 — CORNER/EDGE ROBUSTNESS: a partial-overlap case NTS misses, SOFT=2 catches.
# A LARGE slave facet centred on the small master so EVERY slave node projects OUTSIDE
# the master bounds (NTS node-to-segment ⇒ zero force ⇒ pass-through) while the facets
# DO overlap+penetrate over the whole master (SOFT=2 segment-based ⇒ restrained).
# ===================================================================================
def _big_slave_nodes(z):
    """4 corner nodes of a [-0.5,1.5]² slave facet at height z (all OUTSIDE master [0,1]²)."""
    for t, (x, y) in zip((5, 6, 7, 8), [(-0.5, -0.5), (1.5, -0.5), (1.5, 1.5), (-0.5, 1.5)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)


def _run_explicit(v=-1.0, z0=-0.005, nsteps=1500, soft=0.0, nts=False, kn=1.0e6):
    """Drive the slave facet DOWN into the master at velocity v (starts z0 penetrating), explicit
    CentralDifferenceLadruno. nts=True ⇒ the NTS node-to-segment lane (slave NODES); else SOFT=2
    (-mortar -soft). Returns (min_z, last_z) of node 5 — min_z is the deepest penetration reached."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _big_slave_nodes(z0)
    if nts:
        ops.contactSurface(2, "-slave", 5, 6, 7, 8)
        ops.contact(1, 1, 2, kn, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, v, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    min_z, last_z = z0, z0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            break
        last_z = ops.nodeDisp(5, 3) + z0
        min_z = min(min_z, last_z)
    return min_z, last_z


def test_soft2_catches_edge_overlap_that_nts_misses():
    """GATE 1 — the headline. A slave facet whose every node is OUTSIDE the master footprint but
    whose interior OVERLAPS the master: NTS (node-to-segment) finds no in-bounds slave node ⇒ zero
    force ⇒ the facet falls THROUGH the master (deep, growing penetration). SOFT=2 detects the
    segment overlap ⇒ the facet is restrained/bounced (bounded penetration). Same geometry + dt."""
    nts_min, nts_last = _run_explicit(nts=True)
    s2_min, s2_last = _run_explicit(soft=0.1)
    # NTS passes through: node 5 descends far past the master (no contact restraint).
    assert nts_min < -0.5, f"NTS should pass through the un-noded overlap (min z={nts_min:.3e})"
    # SOFT=2 catches it: penetration stays bounded and the facet rebounds (last z above the deepest).
    assert s2_min > -0.05, f"SOFT=2 should restrain the edge overlap (min z={s2_min:.3e})"
    assert s2_last > s2_min + 1e-4, f"SOFT=2 should rebound (last {s2_last:.3e} > min {s2_min:.3e})"
    assert s2_min - nts_min > 0.4, "SOFT=2 must catch contact the NTS lane misses (much less pen.)"


# ===================================================================================
# GATE 2 — EXPLICIT STABILITY + ENERGY BALANCE (matched footprint ⇒ all nodes supported)
# ===================================================================================
def _matched_slave_nodes(z):
    """4 nodes of a slave facet on the SAME [0,1]² footprint as the master, at height z."""
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)


def _matched_loaded(kn, soft, dt=DT, P=100.0, v=-0.5, nsteps=3000):
    """A matched slave facet pressed onto the fixed master by a sustained load + a velocity kick,
    explicit. Returns (diverged, max|z|). The stiff penalty (ω·dt ≫ 2) diverges; SOFT=2 is bounded."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave_nodes(0.0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if soft > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -P)
        ops.setNodeVel(t, 3, v, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs
        maxabs = max(maxabs, abs(ops.nodeDisp(5, 3)))
        if maxabs > RUNAWAY:
            return True, maxabs
    return False, maxabs


def test_soft2_stabilizes_explicit_segment_contact():
    """GATE 2a — a stiff segment penalty (ω·dt ≫ 2) DIVERGES under central difference at the
    structural dt; the SAME base epsN with `-soft` sizes k_soft = SOFSCL·4·m_eff/dt²
    (ω·dt = 2√SOFSCL < 2) and stays BOUNDED. dt identical across both — only SOFT differs."""
    div_stiff, _ = _matched_loaded(KN_STIFF, 0.0)
    div_soft, maxabs = _matched_loaded(KN_STIFF, 0.1)
    assert div_stiff, "the stiff segment penalty should DIVERGE at the structural dt (ω·dt ≫ 2)"
    assert not div_soft, f"SOFT=2 must be STABLE at the same dt (max|z|={maxabs:.3e})"
    assert maxabs < 0.05, f"SOFT=2 response should be bounded near equilibrium (max|z|={maxabs:.3e})"


def _matched_impact(kn, soft, dt=DT, v_mag=2.0, gap0=0.02, nsteps=3000):
    """A free matched slave facet gap0 above the master, flung down at v_mag (rebounds + flies off).
    Returns (blew, restitution): a conservative contact gives e ≈ 1; a stiff penalty INJECTS (e ≫ 1)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave_nodes(gap0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if soft > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, -v_mag, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    v_out, blew = 0.0, False
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            blew = True
            break
        z = gap0 + ops.nodeDisp(5, 3)
        if z > 0.5 * gap0 and ops.nodeVel(5, 3) > 0.0:
            v_out = ops.nodeVel(5, 3)
    return blew, abs(v_out / v_mag)


def test_soft2_energy_balance_restitution():
    """GATE 2b — a conservative segment penalty must not inject/leak energy. The SOFT=2 impact
    rebounds with restitution ≈ 1 (energy-balanced); the equivalent stiff penalty at the same dt
    INJECTS energy (e ≫ 1) or diverges — the instability SOFT=2 removes."""
    _, e_soft = _matched_impact(KN_STIFF, 0.1)
    blew_stiff, e_stiff = _matched_impact(KN_STIFF, 0.0)
    assert abs(1.0 - e_soft) < 0.05, f"SOFT=2 not energy-conserving (e={e_soft:.4f}, expected ≈1)"
    assert blew_stiff or e_stiff > 2.0, f"stiff penalty should inject/diverge (e={e_stiff:.2f})"


# ===================================================================================
# GATE 3 — SOFT=2 is EXPLICIT-ONLY: implicit is inert (configured epsN ⇒ byte-identical)
# ===================================================================================
def _matched_static(kn, soft):
    """The matched slave facet pressed onto the fixed master by a static load (implicit Newton).
    SOFT=2 is explicit-only ⇒ under implicit it falls through to the shipped mortar penalty with the
    configured epsN. Returns node-5 penetration."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave_nodes(0.0)
    # ADR-76 LAPACK-fix fallout: each slave z is held by NOTHING until the
    # mortar contact engages, so the FIRST tangent is singular. Pre-fix,
    # FullGeneral swallowed that solve (rc 0, X = B) and the garbage first
    # "displacement" is what seated the contact. Ground each slave z with a
    # spring far softer than epsN (0.1 vs 1e6); both compared runs get the
    # identical springs, so the bit-identity assertion is untouched.
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(90 + t, float(x), float(y), 0.0)
        ops.fix(90 + t, 1, 1, 1)
        ops.element("zeroLength", 900 + t, 90 + t, t, "-mat", 90, "-dir", 3)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if soft > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -100.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "static mortar contact did not converge"
    return ops.nodeDisp(5, 3)


def test_soft2_implicit_byte_identical():
    """GATE 3 — SOFT=2 is EXPLICIT-ONLY. Under an implicit integrator the dynamic_cast to
    CentralDifferenceLadruno fails ⇒ -soft is inert and the shipped mortar penalty (configured epsN)
    runs ⇒ a static solve with `-soft` is BIT-IDENTICAL to one without it."""
    KN = 1.0e6
    z_off = _matched_static(KN, 0.0)
    z_on = _matched_static(KN, 0.1)
    assert z_off == z_on, f"-soft perturbed the implicit (explicit-only) solution: {z_off!r} vs {z_on!r}"
    assert z_off < 0.0, f"the static mortar contact should penetrate under load (z={z_off:.3e})"
    # Magnitude gate (guards the grounding spring added for the ADR-76 LAPACK
    # fix): an ENGAGED mortar contact penetrates ~P/epsN ~ 1e-4; if the broad
    # phase ever failed to find the pair after the big spring-only first
    # iterate, the rig would converge on the 0.1 spring alone at z ~ -1e3 and
    # the two asserts above would pass MEANINGLESSLY. 7 orders separate the
    # branches; -1.0 cleanly splits them.
    assert z_off > -1.0, (
        f"penetration magnitude {z_off:.3e} says the mortar contact never "
        "engaged — the rig converged on the grounding spring alone")


# ===================================================================================
# GATE 5 — command-surface refusals + the SOFT=1/SOFT=2 selector
# ===================================================================================
def test_soft2_command_refusals():
    """GATE 5 — `-mortar -soft` is the segment-based explicit penalty: only `-tie` is refused (a
    permanent bond is not a soft contact penalty), and it needs a base penalty. `-visc` (D2.2 normal
    damper) AND `-mu/-cohesion/-tauMax` (Courant-stable Coulomb/Tresca on the SOFT=2 overlap) are now
    allowed. SOFT=1 (NTS) vs SOFT=2 (mortar) is selected by the formulation. OpenSeesPy returns None on
    success, raises on a parse error."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(4 + t, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    # -mortar -soft + -tie ⇒ refused (a permanent bond is not a soft contact penalty)
    with pytest.raises(Exception):
        ops.contact(2, 1, 2, "-mortar", "-soft", 0.1, "-tie", "-epsTie", 1.0e6)
    # -mortar -soft with NO base penalty ⇒ refused
    with pytest.raises(Exception):
        ops.contact(4, 1, 2, "-mortar", "-soft", 0.1)
    # -mortar -soft + -visc ⇒ ACCEPTED (the D2.2 viscous damper on the SOFT=2 active set)
    ops.contact(3, 1, 2, "-mortar", "-epsN", 1.0e6, "-soft", 0.1, "-visc", 1.0,
                "-outward", 0.0, 0.0, 1.0)
    # -mortar -soft + friction ⇒ ACCEPTED (Courant-stable Coulomb/Tresca on the SOFT=2 overlap)
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-soft", 0.1, "-mu", 0.3,
                "-outward", 0.0, 0.0, 1.0)
    # -mortar -soft WITH a base penalty (-epsN auto) ⇒ accepted (SOFT=2)
    ops.contact(5, 1, 2, "-mortar", "-epsN", "auto", "-soft", 0.1, "-outward", 0.0, 0.0, 1.0)


def test_soft2_default_sofscl_accepted():
    """A bare `-mortar ... -soft` (no value) defaults to SOFSCL=0.10 (the LS-DYNA SLSFAC default)
    and runs stably through the real contact command — the value-less flag is a valid form."""
    div, maxabs = None, None
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave_nodes(0.0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", KN_STIFF, "-soft", "-outward", 0.0, 0.0, 1.0)  # bare -soft
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -100.0)
        ops.setNodeVel(t, 3, -0.5, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(2000):
        assert ops.analyze(1, DT) == 0, "bare -soft (default SOFSCL) should be stable"
        assert abs(ops.nodeDisp(5, 3)) < RUNAWAY, "bare -soft should stay bounded"
