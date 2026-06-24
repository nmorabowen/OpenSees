"""ADR-39 P4 / capstone B1 — SOFT=1 Courant-stable contact penalty.

The shipped NTS penalty sizes kn from MATERIAL/element stiffness (fixed -kn or -kn auto). Under the
explicit CentralDifferenceLadruno the contact spring kn lives in the RESIDUAL and adds to the
effective stiffness, so a stiff kn raises the contact mode ω = √(kn/m_eff) and THROTTLES the
central-difference critical step dt_cr = 2/ω: a stiff/auto kn ⇒ a tiny stable dt ⇒ explicit contact
is impractically slow (and a run at the structural dt DIVERGES — the CDL CFL estimate eigensolves
only the structural elements, not the contact spring).

B1 adds the LS-DYNA §26.15 SOFT=1 option: `-soft <SOFSCL>` sizes the contact stiffness from the
nodal MASS and the current timestep instead — k_soft = SOFSCL·4·m_eff/dt² (m_eff = the gap-mode
effective contact mass; SOFSCL default 0.1) — so the contact's own ω·dt = 2√SOFSCL ≤ 2 and it never
shrinks dt_cr below the structural step. Derivation + FD/energy checks: the numpy oracle
contact_prototypes/proto_b1_soft_penalty.py.

Gates (capstone B1 row):
  (1) EXPLICIT STABILITY — a stiff kn that DIVERGES at the structural dt runs STABLY (bounded) with -soft;
  (2) ENERGY BALANCE — the soft-penalty impact conserves energy (restitution ≈ 1; stiff injects, e ≫ 1);
  (3) -soft ABSENT byte-identical; IMPLICIT byte-identical (SOFT is explicit-only);
  (4) the SEGMENT (faceted NTS) path is stabilized too;
  (5) command-surface refusals (-soft needs a base kn; -soft ⊥ -mortar).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

M = 1.0
DT = 1.0e-3          # the "structural" step; ω(kn=1e8)·dt = 10 ≫ 2 ⇒ a stiff penalty is unstable here
KN_STIFF = 1.0e8     # a stiff (accuracy-driven) penalty that throttles the explicit dt
RUNAWAY = 1.0        # |z| beyond this (≫ any physical penetration) ⇒ the integration diverged


def _plane_loaded(kn, soft, dt=DT, P=1.0e3, v=-0.5, nsteps=3000):
    """A mass at the rigid plane (normal +z), pressed by a SUSTAINED downward load P (so the contact
    stays engaged ⇒ the unstable contact mode compounds) and kicked with velocity v. Explicit
    CentralDifferenceLadruno. Returns (diverged, max|z|): diverged = a NaN abort OR a runaway |z|."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)                            # starts AT the plane (immediately in contact)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    if soft > 0.0:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn, "-soft", soft)
    else:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)                             # holds it in contact
    ops.setNodeVel(1, 3, v, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs                          # NaN/Inf abort
        maxabs = max(maxabs, abs(ops.nodeDisp(1, 3)))
        if maxabs > RUNAWAY:
            return True, maxabs                          # runaway ⇒ unstable
    return False, maxabs


def test_soft_stabilizes_explicit_impact():
    """GATE 1 — the headline. At a structural-scale dt a STIFF penalty (ω·dt ≫ 2) DIVERGES under
    central difference (runaway/NaN); the SAME base kn with `-soft` sizes k_soft = SOFSCL·4·m/dt²
    (ω·dt = 2√SOFSCL < 2) and stays BOUNDED. dt is identical across both runs — only SOFT differs."""
    div_stiff, _ = _plane_loaded(KN_STIFF, 0.0)
    div_soft, maxabs_soft = _plane_loaded(KN_STIFF, 0.1)
    assert div_stiff, "the stiff penalty should DIVERGE at the structural dt (ω·dt ≫ 2)"
    assert not div_soft, f"the SOFT=1 penalty must be STABLE at the same dt (max|z|={maxabs_soft:.3e})"
    # the soft penetration is the bounded equilibrium δ = P/k_soft + a small kick amplitude
    assert maxabs_soft < 0.05, f"soft response should be bounded near the equilibrium (max|z|={maxabs_soft:.3e})"


def _plane_impact(kn, soft, dt=DT, v_mag=2.0, gap0=0.02, nsteps=3000):
    """Free mass gap0 above the rigid plane, flung down at v_mag (no sustained load ⇒ it rebounds and
    flies away). Returns (blew, restitution): a CONSERVATIVE contact gives e ≈ 1, an unstable stiff
    penalty INJECTS energy ⇒ e ≫ 1."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    if soft > 0.0:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn, "-soft", soft)
    else:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn)
    ops.setNodeVel(1, 3, -v_mag, "-commit")
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
        z = gap0 + ops.nodeDisp(1, 3)
        if z > 0.5 * gap0 and ops.nodeVel(1, 3) > 0.0:   # separated, moving away
            v_out = ops.nodeVel(1, 3)
    return blew, abs(v_out / v_mag)


def test_soft_uses_element_assembled_mass():
    """GATE 1 (element-mass form) — the explicit headline. A contact node whose mass comes from
    ELEMENT density (no nodal `mass` command) must still be soft-stabilized. Node::getMass() holds
    ONLY the nodal mass, so the handler assembles nodal + Σ element-diagonal mass and the adapter
    reads that. Here a truss supplies node 1's mass (ρ) with NEGLIGIBLE stiffness; with a stiff base
    kn the soft sizing would see m=0 and fall back to the stiff kn ⇒ DIVERGE — so STAYING BOUNDED
    proves the element-assembled mass is what k_soft = SOFSCL·4·m_eff/dt² uses."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)                            # slave on the plane (NO ops.mass on it)
    ops.node(2, 0.0, 0.0, 1.0)
    ops.fix(2, 1, 1, 1)
    ops.fix(1, 1, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e-3)           # negligible truss stiffness
    ops.element("truss", 1, 1, 2, 1.0, 1, "-rho", 2.0)   # lumped node-1 mass = 0.5·ρ·A·L = 1.0
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, "-soft", 0.1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.setNodeVel(1, 3, -0.5, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    maxabs = 0.0
    for _ in range(3000):
        assert ops.analyze(1, DT) == 0, "element-mass soft contact aborted (mass not assembled?)"
        maxabs = max(maxabs, abs(ops.nodeDisp(1, 3)))
        assert maxabs < RUNAWAY, "soft must size from ELEMENT-assembled mass (else falls back ⇒ diverge)"


def test_soft_energy_balance_restitution():
    """GATE 2 — a conservative penalty must not inject/leak energy. The soft impact rebounds with
    restitution ≈ 1 (energy-balanced), whereas the equivalent STIFF penalty at the same dt INJECTS
    energy (e ≫ 1) — the instability. (The residual |1−e| of the soft case is the bounded discrete
    one-sided-contact engagement error, SOFSCL-convergent in the numpy oracle.)"""
    _, e_soft = _plane_impact(KN_STIFF, 0.1)
    _, e_stiff = _plane_impact(KN_STIFF, 0.0)
    assert abs(1.0 - e_soft) < 0.05, f"soft penalty not energy-conserving (e={e_soft:.4f}, expected ≈1)"
    assert e_stiff > 2.0, f"the stiff penalty should INJECT energy at this dt (e={e_stiff:.2f} ≫ 1)"


def _plane_static(kn, soft):
    """A node pressed into the rigid plane by a static downward load (implicit Newton). SOFT is
    explicit-only ⇒ under implicit it is INERT and the configured kn is used (penetration = −P/kn)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    if soft > 0.0:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn, "-soft", soft)
    else:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "static contact did not converge"
    return ops.nodeDisp(1, 3)


def test_soft_implicit_byte_identical():
    """GATE 3 — SOFT is EXPLICIT-ONLY. Under an implicit integrator the dynamic_cast to
    CentralDifferenceLadruno fails ⇒ the configured kn is used ⇒ a static solve with `-soft` is
    BIT-IDENTICAL to one without it (penetration = −P/kn either way)."""
    KN = 1.0e6
    z_off = _plane_static(KN, 0.0)
    z_on = _plane_static(KN, 0.1)
    assert z_off == z_on, f"-soft perturbed the implicit (explicit-only) solution: {z_off!r} vs {z_on!r}"
    assert abs(z_off - (-1.0e3 / KN)) < 1.0e-9, f"penetration {z_off:.3e} != −P/kn (base kn used)"


def _segment_loaded(kn, soft, dt=DT, P=1.0e3, v=-0.5, nsteps=3000):
    """The SEGMENT (faceted NTS) path: a slave node pressed onto a FIXED master quad by a sustained
    load + a velocity kick. m_eff = m_s (the fixed master is infinite-mass ⇒ 0 inverse-mass). Returns
    (diverged, max|z|). Exercises B = [n | −Nᵢ n], distinct from the rigid-plane n-only path."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)                               # fixed master quad (no mass ⇒ ∞ mass)
    ops.node(5, 0.0, 0.0, 0.0)                            # slave at the quad centre, on the surface
    ops.mass(5, M, M, M)
    ops.fix(5, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, kn, 0.0, 0.0, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, 0.0, -P)
    ops.setNodeVel(5, 3, v, "-commit")
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


def test_soft_segment_path_stabilized():
    """GATE 4 — the faceted SEGMENT lane is soft-stabilized too: a stiff kn diverges at the
    structural dt; -soft (m_eff from the slave + shape-projected, here fixed, master masses) is
    bounded. Confirms the fixed-master ⇒ infinite-mass handling in the gap-mode mass."""
    div_stiff, _ = _segment_loaded(KN_STIFF, 0.0)
    div_soft, maxabs = _segment_loaded(KN_STIFF, 0.1)
    assert div_stiff, "stiff SEGMENT penalty should diverge at the structural dt"
    assert not div_soft, f"SOFT=1 SEGMENT penalty must be stable at the same dt (max|z|={maxabs:.3e})"
    assert maxabs < 0.05, f"soft SEGMENT response should be bounded (max|z|={maxabs:.3e})"


def test_soft_command_refusals():
    """GATE 5 — the command surface: -soft needs a base penalty (a modifier, not standalone), and is
    NTS-only (refused on -mortar). OpenSeesPy returns None on success and raises on a parse error."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(4 + t, float(x), float(y), 0.0)
    ops.node(9, 0.5, 0.5, 0.1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contactSurface(3, "-slave", 9)
    # -soft on a mortar contact ⇒ refused (NTS-only)
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-soft", 0.1)
    # -soft with NO base penalty ⇒ refused (needs a positional auto/kn)
    with pytest.raises(Exception):
        ops.contact(2, 1, 3, "-soft", 0.1)               # no kn, no auto
    # -soft WITH a base penalty (positional kn) ⇒ accepted (must NOT raise)
    ops.contact(3, 1, 3, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-soft", 0.1)


def test_soft_default_sofscl_accepted():
    """A bare `-soft` (no value) defaults to SOFSCL=0.10 (the LS-DYNA SLSFAC default) and runs
    stably through the real contactPlane command — the value-less flag is a valid form."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, "-soft")   # bare ⇒ SOFSCL=0.10
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.setNodeVel(1, 3, -0.5, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(2000):
        assert ops.analyze(1, DT) == 0, "bare -soft (default SOFSCL) should be stable"
        assert abs(ops.nodeDisp(1, 3)) < RUNAWAY, "bare -soft should stay bounded"
