"""ADR-57 E7 — edge-edge INTEGRATION + regression (the trust gate).

E0–E6 proved the edge-edge mechanics on near-unit geometries (bare crossed-bar facets + artificial
regularizing springs, numpy oracles). E7 is the first time the lane meets REAL multi-element bodies:
the slave is a deformable `LadrunoBrick` whose mass (material ρ) and stiffness (E) come from the
element — so the point-pair differential mode is regularized by the body itself (NO artificial
springs), and the explicit SOFT m_eff is sized from the ASSEMBLED element-density mass (the full
ladrunoBuildNodalMass path), not a hand-placed nodal mass.

Geometry (cross-stacked bars, the real-body version of test_adr57_edge_edge_1): a fixed master facet
in the y-z plane (top edge ∥y at z=0) crossed perpendicular by a slave BRICK bar whose bottom face
(x-z plane, y=0) has its bottom edge ∥x at z=ZC. The two edges cross over the origin, n=ẑ — exactly
the cos_t→0 case face-mortar AND NTS go blind on, now between two MESHED bodies.

Gates (the ADR-57 E7 row):
  - REAL-BODY restraint: the slave brick, pressed down, is held by the edge-edge contact (Newton
    converges with the brick supplying the regularization — NO springs);
  - the NTS falsifier survives the move to real bodies (NTS transmits ZERO force ⇒ the brick passes
    through, edge-edge restrains);
  - E5 SOFT on REAL element mass: an explicit drop runs at the structural dt (the configured stiff εN
    would throttle dt_cr), m_eff from the assembled element-density mass;
  - E6 ALM on a REAL elastic body: held-load augmentation closes the edge gap to augTol εN-independently;
  - byte-identity at integration scale: bars far apart (no pair routes) ⇒ -edgeedge present == absent;
  - the reference-config / large-sliding LIMIT, documented: a pair driven to slide off its crossing
    releases gracefully (force→0, no crash) — the §7 scope fence + the E6 stale-gap fold, validated.
"""
import math
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

E_MOD = 2.0e4    # brick Young's modulus
NU = 0.0         # uniaxial (no Poisson coupling — the slave nodes are x,y fixed)
W = 0.4          # brick y-extent (the bar's out-of-crossing-plane thickness)
H = 0.4          # brick z-height
EPS = 1.0e6      # edge-edge normal penalty εN


def _master_facet():
    """Fixed master quad in the y-z plane at x=0; top edge m1-m2 (1,2) ∥ y at z=0."""
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _slave_brick(zc, rho=0.0, xoff=0.0, matTag=1, fix=(1, 1, 0)):
    """A deformable LadrunoBrick bar crossing the master perpendicular. Bottom z-face at z=zc, top at
    z=zc+H, spanning x∈[-1,1]+xoff (the crossing axis) and y∈[0,W]. The y=0 face (x-z plane) is the
    slave contact facet; its bottom edge (5,6 ∥x at z=zc) crosses the master edge (1,2 ∥y at z=0).
    `fix` is the per-node SP pattern (default (1,1,0) ⇒ x,y fixed, z free — uniaxial; the brick couples
    the bottom-edge z-DOFs ⇒ regularizes the point-pair differential mode, NO springs). Returns
    (contact_facet, top_face_nodes)."""
    n = [(5, -1.0 + xoff, 0.0, zc), (6, 1.0 + xoff, 0.0, zc), (7, 1.0 + xoff, W, zc), (8, -1.0 + xoff, W, zc),
         (15, -1.0 + xoff, 0.0, zc + H), (16, 1.0 + xoff, 0.0, zc + H),
         (17, 1.0 + xoff, W, zc + H), (18, -1.0 + xoff, W, zc + H)]
    for t, x, y, z in n:
        ops.node(t, x, y, z)
        ops.fix(t, *fix)
    ops.nDMaterial("ElasticIsotropic", matTag, E_MOD, NU, rho)
    ops.element("LadrunoBrick", 1, 5, 6, 7, 8, 15, 16, 17, 18, matTag)
    contact_facet = [5, 6, 16, 15]               # y=0 face quad (x-z plane); bottom edge 5-6 ∥ x
    top_face = [15, 16, 17, 18]                  # z=zc+H face (load here to press the bar down)
    return contact_facet, top_face


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


# ------------------------------------------------------------------ real-body normal restraint
def test_e7_cross_bars_real_brick_restrained():
    """REAL-BODY headline: a deformable brick bar, pressed down, is held by the edge-edge contact at
    the penalty depth δ≈P/εN — and Newton CONVERGES with the BRICK (not a spring) regularizing the
    point-pair differential mode. The bottom edge does not pass through the master edge."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_facet()
    facet, top = _slave_brick(zc=-1.0e-4)        # seeded just penetrating ⇒ contact active step 1
    ops.contactSurface(2, "-slave-segments", 4, *facet)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 200.0
    for t in top:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    _static()
    assert ops.analyze(1) == 0, "real-brick edge-edge Newton did not converge (brick regularization?)"
    z_edge = -1.0e-4 + ops.nodeDisp(5, 3)        # absolute z of the bottom edge
    # MECHANISM-SPECIFIC: by static equilibrium of the brick (top load P down = bottom-edge contact
    # reaction P up), the edge-edge penalty depth is δ = P/εN — INDEPENDENT of the brick stiffness (the
    # brick only sets the top displacement). So z_edge = −P/εN pins the EDGE-EDGE restraint specifically,
    # not "small from the brick." Without edge-edge the bottom edge has NO z-support (a rigid z-mode) ⇒
    # singular / unbounded, so this exact-depth match cannot pass on the brick alone.
    delta = P / EPS
    assert z_edge < 0.0, "bottom edge did not stay in contact"
    assert abs(z_edge - (-delta)) < 0.2 * delta, \
        f"bottom edge not at the edge-edge penalty depth −P/εN={-delta:.2e} (got {z_edge:.2e})"


def test_e7_real_brick_nts_passes_through():
    """THE FALSIFIER survives real bodies: the SAME brick bar with an NTS contact (its bottom-edge
    NODES vs the master facet) — the nodes are 1 unit off the master face in x, so NONE projects
    in-bounds ⇒ NTS transmits ZERO force. A soft z-spring per bottom node makes the problem well-posed
    once the contact contributes nothing (else the rigid z-mode is singular); the brick then sinks to
    the pure-spring depth — ~10³× past the edge-edge restraint. NTS passes through where edge-edge holds."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_facet()
    facet, top = _slave_brick(zc=-1.0e-4)
    ks = 1.0e3
    ops.uniaxialMaterial("Elastic", 9, ks)
    for gnd, sn, x in [(905, 5, -1.0), (906, 6, 1.0)]:
        ops.node(gnd, x, 0.0, -1.0e-4)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", gnd, gnd, sn, "-mat", 9, "-dir", 3)
    ops.contactSurface(3, "-slave", 5, 6)        # NTS slave NODE-set = the bottom-edge nodes
    ops.contact(1, 1, 3, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 200.0
    for t in top:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    _static()
    assert ops.analyze(1) == 0, "NTS real-brick Newton did not converge"
    z_edge = -1.0e-4 + ops.nodeDisp(5, 3)
    # NTS net force ZERO ⇒ the load is carried by the soft springs alone ⇒ the bottom edge sinks to the
    # pure-spring depth ~P/(2ks) ≈ 0.1 (P=200, 2ks=2e3), ~10³× past the edge-edge restraint (P/εN=2e-4).
    spring_depth = P / (2.0 * ks)
    assert z_edge < -0.5 * spring_depth, f"NTS unexpectedly restrained the edge (z_edge={z_edge})"
    assert abs(z_edge) > 100.0 * (P / EPS), "NTS did not pass through (it restrained like edge-edge?)"


# ------------------------------------------------------------------ E5 SOFT on real element mass
def test_e7_soft_real_element_mass():
    """E5 SOFT with REAL element-density mass: a brick bar (material ρ>0) dropped onto the master edge
    under explicit central difference. With -edgeSoft the contact mode ω·dt=2√SOFSCL≤2 ⇒ stable at the
    structural dt (the brick's own dt_cr); the configured stiff εN would throttle dt_cr and DIVERGE.
    m_eff is sized from the ASSEMBLED element-density mass (ladrunoBuildNodalMass element pass)."""
    def _drop(soft, eps=1.0e11, dt=2.0e-4, nsteps=40, v0=-2.0):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _master_facet()
        facet, top = _slave_brick(zc=0.01, rho=10.0)   # ABOVE the master; real element mass (ρ=10)
        ops.contactSurface(2, "-slave-segments", 4, *facet)
        args = [1, 1, 2, "-mortar", "-epsN", eps, "-edgeedge", "-edgeBand", 0.15,
                "-outward", 0.0, 0.0, 1.0]
        if soft:
            args += ["-edgeSoft", 0.1]
        ops.contact(*args)
        for t in (5, 6, 7, 8, 15, 16, 17, 18):
            ops.setNodeVel(t, 3, v0, "-commit")        # the whole bar approaches downward
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.integrator("CentralDifferenceLadruno")
        ops.algorithm("Linear")
        ops.analysis("Transient")
        zmin = 0.01
        for _ in range(nsteps):
            ops.analyze(1, dt)
            z = 0.01 + ops.nodeDisp(5, 3)
            if not math.isfinite(z) or abs(z) > 1.0e4:
                return None, zmin                       # diverged
            zmin = min(zmin, z)
        return 0.01 + ops.nodeDisp(5, 3), zmin
    z_soft, pen_soft = _drop(soft=True)
    z_stiff, _ = _drop(soft=False)
    # SOFT: stable ⇒ finite, bounded penetration AND bounded total excursion (the bar engages then is
    # held/rebounds within the drop scale ~0.01, never ejected). m_eff from the assembled element mass.
    assert z_soft is not None, "SOFT drop diverged (should be Courant-stable at the structural dt)"
    assert pen_soft < 0.0, "SOFT bar never engaged the master edge"
    assert pen_soft > -0.05, f"SOFT penetration unbounded ({pen_soft}) — m_eff sizing wrong?"
    assert abs(z_soft) < 0.1, f"SOFT excursion not bounded near the drop scale (z_soft={z_soft})"
    # stiff εN at the SAME dt: ω_contact·dt ≫ 2 ⇒ HARD divergence (None = |z|>1e4/NaN) OR a gross
    # non-physical ejection ≫ the bounded SOFT excursion. Both are the energy-injecting CFL violation.
    assert z_stiff is None or abs(z_stiff) > 100.0 * abs(z_soft), \
        f"stiff εN did NOT diverge vs the bounded SOFT run (stiff={z_stiff}, soft={z_soft})"


# ------------------------------------------------------------------ E6 ALM on a real elastic body
def test_e7_alm_real_elastic_body():
    """E6 ALM on a deformable brick: held-load augmentation closes the edge gap to augTol at FINITE εN,
    εN-INDEPENDENTLY — the elastic coupling the oracle could not model. Penalty alone leaves an
    εN-dependent residual penetration."""
    AUGTOL = 1.0e-9

    def _push(eps):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _master_facet()
        facet, top = _slave_brick(zc=-1.0e-4)
        ops.contactSurface(2, "-slave-segments", 4, *facet)
        ops.contact(1, 1, 2, "-mortar", "-epsN", eps, "-edgeedge", "-edgeBand", 0.15,
                    "-edgeAlm", "-edgeAugTol", AUGTOL, "-outward", 0.0, 0.0, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, 0.0, -200.0 / 4.0)
        _static()
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, "ALM real-brick push did not converge"
        pen_penalty = ops.ladrunoEdgePenetration()     # λ_N=0 on the first solve ⇒ penalty penetration
        status, _, pen = analyze_augmented(ops, maxAug=60, augTol=AUGTOL,
                                           query=ops.ladrunoEdgePenetration)
        assert status == 0, f"εN={eps:.0e}: real-body augmentation failed (pen {pen})"
        return pen_penalty, pen
    lo = _push(1.0e5)
    hi = _push(1.0e7)
    assert lo[0] > 10.0 * hi[0], f"penalty penetration not εN-dependent: {lo[0]} vs {hi[0]}"
    assert lo[1] < AUGTOL and hi[1] < AUGTOL, f"ALM did not close the gap εN-independently: {lo[1]}, {hi[1]}"


# ------------------------------------------------------------------ byte-identity at integration scale
def test_e7_byte_identical_no_route():
    """Byte-identity at integration scale: a brick bar 10 units away in +x (no edge pair routes) evolves
    BITWISE identically with -edgeedge present vs absent — the off-by-default contract on a REAL body."""
    def _run(edgeedge):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _master_facet()
        facet, top = _slave_brick(zc=0.01, rho=5.0, xoff=10.0)   # far away ⇒ no proximity
        ops.contactSurface(2, "-slave-segments", 4, *facet)
        args = [1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0]
        if edgeedge:
            args += ["-edgeedge", "-edgeBand", 0.15]
        ops.contact(*args)
        for t in (5, 6, 7, 8, 15, 16, 17, 18):
            ops.setNodeVel(t, 3, 0.1, "-commit")
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        out = []
        for _ in range(6):
            assert ops.analyze(1, 1.0e-3) == 0
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8, 15, 16, 17, 18)))   # ALL 8 brick nodes
        return out
    assert _run(True) == _run(False)


# ------------------------------------------------------------------ non-axis-aligned crossing (n ≠ ẑ)
def _inclined_master(alpha):
    """The master facet, the y-z quad ROTATED about the x-axis by α (it stays in the x=0 plane ⇒ its
    normal stays ∥x ⇒ edge-edge still routes), so its top edge runs along ê_m=(0,cosα,sinα) and the
    common-perpendicular normal n=ê_s×ê_m=(0,−sinα,cosα) is NON-vertical. Stresses the kernel normal,
    the body-fixed sign, the B-operator and routing on a non-canonical orientation."""
    ca, sa = math.cos(alpha), math.sin(alpha)
    for t, y0, z0 in [(1, -1.0, 0.0), (2, 1.0, 0.0), (3, 1.0, -0.3), (4, -1.0, -0.3)]:
        ops.node(t, 0.0, y0 * ca - z0 * sa, y0 * sa + z0 * ca)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def test_e7_non_axis_aligned_crossing():
    """A genuinely off-axis crossing (n=(0,−sin30°,cos30°), NOT ẑ) on a real brick: the lane must route
    edge-edge (facets ⊥), restrain the bar along the inclined normal, and the E6 ALM must drive the
    penetration → augTol — exercising the normal/sign/B-operator/routing away from the canonical
    axis-aligned case the other tests use. Guards against an n-orientation or sign-anchor artifact."""
    AUGTOL = 1.0e-9
    alpha = math.radians(30.0)

    def _push(eps):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _inclined_master(alpha)
        facet, top = _slave_brick(zc=-1.0e-4)        # bottom edge ∥x; crosses the inclined master edge
        ops.contactSurface(2, "-slave-segments", 4, *facet)
        ops.contact(1, 1, 2, "-mortar", "-epsN", eps, "-edgeedge", "-edgeBand", 0.15,
                    "-edgeAlm", "-edgeAugTol", AUGTOL, "-outward", 0.0, 0.0, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, 0.0, -200.0 / 4.0)      # push down; n has a z-component cosα ⇒ engages
        _static()
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, "inclined-crossing Newton did not converge (off-axis normal?)"
        pen_penalty = ops.ladrunoEdgePenetration()   # >0 ⇒ the off-axis pair routed + is KKT-active
        status, _, pen = analyze_augmented(ops, maxAug=60, augTol=AUGTOL,
                                           query=ops.ladrunoEdgePenetration)
        assert status == 0, f"εN={eps:.0e}: off-axis augmentation failed (pen {pen})"
        return pen_penalty, pen
    lo = _push(1.0e5)
    hi = _push(1.0e7)
    # the off-axis pair genuinely routed + engaged (a non-zero, εN-dependent penalty penetration)…
    assert lo[0] > 0.0 and hi[0] > 0.0, "the off-axis crossing did not route edge-edge / engage"
    assert lo[0] > 10.0 * hi[0], f"off-axis penalty penetration not εN-dependent: {lo[0]} vs {hi[0]}"
    # …and the ALM closed it to augTol εN-independently on the inclined normal.
    assert lo[1] < AUGTOL and hi[1] < AUGTOL, f"off-axis ALM did not close the gap: {lo[1]}, {hi[1]}"


# ------------------------------------------------------------------ E3 friction on a real body
def test_e7_friction_real_brick():
    """E3 Coulomb/Tresca friction on a REAL deformable brick (the oracle + test_2 used bare facets): the
    brick is held at a fixed penetration (z fixed ⇒ a CONSTANT normal pressure N=εN|ZC|) and dragged
    tangentially in y by a SUB-CONE load. A soft y-spring per bottom node makes the no-friction baseline
    well-posed (the point-pair shares ONE tangential gap ⇒ the differential y-mode is a zero-energy mode
    of the point contact — inherent to point-like penalty contact, as in E2/E3); the BRICK supplies the
    real-body coupling that turns the bottom-edge y into the slip the edgeSlip measure reads. With the
    cone (μN+c) ABOVE the drag the bar STICKS (small y, Newton converges); with μ=0,c=0 it slides to the
    spring limit. The contrast is the assembled friction traction on a meshed body."""
    ZCF = -1.0e-3                                              # held penetration ⇒ N = εN·|ZCF|
    KY = 100.0                                                 # soft y-spring (the μ=0 baseline)

    def _stick(mu, coh):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _master_facet()
        facet, top = _slave_brick(zc=ZCF, fix=(1, 0, 1))      # x,z fixed (held penetration); y = the rake DOF
        ops.contactSurface(2, "-slave-segments", 4, *facet)
        ops.uniaxialMaterial("Elastic", 9, KY)
        for gnd, sn, x in [(905, 5, -1.0), (906, 6, 1.0)]:
            ops.node(gnd, x, 0.0, ZCF)
            ops.fix(gnd, 1, 1, 1)
            ops.element("zeroLength", gnd, gnd, sn, "-mat", 9, "-dir", 2)
        args = [1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
                "-outward", 0.0, 0.0, 1.0]
        if mu > 0.0 or coh > 0.0:
            args += ["-edgeMu", mu, "-edgeCohesion", coh, "-edgeKt", 1.0e4]
        ops.contact(*args)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in (5, 6):
            ops.load(t, 0.0, 5.0, 0.0)                        # +y drag, total 10 (sub-cone if friction on)
        _static()
        assert ops.analyze(1) == 0, f"friction real-brick (μ={mu}) Newton did not converge"
        return ops.nodeDisp(5, 2)                              # the y (slip) displacement of a bottom node
    # cone = μN + c = 0.5·(εN·1e-3) + 50 = 500+50 = 550 ≫ the 10 N drag ⇒ STICK; μ=0,c=0 ⇒ slide.
    y_fric = _stick(0.5, 50.0)
    y_free = _stick(0.0, 0.0)
    # μ=0: contact gives NO tangential resistance ⇒ the bar drags to the spring limit (P/2KY=10/200=0.05).
    # friction on: the cone holds the drag ⇒ STICK ⇒ the y-slip is far smaller (kt-elastic, ≪ the slide).
    assert abs(y_free) > 0.02, f"the μ=0 baseline did not slide (y={y_free}); cannot show the contrast"
    assert abs(y_fric) < 0.2 * abs(y_free), \
        f"Coulomb friction did not resist the drag on the real brick (μ>0 y={y_fric}, μ=0 y={y_free})"


# ------------------------------------------------------------------ the large-sliding limit (documented)
def test_e7_large_sliding_releases_gracefully():
    """The reference-config / large-sliding LIMIT (§7 scope fence), validated as a GRACEFUL release, not
    a crash. The edge pair is enumerated once at the reference config (handle()-time). Under explicit
    central difference the slave bar TRANSLATES along its own axis (+x) at constant velocity; its bottom
    edge slides until the closest point leaves the margin-interior band ⇒ edgeGeom returns false ⇒ the
    adapter goes inert (force→0), and the E6 stale-gap fold (the reset-on-inactive + edgeGCEnd reset)
    keeps gN_committed/εN from corrupting the query. The contact does NOT track the large slide (it needs
    re-emission — the documented hand-off), but it releases CLEANLY: the run stays FINITE (a stale/
    degenerate pair would inject energy and blow up) and the penetration query drops to 0. This pins the
    boundary of the supported regime — the honest E7 finding."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_facet()
    # x FREE (the slide axis), y fixed, z free with real mass; -edgeSoft keeps the z-holding Courant-
    # stable while sliding, -edgeAlm makes the penetration query + the stale-gap reset live.
    facet, top = _slave_brick(zc=-1.0e-4, rho=10.0, fix=(0, 1, 0))
    ops.contactSurface(2, "-slave-segments", 4, *facet)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
                "-edgeSoft", 0.1, "-edgeAlm", "-edgeAugTol", 1.0e-9, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8, 15, 16, 17, 18):
        ops.setNodeVel(t, 1, 500.0, "-commit")       # rigid +x translation (no internal strain ⇒ no waves)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 2.0e-4
    pen_early = 0.0
    for k in range(60):
        ops.analyze(1, dt)                           # x advances ~0.1/step ⇒ slides off the crossing (x>1)
        z = -1.0e-4 + ops.nodeDisp(5, 3)
        assert math.isfinite(z) and abs(z) < 1.0e3, f"sliding-off injected energy / blew up (z={z})"
        if k < 5:
            pen_early = max(pen_early, ops.ladrunoEdgePenetration())   # while still crossing (x<1)
    # ENGAGED → RELEASED (not "never crossed"): the pair was a live, KKT-active edge-edge contact early
    # on (seeded penetrating, still margin-interior ⇒ query > 0), THEN slid past its crossing.
    assert pen_early > 0.0, "the pair never engaged (the slide-off path is not actually exercised)"
    # the bottom edge has slid well past the crossing (x ≫ 1) ⇒ the pair is off its margin-interior band
    x5 = ops.nodeCoord(5, 1) + ops.nodeDisp(5, 1)
    assert x5 > 2.0, f"the bar did not slide off the crossing (x5={x5})"
    # GRACEFUL: the edge-edge query now reports NO active penetration (the adapter went inert as the
    # closest point left the margin-interior band + the stale-gap reset cleared gN_committed/εN ⇒
    # getMaxEdgePenetration skips the dead slot). Engaged → released, no crash, no stale-state injection.
    assert ops.ladrunoEdgePenetration() < 1.0e-9, \
        "edge-edge still reports active penetration after sliding off (stale gap not released)"
