"""LadrunoJ2 + Lemaitre damage — meshed notched-bar fracture (Zone-B, gmsh).

The single-element batteries (test_lemaitre_damage.py, test_lemaitre_validation.py)
prove the local constitutive law and the lch handshake by element-size scaling. This
case confirms the headline property in a *full structural fracture* on a gmsh-meshed
double-edge-notched (DEN) flat bar pulled in tension — the classic simple-geometry
Lemaitre ductile-fracture coupon (the textbook V6/V7 cells of ADR §6,
Ladruno_implementation/15_lemaitre_ductile_damage_adr.md).

Mesh-objectivity / fracture energy:  with `-autoRegularization` the dissipated energy
regularizes by the element characteristic length, so the structural force-elongation
response and the total dissipated energy are MESH-INDEPENDENT even though the damage
localizes into a one-element-wide band across the notch ligament. Two meshes (coarse /
fine, 2x refinement) must give the same peak load and the same dissipated energy. A
broken lch handshake (energy proportional to band volume) would make the fine mesh
markedly more brittle and the curves diverge.

The notch ligament localizes the band at a known place (so the test is robust); the
material runs with IMPL-EX for an SPD tangent through the softening branch.

Honest scope: lch regularizes the dissipated ENERGY, not the localization WIDTH (ADR
D4). The fully triaxial round (Bridgman) notched bar — whose necking triaxiality field
needs a graded 3D mesh — is deferred to the apeGmsh validation campaign (see the report
Ladruno_implementation/16_lemaitre_validation.md).

Zone-B: needs gmsh (raw gmsh API, not apeGmsh); runs nightly.
"""
import pytest

from _testbed import ops

gmsh = pytest.importorskip("gmsh")
pytestmark = [pytest.mark.zone_b]

# --- DEN flat bar geometry (N, mm) ---
LX, HY, TZ = 30.0, 20.0, 4.0       # length, height, thickness (one z-layer)
DN, WN = 5.0, 4.0                  # notch depth (each side), notch half-width in x
LIG = HY - 2.0 * DN               # ligament height = 10

# --- steel: J2 + Voce/linear hardening, Lemaitre ductile damage ---
E, NU, SY = 200000.0, 0.3, 400.0
K, G = E / (3.0 * (1.0 - 2.0 * NU)), E / (2.0 * (1.0 + NU))
ISO = ("-iso", "voce", SY, 80.0, 40.0, 50.0)        # s0,Qinf,b,Hiso (gentle => localizes)
DMG = (0.35, 1.0, 0.01, 0.92)                       # r,s,pD,Dc — drives ligament damage
# FIXED reference char length (NOT the element size): lchScale = lch/LCH_REF must
# DIFFER between meshes for the crack-band regularization to equalize the dissipated
# energy. Passing the mesh size here would make lchScale==1 on both and defeat it.
# Chosen >= the coarse element size so lchScale <= 1 on both meshes: a steeper-than-
# natural local softening (lchScale>1) snap-backs under displacement control and stalls
# the coarse mesh before the bands can be compared. Objectivity holds for any fixed
# LCH_REF; this value just keeps both meshes numerically able to traverse the softening.
LCH_REF = 5.0
# Implicit (NOT IMPL-EX) for this structural run: IMPL-EX's bounded-lag extrapolation
# overshoots and destabilizes the fine mesh well before damage develops; the implicit
# return map + adaptive step-cut traverses the softening branch cleanly on both meshes.
# (IMPL-EX robustness itself is validated at the material point in test_lemaitre_damage.py V9.)


def _mesh(h):
    """Structured quad mesh of the LX x HY face (size ~h), recombined + extruded one
    layer in z -> hexes. Returns (nodes{tag:(x,y,z)}, hexes[[8 tags]])."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("denbar")
        occ = gmsh.model.occ
        occ.addRectangle(0, 0, 0, LX, HY)
        occ.synchronize()
        nx, ny = max(4, round(LX / h)), max(4, round(HY / h))
        for c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c[1])
            horiz = (bb[3] - bb[0]) > (bb[4] - bb[1])
            gmsh.model.mesh.setTransfiniteCurve(c[1], (nx + 1) if horiz else (ny + 1))
        gmsh.model.mesh.setTransfiniteSurface(1)
        gmsh.model.mesh.setRecombine(2, 1)
        occ.synchronize()
        occ.extrude([(2, 1)], 0, 0, TZ, numElements=[1], recombine=True)
        occ.synchronize()
        gmsh.model.mesh.generate(3)
        nt, nc, _ = gmsh.model.mesh.getNodes()
        nodes = {int(t): (nc[3 * i], nc[3 * i + 1], nc[3 * i + 2]) for i, t in enumerate(nt)}
        hexes = []
        for et, _, en in zip(*gmsh.model.mesh.getElements(3)):
            if et == 5:  # 8-node hex
                en = [int(x) for x in en]
                hexes += [en[8 * k:8 * k + 8] for k in range(len(en) // 8)]
        return nodes, hexes
    finally:
        gmsh.finalize()


def _in_notch(cx, cy):
    """True if a hex centroid lies in either symmetric notch bite at mid-length."""
    near_mid = abs(cx - LX / 2.0) < WN
    in_top = cy > HY - DN
    in_bot = cy < DN
    return near_mid and (in_top or in_bot)


def _solve(h, umax=1.5, nsteps=150, regularize=True):
    """Build + displacement-control the DEN bar in tension to elongation umax
    (adaptive multi-level step-cut). `regularize` toggles the crack-band lch handshake
    (the negative control). Returns dict(peak, W, lig_damage, reached, n_elem)."""
    nodes, hexes = _mesh(h)
    keep = [hx for hx in hexes
            if not _in_notch(sum(nodes[n][0] for n in hx) / 8.0,
                             sum(nodes[n][1] for n in hx) / 8.0)]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    used = set(n for hx in keep for n in hx)
    for t in used:
        ops.node(int(t), *nodes[t])
    reg = ["-autoRegularization", LCH_REF] if regularize else []
    ops.nDMaterial("LadrunoJ2", 1, K, G, *ISO, "-damage", "lemaitre", *DMG, *reg)
    lig_eid, eids = None, []
    for i, hx in enumerate(keep, 1):
        ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1, "-formulation", "bbar")
        eids.append(i)
        cx = sum(nodes[n][0] for n in hx) / 8.0
        cy = sum(nodes[n][1] for n in hx) / 8.0
        if abs(cx - LX / 2.0) < h and abs(cy - HY / 2.0) < h:
            lig_eid = i
    xmin, xmax = min(nodes[t][0] for t in used), max(nodes[t][0] for t in used)
    # supports: x=0 face fixed in x; z=0 slice fixed in z; one node pins y (rigid body)
    left = sorted(int(t) for t in used if abs(nodes[t][0] - xmin) < 1e-6)
    for t in used:
        if abs(nodes[t][2]) < 1e-6:
            ops.fix(int(t), 0, 0, 1)
    for t in left:
        ops.fix(int(t), 1, 1 if t == left[0] else 0, 0)   # x on all; pin y at one node
    pull = sorted(int(t) for t in used if abs(nodes[t][0] - xmax) < 1e-6)
    master = pull[0]
    for n in pull:                                  # right face moves together in x
        if n != master:
            ops.equalDOF(master, n, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(master, 1.0, 0.0, 0.0)                # reference tension (DisplacementControl)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 60, 0)
    ops.analysis("Static")

    base = umax / nsteps
    W = uprev = Pprev = Pmax = u = 0.0
    calls = 0
    while u < umax - 1e-12 and calls < 8 * nsteps:
        ok = -1
        for cut, alg in ((1.0, "KrylovNewton"), (0.5, "KrylovNewton"),
                         (0.25, "NewtonLineSearch"), (0.1, "NewtonLineSearch"),
                         (0.03, "NewtonLineSearch")):
            ops.integrator("DisplacementControl", master, 1, base * cut)
            ops.algorithm(alg)
            ok = ops.analyze(1); calls += 1
            if ok == 0:
                break
        if ok != 0:
            break                                  # graceful stop (deep softening)
        ops.reactions()
        P = -sum(ops.nodeReaction(n)[0] for n in left)   # tensile reaction at the fixed end
        u = ops.nodeDisp(master)[0]
        W += 0.5 * (P + Pprev) * (u - uprev)
        Pprev, uprev, Pmax = P, u, max(Pmax, P)
    ligd = max(ops.eleResponse(lig_eid, "material", 1, "damage")) if lig_eid else 0.0
    return {"peak": Pmax, "W": W, "lig_damage": ligd, "reached": u, "n_elem": len(keep)}


def _rel(a, b):
    return abs(a - b) / max(abs(a), 1e-9)


def test_notched_bar_mesh_objectivity():
    """DEN flat bar pulled into the ductile-damage regime at two mesh densities (2x
    refine). With auto-regularization the structural response is mesh-objective: the
    peak load and the dissipated energy agree across meshes. The damage localizes into a
    band ~1-3 element columns wide across the ligament (NOT exactly one element — see the
    geometry figure); lch regularization makes the *dissipated energy* mesh-objective, not
    the band width.

    SCOPE / honesty (per adversarial review): the bar stays in the MILD-damage regime
    (D ~ 0.3-0.4, no full rupture), so the dissipated energy is partly bulk plastic work
    (mesh-objective regardless) plus the localized-damage part. This test gates the
    structural response objectivity; the dedicated negative control
    (test_regularization_discriminates) is what proves the lch handshake itself matters."""
    coarse = _solve(4.0)
    fine = _solve(2.0)

    for tag, r in (("coarse", coarse), ("fine", fine)):
        assert r["reached"] > 1.4, f"{tag}: only reached u={r['reached']:.3f} (<1.4)"
        assert r["lig_damage"] > 0.2, \
            f"{tag}: ligament damage {r['lig_damage']:.2f} — band did not localize"
    assert fine["n_elem"] > 3 * coarse["n_elem"], "fine mesh not actually refined"

    rel_P = _rel(coarse["peak"], fine["peak"])
    assert rel_P < 0.10, (
        f"peak load not mesh-objective: coarse {coarse['peak']:.0f} N vs "
        f"fine {fine['peak']:.0f} N ({rel_P*100:.1f}%)")

    rel_W = _rel(coarse["W"], fine["W"])
    assert rel_W < 0.20, (
        f"dissipated energy not mesh-objective: coarse {coarse['W']:.0f} vs "
        f"fine {fine['W']:.0f} N.mm ({rel_W*100:.1f}%)")


def test_regularization_discriminates():
    """NEGATIVE CONTROL (added after adversarial review): prove the lch handshake actually
    does something — the same coarse/fine pair must be MORE mesh-dependent WITHOUT
    -autoRegularization than WITH it. A test that only ever runs the regularized case
    cannot show the feature works; the earlier '>60% divergence' claim was never asserted
    and is retracted (the real local spread is modest here because damage stays mild and
    bulk plastic work dominates — see test_notched_bar_mesh_objectivity scope note)."""
    reg_c = _solve(4.0, regularize=True)
    reg_f = _solve(2.0, regularize=True)
    loc_c = _solve(4.0, regularize=False)
    loc_f = _solve(2.0, regularize=False)

    spread_reg = _rel(reg_c["W"], reg_f["W"])
    spread_loc = _rel(loc_c["W"], loc_f["W"])
    # regularization must REDUCE the coarse->fine dissipated-energy spread
    assert spread_loc > spread_reg + 0.03, (
        f"regularization did not reduce mesh-dependence: local spread {spread_loc:.1%} "
        f"!> regularized spread {spread_reg:.1%} (+3%) — the lch handshake is not discriminating")
