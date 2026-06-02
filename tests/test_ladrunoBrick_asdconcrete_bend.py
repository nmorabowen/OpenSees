"""LadrunoBrick × ASDConcrete3D — meshed crack-band validation (Zone-B, gmsh).

The single-element Zone-A suite (test_ladrunoBrick_asdconcrete.py) proves the lch
handshake by element-size scaling and the exact Tier-A damage-scale formula. These
two cases confirm both properties in a *full structural crack* on a gmsh-meshed
single-edge-notched 3-point-bend (SENB) beam — the §5 cases of
Ladruno_implementation/11_brick_asdconcrete_integration.md that need real meshing:

1. **Fracture-energy / mesh objectivity** — with `-autoRegularization` the dissipated
   energy regularizes by the element characteristic length, so the structural
   load–deflection response and the dissipated energy are MESH-INDEPENDENT. Two
   meshes (coarse / fine, 2× refinement) must give the same load at a common
   deflection and the same dissipated energy. A broken handshake (energy ∝ volume)
   would make the fine mesh markedly more brittle and the curves diverge.

2. **Tier-A hourglass-energy fraction in the cracked band** — running the same beam
   with the single-point `eas` formulation, the stored hourglass/stabilization
   energy must stay a small fraction (< ~10%) of the internal energy even as the
   crack fully forms. This is the §5 acceptance dial for the damage-scaled `Kstab`:
   it confirms the degraded stabilization neither props the cracked load up nor
   lets the spurious modes dominate.

Zone-B: needs gmsh for meshing (raw gmsh API, not apeGmsh); runs nightly.
"""
import pytest

from _testbed import ops

gmsh = pytest.importorskip("gmsh")          # belt-and-suspenders (conftest also gates zone_b on gmsh)
pytestmark = [pytest.mark.zone_b]

# --- SENB geometry (N, mm) ---
L, S, D, T_, A0, WN = 220.0, 200.0, 50.0, 25.0, 10.0, 6.0
# --- concrete: gentle (ductile) regularized tension backbone -> shallow softening,
#     so load-point displacement control + step-cutting traverses the full crack ---
E, NU, FT = 30000.0, 0.2, 3.0
ET0 = FT / E
TE = [ET0, 2.0e-3, 1.5e-2]
TS = [FT, 1.8, 0.15]
TD = [0.0, 1.0 - 1.8 / (E * 2.0e-3), 1.0 - 0.15 / (E * 1.5e-2)]
CE, CS, CD = [-ET0, -2.0e-3], [-FT, -25.0], [0.0, 0.0]
LCH_REF = 10.0


def _mesh(h):
    """Structured quad mesh of the L×D face (size ~h), recombined + extruded one
    layer in z -> hexes (gmsh hex node order == OpenSees Brick order). Returns
    (nodes{tag:(x,y,z)}, hexes[[8 tags]])."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("senb")
        occ = gmsh.model.occ
        occ.addRectangle(0, 0, 0, L, D)
        occ.synchronize()
        nx, ny = max(2, round(L / h)), max(2, round(D / h))
        for c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c[1])
            horiz = (bb[3] - bb[0]) > (bb[4] - bb[1])
            gmsh.model.mesh.setTransfiniteCurve(c[1], (nx + 1) if horiz else (ny + 1))
        gmsh.model.mesh.setTransfiniteSurface(1)
        gmsh.model.mesh.setRecombine(2, 1)
        occ.synchronize()
        occ.extrude([(2, 1)], 0, 0, T_, numElements=[1], recombine=True)
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


def _solve(h, form, dmax=0.06, nsteps=80, monitor_hg=False):
    """Build + displacement-control the SENB beam to deflection dmax (adaptive
    step-cut through softening). Returns dict(load_at_dmax, dissipated, tip_damage,
    reached, hg_frac_max)."""
    nodes, hexes = _mesh(h)
    keep = [hx for hx in hexes
            if not (sum(nodes[n][1] for n in hx) / 8.0 < A0
                    and abs(sum(nodes[n][0] for n in hx) / 8.0 - L / 2.0) < WN)]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    used = set(n for hx in keep for n in hx)
    for t in used:
        ops.node(int(t), *nodes[t])
    ops.nDMaterial("ASDConcrete3D", 1, E, NU, "-rho", 0.0,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-autoRegularization", LCH_REF, "-implex", "-implexControl", 0.05, 0.01)
    tip_eid, eids = None, []
    for i, hx in enumerate(keep, 1):
        ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1, "-formulation", form)
        eids.append(i)
        cx = sum(nodes[n][0] for n in hx) / 8.0
        cy = sum(nodes[n][1] for n in hx) / 8.0
        if abs(cx - L / 2.0) < WN and A0 <= cy < A0 + 1.5 * h:
            tip_eid = i
    ov = (L - S) / 2.0
    botY, topY = min(nodes[t][1] for t in used), max(nodes[t][1] for t in used)
    for t in used:                                   # dof3 on z=0 face (1-layer slice)
        if abs(nodes[t][2]) < 1e-6:
            ops.fix(int(t), 0, 0, 1)
    sup = []
    for t in used:                                   # supports touch only dof1/dof2
        x, y, _ = nodes[t]
        if abs(y - botY) < 1e-6 and abs(x - ov) < h * 0.6:
            ops.fix(int(t), 1, 1, 0); sup.append(int(t))
        elif abs(y - botY) < 1e-6 and abs(x - (L - ov)) < h * 0.6:
            ops.fix(int(t), 0, 1, 0); sup.append(int(t))
    ld = sorted(int(t) for t in used
                if abs(nodes[t][1] - topY) < 1e-6 and abs(nodes[t][0] - L / 2.0) < h * 0.6)
    master = ld[0]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in ld:
        ops.load(n, 0.0, -1.0 / len(ld), 0.0)        # downward reference load (DisplacementControl)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 80, 0)
    ops.analysis("Static")

    base = dmax / nsteps
    W = dprev = Pprev = Pmax = d = hg_frac_max = 0.0
    calls = 0
    while d < dmax - 1e-12 and calls < 4000:
        ops.integrator("DisplacementControl", master, 2, -base)
        ops.algorithm("KrylovNewton")
        ok = ops.analyze(1); calls += 1
        if ok != 0:                                  # adaptive step-cut through softening
            for cut in (0.5, 0.25, 0.1):
                ops.integrator("DisplacementControl", master, 2, -base * cut)
                ops.algorithm("NewtonLineSearch")
                ok = ops.analyze(1); calls += 1
                if ok == 0:
                    break
            if ok != 0:
                break                                # graceful stop (full fracture / limit)
        ops.reactions()
        R = sum(ops.nodeReaction(n)[1] for n in sup)
        d = -ops.nodeDisp(master)[1]
        W += 0.5 * (R + Pprev) * (d - dprev)
        Pprev, dprev, Pmax = R, d, max(Pmax, R)
        if monitor_hg and W > 1e-6 and d > 0.01:
            ehg = sum((ops.eleResponse(e, "hourglassEnergy") or [0.0])[0] for e in eids)
            hg_frac_max = max(hg_frac_max, ehg / W)
    tipd = max(ops.eleResponse(tip_eid, "material", 1, "damage")) if tip_eid else 0.0
    return {"load": Pprev, "W": W, "tip_damage": tipd, "reached": d,
            "hg_frac_max": hg_frac_max, "n_elem": len(keep)}


# ===========================================================================
# 1. fracture-energy / mesh objectivity (load-deflection + dissipation)
# ===========================================================================
def test_notched_bend_mesh_objectivity():
    """SENB beam cracked under 3-point bending at two mesh densities (2× refine).
    With auto-regularization the structural response is mesh-objective: the load at
    a common deflection and the dissipated energy must agree across meshes (a broken
    lch handshake would make the fine mesh brittle and the curves diverge)."""
    coarse = _solve(5.0, "bbar")
    fine = _solve(2.5, "bbar")

    # both must crack (localized band at the notch tip) and reach the target deflection
    for tag, r in (("coarse", coarse), ("fine", fine)):
        assert r["reached"] > 0.05, f"{tag}: only reached d={r['reached']:.3f} (<0.05)"
        assert r["tip_damage"] > 0.8, f"{tag}: notch-tip damage {r['tip_damage']:.2f} — no crack"
    assert fine["n_elem"] > 3 * coarse["n_elem"]      # genuine refinement

    # load–deflection objectivity at the common deflection
    rel_P = abs(coarse["load"] - fine["load"]) / max(abs(coarse["load"]), 1.0)
    assert rel_P < 0.10, (
        f"load not mesh-objective: coarse {coarse['load']:.1f} N vs fine {fine['load']:.1f} N "
        f"({rel_P*100:.1f}%)")
    # dissipated-energy (fracture-energy) objectivity
    rel_W = abs(coarse["W"] - fine["W"]) / max(abs(coarse["W"]), 1e-9)
    assert rel_W < 0.20, (
        f"dissipated energy not mesh-objective: coarse {coarse['W']:.1f} vs fine {fine['W']:.1f} "
        f"N·mm ({rel_W*100:.1f}%)")


# ===========================================================================
# 2. Tier-A: hourglass energy stays a small fraction of internal energy in the
#    cracked band (the §5 acceptance dial for the damage-scaled Kstab)
# ===========================================================================
def test_eas_hourglass_energy_fraction_in_cracked_band():
    """Same beam with the single-point `eas` formulation. As the crack forms, the
    damage-scaled stabilization energy must stay a small fraction of the internal
    energy — confirming the degraded Kstab neither props the cracked load up nor
    lets spurious hourglass modes dominate."""
    r = _solve(5.0, "eas", monitor_hg=True)
    assert r["reached"] > 0.05, f"eas run only reached d={r['reached']:.3f}"
    assert r["tip_damage"] > 0.5, (
        f"notch-tip damage {r['tip_damage']:.2f} — not actually in the cracked regime")
    assert r["hg_frac_max"] > 0.0, "hourglass energy never sampled (bad probe)"
    assert r["hg_frac_max"] < 0.10, (
        f"hourglass energy reached {r['hg_frac_max']*100:.1f}% of internal energy in the "
        f"cracked band (> 10% ceiling): stabilization out of control")
