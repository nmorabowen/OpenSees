"""LadrunoBrick × ASDConcrete3D — softening (crack-band) validation.

Two things are checked, pairing the element with Petracca's plastic-damage
ASDConcrete3D (`SRC/material/nD/ASDConcrete3DMaterial.{h,cpp}`):

1. **Characteristic-length handshake / mesh-objectivity** (the lch handshake of
   Ladruno_implementation/11_brick_asdconcrete_integration.md §1). With
   ``-autoRegularization`` the fracture energy is regularized by the element's
   characteristic length (LadrunoBrick = min edge), so the dissipated energy is
   ``G_f/lch · V = G_f · A`` — i.e. it scales with crack AREA, not VOLUME. A
   single cube pulled to failure at sizes L and 2L must dissipate energy in the
   ratio ``(2L)²/L² = 4`` (area), NOT ``8`` (volume). At L = lch_ref the curve is
   un-rescaled, so the absolute dissipation equals the backbone's specific area ×
   the volume.

2. **Tier-A damage-scaled hourglass stabilization** (item (a) / §3). The single-
   point STIFFNESS-stabilized formulations (ssp, uri+stiffness) degrade their
   constant elastic ``Kstab`` with the material damage:
   ``Kstab ← max(floor, 1 - max(d_t, d_c)) · Kstab_elastic`` (floor = 1%). The
   "hourglassEnergy" report carries that scale, so for an IDENTICAL prescribed
   deformation the stored stabilization energy ratio between an ASDConcrete3D
   element and an elastic element (same E, ν → same elastic Kstab) equals exactly
   ``max(floor, 1 - max(d_t, d_c))`` — degrading toward the floor as the element
   cracks, never to zero (residual hourglass control in the crack band).

Both stay Zone-A (OpenSeesPy only, no apeGmsh, numpy-free). The notched 3-point
bend with mesh refinement (§5) needs meshing and lives in Zone-B.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# ---- a pure-damage uniaxial backbone (E·x rise to f_t, then softening) -------
E, NU, FT = 30000.0, 0.2, 3.0
ET0 = FT / E                                   # peak (elastic-limit) strain
TE = [ET0, 1.0e-3, 8.0e-3]                     # tension strain points
TS = [FT, 1.0, 0.05]                           # tension stress points
TD = [0.0,                                     # damage points: d = 1 - y/(E·x)
      1.0 - 1.0 / (E * 1.0e-3),
      1.0 - 0.05 / (E * 8.0e-3)]
CE, CS, CD = [-ET0, -1.0e-3], [-FT, -1.0], [0.0, 0.0]   # mild compression (unused)
LCH_REF = 1.0
HG_FLOOR = 0.01                                # must match HG_DAMAGE_FLOOR in the C++

_CONN = [1, 2, 3, 4, 5, 6, 7, 8]


def _g_specific():
    """Specific fracture energy (area under the uniaxial backbone, incl. the
    prepended (0,0) point) = energy per unit volume = G_f / lch_ref."""
    xs = [0.0] + TE
    ys = [0.0] + TS
    return sum(0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1))


def _cube(L):
    return {1: (0, 0, 0), 2: (L, 0, 0), 3: (L, L, 0), 4: (0, L, 0),
            5: (0, 0, L), 6: (L, 0, L), 7: (L, L, L), 8: (0, L, L)}


def _concrete(tag):
    ops.nDMaterial("ASDConcrete3D", tag, E, NU, "-rho", 0.0,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-autoRegularization", LCH_REF)


# ===========================================================================
# 1. lch handshake — dissipated energy scales with crack AREA (mesh-objective)
# ===========================================================================
def _dissipation_uniaxial(L, form_args, nsteps=400):
    """Single cube of side L pulled in uniaxial STRESS to failure (symmetry
    rollers, lateral free). Returns the reaction work = total dissipated energy."""
    C = _cube(L)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in C.items():
        ops.node(t, *c)
    _concrete(1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    for n in (1, 4, 5, 8):      # x = 0 face: u_x = 0
        ops.fix(n, 1, 0, 0)
    for n in (1, 2, 5, 6):      # y = 0 face: u_y = 0
        ops.fix(n, 0, 1, 0)
    for n in (1, 2, 3, 4):      # z = 0 face: u_z = 0
        ops.fix(n, 0, 0, 1)
    ops.equalDOF(2, 3, 1); ops.equalDOF(2, 6, 1); ops.equalDOF(2, 7, 1)  # uniform x face
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-8, 50, 0); ops.algorithm("Newton")
    umax = 8.0e-3 * L
    ops.integrator("DisplacementControl", 2, 1, umax / nsteps)
    ops.analysis("Static")
    W, u_prev, F_prev, sig_peak = 0.0, 0.0, 0.0, 0.0
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"{form_args}: pull diverged (L={L})"
        ops.reactions()
        F = -sum(ops.nodeReaction(n)[0] for n in (1, 4, 5, 8))
        u = ops.nodeDisp(2)[0]
        W += 0.5 * (F + F_prev) * (u - u_prev)
        u_prev, F_prev = u, F
        sig = ops.eleResponse(1, "material", 1, "stress")[0]
        sig_peak = max(sig_peak, sig)
    return W, sig_peak


@pytest.mark.parametrize("form_args", [
    ["-formulation", "std"],
    ["-formulation", "bbar"],
    ["-formulation", "ssp"],
])
def test_lch_handshake_mesh_objectivity(form_args):
    """Dissipated energy regularizes by lch: it scales with crack AREA (∝ L²),
    not volume (∝ L³). Pull a cube at L and 2L; the energy ratio must be ≈ 4,
    decisively not 8 — proof the element feeds its characteristic length to
    ASDConcrete3D's auto-regularization. At L = lch_ref the absolute dissipation
    equals the backbone specific energy × the volume."""
    g = _g_specific()
    W1, peak1 = _dissipation_uniaxial(1.0, form_args)
    W2, _ = _dissipation_uniaxial(2.0, form_args)

    # softening actually occurred (peak ≈ f_t, then the curve came down)
    assert peak1 == pytest.approx(FT, rel=0.05), f"{form_args}: peak {peak1:.3f} != f_t"

    # at L = lch_ref the curve is un-rescaled -> W = g_specific · V (V = 1)
    assert W1 == pytest.approx(g, rel=0.05), (
        f"{form_args}: W(L=1) {W1:.4e} != specific fracture energy {g:.4e}")

    # AREA scaling (objective) gives ratio ≈ 4; VOLUME scaling (no handshake) = 8
    ratio = W2 / W1
    assert 3.5 < ratio < 5.0, (
        f"{form_args}: dissipation ratio {ratio:.3f} not ≈4 (area-objective); "
        f"≈8 would mean lch is NOT being regularized (volume scaling)")
    # dissipation per crack area is mesh-objective (≈ constant = G_f)
    Gf1, Gf2 = W1 / 1.0, W2 / 4.0          # A = L² = 1 and 4
    assert Gf2 == pytest.approx(Gf1, rel=0.12), (
        f"{form_args}: G_f not mesh-objective: {Gf1:.4e} vs {Gf2:.4e}")


# ===========================================================================
# 2. Tier-A damage-scaled Kstab — stabilization energy degrades with damage,
#    floored so a fully-cracked element keeps residual hourglass control
# ===========================================================================
# raw hourglass base vectors (LadrunoBrick.cpp hg0[0], hg0[1]): nodal patterns
# orthogonal to the linear field -> they drive hourglass energy but add NO
# centroid strain, so the centroid material (slot 0) damages from the uniaxial
# part alone.
_HG1 = [1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0]
_HG2 = [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0]


def _prescribed_run(form_args, material, nsteps=120, eps_max=5.0e-3, beta=1.0e-3):
    """Single unit cube with ALL 24 dofs prescribed (Penalty): a uniform uniaxial
    stretch (ramps centroid strain -> damage) plus two pure hourglass modes
    (ramp -> stored stabilization energy). Returns per-step (Ehg, max damage)."""
    C = _cube(1.0)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in C.items():
        ops.node(t, *c)
    if material == "concrete":
        _concrete(1)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for j, (nid, (X, Y, Z)) in enumerate(C.items()):
        ops.sp(nid, 1, eps_max * X)        # uniaxial x  -> centroid eps_xx
        ops.sp(nid, 2, beta * _HG1[j])     # hourglass in y
        ops.sp(nid, 3, beta * _HG2[j])     # hourglass in z
    ops.constraints("Penalty", 1e14, 1e14); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0 / nsteps); ops.algorithm("Linear")
    ops.analysis("Static")
    out = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"{form_args}/{material}: prescribed step diverged"
        ops.eleResponse(1, "forces")
        Ehg = ops.eleResponse(1, "hourglassEnergy")[0]
        if material == "concrete":
            dmg = max(ops.eleResponse(1, "material", 1, "damage"))
        else:
            dmg = 0.0
        out.append((Ehg, dmg))
    return out


@pytest.mark.parametrize("form_args", [
    ["-formulation", "ssp"],
    ["-formulation", "uri", "-hourglass", "stiffness"],
])
def test_tier_a_damage_scaled_kstab(form_args):
    """For an IDENTICAL prescribed deformation, the stored stabilization energy of
    an ASDConcrete3D element relative to an elastic element (same E, ν → identical
    elastic Kstab) must equal exactly max(floor, 1 - max(d_t, d_c)) at every step:
    full elastic value when intact, degrading with damage, never below the floor."""
    elastic = _prescribed_run(form_args, "elastic")
    concrete = _prescribed_run(form_args, "concrete")

    saw_intact = saw_floor = False
    for (Ehg_e, _), (Ehg_c, dmax) in zip(elastic, concrete):
        assert Ehg_e > 1e-12, f"{form_args}: elastic stabilization energy ~0 (bad probe)"
        s_meas = Ehg_c / Ehg_e
        s_exp = max(HG_FLOOR, 1.0 - dmax)
        assert s_meas == pytest.approx(s_exp, abs=2e-3), (
            f"{form_args}: damage scale {s_meas:.4f} != max(floor,1-d)={s_exp:.4f} "
            f"(dmax={dmax:.4f})")
        assert Ehg_c > 0.0, f"{form_args}: cracked element fully unstabilized (no floor)"
        if dmax < 1e-3:
            saw_intact = True
            assert s_meas == pytest.approx(1.0, abs=2e-3)   # intact -> elastic Kstab
        if 1.0 - dmax < HG_FLOOR:
            saw_floor = True
            assert s_meas == pytest.approx(HG_FLOOR, abs=2e-3)  # floored, not zero

    assert saw_intact, f"{form_args}: never sampled the intact (elastic) regime"
    assert saw_floor, f"{form_args}: never reached the damage floor"
