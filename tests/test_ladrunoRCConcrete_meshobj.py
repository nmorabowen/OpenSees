"""LadrunoRCConcrete Phase-3b crack-band regularization — STRUCTURAL mesh-objectivity gate
(`-autoRegularization`, classTag 33015) — Zone-A, the material AS INTEGRATED IN ASDShellQ4.

The Phase-3b unit battery (`test_ladrunoRCConcrete_reg.py`) proves the regularize() *mechanism*
at a single material point: the post-peak dissipated energy DENSITY scales as g_reg = G_f0 *
(lch_ref/lch), so g_density * lch is constant. This file closes the loop the 3b review staged: the
*structural* proof that an actual softening FE response — a row of `ASDShellQ4` driven to full
tensile softening — dissipates a MESH-OBJECTIVE total energy, because each element latches its own
`getCharacteristicLength()` and the kernel rescales the softening accordingly.

The classical Bazant tension-bar set-up:

  * A short uniaxial bar built from square `ASDShellQ4` membrane elements, every transverse and
    out-of-plane DOF locked (a clean 1-D X-chain), the right edge pulled in displacement control.
  * The LEFTMOST element column is given a 10%-thinner section (`weak` band) so the crack
    LOCALIZES there deterministically — the rest unloads elastically. (Thickness, not the material:
    the backbone + regularization are identical everywhere, only the cross-section triggers the
    band, so the band-dissipation factor cancels cleanly between meshes.)
  * The SAME 100x50 specimen is meshed three ways — 2x1 (s=50), 4x2 (s=25), 8x4 (s=12.5) — so the
    element `lch` halves twice while the geometry is unchanged.

Without regularization the localized band is one element wide, so as the mesh refines the
dissipated energy spuriously SHRINKS with the element size (the well-known pathological mesh
dependence — here fine/coarse ~ 0.5). WITH `-autoRegularization` the softening density is scaled
up by lch_ref/lch, so band_density * band_width stays constant and the total dissipated energy is
mesh-OBJECTIVE. Both halves are asserted: the regularized sequence collapses to one value
(objective), and the un-regularized negative control demonstrably does NOT (the gate has teeth).

Scope: an axis-aligned (crack normal to the element edges) bar — the cleanest first cut and
gmsh-free, hence Zone-A. The inclined-crack / notched panel (which additionally exercises the
documented single-scalar-`lch` sqrt(2) mis-regularization of ~45-deg struts, an Option-A residual,
and is snap-back-prone) is the harder Zone-B follow-up.

Plan: Ladruno_implementation/19_ladruno_rc_shell_adr.md (Phase 3b) and the forward handout
rc_shell_phase3plus_handout.md (3b-struct).
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# same softening backbone family as the reg unit battery; tension branch decays to zero stress
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010, 0.0040]
TS = [0.0, 3.0,    0.5,    0.0]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0, 0.999]
FT = max(TS)                       # tensile peak
LCH_REF = 50.0                     # author the backbone at the coarse element size (coarse = no-op)
T = 100.0                          # nominal section thickness
NLAY = 4
WEAK = 0.90                        # left band 10% thinner -> deterministic localization

# the SAME 100 x 50 specimen, three meshes (square elements; lch halves twice)
MESHES = [(2, 1, 50.0), (4, 2, 25.0), (8, 4, 12.5)]

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # np.trapz removed in numpy 2.0


def _mat(tag, auto_reg):
    a = ["LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
         "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if auto_reg:
        a += ["-autoRegularization", LCH_REF]      # NO -lch: each element latches its own size
    ops.nDMaterial(*a)


def _sec(tag, thick):
    layers = []
    for _ in range(NLAY):
        layers += [1, thick / NLAY]
    ops.section("LayeredShell", tag, NLAY, *layers)


def _run_bar(ncols, nrows, s, auto_reg, dmax=0.6, nsteps=400):
    """Build the SAME 100x50 RC membrane bar at the given mesh, pull the right edge to dmax in
    displacement control, return (peak_force, dissipated_energy, n_steps_completed)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nid, k = {}, 1
    for i in range(ncols + 1):
        for j in range(nrows + 1):
            nid[(i, j)] = k
            ops.node(k, i * s, j * s, 0.0)
            k += 1
    _mat(1, auto_reg)
    _sec(1, T)            # normal section
    _sec(2, T * WEAK)     # weak (thinner) left band
    eid = 1
    for i in range(ncols):
        for j in range(nrows):
            sec = 2 if i == 0 else 1
            ops.element("ASDShellQ4", eid, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], sec)
            eid += 1
    probe_elem = eid - 1                           # rightmost column, normal section
    # pure 1-D X-chain: lock Y,Z + all rotations everywhere; left edge fixed in X
    for n in nid.values():
        ops.fix(n, 0, 1, 1, 1, 1, 1)
    for j in range(nrows + 1):
        ops.fix(nid[(0, j)], 1, 0, 0, 0, 0, 0)     # dof2..6 already fixed above
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(nrows + 1):                     # prescribe the right edge (scaled by LoadControl)
        ops.sp(nid[(ncols, j)], 1, dmax)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    Wy = nrows * s
    right_node = nid[(ncols, 0)]
    F, d = [], []
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break                                  # un-regularized fine meshes snap-back (by design)
        ops.eleResponse(probe_elem, "forces")      # prime lazy state
        nxx = list(ops.eleResponse(probe_elem, "stresses"))[0]
        F.append(nxx * Wy)                          # series bar: axial force uniform => F = Nxx * Wy
        d.append(ops.nodeDisp(right_node, 1))
    F, d = np.array(F), np.array(d)
    peak = float(F.max()) if len(F) else 0.0
    energy = float(_trapz(F, d)) if len(F) > 1 else 0.0   # external work ~ dissipation once F -> 0
    return peak, energy, len(F)


@pytest.fixture(scope="module")
def bars():
    """Run every (regularization, mesh) once. The un-regularized finest mesh is brittle and
    snap-backs (it is the disease, not asserted on) so it is omitted."""
    out = {}
    for reg in (False, True):
        for nc, nr, s in MESHES:
            if (not reg) and s < 20.0:
                continue                            # reg-OFF 8x4 diverges by design
            out[(reg, s)] = _run_bar(nc, nr, s, reg)
    return out


# --------------------------------------------------------------------------
def test_meshobj_regon_energy_is_mesh_objective(bars):
    """The headline: WITH -autoRegularization the total dissipated energy is mesh-objective —
    the 2x1 / 4x2 / 8x4 sequence (lch 50/25/12.5, a 4x refinement) collapses to one value."""
    es = [bars[(True, s)][1] for _, _, s in MESHES]
    ns = [bars[(True, s)][2] for _, _, s in MESHES]
    assert all(n == 400 for n in ns), f"a regularized run did not complete: steps {ns}"
    ref = es[0]
    for s, e in zip((50.0, 25.0, 12.5), es):
        assert e == pytest.approx(ref, rel=0.05), f"lch={s}: energy {e} not objective vs {ref} ({es})"


def test_meshobj_regon_peak_is_mesh_objective(bars):
    """Peak capacity (governed by ft x the weak cross-section) is mesh-objective with reg on."""
    ps = [bars[(True, s)][0] for _, _, s in MESHES]
    # ft * (weak thickness) * width
    assert ps[0] == pytest.approx(FT * (T * WEAK) * 50.0, rel=0.10), ps
    for p in ps:
        assert p == pytest.approx(ps[0], rel=0.03), f"peak not objective: {ps}"


def test_meshobj_unregularized_is_mesh_dependent(bars):
    """Negative control — the gate has teeth: WITHOUT regularization the dissipated energy
    spuriously SHRINKS with the element size. The band halves (s: 50 -> 25), so the dissipated
    energy roughly halves — a gross mesh dependence the regularized run does not have."""
    e_coarse = bars[(False, 50.0)][1]
    e_fine = bars[(False, 25.0)][1]
    assert e_fine < 0.65 * e_coarse, (
        f"un-regularized energy not mesh-dependent (coarse {e_coarse}, fine {e_fine}, "
        f"ratio {e_fine / e_coarse:.3f}) — the objectivity gate would be vacuous")


def test_meshobj_regularization_restores_ductility(bars):
    """At the fine mesh the regularization must ADD dissipation (the rescaled longer tail): the
    regularized 4x2 bar dissipates substantially more than the un-regularized one — the same
    energy the coarse mesh already had, i.e. regularization is actively compensating."""
    e_reg = bars[(True, 25.0)][1]
    e_off = bars[(False, 25.0)][1]
    assert e_reg > 1.5 * e_off, f"regularization did not restore fine-mesh dissipation: reg {e_reg}, off {e_off}"
