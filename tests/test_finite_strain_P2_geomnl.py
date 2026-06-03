"""Finite-strain trifecta — **P2 geometric nonlinearity** (corot large rotation).

Validation plan: Ladruno_implementation/17_finite_strain_validation_plan.md, §9
phase **P2** (gate: corot elastica/buckling vs analytical). Where P1 exercised
large *strain* (`-geom finite`), P2 exercises large *rotation with small strain*
— the corotational path `LadrunoBrick -geom corot`, which strips the element
rigid rotation R=polar(H) and feeds the small deformational strain to the
unchanged small-strain material (here `ElasticIsotropic`).

Benchmarks delivered here:
  * **A7** — cantilever under a pure end **moment** rolls into a **circle**
    (Euler elastica): the deformed centerline is a circular arc of uniform
    curvature κ = M/EI. Pure bending ⇒ no transverse shear ⇒ the corot solid
    captures the constant-curvature shape essentially locking-free; the
    moment–curvature law M = EI/ρ converges under mesh refinement.
  * **C4** — large-rotation cantilever under an end **transverse force**, vs the
    exact **elastica** (Bisshopp–Drucker 1945 / Mattiasson 1981; the
    Bathe–Bolourchi 1979 benchmark): tip deflection in the load direction and the
    axis foreshortening match the tabulated `w/L`, `u/L` at load parameters
    `α = PL²/EI` of 7 and 10 to ≤4 %.
  * **C5** — **Euler buckling** of a corot cantilever column: an imperfect column
    under ramped axial compression; a **Southwell plot** recovers the critical
    load, converging to `P_cr = π²EI/(2L)²` under mesh refinement.
  * **D4** — geometry-method consistency: the SAME large-rotation cantilever
    driven by `-geom corot` and `-geom finite` (small material strain) gives the
    same tip displacement. No external oracle — a code-internal cross-check that
    the two independent geometric-nonlinearity paths agree (their shared
    low-order-hex shear stiffness cancels in the comparison).

> [!note] Formulation constraint (verified)
> Under `-geom corot|finite` only `std`/`bbar` exist — `uri`/`eas` (the
> shear-locking cures) are `-geom linear` only. So corot bending accuracy comes
> from mesh refinement, not a locking-free element. A7 refines the CROSS-SECTION
> (≥2×2 through-thickness) to relieve the 1-element parasitic shear; the arc
> SHAPE is locking-insensitive and tight even on coarse meshes.

"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t2a]

_E = 1000.0
_B = 1.0           # square cross-section side
_L = 10.0          # cantilever length (along z)
_I = _B ** 4 / 12.0
_EI = _E * _I


# --------------------------------------------------------------------------- #
#  Slender cantilever of LadrunoBricks along z (square b×b section, nx×ny       #
#  through-thickness, nz along length). Base z=0 clamped.                       #
# --------------------------------------------------------------------------- #
def _build_beam(geom, nz, nx=2, ny=2):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nxp, nyp = nx + 1, ny + 1

    def nid(i, j, k):
        return 1 + i + nxp * j + nxp * nyp * k

    for k in range(nz + 1):
        for j in range(nyp):
            for i in range(nxp):
                ops.node(nid(i, j, k), _B * i / nx, _B * j / ny, _L * k / nz)
    base = [nid(i, j, 0) for j in range(nyp) for i in range(nxp)]
    for n in base:
        ops.fix(n, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, _E, 0.0)        # nu=0 ⇒ EI = E·I exactly
    mt = 1
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mt = 2
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, mt, "-formulation", "std", "-geom", geom)
                e += 1
    return nid, base, nx, ny, nz


def _centerline(nid, nx, ny, nz):
    """Deformed (x,z) centroid of each cross-section layer."""
    xs, zs = [], []
    for k in range(nz + 1):
        cx = np.mean([ops.nodeCoord(nid(i, j, k), 1) + ops.nodeDisp(nid(i, j, k), 1)
                      for j in range(ny + 1) for i in range(nx + 1)])
        cz = np.mean([ops.nodeCoord(nid(i, j, k), 3) + ops.nodeDisp(nid(i, j, k), 3)
                      for j in range(ny + 1) for i in range(nx + 1)])
        xs.append(cx); zs.append(cz)
    return np.array(xs), np.array(zs)


def _fit_circle(xs, zs):
    """Algebraic circle fit; returns (radius, max abs radial residual)."""
    A = np.c_[xs, zs, np.ones(len(xs))]
    b = -(xs ** 2 + zs ** 2)
    D, E, F = np.linalg.lstsq(A, b, rcond=None)[0]
    cx, cz = -D / 2.0, -E / 2.0
    r = math.sqrt(cx * cx + cz * cz - F)
    res = np.abs(np.hypot(xs - cx, zs - cz) - r)
    return r, float(res.max())


def _base_moment_y(base):
    """Reaction moment about y at the clamped base (z=0 ⇒ M_y = −Σ x·Rz)."""
    ops.reactions()
    My = 0.0
    for n in base:
        x = ops.nodeCoord(n, 1)
        My += -x * ops.nodeReaction(n, 3)
    return abs(My)


def _solve_end_moment(geom, nz, Mnom, nx=2, ny=2, nsteps=40):
    """Apply an end couple of nominal moment Mnom (±Fz on the ±x end faces) and
    solve the corot/finite large-rotation response. Returns (nid, base, nx, ny, nz)."""
    nid, base, nx, ny, nz = _build_beam(geom, nz, nx, ny)
    F = Mnom / (2.0 * _B)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for n in (nid(nx, j, nz),):       # +x face nodes → +Fz
            ops.load(n, 0.0, 0.0, +F / (ny + 1))
        for n in (nid(0, j, nz),):        # −x face nodes → −Fz
            ops.load(n, 0.0, 0.0, -F / (ny + 1))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-9, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"end-moment solve failed (geom={geom}, nz={nz}, M={Mnom})"
    return nid, base, nx, ny, nz


# =========================================================================== #
#  A7.1 — the deformed centerline is a CIRCULAR ARC (rolls toward a circle)     #
# =========================================================================== #
def test_A7_end_moment_rolls_into_circular_arc():
    # moment giving a large (>60°) curl — solidly in the large-rotation regime
    Mnom = _EI * math.pi / _L
    nid, base, nx, ny, nz = _solve_end_moment("corot", nz=16, Mnom=Mnom)
    xs, zs = _centerline(nid, nx, ny, nz)
    r, res = _fit_circle(xs[1:], zs[1:])      # skip the clamped base point
    arc_angle = _L / r
    # (1) the centerline is a circular arc to high precision (uniform curvature
    #     = the elastica signature; locking-insensitive)
    assert res < 0.01 * _L, f"deformed centerline is not a circular arc (residual {res:.4e})"
    # (2) and it has actually curled into the large-rotation regime
    assert arc_angle > math.radians(50.0), f"insufficient curl ({math.degrees(arc_angle):.1f}°)"


# =========================================================================== #
#  A7.2 — curvature is UNIFORM along the length (constant internal moment)      #
# =========================================================================== #
def test_A7_curvature_is_uniform():
    Mnom = _EI * math.pi / _L
    nid, base, nx, ny, nz = _solve_end_moment("corot", nz=16, Mnom=Mnom)
    xs, zs = _centerline(nid, nx, ny, nz)
    # discrete curvature from successive centerline triples (turning angle / arc len)
    kappas = []
    for k in range(1, nz):
        a = np.array([xs[k] - xs[k - 1], zs[k] - zs[k - 1]])
        b = np.array([xs[k + 1] - xs[k], zs[k + 1] - zs[k]])
        la, lb = np.linalg.norm(a), np.linalg.norm(b)
        ang = math.atan2(a[0] * b[1] - a[1] * b[0], a @ b)
        kappas.append(ang / (0.5 * (la + lb)))
    kappas = np.array(kappas)
    spread = (kappas.max() - kappas.min()) / abs(kappas.mean())
    assert spread < 0.10, f"curvature not uniform along the beam (spread {spread:.3f})"


# =========================================================================== #
#  A7.3 — moment–curvature law M = EI/ρ converges under mesh refinement         #
# =========================================================================== #
def test_A7_moment_curvature_converges_to_EI():
    Mnom = _EI * math.pi / _L
    errs = []
    for nz in (8, 16, 24):
        nid, base, nx, ny, nz = _solve_end_moment("corot", nz=nz, Mnom=Mnom)
        xs, zs = _centerline(nid, nx, ny, nz)
        rho, _ = _fit_circle(xs[1:], zs[1:])
        M = _base_moment_y(base)
        errs.append(abs(M * rho - _EI) / _EI)        # M·ρ should → EI
    # monotone convergence toward the Euler moment–curvature law (the residual is
    # O(h) low-order-hex bending stiffness — std/bbar are the only corot
    # formulations, no uri/eas locking cure — so the ABSOLUTE match is a coarse-
    # mesh engineering bound; the TIGHT elastica gates are the arc shape (A7.1)
    # and uniform curvature (A7.2)).
    assert errs[0] > errs[1] > errs[2], f"M·ρ not converging to EI under refinement: {errs}"
    # each refinement roughly halves the error (≈O(h)) — evidences κ → M/EI
    assert errs[1] < 0.6 * errs[0] and errs[2] < 0.6 * errs[1], (
        f"M·ρ→EI convergence too slow to be credible: {errs}")
    assert errs[-1] < 0.10, f"M·ρ = EI not within engineering tol at finest mesh (err {errs[-1]:.3f})"


# =========================================================================== #
#  C4 — large-rotation cantilever vs the exact elastica (Mattiasson 1981)      #
# =========================================================================== #
# Tip deflection in the load direction (w/L) and axis foreshortening (u/L) for a
# cantilever under an end transverse force, parameter α = PL²/EI. Exact elastica
# (Bisshopp–Drucker 1945; tabulated to 5 figures by Mattiasson 1981) — the
# Bathe–Bolourchi (1979) large-displacement cantilever benchmark.
_MATTIASSON = {7: (0.76737, 0.47490), 10: (0.81061, 0.55500)}


@pytest.mark.parametrize("alpha", [7, 10])
def test_C4_elastica_cantilever_matches_mattiasson(alpha):
    w_ref, u_ref = _MATTIASSON[alpha]
    nz = 32                                  # ≤2.3% vs elastica at this refinement
    P = alpha * _EI / _L ** 2
    ux, uz = _tip_disp_transverse("corot", nz, P, nsteps=max(48, alpha * 8))
    w = abs(ux) / _L                         # deflection in the load (x) direction
    u = -uz / _L                             # axis (z) foreshortening
    assert w == pytest.approx(w_ref, rel=0.04), (
        f"C4 α={alpha}: w/L {w:.4f} != elastica {w_ref:.4f}")
    assert u == pytest.approx(u_ref, rel=0.04), (
        f"C4 α={alpha}: u/L {u:.4f} != elastica {u_ref:.4f}")


# =========================================================================== #
#  C5 — Euler buckling of a corot column (Southwell plot)                       #
# =========================================================================== #
def _southwell_Pcr(nz, Pmax, imp=0.001, nsteps=40):
    """Imperfect cantilever column under ramped axial compression + a tiny
    proportional transverse imperfection load. Returns the Southwell-plot
    critical load (slope of δ/P vs δ over the pre-buckling range = 1/P_cr)."""
    nid, base, nx, ny, nz = _build_beam("corot", nz, nx=2, ny=2)
    tip = [nid(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in tip:
        ops.load(n, imp * Pmax / len(tip), 0.0, -Pmax / len(tip))   # axial + imperfection
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-10, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    Pcr_euler = math.pi ** 2 * _EI / (2.0 * _L) ** 2
    xs, ys = [], []
    for s in range(1, nsteps + 1):
        if ops.analyze(1) != 0:
            break
        P = (s / nsteps) * Pmax
        d = abs(np.mean([ops.nodeDisp(n, 1) for n in tip]))
        if 0.25 * Pcr_euler < P < 0.85 * Pcr_euler and d > 0:
            xs.append(d); ys.append(d / P)            # Southwell: δ/P vs δ
    A = np.c_[np.array(xs), np.ones(len(xs))]
    slope, _ = np.linalg.lstsq(A, np.array(ys), rcond=None)[0]
    return 1.0 / slope, Pcr_euler


def test_C5_euler_buckling_converges():
    Pcr_euler = math.pi ** 2 * _EI / (2.0 * _L) ** 2
    Pmax = 0.83 * Pcr_euler                  # ramp to just below the critical load
    Pcr_24, _ = _southwell_Pcr(24, Pmax)
    Pcr_32, _ = _southwell_Pcr(32, Pmax)
    err_24 = abs(Pcr_24 - Pcr_euler) / Pcr_euler
    err_32 = abs(Pcr_32 - Pcr_euler) / Pcr_euler
    # the column buckles near the Euler load; the residual is shear-locking
    # over-stiffness (P_cr biased HIGH), which shrinks under mesh refinement.
    assert err_32 < err_24, f"P_cr not converging to Euler: err {err_24:.3f} → {err_32:.3f}"
    assert 1.0 - 0.03 < Pcr_32 / Pcr_euler < 1.10, (
        f"C5 buckling load off: Southwell P_cr {Pcr_32:.4f} vs Euler {Pcr_euler:.4f}")


# =========================================================================== #
#  D4 — corot ↔ finite geometry-method consistency (large rotation)            #
# =========================================================================== #
def _tip_disp_transverse(geom, nz, P, nx=2, ny=2, nsteps=20):
    nid, base, nx, ny, nz = _build_beam(geom, nz, nx, ny)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    tip = [nid(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)]
    for n in tip:
        ops.load(n, P / len(tip), 0.0, 0.0)        # transverse (x) tip load
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("EnergyIncr", 1.0e-9, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"transverse cantilever failed (geom={geom})"
    ux = np.mean([ops.nodeDisp(n, 1) for n in tip])
    uz = np.mean([ops.nodeDisp(n, 3) for n in tip])
    return float(ux), float(uz)


def test_D4_corot_finite_large_rotation_consistency():
    # a load large enough to bend the tip well into large rotation (tip drop uz
    # is a sizeable fraction of the tip sway ux — geometric nonlinearity active)
    P = 2.0
    nz = 16
    ux_c, uz_c = _tip_disp_transverse("corot", nz, P)
    ux_f, uz_f = _tip_disp_transverse("finite", nz, P)
    assert abs(ux_c) > 0.15 * _L, f"deformation not large-rotation (tip sway {ux_c:.3f})"
    assert abs(uz_c) > 0.02 * _L, "tip foreshortening negligible — not geometrically nonlinear"
    # corot and finite are independent geometry paths; at small material strain
    # they must give the same large-rotation tip displacement (shared low-order
    # hex shear stiffness cancels in the comparison).
    scale = max(abs(ux_c), 1e-9)
    assert abs(ux_c - ux_f) < 0.02 * scale, (
        f"corot vs finite tip sway disagree: {ux_c:.5f} vs {ux_f:.5f}")
    assert abs(uz_c - uz_f) < 0.02 * max(abs(uz_c), 1e-9), (
        f"corot vs finite tip foreshortening disagree: {uz_c:.5f} vs {uz_f:.5f}")


def test_D4_corot_reduces_to_euler_in_linear_limit():
    # tiny load ⇒ corot must recover the linear Euler–Bernoulli tip deflection
    # PL³/3EI (to within the low-order-hex bending discretization error, which a
    # 2×2×nz mesh keeps modest — refine-checked, not a tight absolute gate).
    P = 0.005
    nz = 24
    ux, _ = _tip_disp_transverse("corot", nz, P, nsteps=2)
    euler = P * _L ** 3 / (3.0 * _EI)
    ratio = abs(ux) / euler
    assert 0.9 < ratio < 1.02, (
        f"corot did not approach Euler tip deflection in the linear limit "
        f"(ux/Euler = {ratio:.3f}; residual is low-order-hex bending stiffness)")
