"""LadrunoSolidShell (ADR 64 P5.1, ELE_TAG 33020) — Zone-A element gates.

8-node ANS/EAS solid-shell: Dvorkin-Bathe ANS transverse-shear tying +
Betsch-Stein ANS thickness-strain (E33) corner tying + state-dependent
Simo-Rifai EAS-on-E33 (single zeta-mode, center transform, committed alpha).
'-formulation std' is the plain-displacement negative control (no ANS, no
EAS) used by the locking gates as the A/B disease exhibit.

  G1a  test_affine_patch_std_general      covariant->Cartesian pipeline exact on an
                                          ARBITRARILY distorted hex (std: plain B is
                                          exact for affine on any trilinear geometry)
  G1a  test_affine_patch_ans              ANS tying exact for the flat-face,
                                          in-plane-distorted geometry class
  G1a  test_thickness_patch_sigma33       the ADR 64 G1 headline: uniform sigma_33
                                          through a distorted mesh (catches a broken
                                          E33 tying / EAS orthogonality)
  G1b  test_interior_node_patch           2x2x1 mesh, boundary-prescribed affine,
                                          FREE interior nodes land on the affine field
  G2b  test_pure_bending_poisson_eas      EAS Poisson-thickness cure: single-element
                                          pure bending at nu=0.3 matches the exact
                                          sigma_zz=0 energy; std shows the parasitic
                                          sigma_33 plane-strain stiffening (+22%)
  G2a  test_slenderness_sweep_no_locking  cantilever t/L = 1e-1, 1e-2, 1e-3 vs
                                          Timoshenko: ratio flat in t (no locking)
  G2a  test_std_locks_thin                negative control: std collapses at t/L=1e-3
  G4   test_trapezoidal_mesh_bending      skewed zeta-fibers (trapezoidal x-z mesh):
                                          Betsch-Stein tying holds the answer
  --   test_parser_guards                 -nz/-quadz/-geom rejection surface
  --   test_database_roundtrip            sendSelf/recvSelf (committed alpha + mats)
  --   test_characteristic_length         in-plane-projected lch (never the thickness)

Deferred (ADR 64): pinched cylinder + Scordelis-Lo (G3, curved-mesh benchmark
PR), softening EAS stability G6 (P5.2), punching headline G8 (P5.3).
"""
import math

import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

E0 = 1000.0


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)


def _elastic(tag, E=E0, nu=0.3, rho=0.0):
    ops.nDMaterial('ElasticIsotropic', tag, E, nu, rho)


def _exact_stress(A, E, nu):
    """sigma (Voigt xx,yy,zz,xy,yz,zx, engineering shear) for eps = sym(A)."""
    eps = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    tr = eps[0][0] + eps[1][1] + eps[2][2]
    sig = [[lam * tr * (1 if i == j else 0) + 2 * mu * eps[i][j]
            for j in range(3)] for i in range(3)]
    return [sig[0][0], sig[1][1], sig[2][2], sig[0][1], sig[1][2], sig[0][2]]


def _affine_u(A, x, y, z):
    return (A[0][0] * x + A[0][1] * y + A[0][2] * z,
            A[1][0] * x + A[1][1] * y + A[1][2] * z,
            A[2][0] * x + A[2][1] * y + A[2][2] * z)


# an all-components-active affine gradient (exercises every strain channel,
# including eps_33 — the thickness patch rides inside every affine gate)
A_FULL = [[1.0e-3,  2.0e-4,  3.0e-4],
          [-1.0e-4, 5.0e-4,  1.0e-4],
          [2.0e-4, -3.0e-4,  4.0e-4]]


def _impose_affine_and_solve(node_coords, prescribed, A, nsteps=1):
    """sp-prescribe u = A*x on `prescribed` node tags; solve the rest (Penalty)."""
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n in prescribed:
        x, y, z = node_coords[n]
        ux, uy, uz = _affine_u(A, x, y, z)
        ops.sp(n, 1, ux)
        ops.sp(n, 2, uy)
        ops.sp(n, 3, uz)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    assert ops.analyze(nsteps) == 0


def _gp_stresses(tag, ngp):
    data = list(ops.eleResponse(tag, 'stresses'))
    assert len(data) == 6 * ngp
    return [data[6 * i:6 * i + 6] for i in range(ngp)]


# a single distorted hex: flat parallel faces (z = -t/2, +t/2), in-plane
# distorted, straight parallel zeta-fibers — the geometry class the ANS tying
# is exact on (file-header derivation)
def _flat_distorted_coords(t=0.4):
    xy = {1: (0.0, 0.0), 2: (1.1, -0.1), 3: (1.2, 0.9), 4: (-0.15, 1.05)}
    coords = {}
    for i in range(1, 5):
        x, y = xy[i]
        coords[i] = (x, y, -t / 2)
        coords[i + 4] = (x, y, t / 2)
    return coords


# an ARBITRARILY distorted hex (warped faces, tilted fibers) — std-only patch
def _general_distorted_coords():
    return {1: (0.0, 0.0, 0.0), 2: (1.0, 0.05, -0.05), 3: (1.1, 1.0, 0.08),
            4: (-0.05, 0.95, 0.03), 5: (0.05, -0.08, 0.5), 6: (1.05, 0.0, 0.55),
            7: (0.95, 1.05, 0.6), 8: (0.0, 1.0, 0.45)}


def _single_element(coords, mat_kw=None, ele_opts=()):
    _model()
    _elastic(1, **(mat_kw or {}))
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1, *ele_opts)


# ---------------------------------------------------------------------------
# G1a — affine consistency (patch) tests
# ---------------------------------------------------------------------------

def test_affine_patch_std_general():
    """std (plain covariant B) reproduces an affine field EXACTLY on an
    arbitrarily distorted hex — certifies the covariant->Cartesian pipeline
    (jacobian, buildT, covB) independent of the ANS substitutions."""
    coords = _general_distorted_coords()
    _single_element(coords, ele_opts=('-formulation', 'std'))
    _impose_affine_and_solve(coords, list(coords), A_FULL)
    exact = _exact_stress(A_FULL, E0, 0.3)
    for gp, sig in enumerate(_gp_stresses(1, 8)):
        for r in range(6):
            assert sig[r] == pytest.approx(exact[r], rel=1e-6, abs=1e-9), \
                f"std affine patch: GP {gp} component {r}"


def test_affine_patch_ans():
    """The full ANS+EAS recipe reproduces an affine field exactly on the
    flat-face in-plane-distorted geometry class (tying exact there)."""
    coords = _flat_distorted_coords()
    _single_element(coords)
    _impose_affine_and_solve(coords, list(coords), A_FULL)
    exact = _exact_stress(A_FULL, E0, 0.3)
    for gp, sig in enumerate(_gp_stresses(1, 8)):
        for r in range(6):
            assert sig[r] == pytest.approx(exact[r], rel=1e-6, abs=1e-9), \
                f"ans affine patch: GP {gp} component {r}"


def test_thickness_patch_sigma33():
    """ADR 64 G1 headline: a pure through-thickness stretch on the distorted
    flat-face element produces the exact UNIFORM sigma_33 at every GP — the
    test that catches a broken E33 tying or a non-orthogonal EAS mode."""
    A = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0e-3]]
    coords = _flat_distorted_coords()
    _single_element(coords)
    _impose_affine_and_solve(coords, list(coords), A)
    exact = _exact_stress(A, E0, 0.3)
    assert abs(exact[2]) > 0.0
    for gp, sig in enumerate(_gp_stresses(1, 8)):
        assert sig[2] == pytest.approx(exact[2], rel=1e-6), \
            f"thickness patch: sigma_33 at GP {gp}"
        assert sig[0] == pytest.approx(exact[0], rel=1e-6)  # lateral lam*eps33


def test_interior_node_patch():
    """2x2x1 flat-face in-plane-distorted mesh, affine prescribed on the
    BOUNDARY only: the free interior nodes must land on the affine field and
    every GP stress must be the exact constant.

    The patch field must be TRACTION-CONSISTENT with the free top/bottom
    faces (every node of a 1-element-thick patch sits on them): sigma.e_z = 0
    => eps_13 = eps_23 = 0 and eps_33 = -lam*(eps_11+eps_22)/(lam+2mu). The
    E33 channel is still exercised (eps_33 != 0); the full-traction GP-level
    exactness is owned by the fully-prescribed single-element patches above.
    (A full-A interior patch is ill-posed: the unloaded z-faces legitimately
    shed the sigma_33/sigma_i3 the affine field carries — verified against
    the numpy replica, where the traction-consistent field is exact to 1e-16
    for ans, ans-without-eas, and std alike.)"""
    t = 0.4
    nu = 0.3
    lam = E0 * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E0 / (2 * (1 + nu))
    e11, e22 = 1.0e-3, 5.0e-4
    e33 = -lam * (e11 + e22) / (lam + 2 * mu)
    A_PATCH = [[e11, 4.0e-4, 0.0],
               [-2.0e-4, e22, 0.0],
               [0.0, 0.0, e33]]
    # 3x3 in-plane grid, mid-edge + center nodes perturbed in-plane (same top
    # and bottom => flat faces, straight parallel fibers)
    grid = {(0, 0): (0.0, 0.0), (1, 0): (1.05, -0.08), (2, 0): (2.0, 0.0),
            (0, 1): (-0.06, 1.1), (1, 1): (0.93, 1.02), (2, 1): (2.08, 0.94),
            (0, 2): (0.0, 2.0), (1, 2): (1.1, 2.06), (2, 2): (2.0, 2.0)}
    coords, tags = {}, {}
    n = 0
    for (i, j), (x, y) in sorted(grid.items()):
        for k, z in enumerate((-t / 2, t / 2)):
            n += 1
            tags[(i, j, k)] = n
            coords[n] = (x, y, z)
    _model()
    _elastic(1)
    for tag, c in coords.items():
        ops.node(tag, *c)
    e = 0
    for i in range(2):
        for j in range(2):
            e += 1
            ops.element('LadrunoSolidShell', e,
                        tags[(i, j, 0)], tags[(i + 1, j, 0)],
                        tags[(i + 1, j + 1, 0)], tags[(i, j + 1, 0)],
                        tags[(i, j, 1)], tags[(i + 1, j, 1)],
                        tags[(i + 1, j + 1, 1)], tags[(i, j + 1, 1)], 1)
    interior = [tags[(1, 1, 0)], tags[(1, 1, 1)]]
    boundary = [tag for tag in coords if tag not in interior]
    _impose_affine_and_solve(coords, boundary, A_PATCH)

    for n in interior:
        x, y, z = coords[n]
        ux, uy, uz = _affine_u(A_PATCH, x, y, z)
        d = ops.nodeDisp(n)
        assert d[0] == pytest.approx(ux, rel=1e-6, abs=1e-12)
        assert d[1] == pytest.approx(uy, rel=1e-6, abs=1e-12)
        assert d[2] == pytest.approx(uz, rel=1e-6, abs=1e-12)

    exact = _exact_stress(A_PATCH, E0, 0.3)
    for etag in range(1, 5):
        for gp, sig in enumerate(_gp_stresses(etag, 8)):
            for r in range(6):
                assert sig[r] == pytest.approx(exact[r], rel=1e-5, abs=1e-8), \
                    f"interior patch: ele {etag} GP {gp} comp {r}"


# ---------------------------------------------------------------------------
# G2b — EAS Poisson-thickness cure (single element, closed form)
# ---------------------------------------------------------------------------

def _bending_energy(formulation, nu=0.3, t=0.2, kappa=1e-2):
    """Prescribe the pure-bending field u_x = k*x*z, u_y = 0, u_z = -k*x^2/2
    on a regular element (x,y in [-1,1]x[0,1], z in [-t/2,t/2]) and return the
    stored energy 1/2 f.u. Compatible E33 = 0 pointwise; eps_yy = 0; so:
      with EAS  : element relaxes sigma_zz -> 0 (the enhanced zeta-mode is the
                  EXACT optimal eps_33 = -nu/(1-nu)*k*z) => E* = E/(1-nu^2)
      without   : eps_33 stuck at 0 => plane strain in BOTH y,z
                  => E* = E(1-nu)/((1+nu)(1-2nu))  (~22% stiffer at nu=0.3)
    U_exact = 1/2 * E* * k^2 * Lx*Ly*t^3/12.
    """
    coords = {1: (-1.0, 0.0, -t / 2), 2: (1.0, 0.0, -t / 2),
              3: (1.0, 1.0, -t / 2), 4: (-1.0, 1.0, -t / 2),
              5: (-1.0, 0.0, t / 2), 6: (1.0, 0.0, t / 2),
              7: (1.0, 1.0, t / 2), 8: (-1.0, 1.0, t / 2)}
    _single_element(coords, mat_kw={'nu': nu},
                    ele_opts=('-formulation', formulation))
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, kappa * x * z)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, -kappa * x * x / 2.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    U = 0.0
    f = list(ops.eleForce(1))   # 24 entries, connectivity order (nodes 1..8)
    for a, n in enumerate(sorted(coords)):
        d = ops.nodeDisp(n)
        for k in range(3):
            U += 0.5 * f[3 * a + k] * d[k]
    return U


def test_pure_bending_poisson_eas():
    nu, t, kappa = 0.3, 0.2, 1e-2
    vol_factor = kappa ** 2 * (2.0 * 1.0 * t ** 3 / 12.0) / 2.0
    U_eas_exact = (E0 / (1 - nu ** 2)) * vol_factor
    mu = E0 / (2 * (1 + nu))
    # std carries TWO parasitic terms under the nodally-prescribed bending
    # field: (a) plane-strain-in-both sigma_33 stiffening on the bending term
    # (eps_33 stuck at 0), and (b) the classic parasitic shear gamma_xz = k*x
    # (the trilinear u_z can't follow -k*x^2/2 between nodes; the ANS tying at
    # xi=0 kills it exactly). U_par = 1/2 * mu * k^2 * int(x^2) * Ly * t.
    U_std_exact = (E0 * (1 - nu) / ((1 + nu) * (1 - 2 * nu))) * vol_factor \
        + 0.5 * mu * kappa ** 2 * (2.0 / 3.0) * 1.0 * t

    U_ans = _bending_energy('ans', nu=nu, t=t, kappa=kappa)
    # the EAS zeta-mode is the exact optimal thickness strain here — and the
    # ANS shear tying kills the parasitic gamma_xz — so the closed form is exact
    assert U_ans == pytest.approx(U_eas_exact, rel=2e-3), \
        "EAS Poisson-thickness cure: energy should match the sigma_zz=0 closed form"
    # sigma_33 ~ 0 at every GP (the physical signature of the cure)
    smax = max(abs(s[2]) for s in _gp_stresses(1, 8))
    sxx = max(abs(s[0]) for s in _gp_stresses(1, 8))
    assert smax < 1e-3 * sxx, "EAS should drive parasitic sigma_33 to ~0"
    # the enhanced parameter is alive
    a = list(ops.eleResponse(1, 'alpha'))
    assert abs(a[0]) > 0.0

    U_std = _bending_energy('std', nu=nu, t=t, kappa=kappa)
    assert U_std == pytest.approx(U_std_exact, rel=2e-3), \
        "std control: plane-strain bending + parasitic-shear closed form"
    assert U_std > 2.0 * U_ans, \
        "the disease must be visible: std locks (shear + Poisson-thickness)"


# ---------------------------------------------------------------------------
# G2a — slenderness sweep (transverse-shear locking) + G4 trapezoidal
# ---------------------------------------------------------------------------

def _cantilever(t, formulation='ans', nx=10, L=100.0, b=10.0, nu=0.0,
                skew=0.0, P=1.0e-3, axis='x'):
    """nx x 1 x 1 solid-shell cantilever spanning `axis` (x or y), clamped at
    the near end, tip shear P in z at the far end. skew != 0 tilts the
    zeta-fibers alternately (trapezoidal span-z mesh, the Betsch-Stein
    exhibit). Returns tip deflection / Timoshenko.

    axis='x' bends in the x-z plane -> transverse shear gamma_xz = gamma_13
    (the eta-tied row); axis='y' bends in the y-z plane -> gamma_yz = gamma_23
    (the xi-tied row). Running BOTH exercises both Dvorkin-Bathe tying rows."""
    _model()
    _elastic(1, nu=nu)
    ds = L / nx
    tags = {}
    n = 0
    for i in range(nx + 1):
        # alternate the fiber tilt on INTERIOR grid lines only
        s = 0.0 if (i == 0 or i == nx or skew == 0.0) else ((-1) ** i) * skew
        for j in range(2):
            for k, z in enumerate((-t / 2, t / 2)):
                n += 1
                tags[(i, j, k)] = n
                off = s if k == 1 else -s
                if axis == 'x':
                    # span x, width y in {0,b}
                    ops.node(n, i * ds + off, 0.0 if j == 0 else b, z)
                else:
                    # span y, width x in {b,0} (swapped so bottom face stays CCW
                    # with +z up -> detJ > 0)
                    ops.node(n, b if j == 0 else 0.0, i * ds + off, z)
    for i in range(nx):
        ops.element('LadrunoSolidShell', i + 1,
                    tags[(i, 0, 0)], tags[(i + 1, 0, 0)],
                    tags[(i + 1, 1, 0)], tags[(i, 1, 0)],
                    tags[(i, 0, 1)], tags[(i + 1, 0, 1)],
                    tags[(i + 1, 1, 1)], tags[(i, 1, 1)],
                    1, '-formulation', formulation)
    for j in range(2):
        for k in range(2):
            ops.fix(tags[(0, j, k)], 1, 1, 1)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    tip = [tags[(nx, j, k)] for j in range(2) for k in range(2)]
    for nt in tip:
        ops.load(nt, 0.0, 0.0, P / 4.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    dz = sum(ops.nodeDisp(nt)[2] for nt in tip) / 4.0
    I = b * t ** 3 / 12.0
    G = E0 / (2 * (1 + nu))
    delta = P * L ** 3 / (3 * E0 * I) + 6.0 * P * L / (5.0 * G * b * t)
    return dz / delta


def test_transverse_shear_both_axes():
    """The x-spanning slenderness gate only activates gamma_13 (the eta-tied
    Dvorkin-Bathe row); a y-spanning cantilever activates gamma_23 (the xi-tied
    row). A wrong/degenerate gamma_23 tying would locally shear-lock only the
    y-cantilever, so require the SAME thin-shell no-lock behaviour on both."""
    for ax in ('x', 'y'):
        r = _cantilever(0.1, axis=ax)
        assert 0.90 < r < 1.05, \
            f"ans {ax}-cantilever t/L=1e-3: ratio {r:.4f} outside [0.90, 1.05]"
        r_std = _cantilever(0.1, formulation='std', axis=ax)
        assert r_std < 0.5, \
            f"std {ax}-cantilever should lock (control, got {r_std:.4f})"


def test_slenderness_sweep_no_locking():
    """The ADR 64 G2 gate: tip-deflection ratio vs Timoshenko stays flat over
    t/L = 1e-1 -> 1e-3 (the locking failure mode collapses the thin ratios)."""
    ratios = {t: _cantilever(t) for t in (10.0, 1.0, 0.1)}
    for t, r in ratios.items():
        assert 0.90 < r < 1.05, f"ans t={t}: ratio {r:.4f} outside [0.90, 1.05]"
    vals = list(ratios.values())
    spread = max(vals) - min(vals)
    assert spread < 0.05, \
        f"ratio must be flat in t (spread {spread:.4f}) — locking collapses it"


def test_std_locks_thin():
    """Negative control (the disease exhibit): the plain-displacement std
    formulation shear-locks at t/L=1e-3."""
    r = _cantilever(0.1, formulation='std')
    assert r < 0.5, f"std at t/L=1e-3 should lock hard (got ratio {r:.4f})"


def test_trapezoidal_mesh_bending():
    """G4 (Betsch-Stein): alternately tilted zeta-fibers (trapezoidal x-z
    mesh, ~27 deg fiber tilt). The ANS thickness tying holds the bending
    answer to the known distortion residual (~15% stiffening at this tilt —
    the in-plane trapezoid coupling is not fully curable; Betsch-Stein removes
    the curvature-thickness part); std collapses outright."""
    t = 1.0
    r_ans = _cantilever(t, skew=0.25 * t)
    assert 0.80 < r_ans < 1.05, \
        f"ans on the trapezoidal mesh: ratio {r_ans:.4f} outside [0.80, 1.05]"
    r_std = _cantilever(t, formulation='std', skew=0.25 * t)
    assert r_std < 0.6, f"std trapezoidal control should lock (got {r_std:.4f})"
    assert r_ans > 2.0 * r_std, \
        f"discrimination: ans ({r_ans:.3f}) must beat std ({r_std:.3f}) by 2x"


# ---------------------------------------------------------------------------
# parser guards / serialization / lch
# ---------------------------------------------------------------------------

def _try_element(*opts):
    coords = _flat_distorted_coords()
    _model()
    _elastic(1)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1, *opts)


def test_parser_guards():
    # good: defaults, explicit lobatto nz=5, std
    _try_element()
    _try_element('-nz', 5)
    _try_element('-nz', 3, '-quadz', 'lobatto')
    _try_element('-formulation', 'std')
    # bad: nz out of range / lobatto nz=2 (trapezoid) / gauss nz=7 / -geom corot
    for bad in (('-nz', 9), ('-nz', 1),
                ('-nz', 2, '-quadz', 'lobatto'),
                ('-nz', 7, '-quadz', 'gauss'),
                ('-geom', 'corot'), ('-geom', 'finite')):
        with pytest.raises(Exception):
            _try_element(*bad)


def test_database_roundtrip():
    """PLUMBING gate: nz=3 (numGP=12) send/recv wire layout + broker dispatch
    (classTag -> constructor) survive a FE_Datastore save/restore without
    crashing and reproduce the per-GP material stresses. NOTE this uses the
    shared helper, which RE-SOLVES on the restore side and an affine state has
    alpha==0 (EAS orthogonality) — so it does NOT prove alphaCommit survives.
    That is what test_database_roundtrip_alpha_survives below does."""
    def build():
        coords = _flat_distorted_coords()
        _single_element(coords, ele_opts=('-nz', 3))
        _impose_affine_and_solve(coords, list(coords), A_FULL)
    database_roundtrip(build, probe_nodes=[6], ndf=3,
                       probe_fn=lambda: list(ops.eleResponse(1, 'stresses')))


def _build_bending_committed(nz=3):
    """Regular element under the pure-bending prescribed field (u_x=k*x*z,
    u_z=-k*x^2/2), analyzed+committed. Leaves element tag 1 with a NON-ZERO
    committed alpha (the EAS zeta-mode activates under bending). Returns the
    node coords dict."""
    t, kappa = 0.2, 1e-2
    coords = {1: (-1.0, 0.0, -t / 2), 2: (1.0, 0.0, -t / 2),
              3: (1.0, 1.0, -t / 2), 4: (-1.0, 1.0, -t / 2),
              5: (-1.0, 0.0, t / 2), 6: (1.0, 0.0, t / 2),
              7: (1.0, 1.0, t / 2), 8: (-1.0, 1.0, t / 2)}
    _single_element(coords, ele_opts=('-nz', nz))
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, kappa * x * z)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, -kappa * x * x / 2.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    return coords


def test_database_roundtrip_alpha_survives():
    """RIGOROUS serialization gate: committed alpha must be reconstructed by
    recvSelf on a skeleton that is NEVER re-solved. Build a non-zero-alpha
    bending state, save, wipe, rebuild ONLY the model+nodes+element (no loads,
    no analyze), restore, and read the element state directly. Because no
    analysis runs on the restore side, a recvSelf that dropped alphaCommit (or
    the per-GP material state) would read back alpha=0 / zero stress and fail —
    the re-solve masking of the shared helper is removed."""
    import os
    import tempfile
    with tempfile.TemporaryDirectory(prefix='ladruno_rt_alpha_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'ss_alpha')
        _build_bending_committed(nz=3)
        alpha_before = list(ops.eleResponse(1, 'alpha'))
        stress_before = list(ops.eleResponse(1, 'stresses'))
        assert abs(alpha_before[0]) > 1e-8, "precondition: alpha must be non-zero"

        try:
            ops.database('File', dbpath)
        except Exception as exc:  # noqa: BLE001
            pytest.skip(f"database() unsupported: {exc}")
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip("database save failed on this build")

        ops.wipe()
        # SKELETON only — no sp, no pattern, no analysis, no analyze()
        coords = {1: (-1.0, 0.0, -0.1), 2: (1.0, 0.0, -0.1),
                  3: (1.0, 1.0, -0.1), 4: (-1.0, 1.0, -0.1),
                  5: (-1.0, 0.0, 0.1), 6: (1.0, 0.0, 0.1),
                  7: (1.0, 1.0, 0.1), 8: (-1.0, 1.0, 0.1)}
        _single_element(coords, ele_opts=('-nz', 3))
        ops.database('File', dbpath)
        ops.restore(1)
        alpha_after = list(ops.eleResponse(1, 'alpha'))
        stress_after = list(ops.eleResponse(1, 'stresses'))
        ops.wipe()

    assert alpha_after[0] == pytest.approx(alpha_before[0], rel=1e-10, abs=1e-14), \
        f"alphaCommit not restored: {alpha_before[0]} -> {alpha_after[0]}"
    assert len(stress_after) == len(stress_before)
    for i, (b, a) in enumerate(zip(stress_before, stress_after)):
        assert a == pytest.approx(b, rel=1e-9, abs=1e-9), \
            f"per-GP material stress[{i}] not restored: {b} -> {a}"


def _stiff_matrix(resp):
    flat = list(ops.eleResponse(1, resp))
    assert len(flat) == 24 * 24, f"{resp} returned {len(flat)} entries"
    return [flat[24 * i:24 * i + 24] for i in range(24)]


def test_initial_equals_tangent_stiff_elastic():
    """getInitialStiff is a SEPARATE hand-written condensed-initial-tangent
    branch (formANS useInitialTangent, its own Kaa/Kda/Kad accumulation and
    its own Kaa<=0 skip). For a linear-elastic material at the undeformed
    (alpha=0) state it MUST equal getTangentStiff entry-for-entry — otherwise
    the two EAS condensations disagree and Newton starts from a wrong Jacobian.
    Compared directly via the diagnostic stiff/stiffInitial responses."""
    coords = _flat_distorted_coords()
    _single_element(coords)              # ans, elastic, undeformed
    Kt = _stiff_matrix('stiff')
    Ki = _stiff_matrix('stiffInitial')
    scale = max(abs(Kt[i][j]) for i in range(24) for j in range(24))
    assert scale > 0.0
    for i in range(24):
        for j in range(24):
            assert abs(Kt[i][j] - Ki[i][j]) <= 1e-9 * scale, \
                f"K_tangent != K_initial at ({i},{j}): {Kt[i][j]} vs {Ki[i][j]}"
    # the elastic condensed tangent must be symmetric (Simo-Rifai EAS with a
    # symmetric material is symmetric; a broken Kda/Kad transpose would show)
    asym = max(abs(Kt[i][j] - Kt[j][i])
               for i in range(24) for j in range(24))
    assert asym <= 1e-9 * scale, f"elastic tangent not symmetric (asym {asym:.3e})"
    # and the std formulation likewise
    _single_element(coords, ele_opts=('-formulation', 'std'))
    Kt2 = _stiff_matrix('stiff')
    Ki2 = _stiff_matrix('stiffInitial')
    s2 = max(abs(Kt2[i][j]) for i in range(24) for j in range(24))
    for i in range(24):
        for j in range(24):
            assert abs(Kt2[i][j] - Ki2[i][j]) <= 1e-9 * s2, \
                f"std K_tangent != K_initial at ({i},{j})"


def _column(rho=0.0):
    """4-element solid-shell column along z (height 10, 1x1 section), base
    fixed. Returns the top-node tag list (loaded face). Leaves the model built
    with a Constant step lateral load in x on the top face."""
    _model()
    _elastic(1, E=1000.0, nu=0.0, rho=rho)
    nz_ele = 4
    h = 10.0 / nz_ele
    tags = {}
    n = 0
    for k in range(nz_ele + 1):
        for (ii, jj, xx, yy) in ((0, 0, 0.0, 0.0), (1, 0, 1.0, 0.0),
                                 (1, 1, 1.0, 1.0), (0, 1, 0.0, 1.0)):
            n += 1
            tags[(ii, jj, k)] = n
            ops.node(n, xx, yy, k * h)
    for e in range(nz_ele):
        ops.element('LadrunoSolidShell', e + 1,
                    tags[(0, 0, e)], tags[(1, 0, e)], tags[(1, 1, e)], tags[(0, 1, e)],
                    tags[(0, 0, e + 1)], tags[(1, 0, e + 1)],
                    tags[(1, 1, e + 1)], tags[(0, 1, e + 1)], 1)
    for corner in ((0, 0), (1, 0), (1, 1), (0, 1)):
        ops.fix(tags[(corner[0], corner[1], 0)], 1, 1, 1)
    top = [tags[(a, b, nz_ele)] for (a, b) in ((0, 0), (1, 0), (1, 1), (0, 1))]
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for nt in top:
        ops.load(nt, 0.25, 0.0, 0.0)
    return top


def _column_tip_static():
    top = _column(rho=0.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-10, 30, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    return abs(ops.nodeDisp(top[0])[0])


def _column_tip_dynamic(betaK, nsteps=240, dt=0.05):
    """Newmark transient (avg-accel) of the column with ELEMENT mass (rho) and
    stiffness-proportional Rayleigh betaK. Returns peak |tip x-disp| over the
    window (tuned so ~1 period is resolved -> a clean ~2x step-load overshoot)."""
    top = _column(rho=2.0)
    if betaK != 0.0:
        ops.rayleigh(0.0, betaK, 0.0, 0.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-10, 30, 0)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    peak = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        peak = max(peak, abs(ops.nodeDisp(top[0])[0]))
    return peak


def test_dynamic_rayleigh_preserves_inertia():
    """REGRESSION for the critical resid-clobber bug: getResistingForceIncInertia
    accumulates f_int + M*a into the shared static resid, then adds the Rayleigh
    force — but getRayleighDampingForces() re-enters getTangentStiff()/
    getInitialStiff() -> formANS() -> resid.Zero(), silently wiping the inertia
    term (without the local snapshot). A transient with stiffness-proportional
    betaK would then drop element inertia and collapse toward the quasi-static
    response (no overshoot) while Newton still converges.

    Self-validating: (1) the undamped run MUST show the ~2x step-load dynamic
    overshoot vs the static deflection — proving real dynamics are resolved and
    we are in the regime where the bug is visible; (2) adding tiny betaK must
    change the peak by <5% (light damping). Under the bug, betaK!=0 loses
    inertia so its peak collapses to ~1x static (~half the undamped peak), a
    >40% divergence that trips assertion (2)."""
    static = _column_tip_static()
    peak_undamped = _column_tip_dynamic(0.0)
    overshoot = peak_undamped / static
    assert 1.5 < overshoot < 2.1, (
        f"undamped peak/static = {overshoot:.3f}; expected the ~2x step-load "
        "overshoot — the test is not resolving real dynamics, so it cannot "
        "discriminate the inertia bug")
    peak_lightdamp = _column_tip_dynamic(1.0e-5)
    rel = abs(peak_lightdamp - peak_undamped) / peak_undamped
    assert rel < 0.05, (
        f"tiny betaK changed the peak dynamic response by {rel:.1%} "
        f"({peak_undamped:.6g} -> {peak_lightdamp:.6g}, static {static:.6g}) — "
        "inertia is being dropped from getResistingForceIncInertia "
        "(the resid-clobber bug)")


def test_setresponse_forwarding():
    """The 'material $gp' / 'stresses' / 'strains' / 'alpha' setResponse
    branches must all return non-empty (a missing branch = silently empty
    recorder). Covers the D4 defect class via the _testbed smoke helper."""
    from _testbed.roundtrip import setresponse_smoke
    coords = _flat_distorted_coords()
    _single_element(coords, ele_opts=('-nz', 3))
    _impose_affine_and_solve(coords, list(coords), A_FULL)
    setresponse_smoke(1, [('stresses',), ('strains',), ('alpha',),
                          ('material', 1, 'stress'),
                          ('material', 12, 'stress')])


def _tension_softening_curve(mat_flags):
    """Fully-prescribed uniform tension (eps_xx ramped, lateral held) on a
    2 x 1 x 0.05 solid-shell with LadrunoConcrete3D; returns the sigma_xx
    history at GP 1. Fully prescribed => no limit point, plain LoadControl."""
    coords = {1: (0.0, 0.0, 0.0), 2: (2.0, 0.0, 0.0), 3: (2.0, 1.0, 0.0),
              4: (0.0, 1.0, 0.0), 5: (0.0, 0.0, 0.05), 6: (2.0, 0.0, 0.05),
              7: (2.0, 1.0, 0.05), 8: (0.0, 1.0, 0.05)}
    _model()
    ops.nDMaterial('LadrunoConcrete3D', 1, 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0,
                   *mat_flags)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1)
    eps_target = 2.0e-3   # ~20x the tension onset eps0 = ft/E = 1e-4
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, eps_target * x)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, 0.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-10, 40, 0)
    ops.algorithm('Newton')
    nsteps = 40
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    curve = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        curve.append(list(ops.eleResponse(1, 'stresses'))[0])
    return curve


def test_characteristic_length_in_plane():
    """ADR 64 D6: getCharacteristicLength() = sqrt(midsurface area) — the
    IN-PLANE projected size, NEVER the thickness. Probed quantitatively via
    LadrunoConcrete3D's crack-band latch: with -autoRegularization the
    softening backbone encodes the element-supplied lch, so the response must
    MATCH an explicit '-lch sqrt(2)' run and DIFFER from an '-lch 0.05'
    (thickness) run. Element: 2x1 in-plane, t=0.05 => sqrt(A)=sqrt(2)."""
    auto = _tension_softening_curve(['-autoRegularization'])
    exact = _tension_softening_curve(['-lch', math.sqrt(2.0)])
    wrong = _tension_softening_curve(['-lch', 0.05])

    for a, b in zip(auto, exact):
        assert a == pytest.approx(b, rel=1e-9, abs=1e-12), \
            "auto-latched lch must equal sqrt(midsurface area) = sqrt(2)"
    # non-tautology: the thickness-lch curve must be clearly different
    gap = max(abs(a - w) for a, w in zip(auto, wrong))
    peak = max(abs(v) for v in auto)
    assert gap > 0.05 * peak, \
        "an element reporting the THICKNESS as lch would have matched 'wrong'"
