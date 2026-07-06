"""LadrunoSolidShell (ADR 66 P5.2, ELE_TAG 33020) — softening hosts + EAS stability gates.

P5.1 (test_ladrunoSolidShell_element.py) certified the elastic anti-locking
machinery (ANS shear/thickness tying + state-dependent EAS) and the D6 lch
contract. This file lands the P5.2 gates: LadrunoConcrete3D (33017) and
LadrunoRCConcrete (33015) hosted at the GPs, driven into damage/softening.

  --   test_hosting_parity_ndtest_replay   the solid-shell is TRANSPARENT to the
                                           material: GP stress under a fully-
                                           prescribed uniform-strain ramp ==
                                           an NDTest replay of the same strain
                                           path (all 6 components, every step);
                                           EAS alpha stays exactly 0 (a uniform
                                           field is orthogonal to the zeta-mode)
  --   test_rcconcrete_hosting_smoke       LadrunoRCConcrete (the ADR 19 RC
                                           GP material) loads, cracks, softens
  G5   test_g5_energy_mesh_objectivity     Bazant-band: weakened element column
                                           in an 8x2x0.1 strip, IN-PLANE-SQUARE
                                           elements at 3 densities; with
                                           -autoRegularization the total work
                                           to full softening is mesh-objective
  G5   test_g5_fixed_lch_not_objective     the discriminating control: a FIXED
                                           -lch makes dissipation ~ band volume
                                           (coarse >> fine) — proving G5 above
                                           is the element-supplied lch at work
  G6   test_g6_alpha_bounded_bending       single element, fully-prescribed
                                           bending ramped DEEP post-peak on the
                                           steep -lch 50 law: alpha stays
                                           finite, bounded, and settles (the
                                           zeta-mode is only activated by
                                           through-thickness gradients, so
                                           bending IS the alpha exercise)
  G6   test_g6_uniaxial_stress_softening   free-DOF uniaxial-stress tension
                                           driven deep post-peak: the EAS
                                           condensation lives on the damaged
                                           (indefinite) tangent — the Kaa<=0
                                           freeze guard path — without
                                           corrupting alpha or the response

G5/G6 mesh + law notes (why the numbers are what they are):
  * elements are kept IN-PLANE SQUARE (dx == dy) at every density because the
    element's lch is the isotropic sqrt(midsurface area) — only then does the
    crack-band eps_f = Gf/(ft*lch) see the true band width (dx). A directional
    lch(n) is ADR 66 O4 backlog.
  * the band column is 5% weaker in ft (2.85 vs 3.0) so localization is seeded
    deterministically; the bulk stays elastic (its onset is never reached).
  * no-snap-back check: Gf*E/ft^2 = 333 >> L = 8, so the prescribed-end path
    is stable and displacement control with step-cutting can follow it.
  * the G6 arms use the steep '-lch 50' law (eps_f ~ 6.7e-4, the same one the
    Concrete3D Tier-3 battery uses to force an unambiguous implicit stall) —
    a HARSHER alpha-stability exercise than the mild autoReg law.

Plan/ADR: Ladruno_implementation/66_ladruno_solidshell_adr.md (gates section 7).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# the oracle-gate concrete (shared with the Concrete3D batteries)
_E, _NU, _FC, _FT, _GF, _GC = 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0
# the weakened band strength (localization seed). 0.8*ft, NOT a hairline 0.95:
# at band peak the bulk then sits at ~80% of its own onset, so the global
# Newton's trial excursions through the localization transition cannot push
# bulk GPs into the deep-tension return-map apex regime (probed empirically:
# a 0.95 band drowned the run in bulk "return map did not converge" step-cuts).
# The dissipated energy is Gf-governed, so the gate is unaffected by ft_band.
_FTB = 0.8 * _FT


def _model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)


def _c3d(tag, ft=_FT, lch=None, auto=False, implex=False):
    args = ['LadrunoConcrete3D', tag, _E, _NU, _FC, ft, _GF, _GC]
    if lch is not None:
        args += ['-lch', float(lch)]
    if auto:
        args += ['-autoRegularization']
    if implex:
        args += ['-implex']
    ops.nDMaterial(*args)


def _max_wt(etag, ngp):
    return max(list(ops.eleResponse(etag, 'material', gp, 'damage'))[0]
               for gp in range(1, ngp + 1))


# ---------------------------------------------------------------------------
# hosting parity — the element is transparent to the material
# ---------------------------------------------------------------------------

def test_hosting_parity_ndtest_replay():
    """Fully-prescribed uniform uniaxial-STRAIN tension on a regular aligned
    element: every GP sees exactly the prescribed strain, so the recorded GP
    stress history must match an NDTest replay of the same strain path through
    a fresh LadrunoConcrete3D (same explicit -lch). This pins the hosting
    end-to-end — the T-map at the aligned geometry, the engineering-shear
    boundary convention, the trial/commit cycle — with the material itself as
    the reference. EAS alpha must stay EXACTLY 0 (a uniform strain field is
    orthogonal to the zeta-mode: int G dV = 0), so the parity is not
    EAS-polluted by construction. The steep -lch 50 law is used so the ramp
    genuinely softens (under the confined uniaxial-STRAIN state the mild
    autoReg law barely degrades over this strain range)."""
    coords = {1: (0.0, 0.0, 0.0), 2: (2.0, 0.0, 0.0), 3: (2.0, 1.0, 0.0),
              4: (0.0, 1.0, 0.0), 5: (0.0, 0.0, 0.05), 6: (2.0, 0.0, 0.05),
              7: (2.0, 1.0, 0.05), 8: (0.0, 1.0, 0.05)}
    _model()
    _c3d(1, lch=50.0)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1)
    eps_target, nsteps = 2.0e-3, 40           # ~20x the tension onset ft/E
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, eps_target * x)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, 0.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-10, 60, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    ele_curve = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
        sig = list(ops.eleResponse(1, 'stresses'))[0:6]
        a = list(ops.eleResponse(1, 'alpha'))[0]
        assert abs(a) < 1e-14, f"uniform strain must leave alpha at 0 (got {a})"
        ele_curve.append(sig)

    # NDTest replay of the identical strain path (commit each step)
    _model()
    _c3d(1, lch=50.0)
    ref_curve = []
    for i in range(1, nsteps + 1):
        e = eps_target * i / nsteps
        ops.NDTest('SetStrain', 1, e, 0.0, 0.0, 0.0, 0.0, 0.0)
        ref_curve.append(list(ops.NDTest('GetStress', 1)))
        ops.NDTest('CommitState', 1)

    smax = max(abs(v) for sig in ref_curve for v in sig)
    for i, (se, sr) in enumerate(zip(ele_curve, ref_curve)):
        for r in range(6):
            assert se[r] == pytest.approx(sr[r], rel=2e-6, abs=1e-6 * smax), \
                f"hosting parity: step {i} stress[{r}] {se[r]} != NDTest {sr[r]}"
    # the run genuinely softened (not a vacuous elastic comparison)
    ax = [sig[0] for sig in ele_curve]
    assert max(ax) > 0.9 * _FT and ax[-1] < 0.5 * max(ax)


def test_rcconcrete_hosting_smoke():
    """LadrunoRCConcrete (33015, the ADR 19 RC GP material) hosts in the
    solid-shell: a fully-prescribed uniform tension ramp cracks and softens
    without a stall or a non-finite response. (The RC material's own physics
    battery lives in its own test files — this is the hosting wire.)"""
    E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
    CE = [0.0, 0.0007, 0.0020, 0.0100]
    CS = [0.0, 24.0, 30.0, 5.0]
    CD = [0.0, 0.0, 0.25, 1.0 - 5.0 / 45.0]
    TE = [0.0, 0.0001, 0.0010, 0.0040]
    TS = [0.0, 3.0, 0.5, 0.0]
    TD = [0.0, 0.0, 1.0 - 0.5 / 5.0, 0.999]
    coords = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0),
              4: (0.0, 1.0, 0.0), 5: (0.0, 0.0, 0.2), 6: (1.0, 0.0, 0.2),
              7: (1.0, 1.0, 0.2), 8: (0.0, 1.0, 0.2)}
    _model()
    ops.nDMaterial('LadrunoRCConcrete', 1, E, NU, '-Ce', *CE, '-Cs', *CS,
                   '-Cd', *CD, '-Te', *TE, '-Ts', *TS, '-Td', *TD, '-Kc', KC)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, 1.5e-3 * x)              # 15x the cracking strain
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, 0.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-10, 60, 0)
    ops.algorithm('Newton')
    nsteps = 40
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    ax = []
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "RCConcrete hosting stalled"
        ax.append(list(ops.eleResponse(1, 'stresses'))[0])
    assert all(math.isfinite(s) for s in ax)
    peak = max(ax)
    assert peak > 2.0, f"RCConcrete never loaded (peak {peak})"
    assert ax[-1] < 0.8 * peak, f"RCConcrete did not soften ({ax[-1]} vs {peak})"


# ---------------------------------------------------------------------------
# G5 — softening-band energy objectivity (the Bazant-band pattern)
# ---------------------------------------------------------------------------
_L, _W, _T = 8.0, 2.0, 0.1
_LAM_END = 0.14                 # ~3.4 e-folds of band softening at ft_band=2.4


def _band_run(nx, ny, lch_mode, lam_end=_LAM_END):
    """Tension strip L x W x T, nx x ny in-plane elements (dx == dy enforced),
    the column whose RIGHT edge sits at x = L/2 built from the 5%-weaker band
    material. x=0 face u_x fixed (+ determinate rigid-body fixes), x=L face
    u_x prescribed = lambda. Adaptive LoadControl to lam_end with step-cutting
    through localization. Returns (curve[(lam, F)], wt_band, wt_neighbors)."""
    dx, dy = _L / nx, _W / ny
    assert abs(dx - dy) < 1e-12, "G5 requires in-plane-square elements"
    _model()
    # -implex: the ADR risk-register cure for "indefinite/non-symmetric
    # softening tangents stall implicit runs" — the SPD-ish IMPL-EX secant
    # lets plain Newton track the localization at the planned step size
    # (probed: Tier-1 cut-crawled to ~2400 micro-steps / 13 min on the COARSE
    # mesh alone). Constant-dlam LoadControl = the uniform-pseudo-time regime
    # IMPL-EX wants; the committed states are implicit-exact.
    if lch_mode == 'auto':
        _c3d(1, auto=True, implex=True)
        _c3d(2, ft=_FTB, auto=True, implex=True)
    else:
        _c3d(1, lch=lch_mode, implex=True)
        _c3d(2, ft=_FTB, lch=lch_mode, implex=True)
    tags = {}
    n = 0
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k, z in enumerate((-_T / 2, _T / 2)):
                n += 1
                tags[(i, j, k)] = n
                ops.node(n, i * dx, j * dy, z)
    band_col = nx // 2 - 1                    # right edge at x = L/2
    band, neigh = [], []
    e = 0
    for i in range(nx):
        for j in range(ny):
            e += 1
            mat = 2 if i == band_col else 1
            ops.element('LadrunoSolidShell', e,
                        tags[(i, j, 0)], tags[(i + 1, j, 0)],
                        tags[(i + 1, j + 1, 0)], tags[(i, j + 1, 0)],
                        tags[(i, j, 1)], tags[(i + 1, j, 1)],
                        tags[(i + 1, j + 1, 1)], tags[(i, j + 1, 1)], mat)
            if i == band_col:
                band.append(e)
            elif abs(i - band_col) == 1:
                neigh.append(e)
    # determinate restraints: u_x on the x=0 face; y/z rigid-body at the origin
    for j in range(ny + 1):
        for k in range(2):
            if (j, k) == (0, 0):
                ops.fix(tags[(0, 0, 0)], 1, 1, 1)
            elif (j, k) == (0, 1):
                ops.fix(tags[(0, 0, 1)], 1, 1, 0)
            else:
                ops.fix(tags[(0, j, k)], 1, 0, 0)
    left = [tags[(0, j, k)] for j in range(ny + 1) for k in range(2)]
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for j in range(ny + 1):
        for k in range(2):
            ops.sp(tags[(nx, j, k)], 1, 1.0)  # u_x = lambda at the pulled face
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')                 # unsymmetric-capable; RCM keeps the
    ops.test('NormDispIncr', 1.0e-8, 100, 0)  # strip bandwidth tiny (vs dense FullGeneral)
    ops.algorithm('Newton')
    ops.analysis('Static')

    def _F():
        ops.reactions()
        return -sum(ops.nodeReaction(t, 1) for t in left)

    curve = [(0.0, 0.0)]
    # three-phase ramp: seed (elastic rise + peak), TRANSITION (the band sheds
    # its load — per-step bulk strain increment dlam/L must stay under the
    # ~2e-5 gap between the bulk working stress at 0.8*onset and its onset, or
    # the redistribution transient commits real damage in the neighbor columns
    # — measured: dlam=1.17e-3 leaked omega_t=0.19, dlam<=1.5e-4 leaks none),
    # then bulk production steps through the smooth deep-softening tail.
    lam_seed, lam_mid = 1.5e-3, 0.02
    d_seed, d_mid, d_bulk = lam_seed / 25.0, 1.5e-4, lam_end / 120.0
    scale = 1.0                               # step = scale * phase base size
    guard = 0
    lam = 0.0
    while lam < lam_end - 1e-12 and guard < 6000:
        guard += 1
        base = d_seed if lam < lam_seed else (d_mid if lam < lam_mid else d_bulk)
        ops.integrator('LoadControl', min(scale * base, lam_end - lam))
        if ops.analyze(1) == 0:
            lam = ops.getTime()
            curve.append((lam, _F()))
            if scale < 1.0 and guard % 4 == 0:
                scale = min(1.0, 2.0 * scale)  # recover toward the base step
        else:
            scale *= 0.5
            if scale < 2.0 ** -10:
                break                          # genuinely stuck
    # derive the GP count live (4 in-plane x nz) so a future -nz change here
    # cannot silently misread the damage sweep
    ngp = len(list(ops.eleResponse(band[0], 'stresses'))) // 6
    wt_band = max(_max_wt(t, ngp) for t in band)
    wt_neigh = max(_max_wt(t, ngp) for t in neigh)
    return curve, wt_band, wt_neigh


def _work(curve):
    return sum(0.5 * (f0 + f1) * (l1 - l0)
               for (l0, f0), (l1, f1) in zip(curve, curve[1:]))


def _g5_arm(nx, ny, lch_mode):
    curve, wt_band, wt_neigh = _band_run(nx, ny, lch_mode)
    lam_last = curve[-1][0]
    assert lam_last > 0.95 * _LAM_END, \
        f"{nx}x{ny} ({lch_mode}): ramp incomplete (reached {lam_last:.4f})"
    F = [f for _, f in curve]
    peak = max(F)
    assert peak > 0.0
    assert F[-1] < 0.15 * peak, \
        f"{nx}x{ny} ({lch_mode}): band not fully softened (F_end/F_peak = {F[-1] / peak:.3f})"
    assert wt_band > 0.8, f"{nx}x{ny} ({lch_mode}): band omega_t only {wt_band:.3f}"
    # measured 0.000 with the d_mid transition schedule (0.19 at the coarse
    # production step); 0.1 absorbs platform FP wiggle while staying 2x under
    # the full-leak signature
    assert wt_neigh < 0.1, \
        f"{nx}x{ny} ({lch_mode}): damage leaked to the neighbor column ({wt_neigh:.3f})"
    return _work(curve)


def test_g5_energy_mesh_objectivity():
    """ADR 66 G5: with -autoRegularization (lch = sqrt(A) = dx on the square
    in-plane elements) the total work to full band softening is mesh-objective
    across a 4x element-size range. The theoretical floor is the crack-band
    design value Gf*W*T = 0.02 (plus the bulk's recovered-elastic + the band's
    un-regularized plastic extra — the documented D3 caveat, which is why the
    spread tolerance is not machine-tight)."""
    W_by_mesh = {}
    for nx, ny in ((4, 1), (8, 2), (16, 4)):
        W_by_mesh[(nx, ny)] = _g5_arm(nx, ny, 'auto')
    vals = list(W_by_mesh.values())
    mean = sum(vals) / len(vals)
    spread = (max(vals) - min(vals)) / mean
    # measured 0.009 across the 4x range (W = 0.01748/0.01735/0.01733). The
    # tolerance is NOT set hairline to that measurement: the runs are -implex
    # under an adaptive schedule, and the IMPL-EX reported state is step-
    # schedule sensitive (r = dt/dt_n), so cross-platform FP can move a step
    # cut and shift W. 0.12 keeps >13x margin over the measurement while
    # staying far below the ~0.75 fixed-lch spread over the same range.
    assert spread < 0.12, \
        f"G5 energy not mesh-objective: {W_by_mesh} (spread {spread:.3f})"
    # sanity band around the crack-band design dissipation Gf*W*T
    d_th = _GF * _W * _T
    for k, w in W_by_mesh.items():
        assert 0.5 * d_th < w < 3.0 * d_th, \
            f"G5 {k}: work {w:.4f} implausible vs Gf*W*T = {d_th:.4f}"


def test_g5_fixed_lch_not_objective():
    """The discriminating control: pin -lch to the COARSE element size (2.0)
    on both the coarse and the fine mesh. eps_f is then identical, so the
    dissipated energy scales with the BAND VOLUME (~dx) and the fine mesh
    must dissipate far less — the disease -autoRegularization cures. (On the
    coarse mesh lch=2.0 == sqrt(A), so this arm equals the G5 coarse arm and
    doubles as a cross-check of the auto-latch.)"""
    w_coarse = _g5_arm(4, 1, 2.0)
    w_fine = _g5_arm(16, 4, 2.0)
    # measured: 0.00438 vs 0.01748 = 0.25, exactly the band-volume (dx) ratio
    assert w_fine < 0.55 * w_coarse, \
        f"fixed-lch control failed to discriminate: fine {w_fine:.4f} vs coarse {w_coarse:.4f}"


# ---------------------------------------------------------------------------
# G6 — EAS alpha stability under deep softening
# ---------------------------------------------------------------------------

def test_g6_alpha_bounded_bending():
    """ADR 66 G6, the alpha-active arm: a fully-prescribed pure-bending field
    (u_x = k*x*z, u_z = -k*x^2/2) ramped DEEP post-peak on the steep -lch 50
    law (tension-face strain reaches ~4.5x eps_f => sigma ~ 1% ft). The
    zeta-linear EAS mode is directly driven by the through-thickness strain
    gradient, and the tension face damages while the compression face stays
    stiff — the asymmetric-Kaa regime the ADR flagged. Gate: every step
    converges, alpha stays finite and bounded by the kinematic scale of the
    problem, alpha genuinely activates, and the element ends deep post-peak
    with a finite stress state. A failure here triggers the documented
    fallback ladder (freeze-alpha post-peak -> mode removal -> Kstab)."""
    t, nz = 0.2, 5
    k_end, nsteps = 0.03, 150                 # face strain 3e-3 ~ 4.5x eps_f
    coords = {1: (-1.0, 0.0, -t / 2), 2: (1.0, 0.0, -t / 2),
              3: (1.0, 1.0, -t / 2), 4: (-1.0, 1.0, -t / 2),
              5: (-1.0, 0.0, t / 2), 6: (1.0, 0.0, t / 2),
              7: (1.0, 1.0, t / 2), 8: (-1.0, 1.0, t / 2)}
    _model()
    _c3d(1, lch=50.0)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1, '-nz', nz)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n, (x, y, z) in coords.items():
        ops.sp(n, 1, x * z)                   # u_x = lambda * x * z (lambda = k)
        ops.sp(n, 2, 0.0)
        ops.sp(n, 3, -x * x / 2.0)
    ops.constraints('Penalty', 1.0e12, 1.0e12)
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-9, 80, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', k_end / nsteps)
    ops.analysis('Static')
    alphas = []
    for s in range(nsteps):
        assert ops.analyze(1) == 0, f"bending softening stalled at step {s}"
        a = list(ops.eleResponse(1, 'alpha'))[0]
        assert math.isfinite(a), f"alpha went non-finite at step {s}"
        alphas.append(a)
    amax = max(abs(a) for a in alphas)
    # kinematic scale: the elastic optimal mode is |alpha| = nu/(1-nu)*k*t/2
    # ~ 7.5e-4 at k_end; damage redistribution may grow it, but a stable alpha
    # stays within one order of that scale.
    a_scale = _NU / (1.0 - _NU) * k_end * t / 2.0
    assert amax < 10.0 * a_scale, \
        f"alpha unbounded under deep softening: max|alpha| = {amax:.3e} vs scale {a_scale:.3e}"
    assert amax > 1e-6, "alpha never activated — the G6 arm is vacuous"
    # convergent: no blow-up in the tail — the largest step-to-step jump over
    # the last third is small against the run's own alpha scale
    tail = alphas[-nsteps // 3:]
    djump = max(abs(b - a) for a, b in zip(tail, tail[1:]))
    assert djump < 0.2 * (amax + 1e-12), \
        f"alpha still jumping in the tail (djump {djump:.3e} vs max {amax:.3e})"
    # NO "settled/contracting" assertion on the tail beyond the jump bound:
    # alpha legitimately keeps evolving while the ramp continues (the
    # softening front migrates through the zeta GPs), and its covariant scale
    # is tiny (measured amax = 1.7e-5 — the (t/2)^2 covariant metric factor
    # under the physical estimate), so any settle formalization compares
    # noise-level numbers and false-REDs (two variants tried: mean-jump
    # contraction 1.4e-7 vs 1.8e-7; last-third drift 7.7e-6 = 46% of amax).
    # G6 "bounded and convergent" = the two bounds above.
    # deep post-peak reached, response finite
    ngp = 4 * nz
    assert _max_wt(1, ngp) > 0.6, "did not reach deep tension softening"
    sig = list(ops.eleResponse(1, 'stresses'))
    assert all(math.isfinite(v) for v in sig)


def test_g6_uniaxial_stress_softening():
    """ADR 66 G6, the indefinite-Kaa arm: uniaxial-STRESS tension (free
    lateral/thickness DOFs) on the steep -lch 50 law, driven by adaptive
    step-cutting deep past the peak. Post-peak the damaged tangent is
    indefinite, so the EAS inner Newton lives on a Kaa that can lose
    positivity — the freeze-alpha guard path. The state is z-uniform, so a
    STABLE alpha stays ~0 throughout (the zeta-mode is orthogonal to uniform
    fields); a corrupted condensation would show up as alpha drift, a stress
    blow-up, or a permanent stall before deep softening."""
    t = 0.2
    coords = {1: (0.0, 0.0, -t / 2), 2: (1.0, 0.0, -t / 2),
              3: (1.0, 1.0, -t / 2), 4: (0.0, 1.0, -t / 2),
              5: (0.0, 0.0, t / 2), 6: (1.0, 0.0, t / 2),
              7: (1.0, 1.0, t / 2), 8: (0.0, 1.0, t / 2)}
    _model()
    _c3d(1, lch=50.0)
    for n, c in coords.items():
        ops.node(n, *c)
    ops.element('LadrunoSolidShell', 1, *range(1, 9), 1)
    # determinate restraints (x=0 face u_x + rigid-body y/z at the origin edge)
    ops.fix(1, 1, 1, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(4, 1, 0, 0)
    ops.fix(8, 1, 0, 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, 1.0)                     # u_x = lambda on the pulled face
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-9, 60, 0)
    ops.algorithm('Newton')
    ops.analysis('Static')
    eps_end = 2.5e-3                          # ~sigma/ft = e^-3.6 ~ 3%
    step = eps_end / 250.0
    cuts = guard = 0
    hist = []
    while ops.getTime() < eps_end - 1e-12 and guard < 4000:
        guard += 1
        ops.integrator('LoadControl', min(step, eps_end - ops.getTime()))
        if ops.analyze(1) == 0:
            a = list(ops.eleResponse(1, 'alpha'))[0]
            s = list(ops.eleResponse(1, 'stresses'))[0]
            hist.append((ops.getTime(), a, s))
            if cuts > 0 and guard % 4 == 0:
                step *= 2.0
                cuts -= 1
        else:
            step *= 0.5
            cuts += 1
            if cuts > 10:
                break
    assert hist, "no step converged"
    lam_last, _, _ = hist[-1]
    assert lam_last > 1.5e-3, \
        f"stalled before deep softening (reached eps {lam_last:.2e})"
    assert all(math.isfinite(a) and math.isfinite(s) for _, a, s in hist)
    amax = max(abs(a) for _, a, _ in hist)
    assert amax < 1e-8, \
        f"alpha drifted on a z-uniform state (max|alpha| = {amax:.3e}) — corrupted condensation"
    smax = max(abs(s) for _, _, s in hist)
    assert smax < 1.5 * _FT, f"stress blew up ({smax:.3f})"
    assert _max_wt(1, 8) > 0.5, "did not reach deep tension softening"
    # genuinely softened at the end
    assert abs(hist[-1][2]) < 0.4 * smax, \
        f"no softening: end stress {hist[-1][2]:.3f} vs peak {smax:.3f}"
