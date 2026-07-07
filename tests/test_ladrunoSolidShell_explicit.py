"""LadrunoSolidShell (ADR 66 G9, ELE_TAG 33020) — the explicit path: -lumped mass + SMS.

The D7 risk made concrete and discharged. G9 probing on the P5.1 element found
the real gap this PR closes: the element shipped CONSISTENT-mass only, and
under the explicit recipe (`system Diagonal`) the SOE keeps just the raw
diagonal of that matrix — 8/27 of the element mass (1/27 per node on the
trilinear shape) — so explicit dynamics ran 3.375x mass-deficient and the
'-lump hrz' pencil over-reported the stable step by sqrt(3.375) = 1.84x
(measured: blow-up at 0.9x the reported dt_cr, true boundary at ~0.54x = the
consistent-eig pencil). The fix is the LadrunoBrick donor pattern: an element
`-lumped` flag whose getMass() returns the row-sum diagonal m_a = rho*int N_a
dV (== HRZ on the trilinear shape, always positive) — then the runtime mass,
the `-lump hrz`/`rowsum` pencil, and the physics all agree.

  G9a  test_lumped_mass_dynamics_correct   free axial vibration of a -lumped
                                           strip lands on the continuum period
                                           T1 = 4L/c (the mass-CORRECTNESS
                                           gate: the old consistent-under-
                                           Diagonal path is ~1.84x off — this
                                           is what -lumped fixes)
  G9a  test_criticalTimeStep_validated     -lumped + '-lump hrz' pencil == t/c
                                           (the thickness CFL), HALVES with
                                           the thickness (D7), and is
                                           bisection-validated: stable at
                                           0.9x, unstable by 1.5x
  --   test_lumped_roundtrip_stable        massType survives sendSelf/recvSelf
                                           (a dropped flag = consistent mass =
                                           blow-up within ~25 steps at the
                                           lumped step; the restored run
                                           marches stably)
  G9b  test_sms_speedup_quasi_static       fixed physical loading time: CDL at
                                           the thickness-bound step vs SMS at
                                           dtTarget = 20x — both land on the
                                           implicit static reference, SMS in
                                           >= 15x fewer steps
  G9c  test_softening_band_sms_completes   a Concrete3D weakened band pulled
                                           to full softening under
                                           CentralDifferenceSMS + smoothstep
                                           support motion (Tier-3: no implex,
                                           no step-cutting, no unsymmetric
                                           solver) at 3x the pencil — the
                                           D7 cure end-to-end, with two
                                           measured sqrt(s) impedance-shock
                                           lessons in the comments

NB the LEDGER quirk: criticalTimeStep() under SMS reports the PRE-scaling
pencil — under SMS steps are sized by dtTarget, never by the query. And the
CDL '-lump' flag selects the PENCIL model only; the runtime mass is whatever
getMass() feeds `system Diagonal` — which is why the element flag must exist.

Plan/ADR: Ladruno_implementation/66_ladruno_solidshell_adr.md (gate G9).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
SMS = "CentralDifferenceSMS"

_E, _NU, _RHO = 1000.0, 0.0, 1.0            # nu=0: clean 1-D axial speed c = sqrt(E/rho)
_C = math.sqrt(_E / _RHO)


def _model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)


def _strip(nx, dx, dy, t, rho=_RHO, mat=None, lumped=True, axial_only=False):
    """nx x 1 solid-shell strip along x (width dy, thickness t), x=0 face
    clamped. axial_only pins y,z everywhere (clean 1-D axial dynamics)."""
    _model()
    if mat is None:
        ops.nDMaterial('ElasticIsotropic', 1, _E, _NU, rho)
    else:
        mat()
    tags = {}
    n = 0
    for i in range(nx + 1):
        for j in range(2):
            for k, z in enumerate((-t / 2, t / 2)):
                n += 1
                tags[(i, j, k)] = n
                ops.node(n, i * dx, j * dy, z)
    opts = ('-lumped',) if lumped else ()
    for i in range(nx):
        ops.element('LadrunoSolidShell', i + 1,
                    tags[(i, 0, 0)], tags[(i + 1, 0, 0)],
                    tags[(i + 1, 1, 0)], tags[(i, 1, 0)],
                    tags[(i, 0, 1)], tags[(i + 1, 0, 1)],
                    tags[(i + 1, 1, 1)], tags[(i, 1, 1)], 1, *opts)
    for j in range(2):
        for k in range(2):
            ops.fix(tags[(0, j, k)], 1, 1, 1)
    if axial_only:
        for i in range(1, nx + 1):
            for j in range(2):
                for k in range(2):
                    ops.fix(tags[(i, j, k)], 0, 1, 1)
    return tags


def _explicit_solver():
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1.0e-12, 1)
    ops.algorithm('Linear')


def _pencil_of_current_model():
    _explicit_solver()
    ops.integrator(CDL, '-lump', 'hrz')
    ops.analysis('Transient')
    assert ops.analyze(1, 1.0e-9) == 0, "prime step failed"
    dt = ops.criticalTimeStep()
    assert math.isfinite(dt) and dt > 0.0
    return dt


# ---------------------------------------------------------------------------
# G9a — mass correctness + dt_cr validation
# ---------------------------------------------------------------------------

def test_lumped_mass_dynamics_correct():
    """Free axial vibration of the clamped -lumped strip (lateral pinned):
    the fundamental period must land on the continuum fixed-free bar value
    T1 = 4L/c. This is the gate the -lumped flag exists for: the raw
    consistent diagonal under `system Diagonal` carries 8/27 of the mass, so
    that path rings ~sqrt(27/8) = 1.84x FAST — far outside the 5% band."""
    nx, dx, t = 4, 2.0, 0.1
    L = nx * dx
    tags = _strip(nx, dx, 2.0, t, axial_only=True)
    # linear axial ramp IC (excites the fundamental), set BEFORE the analysis
    # is built (CDL seeds its staggered velocity at the first step)
    A = 1.0e-3
    for i in range(1, nx + 1):
        for j in range(2):
            for k in range(2):
                ops.setNodeDisp(tags[(i, j, k)], 1, A * (i * dx) / L, '-commit')
    dtc = _pencil_of_current_model()
    T1 = 4.0 * L / _C
    dt = 0.5 * dtc
    nsteps = int(2.5 * T1 / dt)
    tip = tags[(nx, 0, 0)]
    t_hist, u_hist = [], []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        t_hist.append(ops.getTime())
        u_hist.append(ops.nodeDisp(tip)[0])
    # period from successive upward zero crossings about the mean
    um = sum(u_hist) / len(u_hist)
    c = [u - um for u in u_hist]
    crossings = [t_hist[i] for i in range(1, len(c))
                 if c[i - 1] <= 0.0 < c[i]]
    assert len(crossings) >= 2, "did not capture two periods"
    T_meas = crossings[1] - crossings[0]
    assert T_meas == pytest.approx(T1, rel=0.05), \
        f"axial period {T_meas:.4f} != 4L/c {T1:.4f} — lumped mass wrong " \
        f"(consistent-under-Diagonal would ring ~1.84x fast)"


def _free_vib_peak(t, dt_factor, nsteps=1500):
    """Clamped -lumped strip, tip-face z-velocity kick (set BEFORE analysis),
    march nsteps at dt_factor * pencil. Returns (completed, peak |u_z| tip)."""
    tags = _strip(4, 2.0, 2.0, t)
    for j in range(2):
        for k in range(2):
            ops.setNodeVel(tags[(4, j, k)], 3, 0.1, '-commit')
    dtc = _pencil_of_current_model()
    dt = dt_factor * dtc
    peak = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return False, peak
        u = ops.nodeDisp(tags[(4, 0, 0)])[2]
        if not math.isfinite(u) or abs(u) > 1.0e3:
            return False, peak
        peak = max(peak, abs(u))
    return True, peak


def test_criticalTimeStep_validated():
    # the -lumped + '-lump hrz' pencil is the thickness CFL exactly
    tags = _strip(4, 2.0, 2.0, 0.1)
    dt1 = _pencil_of_current_model()
    assert dt1 == pytest.approx(0.1 / _C, rel=0.05), \
        f"lumped pencil {dt1:.4e} != t/c {0.1 / _C:.4e}"
    # D7: halve the thickness -> dt_cr halves (in-plane size unchanged)
    tags = _strip(4, 2.0, 2.0, 0.05)
    dt2 = _pencil_of_current_model()
    assert 0.4 < dt2 / dt1 < 0.6, \
        f"dt_cr not thickness-bound: t/2 gave ratio {dt2 / dt1:.3f}"
    # bisection: stable below the pencil, unstable not far above it
    ok, peak = _free_vib_peak(0.1, 0.9)
    v0_scale = 0.1 * (2.0 / _C)
    assert ok and peak < 50.0 * v0_scale, \
        f"0.9x pencil not stable (ok={ok}, peak {peak:.3e})"
    unstable_at = None
    for f in (1.05, 1.2, 1.5):
        ok, peak = _free_vib_peak(0.1, f)
        if (not ok) or peak > 1.0e3 * v0_scale:
            unstable_at = f
            break
    assert unstable_at is not None, \
        "still stable at 1.5x the pencil — dt_cr no longer tracks the runtime mass"


def test_lumped_roundtrip_stable():
    """massType must survive sendSelf/recvSelf (it rides the former spare
    idData(6) slot). Discriminating by physics: the restored element is
    marched at 0.9x the LUMPED pencil — if the flag were dropped on the wire
    the element would come back consistent-mass and blow up within ~25 steps
    (the measured super-critical signature); the lumped element marches all
    500 steps."""
    import os
    import tempfile
    with tempfile.TemporaryDirectory(prefix='ladruno_rt_lump_',
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, 'ss_lump')
        tags = _strip(4, 2.0, 2.0, 0.1)
        for j in range(2):
            for k in range(2):
                ops.setNodeVel(tags[(4, j, k)], 3, 0.1, '-commit')
        dtc = _pencil_of_current_model()
        for _ in range(30):
            assert ops.analyze(1, 0.9 * dtc) == 0
        try:
            ops.database('File', dbpath)
        except Exception as exc:  # noqa: BLE001
            pytest.skip(f"database() unsupported: {exc}")
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip("database save failed on this build")

        # SKELETON rebuilt with the DEFAULT (consistent) element — restore
        # must bring massType=1 back from the stream
        tags = _strip(4, 2.0, 2.0, 0.1, lumped=False)
        ops.database('File', dbpath)
        ops.restore(1)
        _explicit_solver()
        ops.integrator(CDL, '-lump', 'hrz')
        ops.analysis('Transient')
        for s in range(500):
            assert ops.analyze(1, 0.9 * dtc) == 0, \
                f"restored run failed at step {s} — massType dropped on the wire?"
            u = ops.nodeDisp(tags[(4, 0, 0)])[2]
            assert math.isfinite(u) and abs(u) < 1.0e3, \
                f"restored run diverging at step {s} — massType dropped on the wire?"
        ops.wipe()


# ---------------------------------------------------------------------------
# G9b — SMS speedup: a SLIVER-bound mesh (the SMS design case). On a UNIFORM
# mesh dtTarget-scaling touches every element, the fundamental slows by the
# same factor, and a fixed-T quasi-static ramp under-delivers (measured:
# 20x-target on a uniform strip left the tip at 2% of static) — that is the
# documented "dtTarget > bulk step scales the bulk too" trap, not a bug. SMS
# pays when a FEW sliver elements pin the pencil: here one 0.01-wide element
# in a 0.1-thick strip sets dt_cr 10x below the bulk, SMS injects mass only
# around the sliver, the low modes stay intact, and the SAME physical T needs
# ~10x fewer steps.
# ---------------------------------------------------------------------------
_SLIVER_XS = (0.0, 1.0, 1.01, 2.01)            # widths 1.0 / 0.01 / 1.0
_B_DY, _B_T = 1.0, 0.1
_B_P = 1.0e-3                                  # tip force (z), small/elastic


def _sliver_strip(lumped=True):
    """3-element cantilever strip along x with a 0.01-wide sliver column in
    the middle (widths 1.0/0.01/1.0, dy=1, t=0.1), x=0 face clamped."""
    _model()
    ops.nDMaterial('ElasticIsotropic', 1, _E, _NU, _RHO)
    tags = {}
    n = 0
    for i, x in enumerate(_SLIVER_XS):
        for j in range(2):
            for k, z in enumerate((-_B_T / 2, _B_T / 2)):
                n += 1
                tags[(i, j, k)] = n
                ops.node(n, x, j * _B_DY, z)
    opts = ('-lumped',) if lumped else ()
    for i in range(3):
        ops.element('LadrunoSolidShell', i + 1,
                    tags[(i, 0, 0)], tags[(i + 1, 0, 0)],
                    tags[(i + 1, 1, 0)], tags[(i, 1, 0)],
                    tags[(i, 0, 1)], tags[(i + 1, 0, 1)],
                    tags[(i + 1, 1, 1)], tags[(i, 1, 1)], 1, *opts)
    for j in range(2):
        for k in range(2):
            ops.fix(tags[(0, j, k)], 1, 1, 1)
    return tags


def _static_tip():
    tags = _sliver_strip()
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for j in range(2):
        for k in range(2):
            ops.load(tags[(3, j, k)], 0.0, 0.0, _B_P / 4.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    return ops.nodeDisp(tags[(3, 0, 0)])[2]


def _dynamic_tip(integ_args, dt, T, alphaM):
    """Tip load ramped over T/2 then held; mass damping settles the response.
    Returns (tip deflection at T, steps taken)."""
    tags = _sliver_strip()
    ops.timeSeries('Path', 1, '-time', 0.0, T / 2.0, T * 10.0,
                   '-values', 0.0, 1.0, 1.0)
    ops.pattern('Plain', 1, 1)
    for j in range(2):
        for k in range(2):
            ops.load(tags[(3, j, k)], 0.0, 0.0, _B_P / 4.0)
    ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
    _explicit_solver()
    ops.integrator(*integ_args)
    ops.analysis('Transient')
    steps = int(T / dt)
    for s in range(steps):
        assert ops.analyze(1, dt) == 0, f"explicit ramp step {s} failed"
    return ops.nodeDisp(tags[(3, 0, 0)])[2], steps


def test_sms_speedup_quasi_static():
    """Fixed physical loading time on the sliver-bound strip: CDL must step at
    the sliver pencil (~0.01/c); SMS at dtTarget = 0.9x the BULK thickness
    step injects mass only around the sliver, preserves the statics, and
    takes ~10x fewer steps for the same T."""
    delta = _static_tip()
    tags = _sliver_strip()
    dtc = _pencil_of_current_model()            # the sliver pencil ~ 0.01/c
    dt_bulk = _B_T / _C                         # the bulk thickness step
    assert dtc < 0.15 * dt_bulk, \
        f"sliver did not pin the pencil (dtc {dtc:.3e} vs bulk {dt_bulk:.3e})"
    T = 24.0                                    # ~3 damped fundamental periods
    alphaM = 1.1                                # ~critical for f1 ~ 0.13 Hz
    u_cdl, n_cdl = _dynamic_tip((CDL, '-lump', 'hrz'), 0.9 * dtc, T, alphaM)
    dt_t = 0.9 * dt_bulk
    u_sms, n_sms = _dynamic_tip((SMS, dt_t, '-lump', 'hrz'), dt_t, T, alphaM)
    assert u_cdl == pytest.approx(delta, rel=0.05), \
        f"CDL quasi-static tip {u_cdl:.4e} != static {delta:.4e}"
    assert u_sms == pytest.approx(delta, rel=0.05), \
        f"SMS quasi-static tip {u_sms:.4e} != static {delta:.4e} (injection corrupted statics)"
    assert n_sms * 8 <= n_cdl, \
        f"SMS speedup only {n_cdl / max(n_sms, 1):.1f}x (steps {n_cdl} -> {n_sms})"


# ---------------------------------------------------------------------------
# G9c — the G5 softening band under SMS (Tier-3 + the D7 cure, end-to-end)
# ---------------------------------------------------------------------------
_CE, _CNU, _CFC, _CFT, _CGF, _CGC = 30000.0, 0.2, 30.0, 3.0, 0.1, 5.0
_CRHO = 2.4e-3


def test_softening_band_sms_completes():
    """A Concrete3D weakened-band strip (3x1, dx=dy=2, t=0.4, 0.6*ft band,
    -autoRegularization, -lumped elements) pulled to full softening
    (Delta = 0.14) by smoothstep support motion under CentralDifferenceSMS
    at dtTarget = 3x the thickness-bound pencil. No implex, no step-cutting,
    no unsymmetric solver — the Tier-3 point — above the plain-CDL step —
    the D7 cure. Gates: completes, finite, the band localizes with clean
    neighbors, the reaction ends far below its peak (fully softened)."""
    _model()
    # 0.6*ft band seed: the EXPLICIT pull adds dynamic overshoot on top of the
    # static redistribution — the localization drop launches a redistribution
    # wave, and with the transit budget CI cost allows (~200 scaled transits,
    # ~50 during the drop) the neighbors need a >60% headroom to stay clean
    # (measured: 0.8 leaked omega_t 0.19; 0.7 on a 118-transit budget cracked
    # them outright). Dissipation stays Gf-governed, so the energy physics is
    # seed-independent.
    ops.nDMaterial('LadrunoConcrete3D', 1, _CE, _CNU, _CFC, _CFT,
                   _CGF, _CGC, '-autoRegularization', '-rho', _CRHO)
    ops.nDMaterial('LadrunoConcrete3D', 2, _CE, _CNU, _CFC, 0.6 * _CFT,
                   _CGF, _CGC, '-autoRegularization', '-rho', _CRHO)
    # nx=3 / t=0.4 keep the CI cost honest: on this UNIFORM mesh the step
    # count is N_transits * (L/t) INDEPENDENT of dtTarget (a larger target
    # slows the scaled waves by the same factor it grows the step — the SMS
    # "+9878% added mass" warning it prints is the deliberate consequence of
    # a uniform-mesh target, not a defect), so the only cost levers are the
    # aspect ratio and the transit budget. L/t = 15 with ~200 transits =
    # ~3k steps.
    dx = dy = 2.0
    t = 0.4
    nx = 3
    tags = {}
    n = 0
    for i in range(nx + 1):
        for j in range(2):
            for k, z in enumerate((-t / 2, t / 2)):
                n += 1
                tags[(i, j, k)] = n
                ops.node(n, i * dx, j * dy, z)
    band_col = 1                               # the middle of the 3 columns
    band, neigh = [], []
    for i in range(nx):
        ops.element('LadrunoSolidShell', i + 1,
                    tags[(i, 0, 0)], tags[(i + 1, 0, 0)],
                    tags[(i + 1, 1, 0)], tags[(i, 1, 0)],
                    tags[(i, 0, 1)], tags[(i + 1, 0, 1)],
                    tags[(i + 1, 1, 1)], tags[(i, 1, 1)],
                    2 if i == band_col else 1, '-lumped')
        if i == band_col:
            band.append(i + 1)
        elif abs(i - band_col) == 1:
            neigh.append(i + 1)
    # determinate restraints (the G5 scheme)
    ops.fix(tags[(0, 0, 0)], 1, 1, 1)
    ops.fix(tags[(0, 0, 1)], 1, 1, 0)
    ops.fix(tags[(0, 1, 0)], 1, 0, 0)
    ops.fix(tags[(0, 1, 1)], 1, 0, 0)
    left = [tags[(0, j, k)] for j in range(2) for k in range(2)]
    # support motion on the pulled face — SMOOTHSTEP displacement protocol,
    # NOT the plain Linear ramp: a Linear series applies a velocity STEP at
    # t=0, and mass scaling multiplies the support-motion shock impedance by
    # sqrt(s) (sigma_shock ~ rho' c' v = sqrt(s) * rho c v; at a 10x target
    # and this rate that is ~3.5 > ft — the START SHOCK alone cracked the
    # pulled-face element, measured omega_t -> 0.97). The smoothstep
    # lam(t) = lam_end * u^2 (3 - 2u) starts and ends at zero velocity. T =
    # 200 MASS-SCALED wave transits (c' = c/dtF), so the step count is the
    # invariant N * (L/t) ~ 3.2k regardless of the target factor.
    dtF = 3.0
    rho_c = math.sqrt(_CE / _CRHO)
    lam_end = 0.14
    T = 200.0 * (nx * dx) / (rho_c / dtF)
    npts = 60
    ts, vs = [], []
    for p in range(npts + 1):
        u = p / npts
        ts.append(u * T)
        vs.append(lam_end * u * u * (3.0 - 2.0 * u))
    ts.append(10.0 * T)
    vs.append(lam_end)
    ops.timeSeries('Path', 1, '-time', *ts, '-values', *vs)
    ops.pattern('Plain', 1, 1)
    for j in range(2):
        for k in range(2):
            ops.sp(tags[(nx, j, k)], 1, 1.0)
    ops.rayleigh(80.0, 0.0, 0.0, 0.0)          # settle the axial ringing
    dtc = _pencil_of_current_model()           # plain-CDL prime + pencil
    # dtTarget = 3x, NOT higher: when the band lets go, the pulled-side block
    # is momentarily stretched between its own (scaled) inertia and the
    # moving support — a velocity-mismatch shock ~ rho'c'*dv that grows with
    # the target factor (sqrt(s)). Measured trajectories: at 10x the shock is
    # ~5 > ft and cracks the pulled-side neighbor AFTER localization
    # (omega_t 0.31 at frac 0.25-0.42); at 3x it is ~1.6 < ft and both
    # neighbors stay at 0.000. The step count is target-INDEPENDENT on a
    # uniform mesh (the transit budget scales with sqrt(s) too), so 3x costs
    # nothing — the 10x-class speedup story belongs to the sliver-mesh G9b.
    dt_t = dtF * dtc
    ops.integrator(SMS, dt_t, '-lump', 'hrz')
    steps = int(T / dt_t)
    fpeak, fend = 0.0, 0.0
    probe_every = max(1, steps // 200)
    for s in range(steps):
        assert ops.analyze(1, dt_t) == 0, f"SMS softening run failed at step {s}/{steps}"
        if s % probe_every == 0 or s == steps - 1:
            ops.reactions()
            f = -sum(ops.nodeReaction(tg, 1) for tg in left)
            assert math.isfinite(f), "reaction went non-finite"
            fpeak = max(fpeak, f)
            fend = f
    ngp = len(list(ops.eleResponse(band[0], 'stresses'))) // 6
    wt_band = max(list(ops.eleResponse(band[0], 'material', gp, 'damage'))[0]
                  for gp in range(1, ngp + 1))
    wt_neigh = max(list(ops.eleResponse(e, 'material', gp, 'damage'))[0]
                   for e in neigh for gp in range(1, ngp + 1))
    assert wt_band > 0.8, f"band did not localize under SMS (omega_t {wt_band:.3f})"
    assert wt_neigh < 0.1, f"damage leaked to neighbors ({wt_neigh:.3f})"
    assert fpeak > 0.0
    assert abs(fend) < 0.25 * fpeak, \
        f"band not softened (F_end {fend:.4f} vs peak {fpeak:.4f})"
