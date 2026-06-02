"""
Damping-channel verification battery.

Confirms the *behavioural / numeric* claims from the OpenSees damping deep-dive
that source-reading alone does not settle (the mechanism is visible in the code,
but not the realised damping ratio). The pure code-level facts (e.g. that
`setRayleighDampingFactors` is plain assignment) are NOT re-tested here -- they
are transcription, not inference.

What IS tested (all via free-vibration log-decrement, implicit Newmark
average-acceleration so there is no algorithmic damping to contaminate zeta):

  1. Element stiffness-proportional Rayleigh (betaK0) realises  zeta = betaK0*w/2.
  2. Nodal mass-proportional Rayleigh (alphaM) realises          zeta = alphaM/(2w).
  3. ZeroLength QUIRK: a zeroLength element does NOT contribute stiffness-
     proportional Rayleigh damping unless built with `-doRayleigh 1`. The
     default (`-doRayleigh 0`) yields ZERO stiffness damping -- a classic
     OpenSees footgun. Both the default-off and the opted-in behaviour are
     pinned here.
  4. OVERWRITE: domain `rayleigh` then per-element `rayleigh -ele` does NOT sum;
     the element keeps only the last value (last-writer-wins).
  5. ADDITIVE: `rayleigh` + `modalDamping` stack (the two channels add).
  6. modalDamping realises the target zeta in the targeted mode (2-DOF).
  7. modalDamping SCOPE: a mode given factor 0 (or left unspecified) is undamped,
     even though it is a computed eigen-mode -- no spillover from the damped mode.
  8. modalDamping (matrix-in-tangent) realises the target zeta on a linear SDOF.
     modalDampingQ (force-only) is xfail-documented here: in a free-vibration
     SDOF it produced ANTI-damping of the same magnitude regardless of
     Newton/Linear -- flagged for a closer look (only the M-weighted
     addModalDampingForce at IncrementalIntegrator.cpp:502 is live; the two
     earlier variants are commented out).

Run:  <py3.12> -m pytest tests/test_damping_channels.py -v
(the dist/bin opensees.pyd is built for CPython 3.12)
"""
import os
import sys
import math

# --- bootstrap the locally-built engine -----------------------------------
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
if os.path.isdir(DIST):
    try:
        os.add_dll_directory(DIST)
    except (FileNotFoundError, OSError):
        pass
    if DIST not in sys.path:
        sys.path.insert(0, DIST)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

import numpy as np


# --- helpers ----------------------------------------------------------------
def _zeta_from_decay(t, u, n_skip=1, n_use=5):
    """Estimate the damping ratio of a free-decay record by log decrement.

    Uses positive local maxima (peaks). Skips the first `n_skip` peaks
    (initial transient) and measures over the next `n_use` cycles.
    """
    u = np.asarray(u)
    # positive local maxima
    peaks = []
    for i in range(1, len(u) - 1):
        if u[i] > u[i - 1] and u[i] >= u[i + 1] and u[i] > 0.0:
            peaks.append((t[i], u[i]))
    assert len(peaks) >= n_skip + n_use + 1, (
        f"not enough peaks ({len(peaks)}) to measure damping"
    )
    p0 = peaks[n_skip][1]
    pn = peaks[n_skip + n_use][1]
    delta = math.log(p0 / pn) / n_use          # log decrement per cycle
    return delta / math.sqrt(4.0 * math.pi ** 2 + delta ** 2)


def _run_sdof_free(m, k, dt, nsteps, d0=1.0,
                   rayleigh_domain=None, rayleigh_ele=None,
                   modal=None, modalQ=None, do_rayleigh=1):
    """Build a 1-DOF zeroLength oscillator (mass on the node, stiffness on the
    element) and integrate free vibration with Newmark (avg accel).

    rayleigh_domain / rayleigh_ele : (alphaM, betaK, betaK0, betaKc) or None
    modal / modalQ                 : zeta (float) or None  -> modalDamping(Q)
    do_rayleigh                    : zeroLength -doRayleigh flag (1 = include
                                     element stiffness in Rayleigh damping).
    Returns (t[], u[]).
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1,
                '-doRayleigh', do_rayleigh)

    if modal is not None or modalQ is not None:
        ops.constraints('Transformation'); ops.numberer('Plain')
        ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-12, 10)
        ops.algorithm('Linear'); ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')
        ops.eigen('-fullGenLapack', 1)
        if modal is not None:
            ops.modalDamping(modal)
        else:
            ops.modalDampingQ(modalQ)
    else:
        ops.constraints('Transformation'); ops.numberer('Plain')
        ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-12, 10)
        ops.algorithm('Linear'); ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')

    if rayleigh_domain is not None:
        ops.rayleigh(*rayleigh_domain)
    if rayleigh_ele is not None:
        ops.rayleigh(*rayleigh_ele, '-ele', 1)

    ops.setNodeDisp(2, 1, d0, '-commit')

    t, u = [], []
    for _ in range(nsteps):
        ops.analyze(1, dt)
        t.append(ops.getTime())
        u.append(ops.nodeDisp(2, 1))
    return np.array(t), np.array(u)


# Common SDOF: f = 1 Hz, w = 2*pi, T = 1 s
M = 1.0
W = 2.0 * math.pi
K = W * W * M
DT = 0.005          # 200 steps / period
NSTEPS = 1600       # 8 periods


# --- 1. element stiffness-proportional Rayleigh (needs -doRayleigh 1) --------
def test_betaK0_realises_target_zeta():
    zeta_target = 0.05
    betaK0 = 2.0 * zeta_target / W                       # zeta = betaK0*w/2
    t, u = _run_sdof_free(M, K, DT, NSTEPS, do_rayleigh=1,
                          rayleigh_domain=(0.0, 0.0, betaK0, 0.0))
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, zeta_target, rel_tol=0.10), \
        f"betaK0: measured zeta={zeta:.4f}, target={zeta_target:.4f}"


# --- 1b. ZeroLength quirk: default -doRayleigh 0 => ZERO stiffness damping ----
def test_zeroLength_doRayleigh_default_off():
    betaK0 = 2.0 * 0.05 / W
    t, u = _run_sdof_free(M, K, DT, NSTEPS, do_rayleigh=0,
                          rayleigh_domain=(0.0, 0.0, betaK0, 0.0))
    zeta = _zeta_from_decay(t, u)
    assert abs(zeta) < 0.005, (
        f"zeroLength with -doRayleigh 0 should give ~zero stiffness damping, "
        f"got zeta={zeta:.4f}")


# --- 2. nodal mass-proportional Rayleigh ------------------------------------
def test_alphaM_realises_target_zeta():
    zeta_target = 0.05
    alphaM = 2.0 * zeta_target * W                       # zeta = alphaM/(2w)
    t, u = _run_sdof_free(M, K, DT, NSTEPS,
                          rayleigh_domain=(alphaM, 0.0, 0.0, 0.0))
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, zeta_target, rel_tol=0.10), \
        f"alphaM: measured zeta={zeta:.4f}, target={zeta_target:.4f}"


# --- 3. OVERWRITE: per-element rayleigh replaces, does not sum ---------------
def test_element_rayleigh_overwrites_domain():
    z_dom, z_ele = 0.02, 0.06
    beta_dom = 2.0 * z_dom / W
    beta_ele = 2.0 * z_ele / W
    # domain sets betaK0 on every element, then element-1 is overwritten
    t, u = _run_sdof_free(
        M, K, DT, NSTEPS, do_rayleigh=1,
        rayleigh_domain=(0.0, 0.0, beta_dom, 0.0),
        rayleigh_ele=(0.0, 0.0, beta_ele, 0.0),
    )
    zeta = _zeta_from_decay(t, u)
    # If it OVERWRITES -> zeta ~ z_ele (0.06). If it SUMMED -> ~0.08.
    assert math.isclose(zeta, z_ele, rel_tol=0.12), \
        f"overwrite: measured zeta={zeta:.4f}, expected ~{z_ele:.4f} (NOT {z_dom+z_ele:.4f})"
    assert abs(zeta - (z_dom + z_ele)) > 0.012, \
        f"overwrite: zeta={zeta:.4f} looks like the SUM -> channels wrongly additive"


# --- 4. ADDITIVE: rayleigh + modalDamping stack -----------------------------
def test_rayleigh_plus_modal_is_additive():
    z_ray = 0.03
    z_mod = 0.04
    beta = 2.0 * z_ray / W
    t, u = _run_sdof_free(M, K, DT, NSTEPS, do_rayleigh=1,
                          rayleigh_domain=(0.0, 0.0, beta, 0.0),
                          modal=z_mod)
    zeta = _zeta_from_decay(t, u)
    assert math.isclose(zeta, z_ray + z_mod, rel_tol=0.12), \
        f"additive: measured zeta={zeta:.4f}, expected ~{z_ray+z_mod:.4f}"


# --- 5/6. modal damping: per-mode value + scope (zero on undamped mode) ------
def _run_2dof_modal_free(zetas, exc_mode, dt, nsteps, d_scale=1.0):
    """2-DOF shear chain. eigen 2; modalDamping with per-mode `zetas`
    (list of len 2; a 0.0 entry leaves that mode undamped). Excite `exc_mode`
    by seeding its mode shape as the initial displacement. Free vibration.
    Returns (t[], roof-disp[]).
    """
    m, k = 1.0, (2.0 * math.pi) ** 2
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 1.0); ops.mass(2, m)
    ops.node(3, 2.0); ops.mass(3, m)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    ops.element('zeroLength', 2, 2, 3, '-mat', 1, '-dir', 1)

    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-12, 10)
    ops.algorithm('Linear'); ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')

    ops.eigen('-fullGenLapack', 2)
    ops.modalDamping(*zetas)

    # seed the requested mode shape
    phi2 = ops.nodeEigenvector(2, exc_mode, 1)
    phi3 = ops.nodeEigenvector(3, exc_mode, 1)
    nrm = max(abs(phi2), abs(phi3))
    ops.setNodeDisp(2, 1, d_scale * phi2 / nrm, '-commit')
    ops.setNodeDisp(3, 1, d_scale * phi3 / nrm, '-commit')

    t, u = [], []
    for _ in range(nsteps):
        ops.analyze(1, dt)
        t.append(ops.getTime())
        u.append(ops.nodeDisp(3, 1))
    return np.array(t), np.array(u)


def test_modal_per_mode_value():
    # damp mode 1 at 0.05, mode 2 at 0.03; excite each, measure its own zeta
    t1, u1 = _run_2dof_modal_free([0.05, 0.03], exc_mode=1, dt=0.004, nsteps=2500)
    z1 = _zeta_from_decay(t1, u1, n_use=4)
    assert math.isclose(z1, 0.05, rel_tol=0.15), f"mode1 zeta={z1:.4f}, target 0.05"

    t2, u2 = _run_2dof_modal_free([0.05, 0.03], exc_mode=2, dt=0.002, nsteps=2500)
    z2 = _zeta_from_decay(t2, u2, n_use=4)
    assert math.isclose(z2, 0.03, rel_tol=0.20), f"mode2 zeta={z2:.4f}, target 0.03"


def test_modal_scope_zero_mode_is_undamped():
    # damp mode 1 only; mode 2 factor 0.0 -> mode 2 must stay (nearly) undamped
    t2, u2 = _run_2dof_modal_free([0.05, 0.0], exc_mode=2, dt=0.002, nsteps=2500)
    z2 = _zeta_from_decay(t2, u2, n_use=5)
    assert z2 < 0.008, f"mode2 should be ~undamped, measured zeta={z2:.4f}"


# --- 8a. modalDamping (matrix form) realises target zeta on a linear SDOF ----
def test_modalDamping_matrix_realises_target():
    z = 0.05
    t_m, u_m = _run_sdof_free(M, K, DT, NSTEPS, modal=z)
    zm = _zeta_from_decay(t_m, u_m)
    assert math.isclose(zm, z, rel_tol=0.12), f"modalDamping zeta={zm:.4f}"


# --- 8b. modalDampingQ (force-only) anomaly: DOCUMENTED, expected to fail -----
import pytest


@pytest.mark.xfail(reason="modalDampingQ (force-only) anti-damps in a free-vib "
                          "SDOF: zeta comes out ~ -target regardless of "
                          "Newton/Linear. Flagged for investigation.",
                   strict=True)
def test_modalDampingQ_force_only_matches_matrix():
    z = 0.05
    t_q, u_q = _run_sdof_free(M, K, DT, NSTEPS, modalQ=z)
    zq = _zeta_from_decay(t_q, u_q)
    # If modalDampingQ behaved like modalDamping this would hold:
    assert math.isclose(zq, z, rel_tol=0.12), f"modalDampingQ zeta={zq:.4f}"


if __name__ == "__main__":
    # allow running as a plain script too
    import traceback
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    npass = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
            npass += 1
        except Exception as e:
            print(f"  FAIL  {fn.__name__}: {e}")
            traceback.print_exc()
    print(f"\n{npass}/{len(fns)} passed")
