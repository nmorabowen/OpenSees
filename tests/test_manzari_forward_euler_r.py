"""WP-158 -- ``ManzariDafalias::ForwardEuler`` (IntScheme 5) had a shadowed ``r``.

THE DEFECT (vanilla; upstream OpenSees master has the same lines). ``ForwardEuler`` read

    Vector r(6);
    if (p > small)
        Vector r = GetDevPart(CurStress) / p;     // a NEW local, dies at the `;`

so the outer ``r`` stayed identically zero and every ``(n:r)`` term vanished
from the plastic multiplier: the ``-K D (n:r)`` of its denominator ``temp4``
and the ``-K de_v (n:r)`` of its numerator. Two tangent-only defects sat in
the same function. ``temp2 = 2G n - (n:r) I`` is missing the bulk modulus
(the numerator it differentiates is ``2G n:de - K de_v (n:r)``), invisible
while ``r == 0``. ``temp1 = 2G mIIdevMix + K mIIvol`` puts 2G, not G, on the
shear diagonal (the mixed-variant identity; the stress update answers an
engineering shear strain with G -- the WP-110 family). WP-158 fixes all three.

REACH. ``ForwardEuler`` is the substep integrator of IntScheme 5, of 4
(``MaxEnergyInc``, ``INT_MAXENE_FE``), of 7, 8 and 9 (``MaxStrainInc``'s
switch sends EVERY case to it), and of scheme 2's opt-in WP-130 guess walk.

WHAT THE BUG DID, MEASURED. Not a model error that survives refinement. The
wrong multiplier pushes the stress off the yield surface; when it lands
inside, ``explicit_integrator`` re-intersects it on the next step, which
enforces consistency geometrically, so in the fine-step limit the old scheme
5 converged onto scheme 1 too. What it broke is the ORDER of the step (gate
1) and, through it, the answer at practical step sizes (gate 2).

GATES (all from IntScheme 1 materials switched to 5 by ``setParameter
IntegrationScheme``: the scheme-3/5 warning latch is a process-wide static in
the CONSTRUCTOR, and ``test_manzari_safety_pack.py`` -- collected after this
file -- must be the first to construct scheme 5 to observe it).

  1. ``test_scheme5_step_satisfies_consistency_to_second_order`` -- THE
     LOAD-BEARING GATE. From a committed plastic state ON the yield surface
     (history under scheme 1), one scheme-5 step of size h, h/4: a multiplier
     that solves the linearised consistency condition leaves f unchanged to
     O(h^2), a wrong one to O(h). Drift ratio for 4x shorter step: 16 vs 4.
  2. ``test_scheme5_agrees_with_scheme1_drained_triaxial`` -- the requested
     comparison. Drained triaxial compression (cell pressure 100 kPa held on
     free lateral faces, 2 % axial strain under DisplacementControl) on one
     SSPbrick; IntScheme 1 at 1600 steps is the reference (converged: 6400
     steps move q/p by 7e-5). At 800 steps scheme 5 must sit within
     ``GAP_TOL`` of it in q/p and in eps_v, and the gap must halve by 3200
     steps (first order). The 800-step bound is the discriminating one; the
     old scheme 5 passes the halving test by accident (its error is erratic).
  3. ``test_scheme5_tangent_is_derivative_of_its_stress_update`` -- the
     ``temp1`` / ``temp2`` part. From a committed plastic state a forward-Euler
     step is LINEAR in the strain increment, so its TanType-1 tangent must
     equal a one-sided finite difference of its own stress, to round-off.

MEASURED 2026-10-01, Windows, Ladruno_scripts\\build.bat on origin/ladruno
d63f49750 unmodified ("pre") and with the WP-158 source ("fixed"), this
file's __main__:

  gate 1  Delta f/p at h, h/4, h/16   pre   -3.78e-3, -1.01e-3, -2.54e-4  (ratio 3.75, 3.97)
                                      fixed  3.55e-5,  2.24e-6,  1.40e-7  (ratio 15.9, 16.0)
  gate 2  gap (q/p, eps_v) at 800     pre    1.72e-2, 5.32e-2   fixed 3.35e-3, 1.77e-3
          ... at 3200                 pre    6.48e-4, 5.33e-4   fixed 9.11e-4, 4.57e-4
          (pre-fix at 50 steps: Newton diverges at eps_a 0.32 %, q/p 1.03 vs 0.75)
  gate 3  max |Ct - Cfd| / max |Cfd|  pre 0.56   r+K fixed, temp1 not 0.64   fixed 1.2e-11

Cost: < 2 s.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

# Toyoura reference set (the WP-110 tangent gate's MATP), medium-dense.
MATP = dict(G0=125.0, nu=0.05, e_init=0.75, Mc=1.25, c=0.712, lambda_c=0.019,
            e0=0.934, ksi=0.7, P_atm=101.3, m=0.01, h0=7.05, Ch=0.968, nb=1.1,
            A0=0.704, nd=3.5, z_max=4.0, cz=600.0, Rho=1.7)
_ORDER = ["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi", "P_atm",
          "m", "h0", "Ch", "nb", "A0", "nd", "z_max", "cz", "Rho"]
P_RESIDUAL = 1.0e-2 * MATP["P_atm"]      # vanilla's hardcoded m_Presidual

COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]

# gate 1
CONSIST_RATIO = 10.0     # second order gives 16, first order 4
# gate 2
P0 = 100.0               # cell pressure (kPa)
N_CONF = 10              # elastic confinement steps
E_AX = 0.02              # total axial compressive strain
N_REF = 1600             # IntScheme 1 reference
N_FE = 800               # IntScheme 5 at a practical step (2.5e-5 axial strain)
N_FE_FINE = 3200
GAP_TOL = 1.0e-2         # fixed 3.4e-3 / 1.8e-3; pre-fix 1.7e-2 / 5.3e-2
CHECK = (0.0025, 0.005, 0.01, 0.015, 0.02)   # axial-strain checkpoints
# gate 3
H_FD = 1e-7
TAN_RTOL = 1e-6          # fixed 1.2e-11; the 2G shear diagonal alone is ~0.5


def _mat(tag, tan_type=1):
    """Always constructed as IntScheme 1 (see the module docstring)."""
    p = MATP
    ops.nDMaterial("ManzariDafalias", tag, *[p[k] for k in _ORDER],
                   1, tan_type, 1, 1e-7, 1e-7)


def _set_scheme(scheme, ptag=1):
    """Runtime switch; argv[1] of ManzariDafalias::setParameter is the material tag."""
    ops.parameter(ptag, "element", 1, "IntegrationScheme", 1)
    ops.updateParameter(ptag, float(scheme))
    ops.remove("parameter", ptag)


def _brick():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in COORDS.items():
        ops.node(tag, float(x), float(y), float(z))
    _mat(1)
    ops.element("SSPbrick", 1, *range(1, 9), 1)


# ------------------------------------------- prescribed strain (gates 1, 3) ---
def _unit_E(j):
    a, b = VOIGT[j]
    E = np.zeros((3, 3))
    if a == b:
        E[a, b] = 1.0
    else:                       # engineering shear: gamma = 1
        E[a, b] = E[b, a] = 0.5
    return E


def _prescribed(strains, scheme, last_scheme=None):
    """Drive one SSPbrick through 3x3 strains (zero free DOFs, so the Gauss
    point sees exactly these): 5 elastic steps, then stage 1. Returns
    (stress, tangent, alpha) after the LAST step; `last_scheme` switches the
    integrator for the last step only."""
    _brick()
    ops.fix(1, 1, 1, 1)
    ops.constraints("Transformation")    # 'Plain' drops non-homogeneous SPs
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-9, 5, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    times = [float(t) for t in range(len(strains) + 1)]
    pts = [np.zeros((3, 3))] + list(strains)
    tag = 100
    for n in range(2, 9):
        X = np.array(COORDS[n], dtype=float)
        for d in range(3):
            tag += 1
            ops.timeSeries("Path", tag, "-time", *times,
                           "-values", *[float((E @ X)[d]) for E in pts])
            ops.pattern("Plain", tag, tag)
            ops.sp(n, d + 1, 1.0)
    if scheme != 1:
        _set_scheme(scheme)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for i in range(len(strains)):
        if i == 5:
            ops.updateMaterialStage("-material", 1, "-stage", 1)
        if last_scheme is not None and i == len(strains) - 1:
            _set_scheme(last_scheme)
        assert ops.analyze(1) == 0, "prescribed step failed"
    return (np.array(ops.eleResponse(1, "stress")),
            np.array(ops.eleResponse(1, "tangent")).reshape(6, 6),
            np.array(ops.eleResponse(1, "alpha")))


def _history():
    """Isotropic compression (elastic), then a deviatoric push with shear, so
    the end state is plastic and n carries shear components."""
    pts = [np.eye(3) * -1.5e-3 * (i + 1) / 5 for i in range(5)]
    dev = np.array([[0.5, 0.2, 0.0], [0.2, 0.5, 0.1], [0.0, 0.1, -1.0]]) * 4e-3
    for i in range(1, 41):
        pts.append(pts[4] + dev * i / 40)
    return pts, dev / 40


def yield_f(stress_ops, alpha):
    """(f, p), f = |s - p alpha| - sqrt(2/3) m p, in the material's convention
    (compression +; the 3D wrapper flips the sign of stress, not of alpha)."""
    sig = -np.asarray(stress_ops)
    p = (sig[0] + sig[1] + sig[2]) / 3.0 + P_RESIDUAL
    t = sig.copy()
    t[:3] -= p - P_RESIDUAL
    t -= p * np.asarray(alpha)
    w = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])     # tensor norm: Voigt shear twice
    return math.sqrt(float(np.sum(w * t * t))) - math.sqrt(2.0 / 3.0) * MATP["m"] * p, p


def consistency_drift(scheme=5, fracs=(1.0, 0.25)):
    """(f0/p, [Delta f / p]) over ONE `scheme` step of frac * (the history's
    step) from the end of _history(), whose path runs under IntScheme 1 so the
    start is ON the yield surface. (Scheme 5 has no drift correction and
    wanders ~10x the surface's width outside it, where its (n:r) is no longer
    the gradient term: the probe would measure that offset, not the fix.)"""
    pts, step = _history()
    s0, _, a0 = _prescribed(pts[:-1], 1)
    f0, p0 = yield_f(s0, a0)
    out = []
    for fr in fracs:
        s1, _, a1 = _prescribed(pts[:-1] + [pts[-2] + fr * step], 1, last_scheme=scheme)
        out.append((yield_f(s1, a1)[0] - f0) / p0)
    return f0 / p0, out


def tangent_vs_fd(scheme=5):
    pts, _ = _history()
    s0, Ct, _ = _prescribed(pts, scheme)
    Cfd = np.zeros((6, 6))
    for j in range(6):
        sj, _, _ = _prescribed(pts[:-1] + [pts[-1] + H_FD * _unit_E(j)], scheme)
        Cfd[:, j] = (sj - s0) / H_FD
    return Ct, Cfd


# ------------------------------------------------- drained triaxial (gate 2) ---
def drained_triaxial(scheme, n_shear, strict=True):
    """{eps_a: (q/p, eps_v)} at CHECK, compression +. Lateral faces x=1, y=1
    are free and carry the cell pressure as nodal forces (unit cube: P0/4 per
    node); the top face is tied in z to node 5 and driven by
    DisplacementControl once the confining load is constant."""
    _brick()
    for tag, (x, y, z) in COORDS.items():
        ops.fix(tag, 1 if x == 0 else 0, 1 if y == 0 else 0, 1 if z == 0 else 0)
    for n in (6, 7, 8):
        ops.equalDOF(5, n, 3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    f = -P0 / 4.0
    for tag, (x, y, z) in COORDS.items():
        ops.load(tag, f if x == 1 else 0.0, f if y == 1 else 0.0, f if z == 1 else 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    # ModifiedEuler's substep count changes between Newton iterates, so its
    # residual floors at ~6e-6 q; 1e-2 kPa is ~5e-5 of q, far below GAP_TOL.
    ops.test("NormUnbalance", 1e-2, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / N_CONF)
    ops.analysis("Static")
    if scheme != 1:
        _set_scheme(scheme)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(N_CONF):
        assert ops.analyze(1) == 0, "confinement step failed"

    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(5, 0.0, 0.0, -1.0)
    ops.integrator("DisplacementControl", 5, 3, -E_AX / n_shear)
    ops.analysis("Static")

    uz0 = ops.nodeDisp(5, 3)
    out, want = {}, list(CHECK)
    for i in range(1, n_shear + 1):
        rc = ops.analyze(1)
        if rc != 0 and not strict:          # __main__ measurement: keep what converged
            out["failed_at_eps_a"] = -(ops.nodeDisp(5, 3) - uz0)
            break
        assert rc == 0, f"scheme {scheme}: shear step {i}/{n_shear} failed"
        ea = -(ops.nodeDisp(5, 3) - uz0)
        if want and ea >= want[0] - 1e-12:
            s = -np.array(ops.eleResponse(1, "stress"))          # compression +
            p = (s[0] + s[1] + s[2]) / 3.0
            q = s[2] - 0.5 * (s[0] + s[1])
            ev = -(ops.nodeDisp(5, 3) + ops.nodeDisp(2, 1) + ops.nodeDisp(4, 2))
            assert abs(s[0] - P0) < 5e-4 * P0 and abs(s[1] - P0) < 5e-4 * P0, (
                "cell pressure not held", s[:3])                 # 4 nodes x the 1e-2 test
            out[want.pop(0)] = (q / p, ev)
    return out


def gaps(test, ref):
    """Max relative gap in q/p, and in eps_v over max|eps_v| of ref, at the
    checkpoints both runs reached."""
    keys = [k for k in CHECK if k in test and k in ref]
    ev_scale = max(abs(ref[k][1]) for k in keys)
    return (max(abs(test[k][0] - ref[k][0]) / abs(ref[k][0]) for k in keys),
            max(abs(test[k][1] - ref[k][1]) / ev_scale for k in keys))


# ------------------------------------------------------------------ gates ---
def test_scheme5_step_satisfies_consistency_to_second_order():
    f0, (d1, d4) = consistency_drift(5)
    assert abs(f0) < 1e-8, ("start state is not on the yield surface", f0)
    ratio = abs(d1) / max(abs(d4), 1e-300)
    assert ratio > CONSIST_RATIO, (
        "one IntScheme-5 step drifts off the yield surface at FIRST order "
        "(Delta f/p at h, h/4, ratio) -- the multiplier is not the linearised "
        "consistency solution", d1, d4, ratio)


def test_scheme5_agrees_with_scheme1_drained_triaxial():
    ref = drained_triaxial(1, N_REF)
    g = gaps(drained_triaxial(5, N_FE), ref)
    gf = gaps(drained_triaxial(5, N_FE_FINE), ref)
    assert g[0] < GAP_TOL and g[1] < GAP_TOL, (
        f"IntScheme 5 at {N_FE} steps is not within {GAP_TOL} of IntScheme 1 "
        "(q/p, eps_v)", g)
    assert gf[0] < 0.5 * g[0] and gf[1] < 0.5 * g[1], (
        "IntScheme 5 does not converge at first order onto IntScheme 1", g, gf)


def test_scheme5_tangent_is_derivative_of_its_stress_update():
    Ct, Cfd = tangent_vs_fd(5)
    scale = np.abs(Cfd).max()
    assert np.abs(Ct - Ct.T).max() > 1e-6 * scale, (
        "last step is elastic (Ce is symmetric) -- the gate would be vacuous")
    err = np.abs(Ct - Cfd).max() / scale
    assert err < TAN_RTOL, ("scheme-5 tangent != FD of its own stress update", err, Ct, Cfd)


if __name__ == "__main__":
    # Measurement mode: the numbers in the module docstring.
    f0, d = consistency_drift(5, (1.0, 0.25, 0.0625))
    print("gate 1: f0/p %.2e  Delta f/p at h, h/4, h/16: %s  ratios %.2f %.2f" % (
        f0, ", ".join("%.3e" % x for x in d), abs(d[0] / d[1]), abs(d[1] / d[2])))
    ref = drained_triaxial(1, N_REF)
    for n in (50, 200, 400, N_FE, 1600, N_FE_FINE):
        r = drained_triaxial(5, n, strict=False)
        tail = "  FAILED at eps_a %.4f" % r["failed_at_eps_a"] if "failed_at_eps_a" in r else ""
        print("gate 2: IntScheme 5 %5d steps  gap (q/p, eps_v) = (%.3e, %.3e)%s" % (
            (n,) + gaps(r, ref) + (tail,)))
    Ct, Cfd = tangent_vs_fd(5)
    print("gate 3: max |Ct - Cfd| / max |Cfd| = %.3e" % (np.abs(Ct - Cfd).max() / np.abs(Cfd).max()))
