"""WP-110 / F15 -- gate for ManzariDafalias's continuum elastoplastic tangent.

WHAT IS GATED. ``ManzariDafalias::GetElastoPlasticTangent`` had two
independent Voigt-convention defects, measured by the F15a probe (this file's
predecessor, ``test_manzari_ep_tangent_probe.py``; results in
``Ladruno_implementation/_wp110_f15a_probe_results.md``):

  1. numerator: ``DoubleDot4_2(aC, ToCovariant(R))`` with ``R`` already
     covariant -- an exact 2x error in aCep's shear rows;
  2. denominator: ``DoubleDot2_2_Contr(temp2, R)`` contracting a contravariant
     ``Q:Ce`` against a covariant ``R`` -- the denominator ~7.6 % too large.

And ``ModifiedEuler`` (IntScheme 1) never wrote ``aCep`` at all, so TanType 1
there handed the element a stale matrix (usually ``Ce``). WP-110 fixes all
three; this file gates the ENGINE's own tangent -- read back through the new
``eleResponse(ele, 'tangent')`` material response -- against:

  * a finite-difference (FD) tangent of the engine's own stress update, which
    is what a tangent must match; tolerance 1.5 % (measured worst 0.496 %) on every entry that carries a
    real signal;
  * the numpy "workbench" formula (full 3x3 tensor algebra, no Voigt
    bookkeeping at all: ``Dep = Ce - (Ce:R) (x) (Q:Ce) / (Kp + Q:Ce:R)``) --
    an independent transcription, tolerance 1e-3.

THE FD MUST BE ONE-SIDED AND REVERSAL-CONSISTENT. ``integrate()`` resets
``alpha_in := alpha_n`` when ``(alpha_n - alpha_in):Ce:d_eps < 0``. The sign
of that test depends on the probe direction, so a central difference
``(F(+h) - F(-h)) / 2h`` differences two DIFFERENT internal states and
straddles a real kink (measured: off-diagonal FD columns 40-100 % wrong). The
harness therefore probes each direction with whichever sign leaves ``alpha_in``
equal to the committed base state, and takes a one-sided difference.

THE STATE MUST CARRY SHEAR. A coaxial (no-shear) triaxial path makes ``R``'s
shear entries zero and hides BOTH defects identically. The two states here
(p' ~ 20 and ~ 214 kPa, eta ~ 1.3) have ``n_xy ~ 0.21``, asserted.

THE LAST STEP IS TINY ON PURPOSE. Both integrators evaluate G and K at the
START of an increment, while the FD samples the rate response AT the committed
state. A final continuation segment of 1e-3 of the main loading segment makes
the last increment start (to ~1e-5) at the state the FD probes, so the two are
comparing the same point, not two points one load step apart.

HARNESS. One 8-node ``SSPbrick`` (single shared integration point) with the 3D
``ManzariDafalias``, every non-pinned DOF driven by an SP under a ``Path``
series (``constraints Transformation`` -- ``Plain`` drops non-homogeneous SPs),
so the material sees an exact affine strain history with zero free DOFs.
``updateMaterialStage 1`` is REQUIRED: the full constructor starts the
material elastic (``mElastFlag = 0``) and without the flip every tangent is Ce.
The 3D wrapper flips the sign of stress AND strain, so the tangent is
invariant and no sign correction is needed.

COST: 2 states x 2 schemes x (1 + <= 12) single-element runs plus the sanity
and TanType-1 legs -- about 60 builds of ~120 steps each. Budget < 3 min;
measured 44-60 s on the dev box (build dee04dbe3), so not marked slow.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

# ---------------------------------------------------------------------------
# Material parameters: Toyoura-sand reference set (Ghofrani/Arduino UW example).
# ---------------------------------------------------------------------------
MATP = dict(G0=125.0, nu=0.05, e_init=0.7, Mc=1.25, c=0.712, lambda_c=0.019,
            e0=0.934, ksi=0.7, P_atm=101.3, m=0.01, h0=7.05, Ch=0.968, nb=1.1,
            A0=0.704, nd=3.5, z_max=4.0, cz=600.0, Rho=1.7)
M_PMIN = 1.0e-4 * MATP['P_atm']
M_PRESIDUAL = 1.0e-2 * MATP['P_atm']
SMALL = 1e-10
ONE3 = 1.0 / 3.0
TWO3 = 2.0 / 3.0
ROOT23 = math.sqrt(2.0 / 3.0)
I1 = np.array([1., 1., 1., 0., 0., 0.])

COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
FREE_DOFS = [(n, d) for n in range(2, 9) for d in (1, 2, 3)]
VOIGT_COMPS = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]

H_FD = 1.0e-8          # the one FD step size (probe measured 0.02-0.17 % at it)
FD_RTOL = 1.5e-2       # engine tangent vs one-sided FD. Measured worst 0.496 %
                       # (p'~20, IntScheme 2) -- FD truncation plus BE's algorithmic
                       # step, not tangent error. 1.5 % keeps cross-platform margin and
                       # still catches either defect alone (7.6 % denominator, 2x shear);
                       # ORACLE_RTOL below carries the precision.
ORACLE_RTOL = 1.0e-3   # engine tangent vs the numpy workbench formula
N_MAIN = 40            # steps per main loading segment
N_TAIL = 4             # steps in the tiny continuation segment
TAIL_FRAC = 1.0e-3     # continuation length, as a fraction of the last segment

# (label, iso, yy, xy): isotropic compression, then a yy + xy shear push.
STATES = [
    ("low_p20", -1.5e-5, -2.7e-3, -5e-4),
    ("high_p214", -4.8e-5, -8.64e-3, -1.6e-3),
]


# --------------------------------------------------- numpy tensor primitives -
def trace(v):
    return v[0] + v[1] + v[2]


def devpart(v):
    p = trace(v)
    r = v.copy()
    r[0] -= ONE3 * p; r[1] -= ONE3 * p; r[2] -= ONE3 * p
    return r


def single_dot(v1, v2):
    r = np.zeros(6)
    r[0] = v1[0]*v2[0] + v1[3]*v2[3] + v1[5]*v2[5]
    r[1] = v1[3]*v2[3] + v1[1]*v2[1] + v1[4]*v2[4]
    r[2] = v1[5]*v2[5] + v1[4]*v2[4] + v1[2]*v2[2]
    r[3] = 0.5*(v1[0]*v2[3]+v1[3]*v2[0]+v1[3]*v2[1]+v1[1]*v2[3]+v1[5]*v2[4]+v1[4]*v2[5])
    r[4] = 0.5*(v1[3]*v2[5]+v1[5]*v2[3]+v1[1]*v2[4]+v1[4]*v2[1]+v1[4]*v2[2]+v1[2]*v2[4])
    r[5] = 0.5*(v1[0]*v2[5]+v1[5]*v2[0]+v1[3]*v2[4]+v1[4]*v2[3]+v1[5]*v2[2]+v1[2]*v2[5])
    return r


def dd22_contr(v1, v2):
    s = 0.0
    for i in range(6):
        s += v1[i]*v2[i] + (1.0 if i > 2 else 0.0)*v1[i]*v2[i]
    return s


def norm_contr(v):
    return math.sqrt(dd22_contr(v, v))


def macauley(x):
    return x if x > 0 else 0.0


def macauley_index(x):
    return 1.0 if x > 0 else 0.0


def g_fun(cos3theta, c):
    return 2*c/((1+c)-(1-c)*cos3theta)


def get_stiffness(K, G):
    a = K + 4.0*ONE3*G; b = K - 2.0*ONE3*G
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = a
    C[3, 3] = C[4, 4] = C[5, 5] = G
    C[0, 1] = C[0, 2] = C[1, 2] = b
    C[1, 0] = C[2, 0] = C[2, 1] = b
    return C


def get_normal_to_yield(stress, alpha):
    p = ONE3*trace(stress) + M_PRESIDUAL
    if abs(p) < SMALL:
        return np.zeros(6)
    n = -p*alpha + devpart(stress)
    nn = norm_contr(n)
    if nn < SMALL:
        nn = 1.0
    return n/nn


def get_elastic_moduli(sigma):
    # stage 1 keeps the sqrt(pn/Patm) pressure factor (ManzariDafalias.cpp GetElasticModuli)
    pn = ONE3*trace(sigma)
    pn = M_PMIN if pn <= M_PMIN else pn
    eG = MATP['e_init']
    G = MATP['G0']*MATP['P_atm']*(2.97-eG)**2/(1+eG)*math.sqrt(pn/MATP['P_atm'])
    K = TWO3*(1+MATP['nu'])/(1-2*MATP['nu'])*G
    return K, G


def get_state_dependent(stress, alpha, fabric, e, alpha_in):
    p = ONE3*trace(stress) + M_PRESIDUAL
    p = SMALL if p < SMALL else p
    n = get_normal_to_yield(stress, alpha)
    alpha_alpha_in_dot_n = dd22_contr(alpha - alpha_in, n)
    psi = e - (MATP['e0'] - MATP['lambda_c']*(p/MATP['P_atm'])**MATP['ksi'])
    cos3theta = math.sqrt(6.0)*trace(single_dot(n, single_dot(n, n)))
    cos3theta = min(1.0, max(-1.0, cos3theta))
    alpha_b = g_fun(cos3theta, MATP['c'])*MATP['Mc']*math.exp(-MATP['nb']*psi) - MATP['m']
    alpha_d = g_fun(cos3theta, MATP['c'])*MATP['Mc']*math.exp(MATP['nd']*psi) - MATP['m']
    b0 = MATP['G0']*MATP['h0']*(1.0-MATP['Ch']*e)/math.sqrt(p/MATP['P_atm'])
    d = ROOT23*alpha_d*n - alpha
    b = ROOT23*alpha_b*n - alpha
    h = 1.0e10 if abs(alpha_alpha_in_dot_n) < SMALL else b0/alpha_alpha_in_dot_n
    A = MATP['A0']*(1+macauley(dd22_contr(fabric, n)))
    D = A*dd22_contr(d, n)
    if p < 0.05*MATP['P_atm']:
        D *= 1.0/(1.0+math.exp(7.6349 - 7.2713*101.0/MATP['P_atm']*p))
    B = 1.0+1.5*(1-MATP['c'])/MATP['c']*g_fun(cos3theta, MATP['c'])*cos3theta
    C = 3.0*math.sqrt(1.5)*(1-MATP['c'])/MATP['c']*g_fun(cos3theta, MATP['c'])
    return dict(n=n, d=d, b=b, h=h, B=B, C=C, D=D)


def _vec_to_ten(v):
    T = np.zeros((3, 3))
    for a, (i, j) in enumerate(VOIGT_COMPS):
        T[i, j] = T[j, i] = v[a]
    return T


def _ten_to_vec(T):
    return np.array([T[i, j] for (i, j) in VOIGT_COMPS])


def _dev3(T):
    return T - np.trace(T)/3.0*np.eye(3)


def _fro(A, B):
    return float(np.tensordot(A, B, axes=2))


def workbench_tangent(sigma, dgamma, K, G, sd):
    """Dep = Ce - (Ce:R) (x) (Q:Ce) / (Kp + Q:Ce:R), in full 3x3 tensor algebra
    (the Workbench's sanisand_cep.py ``Dep6``). Returned in the engine's Voigt
    convention: rows are stress (tensor) components, columns multiply
    engineering (covariant) strain -- hence ``Q:Ce`` enters as tensor
    components, which dotted with engineering strain is the true contraction."""
    S = _vec_to_ten(sigma)
    p = np.trace(S)/3.0 + M_PRESIDUAL
    p = (SMALL+M_PRESIDUAL) if p < (SMALL+M_PRESIDUAL) else p
    r_t = _dev3(S)/p
    n_t = _vec_to_ten(sd['n'])
    Kp = TWO3*p*sd['h']*_fro(_vec_to_ten(sd['b']), n_t)
    R_t = sd['B']*n_t - sd['C']*(n_t@n_t - np.eye(3)/3.0) + sd['D']/3.0*np.eye(3)
    Q_t = n_t - _fro(n_t, r_t)/3.0*np.eye(3)
    CeR = 2.0*G*_dev3(R_t) + K*np.trace(R_t)*np.eye(3)
    QCe = 2.0*G*_dev3(Q_t) + K*np.trace(Q_t)*np.eye(3)
    denom = Kp + _fro(QCe, R_t)
    aC = get_stiffness(K, G)
    if abs(denom) < SMALL:
        return aC
    return aC - macauley_index(dgamma)*np.outer(_ten_to_vec(CeR), _ten_to_vec(QCe))/denom


# --------------------------------------------------------------- FE harness --
_tagctr = [9000]


def build(int_scheme, tan_type, mat_tag=1):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for tag, (x, y, z) in COORDS.items():
        ops.node(tag, float(x), float(y), float(z))
    ops.fix(1, 1, 1, 1)
    p = MATP
    ops.nDMaterial('ManzariDafalias', mat_tag, p['G0'], p['nu'], p['e_init'], p['Mc'], p['c'],
                   p['lambda_c'], p['e0'], p['ksi'], p['P_atm'], p['m'], p['h0'], p['Ch'],
                   p['nb'], p['A0'], p['nd'], p['z_max'], p['cz'], p['Rho'],
                   int_scheme, tan_type, 1, 1e-7, 1e-7)
    ops.element('SSPbrick', 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)
    # the full constructor starts ELASTIC (mElastFlag = 0); without this every
    # tangent below would be Ce and the gate would be vacuous
    ops.updateMaterialStage('-material', mat_tag, '-stage', 1)
    ops.constraints('Transformation')   # 'Plain' silently drops non-homogeneous SPs
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-9, 30)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    _tagctr[0] = 9000


def strain_to_disps(E):
    d = {}
    for n in range(2, 9):
        u = E @ np.array(COORDS[n], dtype=float)
        d[(n, 1)] = u[0]; d[(n, 2)] = u[1]; d[(n, 3)] = u[2]
    return d


def unit_dE(idx):
    a, b = VOIGT_COMPS[idx]
    dE = np.zeros((3, 3))
    if a == b:
        dE[a, b] = 1.0
    else:   # engineering shear strain gamma = 1 -> eps_ab = eps_ba = 0.5
        dE[a, b] = 0.5; dE[b, a] = 0.5
    return dE


def path_for(state):
    """Waypoints and per-segment step counts: isotropic, shear push, tiny tail."""
    _, iso, yy, xy = state
    e1 = np.diag([iso, iso, iso])
    e2 = e1.copy()
    e2[1, 1] += yy; e2[0, 1] += xy; e2[1, 0] += xy
    e3 = e2 + TAIL_FRAC*(e2 - e1)
    return [e1, e2, e3], [N_MAIN, N_MAIN, N_TAIL]


def checkpoints(waypoints, nsteps, extra=None):
    pts = [np.zeros((3, 3))]
    prev = pts[0]
    for w, n in zip(waypoints, nsteps):
        for i in range(1, n+1):
            pts.append(prev + (w-prev)*i/n)
        prev = w
    if extra is not None:
        pts.append(prev + extra)
    vals = {k: [] for k in FREE_DOFS}
    for E in pts:
        dd = strain_to_disps(E)
        for k in FREE_DOFS:
            vals[k].append(dd[k])
    return vals


def run(vals):
    n_points = len(next(iter(vals.values())))
    times = [float(t) for t in range(n_points)]
    for (node, dof), series in vals.items():
        _tagctr[0] += 1
        tag = _tagctr[0]
        ops.timeSeries('Path', tag, '-time', *times, '-values', *[float(v) for v in series])
        ops.pattern('Plain', tag, tag)
        ops.sp(node, dof, 1.0)
    for _ in range(n_points - 1):
        assert ops.analyze(1) == 0, "prescribed-strain step failed"


def wrapper_stress():
    return np.array(ops.eleResponse(1, 'stress'))   # = -mSigma (3D sign flip)


def engine_tangent():
    t = ops.eleResponse(1, 'tangent')
    assert t is not None and len(t) == 36, (
        "eleResponse(1, 'tangent') returned %r -- the WP-110 'tangent' material "
        "response is missing, i.e. this is a STALE opensees binary" % (t,))
    return np.array(t).reshape(6, 6)   # Information::getData flattens row-major


def internal_state():
    state = np.array(ops.eleResponse(1, 'state'))
    return dict(sigma=-1.0*wrapper_stress(),
                alpha=np.array(ops.eleResponse(1, 'alpha')),
                fabric=np.array(ops.eleResponse(1, 'fabric')),
                alpha_in=np.array(ops.eleResponse(1, 'alpha_in')),
                e=state[24], dgamma=state[25])


def base_run(state, int_scheme, tan_type):
    wp, ns = path_for(state)
    build(int_scheme, tan_type)
    run(checkpoints(wp, ns))
    return engine_tangent(), wrapper_stress(), internal_state()


def fd_one_sided(state, int_scheme, tan_type, F0, base_alpha_in):
    """One-sided FD, each direction probed with the sign that does NOT trip the
    alpha_in reversal reset. Returns (C_fd, signs)."""
    wp, ns = path_for(state)
    C = np.zeros((6, 6))
    signs = []
    for j in range(6):
        for sign in (+1, -1):
            build(int_scheme, tan_type)
            run(checkpoints(wp, ns, extra=sign*H_FD*unit_dE(j)))
            if np.allclose(internal_state()['alpha_in'], base_alpha_in, rtol=0.0, atol=1e-9):
                C[:, j] = (wrapper_stress() - F0) / (sign*H_FD)
                signs.append(sign)
                break
        else:
            pytest.fail("direction %d: both probe signs reset alpha_in -- no "
                        "reversal-consistent one-sided FD exists at this state" % j)
    return C, signs


def rel_err_masked(A, ref):
    """Max relative error over entries of `ref` carrying a real signal (> 1 % of
    its largest entry). Entries below that are ~0 analytically and FD round-off
    (a few units at h = 1e-8) would divide noise by noise."""
    mask = np.abs(ref) > 1e-2*np.max(np.abs(ref))
    return float(np.max(np.abs(A[mask]-ref[mask])/np.abs(ref[mask]))), mask


def _build_hash():
    try:
        return ops.ladrunoBuild()
    except Exception:   # pragma: no cover - vanilla OpenSeesPy has no stamp
        return "n/a"


# --------------------------------------------------------------------- tests --
def test_elastic_state_tangent_is_ce():
    """Harness sanity: at an elastic state the engine tangent is Ce, and the
    numpy Ce transcription reproduces it. If this fails the harness (or the
    'tangent' response's row-major reshape) is wrong, not the fix.

    Same tiny-tail rule as the plastic legs: the elastic tangent is built from
    G(p) at the START of the last increment (measured on dee04dbe3: without the
    tail, 10 equal steps from zero leave Ct 9.5 % off Ce(p_final) -- G(0.9 p) =
    sqrt(0.9) G(p)). A 20-step tail of 1e-3 of the path puts the last start
    close to p_final; measured on dee04dbe3: Ct vs Ce(p_final) 4.97e-5."""
    build(2, 1)
    e_iso = np.diag([-2e-4, -2e-4, -2e-4])
    run(checkpoints([e_iso, e_iso*(1.0 + TAIL_FRAC)], [10, 20]))
    st = internal_state()
    K, G = get_elastic_moduli(st['sigma'])
    Ce = get_stiffness(K, G)
    Ct = engine_tangent()
    err, _ = rel_err_masked(Ct, Ce)
    print("\nladrunoBuild = %s; elastic Ct vs numpy Ce: max rel err %.2e" % (_build_hash(), err))
    assert err < 1e-4, (Ct, Ce)


@pytest.mark.parametrize("int_scheme", [2, 1], ids=["scheme2_BE", "scheme1_ME"])
@pytest.mark.parametrize("state", STATES, ids=[s[0] for s in STATES])
def test_tantype1_matches_one_sided_fd(state, int_scheme):
    """The engine's TanType 1 tangent at a plastic shear state matches the
    reversal-consistent one-sided FD of its own stress update to 1.5 %, and the
    workbench formula to 0.1 %. Scheme 2 reaches GetElastoPlasticTangent via
    BackwardEuler_CPPM's end-of-increment call; scheme 1 via the WP-110 (F15c)
    end-of-increment write in ModifiedEuler."""
    Ct, F0, st = base_run(state, int_scheme, tan_type=1)
    p = ONE3*trace(st['sigma'])   # internal convention, compression positive
    sd = get_state_dependent(st['sigma'], st['alpha'], st['fabric'], st['e'], st['alpha_in'])
    assert st['dgamma'] > 0, "%s: never went plastic -- retarget the path" % state[0]
    assert abs(sd['n'][3]) > 0.05, "%s: n carries no shear -- both defects would be invisible" % state[0]

    K, G = get_elastic_moduli(st['sigma'])
    C_wb = workbench_tangent(st['sigma'], st['dgamma'], K, G, sd)
    C_fd, signs = fd_one_sided(state, int_scheme, 1, F0, st['alpha_in'])

    err_fd, mask = rel_err_masked(Ct, C_fd)
    err_wb, _ = rel_err_masked(Ct, C_wb)
    np.set_printoptions(suppress=True, linewidth=140, precision=5)
    print("\n[%s, IntScheme %d] ladrunoBuild=%s p'=%.4g kPa n_xy=%.3f dGamma=%.3g signs=%s"
          % (state[0], int_scheme, _build_hash(), p, sd['n'][3], st['dgamma'], signs))
    print("  engine row3: %s" % Ct[3])
    print("  FD     row3: %s" % C_fd[3])
    print("  wb     row3: %s" % C_wb[3])
    print("  max rel err engine vs FD = %.3e (%d entries), vs workbench = %.3e"
          % (err_fd, int(mask.sum()), err_wb))

    assert err_wb < ORACLE_RTOL, (
        "engine tangent does not match the workbench formula -- the C++ fix and "
        "the tensor-algebra oracle disagree", err_wb, Ct, C_wb)
    assert err_fd < FD_RTOL, (
        "engine tangent is not the derivative of its own stress update to 1.5 %",
        err_fd, Ct, C_fd)


def test_tantype1_under_modified_euler_is_not_ce():
    """F15(c): before WP-110, ModifiedEuler never assigned aCep, so TanType 1 on
    IntScheme 1 handed the element a stale matrix -- Ce, at any state reached
    through ordinary substeps. At a plastic shear state it must now carry the
    plastic correction (here ~10 % of Ce's largest entry, measured by the
    workbench), not be Ce."""
    Ct, _, st = base_run(STATES[0], int_scheme=1, tan_type=1)
    assert st['dgamma'] > 0
    K, G = get_elastic_moduli(st['sigma'])
    Ce = get_stiffness(K, G)
    dev = float(np.max(np.abs(Ct - Ce)) / np.max(np.abs(Ce)))
    print("\nIntScheme 1 TanType 1: max|Ct - Ce| / max|Ce| = %.3e" % dev)
    assert dev > 0.05, ("TanType 1 under IntScheme 1 is still (close to) Ce -- "
                        "ModifiedEuler is not writing aCep", dev)
