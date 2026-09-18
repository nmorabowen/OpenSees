"""WP-110 / F15a -- measurement probe for ManzariDafalias::GetElastoPlasticTangent.

REVISION 2: a coordinator/Workbench cross-check found the first revision's
"numerator-only fix" still ~5-7% off finite-difference (FD) at realistic
states. That review surfaced TWO more findings, both folded in below:

  (a) a SECOND, independent defect in the DENOMINATOR
      (``temp3 = DoubleDot2_2_Contr(temp2, R) + Kp;`` at :5139) -- see
      ``get_ep_tangent_workbench`` and the "second defect" note further down;
  (b) the ORIGINAL central-difference FD was itself contaminated: it straddles
      a real branch in ``integrate()`` (the alpha_in reversal reset) whose sign
      depends on the probe direction, so ``(F(+h)-F(-h))/(2h)`` silently
      differenced two DIFFERENT internal states. ``safe_sign``/``fd_one_sided``
      below fix this by picking, per probe direction, whichever sign leaves
      ``alpha_in`` equal to the committed base state, then using a one-sided
      difference. With both the harness fixed and both Voigt defects corrected
      (the ``workbench`` formula), FD now matches to 0.02%-0.17% at realistic
      p'~20/~214 kPa, eta~1.3 states with genuine shear demand -- see
      Ladruno_implementation/_wp110_f15a_probe_results.md for the full table.

GOAL (measurement only, no C++ changes): confirm or refute that
``GetElastoPlasticTangent`` (SRC/material/nD/UWmaterials/ManzariDafalias.cpp:5110)
computes a wrong continuum elastoplastic tangent at a plastic (loading) state,
because the flow-direction vector ``R`` is converted contravariant -> covariant
TWICE:

    R = ToCovariant(temp0);                          // :5127 -- single, correct
    temp1 = DoubleDot4_2(aC, ToCovariant(R));         // :5129 -- R is covariant
                                                       //          ALREADY; this
                                                       //          re-wraps it,
                                                       //          doubling its
                                                       //          shear entries
                                                       //          a SECOND time.

``ToCovariant`` (:5651) only multiplies the shear entries (Voigt indices 3,4,5 =
xy,yz,xz) by 2 -- it leaves the normal entries (0,1,2 = xx,yy,zz) untouched. So
the bug can ONLY ever corrupt ``temp1``'s shear components; ``aCep = aC -
(Macauley/temp3) * Dyadic2_2(temp1, temp2)`` then has its correction term
row-scaled by ``temp1``, so the bug can ONLY corrupt aCep's SHEAR ROWS (Voigt
rows 3,4,5 -- dSigma_xy/dEps_*, dSigma_yz/dEps_*, dSigma_xz/dEps_*), leaving the
NORMAL rows (0,1,2) bit-identical between the as-written and fixed formulas.
That is the falsifiable, mechanism-specific prediction this probe checks.

METHOD. No C++ edit, no build (project rule -- this worktree is measurement
only). Two independent legs, cross-checked against each other:

  1. FE leg: a single 8-node SSPbrick (1-point stabilized hex, so the whole
     element shares one Gauss-point/material state) with ``ManzariDafalias``
     wrapped 3D (no plane-strain reduction, so the full 6-component internal
     state -- INCLUDING sigma_zz -- is directly readable via ``eleResponse``,
     which the 2D PlaneStrain wrapper does not expose). All 21 non-pinned DOFs
     are driven by SP constraints under a per-DOF ``Path`` TimeSeries
     (``constraints('Transformation')`` -- ``Plain`` silently drops any
     non-homogeneous SP value, discovered the hard way while building this
     harness) so the element sees an EXACT, fully prescribed affine
     (uniform-strain) history with zero free DOFs: no equilibrium solve can
     contaminate the probe. A finite-difference tangent of the wrapper's
     3D stress (6-comp) w.r.t. each of the 6 independent strain components is
     computed by rebuilding the model from scratch per probe direction/sign
     (revert-by-rebuild, replaying the identical committed path then bumping
     ONE more Path checkpoint) and taking a central difference; this converges,
     as h -> 0, to the TRUE continuum (rate-form) tangent regardless of which
     IntScheme produced the probe steps -- a small enough FD step samples the
     local rate response, not the scheme's own discretization.

  2. Numpy-oracle leg: ``GetElastoPlasticTangent`` (and the ``GetStateDependent``
     / ``GetElasticModuli`` / ``GetStiffness`` helpers it needs) is transcribed
     verbatim into plain numpy from the C++ source, in two variants -- AS
     WRITTEN (the double ``ToCovariant``) and FIXED (single ``ToCovariant``) --
     fed with the EXACT internal state (mSigma, mAlpha, mFabric, mAlpha_in,
     void ratio, dGamma) read back from the same committed FE state via
     ``eleResponse('stress'|'alpha'|'fabric'|'alpha_in'|'state')``. The 3D
     wrapper's ``getStress()``/``setTrialStrain()`` apply a UNIFORM ``-1.0``
     "geotechnical sign convention" flip to every stress/strain component
     (ManzariDafalias3D.cpp:80,122) -- since both sides of any dSigma/dEps
     ratio flip together, the tangent MATRIX is invariant under that flip, so
     no sign correction is needed to compare the numpy oracle (internal
     convention) against the wrapper-convention FD.

GOTCHA recorded for the ledger: the full ``ManzariDafalias`` constructor sets
``mElastFlag = 0`` ("stage 0") at construction, which makes ``integrate()``
call ``elastic_integrator`` UNCONDITIONALLY -- ``GetElastoPlasticTangent`` is
never reached, ``mDGamma`` stays 0, and TanType 0/1/2 are all silently
identical (aCep == aC). ``ops.updateMaterialStage('-material', tag, '-stage',
1)`` must be called before loading past the elastic nucleus; this cost real
time to find and is repeated here as its own regression-style assertion.

This is a MEASUREMENT probe, not a merge gate: it is intentionally UNMARKED
(no zone_a/zone_b) and skips (never fails/errors) when the installed Ladruno
binary this probe needs is not present, and it does not assert anything about
whether the bug is present -- it prints/records the comparison tables and a
verdict line so a human (or the next agent) can read them off directly.
Results are duplicated in Ladruno_implementation/_wp110_f15a_probe_results.md.
"""
import math
import os
import sys

import numpy as np
import pytest

BIN_DIR = r"C:\Program Files\Ladruno\OpenSees\bin"
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

ops = None
BUILD_HASH = None
_SKIP_REASON = None
if os.path.isdir(BIN_DIR):
    if BIN_DIR not in sys.path:
        sys.path.insert(0, BIN_DIR)
    try:
        import opensees as ops  # noqa: E402

        eng = os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__)))
        if eng != os.path.normcase(BIN_DIR):
            _SKIP_REASON = f"wrong engine bound: {ops.__file__}"
            ops = None
        else:
            BUILD_HASH = ops.ladrunoBuild()
    except Exception as exc:  # pragma: no cover - environment dependent
        _SKIP_REASON = f"could not import installed opensees: {exc!r}"
        ops = None
else:
    _SKIP_REASON = f"installed Ladruno binary directory not found: {BIN_DIR}"

pytestmark = pytest.mark.skipif(
    ops is None, reason=_SKIP_REASON or "installed Ladruno opensees.pyd unavailable"
)

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

# unit-cube 8-node hex, standard bottom-then-top node ordering
COORDS = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
FREE_DOFS = [(n, d) for n in range(2, 9) for d in (1, 2, 3)]
VOIGT_COMPS = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]


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


def to_covariant(v):
    r = v.copy(); r[3] *= 2.0; r[4] *= 2.0; r[5] *= 2.0
    return r


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
    # mElastFlag==1 ("stage 1", see module docstring GOTCHA) keeps the
    # sqrt(pn/Patm) pressure-dependence factor (ManzariDafalias.cpp:5028).
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
        D_factor = 1.0/(1.0+math.exp(7.6349 - 7.2713*101.0/MATP['P_atm']*p))
    else:
        D_factor = 1.0
    D *= D_factor
    B = 1.0+1.5*(1-MATP['c'])/MATP['c']*g_fun(cos3theta, MATP['c'])*cos3theta
    C = 3.0*math.sqrt(1.5)*(1-MATP['c'])/MATP['c']*g_fun(cos3theta, MATP['c'])
    R = B*n - C*(single_dot(n, n)-ONE3*I1) + ONE3*D*I1
    return dict(n=n, d=d, b=b, h=h, B=B, C=C, D=D, R=R)


def _get_ep_tangent(next_stress, next_dgamma, G, K, B, C, D, h, n, d, b, buggy):
    p = ONE3*trace(next_stress) + M_PRESIDUAL
    p = (SMALL+M_PRESIDUAL) if p < (SMALL+M_PRESIDUAL) else p
    r = devpart(next_stress)/p
    Kp = TWO3*p*h*dd22_contr(b, n)
    aC = get_stiffness(K, G)
    temp0 = B*n - C*(single_dot(n, n)-ONE3*I1) + ONE3*D*I1
    R = to_covariant(temp0)
    if buggy:
        temp1 = aC @ to_covariant(R)   # AS-WRITTEN (ManzariDafalias.cpp:5129): double ToCovariant
    else:
        temp1 = aC @ R                 # FIXED: single ToCovariant (R is already covariant)
    temp0b = to_covariant(n - ONE3*dd22_contr(n, r)*I1)
    temp2 = temp0b @ aC                # DoubleDot2_4(v, m) == v^T . m
    temp3 = dd22_contr(temp2, R) + Kp
    if abs(temp3) < SMALL:
        return aC
    return aC - (macauley_index(next_dgamma)/temp3) * np.outer(temp1, temp2)


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


def get_ep_tangent_workbench(next_stress, next_dgamma, G, K, B, C, D, h, n_vec, b_vec):
    """SECOND defect check: the Workbench's (C4-sanisand/sanisand_cep.py,
    ``Dep6``) formula, using proper full-tensor contractions throughout (no
    Voigt ToCovariant/DoubleDot2_2_Contr detours at all):

        CeR   = Ce:R    = 2G dev(R) + K tr(R) I
        QCe   = Q:Ce    = 2G dev(Q) + K tr(Q) I,   Q = n - (n:r)/3 I
        denom = Kp + Q:Ce:R   (a TRUE Frobenius double-dot)
        Dep   = Ce - CeR (x) QCe / denom

    Measured: the engine's own ``temp3`` (DoubleDot2_2_Contr(temp2, R) + Kp,
    :5139) is a near-constant ~1.076x too LARGE relative to `denom` here
    whenever R/n carry nonzero shear (Voigt 3-5) components -- a second,
    independent Voigt-convention mismatch from the numerator's double
    ToCovariant, surviving even after that first fix. It rescales the WHOLE
    plastic correction term uniformly (same ratio on every entry, diagonal
    and off-diagonal alike), because it is a single scalar denominator."""
    S = _vec_to_ten(next_stress)
    p = np.trace(S)/3.0 + M_PRESIDUAL
    p = (SMALL+M_PRESIDUAL) if p < (SMALL+M_PRESIDUAL) else p
    r_t = _dev3(S)/p
    n_t = _vec_to_ten(n_vec)
    b_t = _vec_to_ten(b_vec)
    Kp = TWO3*p*h*_fro(b_t, n_t)
    R_t = B*n_t - C*(n_t@n_t - np.eye(3)/3.0) + D/3.0*np.eye(3)
    Q_t = n_t - _fro(n_t, r_t)/3.0*np.eye(3)
    CeR = 2.0*G*_dev3(R_t) + K*np.trace(R_t)*np.eye(3)
    QCe = 2.0*G*_dev3(Q_t) + K*np.trace(Q_t)*np.eye(3)
    denom = Kp + _fro(QCe, R_t)
    aC = get_stiffness(K, G)
    if abs(denom) < SMALL:
        return aC, denom
    aCep = aC - macauley_index(next_dgamma)*np.outer(_ten_to_vec(CeR), _ten_to_vec(QCe))/denom
    return aCep, denom


# --------------------------------------------------------------- FE harness --
_tagctr = [9000]


def _fresh_tag():
    _tagctr[0] += 1
    return _tagctr[0]


def build(mat_tag=1, int_scheme=2, tan_type=1):
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
    # See module docstring GOTCHA: mElastFlag defaults to 0 (force-elastic) at
    # construction; without this the whole probe would silently measure Ce.
    ops.updateMaterialStage('-material', mat_tag, '-stage', 1)
    ops.constraints('Transformation')  # 'Plain' silently drops non-homogeneous SPs
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-9, 30)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    _tagctr[0] = 9000


def affine_disp(node, E):
    x, y, z = COORDS[node]
    return E @ np.array([x, y, z])


def strain_to_disps(E):
    d = {}
    for n in range(2, 9):
        u = affine_disp(n, E)
        d[(n, 1)] = u[0]; d[(n, 2)] = u[1]; d[(n, 3)] = u[2]
    return d


def make_checkpoints(e_targets, nsteps):
    waypoints = [np.zeros((3, 3))] + list(e_targets)
    vals = {k: [0.0] for k in FREE_DOFS}
    for seg in range(len(waypoints)-1):
        a, b = waypoints[seg], waypoints[seg+1]
        for i in range(1, nsteps+1):
            s = a + (b-a)*i/nsteps
            dd = strain_to_disps(s)
            for k in FREE_DOFS:
                vals[k].append(dd[k])
    return vals


def unit_dE(idx):
    a, b = VOIGT_COMPS[idx]
    dE = np.zeros((3, 3))
    if a == b:
        dE[a, b] = 1.0
    else:
        dE[a, b] = 0.5; dE[b, a] = 0.5
    return dE


def make_checkpoints_strain(e_targets, nsteps, probe_idx, h):
    vals = make_checkpoints(e_targets, nsteps)
    e_final = list(e_targets)[-1] + h*unit_dE(probe_idx)
    dd = strain_to_disps(e_final)
    for k in FREE_DOFS:
        vals[k].append(dd[k])
    return vals


def apply_path_and_run(vals):
    n_points = len(next(iter(vals.values())))
    times = list(range(n_points))
    for (node, dof), series in vals.items():
        tag = _fresh_tag()
        ops.timeSeries('Path', tag, '-time', *[float(t) for t in times],
                        '-values', *[float(v) for v in series])
        ops.pattern('Plain', tag, tag)
        ops.sp(node, dof, 1.0)
    for _ in range(n_points - 1):
        ret = ops.analyze(1)
        if ret != 0:
            raise RuntimeError(f"analyze failed, ret={ret}")


def wrapper_stress():
    return np.array(ops.eleResponse(1, 'stress'))  # 6-comp, = -mSigma


def internal_state():
    m_sigma = -1.0*wrapper_stress()
    m_alpha = np.array(ops.eleResponse(1, 'alpha'))
    m_fabric = np.array(ops.eleResponse(1, 'fabric'))
    m_alpha_in = np.array(ops.eleResponse(1, 'alpha_in'))
    state = np.array(ops.eleResponse(1, 'state'))
    return m_sigma, m_alpha, m_fabric, m_alpha_in, state[24], state[25]


def p_q_eta(sig_w):
    p = -ONE3*trace(sig_w)
    s = devpart(sig_w)
    q = math.sqrt(1.5*dd22_contr(s, s))
    return p, q, (q/p if p != 0 else float('nan'))


def fd_material_tangent(e_targets, nsteps, h, int_scheme, tan_type, mat_tag=1):
    C = np.zeros((6, 6))
    for j in range(6):
        fs = {}
        for sign in (+1, -1):
            build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
            vals = make_checkpoints_strain(e_targets, nsteps, probe_idx=j, h=sign*h)
            apply_path_and_run(vals)
            fs[sign] = wrapper_stress()
        C[:, j] = (fs[+1] - fs[-1]) / (2*h)
    return C


def committed_alpha_in(e_targets, nsteps, int_scheme, tan_type, mat_tag=1):
    build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
    apply_path_and_run(make_checkpoints(e_targets, nsteps))
    return internal_state()[3]


def safe_probe_sign(e_targets, nsteps, h, int_scheme, tan_type, probe_idx, base_alpha_in, mat_tag=1):
    """HARNESS FIX for the reversal-branch contamination (see module
    docstring REVISION 2 note b): return whichever of +1/-1 leaves alpha_in
    equal to the committed base state's, i.e. does NOT trip integrate()'s
    reversal reset (~ManzariDafalias.cpp:1005-1013). A one-sided FD in that
    direction samples one smooth branch instead of straddling the kink."""
    for sign in (+1, -1):
        build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
        vals = make_checkpoints_strain(e_targets, nsteps, probe_idx=probe_idx, h=sign*h)
        apply_path_and_run(vals)
        if np.allclose(internal_state()[3], base_alpha_in, atol=1e-9):
            return sign
    return None


def fd_material_tangent_one_sided(e_targets, nsteps, h, int_scheme, tan_type, signs, mat_tag=1):
    build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
    apply_path_and_run(make_checkpoints(e_targets, nsteps))
    F0 = wrapper_stress()
    C = np.zeros((6, 6))
    for j in range(6):
        build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
        vals = make_checkpoints_strain(e_targets, nsteps, probe_idx=j, h=signs[j]*h)
        apply_path_and_run(vals)
        F1 = wrapper_stress()
        C[:, j] = (F1 - F0) / (signs[j]*h)
    return C


def analytic_tangent_from_state(e_targets, nsteps, int_scheme, tan_type, mat_tag=1, buggy=True):
    build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
    vals = make_checkpoints(e_targets, nsteps)
    apply_path_and_run(vals)
    m_sigma, m_alpha, m_fabric, m_alpha_in, e, dgamma = internal_state()
    K, G = get_elastic_moduli(m_sigma)
    sd = get_state_dependent(m_sigma, m_alpha, m_fabric, e, m_alpha_in)
    aCep = _get_ep_tangent(m_sigma, dgamma, G, K, sd['B'], sd['C'], sd['D'], sd['h'],
                            sd['n'], sd['d'], sd['b'], buggy=buggy)
    return aCep, dict(mSigma=m_sigma, e=e, dGamma=dgamma, K=K, G=G, **sd)


def analytic_tangents_all(e_targets, nsteps, int_scheme, tan_type, mat_tag=1):
    """buggy (as-written) / fixed (numerator only) / workbench (both defects
    corrected) in one committed build, plus the state dict."""
    build(mat_tag=mat_tag, int_scheme=int_scheme, tan_type=tan_type)
    apply_path_and_run(make_checkpoints(e_targets, nsteps))
    m_sigma, m_alpha, m_fabric, m_alpha_in, e, dgamma = internal_state()
    K, G = get_elastic_moduli(m_sigma)
    sd = get_state_dependent(m_sigma, m_alpha, m_fabric, e, m_alpha_in)
    Ca_bug = _get_ep_tangent(m_sigma, dgamma, G, K, sd['B'], sd['C'], sd['D'], sd['h'],
                              sd['n'], sd['d'], sd['b'], buggy=True)
    Ca_fix = _get_ep_tangent(m_sigma, dgamma, G, K, sd['B'], sd['C'], sd['D'], sd['h'],
                              sd['n'], sd['d'], sd['b'], buggy=False)
    Ca_wb, denom_wb = get_ep_tangent_workbench(m_sigma, dgamma, G, K, sd['B'], sd['C'], sd['D'],
                                                sd['h'], sd['n'], sd['b'])
    return dict(Ca_bug=Ca_bug, Ca_fix=Ca_fix, Ca_wb=Ca_wb, mSigma=m_sigma, e=e, dGamma=dgamma,
                K=K, G=G, alpha_in=m_alpha_in, denom_wb=denom_wb, **sd)


def _fmt(M):
    return "\n".join(" ".join(f"{x:11.4g}" for x in row) for row in M)


# --------------------------------------------------------------------- test --
def test_manzari_ep_tangent_probe():
    """Records (does not gate on) the ManzariDafalias GetElastoPlasticTangent
    measurement. See module docstring for method and the companion results
    note Ladruno_implementation/_wp110_f15a_probe_results.md for the table."""
    np.set_printoptions(suppress=True, linewidth=140)
    lines = []

    def out(s=""):
        lines.append(str(s))

    out(f"ladrunoBuild = {BUILD_HASH}")
    assert BUILD_HASH is not None and len(BUILD_HASH) >= 7, "no build hash read back"

    # --- elastic sanity gate: this one DOES assert, since it validates the
    # harness itself (FE kinematics + the numpy get_elastic_moduli/get_stiffness
    # transcription), not the suspected bug.
    out("="*70); out("ELASTIC SANITY (TanType=0)"); out("="*70)
    e_elastic = [np.diag([-2e-4, -2e-4, -2e-4])]
    Ce_analytic, _ = analytic_tangent_from_state(e_elastic, 10, int_scheme=2, tan_type=0, buggy=True)
    Ce_fd = fd_material_tangent(e_elastic, 10, 1e-7, int_scheme=2, tan_type=0)
    out("analytic Ce:\n" + _fmt(Ce_analytic))
    out("FD:\n" + _fmt(Ce_fd))
    max_rel = np.max(np.abs(Ce_fd-Ce_analytic)/np.where(np.abs(Ce_analytic) > 1e-6, np.abs(Ce_analytic), 1.0))
    out(f"elastic sanity max rel err = {max_rel:.3e}")
    assert max_rel < 1e-4, "harness sanity check FAILED: FD does not reproduce Ce -- fix the harness, not the bug"

    # Realistic states (p'~20 and ~214 kPa, eta~1.3) WITH genuine shear demand
    # (n's Voigt-shear component ~0.21) -- a pure coaxial/no-shear triaxial
    # path makes R's shear component identically zero and hides BOTH defects,
    # which is what an earlier (unshared) attempt at "realistic states"
    # accidentally measured. iso/yy/xy tuned empirically against this harness.
    cases = [
        ("LOW p'~20, eta~1.3", -1.5e-5, -2.7e-3, -5e-4),
        ("HIGH p'~214, eta~1.3", -4.8e-5, -8.64e-3, -1.6e-3),
    ]
    for label, iso, yy, xy in cases:
        out("="*70); out(f"{label} (IntScheme=2 BackwardEuler_CPPM, TanType=1)"); out("="*70)
        e_targets = [np.diag([iso, iso, iso])]
        e2 = np.diag([iso, iso, iso]).copy()
        e2[1, 1] += yy; e2[0, 1] += xy; e2[1, 0] += xy
        e_targets.append(e2)
        nsteps = 40

        st = analytic_tangents_all(e_targets, nsteps, int_scheme=2, tan_type=1)
        p, q, eta = p_q_eta(-1.0*st['mSigma'])
        out(f"p'={p:.6g} kPa  q={q:.6g}  eta={eta:.4g}  dGamma={st['dGamma']:.4g}"
            f"  n_shear={st['n'][3]:.4g}")
        assert st['dGamma'] > 0, f"{label}: material never went plastic -- retarget the strain path"
        assert abs(st['n'][3]) > 0.05, f"{label}: n has ~no shear component -- both defects would be invisible"

        out("Analytic aCep AS-WRITTEN (buggy, double ToCovariant) row3:\n  " + str(st['Ca_bug'][3]))
        out("Analytic aCep FIXED (numerator only) row3:\n  " + str(st['Ca_fix'][3]))
        out("Analytic aCep WORKBENCH (numerator + denominator) row3:\n  " + str(st['Ca_wb'][3]))

        # Mechanism 1 (numerator, exact 2x): off-diagonal shear-row entries.
        offdiag_ratio = st['Ca_bug'][3, 0:3] / st['Ca_fix'][3, 0:3]
        out(f"buggy/fixed ratio, row3 cols0-2 (expect exactly 2.0): {offdiag_ratio}")
        assert np.allclose(offdiag_ratio, 2.0, rtol=1e-6), "numerator defect signature (exact 2x) not reproduced"

        # Mechanism 2 (denominator, ~1.076x): fixed and workbench share the
        # SAME numerator (temp1/temp2 equivalent to Ce:R / Q:Ce) and differ
        # ONLY in the denominator scalar, so (aC-Ca_wb)/(aC-Ca_fix) must be
        # the SAME constant across every entry of the correction term.
        aC = get_stiffness(st['K'], st['G'])
        corr_fix = aC - st['Ca_fix']
        corr_wb = aC - st['Ca_wb']
        mask = np.abs(corr_fix) > 1.0
        denom_ratio = corr_wb[mask] / corr_fix[mask]
        out(f"denominator-only ratio (correction_wb/correction_fixed), should be one "
            f"constant across all entries: min={denom_ratio.min():.5f} max={denom_ratio.max():.5f}")
        assert (denom_ratio.max() - denom_ratio.min()) < 1e-3 * abs(denom_ratio.mean()), (
            "workbench and fixed differ by more than a single denominator rescale -- "
            "second-defect isolation broken")

        # Harness fix: one-sided, reversal-consistent FD. This is the
        # regression-style, tightly-toleranced check.
        base_ain = st['alpha_in']
        signs = [safe_probe_sign(e_targets, nsteps, 1e-6, 2, 1, j, base_ain) for j in range(6)]
        out(f"reversal-safe probe signs per direction: {signs}")
        Cfd_1sided = fd_material_tangent_one_sided(e_targets, nsteps, 1e-8, 2, 1, signs)
        out("FD one-sided (h=1e-8, reversal-safe) row3:\n  " + str(Cfd_1sided[3]))

        # Two entries in this row (yz, xz -- indices 4,5) carry NO real signal
        # (Ca_wb ~ 0 exactly -- no yz/xz loading imposed) and h=1e-8 finite
        # differencing of double-precision stress leaves ~1-3 units of pure
        # round-off there (Cfd_1sided ~ 3, not ~0), which is noise, not a
        # mismatch -- checking them against ANY nonzero floor divides noise by
        # noise. Only assert on the entries carrying a real analytic signal.
        real_signal = np.abs(st['Ca_wb'][3]) > 1.0
        rel_vs_wb = np.abs(Cfd_1sided[3] - st['Ca_wb'][3]) / np.where(real_signal, np.abs(st['Ca_wb'][3]), 1.0)
        rel_vs_fix = np.abs(Cfd_1sided[3] - st['Ca_fix'][3]) / np.where(real_signal, np.abs(st['Ca_fix'][3]), 1.0)
        out(f"row3 rel err, one-sided FD vs WORKBENCH: {rel_vs_wb}  (real-signal mask: {real_signal})")
        out(f"row3 rel err, one-sided FD vs FIXED (numerator only): {rel_vs_fix}")
        out(f"no-signal entries (yz,xz) absolute FD noise: {Cfd_1sided[3][~real_signal]}")
        assert np.max(rel_vs_wb[real_signal]) < 0.01, (
            f"{label}: one-sided FD does not match the workbench (both-defects-fixed) "
            f"formula to 1% -- rel err {rel_vs_wb}")

        # Context only (not asserted): the ORIGINAL central-difference FD,
        # which straddles the alpha_in reversal branch and is expected to be
        # noisy/wrong on the off-diagonal columns.
        Cfd_central = fd_material_tangent(e_targets, nsteps, 1e-8, int_scheme=2, tan_type=1)
        out("FD central-difference (h=1e-8, KNOWN CONTAMINATED, for context only) row3:\n  "
            + str(Cfd_central[3]))

    report = "\n".join(lines)
    print(report)
    log_path = os.path.join(os.path.dirname(__file__), "..",
                             "Ladruno_implementation", "_wp110_f15a_probe_raw_output.txt")
    try:
        with open(log_path, "w") as fh:
            fh.write(report + "\n")
    except OSError:
        pass
