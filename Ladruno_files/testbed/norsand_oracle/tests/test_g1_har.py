"""WP-144 round 3b, G1 (Zone B): the HAR energy option `-energy HAR` on both oracles.  Equation sheet 144a 2.3, 2.4, 13.1h, 13.11.

TEST RULE (plan 5): every expected value and tolerance in this file comes from (a) a closed form of the sheet (coded in
g1_hf_common.py from the sheet's formulas, section cited at each use), (b) a number PRINTED in the sheet (compared to the
printed digits), or (c) a stated convergence / truncation argument; all of it was written BEFORE the oracles were run on
these paths.  No oracle output is an expected value.  A failing assertion is a finding (against the sheet or an oracle),
never a reason to edit an oracle or loosen a bound.

Gates (ids follow the sheet's K1 numbering, plus H for the FD / census re-runs sheet 2.4 lists)
  K1.1h   isotropic compression p(eps_v) = -p_a [1 - k(1-n) eps_v]^(1/(1-n)), TIMs printed values; general n incl. n = 0;
          the domain edge eps_v = 1/(k(1-n)) is the refusal when p_min = 0                         sheet 2.3, 13.1h
  K1.11   constant-volume elastic shear from eps^e = 0: p(eps_s), q(eps_s), eta = 3 g eps_s EXACTLY, printed values;
          BA06 control (p constant) so the discriminator has power                                  sheet 13.11
  K1.11i  inverse map (S.5h''): the printed spot values (p = -3.5, eta = 2.1) and a seeded round trip  sheet 2.3, 13.11
  K1.2h   closed non-coaxial elastic loop under HAR: zero net work (1e-12), state returns; control: the open segment
          stores Psi(E1) - Psi(0) of (S.4h)                                                         sheet 13.2, 2.3
  H.hess  the 6 x 6 elastic tangent against the sheet's gate table (lambda_min(6-D)/K_iso, 9 values of eta), the 4-fold
          2G eigenvalue and det D = 3 k g p_a^2 (varpi/p_a)^(2n) (S.5c), through the public tangent  sheet 2.3
  H.jac   (S.30) Jacobian vs central FD of (S.29), HAR, fork and paper CSL                           sheet 2.4 FD list
  H.cto   (S.33) consistent tangent vs central FD of the O2 stress update, HAR, off the corners, three shears  sheet 2.4
  H.fin   (S.34) finite-strain tangent vs central FD of P = tau (1 + hE)^-T, HAR; the half-spin and no-tau(+)1 variants are
          negative controls that must miss the bound by >= 10x                                       sheet 2.4, 9.5
  H.chain (S.45)-(S.47) chained tangent of substepped increments vs central FD, HAR, paper + fork CSL, m = 2..8, plus
          non-uniform fractions                                                                       sheet 2.4, 9.6
  H.conv  O2 -> O1 first-order convergence under HAR                                                  plan 5.1
  H.diss  D >= 0 census under HAR on the four convergence paths (both oracles) and 200 random increments (O2)  sheet 2.4, 11
  H.ref   the parser refusals of sheet 2.4 and the defaults of 1.3 / 9.7

MUTANTS (named for the mutation gate, plan 5.3).  "HAR silently replaced by BA06" is killed by K1.1h, K1.11 (p is constant
under BA06), K1.11i, H.hess and by the O2-parity of the G2 files, and is NOT killed by any FD test here (an FD test compares
a code with itself: a consistent swap passes its own FD check; sheet 2.4, round 3b A4).  The (S.3) t4 term dropped / AB06 eq 64
used: H.hess (eigenvalue table), H.cto (a wrong a^e is FD-inconsistent with the stress).  D12 dropped: H.hess, K1.11 (eta = 3 g eps_s
needs the stress-induced |p| growth).  n hard-coded to 1/2: K1.1h general-n.  p_a (HAR) taken from the CSL default
(101.325) instead of the one shared flag: every p_a = 101 closed form here.  p_ref = |p0| instead of p_a in the F_tol / r4 scalings
or the p_min default: H.ref default test and the p_min = 0.505 floor tests of test_g1_floor.py.
"""
import functools
import math
import warnings

import numpy as np
import pytest

import g1_hf_common as C
import test_g1_convergence_tangents as CT
import test_g1_finite_strain_k2 as FS
import test_g1_substep_tangent as ST
from conftest import ORACLES, make_params

O1, O2 = ORACLES["O1"], ORACLES["O2"]
KER = O2.kernel
I3 = np.eye(3)

# a plastic-free yield surface far from every state of the elastic tests (N = 0.4: eta_surface -> M/N = 3.33)
PI_FAR = -5000.0


def hf_state(oname, kw, sigma0, pi0, psi0=-0.05, v0=None):
    """(Params, State) at sigma0 with the fork/paper CSL v0 built for the target psi_i0 at pi_i0 (sheet 6, S.22)."""
    P = C.make_hf(oname, **kw)
    v0 = C.v_for_psi(kw, pi0, psi0) if v0 is None else v0
    return P, ORACLES[oname].initial_state(P, sigma0, v0, pi0)


def run_deps(oname, P, st0, deps):
    ora = ORACLES[oname]
    deps = np.asarray(deps, float)
    return ora.run_path(P, st0, deps, rtol=1e-10) if oname == "O1" else ora.run_path(P, st0, deps)


def pressure(st):
    return float(np.trace(st.sigma)) / 3.0


def eps_vs(st):
    """(eps_v, eps_s) of the elastic strain tensor: eps_v = tr, eps_s = sqrt(2/3) ||dev||."""
    ev = float(np.trace(st.eps_e))
    return ev, C.SQ23 * float(np.linalg.norm(st.eps_e - ev / 3.0 * I3))


def bulk_modulus(C4):
    """K = C_iijj / 9 (the hydrostatic modulus of a 4th-order tangent)."""
    return float(np.einsum("iijj->", C4)) / 9.0


def is_plastic(st):
    return bool(st.flags.get("plastic", False))


# ==============================================================================================
# K1.1h  isotropic compression under HAR (sheet 2.3, 13.1h)
# ==============================================================================================
@pytest.mark.parametrize("pmin", [0.0, "default"])
def test_k1_1h_isotropic_compression_matches_closed_form_and_printed_values(oracle_name, pmin):
    """p(eps_v) = -p_a [1 - k(1-n) eps_v]^(1/(1-n)) (13.1h), K = k p_a (|p|/p_a)^n, TIMs set.  Path eps_v: 0 -> -1e-3 -> +5e-4
    from eps^e = 0 (p = -p_a, the sheet's origin shift).  The sheet prints p(-1e-3) = -381.983585, K = 371129.585,
    p(+5e-4) = -28.117707 (compared to the printed digits: 2e-9, 5e-9, 5e-8 relative).  Closed form on every step to 1e-10
    (elastic exact algebra; O1 rtol 1e-10).  The path never comes near the floor (|p| >= 28 >> 0.505): the result with the
    default floor must be identical to p_min = 0 (a floor that acts above p_min, or a default that breaks the
    option, is caught).  KILLS: HAR replaced by BA06; n hard-coded; p_a left at the CSL default 101.325; k(1-n) mis-scaled;
    a floor active away from p_min."""
    kw = C.har_kw(p_min=pmin)
    P, st0 = hf_state(oracle_name, kw, -C.PA_T * I3, PI_FAR)
    assert abs(np.trace(st0.eps_e)) <= 1e-12, "premise: eps^e = 0 at p = -p_a (sheet 2.3 origin shift)"
    steps = [-1e-4] * 10 + [1e-4] * 15
    deps = np.array([np.eye(3) * s / 3.0 for s in steps])
    sts = run_deps(oracle_name, P, st0, deps)
    cum = np.cumsum(steps)
    for k, (s, c) in enumerate(zip(sts, cum)):
        assert C.ok_state(oracle_name, s), (k, C.why_state(oracle_name, s))
        assert not is_plastic(s), f"step {k} yielded: not an elastic path"
        p_ref, _ = C.har_pq(c, 0.0)
        assert abs(pressure(s) - p_ref) <= 1e-10 * abs(p_ref), (k, pressure(s), p_ref)
        assert np.linalg.norm(s.sigma - pressure(s) * I3) <= 1e-9 * abs(p_ref), "isotropic strain -> no deviatoric stress"
        if oracle_name == "O2":
            assert (s.n_f_tr, s.n_f_post, s.n_f_init) == (0, 0, 0), "floor event on a path far above p_min"
    # printed values of 13.1h
    s10, s25 = sts[9], sts[24]
    assert abs(pressure(s10) - (-381.983585)) <= 2e-9 * 381.983585
    assert abs(pressure(s25) - (-28.117707)) <= 5e-8 * 28.117707
    K10 = bulk_modulus(ORACLES[oracle_name].tangent(P, s10))
    assert abs(K10 - 371129.585) <= 5e-9 * 371129.585, K10
    assert abs(K10 - C.har_Kiso(pressure(s10))) <= 1e-9 * K10          # K = k p_a (|p|/p_a)^n at q = 0 (S.5h')


@pytest.mark.parametrize("n", [0.0, 0.3, 0.7])
def test_k1_1h_closed_form_for_general_n_including_n_zero(oracle_name, n):
    """13.1h for n in {0, 0.3, 0.7}, k = 1200, g = 600, p_a = 101 (not the TIMs numbers).  n = 0 is HAR05 eq 22 with the
    shift: p = -p_a (1 - k eps_v) EXACT (sheet 2.3), asserted against the oracle AND against the general formula.  Path
    eps_v: 0 -> -1e-3 -> 0.7 x edge, edge = 1/(k(1-n)).  p_min = 0 (the n = 0.7 end point is at |p| = 1.8 kPa).
    KILLS: n hard-coded to 1/2 (the TIMs value); an exponent 1/n or n/(1-n) swapped in p; k(1-n) scaled by k only."""
    k, g = 1200.0, 600.0
    kw = C.har_kw(k=k, g=g, n_e=n, p_min=0.0)
    edge = 1.0 / (k * (1.0 - n))
    P, st0 = hf_state(oracle_name, kw, -C.PA_T * I3, PI_FAR)
    top = 0.7 * edge
    steps = [-2.5e-4] * 4 + [(top + 1e-3) / 16.0] * 16
    deps = np.array([np.eye(3) * s / 3.0 for s in steps])
    sts = run_deps(oracle_name, P, st0, deps)
    cum = np.cumsum(steps)
    for kk, (s, c) in enumerate(zip(sts, cum)):
        assert C.ok_state(oracle_name, s), (kk, C.why_state(oracle_name, s))
        p_ref, _ = C.har_pq(c, 0.0, k, g, n, C.PA_T)
        if n == 0.0:
            assert abs(p_ref - (-C.PA_T * (1.0 - k * c))) <= 1e-12 * abs(p_ref)       # sheet 2.3: n = 0 reduces to eq 22
        assert abs(pressure(s) - p_ref) <= 1e-10 * abs(p_ref), (kk, pressure(s), p_ref)
        assert not is_plastic(s)
    assert abs(cum[-1] - top) <= 1e-12


def test_k1_1h_domain_edge_is_the_refusal_when_pmin_is_zero(oracle_name):
    """13.1h / 2.4: p = 0 at eps_v = 1/(k(1-n)) = 1.058491699e-3 (printed), 'the domain edge, refused' when p_min = 0.  A
    step to edge - 1e-5 (|p| = 9e-3 kPa, above O1's p -> 0 event stop at 1e-6 p_ref = 1e-4 kPa) is accepted and p equals the
    closed form; a step past the edge, edge + 1e-5, is refused (O2 flags refused with a domain reason; O1 reports a non-ok
    status).  The edge value is the printed one.  No assertion on the stress carried by a refused state: the sheet fixes only
    that the elastic() evaluation reports a failure code, never NaN, inside the kernel (2.4 table); what a driver freezes is the
    oracle's own convention (O2 carries NaN there, reported in the round-3b test report, not asserted).
    KILLS: a missing domain test (the HAR formula past the edge is the mirror image, p > 0, and the trial would be accepted
    with a tensile pressure); a refusal threshold that fires inside the domain."""
    assert abs(C.EDGE_T - 1.058491699e-3) <= 1e-12
    kw = C.har_kw(p_min=0.0)
    P, st0 = hf_state(oracle_name, kw, -C.PA_T * I3, PI_FAR)
    ok_step = run_deps(oracle_name, P, st0, [np.eye(3) * (C.EDGE_T - 1e-5) / 3.0])[-1]
    assert C.ok_state(oracle_name, ok_step), C.why_state(oracle_name, ok_step)
    p_cf = C.har_pq(C.EDGE_T - 1e-5, 0.0)[0]
    assert abs(pressure(ok_step) - p_cf) <= 1e-9 * abs(p_cf), (pressure(ok_step), p_cf)
    P, st0 = hf_state(oracle_name, kw, -C.PA_T * I3, PI_FAR)
    bad = run_deps(oracle_name, P, st0, [np.eye(3) * (C.EDGE_T + 1e-5) / 3.0])[-1]
    assert not C.ok_state(oracle_name, bad), "a trial outside the HAR domain was accepted with p_min = 0"
    if oracle_name == "O2":
        assert "domain" in str(bad.flags["reason"]), bad.flags["reason"]


# ==============================================================================================
# K1.11  constant-volume elastic shear under HAR (sheet 13.11)
# ==============================================================================================
DEV_TXC = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)     # compression corner (theta = pi/3), zeta = 1: the lowest F


def shear_increments(es_targets):
    """Diagonal tensor increments with tr = 0 and cumulative eps_s = es_targets[j]: e = sqrt(3/2) eps_s n_hat."""
    out, prev = [], 0.0
    for t in es_targets:
        out.append(np.diag(C.SQ32 * (t - prev) * DEV_TXC))
        prev = t
    return out


ES_STEPS = [2e-4, 4e-4, 6e-4, 8e-4, 8.66546968e-4, 1e-3]


def test_k1_11_constant_volume_shear_stress_growth_and_eta_equals_3g_eps_s(oracle_name):
    """13.11: at eps_v = 0 from eps^e = 0, p(eps_s) = -p_a [1 + 3 g k(1-n) eps_s^2]^(n/(2(1-n))), q = 3 g p_a eps_s [.]^same,
    eta = q/|p| = 3 g eps_s EXACTLY.  Printed (TIMs): eps_s = 8.66546968e-4 -> eta = 2.1, p = -166.548671, q = 349.752210;
    eps_s = 1e-3 -> p = -183.183351, q = 443.928664, eta = 2.42341163.  Path: six elastic steps (eps_s 2e-4 ... 1e-3) of
    pure (tr = 0) shear along the compression meridian; closed form at every step to 1e-10; printed values to their
    digits (3e-9, 3e-9, 5e-9 relative).  The surface is made inaccessible (pi_i0 = -5000), asserted by 'never plastic'.
    The control: the same path under BA06 (K2, alpha0 = 0) has p = p0 exp(-eps_v/kappa) = -100 CONSTANT and eta = 3 mu0 eps_s/|p|:
    the discriminating check between the two energies (sheet 13.11).
    KILLS: HAR replaced by BA06 (p constant); D12 / the stress-induced |p| growth dropped; q computed from eps_s without the
    w factor; n = 1/2 hard-coded exponents (this is also a n = 1/2 closed form, so the general-n gate is K1.1h)."""
    kw = C.har_kw(p_min=0.0)
    P, st0 = hf_state(oracle_name, kw, -C.PA_T * I3, PI_FAR)
    sts = run_deps(oracle_name, P, st0, shear_increments(ES_STEPS))
    for t, s in zip(ES_STEPS, sts):
        assert C.ok_state(oracle_name, s) and not is_plastic(s)
        p_ref, q_ref = C.har_pq(0.0, t)
        w = np.linalg.eigvalsh(s.sigma)
        p, q = float(w.sum()) / 3.0, C.SQ32 * float(np.linalg.norm(w - w.mean()))
        assert abs(p - p_ref) <= 1e-10 * abs(p_ref), (t, p, p_ref)
        assert abs(q - q_ref) <= 1e-10 * q_ref, (t, q, q_ref)
        assert abs(q / abs(p) - 3.0 * C.G_T * t) <= 1e-10 * 3.0 * C.G_T * t          # eta = 3 g eps_s exactly
        assert abs(np.trace(s.eps_e)) <= 1e-12, "constant volume"
    s5, s6 = sts[4], sts[5]
    for s, (eta, p_pr, q_pr) in ((s5, (2.1, -166.548671, 349.752210)), (s6, (2.42341163, -183.183351, 443.928664))):
        w = np.linalg.eigvalsh(s.sigma)
        p, q = float(w.sum()) / 3.0, C.SQ32 * float(np.linalg.norm(w - w.mean()))
        assert abs(q / abs(p) - eta) <= 3e-9 * eta
        assert abs(p - p_pr) <= 3e-9 * abs(p_pr) + 5e-7 and abs(q - q_pr) <= 5e-9 * q_pr + 5e-7


def test_k1_11_control_under_ba06_p_is_constant_so_the_check_discriminates(oracle_name):
    """The BA06 control of 13.11: same strain path under the K2 energy (alpha0 = 0): p = p0 e^omega = -100 exactly, q = 3 mu0 eps_s.
    This is not a HAR test: it proves the HAR assertions above would FAIL under a silent replacement (|p| grows by 83 % over
    the path under HAR; it is constant here)."""
    kw = C.ba06_kw(p_min=0.0)
    kw.update(rho_bar=0.8)
    P, st0 = hf_state(oracle_name, kw, -100.0 * I3, PI_FAR, v0=1.7)
    sts = run_deps(oracle_name, P, st0, shear_increments(ES_STEPS))
    for t, s in zip(ES_STEPS, sts):
        w = np.linalg.eigvalsh(s.sigma)
        p, q = float(w.sum()) / 3.0, C.SQ32 * float(np.linalg.norm(w - w.mean()))
        assert abs(p - (-100.0)) <= 1e-10 * 100.0, (t, p)
        assert abs(q - 3.0 * 5400.0 * t) <= 1e-10 * q


def _tims_sigma(p, eta, theta=None):
    q = eta * abs(p)
    nh = DEV_TXC if theta is None else C.dir_for_theta(theta)
    return np.diag(C.sigma_pq(p, q, nh))


def test_k1_11_inverse_map_printed_spot_values(oracle_name):
    """(S.5h'') / 13.11 printed spot value for the ring envelope: p = -3.5 kPa, eta = 2.1 -> varpi = 5.7714886,
    eps_v^e = 9.0504735e-4, eps_s^e = 1.2561906e-4.  Through initial_state (the oracle's own inverse map): the elastic strain
    it builds has these invariants (eps_v to 1e-7, eps_s to 5e-7 relative: printed 8 digits) and reproduces the stress to
    1e-12.  The state is above the floor (|p| = 3.5 > 0.505), so no floor event.  Plus the forward map closure: (S.5h) of the
    returned strain gives back (p, q) to 1e-12 (the inverse is the inverse of the sheet's forward map).
    KILLS: the BA06 inverse (eps_v = -kappa ln(p/p0)) used under HAR (orders of magnitude off); (S.5h'') with the exponent
    of (|p|/varpi) dropped; eps_s computed without the (varpi/p_a)^n factor (n = 0 stand-in)."""
    kw = C.har_kw(p_min="default")
    sig = _tims_sigma(-3.5, 2.1)
    P, st = hf_state(oracle_name, kw, sig, PI_FAR)
    ev, es = eps_vs(st)
    assert abs(ev - 9.0504735e-4) <= 1e-7 * 9.0504735e-4, ev
    assert abs(es - 1.2561906e-4) <= 5e-7 * 1.2561906e-4, es
    ev_f, es_f, vp = C.har_inverse(-3.5, 2.1 * 3.5)
    assert abs(vp - 5.7714886) <= 1e-7 * 5.7714886
    assert abs(ev - ev_f) <= 1e-12 and abs(es - es_f) <= 1e-12
    assert np.linalg.norm(st.sigma - sig) <= 1e-12 * 3.5
    p_fwd, q_fwd = C.har_pq(ev, es)
    assert abs(p_fwd + 3.5) <= 1e-12 * 3.5 and abs(q_fwd - 2.1 * 3.5) <= 1e-12 * 7.35
    if oracle_name == "O2":
        assert st.n_f_init == 0


def test_k1_11_inverse_map_round_trip_seeded(oracle_name):
    """Seeded round trip over 12 (p, eta, theta) states, p in [-300, -0.8] (above the floor), eta in [0, 12.9], theta in
    [0.05, pi/3 - 0.05]: initial_state's strain satisfies the forward map (S.5h) to 1e-11 relative (the closed form needs no
    Newton), its stress equals the input, and eps^e is co-axial with sigma (the principal strains are the stress's)."""
    rng = np.random.default_rng(144)
    for _ in range(12):
        p = -float(np.exp(rng.uniform(math.log(0.8), math.log(300.0))))
        eta = float(rng.uniform(0.0, 12.9))
        th = float(rng.uniform(0.05, math.pi / 3 - 0.05))
        sig = _tims_sigma(p, eta, th)
        P, st = hf_state(oracle_name, C.har_kw(p_min="default"), sig, PI_FAR)
        ev, es = eps_vs(st)
        p_f, q_f = C.har_pq(ev, es)
        assert abs(p_f - p) <= 1e-11 * abs(p), (p, p_f)
        assert abs(q_f - eta * abs(p)) <= 1e-11 * max(abs(p) * eta, abs(p) * 1e-3), (p, eta, q_f)
        assert np.linalg.norm(st.sigma - sig) <= 1e-11 * abs(p) + 1e-13
        assert np.linalg.norm(st.eps_e @ st.sigma - st.sigma @ st.eps_e) <= 1e-12 * np.linalg.norm(st.eps_e) * np.linalg.norm(st.sigma)


# ==============================================================================================
# K1.2h  closed non-coaxial elastic loop under HAR (sheet 13.2)
# ==============================================================================================
def _sym(a):
    return np.array(a, float)


SCALE_LOOP = 0.2        # the K1.2 loop (amplitude 2e-3) x 0.2: |eps_ij| <= 4e-4 << the HAR domain scale 1.06e-3
LOOP = [SCALE_LOOP * m for m in (
    np.zeros((3, 3)),
    _sym([[0.0020, 0.0015, 0.0], [0.0015, -0.0010, 0.0010], [0.0, 0.0010, 0.0005]]),
    _sym([[-0.0010, -0.0020, 0.0015], [-0.0020, 0.0015, 0.0], [0.0015, 0.0, -0.0015]]),
    _sym([[0.0005, 0.0, -0.0020], [0.0, 0.0020, 0.0015], [-0.0020, 0.0015, -0.0005]]),
    np.zeros((3, 3)))]


def test_k1_2h_closed_nonconaxial_loop_zero_work_and_state_returns(oracle_name):
    """13.2 under HAR: W = oint sigma:d eps = 0 to 1e-12 of oint |sigma:d eps| (sigma = dPsi/d eps^e, (S.4h)), the state returns
    (sigma, eps^e to 1e-12; pi_i exactly).  Start: isotropic p = -400 (eps_v = -1.05e-3: room to the edge).  Quadrature: along
    each segment sigma(E_j + t D_j) is obtained by one oracle step of t D_j from the segment start (elastic, path independent);
    sigma:D_j is sqrt(quadratic) x linear in t for n = 1/2 with the nearest complex singularity at distance
    >= 1.7 segment-lengths-units from [0,1] (eps* >= 1.1e-3 against a strain variation <= 1.2e-3 x 0.2 ... ), so 16-point
    Gauss-Legendre is exact to ~1e-18.  NON-VACUITY control: the open first segment stores W_1 = Psi(E_1) - Psi(0) with Psi of
    (S.4h) at the segment ends, to 1e-10 relative (the HAR energy is the potential of the HAR stress).
    KILLS: a non-conservative HAR tangent / stress pair (D12 sign or the t2 term: sigma != dPsi/d eps), a state that does not
    return (hidden plastic or floor work), a Psi-stress inconsistency."""
    kw = C.har_kw(p_min="default")
    P, st0 = hf_state(oracle_name, kw, -400.0 * I3, PI_FAR)
    ev0, es0 = eps_vs(st0)
    xg, wg = np.polynomial.legendre.leggauss(16)
    tg, wg = 0.5 * (xg + 1.0), 0.5 * wg
    W = Wabs = 0.0
    st = st0
    W1 = None
    for j in range(len(LOOP) - 1):
        D = LOOP[j + 1] - LOOP[j]
        w_seg = 0.0
        for t, w in zip(tg, wg):
            s = run_deps(oracle_name, P, st, [t * D])[-1]
            assert C.ok_state(oracle_name, s) and not is_plastic(s), "loop left the elastic domain"
            d = float(np.sum(s.sigma * D))
            W += w * d
            w_seg += w * d
            Wabs += w * abs(d)
        if j == 0:
            W1 = w_seg
        st = run_deps(oracle_name, P, st, [D])[-1]
        assert C.ok_state(oracle_name, st) and not is_plastic(st)
    assert Wabs > 0.0
    assert abs(W) <= 1e-12 * Wabs, f"net work {W:.3e}, sum|.| {Wabs:.3e}"
    assert np.linalg.norm(st.sigma - st0.sigma) <= 1e-12 * 400.0
    assert np.linalg.norm(st.eps_e - st0.eps_e) <= 1e-12
    assert st.pi_i == st0.pi_i
    # control: W_1 = Psi(eps^e_0 + E_1) - Psi(eps^e_0)
    e1 = st0.eps_e + LOOP[1]
    ev1 = float(np.trace(e1))
    es1 = C.SQ23 * float(np.linalg.norm(e1 - ev1 / 3.0 * I3))
    dPsi = C.har_psi(ev1, es1) - C.har_psi(ev0, es0)
    assert abs(W1 - dPsi) <= 1e-10 * abs(dPsi), (W1, dPsi)


# ==============================================================================================
# H.hess  the elastic tangent against the sheet's gate table (2.3), (S.3) in full
# ==============================================================================================
GATE_ETA = [0.0, 0.5, 1.0, 1.331, 1.71, 2.1, 3.0, 6.0, 12.87]
GATE_LMIN6 = ["0.8551", "0.8752", "0.9284", "0.9750", "1.034", "1.0980", "1.2460", "1.6837", "2.4332"]   # sheet 2.3 table, row 2


def _printed_tol(s):
    """Half a unit in the last printed decimal (+1e-9)."""
    dec = len(s.split(".")[1]) if "." in s else 0
    return 0.5 * 10.0 ** (-dec) * 1.0001 + 1e-9


@pytest.mark.parametrize("p", [-10.0, -3.5])
@pytest.mark.parametrize("eta,printed", list(zip(GATE_ETA, GATE_LMIN6)))
def test_har_elastic_tangent_matches_the_sheet_gate_table_and_eigenstructure(oracle_name, eta, printed, p):
    """Sheet 2.3 'Gate table': lambda_min of the 6-D Mandel Hessian over K_iso(p) = k p_a (|p|/p_a)^n for the 9 stress ratios
    of the table (printed: 0.8551 0.8752 0.9284 0.9750 1.034 1.0980 1.2460 1.6837 2.4332), a function of eta ALONE (p = -10 and
    -3.5 give the same row).  Also (S.5c) and the structure line of 2.3: the 6 eigenvalues are {eig of [[3 D11, sqrt2 D12],
    [sqrt2 D12, 2 D22/3]]} u {2q/(3 eps_s) = 2 g p_a (varpi/p_a)^n x 4}, with D11, D12, D22 from (S.5h') coded from the
    sheet, and the product of the two non-shear eigenvalues is 2 det D = 6 k g p_a^2 (varpi/p_a)^(2n).  The oracle supplies
    only its public elastic tangent at the state (initial_state: elastic, a^e); direction theta = 0.5 (distinct principal
    strains: no repeated-eigenvalue limit).
    KILLS: HAR replaced by BA06 (the eigenvalues are 2 mu0 / 3 kappa-type numbers, ~1e3 times off); the t4 term of (S.3)
    dropped or AB06 eq 64 used (the n_hat n_hat block moves: the table row at eta >= 1 misses by 1e-2 .. 1e-1); D12 = 0; a
    G(p) taken at q = 0 instead of at the state."""
    sig = _tims_sigma(p, eta, theta=0.5)
    P, st = hf_state(oracle_name, C.har_kw(p_min="default"), sig, PI_FAR)
    M6 = C.mandel6(ORACLES[oracle_name].tangent(P, st))
    assert np.abs(M6 - M6.T).max() <= 1e-9 * np.abs(M6).max(), "an elastic HAR tangent must be symmetric"
    eig = np.linalg.eigvalsh(0.5 * (M6 + M6.T))
    Kiso = C.har_Kiso(p)
    assert abs(eig[0] / Kiso - float(printed)) <= _printed_tol(printed) / float(printed) * float(printed), (eig[0] / Kiso, printed)
    q = eta * abs(p)
    D11, D12, D22, _ = C.har_D(p, q)
    kn = C.K_T * (1.0 - C.N_T)
    vp = math.sqrt(p * p + kn * q * q / (3.0 * C.G_T))
    two_G = 2.0 * C.G_T * C.PA_T * (vp / C.PA_T) ** C.N_T
    blk = np.array([[3.0 * D11, math.sqrt(2.0) * D12], [math.sqrt(2.0) * D12, 2.0 * D22 / 3.0]])
    want = np.sort(np.concatenate([np.linalg.eigvalsh(blk), [two_G] * 4]))
    assert np.abs(eig - want).max() <= 1e-9 * want.max(), (eig, want)
    det_blk = np.linalg.det(blk)
    assert abs(det_blk - 6.0 * C.K_T * C.G_T * C.PA_T ** 2 * (vp / C.PA_T) ** (2 * C.N_T)) <= 1e-10 * det_blk     # (S.5c)
    assert np.all(eig > 0.0), "positive definiteness (S.5c) at every eta"


# ==============================================================================================
# H.jac / H.cto / H.fin / H.chain : the FD re-runs sheet 2.4 lists, under HAR
# ==============================================================================================
S_HAR = 1.0e-3          # the HAR short strain length: eps* = 1/(k(1-n)) (|p|/p_a)^(1/2) ~ 1.05e-3 at p = -100


def fd_bound(h, delta=1.0e-12):
    """Truncation + noise bound of a central difference of the stress map: 10 (h^2/(6 s^2) + delta/h), s = S_HAR (the
    analogue of test_g1_finite_strain_k2.fd_bound with kappa_hat -> the HAR strain length)."""
    return 10.0 * (h * h / (6.0 * S_HAR ** 2) + delta / h)


# HAR preload / increment of the CTO and Jacobian gates: the BA06 pair (CT.E_PRE, CT.DEPS_FD) drives theta to 1.03 (0.016 from the
# compression corner) under the stiffer HAR law (measured with the oracle, scan on Esmeralda 2026-10-03), so the gates' own premise
# 'theta off both corners' fails; this pair keeps theta = 0.564 after the preload and 0.793 after the gated increment
E_PRE_H = np.diag([4.0e-4, 0.0, 0.0]) + CT.SHEAR_PRE
DEPS_H = 0.3 * CT.DEPS_FD


def har_pre_state(csl="fork", e_pre=E_PRE_H, npre=CT.NPRE, **over):
    kw = C.har_kw(p_min="default", **over)
    if csl == "paper":
        kw.update(csl_mode="paper", lambda_tilde=0.0135, v_c0=1.81)
    P = C.make_hf("O2", **kw)
    st0 = O2.initial_state(P, -100.0 * I3, C.v_for_psi(kw, -80.0, -0.05), -80.0)
    sts = O2.run_path(P, st0, np.array([e_pre] * npre))
    assert not any(s.flags["refused"] for s in sts), next(s.flags["reason"] for s in sts if s.flags["refused"])
    return P, sts[-1], kw


@pytest.mark.parametrize("csl", ["fork", "paper"])
def test_har_jacobian_S30_matches_central_fd_of_the_residual(csl):
    """Sheet 2.4 FD list: the 4 x 4 Jacobian (S.30) against the central FD of the residual (S.29) r(x), x = (eps^e, dlam), the
    nested pi_i solve inside, with the HAR law (D12 != 0, D22 != q/eps_s make the t2, t3/t4 terms of (S.3) live).  Evaluated at
    the converged iterate AND at a perturbed one (eps^e += 3e-6 on every component, dlam x 1.1: still a valid evaluation).
    Rows 0-2 (strain units) and row 3 (kPa) are scored separately, per-column relative to the row-block maximum; gate = best h of
    {1e-6, 1e-7, 1e-8} <= 1e-6 (the O(h^2) truncation (h/s)^2/6 with s = 1e-3 is 1.7e-7 at 1e-6, 1.7e-9 at 1e-7).
    KILLS: a (S.30) entry built from the BA06 D (D12 = 0, q/eps_s = D22); the Pi_b = P_a a^e product with a symmetrised a^e;
    q_ab / q_api / f_a wrong sign (found at G0 for BA06: still must hold for HAR)."""
    P, st, kw = har_pre_state(csl)
    deps = DEPS_H
    w, V = np.linalg.eigh(st.eps_e + 0.5 * (deps + deps.T))
    v = st.v * math.exp(float(np.trace(deps)))
    res = KER.return_map(P, w, st.pi_i, v, v)
    assert res.plastic and not res.refused
    for label, x in (("converged", np.concatenate([res.eps_e, [res.dlam]])),
                     ("perturbed", np.concatenate([res.eps_e + 3e-6, [1.1 * res.dlam]]))):
        def r_of(xx):
            return KER.evaluate(P, xx[:3], xx[3], w, v, st.pi_i).r
        pe = KER.evaluate(P, x[:3], x[3], w, v, st.pi_i)
        J = KER.jacobian(P, pe)
        best = None
        for h in (1e-6, 1e-7, 1e-8):
            fd = np.zeros((4, 4))
            for j in range(4):
                e = np.zeros(4)
                e[j] = h * (max(abs(x[3]), 1e-6) if j == 3 else 1.0)
                fd[:, j] = (r_of(x + e) - r_of(x - e)) / (2.0 * e[j])
            err = 0.0
            for rows in (slice(0, 3), slice(3, 4)):
                sc = np.abs(J[rows, :]).max()
                err = max(err, np.abs(J[rows, :] - fd[rows, :]).max() / sc)
            best = err if best is None else min(best, err)
        assert best <= 1e-6, f"{csl}/{label}: Jacobian (S.30) vs FD best-h error {best:.3e}"


@pytest.mark.parametrize("csl", ["fork", "paper"])
def test_har_cto_S33_matches_central_fd_off_the_corners_all_six_columns(csl):
    """Sheet 2.4: (S.33) non-coaxial CTO off the corners under HAR vs the central FD of O2's own stress update (the G1.fd protocol
    of test_g1_convergence_tangents, six independent columns, three shears).  State: 10 non-coaxial plastic preload steps
    (E_PRE, three distinct principal strains, theta mid-range asserted); increment DEPS_FD continues the loading.  Gate: best h
    of {1e-6, 1e-7, 1e-8} <= 1e-6 per column, shear columns named (truncation (h/s)^2/6, s = 1e-3: 1.7e-7 at 1e-6).  Negative
    control: the ELASTIC tangent a^e at the same state misses the plastic CTO by >= 1e-2 relative (the gate is not vacuous).
    KILLS: a^ep built with the BA06 D; (S.3) t4 / t2 terms dropped; the spin g_ab = (sigma_a - sigma_b)/(eps~_a - eps~_b) wrong."""
    P, st, kw = har_pre_state(csl)
    stn = O2.step(P, st, DEPS_H)
    th = CT.theta_of(stn.sigma)
    assert 0.1 < th < math.pi / 3 - 0.1 and stn.flags["plastic"] and stn.flags["substeps"] == 1
    best = {}
    for h in (1e-6, 1e-7, 1e-8):
        _, _, errs = CT._fd_and_cto(P, st, DEPS_H, h)
        best[h] = (max(errs.values()), max(errs[kl] for kl in ((0, 1), (0, 2), (1, 2))))
    hb = min(best, key=lambda h: best[h][0])
    print(f"\n[HAR {csl}] CTO vs FD best-h {hb}: max col {best[hb][0]:.3e}, shear {best[hb][1]:.3e}; all {best}")
    assert best[hb][0] <= 1e-6 and best[hb][1] <= 1e-6, best
    Ct = O2.tangent(P, stn)
    ae = O2.tangent(P, O2.step(P, st, 1e-12 * DEPS_H))        # a sub-F_tol increment: the elastic branch, a^e
    assert not O2.step(P, st, 1e-12 * DEPS_H).flags["plastic"]
    assert np.linalg.norm(Ct - ae) / np.linalg.norm(Ct) >= 1e-2, "elastic and plastic tangents indistinguishable: the FD gate has no power"


F1H = 0.25 * FS.F1_LOG          # the (S.43) f1 increment scaled for the stiffer HAR law (3 steps -> a plastic state)


@functools.lru_cache(maxsize=None)
def _fin_state(csl):
    kw = C.har_kw(p_min="default")
    if csl == "paper":
        kw.update(csl_mode="paper", lambda_tilde=0.0135, v_c0=1.81)
    P = C.make_hf("O2", **kw)
    st = O2.initial_state(P, -100.0 * I3, C.v_for_psi(kw, -80.0, -0.05), -80.0, finite=True)
    prev = None
    for _ in range(3):
        prev = st
        st = O2.step(P, st, np.diag(F1H))
        assert not st.flags["refused"], st.flags["reason"]
    return P, prev, st


@pytest.mark.parametrize("h", [1e-6, 1e-7])
@pytest.mark.parametrize("csl", ["fork", "paper"])
def test_har_finite_tangent_S34_matches_fd_of_the_nominal_stress(csl, h):
    """Sheet 2.4 / 9.5: (S.34) under HAR against the central FD of P(E) = tau((1 + hE) f0)(1 + hE)^-T over the nine unit E_kl
    and a dense non-symmetric g (the G1 gate FS.fd_o2 on the HAR law; the log-strain protocol of the (S.43) f1 steps, 3 steps of
    0.25 F1 from p = -100).  Bound fd_bound(h) = 10 (h^2/(6 s^2) + 1e-12/h), s = 1e-3.  NEGATIVE CONTROLS from the same blocks (must
    miss the bound by >= 10x): the half-spin variant (BA06 3.48) and the missing tau(+)1 term.
    KILLS: the half-spin sum; tau(+)1 dropped; gamma~_ab wrong; a v-term of (S.31) with v0 or v_n under HAR."""
    P, prev, st = _fin_state(csl)
    assert st.flags["plastic"]
    f0 = np.diag(np.exp(F1H))
    w_prev, V_prev = np.linalg.eigh(prev.eps_e)
    b_n = (V_prev * np.exp(2.0 * w_prev)) @ V_prev.T

    def tau_of(f):
        b_tr = f @ b_n @ f.T
        wb, Vb = np.linalg.eigh(0.5 * (b_tr + b_tr.T))
        vv = prev.v * np.linalg.det(f)
        res = KER.return_map(P, 0.5 * np.log(wb), prev.pi_i, vv, vv)
        assert res.plastic and not res.refused, f"FD map left the plastic branch: {res.reason}"
        return (Vb * res.sig) @ Vb.T

    assert np.abs(tau_of(f0) - st.sigma).max() <= 1e-10 * 100.0, "FD map is not the O2 step map"
    V, ct, gam, tau = FS._blocks_from_o2_cache(st)
    a_o2 = O2.tangent_finite(P, st)
    assert np.abs(FS._assemble(V, ct, gam, tau) - a_o2).max() <= 1e-9 * np.abs(a_o2).max()
    bound = fd_bound(h)
    _, g = FS._fd_directions()
    e_o2 = FS._fd_errors(a_o2, tau_of, f0, h, g)
    assert max(e_o2) <= bound, f"HAR {csl} h={h}: (S.34) FD error unit {e_o2[0]:.3e}, dense {e_o2[1]:.3e} > bound {bound:.3e}"
    for name, kw in (("half-spin", dict(spin=0.5)), ("no-tau(+)1", dict(tau_plus_one=False))):
        e_v = FS._fd_errors(FS._assemble(V, ct, gam, tau, **kw), tau_of, f0, h, g)
        assert max(e_v) >= 10.0 * bound, f"negative control '{name}' error {max(e_v):.3e} < 10 x {bound:.3e}: the gate does not discriminate"


# ---- chained tangent (S.45)-(S.47) under HAR ------------------------------------------------------
def preloaded_har(csl, e_pre, npre):
    kw = C.har_kw(p_min="default")
    if csl == "paper":
        kw.update(csl_mode="paper", lambda_tilde=0.0135, v_c0=1.81)
    P = C.make_hf("O2", **kw)
    st0 = O2.initial_state(P, ST.SIG0, C.v_for_psi(kw, -80.0, -0.05), -80.0)
    sts = O2.run_path(P, st0, np.array([e_pre] * npre))
    assert not any(s.flags["refused"] for s in sts), next(s.flags["reason"] for s in sts if s.flags["refused"])
    st = sts[-1]
    off = np.linalg.norm(st.sigma - np.diag(np.diag(st.sigma)))
    assert off > 1e-3 * np.linalg.norm(st.sigma), "preloaded state is coaxial"
    return P, st


# forced-subdivision scales for the stiffer HAR law: the whole-increment return map is refused beyond ~8 D (probe: m = 2, 4, 8 at
# s = 16, 32, 64), so these three scales give ladder levels m = 2, 4, 8 (premise asserted by gate_increment: substeps >= 2)
SCALES_H = (16, 32, 64)
HAR_GENERIC = {
    "har_fork_A": ("fork", ST.E_PRE_A, 10, ST.D_A),
    "har_fork_B": ("fork", ST.E_PRE_B, 12, ST.D_B),
    "har_fork_C": ("fork", ST.E_PRE_C, 10, ST.D_C),
    "har_paper_A": ("paper", ST.E_PRE_A, 10, ST.D_A),
    "har_paper_C": ("paper", ST.E_PRE_C, 10, ST.D_C),
}


def _har_rows(names, scales):
    rows = []
    for name in names:
        csl, e_pre, npre, d = HAR_GENERIC[name]
        P, st = preloaded_har(csl, e_pre, npre)
        for s in scales:
            rows.append(ST.gate_increment(P, st, s * d, f"{name} s={s}"))
    return rows


def test_har_chained_tangent_matches_fd_on_substepped_increments_default_tier():
    """Sheet 2.4 / 9.6: the chained tangent (S.45)-(S.47) under HAR on forced-substep increments, one per (state, scale): fork
    A/B/C at s = 32 (m = 4), paper A at s = 16 (m = 2), paper C at s = 64 (m = 8).  Gate (G1 substep protocol): best h of
    {1e-6, 1e-7}: max column error <= 1e-6, shear columns included; every FD point plastic, accepted, with the same substep count
    and elastic/plastic pattern; the base increment off the theta corners and not on the vertex.
    KILLS: the last-sub-increment tangent (3.6e-2 .. 0.9 off); a chain that drops the pi_i / v columns; any (S.45) entry built from
    the BA06 D (D12 = 0) -- the chain consumes u = b t, w, kappa that depend on a^e through (S.3) (t2, t4 live under HAR)."""
    rows = []
    for name, s in (("har_fork_A", 32), ("har_fork_B", 32), ("har_fork_C", 32), ("har_paper_A", 16), ("har_paper_C", 64)):
        csl, e_pre, npre, d = HAR_GENERIC[name]
        P, st = preloaded_har(csl, e_pre, npre)
        rows.append(ST.gate_increment(P, st, s * d, f"{name} s={s}"))
    ST.report(rows, "H.chain default tier")
    assert len(rows) == 5
    ST.assert_rows(rows, "H.chain", min_valid_frac=1.0)


@pytest.mark.slow
def test_har_chained_tangent_matches_fd_full_scale_grid():
    """H.chain, the full grid: 5 HAR states (fork A/B/C, paper A/C) x scales {16, 32, 64} (ladder levels m = 2, 4, 8): 15
    increments, all FD-valid and within the gate.  Wall time ~ 8 min (the FD evaluates a 2^k ladder at 24 points per increment)."""
    rows = _har_rows(HAR_GENERIC, SCALES_H)
    ST.report(rows, "H.chain full grid")
    assert len(rows) == 15
    ST.assert_rows(rows, "H.chain full", min_valid_frac=1.0)


def test_har_chained_tangent_non_uniform_fractions():
    """Sheet 9.6 (E) under HAR: recursive-halving fractions through step_fractions (no ladder), fork A at s = 16 with
    (1/2, 1/4, 1/4) and fork B at s = 32 with (1/4, 1/4, 1/4, 1/8, 1/8) (largest sub-increment 8 D, inside the basin of a
    single backward-Euler step under HAR: probe m = 1 at s <= 8).  Same gate and protocol as the substep-tangent file.
    KILLS: a chain that assumes alpha_k = 1/m (T_k = S + alpha_k E_J, cumulative fraction in S^v)."""
    rows = []
    for name, s, frs in (("har_fork_A", 16, (0.5, 0.25, 0.25)), ("har_fork_B", 32, (0.25, 0.25, 0.25, 0.125, 0.125))):
        csl, e_pre, npre, d = HAR_GENERIC[name]
        P, st = preloaded_har(csl, e_pre, npre)
        rows.append(ST.gate_increment(P, st, s * d, f"{name} s={s} {frs}",
                                      stepper=lambda P_, A_, d_, frs=frs: O2.step_fractions(P_, A_, d_, frs)))
    ST.report(rows, "H.chain fractions")
    ST.assert_rows(rows, "H.chain fractions", min_valid_frac=1.0)


# ==============================================================================================
# H.conv  O2 -> O1 under HAR (plan 5.1)
# ==============================================================================================
NS_H = (25, 50, 100, 200, 400)
HAR_PATHS = ("TXC_drained", "TXE_drained", "NONCOAXIAL")


def _har_start(oname, spec, **over):
    kw = C.har_kw(p_min="default", **over)
    P = C.make_hf(oname, **kw)
    v0 = C.v_for_psi(kw, spec["pi0"], spec["psi0"])
    return P, ORACLES[oname].initial_state(P, -100.0 * I3, v0, spec["pi0"])


def _har_run(oname, path, n):
    spec = CT.PATHS[path]
    P, st0 = _har_start(oname, spec)
    ora = ORACLES[oname]
    if spec["kind"] == "tx":
        sts = ora.triaxial(P, st0, spec["tx"], spec["ax"], n, **({"rtol": 1e-10} if oname == "O1" else {}))
    else:
        d = CT.noncoax_deps(n)
        sts = ora.run_path(P, st0, d, rtol=1e-10) if oname == "O1" else ora.run_path(P, st0, d)
    return P, st0, sts


@functools.lru_cache(maxsize=None)
def _har_truth(path):
    return _har_run("O1", path, CT.PATHS[path]["n_truth"])


@pytest.mark.parametrize("path", HAR_PATHS)
def test_har_o2_converges_to_o1_first_order(path):
    """Plan 5.1 under HAR: backward Euler (O2) converges to the rate integrator (O1) at first order.  Stated argument: err(n) ~ C/n,
    so err strictly decreases over n = 25 .. 400, the log-log slope is in [0.8, 1.3], and err(400) <= err(25)/8 (first order
    predicts 1/16; the G1 file's a-priori size bound uses the BA06 strain scale kappa_hat and is replaced here by this ratio,
    since the HAR short strain length 1/(k(1-n)) = 1.06e-3 is 10x shorter than kappa_hat).  Quantities: Cauchy-axis stress (axial,
    lateral for triaxial; the full tensor for NONCOAXIAL) and pi_i.  The O1 truth must complete.
    KILLS: an O1 or O2 HAR law that is not the same (D12, q/eps_s, p_a): the gap does not close; a non-first-order scheme."""
    P1, st01, sts1 = _har_truth(path)
    n_truth = CT.PATHS[path]["n_truth"]
    assert len(sts1) == n_truth and all(C.ok_state("O1", s) for s in sts1), "O1 truth did not complete"
    sig1, pi1, _ = CT.endpoint("O1", path, sts1[-1])
    e_sig, e_pi = [], []
    for n in NS_H:
        P2, st02, sts2 = _har_run("O2", path, n)
        assert len(sts2) == n and not any(s.flags["refused"] for s in sts2), f"O2 n={n} refused"
        sig2, pi2, _ = CT.endpoint("O2", path, sts2[-1])
        e_sig.append(np.linalg.norm(sig2 - sig1) / np.linalg.norm(sig1))
        e_pi.append(abs(pi2 - pi1) / abs(pi1))
    print(f"\n[HAR {path}] n = {NS_H}\n" + CT.table(["sigma", "pi_i"], [e_sig, e_pi]))
    for name, errs in (("sigma", e_sig), ("pi_i", e_pi)):
        assert all(e > 0.0 for e in errs)
        assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)), f"{name}: not strictly decreasing {errs}"
        p = CT.observed_order(NS_H, errs)
        assert 0.8 <= p <= 1.3, f"{name}: observed order {p:.3f}; errs {errs}"
        assert errs[-1] <= errs[0] / 8.0, f"{name}: err(400)/err(25) = {errs[-1] / errs[0]:.3f} > 1/8"


# ==============================================================================================
# H.diss  dissipation census under HAR (K1.9, S.38 / S.39)
# ==============================================================================================
@pytest.mark.parametrize("path", list(CT.PATHS))
def test_har_dissipation_census_every_step(oracle_name, path):
    """K1.9 under HAR (the dissipation proof never uses the energy: sheet 2.3 item (8), 11.4): D >= 0 at every plastic step and D = 0
    at every elastic one, both oracles, the four convergence paths (O2 at n = 100, O1 at its truth increments), condition A holding
    (N_bar 0.2 <= N 0.4; rho/rho_bar = 1 >= 0.75).  O2 exact; O1 the integration tolerance (G1: -1e-9 |p0| eps_max).
    KILLS: a HAR stress that is not the derivative of the HAR Psi feeding a negative D through the plastic update; a floor leaking
    work into D (the floor is inert on these paths: asserted by the zero counters)."""
    if oracle_name == "O1":
        P, st0, sts = _har_truth_any(path)
        n = CT.PATHS[path]["n_truth"]
        emax = abs(CT.PATHS[path]["ax"]) / n if CT.PATHS[path]["kind"] == "tx" else np.abs(CT.noncoax_deps(n)).max()
        tolD = -1e-9 * 100.0 * emax
    else:
        P, st0, sts = _har_run("O2", path, 100)
        n, tolD = 100, 0.0
    assert len(sts) == n and all(C.ok_state(oracle_name, s) for s in sts), "path did not complete"
    D = np.array([s.D for s in sts])
    plastic = np.array([is_plastic(s) for s in sts])
    assert plastic.any(), "the path never yields: not a dissipation test"
    assert np.all(D >= tolD), f"negative dissipation under HAR: min D = {D.min():.3e}"
    assert np.all(D[~plastic] == 0.0), "D must be exactly 0 on elastic steps"
    if oracle_name == "O2":
        assert all((s.n_f_tr, s.n_f_post) == (0, 0) for s in sts), "floor event on a p ~ 100 kPa path"


def _har_truth_any(path):
    return _har_truth(path)


def test_har_dissipation_census_random_increments_o2():
    """D >= 0 on 200 seeded random increments (mixed normal + three shears, |d eps| <= 6e-4 per component) from the dense HAR start,
    O2; D = 0 on elastic steps; refused steps are skipped (counted: a census needs >= 150 accepted plastic steps).
    KILLS: as the census above, on strain directions no triaxial path visits (all three shears, extension, compression)."""
    rng = np.random.default_rng(1441)
    spec = CT.PATHS["TXC_drained"]
    P, st = _har_start("O2", spec)
    n_pl = n_acc = 0
    for _ in range(200):
        a = rng.uniform(-6e-4, 6e-4, size=(3, 3))
        d = 0.5 * (a + a.T) - 3e-5 * I3          # slight compressive bias keeps p away from the edge
        nxt = O2.step(P, st, d)
        if nxt.flags["refused"]:
            continue
        n_acc += 1
        if is_plastic(nxt):
            n_pl += 1
            assert nxt.D >= 0.0, f"negative dissipation {nxt.D:.3e}"
        else:
            assert nxt.D == 0.0
        st = nxt
    assert n_acc >= 150 and n_pl >= 50, (n_acc, n_pl)


# ==============================================================================================
# H.ref  parser refusals (sheet 2.4) and defaults (1.3 / 9.7), both oracles
# ==============================================================================================
def _validated(oname, **kw):
    P = C.make_hf(oname, **kw)
    return P.validate()


@pytest.mark.parametrize("field,val", [("p0", -100.0), ("kappa_hat", 0.01), ("eps_v0", 0.0), ("mu0", 5400.0), ("alpha0", 0.0)])
def test_har_refuses_every_ba06_energy_parameter_even_at_its_default_value(oracle_name, field, val):
    """Sheet 2.4: any of -p0 -kappa_hat -eps_v0 -mu0 -alpha0 given together with -energy HAR is REFUSED, never ignored (HAR replaces
    alpha0 and the other four are not read); a GIVEN value equal to the BA06 default (0.0 for alpha0, eps_v0) is still given.
    KILLS: a parser that silently ignores the BA06 parameters under HAR (the user thinks alpha0 = 5 is active)."""
    with pytest.raises(ValueError):
        _validated(oracle_name, **C.har_kw(p_min=0.0, **{field: val}))
    _validated(oracle_name, **C.har_kw(p_min=0.0))            # positive control: the same set without the field


@pytest.mark.parametrize("field,val", [("k", 1000.0), ("g", 500.0), ("n_e", 0.5)])
def test_ba06_refuses_every_har_parameter(oracle_name, field, val):
    """Sheet 2.4: -k -g -n given under -energy BA06 (the default) are refused.  `-p_a` is NOT in this list: it is the fork-CSL
    parameter too (positive controls below).
    KILLS: HAR parameters silently dropped on a BA06 deck (the user believes the energy is HAR)."""
    with pytest.raises(ValueError):
        _validated(oracle_name, **C.ba06_kw(p_min=0.0, **{field: val}))


def test_p_a_is_one_flag_accepted_under_ba06_and_har_and_both_csl(oracle_name):
    """Round 3b A4: p_a is ONE flag shared by the HAR energy and the fork CSL.  Accepted: BA06 + paper CSL + p_a (only echoed);
    BA06 + fork CSL + p_a; HAR + fork CSL + p_a; HAR + paper CSL + p_a.  KILLS: p_a placed in the BA06 refusal list (a BA06 fork deck
    would be refused); a separate HAR p_a and CSL p_a."""
    _validated(oracle_name, **C.ba06_kw(p_min=0.0, p_a=101.0))
    _validated(oracle_name, **C.ba06_kw(p_min=0.0, csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45, p_a=101.0))
    _validated(oracle_name, **C.har_kw(p_min=0.0))
    _validated(oracle_name, **C.har_kw(p_min=0.0, csl_mode="paper", lambda_tilde=0.0135, v_c0=1.81))


@pytest.mark.parametrize("label,over", [
    ("k=0", dict(k=0.0)), ("k<0", dict(k=-1.0)), ("g=0", dict(g=0.0)), ("g<0", dict(g=-5.0)),
    ("n<0", dict(n_e=-0.1)), ("n=1", dict(n_e=1.0)), ("n>1", dict(n_e=1.5)), ("p_a=0", dict(p_a=0.0)), ("p_a<0", dict(p_a=-101.0)),
    ("p_min<0", dict(p_min=-0.1)), ("energy typo", dict(energy="har")), ("energy unknown", dict(energy="HAR05"))])
def test_har_parameter_range_refusals(oracle_name, label, over):
    """Sheet 2.4 hard refusals: k <= 0, g <= 0, n < 0, n >= 1 (n = 1 is HAR05 eq 47-48, another closed form, not shipped), p_a <= 0,
    p_min < 0, an unknown energy name.  KILLS: a missing range check (n = 1 divides by zero at k(1-n); p_a <= 0 inverts p)."""
    kw = C.har_kw(p_min=0.0)
    kw.update(over)
    with pytest.raises(ValueError):
        _validated(oracle_name, **kw)


@pytest.mark.parametrize("missing", ["k", "g", "n_e"])
def test_har_needs_k_g_and_n(oracle_name, missing):
    """Sheet 2.4: HAR parameters are k, g, n (p_a is shared and has its own default).  A missing one is a refusal."""
    kw = C.har_kw(p_min=0.0)
    del kw[missing]
    with pytest.raises(ValueError):
        _validated(oracle_name, **kw)


def test_n_zero_is_accepted_and_n_one_is_refused(oracle_name):
    """0 <= n < 1 (sheet 2.3): n = 0 (HAR05 eq 22) is shipped, n = 1 is not."""
    _validated(oracle_name, **C.har_kw(p_min=0.0, n_e=0.0))
    with pytest.raises(ValueError):
        _validated(oracle_name, **C.har_kw(p_min=0.0, n_e=1.0))


@pytest.mark.parametrize("oname", ["O1", "O2"])
def test_default_pmin_is_five_thousandths_of_the_reference_pressure(oname):
    """Sheet 1.3 / 9.7 / 2.4 table: the p_min default is 5e-3 p_ref with p_ref = |p0| (BA06) or p_a (HAR): 0.5 kPa on the K2 set,
    0.505 kPa on the TIMs set, and it follows the reference (|p0| = 50 -> 0.25; p_a = 200 -> 1.0).
    KILLS (M-F8): a default taken from the wrong reference (|p0| under HAR: p0 is None), a unit-blind constant (0.5 kPa for every
    deck)."""
    def pm(P):
        return P.p_min if oname == "O2" else P.pmin
    assert abs(pm(C.make_hf(oname, **C.ba06_kw(p_min="default"))) - 0.5) <= 1e-15
    assert abs(pm(C.make_hf(oname, **C.har_kw(p_min="default"))) - 0.505) <= 1e-15
    assert abs(pm(C.make_hf(oname, **C.ba06_kw(p_min="default", p0=-50.0))) - 0.25) <= 1e-15
    assert abs(pm(C.make_hf(oname, **C.har_kw(p_min="default", p_a=200.0))) - 1.0) <= 1e-15
    assert abs(pm(C.make_hf(oname, **C.har_kw(p_min=0.3))) - 0.3) <= 1e-15           # an explicit value is used as given
    assert pm(C.make_hf(oname, **C.har_kw(p_min=0.0))) == 0.0                         # 0 = off


@pytest.mark.parametrize("oname", ["O1", "O2"])
def test_initial_state_with_positive_pressure_is_refused(oname):
    """Sheet 2.4: an initial state with p >= 0 after the floor (only possible with p_min = 0) is refused, in both energies."""
    for kw in (C.har_kw(p_min=0.0), C.ba06_kw(p_min=0.0)):
        P = C.make_hf(oname, **kw)
        with pytest.raises(ValueError):
            ORACLES[oname].initial_state(P, +1.0 * I3, 1.7, -50.0)
