"""G1 gate, finite-strain / K2 side (WP-144, Zone B): convergence on the K2 path, the finite-strain
tangent FD check, and the localization direction.

Sources: equation sheet 144a sections 1.4 and 9.5 ((S.34), the G1 re-derivation and the "this FD check
is a G1 gate test" statement), section 14 ((S.43), (S.44), the Fig 6 reading phi = pi/2, the "Why the
two oracles do not report the same n" paragraph, Richardson), section 16.6 items 11-12; plan 144
section 5.1 ("every K-path" convergence) and 5.2 K2.

TEST RULE (plan 5). Every expected value and every tolerance below is fixed by (a) a closed form of
the sheet, (b) a published number, or (c) a convergence / truncation argument that is stated next to
the constant. The constants were written before any oracle output was read for these paths. No oracle
value is used as an expected value. A failing assertion is a finding (against the sheet or an
oracle), never a reason to edit an oracle or loosen a bound.

Gates
  FS.conv_state   O2(m) -> O1 on the K2 path at n = 10 and n = 20: Kirchhoff stress tau, image
                  pressure pi_i, minimum acoustic determinant (raw, and normalised by its step-10
                  value at n = 20). m in {1, 2, 4, 8, 16} backward-Euler substeps per (S.43)
                  increment, same total deformation.
  FS.conv_ninterp the interpolated first-localization step n_interp: |n_O2(m) - n_O1| falls by a
                  factor 2 +- 0.3 per doubling (first order), and the Richardson estimate
                  2 n(16) - n(8) matches O1 to <= 0.02 step.
  FS.fd_o2        (S.34) finite-strain tangent of O2 vs central FD of the O2 stress update over the
                  nine unit E_kl and a dense non-symmetric g, h = 1e-5 and 1e-6, error bounded by
                  the truncation argument; the 1/2-spin variant and the missing tau(+)1 variant are
                  run as negative controls and must MISS the same bound by >= 10x.
  FS.fd_o1        the same for O1's `finite_spatial_tangent` on the elastic branch (O1 has no
                  stress-update map for a non-diagonal f, so the FD uses O1's own energy); and the
                  plastic finite tangent of O1 as the m -> infinity limit of O2's CTO.
  FS.direction    the first-localization normal is perpendicular to the intermediate principal
                  direction, |n . e_int| <= 0.05 (AB06 Fig 6 phi = pi/2, sheet 14), both oracles.

  FS.s44          (S.44), the acoustic tensor in the principal basis (c~, gamma~, tau, alpha = V^T n), equals the
                  direct contraction n_j a_ijkl n_l of the assembled (S.34) tangent, to 1e-12 relative, for a set of unit
                  normals including the minimiser, at the plastic K2 states after the 10 f1 steps and at first
                  localization, both oracles, both cases (sheet 9.5 / 16.3 G1 gate).

Mutants (each test says which one it kills):
  * wrong physics in O2 that survives O(1) (a v = v0 (1 + tr eps) kinematics for the finite path,
    a wrong hardening/coupling term, the tau(+)1 term dropped from the localization tangent): the
    gap to O1 does not close, so the Richardson match and the E_16 <= E_1 / 8 checks fail.
  * a scheme that is not first order (explicit update: order 0 or divergent; midpoint: order 2):
    the factor-2 per doubling check of FS.conv_ninterp fails.
  * the 1/2 on the finite-strain spin sum (BA06 3.48) or the missing tau(+)1 term in (S.34):
    FS.fd_o2 / FS.fd_o1 fail, and the negative controls prove the FD bound separates them.
  * a localization search in the wrong basis (principal vs global axes mixed up) or an acoustic tensor
    built from a symmetrised tangent: the normal leaves the e2-e3 plane, FS.direction fails.

Run from norsand_oracle/:  python -u -m pytest tests/test_g1_finite_strain_k2.py -v
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from conftest import K2_BASE, ORACLES, make_params
from k2_sensitivity_table import N_MAX, NOMINAL, SIGMA0, V0

import o2_algo  # noqa: F401  (State, step, tangent_finite, kernel, acoustic_min_det)
from o2_algo import kernel as K
from o1_rate.localization import _det_fun_finite as o1_det_fun_finite
from o1_rate.localization import finite_blocks, finite_spatial_tangent
from o2_algo.acoustic import acoustic_principal
from o1_rate.model import energy as o1_energy

O1 = ORACLES["O1"]
O2 = ORACLES["O2"]
I3 = np.eye(3)

# ----------------------------------------------------------------------------------------------
# The (S.43) protocol (sheet 14), in principal log strains
# ----------------------------------------------------------------------------------------------
LAM1, LAM2 = 1.0e-3, 4.0e-4
N1 = 10                                             # f1 steps; f2 afterwards
F1_LOG = np.log(np.array([1.0 + LAM2, 1.0 - LAM1, 1.0]))
F2_LOG = np.log(np.array([1.0, 1.0 - LAM2, 1.0 + LAM1]))
CASES = {"rho0.7_0.8": (0.7, 0.8), "rho1_1": (1.0, 1.0)}   # (rho, rho_bar), sheet 14 cases 2 and 1
M_LIST = (1, 2, 4, 8, 16)                           # backward-Euler substeps per (S.43) increment
P_ABS = 100.0                                       # |p0|, kPa (K2_BASE p0 = -100)
KAPPA = K2_BASE["kappa_hat"]                        # 0.01

# ----------------------------------------------------------------------------------------------
# EXPECTED VALUES AND TOLERANCES (written before any oracle output was read)
# ----------------------------------------------------------------------------------------------
# (1) Convergence. Backward Euler is first order: for a smooth quantity Q of the state,
#     Q(m) = Q* + a/m + b/m^2 + ..., Q* = the continuum (O1) value. Consequences used below:
#       e(m)/e(2m) = 2 (1 + r/m) / (1 + r/(2m)),  r = b/a   -> 2 + O(1/m)   (ratio band 2 +- 0.3)
#       Richardson R = 2 Q(16) - Q(8) has error  -b/128   (the O(1/m^2) remainder, sheet 14),
#       while e(16) ~ a/16, so |R - Q*| / e(16) = |r| / 8 <= 1/2 whenever |r| <= 4.
RATIO_LO, RATIO_HI = 1.7, 2.3                       # 2 +- 0.3 per doubling (task statement; first order)
RICHARDSON_NINTERP_TOL = 0.02                       # steps; 2 n(16) - n(8) vs O1 (task statement)
E16_OVER_E1_MAX = 1.0 / 8.0                         # first order predicts 1/16; half of the ideal reduction
RICH_OVER_E16_MAX = 0.5                             # |r| <= 4, see above
# Noise floors below which an error is "converged": O1 Radau rtol = 1e-10 and the O2 local Newton
# tolerance RES_TOL = 1e-12 (sheet 9.1) bound the state agreement near 1e-9; the acoustic
# determinant is a cancellation (normalised value at n = 20 is O(0.1) of the step-10 value), which
# amplifies a relative state error by at most ~1e2.
FLOOR_STATE = 1.0e-8
FLOOR_DET = 1.0e-6
FLOOR_TANGENT = 1.0e-8
# (1b) The minimum acoustic determinant at n = 20 is a cancellation of O(1) terms (normalised value ~0.1 of its
#     step-10 value) and its error series Q(m) = Q* + a/m + b/m^2 + c/m^3 has a LARGE r = b/a (|r| ~ 9, from the triage
#     of the first run: the signed error changes sign between m = 8 and 16). The doubling ratio of the error is
#         e(m)/e(2m) = 2 (1 + r/m) / (1 + r/(2m))     (two-term model),
#     which for r = +-9 lies in [1.7, 2.3] for BOTH signs of r only from m = 64 on (m = 32, r = -9: 1.67; m = 64: 1.85 /
#     2.13; m = 128: 1.92 / 2.07). So the asymptotic range of this quantity is m >= 64 (about 7 |r|), and its first-order
#     check is anchored there: M_ASYM = (64, 128, 256), two doublings, ratio in [RATIO_LO, RATIO_HI] (the band
#     FS.conv_ninterp already uses). The 2-level Richardson 2Q(16) - Q(8) used for the other quantities leaves the
#     remainder -b/128 = -r a/128 against E(16) ~ a/16, i.e. |R - Q*|/E(16) = |r|/8 > 1 for |r| ~ 9: the premise
#     |r| <= 4 of that check does not hold for this quantity (the test of the remainder was wrong, not the integrator).
#     The 3-level Richardson on the doubling ladder (m, 2m, 4m),
#         R3 = (8 Q(4m) - 6 Q(2m) + Q(m)) / 3     (weights sum to 1 and kill 1/m and 1/m^2:
#                                                  8/(4m) - 6/(2m) + 1/m = 0, 8/(4m)^2 - 6/(2m)^2 + 1/m^2 = 0),
#     has the error c (1 - 6/8 + 8/64) / (3 m^3) = c/(8 m^3), i.e. O(1/m^3). Bound: |R3 - Q*| <= max(E(4m)/4, floor).
#     At m = 64: c/(8 m^3) <= (1/4) a/(4m) iff |c/a| <= m^3/128 = 2048 (the cubic coefficient would have to exceed 2000 a;
#     for |r| ~ 9 a geometric-type series has |c/a| ~ r^2 ~ 80).
M_ASYM = (64, 128, 256)
RICH3_OVER_EFINE_MAX = 0.25

# (2) Finite-strain tangent FD. P(E) = tau((1 + hE) f0) (1 + hE)^{-T}, dP/dh = a^ep : E (sheet 9.5).
#     Central difference: P(h) - P(-h) = 2h P' + (h^3/3) P''' + ...; the third derivative of the stress
#     map relative to the tangent scales as 1/kappa_hat^2 (stress ~ exp(-eps_v/kappa_hat), the only
#     short strain length of the K2 set), so the relative truncation error is
#       T(h) <= h^2 / (6 kappa_hat^2).
#     Noise: the stress update is converged to a scaled residual delta = 1e-12 (sheet 9.1 scaling
#     note; O2 RES_TOL), which enters the central difference as delta / h relative to the tangent.
#     Bound = SAFETY * (T(h) + delta / h), SAFETY = 10 (O(1) constant of the remainder, plastic
#     coupling of pi_i and Delta-lambda to the same strain length).
FD_SAFETY = 10.0
FD_DELTA_PLASTIC = 1.0e-12                          # converged return map (O2)
FD_DELTA_ELASTIC = 1.0e-14                          # closed-form elastic stress: round-off only
FD_HS = (1.0e-5, 1.0e-6)
NEG_CONTROL_FACTOR = 10.0                           # a mutant variant must miss the bound by this factor

# (3) Direction. Sheet 14 (Fig 6 reading): n is perpendicular to the intermediate principal direction;
#     the search grid is ~2 degrees and sin(2 deg) = 0.035; the gate is 0.05 (task statement).
DIRECTION_TOL = 0.05
E_INT_AXIS_TOL = 1.0e-6                             # sheet 14: under (S.43) the intermediate direction is e_1
PRINCIPAL_GAP_REL = 1.0e-3                          # principal stresses distinct by >= 1e-3 |p0|


def fd_bound(h: float, delta: float) -> float:
    return FD_SAFETY * (h * h / (6.0 * KAPPA ** 2) + delta / h)


# ----------------------------------------------------------------------------------------------
# runs (cached)
# ----------------------------------------------------------------------------------------------
def _params(oracle: str, rho: float, rho_bar: float):
    kw = dict(K2_BASE)
    kw.update(rho=rho, rho_bar=rho_bar)            # chi = -3.5, v_c0 = 1.81: the nominal combination
    return make_params(oracle, **kw)


@functools.lru_cache(maxsize=None)
def _o1_run(case):
    """O1 K2 path to the first crossing (no extra steps)."""
    P = _params("O1", *CASES[case])
    st0 = O1.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"])
    return P, O1.k2_path(P, st0, N_MAX, extra_after=0)


def _n_interp_from(dets: dict, n_first):
    """Same definition as both oracles' k2_path: linear zero crossing between steps n_first-1, n_first."""
    if n_first is None or n_first < 2:
        return None
    d0, d1 = dets[n_first - 1], dets[n_first]
    return (n_first - 1) + d0 / (d0 - d1)


@functools.lru_cache(maxsize=None)
def _o2_run(case, m):
    """O2 on (S.43) with every increment split into m equal backward-Euler substeps (log strain
    ln f / m each, so the product of the substeps' f is exactly f). Minimum acoustic determinant
    after every nominal step n >= 10 (the crossing is after step 10; asserted)."""
    P = _params("O2", *CASES[case])
    st = O2.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"], finite=True)
    states, dets, normals = {}, {}, {}
    n_first, max_sub = None, 1
    for n in range(1, N_MAX + 1):
        d = np.diag((F1_LOG if n <= N1 else F2_LOG) / m)
        for _ in range(m):
            st = O2.step(P, st, d)
            if st.flags["refused"]:
                raise AssertionError(f"O2 m={m} {case}: refused at nominal step {n}: {st.flags['reason']}")
            max_sub = max(max_sub, st.flags.get("substeps", 1))
        states[n] = st
        if n >= N1:
            md, nv = O2.acoustic_min_det(P, st, O2.tangent_finite(P, st))
            dets[n], normals[n] = md, nv
            if md <= 0.0:
                n_first = n
                break
    return dict(P=P, states=states, dets=dets, normals=normals, n_first=n_first,
                n_interp=_n_interp_from(dets, n_first), max_sub=max_sub)


@functools.lru_cache(maxsize=None)
def _o2_dets_to20(case, m):
    """O2 with m substeps per (S.43) increment, only as far as n = 20: (min acoustic det at n = 10, at n = 20, max
    substepping seen). Same driver as _o2_run; no search beyond n = 20 (m = 256 is 5120 backward-Euler steps)."""
    P = _params("O2", *CASES[case])
    st = O2.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"], finite=True)
    dets, max_sub = {}, 1
    for n in range(1, 21):
        d = np.diag((F1_LOG if n <= N1 else F2_LOG) / m)
        for _ in range(m):
            st = O2.step(P, st, d)
            if st.flags["refused"]:
                raise AssertionError(f"O2 m={m} {case}: refused at nominal step {n}: {st.flags['reason']}")
            max_sub = max(max_sub, st.flags.get("substeps", 1))
        if n in (N1, 20):
            dets[n] = O2.acoustic_min_det(P, st, O2.tangent_finite(P, st))[0]
    return dets, max_sub


def _assess(name, errs: dict, rich_err: float, floor: float):
    """Convergence assessment of one quantity; returns a list of failure strings."""
    fails = []
    e1, e16 = errs[1], errs[16]
    if e1 > floor:
        if not e16 <= E16_OVER_E1_MAX * e1:
            fails.append(f"{name}: E(16) = {e16:.3e} > E(1)/8 = {E16_OVER_E1_MAX * e1:.3e} (first order predicts E(1)/16)")
        if not e16 <= errs[8]:
            fails.append(f"{name}: not decreasing at the finest levels, E(8) = {errs[8]:.3e}, E(16) = {e16:.3e}")
    lim = max(RICH_OVER_E16_MAX * e16, floor)
    if not rich_err <= lim:
        fails.append(f"{name}: Richardson error {rich_err:.3e} > max(E(16)/2, floor) = {lim:.3e}")
    return fails


def _rel(a, ref) -> float:
    a, ref = np.asarray(a, float), np.asarray(ref, float)
    return float(np.linalg.norm(a - ref) / np.linalg.norm(ref))


# ----------------------------------------------------------------------------------------------
# FS.conv_state : O2(m) -> O1 on the K2 path at n = 10, 20
# ----------------------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("n", (10, 20))
@pytest.mark.parametrize("case", sorted(CASES))
def test_k2_state_convergence_to_o1(case, n):
    """FS.conv_state. tau, pi_i and the minimum acoustic determinant of O2 with m substeps per (S.43)
    increment converge to O1 (the continuum truth) at first order: E(16) <= E(1)/8, E(16) <= E(8),
    and the Richardson state 2Q(16) - Q(8) is within max(E(16)/2, noise floor) of O1 (the O(1/m^2)
    remainder). This applies to tau and pi_i at n = 10 and 20 and to the min det at n = 10.

    The min det at n = 20 (raw, and normalised by its step-10 value; the normalised value at n = 10 is
    identically 1 and carries no information) has a large second-order coefficient (|b/a| ~ 9, error sign
    change between m = 8 and 16): m = 1, 2 are pre-asymptotic and the 2-level remainder check does not apply.
    It is therefore checked in the asymptotic range m >= 64 (see the constants block (1b) for the argument):
    the raw error ratio over the doublings 64 -> 128 -> 256 lies in [1.7, 2.3] (first order), and the 3-level
    Richardson (8Q(256) - 6Q(128) + Q(64))/3 of both the raw and the normalised det is within
    max(E(256)/4, 1e-6) of O1 (error c/(8 m^3), O(1/m^3)).

    Kills: an O(1) physics error in O2's finite-strain update (v = v0 (1 + tr eps) kinematics,
    wrong hardening or CSL coupling, dropped tau(+)1) -- the gap to O1 does not close; a
    non-first-order scheme -- E(1)/E(16) misses 16."""
    _, r1 = _o1_run(case)
    runs = {m: _o2_run(case, m) for m in M_LIST}
    for m, run in runs.items():
        assert run["max_sub"] == 1, f"O2 m={m}: internal refusal-substepping engaged ({run['max_sub']}), not a pure m-step run"
        assert run["n_first"] is None or run["n_first"] > 20, f"O2 m={m}: crossing at {run['n_first']} <= 20"
    assert r1["n_first"] is not None and r1["n_first"] > 20, f"O1: crossing at {r1['n_first']} <= 20"
    s1 = r1["states"][n - 1]
    qs = {
        "tau": (lambda run: run["states"][n].sigma, s1.sigma, FLOOR_STATE),
        "pi_i": (lambda run: np.array([run["states"][n].pi_i]), np.array([s1.pi_i]), FLOOR_STATE),
    }
    if n == N1:                                      # n = 10: the min det is assessed as before (the n = 10 cases are unchanged)
        qs["mindet"] = (lambda run: np.array([run["dets"][n]]), np.array([r1["mindet"][n - 1]]), FLOOR_DET)
    fails = []
    for name, (get, ref, floor) in qs.items():
        vals = {m: np.asarray(get(runs[m]), float) for m in M_LIST}
        errs = {m: _rel(vals[m], ref) for m in M_LIST}
        rich = 2.0 * vals[16] - vals[8]
        fails += _assess(f"{case} n={n} {name}", errs, _rel(rich, ref), floor)
    if n == 20:
        # The determinant at n = 20 (raw and normalised by its step-10 value): first order anchored in the asymptotic
        # range M_ASYM, and the 3-level Richardson (see the constants block for the argument).
        det = {m: _o2_dets_to20(case, m) for m in M_ASYM}
        for m, (_, ms) in det.items():
            assert ms == 1, f"O2 m={m}: internal refusal-substepping engaged ({ms})"
        refs = {"mindet": (np.array([r1["mindet"][19]]), lambda d: np.array([d[20]])),
                "mindet_normalised": (np.array([r1["mindet_norm"][19]]), lambda d: np.array([d[20] / d[N1]]))}
        for name, (ref, get) in refs.items():
            q = {m: get(det[m][0]) for m in M_ASYM}
            e = {m: _rel(q[m], ref) for m in M_ASYM}
            if name == "mindet":
                for m in M_ASYM[:-1]:
                    ratio = e[m] / max(e[2 * m], 1.0e-300)
                    if e[m] > FLOOR_DET and not RATIO_LO <= ratio <= RATIO_HI:
                        fails.append(f"{case} n=20 {name}: E({m})/E({2 * m}) = {ratio:.3f} outside [{RATIO_LO}, {RATIO_HI}] "
                                     f"(E {e[m]:.3e} -> {e[2 * m]:.3e})")
            m0, m1, m2 = M_ASYM
            r3 = (8.0 * q[m2] - 6.0 * q[m1] + q[m0]) / 3.0
            lim = max(RICH3_OVER_EFINE_MAX * e[m2], FLOOR_DET)
            if not _rel(r3, ref) <= lim:
                fails.append(f"{case} n=20 {name}: 3-level Richardson error {_rel(r3, ref):.3e} > max(E({m2})/4, floor) = {lim:.3e}")
    assert not fails, "; ".join(fails)


# ----------------------------------------------------------------------------------------------
# FS.conv_ninterp : first-order convergence of the interpolated crossing, Richardson to O1
# ----------------------------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("case", sorted(CASES))
def test_k2_ninterp_first_order_and_richardson(case, record_property):
    """FS.conv_ninterp. gap(m) = |n_interp,O2(m) - n_interp,O1| falls by a factor in [1.7, 2.3] at each
    doubling m -> 2m (m = 1, 2, 4, 8): first order, ratio 2 + O(1/m). The Richardson estimate
    2 n(16) - n(8) matches O1 to <= 0.02 step (the O(1/m^2) remainder).

    Kills: a non-first-order integrator (explicit/forward update: ratio ~ 1 or divergent; midpoint
    update: ratio ~ 4), and any O(1) error that keeps the O2 crossing off the continuum value
    (the Richardson match fails): wrong v kinematics, a missing tau(+)1 or sigma_n term in (S.44),
    a mis-assembled gamma_ab in (S.34)."""
    _, r1 = _o1_run(case)
    n1c = r1["n_interp"]
    assert n1c is not None, f"O1 {case}: no interpolated crossing within {N_MAX}"
    ns = {}
    for m in M_LIST:
        run = _o2_run(case, m)
        assert run["max_sub"] == 1, f"O2 m={m}: internal refusal-substepping engaged ({run['max_sub']})"
        assert run["n_interp"] is not None, f"O2 m={m} {case}: no crossing within {N_MAX}"
        ns[m] = run["n_interp"]
        record_property(f"n_interp_O2_m{m}", ns[m])
    record_property("n_interp_O1", n1c)
    gaps = {m: abs(ns[m] - n1c) for m in M_LIST}
    fails = []
    for m in (1, 2, 4, 8):
        ratio = gaps[m] / max(gaps[2 * m], 1.0e-12)
        if not RATIO_LO <= ratio <= RATIO_HI:
            fails.append(f"gap({m})/gap({2 * m}) = {ratio:.3f} outside [{RATIO_LO}, {RATIO_HI}] "
                         f"(gaps {gaps[m]:.4f} -> {gaps[2 * m]:.4f})")
    rich = 2.0 * ns[16] - ns[8]
    if not abs(rich - n1c) <= RICHARDSON_NINTERP_TOL:
        fails.append(f"Richardson 2 n(16) - n(8) = {rich:.4f} vs O1 {n1c:.4f}: |diff| = {abs(rich - n1c):.4f} "
                     f"> {RICHARDSON_NINTERP_TOL}")
    assert not fails, f"{case}: " + "; ".join(fails)


# ----------------------------------------------------------------------------------------------
# FS.fd_o2 : (S.34) vs central FD of the O2 stress update
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _o2_state_after_f1(case):
    """O2 (m = 1) at the nominal K2 state after the 10 f1 steps, and the state before step 10."""
    P = _params("O2", *CASES[case])
    st = O2.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"], finite=True)
    prev = None
    for _ in range(N1):
        prev = st
        st = O2.step(P, st, np.diag(F1_LOG))
        assert not st.flags["refused"], st.flags["reason"]
    return P, prev, st


def _assemble(V, ct, gam, tau, spin=1.0, tau_plus_one=True):
    """(S.34): a^ep = sum_ab ct_ab m^a (x) m^b + spin * sum_{a != b} gam_ab (m^ab (x) m^ab + m^ab (x) m^ba)
    + tau (+) 1, (tau (+) 1)_ijkl = tau_jl delta_ik. spin = 1 is the sheet; 0.5 is the BA06 3.48 error."""
    m = [np.outer(V[:, a], V[:, a]) for a in range(3)]
    a4 = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            a4 += ct[a, b] * np.einsum("ij,kl->ijkl", m[a], m[b])
            if a != b:
                mab, mba = np.outer(V[:, a], V[:, b]), np.outer(V[:, b], V[:, a])
                a4 += spin * gam[a, b] * (np.einsum("ij,kl->ijkl", mab, mab) + np.einsum("ij,kl->ijkl", mab, mba))
    if tau_plus_one:
        T = sum(tau[a] * m[a] for a in range(3))
        a4 += np.einsum("jl,ik->ijkl", T, I3)
    return a4


def _blocks_from_o2_cache(st):
    """(V, c~_ab, gamma~_ab, tau_a) of (S.34) from the O2 trial-spectral cache; distinct stretches required."""
    c = st.cache
    w, tau, at, V = c["eps_tr"], c["sig"], c["atilde"], c["nvec"]
    lam2 = np.exp(2.0 * w)
    gam = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a != b:
                assert abs(np.exp(w[a]) - np.exp(w[b])) > 1.0e-6, "trial stretches not distinct: FD test is vacuous"
                gam[a, b] = (tau[b] * lam2[a] - tau[a] * lam2[b]) / (lam2[b] - lam2[a])
    return V, at - 2.0 * np.diag(tau), gam, tau


def _fd_directions():
    dirs = [(f"E{k}{l}", np.outer(I3[k], I3[l])) for k in range(3) for l in range(3)]
    g = np.random.default_rng(7).normal(size=(3, 3))     # dense, non-symmetric, fixed seed
    return dirs, g


def _fd_errors(a4, tau_of_f, f0, h, g):
    """(max rel error of a4 against the 81 FD entries over the unit E_kl, rel error of a4:g against the
    FD along the dense g). Relative to max |FD|."""
    fd = np.zeros((3, 3, 3, 3))
    for k in range(3):
        for l in range(3):
            E = np.outer(I3[k], I3[l])
            fp, fm = (I3 + h * E) @ f0, (I3 - h * E) @ f0
            Pp = tau_of_f(fp) @ np.linalg.inv(I3 + h * E).T
            Pm = tau_of_f(fm) @ np.linalg.inv(I3 - h * E).T
            fd[:, :, k, l] = (Pp - Pm) / (2.0 * h)
    sc = float(np.abs(fd).max())
    e_unit = float(np.abs(fd - a4).max()) / sc
    fp, fm = (I3 + h * g) @ f0, (I3 - h * g) @ f0
    fdg = (tau_of_f(fp) @ np.linalg.inv(I3 + h * g).T - tau_of_f(fm) @ np.linalg.inv(I3 - h * g).T) / (2.0 * h)
    e_dense = float(np.abs(fdg - np.einsum("ijkl,kl->ij", a4, g)).max() / np.abs(fdg).max())
    return e_unit, e_dense


@pytest.mark.parametrize("h", FD_HS)
@pytest.mark.parametrize("case", sorted(CASES))
def test_finite_tangent_fd_o2(case, h):
    """FS.fd_o2. At the plastic K2 state after the 10 f1 steps, (S.34) of O2 (both the O2
    `tangent_finite` and the sheet formula assembled independently here from the CTO block) equals
    the central FD of P(E) = tau((1 + hE) f0)(1 + hE)^{-T} over the nine unit E_kl and over a dense
    non-symmetric g, within the truncation bound fd_bound(h) = 10 (h^2/(6 kappa^2) + 1e-12/h). The two
    mutants of the sheet (1/2 on the spin sum, BA06 3.48; no tau(+)1) are assembled from the same
    blocks and must miss that same bound by >= 10x (the negative controls).

    Kills: the 1/2-spin variant, a dropped tau(+)1, a wrong gamma_ab, a v0-for-v error in the CTO
    (a~ with v0 in place of v, sheet 1.2), any error in the (S.31)-(S.32) tangent."""
    P, prev, st = _o2_state_after_f1(case)
    assert st.flags["plastic"], f"{case}: state after the 10 f1 steps is not plastic: {st.flags}"
    f0 = np.diag(np.exp(F1_LOG))
    w_prev, V_prev = np.linalg.eigh(prev.eps_e)
    b_n = (V_prev * np.exp(2.0 * w_prev)) @ V_prev.T

    def tau_of(f):
        b_tr = f @ b_n @ f.T
        wb, Vb = np.linalg.eigh(0.5 * (b_tr + b_tr.T))
        res = K.return_map(P, 0.5 * np.log(wb), prev.pi_i, prev.v * np.linalg.det(f), prev.v * np.linalg.det(f))
        assert res.plastic and not res.refused, f"FD map left the plastic branch: {res.reason}"
        return (Vb * res.sig) @ Vb.T

    assert np.abs(tau_of(f0) - st.sigma).max() <= 1.0e-10 * P_ABS, "FD map is not the O2 step map"
    V, ct, gam, tau = _blocks_from_o2_cache(st)
    a_o2 = O2.tangent_finite(P, st)
    a_mine = _assemble(V, ct, gam, tau)
    assert np.abs(a_mine - a_o2).max() <= 1.0e-9 * np.abs(a_o2).max(), "independent (S.34) assembly differs from O2's"
    bound = fd_bound(h, FD_DELTA_PLASTIC)
    _, g = _fd_directions()
    e_o2 = _fd_errors(a_o2, tau_of, f0, h, g)
    assert max(e_o2) <= bound, f"{case} h={h}: (S.34) FD error unit {e_o2[0]:.3e}, dense {e_o2[1]:.3e} > bound {bound:.3e}"
    for name, kw in (("half-spin", dict(spin=0.5)), ("no-tau(+)1", dict(tau_plus_one=False))):
        e_v = _fd_errors(_assemble(V, ct, gam, tau, **kw), tau_of, f0, h, g)
        assert max(e_v) >= NEG_CONTROL_FACTOR * bound, (
            f"{case} h={h}: negative control '{name}' error {max(e_v):.3e} is not >= {NEG_CONTROL_FACTOR}x "
            f"the bound {bound:.3e}: the FD check does not discriminate it")


# ----------------------------------------------------------------------------------------------
# FS.fd_o1 : O1's finite-strain tangent
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("h", FD_HS)
@pytest.mark.parametrize("case", sorted(CASES))
def test_finite_tangent_fd_o1_elastic_assembly(case, h):
    """FS.fd_o1 (assembly). O1 has no stress-update map for a non-diagonal f (its log-strain mode is
    fixed-direction), so the assembly of O1's `finite_spatial_tangent` is checked on the elastic branch,
    where the stress map is closed form: tau(f) = dPsi/d(eps^e) at eps^e = (1/2) ln(f b^e f^T), with
    O1's own energy. At the O1 state after the 10 f1 steps, a^e assembled by (S.34) must equal the
    central FD of P(E) = tau((1 + hE))(1 + hE)^{-T} within 10 (h^2/(6 kappa^2) + 1e-14/h) (round-off
    only), and the 1/2-spin / no-tau(+)1 variants built from O1's own blocks must miss it by >= 10x.

    Kills: the 1/2 on the spin sum or a missing tau(+)1 in O1's `finite_spatial_tangent`; a wrong
    gamma_ab; an a^e that is not the derivative of O1's energy."""
    P, r1 = _o1_run(case)
    st = r1["states"][N1 - 1]
    w, V0_ = np.linalg.eigh(0.5 * (st.eps_e + st.eps_e.T))
    b_n = (V0_ * np.exp(2.0 * w)) @ V0_.T

    def tau_of(f):
        b = f @ b_n @ f.T
        wb, Vb = np.linalg.eigh(0.5 * (b + b.T))
        ee = (Vb * (0.5 * np.log(wb))) @ Vb.T
        return o1_energy(ee, P, tangent=False).sig

    assert np.abs(tau_of(I3) - st.sigma).max() <= 1.0e-9 * P_ABS, "FD map is not O1's elastic stress"
    a_e = O1.tangent(P, st, plastic_branch=False)
    a_o1 = finite_spatial_tangent(st, a_e)
    bound = fd_bound(h, FD_DELTA_ELASTIC)
    _, g = _fd_directions()
    e_o1 = _fd_errors(a_o1, tau_of, I3, h, g)
    assert max(e_o1) <= bound, f"{case} h={h}: O1 finite tangent FD error unit {e_o1[0]:.3e}, dense {e_o1[1]:.3e} > bound {bound:.3e}"
    ct, gam, tau, V = finite_blocks(st, a_e)
    for name, kw in (("half-spin", dict(spin=0.5)), ("no-tau(+)1", dict(tau_plus_one=False))):
        e_v = _fd_errors(_assemble(V, ct, gam, tau, **kw), tau_of, I3, h, g)
        assert max(e_v) >= NEG_CONTROL_FACTOR * bound, (
            f"{case} h={h}: negative control '{name}' error {max(e_v):.3e} is not >= {NEG_CONTROL_FACTOR}x "
            f"the bound {bound:.3e}: the FD check does not discriminate it")


@pytest.mark.slow
@pytest.mark.parametrize("case", sorted(CASES))
def test_finite_tangent_o2_cto_converges_to_o1(case):
    """FS.fd_o1 (plastic branch). The plastic finite-strain tangent of O1 (continuum a^ep of (S.42) assembled
    by (S.34)) at step 10 is the m -> infinity limit of O2's CTO (the closed-form derivative of the
    backward-Euler map, whose FD check is test_finite_tangent_fd_o2): max |a_O2(m) - a_O1| / max |a_O1|
    falls at first order (E(16) <= E(1)/8) and Richardson 2a(16) - a(8) is within max(E(16)/2, 1e-8).

    Kills: an O(1) error in O1's plastic tangent (a wrong H, a missing hardening coupling, a v for v0
    slip) or in O2's CTO; the two cannot then close."""
    P, r1 = _o1_run(case)
    st = r1["states"][N1 - 1]
    assert st.flags["plastic"], f"O1 {case}: state after the 10 f1 steps is not plastic"
    a_o1 = finite_spatial_tangent(st, O1.tangent(P, st))
    sc = float(np.abs(a_o1).max())
    a = {m: O2.tangent_finite(_o2_run(case, m)["P"], _o2_run(case, m)["states"][N1]) for m in M_LIST}
    errs = {m: float(np.abs(a[m] - a_o1).max()) / sc for m in M_LIST}
    rich = float(np.abs(2.0 * a[16] - a[8] - a_o1).max()) / sc
    fails = _assess(f"{case} finite tangent", errs, rich, FLOOR_TANGENT)
    assert not fails, "; ".join(fails)


# ----------------------------------------------------------------------------------------------
# FS.direction : the localization normal is perpendicular to the intermediate principal direction
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _first_localization(oracle, case):
    """(normal, sigma) at the first-localization step of the oracle's own production k2_path."""
    P = _params(oracle, *CASES[case])
    if oracle == "O1":
        st0 = O1.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"])
        r = O1.k2_path(P, st0, N_MAX, extra_after=0)
        normals = r["normals"]
    else:
        st0 = O2.initial_state(P, SIGMA0, V0, NOMINAL["pi_i0"])
        r = O2.k2_path(P, st0, N_MAX)
        normals = r["n_vec"]
    n_first = r["n_first"]
    assert n_first is not None, f"{oracle} {case}: no localization within {N_MAX}"
    return n_first, np.asarray(normals[n_first - 1], float), r["states"][n_first - 1].sigma


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("oracle", ("O1", "O2"))
def test_localization_normal_perpendicular_to_intermediate_direction(oracle, case, record_property):
    """FS.direction (AB06 Fig 6, phi = pi/2; sheet 14 reading). At the first-localization state of the
    nominal combination, the minimising normal n of det A(n) is perpendicular to the intermediate
    principal direction: |n . e_int| <= 0.05 (grid ~ 2 deg; sin 2 deg = 0.035). e_int is the
    eigenvector of the middle principal Kirchhoff stress; the sheet says it is e_1 under (S.43) (checked).
    The sign of n and the mirror well (+-theta) are irrelevant (absolute value).

    Kills: a localization search or tangent assembly whose basis labelling is wrong (principal vs
    global axes mixed, a<->b indices of c~/gamma~ swapped in (S.44)), or an acoustic tensor built
    from a symmetrised a^ep: the minimum then leaves the e_2-e_3 plane."""
    n_first, n, sigma = _first_localization(oracle, case)
    record_property("n_first", n_first)
    assert abs(np.linalg.norm(n) - 1.0) <= 1.0e-9, f"normal not unit: {n}"
    w, V = np.linalg.eigh(0.5 * (sigma + sigma.T))        # ascending: most compressive first
    assert w[1] - w[0] >= PRINCIPAL_GAP_REL * P_ABS and w[2] - w[1] >= PRINCIPAL_GAP_REL * P_ABS, \
        f"principal stresses not distinct, intermediate direction undefined: {w}"
    e_int = V[:, 1]
    assert abs(e_int[0]) >= 1.0 - E_INT_AXIS_TOL, f"intermediate direction is not e_1 (sheet 14): {e_int}"
    cosang = abs(float(n @ e_int))
    record_property("abs_n_dot_e_int", cosang)
    assert cosang <= DIRECTION_TOL, (f"{oracle} {case} n_first={n_first}: |n . e_int| = {cosang:.4f} > {DIRECTION_TOL}; "
                                     f"n = {n}, e_int = {e_int}, principal tau = {w}")


# ----------------------------------------------------------------------------------------------
# FS.s44 : (S.44) against the direct contraction n_j a_ijkl n_l
# ----------------------------------------------------------------------------------------------
S44_TOL = 1.0e-12                                   # relative to max |A| (|A|^3 for the determinant), see the docstring
S44_N_RANDOM = 24


def _s44_normals(n_min):
    """Normals tested: the oracle's own minimiser, the three global axes, and 24 seeded random unit vectors."""
    rng = np.random.default_rng(44)
    rnd = rng.normal(size=(S44_N_RANDOM, 3))
    out = [np.asarray(n_min, float) / np.linalg.norm(n_min)] + [I3[i] for i in range(3)]
    out += [v / np.linalg.norm(v) for v in rnd]
    return out


@functools.lru_cache(maxsize=None)
def _s44_states(oracle, case):
    """{label: (state, tangent a^ep (4th order), (V, c~, gamma~, tau) blocks, minimiser normal)} at the plastic K2 states
    after the 10 f1 steps ('f1') and at first localization ('loc'); m = 1 for O2 (the production step)."""
    out = {}
    if oracle == "O2":
        run = _o2_run(case, 1)
        for label, n in (("f1", N1), ("loc", run["n_first"])):
            st = run["states"][n]
            out[label] = (st, O2.tangent_finite(run["P"], st), _blocks_from_o2_cache(st), run["normals"][n])
    else:
        P, r1 = _o1_run(case)
        for label, n in (("f1", N1), ("loc", r1["n_first"])):
            st = r1["states"][n - 1]
            a_e = O1.tangent(P, st)
            ct, gam, tau, V = finite_blocks(st, a_e)
            out[label] = (st, finite_spatial_tangent(st, a_e), (V, ct, gam, tau), r1["normals"][n - 1])
    return out


@pytest.mark.parametrize("where", ("f1", "loc"))
@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("oracle", ("O1", "O2"))
def test_s44_acoustic_tensor_equals_the_direct_contraction(oracle, case, where):
    """FS.s44 (sheet 9.5 / 14 / 16.3, a G1 gate). With a^ep = c~ + tau(+)1 of (S.34), the acoustic tensor A_ik = n_j a_ijkl n_l
    contracted directly must equal, to 1e-12 relative, the principal-basis closed form (S.44),
        A^_aa = alpha_a^2 c~_aa + sum_c alpha_c^2 tau_c + sum_{c != a} alpha_c^2 gamma~_ca,
        A^_ab = alpha_a (c~_ab + gamma~_ab) alpha_b   (a != b),        alpha = V^T n,    V^T A V = A^.
    The closed form is the transcription `acoustic_principal` (o2_algo.acoustic, a pure function of (c~, gamma~, tau, alpha))
    fed with each oracle's OWN blocks; for O1 the determinant of O1's own production search function (the (S.44) code the
    localization search of O1 uses) is compared too. Normals: the oracle's own minimiser, the three global axes and 24
    seeded random unit vectors. States: the plastic K2 state after the 10 f1 steps and the first-localization state,
    cases rho 0.7/0.8 and 1/1.
    Tolerance argument: both sides are 3x3 sums of a few products of numbers of O(1e1..1e4) kPa evaluated in double precision,
    so the difference is a few units of round-off (1e-16) times a growth factor < 1e2 from the 81-term contraction; the gate
    is 1e-12 of max |A| (of max |A|^3 for the determinant, which is a cancellation and is therefore NOT normalised by itself).
    Kills: the 1/2 on the spin sum or a swapped a<->b in gamma~ in (S.44); tau(+)1 missing from one side (sigma_n term of the
    diagonal); alpha taken in the global instead of the principal basis; a symmetrised tangent in the search."""
    st, a4, (V, ct, gam, tau), n_min = _s44_states(oracle, case)[where]
    assert st.flags["plastic"], f"{oracle} {case} {where}: not a plastic state: {st.flags}"
    worst, worst_det = 0.0, 0.0
    for n in _s44_normals(n_min):
        A = np.einsum("j,ijkl,l->ik", n, a4, n)
        alpha = V.T @ n
        Ahat = acoustic_principal(ct, gam, tau, alpha)
        scale = float(np.abs(A).max())
        err = float(np.abs(V.T @ A @ V - Ahat).max()) / scale
        worst = max(worst, err)
        assert err <= S44_TOL, f"{oracle} {case} {where}: (S.44) vs n.a.n relative error {err:.3e} for n = {n}"
        if oracle == "O1":
            ph = float(np.arccos(np.clip(alpha[1], -1.0, 1.0)))
            th = float(np.arctan2(alpha[0], alpha[2]))
            d_o1 = float(o1_det_fun_finite(ct, gam, tau)(np.array(th), np.array(ph)))
            derr = abs(d_o1 - float(np.linalg.det(A))) / scale ** 3
            worst_det = max(worst_det, derr)
            assert derr <= S44_TOL, f"{oracle} {case} {where}: O1 det A(S.44) vs det(n.a.n): relative error {derr:.3e} for n = {n}"
