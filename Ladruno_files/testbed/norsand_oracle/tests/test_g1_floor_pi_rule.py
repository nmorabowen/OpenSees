"""WP-144 round 4, G1 (Zone B): the image-pressure rule at a post floor, (S.57)-(S.60), on both oracles.  Equation sheet 144a 9.1 step 5g, 9.7
(S.57)-(S.60) and the wet-floor table, 13.12 (C), 13.14b, 13.16 (K1.16).

TEST RULE (plan 5): every expected value and tolerance in this file is (a) a number PRINTED in the sheet (round 4b, commit 6ae9fe5ab), compared to
its printed digits, (b) a closed form of the sheet coded in g1_hf_common.py or here (the surface through a stress (S.12)/(S.57), c_N (S.58), the
(S.50) floor), or (c) a stated convergence argument.  All of it was written before the oracles were run on these paths; no oracle output is an
expected value.  A failing assertion is a finding against the sheet or an oracle, never a reason to edit an oracle or loosen a bound.

Gates (ids follow the sheet)
  K1.12 (C)  BA06 wet-side post floor with the rule: pi_i = -40.516314098 (was -40.198569107), sigma_f, eta_f, F(sigma_f, pi_BE) = +0.010899 p_min,
             pattern -Pf*, c_N p_min = 23.2379; cases (A), (B) unchanged (no post floor: nothing counted); p_min = 0 turns the rule off
  wet table  the 12 printed cells of the wet-floor table (9.7 item 3): F(sigma_f, pi_BE)/p_min, pi_BE, eta_f, pi_i by (S.57), HAR TIMs, pure shear
  K1.16      the H19 states `s185.txt` / `stall.txt` (copied to tests/data_h19/): the printed table (steps 1, 2, 60), the step-1 detail, the stall row
             asserted by digits only when the step ran at m = 1 and by the invariants otherwise (9.7, round 4b), F <= F_tol and |pi_i| >= c_N p_min
  S.58       seeded random census of post-floor events: F <= F_tol, |pi_i| >= c_N |p|(1 + ...), never shrunk, pi_i = pi_f(sigma_f) when it fired
  S.59       the forced (1/2, 1/2) HAR fixture: the chained tangent against central FD at three h, and the same chain with the replaced row (monkeypatched
             kernel.chain_propagate: M-R1 row not replaced, M-R3 row added) far from FD, so the row matters
  counter    n_f_pi: +1 per reset, sums over sub-increments, nothing on an elastic / non-post-floored step, nothing on a refused increment
  conv       O2 -> O1 at first order on a path where the rule fires (the sheet's rate table, O1 with the rule on); the gap to the Koiter form
             (O1 floor_rule=False) is the sheet's plateau; the round-3b split (rule off) converges to the Koiter form (sheet table, third column)

MUTANTS (sheet 9.7 "Tests and mutants", round 4; each test names the ones it kills).
  M-R1 chain row not replaced   M-R2 zeta_y term dropped from f_a   M-R3 row added to the (S.54) row   M-R4 rule absent / evaluated at the pre-floor stress
  M-R5 a^e of the pre-floor state in gamma   M-R6 row built before Phi^post   M-R7 rule applied when F <= 0 (surface shrunk)   M-R8 rule on a step
  without a post floor / elastic step   M-R9 single-step tangent altered by the rule
"""
import dataclasses
import math
import os

import numpy as np
import pytest

import g1_hf_common as C
from conftest import ORACLES

O1, O2 = ORACLES["O1"], ORACLES["O2"]
K = O2.kernel
ONES = np.ones(3)
SQ23 = C.SQ23
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data_h19")

# ----------------------------------------------------------------------------------------------------------------------
# the K3 set of sheet 13.16 (printed): HAR, WW, fork CSL, smooth cap, p_min = 0.505 kPa
# ----------------------------------------------------------------------------------------------------------------------
K3 = dict(energy="HAR", k=298.34244992965336, g=383.5831499095543, n_e=0.5, p_a=101.0, M=1.25, N=0.4686136860640983,
          N_bar=0.16034970350694513, rho=0.5828981510111576, rho_bar=0.7213710496292179, zeta="WW", chi=-2.8923322243122183,
          h=308.3223484879116, csl_mode="fork", e0=0.934, lambda_c=0.019, xi=0.7, cap="smooth", c1=0.05, c2=0.15, p_min=0.505)
PMIN3 = 0.505
FTOL3 = 1.0e-10 * 101.0                     # F_tol = 1e-10 p_ref, p_ref = p_a (9.1 step 2, 9.7)
DELTA_F = FTOL3 / PMIN3                     # delta_F of (S.58), 2e-8


def c_N(N):
    """(S.58): c_N = (1-N)^((1-N)/N) (N > 0), e^-1 (N = 0): |pi_i|/|p| at the tip eta = 0."""
    return (1.0 - N) ** ((1.0 - N) / N) if N > 0.0 else math.exp(-1.0)


def F_of(kw, sig, pi):
    return C.F_yield(kw["M"], kw["N"], kw["rho"], sig, pi)


def pi_f_of(kw, sig):
    """(S.57): pi_f = the inverse of (S.12) at (p_f, eta_f = zeta(theta_f) q_f/|p_f|), the surface through the stress `sig`; returns (pi_f, eta_f)."""
    p, q, th = C.invariants_pq(sig)
    eta = C.zeta_ww(th, kw["rho"]) * q / abs(p)
    return C.pi_S53(kw["M"], kw["N"], p, eta), eta


def q_of(sig):
    return C.invariants_pq(sig)[1]


# ----------------------------------------------------------------------------------------------------------------------
# the H19 states (data copied from norsand_calib/out/h19/; format at the top of out/h19/h19_replay.cpp):
#   id ee11 ee22 eg12 s11 s22 s12 s33 pi_i v v0 epv eps at_floor n_f_tr n_f_post eps_f_v W_f d11 d22 dg12   (plane strain; eg12, dg12 engineering)
# eps^e_33 is NOT in the file: it is recovered from s33 by the oracle's inverse map and compared to the printed value of 13.16.
# ----------------------------------------------------------------------------------------------------------------------
def h19_records():
    recs = {}
    for fn in ("s185.txt", "stall.txt"):
        with open(os.path.join(DATA, fn)) as fh:
            for line in fh:
                t = line.split()
                if not t or t[0].startswith("#"):
                    continue
                x = [float(a) for a in t[1:]]
                recs[t[0]] = dict(ee=x[0:3], s=x[3:7], pi=x[7], v=x[8], d=x[17:20])
    return recs


def h19_state(oname, P, rec, norm=2.0e-4):
    """(initial state, increment tensor scaled to |d eps| = (d11^2 + d22^2 + dg12^2/2)^(1/2) = norm (13.16))."""
    s11, s22, s12, s33 = rec["s"]
    sig = np.array([[s11, s12, 0.0], [s12, s22, 0.0], [0.0, 0.0, s33]])
    st = ORACLES[oname].initial_state(P, sig, rec["v"], rec["pi"])
    d11, d22, dg = rec["d"]
    sc = norm / math.sqrt(d11 * d11 + d22 * d22 + 0.5 * dg * dg)
    return st, sc * np.array([[d11, 0.5 * dg, 0.0], [0.5 * dg, d22, 0.0], [0.0, 0.0, 0.0]])


# the printed table of 13.16: start data, step-1 detail, and the rows (sigma11, sigma22, sigma33, sigma12, pi_i) at the listed steps
K116 = {
    "e7367g1_s185": dict(
        e33=2.099254865382e-3, F0=0.2315, F0_pmin=0.458, pi0=-0.3224962955,
        rows={1: (-0.5101907038, -0.5976456953, -0.4071636009, -0.2151586912, -0.3989652918),
              2: (-0.5082384097, -0.6104039680, -0.3963576223, -0.2510036052, -0.4461726222),
              60: (-0.5114970385, -0.6473686614, -0.3561343001, -0.3336845167, -0.5923861407)},
        detail=dict(dlam=1.427733871e-4, p_c=-0.365744965, eta_c=1.092426, pi_be=-0.3245384951, F_post=0.175979368,
                    eta_f=0.922543205, dfv=7.223598e-5)),
    "e7367g2_s185": dict(
        e33=2.100404994193e-3, F0=0.1472, F0_pmin=0.291, pi0=-0.3666004637,
        rows={1: (-0.4889744155, -0.6273143992, -0.3987111852, -0.2332250428, -0.4276028897),
              2: (-0.4850426738, -0.6393433615, -0.3906139647, -0.2596017804, -0.4672749361),
              60: (-0.4814251232, -0.6764194559, -0.3571554210, -0.3278584883, -0.5976134600)},
        detail=dict(dlam=1.401240162e-4, p_c=-0.388932212, eta_c=1.177716, pi_be=-0.3676032012, F_post=0.118226999,
                    eta_f=1.026016850, dfv=5.957758e-5)),
    "e7367g1_m53": dict(
        e33=2.088911465473e-3, F0=1.6596, F0_pmin=3.29, pi0=-0.1116066186,
        rows={1: (-0.5071972881, -0.5355837599, -0.4722189520, -0.0691517825, -0.2815026983),
              2: (-0.5031271070, -0.5765889736, -0.4352839194, -0.1789591553, -0.3671148499)},
        detail=dict(dlam=1.528227673e-4, p_c=-0.187753246, eta_c=0.484384, pi_be=-0.1150543895, F_post=1.439692306,
                    eta_f=0.294256201, dfv=1.856529e-4)),
}
ROW_TOL = 5.0e-10        # printed to 10 decimals (rounding 5e-11) + the two solvers' agreement (the sheet: <= 1e-11 over 60 steps)
F_DIGITS_TOL = 1.0e-12   # the printed |F| <= 6e-15 on the digit rows, 'to round-off'


def march(P, st, d, n, stepper=None):
    """n repeated steps of d; returns the list of (state_after) (stops at a refusal)."""
    stepper = stepper or O2.step
    out = []
    for _ in range(n):
        st = stepper(P, st, d)
        out.append(st)
        if st.flags["refused"]:
            break
    return out


class rule_off:
    """Context: the round-3b return map (step 5g absent: pi_i := pi_i^BE whatever F(sigma_f)); the M-R4 mutant."""

    def __enter__(self):
        self._orig = K.floor_pi_rule
        K.floor_pi_rule = lambda P, sig_f, ae_f, pi_be: dict(pi=pi_be, pi_be=pi_be)

    def __exit__(self, *a):
        K.floor_pi_rule = self._orig


# ======================================================================================================================
# sheet checks that need no oracle
# ======================================================================================================================
@pytest.mark.sheet_check
def test_printed_cN_pmin_and_the_k112c_fixture_values_reproduce_from_the_closed_forms():
    """(S.58c) / 13.12 (C) / 13.16 printed digits against the closed forms: c_N p_min = 0.246558167 (K3: N = 0.46861, p_min 0.505) and 23.2379 (K2,
    N = 0.4, p_min 50); the K1.12 (C) start pi_i = -33.125471695 is the inverse of (S.12) at (p = -51, eta = 0.6); the wet-table start (p = -0.505,
    eta = 0.5 M = 0.66545) has pi_i = -0.328007; the K1.16 stall fixtures' F0 = +0.2315 / +0.1472 / +1.6596 kPa.  Pure sheet check.
    KILLS (as a reference for the tests below): a wrong c_N, a wrong (S.12) inverse, a wrong F evaluated on the H19 states."""
    assert abs(c_N(K3["N"]) * PMIN3 - 0.246558167) <= 6e-10
    assert abs(c_N(0.4) * 50.0 - 23.2379) <= 6e-5
    assert abs(C.pi_S53(1.2, 0.4, -51.0, 0.6) - (-33.125471695)) <= 6e-10
    assert abs(C.pi_S53(1.3309, 0.4, -0.505, 0.5 * 1.3309) - (-0.328007)) <= 6e-7
    recs = h19_records()
    for name, X in K116.items():
        r = recs[name]
        s11, s22, s12, s33 = r["s"]
        sig = np.array([[s11, s12, 0.0], [s12, s22, 0.0], [0.0, 0.0, s33]])
        assert abs(r["pi"] - X["pi0"]) <= 6e-11, name
        assert abs(C.invariants_pq(sig)[0] + PMIN3) <= 1e-12, name
        F0 = F_of(K3, sig, r["pi"])
        assert abs(F0 - X["F0"]) <= 6e-5 and abs(F0 / PMIN3 - X["F0_pmin"]) <= (6e-4 if X["F0_pmin"] < 1.0 else 6e-3), (name, F0)


# ======================================================================================================================
# K1.12 (C): BA06 wet-side post floor with the rule
# ======================================================================================================================
def k112c(oname, pmin=50.0):
    """13.12 (C) 'Inputs to rebuild': BA06 p0 -100, kappa 0.01, mu0 5400, alpha0 0; M 1.2, N 0.4, N_bar 0.2, chi -3.5, h 280; WW rho 0.7, rho_bar 0.8; paper
    CSL lambda_tilde 0.0135, v_c0 1.81; cap none; p_min 50.  State p = -51, eta = 0.6, n_hat = dev(-1.2, -0.1, 1.3), pi_i from the inverse of (S.12),
    v_n = 1.70; increment (principal, co-axial) 6e-3 (0.55, 0.40, -0.95) - 1e-5."""
    kw = C.ba06_kw(p_min=pmin)
    P = C.make_hf(oname, **kw)
    sig, th, q = C.sigma_on_surface(1.2, 0.4, 0.7, -51.0, 0.6, C.WET_DIR)
    pi0 = C.pi_S53(1.2, 0.4, -51.0, 0.6)
    st0 = ORACLES[oname].initial_state(P, np.diag(sig), 1.70, pi0)
    deps = np.diag(6e-3 * np.array([0.55, 0.40, -0.95]) - 1e-5)
    return kw, P, st0, deps, th, pi0


def test_k1_12c_o2_post_floor_with_the_rule_printed_values():
    """13.12 (C) with the round-4 rule: theta = 0.454427, pi_i start -33.125471695; return: d lambda = 1.546466e-3, p_c = -48.764495, q_c = 46.421589,
    eta_c = 0.952607, pi_i^BE = -40.198569107; post floor: sigma_f = (-34.8599764962, -34.1946812138, -80.9453422900), eta_f = 0.929068, F(sigma_f,
    pi^BE) = +0.010899 p_min > 0 -> pi_i = -40.516314098 (was -40.198569107), F(sigma_f, pi_i) = 0 to 3e-16 p_min, pattern -Pf*; v, the floor event
    and the tangent (S.32f) are those of round 3b; |pi_i| >= c_N p_min = 23.2379; counters (n_f_tr, n_f_post, n_f_pi) = (0, 1, 1).
    KILLS: M-R4 (rule absent or evaluated at sigma_c, where F = 0: pi_i stays -40.198569107, F = +0.010899 p_min), M-R7 (always reset), M-R8 (reset
    without the post floor), a wrong inverse (S.12) in pi_f, pi_f evaluated at p_c instead of the floored p_f."""
    kw, P, st0, deps, th, pi0 = k112c("O2")
    assert abs(th - 0.454427) <= 6e-7
    st = O2.step(P, st0, deps)
    assert not st.flags["refused"] and st.flags["substeps"] == 1
    assert st.flags["rpattern"] == "-Pf*" and st.flags["floor_pi"] == 1 and st.flags["floor_post"] == 1 and st.flags["floor_tr"] == 0
    assert (st.n_f_tr, st.n_f_post, st.n_f_pi) == (0, 1, 1)
    res = st.cache["res"]
    assert abs(st.dlam - 1.546466e-3) <= 6e-10
    pre = st.cache["floor_events"][-1][0]
    ev, es, _ = C.invariants_eps(pre)
    p_c, q_c = C.ba06_pq(ev, es)
    assert abs(p_c + 48.764495) <= 6e-7 and abs(q_c - 46.421589) <= 6e-7 and abs(st.eta - 0.952607) <= 6e-7
    assert abs(res.pi_be - (-40.198569107)) <= 6e-9
    # sigma_f is printed in the principal order of the start (1-2-3 of d = (-1.2, -0.1, 1.3)): compare as ordered by the direction
    want = np.array([-34.8599764962, -34.1946812138, -80.9453422900])
    assert np.abs(np.diag(st.sigma) - want).max() <= 6e-9, np.diag(st.sigma)
    assert abs(res.eta_f - 0.929068) <= 6e-7
    assert abs(res.F_post / 50.0 - 0.010899) <= 6e-7
    assert abs(st.pi_i - (-40.516314098)) <= 6e-9, st.pi_i
    assert abs(F_of(kw, st.sigma, st.pi_i) / 50.0) <= 1e-12                                # 3e-16 p_min printed
    assert abs(st.pi_i - pi_f_of(kw, st.sigma)[0]) <= 1e-12 * abs(st.pi_i)                 # (S.57): the surface through sigma_f
    assert abs(float(np.trace(st.sigma)) / 3.0 + 50.0) <= 1e-12 * 50.0
    assert abs(st.pi_i) >= c_N(0.4) * 50.0 and abs(c_N(0.4) * 50.0 - 23.2379) <= 6e-5
    assert abs(st.v - 1.70 * math.exp(float(np.trace(deps)))) <= 1e-13 * st.v              # v untouched by the rule


def test_k1_12c_o1_ends_on_the_floor_with_the_rule_in_its_rate_form():
    """O1 (rate oracle) on the K1.12 (C) increment: the wet side, f:a:P < 0, so the split's rate limit is (S.60): the path runs in the 'floor_rule' mode
    on at least one segment, ends at p = -p_min, on the surface through the committed stress (F ~ 0 within the integration tolerance) and with |pi_i| >=
    c_N p_min = 23.2379 (S.58c).  The digits of the one-step BE split (-40.516314098) are NOT asserted for O1: the split and its rate limit differ at O(increment)
    (the convergence is the test below).
    KILLS: M-R4 on O1 (the Koiter form (S.55) everywhere: no floor_rule segment), a rate form that leaves the surface, a lost floor."""
    kw, P, st0, deps, th, pi0 = k112c("O1")
    st = O1.run_path(P, st0, deps[None], rtol=1e-10)[-1]
    assert st.flags["status"] == "ok", st.flags["status"]
    assert st.flags["rule_segments"] >= 1 and st.flags["mode"] in ("floor_rule", "floor", "plastic_floor")
    assert abs(float(np.trace(st.sigma)) / 3.0 + 50.0) <= 1e-6 * 50.0
    assert abs(F_of(kw, st.sigma, st.pi_i)) <= 1e-6 * 50.0
    assert abs(st.pi_i) >= c_N(0.4) * 50.0
    assert abs(st.pi_i - pi_f_of(kw, st.sigma)[0]) <= 1e-5 * abs(st.pi_i)


def test_k1_12c_pmin_zero_turns_the_rule_off_and_the_return_is_the_round_3b_one():
    """9.7: with p_min = 0 the post floor never acts and the model is the round-3b one: the K1.12 (C) increment is an ordinary plastic step ending ABOVE the
    would-be floor at p_c = -48.764495 with pi_i = pi_i^BE = -40.198569107 (printed), d lambda 1.546466e-3, no event, counters zero, pattern -P-, no star.
    KILLS: M-R8 (a reset without a post floor), a rule that does not look at p_min, an event counted with the floor off."""
    kw, P, st0, deps, th, pi0 = k112c("O2", pmin=0.0)
    st = O2.step(P, st0, deps)
    assert not st.flags["refused"] and st.flags["rpattern"] == "-P-"
    assert (st.n_f_tr, st.n_f_post, st.n_f_pi, st.flags["floor_pi"]) == (0, 0, 0, 0)
    assert abs(float(np.trace(st.sigma)) / 3.0 - (-48.764495)) <= 6e-7
    assert abs(st.pi_i - (-40.198569107)) <= 6e-9 and abs(st.dlam - 1.546466e-3) <= 6e-10
    kw1, P1, st01, deps1, _, _ = k112c("O1", pmin=0.0)
    s1 = O1.run_path(P1, st01, deps1[None], rtol=1e-10)[-1]
    assert s1.flags["status"] == "ok" and s1.flags.get("rule_segments", 0) == 0 and not s1.flags.get("floor_rule", False)


def test_k1_12_cases_a_and_b_have_no_post_floor_and_nothing_is_reset():
    """13.12 / 9.7 (round 4: 'Cases (A), (B) do not change (no post floor)'): (A) elastic trial-floored step 'FE-': pi_i unchanged bit for bit; (B) plastic
    from a floored trial, p_c = -50.70 'FP-': no post floor, F = 0 by the return, no reset.  Both: floor_pi = 0, n_f_pi = 0, no star in the pattern.
    The cases are those of test_g1_floor.py (K2 BA06, p_min = 50, alpha0 = 0 and 5).
    KILLS: M-R8 (the rule applied to an elastic step or to a plastic step without a post floor)."""
    from test_g1_floor import case_A, case_B
    for a0 in (0.0, 5.0):
        P, st0, deps = case_A(a0)
        sa = O2.step(P, st0, deps)
        assert sa.flags["rpattern"] == "FE-" and sa.pi_i == st0.pi_i and (sa.n_f_pi, sa.flags["floor_pi"]) == (0, 0)
        P, st0, deps = case_B(a0)
        sb = O2.step(P, st0, deps)
        assert sb.flags["rpattern"] == "FP-" and (sb.n_f_pi, sb.flags["floor_pi"], sb.n_f_post) == (0, 0, 0)


# ======================================================================================================================
# the wet-floor table (9.7 item 3, round 4)
# ======================================================================================================================
WET_TABLE = {   # psi_i -> {dgamma: (F(sigma_f, pi_BE)/p_min, pi_BE, eta_f, pi_i by (S.57))}
    0.0: {1e-4: (0.2699, -0.328889, 0.940109, -0.386215057), 1e-3: (0.3496, -0.329228, 1.021599, -0.406870322),
          1e-2: (0.3572, -0.329262, 1.029388, -0.408940969)},
    0.04: {1e-4: (0.2717, -0.328346, 0.938960, -0.385936497), 1e-3: (0.3670, -0.324406, 1.012743, -0.404537473),
           1e-2: (0.5322, -0.287273, 0.951647, -0.389032266)},
    0.08: {1e-4: (0.2732, -0.327887, 0.937988, -0.385701048), 1e-3: (0.3820, -0.320380, 1.005324, -0.402600177),
           1e-2: (0.6866, -0.258876, 0.897090, -0.376005318)},
    0.12: {1e-4: (0.2745, -0.327494, 0.937157, -0.385499896), 1e-3: (0.3950, -0.316976, 0.999032, -0.400969027),
           1e-2: (0.8217, -0.238442, 0.856623, -0.366805331)},
}


def wet_state(P, kw, psi):
    """9.7 item 3 'Inputs to rebuild': HAR TIMs k, g, n = 1/2, p_a 101; M 1.3309, N 0.4, N_bar 0.2, chi -3.5, h 280; WW rho = rho_bar = 0.71; fork CSL e0 0.83,
    lambda_c 0.027, xi 0.45; cap none; p_min 0.505.  State p = -0.505, eta = 0.5 M, n_hat = dev(-1.2, -0.1, 1.3), pi_i from (S.12), v_n = 1 + e0 - lambda_c
    (-pi_i/p_a)^xi + psi_i."""
    M, N, rho = kw["M"], kw["N"], kw["rho"]
    sig, th, q = C.sigma_on_surface(M, N, rho, -PMIN3, 0.5 * M, C.WET_DIR)
    pi = C.pi_S53(M, N, -PMIN3, 0.5 * M)
    v = C.v_for_psi(kw, pi, psi)
    return O2.initial_state(P, np.diag(sig), v, pi), th, q, pi


def test_wet_floor_table_twelve_printed_cells():
    """9.7 item 3 (round 4) table, 'every case -Pf*, m = 1': start p = -0.505, q = 0.2583 kPa, theta = 0.4544, pi_i = -0.328007 kPa; for psi_i in {0, +0.04,
    +0.08, +0.12} and d gamma in {1e-4, 1e-3, 1e-2} (increment (d gamma/sqrt 2) n_hat, tr = 0, one BE step) the printed F(sigma_f, pi^BE)/p_min, pi_i^BE,
    eta_f and the pi_i by (S.57) to 9 decimals; F(sigma_f, pi_i)/p_min = 0 to 7e-16 after the rule.
    KILLS: M-R4 (pi_i stays pi^BE: F = +0.27 ... +0.82 p_min), M-R7, M-R2 (the zeta_y term does not enter pi_f, only the chain: caught by the digits of pi_f
    through eta_f), pi_f from the wrong (p_f, eta_f), a rule that shrinks."""
    kw = C.har_kw(p_min="default")
    P = C.make_hf("O2", **kw)
    for psi, row in WET_TABLE.items():
        st0, th, q, pi0 = wet_state(P, kw, psi)
        assert abs(q - 0.2583) <= 6e-5 and abs(th - 0.4544) <= 6e-5 and abs(pi0 + 0.328007) <= 6e-7
        assert abs(st0.psi_i - psi) <= 1e-12
        nhat = (C.WET_DIR - C.WET_DIR.mean())
        nhat = nhat / np.linalg.norm(nhat)
        for dg, (Fp, pibe, etaf, pii) in row.items():
            st = O2.step(P, st0, np.diag(dg / math.sqrt(2.0) * nhat))
            assert not st.flags["refused"] and st.flags["substeps"] == 1, (psi, dg, st.flags["reason"])
            assert st.flags["rpattern"] == "-Pf*", (psi, dg, st.flags["rpattern"])
            res = st.cache["res"]
            assert abs(res.F_post / PMIN3 - Fp) <= 6e-5, (psi, dg, res.F_post / PMIN3)
            assert abs(res.pi_be - pibe) <= 6e-7, (psi, dg, res.pi_be)
            assert abs(res.eta_f - etaf) <= 6e-7, (psi, dg, res.eta_f)
            assert abs(st.pi_i - pii) <= 6e-10, (psi, dg, st.pi_i)
            assert abs(F_of(kw, st.sigma, st.pi_i) / PMIN3) <= 1e-12 and abs(st.pi_i) >= c_N(0.4) * PMIN3


# ======================================================================================================================
# K1.16: the H19 regression points
# ======================================================================================================================
@pytest.mark.parametrize("name", ["e7367g1_s185", "e7367g2_s185"])
def test_k1_16_s185_states_printed_rows_and_invariants(name):
    """13.16, K3 set, state `s185`, increment scaled to |d eps| = 2e-4 and repeated: every step m = 1, pattern -Pf* (floor_tr 0, plastic, floor_post 1, rule fired),
    p = -p_min, F <= F_tol and |pi_i| >= c_N p_min after every accepted step (S.58), pi_i = pi_f(sigma_f) to 1e-12 relative (S.57); the printed rows (sigma11,
    sigma22, sigma33, sigma12, pi_i) at steps 1, 2 and 60 to 5e-10, |F| <= 1e-12 (printed <= 6e-15); the start eps^e_33 recovered by the inverse map equals the printed
    value; the start is OUTSIDE F by the printed +0.458 / +0.291 p_min (a legacy state); the step-1 detail (d lambda, p_c, eta_c, pi^BE, F(sigma_f, pi^BE), eta_f,
    d eps^f_v).  Counter: n_f_pi = k after k steps.
    KILLS: M-R4 (round 3b leaves F = +0.176 / +0.118 kPa and pi_i = -0.3245 / -0.3676 at step 1), M-R7, M-R8, M-R2 / M-R5 only through the chain (see the (S.59) test),
    a rule applied at entry to the legacy state, an uncounted reset."""
    X = K116[name]
    P = C.make_hf("O2", **K3)
    st0, d = h19_state("O2", P, h19_records()[name])
    assert abs(st0.eps_e[2, 2] - X["e33"]) <= 1e-12
    assert st0.n_f_pi == 0 and abs(F_of(K3, st0.sigma, st0.pi_i) - X["F0"]) <= 6e-5
    cN = c_N(K3["N"])
    st, prev = st0, st0
    for k in range(1, 61):
        st = O2.step(P, prev, d)
        assert not st.flags["refused"], (k, st.flags["reason"])
        assert st.flags["substeps"] == 1, f"step {k} not at m = 1"
        assert st.flags["rpattern"] == "-Pf*" and (st.flags["floor_tr"], st.flags["floor_post"], st.flags["floor_pi"]) == (0, 1, 1), (k, st.flags["rpattern"])
        assert st.n_f_pi == k and st.n_f_post == k
        assert abs(float(np.trace(st.sigma)) / 3.0 + PMIN3) <= 1e-12 * PMIN3
        F = F_of(K3, st.sigma, st.pi_i)
        assert F <= FTOL3 and abs(st.pi_i) >= cN * PMIN3 * (1.0 - 1e-12), (k, F, st.pi_i)
        assert abs(st.pi_i - pi_f_of(K3, st.sigma)[0]) <= 1e-12 * abs(st.pi_i), k
        assert st.pi_i <= st.cache["res"].pi_be                      # never shrunk (|pi_i| >= |pi^BE|)
        if k in X["rows"]:
            s11, s22, s33, s12, pi = X["rows"][k]
            got = (st.sigma[0, 0], st.sigma[1, 1], st.sigma[2, 2], st.sigma[0, 1], st.pi_i)
            assert np.abs(np.array(got) - np.array([s11, s22, s33, s12, pi])).max() <= ROW_TOL, (k, got)
            assert abs(F) <= F_DIGITS_TOL, (k, F)
        if k == 1:
            res, D = st.cache["res"], X["detail"]
            assert abs(st.dlam - D["dlam"]) <= 6e-12 * 1e3                    # 10 printed digits of 1.4e-4: 5e-14 rounding, solver 1e-12 relative
            pre = st.cache["floor_events"][-1][0]
            ev, es, _ = C.invariants_eps(pre)
            p_c, _ = C.har_pq(ev, es, k=K3["k"], g=K3["g"], n=0.5)
            assert abs(p_c - D["p_c"]) <= 6e-10 and abs(st.eta - D["eta_c"]) <= 6e-7
            assert abs(res.pi_be - D["pi_be"]) <= 6e-10 and abs(res.F_post - D["F_post"]) <= 6e-10
            assert abs(res.eta_f - D["eta_f"]) <= 6e-10 and abs(st.cache["floor_events"][-1][2] - D["dfv"]) <= 6e-12
        prev = st


def test_k1_16_stall_row_by_digits_at_m1_else_by_the_invariants():
    """13.16, `e7367g1_m53` (the stall state, F = +1.6596 kPa = +3.29 p_min, a legacy round-3b state): the first step of 2e-4 is attempted by the ordinary return;
    it converges, the post floor acts and 5g commits F = 0.  The printed rows (steps 1, 2) are asserted by digits ONLY if the step ran at m = 1 (the return from
    3 p_min outside F is far from its start and its path depends on the solver's globalisation); otherwise by the invariants of the sheet: accepted, pattern
    -Pf*, F <= F_tol, pi_i = pi_f(sigma_f) to 1e-12 relative, p = -p_min, |pi_i| >= c_N p_min.  The step-1 detail (d lambda, p_c, q_c, eta_c, pi^BE, F(sigma_f, pi^BE),
    eta_f, d eps^f_v) is a digit assertion under the same condition.  A refusal is NOT admissible at this increment (9.7: 'at |d eps| = 2e-4 it converges').
    KILLS: M-R4 (F = +1.4397 kPa after step 1), the rule applied at entry (dead code: the first step would not be -Pf*), a committed state with F > F_tol."""
    name = "e7367g1_m53"
    X = K116[name]
    P = C.make_hf("O2", **K3)
    st0, d = h19_state("O2", P, h19_records()[name])
    assert abs(st0.eps_e[2, 2] - X["e33"]) <= 1e-12 and abs(F_of(K3, st0.sigma, st0.pi_i) / PMIN3 - X["F0_pmin"]) <= 6e-3
    cN = c_N(K3["N"])
    st = st0
    for k in (1, 2):
        prev = st
        st = O2.step(P, prev, d)
        assert not st.flags["refused"], (k, st.flags["reason"])
        assert st.flags["rpattern"].endswith("-Pf*") and st.flags["floor_post"] >= 1 and st.flags["floor_pi"] >= 1, (k, st.flags["rpattern"])
        assert abs(float(np.trace(st.sigma)) / 3.0 + PMIN3) <= 1e-12 * PMIN3
        assert F_of(K3, st.sigma, st.pi_i) <= FTOL3 and abs(st.pi_i) >= cN * PMIN3 * (1.0 - 1e-12), k
        assert abs(st.pi_i - pi_f_of(K3, st.sigma)[0]) <= 1e-12 * abs(st.pi_i), k
        assert st.n_f_pi == prev.n_f_pi + st.flags["floor_pi"] and st.n_f_pi >= k
        if st.flags["substeps"] == 1:
            s11, s22, s33, s12, pi = X["rows"][k]
            got = (st.sigma[0, 0], st.sigma[1, 1], st.sigma[2, 2], st.sigma[0, 1], st.pi_i)
            assert np.abs(np.array(got) - np.array([s11, s22, s33, s12, pi])).max() <= ROW_TOL, (k, got)
            assert abs(F_of(K3, st.sigma, st.pi_i)) <= F_DIGITS_TOL
            if k == 1:
                res, D = st.cache["res"], X["detail"]
                assert abs(st.dlam - D["dlam"]) <= 6e-9
                assert abs(res.pi_be - D["pi_be"]) <= 6e-10 and abs(res.F_post - D["F_post"]) <= 6e-10 and abs(res.eta_f - D["eta_f"]) <= 6e-10
                assert abs(st.cache["floor_events"][-1][2] - D["dfv"]) <= 6e-7 * D["dfv"]     # printed to 7 digits
                pre = st.cache["floor_events"][-1][0]
                ev, es, _ = C.invariants_eps(pre)
                p_c, q_c = C.har_pq(ev, es, k=K3["k"], g=K3["g"], n=0.5)
                assert abs(p_c - D["p_c"]) <= 6e-10 and abs(q_c - 0.080655979) <= 6e-10 and abs(st.eta - D["eta_c"]) <= 6e-7


def test_k1_16_stall_state_small_increment_is_refused_frozen_or_accepted_with_the_invariants():
    """9.7 'A state the rule never commits': from the stall state at |d eps| = 1e-6 a return from 3 p_min outside F is not continuous with the start; the
    admissible outcomes are (i) a REFUSAL with the state frozen (pi_i = -0.1116066186 unchanged) and NOTHING counted, or (ii) an accepted step that satisfies
    the invariants.  'It must never commit a state with F > F_tol.'
    KILLS: a counter incremented on a refused increment, a refused increment that moves pi_i, a committed state outside F."""
    P = C.make_hf("O2", **K3)
    st0, d = h19_state("O2", P, h19_records()["e7367g1_m53"], norm=1.0e-6)
    st = O2.step(P, st0, d)
    if st.flags["refused"]:
        assert st.pi_i == st0.pi_i and abs(st.pi_i + 0.1116066186) <= 6e-11
        assert (st.n_f_tr, st.n_f_post, st.n_f_pi) == (st0.n_f_tr, st0.n_f_post, st0.n_f_pi)
    else:
        assert F_of(K3, st.sigma, st.pi_i) <= FTOL3 and abs(st.pi_i) >= c_N(K3["N"]) * PMIN3 * (1.0 - 1e-12)
        assert st.n_f_pi == st0.n_f_pi + st.flags["floor_pi"]


def test_k1_16_o1_marches_keep_the_invariants_after_the_first_accepted_step():
    """13.16 / (S.58) on the rate oracle (O1, the rule in its rate form (S.60)): the printed start states are legacy states OUTSIDE F (+0.458 p_min), from which O1 (an
    event integrator) cannot start (status 'outside'; 'the assertions of K1.16 apply after the first accepted step'), so O1 starts from the admissible state after the first
    accepted step of the split (its sigma, v, pi_i: F = 0, p = -p_min) and takes three increments of 2e-4: every step status ok, p = -p_min, |pi_i| >= c_N p_min and F within the
    integration tolerance of the surface (O1 integrates to rtol 1e-10; no printed O1 tolerance: |F| <= 1e-6 p_min).  The 60-step digits are the split's (O2 / kernel), not O1's.
    KILLS: an O1 that runs the Koiter form (S.55) on this path (round-3b collapse), a rate form that leaves the surface."""
    P2 = C.make_hf("O2", **K3)
    st20, d = h19_state("O2", P2, h19_records()["e7367g1_s185"])
    s1 = O2.step(P2, st20, d)
    P = C.make_hf("O1", **K3)
    st = O1.initial_state(P, s1.sigma, s1.v, s1.pi_i)
    cN = c_N(K3["N"])
    for k in range(3):
        st = O1.run_path(P, st, d[None], rtol=1e-10)[-1]
        assert st.flags["status"] == "ok", (k, st.flags["status"])
        assert abs(float(np.trace(st.sigma)) / 3.0 + PMIN3) <= 1e-6 * PMIN3
        assert abs(F_of(K3, st.sigma, st.pi_i)) <= 1e-6 * PMIN3 and abs(st.pi_i) >= cN * PMIN3
    assert st.flags.get("n_rule_increments", 0) >= 1


def test_the_rule_absent_mutant_m_r4_is_far_from_every_printed_row():
    """Power of the gates above (M-R4: the rule absent): with step 5g patched out the very first K1.16 step from `e7367g1_s185` leaves F = +0.1760 kPa and commits
    pi_i = pi^BE = -0.3245384951 (printed, round 3b), 3e-1 away from the printed -0.3989652918; after 60 steps the round-3b state is pi_i = -0.3894275650,
    F = +0.0929 kPa, sigma = (-0.5098218813, -0.6109439554, -0.3942341633; sigma12 -0.2483437655), q = 0.469 against 0.631 kPa with the rule (all printed).
    KILLS: nothing in the code; it proves that test_k1_16_s185_states_printed_rows_and_invariants kills M-R4 and that the printed round-3b figures are the rule-off split's."""
    P = C.make_hf("O2", **K3)
    st0, d = h19_state("O2", P, h19_records()["e7367g1_s185"])
    with rule_off():
        out = march(P, st0, d, 60)
    s1, s60 = out[0], out[-1]
    assert abs(s1.pi_i - (-0.3245384951)) <= 6e-10 and abs(F_of(K3, s1.sigma, s1.pi_i) - 0.1760) <= 6e-5
    assert abs(s1.pi_i - K116["e7367g1_s185"]["rows"][1][4]) >= 6e-2
    assert abs(s60.pi_i - (-0.3894275650)) <= 6e-9 and abs(F_of(K3, s60.sigma, s60.pi_i) - 0.0929) <= 6e-5
    want = np.array([-0.5098218813, -0.6109439554, -0.3942341633])
    assert np.abs(np.diag(s60.sigma) - want).max() <= 6e-9 and abs(s60.sigma[0, 1] - (-0.2483437655)) <= 6e-9
    assert abs(q_of(s60.sigma) - 0.469) <= 6e-4
    with_rule = march(P, st0, d, 60)[-1]
    assert abs(q_of(with_rule.sigma) - 0.631) <= 6e-4


# ======================================================================================================================
# K1.14 / K1.14b: the rule never fires, and "never shrink" (M-R7)
# ======================================================================================================================
def fpf_build():
    """The state of 13.14b (round 4: 'no value changes'), inputs of 9.7 item 3 / 13.14b: HAR TIMs + M 1.3309, rho = rho_bar 0.71, fork CSL, K2 plastic constants,
    p_min 0.505; state p = -0.6, eta = 1.2 M, n_hat = dev(-1, -0.35, 1.35), pi_i from (S.12) (-0.743658), psi_i = -0.10 (v = 1.72704)."""
    kw = C.har_kw(p_min="default")
    M, N, rho = kw["M"], kw["N"], kw["rho"]
    sig, th, q = C.sigma_on_surface(M, N, rho, -0.6, 1.2 * M, C.OFFC_DIR)
    pi = C.pi_S53(M, N, -0.6, 1.2 * M)
    v = C.v_for_psi(kw, pi, -0.10)
    P = C.make_hf("O2", **kw)
    return kw, P, O2.initial_state(P, np.diag(sig), v, pi)


def fpf_inc(tr, sh):
    return np.diag(tr / 3.0 * ONES + sh * C.OFFC_DIR)


def test_k1_14b_the_rule_never_fires_and_the_surface_is_never_shrunk():
    """13.14b (round 4: 'every post floor has F(sigma_f, pi^BE) <= 0, so (S.57) does not fire'): single FPf: F(sigma_f)/p_min = -0.004098, pi_i = -0.743350807,
    sigma_f = (-0.8497174515, -0.5790682544, -0.0862142940), pattern FPf (no star), n_f_pi = 0; chain FPf,FPf (twice the increment, fractions (1/2, 1/2)): F -0.00410, -0.00413
    in the two sub-increments, final pi_i = -0.743038708; chain -P-,-Pf (tr 1.4e-5 + shear 2e-5): F -0.00066 in the second, final pi_i = -0.743573074.
    The always-reset mutant commits -0.740398976 (printed): shown here by patching the rule to ignore the F > 0 test.
    KILLS: M-R7 (the surface shrunk: -0.740398976, 4.0e-3 away), M-R8 (a reset on a post floor with F <= 0), a rule that fires on a dry-side FPf."""
    kw, P, st0 = fpf_build()
    s = O2.step(P, st0, fpf_inc(2e-5, 2e-5))
    assert s.flags["rpattern"] == "FPf" and (s.n_f_pi, s.flags["floor_pi"], s.n_f_post) == (0, 0, 1)
    assert abs(s.cache["res"].F_post / PMIN3 - (-0.004098)) <= 6e-7
    assert abs(s.pi_i - (-0.743350807)) <= 6e-10 and s.pi_i == s.cache["res"].pi_be
    assert np.abs(np.sort(np.diag(s.sigma)) - np.sort(np.array([-0.8497174515, -0.5790682544, -0.0862142940]))).max() <= 6e-10
    c = O2.step_fractions(P, st0, fpf_inc(4e-5, 4e-5), (0.5, 0.5))
    assert c.flags["rpattern"] == "FPf,FPf" and c.n_f_pi == 0 and c.flags["floor_pi"] == 0
    assert abs(c.pi_i - (-0.743038708)) <= 6e-10
    h1 = O2.step_fractions(P, st0, 0.5 * fpf_inc(4e-5, 4e-5), [1.0])
    h2 = O2.step_fractions(P, h1, 0.5 * fpf_inc(4e-5, 4e-5), [1.0])
    assert abs(h1.cache["res"].F_post / PMIN3 - (-0.00410)) <= 6e-6 and abs(h2.cache["res"].F_post / PMIN3 - (-0.00413)) <= 6e-6
    d2 = O2.step_fractions(P, st0, fpf_inc(1.4e-5, 2e-5), (0.5, 0.5))
    assert d2.flags["rpattern"] == "-P-,-Pf" and d2.n_f_pi == 0 and abs(d2.pi_i - (-0.743573074)) <= 6e-10
    g1 = O2.step_fractions(P, st0, 0.5 * fpf_inc(1.4e-5, 2e-5), [1.0])
    g2 = O2.step_fractions(P, g1, 0.5 * fpf_inc(1.4e-5, 2e-5), [1.0])
    assert abs(g2.cache["res"].F_post / PMIN3 - (-0.00066)) <= 6e-6
    # the mutant: reset whenever the post floor acted, whatever the sign of F (M-R7)
    orig = K.floor_pi_rule

    def always_reset(P_, sig_f, ae_f, pi_be):
        inv = K.invariants(sig_f)
        eta_f = K.zeta_y(inv.theta, P_.rho, P_.zeta)[0] * inv.q / abs(inv.p)
        return dict(pi=K.pi_of_eta(P_, inv.p, eta_f), pi_be=pi_be, floor_pi=True, gamma_pi=np.zeros(3), F_post=0.0,
                    floor_pi_band=False)
    K.floor_pi_rule = always_reset
    try:
        m = O2.step(P, st0, fpf_inc(2e-5, 2e-5))
    finally:
        K.floor_pi_rule = orig
    assert abs(m.pi_i - (-0.740398976)) <= 6e-9 and abs(m.pi_i - s.pi_i) >= 2e-3


def test_k1_14_elastic_floor_under_shear_resets_nothing():
    """13.14 (round 4: 14 is the elastic operator): the elastic trial-floored HAR step under shear (M = 50 keeps it elastic) leaves pi_i = -5000 bit for bit, pattern
    'FE-' without a star, n_f_pi = 0.  KILLS: M-R8 (a reset on an elastic step)."""
    kw = C.har_kw(p_min="default", M=50.0)
    P = C.make_hf("O2", **kw)
    st0 = O2.initial_state(P, -100.0 * np.eye(3), 1.70, -5000.0)
    ev0 = C.invariants_eps(st0.eps_e)[0]
    nh = np.sort(C.dir_for_theta(0.5))
    deps = np.diag((1.050e-3 - ev0) / 3.0 * ONES + math.sqrt(1.5) * 2e-4 * nh)
    st = O2.step(P, st0, deps)
    assert st.flags["rpattern"] == "FE-" and st.pi_i == -5000.0 and (st.n_f_pi, st.flags["floor_pi"]) == (0, 0)


# ======================================================================================================================
# (S.58): the seeded random census
# ======================================================================================================================
def census(P, kw, st, rng, n_steps, amp, bias_dir):
    """Random symmetric increments (amplitude `amp`, plus `bias_dir` * amp so that shear and expansion persist); every accepted step is checked against (S.58) and
    the rule's contract; returns the counts."""
    cN = c_N(kw["N"])
    pmin = P.p_min
    stats = dict(acc=0, ref=0, post=0, fired=0, post_inside=0, band=0)
    for _ in range(n_steps):
        a = rng.uniform(-amp, amp, size=(3, 3))
        d = 0.5 * (a + a.T) + amp * bias_dir
        nxt = O2.step(P, st, d)
        if nxt.flags["refused"]:
            stats["ref"] += 1
            assert (nxt.n_f_tr, nxt.n_f_post, nxt.n_f_pi) == (st.n_f_tr, st.n_f_post, st.n_f_pi), "a refused increment counted something"
            assert nxt.pi_i == st.pi_i
            continue
        stats["acc"] += 1
        F = F_of(kw, nxt.sigma, nxt.pi_i)
        p = float(np.trace(nxt.sigma)) / 3.0
        assert -p >= pmin * (1.0 - 1e-12), "(S.58a)"
        assert F <= FTOL3, f"(S.58b): F = {F}"
        assert abs(nxt.pi_i) >= cN * (-p) * (1.0 + DELTA_F * kw["N"] / kw["M"]) ** (-(1.0 - kw["N"]) / kw["N"]) * (1.0 - 1e-9), "(S.58c) tip bound"
        assert abs(nxt.pi_i) >= cN * pmin * (1.0 - 1e-12 - DELTA_F * (1.0 - kw["N"]) / kw["M"]), "(S.58c)"
        assert nxt.n_f_pi == st.n_f_pi + nxt.flags["floor_pi"]
        if nxt.flags["floor_post"] and nxt.flags["plastic"] and nxt.flags["substeps"] == 1:
            res = nxt.cache["res"]
            stats["post"] += 1
            assert nxt.pi_i <= res.pi_be, "the surface was shrunk"
            if res.F_post > FTOL3:
                assert nxt.flags["floor_pi"] == 1 and abs(nxt.pi_i - pi_f_of(kw, nxt.sigma)[0]) <= 1e-12 * abs(nxt.pi_i)
                stats["fired"] += 1
            elif res.F_post < -FTOL3:
                assert nxt.flags["floor_pi"] == 0 and nxt.pi_i == res.pi_be
                stats["post_inside"] += 1
            else:
                stats["band"] += 1
        elif nxt.flags["floor_post"] == 0:
            assert nxt.flags["floor_pi"] == 0                              # no post floor, no reset (M-R8)
        st = nxt
    return stats


def test_s58_random_census_of_post_floor_events_k3_har():
    """(S.58) on a seeded random census (K3 HAR, p_min 0.505) from an admissible start (`e7367g1_s185` after its first accepted step: 'the assertions of K1.16 apply
    after the first accepted step'): at every accepted committed state (a) p <= -p_min (1 - 1e-12), (b) F <= F_tol = 1e-10 |p_a|, (c) |pi_i| >= c_N |p| (1 + delta_F N/M)^(-(1-N)/N)
    and >= c_N p_min (1 - 1e-12 - delta_F (1-N)/M); on a post-floored plastic single-increment step: never shrunk (pi_i <= pi^BE), F(sigma_f, pi^BE) > F_tol => the rule fired and
    pi_i = pi_f(sigma_f) to 1e-12 relative, F(sigma_f, pi^BE) < -F_tol => pi_i = pi^BE; no post floor => no reset (M-R8); a refused increment counts nothing.  The census engages the
    rule at least 25 times and meets at least 5 post floors with F(sigma_f, pi^BE) < -F_tol (no reset).
    KILLS: M-R4 (F > F_tol after a post floor), M-R7 (a shrunk surface), M-R8, a pi_i that collapses with the pre-floor pressure (the H19 failure: |pi_i| below c_N p_min)."""
    kw = K3
    P = C.make_hf("O2", **kw)
    st0, d = h19_state("O2", P, h19_records()["e7367g1_s185"])
    st = O2.step(P, st0, d)
    assert not st.flags["refused"]
    rng = np.random.default_rng(5801)
    bias = np.array([[0.3, -0.3, 0.0], [-0.3, 0.1, 0.0], [0.0, 0.0, 0.0]])
    stats = {}
    for s in census(P, kw, st, rng, 300, 4e-4, bias).items():
        stats[s[0]] = s[1]
    assert stats["fired"] >= 25 and stats["post_inside"] >= 5, stats
    assert stats["acc"] >= 100, stats


def test_s58_random_census_k1_14b_har_dry_side_never_shrinks():
    """(S.58) / M-R7 on a second seeded census: the TIMs HAR dry-side state of 13.14b (where every printed post floor has F(sigma_f, pi^BE) <= 0), random increments around the printed
    expansion + shear: every accepted state satisfies (a)-(c); every post-floored single step with F(sigma_f, pi^BE) < -F_tol keeps pi_i = pi^BE (never shrunk).  The census meets at
    least 10 post-floored plastic steps with F(sigma_f, pi^BE) < -F_tol.
    KILLS: M-R7 (an always-reset rule would move pi_i on every such step), M-R8."""
    kw, P, st0 = fpf_build()
    rng = np.random.default_rng(5802)
    bias = np.diag(0.5 * ONES / 3.0 + 0.5 * C.OFFC_DIR)
    stats = census(P, kw, st0, rng, 300, 3e-5, bias)
    assert stats["post"] >= 10 and stats["post_inside"] >= 10, stats


# ======================================================================================================================
# (S.59): the chain row
# ======================================================================================================================
S59 = dict(pi_start=-0.3989652918, F_post_pmin=(0.162, 0.099), pi=-0.479306004, C_fd_printed=(1.0e-5, 1.0e-7, 1.3e-9))
CHAIN_IJ = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))


def chain_E():
    out = []
    for (k, l) in CHAIN_IJ:
        E = np.zeros((3, 3))
        E[k, l] += 1.0
        E[l, k] += 1.0 if k != l else 0.0
        out.append(E)
    return out


def s59_fixture():
    """9.7 'Gate that makes (S.59) testable': K3 set; start = `e7367g1_s185` after ONE K1.16 step; increment = 2x the K1.16 increment (|d eps| = 4e-4); fractions (1/2, 1/2)."""
    P = C.make_hf("O2", **K3)
    st0, d1 = h19_state("O2", P, h19_records()["e7367g1_s185"])
    st1 = O2.step(P, st0, d1)
    return P, st1, 2.0 * d1, d1, (0.5, 0.5)


def chain_err(P, st1, deps, fr, h, tangent_fn=None):
    """max_J ||C : E_J - FD_J|| / ||C : E_J|| (the measure of 9.6), central FD of the final stress over the six columns, fractions fixed, the SAME reset pattern at every FD
    point and no FD point in the noise band."""
    base = O2.step_fractions(P, st1, deps, fr)
    Cc = O2.tangent(P, base) if tangent_fn is None else tangent_fn(P, base)
    e = 0.0
    for E in chain_E():
        sp = O2.step_fractions(P, st1, deps + h * E, fr)
        sm = O2.step_fractions(P, st1, deps - h * E, fr)
        for s in (sp, sm):
            assert not s.flags["refused"] and s.flags["rpattern"] == base.flags["rpattern"], (s.flags["rpattern"], base.flags["rpattern"])
            assert s.flags["floor_pi_band"] == 0, "an FD point inside the noise band of the strict activation"
        fd = (sp.sigma - sm.sigma) / (2.0 * h)
        col = np.einsum("ijkl,kl->ij", Cc, E)
        e = max(e, np.linalg.norm(col - fd) / np.linalg.norm(col))
    return e


def test_s59_forced_ladder_fixture_values_and_chained_tangent_vs_fd():
    """9.7 gate (round 4b), case (b): start `e7367g1_s185` after ONE K1.16 step (pi_i = -0.3989652918), increment 2x (|d eps| = 4e-4), fractions (1/2, 1/2): pattern -Pf*,-Pf* with the
    rule firing in the NON-FINAL sub-increment (F(sigma_f, pi^BE) = +0.162 and +0.099 p_min, outside the noise band), committed pi_i = -0.479306004, p = -0.505.  The chained tangent
    (S.47) with the (S.59) row against the central FD of the whole increment over the six columns: printed 1.0e-5 / 1.0e-7 / 1.3e-9 at h = 1e-6 / 1e-7 / 1e-8 (O(h^2)); pass <= 3e-7 at
    h = 1e-7 (the sheet's threshold), 3x the printed value at the other two; delta:C = 0 (the last sub-increment is post-floored).
    KILLS: M-R1 / M-R3 (see the next test), M-R2 (zeta_y term dropped from f_a: 0.19), M-R5 (a^e of the pre-floor state in gamma: 0.016), M-R6 (the row built before Phi^post: 0.26),
    a chain that carries S^pi_k through the reset sub-increment."""
    P, st1, deps, d1, fr = s59_fixture()
    assert abs(st1.pi_i - S59["pi_start"]) <= 6e-10
    base = O2.step_fractions(P, st1, deps, fr)
    assert not base.flags["refused"] and base.flags["rpattern"] == "-Pf*,-Pf*" and base.flags["floor_pi"] == 2 and base.flags["floor_pi_band"] == 0
    assert base.n_f_pi == st1.n_f_pi + 2
    h1 = O2.step_fractions(P, st1, fr[0] * deps, [1.0])
    h2 = O2.step_fractions(P, h1, fr[1] * deps, [1.0])
    assert h1.flags["floor_pi"] == 1, "the rule must fire in the NON-FINAL sub-increment"
    assert abs(h1.cache["res"].F_post / PMIN3 - S59["F_post_pmin"][0]) <= 6e-4 and abs(h2.cache["res"].F_post / PMIN3 - S59["F_post_pmin"][1]) <= 6e-4
    assert abs(base.pi_i - S59["pi"]) <= 2e-9 and abs(float(np.trace(base.sigma)) / 3.0 + PMIN3) <= 1e-12 * PMIN3
    assert np.abs(base.sigma - h2.sigma).max() <= 1e-12 and abs(base.pi_i - h2.pi_i) <= 1e-12     # the chain driver = two single steps
    errs = {h: chain_err(P, st1, deps, fr, h) for h in (1e-6, 1e-7, 1e-8)}
    assert errs[1e-6] <= 3.0 * S59["C_fd_printed"][0] and errs[1e-7] <= 3e-7 and errs[1e-8] <= 5e-9, errs
    assert 30.0 <= errs[1e-6] / errs[1e-7] <= 300.0, f"not O(h^2): {errs}"
    Ct = O2.tangent(P, base)
    assert np.abs(np.einsum("iikl->kl", Ct)).max() <= 1e-13 * np.abs(Ct).max()


@pytest.mark.parametrize("mutant,printed,bound", [("M-R1", 0.20, 0.05), ("M-R3", 0.020, 5e-3)])
def test_s59_the_replaced_row_matters_the_chain_without_it_is_far_from_fd(monkeypatch, mutant, printed, bound):
    """The same fixture with kernel.chain_propagate monkeypatched on every sub-increment where the rule fired: M-R1 (the (S.54) S^pi row kept after the reset; printed 0.20, the
    Adversary 0.15) and M-R3 (the (S.59) row ADDED to the (S.54) row; printed 0.020): the chained C against the central FD at h = 1e-7 is >= 0.05 / >= 5e-3, i.e. at least
    1.6e4 x / 1.6e4 x the 3e-7 threshold, while the committed stress and pi_i (the state) are those of the shipped chain: a mutant of the ROW only.
    KILLS: (the gate's power) M-R1, M-R3; the previous test would not otherwise see a missing or additive row."""
    P, st1, deps, d1, fr = s59_fixture()
    good = O2.step_fractions(P, st1, deps, fr)
    orig = K.chain_propagate

    def patched(S_eps, S_pi, cum, alpha, v_new, res, nvec):
        if not res.floor_pi:
            return orig(S_eps, S_pi, cum, alpha, v_new, res, nvec)
        res54 = dataclasses.replace(res, floor_pi=False)
        if mutant == "M-R1":
            return orig(S_eps, S_pi, cum, alpha, v_new, res54, nvec)
        e54 = orig(S_eps, S_pi, cum, alpha, v_new, res54, nvec)
        e59 = orig(S_eps, S_pi, cum, alpha, v_new, res, nvec)
        return e59[0], e54[1] + e59[1], e59[2]
    monkeypatch.setattr(K, "chain_propagate", patched)
    bad = O2.step_fractions(P, st1, deps, fr)
    assert np.array_equal(bad.sigma, good.sigma) and bad.pi_i == good.pi_i
    e = chain_err(P, st1, deps, fr, 1e-7)
    assert e >= bound, (mutant, e, printed)
    assert e <= 3.0 * printed + 0.2, "the mutant is not the sheet's mutant (error far above the printed figure)"


def test_s59_single_step_tangent_is_unchanged_by_the_rule():
    """(S.32f) / M-R9: at m = 1 the tangent of the step where the rule fired is the closed form (S.32f)/(S.33): it equals, bit for bit, the tangent of the same step with
    step 5g patched out (sigma_f does not depend on the replacement); the chain at m = 1 reduces to it; the stress is the same, pi_i differs.
    KILLS: M-R9 (the single-step tangent altered by the rule)."""
    P, st1, deps, d1, fr = s59_fixture()
    s = O2.step_fractions(P, st1, d1, [1.0])
    assert s.flags["rpattern"] == "-Pf*"
    with rule_off():
        s_off = O2.step_fractions(P, st1, d1, [1.0])
    assert np.array_equal(O2.tangent(P, s), O2.tangent(P, s_off))
    assert np.array_equal(s.sigma, s_off.sigma) and s.pi_i != s_off.pi_i
    c33 = O2.tangent_last_substep(P, s)
    assert np.abs(O2.tangent(P, s) - c33).max() <= 1e-9 * np.abs(c33).max()


# ======================================================================================================================
# counters
# ======================================================================================================================
def test_n_f_pi_counts_resets_sums_over_sub_increments_and_nothing_else():
    """9.7 'Counted' (round 4): n_f_pi is the committed total of resets: +1 per reset step (K1.16: n_f_pi = k after k steps), a substepped increment SUMS its sub-increments (the
    forced (1/2, 1/2) increment adds 2, floor_pi = 2), an elastic or non-post-floored step adds 0, a refused increment counts nothing.  floor_pi of a step is its own count, n_f_pi the
    cumulative one; O2's State.floor dict carries n_f_pi last and the earlier keys where they were.
    KILLS: an uncounted reset (M-R2 of the plan: counted-or-not), a count not summed over sub-increments, a counter incremented by an elastic step, a refusal that counts."""
    P, st1, deps, d1, fr = s59_fixture()
    assert st1.n_f_pi == 1
    b = O2.step_fractions(P, st1, deps, fr)
    assert (b.flags["floor_pi"], b.n_f_pi) == (2, 3)
    e = O2.step(P, b, np.zeros((3, 3)))                                    # F <= F_tol: an elastic step, nothing counted
    assert e.flags["rpattern"] == "-E-" and e.n_f_pi == 3 and e.flags["floor_pi"] == 0
    keys = list(st1.floor.keys())
    assert keys[:5] == ["at_floor", "n_f_tr", "n_f_post", "n_f_init", "eps_f_v"] and keys[-1] == "n_f_pi" and st1.floor["n_f_pi"] == 1
    # a refused increment: freeze (wild increment from the stall state)
    P2 = C.make_hf("O2", **K3)
    st0, d = h19_state("O2", P2, h19_records()["e7367g1_m53"], norm=1.0e-6)
    r = O2.step(P2, st0, d)
    if r.flags["refused"]:
        assert r.n_f_pi == st0.n_f_pi == 0 and r.flags["floor_pi"] == 0


# ======================================================================================================================
# O2 -> O1 at first order where the rule fires; the gap to the Koiter form; the round-3b split -> Koiter
# ======================================================================================================================
RATE = {   # the rate table of 9.7 (round 4): printed first and last entries, observed orders, the plateau against (S.55)
    "ba06_k2": dict(split_vs_s60_sigma=(1.0e-2, 4.5e-5), split_vs_s60_pi=(2.2e-2, 9.9e-5), plateau=(5.1e-3, 1.09e-2),
                    r3b_vs_s55_sigma=(6.5e-3, 3.1e-5), r3b_vs_s55_pi=(2.4e-2, 1.1e-4)),
    "har_k3": dict(split_vs_s60_sigma=(4.6e-2, 2.2e-4), split_vs_s60_pi=(5.1e-2, 2.8e-4), plateau=(9.4e-2, 0.18),
                   r3b_vs_s55_sigma=(4.9e-2, 7.4e-4), r3b_vs_s55_pi=(1.2e-2, 1.1e-4)),
}
LADDER = (1, 2, 4, 8, 16, 32, 64, 128, 256)


def rate_setup(case, oname):
    """BA06 K2 set, p_min 50, state on the floor at eta = 0.5 M, shear 2e-3 along n_hat (WET_DIR); HAR K3 set, p_min 0.505, eta = 0.5 M, psi_i +0.02, shear 2e-4 along n_hat.
    (The sheet gives no v for the BA06 row: v = 1.70, as K1.12; and no direction for either: the wet-table direction dev(-1.2, -0.1, 1.3).)"""
    nh = C.WET_DIR - C.WET_DIR.mean()
    nh = nh / np.linalg.norm(nh)
    if case == "ba06_k2":
        kw = C.ba06_kw(p_min=50.0)
        P = C.make_hf(oname, **kw)
        sig, _, _ = C.sigma_on_surface(1.2, 0.4, 0.7, -50.0, 0.6, C.WET_DIR)
        pi0 = C.pi_S53(1.2, 0.4, -50.0, 0.6)
        v, amp = 1.70, 2e-3
    else:
        kw = dict(K3)
        P = C.make_hf(oname, **kw)
        sig, _, _ = C.sigma_on_surface(kw["M"], kw["N"], kw["rho"], -PMIN3, 0.5 * kw["M"], C.WET_DIR)
        pi0 = C.pi_S53(kw["M"], kw["N"], -PMIN3, 0.5 * kw["M"])
        v, amp = C.v_for_psi(kw, pi0, 0.02), 2e-4
    st0 = ORACLES[oname].initial_state(P, np.diag(sig), v, pi0)
    return kw, P, st0, np.diag(amp * nh)


def errs_vs(ref, outs, norm):
    """Sheet 9.7 (round 4c) 'sigma-norm of the table': max_a|sigma_a - sigma_a^ref| / max_a|sigma_a^(S.60)|, EVERY column normalised by the (S.60) end state `norm`;
    the pi_i columns are relative to the reference of the column."""
    sig_e = [float(np.abs(s.sigma - ref.sigma).max() / np.abs(norm.sigma).max()) for s in outs]
    pi_e = [abs(s.pi_i - ref.pi_i) / abs(ref.pi_i) for s in outs]
    return sig_e, pi_e


def ladder(P, st0, d, rule=True):
    outs = []
    for m in LADDER:
        if rule:
            path = O2.run_path(P, st0, np.array([d / m] * m))
        else:
            with rule_off():
                path = O2.run_path(P, st0, np.array([d / m] * m))
        assert not any(s.flags["refused"] for s in path), f"m = {m} refused"
        outs.append(path[-1])
    return outs


def orders(e):
    return [math.log(e[i] / e[i + 1]) / math.log(2.0) for i in range(len(e) - 1)]


@pytest.mark.parametrize("case", ["ba06_k2", "har_k3"])
def test_o2_converges_to_o1_rate_form_s60_at_first_order_and_misses_koiter_by_the_plateau(case):
    """9.7 rate table (round 4): the split WITH the rule converges at first order to (S.60) (O1 with floor_rule=True), NOT to the Koiter form (S.55) (O1 floor_rule=False).  Ladder m = 1 ... 256.
    Printed (sigma-norm max|dsigma|/max|sigma| and pi_i relative): split vs (S.60) 1.0e-2 -> 4.5e-5 and 2.2e-2 -> 9.9e-5 (BA06 K2), 4.6e-2 -> 2.2e-4 and 5.1e-2 -> 2.8e-4 (HAR K3), observed orders
    0.90 ... 1.00 / 0.79 ... 1.00; the same split vs (S.55) plateaus at 5.1e-3, 1.09e-2 (BA06) and 9.4e-2, 0.18 (HAR).  Gates: err strictly decreasing, orders of the last four steps in
    [0.9, 1.1], err(m = 1) and err(m = 256) within 1.5x the printed values, the plateau within 3 % of the printed one at m = 256 (two printed digits; definitions of sheet round 4c: Delta eps = X n_hat, v = 1.70 on the BA06 row, every sigma column normalised by the (S.60) end state).
    KILLS: an O1 that integrates (S.55) where the rule fires (the gap to the split does not close), a split without the rule (it would converge to (S.55): the plateau vanishes), a floor
    applied inside the Newton (order lost)."""
    X = RATE[case]
    kw, P1, s01, d = rate_setup(case, "O1")
    _, P2, s02, _ = rate_setup(case, "O2")
    ref60 = O1.run_path(P1, s01, d[None], rtol=1e-10, floor_rule=True)[-1]
    ref55 = O1.run_path(P1, s01, d[None], rtol=1e-10, floor_rule=False)[-1]
    assert ref60.flags["status"] == "ok" and ref55.flags["status"] == "ok"
    assert ref60.flags["rule_segments"] >= 1 and ref55.flags.get("rule_segments", 0) == 0
    outs = ladder(P2, s02, d)
    se, pe = errs_vs(ref60, outs, ref60)
    print(f"\n[{case}] split vs (S.60): sigma {['%.2e' % e for e in se]}\n  pi {['%.2e' % e for e in pe]}")
    for e in (se, pe):
        assert all(e[i + 1] < e[i] for i in range(len(e) - 1)), e
        assert all(0.9 <= o <= 1.1 for o in orders(e)[-4:]), orders(e)
    assert se[0] <= 1.5 * X["split_vs_s60_sigma"][0] and se[-1] <= 1.5 * X["split_vs_s60_sigma"][1], (se[0], se[-1])
    assert pe[0] <= 1.5 * X["split_vs_s60_pi"][0] and pe[-1] <= 1.5 * X["split_vs_s60_pi"][1], (pe[0], pe[-1])
    se55, pe55 = errs_vs(ref55, outs, ref60)
    print(f"  split vs (S.55): sigma {se55[-1]:.3e} pi {pe55[-1]:.3e}  (plateau printed {X['plateau']})")
    assert abs(se55[-1] / X["plateau"][0] - 1.0) <= 0.03, (se55[-1], X["plateau"][0])
    assert abs(pe55[-1] / X["plateau"][1] - 1.0) <= 0.03, (pe55[-1], X["plateau"][1])
    assert abs(se55[-1] - se55[-2]) <= 0.1 * se55[-1], "the gap to (S.55) is not a plateau"


@pytest.mark.parametrize("case", ["ba06_k2", "har_k3"])
def test_round_3b_split_rule_off_converges_to_koiter_at_first_order(case):
    """9.7 rate table, third column: the round-3b split (the rule patched out) against (S.55) (O1 floor_rule=False) converges at first order: BA06 6.5e-3 -> 3.1e-5 (sigma), 2.4e-2 -> 1.1e-4 (pi_i);
    HAR 4.9e-2 -> 7.4e-4, 1.2e-2 -> 1.1e-4.  Same gates as the previous test (orders in [0.9, 1.1] over the last four, ends within 1.5x).
    KILLS: (the gate's power for the previous test) a rule-off split that does not converge to (S.55); an O1 floor_rule=False that is not the Koiter form."""
    X = RATE[case]
    kw, P1, s01, d = rate_setup(case, "O1")
    _, P2, s02, _ = rate_setup(case, "O2")
    ref55 = O1.run_path(P1, s01, d[None], rtol=1e-10, floor_rule=False)[-1]
    ref60 = O1.run_path(P1, s01, d[None], rtol=1e-10, floor_rule=True)[-1]
    outs = ladder(P2, s02, d, rule=False)
    se, pe = errs_vs(ref55, outs, ref60)
    print(f"\n[{case}] round-3b split vs (S.55): sigma {['%.2e' % e for e in se]}\n  pi {['%.2e' % e for e in pe]}")
    for e in (se, pe):
        assert all(e[i + 1] < e[i] for i in range(len(e) - 1)), e
        assert all(0.9 <= o <= 1.1 for o in orders(e)[-4:]), orders(e)
    assert se[0] <= 1.5 * X["r3b_vs_s55_sigma"][0] and se[-1] <= 1.5 * X["r3b_vs_s55_sigma"][1], (se[0], se[-1])
    assert pe[0] <= 1.5 * X["r3b_vs_s55_pi"][0] and pe[-1] <= 1.5 * X["r3b_vs_s55_pi"][1], (pe[0], pe[-1])


def test_the_model_change_on_the_floor_q_with_the_rule_against_koiter_k3():
    """9.7 'What the continuum model is on the floor with the rule': on the K3 path of the rate table, q = 0.405 kPa with (S.60) against 0.293 kPa with (S.55) (+38 %) at a total shear of
    2e-4, and 0.535 against 0.349 kPa (+53 %) at 1e-3 (O1, floor_rule True / False).  Printed to 3 digits: 0.002 abs, the ratio's gap +38 % / +53 % to 3 points.
    KILLS: an O1 rate form (S.60) that reduces to (S.55) (no gap), a gap that is not the sheet's figure."""
    for amp, qs, qk, gap in ((2e-4, 0.405, 0.293, 0.38), (1e-3, 0.535, 0.349, 0.53)):
        kw, P1, s01, d = rate_setup("har_k3", "O1")
        nh = d / np.linalg.norm(d)
        dd = nh * amp
        a = O1.run_path(P1, s01, dd[None], rtol=1e-10, floor_rule=True)[-1]
        b = O1.run_path(P1, s01, dd[None], rtol=1e-10, floor_rule=False)[-1]
        qa, qb = q_of(a.sigma), q_of(b.sigma)
        print(f"\n[q gap amp={amp}] (S.60) {qa:.4f}  (S.55) {qb:.4f}  gap {qa / qb - 1:.3f}")
        assert abs(qa - qs) <= 2e-3 and abs(qb - qk) <= 2e-3 and abs((qa / qb - 1.0) - gap) <= 0.03, (qa, qb)
