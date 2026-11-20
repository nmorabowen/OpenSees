"""WP-144 P1a (Zone B): parity of the C++ kernel SRC/material/nD/LadrunoNorSandKernel.h with O2.

O2 (o2_algo) is the kernel's contract: the same algorithm, branches and constants. This file
drives the SAME strain increments through O2 (o2_algo.step / o2_algo.tangent) and through the
C++ kernel (ctypes shim ns_shim.cpp, built here with g++) and compares, step by step:
    sigma, eps_e, pi_i, v, eps_p_v, eps_p_s, D (the step's dissipation), the 6x6 tangent,
    and the flags: refused (and its finest reason), plastic, vertex, cap_active, substeps,
    and the local-Newton / nested-solve iteration counts.

ITERATION-COUNT EQUALITY is gated on the FIXED paths of this file only (CASES and FRACTION_CASES
below: equal on every step -- the same algorithm with the same constants; a constant changed
without changing the answer shows only there). It is NOT a gate on randomised paths: the
independent census (scratch_indep/census.py, 2 x 200 random paths) found 10 and 19 steps whose
counts differ where a Newton / scan stopping test sits within round-off of its threshold (a
1e-16 difference flips one comparison; the answers still agree at 1e-10). A fixed path that
starts doing that after a change is a signal to re-check, not noise to tolerate.

TANGENT OF A SUBSTEPPED INCREMENT (sheet 9.6, owner decision 2026-10-01): the CHAINED consistent
tangent (S.46)-(S.47), the exact derivative of the final stress with respect to the TOTAL strain
increment. Gated against O2.tangent (which returns cache['C_chain']) at the same 1e-10 as every
other tangent, with the band rules below applied to EVERY sub-increment of the increment (a
sub-increment in a band puts the whole chained tangent in that band). The smooth-cap AMP_STOP
paths substep 30 increments each, and O2's own last-sub-increment CTO is 0.5-0.9 off the chain
there (asserted, so the comparison is not vacuous). FRACTION_CASES drive detail::step_fractions
against O2 api.step_fractions: uniform m = 2, 4, 8, the recursive-halving shapes
(1/2, 1/4, 1/8, 1/8), (1/4, 1/4, 1/4, 1/8, 1/8), (1/2, 1/4, 1/4), the m = 1 chain,
chain = False (last-sub CTO) and the vertex.

GATE (relative, every quantity, every step, every path): 1e-10.
Argument: the kernel IS O2's algorithm in IEEE double precision, so the two can differ only
by the order of floating-point operations (numpy's BLAS dot/gemm and LAPACK LU / eigh against
hand-written sums, an LU and a cyclic Jacobi eigen-solver; numpy's vectorised power against
libm pow). Each is O(eps) = 1e-16 relative per operation. Those differences pass through
(a) the local 4x4 Newton, whose iterates converge quadratically to a root that both sides
compute to round-off (RES_TOL = 1e-12 bounds only where the iteration STOPS, and the two
iteration sequences coincide to round-off, so they stop at the same iterate), and (b) at
most 40 chained history steps of a dissipative (contractive) map. A propagated O(1e-13)
difference is the expectation; 1e-10 leaves three decades of margin and still catches every
formula or branch error, which shows up at >= 1e-8 (the corner switch of (S.9) is O(1e-8)
relative, a dropped term is O(dlambda) ~ 1e-3). Branch flags and the refusal step must
agree exactly: a decision that differs is a contract failure, not a tolerance issue.
Measure: per step, max|kernel - O2| / max(max|O2|, 1e-6 x natural scale), the floor only
for quantities that pass through zero (natural scales: |p0| for stresses and pi_i,
kappa_hat for strains, 1 for v); the tangent normwise (max |dC| / max |C|). D, the
dissipation INCREMENT of a step, is measured against the path's dissipation scale,
max_n |D_n - D_n(O2)| / max_n |D_n(O2)|: on a near-neutral plastic step (the NEUTRAL path)
dlambda is fixed by F_tr ~ 1e-7 kPa, a difference of O(100 kPa) terms, so the step's own
D is conditioned at ~1e-7 relative in O2 itself.

THE EXCEPTION -- the tangent in three round-off conditioning bands (state quantities stay at
1e-10 inside them too). O2 itself does not determine its tangent to 1e-10 there; the
evidence is O2's own sensitivity to a 1-ulp perturbation of its committed elastic strain,
measured at every band step and printed next to the kernel error ("O2 1-ulp").
 (1) CORNER band, |sin 3theta| < 1e-7 at O2's converged stress (not at the vertex). On an
     exactly axisymmetric path (TXC) theta sits on a corner and is set by the last ulp of
     sqrt(6) y fed to arccos: |sin 3theta| takes the quantised values 0, 1.5e-8, 2.1e-8,
     2.6e-8, ... There zeta_y switches between the corner value (S.9) and the quotient (S.8),
     which agree only to O(|theta - theta_c|) ~ 1e-8, and (S.8) loses ~8 digits to
     cancellation (sheet 3.1). zeta_y multiplies y_ab = O(1/R) in q_ab, so the TANGENT moves
     at 1e-9..1e-8 (O2 1-ulp: up to 1.8e-8) while the state does not (y_a = O(theta -
     theta_c) there, so q_a and the converged state are insensitive).
 (2) NEAR-VERTEX band, R < 1e-3 |p|: n_hat_ab ~ 1/R and y_ab ~ 1/R^2 carry the absolute
     round-off of xi (~eps |p|) as a relative error ~eps |p| / R, amplified through J^-1.
 (3) NEAR-COALESCENT trial eigenvalues, min_{a<b} |eps~_a - eps~_b| < 1e-6 max_a |eps~_a| (any
     sub-increment, elastic or plastic; also the final converged eps^e of the chain assembly):
     the spin terms (sigma_a - sigma_b)/(eps~_a - eps~_b) and the eigenvectors themselves are
     conditioned as |eps~| / gap, so the ulp-level difference between numpy eigh and the kernel's
     Jacobi is amplified to ~eps |eps~| / gap; and below REPEATED_EIG_TOL = 1e-10 the (S.33)
     quotient switches to its limit, the two agreeing only to O(gap). Found by the independent
     census (scratch_indep/, tangent outliers 1e-10..2e-8 relative at trial gaps 1e-18..3e-7).
     An EXACTLY zero gap (axisymmetric / isotropic) is not in it: gated at the full 1e-10.
Inside the three bands the tangent gate is 1e-7: the O(|theta - theta_c|) accuracy of (S.8)/(S.9)
at the corner-band edge. A formula error in the tangent is O(dlambda) ~ 1e-3 and is still
caught; the same tangent code is gated at 1e-10 on every off-band step (generic theta: the
non-coaxial, GA, cap and TXE paths).

On a refused step the kernel returns the COMMITTED state n (API contract: frozen) while O2
returns the trial-elastic state of the whole increment; for refused steps only the decision,
the reason and the frozen state are compared.

Paths: the G1 set (tests/test_g1_convergence_tangents.py, test_g1_cap.py) -- drained /
undrained TXC, TXE, the NONCOAXIAL path; paper + fork CSL; WW + GA; cap none / planar /
smooth including the AMP_STOP smooth-cap path (n = 40, where substepping is normal) and the
planar / no-cap AMP_STOP paths that O2 refuses; the vertex / isotropic paths; plus branch
coverage cases (N = 0, alpha0 != 0 with an on-surface non-hydrostatic start); and the
parameter refusals (validate must refuse the counterexample).

ROUND 3B (energy option, p' floor, unified pi_i0; sheet §2.3-§2.4, §5.4, §9.7, §10.2 (S.56), §13 K1.12-K1.15):
every path above now runs with O2's DEFAULT floor p_min = 5e-3 p_ref (0.5 kPa on K2), which none of them reaches
(the floor counters are compared and stay 0: the defaults leave every result unchanged). New paths: HAR (TIMs set,
n = 1/2, p_a 101) drained / undrained TXC, TXE, non-coaxial, GA and the smooth-cap AMP_STOP path, all at p ~ 100
kPa (the HAR->BA06 mutant dies here, A4); floor paths: K1.12 BA06 (FE-), the §9.7 FD-record cases (A) FE-, (B) FP-,
(C) -Pf under BA06 at p_min = 50 kPa with alpha0 = 0 and 5, K1.13 HAR in-domain and OUT-OF-DOMAIN trials (and the
p_min = 0 refusal 'trial_elastic_domain'), K1.14b HAR dry-side FPf x 4, the HAR wet-side -Pf, the near-floor -P-,
an elastic HAR FE-, and an initial state above the floor (n_f_init = 1); the unified pi_i0 rule (S.53) on a
smooth-cap drained TXC start. Per step the floor counters are compared EXACTLY (floor_tr, floor_post per increment,
n_f_tr, n_f_post, n_f_init, at_floor) and eps_f_v, W_f, E_f (S.52) at the 1e-10 gate; at every step whose last
operator is an active Pi_f the kernel's delta : C must vanish (<= 1e-12 max|C|). FRACTION_CASES add the HAR chains
FPf,FPf and -P-,-Pf (S.54), the m = 1 chain at FPf, and the BA06 m = 2 FP-,FP- chain. Closed-form unit tests:
K1.1h / K1.11 (HAR elastic), K1.12-K1.14 floor closed forms (S.49)/(S.50) incl. general n, Pi_f against O2, the
(S.53) values of K1.15 (unified / legacy / the B guard), p_ref and the p_min default, and the new validate()
refusals (HAR constants, p_min < 0, (S.56) gated to cap = smooth).

Run (Esmeralda, WP-144 venv), from Ladruno_files/testbed/norsand_oracle/:
    python -m pytest kernel_parity -q -p no:cacheprovider -s
"""
from __future__ import annotations

import functools
import math
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))   # norsand_oracle/ (o2_algo)
sys.path.insert(0, HERE)

import o2_algo as O2                     # noqa: E402
from o2_algo import api as O2api         # noqa: E402
from o2_algo import kernel as O2K        # noqa: E402
import ns_kernel as NK                   # noqa: E402

GATE = 1.0e-10
ZERO_FLOOR = 1.0e-6          # x natural scale, for quantities that pass through zero
CORNER_BAND = 1.0e-7         # |sin 3theta| below: round-off corner band
VERTEX_BAND = 1.0e-3         # R/|p| below: near-vertex band
# Near-coalescent band: 0 < min trial-eigenvalue gap / max|eps~| < COALESCENT_BAND (a NON-ZERO gap). An
# EXACTLY zero gap (axisymmetric / isotropic trial: both sides take the same repeated-eigenvalue limit) is NOT
# in the band and is gated at the full GATE: the ~1e-8 effect is the (a,b)/(b,a) row-convention difference
# of the limit, which needs a non-zero gap (the rows agree only for EXACTLY repeated eigenvalues).
COALESCENT_BAND = 1.0e-6
BAND_TANGENT_GATE = 1.0e-7   # tangent gate inside the three bands
LASTSUB_MIN_GAP = 0.1        # O2's last-sub-increment CTO must be at least this far from the chain somewhere
I3 = np.eye(3)
SIG0 = -100.0 * I3

QUANTITIES = ("sigma", "eps_e", "pi_i", "v", "eps_p_v", "eps_p_s", "D", "tangent", "eps_f_v", "W_f", "E_f")
REPORTED = QUANTITIES + ("tangent_band", "tangent_coalescent", "o2_1ulp_band", "tangent_substepped", "dC_floor")
DC_FLOOR_GATE = 1.0e-12      # delta : C / max|C| at a step whose last operator is an active Pi_f (sheet 9.7)

# ----------------------------------------------------------------------------------------------
# kernel build (once per session)
# ----------------------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def kern(tmp_path_factory):
    if NK.compiler() is None:
        pytest.skip("no C++ compiler (g++) on PATH: kernel parity is a Zone B test for Esmeralda/Linux")
    lib = NK.build(str(tmp_path_factory.mktemp("ns_kernel")))
    return NK.Kernel(lib)


# ----------------------------------------------------------------------------------------------
# parameter sets and cases (O2 field names)
# ----------------------------------------------------------------------------------------------
K2 = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0, M=1.2, N=0.4, N_bar=0.2,
          chi=-3.5, h=280.0, lam_tilde=0.0135, v_c0=1.81, csl_mode="paper", zeta="WW", cap="none")
MODES = {
    "paper": dict(rho=0.7, rho_bar=0.8),
    "fork": dict(csl_mode="fork", M=1.3309, N=0.3, N_bar=0.2, rho=0.71, rho_bar=0.75,
                 e0=0.83, lam_c=0.027, xi=0.45, p_a=101.325),
}
# HAR TIMs set (sheet §2.3, §13.14b; O2 selfcheck tims_params): n = 1/2, p_a = 101 (round 3b, A4), fork CSL, the K2
# plastic constants; under HAR the five BA06 values are NOT given (O2 refuses them, sheet §2.4)
TIMS = dict(energy="HAR", k=1889.48104361, g=807.80387674, n_e=0.5, p_a=101.0, M=1.3309, N=0.4, N_bar=0.2,
            rho=0.71, rho_bar=0.71, zeta="WW", chi=-3.5, h=280.0, csl_mode="fork", e0=0.83, lam_c=0.027, xi=0.45,
            cap="none")
GA = dict(zeta="GA", rho=0.8, rho_bar=0.9)
SMOOTH = dict(cap="smooth", c1=0.05, c2=0.15)
PLANAR = dict(cap="planar", c1=0.10, c2=0.10)
DENSE = (-80.0, -0.05)       # (pi_i0, psi_i0): dense start for drained / cap paths
LOOSE = (-105.0, 0.01)       # nearly critical loose start for undrained / non-coaxial paths

E1_NC = np.diag([0.003, 0.003, -0.010])
E2_NC = np.array([[0.0005, 0.006, 0.0], [0.006, 0.0005, 0.004], [0.0, 0.004, -0.002]])
KNOT = 0.6


def noncoax_deps(n):
    def eps(s):
        return E1_NC * min(s, KNOT) / KNOT + E2_NC * max(s - KNOT, 0.0) / (1.0 - KNOT)
    t = np.arange(n + 1) / n
    return [eps(t[k + 1]) - eps(t[k]) for k in range(n)]


def e_iso_total(amp):
    return -0.01 * I3 + amp * np.diag([1.0, 0.0, -1.0])


def case(mode, path, start, n, over=None, sigma0=None, **kw):
    return dict(mode=mode, path=path, start=start, n=n, over=over or {},
                sigma0=SIG0 if sigma0 is None else sigma0, **kw)


CASES = {
    # G1 convergence paths, paper + fork
    "TXC_drained_paper": case("paper", "drained", DENSE, 25, ax=-0.06),
    "TXC_drained_fork": case("fork", "drained", DENSE, 25, ax=-0.06),
    "TXC_undrained_paper": case("paper", "undrained", LOOSE, 25, ax=-0.05),
    "TXC_undrained_fork": case("fork", "undrained", LOOSE, 25, ax=-0.05),
    "TXE_drained_paper": case("paper", "drained", DENSE, 25, ax=+0.06),
    "TXE_drained_fork": case("fork", "drained", DENSE, 25, ax=+0.06),
    "TXE_undrained_paper": case("paper", "undrained", LOOSE, 25, ax=+0.05),
    "NONCOAXIAL_paper": case("paper", "noncoax", LOOSE, 25),
    "NONCOAXIAL_fork": case("fork", "noncoax", LOOSE, 25),
    # Gudehus-Argyris
    "TXC_undrained_paper_GA": case("paper", "undrained", LOOSE, 25, over=GA, ax=-0.05),
    "NONCOAXIAL_paper_GA": case("paper", "noncoax", LOOSE, 25, over=GA),
    "NONCOAXIAL_fork_GA": case("fork", "noncoax", LOOSE, 25, over=dict(GA, rho=0.85, rho_bar=0.9)),
    # caps on the near-isotropic G1.cap path (n = 40)
    "CAP_smooth_AMPSTOP_n40": case("paper", "iso", DENSE, 40, over=SMOOTH, amp=2.0e-3),
    "CAP_smooth_dev1e-4_n40": case("paper", "iso", DENSE, 40, over=SMOOTH, amp=1.0e-4),
    "CAP_smooth_AMPSTOP_n40_fork": case("fork", "iso", DENSE, 40, over=SMOOTH, amp=2.0e-3),
    "CAP_planar_AMPSTOP_n40": case("paper", "iso", DENSE, 40, over=PLANAR, amp=2.0e-3),
    "CAP_none_AMPSTOP_n40": case("paper", "iso", DENSE, 40, amp=2.0e-3),
    "CAP_smooth_NONCOAXIAL": case("paper", "noncoax", LOOSE, 25, over=SMOOTH),
    # vertex rule (sheet 3.2): apex hydrostatic steps; isotropic compression past pi_c
    "VERTEX_apex_hydrostatic": case("paper", "hydro", (None, None), 3, dv=-3.0e-3),
    "VERTEX_iso_compression_n20": case("paper", "iso", DENSE, 20, amp=0.0),
    "VERTEX_iso_compression_smooth_n20": case("paper", "iso", DENSE, 20, over=SMOOTH, amp=0.0),
    # branch coverage: N = 0 (log branches of eta, pi*), alpha0 != 0 with an on-surface start
    "N0_TXC_undrained": case("paper", "undrained", LOOSE, 25, over=dict(N=0.0, N_bar=0.0, rho=0.8, rho_bar=0.8),
                             ax=-0.05),
    "ALPHA0_NONCOAXIAL_onsurface": case("paper", "noncoax", (None, None), 25, over=dict(alpha0=2.0),
                                        sigma0=np.array([[-90.0, 6.0, 0.0], [6.0, -100.0, 3.0], [0.0, 3.0, -125.0]]),
                                        v0=1.75),
    # trial contract (sheet 9.1): a tensor-shear increment h on a coaxial yielded TXC state gives
    # F_tr - F_n = O(h^2) (~1e7 h^2 kPa): h = 1e-8 stays below F_TRIAL_TOL = 1e-10 |p0| (elastic), h = 1e-7
    # is above it (plastic). Both decisions must be O2's.
    "NEUTRAL_shear_on_yielded_TXC": case("paper", "custom", LOOSE, 15, deps=(
        [np.diag([0.001, 0.001, -0.002])] * 10
        + [np.array([[0.0, 1e-8, 0.0], [1e-8, 0.0, 0.0], [0.0, 0.0, 0.0]])] * 2
        + [np.array([[0.0, 1e-7, 0.0], [1e-7, 0.0, 0.0], [0.0, 0.0, 0.0]])]
        + [np.diag([0.001, 0.001, -0.002])] * 2)),
    # ---- round 3b: HAR energy (TIMs set) at p ~ 100 kPa (floor inactive; kills the HAR->BA06 mutant, A4) ----
    "HAR_TXC_drained": case("tims", "drained", DENSE, 25, ax=-0.03),
    "HAR_TXC_undrained": case("tims", "undrained", LOOSE, 25, ax=-0.03),
    "HAR_TXE_drained": case("tims", "drained", DENSE, 25, ax=+0.03),
    "HAR_NONCOAXIAL": case("tims", "noncoax", LOOSE, 25),
    "HAR_NONCOAXIAL_GA": case("tims", "noncoax", LOOSE, 25, over=dict(zeta="GA", rho=0.85, rho_bar=0.9)),
    "HAR_CAP_smooth_AMPSTOP_n40": case("tims", "iso", DENSE, 40, over=SMOOTH, amp=2.0e-3),
    # ---- round 3b: the unified pi_i0 rule (S.53) on a smooth-cap start (pi_i0 = None: ramp_end, K1.15) ----
    "PI0_unified_smooth_TXC_drained": case("paper", "drained", (None, None), 25, over=SMOOTH, ax=-0.06, v0=1.70),
    "PI0_unified_smooth_TXC_undrained_fork": case("fork", "undrained", (None, None), 25, over=SMOOTH, ax=-0.05,
                                                  v0=1.75),
}


def o2_params(c):
    if c["mode"] == "tims":
        kw = dict(TIMS)
    else:
        kw = dict(K2)
        kw.update(MODES[c["mode"]])
    kw.update(c["over"])
    return O2.Params(**kw)


# ---- round 3b floor paths (sheet 9.7, K1.12-K1.14b): built from O2 (the selfcheck constructions) ----------------
def _tims(**kw):
    d = dict(TIMS)
    d.update(kw)
    return O2.Params(**d).validate()


def _k2(**kw):
    d = dict(K2, **MODES["paper"])
    d.update(kw)
    return O2.Params(**d).validate()


def _off_corner(P, p_s, eta_s, direction=(-1.0, -0.35, 1.35)):
    from o2_algo.selfcheck import off_corner_sig
    return off_corner_sig(P, p_s, eta_s, direction)


def _v_for_psi(P, pi, psi):
    from o2_algo.selfcheck import v_for_psi
    return v_for_psi(P, pi, psi)


def _first(P, st, cands, pattern, fractions=None):
    """the first candidate increment whose O2 floor pattern is `pattern` (the selfcheck searches)."""
    for d in cands:
        o = O2.step(P, st, d) if fractions is None else O2.step_fractions(P, st, d, fractions)
        if not o.flags["refused"] and o.flags["fpattern"] == pattern:
            return d
    raise AssertionError(f"no candidate gives the floor pattern {pattern!r}")


def setup_k112():
    """K1.12: BA06 K2 set, default p_min 0.5; a trial at p = -p_min/2 floors elastically (FE-), then a shear step at
    the floor (stays there), then compression back off the floor."""
    P = _k2()
    ev_tr = O2K.floor_ev(P, 0.0)[0] + P.kappa_hat * math.log(2.0)
    shear = 1e-6 * np.array([[0, 1, 0.5], [1, 0, 0], [0.5, 0, 0]])
    d1 = ev_tr / 3 * I3 + shear
    return P, -100.0 * I3, 1.59, -2.0e4, [d1, 2e-5 * np.diag([1.0, -1.0, 0.0]), -0.5 * ev_tr / 3 * I3]


def setup_ba06_50(which, alpha0):
    """§9.7 FD record, BA06 K2 set at p_min = 50 kPa: (A) elastic + trial floor, (B) plastic + trial floor (FP-),
    (C) plastic + post floor (-Pf), alpha0 = 0 / 5 (D12 != 0)."""
    P = _k2(alpha0=alpha0, p_min=50.0)
    if which == "A":
        dA = np.diag([2.0e-3, 1.5e-3, 2.5e-3]) + 2e-4 * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        return P, np.diag([-57.0, -61.0, -64.0]), 1.70, -300.0, [dA, -0.5 * dA]
    if which == "B":
        sig, _, nh = _off_corner(P, -55.0, 1.4)
        pi0 = O2K.pi_of_eta(P, -55.0, 1.4)
        st = O2.initial_state(P, np.diag(sig), 1.70, pi0)
        d = _first(P, st, [np.diag(1.0e-3 * O2K.ONES + s * nh) + 1e-4 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
                           for s in (0.5e-3, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3)], "FP-")
        return P, np.diag(sig), 1.70, pi0, [d, 0.5 * d]
    sig, _, _ = _off_corner(P, -51.0, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
    pi0 = O2K.pi_of_eta(P, -51.0, 0.5 * P.M)
    st = O2.initial_state(P, np.diag(sig), 1.70, pi0)
    d = _first(P, st, [np.diag(a * np.array([0.55, 0.40, -0.95]) - 1e-5) + 1e-4 * np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
                       for a in (2e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1.2e-2)], "-Pf")
    return P, np.diag(sig), 1.70, pi0, [d, 0.25 * d]


def setup_k113(dv, p_min=None):
    """K1.13: HAR TIMs set, isotropic p = -1 kPa, d eps_v = +1e-4 (in domain, floored) or +1.1e-4 (OUT of the domain:
    a floor event, no refusal); p_min = 0: the pre-round-3 refusal 'trial_elastic_domain' (M-F5)."""
    P = _tims() if p_min is None else _tims(p_min=p_min)
    return P, -1.0 * I3, 1.70, -50.0, [dv / 3 * I3, -0.5e-4 / 3 * I3 + 1e-6 * np.diag([1.0, -1.0, 0.0])]


def _fpf_state():
    """K1.14b: TIMs set, p_min 0.505, surface state p = -0.6 kPa, eta = 1.2 M, theta 0.271, psi_i = -0.10."""
    P = _tims(p_min=0.505)
    sig, _, nh = _off_corner(P, -0.6, 1.2 * P.M)
    pi0 = O2K.pi_of_eta(P, -0.6, 1.2 * P.M)
    v0 = _v_for_psi(P, pi0, -0.10)
    return P, sig, nh, pi0, v0


def setup_fpf():
    P, sig, nh, pi0, v0 = _fpf_state()
    dF = np.diag((2e-5 / 3) * O2K.ONES + 2e-5 * nh)
    return P, np.diag(sig), v0, pi0, [dF] * 4


def setup_nearfloor():
    P, sig, nh, pi0, v0 = _fpf_state()
    dN = np.diag((-1.5e-5 / 3) * O2K.ONES + 1.5e-5 * nh) + 2e-6 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    return P, np.diag(sig), v0, pi0, [dN]


def setup_wet():
    """§9.7 item 3 (wet side, HAR): a state on the surface near the floor at eta = 0.5 M (theta 0.454), psi_i = +0.04,
    a pure-shear step along n^: the return relaxes p below the floor, the post floor acts (-Pf)."""
    P = _tims(p_min=0.505)
    sig, _, nh = _off_corner(P, -0.6, 0.5 * P.M, direction=(-1.2, -0.1, 1.3))
    pi0 = O2K.pi_of_eta(P, -0.6, 0.5 * P.M)
    v0 = _v_for_psi(P, pi0, 0.04)
    st = O2.initial_state(P, np.diag(sig), v0, pi0)
    d = _first(P, st, [np.diag(a * nh / math.sqrt(2.0)) for a in (1e-4, 3e-4, 1e-3, 3e-3)], "-Pf")
    return P, np.diag(sig), v0, pi0, [d, 0.5 * d]


def setup_har_fe():
    """an elastic trial-floored HAR step (FE-), O2 selfcheck 9 (c): low eta, pi_i far inside."""
    P = _tims(p_min=0.505)
    es2 = 2e-6
    nh = np.array([1, 0, -1]) / math.sqrt(2)
    eps0 = (O2K.floor_ev(P, es2)[0] - 2e-5) / 3 * O2K.ONES + math.sqrt(1.5) * es2 * nh
    sig0 = np.diag(O2K.elastic(P, eps0).sig)
    return P, sig0, 1.70, -2.0e4, [3e-5 / 3 * I3 + 1e-7 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
                                   -3e-5 / 3 * I3]


def setup_har_init(rule_start):
    """an initial state ABOVE the floor (p = -0.2 kPa): projected and counted (n_f_init = 1). rule_start: pi_i0 = None
    -> (S.53) through the FLOORED p (the apex, cap none), then hydrostatic compression at the vertex (sheet 3.2);
    else pi_i0 = -50 (inside), then compressive + shear elastic steps off the floor."""
    P = _tims(p_min=0.505)
    if rule_start:
        return P, -0.2 * I3, 1.70, None, [-2e-5 / 3 * I3, -2e-5 / 3 * I3]
    d = -2e-5 / 3 * I3 + 1e-6 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    return P, -0.2 * I3, 1.70, -50.0, [d, d]


SETUP_CASES = {
    "FLOOR_BA06_K112": setup_k112,
    "FLOOR_BA06_p50_A_alpha0": lambda: setup_ba06_50("A", 0.0),
    "FLOOR_BA06_p50_A_alpha5": lambda: setup_ba06_50("A", 5.0),
    "FLOOR_BA06_p50_B_alpha0": lambda: setup_ba06_50("B", 0.0),
    "FLOOR_BA06_p50_B_alpha5": lambda: setup_ba06_50("B", 5.0),
    "FLOOR_BA06_p50_C_alpha0": lambda: setup_ba06_50("C", 0.0),
    "FLOOR_BA06_p50_C_alpha5": lambda: setup_ba06_50("C", 5.0),
    "FLOOR_HAR_K113_in_domain": lambda: setup_k113(1.0e-4),
    "FLOOR_HAR_K113_out_of_domain": lambda: setup_k113(1.1e-4),
    "FLOOR_HAR_K113_out_pmin0_refused": lambda: setup_k113(1.1e-4, p_min=0.0),
    "FLOOR_HAR_FPf_x4": setup_fpf,
    "FLOOR_HAR_nearfloor_mPm": setup_nearfloor,
    "FLOOR_HAR_wet_mPf": setup_wet,
    "FLOOR_HAR_FE": setup_har_fe,
    "FLOOR_HAR_init_above_floor": lambda: setup_har_init(False),
    "FLOOR_HAR_init_above_floor_S53": lambda: setup_har_init(True),
}
for _nm in SETUP_CASES:
    CASES[_nm] = dict(setup=_nm)


def v0_of(P, c):
    if "v0" in c:
        return c["v0"]
    pi0, psi0 = c["start"]
    if pi0 is None:
        return 1.7
    if P.csl_mode == "paper":
        return psi0 + P.v_c0 - P.lam_tilde * math.log(-pi0)
    return 1.0 + P.e0 + psi0 - P.lam_c * (-pi0 / P.p_a) ** P.xi


def drained_deps(P, st0, ax, n):
    """O2's drained triaxial driver (api.triaxial, kind='drained'), recording the accepted increments."""
    da = ax / n
    sig_lat0 = st0.sigma[0, 0]
    st, out = st0, []
    for _ in range(n):
        dl = 0.0
        if st.cache:
            Ct = O2.tangent(P, st)
            dl = -(Ct[0, 0, 2, 2] * da) / (Ct[0, 0, 0, 0] + Ct[0, 0, 1, 1])
        ok = False
        for _it in range(O2api.TRIAX_LAT_MAX_ITERS):
            d = np.diag([dl, dl, da])
            trial = O2.step(P, st, d)
            if trial.flags["refused"]:
                break
            res = trial.sigma[0, 0] - sig_lat0
            if abs(res) <= O2api.TRIAX_LAT_TOL_REL * P.p_ref:
                ok = True
                break
            Ct = O2.tangent(P, trial)
            dl -= res / (Ct[0, 0, 0, 0] + Ct[0, 0, 1, 1])
        out.append(d)
        st = trial
        if not ok:
            break
    return out


def increments(P, c, st0):
    n = c["n"]
    if c["path"] == "drained":
        return drained_deps(P, st0, c["ax"], n)
    if c["path"] == "undrained":
        da = c["ax"] / n
        return [np.diag([-0.5 * da, -0.5 * da, da])] * n
    if c["path"] == "noncoax":
        return noncoax_deps(n)
    if c["path"] == "iso":
        return [e_iso_total(c["amp"]) / n] * n
    if c["path"] == "hydro":
        return [c["dv"] / 3.0 * I3] * n
    if c["path"] == "custom":
        return list(c["deps"])
    raise KeyError(c["path"])


# ----------------------------------------------------------------------------------------------
# the comparison
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _o2_setup(name):
    """(P, deck sigma0, v0, pi0, increments) of a CASES entry."""
    c = CASES[name]
    if "setup" in c:
        return SETUP_CASES[c["setup"]]()
    P = o2_params(c)
    v0 = v0_of(P, c)
    pi0 = c["start"][0]
    st0 = O2.initial_state(P, c["sigma0"], v0, pi0)
    return P, c["sigma0"], v0, pi0, increments(P, c, st0)


@functools.lru_cache(maxsize=None)
def _o2_run(name):
    P, sigma0, v0, pi0, deps = _o2_setup(name)
    st0 = O2.initial_state(P, sigma0, v0, pi0)
    sts, tans = [], []
    st = st0
    for d in deps:
        st = O2.step(P, st, d)
        sts.append(st)
        tans.append(O2.tangent(P, st))
        if st.flags["refused"]:
            break
    return P, v0, pi0, st0, deps, sts, tans


def _strain_scale(P):
    """natural elastic-strain scale: kappa_hat (BA06), 1/(k(1-n)) (HAR: the domain edge, 1.06e-3 on TIMs)."""
    return P.kappa_hat if P.energy == "BA06" else 1.0 / (P.k * (1.0 - P.n_e))


def _natural_scales(P):
    e = _strain_scale(P)
    pr = P.p_ref
    return dict(sigma=pr, eps_e=e, pi_i=pr, v=1.0, eps_p_v=e, eps_p_s=e, D=pr * e, eps_f_v=e,
                W_f=max(P.p_min, 1e-300) * e)


def _last_op_floor(o):
    """True when the last operator of the increment is an active Pi_f (pattern of the last sub-increment ends in 'f',
    or is 'FE-'): delta : C = 0 there (sheet 9.7)."""
    last = o.flags.get("fpattern", "").split(",")[-1]
    return last.endswith("f") or last == "FE-"


def _o2_kstate(st):
    """O2 State -> the kernel state[18]."""
    return np.array([*NK.t6(st.eps_e), st.pi_i, st.v, st.v0, st.eps_p_v, st.eps_p_s, st.D, st.eps_f_v, st.W_f,
                     st.n_f_tr, st.n_f_post, st.n_f_init, float(bool(st.flags.get("at_floor", False)))], float)


def _floor_mismatch(k, ns, inf, o, prev):
    """exact floor counters / flags of one accepted increment (sheet 9.7 'Counted'): list of (step, what, O2, kernel)."""
    out = []
    for what, ko, oo in (("floor_tr", inf["floor_tr"], o.flags["floor_tr"]),
                         ("floor_post", inf["floor_post"], o.flags["floor_post"]),
                         ("at_floor", bool(inf["at_floor"]), bool(o.flags["at_floor"])),
                         ("state.at_floor", bool(ns[17]), bool(o.flags["at_floor"])),
                         ("n_f_tr", int(ns[14]), o.n_f_tr), ("n_f_post", int(ns[15]), o.n_f_post),
                         ("n_f_init", int(ns[16]), o.n_f_init)):
        if ko != oo:
            out.append((k, what, oo, ko))
    if inf["W_f"] != ns[13]:
        out.append((k, "info.W_f != state.W_f", ns[13], inf["W_f"]))
    return out


def _rel(a, b, natural):
    """max|a - b| / max(max|b|, ZERO_FLOOR * natural) (see module docstring)."""
    a = np.atleast_1d(np.asarray(a, float))
    b = np.atleast_1d(np.asarray(b, float))
    num = float(np.abs(a - b).max())
    if num == 0.0:
        return 0.0
    return num / max(float(np.abs(b).max()), ZERO_FLOOR * natural)


def _o2_1ulp(P, ost, d, Co):
    """O2's own tangent sensitivity to a 1-ulp perturbation of its committed elastic strain."""
    out = [0.0]
    for f in (1.0 + 2.0 ** -52, 1.0 - 2.0 ** -53):
        s2 = ost.copy()
        s2.eps_e = ost.eps_e * f
        ob = O2.step(P, s2, d)
        if not ob.flags["refused"]:
            out.append(np.abs(NK.c4_to_c6(O2.tangent(P, ob)) - Co).max() / np.abs(Co).max())
    return max(out)


def _o2_subincrements(P, ost, d, fractions):
    """Replay O2's sub-increments (api._step_once, the same calls api._run_fractions makes) and return
    their states (each with cache eps_tr / sig / res and flags.plastic)."""
    out, cur = [], ost
    for a in fractions:
        cur = O2api._step_once(P, cur, a * d)
        out.append(cur)
        if cur.flags["refused"]:
            break
    return out


def _min_rel_gap(w):
    w = np.asarray(w, float)
    gap = min(abs(w[a] - w[b]) for a in range(3) for b in range(a + 1, 3))
    return gap / max(float(np.abs(w).max()), 1e-300)


def _coalescent(w):
    g = _min_rel_gap(w)
    return 0.0 < g < COALESCENT_BAND          # exactly 0: gated at the full GATE (see COALESCENT_BAND)


def _bands(subs):
    """(corner_or_vertex, coalescent) over every sub-increment (see the module docstring)."""
    cv, coal = False, False
    for s in subs:
        coal = coal or _coalescent(s.cache["eps_tr"])
        if s.flags["plastic"]:
            # at the reconstructed sigma tensor (eigvalsh), as the P1a gate always measured it
            inv = O2.kernel.invariants(np.linalg.eigvalsh(0.5 * (s.sigma + s.sigma.T)))
            cv = cv or ((not inv.vertex) and (abs(math.sin(3.0 * inv.theta)) < CORNER_BAND
                                              or inv.R < VERTEX_BAND * abs(inv.p)))
    res = subs[-1].cache.get("res")
    if res is not None:                       # the (S.47) assembly runs on the final converged eps^e
        coal = coal or _coalescent(res.eps_e)
    return cv, coal


def _gate_tangent(errs, rep, et, cv, coal, one_ulp):
    """Route one tangent error to its band (corner/vertex first, then coalescent, else the 1e-10 gate)."""
    if cv or coal:
        rep["band_steps"] += 1
        key = "tangent_band" if cv else "tangent_coalescent"
        rep["coalescent_steps"] = rep.get("coalescent_steps", 0) + (0 if cv else 1)
        errs[key] = max(errs[key], et)
        errs["o2_1ulp_band"] = max(errs["o2_1ulp_band"], one_ulp())
    else:
        errs["tangent"] = max(errs["tangent"], et)


def compare(kern, name):
    P, v0, pi0, st0, deps, sts, tans = _o2_run(name)
    nat = _natural_scales(P)
    rep = dict(name=name, n_o2=len(sts), errs={q: 0.0 for q in REPORTED}, flag_mismatch=[], iter_mismatch=0,
               refusal=None, substepped=0, plastic=0, band_steps=0, lastsub_gap=0.0, floor_mismatch=[],
               floor_events=0, patterns=[])
    errs = rep["errs"]
    # initial state (sigma0 is recovered through the elastic inversion, whose Newton stops at 1e-13 |p|). A deck
    # sigma0 above the floor is passed as given (O2's st0.sigma is the FLOORED one, which would not re-project).
    sig_in = _o2_setup(name)[1] if st0.n_f_init else st0.sigma
    rc, ks, msg = kern.initial_state(P, sig_in, v0, pi0)
    assert rc == 0, f"kernel initialState refused: {rc} {msg}"
    rep["init"] = max(_rel(ks[:6], NK.t6(st0.eps_e), nat["eps_e"]), _rel(ks[6], st0.pi_i, nat["pi_i"]),
                      _rel(ks[7], st0.v, 1.0), _rel(ks[8], v0, 1.0), _rel(ks[12], st0.eps_f_v, nat["eps_f_v"]),
                      _rel(ks[13], st0.W_f, nat["W_f"]))
    if (int(ks[16]), bool(ks[17]), int(ks[14]), int(ks[15])) != (st0.n_f_init, bool(st0.flags["at_floor"]), 0, 0):
        rep["floor_mismatch"].append(("init", "n_f_init/at_floor", (st0.n_f_init, st0.flags["at_floor"]),
                                      (ks[16], ks[17])))
    rep["init_sigma"] = _rel(kern.stress(P, ks), NK.t6(st0.sigma), nat["sigma"])
    # step by step
    D_k, D_o = [], []
    st, ost = ks, st0
    for k, (d, o) in enumerate(zip(deps, sts)):
        r = kern.step(P, st, d)
        inf = r["info"]
        if o.flags["refused"] or inf["refusal"] != "OK":
            assert o.flags["refused"] and inf["refusal"] == "SUBSTEPS_EXHAUSTED", \
                f"{name} step {k}: refusal decision differs: O2 {o.flags['reason']!r}, kernel {inf}"
            f2, s2 = NK.parse_o2_reason(o.flags["reason"])
            assert (inf["finest"], inf["finest_sub"]) == (f2, s2), \
                f"{name} step {k}: refusal reason differs: O2 {o.flags['reason']!r} kernel {inf}"
            assert inf["substeps"] == 256 == o.flags["substeps"]
            assert np.array_equal(r["state"], st), f"{name} step {k}: refused step did not freeze the state"
            assert np.array_equal(r["sigma"], kern.stress(P, st))
            rep["refusal"] = (k, o.flags["reason"])
            break
        for f in ("plastic", "vertex", "cap_active", "substeps"):
            if inf[f] != o.flags[f]:
                rep["flag_mismatch"].append((k, f, o.flags[f], inf[f]))
        if (inf["local_iters"], inf["pi_iters"]) != (o.flags["local_iters"], o.flags["pi_iters"]):
            rep["iter_mismatch"] += 1
        rep["substepped"] += o.flags["substeps"] > 1
        rep["plastic"] += bool(o.flags["plastic"])
        ns = r["state"]
        assert ns[8] == v0, "v0 must be carried unchanged"
        for q, x, y in (("sigma", r["sigma"], NK.t6(o.sigma)), ("eps_e", ns[:6], NK.t6(o.eps_e)),
                        ("pi_i", ns[6], o.pi_i), ("v", ns[7], o.v), ("eps_p_v", ns[9], o.eps_p_v),
                        ("eps_p_s", ns[10], o.eps_p_s), ("eps_f_v", ns[12], o.eps_f_v), ("W_f", ns[13], o.W_f)):
            errs[q] = max(errs[q], _rel(x, y, nat[q]))
        # floor (sheet 9.7): exact counters; the increment's d eps^f_v and E_f (S.52) against O2's floor events
        oprev = ost
        rep["floor_mismatch"] += _floor_mismatch(k, ns, inf, o, oprev)
        rep["floor_events"] += o.flags["floor_tr"] + o.flags["floor_post"]
        rep["patterns"].append(o.flags["fpattern"])
        dfv_o = o.eps_f_v - oprev.eps_f_v
        errs["eps_f_v"] = max(errs["eps_f_v"], _rel(inf["deps_f_v"], dfv_o, nat["eps_f_v"]))
        Ef_o = O2api.floor_energy(P, o)
        errs["E_f"] = max(errs["E_f"], _rel(inf["E_f"], Ef_o, max(P.p_min * abs(dfv_o), 1e-300)))
        if inf["E_f"] > P.p_min * inf["deps_f_v"] * (1.0 + 1e-12) + 1e-300:
            rep["floor_mismatch"].append((k, "E_f > p_min deps_f_v (S.52)", P.p_min * inf["deps_f_v"], inf["E_f"]))
        D_k.append(ns[11])
        D_o.append(o.D)
        Co = NK.c4_to_c6(tans[k])
        et = np.abs(r["C"] - Co).max() / np.abs(Co).max()
        if _last_op_floor(o):       # delta : C = 0 when the last operator is an active Pi_f (no bulk stiffness)
            errs["dC_floor"] = max(errs["dC_floor"], np.abs(r["C"][:3, :].sum(axis=0)).max() / np.abs(r["C"]).max())
        m = o.flags["substeps"]
        subs = [o] if m == 1 else _o2_subincrements(P, ost, d, [1.0 / m] * m)
        if m > 1:
            assert np.array_equal(subs[-1].sigma, o.sigma), f"{name} step {k}: O2 sub-increment replay drifted"
            Cl = NK.c4_to_c6(O2.tangent_last_substep(P, o))
            rep["lastsub_gap"] = max(rep["lastsub_gap"], np.abs(Cl - Co).max() / np.abs(Co).max())
        cv, coal = _bands(subs)
        if m > 1 and not (cv or coal):
            errs["tangent_substepped"] = max(errs["tangent_substepped"], et)
        _gate_tangent(errs, rep, et, cv, coal, lambda: _o2_1ulp(P, ost, d, Co))
        st, ost = ns, o
    if D_o and max(abs(x) for x in D_o) > 0.0:
        errs["D"] = max(abs(x - y) for x, y in zip(D_k, D_o)) / max(abs(x) for x in D_o)
    else:
        errs["D"] = max([abs(x) for x in D_k] + [0.0])
    return rep


_REPORTS = {}


def report(kern, name):
    if name not in _REPORTS:
        _REPORTS[name] = compare(kern, name)
    return _REPORTS[name]


# ----------------------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(CASES))
def test_kernel_matches_o2_on_path(kern, name):
    rep = report(kern, name)
    errs = rep["errs"]
    line = "  ".join(f"{q}={errs[q]:.2e}" for q in REPORTED)
    print(f"\n[{name}] steps {rep['n_o2']} plastic {rep['plastic']} substepped {rep['substepped']} "
          f"band {rep['band_steps']} (coalescent only {rep.get('coalescent_steps', 0)}) refusal {rep['refusal']} "
          f"iter-count mismatches {rep['iter_mismatch']}"
          f"\n  init={rep['init']:.2e} init_sigma={rep['init_sigma']:.2e}  {line}"
          + (f"\n  floor events {rep['floor_events']}, patterns {rep['patterns']}" if rep["floor_events"] else ""))
    assert not rep["flag_mismatch"], f"branch flags differ (step, flag, O2, kernel): {rep['flag_mismatch']}"
    assert not rep["floor_mismatch"], f"floor counters differ (step, what, O2, kernel): {rep['floor_mismatch']}"
    assert errs["dC_floor"] <= DC_FLOOR_GATE, f"delta : C = {errs['dC_floor']:.3e} at a floored step (must vanish)"
    # same algorithm, same constants => the same number of local Newton iterations and of nested
    # r(pi_i) evaluations on every step (a constant changed without changing the answer shows here)
    assert rep["iter_mismatch"] == 0, f"{rep['iter_mismatch']} steps with different local/nested iteration counts"
    assert rep["init"] <= GATE, f"initial state differs: {rep['init']:.3e}"
    assert rep["init_sigma"] <= 1e-12, f"stress(initial state) != sigma0: {rep['init_sigma']:.3e}"
    bad = {q: errs[q] for q in QUANTITIES if not errs[q] <= GATE}
    assert not bad, f"parity gate {GATE:.0e} exceeded: {bad}"
    for key in ("tangent_band", "tangent_coalescent"):
        assert errs[key] <= BAND_TANGENT_GATE, \
            f"{key} {errs[key]:.3e} > {BAND_TANGENT_GATE:.0e} (O2 1-ulp {errs['o2_1ulp_band']:.3e})"


def test_v_law_is_exponential(kern):
    """G2 owner decision 2026-10-01 (sheet 1.2 (S.26)): v_{n+1} = v_n exp(tr deps), so v = v0 exp(tr eps) along
    every path, in the kernel (checked here directly) and in O2 (the per-step v parity above). Non-vacuity: on
    the fixed paths the superseded linear law v0 (1 + tr eps) differs from O2's v by far more than the 1e-10
    parity gate, so the v comparison discriminates the two laws (mutate_kernel.sh linear_v_update)."""
    worst_step = worst_path = lin_gap = 0.0
    for name in CASES:
        P, v0, pi0, st0, deps, sts, tans = _o2_run(name)
        rc, st, msg = kern.initial_state(P, st0.sigma, v0, pi0)
        assert rc == 0, msg
        trsum = 0.0
        for d, o in zip(deps, sts):
            if o.flags["refused"]:
                break
            r = kern.step(P, st, d)
            assert r["info"]["refusal"] == "OK", name
            tr = float(np.trace(np.asarray(d, float)))
            trsum += tr
            v_n, v_np1 = st[7], r["state"][7]
            worst_step = max(worst_step, abs(v_np1 - v_n * math.exp(tr)) / v_n)
            worst_path = max(worst_path, abs(v_np1 - v0 * math.exp(trsum)) / v0)
            lin_gap = max(lin_gap, abs(o.v - v0 * (1.0 + trsum)) / v0)
            st = r["state"]
    print(f"\nv-law: kernel per-step |v - v_n exp(tr)|/v {worst_step:.3e}, path |v - v0 exp(tr eps)|/v0 "
          f"{worst_path:.3e}; O2 v vs linear law v0 (1 + tr eps): {lin_gap:.3e}")
    assert worst_step <= 1e-14 and worst_path <= 1e-13
    assert lin_gap > 1e3 * GATE, "the fixed paths do not discriminate the exponential from the linear v-law"


EXPECTED_REFUSALS = {          # sheet 3.2 / 10.1 / 16.6: O2 refuses these after 2^8 substeps
    "CAP_planar_AMPSTOP_n40": True,
    "CAP_none_AMPSTOP_n40": True,
    # sheet 9.7 / 2.4: with p_min = 0 a HAR trial outside dom Psi is the refusal 'trial_elastic_domain' (M-F5) ...
    "FLOOR_HAR_K113_out_pmin0_refused": True,
    # ... and with the floor on it is a floor event, never a refusal (M-F1)
    "FLOOR_HAR_K113_out_of_domain": False,
    "FLOOR_HAR_K113_in_domain": False,
    "FLOOR_HAR_init_above_floor": False,
    "FLOOR_HAR_init_above_floor_S53": False,
    "PI0_unified_smooth_TXC_drained": False,
    "PI0_unified_smooth_TXC_undrained_fork": False,
}

# the floor patterns each floor path must exercise (sheet 9.7, K1.12-K1.14b; O2's fpattern of the increment)
EXPECTED_PATTERNS = {
    "FLOOR_BA06_K112": "FE-",
    "FLOOR_BA06_p50_A_alpha0": "FE-", "FLOOR_BA06_p50_A_alpha5": "FE-",
    "FLOOR_BA06_p50_B_alpha0": "FP-", "FLOOR_BA06_p50_B_alpha5": "FP-",
    "FLOOR_BA06_p50_C_alpha0": "-Pf", "FLOOR_BA06_p50_C_alpha5": "-Pf",
    "FLOOR_HAR_K113_in_domain": "FE-", "FLOOR_HAR_K113_out_of_domain": "FE-",
    "FLOOR_HAR_FPf_x4": "FPf", "FLOOR_HAR_nearfloor_mPm": "-P-", "FLOOR_HAR_wet_mPf": "-Pf", "FLOOR_HAR_FE": "FE-",
}


def test_refusal_reason_elastic_domain(kern):
    """M-F5 / §2.4: the p_min = 0 HAR out-of-domain trial is refused with the finest reason trial_elastic_domain."""
    rep = report(kern, "FLOOR_HAR_K113_out_pmin0_refused")
    assert rep["refusal"] is not None and "trial_elastic_domain" in rep["refusal"][1], rep["refusal"]


@pytest.mark.parametrize("name", list(EXPECTED_PATTERNS))
def test_floor_path_exercises_its_pattern(kern, name):
    """non-vacuity: each floor path really produces its pattern on the first increment (the parity above then
    compares that pattern's numbers, counters and tangent), and the FPf path stays FPf for 4 increments (A1)."""
    rep = report(kern, name)
    assert rep["patterns"] and rep["patterns"][0] == EXPECTED_PATTERNS[name], (name, rep["patterns"])
    if name == "FLOOR_HAR_FPf_x4":
        assert rep["patterns"] == ["FPf"] * 4, rep["patterns"]


def test_refusing_paths_refuse_in_both(kern):
    """The planar / no-cap AMP_STOP paths are refused by O2 (sheet 3.2, 10.1); the kernel must refuse at the
    same step with the same finest reason (asserted inside compare) -- here: that the refusal is actually
    exercised, so the comparison above is not vacuous."""
    for name, want in EXPECTED_REFUSALS.items():
        rep = report(kern, name)
        assert (rep["refusal"] is not None) == want, f"{name}: refusal {rep['refusal']}"


def test_smooth_cap_path_exercises_substepping(kern):
    """O2 README: the AMP_STOP smooth-cap path at n = 40 substeps 30 of its 34 plastic increments. The kernel
    matches the substep counts exactly (flags compared per step) -- check the path really substeps."""
    for name in ("CAP_smooth_AMPSTOP_n40", "CAP_smooth_AMPSTOP_n40_fork"):
        rep = report(kern, name)
        assert rep["refusal"] is None
        assert rep["substepped"] >= 20 and rep["plastic"] >= 1, (name, rep["substepped"])
        # the chained tangent is gated (off-band at 1e-10) and O2's last-sub-increment CTO is far from it:
        # a kernel returning the last-sub CTO cannot pass (mutate_kernel.sh "last_substep_tangent")
        assert rep["errs"]["tangent_substepped"] > 0.0, "no off-band substepped tangent was compared"
        assert rep["lastsub_gap"] > LASTSUB_MIN_GAP, (name, rep["lastsub_gap"])
        print(f"\n[{name}] substepped {rep['substepped']}: chained tangent vs O2 {rep['errs']['tangent_substepped']:.2e}"
              f" (off-band); O2 last-sub CTO vs O2 chain, max {rep['lastsub_gap']:.2f}")


# ----------------------------------------------------------------------------------------------
# detail::step_fractions vs O2 api.step_fractions (sheet 9.6: arbitrary fractions, the m = 1 chain)
# ----------------------------------------------------------------------------------------------
SHEAR3 = np.array([[0.0, 3e-4, 1e-4], [3e-4, 0.0, 2e-4], [1e-4, 2e-4, 0.0]])


def _fork_b():
    """O2 selfcheck fork_params(): fork CSL, WW, rho = rho_bar = 0.71, no cap."""
    return O2.Params(**dict(K2, csl_mode="fork", M=1.3309, N=0.3, N_bar=0.2, rho=0.71, rho_bar=0.71,
                            e0=0.83, lam_c=0.027, xi=0.45, p_a=101.325)).validate()


@functools.lru_cache(maxsize=None)
def _frac_setup(which):
    """(P, O2 state entering the increment, increment) of the O2 selfcheck chain group."""
    if which == "generic":          # (B): fork WW, all three shears, 5 pre-steps
        P = _fork_b()
        s0 = O2.initial_state(P, SIG0, 1.65, -60.4)
        st = O2.run_path(P, s0, np.array([np.diag([4e-4, -1e-3, 0.0]) + SHEAR3] * 5))[-1]
        return P, st, np.diag([1e-4, -6e-4, 2e-4]) + 0.5 * SHEAR3
    if which.startswith("amp"):     # (A)/(E): AMP_STOP smooth cap n = 40, the state entering step n (1-based)
        P = O2.Params(**dict(K2, **MODES["paper"], **SMOOTH)).validate()
        v0 = -0.05 + P.v_c0 - P.lam_tilde * math.log(80.0)
        s0 = O2.initial_state(P, SIG0, v0, -80.0)
        d = (-0.01 * I3 + 2e-3 * np.diag([1.0, 0.0, -1.0])) / 40
        n = int(which[3:])
        st = O2.run_path(P, s0, np.array([d] * (n - 1)))[-1] if n > 1 else s0
        return P, st, d
    if which == "vertex":           # (F): hydrostatic plastic step from the apex, no cap
        P = O2.Params(**dict(K2, **MODES["paper"])).validate()
        return P, O2.initial_state(P, SIG0, 1.59, None), -1e-3 * I3
    # ---- round 3b (sheet 9.7 (S.54), K1.14b): the floored chains under HAR and BA06 ----
    if which.startswith("har_"):
        P, sig, nh, pi0, v0 = _fpf_state()
        st = O2.initial_state(P, np.diag(sig), v0, pi0)
        dF = np.diag((2e-5 / 3) * O2K.ONES + 2e-5 * nh)
        if which == "har_fpf":          # the K1.14b increment itself (FPf; halved: -P-,FPf)
            return P, st, dF
        if which == "har_2fpf":         # twice the increment: FPf,FPf when halved
            return P, st, 2.0 * dF
        if which == "har_mP_mPf":       # -P-,-Pf (the O2 selfcheck search, sheet K1.14b)
            cands = [np.diag((av / 3) * O2K.ONES + a_s * nh) for av in np.linspace(1.0e-5, 2.2e-5, 13)
                     for a_s in np.linspace(1.5e-5, 9.0e-5, 16)]
            return P, st, _first(P, st, cands, "-P-,-Pf", fractions=[0.5, 0.5])
        if which == "har_nearfloor":    # -P- next to the floor, with a shear (generic eigen-data)
            return P, st, np.diag((-1.5e-5 / 3) * O2K.ONES + 1.5e-5 * nh) + 2e-6 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    if which.startswith("ba06_p50_"):   # §9.7 (D): the (B) / (C) increments at p_min = 50 kPa, alpha0 = 0 / 5
        kind, alpha0 = which[len("ba06_p50_")], float(which[-1])
        P, sig0, v0, pi0, deps = setup_ba06_50(kind, alpha0)
        return P, O2.initial_state(P, sig0, v0, pi0), deps[0]
    raise KeyError(which)


FRACTION_CASES = {
    "generic_m8": ("generic", (0.125,) * 8, True),
    "generic_m2": ("generic", (0.5, 0.5), True),
    "generic_halving_2_4_8_8": ("generic", (0.5, 0.25, 0.125, 0.125), True),
    "generic_m1_chain": ("generic", (1.0,), True),
    "generic_m4_nochain": ("generic", (0.25,) * 4, False),
    "amp20_halving_4_4_4_8_8": ("amp20", (0.25, 0.25, 0.25, 0.125, 0.125), True),
    "amp11_halving_2_4_4": ("amp11", (0.5, 0.25, 0.25), True),
    "amp11_m2": ("amp11", (0.5, 0.5), True),
    "amp30_m4": ("amp30", (0.25,) * 4, True),
    "vertex_m1_chain": ("vertex", (1.0,), True),
    "vertex_halving_2_4_4": ("vertex", (0.5, 0.25, 0.25), True),
    # round 3b, sheet 9.7 (S.54) / K1.14b: floored chains (HAR TIMs set, p_min 0.505; BA06 at p_min = 50 kPa)
    "har_FPf_FPf_m2": ("har_2fpf", (0.5, 0.5), True),
    "har_mP_mPf_m2": ("har_mP_mPf", (0.5, 0.5), True),
    "har_FPf_m1_chain": ("har_fpf", (1.0,), True),
    "har_FPf_halving_2_4_4": ("har_2fpf", (0.5, 0.25, 0.25), True),
    "har_FPf_m4_nochain": ("har_2fpf", (0.25,) * 4, False),
    "har_nearfloor_m2": ("har_nearfloor", (0.5, 0.5), True),
    "ba06_p50_FP_FP_m2_alpha0": ("ba06_p50_B0", (0.5, 0.5), True),
    "ba06_p50_FP_FP_m2_alpha5": ("ba06_p50_B5", (0.5, 0.5), True),
    "ba06_p50_mPf_m4_alpha5": ("ba06_p50_C5", (0.25,) * 4, True),
}
FRACTION_PATTERNS = {"har_FPf_FPf_m2": "FPf,FPf", "har_mP_mPf_m2": "-P-,-Pf", "har_FPf_m1_chain": "FPf"}


def _o2_to_kstate(st):
    return _o2_kstate(st)


@pytest.mark.parametrize("label", list(FRACTION_CASES))
def test_kernel_step_fractions_matches_o2(kern, label):
    which, fr, chain = FRACTION_CASES[label]
    P, ost, d = _frac_setup(which)
    o = O2.step_fractions(P, ost, d, list(fr), chain=chain)
    assert not o.flags["refused"], o.flags["reason"]
    r = kern.step_fractions(P, _o2_to_kstate(ost), d, fr, chain=chain)
    inf = r["info"]
    nat = _natural_scales(P)
    assert inf["refusal"] == "OK", inf
    for f in ("plastic", "vertex", "cap_active", "substeps"):
        assert inf[f] == o.flags[f], (f, inf[f], o.flags[f])
    assert (inf["local_iters"], inf["pi_iters"]) == (o.flags["local_iters"], o.flags["pi_iters"]), \
        f"iteration counts differ: kernel {(inf['local_iters'], inf['pi_iters'])} O2 " \
        f"{(o.flags['local_iters'], o.flags['pi_iters'])}"
    assert ("C_chain" in o.cache) == chain
    ns = r["state"]
    errs = {q: _rel(x, y, nat[q]) for q, x, y in (
        ("sigma", r["sigma"], NK.t6(o.sigma)), ("eps_e", ns[:6], NK.t6(o.eps_e)), ("pi_i", ns[6], o.pi_i),
        ("v", ns[7], o.v), ("eps_p_v", ns[9], o.eps_p_v), ("eps_p_s", ns[10], o.eps_p_s),
        ("eps_f_v", ns[12], o.eps_f_v), ("W_f", ns[13], o.W_f))}
    errs["D"] = abs(ns[11] - o.D) / max(abs(o.D), ZERO_FLOOR * nat["D"])
    dfv_o = o.eps_f_v - ost.eps_f_v
    errs["E_f"] = _rel(inf["E_f"], O2api.floor_energy(P, o), max(P.p_min * abs(dfv_o), 1e-300))
    fm = _floor_mismatch(0, ns, inf, o, ost)
    assert not fm, f"floor counters differ: {fm}"
    if label in FRACTION_PATTERNS:
        assert o.flags["fpattern"] == FRACTION_PATTERNS[label], (label, o.flags["fpattern"])
    Co = NK.c4_to_c6(O2.tangent(P, o))
    et = np.abs(r["C"] - Co).max() / np.abs(Co).max()
    cv, coal = _bands(_o2_subincrements(P, ost, d, fr))
    gate = BAND_TANGENT_GATE if (cv or coal) else GATE
    Cl = NK.c4_to_c6(O2.tangent_last_substep(P, o))
    lastsub = np.abs(Cl - Co).max() / np.abs(Co).max()
    dC = np.abs(r["C"][:3, :].sum(axis=0)).max() / np.abs(r["C"]).max()
    print(f"\n[{label}] pattern {o.flags.get('pattern')} fpattern {o.flags.get('fpattern')} tangent {et:.2e} (gate {gate:.0e}"
          f"{', band' if gate > GATE else ''})  O2 last-sub CTO vs O2 tangent {lastsub:.2e}  delta:C {dC:.1e}  "
          + "  ".join(f"{q}={v:.1e}" for q, v in errs.items()))
    bad = {q: v for q, v in errs.items() if not v <= GATE}
    assert not bad, f"parity gate {GATE:.0e} exceeded: {bad}"
    assert et <= gate, f"tangent {et:.3e} > {gate:.0e}"
    if _last_op_floor(o):
        assert dC <= DC_FLOOR_GATE, f"delta : C = {dC:.3e} at a floored increment"


# validate(): owner decisions in force
REFUSE = {
    "counterexample_Nbar_eq_N": dict(MODES["paper"], N_bar=0.4),       # rho < rho_bar, yet beta = 1 > rho/rho_bar
    "Nbar_gt_N": dict(MODES["paper"], N_bar=0.5),
    "WW_rho_half": dict(MODES["paper"], rho=0.5),
    "WW_rhobar_half": dict(MODES["paper"], rho=0.5, rho_bar=0.5),
    "GA_rho_below_7_9": dict(MODES["paper"], zeta="GA", rho=0.75, rho_bar=0.9),
    "GA_rhobar_above_1": dict(MODES["paper"], zeta="GA", rho=0.8, rho_bar=1.01),
    "planar_c1_ne_c2": dict(MODES["paper"], cap="planar", c1=0.05, c2=0.15),
    "p0_positive": dict(MODES["paper"], p0=100.0),
}


@pytest.mark.parametrize("label", list(REFUSE))
def test_validate_refuses_like_o2(kern, label):
    kw = dict(K2)
    kw.update(REFUSE[label])
    P = O2.Params(**kw)          # the dataclass constructor does not validate
    with pytest.raises(ValueError):
        P.validate()
    rc, msg, _ = kern.validate(P)
    assert rc != 0, f"kernel accepted a parameter set O2 refuses ({label})"
    print(f"\n[{label}] kernel rc {rc}: {msg}")


def test_validate_warns_only_for_rho_gt_rhobar(kern):
    kw = dict(K2, rho=0.9, rho_bar=0.8)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        P = O2.Params(**kw).validate()
    assert any("rho_bar" in str(x.message) for x in w)
    rc, msg, warn = kern.validate(P)
    assert rc == 0 and warn, (rc, msg, warn)


# ----------------------------------------------------------------------------------------------
# round 3b: closed forms of the kernel's energy / floor / pi_i0 pieces (sheet §13 K1.1h, K1.11-K1.15) and validate()
# ----------------------------------------------------------------------------------------------
def test_har_elastic_closed_forms(kern):
    """K1.1h / K1.11 (sheet 13, TIMs set, p_a 101): the HAR law in the kernel against the sheet's values (kills
    'HAR silently replaced by BA06', A4: an FD test cannot), and against O2's elastic (p, q, D, a^e of (S.3) in FULL)
    on random states; the domain edge eps* <= 0 is EE_ELASTIC_DOMAIN."""
    P = _tims()
    rc, el = kern.elastic(P, np.full(3, -1e-3 / 3))
    assert rc == 0 and abs(el["p"] + 381.983585) < 1e-6, el["p"]                 # K1.1h
    rc, el = kern.elastic(P, np.full(3, 5e-4 / 3))
    assert rc == 0 and abs(el["p"] + 28.117707) < 1e-6, el["p"]
    nh = np.array([1.0, 0.0, -1.0]) / math.sqrt(2.0)
    rc, el = kern.elastic(P, math.sqrt(1.5) * 8.66546968e-4 * nh)                  # K1.11: eta = 3 g eps_s = 2.1
    assert rc == 0 and abs(el["p"] + 166.548671) < 1e-6 and abs(el["q"] - 349.752210) < 1e-6, el
    rc, el = kern.elastic(P, math.sqrt(1.5) * 1e-3 * nh)
    assert abs(el["q"] / abs(el["p"]) - 2.42341163) < 1e-8
    edge = 1.0 / (P.k * (1.0 - P.n_e))
    assert abs(edge - 1.058491699e-3) < 1e-12
    assert kern.elastic(P, np.full(3, 1.0001 * edge / 3))[0] == NK.EVALERR.index("elastic_domain")
    rng = np.random.default_rng(144)
    worst = 0.0
    for n_e in (0.5, 0.3, 0.0, 0.7):
        Pn = _tims(n_e=n_e)
        e_edge = 1.0 / (Pn.k * (1.0 - n_e))
        for _ in range(20):
            eps = rng.uniform(-3e-3, 0.9 * e_edge / 3, 3) + rng.normal(0, 3e-4, 3)
            if eps.sum() >= 0.95 * e_edge:
                continue
            rc, el = kern.elastic(Pn, eps)
            o = O2K.elastic(Pn, eps)
            assert rc == 0
            worst = max(worst, _rel(el["sig"], o.sig, Pn.p_ref), _rel(el["ae"], o.ae, 1.0),
                        _rel([el["p"], el["q"]], [o.p, o.q], Pn.p_ref),
                        _rel([el["D11"], el["D12"], el["D22"]], [o.D11, o.D12, o.D22], 1.0))
            rcp, psi = kern.energy_psi(Pn, eps)
            worst = max(worst, _rel(psi, O2K.energy_psi(Pn, eps), 1e-12))
    print(f"\nHAR elastic / Psi kernel vs O2 (n = 0.5, 0.3, 0, 0.7; random states): {worst:.2e}")
    assert worst <= 1e-13


def test_floor_closed_forms(kern):
    """K1.12-K1.14 (sheet 13; (S.49)/(S.50)): eps_v,f, eps'_f, Pi_f and its block (S.51a) against the sheet values and
    O2 (BA06 alpha0 0 / 5; HAR n = 1/2 closed form and the general-n bracketed root); p_ref and the p_min default."""
    Pb = _k2()
    assert abs(kern.floor_ev(Pb, 0.0)[0] - 0.0529831737) < 1e-10                    # K1.12
    assert kern.pref(Pb) == (Pb.p_ref, Pb.p_min) and Pb.p_min == 0.5                # defaultPmin = O2's default
    Ph = _tims()
    assert kern.pref(Ph) == (Ph.p_ref, Ph.p_min) and Ph.p_ref == 101.0               # p_ref = p_a under HAR (M-F8)
    assert abs(Ph.p_min - 0.505) < 1e-15
    assert abs(kern.floor_ev(Ph, 0.0)[0] - 9.83645033e-4) < 1e-12                  # K1.13
    evf, epsp = kern.floor_ev(Ph, 2e-4)                                              # K1.14
    assert abs(evf - 1.04102893e-3) < 6e-12 and abs(epsp - 0.08679793) < 6e-9, (evf, epsp)   # the sheet's digits
    worst = 0.0
    for P in (Pb, _k2(alpha0=5.0, p_min=50.0), Ph, _tims(n_e=0.3), _tims(n_e=0.7), _tims(n_e=0.0)):
        for es in (0.0, 2e-5, 2e-4, 3e-3):
            ko, oo = kern.floor_ev(P, es), O2K.floor_ev(P, es)
            worst = max(worst, _rel(ko[0], oo[0], 1e-3), _rel(ko[1], oo[1], 1e-3))
        rng = np.random.default_rng(7)
        for _ in range(12):
            ev = O2K.floor_ev(P, 0.0)[0] + rng.uniform(-2e-4, 4e-4)
            eps = ev / 3 + rng.normal(0, 1e-4, 3)
            eps -= (eps.sum() - ev) / 3
            kf, of = kern.floor_project(P, eps), O2K.floor_project(P, eps)
            assert kf["active"] == of.active and kf["in_domain"] == of.in_domain, (kf, of)
            worst = max(worst, _rel(kf["eps_f"], of.eps_f, 1e-3), _rel(kf["Phi"], of.Phi, 1.0), _rel(kf["dfv"], of.dfv, 1e-3))
            if kf["active"]:          # p = -p_min after the projection, eps_s and the eigenvalue differences kept (M-F6)
                rc, el = kern.elastic(P, kf["eps_f"])
                assert rc == 0 and abs(el["p"] + P.p_min) <= 1e-12 * P.p_min
                d0, d1 = eps - eps.mean(), kf["eps_f"] - kf["eps_f"].mean()
                assert np.abs(d0 - d1).max() <= 1e-15
    print(f"\nfloor closed forms / Pi_f kernel vs O2: {worst:.2e}")
    assert worst <= 1e-12


def test_initial_state_pi0_rules(kern):
    """K1.15 (S.53): unified rule (c2 of the cap; c2 := 0 for none) and the legacy rule, against the sheet values and
    O2; the §7 guard refusal (107) at a very dense non-hydrostatic start; an explicit pi_i0 overrides."""
    Pc = _k2(cap="smooth", c1=0.05, c2=0.15)
    Pn = _k2()
    v0 = 1.59
    rc, s, _ = kern.initial_state(Pc, -100 * I3, v0, None)
    assert rc == 0 and abs(s[6] + 50.995881) < 1e-6, s[6]                          # ramp_end
    rc, s, _ = kern.initial_state(Pc, -100 * I3, v0, None, pi0_rule="legacy")
    assert rc == 0 and abs(s[6] + 46.475800) < 1e-6, s[6]                          # apex
    rc, s, _ = kern.initial_state(Pn, -100 * I3, v0, None)
    assert rc == 0 and abs(s[6] + 46.475800) < 1e-6
    sig75 = np.diag(-100 * np.ones(3) + np.array([1, 1, -2]) * 75.0 / 3)
    rc, s, _ = kern.initial_state(Pn, sig75, v0, None)
    assert rc == 0 and abs(s[6] + 71.554175) < 1e-6
    rc, s, _ = kern.initial_state(Pc, -100 * I3, v0, -60.0)
    assert rc == 0 and s[6] == -60.0
    worst = 0.0
    for P in (Pc, _k2(cap="planar", c1=0.1, c2=0.1), Pn, _tims(cap="smooth", c1=0.05, c2=0.15), _tims()):
        for sig in (-100 * I3, sig75, np.array([[-90.0, 6.0, 0.0], [6.0, -100.0, 3.0], [0.0, 3.0, -125.0]]), -0.2 * I3):
            for rule in ("unified", "legacy"):
                o = O2.initial_state(P, sig, 1.70, None, pi0_rule=rule)
                rc, s, msg = kern.initial_state(P, sig, 1.70, None, pi0_rule=rule)
                assert rc == 0, msg
                worst = max(worst, _rel(s[6], o.pi_i, P.p_ref), _rel(s[:6], NK.t6(o.eps_e), _strain_scale(P)),
                            _rel(s[12], o.eps_f_v, _strain_scale(P)))
                assert int(s[16]) == o.n_f_init and bool(s[17]) == bool(o.flags["at_floor"])
    print(f"\ninitial_state (unified / legacy, floor at init) kernel vs O2: {worst:.2e}")
    assert worst <= GATE
    # the §7 guard B <= 0 at (pi_i0, psi_i0) (unified rule): very dense (psi_i0 ~ -1), non-hydrostatic start
    sigd = np.diag([-90.0, -100.0, -110.0])
    with pytest.raises(ValueError):
        O2.initial_state(Pc, sigd, 0.75, None)
    rc, _, msg = kern.initial_state(Pc, sigd, 0.75, None)
    assert rc == 107, (rc, msg)


HAR_REFUSE = {"k_zero": dict(k=0.0), "g_neg": dict(g=-1.0), "n_e_one": dict(n_e=1.0), "n_e_neg": dict(n_e=-0.1),
              "p_a_zero": dict(p_a=0.0), "p_min_neg": dict(p_min=-1.0)}


@pytest.mark.parametrize("label", list(HAR_REFUSE))
def test_validate_refuses_har_and_floor_like_o2(kern, label):
    P = O2.Params(**dict(TIMS, **HAR_REFUSE[label]))
    with pytest.raises(ValueError):
        P.validate()
    rc, msg, _ = kern.validate(P)
    assert rc != 0, f"kernel accepted a set O2 refuses ({label})"


def test_validate_s56_gated_to_smooth(kern):
    """(S.56) PI_SCAN_REL <= W_ramp/10, cap = smooth ONLY (round 3b, A3): c2 = 0.06 refused (W_ramp 0.0061), c2 = 0.07
    accepted (0.0122), planar (W_ramp = 0) and none accepted; BA06 p_min < 0 refused."""
    bad = O2.Params(**dict(K2, **MODES["paper"], cap="smooth", c1=0.05, c2=0.06))
    with pytest.raises(ValueError):
        bad.validate()
    assert kern.validate(bad)[0] == 26
    for over in (dict(cap="smooth", c1=0.05, c2=0.07), dict(cap="planar", c1=0.1, c2=0.1), dict(cap="none"),
                 dict(cap="planar", c1=0.0, c2=0.0)):
        P = O2.Params(**dict(K2, **MODES["paper"], **over)).validate()
        rc, msg, _ = kern.validate(P)
        assert rc == 0, (over, rc, msg)
    P = O2.Params(**dict(K2, **MODES["paper"], p_min=-0.1))
    assert kern.validate(P)[0] == 25


def test_floor_init_counted(kern):
    """sheet 9.7 'initialState': a deck sigma0 above the floor is projected and counted (n_f_init = 1, at_floor), in
    the kernel exactly as in O2 (compared inside the path parity); non-vacuity here."""
    for name in ("FLOOR_HAR_init_above_floor", "FLOOR_HAR_init_above_floor_S53"):
        st0 = _o2_run(name)[3]
        assert st0.n_f_init == 1 and st0.flags["at_floor"] and st0.eps_f_v > 0.0
        assert report(kern, name)["init"] <= GATE


PRE_ROUND3 = [nm for nm in CASES if not nm.startswith(("HAR_", "PI0_", "FLOOR_"))]


def test_default_floor_leaves_existing_paths_bit_identical(kern):
    """owner decision (e): on every pre-round-3 path the kernel with the DEFAULT floor (p_min = 5e-3 p_ref, as the
    parity above runs it) is BIT-IDENTICAL to the kernel with the floor off (p_min = 0): state, stress, tangent,
    flags and iteration counts, step by step (the pi_i0 = None starts there have cap = none: the (S.53) rule is the
    pre-round-3 apex)."""
    import dataclasses
    for name in PRE_ROUND3:
        P, v0, pi0, st0, deps, sts, tans = _o2_run(name)
        P0 = dataclasses.replace(P, p_min=0.0)
        P0._ba06_given, P0._har_given = P._ba06_given, P._har_given
        assert P.p_min == 5e-3 * P.p_ref and P0.p_min == 0.0
        rc, s, _ = kern.initial_state(P, st0.sigma, v0, pi0)
        rc0, s0, _ = kern.initial_state(P0, st0.sigma, v0, pi0)
        assert rc == rc0 == 0 and np.array_equal(s, s0), name
        for d in deps:
            r, r0 = kern.step(P, s, d), kern.step(P0, s0, d)
            assert np.array_equal(r["state"], r0["state"]) and np.array_equal(r["sigma"], r0["sigma"]) \
                and np.array_equal(r["C"], r0["C"]) and r["info"] == r0["info"], (name, r["info"], r0["info"])
            if r["info"]["refusal"] != "OK":
                break
            s, s0 = r["state"], r0["state"]
    print(f"\ndefault floor vs p_min = 0, kernel, {len(PRE_ROUND3)} pre-round-3 paths: bit-identical")


def test_kernel_self_check_driver():
    """tests/ladrunonorsand_kernel_check.cpp compiled with -Wall -Wextra -Werror: every CHECK passes."""
    cxx = NK.compiler()
    if cxx is None:
        pytest.skip("no C++ compiler")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "nsk" + (".exe" if sys.platform == "win32" else ""))
        b = subprocess.run([cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-I", NK.KERNEL_DIR,
                            NK.DRIVER, "-o", exe], capture_output=True, text=True)
        assert b.returncode == 0, b.stderr
        r = subprocess.run([exe], capture_output=True, text=True, timeout=600)
    checks = [ln for ln in r.stdout.splitlines() if ln.startswith(("CHECK", "INFO", "SUMMARY"))]
    print("\n" + "\n".join(checks))
    assert r.returncode == 0 and "SUMMARY PASS" in r.stdout, "\n".join(checks)


def test_parity_summary(kern):
    """Max relative error per quantity over every path (the P1a parity numbers)."""
    rows = {name: report(kern, name) for name in CASES}
    worst = {q: max(r["errs"][q] for r in rows.values()) for q in REPORTED}
    worst_init = max(r["init"] for r in rows.values())
    print("\nkernel vs O2, max relative error per quantity over all paths (gate %.0e; band tangent gate %.0e):"
          % (GATE, BAND_TANGENT_GATE))
    for q in REPORTED:
        arg = max(rows, key=lambda nm: rows[nm]["errs"][q])
        print(f"  {q:>13s}: {worst[q]:.3e}   (worst path {arg})")
    print(f"  {'init':>13s}: {worst_init:.3e}")
    print(f"  {'init_sigma':>13s}: {max(r['init_sigma'] for r in rows.values()):.3e}")
    print(f"  {'lastsub_gap':>13s}: {max(r['lastsub_gap'] for r in rows.values()):.3e}   "
          "(O2 last-sub-increment CTO vs O2 chain on substepped steps; must be large)")
    tot = {key: sum(r[key] for r in rows.values())
           for key in ("n_o2", "plastic", "substepped", "band_steps", "iter_mismatch")}
    print("  paths %d, steps compared %d, plastic %d, substepped %d, band %d, refusal paths %d, "
          "iter-count mismatches %d" % (len(rows), tot["n_o2"], tot["plastic"], tot["substepped"], tot["band_steps"],
                                        sum(r["refusal"] is not None for r in rows.values()), tot["iter_mismatch"]))
    assert all(worst[q] <= GATE for q in QUANTITIES) and worst_init <= GATE
    assert worst["tangent_band"] <= BAND_TANGENT_GATE and worst["tangent_coalescent"] <= BAND_TANGENT_GATE
