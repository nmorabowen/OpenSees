"""WP-144 P1a (Zone B): parity of the C++ kernel SRC/material/nD/LadrunoNorSandKernel.h with O2.

O2 (o2_algo) is the kernel's contract: the same algorithm, branches and constants. This file
drives the SAME strain increments through O2 (o2_algo.step / o2_algo.tangent) and through the
C++ kernel (ctypes shim ns_shim.cpp, built here with g++) and compares, step by step:
    sigma, eps_e, pi_i, v, eps_p_v, eps_p_s, D (the step's dissipation), the 6x6 tangent,
    and the flags: refused (and its finest reason), plastic, vertex, cap_active, substeps,
    and the local-Newton / nested-solve iteration counts (equal on every step: the same
    algorithm with the same constants; a constant changed without changing the answer shows
    only there).

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

THE EXCEPTION -- the tangent in two round-off conditioning bands (state quantities stay at
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
Inside both bands the tangent gate is 1e-7: the O(|theta - theta_c|) accuracy of (S.8)/(S.9)
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
import ns_kernel as NK                   # noqa: E402

GATE = 1.0e-10
ZERO_FLOOR = 1.0e-6          # x natural scale, for quantities that pass through zero
CORNER_BAND = 1.0e-7         # |sin 3theta| below: round-off corner band
VERTEX_BAND = 1.0e-3         # R/|p| below: near-vertex band
BAND_TANGENT_GATE = 1.0e-7   # tangent gate inside the two bands
I3 = np.eye(3)
SIG0 = -100.0 * I3

QUANTITIES = ("sigma", "eps_e", "pi_i", "v", "eps_p_v", "eps_p_s", "D", "tangent")
REPORTED = QUANTITIES + ("tangent_band", "o2_1ulp_band")

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
}


def o2_params(c):
    kw = dict(K2)
    kw.update(MODES[c["mode"]])
    kw.update(c["over"])
    return O2.Params(**kw)


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
            if abs(res) <= O2api.TRIAX_LAT_TOL_REL * abs(P.p0):
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
def _o2_run(name):
    c = CASES[name]
    P = o2_params(c)
    v0 = v0_of(P, c)
    pi0 = c["start"][0]
    st0 = O2.initial_state(P, c["sigma0"], v0, pi0)
    deps = increments(P, c, st0)
    sts, tans = [], []
    st = st0
    for d in deps:
        st = O2.step(P, st, d)
        sts.append(st)
        tans.append(O2.tangent(P, st))
        if st.flags["refused"]:
            break
    return P, v0, pi0, st0, deps, sts, tans


def _natural_scales(P):
    return dict(sigma=abs(P.p0), eps_e=P.kappa_hat, pi_i=abs(P.p0), v=1.0, eps_p_v=P.kappa_hat,
                eps_p_s=P.kappa_hat, D=abs(P.p0) * P.kappa_hat)


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


def compare(kern, name):
    P, v0, pi0, st0, deps, sts, tans = _o2_run(name)
    nat = _natural_scales(P)
    rep = dict(name=name, n_o2=len(sts), errs={q: 0.0 for q in REPORTED}, flag_mismatch=[], iter_mismatch=0,
               refusal=None, substepped=0, plastic=0, band_steps=0)
    errs = rep["errs"]
    # initial state (sigma0 is recovered through the elastic inversion, whose Newton stops at 1e-13 |p|)
    rc, ks, msg = kern.initial_state(P, st0.sigma, v0, pi0)
    assert rc == 0, f"kernel initialState refused: {rc} {msg}"
    rep["init"] = max(_rel(ks[:6], NK.t6(st0.eps_e), nat["eps_e"]), _rel(ks[6], st0.pi_i, nat["pi_i"]),
                      _rel(ks[7], st0.v, 1.0), _rel(ks[8], v0, 1.0))
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
                        ("eps_p_s", ns[10], o.eps_p_s)):
            errs[q] = max(errs[q], _rel(x, y, nat[q]))
        D_k.append(ns[11])
        D_o.append(o.D)
        Co = NK.c4_to_c6(tans[k])
        et = np.abs(r["C"] - Co).max() / np.abs(Co).max()
        inv = O2.kernel.invariants(np.linalg.eigvalsh(0.5 * (o.sigma + o.sigma.T)))
        in_band = (not inv.vertex) and (abs(math.sin(3.0 * inv.theta)) < CORNER_BAND
                                        or inv.R < VERTEX_BAND * abs(inv.p))
        if in_band and o.flags["plastic"]:
            rep["band_steps"] += 1
            errs["tangent_band"] = max(errs["tangent_band"], et)
            errs["o2_1ulp_band"] = max(errs["o2_1ulp_band"], _o2_1ulp(P, ost, d, Co))
        else:
            errs["tangent"] = max(errs["tangent"], et)
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
          f"band {rep['band_steps']} refusal {rep['refusal']} iter-count mismatches {rep['iter_mismatch']}"
          f"\n  init={rep['init']:.2e} init_sigma={rep['init_sigma']:.2e}  {line}")
    assert not rep["flag_mismatch"], f"branch flags differ (step, flag, O2, kernel): {rep['flag_mismatch']}"
    # same algorithm, same constants => the same number of local Newton iterations and of nested
    # r(pi_i) evaluations on every step (a constant changed without changing the answer shows here)
    assert rep["iter_mismatch"] == 0, f"{rep['iter_mismatch']} steps with different local/nested iteration counts"
    assert rep["init"] <= GATE, f"initial state differs: {rep['init']:.3e}"
    assert rep["init_sigma"] <= 1e-12, f"stress(initial state) != sigma0: {rep['init_sigma']:.3e}"
    bad = {q: errs[q] for q in QUANTITIES if not errs[q] <= GATE}
    assert not bad, f"parity gate {GATE:.0e} exceeded: {bad}"
    assert errs["tangent_band"] <= BAND_TANGENT_GATE, \
        f"band tangent {errs['tangent_band']:.3e} > {BAND_TANGENT_GATE:.0e} (O2 1-ulp {errs['o2_1ulp_band']:.3e})"


EXPECTED_REFUSALS = {          # sheet 3.2 / 10.1 / 16.6: O2 refuses these after 2^8 substeps
    "CAP_planar_AMPSTOP_n40": True,
    "CAP_none_AMPSTOP_n40": True,
}


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
    rep = report(kern, "CAP_smooth_AMPSTOP_n40")
    assert rep["refusal"] is None
    assert rep["substepped"] >= 1 and rep["plastic"] >= 1


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
    tot = {key: sum(r[key] for r in rows.values())
           for key in ("n_o2", "plastic", "substepped", "band_steps", "iter_mismatch")}
    print("  paths %d, steps compared %d, plastic %d, substepped %d, band %d, refusal paths %d, "
          "iter-count mismatches %d" % (len(rows), tot["n_o2"], tot["plastic"], tot["substepped"], tot["band_steps"],
                                        sum(r["refusal"] is not None for r in rows.values()), tot["iter_mismatch"]))
    assert all(worst[q] <= GATE for q in QUANTITIES) and worst_init <= GATE
    assert worst["tangent_band"] <= BAND_TANGENT_GATE
