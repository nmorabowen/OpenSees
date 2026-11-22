"""WP-144 (Zone B): the chained consistent tangent across substeps (equation sheet 144a 9.6, S.45-S.47; plan 2.8; owner decision 2026-10-01).

TEST RULE (plan 5): every expected value and tolerance in this file comes from the sheet (section cited), a closed form or a stated
truncation / round-off argument, and was written BEFORE the new O2 tangent was run on any of these paths.  No oracle tangent is an
expected value.  A tolerance is never loosened to make a test pass; a failing assertion is a finding (against the sheet or against
O2), not a reason to edit an oracle or a bound.

CONTRACT TESTED (sheet 9.6, first paragraph).  On a substepped increment O2.tangent(P, step(P, A, deps)) is the exact derivative of the
increment's FINAL stress with respect to the TOTAL strain increment deps, propagated through every sub-increment (any recursive-halving
sequence included):
        C[:, :, k, l] = d sigma_{n+1} / d deps_(kl)         (tensor shear component: eps_kl and eps_lk moved together)
and the tangent of the last sub-increment alone is NOT that derivative (the sheet measures 0.51 at m = 2 and 0.81-0.90 at m = 4 on the
AMP_STOP cap path, 3.6e-2 at m = 8 on a generic step).  The expected value is therefore the central finite difference of O2's own
stress update over the whole increment, which is exactly the object the contract names.

FD PROTOCOL (all tests).  Central difference over the six independent components E_J = e_k (x) e_k (normal), e_k (x) e_l + e_l (x) e_k
over two (shear: each tensor component at 1/2 each -> d sigma/dh = C_ij(kl) by the minor symmetry), at h in {1e-6, 1e-7}; the gate is
the best h (max over the six columns of ||C_J - FD_J|| / ||C_J|| at each h, then the smaller of the two), <= 1e-6, shear columns
included.  Argument for the tolerance: the FD error is O(h^2) truncation plus O(eps sigma / h) round-off; the sheet's own measured
chained-vs-FD errors are 1e-5 -> 9e-8 -> 9e-9 for h = 1e-6, 1e-7, 1e-8 (AMP step 11, the worst), 1e-7 -> 1e-9 for the generic m = 8
step; round-off at h = 1e-7 is ~1e-9 relative.  1e-6 is therefore >= 10x above the expected best-h value and ~1e4x below the
last-sub-increment error that the mutant leaves.  An h counts only if EVERY FD point is a plastic, non-refused increment with the SAME
substep count as the base increment (the ladder is a decision, not differentiated, sheet 9.6: an FD point that moves the ladder or the
branch differentiates a different map); the base increment must be plastic, not on the vertex branch, with the final stress off the
theta corners (theta in [0.1, pi/3 - 0.1], S.1; sheet 4.3).

PATHS (task WP-144 substep-tangent gate).
  P1  AMP_STOP smooth cap, K2 paper set (rho 0.7 / rho_bar 0.8), n = 40: sheet 9.6 (A) states that all 30 plastic increments 11-40 are
      substepped; every substepped plastic increment is gated, >= 30 of them required (the count is a path-premise guard, not an
      oracle value).
  P2  generic NON-COAXIAL plastic states forced to substep.  HOW: O2.step first attempts the whole increment and subdivides (2, 4, 8)
      only when the return map refuses it (api.step contract).  A preload of NPRE equal non-coaxial increments (pure-shear deviator plus
      a full shear tensor) from a dense isotropic start gives a plastic state with off-diagonal stress; a step of s * D with
      D = -1e-4 I + 2e-4 dev + SH/2 (compression plus the preload's deviator and shear) and s in {32, 64, 128} (total volumetric strain
      -9.6e-3 .. -3.8e-2) is far outside the basin of the whole-increment return map and is subdivided (premise asserted:
      flags['substeps'] >= 2, plastic, theta of the final stress mid-range).  Six states (three shear/deviator mixes x paper / fork CSL),
      each at the three scales: eighteen gated increments.
  P3  fork CSL + smooth cap: the shipped fork set of test_g1_cap.FORK_CAP_CASES['fork_default'] (M = 1.3309, rho 0.71 / rho_bar 0.75,
      csl_mode 'fork', cap 'smooth' c1 = 0.05, c2 = 0.15) on the near-isotropic path of deviator 2e-3 and 6e-3, n = 40; every substepped
      plastic increment is gated.
  P4  Gudehus-Argyris (rho 0.8, rho_bar 0.85, the GA census set), the P2 construction (paper CSL mixes A and B, fork CSL mix C), three scales each.
  P5  non-uniform fractions (recursive-halving shapes, sheet 9.6 (E)) through o2_algo.step_fractions on P1/P2/P4 states.
  R   regression: on every NON-substepped increment the tangent is unchanged, i.e. equals the (S.33) closed form (sheet 9.6 'm = 1
      reduces to (S.33)': exact to round-off at distinct trial eigenvalues, measured 1e-15) assembled by o2_algo.kernel.tangent_small
      from the step's own spectral data, to 1e-11 relative to max|C|; elastic increments give a^e (same function).  A central FD (the
      protocol above, h = 1e-6, 1e-7, 1e-6 gate) is added on a sample so the regression is not coupled to one assembly route.

MUTANT.  "Returns the last-substep tangent" (the behaviour this test replaces): killed by P1 (0.51 .. 0.90), P2/P3/P4 at m >= 4.  Note a
two-level ladder whose first half is elastic (pattern EP) has T_1 = identity and the last-sub tangent coincides with the chain there
(sheet 9.6 (B)), so the kill rests on the m >= 4 increments, which P1-P4 contain in number.
"""
import math

import numpy as np
import pytest

from conftest import K2_BASE, ORACLES, make_params

pytestmark = pytest.mark.slow          # ~8 min wall time (the forced-substep ladder refuses whole increments first)
O2 = ORACLES["O2"]
KER = O2.kernel

I3 = np.eye(3)
SIG0 = -100.0 * I3                 # isotropic start, p = p0 = -100 kPa, eps^e = 0
H_SET = (1.0e-6, 1.0e-7)
GATE = 1.0e-6
REG_TOL = 1.0e-11                  # round-off argument: sheet 9.6, measured 1e-15 / 8e-17 at m = 1
COLS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
SHEAR = [(0, 1), (0, 2), (1, 2)]


# ----------------------------------------------------------------------------------------------
# parameter sets (sheet 14 K2 constants; fork CSL mapping of sheet 15), restated here
# ----------------------------------------------------------------------------------------------
PAPER = dict(K2_BASE, rho=0.7, rho_bar=0.8)                                       # condition A: 0.875 >= beta = 0.75
FORK = dict(K2_BASE, csl_mode="fork", M=1.3309, N=0.3, N_bar=0.2, rho=0.71, rho_bar=0.75,
            e0=0.83, lambda_c=0.027, xi=0.45, p_a=101.325)                        # 0.9467 >= 0.875
GA_RHO, GA_RHOB = 0.8, 0.85                                                       # GA range [7/9, 1]; 0.941 >= 0.75
PAPER_GA = dict(PAPER, zeta="GA", rho=GA_RHO, rho_bar=GA_RHOB)
FORK_GA = dict(FORK, zeta="GA", rho=GA_RHO, rho_bar=GA_RHOB)
SMOOTH = dict(cap="smooth", c1=0.05, c2=0.15)                                     # sheet 10.2 defaults
PAPER_CAP = dict(PAPER, **SMOOTH)
FORK_CAP = dict(K2_BASE, M=1.3309, rho=0.71, rho_bar=0.75, csl_mode="fork", e0=0.83, lambda_c=0.027, xi=0.45,
                p_a=101.325, **SMOOTH)                                            # test_g1_cap FORK_CAP_CASES['fork_default']


def v0_for(kw, pi0, psi0):
    """Specific volume that gives image state parameter psi0 at pi_i0 (sheet 6, S.22)."""
    if kw["csl_mode"] == "paper":
        return psi0 + kw["v_c0"] - kw["lambda_tilde"] * math.log(-pi0)
    return 1.0 + kw["e0"] + psi0 - kw["lambda_c"] * (-pi0 / kw["p_a"]) ** kw["xi"]


def theta_of(sig):
    """Lode angle from (S.1): y = tr(s^3)/R^3 = cos3theta/sqrt6, theta = acos(sqrt6 y)/3 in [0, pi/3]."""
    p = np.trace(sig) / 3.0
    s = sig - p * I3
    R = np.linalg.norm(s)
    y = np.trace(s @ s @ s) / R ** 3
    return math.acos(max(-1.0, min(1.0, math.sqrt(6.0) * y))) / 3.0


def e_iso_total(amp):
    """Near-isotropic compression of the cap tests: -0.01 on every axis plus the deviator amp*diag(1, 0, -1)."""
    return -0.01 * I3 + amp * np.diag([1.0, 0.0, -1.0])


AMP_STOP = 2.0e-3
N_CAP = 40
PI0_CAP, PSI0_CAP = -80.0, -0.05          # dense start of the cap tests (pi_c = -172 kPa)

# generic non-coaxial preloads (pure-shear deviator + a full shear tensor: the stress is non-coaxial and theta stays mid-range) and
# loading directions D = -a I + b dev + SH/2 (compression keeps the increments outside the basin of the whole-increment return map
# at large scale; the deviator keeps theta mid-range: a plain triaxial-compression direction drives theta to the pi/3 corner)
SH_A = np.array([[0.0, 2.0e-4, 1.0e-4], [2.0e-4, 0.0, 0.0], [1.0e-4, 0.0, 0.0]])
SH_B = np.array([[0.0, 1.5e-4, -2.0e-4], [1.5e-4, 0.0, 3.0e-4], [-2.0e-4, 3.0e-4, 0.0]])
SH_C = np.array([[0.0, -3.0e-4, 0.0], [-3.0e-4, 0.0, 2.0e-4], [0.0, 2.0e-4, 0.0]])
DEV_A, DEV_B, DEV_C = np.diag([1.0, 0.0, -1.0]), np.diag([-1.0, 0.0, 1.0]), np.diag([0.0, 1.0, -1.0])
E_PRE_A, E_PRE_B, E_PRE_C = 4.0e-4 * DEV_A + SH_A, 4.0e-4 * DEV_B + SH_B, 4.0e-4 * DEV_C + SH_C
A_ISO, B_DEV = 1.0e-4, 2.0e-4


def _dirn(dev, sh):
    return -A_ISO * I3 + B_DEV * dev + 0.5 * sh


D_A, D_B, D_C = _dirn(DEV_A, SH_A), _dirn(DEV_B, SH_B), _dirn(DEV_C, SH_C)
SCALES = (32, 64, 128)


# ----------------------------------------------------------------------------------------------
# the FD machinery
# ----------------------------------------------------------------------------------------------
def unit_tensor(k, l):
    E = np.zeros((3, 3))
    if k == l:
        E[k, k] = 1.0
    else:
        E[k, l] = E[l, k] = 0.5
    return E


def _single_branch(s, m, pattern=None):
    """An FD point on the base map: plastic overall, accepted, same substep count m and (when O2 reports it) the same per-sub-increment
    elastic/plastic pattern as the base increment."""
    f = s.flags
    return (bool(f["plastic"]) and not f["refused"] and f.get("substeps", 1) == m
            and (pattern is None or f.get("pattern") == pattern))


def fd_columns(P, A, deps, h, m, stepper=None, pattern=None):
    """Central FD of the whole-increment final stress per column; None if any FD point is not a plastic, non-refused increment with
    the same substep count m (and pattern) as the base increment (a different map)."""
    stepper = stepper or O2.step
    out = {}
    for (k, l) in COLS:
        E = unit_tensor(k, l)
        sp_, sm_ = stepper(P, A, deps + h * E), stepper(P, A, deps - h * E)
        if not (_single_branch(sp_, m, pattern) and _single_branch(sm_, m, pattern)):
            return None
        out[(k, l)] = (sp_.sigma - sm_.sigma) / (2.0 * h)
    return out


def col_err(C, fd):
    return {kl: float(np.linalg.norm(C[:, :, kl[0], kl[1]] - fd[kl]) / np.linalg.norm(C[:, :, kl[0], kl[1]])) for kl in COLS}


def gate_increment(P, A, deps, label, stepper=None):
    """Substepped increment: tangent vs central FD of the whole increment.  Returns a row dict (best-h max column error, per-h values,
    substeps, theta of the final stress, valid h) or raises AssertionError on a broken premise.  `stepper(P, A, deps)` defaults to
    O2.step (the ladder); P5 passes O2.step_fractions with fixed fractions (a ladder-free sequence of sub-increments)."""
    stepper = stepper or O2.step
    stn = stepper(P, A, deps)
    f = stn.flags
    assert f["plastic"] and not f["refused"], f"{label}: premise broken, increment not plastic/accepted ({f})"
    m = f.get("substeps", 1)
    assert m >= 2, f"{label}: premise broken, increment was not substepped (substeps = {m})"
    assert not f["vertex"], f"{label}: premise broken, base increment is on the vertex branch"
    th = theta_of(stn.sigma)
    assert 0.1 < th < math.pi / 3.0 - 0.1, f"{label}: premise broken, final stress at the theta corner (theta = {th:.4f})"
    C = O2.tangent(P, stn)
    pattern = f.get("pattern")
    per_h, shear = {}, {}
    for h in H_SET:
        fd = fd_columns(P, A, deps, h, m, stepper, pattern)
        if fd is None:
            continue
        e = col_err(C, fd)
        per_h[h] = max(e.values())
        shear[h] = max(e[kl] for kl in SHEAR)
    best = min(per_h.values()) if per_h else None
    return dict(label=label, m=m, theta=th, per_h=per_h, shear=shear, best=best, pattern=pattern)


def report(rows, title):
    print(f"\n[{title}] substepped increments gated: {len(rows)}")
    for r in rows:
        ph = ", ".join(f"h={h:.0e}: {e:.2e}" for h, e in r["per_h"].items()) or "no valid h"
        print(f"  {r['label']:<34s} m={r['m']} {r.get('pattern') or '':<6s} theta={r['theta']:.3f}  {ph}")


def assert_rows(rows, title, min_valid_frac=1.0):
    valid = [r for r in rows if r["best"] is not None]
    assert len(valid) >= min_valid_frac * len(rows) and valid, \
        f"{title}: {len(rows) - len(valid)} of {len(rows)} increments have no FD-valid h (FD points leave the base ladder/branch)"
    bad = [r for r in valid if r["best"] > GATE]
    assert not bad, (f"{title}: {len(bad)} of {len(valid)} substepped increments miss the chained tangent by more than {GATE:.0e} "
                     f"(best-h max column error): "
                     + "; ".join(f"{r['label']} m={r['m']} best {r['best']:.3e} ({r['per_h']})" for r in bad[:8]))
    # shear columns are inside the max; restate them at the best h so a shear-only failure is named
    for r in valid:
        hb = min(r["per_h"], key=r["per_h"].get)
        assert r["shear"][hb] <= GATE, f"{title}: {r['label']} shear columns {r['shear'][hb]:.3e} at h = {hb:.0e}"


# ----------------------------------------------------------------------------------------------
# P1 : AMP_STOP smooth-cap path, n = 40
# ----------------------------------------------------------------------------------------------
def _cap_path(kw, amp, n):
    P = make_params("O2", **kw)
    st0 = O2.initial_state(P, SIG0, v0_for(kw, PI0_CAP, PSI0_CAP), PI0_CAP)
    deps = e_iso_total(amp) / n
    sts = O2.run_path(P, st0, np.array([deps] * n))
    assert not any(s.flags["refused"] for s in sts), "the cap path was refused: " + next(
        s.flags["reason"] for s in sts if s.flags["refused"])
    return P, st0, sts, deps


def test_substep_tangent_amp_stop_smooth_cap_n40_every_substepped_increment():
    """P1.  Sheet 9.6 (A): AMP_STOP smooth cap, n = 40, plastic increments 11-40 all substepped (>= 30 asserted).  Each gated against the
    whole-increment central FD (best h in {1e-6, 1e-7} <= 1e-6, shear columns included; theta = pi/6 off the corners).  >= 90% of the
    substepped increments must have an FD-valid h (the exclusion rule is stated here, before the run).
    Kills: last-sub-increment tangent (sheet: 0.51 at m = 2, 0.81-0.90 at m = 4); a chain that drops the pi_i column (S.45
    d pi/d pi_n = (1 - kappa)/c) or the v column; wrong cumulative fraction in S^v; a missing 1/c in d x/d pi_n."""
    P, st0, sts, deps = _cap_path(PAPER_CAP, AMP_STOP, N_CAP)
    states = [st0] + list(sts)
    rows = []
    for i, s in enumerate(sts):
        if s.flags["plastic"] and s.flags.get("substeps", 1) >= 2:
            rows.append(gate_increment(P, states[i], deps, f"AMP step {i + 1}"))
    report(rows, "P1 AMP_STOP smooth cap n=40")
    assert len(rows) >= 30, f"only {len(rows)} substepped plastic increments on the AMP_STOP n = 40 path (sheet 9.6 (A): 30)"
    assert_rows(rows, "P1", min_valid_frac=0.9)


# ----------------------------------------------------------------------------------------------
# P2 / P4 : generic non-coaxial states forced to substep
# ----------------------------------------------------------------------------------------------
def preloaded(kw, e_pre, npre):
    P = make_params("O2", **kw)
    st0 = O2.initial_state(P, SIG0, v0_for(kw, -80.0, -0.05), -80.0)
    sts = O2.run_path(P, st0, np.array([e_pre] * npre))
    assert not any(s.flags["refused"] for s in sts), "preload refused: " + next(s.flags["reason"] for s in sts if s.flags["refused"])
    st = sts[-1]
    off = np.linalg.norm(st.sigma - np.diag(np.diag(st.sigma)))
    assert off > 1e-3 * np.linalg.norm(st.sigma), f"preloaded state is coaxial (offdiag/|sigma| = {off / np.linalg.norm(st.sigma):.2e})"
    return P, st


GENERIC = {
    "paper_A": (PAPER, E_PRE_A, 10, D_A),
    "paper_B": (PAPER, E_PRE_B, 12, D_B),
    "paper_C": (PAPER, E_PRE_C, 10, D_C),
    "fork_A": (FORK, E_PRE_A, 10, D_A),
    "fork_B": (FORK, E_PRE_B, 12, D_B),
    "fork_C": (FORK, E_PRE_C, 10, D_C),
}
GENERIC_GA = {
    "GA_paper_A": (PAPER_GA, E_PRE_A, 10, D_A),
    "GA_paper_B": (PAPER_GA, E_PRE_B, 12, D_B),
    "GA_fork_C": (FORK_GA, E_PRE_C, 10, D_C),
}


def _generic_rows(cases):
    rows = []
    for name, (kw, e_pre, npre, d) in cases.items():
        P, st = preloaded(kw, e_pre, npre)
        for s in SCALES:
            rows.append(gate_increment(P, st, s * d, f"{name} s={s}"))
    return rows


def test_substep_tangent_generic_noncoaxial_states_forced_to_substep():
    """P2.  Six distinct non-coaxial plastic states (three shear/deviator mixes x paper / fork CSL), each stepped by s*D, s in {32, 64, 128}
    (forced subdivision, see module docstring); eighteen increments, every one must be FD-valid (generic constructed cases carry no
    exclusion).  theta off the corners.  Kills: last-sub tangent (3.6e-2 at m = 8 on this kind of increment); missing spin
    term g^Phi_ab in S^eps (S.46: the (eps^e_a - eps^e_b)/(eps~_a - eps~_b) of the map eps~ -> eps^e) or a^e spin; chain through the
    trial basis of the wrong sub-increment; E_J rotated with the committed rather than the trial eigenvectors."""
    rows = _generic_rows(GENERIC)
    report(rows, "P2 generic non-coaxial, forced substep")
    assert len(rows) >= 18 and len({r["label"].split()[0] for r in rows}) >= 6
    assert_rows(rows, "P2", min_valid_frac=1.0)


def test_substep_tangent_gudehus_argyris_case():
    """P4.  Same construction with zeta = 'GA' (rho 0.8, rho_bar 0.85) on the paper and fork CSL (GA has no compression-corner
    singularity of the flow Hessian, so its q_ab differs from WW's: a chain that hard-codes a WW q_ab or Omega_a is caught)."""
    rows = _generic_rows(GENERIC_GA)
    report(rows, "P4 GA generic, forced substep")
    assert_rows(rows, "P4", min_valid_frac=1.0)


# ----------------------------------------------------------------------------------------------
# P3 : fork CSL + smooth cap
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("amp", [AMP_STOP, 6.0e-3], ids=["dev2e-3", "dev6e-3"])
def test_substep_tangent_fork_csl_smooth_cap(amp):
    """P3.  Fork CSL (e_c = e0 - lambda_c (-p/p_a)^xi; d psi/d pi differs from the paper log law) + smooth cap on the near-isotropic
    path, n = 40; every substepped plastic increment gated (>= 20 required as a premise, >= 90% FD-valid).  Kills: a chain that
    differentiates the paper CSL inside psi_i (Pi_v, c = r'(pi_i) of S.37 built from the log law), the cap terms (q_a, Omega_pi of
    S.36) dropped from A or t of (S.45), v column missing."""
    P, st0, sts, deps = _cap_path(FORK_CAP, amp, N_CAP)
    states = [st0] + list(sts)
    rows = []
    for i, s in enumerate(sts):
        if s.flags["plastic"] and s.flags.get("substeps", 1) >= 2:
            rows.append(gate_increment(P, states[i], deps, f"forkcap amp={amp:g} step {i + 1}"))
    report(rows, f"P3 fork CSL + smooth cap amp={amp:g}")
    assert len(rows) >= 20, f"only {len(rows)} substepped plastic increments on the fork-cap path"
    assert_rows(rows, f"P3 amp={amp:g}", min_valid_frac=0.9)


# ----------------------------------------------------------------------------------------------
# P5 : non-uniform fractions (a recursive-halving shape), sheet 9.6 (E)
# ----------------------------------------------------------------------------------------------
def test_substep_tangent_non_uniform_fractions_recursive_halving_shapes():
    """P5.  Sheet 9.6: 'any recursive-halving sequence is covered by the same algebra' (state map for alpha_k > 0, sum alpha_k = 1).  The
    increment is taken as fixed fractions alpha_k * deps via o2_algo.step_fractions (no ladder: the FD points cannot move the
    decomposition), and the returned tangent is gated against the central FD of the same fixed-fraction map, same protocol and
    tolerance.  Fractions (a recursive-halving shape each): (1/2, 1/4, 1/4) on AMP step 11, (1/4, 1/4, 1/4, 1/8, 1/8) on AMP step 20,
    (1/2, 1/4, 1/8, 1/8) on the generic paper_A state at s = 32, (3/8, 3/8, 1/4) on the generic fork_B state at s = 40, and the GA
    state GA_paper_B at s = 32 with (1/2, 1/4, 1/8, 1/8).  Every sub-increment must be accepted (premise; the scales are chosen so the largest sub-increment is <= 16 D, inside the basin in which the return map accepts a single backward-Euler step: the s = 64 ladder level is m = 4, so a half step of 32 D is refused); the per-sub-increment
    elastic/plastic pattern is held fixed across FD points.
    Kills: a chain that assumes alpha_k = 1/m (the recursion uses alpha_k in T_k = S^eps_k + alpha_k E_J and in the cumulative fraction
    of S^v); the last-sub-increment tangent."""
    P, st0, sts, deps = _cap_path(PAPER_CAP, AMP_STOP, N_CAP)
    states = [st0] + list(sts)
    cases = [
        ("AMP step 11 (1/2,1/4,1/4)", P, states[10], deps, (0.5, 0.25, 0.25)),
        ("AMP step 20 (1/4x3,1/8x2)", P, states[19], deps, (0.25, 0.25, 0.25, 0.125, 0.125)),
    ]
    for name, scale, frs in (("paper_A", 32, (0.5, 0.25, 0.125, 0.125)), ("fork_B", 40, (0.375, 0.375, 0.25))):
        kw, e_pre, npre, d = GENERIC[name]
        Pg, st = preloaded(kw, e_pre, npre)
        cases.append((f"{name} s={scale} {frs}", Pg, st, scale * d, frs))
    kw, e_pre, npre, d = GENERIC_GA["GA_paper_B"]
    Pg, st = preloaded(kw, e_pre, npre)
    cases.append(("GA_paper_B s=32 (1/2,1/4,1/8,1/8)", Pg, st, 32 * d, (0.5, 0.25, 0.125, 0.125)))
    rows = []
    for label, Pc, A, dp, frs in cases:
        assert abs(sum(frs) - 1.0) < 1e-15
        rows.append(gate_increment(Pc, A, dp, label, stepper=lambda P_, A_, d_, frs=frs: O2.step_fractions(P_, A_, d_, frs)))
    report(rows, "P5 non-uniform fractions")
    assert_rows(rows, "P5", min_valid_frac=1.0)


# ----------------------------------------------------------------------------------------------
# V : the v column of the chain follows the EXPONENTIAL update (sheet 1.2, 9.6; G2 owner decision 2026-10-01)
# ----------------------------------------------------------------------------------------------
D2_ALL_SHEARS = np.array([[-1.5e-3, 0.6e-3, -0.4e-3],
                          [0.6e-3, -0.2e-3, 0.5e-3],
                          [-0.4e-3, 0.5e-3, 0.9e-3]])             # the sheet 9.6 / chain_fd (B) increment: all three shears
V_MID_PATH = 1.45                                                  # a state "built mid-path": v far from v0 = 1.70


def _state_B(v_forced=None):
    P = make_params("O2", **PAPER)
    st0 = O2.initial_state(P, SIG0, v0_for(PAPER, PI0_CAP, PSI0_CAP), PI0_CAP)
    d1 = -2.0e-3 * I3 + 2.0e-3 * np.diag([1.0, 0.0, -1.0])
    st1 = O2.step(P, st0, d1)
    assert st1.flags["plastic"] and not st1.flags["refused"]
    if v_forced is not None:
        st1 = st1.copy()
        st1.v = v_forced
    return P, st1


V_FRACTIONS = {"m=8": (0.125,) * 8, "m=2": (0.5, 0.5), "alpha=(1/2,1/4,1/8,1/8)": (0.5, 0.25, 0.125, 0.125)}


def _v_rows(P, A):
    return [gate_increment(P, A, D2_ALL_SHEARS, name,
                           stepper=lambda P_, A_, d_, frs=frs: O2.step_fractions(P_, A_, d_, frs))
            for name, frs in V_FRACTIONS.items()]


def test_substep_chain_v_column_is_v_k_plus_1_not_v0(monkeypatch):
    """Sheet 1.2 / 9.6 (G2 owner decision 2026-10-01): S^v_{k+1} = v_{k+1} (sum_{j<=k} alpha_j) tr E_J, the converged v of the
    sub-increment, NOT v0 (the superseded G0/G1 closed form S^v = v0 cum tr E_J).  Expected values from the sheet 9.6
    FD record on this increment (chain_fd (B), all three shears) at m = 8, 2 and the non-uniform (1/2, 1/4, 1/8, 1/8):
      * chain with v_{k+1}: chained tangent vs central FD of the whole increment 2.7e-8 .. 4.2e-10 (v forced to 1.45) and
        1e-9 .. 1.8e-9 (the real v), i.e. <= GATE = 1e-6 at the best h (the protocol of this file);
      * the v0 variant of S^v: 2.2e-6 (real v), 2.2e-5 .. 2.6e-5 (v = 1.45), for every h: the gate must REJECT it, i.e.
        the best-h error is > GATE on the real-v increments and >= 10 GATE on the v = 1.45 ones.
    The mutant is injected by replacing o2_algo.kernel.chain_propagate with a wrapper that hands it the committed v0 instead
    of the sub-increment's v_{k+1} (the exact change of the superseded rule).
    Kills: a chain (O2 or the C++ kernel it is the contract for) whose S^v still reads v0 or v_n; the gate's blindness to it."""
    for label, vf, lo in (("real v", None, GATE), ("v = 1.45", V_MID_PATH, 10.0 * GATE)):
        P, A = _state_B(vf)
        assert abs(A.v - A.v0) >= (1e-3 if vf is None else 0.2), "premise: v must differ from v0 for the discrimination"
        good = _v_rows(P, A)
        report(good, f"V chain with v_(k+1), {label}")
        assert_rows(good, f"V correct chain, {label}", min_valid_frac=1.0)

        orig = KER.chain_propagate

        def v0_variant(S_eps, S_pi, cum, alpha, v_new, res, nvec, _v0=A.v0, _orig=orig):
            return _orig(S_eps, S_pi, cum, alpha, _v0, res, nvec)

        with monkeypatch.context() as mp:
            mp.setattr(KER, "chain_propagate", v0_variant)
            bad = _v_rows(P, A)
        report(bad, f"V chain with v0 (the superseded rule), {label}")
        for r in bad:
            assert r["best"] is not None, f"{label} {r['label']}: no FD-valid h for the mutant"
            assert r["best"] > lo, (f"{label} {r['label']}: the v0 variant is NOT rejected: best-h error "
                                    f"{r['best']:.3e} <= {lo:.1e} (sheet 9.6: 2.2e-6 real v, 2.2e-5..2.6e-5 at v = 1.45)")


# ----------------------------------------------------------------------------------------------
# R : regression on non-substepped increments
# ----------------------------------------------------------------------------------------------
def s33_reference(st):
    """(S.33) closed form from the step's own spectral data (sheet 9.4): the pre-chain assembly, unchanged by the contract."""
    c = st.cache
    return KER.tangent_small(c["atilde"], c["sig"], c["eps_tr"], c["nvec"])


def assert_equals_s33(P, st, label):
    C = O2.tangent(P, st)
    ref = s33_reference(st)
    err = float(np.max(np.abs(C - ref)) / np.max(np.abs(ref)))
    assert err <= REG_TOL, f"{label}: tangent differs from the (S.33) closed form by {err:.3e} of max|C| (> {REG_TOL:.0e})"
    return err


def fd_gate_single(P, A, deps, label):
    """Central FD gate for a NON-substepped increment (same protocol, any branch, same substep count 1)."""
    stn = O2.step(P, A, deps)
    assert stn.flags.get("substeps", 1) == 1
    C = O2.tangent(P, stn)
    plastic = bool(stn.flags["plastic"])
    best = None
    for h in H_SET:
        pts = {}
        ok = True
        for (k, l) in COLS:
            E = unit_tensor(k, l)
            sp_, sm_ = O2.step(P, A, deps + h * E), O2.step(P, A, deps - h * E)
            for s in (sp_, sm_):
                if s.flags["refused"] or s.flags.get("substeps", 1) != 1 or bool(s.flags["plastic"]) != plastic:
                    ok = False
            if not ok:
                break
            pts[(k, l)] = (sp_.sigma - sm_.sigma) / (2.0 * h)
        if ok:
            e = max(col_err(C, pts).values())
            best = e if best is None else min(best, e)
    assert best is not None, f"{label}: no FD-valid h"
    assert best <= GATE, f"{label}: best-h FD error {best:.3e} > {GATE:.0e}"
    return best


def test_regression_non_substepped_tangent_is_unchanged():
    """R.  (i) AMP smooth-cap path: every increment with substeps == 1 (elastic and plastic) returns the (S.33) tangent, 1e-11 of max|C|
    (sheet 9.6 'm = 1 reduces to (S.33)'; at distinct trial eigenvalues exact to round-off, 1e-15 measured); (ii) the three generic
    states and the two GA states at s = 1 (plastic, not substepped), same equality; (iii) an elastic unloading increment (-D, one preload-sized step back) of a generic
    state returns a^e (same function; flag plastic == False asserted); (iv) FD gate (best h <= 1e-6) on a sample of (i)-(iii), so the
    regression is not coupled to one assembly route; (v) the chain forced on one sub-increment (step_fractions, alpha = (1,)) equals
    (S.33) to 1e-11 of max|C| on the generic states (sheet 9.6 'm = 1 reduces to (S.33)').  Kills: a chained assembly that is applied (and wrong) at m = 1; a chain that
    perturbs the elastic or single-step return path; a substep-count test that mis-classifies m = 1."""
    P, st0, sts, deps = _cap_path(PAPER_CAP, AMP_STOP, N_CAP)
    states = [st0] + list(sts)
    n_checked, n_plastic, n_elastic, worst = 0, 0, 0, 0.0
    fd_idx = []
    for i, s in enumerate(sts):
        if s.flags.get("substeps", 1) == 1:
            worst = max(worst, assert_equals_s33(P, s, f"AMP step {i + 1}"))
            n_checked += 1
            if s.flags["plastic"]:
                n_plastic += 1
            else:
                n_elastic += 1
            fd_idx.append(i)
    assert n_plastic >= 3 and n_elastic >= 3, f"AMP path gives {n_plastic} plastic / {n_elastic} elastic unsubstepped steps"
    # FD sample: the first elastic, the first and last plastic unsubstepped steps
    pick = [next(i for i in fd_idx if not sts[i].flags["plastic"])]
    pl = [i for i in fd_idx if sts[i].flags["plastic"]]
    pick += [pl[0], pl[-1]]
    fd_errs = {f"AMP step {i + 1}": fd_gate_single(P, states[i], deps, f"AMP step {i + 1}") for i in pick}
    for name, (kw, e_pre, npre, d) in {**GENERIC, **GENERIC_GA}.items():
        Pg, st = preloaded(kw, e_pre, npre)
        n1 = O2.step(Pg, st, d)
        assert n1.flags["plastic"] and n1.flags.get("substeps", 1) == 1, f"{name}: s = 1 premise broken ({n1.flags})"
        worst = max(worst, assert_equals_s33(Pg, n1, f"{name} s=1"))
        fd_errs[f"{name} s=1"] = fd_gate_single(Pg, st, d, f"{name} s=1")
        n_checked += 1
        unload = O2.step(Pg, st, -1.0 * d)
        assert not unload.flags["plastic"] and unload.flags.get("substeps", 1) == 1, f"{name}: unloading is not elastic ({unload.flags})"
        worst = max(worst, assert_equals_s33(Pg, unload, f"{name} unloading"))
        fd_errs[f"{name} unloading"] = fd_gate_single(Pg, st, -1.0 * d, f"{name} unloading")
        n_checked += 1
    # sheet 9.6 'm = 1 reduces to (S.33)': the chain forced on a single sub-increment (alpha = (1,)) equals (S.33) at distinct trial
    # eigenvalues, exact to round-off (measured 1e-15 and 8e-17), 1e-11 of max|C|
    chain1 = {}
    for name, (kw, e_pre, npre, d) in {**GENERIC, **GENERIC_GA}.items():
        Pg, st = preloaded(kw, e_pre, npre)
        s1 = O2.step_fractions(Pg, st, d, (1.0,))
        assert s1.flags["plastic"] and not s1.flags["refused"]
        C1 = O2.tangent(Pg, s1)
        ref = s33_reference(s1)
        chain1[name] = float(np.max(np.abs(C1 - ref)) / np.max(np.abs(ref)))
        assert chain1[name] <= REG_TOL, f"{name}: chain(m = 1) differs from (S.33) by {chain1[name]:.3e} of max|C|"
    print(f"\n[R] chain(m=1) vs (S.33), max over generic states: {max(chain1.values()):.2e}")
    print(f"\n[R] non-substepped increments checked {n_checked}; max |C - (S.33)|/max|C| = {worst:.2e}; FD best-h errors:")
    for k, v in fd_errs.items():
        print(f"   {k:<24s} {v:.2e}")
