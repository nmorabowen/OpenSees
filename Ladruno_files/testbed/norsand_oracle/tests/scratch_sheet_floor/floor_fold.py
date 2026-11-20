"""Sheet round 2026-10-03, item 4: which nested-scan step is right, 1e-4 (sheet 8/10.2) or 1e-3 (O2, kernel)?

Geometry of r(pi_i) on the AMP_STOP smooth-cap path (K2 paper set, rho 0.7 / rho_bar 0.8, c1 0.05, c2 0.15), with
the shipped O2 (exponential v-update, PI_SCAN_REL = 1e-3), imported read-only from the WP-144 worktree:
  (A) over the whole path at n = 20, 40, 80, 160: at the CONVERGED iterate of every non-substepped plastic step a
      dense scan of r(pi_i) gives the roots pi_1 (near), pi_2, pi_3 and the first local extremum of |r| beyond pi_2;
      d2 = |pi_2 - pi_n|/|pi_n| is the distance a scan step must not exceed to bracket the near root at the first
      point, d_x = |pi_x - pi_n|/|pi_n| the distance it must not exceed for a skipped pair to still be caught by the
      |r|-growth rule (beyond d_x a scan would walk to the far root in silence). Nested EvalErrors raised inside
      the line search are counted by kind.
  (B) at the n = 40 step-8 iterate (sheet 10.2): roots and dip width at dlam = 7.5e-4 ... 8.3e-4 through the fold,
      and what solve_pi does at PI_SCAN_REL = 1e-3 and 1e-4 (scan range held at |pi_n| for both).
  (C) the ramp width in pi_i at fixed p (the scale of d_x) and its closed-form estimate.
Run:  python -u floor_fold.py
"""
import math
import os
import sys
import warnings
from collections import Counter

import numpy as np

ORACLE = r"C:/Users/nmora/Documents/Github/OpenSees/.claude/worktrees/ladrunonorsand-implementation-review-7bbd75/Ladruno_files/testbed/norsand_oracle"
sys.path.insert(0, ORACLE)
sys.path.insert(0, os.path.join(ORACLE, "tests"))
sys.path.insert(0, os.path.join(ORACLE, "tests", "scratch_triage_conv"))
warnings.simplefilter("ignore")
from conftest import make_params                          # noqa: E402
import o2_algo                                            # noqa: E402
from o2_algo import kernel as K                           # noqa: E402
from diag_cap import kw_paper, v0_for, e_iso_total, SIG0  # noqa: E402

kw = kw_paper(cap="smooth", c1=0.05, c2=0.15)
v0 = v0_for(kw, -80.0, -0.05)
P = make_params("O2", **kw)

events = Counter()
ncalls = [0]
_orig_solve_pi = K.solve_pi


def solve_pi_logged(P_, inv, dlam, v, pi_n):
    ncalls[0] += 1
    try:
        return _orig_solve_pi(P_, inv, dlam, v, pi_n)
    except K.EvalError as e:
        events[str(e)] += 1
        raise


K.solve_pi = solve_pi_logged


def root_structure(inv, dlam, v, pi_n, rel_max=0.5, npts=20001):
    s = np.linspace(0.0, rel_max, npts)
    pis = pi_n * (1.0 + s)
    rs = np.full(npts, np.nan)
    for i, pi in enumerate(pis):
        try:
            rs[i], _ = K._pi_residual(P, inv, dlam, v, pi_n, pi)
        except K.EvalError:
            break
    roots = [0.5 * (s[i] + s[i + 1]) for i in range(npts - 1)
             if np.isfinite(rs[i]) and np.isfinite(rs[i + 1]) and rs[i] * rs[i + 1] < 0]
    dx = None
    if len(roots) >= 2:
        i2 = int(np.searchsorted(s, roots[1]))
        j = i2 + 1
        while j + 1 < npts and np.isfinite(rs[j + 1]) and abs(rs[j + 1]) >= abs(rs[j]):
            j += 1
        dx = s[j]
    return roots, dx, rs[0]


print("=== (A) fold geometry along the AMP_STOP path, shipped O2 (PI_SCAN_REL = 1e-3) ===")
for n in (20, 40, 80, 160):
    events.clear()
    ncalls[0] = 0
    st = o2_algo.initial_state(P, SIG0, v0, -80.0)
    deps = e_iso_total(2e-3) / n
    sts = o2_algo.run_path(P, st, np.array([deps] * n))
    refused = [i for i, s_ in enumerate(sts) if s_.flags["refused"]]
    nsub = Counter(s_.flags.get("substeps", 1) for s_ in sts)
    ev = dict(events)
    d2s, dxs, nroots = [], [], Counter()
    for i, s_ in enumerate(sts):
        if not s_.flags["plastic"] or s_.flags["refused"] or s_.flags.get("substeps", 1) > 1:
            continue
        w_ = np.linalg.eigh(s_.eps_e)[0]
        inv = K.invariants(K.elastic(P, w_).sig)
        pi_prev = sts[i - 1].pi_i if i > 0 else -80.0
        roots, dx, r0 = root_structure(inv, s_.dlam, s_.v, pi_prev)
        nroots[len(roots)] += 1
        if len(roots) >= 2:
            d2s.append(roots[1])
            dxs.append(dx)
    print(f"  n={n:4d}: refused={refused[:1] or 'none'}, substep levels={dict(nsub)}, nested EvalErrors={ev or 'none'}, "
          f"nested solves={ncalls[0]}; converged-iterate root counts {dict(nroots)}; "
          f"min d2 = {min(d2s) if d2s else float('nan'):.3e}|pi_n|, min d_x = {min(dxs) if dxs else float('nan'):.3e}|pi_n|")

print("\n=== (B) n = 40, step-8 iterate (sheet 10.2): roots through the fold; what each scan step selects ===")
n = 40
st = o2_algo.initial_state(P, SIG0, v0, -80.0)
deps = e_iso_total(2e-3) / n
for _ in range(8):
    st = o2_algo.step(P, st, deps)
    assert not st.flags["refused"]
v = st.v * math.exp(float(np.trace(deps)))          # exponential v-update (G2); r2_fold used the linear one
eps_tr = np.linalg.eigh(st.eps_e + deps)[0]
pi_n = st.pi_i
print(f"  pi_n = {pi_n:.4f}, p_tr = {K.elastic(P, eps_tr).p:.4f}")
for dlam in (7.5e-4, 7.7e-4, 7.877e-4, 8.0e-4, 8.05e-4, 8.1e-4, 8.3e-4):
    eps_e = eps_tr + dlam / 3.0
    inv = K.invariants(K.elastic(P, eps_e).sig)
    roots, dx, r0 = root_structure(inv, dlam, v, pi_n, rel_max=0.4, npts=40001)
    desc = ", ".join(f"{abs(pi_n)*r:.4f} kPa ({r:.2e}|pi_n|)" for r in roots)
    picks = {}
    for rel in (1e-3, 1e-4):
        K.PI_SCAN_REL = rel
        K.PI_SCAN_MAX = int(round(1.0 / rel))
        try:
            pi_sel, c, its = _orig_solve_pi(P, inv, dlam, v, pi_n)
            d_sel = abs(pi_sel - pi_n) / abs(pi_n)
            if roots:
                j = min(range(len(roots)), key=lambda j: abs(d_sel - roots[j]))
                which = "near root" if j == 0 else f"ROOT {j+1} (NOT the near one)"
            else:
                which = "?"
            picks[rel] = f"{which} at {abs(pi_sel-pi_n):.4f} kPa, {its} evals, c = {c:+.3f}"
        except K.EvalError as e:
            picks[rel] = f"EvalError '{e}'"
    K.PI_SCAN_REL, K.PI_SCAN_MAX = 1e-3, 1000
    dip = f"{abs(pi_n)*(roots[1]-roots[0]):.4f} kPa = {roots[1]-roots[0]:.2e}|pi_n|" if len(roots) >= 2 else "-"
    print(f"  dlam={dlam:.4e}: r(pi_n)={r0:+.3e}; {len(roots)} root(s) at [{desc}]; dip root1->root2 = {dip}; "
          f"|r| extremum beyond root 2 at {abs(pi_n)*dx if dx else float('nan'):.3f} kPa ({dx if dx else float('nan'):.2e}|pi_n|)")
    print(f"      scan 1e-3 -> {picks[1e-3]}\n      scan 1e-4 -> {picks[1e-4]}")

print("\n=== (C) ramp width in pi_i at fixed p (the scale of d_x) ===")
for p in (-100.0, -160.0, -172.0):
    pi1 = K.pi_of_eta(P, p, P.c1 * P.M)
    pi2 = K.pi_of_eta(P, p, P.c2 * P.M)
    est = (P.c2 - P.c1) * abs(pi2) * (pi2 / p) ** (P.N / (1 - P.N))
    print(f"  p = {p:.1f}: pi(eta1) = {pi1:.3f}, pi(eta2) = {pi2:.3f}, ramp width = {abs(pi1-pi2):.3f} kPa = "
          f"{abs(pi1-pi2)/abs(pi2):.3e}|pi(eta2)|; estimate (c2-c1)|pi|(pi/p)^(N/(1-N)) = {est:.3f}")
print("\ndone")
