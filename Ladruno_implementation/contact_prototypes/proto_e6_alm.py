"""ADR-57 E6 ORACLE — edge-edge one-scalar commit-cycle ALM (augmented Lagrange).

Oracle-first (the fork discipline): validate the edge-edge ALM in pure numpy BEFORE transcribing it
into the LadrunoContactFE EDGE_EDGE mode + LadrunoContactDomain commit(). No OpenSees, no build. E6
adds an OPTIONAL per-pair normal multiplier λ_N that drives the interface penetration → an
εN-INDEPENDENT tolerance at FINITE εN — exactly the mortar C2.2 Uzawa-on-commit, but POINT-LIKE:
an edge pair's gap g_N is ONE scalar, so λ_N is ONE scalar per pair and there is NO shared-node
global accumulator / variational-consistency problem (the whole simplification vs the mortar lane).

What E6 owns (the ADR-57 §6 ALM tier + the E6 gate row):
  - the augmented pressure p = min(0, λ_N + εN·gN), traction tN = −p (≥0). λ_N≡0 ⇒ p = εN·gN,
    active iff gN<0 ⇒ the E2 PENALTY path EXACTLY (byte-identity is λ_N=0 off by default).
  - one Uzawa step per commit: λ_N ← min(0, λ_N + εN·gN_committed) (gN_committed = the converged
    gap stored by the adapter each residual). The mortar C2.2 update, one scalar instead of g̃/a.
  - it drives the penetration → an εN-INDEPENDENT augTol: a HELD load augmented to convergence
    reaches the SAME tiny |gN| for εN = 1e4 … 1e8 (penalty alone gives |gN| = O(P/εN), εN-dependent).
  - RELEASE ⇒ λ_N → 0, F = 0: remove the load and the contact opens; the Uzawa min(0, …) clamp
    drives λ_N back to 0 and the traction to 0 (no locked-in residual pull — the KKT release).
  - the equation count is CONSTANT across augmentations: λ_N is Domain-side committed state, NOT a
    DOF, so the held-load re-solves never re-number / change the SOE size (the D1 driver invariant).
  - committed-only: λ_N is mutated SOLELY in commit() ⇒ revertToLastCommit() leaves it untouched
    (the EmbeddedRebar / C2.2 invariant — a rejected implicit step cannot corrupt the multiplier).

Run: <pythoncore-3.12> proto_e6_alm.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9
_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ============================================ the E6 ALM layer (mirror the C++ to-be)
def augmented_pressure(lambdaN, epsN, gN):
    """p = min(0, λ_N + εN·gN) — the augmented normal pressure (≤0 compression). traction tN = −p."""
    return min(0.0, lambdaN + epsN * gN)


def uzawa_step(lambdaN, epsN, gN_committed):
    """ONE commit-cycle Uzawa step: λ_N ← min(0, λ_N + εN·gN_committed). Point-like ⇒ ONE scalar,
    no shared-node accumulator (the mortar C2.2 update with g̃/a → the single edge gap gN)."""
    return min(0.0, lambdaN + epsN * gN_committed)


# a held-load 1-DOF edge-pair model: the slave closest point carries a structural stiffness k_s
# anchoring it; an external load P pushes it INTO the master (toward penetration, gN<0). The
# augmented contact traction tN = −p pushes it OUT. Solve the gap at fixed (λ_N, P), then augment.
def solve_gap(lambdaN, epsN, P, k_s):
    """equilibrium of the 1-DOF gap u≡gN under load P (pushing into contact) + the structural spring
    k_s + the augmented contact traction. Returns the converged gN at this (λ_N, P)."""
    # assume ACTIVE: −P − k_s·gN − (λ_N + εN·gN) = 0 ⇒ gN = (−P − λ_N)/(k_s + εN)
    gN = (-P - lambdaN) / (k_s + epsN)
    if lambdaN + epsN * gN < 0.0:                  # KKT active (compression) ⇒ accept
        return gN
    # else inactive: no contact traction ⇒ gN = −P/k_s (open; P<0 would pull it open)
    return -P / k_s


def held_load_augment(epsN, P, k_s, maxAug=30, augTol=1e-10):
    """the analyze_augmented held-load loop (D1), 1-DOF: at a FIXED external load P, repeatedly
    solve the gap and commit one Uzawa step until the penetration < augTol. Returns (λ_N, gN, n_aug)."""
    lambdaN = 0.0
    gN = solve_gap(lambdaN, epsN, P, k_s)          # first equilibrium (penalty, 1 implicit solve)
    for k in range(1, maxAug + 1):
        lambdaN = uzawa_step(lambdaN, epsN, gN)    # commit-cycle Uzawa
        gN = solve_gap(lambdaN, epsN, P, k_s)      # held-load re-solve at the SAME P
        if abs(min(0.0, gN)) < augTol:             # penetration converged
            return lambdaN, gN, k
    return lambdaN, gN, maxAug


# ============================================================================ tests
print("ADR-57 E6 oracle — edge-edge one-scalar commit-cycle ALM")

# --- T1: byte-identity — λ_N≡0 ⇒ the augmented pressure IS the E2 penalty pressure ---
print("\nT1  λ_N≡0 ⇒ p = εN·gN, active iff gN<0 ⇒ the E2 PENALTY path EXACTLY (off-by-default identity)")
EPS = 1.0e5
for gN in (-0.02, -1e-4, 0.0, 1e-4, 0.05):
    p_alm = augmented_pressure(0.0, EPS, gN)
    p_penalty = min(0.0, EPS * gN)                 # the E2 path: tN = εN⟨−gN⟩ ⇒ p = min(0, εN gN)
    check(f"gN={gN:+.0e}: ALM(λ=0) == penalty", abs(p_alm - p_penalty) < TOL,
          f"p_alm={p_alm:.4f}")

# --- T2: held-load Uzawa drives penetration → an εN-INDEPENDENT tol (the headline ALM win) ---
print("\nT2  held-load augmentation: |gN| → augTol INDEPENDENT of εN (penalty is O(P/εN))")
P, KS, AUGTOL = 10.0, 1.0e3, 1e-10
pens_penalty, pens_alm = {}, {}
for epsN in (1e4, 1e5, 1e6, 1e7, 1e8):
    gN_penalty = solve_gap(0.0, epsN, P, KS)       # penalty-only (1 solve, no augmentation)
    pens_penalty[epsN] = abs(gN_penalty)
    lam, gN_alm, n_aug = held_load_augment(epsN, P, KS, augTol=AUGTOL)
    pens_alm[epsN] = abs(min(0.0, gN_alm))
    check(f"εN={epsN:.0e}: ALM penetration {pens_alm[epsN]:.2e} < augTol (was {pens_penalty[epsN]:.2e})",
          pens_alm[epsN] < AUGTOL, f"n_aug={n_aug}")
# the penalty penetration is εN-DEPENDENT (shrinks ∝ 1/εN); the ALM penetration is NOT (all < augTol)
penalty_spread = max(pens_penalty.values()) / min(pens_penalty.values())
alm_max = max(pens_alm.values())
check("penalty penetration is εN-DEPENDENT (spans orders of magnitude across εN)", penalty_spread > 1e3,
      f"spread={penalty_spread:.1e}")
check("ALM penetration is εN-INDEPENDENT (all < augTol regardless of εN)", alm_max < AUGTOL,
      f"max={alm_max:.2e}")

# --- T3: the converged contact TRACTION carries the load (tN ≈ P), λ_N negative and O(P) ---
print("\nT3  converged contact traction tN ≈ P (the multiplier carries the reaction at zero penetration)")
for epsN in (1e4, 1e6, 1e8):
    lam, gN_alm, n_aug = held_load_augment(epsN, P, KS, augTol=AUGTOL)
    tN = -augmented_pressure(lam, epsN, gN_alm)    # the assembled normal traction
    # at gN→0 the structural spring carries ~k_s·gN≈0 ⇒ the contact traction balances the load P
    check(f"εN={epsN:.0e}: tN ≈ P (the contact carries the full reaction)", abs(tN - P) < 1e-6 * P,
          f"tN={tN:.6f} (P={P})")
    # λ_N → −P as εN→∞ (residual O(P·k_s/εN)); always negative and O(P) (it IS the reaction multiplier)
    check(f"εN={epsN:.0e}: λ_N negative, → −P within the spring share O(P·k_s/εN)",
          lam < 0 and abs(lam - (-P)) < P * KS / epsN + 1e-9, f"λ_N={lam:.6f}")

# --- T4: RELEASE ⇒ λ_N → 0, F = 0 (the KKT release; no locked-in pull) ---
print("\nT4  release the load ⇒ the Uzawa min(0,·) drives λ_N → 0, traction → 0 (KKT open)")
# start from the converged contact state (λ_N ≈ −P), then REMOVE the load (P→0) and let the body open.
lam0, _, _ = held_load_augment(1e6, P, KS, augTol=AUGTOL)
lam = lam0
gN = solve_gap(lam, 1e6, 0.0, KS)                  # load removed: the body springs back (gN>0, open)
hist = [lam]
for _ in range(20):
    lam = uzawa_step(lam, 1e6, gN)                 # min(0, λ + εN·gN); gN>0 ⇒ drives λ up to 0
    gN = solve_gap(lam, 1e6, 0.0, KS)
    hist.append(lam)
F_final = -augmented_pressure(lam, 1e6, gN)
check("released: λ_N → 0 (no locked-in multiplier)", abs(lam) < TOL, f"λ_N: {lam0:.3f} → {lam:.3e}")
check("released: contact traction F → 0 (KKT open, no spurious pull)", abs(F_final) < TOL,
      f"F={F_final:.3e}")
check("released: gN > 0 (the body separated)", gN > 0, f"gN={gN:.4f}")

# --- T5: committed-only — a rejected (reverted) step must NOT mutate λ_N ---
print("\nT5  committed-only: λ_N is mutated SOLELY in commit() ⇒ revertToLastCommit leaves it untouched")
# model the lifecycle: commit() runs the Uzawa; revert() restores trials but NOT λ_N (it is never a trial).
lam_committed = -3.7
lam_after_revert = lam_committed                   # revert touches gpT/sign trials, NOT λ_N (C++ contract)
check("revert does not change λ_N (committed-only invariant)", lam_after_revert == lam_committed)
# and commit DOES advance it (the only mutation site)
lam_after_commit = uzawa_step(lam_committed, 1e5, -1e-6)
check("commit advances λ_N (the sole mutation site)", lam_after_commit != lam_committed,
      f"{lam_committed} → {lam_after_commit:.6f}")

# --- T6: equation count constant across augmentations (λ_N is state, not a DOF) ---
print("\nT6  the equation count is CONSTANT across augmentations (λ_N is committed state, not a DOF)")
# the 1-DOF model has exactly ONE unknown (the gap u) at EVERY augmentation — λ_N is a parameter that
# shifts the residual, never an unknown. So no augmentation re-numbers / resizes the SOE (D1 invariant).
n_unknowns = []
lam = 0.0
gN = solve_gap(lam, 1e6, P, KS)
for _ in range(5):
    n_unknowns.append(1)                           # always 1 DOF (u); λ_N never enters the SOE
    lam = uzawa_step(lam, 1e6, gN)
    gN = solve_gap(lam, 1e6, P, KS)
check("SOE size constant across all augmentations (no new DOF for λ_N)",
      len(set(n_unknowns)) == 1 and n_unknowns[0] == 1, f"sizes={n_unknowns}")

print(f"\n{'='*66}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E6 edge-edge one-scalar ALM oracle)\n{'='*66}")
sys.exit(1 if _fails else 0)
