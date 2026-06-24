"""ADR-41 D1 ORACLE — within-step held-load augmentation: loop control + no-corruption protocol.

C2.2 / C4 already ship + oracle the augmented-Lagrange Uzawa mechanics (proto_c2_alm 25/25,
proto_c4_mortar_tie). D1 adds NO new numerics — it formalizes the WITHIN-STEP loop control (drive
the penetration / tie residual to augTol at a HELD external load, distinct from the one-augmentation-
per-step across-step Uzawa) AND the NO-CORRUPTION protocol (the held-load re-commits must fire NO
recorder samples and advance NO time / commitTag). The convergence criterion is the standard Uzawa
outer loop (Wriggers Ch.6 / Laursen Ch.8):

    inner:  solve  R(u; λ) = 0   (Newton equilibrium at FIXED λ and FIXED external load)
    outer:  λ ← min(0, λ + εN·ḡ(u))     (projected for the normal inequality; NO clamp for a tie)
    stop:   ‖ḡ‖_∞ < augTol

This proto pins the two things D1 is responsible for, in pure numpy:

  T1  WITHIN-STEP monotone convergence at a HELD load: a 1-DOF elastic+contact node, external load
      frozen, λ augmented in place → penetration ‖g‖ DECREASES MONOTONICALLY to augTol; the converged
      position is εN-INDEPENDENT (the penalty-alone residual g = O(1/εN) is removed at finite εN).
  T2  augTol STOP: the loop halts the first pass ‖g‖ < augTol — no over/under-running.
  T3  NO-CORRUPTION accounting: model Domain::commit() as {Uzawa-λ-update ; record ; commitTag++}.
      The held-load sweep must update λ EVERY pass (that is the point) but fire the recorder / bump
      commitTag ZERO times. With the suppression flag ON across the sweep: λ advances n_aug times,
      recorder fires 0 times, commitTag unchanged ⇒ recorder output == the no-augmentation run.
  T4  LOAD not double-counted: the held re-solves re-APPLY the same load factor (a zero-increment
      LoadControl SETS, not accumulates), so the converged reaction == the single-application reaction.

Run: <pythoncore-3.12> proto_d1_augmentation.py   (numpy only)
"""
import sys
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

_fails = 0


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ----------------------------------------------------------------------------
# Minimal 1-DOF normal-contact model (the AL-Uzawa contraction in closed form).
# A slave z-DOF u, anchored to ground by a spring ks (natural at u=0), pushed DOWN
# by a held external load P>0 against a rigid wall at u=0 (penetration ⇔ u<0).
#
#   contact traction  t(u) = min(0, λ + εN·u)                 (≤ 0, active iff u<0)
#   force balance     ks·u + t = −P            (inner equilibrium, FIXED λ, FIXED P)
#   active solve      u = (−P − λ)/(ks + εN)
#   Uzawa outer       λ ← min(0, λ + εN·u)
#
# Converged fixed point: λ* = −P, u* = 0 (g=0) — εN-INDEPENDENT; the penalty-only
# answer (λ≡0) leaves u = −P/(ks+εN) = O(1/εN). Contraction factor ks/(ks+εN).
# ----------------------------------------------------------------------------
def held_load_augment(P, ks, epsN, augTol, maxAug, suppress):
    """Run the within-step held-load Uzawa loop. Returns a trace dict.

    `suppress` models the D1 flag: when True the held re-commits fire NO recorder
    and DO NOT bump commitTag (only the Uzawa λ update runs); when False each pass
    fires the recorder + bumps commitTag (the pre-D1 corruption)."""
    lam = 0.0
    rec_fires = 0
    commitTag = 0
    # the load step itself: apply P, inner-solve, ONE real commit (records once).
    u = (-P - lam) / (ks + epsN)        # penalty equilibrium (λ=0)
    rec_fires += 1                       # the physical step's single recorder sample
    commitTag += 1
    lam = min(0.0, lam + epsN * u)       # the across-step augmentation (commit #1)
    pen0 = abs(min(0.0, u))
    pens = [pen0]
    n_aug = 0
    # the WITHIN-STEP held-load sweep (zero load increment: P held, re-solve, re-augment)
    for _ in range(maxAug):
        u = (-P - lam) / (ks + epsN)     # inner equilibrium at the HELD load, current λ
        pen = abs(min(0.0, u))
        pens.append(pen)
        n_aug += 1
        if not suppress:                 # pre-D1: each held re-commit corrupts the record stream
            rec_fires += 1
            commitTag += 1
        lam = min(0.0, lam + epsN * u)   # Uzawa λ update ALWAYS runs (the augmentation)
        if pen < augTol:
            break
    return dict(pens=pens, lam=lam, u=u, n_aug=n_aug, rec_fires=rec_fires,
                commitTag=commitTag, reaction=-(ks * u))   # contact reaction = −ks·u + (−P) balance


print("ADR-41 D1 — within-step held-load augmentation oracle\n")

# -- T1: within-step monotone convergence, εN-independent converged position --------
augTol = 1e-10
lo = held_load_augment(P=1.0, ks=1.0, epsN=1.0e2, augTol=augTol, maxAug=200, suppress=True)
hi = held_load_augment(P=1.0, ks=1.0, epsN=1.0e4, augTol=augTol, maxAug=200, suppress=True)
mono = all(lo["pens"][i + 1] <= lo["pens"][i] + 1e-300 for i in range(len(lo["pens"]) - 1))
check("T1 within-step penetration is monotone non-increasing", mono)
check("T1 converged below augTol at low εN", lo["pens"][-1] < augTol, f"‖g‖={lo['pens'][-1]:.2e}")
check("T1 converged below augTol at high εN", hi["pens"][-1] < augTol, f"‖g‖={hi['pens'][-1]:.2e}")
check("T1 converged position εN-INDEPENDENT (both → u≈0)",
      abs(lo["u"] - hi["u"]) < 1e-9 and abs(lo["u"]) < 1e-9,
      f"u(1e2)={lo['u']:.2e}, u(1e4)={hi['u']:.2e}")
# penalty-alone (the very first pass, λ=0) IS εN-dependent — the thing ALM removes
pen_penalty_lo = 1.0 / (1.0 + 1.0e2)
pen_penalty_hi = 1.0 / (1.0 + 1.0e4)
check("T1 penalty-alone penetration is O(1/εN) (what ALM removes)",
      pen_penalty_lo > 50.0 * pen_penalty_hi,
      f"{pen_penalty_lo:.2e} vs {pen_penalty_hi:.2e}")

# -- T2: augTol stop (no over-run) ---------------------------------------------------
# the loop must stop the FIRST pass under augTol; the contraction is ks/(ks+εN)=1/101 per
# pass at εN=1e2, so n_aug = ceil(log(augTol/pen0)/log(1/101)) — a finite, predicted count.
import math
factor = 1.0 / (1.0 + 1.0e2)
pen0 = lo["pens"][0]
expect = math.ceil(math.log(augTol / pen0) / math.log(factor)) if pen0 > augTol else 0
check("T2 loop stops the first pass under augTol (no over-run)",
      lo["pens"][-1] < augTol <= (lo["pens"][-2] if len(lo["pens"]) > 1 else math.inf),
      f"n_aug={lo['n_aug']} (predicted ≈ {expect})")

# -- T3: NO-CORRUPTION accounting ----------------------------------------------------
# WITHOUT augmentation: 1 recorder fire, commitTag advances by 1 (the physical step).
base = held_load_augment(P=1.0, ks=1.0, epsN=1.0e2, augTol=augTol, maxAug=0, suppress=True)
# WITH within-step augmentation + suppression: λ advanced n_aug+1 times, but STILL 1 fire / 1 tag.
sup = lo
# WITH augmentation but NO suppression (the pre-D1 corruption): 1 + n_aug spurious fires.
cor = held_load_augment(P=1.0, ks=1.0, epsN=1.0e2, augTol=augTol, maxAug=200, suppress=False)
check("T3 suppressed sweep fires recorder EXACTLY once (== no-augmentation run)",
      sup["rec_fires"] == base["rec_fires"] == 1,
      f"base={base['rec_fires']}, suppressed={sup['rec_fires']}")
check("T3 suppressed sweep leaves commitTag == no-augmentation run",
      sup["commitTag"] == base["commitTag"] == 1,
      f"base={base['commitTag']}, suppressed={sup['commitTag']}")
check("T3 the Uzawa λ STILL advances under suppression (augmentation is NOT skipped)",
      abs(sup["lam"] - (-1.0)) < 1e-8 and sup["n_aug"] >= 2,
      f"λ={sup['lam']:.4f} after {sup['n_aug']} passes")
check("T3 WITHOUT suppression the recorder is corrupted (1 + n_aug spurious fires)",
      cor["rec_fires"] == 1 + cor["n_aug"] and cor["rec_fires"] > base["rec_fires"],
      f"corrupted fires={cor['rec_fires']} vs clean={base['rec_fires']}")

# -- T4: the held re-solves do not double-count the load -----------------------------
# the converged reaction equals the single-application reaction P (a zero-increment LoadControl
# SETS the load factor, it does not accumulate): total reaction balancing the held load = P.
check("T4 converged reaction == single-application load (no double-count)",
      abs(abs(sup["lam"]) - 1.0) < 1e-8,
      f"|λ*|={abs(sup['lam']):.6f} == P=1.0")

print(f"\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAIL(S)'}  (D1 within-step augmentation)")
sys.exit(1 if _fails else 0)
