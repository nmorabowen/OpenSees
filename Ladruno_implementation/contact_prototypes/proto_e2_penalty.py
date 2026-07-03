"""ADR-57 E2 ORACLE — edge-edge penalty NORMAL force + main tangent.

Oracle-first (the fork discipline): validate the penalty enforcement in pure numpy BEFORE
transcribing it into the LadrunoContactFE EDGE_EDGE mode. No OpenSees, no build. E2 turns
the E0 geometry (segment-segment closest point + common-perpendicular normal/gap + the
B-operator) into FORCE: the penalty normal traction, its self-equilibrating scatter, and the
main (penalty-Gram) tangent — the cos_t→0 contact the mortar clip degenerates on now produces
a restoring force.

What E2 owns (the ADR-57 E2 gate row + the 3-reviewer gate folds):
  - penalty normal traction tN = εN·⟨−gN⟩ (Macaulay; active only in penetration gN<0).
  - residual f = tN·B over the 4-node edge pair [a0,a1,b0,b1] (B from LadrunoEdgeKernel::bOperator
    = [(1−s)n | s n | −(1−t)n | −t n]) — the NTS B-operator shape with the master shape fns
    replaced by the edge linear weights (1−t,t) and the slave node by (1−s,s).
  - SELF-EQUILIBRIUM Σf = 0 (slave weights sum to 1, master weights sum to 1, opposite signs).
  - the force is ⟂ BOTH edges (along n = ê_s×ê_m) — the common-perpendicular line of action.
  - MAIN TANGENT K = εN·BᵀB (12×12, symmetric, rank-1, PSD) — the penalty-Gram / frozen-
    direction block (the geometric ∂n/∂u, ∂s/∂u, ∂t/∂u curvature terms are E4, gated off). FD-
    checked against the frozen-direction residual Jacobian (the εN B(∂gN/∂u)ᵀ form; ∂gN/∂u=B is
    the E0 envelope identity), and verified symmetric + PSD (eigenvalues ≥ 0).
  - penetration δ = P/εN: a unit edge-edge contact loaded by P penetrates δ = P/εN (the penalty
    compliance — the εN-dependence the E6 ALM removes).
  - MONOTONE gN THROUGH CONTACT driven through a REAL penalty step: as the slave edge is driven
    from gap>0 through 0 into penetration, the BODY-FIXED committed sign (E0 A-4 fix) keeps gN
    strictly DECREASING and the penalty force tN strictly INCREASING from 0 — no sign flip, no
    spurious release. The self-referential w·n rule MASKS the penetration (force stays 0).

Run: <pythoncore-3.12> proto_e2_penalty.py   (numpy only)
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

# status codes + gauges (mirror the E0 kernel / LadrunoEdgeKernel)
OK, PARALLEL, DEGENERATE = 0, -1, -2
TAU_PAR = np.sin(np.deg2rad(15.0)) ** 2
TAU_LEN_REL = 1e-12   # relative zero-edge floor (squared-length ratio; contact-review P5)
MARGIN = 1e-6


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# ============================================ E0 kernel (compact, mirrors LadrunoEdgeKernel)
def closest_pt_segseg(p1, q1, p2, q2, tau_par=TAU_PAR):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    out = {"status": OK, "s": 0.0, "t": 0.0, "c1": p1.copy(), "c2": p2.copy(),
           "denom": 0.0, "interior": False}
    len_ref = max(a, e)
    if len_ref <= 0.0 or a <= TAU_LEN_REL * len_ref or e <= TAU_LEN_REL * len_ref:
        out["status"] = DEGENERATE
        return out
    c, b, f = d1 @ r, d1 @ d2, d2 @ r
    denom = a * e - b * b
    out["denom"] = denom
    s = np.clip((b * f - c * e) / denom, 0.0, 1.0) if denom > 1e-300 else 0.0
    t = (b * s + f) / e
    if t < 0.0:
        t, s = 0.0, np.clip(-c / a, 0.0, 1.0)
    elif t > 1.0:
        t, s = 1.0, np.clip((b - c) / a, 0.0, 1.0)
    out["s"], out["t"] = s, t
    out["c1"], out["c2"] = p1 + s * d1, p2 + t * d2
    out["interior"] = (MARGIN <= s <= 1.0 - MARGIN) and (MARGIN <= t <= 1.0 - MARGIN)
    if denom < tau_par * a * e:
        out["status"] = PARALLEL
    return out


def edge_normal(d1, d2):
    m = np.cross(d1, d2)
    nrm = np.linalg.norm(m)
    return None if nrm < 1e-300 else m / nrm


def edge_gap(c1, c2, d1, d2, sign_ref=None, committed_sign=0):
    """signed gap gN = w.n, w = c1 - c2 (slave - master); body-fixed sign (E0 A-4).
    committed_sign in {+1,-1} ⇒ verbatim; else sign_ref ⇒ orient n.sign_ref ≥ 0 (return it);
    else the buggy self-referential w.n rule (kept only to DEMO the masking)."""
    n = edge_normal(d1, d2)
    if n is None:
        return 0.0, None, 0
    w = c1 - c2
    if committed_sign in (1, -1):
        sign = committed_sign
    elif sign_ref is not None:
        sign = 1 if (n @ sign_ref) >= 0.0 else -1
    else:
        sign = 1 if (n @ w) >= 0.0 else -1
    n = sign * n
    return w @ n, n, sign


def b_operator(n, s, t):
    """first variation dgN = B.[dp1,dq1,dp2,dq2]; rows = [a0,a1,b0,b1] (slave edge, master edge)."""
    return np.array([(1 - s) * n, s * n, -(1 - t) * n, -t * n])


# ============================================================ the E2 penalty layer
def penalty_residual(p1, q1, p2, q2, epsN, sign_ref=None, committed_sign=0):
    """edge-edge penalty: tN = εN⟨−gN⟩, residual f = tN·B over [a0,a1,b0,b1] (12-vec).
    returns dict: active, gN, tN, n, s, t, B(4x3), f(12). f≡0 unless penetrating (gN<0)."""
    k = closest_pt_segseg(p1, q1, p2, q2)
    out = {"active": False, "gN": 0.0, "tN": 0.0, "n": None, "s": k["s"], "t": k["t"],
           "B": np.zeros((4, 3)), "f": np.zeros(12), "status": k["status"], "interior": k["interior"]}
    if k["status"] != OK or not k["interior"]:
        return out                                   # parallel/degenerate/near-vertex ⇒ inert
    d1, d2 = q1 - p1, q2 - p2
    gN, n, sgn = edge_gap(k["c1"], k["c2"], d1, d2, sign_ref, committed_sign)
    out["gN"], out["n"], out["sign"] = gN, n, sgn
    B = b_operator(n, k["s"], k["t"])
    out["B"] = B
    if gN < 0.0:                                     # Macaulay ⟨−gN⟩: active only in penetration
        tN = -epsN * gN                              # = εN|gN| > 0
        out["active"], out["tN"] = True, tN
        out["f"] = (tN * B).reshape(-1)              # f[3i+d] = tN·B[i][d]
    return out


def penalty_main_tangent(B, epsN):
    """main (penalty-Gram) tangent K = εN·BᵀB (12×12). Brow = B flattened to the 1×12 gap operator."""
    Brow = B.reshape(-1)
    return epsN * np.outer(Brow, Brow)


# ============================================================================ tests
print("ADR-57 E2 oracle — edge-edge penalty normal force + main tangent")

# canonical two crossed bars, edge-on: slave edge ∥ x at height zc, master edge ∥ y at z=0.
# they cross over the origin; zc>0 ⇒ separated, zc<0 ⇒ penetration. orientDir = +z (toward slave).
def crossed_bars(zc):
    a0, a1 = np.array([-1, 0, zc]), np.array([1, 0, zc])      # slave edge (∥ x)
    b0, b1 = np.array([0, -1, 0.]), np.array([0, 1, 0.])      # master edge (∥ y)
    return a0, a1, b0, b1
ORIENT = np.array([0, 0, 1.])     # toward the slave's allowed half-space
EPS = 1.0e3

# --- T1: penetrating crossed bars -> active penalty, force along n, ⟂ both edges ---
print("\nT1  penetrating crossed bars: tN=εN⟨−gN⟩, force along n ⟂ both edges")
a0, a1, b0, b1 = crossed_bars(-0.02)
res = penalty_residual(a0, a1, b0, b1, EPS, sign_ref=ORIENT)
check("penetration ⇒ active", res["active"], f"gN={res['gN']:+.4f} tN={res['tN']:.4f}")
check("gN < 0 (penetration registered, body-fixed sign)", res["gN"] < 0)
check("tN == εN·|gN|", abs(res["tN"] - EPS * abs(res["gN"])) < TOL)
# slave total force = +tN n, master total = −tN n (the 4-node scatter)
f = res["f"].reshape(4, 3)
fs = f[0] + f[1]                                  # slave edge nodes a0,a1
fm = f[2] + f[3]                                  # master edge nodes b0,b1
n = res["n"]
check("slave total force ∥ +n", np.allclose(fs, res["tN"] * n, atol=TOL),
      f"fs={fs} tN·n={res['tN']*n}")
d1, d2 = a1 - a0, b1 - b0
check("force ⟂ slave edge (f_s·ê_s == 0)", abs(fs @ d1) < TOL)
check("force ⟂ master edge (f_s·ê_m == 0)", abs(fs @ d2) < TOL)

# --- T2: self-equilibrium Σf = 0 (no net force/moment injected) ---
print("\nT2  self-equilibrium: Σf = 0 over the 4 edge nodes")
check("Σf == 0", np.allclose(f.sum(axis=0), 0.0, atol=TOL), f"Σf={f.sum(axis=0)}")
check("slave + master totals cancel", np.allclose(fs + fm, 0.0, atol=TOL))

# --- T3: separated bars (gN>0) -> ZERO force (Macaulay inactive), no spurious traction ---
print("\nT3  separated bars (gN>0) -> inert (Macaulay ⟨−gN⟩ = 0)")
a0, a1, b0, b1 = crossed_bars(+0.05)
resS = penalty_residual(a0, a1, b0, b1, EPS, sign_ref=ORIENT)
check("separated ⇒ NOT active", not resS["active"], f"gN={resS['gN']:+.4f}")
check("separated ⇒ f == 0", np.allclose(resS["f"], 0.0, atol=TOL))
check("separated gN > 0 (body-fixed sign, open)", resS["gN"] > 0)

# --- T4: main tangent K = εN BᵀB — symmetric, PSD, rank-1 ---
print("\nT4  main tangent K = εN BᵀB: symmetric, PSD, rank-1 (the gap mode)")
a0, a1, b0, b1 = crossed_bars(-0.02)
res = penalty_residual(a0, a1, b0, b1, EPS, sign_ref=ORIENT)
K = penalty_main_tangent(res["B"], EPS)
check("K symmetric", np.allclose(K, K.T, atol=TOL), f"max asym={np.abs(K-K.T).max():.2e}")
evals = np.linalg.eigvalsh(0.5 * (K + K.T))
check("K PSD (all eigenvalues ≥ 0)", evals.min() > -1e-9, f"min eig={evals.min():.3e}")
check("K rank-1 (single nonzero eigenvalue = the gap mode)",
      np.sum(evals > 1e-6 * evals.max()) == 1, f"eigs>0: {int(np.sum(evals>1e-6*evals.max()))}")

# --- T5: FD-check the main (frozen-direction) tangent vs the residual Jacobian ---
print("\nT5  K == εN BᵀB FD-checked vs the frozen-direction residual Jacobian")
# the MAIN term is the penalty-Gram block: freeze n,s,t at the base config and differentiate
# the residual f = tN·B_base with gN_frozen(u) = n_base·(Ps(u)−Pm(u)) (s,t,n held). Then
# ∂f/∂u = εN B_base (∂gN_frozen/∂u)ᵀ = εN B_baseᵀB_base (since ∂gN_frozen/∂u = B_base, the E0
# envelope identity). The geometric ∂{n,s,t}/∂u terms are E4 (gated off) — excluded here.
base = penalty_residual(a0, a1, b0, b1, EPS, sign_ref=ORIENT)
nB, sB, tB, Bb = base["n"], base["s"], base["t"], base["B"]
X0 = np.concatenate([a0, a1, b0, b1]).astype(float)    # 12 dof
def f_frozen(X):
    P = X.reshape(4, 3)
    Ps = (1 - sB) * P[0] + sB * P[1]                   # closest pt on slave edge (frozen s)
    Pm = (1 - tB) * P[2] + tB * P[3]                   # closest pt on master edge (frozen t)
    gN = nB @ (Ps - Pm)                                # frozen-direction gap
    tN = -EPS * gN if gN < 0 else 0.0
    return (tN * Bb).reshape(-1)
h = 1e-7
Kfd = np.zeros((12, 12))
for j in range(12):
    Xp = X0.copy(); Xp[j] += h
    Xm = X0.copy(); Xm[j] -= h
    Kfd[:, j] = (f_frozen(Xp) - f_frozen(Xm)) / (2 * h)
# residual r = tN·B = −f_int ⇒ the stiffness K = −∂r/∂u (the solver tangent), == εN BᵀB
Kana = penalty_main_tangent(Bb, EPS)
check("max |(−∂f/∂u) − εN BᵀB| < 1e-3", np.abs(-Kfd - Kana).max() < 1e-3,
      f"max err={np.abs(-Kfd - Kana).max():.2e}")

# --- T6: penetration δ = P/εN (the penalty compliance) ---
print("\nT6  penetration δ = P/εN under a normal load P (penalty compliance)")
# single-dof model: drive the slave edge straight down by δ (rigid translation along −z), master
# fixed. Equilibrium of the gap dof: penalty force tN = εN·δ balances the applied P ⇒ δ = P/εN.
P = 5.0
for eps in (1e2, 1e3, 1e4):
    # solve εN·δ = P ⇒ δ = P/εN; verify the penalty residual reproduces it at that penetration
    delta = P / eps
    a0, a1, b0, b1 = crossed_bars(-delta)
    r = penalty_residual(a0, a1, b0, b1, eps, sign_ref=ORIENT)
    # the total normal reaction transmitted = tN (= the slave total |force|)
    check(f"εN={eps:.0e}: tN == P (δ=P/εN={delta:.2e})", abs(r["tN"] - P) < 1e-6 * P,
          f"tN={r['tN']:.6f}")

# --- T7: monotone gN through contact (committed sign) driven through a REAL penalty step ---
print("\nT7  committed sign -> gN monotone, force monotone-from-0; w·n rule MASKS penetration")
zs = np.linspace(0.06, -0.06, 25)                  # drive slave down through the master
# capture the body-fixed sign ONCE at first contact engagement (gap just > 0), ref = +z (toward slave)
a0, a1, b0, b1 = crossed_bars(zs[0])
k0 = closest_pt_segseg(a0, a1, b0, b1)
_, _, SGN = edge_gap(k0["c1"], k0["c2"], a1 - a0, b1 - b0, sign_ref=ORIENT)
gN_hist, tN_hist, gN_buggy = [], [], []
for z in zs:
    a0, a1, b0, b1 = crossed_bars(z)
    r = penalty_residual(a0, a1, b0, b1, EPS, committed_sign=SGN)   # body-fixed committed sign
    gN_hist.append(r["gN"]); tN_hist.append(r["tN"])
    rb = penalty_residual(a0, a1, b0, b1, EPS)                      # buggy self-referential w·n rule
    gN_buggy.append(rb["gN"])
gN_hist, tN_hist, gN_buggy = np.array(gN_hist), np.array(tN_hist), np.array(gN_buggy)
check("committed sign: gN strictly decreasing through contact (no flip)",
      np.all(np.diff(gN_hist) < 1e-12), f"gN: {gN_hist[0]:+.3f} -> {gN_hist[-1]:+.3f}")
check("committed sign: gN goes NEGATIVE (penetration registered)", gN_hist[-1] < 0)
# the penalty force is 0 while open, then rises MONOTONICALLY once penetrating — no spurious release
active = tN_hist > 0
check("force == 0 while gN ≥ 0 (Macaulay inactive when open)",
      np.all(tN_hist[gN_hist >= 0] == 0.0))
check("force rises monotonically once penetrating (no chatter/release)",
      np.all(np.diff(tN_hist[active]) > -1e-12) and active.sum() >= 2,
      f"tN_max={tN_hist.max():.3f}, active steps={int(active.sum())}")
check("buggy w·n rule: gN stays ≥ 0 (penetration MASKED — the A-4 bug)",
      np.all(gN_buggy >= -1e-9), f"buggy min gN={gN_buggy.min():+.4f}")

print(f"\n{'='*66}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E2 edge-edge penalty force + main tangent oracle)\n{'='*66}")
sys.exit(1 if _fails else 0)
