"""ADR-57 E0 ORACLE — LadrunoEdgeKernel: segment-segment closest point + edge normal/gap.

Oracle-first (the fork discipline): validate the NET-NEW edge-edge geometry in pure
numpy BEFORE transcribing it into SRC/domain/contact/LadrunoEdgeKernel.h. No OpenSees,
no build. This is the geometry the perpendicular edge-edge narrow phase (the cos_t->0
case the mortar face-clip degenerates on) stands on.

What E0 owns (and the falsifiers, per the ADR-57 E0 gate row + the 3-reviewer gate folds):
  - segment-segment closest point via ERICSON'S EXACT BRANCH LADDER (RTCD 5.1.9
    ClosestPtSegmentSegment) — the conditional clamp ladder, NOT a single back-substitution
    (Lens-A MAJOR A-1). Validated vs an INDEPENDENT 1D-scan+refine reference across all 8
    Voronoi sub-regions, incl. the vertex-vertex double-clamp.
  - dual-purpose degeneracy gauge denom = a*e - b*b = ||d1||^2 ||d2||^2 sin^2(theta_edge):
    governs BOTH closest-point uniqueness AND normal conditioning (Lens-A A-3). The
    parallel/ill-conditioned REFUSAL is denom < tau_par*a*e <=> sin^2 theta < tau_par;
    tau_par set by CONDITIONING (default sin^2(15deg) ~ 0.067), not the draft's 5deg.
  - zero-length / coincident-node edge REFUSAL (Lens-A NIT; real once elements are removed).
  - common-perpendicular normal n = (d1 x d2)/||.|| with a BODY-FIXED (not w-derived) sign
    anchor (Lens-A A-4: the self-referential w.n rule masks interpenetration through contact).
  - edge gap gN = w.n; the identity gN = +-||w|| holds ONLY in the margin-interior set
    (Lens-A A-2) — demonstrated to FAIL at a clamped near-vertex.
  - the first-variation B-operator B = [(1-s)n | s n | -(1-t)n | -t n] is the EXACT dgN/du
    (envelope theorem: w _|_ tangents and w || n kill the ds/dt and dn terms) — FD-checked.

CONDITIONING NOTE (caught oracle-first; fold into the ADR wording in the E0 PR): the
NORMAL DIRECTION dn/du scales like 1/sin(theta) (normalization of d1 x d2), while the
closest-point PARAMETERS d(s,t)/du scale like 1/sin^2(theta) (the denom in the solve) —
the latter is the dominant gauge and is what the ADR's "1/sin^2" refers to; both are
bounded at the admitted band (theta >= 15deg) and the tau_par gate refuses below.

Run: <pythoncore-3.12> proto_e0_closest_point.py   (numpy only)
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

# status codes (mirror the future LadrunoEdgeKernel)
OK, PARALLEL, DEGENERATE = 0, -1, -2
TAU_PAR = np.sin(np.deg2rad(15.0)) ** 2     # conditioning-justified parallel gate ~0.0670
TAU_LEN = 1e-9                              # relative zero-length floor (on squared length)
MARGIN = 1e-6                               # strict-interior margin delta (parametric)


def check(name, ok, extra=""):
    global _fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not ok:
        _fails += 1


# =================================================== the kernel (mirror C++ to come)
def closest_pt_segseg(p1, q1, p2, q2, tau_par=TAU_PAR):
    """Ericson RTCD 5.1.9 ClosestPtSegmentSegment + a fork parallel/degenerate refusal.

    seg1 = p1->q1 (slave edge), seg2 = p2->q2 (master edge).
    returns dict: status, s, t, c1, c2, denom, a, e, interior(bool, margin'd).
    On PARALLEL we STILL return the Ericson s=0 fallback (s,t,c1,c2) so the handler can
    run the clip cross-check before dropping the pair (Lens-C: no silent lost contact).
    """
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = d1 @ d1            # ||seg1||^2 >= 0
    e = d2 @ d2            # ||seg2||^2 >= 0
    out = {"status": OK, "s": 0.0, "t": 0.0, "c1": p1.copy(), "c2": p2.copy(),
           "denom": 0.0, "a": a, "e": e, "interior": False}

    # --- zero-length / coincident-node guard (Lens-A NIT) ---
    if a <= TAU_LEN or e <= TAU_LEN:
        out["status"] = DEGENERATE
        return out

    c = d1 @ r
    b = d1 @ d2
    f = d2 @ r
    denom = a * e - b * b          # = a*e*sin^2(theta_edge) >= 0
    out["denom"] = denom

    # --- Ericson's EXACT branch ladder (NOT a single pass) ---
    if denom > 1e-300:
        s = np.clip((b * f - c * e) / denom, 0.0, 1.0)
    else:
        s = 0.0                    # parallel: pick an arbitrary point on seg1 (Ericson)
    t = (b * s + f) / e
    if t < 0.0:
        t = 0.0
        s = np.clip(-c / a, 0.0, 1.0)
    elif t > 1.0:
        t = 1.0
        s = np.clip((b - c) / a, 0.0, 1.0)
    out["s"], out["t"] = s, t
    out["c1"] = p1 + s * d1
    out["c2"] = p2 + t * d2
    out["interior"] = (MARGIN <= s <= 1.0 - MARGIN) and (MARGIN <= t <= 1.0 - MARGIN)

    # --- parallel / ill-conditioned-normal refusal (Lens-A A-3) ---
    #     denom < tau_par*a*e  <=>  sin^2(theta_edge) < tau_par
    if denom < tau_par * a * e:
        out["status"] = PARALLEL
    return out


def edge_normal(d1, d2):
    """common-perpendicular UNIT normal n = (d1 x d2)/||.||; None if degenerate."""
    m = np.cross(d1, d2)
    nrm = np.linalg.norm(m)
    if nrm < 1e-300:
        return None
    return m / nrm


def edge_gap(p1, q1, p2, q2, s, t, sign_ref=None, committed_sign=0):
    """signed gap gN = w.n with w = c1 - c2 (slave - master).

    sign discipline (Lens-A A-4):
      - committed_sign in {+1,-1}: BODY-FIXED, applied verbatim (the shipped behaviour).
      - else sign_ref given: orient n so n.sign_ref > 0 and RETURN that sign to commit.
      - else (the BUGGY draft rule, kept only to demonstrate the failure): n.w > 0.
    returns (gN, n, sign).
    """
    d1, d2 = q1 - p1, q2 - p2
    c1, c2 = p1 + s * d1, p2 + t * d2
    n = edge_normal(d1, d2)
    if n is None:
        return 0.0, None, 0
    w = c1 - c2
    if committed_sign in (1, -1):
        n = committed_sign * n
        sign = committed_sign
    elif sign_ref is not None:
        sign = 1 if (n @ sign_ref) >= 0.0 else -1
        n = sign * n
    else:
        sign = 1 if (n @ w) >= 0.0 else -1      # the self-referential rule (BUG demo)
        n = sign * n
    return w @ n, n, sign


def b_operator(n, s, t):
    """first variation dgN = B . [dp1,dq1,dp2,dq2] (each a 3-vec). Returns 4x3."""
    return np.array([(1 - s) * n, s * n, -(1 - t) * n, -t * n])


# ------------------------------------------ independent reference (different algorithm)
def closest_ref(p1, q1, p2, q2, ns=200001):
    """1D scan over s in [0,1] with the inner closed-form optimal t, then a local refine.
    Independent of Ericson's ladder (validates it). Returns (s,t,dist,c1,c2)."""
    d1, d2 = q1 - p1, q2 - p2
    e = d2 @ d2

    def best_t(s):
        c1 = p1 + s * d1
        return float(np.clip(((c1 - p2) @ d2) / e, 0.0, 1.0)) if e > 0 else 0.0

    def dist2(s):
        t = best_t(s)
        w = (p1 + s * d1) - (p2 + t * d2)
        return w @ w

    ss = np.linspace(0.0, 1.0, ns)
    d2s = np.array([dist2(s) for s in ss])
    i = int(np.argmin(d2s))
    lo, hi = ss[max(0, i - 1)], ss[min(ns - 1, i + 1)]
    for _ in range(80):                      # golden-section refine on [lo,hi]
        m1 = lo + 0.382 * (hi - lo)
        m2 = lo + 0.618 * (hi - lo)
        if dist2(m1) < dist2(m2):
            hi = m2
        else:
            lo = m1
    s = 0.5 * (lo + hi)
    t = best_t(s)
    c1, c2 = p1 + s * d1, p2 + t * d2
    return s, t, float(np.linalg.norm(c1 - c2)), c1, c2


def rot(axis, ang):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    c, s = np.cos(ang), np.sin(ang)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) * c + s * K + (1 - c) * np.outer(a, a)


# ============================================================================ tests
print("ADR-57 E0 oracle — LadrunoEdgeKernel: segment-segment closest point + normal/gap")

# --- T1: Ericson ladder vs the independent reference across the 8 Voronoi sub-regions ---
print("\nT1  closest point vs independent scan reference (all sub-regions)")
cases = {
    "interior crossing (skew)":
        (np.array([-1, 0, 0.]), np.array([1, 0, 0.]),
         np.array([0, -1, 1.]), np.array([0, 1, 1.])),
    "edge-region (one clamped: seg1 vertex vs seg2 interior)":
        (np.array([0, 0, 0.]), np.array([1, 0, 0.]),
         np.array([2, -1, 1.]), np.array([2, 1, 1.])),
    "vertex-vertex DOUBLE-CLAMP (both endpoints)":
        (np.array([0, 0, 0.]), np.array([1, 0, 0.]),
         np.array([2, 0, 1.]), np.array([2, 1, 1.])),
    "tilted skew (rotated interior crossing)":
        (rot([0.3, 1, 0.2], 0.6) @ np.array([-1, 0, 0.]),
         rot([0.3, 1, 0.2], 0.6) @ np.array([1, 0, 0.]),
         rot([0.3, 1, 0.2], 0.6) @ np.array([0.2, -1, 1.]),
         rot([0.3, 1, 0.2], 0.6) @ np.array([0.2, 1, 1.])),
}
for name, (p1, q1, p2, q2) in cases.items():
    k = closest_pt_segseg(p1, q1, p2, q2)
    rs, rt, rd, rc1, rc2 = closest_ref(p1, q1, p2, q2)
    kd = float(np.linalg.norm(k["c1"] - k["c2"]))
    check(f"{name}: dist matches ref", abs(kd - rd) < 1e-5, f"ker={kd:.7f} ref={rd:.7f}")
    check(f"{name}: (s,t) match ref", abs(k["s"] - rs) < 1e-4 and abs(k["t"] - rt) < 1e-4,
          f"ker=({k['s']:.5f},{k['t']:.5f}) ref=({rs:.5f},{rt:.5f})")

# explicit closed-form check on the vertex-vertex case (s=1 clamp, t=0 clamp, dist=sqrt2)
kvv = closest_pt_segseg(np.array([0, 0, 0.]), np.array([1, 0, 0.]),
                        np.array([2, 0, 1.]), np.array([2, 1, 1.]))
check("vertex-vertex: s==1 (clamped), t==0 (clamped)",
      abs(kvv["s"] - 1.0) < TOL and abs(kvv["t"]) < TOL, f"s={kvv['s']} t={kvv['t']}")
check("vertex-vertex: dist == sqrt(2)", abs(np.linalg.norm(kvv["c1"] - kvv["c2"]) - np.sqrt(2)) < TOL)
check("vertex-vertex: NOT margin-interior", not kvv["interior"])

# --- T2: parallel + near-parallel-edge-on REFUSAL; fallback point still computable ---
print("\nT2  parallel / near-parallel refusal (denom < tau_par), fallback point usable")
# exactly parallel (both along x, offset in z)
kp = closest_pt_segseg(np.array([0, 0, 0.]), np.array([1, 0, 0.]),
                       np.array([0, 0, 1.]), np.array([1, 0, 1.]))
check("exactly parallel -> PARALLEL", kp["status"] == PARALLEL, f"status={kp['status']}")
check("parallel: fallback (s,t,c1,c2) still finite",
      np.all(np.isfinite(kp["c1"])) and np.all(np.isfinite(kp["c2"])))
# near-parallel: 4deg between edges -> below the 15deg gate -> refused
ang = np.deg2rad(4.0)
kn = closest_pt_segseg(np.array([-1, 0, 0.]), np.array([1, 0, 0.]),
                       np.array([-1, 0, 0.5]) + np.array([0, 0, 0]),
                       np.array([1, 0, 0.5]) + np.array([0, np.sin(ang) * 2, 0]) - np.array([0, np.sin(ang), 0]))
# build cleanly instead: seg2 dir = (cos4,sin4,0), lifted z=0.5, centered
d2dir = np.array([np.cos(ang), np.sin(ang), 0.0])
kn = closest_pt_segseg(np.array([-1, 0, 0.]), np.array([1, 0, 0.]),
                       np.array([0, 0, 0.5]) - d2dir, np.array([0, 0, 0.5]) + d2dir)
check("4deg near-parallel (< 15deg gate) -> PARALLEL", kn["status"] == PARALLEL,
      f"status={kn['status']} sin^2={kn['denom']/(kn['a']*kn['e']):.5f} tau={TAU_PAR:.5f}")
# 30deg: comfortably admitted
d2dir = np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0.0])
k30 = closest_pt_segseg(np.array([-1, 0, 0.]), np.array([1, 0, 0.]),
                        np.array([0, 0, 0.5]) - d2dir, np.array([0, 0, 0.5]) + d2dir)
check("30deg -> OK (admitted), margin-interior", k30["status"] == OK and k30["interior"])

# --- T3: zero-length / coincident-node edge -> DEGENERATE ---
print("\nT3  zero-length / coincident-node edge -> DEGENERATE")
kz = closest_pt_segseg(np.array([1, 1, 1.]), np.array([1, 1, 1.]),     # collapsed seg1
                       np.array([0, -1, 0.]), np.array([0, 1, 0.]))
check("collapsed seg1 -> DEGENERATE", kz["status"] == DEGENERATE, f"status={kz['status']}")
kz2 = closest_pt_segseg(np.array([0, 0, 0.]), np.array([1, 0, 0.]),
                        np.array([2, 2, 2.]), np.array([2, 2, 2.]))   # collapsed seg2
check("collapsed seg2 -> DEGENERATE", kz2["status"] == DEGENERATE)

# --- T4: conditioning — d(s,t)/du ~ 1/sin^2, dn/du ~ 1/sin; bounded at 15deg gate ---
print("\nT4  normal/param conditioning across the angle band (gauge = 1/sin(theta))")
def stn_sensitivity(theta_deg, h=1e-7):
    th = np.deg2rad(theta_deg)
    dd = np.array([np.cos(th), np.sin(th), 0.0])
    p1, q1 = np.array([-1, 0, 0.]), np.array([1, 0, 0.])
    p2, q2 = np.array([0, 0, 0.5]) - dd, np.array([0, 0, 0.5]) + dd
    base = closest_pt_segseg(p1, q1, p2, q2)
    n0 = edge_normal(q1 - p1, q2 - p2)
    # perturb one master node component, measure d(s,t) and dn
    q2p = q2.copy(); q2p[1] += h
    pert = closest_pt_segseg(p1, q1, p2, q2p)
    np_ = edge_normal(q1 - p1, q2p - p2)
    dst = (abs(pert["s"] - base["s"]) + abs(pert["t"] - base["t"])) / h
    dn = np.linalg.norm(np_ - n0) / h
    return dst, dn, np.sin(th)

prev = None
for thd in (60.0, 30.0, 15.0):
    dst, dn, sinth = stn_sensitivity(thd)
    # dn/du should track ~1/sin(theta); param sens ~1/sin^2(theta)
    ok_finite = np.isfinite(dst) and np.isfinite(dn) and dst < 1e4 and dn < 1e3
    check(f"theta={thd:>4}deg: (s,t)&n sensitivity finite/bounded", ok_finite,
          f"d(s,t)/du={dst:.2f} (~1/sin^2={1/sinth**2:.2f})  dn/du={dn:.2f} (~1/sin={1/sinth:.2f})")
check("15deg gate: param sensitivity bounded (~1/sin^2(15)=14.9, not blown up)",
      stn_sensitivity(15.0)[0] < 100.0)

# --- T5: gN = +-||w|| ONLY in the margin-interior set (fails at a clamped near-vertex) ---
print("\nT5  gN == +-||w|| identity holds interior, FAILS clamped (the demote rule)")
# interior crossing
p1, q1 = np.array([-1, 0, 0.]), np.array([1, 0, 0.])
p2, q2 = np.array([0, -1, 1.]), np.array([0, 1, 1.])
ki = closest_pt_segseg(p1, q1, p2, q2)
gN, n, _ = edge_gap(p1, q1, p2, q2, ki["s"], ki["t"], sign_ref=np.array([0, 0, -1.]))
w = (p1 + ki["s"] * (q1 - p1)) - (p2 + ki["t"] * (q2 - p2))
check("interior: |gN| == ||w||", abs(abs(gN) - np.linalg.norm(w)) < TOL,
      f"|gN|={abs(gN):.7f} ||w||={np.linalg.norm(w):.7f}")
check("interior: pair IS margin-interior (fires edge-edge)", ki["interior"])
# clamped (edge-region): seg1 vertex vs seg2 interior
p2c, q2c = np.array([2, -1, 1.]), np.array([2, 1, 1.])
kc = closest_pt_segseg(p1, q1, p2c, q2c)
gNc, nc, _ = edge_gap(p1, q1, p2c, q2c, kc["s"], kc["t"], sign_ref=np.array([0, 0, -1.]))
wc = (p1 + kc["s"] * (q1 - p1)) - (p2c + kc["t"] * (q2c - p2c))
check("clamped: |gN| < ||w|| (cross-product n NOT the connector)",
      abs(gNc) < np.linalg.norm(wc) - 1e-6, f"|gN|={abs(gNc):.5f} ||w||={np.linalg.norm(wc):.5f}")
check("clamped: pair NOT margin-interior -> routing demotes to NTS/vertex", not kc["interior"])

# --- T6: FD dgN/d(coord) vs the analytic B-operator (envelope theorem -> exact) ---
print("\nT6  dgN/du: analytic B-operator vs central-difference (envelope theorem)")
p1, q1 = np.array([-1, 0.1, 0.]), np.array([1, -0.05, 0.])
p2, q2 = np.array([0.1, -1, 0.8]), np.array([-0.1, 1, 1.2])
k = closest_pt_segseg(p1, q1, p2, q2)
gN0, n0, sgn0 = edge_gap(p1, q1, p2, q2, k["s"], k["t"], sign_ref=np.array([0, 0, -1.]))
B = b_operator(n0, k["s"], k["t"])          # 4x3, rows = [p1,q1,p2,q2]
nodes = [p1, q1, p2, q2]
maxerr = 0.0
h = 1e-7
for ni in range(4):
    for d in range(3):
        pert = [x.copy() for x in nodes]
        pert[ni][d] += h
        kp = closest_pt_segseg(pert[0], pert[1], pert[2], pert[3])
        gp, _, _ = edge_gap(pert[0], pert[1], pert[2], pert[3], kp["s"], kp["t"], committed_sign=sgn0)
        pert[ni][d] -= 2 * h
        km = closest_pt_segseg(pert[0], pert[1], pert[2], pert[3])
        gm, _, _ = edge_gap(pert[0], pert[1], pert[2], pert[3], km["s"], km["t"], committed_sign=sgn0)
        fd = (gp - gm) / (2 * h)
        maxerr = max(maxerr, abs(fd - B[ni][d]))
check("max |B_analytic - FD| < 1e-6 (B is the exact first variation)", maxerr < 1e-6,
      f"max err={maxerr:.2e}")

# --- T7: BODY-FIXED committed sign keeps gN monotone through contact; w-rule FLIPS ---
print("\nT7  committed sign -> gN monotone through contact; the w.n rule masks penetration")
# seg1 slides along -z from above (z>1) through seg2 (at z=1) to below (z<1)
p2, q2 = np.array([0, -1, 1.]), np.array([0, 1, 1.])          # master edge along y at z=1
zs = np.linspace(1.6, 0.4, 13)                                 # drive slave down through master
# capture the committed sign ONCE at first contact engagement (gap just > 0), ref = master outward +z
p1a, q1a = np.array([-1, 0, zs[0]]), np.array([1, 0, zs[0]])
k0 = closest_pt_segseg(p1a, q1a, p2, q2)
_, _, sgn_committed = edge_gap(p1a, q1a, p2, q2, k0["s"], k0["t"], sign_ref=np.array([0, 0, 1.]))
gN_committed, gN_buggy = [], []
for z in zs:
    p1, q1 = np.array([-1, 0, z]), np.array([1, 0, z])
    k = closest_pt_segseg(p1, q1, p2, q2)
    gc, _, _ = edge_gap(p1, q1, p2, q2, k["s"], k["t"], committed_sign=sgn_committed)
    gb, _, _ = edge_gap(p1, q1, p2, q2, k["s"], k["t"])     # buggy w.n rule
    gN_committed.append(gc); gN_buggy.append(gb)
gN_committed, gN_buggy = np.array(gN_committed), np.array(gN_buggy)
mono = np.all(np.diff(gN_committed) < 1e-9)                    # strictly decreasing through contact
check("committed sign: gN strictly decreasing through contact (no flip)", mono,
      f"gN: {gN_committed[0]:+.2f} -> {gN_committed[-1]:+.2f}")
check("committed sign: gN actually goes NEGATIVE (penetration registered)", gN_committed[-1] < 0)
check("buggy w.n rule: gN stays >= 0 (penetration MASKED — the A-4 bug)",
      np.all(gN_buggy >= -1e-9), f"buggy min={gN_buggy.min():+.3f}")

print(f"\n{'='*64}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(ADR-57 E0 edge-edge closest-point oracle)\n{'='*64}")
sys.exit(1 if _fails else 0)
