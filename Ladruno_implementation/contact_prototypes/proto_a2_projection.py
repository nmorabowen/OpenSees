"""ADR-41 A2 ORACLE — LadrunoContactProjection.projectFull (mortar GP geometry).

Oracle-first (the fork discipline): validate the RICH projection the mortar
Gauss-point loop needs — covariant tangents t1=g1, t2=g2, the surface metric
g_ab = g_a·g_b, and the master shape fns phi — in pure numpy BEFORE trusting the
transcription in SRC/domain/contact/LadrunoContactProjection.h. No OpenSees, no build.

The base projection (bounded Newton, oriented winding-immune normal, gap) is already
validated by proto_p2b_nts.py (7/7) and is moved here bit-for-bit; this oracle covers
the A2 ADDITIONS:
  - covariant tangents g1=∂x̄/∂ξ, g2=∂x̄/∂η  (NOT unit)
  - surface metric g[a][b] = g_a·g_b   (the slip return uses R=g·r, NOT I₂)
  - master shape fns phi at the projection  (partition of unity Σφ=1)
  - projectFull returns geometry for ANY converged projection (gap as DATA, no
    penetration/in-bounds gate, unlike evalSegment)
  - bounded-Newton sentinel: a degenerate/collinear segment returns status -1 and
    ZEROED geometry (never garbage)
  - projectFull vs evalSegment CONSISTENCY (same gap/n on a penetrating in-bounds point)

Run: <pythoncore-3.12> proto_a2_projection.py   (numpy only)
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


# ---------------------------------------------------------------- geom (mirror C++)
def shape(nps, xi, eta):
    if nps == 4:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dNe = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    else:  # tri-3, N = [1-xi-eta, xi, eta]
        N = np.array([1 - xi - eta, xi, eta])
        dNx = np.array([-1.0, 1.0, 0.0])
        dNe = np.array([-1.0, 0.0, 1.0])
    return N, dNx, dNe


def project(nps, X, xs, tolR=1e-12, maxit=10):
    """returns (xi, eta, status): 0 conv in-bounds, 1 conv oob, -1 no valid proj."""
    xi, eta = (0.0, 0.0) if nps == 4 else (1 / 3, 1 / 3)
    conv = False
    for _ in range(maxit):
        N, dNx, dNe = shape(nps, xi, eta)
        xbar = N @ X
        g1, g2 = dNx @ X, dNe @ X
        d = xs - xbar
        R0, R1 = d @ g1, d @ g2
        if np.hypot(R0, R1) < tolR:
            conv = True
            break
        K = np.array([[g1 @ g1, g1 @ g2], [g2 @ g1, g2 @ g2]])
        detK = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]
        scale = np.linalg.norm(g1) * np.linalg.norm(g2)
        if abs(detK) < 1e-14 * (scale + 1e-300):
            return xi, eta, -1                       # degenerate
        dxi = (K[1, 1] * R0 - K[0, 1] * R1) / detK
        deta = (-K[1, 0] * R0 + K[0, 0] * R1) / detK
        xi += dxi
        eta += deta
    if not conv:
        return xi, eta, -1
    t = 1e-9
    if nps == 4:
        inb = (-1 - t <= xi <= 1 + t) and (-1 - t <= eta <= 1 + t)
    else:
        inb = (xi >= -t) and (eta >= -t) and (xi + eta <= 1 + t)
    return xi, eta, (0 if inb else 1)


def normal_oriented(nps, xi, eta, X, refDir):
    N, dNx, dNe = shape(nps, xi, eta)
    g1, g2 = dNx @ X, dNe @ X
    raw = np.cross(g1, g2)
    j = np.linalg.norm(raw)
    if j < 1e-300:
        return None
    n = raw / j
    p = n @ refDir
    if abs(p) < 1e-12 * (np.linalg.norm(refDir) + 1e-300):
        return None
    return -n if p < 0 else n


def project_full(nps, X, xs, refDir, tolR=1e-12, maxit=10):
    """mirror of LadrunoContactProjection::projectFull."""
    p = {"status": -1, "xi": 0.0, "eta": 0.0, "gap": 0.0,
         "n": np.zeros(3), "t1": np.zeros(3), "t2": np.zeros(3),
         "g": np.zeros((2, 2)), "phi": np.zeros(4)}
    xi, eta, st = project(nps, X, xs, tolR, maxit)
    p["status"], p["xi"], p["eta"] = st, xi, eta
    if st < 0:
        return p
    N, dNx, dNe = shape(nps, xi, eta)
    xbar = N @ X
    g1, g2 = dNx @ X, dNe @ X
    n = normal_oriented(nps, xi, eta, X, refDir)
    if n is None:
        p["status"] = -1
        return p
    p["n"], p["t1"], p["t2"] = n, g1, g2
    p["g"] = np.array([[g1 @ g1, g1 @ g2], [g2 @ g1, g2 @ g2]])
    p["gap"] = n @ (xs - xbar)
    p["phi"][:len(N)] = N
    return p


def eval_segment(nps, X, xs, refDir, tolR=1e-12, maxit=10):
    """mirror of evalSegment (penetrating-only); returns (active, gap, n, N)."""
    xi, eta, st = project(nps, X, xs, tolR, maxit)
    if st != 0:
        return False, 0.0, None, None
    n = normal_oriented(nps, xi, eta, X, refDir)
    if n is None:
        return False, 0.0, None, None
    N, _, _ = shape(nps, xi, eta)
    xbar = N @ X
    gap = n @ (xs - xbar)
    return (gap < 0.0), gap, n, N


def rot_matrix(axis, ang):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    c, s = np.cos(ang), np.sin(ang)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) * c + s * K + (1 - c) * np.outer(a, a)


# ======================================================================= tests
print("ADR-41 A2 oracle — LadrunoContactProjection.projectFull")

# --- T1: flat axis-aligned quad (side 2) — closed form metric=I, n=ẑ, phi PoU ---
print("\nT1  flat unit-mapped quad (X≡parametric): metric=I, n=ẑ, gap=z, Σφ=1")
Xq = np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]], float)
xs = np.array([0.3, -0.2, 0.5])
p = project_full(4, Xq, xs, refDir=np.array([0, 0, 1.0]))
check("status in-bounds", p["status"] == 0, f"status={p['status']}")
check("xi,eta == (0.3,-0.2)", np.allclose([p["xi"], p["eta"]], [0.3, -0.2], atol=TOL))
check("metric g == I", np.allclose(p["g"], np.eye(2), atol=TOL), f"g={p['g'].tolist()}")
check("t1==x̂, t2==ŷ", np.allclose(p["t1"], [1, 0, 0]) and np.allclose(p["t2"], [0, 1, 0]))
check("n == ẑ", np.allclose(p["n"], [0, 0, 1], atol=TOL))
check("gap == +0.5 (separated, data not gated)", abs(p["gap"] - 0.5) < TOL)
check("Σφ == 1 (partition of unity)", abs(p["phi"][:4].sum() - 1.0) < TOL)

# --- T2: SAME quad rigidly tilted — metric is rotation-invariant (still I) ---
print("\nT2  tilted quad: metric rotation-invariant (==I); n rotates; gap preserved")
R = rot_matrix([0.4, 1.0, 0.3], 0.7)
Xt = Xq @ R.T
xst = xs @ R.T
pt = project_full(4, Xt, xst, refDir=R @ np.array([0, 0, 1.0]))
check("status in-bounds", pt["status"] == 0)
check("metric g still == I (rotation-invariant)", np.allclose(pt["g"], np.eye(2), atol=1e-8))
check("n == R·ẑ", np.allclose(pt["n"], R @ np.array([0, 0, 1.0]), atol=1e-8))
check("gap preserved == 0.5", abs(pt["gap"] - 0.5) < 1e-8)

# --- T3: anisotropic scaling — metric reads off the stretch (g00=a², g11=b²) ---
print("\nT3  scaled quad (x×a, y×b): metric diag == (a², b²), off-diag 0")
a, b = 2.5, 0.6
Xs = Xq * np.array([a, b, 1.0])
ps = project_full(4, Xs, np.array([0.0, 0.0, 0.4]), refDir=np.array([0, 0, 1.0]))
check("g00 == a²", abs(ps["g"][0, 0] - a * a) < TOL, f"g00={ps['g'][0,0]:.6f} a²={a*a:.6f}")
check("g11 == b²", abs(ps["g"][1, 1] - b * b) < TOL)
check("g01 == 0 (orthogonal tangents)", abs(ps["g"][0, 1]) < TOL)

# --- T4: WARPED (non-planar) quad — project converges, metric SPD, n unit, Σφ=1 ---
# Slave sits NEAR the facet (small gap): the regime the mortar GP loop uses and the
# bounded Gauss-Newton supports (it drops the O(gap) curvature term by design, so a
# far point over a curved facet would not reach 1e-12 in maxit — that is accepted).
print("\nT4  warped non-planar quad (near point): converges, metric SPD, |n|=1, Σφ=1")
Xw = np.array([[-1, -1, 0.0], [1, -1, 0.1], [1, 1, -0.05], [-1, 1, 0.15]])
pw = project_full(4, Xw, np.array([0.1, 0.2, 0.2]), refDir=np.array([0, 0, 1.0]))
detg = pw["g"][0, 0] * pw["g"][1, 1] - pw["g"][0, 1] ** 2
check("status valid (0/1)", pw["status"] in (0, 1), f"status={pw['status']}")
check("metric SPD (det>0, g00>0)", detg > TOL and pw["g"][0, 0] > TOL, f"det(g)={detg:.6f}")
check("|n| == 1", abs(np.linalg.norm(pw["n"]) - 1.0) < TOL)
check("metric symmetric", abs(pw["g"][0, 1] - pw["g"][1, 0]) < TOL)
check("Σφ == 1", abs(pw["phi"][:4].sum() - 1.0) < TOL)

# --- T5: tri-3 — metric of a known right triangle, φ[3]==0 ---
print("\nT5  tri-3 (legs Lx,Ly): metric==[[Lx²,0],[0,Ly²]], φ[3]==0")
Lx, Ly = 3.0, 2.0
Xtri = np.array([[0, 0, 0], [Lx, 0, 0], [0, Ly, 0]], float)
ptri = project_full(3, Xtri, np.array([0.5, 0.5, 0.7]), refDir=np.array([0, 0, 1.0]))
check("g == diag(Lx²,Ly²)",
      np.allclose(ptri["g"], np.array([[Lx * Lx, 0], [0, Ly * Ly]]), atol=TOL),
      f"g={ptri['g'].tolist()}")
check("tri-3 φ[3] == 0", abs(ptri["phi"][3]) < TOL)
check("Σφ == 1 (first 3)", abs(ptri["phi"][:3].sum() - 1.0) < TOL)

# --- T6: bounded-Newton SENTINEL — degenerate (collinear) segment → status -1, zeroed ---
print("\nT6  degenerate/collinear segment: status==-1, geometry ZEROED (no garbage)")
Xdeg = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)  # all collinear
pdeg = project_full(4, Xdeg, np.array([1.5, 0.0, 1.0]), refDir=np.array([0, 0, 1.0]))
check("status == -1 (sentinel)", pdeg["status"] == -1, f"status={pdeg['status']}")
check("n zeroed", np.allclose(pdeg["n"], 0.0))
check("metric zeroed", np.allclose(pdeg["g"], 0.0))
check("gap zeroed", pdeg["gap"] == 0.0)

# --- T7: projectFull vs evalSegment CONSISTENCY on a penetrating in-bounds point ---
print("\nT7  projectFull vs evalSegment: identical gap & n on a penetrating point")
xpen = np.array([0.2, -0.1, -0.3])  # below the z=0 quad ⇒ penetrating (gap<0)
pf = project_full(4, Xq, xpen, refDir=np.array([0, 0, 1.0]))
active, gE, nE, NE = eval_segment(4, Xq, xpen, refDir=np.array([0, 0, 1.0]))
check("evalSegment active (penetrating)", active and gE < 0)
check("projectFull gap == evalSegment gap", abs(pf["gap"] - gE) < TOL, f"{pf['gap']:.9f} vs {gE:.9f}")
check("projectFull n == evalSegment n", np.allclose(pf["n"], nE, atol=TOL))
check("projectFull φ == evalSegment N", np.allclose(pf["phi"][:4], NE, atol=TOL))

print(f"\n{'='*60}\n{'ALL PASS' if _fails == 0 else str(_fails) + ' FAILURE(S)'}  "
      f"(A2 projection oracle)\n{'='*60}")
sys.exit(1 if _fails else 0)
