"""ADR-85 T1a ORACLE -- LadrunoContact2DKernel (2D NTS geometry + the vertex policy).

Oracle-first (the fork discipline): an INDEPENDENT reference implementation of the 2D
node-to-segment geometry, written from the ADR-85 SS How/1 formulas, validated against
closed-form geometry and against EXACT (60-digit) arithmetic, and only then used as the
bit-parity gate for the transcription in SRC/domain/contact/LadrunoContact2DKernel.h.

What is covered (ADR-85 gate G-T1a):
  projectSegment2D  -- closed-form point-on-segment projection, oriented normal, gap,
                       the four-way refusal status (IN_BOUNDS / OUT_LOW / OUT_HIGH /
                       DEGENERATE), out-of-bounds REFUSAL (never a clamp).
  vertexWedge2D     -- the corner classifier (wedge membership, convex/concave, side sign).
  vertexEval2D      -- the vertex pair (radial normal, signed gap, coincidence refusal).
  sigmaFromRef2D    -- orientation from a reference datum + the tangent-datum refusal.
  bOperatorSegment2D / bOperatorVertex2D -- FD-gated as the EXACT first variation of gap.

Falsifier families (every one exits nonzero on failure):
  F1 seg in-bounds parity      F2 refusal/status parity (EXACT, incl. the tolIn boundary)
  F3 far-from-origin           F4 zero-length + bad-Lref refusal
  F5 sigma flip / winding      F6 vertex parity + coincidence refusal
  F7 wedge parity              F8 B-operator parity + FD
  F9 sigmaFromRef parity       CA/CB/CC/CD/CE the deterministic corner suites
  XA exact-arithmetic gap reference (the absolute noise floor, measured not assumed)

------------------------------------------------------------------------------
RESULT 1 -- THE CONCAVE-VERTEX DECISION (settles the ADR-85 open question
"vertex-policy sign robustness", and it INVERTS the ADR's premise).

  Measured (family CD, 40000 random corners + the deterministic sweeps):
    * CONVEX corner: 20000/20000 random wedge points have gap = +||xs-Xv|| > 0 and lie
      OUTSIDE the master material. A convex vertex pair is a valid C0 extension of its
      two neighbours but is NEVER penetrating, so it can NEVER carry force.
    * CONCAVE corner: 20000/20000 random wedge points have gap = -||xs-Xv|| < 0 and lie
      INSIDE the master material, while BOTH adjacent segments refuse them. That is a
      genuine finite-measure coverage hole at finite penetration.
    * The concave interior sweep at penetration depth R = 0.25 has a 0.3927-long wedge
      arc (19.7% of the swept path, = R*|turn|) that policy "concave never activates"
      leaves owned by nobody -- penetrating by R with ZERO restoring force.

  => DECISION: concave vertices ACTIVATE and need the vertex rule. They are the ONLY
     vertex pairs that can ever carry force. Convex vertex pairs are retained for the
     ownership map and normal continuity, but a handler that instantiates vertex pairs
     only at CONCAVE corners loses no force. (The interior of a convex corner is covered
     TWICE by the two adjacent strips -- the pre-existing NTS corner double-count, not a
     hole, and out of T1a scope.)

RESULT 2 -- THE MEASURED HANDOFF BAND, and a correction to the naive reading of
ADR-85 SS How/1 (family CC).

  "Segment in-bounds" and "in the vertex wedge" are ALGEBRAIC complements at the same
  parametric slack tolIn = 1e-9 -- but they are NOT floating-point complements, so the
  obvious "exactly one owner = segment XOR wedge" reading of the ADR is wrong:

    * projectSegment2D forms xi_prev from X_prev, so at a corner xi_prev ~ 1 and the
      test xi > 1+tolIn can only resolve the overshoot to ulp(1) = 2.220e-16, whereas
      aPrev resolves it exactly. MEASURED disagreement between the two expressions:
          max |(xi_prev - 1) - aPrev| = 1.149e-15  (5.2 ulp(1) = 1.15e-06 of tolIn)
          max |xi_next  -  aNext|     = 0.0        (same expression -> no band)
    * A 16400-station scan across the two seams found the naive XOR giving 80
      DOUBLE-owner and 1 ZERO-owner stations. No tolerance value closes that -- the two
      quantities have different resolutions.

  CONTRACT FOR T1b (implemented and gated here, and written into the header):
      1. previous segment in-bounds  -> the previous segment owns
      2. else next segment in-bounds -> the next segment owns
      3. else aPrev > 0 AND aNext < 0 (UNSLACKED) -> the vertex owns
      4. else -> nobody
  Ordered like ADR-57 E7. Total and unique over all 16400 seam stations, with no
  tolerance coupling. Inside the overlap the two candidates are physically the same
  contact: relative gap difference <= 1.7e-14 and normal disagreement <= 1.06e-07 rad,
  the latter matching its analytic bound atan(tolIn*L/|gap|) to 1.00006.

RESULT 3 -- A CONDITIONING DEFECT THE ORACLE CAUGHT IN THE FIRST DRAFT (family F7).
  The vertex side sign comes from the bisector n_prev + n_next. The first draft gated it
  with |r.bis| > tauPerp*|r|*|bis| -- "relative to everything". But bis is a CANCELLING
  sum of two UNIT vectors: its components carry an ~eps absolute error however small
  |bis| becomes, so scaling the threshold by |bis| sinks it BELOW the noise floor as the
  corner folds back. Measured over 48902 wedge hits (25% near-fold-back): 707 WRONG side
  signs. Gating on |r.bis| > tauPerp*|r| instead gives 0 wrong signs (2412 deferrals).
  Fixed in the kernel; the F7 family is the regression witness.

RESULT 4 -- parity. C++ (g++ -O2 -ffp-contract=off) vs this reference over 190600
  randomized records: every discrete field (status, inWedge, corner, sign) agreed
  EXACTLY, and every continuous field agreed to 0.0 -- bit-for-bit, no ulp of drift,
  including the far-field family. Identical operation order on both sides is deliberate:
  it turns "within tolerance" into "bit-for-bit". The RELATIVE bounds are still asserted
  so a future divergence is caught rather than absorbed.
------------------------------------------------------------------------------

Run: python3.12 proto_t1_nts2d.py        (numpy + a g++ on PATH)
     python3.12 proto_t1_nts2d.py --nocpp    (numpy families only, no compiler needed)
"""
import math
import os
import subprocess
import sys
import tempfile
from decimal import Decimal, getcontext

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

getcontext().prec = 60
SEED = 20260818

# --- kernel status codes (mirror LadrunoContact2DKernel.h) -------------------
SEG_IN, SEG_LOW, SEG_HIGH, SEG_DEG = 0, 1, 2, -1
VTX_OK, VTX_DEG = 0, -1
C_CONVEX, C_STRAIGHT, C_CONCAVE = -1, 0, 1

TOL_IN = 1e-9      # parametric slack on xi (dimensionless)
TAU_SEG = 1e-8     # zero-length refusal, RELATIVE to Lref
TAU_PERP = 1e-12   # ambiguity gate

EPS = np.finfo(np.float64).eps

# --- reporting --------------------------------------------------------------
_fails = 0
_fam = {}          # family -> [checks, fails, notes]


def check(fam, name, ok, extra=""):
    global _fails
    st = _fam.setdefault(fam, [0, 0])
    st[0] += 1
    if not ok:
        st[1] += 1
        _fails += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {fam}: {name}{(' -- ' + extra) if extra else ''}")


# ============================================================================
# REFERENCE IMPLEMENTATION -- written from ADR-85 SS How/1.
#
# The arithmetic is done in explicit scalar double operations (not np.dot) on purpose:
# numpy is free to reassociate/vectorize a reduction, and a reassociated dot product
# would turn an exact parity statement into a "within 1 ulp" statement. Keeping the
# operation order identical on both sides makes the C++ comparison a BIT-LEVEL gate.
# numpy is used for the RNG and the bookkeeping.
# ============================================================================
def shape2d(xi):
    return (1.0 - xi, xi)


def project_segment_2d(X0, X1, xs, sigma, Lref, tolIn=TOL_IN, tauSeg=TAU_SEG):
    """(status, xi, gap, (nx,ny)) -- ADR-85 SS How/1.

      t  = X1 - X0,  xi = (xs-X0).t / L^2,  n = sigma*perp(t)/L,  g = (xs-x(xi)).n
    Refusals: L < tauSeg*Lref (or Lref<=0) -> DEGENERATE with zeroed outputs;
              xi outside [-tolIn, 1+tolIn] -> OUT_LOW / OUT_HIGH (NEVER clamped).
    """
    tx = X1[0] - X0[0]
    ty = X1[1] - X0[1]
    L2 = tx * tx + ty * ty
    L = math.sqrt(L2)
    if not (Lref > 0.0):
        return (SEG_DEG, 0.0, 0.0, (0.0, 0.0))
    if L < tauSeg * Lref:
        return (SEG_DEG, 0.0, 0.0, (0.0, 0.0))

    d0x = xs[0] - X0[0]
    d0y = xs[1] - X0[1]
    xi = (d0x * tx + d0y * ty) / L2

    nx = sigma * (-ty) / L
    ny = sigma * (tx) / L

    N0, N1 = shape2d(xi)
    bx = N0 * X0[0] + N1 * X1[0]
    by = N0 * X0[1] + N1 * X1[1]
    gap = (xs[0] - bx) * nx + (xs[1] - by) * ny

    if xi < -tolIn:
        st = SEG_LOW
    elif xi > 1.0 + tolIn:
        st = SEG_HIGH
    else:
        st = SEG_IN
    return (st, xi, gap, (nx, ny))


def vertex_wedge_2d(tPrev, tNext, sigma, xs, Xv, Lref,
                    tolIn=TOL_IN, tauSeg=TAU_SEG, tauPerp=TAU_PERP):
    """dict(inWedge, corner, sideSign, turnSin, aPrev, aNext) -- ADR-85 SS How/1.

    Wedge  <=>  xi_prev > 1+tolIn AND xi_next < -tolIn, written vertex-relative as
    aPrev > tolIn and aNext < -tolIn (the "toward the shared vertex" qualifier is what
    separates a corner from the far end of an open patch).
    Corner type from the signed turn sin; SIDE SIGN from the bisector of the two outward
    normals, which stays well conditioned where sign(turnSin) coin-flips.
    """
    out = dict(inWedge=0, corner=C_STRAIGHT, sideSign=0.0,
               turnSin=0.0, aPrev=0.0, aNext=0.0)
    LP2 = tPrev[0] * tPrev[0] + tPrev[1] * tPrev[1]
    LN2 = tNext[0] * tNext[0] + tNext[1] * tNext[1]
    LP = math.sqrt(LP2)
    LN = math.sqrt(LN2)
    if not (Lref > 0.0):
        return out
    if LP < tauSeg * Lref or LN < tauSeg * Lref:
        return out

    rx = xs[0] - Xv[0]
    ry = xs[1] - Xv[1]
    out["aPrev"] = (rx * tPrev[0] + ry * tPrev[1]) / LP2
    out["aNext"] = (rx * tNext[0] + ry * tNext[1]) / LN2
    out["inWedge"] = 1 if (out["aPrev"] > tolIn and out["aNext"] < -tolIn) else 0

    out["turnSin"] = sigma * (tPrev[0] * tNext[1] - tPrev[1] * tNext[0]) / (LP * LN)
    if out["turnSin"] < 0.0:
        out["corner"] = C_CONVEX
    elif out["turnSin"] > 0.0:
        out["corner"] = C_CONCAVE
    else:
        out["corner"] = C_STRAIGHT

    nPx = sigma * (-tPrev[1]) / LP
    nPy = sigma * (tPrev[0]) / LP
    nNx = sigma * (-tNext[1]) / LN
    nNy = sigma * (tNext[0]) / LN
    bx = nPx + nNx
    by = nPy + nNy
    bl = math.sqrt(bx * bx + by * by)
    rl = math.sqrt(rx * rx + ry * ry)
    # conditioning gate scaled by |r| ONLY: bis is a cancelling sum of two UNIT vectors,
    # so its components carry an ~eps absolute error however small |bis| is. A threshold
    # that ALSO scales with |bis| sinks under that noise floor as the corner folds back
    # and admits a coin-flip sign (measured: 707 wrong signs in 48902 wedge hits).
    if bl > 2.0 * tauPerp and rl > 0.0:
        p = rx * bx + ry * by
        if math.fabs(p) > tauPerp * rl:
            out["sideSign"] = 1.0 if p > 0.0 else -1.0
    return out


def vertex_eval_2d(Xv, xs, sigmaSide, tolLen):
    """(status, gap, (nx,ny)) -- n = sigmaSide*(xs-Xv)/d, gap = sigmaSide*d."""
    rx = xs[0] - Xv[0]
    ry = xs[1] - Xv[1]
    d = math.sqrt(rx * rx + ry * ry)
    if not (tolLen > 0.0):
        return (VTX_DEG, 0.0, (0.0, 0.0))
    if d < tolLen:
        return (VTX_DEG, 0.0, (0.0, 0.0))
    return (VTX_OK, sigmaSide * d, (sigmaSide * rx / d, sigmaSide * ry / d))


def sigma_from_ref_2d(t, refDir, tauPerp=TAU_PERP):
    L = math.sqrt(t[0] * t[0] + t[1] * t[1])
    R = math.sqrt(refDir[0] * refDir[0] + refDir[1] * refDir[1])
    if not (L > 0.0) or not (R > 0.0):
        return 0.0
    proj = (-t[1] * refDir[0] + t[0] * refDir[1]) / L
    if math.fabs(proj) < tauPerp * R:
        return 0.0
    return 1.0 if proj > 0.0 else -1.0


def b_operator_segment_2d(N0, N1, n):
    return [n[0], n[1], -N0 * n[0], -N0 * n[1], -N1 * n[0], -N1 * n[1]]


def b_operator_vertex_2d(n):
    return [n[0], n[1], -n[0], -n[1]]


# ---------------------------------------------------------------- geometry helpers
def unit(a):
    L = math.hypot(a[0], a[1])
    return (a[0] / L, a[1] / L)


def outward_normal(t, sigma):
    L = math.hypot(t[0], t[1])
    return (sigma * (-t[1]) / L, sigma * t[0] / L)


def point_from_params(tPrev, tNext, Xv, alpha, beta):
    """Slave point with prescribed (aPrev, aNext) = (alpha, beta). Solves the 2x2
    system r.tPrev = alpha*LP^2, r.tNext = beta*LN^2 -- lets the seam be probed at any
    parametric offset, which is how the handoff band is measured."""
    LP2 = tPrev[0] ** 2 + tPrev[1] ** 2
    LN2 = tNext[0] ** 2 + tNext[1] ** 2
    M = np.array([[tPrev[0], tPrev[1]], [tNext[0], tNext[1]]], dtype=float)
    b = np.array([alpha * LP2, beta * LN2], dtype=float)
    r = np.linalg.solve(M, b)
    return (Xv[0] + float(r[0]), Xv[1] + float(r[1]))


def gap_exact(X0, X1, xs, sigma):
    """Exact (60-digit) signed distance to the segment's line, via the algebraically
    equivalent cross-product form  g = sigma*cross(t, xs-X0)/L.  An INDEPENDENT
    derivation (no shape functions), used to measure the double-precision noise floor."""
    tx = Decimal(X1[0]) - Decimal(X0[0])
    ty = Decimal(X1[1]) - Decimal(X0[1])
    dx = Decimal(xs[0]) - Decimal(X0[0])
    dy = Decimal(xs[1]) - Decimal(X0[1])
    L = (tx * tx + ty * ty).sqrt()
    return Decimal(sigma) * (tx * dy - ty * dx) / L


def owner_xor(Xprev, Xv, Xnext, sigma, Lref, xs):
    """The NAIVE reading of ADR-85 SS How/1: "segment in-bounds" XOR "in the slacked
    wedge". Kept ONLY to measure how far it is from being floating-point exact -- see
    family CC. Returns the number of claimants."""
    tP = (Xv[0] - Xprev[0], Xv[1] - Xprev[1])
    tN = (Xnext[0] - Xv[0], Xnext[1] - Xv[1])
    stP = project_segment_2d(Xprev, Xv, xs, sigma, Lref)[0]
    stN = project_segment_2d(Xv, Xnext, xs, sigma, Lref)[0]
    w = vertex_wedge_2d(tP, tN, sigma, xs, Xv, Lref)
    return (1 if stP == SEG_IN else 0) + (1 if stN == SEG_IN else 0) + w["inWedge"]


def owner_of(Xprev, Xv, Xnext, sigma, Lref, xs, committed_side=0.0):
    """The ORDERED ownership rule documented in LadrunoContact2DKernel.h: segment
    precedence, then the UNSLACKED vertex claim. Total and unique with no tolerance
    coupling. Returns (n_candidates, which, gap, n) with which in
    {0 prev-seg, 1 next-seg, 2 vertex, -1 nobody}."""
    tP = (Xv[0] - Xprev[0], Xv[1] - Xprev[1])
    tN = (Xnext[0] - Xv[0], Xnext[1] - Xv[1])
    stP, _, gP, nP = project_segment_2d(Xprev, Xv, xs, sigma, Lref)
    stN, _, gN, nN = project_segment_2d(Xv, Xnext, xs, sigma, Lref)
    w = vertex_wedge_2d(tP, tN, sigma, xs, Xv, Lref)
    vtx_claim = w["aPrev"] > 0.0 and w["aNext"] < 0.0        # UNSLACKED
    ncand = (1 if stP == SEG_IN else 0) + (1 if stN == SEG_IN else 0) + (1 if vtx_claim else 0)
    if stP == SEG_IN:
        return (ncand, 0, gP, nP)
    if stN == SEG_IN:
        return (ncand, 1, gN, nN)
    if vtx_claim:
        side = committed_side if committed_side != 0.0 else w["sideSign"]
        st, g, n = vertex_eval_2d(Xv, xs, side, TAU_SEG * Lref)
        return (ncand, 2 if st == VTX_OK else -1, g, n)
    return (ncand, -1, 0.0, (0.0, 0.0))


def offset_sweep(Xprev, Xv, Xnext, sigma, off, frac, nsta):
    """Stations along the OFFSET (parallel) curve of the corner at signed distance off:
    parallel to the previous segment, an arc of radius |off| about the vertex, parallel
    to the next segment. off > 0 = exterior, off < 0 = into the material. Constant
    physical distance to the polyline by construction, so any gap variation along it is
    a defect. Returns (list of stations, arclength step, arc length of the wedge phase)."""
    tP = (Xv[0] - Xprev[0], Xv[1] - Xprev[1])
    tN = (Xnext[0] - Xv[0], Xnext[1] - Xv[1])
    lp, ln = math.hypot(*tP), math.hypot(*tN)
    nP, nN = outward_normal(tP, sigma), outward_normal(tN, sigma)
    sgn = 1.0 if off >= 0.0 else -1.0
    a0 = math.atan2(sgn * nP[1], sgn * nP[0])
    a1 = math.atan2(sgn * nN[1], sgn * nN[0])
    d = a1 - a0
    while d > math.pi:
        d -= 2 * math.pi
    while d < -math.pi:
        d += 2 * math.pi
    R = abs(off)
    L1, Larc, L2 = frac * lp, R * abs(d), frac * ln
    Ltot = L1 + Larc + L2
    ds = Ltot / (nsta - 1)
    pts = []
    for i in range(nsta):
        s = i * ds
        if s <= L1:
            back = L1 - s
            pts.append((Xv[0] - back * tP[0] / lp + off * nP[0],
                        Xv[1] - back * tP[1] / lp + off * nP[1]))
        elif s <= L1 + Larc:
            u = (s - L1) / Larc if Larc > 0 else 0.0
            a = a0 + u * d
            pts.append((Xv[0] + R * math.cos(a), Xv[1] + R * math.sin(a)))
        else:
            fwd = s - L1 - Larc
            pts.append((Xv[0] + fwd * tN[0] / ln + off * nN[0],
                        Xv[1] + fwd * tN[1] / ln + off * nN[1]))
    return pts, ds, Larc, Ltot


def ang_between(a, b):
    """atan2(|cross|, dot), NOT acos(dot): acos is ill-conditioned near parallel vectors
    (its answer has a sqrt(eps) ~ 1.5e-8 rad noise floor), which would swamp exactly the
    small-angle measurements this oracle makes at the segment/vertex handoff."""
    return math.atan2(abs(a[0] * b[1] - a[1] * b[0]), a[0] * b[0] + a[1] * b[1])


# ============================================================================
print("=" * 78)
print("ADR-85 T1a oracle -- LadrunoContact2DKernel (2D NTS geometry + vertex policy)")
print("=" * 78)

rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------- A. closed form
print("\nA  closed-form geometry (no C++, no parity -- the independent anchor)")
st, xi, gap, n = project_segment_2d((0.0, 0.0), (2.0, 0.0), (0.5, -0.3), 1.0, 2.0)
check("A", "flat segment: status IN_BOUNDS", st == SEG_IN)
check("A", "flat segment: xi == 0.25", abs(xi - 0.25) < 1e-15)
check("A", "flat segment: n == (0,1)", abs(n[0]) < 1e-15 and abs(n[1] - 1) < 1e-15)
check("A", "flat segment: gap == -0.3 (penetration <=> gap<0)", abs(gap + 0.3) < 1e-15)

# a 3-4-5 triangle: the perpendicular distance from (0,0) to the line through
# (4,0)-(0,3) is 12/5 = 2.4 exactly, with the foot at xi = 16/25.
st, xi, gap, n = project_segment_2d((4.0, 0.0), (0.0, 3.0), (0.0, 0.0), 1.0, 5.0)
check("A", "3-4-5 line: xi == 16/25", abs(xi - 16.0 / 25.0) < 1e-15, f"xi={xi!r}")
check("A", "3-4-5 line: |gap| == 12/5 exactly", abs(abs(gap) - 2.4) < 1e-15, f"gap={gap!r}")

# gap must be the SIGNED DISTANCE, so a rigid rotation of the whole configuration
# leaves it invariant while n co-rotates.
th = 0.7
c, s = math.cos(th), math.sin(th)


def rot(p):
    return (c * p[0] - s * p[1], s * p[0] + c * p[1])


X0, X1, xs = (0.3, -1.2), (2.7, 0.9), (1.1, 0.4)
st0, xi0, g0, n0 = project_segment_2d(X0, X1, xs, 1.0, 1.0)
stR, xiR, gR, nR = project_segment_2d(rot(X0), rot(X1), rot(xs), 1.0, 1.0)
check("A", "rigid rotation: xi invariant", abs(xiR - xi0) < 1e-14)
check("A", "rigid rotation: gap invariant", abs(gR - g0) < 1e-14, f"d={abs(gR-g0):.2e}")
check("A", "rigid rotation: n co-rotates", abs(nR[0] - rot(n0)[0]) < 1e-14
      and abs(nR[1] - rot(n0)[1]) < 1e-14)

# uniform scaling: xi invariant (dimensionless), gap scales linearly
k = 1e3
stS, xiS, gS, nS = project_segment_2d((k * X0[0], k * X0[1]), (k * X1[0], k * X1[1]),
                                      (k * xs[0], k * xs[1]), 1.0, k * 1.0)
check("A", "uniform scaling: xi invariant (dimensionless)", abs(xiS - xi0) < 1e-13)
check("A", "uniform scaling: gap scales by k", abs(gS - k * g0) < 1e-11 * k)

check("A", "shape2D partition of unity", abs(sum(shape2d(0.37)) - 1.0) < 1e-16)
check("A", "B-operator rows sum to zero (self-equilibrated force)",
      all(abs(sum(b_operator_segment_2d(*shape2d(0.37), n=(0.6, 0.8))[i::2])) < 1e-15
          for i in (0, 1)))
check("A", "vertex B rows sum to zero", all(
    abs(sum(b_operator_vertex_2d((0.6, 0.8))[i::2])) < 1e-16 for i in (0, 1)))

# --------------------------------------------------- B. B == d(gap)/du by FD (exact)
print("\nB  B-operator == the EXACT first variation of gap (finite differences)")
X0, X1, xs = (0.2, -0.4), (2.3, 0.7), (1.0, 0.05)
_, xi, g, n = project_segment_2d(X0, X1, xs, 1.0, 1.0)
N0, N1 = shape2d(xi)
B = b_operator_segment_2d(N0, N1, n)
q = [xs[0], xs[1], X0[0], X0[1], X1[0], X1[1]]
h, maxe = 1e-7, 0.0
for j in range(6):
    qp, qm = list(q), list(q)
    qp[j] += h
    qm[j] -= h
    gp = project_segment_2d((qp[2], qp[3]), (qp[4], qp[5]), (qp[0], qp[1]), 1.0, 1.0)[2]
    gm = project_segment_2d((qm[2], qm[3]), (qm[4], qm[5]), (qm[0], qm[1]), 1.0, 1.0)[2]
    maxe = max(maxe, abs((gp - gm) / (2 * h) - B[j]))
check("B", "segment B == FD d(gap)/du", maxe < 1e-7, f"max err={maxe:.2e}")

Xv, xp = (0.4, 0.9), (1.3, -0.2)
_, gv, nv = vertex_eval_2d(Xv, xp, -1.0, 1e-9)
BV = b_operator_vertex_2d(nv)
qv = [xp[0], xp[1], Xv[0], Xv[1]]
maxev = 0.0
for j in range(4):
    qp, qm = list(qv), list(qv)
    qp[j] += h
    qm[j] -= h
    gp = vertex_eval_2d((qp[2], qp[3]), (qp[0], qp[1]), -1.0, 1e-9)[1]
    gm = vertex_eval_2d((qm[2], qm[3]), (qm[0], qm[1]), -1.0, 1e-9)[1]
    maxev = max(maxev, abs((gp - gm) / (2 * h) - BV[j]))
check("B", "vertex B == FD d(gap)/du", maxev < 1e-7, f"max err={maxev:.2e}")

# -------------------------------------------- XA. exact-arithmetic noise floor
print("\nXA exact-arithmetic gap reference (the eps*|x| floor MEASURED, not assumed)")
for label, base, seg, npts in (("near origin", 1.0, 1.0, 2000),
                               ("far field (|x|~5e4, facet ~5e2)", 5e4, 5e2, 2000)):
    worst, worst_rel, scale_seen = 0.0, 0.0, 0.0
    for _ in range(npts):
        a = rng.uniform(0, 2 * math.pi)
        c0 = (base * math.cos(a), base * math.sin(a))
        b = rng.uniform(0, 2 * math.pi)
        L = seg * rng.uniform(0.5, 2.0)
        c1 = (c0[0] + L * math.cos(b), c0[1] + L * math.sin(b))
        t = ((c1[0] - c0[0]) / L, (c1[1] - c0[1]) / L)
        u = rng.uniform(0.0, 1.0)
        hgt = seg * rng.uniform(-0.2, 0.2)
        p = (c0[0] + u * L * t[0] - hgt * t[1], c0[1] + u * L * t[1] + hgt * t[0])
        _, _, gd, _ = project_segment_2d(c0, c1, p, 1.0, L)
        ge = gap_exact(c0, c1, p, 1.0)
        err = abs(Decimal(gd) - ge)
        scale = max(abs(c0[0]), abs(c0[1]), abs(c1[0]), abs(c1[1]), abs(p[0]), abs(p[1]))
        scale_seen = max(scale_seen, scale)
        worst = max(worst, float(err))
        worst_rel = max(worst_rel, float(err) / scale)
    floor = EPS * scale_seen
    check("XA", f"{label}: |gap_double - gap_exact| <= 8*eps*|x|",
          worst <= 8.0 * floor,
          f"abs drift={worst:.3e}  eps*|x|={floor:.3e}  ratio={worst/floor:.2f}  "
          f"rel={worst_rel:.3e}")

# the two algebraically identical gap forms must agree to the same floor
worst2 = 0.0
for _ in range(2000):
    a = rng.uniform(0, 2 * math.pi)
    c0 = (5e4 * math.cos(a), 5e4 * math.sin(a))
    b = rng.uniform(0, 2 * math.pi)
    L = 5e2
    c1 = (c0[0] + L * math.cos(b), c0[1] + L * math.sin(b))
    p = (c0[0] + rng.uniform(-1e3, 1e3), c0[1] + rng.uniform(-1e3, 1e3))
    _, _, gd, nn = project_segment_2d(c0, c1, p, 1.0, L)
    # line form: g == (xs - X0).n  (identical because x(xi)-X0 is parallel to t, t.n == 0)
    gl = (p[0] - c0[0]) * nn[0] + (p[1] - c0[1]) * nn[1]
    worst2 = max(worst2, abs(gd - gl))
check("XA", "shape-function gap == line-form gap (independent identity)",
      worst2 <= 8.0 * EPS * 5e4, f"max diff={worst2:.3e}")

# ============================================================== corner falsifiers
# Canonical convex corner: material in the third quadrant, surface (-1,0)->(0,0)->(0,-1)
CVX = ((-1.0, 0.0), (0.0, 0.0), (0.0, -1.0), 1.0)
# Canonical concave (re-entrant) corner: free space is the second quadrant,
# surface (-1,0)->(0,0)->(0,1); the material is everything else nearby.
CCV = ((-1.0, 0.0), (0.0, 0.0), (0.0, 1.0), 1.0)

print("\nCA corner impact -- refused by BOTH segments, owned by the vertex")
for name, (Xp, Xv, Xn, sg), probe, want_corner, want_side in (
        ("convex", CVX, (0.3, 0.4), C_CONVEX, +1.0),
        ("concave", CCV, (0.3, -0.4), C_CONCAVE, -1.0)):
    tP = (Xv[0] - Xp[0], Xv[1] - Xp[1])
    tN = (Xn[0] - Xv[0], Xn[1] - Xv[1])
    Lref = min(math.hypot(*tP), math.hypot(*tN))
    stP = project_segment_2d(Xp, Xv, probe, sg, Lref)[0]
    stN = project_segment_2d(Xv, Xn, probe, sg, Lref)[0]
    w = vertex_wedge_2d(tP, tN, sg, probe, Xv, Lref)
    stv, gv, nv = vertex_eval_2d(Xv, probe, w["sideSign"], TAU_SEG * Lref)
    d = math.hypot(probe[0] - Xv[0], probe[1] - Xv[1])
    check("CA", f"{name}: prev segment REFUSES (OUT_HIGH)", stP == SEG_HIGH)
    check("CA", f"{name}: next segment REFUSES (OUT_LOW)", stN == SEG_LOW)
    check("CA", f"{name}: vertexWedge2D accepts", w["inWedge"] == 1)
    check("CA", f"{name}: corner classified", w["corner"] == want_corner)
    check("CA", f"{name}: sideSign == -corner", w["sideSign"] == want_side
          and w["sideSign"] == -float(w["corner"]))
    check("CA", f"{name}: |gap| == ||xs-Xv||, sign == sideSign",
          stv == VTX_OK and abs(abs(gv) - d) < 1e-15 and math.copysign(1, gv) == want_side,
          f"gap={gv:+.6f} d={d:.6f}")
    # independent geometric cross-check: in the wedge both adjacent outward normals
    # give the SAME signed offset, and that sign IS the side.
    nP, nN = outward_normal(tP, sg), outward_normal(tN, sg)
    dP = (probe[0] - Xv[0]) * nP[0] + (probe[1] - Xv[1]) * nP[1]
    dN = (probe[0] - Xv[0]) * nN[0] + (probe[1] - Xv[1]) * nN[1]
    check("CA", f"{name}: independent half-plane test agrees with sideSign",
          math.copysign(1, dP) == want_side and math.copysign(1, dN) == want_side)
    # brute-force: the vertex really is the closest polyline feature
    best = min(min(math.hypot(probe[0] - (Xp[0] + u * tP[0]), probe[1] - (Xp[1] + u * tP[1]))
                   for u in np.linspace(0, 1, 4001)),
               min(math.hypot(probe[0] - (Xv[0] + u * tN[0]), probe[1] - (Xv[1] + u * tN[1]))
                   for u in np.linspace(0, 1, 4001)))
    check("CA", f"{name}: vertex IS the closest master feature (brute force)",
          abs(best - d) < 1e-12, f"brute={best:.12f} d={d:.12f}")

# out of bounds on BOTH segments but at the FAR ends is the patch boundary, not a corner
Xp, Xv, Xn, sg = CVX
tP = (Xv[0] - Xp[0], Xv[1] - Xp[1])
tN = (Xn[0] - Xv[0], Xn[1] - Xv[1])
far = (-3.0, -3.0)
check("CA", "far-end double refusal is NOT a wedge (owned by nobody, correctly)",
      project_segment_2d(Xp, Xv, far, sg, 1.0)[0] == SEG_LOW
      and project_segment_2d(Xv, Xn, far, sg, 1.0)[0] == SEG_HIGH
      and vertex_wedge_2d(tP, tN, sg, far, Xv, 1.0)["inWedge"] == 0)

print("\nCB slide around a corner at constant physical distance (handoff + continuity)")
for name, (Xp, Xv, Xn, sg), off in (("convex exterior", CVX, +0.25),
                                    ("concave interior", CCV, -0.25)):
    pts, ds, Larc, Ltot = offset_sweep(Xp, Xv, Xn, sg, off, 0.8, 2000)
    Lref = 1.0
    committed = 0.0
    seq, gaps, normals, nown = [], [], [], []
    for p in pts:
        w = vertex_wedge_2d((Xv[0] - Xp[0], Xv[1] - Xp[1]), (Xn[0] - Xv[0], Xn[1] - Xv[1]),
                            sg, p, Xv, Lref)
        if w["inWedge"] and committed == 0.0:
            committed = w["sideSign"]          # COMMIT once, at first capture
        no, which, g, n = owner_of(Xp, Xv, Xn, sg, Lref, p, committed)
        seq.append(which)
        gaps.append(g)
        normals.append(n)
        nown.append(no)
    rank = {0: 0, 2: 1, 1: 2}
    ranks = [rank[w] for w in seq if w >= 0]
    maxgap = max(abs(g - off) for g in gaps)
    maxang = max(ang_between(normals[i - 1], normals[i]) for i in range(1, len(normals)))
    check("CB", f"{name}: the ordered rule resolves ONE owner at all {len(pts)} stations",
          all(w >= 0 for w in seq) and all(v >= 1 for v in nown))
    check("CB", f"{name}: handoff order segment -> vertex -> segment",
          ranks == sorted(ranks) and 1 in ranks)
    check("CB", f"{name}: gap constant (NO jump) along the offset curve",
          maxgap < 1e-14, f"max |gap - {off}| = {maxgap:.3e}")
    check("CB", f"{name}: normal rotation per step <= ds/R (continuous)",
          maxang <= ds / abs(off) * (1 + 1e-9),
          f"max step={maxang:.4e} bound={ds/abs(off):.4e}")
    # SEAM C0 test done properly: evaluate BOTH candidate features at the SAME point
    # sitting exactly on the wedge boundary (aPrev == 0, i.e. xi_prev == 1, and then
    # aNext == 0, i.e. xi_next == 0). A vertex pair is only an extension of its
    # neighbours if the two agree there.
    tP = (Xv[0] - Xp[0], Xv[1] - Xp[1])
    tN = (Xn[0] - Xv[0], Xn[1] - Xv[1])
    nPv, nNv = outward_normal(tP, sg), outward_normal(tN, sg)
    side = -1.0 if off < 0 else 1.0
    for lbl, seam, seg_args in (("prev", (Xv[0] + off * nPv[0], Xv[1] + off * nPv[1]),
                                 (Xp, Xv)),
                                ("next", (Xv[0] + off * nNv[0], Xv[1] + off * nNv[1]),
                                 (Xv, Xn))):
        _, _, gs, ns = project_segment_2d(seg_args[0], seg_args[1], seam, sg, 1.0)
        _, gv, nv = vertex_eval_2d(Xv, seam, side, TAU_SEG)
        check("CB", f"{name}: {lbl} seam -- segment and vertex return the SAME gap",
              abs(gv - gs) <= 1e-15 * max(1.0, abs(gs)), f"|dgap|={abs(gv-gs):.3e}")
        check("CB", f"{name}: {lbl} seam -- segment and vertex return the SAME normal",
              ang_between(ns, nv) <= 1e-12, f"angle={ang_between(ns, nv):.3e} rad")

print("\nCC no double count -- the MEASURED handoff band and the ordered ownership rule")
# The naive reading ("segment in-bounds XOR slacked wedge") is ALGEBRAICALLY exclusive
# and exhaustive, but NOT in floating point: xi_prev is formed from X_prev so it sits at
# ~1 and can only resolve the overshoot to ulp(1) = 2.22e-16, while aPrev resolves it
# exactly. This family MEASURES the resulting band and shows the ordered rule closes it.
band_res = float(np.spacing(1.0))
xor_two, xor_none, rule_bad, nsta = 0, 0, 0, 0
band_prev, band_next = 0.0, 0.0
overlap_dgap, overlap_dang, overlap_hits = 0.0, 0.0, 0
overlap_bad, overlap_bnd = 0.0, 0.0
for _ in range(200):
    Xv = (rng.uniform(-5, 5), rng.uniform(-5, 5))
    a0 = rng.uniform(0, 2 * math.pi)
    turn = float(rng.choice([-1.0, 1.0])) * rng.uniform(0.3, 2.5)
    lp, ln = rng.uniform(0.3, 3.0), rng.uniform(0.3, 3.0)
    tP = (lp * math.cos(a0), lp * math.sin(a0))
    tN = (ln * math.cos(a0 + turn), ln * math.sin(a0 + turn))
    Xp = (Xv[0] - tP[0], Xv[1] - tP[1])
    Xn = (Xv[0] + tN[0], Xv[1] + tN[1])
    tPr = (Xv[0] - Xp[0], Xv[1] - Xp[1])          # as a handler recomputes it
    tNr = (Xn[0] - Xv[0], Xn[1] - Xv[1])
    sg = float(rng.choice([-1.0, 1.0]))
    Lref = min(math.hypot(*tPr), math.hypot(*tNr))
    # scan a window across each seam that is ~4500x wider than ulp(1)
    for alpha0, beta0, which_seam in ((TOL_IN, -0.5, "prev"), (0.5, -TOL_IN, "next")):
        for dd in np.linspace(-1e-12, 1e-12, 41):
            alpha = alpha0 + (dd if which_seam == "prev" else 0.0)
            beta = beta0 + (dd if which_seam == "next" else 0.0)
            p = point_from_params(tPr, tNr, Xv, alpha, beta)
            nsta += 1
            nx = owner_xor(Xp, Xv, Xn, sg, Lref, p)
            if nx > 1:
                xor_two += 1
            if nx == 0:
                xor_none += 1
            _, which, _, _ = owner_of(Xp, Xv, Xn, sg, Lref, p)
            if which < 0:
                rule_bad += 1
            xiP = project_segment_2d(Xp, Xv, p, sg, Lref)[1]
            xiN = project_segment_2d(Xv, Xn, p, sg, Lref)[1]
            w = vertex_wedge_2d(tPr, tNr, sg, p, Xv, Lref)
            band_prev = max(band_prev, abs((xiP - 1.0) - w["aPrev"]))
            band_next = max(band_next, abs(xiN - w["aNext"]))
check("CC", f"seam scan ({nsta} stations): the ordered rule ALWAYS resolves an owner",
      rule_bad == 0, f"{rule_bad} unowned")
check("CC", "the naive slacked XOR is NOT floating-point exact (documented, not assumed)",
      xor_two > 0 or xor_none > 0,
      f"{xor_two} double-owner + {xor_none} zero-owner stations of {nsta}")
check("CC", "the band is a few ulp of the O(1) parametric coordinate, NOT a modelling "
      "tolerance (no tolerance value can close it)",
      max(band_prev, band_next) <= 32.0 * band_res
      and max(band_prev, band_next) <= 1e-5 * TOL_IN,
      f"band_prev={band_prev:.3e} band_next={band_next:.3e} = "
      f"{max(band_prev, band_next)/band_res:.1f} ulp(1) = "
      f"{max(band_prev, band_next)/TOL_IN:.2e} of tolIn")
HANDOFF_BAND = (band_prev, band_next)

# In the overlap band BOTH candidate features are legitimate: measure how much the
# choice costs physically, at a realistic penetration.
for _ in range(4000):
    Xv = (rng.uniform(-5, 5), rng.uniform(-5, 5))
    a0 = rng.uniform(0, 2 * math.pi)
    turn = float(rng.choice([-1.0, 1.0])) * rng.uniform(0.3, 2.5)
    lp, ln = rng.uniform(0.3, 3.0), rng.uniform(0.3, 3.0)
    tP = (lp * math.cos(a0), lp * math.sin(a0))
    tN = (ln * math.cos(a0 + turn), ln * math.sin(a0 + turn))
    Xp = (Xv[0] - tP[0], Xv[1] - tP[1])
    Xn = (Xv[0] + tN[0], Xv[1] + tN[1])
    sg = float(rng.choice([-1.0, 1.0]))
    Lref = min(lp, ln)
    # a slave inside the overlap 0 < aPrev <= tolIn, at a physically meaningful offset
    alpha = TOL_IN * rng.uniform(0.0, 1.0)
    beta = -rng.uniform(0.05, 0.5)
    p = point_from_params(tP, tN, Xv, alpha, beta)
    stP, _, gs, ns = project_segment_2d(Xp, Xv, p, sg, Lref)
    w = vertex_wedge_2d(tP, tN, sg, p, Xv, Lref)
    if stP != SEG_IN or not (w["aPrev"] > 0.0 and w["aNext"] < 0.0):
        continue
    _, gv, nv = vertex_eval_2d(Xv, p, w["sideSign"], TAU_SEG * Lref)
    overlap_hits += 1
    ang = ang_between(ns, nv)
    # ANALYTIC bound: the slave sits gs along n_P and a = aPrev*L_P along t_P past the
    # vertex, so the two normals differ by exactly atan(a/|gs|) -- i.e. the disagreement
    # is O(tolIn * L / |gap|), controlled by the PARAMETRIC slack, and the gap
    # disagreement is second order in it. Assert against that bound, not a magic number.
    bound = math.atan(w["aPrev"] * math.hypot(*tP) / abs(gs))
    if bound > 0.0:
        overlap_bad = max(overlap_bad, ang / bound)      # worst measured/analytic ratio
    overlap_dgap = max(overlap_dgap, abs(gv - gs) / max(1e-30, abs(gs)))
    overlap_dang = max(overlap_dang, ang)
    overlap_bnd = max(overlap_bnd, bound)
check("CC", f"in the overlap band ({overlap_hits} hits) the normal disagreement obeys "
      "the analytic atan(tolIn*L/|gap|) bound", overlap_hits > 500 and overlap_bad <= 1.001,
      f"max angle={overlap_dang:.3e} rad, worst analytic bound={overlap_bnd:.3e} rad, "
      f"worst measured/analytic ratio={overlap_bad:.6f}")
check("CC", "in the overlap band the two candidates agree on the GAP to 2nd order -- "
      "which one wins is physically immaterial", overlap_dgap < 1e-10,
      f"max relative gap difference={overlap_dgap:.3e}")

print("\nCD the CONCAVE decision -- does a concave vertex need its own rule?")
# Independent material test at a corner (NO kernel quantities involved):
#   convex  body = the INTERSECTION of the two half-planes  -> inside <=> dP<0 AND dN<0
#   concave body = the UNION of the two half-planes         -> inside <=> dP<0 OR  dN<0
# (The "both dots share the side sign" shortcut only holds for |turn| <= 90 deg -- the
# wedge fan has angular width |turn|, so for a blunter corner a wedge direction can be
# more than 90 deg from ONE of the two bounding normals. The bisector test in the kernel
# is what stays correct for every |turn| < 180 deg; family CE exercises that.)
n_cvx_pos, n_cvx, n_ccv_neg, n_ccv, n_ccv_refused, n_fan = 0, 0, 0, 0, 0, 0
for _ in range(20000):
    Xv = (rng.uniform(-5, 5), rng.uniform(-5, 5))
    a0 = rng.uniform(0, 2 * math.pi)
    sg = float(rng.choice([-1.0, 1.0]))
    for want in (C_CONVEX, C_CONCAVE):
        # turn sign chosen so that sigma*cross(tP,tN) has the wanted sign
        mag = rng.uniform(0.2, 2.8)
        turn = mag * (1.0 if want == C_CONCAVE else -1.0) * sg
        lp, ln = rng.uniform(0.3, 3.0), rng.uniform(0.3, 3.0)
        tP = (lp * math.cos(a0), lp * math.sin(a0))
        tN = (ln * math.cos(a0 + turn), ln * math.sin(a0 + turn))
        Xp = (Xv[0] - tP[0], Xv[1] - tP[1])
        Xn = (Xv[0] + tN[0], Xv[1] + tN[1])
        Lref = min(lp, ln)
        p = point_from_params(tP, tN, Xv, rng.uniform(0.02, 0.6), -rng.uniform(0.02, 0.6))
        w = vertex_wedge_2d(tP, tN, sg, p, Xv, Lref)
        if not w["inWedge"] or w["corner"] != want:
            continue
        _, g, _ = vertex_eval_2d(Xv, p, w["sideSign"], TAU_SEG * Lref)
        nP, nN = outward_normal(tP, sg), outward_normal(tN, sg)
        dP = (p[0] - Xv[0]) * nP[0] + (p[1] - Xv[1]) * nP[1]
        dN = (p[0] - Xv[0]) * nN[0] + (p[1] - Xv[1]) * nN[1]
        # the wedge fan is the cone of half-width |turn|/2 about sideSign*bisector --
        # an independent characterisation of the region the kernel classifies
        nb = (nP[0] + nN[0], nP[1] + nN[1])
        if math.hypot(*nb) > 1e-9:
            bis = unit(nb)
            rr = unit((p[0] - Xv[0], p[1] - Xv[1]))
            half = ang_between(rr, (w["sideSign"] * bis[0], w["sideSign"] * bis[1]))
            if half <= abs(mag) / 2 + 1e-9:
                n_fan += 1
        if want == C_CONVEX:
            n_cvx += 1
            # OUTSIDE a convex body (the INTERSECTION of the two half-planes)
            if g > 0.0 and not (dP < 0.0 and dN < 0.0):
                n_cvx_pos += 1
        else:
            n_ccv += 1
            # INSIDE a re-entrant body (the UNION of the two half-planes)
            if g < 0.0 and (dP < 0.0 or dN < 0.0):
                n_ccv_neg += 1
            if (project_segment_2d(Xp, Xv, p, sg, Lref)[0] == SEG_HIGH
                    and project_segment_2d(Xv, Xn, p, sg, Lref)[0] == SEG_LOW):
                n_ccv_refused += 1
check("CD", "CONVEX wedge: gap > 0 and the point is OUTSIDE the material, always",
      n_cvx > 5000 and n_cvx_pos == n_cvx, f"{n_cvx_pos}/{n_cvx}")
check("CD", "CONCAVE wedge: gap < 0 and the point is INSIDE the material, always",
      n_ccv > 5000 and n_ccv_neg == n_ccv, f"{n_ccv_neg}/{n_ccv}")
check("CD", "CONCAVE wedge: BOTH adjacent segments refuse it (the coverage hole)",
      n_ccv_refused == n_ccv, f"{n_ccv_refused}/{n_ccv}")
check("CD", "the wedge IS the cone of half-width |turn|/2 about sideSign*bisector",
      n_fan == n_cvx + n_ccv, f"{n_fan}/{n_cvx + n_ccv}")

# quantify the hole policy A ("concave never activates") would leave, at several depths
hole_rows = []
for R in (0.05, 0.1, 0.25, 0.5):
    Xp, Xv, Xn, sg = CCV
    pts, ds, Larc, Ltot = offset_sweep(Xp, Xv, Xn, sg, -R, 0.8, 4000)
    unowned = 0
    for p in pts:
        w = vertex_wedge_2d((Xv[0] - Xp[0], Xv[1] - Xp[1]), (Xn[0] - Xv[0], Xn[1] - Xv[1]),
                            sg, p, Xv, 1.0)
        stP = project_segment_2d(Xp, Xv, p, sg, 1.0)[0]
        stN = project_segment_2d(Xv, Xn, p, sg, 1.0)[0]
        if stP != SEG_IN and stN != SEG_IN and w["inWedge"]:
            unowned += 1          # policy A would drop exactly these
    hole_rows.append((R, Larc, Ltot, 100.0 * unowned / len(pts)))
    check("CD", f"policy 'concave never activates' leaves a hole at depth R={R}",
          unowned > 0,
          f"unowned arc = {Larc:.4f} ({100.0*unowned/len(pts):.1f}% of the sweep), "
          f"penetration {R} with ZERO restoring force")
check("CD", "the hole measure scales with the penetration depth (it is not a sliver)",
      all(abs(r[1] - r[0] * math.pi / 2) < 1e-9 for r in hole_rows),
      "arc == R*|turn| = R*pi/2 exactly at every depth")
CONCAVE_DECISION = "ACTIVATE (concave vertices need the rule; convex ones never carry force)"

print("\nCE side-sign robustness -- the bisector test vs sign(turnSin) at near-straight")
worst_disagree, n_near, n_spike_defer = 0, 0, 0
for _ in range(20000):
    Xv = (rng.uniform(-3, 3), rng.uniform(-3, 3))
    a0 = rng.uniform(0, 2 * math.pi)
    sg = float(rng.choice([-1.0, 1.0]))
    turn = float(rng.choice([-1.0, 1.0])) * 10.0 ** rng.uniform(-11, 0.4)
    lp, ln = rng.uniform(0.3, 3.0), rng.uniform(0.3, 3.0)
    tP = (lp * math.cos(a0), lp * math.sin(a0))
    tN = (ln * math.cos(a0 + turn), ln * math.sin(a0 + turn))
    Lref = min(lp, ln)
    p = point_from_params(tP, tN, Xv, rng.uniform(0.02, 0.6), -rng.uniform(0.02, 0.6))
    w = vertex_wedge_2d(tP, tN, sg, p, Xv, Lref)
    if not w["inWedge"]:
        continue
    n_near += 1
    if w["sideSign"] != -float(w["corner"]):
        worst_disagree += 1
check("CE", "inside the wedge sideSign == -corner (incl. turn angles down to 1e-11)",
      n_near > 3000 and worst_disagree == 0, f"{n_near} wedge hits, {worst_disagree} disagree")
# a fold-back "spike" corner (antiparallel outward normals) must DEFER, never guess
tP = (1.0, 0.0)
tN = (-1.0, 0.0)
w = vertex_wedge_2d(tP, tN, 1.0, (0.5, 0.3), (0.0, 0.0), 1.0)
check("CE", "fold-back spike corner -> sideSign 0 (DEFER, never guess a sign)",
      w["sideSign"] == 0.0)
check("CE", "slave exactly at the vertex -> sideSign 0",
      vertex_wedge_2d((1.0, 0.0), (0.0, 1.0), 1.0, (0.0, 0.0), (0.0, 0.0), 1.0)["sideSign"] == 0.0)

# =========================================================== randomized C++ parity
NOCPP = "--nocpp" in sys.argv


def build_driver(tmp):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, os.pardir, os.pardir))
    inc = os.path.join(root, "SRC", "domain", "contact")
    src = os.path.join(here, "t1a_cpp_check.cpp")
    exe = os.path.join(tmp, "t1a_cpp_check" + (".exe" if os.name == "nt" else ""))
    cmd = ["g++", "-O2", "-ffp-contract=off", "-I", inc, src, "-o", exe]
    subprocess.run(cmd, check=True)
    return exe


def hx(v):
    return float(v).hex()


class Trials(object):
    """Ordered trial tape: (family, kind, args) in, parsed C++ outputs out."""

    def __init__(self):
        self.recs = []

    def add(self, fam, kind, args):
        self.recs.append((fam, kind, [float(a) for a in args]))

    def render(self):
        out = []
        for fam, kind, args in self.recs:
            out.append(kind + " " + " ".join(hx(a) for a in args))
        return "\n".join(out) + "\n"


def parse_out(text, n):
    rows = []
    for line in text.splitlines():
        p = line.split()
        if not p:
            continue
        kind = p[0]
        if kind == "S":
            rows.append((kind, int(p[1]), [float.fromhex(x) for x in p[2:]]))
        elif kind == "W":
            rows.append((kind, (int(p[1]), int(p[2])), [float.fromhex(x) for x in p[3:]]))
        elif kind == "V":
            rows.append((kind, int(p[1]), [float.fromhex(x) for x in p[2:]]))
        elif kind == "R":
            rows.append((kind, 0, [float.fromhex(p[1])]))
    if len(rows) != n:
        raise RuntimeError(f"driver returned {len(rows)} records, expected {n}")
    return rows


def rand_seg(scale=1.0, lo=0.1, hi=5.0):
    X0 = (rng.uniform(-10, 10) * scale, rng.uniform(-10, 10) * scale)
    a = rng.uniform(0, 2 * math.pi)
    L = rng.uniform(lo, hi) * scale
    X1 = (X0[0] + L * math.cos(a), X0[1] + L * math.sin(a))
    return X0, X1, L, (math.cos(a), math.sin(a))


def point_at(X0, t, L, xi, hgt):
    return (X0[0] + xi * L * t[0] - hgt * t[1], X0[1] + xi * L * t[1] + hgt * t[0])


if not NOCPP:
    print("\nP  randomized C++ parity (hexfloat round-trip, relative tolerances)")
    T = Trials()

    # F1 -- in-bounds projections
    for _ in range(30000):
        X0, X1, L, t = rand_seg()
        p = point_at(X0, t, L, rng.uniform(0.0, 1.0), rng.uniform(-0.5, 0.5) * L)
        T.add("F1", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1],
                          float(rng.choice([-1.0, 1.0])), L, TOL_IN, TAU_SEG])
    # F2 -- refusals, including a dense band straddling the tolIn slack boundary
    for _ in range(20000):
        X0, X1, L, t = rand_seg()
        xi = rng.uniform(1.001, 4.0) if rng.random() < 0.5 else rng.uniform(-3.0, -0.001)
        p = point_at(X0, t, L, xi, rng.uniform(-0.5, 0.5) * L)
        T.add("F2", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1],
                          float(rng.choice([-1.0, 1.0])), L, TOL_IN, TAU_SEG])
    for _ in range(10000):
        X0, X1, L, t = rand_seg()
        edge = 1.0 if rng.random() < 0.5 else 0.0
        xi = edge + rng.choice([1.0, -1.0]) * TOL_IN * rng.uniform(0.0, 3.0)
        p = point_at(X0, t, L, xi, rng.uniform(-0.1, 0.1) * L)
        T.add("F2b", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1],
                           float(rng.choice([-1.0, 1.0])), L, TOL_IN, TAU_SEG])
    # F3 -- far from the origin (coordinates ~5e4, facet ~5e2)
    for _ in range(20000):
        a = rng.uniform(0, 2 * math.pi)
        X0 = (5e4 * math.cos(a), 5e4 * math.sin(a))
        b = rng.uniform(0, 2 * math.pi)
        L = 5e2 * rng.uniform(0.5, 2.0)
        t = (math.cos(b), math.sin(b))
        X1 = (X0[0] + L * t[0], X0[1] + L * t[1])
        p = point_at(X0, t, L, rng.uniform(-0.2, 1.2), rng.uniform(-0.2, 0.2) * L)
        T.add("F3", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1],
                          float(rng.choice([-1.0, 1.0])), L, TOL_IN, TAU_SEG])
    # F4 -- zero-length + bad-Lref refusal (factor straddles the tauSeg threshold)
    for _ in range(5000):
        X0 = (rng.uniform(-5, 5), rng.uniform(-5, 5))
        Lref = rng.uniform(0.2, 4.0)
        fac = 10.0 ** rng.uniform(-2, 2)
        L = TAU_SEG * Lref * fac
        a = rng.uniform(0, 2 * math.pi)
        X1 = (X0[0] + L * math.cos(a), X0[1] + L * math.sin(a))
        p = (X0[0] + rng.uniform(-1, 1), X0[1] + rng.uniform(-1, 1))
        T.add("F4", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1], 1.0, Lref, TOL_IN, TAU_SEG])
    for badLref in (0.0, -1.0, float("nan")):
        for _ in range(200):
            X0, X1, L, t = rand_seg()
            p = point_at(X0, t, L, 0.5, 0.1 * L)
            T.add("F4b", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1], 1.0,
                               badLref, TOL_IN, TAU_SEG])
    # F5 -- sigma flip + reversed winding (3 records per configuration, in order)
    for _ in range(15000):
        X0, X1, L, t = rand_seg()
        p = point_at(X0, t, L, rng.uniform(-0.2, 1.2), rng.uniform(-1, 1) * L)
        s = float(rng.choice([-1.0, 1.0]))
        T.add("F5", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1], s, L, TOL_IN, TAU_SEG])
        T.add("F5", "S", [X0[0], X0[1], X1[0], X1[1], p[0], p[1], -s, L, TOL_IN, TAU_SEG])
        T.add("F5", "S", [X1[0], X1[1], X0[0], X0[1], p[0], p[1], -s, L, TOL_IN, TAU_SEG])
    # F6 -- vertex pairs incl. the coincidence refusal band
    for _ in range(20000):
        Xv = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        Lref = rng.uniform(0.2, 4.0)
        tolLen = TAU_SEG * Lref
        d = tolLen * 10.0 ** rng.uniform(-2, 8)
        a = rng.uniform(0, 2 * math.pi)
        p = (Xv[0] + d * math.cos(a), Xv[1] + d * math.sin(a))
        T.add("F6", "V", [Xv[0], Xv[1], p[0], p[1], float(rng.choice([-1.0, 1.0])), tolLen])
    # F7 -- wedge classification over convex / concave / near-straight / spike corners
    for _ in range(30000):
        Xv = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        a0 = rng.uniform(0, 2 * math.pi)
        r = rng.random()
        if r < 0.15:
            turn = float(rng.choice([-1.0, 1.0])) * 10.0 ** rng.uniform(-12, -6)
        elif r < 0.25:
            turn = float(rng.choice([-1.0, 1.0])) * (math.pi - 10.0 ** rng.uniform(-12, -3))
        else:
            turn = float(rng.choice([-1.0, 1.0])) * rng.uniform(0.05, 3.0)
        lp, ln = rng.uniform(0.2, 4.0), rng.uniform(0.2, 4.0)
        tP = (lp * math.cos(a0), lp * math.sin(a0))
        tN = (ln * math.cos(a0 + turn), ln * math.sin(a0 + turn))
        Lref = min(lp, ln)
        if rng.random() < 0.7:
            try:
                p = point_from_params(tP, tN, Xv,
                                      rng.uniform(-0.6, 0.6), rng.uniform(-0.6, 0.6))
            except np.linalg.LinAlgError:
                continue
        else:
            p = (Xv[0] + rng.uniform(-4, 4), Xv[1] + rng.uniform(-4, 4))
        T.add("F7", "W", [tP[0], tP[1], tN[0], tN[1], float(rng.choice([-1.0, 1.0])),
                          p[0], p[1], Xv[0], Xv[1], Lref, TOL_IN, TAU_SEG, TAU_PERP])
    # F9 -- orientation from a datum, incl. near-tangent data
    for _ in range(10000):
        a = rng.uniform(0, 2 * math.pi)
        L = rng.uniform(0.1, 5.0)
        t = (L * math.cos(a), L * math.sin(a))
        if rng.random() < 0.3:
            eps_ang = 10.0 ** rng.uniform(-16, -10)
            b = a + float(rng.choice([-1.0, 1.0])) * eps_ang
        else:
            b = rng.uniform(0, 2 * math.pi)
        Rl = rng.uniform(0.01, 100.0)
        T.add("F9", "R", [t[0], t[1], Rl * math.cos(b), Rl * math.sin(b), TAU_PERP])

    tmp = tempfile.mkdtemp(prefix="adr85_t1a_")
    exe = build_driver(tmp)
    tape = os.path.join(tmp, "trials.txt")
    with open(tape, "w") as fh:
        fh.write(T.render())
    with open(tape, "r") as fh:
        res = subprocess.run([exe, "--trials"], stdin=fh, stdout=subprocess.PIPE,
                             check=True)
    rows = parse_out(res.stdout.decode("ascii"), len(T.recs))
    print(f"   {len(rows)} records through the g++ driver "
          f"(-O2 -ffp-contract=off), tape={os.path.getsize(tape)//1024} KiB")

    # ---- compare, per family, at RELATIVE tolerances
    TOL_PAR, TOL_N, TOL_GAP = 1e-12, 1e-13, 1e-13
    stats = {}
    status_bad = {}
    for (fam, kind, args), row in zip(T.recs, rows):
        st = stats.setdefault(fam, dict(n=0, dxi=0.0, dgap=0.0, dn=0.0, dB=0.0,
                                        dmisc=0.0, dabs=0.0))
        st["n"] += 1
        if kind == "S":
            X0 = (args[0], args[1]); X1 = (args[2], args[3]); p = (args[4], args[5])
            sg, Lref, tolIn, tauSeg = args[6], args[7], args[8], args[9]
            pst, pxi, pgap, pn = project_segment_2d(X0, X1, p, sg, Lref, tolIn, tauSeg)
            pB = b_operator_segment_2d(*shape2d(pxi), n=pn)
            cst, vals = row[1], row[2]
            if cst != pst:
                status_bad[fam] = status_bad.get(fam, 0) + 1
            scale = max(abs(X0[0]), abs(X0[1]), abs(X1[0]), abs(X1[1]),
                        abs(p[0]), abs(p[1]), 1.0)
            st["dxi"] = max(st["dxi"], abs(vals[0] - pxi) / (1.0 + abs(pxi)))
            st["dabs"] = max(st["dabs"], abs(vals[1] - pgap))
            st["dgap"] = max(st["dgap"], abs(vals[1] - pgap) / scale)
            st["dn"] = max(st["dn"], max(abs(vals[2] - pn[0]), abs(vals[3] - pn[1])))
            st["dB"] = max(st["dB"], max(abs(vals[4 + i] - pB[i]) for i in range(6)))
        elif kind == "W":
            tP = (args[0], args[1]); tN = (args[2], args[3]); sg = args[4]
            p = (args[5], args[6]); Xv = (args[7], args[8])
            pw = vertex_wedge_2d(tP, tN, sg, p, Xv, args[9], args[10], args[11], args[12])
            (cw, cc), vals = row[1], row[2]
            if cw != pw["inWedge"] or cc != pw["corner"] or vals[0] != pw["sideSign"]:
                status_bad[fam] = status_bad.get(fam, 0) + 1
            st["dmisc"] = max(st["dmisc"], abs(vals[1] - pw["turnSin"]))
            st["dxi"] = max(st["dxi"],
                            abs(vals[2] - pw["aPrev"]) / (1.0 + abs(pw["aPrev"])),
                            abs(vals[3] - pw["aNext"]) / (1.0 + abs(pw["aNext"])))
        elif kind == "V":
            Xv = (args[0], args[1]); p = (args[2], args[3])
            pst, pgap, pn = vertex_eval_2d(Xv, p, args[4], args[5])
            pB = b_operator_vertex_2d(pn)
            cst, vals = row[1], row[2]
            if cst != pst:
                status_bad[fam] = status_bad.get(fam, 0) + 1
            scale = max(abs(Xv[0]), abs(Xv[1]), abs(p[0]), abs(p[1]), 1.0)
            st["dabs"] = max(st["dabs"], abs(vals[0] - pgap))
            st["dgap"] = max(st["dgap"], abs(vals[0] - pgap) / scale)
            st["dn"] = max(st["dn"], max(abs(vals[1] - pn[0]), abs(vals[2] - pn[1])))
            st["dB"] = max(st["dB"], max(abs(vals[3 + i] - pB[i]) for i in range(4)))
        else:  # R
            ps = sigma_from_ref_2d((args[0], args[1]), (args[2], args[3]), args[4])
            if row[2][0] != ps:
                status_bad[fam] = status_bad.get(fam, 0) + 1

    for fam in sorted(stats):
        s = stats[fam]
        bad = status_bad.get(fam, 0)
        check("P", f"{fam} ({s['n']} trials): discrete fields EXACT", bad == 0,
              f"{bad} mismatches")
        check("P", f"{fam}: continuous fields within relative bounds",
              s["dxi"] <= TOL_PAR and s["dgap"] <= TOL_GAP and s["dn"] <= TOL_N
              and s["dB"] <= TOL_N and s["dmisc"] <= TOL_N,
              f"dxi/1={s['dxi']:.2e} dgap/|x|={s['dgap']:.2e} dn={s['dn']:.2e} "
              f"dB={s['dB']:.2e} (abs gap drift {s['dabs']:.2e})")

    # ---- family-specific invariants on the C++ answers themselves
    # F5: sigma antisymmetry must be EXACT, and reversing the winding with the sign
    # flipped must reproduce the original normal exactly and xi -> 1-xi.
    idx = [i for i, r in enumerate(T.recs) if r[0] == "F5"]
    anti_bad, wind_bad = 0, 0
    for k in range(0, len(idx), 3):
        a, b, cw = rows[idx[k]][2], rows[idx[k + 1]][2], rows[idx[k + 2]][2]
        if not (b[2] == -a[2] and b[3] == -a[3] and b[1] == -a[1] and b[0] == a[0]):
            anti_bad += 1
        if not (abs(cw[0] - (1.0 - a[0])) <= 1e-12 * (1 + abs(a[0]))
                and abs(cw[2] - a[2]) <= 1e-15 and abs(cw[3] - a[3]) <= 1e-15
                and abs(cw[1] - a[1]) <= 1e-11 * max(1.0, abs(a[1]))):
            wind_bad += 1
    check("P", "F5 sigma flip: n(-s)==-n(s) and gap(-s)==-gap(s) EXACTLY (bit-for-bit)",
          anti_bad == 0, f"{anti_bad}/{len(idx)//3} violations")
    check("P", "F5 reversed winding + flipped sigma reproduces n, gap, and xi->1-xi",
          wind_bad == 0, f"{wind_bad}/{len(idx)//3} violations")

    # F4/F4b: everything below the threshold must be DEGENERATE with zeroed outputs
    deg_bad = 0
    for (fam, kind, args), row in zip(T.recs, rows):
        if fam not in ("F4", "F4b"):
            continue
        L = math.hypot(args[2] - args[0], args[3] - args[1])
        expect_deg = (not (args[7] > 0.0)) or (L < TAU_SEG * args[7])
        if expect_deg and (row[1] != SEG_DEG or any(v != 0.0 for v in row[2][:4])):
            deg_bad += 1
        if (not expect_deg) and row[1] == SEG_DEG:
            deg_bad += 1
    check("P", "F4 zero-length / bad-Lref: refusal exactly at L < tauSeg*Lref, outputs zeroed",
          deg_bad == 0, f"{deg_bad} violations")

    # F6: the coincidence refusal, and gap == sigmaSide*d elsewhere
    v_bad = 0
    for (fam, kind, args), row in zip(T.recs, rows):
        if fam != "F6":
            continue
        d = math.hypot(args[2] - args[0], args[3] - args[1])
        if d < args[5]:
            if row[1] != VTX_DEG or any(v != 0.0 for v in row[2][:3]):
                v_bad += 1
        else:
            if row[1] != VTX_OK or abs(abs(row[2][0]) - d) > 1e-12 * max(1.0, d):
                v_bad += 1
    check("P", "F6 vertex: coincidence refusal + |gap| == ||xs-Xv||", v_bad == 0,
          f"{v_bad} violations")

    # F7: inside the wedge the side sign is EITHER the geometric -corner OR a deliberate
    # DEFERRAL (0). A WRONG sign is the only real defect -- a deferral is the designed
    # fold-back-spike fail-safe, and it fires exactly where `corner` itself becomes a
    # coin flip (|n_P + n_N| -> 0 when the turn approaches pi).
    w_bad, w_hits, w_defer, w_defer_bad = 0, 0, 0, 0
    for (fam, kind, args), row in zip(T.recs, rows):
        if fam != "F7":
            continue
        (inw, cor), vals = row[1], row[2]
        if not inw:
            continue
        w_hits += 1
        if not (vals[2] > args[10] and vals[3] < -args[10]):
            w_bad += 1
        if vals[0] == 0.0:
            w_defer += 1
            ct = ((args[0] * args[2] + args[1] * args[3])
                  / (math.hypot(args[0], args[1]) * math.hypot(args[2], args[3])))
            if ct > -0.9999:                 # a deferral on a WELL-conditioned corner
                w_defer_bad += 1
        elif vals[0] != -float(cor):
            w_bad += 1                       # a WRONG sign: the only real defect
    check("P", f"F7 wedge ({w_hits} hits): sideSign is -corner or a deferral, NEVER the "
          "wrong sign; wedge implies both projections out toward the vertex",
          w_hits > 5000 and w_bad == 0, f"{w_bad} violations")
    check("P", f"F7 deferrals ({w_defer}) occur ONLY at fold-back (near-antiparallel) "
          "corners -- the region where sign(turnSin) is itself a coin flip",
          w_defer_bad == 0, f"{w_defer_bad} deferrals on well-conditioned corners")

# ==================================================================== summary
print("\n" + "=" * 78)
print("PER-FAMILY SUMMARY")
for fam in sorted(_fam):
    n, f = _fam[fam]
    print(f"  {fam:4s}  {n - f:4d}/{n:<4d} pass" + ("" if f == 0 else f"   <-- {f} FAIL"))
print("-" * 78)
print(f"  handoff band (parametric): prev {HANDOFF_BAND[0]:.3e}  next {HANDOFF_BAND[1]:.3e}"
      f"   vs tolIn {TOL_IN:.0e}")
print(f"  concave-vertex decision:   {CONCAVE_DECISION}")
for R, Larc, Ltot, pct in hole_rows:
    print(f"    depth R={R:<5} unowned wedge arc {Larc:.4f} of {Ltot:.4f} ({pct:.1f}%)"
          " if concave vertices are dropped")
print("=" * 78)
print(("ALL PASS" if _fails == 0 else f"{_fails} FAILURE(S)") + "   (ADR-85 T1a oracle)")
print("=" * 78)
sys.exit(1 if _fails else 0)
