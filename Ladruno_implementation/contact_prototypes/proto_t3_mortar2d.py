"""ADR-85 T3 ORACLE -- Ladruno2DMortarKernel (2D mortar: interval clip on the
slave-trace measure, ALM/Uzawa, tie, thickness).

Oracle-first (the fork discipline): an INDEPENDENT reference implementation of the 2D
segment-to-segment mortar integrator, written from ADR-85 SS How/3 ("2D mortar: interval
clip on the slave-trace measure") and SS How/7 (thickness), validated against closed-form
slave-trace integrals and against EXACT rational (fractions.Fraction) quadrature, and only
then used as the contract the C++ 2D mortar integrator must match bit-for-bit-testably.

  ADR-85 SS How/3 (the governing sentences):
      [a,b] = [xiM(s0), xiM(s1)]  n  [0,1]        (sorted)
      b - a <= tau_ov  =>  status = -1 (empty)    tau_ov RELATIVE to min(L_s,L_m)/L_m
      dGamma_s = ( L_m / |t_s . t_m| ) dxi_m      <-- THE SLAVE-TRACE MEASURE
  the 2D analogue of the shipped 3D kernel's  A_phys = A_aux / cos_t
  (SRC/domain/contact/LadrunoMortarKernel.h:460) and of its sliver guard
  A_poly < 1e-12*A_min (:444).

  Integrands, on that measure:
      D_ij   = INT N_i^s N_j^s dGamma_s        (slave-slave  mortar matrix)
      M_ik   = INT N_i^s N_k^m dGamma_s        (slave-master mortar matrix)
      gtil_i = INT N_i^s g     dGamma_s        (weighted gap, NORMAL INSIDE the integral)
  with n = sigma * perp(t_m)/L_m (sigma = +-1, the interface orientation vote -- HANDLER
  business, an INPUT here) and the fork sign convention  penetration <=> g < 0.

  EXACTNESS: for straight 2-node segments the map xi_s -> xi_m is affine and both n and
  the Jacobian are constant per pair, so every integrand is a polynomial of degree <= 2 in
  xi_m and 2-POINT GAUSS IS EXACT. CHECK 2 proves it three ways (2pt vs 50pt vs exact
  rational), not by assertion.

WHY THE TILTED CASE EXISTS (ADR-85 SS How/3, verbatim rationale):
  "Omitting the 1/|t_s.t_m| Jacobian would bias D/M/gtil by cos(theta) for every
   non-parallel pair -- and the flat patch test could never catch it (cos = 1 there),
   so the T3 oracle includes a tilted non-parallel-overlap case."
  CHECK 3 therefore runs a DELIBERATELY BROKEN twin (jac_mode="no_align", J = L_m) beside
  the correct kernel at cos(theta) in {0.9, 0.7, 0.5} and asserts the analytic check
  CATCHES it by exactly the cos(theta) factor -- and that at cos(theta) = 1 the two are
  BIT-IDENTICAL, i.e. a flat patch test has ZERO discriminating power against this bug.

Gate coverage (ADR-85 G-T3):
  (a) machine-precision uniform-pressure transfer across a NON-MATCHING 2D interface -- CHECK 4
  (b) augmentation drives penetration below an epsN-INDEPENDENT augTol             -- CHECK 7
  (c) tie patch exact on non-matching meshes incl. a shared-node T-junction
      (the ADR-41 C4 crux)                                                         -- CHECK 5
  (f) a -thickness 2 twin scales tractions by exactly 2 and leaves geometry untouched -- CHECK 8
  ((d) routing / WARN_2D_PERP_NO_NTS and (e) mortar block-on-incline are handler- and
   friction-lane business; this oracle supplies the REFUSE_ALIGN classification (d)
   depends on -- see A2 -- and the h-scaled unified cone (e) depends on, gating both at
   kernel level.)

------------------------------------------------------------------------------
SPEC AMBIGUITIES RESOLVED HERE (these belong in _adr85_t3_design.md).

A1. "tau_ov relative to min(L_s,L_m)/L_m" is read as
        threshold = tau_ov * min(L_s,L_m)/L_m ,  tau_ov a dimensionless factor.
    Justification, not taste: min(L_s,L_m)/L_m IS the largest parametric overlap the pair
    can attain (a contained slave gives L_s|cos|/L_m; a covered master gives 1), so it is
    the exact 2D analogue of the 3D A_min normalizer at LadrunoMortarKernel.h:444.

A2. STATUS vs REFUSAL CLASSIFICATION -- reconciled with binding decision 2b.1.
    The ADR names only "status = -1 (empty)", and the orchestrator's binding decision
    ("status codes mirror the 3D kernel: 0, -1 empty/sliver/perpendicular skip-with-
    routing, no -2 case") keeps it that way.  But decision 2b.4 needs the handler to tell a
    NEAR-PERPENDICULAR refusal apart from an ordinary empty clip, to fire
    WARN_2D_PERP_NO_NTS on exactly one of them -- and the PairResult2D field list in 2b.1
    carries nothing that would let it, so the handler would have to RE-DERIVE cos_t from
    the coordinates: the fork's staleness trap, on the very quantity the gate turns on.
    Resolution adopted here, which costs no status code and no re-derivation:
        status  in {0, -1}                      <- exactly as 2b.1 mandates
        refusal in {NONE, EMPTY, DEGEN, ALIGN}  <- a NEW PairResult2D field, set by the
                                                   kernel that already computed cos_t
    "skip-with-routing" in 2b.1 then has something to route ON.  RECOMMENDATION to the
    T3 implementer: add `refusal` (and `cos_t`) to the C++ PairResult2D; do NOT recover the
    perpendicular case in the handler.  Zero-length keeps status -1 with refusal = DEGEN
    (2b.1: "zero-length is refused by tau_seg as in NTS"), so there is no -2 tier.

A3. The ADR's normalizer omits |cos|; the true attainable maximum is
    min(L_s|cos|, L_m)/L_m.  Keeping the literal (|cos|-free) form makes the guard up to
    1/|cos| stricter.  That is harmless ONLY because tau_align > tau_ov, which is now a
    GATED requirement (CHECK 1), not an assumption: for a contained slave the sliver guard
    bites exactly when |cos| < tau_ov, so the alignment refusal must fire first.

A4. Defaults, measured not guessed: tau_ov = 1e-9 (the far-field TRANSLATION noise floor
    is 9.0e-12 at |x| ~ 5e4, so 1e-10 would leave ~11x; 1e-9 buys ~100x and refuses only
    overlaps below 1e-9 of the attainable one), tau_align = 1e-6 (caps the slave-trace
    Jacobian 1/|cos| at 1e6).

A5. The tie patch is exact for fields whose gradient is TANGENTIAL to the master.  A
    normal-gradient component differs between a slave point and its projection by exactly
    g*n -- physics, not a defect -- so CHECK 5 uses u = c0 + c1*x on a horizontal master.

A6. "lambda_N is NOT separately scaled" (SS How/7) means lambda INHERITS h, because it
    accumulates from epsN*h*gbar.  So the correct -thickness twin's committed multiplier is
    h*lambda_1, and "tractions scale by exactly h" then holds BITWISE for representable h.
    Scaling lambda again in the residual is the h^2 error; CHECK 8 gates it.

A7. SS How/3 does not restate the ALM normalization, so this oracle adopts the shipped 3D
    convention verbatim (proto_c2_alm.py): gbar = gtil/w with w = D row-sum,
    t = min(0, lambda + epsN*gbar) with COMPRESSION NEGATIVE, F^s = D t, F^m = -M^T t,
    and lambda <- min(0, lambda + epsN*gbar) once per commit.

A8. A noise-sized ACCEPTED interval is benign in the kernel (D, M, gtil all scale with
    b-a, so the contribution is noise-sized too) but NOT in the enforcement layer, which
    forms gbar = gtil/w and would divide by a noise-sized tributary.  Flagged for the C++
    handler: the active-set test needs its own tributary floor.  Kernel-level scope here.
------------------------------------------------------------------------------

Style rules honoured (fork discipline): RELATIVE tolerances only (scaled by L_seg, by
|x|_scale, or by the integral magnitude -- never a bare 1e-12 absolute); every random
family seeded; refusals assert the STATUS CODE, never an exception; the "exact" claims are
adjudicated by fractions.Fraction (the proto_t2_reuse_vs_scalar.py precedent), not by a
finer float rule.

Run:  python3.12 proto_t3_mortar2d.py        (raw numpy + fractions; no compiler needed)
"""
import math
import sys
from fractions import Fraction as F

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

EPS = np.finfo(np.float64).eps

# --- kernel status codes (mirror the future Ladruno2DMortarKernel.h) ---------------
# The 3D kernel uses 0 / -1 (empty, skip) / -2 (invalid geometry, caller must hard-error).
# 2D adds -3: the ALIGNMENT refusal, which the ADR-85 SS How/3 ownership protocol must be
# able to tell apart from an ordinary empty clip (it is the case that fires
# WARN_2D_PERP_NO_NTS when no NTS lane is declared on the same surfaces).
STATUS_OK = 0         # accumulated
STATUS_SKIP = -1      # ANY refusal -- skip-with-routing (no -2 tier: straight 2-node
#                       segments cannot be "invalid geometry" the way quadratic facets can)
# The refusal CLASSIFICATION is a separate field, not a separate status code.  See A2.
REFUSE_NONE = 0
REFUSE_EMPTY = 1      # empty / sliver overlap
REFUSE_DEGEN = 2      # zero-length segment (tau_seg)
REFUSE_ALIGN = 3      # |t_s . t_m| below tau_align -> WARN_2D_PERP_NO_NTS / NTS ownership
_RNAME = {REFUSE_NONE: "NONE", REFUSE_EMPTY: "EMPTY", REFUSE_DEGEN: "DEGEN",
          REFUSE_ALIGN: "ALIGN"}


def refused(r, code):
    """A refusal assertion is on the STATUS CODE plus the classification -- never on an
    exception, and never on a bare boolean."""
    return r.status == STATUS_SKIP and r.refusal == code

# --- thresholds (the defaults this oracle VALIDATES) ------------------------------
# TAU_OV is set by the far-field TRANSLATION noise floor measured in CHECK 1.7: the same
# physical pair re-represented at |x| ~ 5e4 moves [a,b] by ~9e-12, so 1e-10 would leave only
# a ~11x margin.  1e-9 buys ~100x and costs nothing measurable -- it refuses only overlaps
# below 1e-9 of the attainable one, i.e. below any force that matters.
TAU_OV = 1e-9         # parametric overlap sliver, RELATIVE to min(L_s,L_m)/L_m
TAU_SEG = 1e-8        # zero-length segment refusal, RELATIVE to Lref (ADR-85 SS How/1)
# TAU_ALIGN caps the slave-trace Jacobian 1/|cos| at 1e6 -- already far past any well-posed
# pair.  The ORDERING TAU_ALIGN > TAU_OV is a hard requirement, not a preference: for a
# fully contained slave, b-a = L_s|cos|/L_m while the sliver threshold is tau_ov*L_s/L_m,
# so b-a > thr  <=>  |cos| > tau_ov.  If tau_ov >= tau_align the sliver guard would swallow
# legitimate strongly-tilted FULL overlaps before the alignment refusal could route them.
TAU_ALIGN = 1e-6      # |t_s . t_m| alignment floor -> REFUSE_ALIGN (gated in CHECK 1)

# =================================================================== reporting
_fails = 0
_sec = {}          # section -> [checks, fails]
_secerr = {}       # section -> max relative error seen


def check(sec, name, ok, extra=""):
    global _fails
    st = _sec.setdefault(sec, [0, 0])
    st[0] += 1
    if not ok:
        st[1] += 1
        _fails += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] C{sec}: {name}{(' -- ' + extra) if extra else ''}")


def rel(sec, x, ref, scale=None):
    """Relative error of x vs ref against an explicit SCALE (never |ref| when ref may
    vanish).  Records the section maximum for the final report."""
    x = np.asarray(x, float)
    ref = np.asarray(ref, float)
    if scale is None:
        scale = max(float(np.max(np.abs(ref))), 0.0)
    scale = float(scale)
    if scale == 0.0:
        e = float(np.max(np.abs(x - ref)))
    else:
        e = float(np.max(np.abs(x - ref))) / scale
    _secerr[sec] = max(_secerr.get(sec, 0.0), e)
    return e


# ================================================================== the 2D kernel
# Scalar operations, written out, NOT np.dot: numpy is free to reassociate a reduction and
# a reassociated dot product is not bit-reproducible against the C++ transcription.  The
# operation order below is the contract (proto_t1_nts2d.py precedent).
def _dot2(ax, ay, bx, by):
    return ax * bx + ay * by


def gauss01(n):
    """n-point Gauss-Legendre mapped to [0,1]:  INT_0^1 f = sum w_i f(u_i)."""
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


class PairResult2D:
    """The 2D analogue of LadrunoMortarKernel::PairResult (2 slave x 2 master nodes)."""

    __slots__ = ("status", "refusal", "ngp", "a", "b", "len", "gapL1", "J", "cos_t",
                 "Ls", "Lm", "al0", "al1", "n", "D", "M", "g")

    def __init__(self):
        self.status = STATUS_SKIP
        self.refusal = REFUSE_EMPTY
        self.ngp = 0
        self.a = self.b = 0.0
        self.len = 0.0
        self.gapL1 = 0.0
        self.J = 0.0
        self.cos_t = 0.0
        self.Ls = self.Lm = 0.0
        self.al0 = self.al1 = 0.0
        self.n = np.zeros(2)
        self.D = np.zeros((2, 2))
        self.M = np.zeros((2, 2))
        self.g = np.zeros(2)


def integrate_pair_2d(Xs, Xm, sigma, order=2, Lref=None,
                      tau_ov=TAU_OV, tau_seg=TAU_SEG, tau_align=TAU_ALIGN,
                      jac_mode="slave_trace"):
    """ONE slave/master 2-node segment pair.  ADR-85 SS How/3.

    Xs = [[x,y],[x,y]] slave segment, Xm likewise master, sigma = +-1 (orientation vote).
    jac_mode:
      "slave_trace" -- the CORRECT measure, J = L_m/|t_s.t_m|          (dGamma_s)
      "no_align"    -- the DELIBERATELY BROKEN twin, J = L_m           (dGamma_m; the bug
                       CHECK 3 exists to catch; identical when cos = 1)
    """
    out = PairResult2D()
    s0x, s0y = float(Xs[0][0]), float(Xs[0][1])
    s1x, s1y = float(Xs[1][0]), float(Xs[1][1])
    m0x, m0y = float(Xm[0][0]), float(Xm[0][1])
    m1x, m1y = float(Xm[1][0]), float(Xm[1][1])

    tsx, tsy = s1x - s0x, s1y - s0y
    tmx, tmy = m1x - m0x, m1y - m0y
    Ls = math.hypot(tsx, tsy)
    Lm = math.hypot(tmx, tmy)
    out.Ls, out.Lm = Ls, Lm
    if Lref is None:
        Lref = Ls if Ls > Lm else Lm
    if Ls <= tau_seg * Lref or Lm <= tau_seg * Lref:
        out.status, out.refusal = STATUS_SKIP, REFUSE_DEGEN
        return out

    # alignment: cos(theta) between the unit tangents.  1/|cos| is the load-bearing factor.
    cos_t = _dot2(tsx / Ls, tsy / Ls, tmx / Lm, tmy / Lm)
    out.cos_t = cos_t
    if abs(cos_t) < tau_align:
        out.status, out.refusal = STATUS_SKIP, REFUSE_ALIGN          # near-perpendicular: ill-posed measure
        return out

    # closed-form projection of the slave ENDPOINTS onto the master parametrization
    Lm2 = _dot2(tmx, tmy, tmx, tmy)
    al0 = _dot2(s0x - m0x, s0y - m0y, tmx, tmy) / Lm2
    al1 = _dot2(s1x - m0x, s1y - m0y, tmx, tmy) / Lm2
    out.al0, out.al1 = al0, al1
    lo = al0 if al0 < al1 else al1
    hi = al1 if al0 < al1 else al0
    a = lo if lo > 0.0 else 0.0
    b = hi if hi < 1.0 else 1.0
    out.a, out.b = a, b

    # sliver guard -- RELATIVE to min(L_s,L_m)/L_m, which is exactly the largest parametric
    # overlap the pair can attain (the 2D analogue of the 3D A_min normalizer).
    thr = tau_ov * (Ls if Ls < Lm else Lm) / Lm
    if b - a <= thr:
        out.status, out.refusal = STATUS_SKIP, REFUSE_EMPTY
        return out

    nx = -sigma * tmy / Lm
    ny = sigma * tmx / Lm
    out.n = np.array([nx, ny])

    if jac_mode == "slave_trace":
        J = Lm / abs(cos_t)                 # dGamma_s = (L_m/|t_s.t_m|) dxi_m
    elif jac_mode == "no_align":
        J = Lm                              # the BUG: master measure, alignment omitted
    else:
        raise ValueError(jac_mode)
    out.J = J

    dal = al1 - al0                          # nonzero: |cos| >= tau_align
    u, w = gauss01(order)
    span = b - a
    for gp in range(order):
        xim = a + span * u[gp]
        wgt = w[gp] * span
        xis = (xim - al0) / dal
        Ns = (1.0 - xis, xis)
        Nm = (1.0 - xim, xim)
        xsx = s0x + xis * tsx
        xsy = s0y + xis * tsy
        xmx = m0x + xim * tmx
        xmy = m0y + xim * tmy
        gN = _dot2(xsx - xmx, xsy - xmy, nx, ny)
        wJ = wgt * J
        out.len += wJ
        out.gapL1 += wJ * abs(gN)
        for i in range(2):
            out.g[i] += wJ * Ns[i] * gN
            for j in range(2):
                out.D[i][j] += wJ * Ns[i] * Ns[j]
            for k in range(2):
                out.M[i][k] += wJ * Ns[i] * Nm[k]
        out.ngp += 1
    out.status, out.refusal = STATUS_OK, REFUSE_NONE
    return out


# ==================================================== EXACT rational referee (Fraction)
# Polynomials in xi_m as Fraction coefficient lists, lowest order first.
def pmul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            r[i + j] += pi * qj
    return r


def pint(p, a, b):
    s = F(0)
    for i, ci in enumerate(p):
        s += ci * (b ** (i + 1) - a ** (i + 1)) / (i + 1)
    return s


def exact_reduced(Xs, Xm, a, b):
    """EXACT (rational) values of the DIMENSIONLESS integrals over [a,b] in xi_m:

        Dred_ij = INT N_i^s N_j^s dxi_m      Mred_ik = INT N_i^s N_k^m dxi_m
        ghat_i  = INT N_i^s ghat  dxi_m ,    ghat = (x_s - m0).perp(t_m)   (= g * L_m/sigma)

    Every input float is converted EXACTLY by Fraction(float), so this adjudicates the
    QUADRATURE alone -- the irrational constants (J, L_m) are supplied by the kernel and
    multiplied back in float.  The clip [a,b] is taken from the kernel too, so CHECK 1
    (projection) and CHECK 2 (quadrature) are adjudicated separately, never entangled.
    """
    s0 = (F(Xs[0][0]), F(Xs[0][1]))
    s1 = (F(Xs[1][0]), F(Xs[1][1]))
    m0 = (F(Xm[0][0]), F(Xm[0][1]))
    m1 = (F(Xm[1][0]), F(Xm[1][1]))
    ts = (s1[0] - s0[0], s1[1] - s0[1])
    tm = (m1[0] - m0[0], m1[1] - m0[1])
    Lm2 = tm[0] * tm[0] + tm[1] * tm[1]
    al0 = ((s0[0] - m0[0]) * tm[0] + (s0[1] - m0[1]) * tm[1]) / Lm2
    al1 = ((s1[0] - m0[0]) * tm[0] + (s1[1] - m0[1]) * tm[1]) / Lm2
    dal = al1 - al0
    # xi_s(xi_m) = (xi_m - al0)/dal
    xis = [-al0 / dal, F(1) / dal]
    Ns = [[F(1) - xis[0], -xis[1]], [xis[0], xis[1]]]
    Nm = [[F(1), F(-1)], [F(0), F(1)]]
    px, py = -tm[1], tm[0]
    tsp = ts[0] * px + ts[1] * py
    d0 = (s0[0] - m0[0]) * px + (s0[1] - m0[1]) * py
    ghat = [d0 + xis[0] * tsp, xis[1] * tsp]
    A, B = F(a), F(b)
    Dred = [[pint(pmul(Ns[i], Ns[j]), A, B) for j in range(2)] for i in range(2)]
    Mred = [[pint(pmul(Ns[i], Nm[k]), A, B) for k in range(2)] for i in range(2)]
    gred = [pint(pmul(Ns[i], ghat), A, B) for i in range(2)]
    return Dred, Mred, gred, (al0, al1)


def exact_reference(res, Xs, Xm, sigma):
    """Float reference D/M/gtil built from the EXACT reduced integrals."""
    Dred, Mred, gred, _ = exact_reduced(Xs, Xm, res.a, res.b)
    D = np.array([[res.J * float(Dred[i][j]) for j in range(2)] for i in range(2)])
    M = np.array([[res.J * float(Mred[i][k]) for k in range(2)] for i in range(2)])
    g = np.array([res.J * sigma / res.Lm * float(gred[i]) for i in range(2)])
    return D, M, g


# ============================================== closed-form full-coverage reference
def analytic_full_cover(Xs, Xm, sigma):
    """Closed-form slave-trace integrals for a slave segment ENTIRELY covered by the
    master's parametric range.  Derived directly from dGamma_s = L_s dxi_s:

        area  = L_s
        D     = L_s * [[1/3, 1/6], [1/6, 1/3]]
        gtil  = L_s * [ga/3 + gb/6 , ga/6 + gb/3]        (g linear along the slave)
        M     = L_s * [[1/2-a0/2-d/6, a0/2+d/6],
                       [1/2-a0/2-d/3, a0/2+d/3]]         a0 = xiM(s0), d = xiM(s1)-xiM(s0)
    Row sums of M equal row sums of D (= L_s/2): partition of unity, independent of tilt.
    """
    s0 = np.asarray(Xs[0], float)
    s1 = np.asarray(Xs[1], float)
    m0 = np.asarray(Xm[0], float)
    m1 = np.asarray(Xm[1], float)
    ts, tm = s1 - s0, m1 - m0
    Ls, Lm = float(np.hypot(*ts)), float(np.hypot(*tm))
    n = sigma * np.array([-tm[1], tm[0]]) / Lm
    a0 = float(np.dot(s0 - m0, tm)) / (Lm * Lm)
    a1 = float(np.dot(s1 - m0, tm)) / (Lm * Lm)
    d = a1 - a0
    ga = float(np.dot(s0 - m0, n))
    gb = float(np.dot(s1 - m0, n))
    D = Ls * np.array([[1.0 / 3, 1.0 / 6], [1.0 / 6, 1.0 / 3]])
    M = Ls * np.array([[0.5 - a0 / 2 - d / 6, a0 / 2 + d / 6],
                       [0.5 - a0 / 2 - d / 3, a0 / 2 + d / 3]])
    g = Ls * np.array([ga / 3 + gb / 6, ga / 6 + gb / 3])
    return Ls, D, M, g


# ================================================================ interface assembly
def assemble_interface(sN, sC, mN, mC, sigma, order=2, jac_mode="slave_trace",
                       tau_ov=TAU_OV):
    """Loop every (slave segment, master segment) pair; accumulate the global mortar
    operators.  Returns D (ns x ns), M (ns x nm), gtil (ns), total area, and the per-pair
    status table (the refusal bookkeeping CHECK 5 gates)."""
    sN = np.asarray(sN, float)
    mN = np.asarray(mN, float)
    ns, nm = len(sN), len(mN)
    D = np.zeros((ns, ns))
    M = np.zeros((ns, nm))
    g = np.zeros(ns)
    area = 0.0
    Lref = max(float(np.hypot(*(sN[e[1]] - sN[e[0]]))) for e in sC)
    table = []
    for es, (i0, i1) in enumerate(sC):
        for em, (k0, k1) in enumerate(mC):
            r = integrate_pair_2d((sN[i0], sN[i1]), (mN[k0], mN[k1]), sigma,
                                  order=order, Lref=Lref, tau_ov=tau_ov,
                                  jac_mode=jac_mode)
            table.append((es, em, (r.status, r.refusal), r.a, r.b))
            if r.status != STATUS_OK:
                continue
            area += r.len
            loc_s, loc_m = (i0, i1), (k0, k1)
            for i in range(2):
                g[loc_s[i]] += r.g[i]
                for j in range(2):
                    D[loc_s[i], loc_s[j]] += r.D[i][j]
                for k in range(2):
                    M[loc_s[i], loc_m[k]] += r.M[i][k]
    return D, M, g, area, table


def chain(x0, xs, y):
    """Node list along a straight horizontal line + a head-to-tail 2-node connectivity."""
    N = np.array([[x0 + t, y] for t in xs], float)
    C = [(i, i + 1) for i in range(len(xs) - 1)]
    return N, C


def chain_dir(origin, direction, arclen, y_off=0.0):
    """Node list along an arbitrary straight line, arc-length parametrized."""
    o = np.asarray(origin, float)
    d = np.asarray(direction, float)
    d = d / np.hypot(*d)
    nrm = np.array([-d[1], d[0]])
    N = np.array([o + s * d + y_off * nrm for s in arclen], float)
    C = [(i, i + 1) for i in range(len(arclen) - 1)]
    return N, C


print("=" * 78)
print("ADR-85 T3 ORACLE -- 2D MORTAR (interval clip on the slave-trace measure)")
print(f"  tau_ov = {TAU_OV:.0e} (RELATIVE to min(L_s,L_m)/L_m)   "
      f"tau_align = {TAU_ALIGN:.0e}   tau_seg = {TAU_SEG:.0e}")
print("=" * 78)

# ============================================================================
# CHECK 1 -- INTERVAL CLIP CORRECTNESS
# ============================================================================
print("\nCHECK 1  interval clip: exact-rational parity, brute-force set parity, the case "
      "table,\n         parametrization invariance, and the tau_ov RELATIVE guard.")

rng = np.random.default_rng(20260818)

# 1.1 -- exact rational parity of the clipped interval over randomized pairs
n11 = 4000
worst_alpha, worst_ulp = 0.0, 0.0
nacc = 0
for _ in range(n11):
    m0 = rng.uniform(-3, 3, 2)
    ang = rng.uniform(0, 2 * math.pi)
    Lm = rng.uniform(0.5, 4.0)
    m1 = m0 + Lm * np.array([math.cos(ang), math.sin(ang)])
    s0 = m0 + rng.uniform(-1.5, 1.5, 2)
    dang = rng.uniform(-1.2, 1.2)
    Ls = rng.uniform(0.2, 5.0)
    s1 = s0 + Ls * np.array([math.cos(ang + dang), math.sin(ang + dang)])
    r = integrate_pair_2d((s0, s1), (m0, m1), 1.0)
    _, _, _, (ea0, ea1) = exact_reduced((s0, s1), (m0, m1), 0.0, 1.0)
    elo, ehi = (ea0, ea1) if ea0 < ea1 else (ea1, ea0)
    ea = max(elo, F(0)); eb = min(ehi, F(1))
    # parametric conditioning: xi = (x - m0).t_m/L_m^2 loses |x - m0|/L_m relative digits
    cond = 1.0 + float(np.hypot(*(s0 - m0)) + np.hypot(*(s1 - m0))) / float(np.hypot(*(m1 - m0)))
    if r.status == STATUS_OK:
        nacc += 1
        da = abs(r.a - float(ea)) + abs(r.b - float(eb))
        worst_alpha = max(worst_alpha, da)
        worst_ulp = max(worst_ulp, da / (EPS * cond))
    else:
        # a refusal must correspond to an exactly-empty / sliver exact interval
        thr = TAU_OV * min(r.Ls, r.Lm) / r.Lm if r.Lm > 0 else 0.0
        if refused(r, REFUSE_EMPTY) and float(eb - ea) > 10.0 * thr:
            worst_alpha = float("inf")
_secerr[1] = max(_secerr.get(1, 0.0), worst_alpha)
check(1, f"exact-rational clip parity over {n11} random pairs ({nacc} accepted)",
      worst_alpha < 1e-13, f"max |dxi| = {worst_alpha:.3e}  ({worst_ulp:.1f} x eps*cond)")

# 1.2 -- brute-force geometric reference: the accepted interval IS the covered set
n12, NS = 800, 20001
grid = np.linspace(0.0, 1.0, NS)
bad12 = 0
for _ in range(n12):
    m0 = rng.uniform(-2, 2, 2)
    ang = rng.uniform(0, 2 * math.pi)
    Lm = rng.uniform(0.7, 3.0)
    tm = Lm * np.array([math.cos(ang), math.sin(ang)])
    m1 = m0 + tm
    s0 = m0 + rng.uniform(-2.0, 2.0, 2)
    dang = rng.uniform(-1.0, 1.0)
    Ls = rng.uniform(0.3, 3.0)
    ts = Ls * np.array([math.cos(ang + dang), math.sin(ang + dang)])
    s1 = s0 + ts
    r = integrate_pair_2d((s0, s1), (m0, m1), 1.0)
    a0 = float(np.dot(s0 - m0, tm)) / (Lm * Lm)
    a1 = float(np.dot(s1 - m0, tm)) / (Lm * Lm)
    xis = (grid - a0) / (a1 - a0)
    cov = (xis >= 0.0) & (xis <= 1.0)
    h = 2.0 / (NS - 1)
    if not cov.any():
        # nothing covered -> either refused, or an accepted interval below grid resolution
        if not (r.status != STATUS_OK or (r.b - r.a) < h):
            bad12 += 1
        continue
    glo, ghi = grid[cov][0], grid[cov][-1]
    if r.status != STATUS_OK:
        if (ghi - glo) > h:          # a real interval was refused: a bug
            bad12 += 1
        continue
    if abs(r.a - glo) > h or abs(r.b - ghi) > h:
        bad12 += 1
check(1, f"brute-force covered-set parity over {n12} pairs ({NS}-station scan)",
      bad12 == 0, f"{bad12} disagreements (grid resolution {2.0 / (NS - 1):.2e})")

# 1.3 -- the deterministic case table (overlap / partial / disjoint / containment /
#        reversed parametrization / endpoint-touch sliver)
M0, M1 = (0.0, 0.0), (10.0, 0.0)
cases = [
    ("full overlap        ", (2.0, 0.0), (8.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.2, 0.8),
    ("partial low         ", (-3.0, 0.0), (4.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.0, 0.4),
    ("partial high        ", (6.0, 0.0), (14.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.6, 1.0),
    ("disjoint left       ", (-9.0, 0.0), (-1.0, 0.0), M0, M1, (STATUS_SKIP, REFUSE_EMPTY), 0.0, 0.0),
    ("disjoint right      ", (11.0, 0.0), (19.0, 0.0), M0, M1, (STATUS_SKIP, REFUSE_EMPTY), 0.0, 0.0),
    ("slave inside master ", (3.0, 0.0), (5.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.3, 0.5),
    ("master inside slave ", (-5.0, 0.0), (15.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.0, 1.0),
    ("endpoint touch      ", (10.0, 0.0), (12.0, 0.0), M0, M1, (STATUS_SKIP, REFUSE_EMPTY), 0.0, 0.0),
    ("reversed slave param", (8.0, 0.0), (2.0, 0.0), M0, M1, (STATUS_OK, REFUSE_NONE), 0.2, 0.8),
    ("reversed mastr param", (2.0, 0.0), (8.0, 0.0), M1, M0, (STATUS_OK, REFUSE_NONE), 0.2, 0.8),
]
bad13, e13 = 0, 0.0
for nm_, s0, s1, m0, m1, st, ea, eb in cases:
    r = integrate_pair_2d((s0, s1), (m0, m1), 1.0)
    ok = ((r.status, r.refusal) == st)
    if st[0] == STATUS_OK:
        e = max(abs(r.a - ea), abs(r.b - eb))
        e13 = max(e13, e)
        ok = ok and e < 1e-15
    if not ok:
        bad13 += 1
        print(f"        case '{nm_.strip()}': status {r.status}/"
              f"{_RNAME[r.refusal]} [a,b]=[{r.a},{r.b}]")
_secerr[1] = max(_secerr.get(1, 0.0), e13)
check(1, f"deterministic case table ({len(cases)} cases: overlap/partial/disjoint/"
         "containment/reversed/touch)", bad13 == 0, f"max |dxi| = {e13:.3e}")

# 1.4 -- parametrization invariance.  Reversing the MASTER listing flips perp(t_m), so the
#        physical normal is preserved only if sigma flips too (handler business); reversing
#        the SLAVE listing permutes the slave rows/cols.
Xs = [(1.5, 0.35), (7.25, 0.35)]
Xm = [(0.0, 0.0), (10.0, 0.0)]
rA = integrate_pair_2d(Xs, Xm, +1.0)
rMrev = integrate_pair_2d(Xs, [Xm[1], Xm[0]], -1.0)
P = np.array([[0.0, 1.0], [1.0, 0.0]])
e_m = max(rel(1, rMrev.D, rA.D, rA.len),
          rel(1, rMrev.M @ P, rA.M, rA.len),
          rel(1, rMrev.g, rA.g, rA.len * abs(Xs[0][1])))
check(1, "master-listing reversal (with sigma flip) leaves D, gtil and column-permuted M "
         "invariant", e_m < 64 * EPS, f"max rel = {e_m:.3e}")
rSrev = integrate_pair_2d([Xs[1], Xs[0]], Xm, +1.0)
e_s = max(rel(1, P @ rSrev.D @ P, rA.D, rA.len),
          rel(1, P @ rSrev.M, rA.M, rA.len),
          rel(1, P @ rSrev.g, rA.g, rA.len * abs(Xs[0][1])))
check(1, "slave-listing reversal permutes slave rows/cols and nothing else",
      e_s < 64 * EPS, f"max rel = {e_s:.3e}")

# 1.5 -- tau_ov is RELATIVE: the fine-slave / coarse-master panel catch
print("        tau_ov sweep (L_s/L_m fully-overlapping fine slave vs a NAIVE ABSOLUTE guard):")
TAU_ABS_NAIVE = 1e-4          # a plausible-looking ABSOLUTE parametric sliver guard
ratios = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
rel_ok, abs_wrong = 0, 0
fine_area_err = 0.0
for q in ratios:
    Ls_ = 10.0 * q
    s0 = (4.0, 0.0)
    s1 = (4.0 + Ls_, 0.0)
    r = integrate_pair_2d((s0, s1), (M0, M1), 1.0)
    acc_rel = (r.status == STATUS_OK)
    acc_abs = (r.b - r.a) > TAU_ABS_NAIVE
    rel_ok += acc_rel
    abs_wrong += (acc_rel and not acc_abs)
    # NOT folded into the section max: this is the DOCUMENTED cost of a fine slave, not a
    # defect -- b-a is an O(eps) difference of O(1) parametric coordinates, so its own
    # relative accuracy degrades like eps/(L_s/L_m).  Reported, not gated.
    ar = abs((r.len if acc_rel else 0.0) - Ls_) / Ls_
    fine_area_err = max(fine_area_err, ar)
    print(f"          L_s/L_m={q:.0e}  b-a={r.b - r.a:.3e}  relative:{'ACCEPT' if acc_rel else 'refuse'}"
          f"  naive-abs({TAU_ABS_NAIVE:.0e}):{'accept' if acc_abs else 'REFUSE'}"
          f"   area rel err {ar:.2e}")
check(1, f"RELATIVE tau_ov accepts every legitimate fine-slave pair ({len(ratios)}/"
         f"{len(ratios)}); the naive ABSOLUTE guard wrongly refuses {abs_wrong} of them",
      rel_ok == len(ratios) and abs_wrong >= 2,
      f"relative accepted {rel_ok}/{len(ratios)}, absolute lost {abs_wrong}")

# a GENUINE sliver (parametric overlap far below the relative threshold) is refused
q_sl = TAU_OV * 1e-2                      # b-a = 1e-12 with min(L_s,L_m)/L_m = 1
r_sl = integrate_pair_2d(((10.0 - 10.0 * q_sl, 0.0), (20.0, 0.0)), (M0, M1), 1.0)
check(1, "a genuine sliver (b-a = tau_ov/100 of the attainable overlap) is REFUSED with "
         "REFUSE_EMPTY", refused(r_sl, REFUSE_EMPTY), f"status {r_sl.status}/{_RNAME[r_sl.refusal]}")

# 1.6 -- the other two refusals, by STATUS CODE (never an exception)
r_dg = integrate_pair_2d(((3.0, 0.0), (3.0 + 1e-12, 0.0)), (M0, M1), 1.0, Lref=10.0)
check(1, "zero-length slave segment -> status -1 / refusal DEGEN, not an exception",
      refused(r_dg, REFUSE_DEGEN), f"status {r_dg.status}/{_RNAME[r_dg.refusal]}")
r_pp = integrate_pair_2d(((5.0, -1.0), (5.0, 1.0)), (M0, M1), 1.0)
check(1, "exactly perpendicular slave -> status -1 / refusal ALIGN, the "
         "classification the ADR-85 SS How/3 ownership protocol and WARN_2D_PERP_NO_NTS "
         "key on (see A2)",
      refused(r_pp, REFUSE_ALIGN), f"status {r_pp.status}/{_RNAME[r_pp.refusal]}")

# tau_ov MUST stay below tau_align: for a contained slave b-a = L_s|cos|/L_m against a
# threshold tau_ov*L_s/L_m, so the sliver guard bites exactly when |cos| < tau_ov.  If the
# ordering were reversed it would swallow legitimate strongly-tilted FULL overlaps before
# the alignment refusal could route them to the NTS lane.
cbar = 3.0 * TAU_ALIGN
th = math.acos(cbar)
r_tt = integrate_pair_2d(((0.0, 0.0), (10.0 * math.cos(th), 10.0 * math.sin(th))),
                         (M0, M1), 1.0)
r_ov = integrate_pair_2d(((0.0, 0.0), (10.0 * TAU_OV * 0.5,
                                       10.0 * math.sqrt(1.0 - (TAU_OV * 0.5) ** 2))),
                         (M0, M1), 1.0)
check(1, f"threshold ORDERING tau_ov ({TAU_OV:.0e}) << tau_align ({TAU_ALIGN:.0e}) holds: a "
         f"full overlap at |cos| = {cbar:.0e} (just above tau_align) is ACCEPTED, while a "
         f"|cos| = {TAU_OV * 0.5:.0e} pair is refused by the ALIGNMENT gate, never silently "
         "by the sliver gate",
      TAU_OV < 1e-2 * TAU_ALIGN and r_tt.status == STATUS_OK
      and refused(r_ov, REFUSE_ALIGN),
      f"b-a = {r_tt.b - r_tt.a:.3e} status {r_tt.status}; deep-tilt "
      f"{r_ov.status}/{_RNAME[r_ov.refusal]}")

# 1.7 -- measured far-from-origin parametric noise floor vs tau_ov.  The master is TILTED
# and of non-unit length on purpose: an axis-aligned unit master makes the projection
# arithmetic exact and would measure nothing.
X0 = 5.0e4
noise = 0.0
for _ in range(3000):
    m0 = np.array([X0, -3.0e4]) + rng.uniform(-1.0, 1.0, 2)
    ang = rng.uniform(0, 2 * math.pi)
    Lm = rng.uniform(0.7, 3.0)
    tm = Lm * np.array([math.cos(ang), math.sin(ang)])
    m1 = m0 + tm
    nrm = np.array([-tm[1], tm[0]]) / Lm
    t = rng.uniform(0.05, 0.95)
    s0 = m0 + t * tm + 1e-3 * nrm
    dang = rng.uniform(-0.8, 0.8)
    s1 = s0 + 0.2 * np.array([math.cos(ang + dang), math.sin(ang + dang)])
    r = integrate_pair_2d((s0, s1), (m0, m1), 1.0)
    _, _, _, (ea0, ea1) = exact_reduced((s0, s1), (m0, m1), 0.0, 1.0)
    noise = max(noise, abs(r.al0 - float(ea0)), abs(r.al1 - float(ea1)))
_secerr[1] = max(_secerr.get(1, 0.0), noise)
check(1, f"far-from-origin (|x| ~ {X0:.0e}) ALGORITHMIC noise: the float clip agrees with "
         "exact rational arithmetic on the same inputs to ~1 ulp, because the coordinate "
         "differences are EXACT (Sterbenz) however far from the origin",
      noise < 64 * EPS, f"measured max |dxi| = {noise:.3e}  ({noise / EPS:.1f} ulp)")

# The quantity the sliver guard actually has to clear is the TRANSLATION noise: the same
# physical pair, its coordinates re-represented at |x| ~ 5e4, moves [a,b] by ~eps*|x|/L_m.
# THAT is the number tau_ov must sit above.
tnoise = 0.0
for _ in range(3000):
    ang = rng.uniform(0, 2 * math.pi)
    Lm = rng.uniform(0.7, 3.0)
    tm = Lm * np.array([math.cos(ang), math.sin(ang)])
    m0 = rng.uniform(-1.0, 1.0, 2)
    m1 = m0 + tm
    nrm = np.array([-tm[1], tm[0]]) / Lm
    s0 = m0 + rng.uniform(0.05, 0.6) * tm + 1e-3 * nrm
    dang = rng.uniform(-0.8, 0.8)
    s1 = s0 + 0.2 * np.array([math.cos(ang + dang), math.sin(ang + dang)])
    T = np.array([X0, -3.0e4])
    rA_ = integrate_pair_2d((s0, s1), (m0, m1), 1.0)
    rB_ = integrate_pair_2d((s0 + T, s1 + T), (m0 + T, m1 + T), 1.0)
    if rA_.status == STATUS_OK and rB_.status == STATUS_OK:
        tnoise = max(tnoise, abs(rA_.a - rB_.a), abs(rA_.b - rB_.b))
check(1, f"far-from-origin TRANSLATION noise floor stays below tau_ov, so the sliver guard "
         f"is a sliver guard and not a noise gate at |x| ~ {X0:.0e}",
      tnoise < TAU_OV, f"measured max |dxi| = {tnoise:.3e}  vs tau_ov = {TAU_OV:.0e} "
                       f"(margin {TAU_OV / max(tnoise, 1e-300):.1f}x)")

# ============================================================================
# CHECK 2 -- 2-POINT GAUSS IS EXACT
# ============================================================================
print("\nCHECK 2  2-pt Gauss exactness: 2pt vs 50pt vs EXACT rational quadrature.")

rng2 = np.random.default_rng(4110085)
n2 = 3000
e_2v50 = {"D": 0.0, "M": 0.0, "g": 0.0, "len": 0.0}
e_2vex = {"D": 0.0, "M": 0.0, "g": 0.0}
e_50vex = {"D": 0.0, "M": 0.0, "g": 0.0}
ntried = 0
while ntried < n2:
    m0 = rng2.uniform(-4, 4, 2)
    ang = rng2.uniform(0, 2 * math.pi)
    Lm = rng2.uniform(0.5, 4.0)
    m1 = m0 + Lm * np.array([math.cos(ang), math.sin(ang)])
    dang = rng2.uniform(-1.2, 1.2)
    Ls = rng2.uniform(0.2, 6.0)
    s0 = m0 + rng2.uniform(-2.0, 2.0, 2)
    s1 = s0 + Ls * np.array([math.cos(ang + dang), math.sin(ang + dang)])
    sg = 1.0 if rng2.random() < 0.5 else -1.0
    r2 = integrate_pair_2d((s0, s1), (m0, m1), sg, order=2)
    if r2.status != STATUS_OK or abs(r2.cos_t) < 0.15:
        continue
    ntried += 1
    r50 = integrate_pair_2d((s0, s1), (m0, m1), sg, order=50)
    De, Me, ge = exact_reference(r2, (s0, s1), (m0, m1), sg)
    # scales: the integral MAGNITUDE, never the entry itself (entries can vanish)
    sc = r2.len
    gsc = r2.gapL1 if r2.gapL1 > 0 else sc
    e_2v50["D"] = max(e_2v50["D"], rel(2, r2.D, r50.D, sc))
    e_2v50["M"] = max(e_2v50["M"], rel(2, r2.M, r50.M, sc))
    e_2v50["g"] = max(e_2v50["g"], rel(2, r2.g, r50.g, gsc))
    e_2v50["len"] = max(e_2v50["len"], rel(2, r2.len, r50.len, sc))
    e_2vex["D"] = max(e_2vex["D"], rel(2, r2.D, De, sc))
    e_2vex["M"] = max(e_2vex["M"], rel(2, r2.M, Me, sc))
    e_2vex["g"] = max(e_2vex["g"], rel(2, r2.g, ge, gsc))
    e_50vex["D"] = max(e_50vex["D"], rel(2, r50.D, De, sc))
    e_50vex["M"] = max(e_50vex["M"], rel(2, r50.M, Me, sc))
    e_50vex["g"] = max(e_50vex["g"], rel(2, r50.g, ge, gsc))

TOL2 = 1e-13
check(2, f"2-pt vs 50-pt Gauss over {n2} random pairs: difference is ROUNDING, not "
         "quadrature",
      max(e_2v50.values()) < TOL2,
      "  ".join(f"{k} {v:.2e}" for k, v in e_2v50.items()))
check(2, "2-pt Gauss vs EXACT RATIONAL quadrature (fractions.Fraction): agreement at "
         "rounding level -- 2-pt IS exact for these degree<=2 integrands",
      max(e_2vex.values()) < TOL2, "  ".join(f"{k} {v:.2e}" for k, v in e_2vex.items()))
check(2, "50-pt Gauss vs EXACT RATIONAL: the same rounding level (the referee is the "
         "referee, not the finer float rule)",
      max(e_50vex.values()) < TOL2, "  ".join(f"{k} {v:.2e}" for k, v in e_50vex.items()))

# a 1-pt rule is NOT exact -- the gate has discriminating power.  (Deliberately NOT folded
# into the section max: this error is the point of the check, not a defect.)
Xs_d = [(1.0, 0.4), (8.0, 0.4)]
Xm_d = [(0.0, 0.0), (10.0, 0.0)]
rr1 = integrate_pair_2d(Xs_d, Xm_d, 1.0, order=1)
rr2 = integrate_pair_2d(Xs_d, Xm_d, 1.0, order=2)
e1 = float(np.max(np.abs(rr1.D - rr2.D))) / rr2.len
check(2, "discriminating power: a 1-pt rule is measurably WRONG on the same integrands "
         "(so 'exact' above is a result, not a tautology)", e1 > 1e-3,
      f"1pt vs 2pt rel = {e1:.3e}")

# ============================================================================
# CHECK 3 -- THE TILTED NON-PARALLEL CASE (the missing-Jacobian falsifier)
# ============================================================================
print("\nCHECK 3  tilted non-parallel pair: analytic slave-trace integrals vs the kernel,")
print("         and the DELIBERATELY BROKEN twin (1/|t_s.t_m| omitted) it must catch.")

Xm_t = [(0.0, 0.0), (20.0, 0.0)]          # master along +x, length 20 (power-of-two-safe)
Ls_t = 8.0
rows3 = []
bad3 = 0
for cth in (1.0, 0.9, 0.7, 0.5):
    sth = math.sqrt(1.0 - cth * cth)
    s0 = np.array([4.0, 0.75])
    s1 = s0 + Ls_t * np.array([cth, sth])
    Xs_t = [tuple(s0), tuple(s1)]
    rok = integrate_pair_2d(Xs_t, Xm_t, +1.0, jac_mode="slave_trace")
    rbad = integrate_pair_2d(Xs_t, Xm_t, +1.0, jac_mode="no_align")
    Ls_a, Da, Ma, ga = analytic_full_cover(Xs_t, Xm_t, +1.0)
    gsc = max(abs(ga[0]), abs(ga[1]))
    e_ok = max(rel(3, rok.len, Ls_a, Ls_a), rel(3, rok.D, Da, Ls_a),
               rel(3, rok.M, Ma, Ls_a), rel(3, rok.g, ga, gsc))
    e_bad = max(abs(rbad.len - Ls_a) / Ls_a,
                float(np.max(np.abs(rbad.D - Da))) / Ls_a,
                float(np.max(np.abs(rbad.M - Ma))) / Ls_a,
                float(np.max(np.abs(rbad.g - ga))) / gsc)
    # the broken twin must be WRONG BY EXACTLY cos(theta) on every quantity
    ratios3 = np.concatenate([[rbad.len / rok.len],
                              (rbad.D / rok.D).ravel(), (rbad.M / rok.M).ravel(),
                              (rbad.g / rok.g).ravel()])
    e_ratio = float(np.max(np.abs(ratios3 - abs(rok.cos_t))))
    bitident = (rbad.len == rok.len and np.array_equal(rbad.D, rok.D)
                and np.array_equal(rbad.M, rok.M) and np.array_equal(rbad.g, rok.g))
    rows3.append((cth, e_ok, e_bad, e_ratio, bitident))
    if e_ok >= 64 * EPS:
        bad3 += 1

TOL3 = 64 * EPS
print("          cos(th)   correct-vs-analytic   broken-vs-analytic   |broken/correct - cos|   bit-identical")
for cth, e_ok, e_bad, e_r, bi in rows3:
    print(f"          {cth:5.2f}      {e_ok:.3e}            {e_bad:.3e}             "
          f"{e_r:.3e}          {bi}")

check(3, "CORRECT kernel matches the closed-form slave-trace integrals (area, D, M, gtil) "
         "at every tilt", bad3 == 0,
      f"max rel = {max(r[1] for r in rows3):.3e}  (tol {TOL3:.2e})")
tilted = [r for r in rows3 if r[0] < 1.0]
check(3, "the analytic check CATCHES the missing-Jacobian twin at every tilt "
         "(cos = 0.9/0.7/0.5)",
      all(r[2] > 0.05 for r in tilted),
      "broken rel err = " + ", ".join(f"{r[2]:.3f}" for r in tilted))
check(3, "and it fails by EXACTLY the predicted cos(theta) factor on every entry of "
         "area/D/M/gtil",
      all(r[3] < 1e-14 for r in tilted),
      f"max |broken/correct - cos| = {max(r[3] for r in tilted):.3e}")
flat = [r for r in rows3 if r[0] == 1.0][0]
check(3, "a FLAT parallel patch (cos = 1) cannot distinguish the two -- BIT-IDENTICAL. "
         "This is why the tilted case exists.", flat[4] is True and flat[2] < TOL3,
      f"broken-vs-analytic = {flat[2]:.3e}, bit-identical = {flat[4]}")

# partition of unity is measure-INDEPENDENT and therefore also blind to the bug -- state it
r_flat = integrate_pair_2d([(4.0, 0.75), (12.0, 0.75)], Xm_t, +1.0)
r_tilt = integrate_pair_2d([(4.0, 0.75), (4.0 + 8.0 * 0.5, 8.0 * math.sqrt(0.75))],
                           Xm_t, +1.0, jac_mode="no_align")
pou = max(rel(3, r_tilt.M.sum(1), r_tilt.D.sum(1), r_tilt.len),
          rel(3, r_flat.M.sum(1), r_flat.D.sum(1), r_flat.len))
check(3, "partition of unity sum_k M_ik == sum_j D_ij holds even for the BROKEN twin "
         "(it is measure-independent -- another gate with no power here)",
      pou < 64 * EPS, f"max rel = {pou:.3e}")

# ============================================================================
# CHECK 4 -- MORTAR PATCH TEST (uniform pressure, NON-MATCHING interface)
# ============================================================================
print("\nCHECK 4  mortar patch test: uniform pressure p across a NON-MATCHING 2D interface")
print("         (slave 3 segments vs master 2 segments over the same span) + far-field twin.")


def patch(x0, y0, g0, tilt=None):
    """Slave 3 segments (nodes at 0,3,7,10) vs master 2 segments (nodes at 0,4.5,10).
    The slave TRACE length is 10 in both the flat and the tilted layout (the tilted chain
    is arc-length parametrized), and in both the whole slave projects inside the master."""
    mN, mC = chain(x0, [0.0, 4.5, 10.0], y0)
    if tilt is None:
        sN, sC = chain(x0, [0.0, 3.0, 7.0, 10.0], y0 + g0)
    else:
        cth = tilt
        sth = math.sqrt(1.0 - cth * cth)
        sN, sC = chain_dir((x0 + 0.5, y0 + g0), (cth, sth), [0.0, 3.0, 7.0, 10.0])
    D, M, g, area, tab = assemble_interface(sN, sC, mN, mC, +1.0)
    return sN, sC, mN, mC, D, M, g, area, 10.0


def patch_report(tag, x0, y0, g0, tilt=None, p=1.0, sec=4):
    sN, sC, mN, mC, D, M, g, area, trace = patch(x0, y0, g0, tilt)
    xscale = max(1.0, float(np.max(np.abs(np.concatenate([sN.ravel(), mN.ravel()])))))
    ctol = 256 * EPS * xscale
    ns = len(sN)
    # (a) total transferred force: slave side and master side both = p * slave trace
    lam = np.full(ns, p)
    fs = D @ lam
    fm = -(M.T @ lam)
    e_fs = abs(float(fs.sum()) - p * trace)
    e_fm = abs(float(fm.sum()) + p * trace)
    _secerr[sec] = max(_secerr.get(sec, 0.0), max(e_fs, e_fm) / max(1.0, p * trace))
    check(sec, f"[{tag}] uniform pressure: sum F_slave = p*trace and sum F_master = "
               f"-p*trace to machine precision",
          e_fs <= ctol * p and e_fm <= ctol * p,
          f"|dF_s| = {e_fs:.3e}, |dF_m| = {e_fm:.3e}  (tol {ctol * p:.3e}, "
          f"|x|_scale = {xscale:.3g})")
    # (b) self-equilibrium of the pair of force vectors
    e_eq = abs(float(fs.sum() + fm.sum()))
    check(sec, f"[{tag}] the transfer is SELF-EQUILIBRATING (sum F_s + sum F_m == 0)",
          e_eq <= ctol * p, f"residual = {e_eq:.3e}")
    # (c) partition of unity, per slave node
    e_pou = rel(sec, M.sum(1), D.sum(1), float(np.max(D.sum(1))))
    check(sec, f"[{tag}] partition of unity sum_k M_ik == sum_j D_ij at every slave node",
          e_pou < 256 * EPS, f"max rel = {e_pou:.3e}")
    # (d) measure: total assembled area == the slave trace length
    e_area = abs(area - trace) / trace
    _secerr[sec] = max(_secerr.get(sec, 0.0), e_area)
    check(sec, f"[{tag}] assembled measure sum(area) == slave trace length",
          abs(area - trace) <= ctol, f"rel = {e_area:.3e}")
    # (e) linear field reproduction: P = D^-1 M maps master nodal values of a LINEAR field
    #     to the slave nodal values EXACTLY (the mortar consistency condition)
    c0, c1 = 0.375, -0.25
    um = c0 + c1 * mN[:, 0]
    us = c0 + c1 * sN[:, 0]
    Pmat = np.linalg.solve(D, M)
    uscale = float(np.max(np.abs(us))) + abs(c1) * xscale
    e_lin = rel(sec, Pmat @ um, us, uscale)
    check(sec, f"[{tag}] linear field reproduced exactly by P = D^-1 M "
               "(mortar consistency)", e_lin < 256 * EPS,
          f"max rel = {e_lin:.3e}  (|u|_scale = {uscale:.3g})")
    # (f) weighted gap reproduces the NODAL gap field: D^-1 gtil == g(node) exactly.  The
    #     master is horizontal here, so n = (0,1) and g(node) = y_node - y_master.  For the
    #     flat patches that field is the CONSTANT g0; for the tilted patch it is LINEAR --
    #     both are in the slave space, so both must be reproduced exactly.
    gnod = sN[:, 1] - y0
    gb = np.linalg.solve(D, g)
    e_gap = float(np.max(np.abs(gb - gnod))) / xscale
    _secerr[sec] = max(_secerr.get(sec, 0.0), e_gap)
    kind = "uniform g0" if tilt is None else "LINEAR (tilted slave)"
    check(sec, f"[{tag}] weighted gap: D^-1 gtil reproduces the nodal gap field "
               f"[{kind}] exactly (relative to |x|_scale)",
          float(np.max(np.abs(gb - gnod))) <= ctol,
          f"max |dg| = {float(np.max(np.abs(gb - gnod))):.3e}  (tol {ctol:.3e})")
    return D, M, g


patch_report("flat, origin ", 0.0, 0.0, 1.0e-3)
patch_report("flat, far-fld", 5.0e4, -3.0e4, 1.0e-3)
patch_report("tilted cos0.8", 0.0, 0.0, 1.0e-3, tilt=0.8)

# the far-field twin under the BROKEN Jacobian: still passes the FLAT patch (the point)
_, _, _, _, Db, Mb, _, areab, _ = patch(0.0, 0.0, 1.0e-3)
sNf, sCf = chain(0.0, [0.0, 3.0, 7.0, 10.0], 1.0e-3)
mNf, mCf = chain(0.0, [0.0, 4.5, 10.0], 0.0)
Dbad, Mbad, gbad, areabad, _ = assemble_interface(sNf, sCf, mNf, mCf, +1.0,
                                                  jac_mode="no_align")
check(4, "the FLAT patch test passes IDENTICALLY under the broken Jacobian -- documented "
         "blindness, and the reason CHECK 3 exists",
      np.array_equal(Dbad, Db) and np.array_equal(Mbad, Mb) and areabad == areab,
      "bit-identical D, M, area")

# ============================================================================
# CHECK 5 -- TIE / T-JUNCTION BOOKKEEPING (the ADR-41 C4 crux)
# ============================================================================
print("\nCHECK 5  shared-node T-junction: a slave node sits exactly on a master "
      "discretization\n         boundary. Refusal bookkeeping, partition of unity, exact "
      "tie patch.")

# slave nodes 0, 5, 8, 10 (3 segments) -- the node at x = 5 is the T-junction, shared by
# slave segments 0 and 1 AND coincident with the master's interior node.
sNt, sCt = chain(0.0, [0.0, 5.0, 8.0, 10.0], 1.0e-3)
mNt, mCt = chain(0.0, [0.0, 5.0, 10.0], 0.0)
Dt, Mt, gt, areat, tabt = assemble_interface(sNt, sCt, mNt, mCt, +1.0)

OKC, EMP = (STATUS_OK, REFUSE_NONE), (STATUS_SKIP, REFUSE_EMPTY)
expected = {(0, 0): OKC, (0, 1): EMP,
            (1, 0): EMP, (1, 1): OKC,
            (2, 0): EMP, (2, 1): OKC}
got = {(es, em): st for es, em, st, _, _ in tabt}
check(5, "T-junction refusal bookkeeping: the four zero-measure pairs across the junction "
         "refuse with REFUSE_EMPTY, the three real pairs accept -- no double count, no "
         "lost slice", got == expected,
      f"{sum(1 for v in got.values() if v == EMP)} refusals / "
      f"{sum(1 for v in got.values() if v == OKC)} accepted")

e_meas = abs(areat - 10.0) / 10.0
_secerr[5] = max(_secerr.get(5, 0.0), e_meas)
check(5, "total assembled measure == the full slave trace (10.0) despite the junction",
      abs(areat - 10.0) < 256 * EPS * 10.0, f"rel = {e_meas:.3e}")

rowsum = Dt.sum(1)
tributary = np.array([2.5, 2.5 + 1.5, 1.5 + 1.0, 1.0])       # nodal tributary lengths
e_row = rel(5, rowsum, tributary, float(np.max(tributary)))
check(5, "D row-sums == the nodal tributary lengths (2.5 / 4.0 / 2.5 / 1.0); the junction "
         "node collects BOTH of its slave segments", e_row < 256 * EPS,
      f"max rel = {e_row:.3e}  rowsum = {np.array2string(rowsum, precision=15)}")

e_pou_t = rel(5, Mt.sum(1), Dt.sum(1), float(np.max(rowsum)))
check(5, "partition of unity at the T-junction node (the ADR-41 C4 crux): "
         "sum_k M_ik == sum_j D_ij", e_pou_t < 256 * EPS, f"max rel = {e_pou_t:.3e}")

e_tot = abs(float(Dt.sum()) - 10.0) / 10.0
check(5, "global partition of unity sum_ij D_ij == sum_ik M_ik == total measure",
      abs(float(Dt.sum()) - float(Mt.sum())) < 256 * EPS * 10.0
      and abs(float(Dt.sum()) - 10.0) < 256 * EPS * 10.0,
      f"|sum D - sum M| = {abs(float(Dt.sum()) - float(Mt.sum())):.3e}, "
      f"rel to 10 = {e_tot:.3e}")

# exact tie patch: a LINEAR displacement field is transferred exactly on the non-matching,
# T-junctioned mesh.  (The field is a function of x, i.e. its gradient is TANGENTIAL to the
# master; a normal-gradient component would legitimately differ by the gap -- see the note
# in the final report.)
for c0, c1 in ((1.0, 0.0), (0.375, -0.25), (-2.5, 0.125)):
    um = c0 + c1 * mNt[:, 0]
    us = c0 + c1 * sNt[:, 0]
    ut = np.linalg.solve(Dt, Mt @ um)
    usc = abs(c0) + abs(c1) * 10.0
    e_tie = rel(5, ut, us, usc)
    check(5, f"exact tie patch on the T-junctioned non-matching mesh, u = {c0} + {c1}*x",
          e_tie < 256 * EPS, f"max rel = {e_tie:.3e}")

# and a far-from-origin T-junction twin, with RELATIVE tolerances
sNt2, sCt2 = chain(5.0e4, [0.0, 5.0, 8.0, 10.0], -3.0e4 + 1.0e-3)
mNt2, mCt2 = chain(5.0e4, [0.0, 5.0, 10.0], -3.0e4)
Dt2, Mt2, gt2, areat2, tabt2 = assemble_interface(sNt2, sCt2, mNt2, mCt2, +1.0)
xsc2 = 5.0e4
um2 = 0.375 - 0.25 * mNt2[:, 0]
us2 = 0.375 - 0.25 * sNt2[:, 0]
ut2 = np.linalg.solve(Dt2, Mt2 @ um2)
e_tie2 = rel(5, ut2, us2, float(np.max(np.abs(us2))))
check(5, f"far-from-origin (|x| ~ {xsc2:.0e}) T-junction tie twin, RELATIVE tolerance",
      e_tie2 < 4096 * EPS, f"max rel = {e_tie2:.3e}")
check(5, "far-from-origin T-junction refusal bookkeeping unchanged",
      {(es, em): st for es, em, st, _, _ in tabt2} == expected)

# ============================================================================
# CHECK 6 -- WEIGHTED-GAP SIGN CONVENTION
# ============================================================================
print("\nCHECK 6  sign convention: penetration <=> g < 0, on both orientations.")

d_pen = 0.02
# (i) master along +x, slave BELOW it, sigma = +1  ->  n = (0,+1), slave penetrates
r_p = integrate_pair_2d([(2.0, -d_pen), (8.0, -d_pen)], [(0.0, 0.0), (10.0, 0.0)], +1.0)
check(6, "uniformly penetrated pair (sigma=+1): gtil_i < 0 for EVERY slave node",
      r_p.status == STATUS_OK and bool(np.all(r_p.g < 0.0)),
      f"gtil = {np.array2string(r_p.g, precision=6)}")
e_s6 = rel(6, r_p.g, -d_pen * r_p.D.sum(1), d_pen * r_p.len)
check(6, "and gtil_i == -d * (row sum of D) exactly (uniform gap factors out)",
      e_s6 < 64 * EPS, f"max rel = {e_s6:.3e}")

# (ii) the mirrored DECLARATION of the SAME physical state: master listed backwards (so
#      perp(t_m) flips) and both segments listed backwards; sigma = -1 restores the same
#      physical normal n = (0,+1).  The gap must be bit-identical, not merely negative.
r_m = integrate_pair_2d([(8.0, -d_pen), (2.0, -d_pen)], [(10.0, 0.0), (0.0, 0.0)], -1.0)
Pp = np.array([[0.0, 1.0], [1.0, 0.0]])
e_mir = rel(6, Pp @ r_m.g, r_p.g, float(np.max(np.abs(r_p.g))))
check(6, "mirrored DECLARATION (master reversed + sigma flip) of the same physical "
         "penetration: gtil_i < 0 still, and identical up to the slave permutation",
      r_m.status == STATUS_OK and bool(np.all(r_m.g < 0.0)) and e_mir < 64 * EPS,
      f"gtil = {np.array2string(r_m.g, precision=6)}, max rel = {e_mir:.3e}")
# and the SIGMA VOTE is what decides it: the same geometry with the wrong sigma reports
# separation.  (This is why sigma is an INPUT here -- the vote is handler business.)
r_bad_sig = integrate_pair_2d([(2.0, -d_pen), (8.0, -d_pen)], [(0.0, 0.0), (10.0, 0.0)], -1.0)
check(6, "sigma is load-bearing: the SAME penetrating geometry with the wrong sigma "
         "reports separation (gtil > 0) -- the orientation vote is the handler's job, and "
         "this oracle takes sigma as an input",
      bool(np.all(r_bad_sig.g > 0.0)) and np.array_equal(r_bad_sig.g, -r_p.g),
      f"gtil = {np.array2string(r_bad_sig.g, precision=6)}")

# (iii) separation is strictly positive, and the sign is not a coin flip near zero
r_sep = integrate_pair_2d([(2.0, +d_pen), (8.0, +d_pen)], [(0.0, 0.0), (10.0, 0.0)], +1.0)
check(6, "separated pair gives gtil_i > 0 for every slave node (no sign ambiguity)",
      bool(np.all(r_sep.g > 0.0)), f"gtil = {np.array2string(r_sep.g, precision=6)}")
check(6, "gapL1 = INT |g| dGamma_s is sign-free and equals |gtil| summed for a "
         "single-signed gap",
      abs(r_p.gapL1 - abs(float(r_p.g.sum()))) < 64 * EPS * r_p.gapL1,
      f"gapL1 = {r_p.gapL1:.12e}")

# a linear gap field that CHANGES SIGN inside the pair: gtil follows the field, not |g|
r_x = integrate_pair_2d([(2.0, -d_pen), (8.0, +d_pen)], [(0.0, 0.0), (10.0, 0.0)], +1.0)
_, _, _, ga_x = analytic_full_cover([(2.0, -d_pen), (8.0, +d_pen)],
                                    [(0.0, 0.0), (10.0, 0.0)], +1.0)
e_x = rel(6, r_x.g, ga_x, float(np.max(np.abs(ga_x))))
check(6, "sign-changing gap: gtil_0 < 0 < gtil_1 and both match the closed form "
         "(the normal is INSIDE the integral)",
      r_x.g[0] < 0.0 < r_x.g[1] and e_x < 64 * EPS, f"max rel = {e_x:.3e}")

# ============================================================================
# CHECK 7 -- ALM / UZAWA MODEL PROBLEM
# ============================================================================
print("\nCHECK 7  ALM/Uzawa: lambda <- min(0, lambda + epsN*gbar) on commit drives "
      "penetration\n         below an epsN-INDEPENDENT augTol; pure penalty does not.")

# --- 7a: the 1-DOF model problem -------------------------------------------------
# elastic body k, external load F pushing into a rigid master; gap g = -u (u = penetration
# direction).  traction t = min(0, lam + epsN*g)   (fork convention: compression NEGATIVE,
# proto_c2_alm.py nodal_pressure).  equilibrium  k*u + t_react = F  with t_react = -t.
k_b, F_ext = 250.0, 1000.0
u_ref = F_ext / k_b
augTol = 1e-9 * u_ref                 # RELATIVE: 1e-9 of the free-body displacement scale


def uzawa_1dof(epsN, augment=True, maxAug=400):
    lam, u = 0.0, 0.0
    for it in range(1, maxAug + 1):
        # solve k*u - t = F with t = min(0, lam + epsN*(-u))   (active while u > 0)
        u = (F_ext + lam) / (k_b + epsN)
        g = -u
        if -g < augTol:
            return u, lam, it
        if not augment:
            return u, lam, it
        lam = min(0.0, lam + epsN * g)
    return u, lam, maxAug


rows7 = []
alm_ok = True
for mult in (1.0, 10.0, 100.0, 1000.0):
    epsN = mult * k_b
    u, lam, it = uzawa_1dof(epsN)
    rows7.append((epsN, u, lam, it))
    if not (u < augTol):
        alm_ok = False
print("          epsN        final penetration     lambda            augmentations")
for epsN, u, lam, it in rows7:
    print(f"          {epsN:9.1f}   {u:.6e}        {lam:+.9f}      {it}")
check(7, f"[1-DOF] augmentation reaches |g| < augTol = {augTol:.2e} at EVERY epsN over 3 "
         "orders of magnitude", alm_ok,
      f"max penetration = {max(r[1] for r in rows7):.3e}")
lam_spread = max(abs(r[2]) for r in rows7) - min(abs(r[2]) for r in rows7)
check(7, "the converged multiplier is epsN-INDEPENDENT (lambda -> -F for every epsN)",
      lam_spread < 1e-6 * F_ext and all(abs(abs(r[2]) - F_ext) < 1e-6 * F_ext for r in rows7),
      f"spread = {lam_spread:.3e}, |lambda| - F = "
      f"{max(abs(abs(r[2]) - F_ext) for r in rows7):.3e}")
u_pen, _, _ = uzawa_1dof(1.0 * k_b, augment=False)
check(7, "DISCRIMINATING POWER: pure penalty (no augmentation) at the SOFTEST epsN does "
         "NOT reach augTol", u_pen > 1e5 * augTol,
      f"penalty penetration = {u_pen:.3e} = {u_pen / augTol:.2e} x augTol")

# --- 7b: the same statement on the REAL mortar operators -------------------------
# 2-node slave on a 1-segment master; gap field g(xi) = g0 - sum N_j u_j, so
# gtil = g0*w - D u  with w = D.rowsum (partition of unity) -- the mortar structure itself.
r_alm = integrate_pair_2d([(2.0, 0.0), (8.0, 0.0)], [(0.0, 0.0), (10.0, 0.0)], +1.0)
Dm = r_alm.D.copy()
w_alm = Dm.sum(1)
K_b = 400.0 * np.array([[2.0, -1.0], [-1.0, 2.0]])
f_b = np.array([900.0, 1500.0])
g0_alm = -4.0e-3          # initial INTERFERENCE, so the converged state is NOT the trivial
#                           u = 0 (D^-1 w == 1 by partition of unity => u_inf = g0*ones)
u_free = float(np.max(np.abs(np.linalg.solve(K_b, f_b))))
augTol2 = 1e-9 * u_free


def uzawa_mortar(epsN, augment=True, maxAug=1200):
    """k u + (contact) = f with gtil = g0*w - D u, gbar = gtil/w,
    t = min(0, lam + epsN*gbar) (compression NEGATIVE -- proto_c2_alm.py convention),
    contact force conjugate to u = -D t.  Uzawa on commit: lam <- min(0, lam + epsN*gbar)."""
    lam = np.zeros(2)
    u = np.zeros(2)
    gbar = g0_alm * np.ones(2)
    A = K_b + epsN * (Dm @ np.diag(1.0 / w_alm) @ Dm)
    for it in range(1, maxAug + 1):
        u = np.linalg.solve(A, f_b + Dm @ (lam + epsN * g0_alm))
        gbar = (g0_alm * w_alm - Dm @ u) / w_alm
        pen = float(np.max(np.maximum(0.0, -gbar)))
        if pen < augTol2:
            return u, gbar, lam, it, True
        if not augment:
            return u, gbar, lam, it, False
        lam = np.minimum(0.0, lam + epsN * gbar)
    return u, gbar, lam, maxAug, False


ok7b, its7b, us7b, ls7b = True, [], [], []
for mult in (1.0, 30.0, 1000.0):
    epsN = mult * 400.0
    u, gbar, lam, it, conv = uzawa_mortar(epsN)
    its7b.append(it)
    us7b.append(u)
    ls7b.append(lam)
    if not conv:
        ok7b = False
check(7, f"[mortar 2-DOF, real D] augmentation reaches ||penetration||_inf < augTol = "
         f"{augTol2:.2e} at every epsN (3 orders)", ok7b,
      f"augmentations = {its7b}")
e7u = rel(7, us7b[0], us7b[-1], abs(g0_alm))
e7l = rel(7, ls7b[0], ls7b[-1], float(np.max(np.abs(ls7b[0]))))
check(7, "the converged mortar state is epsN-INDEPENDENT (u -> g0*ones by partition of "
         "unity; lambda -> the exact reaction)",
      e7u < 1e-5 and e7l < 1e-4, f"max rel |du| = {e7u:.3e}, |dlambda| = {e7l:.3e}")
e7g = rel(7, us7b[0], g0_alm * np.ones(2), abs(g0_alm))
check(7, "and it converges to the CLOSED-FORM limit u_inf = g0 * ones (D^-1 w == 1)",
      e7g < 1e-5, f"max rel = {e7g:.3e}")
u_pen2, gbar_pen2, _, _, _ = uzawa_mortar(400.0, augment=False)
check(7, "DISCRIMINATING POWER: mortar pure penalty at the softest epsN misses augTol",
      float(np.max(np.maximum(0.0, -gbar_pen2))) > 1e5 * augTol2,
      f"penalty penetration = {float(np.max(np.maximum(0.0, -gbar_pen2))):.3e}")

# ============================================================================
# CHECK 8 -- ANISOTROPIC (THICKNESS) SCALING -- ADR-85 SS How/7
# ============================================================================
print("\nCHECK 8  -thickness h: h scales tractions/forces by EXACTLY h and leaves every")
print("         geometric quantity BIT-IDENTICAL. lambda is NOT separately scaled.")


def return_map_1d(epsT, ell, slip, slip_p, N, mu, coh, tau_max, h):
    """ADR-85 SS How/4 scalar return map with the SS How/7 h-scaled unified cone
    min(mu*N + c*h*l, tau_max*h*l).  N already carries h (it comes from epsN*h)."""
    kt = epsT * h * ell
    tr = kt * (slip - slip_p)
    cap = mu * abs(N) + coh * h * ell
    if tau_max > 0.0:
        cap = min(cap, tau_max * h * ell)
    if abs(tr) <= cap:
        return tr, slip_p, True
    t = math.copysign(cap, tr)
    return t, slip - t / kt, False


# geometry: one pair, penetrating.  h enters ONLY at epsN / epsT / the stress clamps.
Xs8 = [(2.0, -0.02), (8.0, -0.02)]
Xm8 = [(0.0, 0.0), (10.0, 0.0)]
r8 = integrate_pair_2d(Xs8, Xm8, +1.0)
w8 = r8.D.sum(1)
gbar8 = r8.g / w8
epsN8, epsT8 = 5.0e5, 3.0e5
lam8_unit = np.array([-120.0, -80.0])     # the committed multiplier of the h = 1 run
ell8 = r8.len
slip8, slipp8 = 3.0e-4, 1.0e-4
mu8, coh8, tmax8 = 0.4, 25.0, 1.0e9

res = {}
for h in (1.0, 2.0):
    # ADR-85 SS How/7: lambda is NOT separately scaled -- it ACCUMULATES from epsN*h*gbar
    # in the Uzawa update and therefore INHERITS h.  The h-twin's committed multiplier is
    # h*lambda_1 by construction; scaling it again in the residual is the h^2 error the
    # rev-2 panel caught, and is the falsifier at the bottom of this section.
    lam_h = h * lam8_unit
    t8 = np.minimum(0.0, lam_h + (epsN8 * h) * gbar8)
    Fs8 = r8.D @ t8
    Fm8 = -(r8.M.T @ t8)
    N8 = abs(float(t8 @ w8))
    tT, spN, stick = return_map_1d(epsT8, ell8, slip8, slipp8, N8, mu8, coh8, tmax8, h)
    res[h] = dict(t=t8, Fs=Fs8, Fm=Fm8, N=N8, tT=tT, spN=spN, stick=stick)

# geometry BIT-identical (h never touches it)
r8h = integrate_pair_2d(Xs8, Xm8, +1.0)
check(8, "geometry is BIT-IDENTICAL under -thickness (interval [a,b], area, D, M, gtil, "
         "gbar): h never enters the integrator",
      r8h.a == r8.a and r8h.b == r8.b and r8h.len == r8.len
      and np.array_equal(r8h.D, r8.D) and np.array_equal(r8h.M, r8.M)
      and np.array_equal(r8h.g, r8.g),
      f"[a,b] = [{r8.a}, {r8.b}], area = {r8.len!r}")

bit2 = (np.array_equal(res[2.0]["t"], 2.0 * res[1.0]["t"])
        and np.array_equal(res[2.0]["Fs"], 2.0 * res[1.0]["Fs"])
        and np.array_equal(res[2.0]["Fm"], 2.0 * res[1.0]["Fm"]))
check(8, "h = 2 scales the normal traction and BOTH force vectors by exactly 2 "
         "(bit-identical, h a power of two)", bit2,
      f"t(h=1) = {np.array2string(res[1.0]['t'], precision=9)}")
check(8, "h = 2 scales the friction traction by exactly 2 and leaves the committed "
         "plastic slip BIT-IDENTICAL (a geometric quantity)",
      res[2.0]["tT"] == 2.0 * res[1.0]["tT"] and res[2.0]["spN"] == res[1.0]["spN"]
      and res[2.0]["stick"] == res[1.0]["stick"],
      f"tT(h=1) = {res[1.0]['tT']:.12e}, stick = {res[1.0]['stick']}, "
      f"slip_p = {res[1.0]['spN']:.12e}")

# the stick/slip DECISION must be h-invariant -- this is the whole point of h-scaling the
# stress clamps (an unscaled clamp makes the threshold wrong by exactly h).
def clamp_unscaled(epsT, ell, slip, slip_p, N, mu, coh, tau_max, h):
    kt = epsT * h * ell
    tr = kt * (slip - slip_p)
    cap = mu * abs(N) + coh * ell           # BUG: cohesion clamp NOT h-scaled
    if tau_max > 0.0:
        cap = min(cap, tau_max * ell)       # BUG: tau_max clamp NOT h-scaled
    if abs(tr) <= cap:
        return tr, slip_p, True
    t = math.copysign(cap, tr)
    return t, slip - t / kt, False


# A COHESION-DOMINATED state, chosen so the unscaled-clamp bug is decisive: the cone at
# h = 1 is mu*N1 + c*l = 400 + 30000; the traction demand is 0.9 of it.  At h = 2 the
# correctly h-scaled cone doubles (still stick), while the unscaled-clamp twin leaves the
# cohesion term at c*l and therefore SLIPS -- wrong by exactly the factor h on the clamp.
N_c, mu_c, coh_c, ell_c = 1000.0, 0.4, 5000.0, 6.0
cap1 = mu_c * N_c + coh_c * ell_c
slip_c = 1.0e-4 + 0.9 * cap1 / (epsT8 * ell_c)
a_ok = return_map_1d(epsT8, ell_c, slip_c, 1.0e-4, N_c, mu_c, coh_c, 0.0, 1.0)
b_ok = return_map_1d(epsT8, ell_c, slip_c, 1.0e-4, 2.0 * N_c, mu_c, coh_c, 0.0, 2.0)
a_bad = clamp_unscaled(epsT8, ell_c, slip_c, 1.0e-4, N_c, mu_c, coh_c, 0.0, 1.0)
b_bad = clamp_unscaled(epsT8, ell_c, slip_c, 1.0e-4, 2.0 * N_c, mu_c, coh_c, 0.0, 2.0)
check(8, "the stick/slip DECISION is h-invariant with h-scaled clamps, and the "
         "unscaled-clamp twin FLIPS it (the SS How/7 panel catch, reproduced)",
      a_ok[2] == b_ok[2] and a_bad[2] != b_bad[2],
      f"h-scaled clamps: stick {a_ok[2]}/{b_ok[2]}   unscaled clamps: {a_bad[2]}/{b_bad[2]}")

# lambda must NOT be scaled a SECOND time in the residual: h*(h*lambda_1) is an h^2 error
# on the lambda part.  The falsifier must both differ from the correct h-twin AND differ by
# exactly (h-1)*h*lambda_1 -- i.e. the error is entirely in the multiplier term.
t_h2_wrong = np.minimum(0.0, 2.0 * (2.0 * lam8_unit) + (epsN8 * 2.0) * gbar8)
delta = t_h2_wrong - res[2.0]["t"]
e_lam2 = float(np.max(np.abs(delta - 2.0 * lam8_unit)))
check(8, "lambda must NOT be separately scaled: the h^2 twin is CAUGHT, and its error is "
         "exactly (h-1)*h*lambda_1 -- entirely in the multiplier term",
      float(np.max(np.abs(delta))) > 1e-6 * float(np.max(np.abs(res[2.0]["t"])))
      and e_lam2 < 1e-9,
      f"|t_wrong - t_correct|_inf = {float(np.max(np.abs(delta))):.6e}, "
      f"residual vs (h-1)*h*lambda_1 = {e_lam2:.3e}")

# a non-power-of-two h still scales to rounding level (h = 2 is bitwise; h = 1.7 is not
# representable-exact, so the claim there is 'to rounding', stated honestly).
h3 = 1.7
t_a = np.minimum(0.0, 1.0 * lam8_unit + (epsN8 * 1.0) * gbar8)
t_b = np.minimum(0.0, h3 * lam8_unit + (epsN8 * h3) * gbar8)
e_h3 = rel(8, t_b, h3 * t_a, float(np.max(np.abs(t_a))) * h3)
check(8, f"a non-power-of-two thickness h = {h3} scales tractions by h to rounding level "
         "(bitwise is only claimed for representable h)",
      e_h3 < 8 * EPS, f"max rel = {e_h3:.3e}")

# ==================================================================== summary
print("\n" + "=" * 78)
print("PER-SECTION SUMMARY")
titles = {
    1: "interval clip / refusals / tau_ov",
    2: "2-pt Gauss exactness (vs 50-pt and vs exact rational)",
    3: "tilted non-parallel case + missing-Jacobian falsifier",
    4: "mortar patch test (non-matching, + far-field twin)",
    5: "tie / T-junction bookkeeping (ADR-41 C4)",
    6: "weighted-gap sign convention",
    7: "ALM / Uzawa augmentation",
    8: "-thickness anisotropic scaling",
}
notes = {7: "  [= the augTol convergence residual, by design -- not a deviation "
            "from an exact reference]"}
tot, totf = 0, 0
for s in sorted(_sec):
    n, f = _sec[s]
    tot += n
    totf += f
    print(f"  CHECK {s}  {n - f:2d}/{n:<2d} pass   max rel err {_secerr.get(s, 0.0):.3e}"
          f"   {titles[s]}" + ("" if f == 0 else f"   <-- {f} FAIL"))
    if s in notes:
        print(f"           {notes[s]}")
print("-" * 78)
print("MEASURED / VALIDATED CONSTANTS")
print(f"  tau_ov   = {TAU_OV:.0e}  (RELATIVE to min(L_s,L_m)/L_m -- the largest parametric")
print("                    overlap a pair can attain; the 2D analogue of the 3D A_min)")
print(f"  tau_align= {TAU_ALIGN:.0e}  (|t_s.t_m| floor -> REFUSE_ALIGN; MUST exceed tau_ov)")
print(f"  tau_seg  = {TAU_SEG:.0e}  (zero-length segment, relative to Lref)")
print(f"  far-field (|x| ~ {X0:.0e}) algorithmic clip noise vs exact rational : {noise:.3e}"
      f"  ({noise / EPS:.1f} ulp)")
print(f"  far-field (|x| ~ {X0:.0e}) TRANSLATION clip noise (what tau_ov clears): {tnoise:.3e}"
      f"  ({TAU_OV / max(tnoise, 1e-300):.0f}x margin)")
print(f"  fine-slave cost: at L_s/L_m = 1e-6 the interval's own relative accuracy is "
      f"{fine_area_err:.2e}")
print("                    (b-a is an O(eps) difference of O(1) parametric coordinates;")
print("                     documented, not a defect -- D/M/gtil all scale with b-a)")
print("=" * 78)
print(f"{tot - totf}/{tot} PASS" + ("" if _fails == 0 else f"   <-- {_fails} FAILURE(S)")
      + "   (ADR-85 T3 oracle, gate G-T3)")
print("=" * 78)
sys.exit(1 if _fails else 0)
