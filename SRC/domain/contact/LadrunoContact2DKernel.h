// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-85 T1a -- LadrunoContact2DKernel: header-only, OpenSees-free geometry for the
// 2D (-ndm 2) node-to-segment (NTS) contact lane. Slave POINT vs 2-node master LINE
// SEGMENT, plus the 2D corner primitive (the VERTEX pair) that the out-of-bounds
// refusal leaves uncovered. NO OpenSees types -- only raw double[] + <cmath>; numpy-
// oracle-testable build-free. Mirrors the LadrunoContactProjection.h /
// LadrunoEdgeKernel.h discipline (raw arrays, refusal status codes, relative gauges,
// every tolerance justified in a comment).
//
// T1a is GEOMETRY ONLY. There is deliberately NO state, NO ownership arbitration and
// NO traction here: the committed orientation sigma, the committed vertex side sign,
// and the segment-vs-vertex ownership dedup all live in the handler (T1b). What the
// kernel DOES supply is every pure PREDICATE the handler needs, so T1b adds no
// geometry of its own:
//   projectSegment2D()  -- closed-form point-on-segment projection + oriented normal +
//                          gap, with the four-way refusal status (ADR-85 SS How/1).
//   vertexWedge2D()     -- the corner CLASSIFIER: is the slave in this vertex's Voronoi
//                          wedge, is the corner convex or concave from the contact side,
//                          and which side of the corner is the slave on.
//   vertexEval2D()      -- the vertex PAIR: radial normal + signed gap from a COMMITTED
//                          side sign (ADR-85 SS How/1 vertex policy).
//   sigmaFromRef2D()    -- the 2D sibling of LadrunoContactProjection::normalOriented's
//                          refDir test (the auto-kn / orientation-vote primitive,
//                          ADR-85 SS Where + SS How/2).
//   bOperatorSegment2D() / bOperatorVertex2D() / shape2D().
//
// The 3D projection's bounded Gauss-Newton, its detK guard and its tolR/tolP escape
// hatch are all STRUCTURALLY ABSENT: in 2D the projection is one division.
//
// -------------------------------------------------------------------------------
// CONCAVE-VERTEX DECISION (ADR-85 "Risks / open questions" -> vertex-policy sign
// robustness). SETTLED IN T1a BY THE ORACLE, and it INVERTS the ADR's premise:
//
//   * A CONVEX corner's Voronoi wedge lies entirely OUTSIDE the master body, so a
//     slave in it has gap = +||xs - Xv|| > 0 STRICTLY. A convex vertex pair is a
//     valid, C0-continuous extension of its two neighbours but it can NEVER be
//     penetrating, hence NEVER carries force. (The interior of a convex corner is
//     covered TWICE -- both adjacent strips contain it -- which is the pre-existing
//     NTS corner double-count, not a hole.)
//   * A CONCAVE (re-entrant) corner is the opposite: its free-side strips OVERLAP, so
//     there is no exterior wedge, while the region behind the corner INSIDE the
//     material is refused by BOTH adjacent segments (aPrev > 0 and aNext < 0) and is
//     ENTIRELY PENETRATING. That is the real coverage hole: a finite-measure wedge, at
//     finite penetration, with zero restoring force -- silent pass-through.
//
//   => POLICY: concave vertices DO activate and DO need the vertex rule; they are the
//      ONLY vertex pairs that can ever carry force. Convex vertex pairs are kept for
//      the ownership map (exactly-one-owner) and for normal continuity while sliding
//      around a corner, but a handler that instantiates vertex pairs ONLY at concave
//      corners loses no force. Evidence: contact_prototypes/proto_t1_nts2d.py families
//      CA/CB/CC/CD (the concave sweep carries gap == -R over a wedge arc that policy
//      "concave never activates" would leave unowned).
// -------------------------------------------------------------------------------
//
// TOLERANCE SEMANTICS (ADR-85 SS How/1 -- every gauge here is RELATIVE or
// dimensionless; the fork has been bitten twice by absolute floors, ADR-57 P5 and the
// projection tolR incident):
//   tolIn  -- PARAMETRIC slack on xi. xi is dimensionless by construction, so 1e-9 is
//             unit-safe at any model scale (same value and same role as
//             LadrunoContactProjection::inBounds).
//   tauSeg -- zero-length refusal, RELATIVE to a caller-supplied Lref. CALLER CONTRACT:
//             Lref is the surface's minimum INITIAL segment length, computed ONCE per
//             surface at handle() and passed unchanged; a segment is refused when
//             L < tauSeg*Lref. Lref <= 0 (or NaN) is a contract violation and is
//             REFUSED, never silently treated as "no guard".
//   tolLen -- vertexEval2D's degenerate-coincidence floor. It is a LENGTH and has no
//             default on purpose: a hardcoded absolute length is exactly the unit trap
//             ADR-57 P5 fixed. Keep it RELATIVE to Lref. T1b AMENDMENT (measured): the
//             original "pass tauSeg*Lref" suggestion collides with the fork's standard
//             1e-8 seeded penetration (tauSeg*Lref = 1.118e-8 on the G-T1b notch deck
//             refused its own seed -> singular slave DOF; boundary probed at exactly
//             tauSeg*Lref). The floor guards DIRECTION conditioning, not segment
//             degeneracy -- pass the conditioning gauge tauPerp*Lref instead, as the
//             T1b adapter does.
//   tauPerp -- dimensionless conditioning gate for the two ambiguity fail-safes
//             (refDir tangent to the segment; a fold-back "spike" corner whose two
//             outward normals are antiparallel). Mirrors normalOriented's 1e-12.
//
// GAP IS NOT SCALE-FREE. xi is a ratio of dot products and is clean, but
// gap = (xs - x(xi)).n is a difference of large coordinates: it carries the
// eps*|x| ABSOLUTE noise floor (~1e-11 at coordinates ~5e4). Oracle tolerances on gap
// are therefore RELATIVE to the coordinate scale, never a bare 1e-12 absolute -- the
// far-from-origin family in proto_t1_nts2d.py measures that floor against exact
// (arbitrary-precision) arithmetic and records the drift.
//
// ORIENTATION. sigma (+1/-1) is an INPUT here (ADR-85 SS How/2): the interface-level
// reference vote, its named degenerate abort and the -outward override all live in the
// handler (T1b). The kernel never guesses a sign; sigmaFromRef2D returns 0 rather than
// pick one when the reference direction is tangent.
//
// Oracle: Ladruno_implementation/contact_prototypes/proto_t1_nts2d.py (numpy reference
// + >=1e5 randomized parity trials against t1a_cpp_check.cpp + the deterministic corner
// falsifiers). See Ladruno_implementation/85_ladruno_contact_2d_adr.md.

#ifndef LadrunoContact2DKernel_h
#define LadrunoContact2DKernel_h

#include <cmath>

namespace LadrunoContact2DKernel {

// --------------------------------------------------------------- status codes
// projectSegment2D: OUT_LOW / OUT_HIGH are REFUSALS, not clamps (ADR-85 SS How/1).
// They are split because the vertex-wedge routing needs to know WHICH end the slave
// ran past: the shared vertex of two segments is the HIGH end of the previous segment
// and the LOW end of the next one.
enum SegStatus {
    SEG2D_IN_BOUNDS  =  0,   // xi in [-tolIn, 1+tolIn]
    SEG2D_OUT_LOW    =  1,   // xi < -tolIn        (past X0)
    SEG2D_OUT_HIGH   =  2,   // xi > 1+tolIn       (past X1)
    SEG2D_DEGENERATE = -1    // zero-length segment, or a bad Lref contract
};

enum VtxStatus {
    VTX2D_OK         =  0,
    VTX2D_DEGENERATE = -1    // slave sits ON the vertex (||xs-Xv|| < tolLen): no direction
};

// Corner type SEEN FROM THE CONTACT (outward) SIDE. The numeric values are chosen so
// that corner == sign(turnSin) and the geometric side sign of the wedge is -corner.
enum CornerType {
    C2D_CONVEX   = -1,   // outward normals fan OUT: wedge is exterior, gap > 0, never active
    C2D_STRAIGHT =  0,   // collinear: the wedge is empty (the predicate cannot fire)
    C2D_CONCAVE  = +1    // re-entrant: wedge is interior, gap < 0, ALWAYS penetrating
};

// ------------------------------------------------------------------- gauges
static const double TOL_IN_BOUNDS_DEFAULT = 1e-9;    // parametric slack on xi (dimensionless)
static const double TAU_SEG_DEFAULT       = 1e-8;    // zero-length refusal, RELATIVE to Lref
static const double TAU_PERP_DEFAULT      = 1e-12;   // ambiguity gate (tangent refDir / spike corner)

// ------------------------------------------------------------------ small vec2
inline double dot2(const double a[2], const double b[2]) { return a[0]*b[0] + a[1]*b[1]; }
inline double cross2(const double a[2], const double b[2]) { return a[0]*b[1] - a[1]*b[0]; }
inline double norm2(const double a[2]) { return std::sqrt(dot2(a, a)); }
// perp(t) = (-t_y, t_x): the +90deg rotation. n = sigma*perp(t)/L (ADR-85 SS How/1).
inline void perp2(const double a[2], double p[2]) { p[0] = -a[1]; p[1] = a[0]; }

// ------------------------------------------------------------ shape functions
// Linear 2-node segment, x(xi) = (1-xi)*X0 + xi*X1, xi in [0,1]. N[0] is the X0 weight.
inline void shape2D(double xi, double N[2]) { N[0] = 1.0 - xi; N[1] = xi; }

// --------------------------------------------------- closed-form 2D projection
// Project the slave point xs onto the master segment X0->X1 (ADR-85 SS How/1):
//
//   t  = X1 - X0,  L2 = t.t,  L = sqrt(L2)
//   xi = (xs - X0).t / L2                       [closed form -- no Newton, no detK guard]
//   n  = sigma * perp(t) / L                    [sigma is an INPUT: SS How/2]
//   g  = (xs - x(xi)) . n                       [penetration <=> g < 0, fork convention]
//
// REFUSALS (never a clamp -- clamping is what re-creates the corner double-count the
// vertex policy exists to avoid):
//   * L < tauSeg*Lref, or Lref <= 0 / NaN  -> SEG2D_DEGENERATE, ALL outputs zeroed
//     (the projectFull precedent: a refused evaluation never leaks stale geometry).
//   * xi outside [-tolIn, 1+tolIn]         -> SEG2D_OUT_LOW / SEG2D_OUT_HIGH.
//
// OUT-OF-BOUNDS OUTPUT CONTRACT (read this before consuming them). On OUT_LOW/OUT_HIGH
// xi, n and gap ARE filled, because the wedge routing and the diagnostics need them.
// The gap returned there is the signed distance to the segment's INFINITE LINE (it is
// identically (xs - X0).n, since x(xi) - X0 is parallel to t and t.n == 0) -- it is
// NOT a contact gap, the closest master point is an ENDPOINT, and the vertex pair owns
// that region. Using it as a contact gap silently re-introduces the clamp.
//
// The gap is computed through the shape functions rather than through the algebraically
// identical cross-product form sigma*cross2(t, xs-X0)/L so that the segment residual and
// the B-operator see the SAME x(xi); the oracle asserts the two forms agree to the
// coordinate-scale noise floor.
inline int projectSegment2D(const double X0[2], const double X1[2], const double xs[2],
                            double sigma, double Lref,
                            double &xi, double &gap, double n[2],
                            double tolIn  = TOL_IN_BOUNDS_DEFAULT,
                            double tauSeg = TAU_SEG_DEFAULT)
{
    xi = 0.0; gap = 0.0; n[0] = 0.0; n[1] = 0.0;

    const double t[2] = { X1[0] - X0[0], X1[1] - X0[1] };
    const double L2 = dot2(t, t);
    const double L  = std::sqrt(L2);

    // Lref is a CALLER CONTRACT (min initial segment length of the surface, computed
    // once at handle()). A missing/garbage Lref would disable the zero-length guard
    // silently, so it is a refusal, not a fallback. The !(x > 0) form also traps NaN.
    if (!(Lref > 0.0))        return SEG2D_DEGENERATE;
    if (L < tauSeg * Lref)    return SEG2D_DEGENERATE;

    const double d0[2] = { xs[0] - X0[0], xs[1] - X0[1] };
    xi = dot2(d0, t) / L2;

    n[0] = sigma * (-t[1]) / L;
    n[1] = sigma * ( t[0]) / L;

    double N[2];
    shape2D(xi, N);
    const double xbar[2] = { N[0]*X0[0] + N[1]*X1[0],
                             N[0]*X0[1] + N[1]*X1[1] };
    const double d[2] = { xs[0] - xbar[0], xs[1] - xbar[1] };
    gap = dot2(d, n);

    if (xi < -tolIn)       return SEG2D_OUT_LOW;
    if (xi > 1.0 + tolIn)  return SEG2D_OUT_HIGH;
    return SEG2D_IN_BOUNDS;
}

// ------------------------------------------------------------- corner classifier
// Result of vertexWedge2D. Pure classification -- no state, no commit; the handler
// owns both (T1b).
struct WedgeResult {
    int    inWedge;    // 1 <=> xs is in this vertex's Voronoi wedge (BOTH adjacent
                       //       projections out of bounds ON THE VERTEX SIDE)
    int    corner;     // CornerType, = sign(turnSin). INFORMATIONAL / policy only --
                       //       do NOT take the side sign from it (see sideSign).
    double sideSign;   // +1 xs on the exterior (separated) side of the corner,
                       // -1 on the interior (penetrating) side,
                       //  0 UNDETERMINED (fold-back spike corner, xs at the vertex, or xs
                       //    inside the conditioning gate of the wedge boundary) -> DEFER.
    double turnSin;    // sigma * cross2(t^_prev, t^_next): signed sine of the turn angle
    double aPrev;      // (xs-Xv).t_prev/L_prev^2  ==  xi_prev - 1   (vertex-side <=> > tolIn)
    double aNext;      // (xs-Xv).t_next/L_next^2  ==  xi_next       (vertex-side <=> < -tolIn)
};

// The 2D corner primitive (ADR-85 SS How/1 vertex policy). tPrev and tNext are the two
// adjacent segments' tangent vectors IN SURFACE ORDER at the shared vertex Xv:
//   tPrev = Xv - X_prev   (the segment that ENDS at Xv)
//   tNext = X_next - Xv   (the segment that STARTS at Xv)
// sigma is the surface's committed orientation (the same sigma both segments use).
//
// WEDGE PREDICATE. The vertex owns xs iff xs ran off BOTH adjacent segments TOWARD the
// shared vertex:  xi_prev > 1+tolIn  AND  xi_next < -tolIn. Written vertex-relative as
// aPrev > tolIn and aNext < -tolIn -- algebraically identical (xi_prev = 1 + aPrev) but
// better conditioned, and it needs no reconstruction of X_prev / X_next.
//   * The "toward the vertex" qualifier is load-bearing: a slave that is out of bounds
//     on both segments at their FAR ends is past the patch boundary and is correctly
//     owned by nobody.
//   * The predicate is the ALGEBRAIC complement of projectSegment2D's in-bounds test at
//     the same tolIn -- but NOT the floating-point complement, and T1b MUST NOT assume
//     it is. projectSegment2D forms xi_prev from X_prev, so at a corner xi_prev ~ 1 and
//     the comparison xi > 1+tolIn can only resolve the overshoot to ulp(1) = 2.22e-16,
//     while aPrev resolves it exactly. MEASURED (oracle family CC, 16400 stations
//     scanned across the seams): a raw XOR of the two predicates produced 80
//     double-owner AND 1 zero-owner stations, inside a band of 1.15e-15 in the
//     parametric coordinate (5.2 ulp(1), i.e. 1.2e-6 of tolIn). No tolerance value fixes
//     that -- the two quantities simply have different resolutions.
//
//     OWNERSHIP RULE (total, unique, and free of any tolerance coupling -- the ADR-57 E7
//     ordered-ownership precedent; T1b implements exactly this and needs no geometry):
//        1. previous segment SEG2D_IN_BOUNDS   -> the PREVIOUS SEGMENT owns
//        2. else next segment SEG2D_IN_BOUNDS  -> the NEXT SEGMENT owns
//        3. else aPrev > 0 AND aNext < 0       -> the VERTEX owns  (UNSLACKED: this is
//           what closes the band the slacked inWedge leaves open)
//        4. else                               -> nobody (genuinely off the patch)
//     Segment precedence removes the double count; the unslacked step 3 removes the
//     hole. Where the two overlap they describe the SAME contact: measured relative gap
//     difference <= 1.7e-14 and normal disagreement <= 1.06e-7 rad, the latter matching
//     its analytic bound atan(tolIn*L/|gap|) to 1.00006 (family CC). Which candidate
//     wins there is physically immaterial -- only that exactly one does.
//
// SIDE SIGN -- derived from the BISECTOR of the two outward normals, NOT from the corner
// type and NOT from either normal alone. The wedge is exactly the cone of half-width
// |turn|/2 about sideSign*(n_prev + n_next) (oracle family CD verifies this on 40000
// random corners), so the bisector dot test is the one criterion that is correct for
// every |turn| < pi. Two tempting alternatives are NOT:
//   * sign((xs-Xv).n_prev) -- agrees with the side only for |turn| <= pi/2. Past that the
//     wedge fan is wider than a right angle, so a legitimate wedge direction can sit more
//     than 90 deg from ONE of the two bounding normals and that single dot flips.
//   * -sign(turnSin) -- correct, but it coin-flips at a near-straight corner where
//     turnSin is at the rounding floor, and again near a fold-back where it is the
//     difference of two nearly antiparallel tangents.
// Inside the wedge sideSign is either -corner or 0 (deferred); the oracle asserts that,
// including on corners with turn angles down to 1e-11 rad.
//   sideSign == 0 (a fold-back "spike" corner whose outward normals are antiparallel, a
//   slave within the conditioning gate of the wedge boundary, or xs exactly at the
//   vertex) means AMBIGUOUS: the handler must DEFER the capture and retry at the next,
//   better-conditioned configuration -- the LadrunoEdgeKernel::edgeGap deferral
//   precedent. Never guess. Measured over 48902 randomized wedge hits (with a
//   deliberately over-sampled 25% near-fold-back population): 0 wrong signs, 2412
//   deferrals, and every deferral on a near-antiparallel corner -- i.e. exactly the
//   region where sign(turnSin), and therefore `corner`, is itself a coin flip.
//
// HANDLER FLOW THIS IS DESIGNED FOR (T1b adds no geometry):
//   1. w = vertexWedge2D(...);            if (!w.inWedge)      -> this vertex does not own xs
//   2. if (w.sideSign == 0.0)             -> DEFER the capture (ambiguous corner)
//   3. at FIRST capture COMMIT sigmaSide = w.sideSign; thereafter pass the COMMITTED
//      value (ADR-57 lesson: a per-step re-derived sign flips exactly while
//      interpenetrating). The committed value and the live w.sideSign agree for as long
//      as the corner keeps its type, so a disagreement is a real diagnostic.
//   4. vertexEval2D(Xv, xs, committedSigmaSide, ...) -> gap < 0 <=> active.
//   w.corner is for POLICY only: only C2D_CONCAVE vertices can ever carry force (see the
//   concave-vertex decision at the top of this file).
inline WedgeResult vertexWedge2D(const double tPrev[2], const double tNext[2],
                                 double sigma, const double xs[2], const double Xv[2],
                                 double Lref,
                                 double tolIn   = TOL_IN_BOUNDS_DEFAULT,
                                 double tauSeg  = TAU_SEG_DEFAULT,
                                 double tauPerp = TAU_PERP_DEFAULT)
{
    WedgeResult w;
    w.inWedge = 0; w.corner = C2D_STRAIGHT; w.sideSign = 0.0;
    w.turnSin = 0.0; w.aPrev = 0.0; w.aNext = 0.0;

    const double LP2 = dot2(tPrev, tPrev), LN2 = dot2(tNext, tNext);
    const double LP  = std::sqrt(LP2),     LN  = std::sqrt(LN2);

    // Same refusal as projectSegment2D, on the same gauge: a vertex whose adjacent
    // segment the projection would refuse must not be classified either, or the
    // ownership map goes inconsistent (a slave owned by nobody, or by two features).
    if (!(Lref > 0.0))                                   return w;
    if (LP < tauSeg * Lref || LN < tauSeg * Lref)         return w;

    const double r[2] = { xs[0] - Xv[0], xs[1] - Xv[1] };
    w.aPrev = dot2(r, tPrev) / LP2;
    w.aNext = dot2(r, tNext) / LN2;
    w.inWedge = (w.aPrev > tolIn && w.aNext < -tolIn) ? 1 : 0;

    // signed turn angle seen from the outward side. No near-straight threshold is
    // applied: a razor-thin wedge is still EXACTLY the region both segments refuse, so
    // suppressing it would open a hole of that measure. The classification is only a
    // policy hint; the load-bearing sign is sideSign below.
    w.turnSin = sigma * (tPrev[0]*tNext[1] - tPrev[1]*tNext[0]) / (LP * LN);
    if (w.turnSin < 0.0)      w.corner = C2D_CONVEX;
    else if (w.turnSin > 0.0) w.corner = C2D_CONCAVE;
    else                      w.corner = C2D_STRAIGHT;

    const double nP[2] = { sigma*(-tPrev[1])/LP, sigma*(tPrev[0])/LP };
    const double nN[2] = { sigma*(-tNext[1])/LN, sigma*(tNext[0])/LN };
    const double bis[2] = { nP[0] + nN[0], nP[1] + nN[1] };
    const double bl = norm2(bis);          // in [0,2]; ~0 only for a fold-back spike
    const double rl = norm2(r);
    // CONDITIONING GATE -- scaled by |r| ONLY, deliberately NOT by |bis|. bis is a
    // CANCELLING sum of two UNIT vectors, so its components carry an absolute error of
    // ~eps no matter how small |bis| gets; a threshold that also scales with |bis| (the
    // naive "make it relative to everything" reading) sinks BELOW that noise floor as
    // the corner folds back and then admits a coin-flip sign. MEASURED on 48902 wedge
    // hits with 25% near-fold-back corners: the |r|*|bis| threshold produced 707 WRONG
    // side signs; this |r|-only threshold produces 0 (2412 deferrals instead).
    if (bl > 2.0*tauPerp && rl > 0.0) {
        const double p = dot2(r, bis);
        if (std::fabs(p) > tauPerp * rl) w.sideSign = (p > 0.0) ? 1.0 : -1.0;
    }
    return w;
}

// ------------------------------------------------------------------ vertex pair
// The 2D corner contact pair (ADR-85 SS How/1): slave point vs a single master vertex.
//
//   n   = sigmaSide * (xs - Xv)/||xs - Xv||
//   gap = sigmaSide * ||xs - Xv||            (penetration <=> gap < 0)
//   B   = [ n | -n ]   (ndof = 4: slave x,y then vertex x,y -- bOperatorVertex2D)
//
// sigmaSide is the COMMITTED side sign from vertexWedge2D (an INPUT -- the commit lives
// in the handler). Note the structural difference from a segment or an edge-edge pair:
// ||xs - Xv|| is a DISTANCE, so no motion of the slave can make the gap change sign by
// itself. The sign MUST come from the side classification; that is why it is committed
// rather than re-derived, and why sideSign == 0 has to defer instead of defaulting.
//
// Consistency with the neighbours (this is what makes the vertex pair an extension and
// not a patch): on the wedge boundary aPrev = 0 the vertex (gap, n) equals the previous
// segment's evaluation at xi = 1 exactly, and on aNext = 0 it equals the next segment's
// at xi = 0. Both seams are gated in the oracle (families CB / CD).
//
// tolLen is REQUIRED and is a LENGTH, kept RELATIVE to the surface scale (a hardcoded
// absolute floor here would be the ADR-57 P5 unit trap -- a model meshed in small units
// silently loses all vertex contact). Pass tauPerp*Lref, the DIRECTION-CONDITIONING
// gauge (see the T1b amendment in the tolerance block at the top of this file: the
// original tauSeg*Lref suggestion sat ON the fork's standard 1e-8 seeded penetration
// and refused its own seed). tolLen <= 0 / NaN is a contract violation -> refuse.
inline int vertexEval2D(const double Xv[2], const double xs[2], double sigmaSide,
                        double &gap, double n[2], double tolLen)
{
    gap = 0.0; n[0] = 0.0; n[1] = 0.0;

    const double r[2] = { xs[0] - Xv[0], xs[1] - Xv[1] };
    const double d = norm2(r);
    if (!(tolLen > 0.0)) return VTX2D_DEGENERATE;   // contract violation (also traps NaN)
    if (d < tolLen)      return VTX2D_DEGENERATE;   // slave sits on the vertex: no direction

    n[0] = sigmaSide * r[0] / d;
    n[1] = sigmaSide * r[1] / d;
    gap  = sigmaSide * d;
    return VTX2D_OK;
}

// ------------------------------------------------------- orientation from a datum
// 2D sibling of LadrunoContactProjection::normalOriented's refDir test (ADR-85 SS Where:
// the auto-kn path calls the tri/quad-only normalOriented and needs a 2D sibling; the
// SS How/2 interface-level vote needs the same primitive). Returns the sigma that makes
// the segment normal point along refDir:
//   proj = perp(t).refDir / L ;  sigma = sign(proj)
// FAIL-SAFE: |proj| < tauPerp*||refDir|| means refDir is (numerically) TANGENT to the
// segment -- the outward sense is genuinely ambiguous (the flush-interface case ADR-85
// SS How/2 calls out). Return 0 and let the caller abort or demand -outward. NEVER guess
// a sign: the shipped 3D default's silent drop in exactly this configuration is the
// behaviour ADR-85 rev 2 replaced with a named abort.
inline double sigmaFromRef2D(const double t[2], const double refDir[2],
                             double tauPerp = TAU_PERP_DEFAULT)
{
    const double L = norm2(t);
    const double R = norm2(refDir);
    if (!(L > 0.0) || !(R > 0.0)) return 0.0;
    double p[2];
    perp2(t, p);
    const double proj = dot2(p, refDir) / L;
    if (std::fabs(proj) < tauPerp * R) return 0.0;
    return (proj > 0.0) ? 1.0 : -1.0;
}

// ---------------------------------------------------------------- B-operators
// FIRST variation of the gap, d(gap) = B . [du_slave, du_master...]. In 2D both
// operators are the EXACT GRADIENT -- which grounds ADR-85 SS How/5 keeping
// K_c = kn*B^T B as the MAIN tangent term:
//   segment: xs - x(xi) is exactly gap*n (the closest-point residual is parallel to n
//            because n is perpendicular to the ONE tangent), and n.dn = 0 for a unit
//            normal, so the dn and dxi FIRST-ORDER contributions cancel identically.
//   vertex:  d||xs-Xv|| = u.(dxs - dXv) with u = (xs-Xv)/||xs-Xv|| = sigmaSide*n.
// Both identities are FD-gated in the oracle (family F8). NOTE (T1b correction of an
// earlier over-read): an exact first variation does NOT make the SECOND variation zero.
// The Newton geometric term kn*g*(d2 g/du2) is nonzero whenever the MASTER moves (a
// rotating segment rotates n; the vertex distance carries the radial (I - u u^T)/||r||
// Hessian) -- it is a NAMED DEFERRAL (the 3D default path omits its analogue too), so
// kn*B^T B is Gauss-Newton, exact for a FIXED master, superlinear otherwise.
//
// DOF ORDER, segment: [slave_x, slave_y, X0_x, X0_y, X1_x, X1_y]  (ndof = 6)
//   B = [ n | -N0*n | -N1*n ],  r = B^T * tn,  tn = kn*<-gap>_+
inline void bOperatorSegment2D(double N0, double N1, const double n[2], double B[6])
{
    B[0] =  n[0];      B[1] =  n[1];
    B[2] = -N0*n[0];   B[3] = -N0*n[1];
    B[4] = -N1*n[0];   B[5] = -N1*n[1];
}

// DOF ORDER, vertex: [slave_x, slave_y, Xv_x, Xv_y]  (ndof = 4);  B = [ n | -n ]
inline void bOperatorVertex2D(const double n[2], double B[4])
{
    B[0] =  n[0];   B[1] =  n[1];
    B[2] = -n[0];   B[3] = -n[1];
}

// =============================================================================
// ADR-85 T3 -- 2D MORTAR: the slave-trace interval integrator (segment-to-segment
// mortar for 2 two-node segments). Header-only, pure append below the T1a/T1b
// NTS geometry above (byte-stable): NO edit above this marker.
//
// Mirrors LadrunoMortarKernel.h's PairResult/integratePair discipline (ADR-41 C1)
// collapsed from the 3D clip -> sub-triangle Gauss to the 1D interval clip -> 2-pt
// Gauss the ADR-85 SS How/3 kinematics prescribe: project the slave ENDPOINTS onto
// the master's affine parametrization, clip [xiM(s0),xiM(s1)] against [0,1], and
// integrate on the SLAVE-TRACE measure dGamma_s = (L_m/|t_s.t_m|)*dxi_m -- the 2D
// analogue of the shipped 3D kernel's A_phys = A_aux/cos_t (LadrunoMortarKernel.h:460).
// Oracle: Ladruno_implementation/contact_prototypes/proto_t3_mortar2d.py (70/70 PASS);
// this function is a line-for-line C++ transcription of integrate_pair_2d() there --
// same operation order (no reassociated dot products), same refusal ordering, same
// thresholds -- so it is testable bit-for-bit against the numpy oracle.
// =============================================================================

// A2 (proto_t3_mortar2d.py) -- status stays the {0,-1} pair binding decision 2b.1
// mandates (0 accumulated, -1 refused/skip; NO -2 tier -- straight 2-node segments
// cannot be "invalid geometry" the way quadratic facets can), but PairResult2D ALSO
// carries a `refusal` CLASSIFICATION and the `cos_t` the ALIGN test already computed
// -- set BEFORE returning on every path that reaches it -- so the HANDLER can route
// WARN_2D_PERP_NO_NTS / the NTS-ownership stand-down (ADR-85 How/3, decision 2b.4)
// WITHOUT re-deriving cos_t from the coordinates (the fork's staleness trap, applied
// to the exact quantity the gate turns on).
enum Mort2DStatus { MORT2D_OK = 0, MORT2D_SKIP = -1 };
enum Mort2DRefusal {
    MORT2D_REF_NONE  = 0,
    MORT2D_REF_EMPTY = 1,   // clipped overlap b-a <= tauOv*min(Ls,Lm)/Lm (sliver / disjoint)
    MORT2D_REF_DEGEN = 2,   // zero-length slave OR master segment (tauSeg*Lref gauge, the
                            // SAME relative gauge projectSegment2D uses)
    MORT2D_REF_ALIGN = 3    // |t_hat_s . t_hat_m| < tauAlign -- near-perpendicular pair;
                            // the slave-trace Jacobian 1/|cos| is ill-conditioned there
};

// TOLERANCE SEMANTICS (both RELATIVE gauges -- never a bare absolute; ADR-85 T3 design
// note SS3, measured not guessed):
//   MORT2D_TAU_OV_DEFAULT    -- the parametric-overlap SLIVER guard, RELATIVE to
//       min(Ls,Lm)/Lm: that ratio IS the largest parametric overlap a pair can attain
//       (a fully-contained slave gives Ls*|cos|/Lm, a fully-covered master gives 1), so
//       it is the exact 2D analogue of the shipped 3D kernel's A_min sliver normalizer
//       (LadrunoMortarKernel.h:444). Default 1e-9: the far-field TRANSLATION-noise floor
//       measured at |x|~5e4 is ~9e-12 (oracle CHECK 1.7), so 1e-9 buys ~100x headroom.
//   MORT2D_TAU_ALIGN_DEFAULT -- the alignment floor on |t_hat_s.t_hat_m|; caps the
//       slave-trace Jacobian 1/|cos| at 1e6, far past any well-posed pair.
// ORDERING IS A HARD, GATED REQUIREMENT (oracle CHECK 1.6, proto_t3_mortar2d.py A3):
// MORT2D_TAU_OV_DEFAULT MUST stay far BELOW MORT2D_TAU_ALIGN_DEFAULT. For a fully
// CONTAINED slave the clipped overlap is exactly b-a = Ls*|cos|/Lm against the sliver
// threshold tauOv*Ls/Lm, so the sliver guard bites EXACTLY when |cos| < tauOv. If the
// ordering were reversed (tauOv >= tauAlign), the EMPTY sliver guard would silently
// swallow a legitimate, strongly-tilted FULL overlap as "no contribution" before the
// ALIGN classification could route it to WARN_2D_PERP_NO_NTS / the NTS-ownership
// stand-down -- a silent hole in exactly the case the routing protocol exists to catch.
// mortarSegmentPair2D therefore tests ALIGN strictly before EMPTY; do not reorder.
static const double MORT2D_TAU_OV_DEFAULT    = 1e-9;
static const double MORT2D_TAU_ALIGN_DEFAULT = 1e-6;

// The 2D mortar pair result: 2 slave nodes x 2 master nodes (straight 2-node segments
// only -- there is no serendipity/quadratic 2D facet). Mirrors
// LadrunoMortarKernel::PairResult's field list (status/ngp/area~len/gapL1/D/M/g) plus
// the A2 refusal+cos_t fields the 2D routing protocol needs.
struct PairResult2D {
    int    status;    // MORT2D_OK (0, accumulated) or MORT2D_SKIP (-1, refused/skip)
    int    refusal;   // Mort2DRefusal; MORT2D_REF_NONE unless status == MORT2D_SKIP
    double cos_t;      // t_hat_s . t_hat_m -- set as soon as computed (before the ALIGN
                       // test), so it is valid on ALIGN and on a later EMPTY refusal;
                       // stays 0.0 on a DEGEN refusal (never reached -> never computed)
    int    ngp;        // Gauss points integrated this pair (0 or 2; 0 <=> no contribution)
    double len;        // INT dGamma_s over the clipped overlap (the slave-trace measure)
    double gapL1;       // INT |gN| dGamma_s -- the sign-free gap guard (mirrors
                        // LadrunoMortarKernel::PairResult::gapL1)
    double D[2][2];    // slave-slave mortar D_ij = INT N_i^s N_j^s dGamma_s
    double M[2][2];    // slave-master mortar M_ik = INT N_i^s N_k^m dGamma_s
    double g[2];       // weighted gap g~_i = INT N_i^s gN dGamma_s -- gN computed INSIDE
                       // the integral (mirrors LadrunoMortarKernel.h:484, the ADR-41 D1
                       // flat-facet fix: factoring gN out of the integral is wrong
                       // whenever the gap is not uniform along the overlap, which a
                       // tilted/rotated pair always is)
};

// Integrate ONE slave/master 2-node segment pair on the slave-trace measure (ADR-85
// SS How/3). Xs/Xm = [[x,y],[x,y]] segment endpoints (CURRENT config -- caller's
// choice of reference vs trial); sigma = +-1 the committed interface orientation vote
// (How/2, a caller INPUT here, never re-derived); Lref = the surface's min INITIAL
// segment length (the SAME relative-gauge caller contract projectSegment2D/
// vertexWedge2D use, computed once per contact at handle()). Zeroes `out` on entry (a
// refused evaluation never leaks stale geometry -- the projectSegment2D precedent).
//
// tauOv/tauAlign/tauSeg default to the oracle-validated constants above; a caller
// passing non-default values is responsible for preserving the tauOv << tauAlign
// ordering (see the block comment above).
//
// Refusal order (never reorder -- see the ordering note above): DEGEN (either segment
// too short) -> ALIGN (near-perpendicular) -> EMPTY (sliver/disjoint clip) -> OK.
inline int mortarSegmentPair2D(const double Xs[2][2], const double Xm[2][2],
                               double sigma, double Lref, PairResult2D &out,
                               double tauOv    = MORT2D_TAU_OV_DEFAULT,
                               double tauAlign = MORT2D_TAU_ALIGN_DEFAULT,
                               double tauSeg   = TAU_SEG_DEFAULT)
{
    out.status = MORT2D_SKIP; out.refusal = MORT2D_REF_EMPTY; out.cos_t = 0.0;
    out.ngp = 0; out.len = 0.0; out.gapL1 = 0.0;
    for (int i = 0; i < 2; i++) {
        out.g[i] = 0.0;
        for (int j = 0; j < 2; j++) { out.D[i][j] = 0.0; out.M[i][j] = 0.0; }
    }

    const double ts[2] = { Xs[1][0] - Xs[0][0], Xs[1][1] - Xs[0][1] };
    const double tm[2] = { Xm[1][0] - Xm[0][0], Xm[1][1] - Xm[0][1] };
    const double Ls = norm2(ts), Lm = norm2(tm);

    // Lref is a CALLER CONTRACT, exactly as in projectSegment2D: a missing/garbage
    // Lref would silently disable the zero-length guard, so it is a refusal, not a
    // fallback. !(x>0) also traps NaN.
    if (!(Lref > 0.0))               { out.refusal = MORT2D_REF_DEGEN; return out.status; }
    if (Ls <= tauSeg * Lref || Lm <= tauSeg * Lref) {
        out.refusal = MORT2D_REF_DEGEN; return out.status;
    }

    // alignment: cos(theta) between the UNIT tangents. 1/|cos| is the load-bearing
    // slave-trace Jacobian factor below. Set cos_t as soon as it is known (A2).
    const double tsHat[2] = { ts[0] / Ls, ts[1] / Ls };
    const double tmHat[2] = { tm[0] / Lm, tm[1] / Lm };
    const double cosT = dot2(tsHat, tmHat);
    out.cos_t = cosT;
    if (std::fabs(cosT) < tauAlign) {           // near-perpendicular: ill-posed measure
        out.refusal = MORT2D_REF_ALIGN;
        return out.status;
    }

    // closed-form projection of the slave ENDPOINTS onto the master parametrization.
    const double Lm2 = dot2(tm, tm);
    const double d0[2] = { Xs[0][0] - Xm[0][0], Xs[0][1] - Xm[0][1] };
    const double d1[2] = { Xs[1][0] - Xm[0][0], Xs[1][1] - Xm[0][1] };
    const double al0 = dot2(d0, tm) / Lm2;
    const double al1 = dot2(d1, tm) / Lm2;
    const double lo = (al0 < al1) ? al0 : al1;
    const double hi = (al0 < al1) ? al1 : al0;
    const double a  = (lo > 0.0) ? lo : 0.0;
    const double b  = (hi < 1.0) ? hi : 1.0;

    // sliver guard -- RELATIVE to min(Ls,Lm)/Lm (see the tolerance-semantics block
    // above); tested AFTER alignment, never before (the tauOv << tauAlign ordering).
    const double thr = tauOv * ((Ls < Lm) ? Ls : Lm) / Lm;
    if (b - a <= thr) {
        out.refusal = MORT2D_REF_EMPTY;
        return out.status;
    }

    const double nx = -sigma * tm[1] / Lm;      // n = sigma*perp(t_m)/L_m (How/1 convention)
    const double ny =  sigma * tm[0] / Lm;

    const double J = Lm / std::fabs(cosT);      // dGamma_s = (L_m/|t_s.t_m|) dxi_m
    const double dal = al1 - al0;                // nonzero: |cosT| >= tauAlign > 0 (see
                                                 // the design note: dal == Ls*cosT/Lm)

    // 2-POINT GAUSS, HARDCODED -- exact for these degree<=2 integrands (oracle CHECK 2,
    // proto_t3_mortar2d.py; a 1-pt rule is measurably wrong there). No order argument:
    // the ADR-85 T3 binding decision (2b.7 / How/3) fixes 2-pt for the 2D mortar lane.
    static const double GP2_U0 = 0.21132486540518711775;  // 0.5 - 1/(2*sqrt(3))
    static const double GP2_U1 = 0.78867513459481288225;  // 0.5 + 1/(2*sqrt(3))
    static const double GP2_W  = 0.5;
    const double u[2] = { GP2_U0, GP2_U1 };
    const double span = b - a;

    for (int gp = 0; gp < 2; gp++) {
        const double xim = a + span * u[gp];
        const double wgt = GP2_W * span;
        const double xis = (xim - al0) / dal;
        const double Ns[2] = { 1.0 - xis, xis };
        const double Nm[2] = { 1.0 - xim, xim };
        const double xsx = Xs[0][0] + xis * ts[0];
        const double xsy = Xs[0][1] + xis * ts[1];
        const double xmx = Xm[0][0] + xim * tm[0];
        const double xmy = Xm[0][1] + xim * tm[1];
        // gN computed INSIDE the loop, from the SAME xis/xim the shape functions use --
        // one consistent contact point (mirrors LadrunoMortarKernel.h:484).
        const double gN = (xsx - xmx) * nx + (xsy - xmy) * ny;
        const double wJ = wgt * J;
        out.len   += wJ;
        out.gapL1 += wJ * std::fabs(gN);
        for (int i = 0; i < 2; i++) {
            out.g[i] += wJ * Ns[i] * gN;
            for (int j = 0; j < 2; j++) out.D[i][j] += wJ * Ns[i] * Ns[j];
            for (int k = 0; k < 2; k++) out.M[i][k] += wJ * Ns[i] * Nm[k];
        }
        out.ngp++;
    }
    out.status = MORT2D_OK;
    out.refusal = MORT2D_REF_NONE;
    return out.status;
}

} // namespace LadrunoContact2DKernel

#endif
