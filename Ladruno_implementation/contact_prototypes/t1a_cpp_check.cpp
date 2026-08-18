// ADR-85 T1a -- standalone C++ check + hexfloat parity driver for
// SRC/domain/contact/LadrunoContact2DKernel.h. No OpenSees, no build:
//
//   g++ -O2 -ffp-contract=off -I SRC/domain/contact t1a_cpp_check.cpp -o t1a_cpp_check
//   ./t1a_cpp_check                 # deterministic self-check (the e0_cpp_check precedent)
//   ./t1a_cpp_check --trials        # read hexfloat trial records on stdin, echo results
//
// -ffp-contract=off is REQUIRED for the parity mode: with GCC's default (=fast) an
// a*b - c*d could contract to an FMA on an FMA-capable target and stop matching the
// numpy oracle bit-for-bit. proto_t1_nts2d.py compiles with that flag.
//
// TRIAL RECORD GRAMMAR (whitespace-separated; every double is a C99 hexfloat so the
// round-trip python <-> C++ is EXACT -- 'S'egment, 'W'edge, 'V'ertex, 'R'ef-orientation):
//   S X0x X0y X1x X1y xsx xsy sigma Lref tolIn tauSeg
//     -> S status xi gap nx ny B0 B1 B2 B3 B4 B5
//   W tPx tPy tNx tNy sigma xsx xsy Xvx Xvy Lref tolIn tauSeg tauPerp
//     -> W inWedge corner sideSign turnSin aPrev aNext
//   V Xvx Xvy xsx xsy sigmaSide tolLen
//     -> V status gap nx ny B0 B1 B2 B3
//   R tx ty refx refy tauPerp
//     -> R sigma
//
// Mirrors Ladruno_implementation/contact_prototypes/proto_t1_nts2d.py.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "LadrunoContact2DKernel.h"
namespace K = LadrunoContact2DKernel;

static const double PI = 3.14159265358979323846;
static int fails = 0;

static void check(const char *name, bool ok, const char *extra = "") {
    std::printf("  [%s] %s%s%s\n", ok ? "PASS" : "FAIL", name, extra[0] ? " -- " : "", extra);
    if (!ok) fails++;
}

// ===================================================================== trial mode
static double rdD() {
    char b[80];
    if (std::scanf("%79s", b) != 1) std::exit(2);
    return std::strtod(b, 0);
}

static int runTrials() {
    char tag[8];
    while (std::scanf("%7s", tag) == 1) {
        if (tag[0] == 'S') {
            double X0[2], X1[2], xs[2], sigma, Lref, tolIn, tauSeg;
            X0[0]=rdD(); X0[1]=rdD(); X1[0]=rdD(); X1[1]=rdD();
            xs[0]=rdD(); xs[1]=rdD(); sigma=rdD(); Lref=rdD(); tolIn=rdD(); tauSeg=rdD();
            double xi, gap, n[2], N[2], B[6];
            int st = K::projectSegment2D(X0, X1, xs, sigma, Lref, xi, gap, n, tolIn, tauSeg);
            K::shape2D(xi, N);
            K::bOperatorSegment2D(N[0], N[1], n, B);
            std::printf("S %d %a %a %a %a %a %a %a %a %a %a\n", st, xi, gap, n[0], n[1],
                        B[0], B[1], B[2], B[3], B[4], B[5]);
        } else if (tag[0] == 'W') {
            double tP[2], tN[2], sigma, xs[2], Xv[2], Lref, tolIn, tauSeg, tauPerp;
            tP[0]=rdD(); tP[1]=rdD(); tN[0]=rdD(); tN[1]=rdD(); sigma=rdD();
            xs[0]=rdD(); xs[1]=rdD(); Xv[0]=rdD(); Xv[1]=rdD();
            Lref=rdD(); tolIn=rdD(); tauSeg=rdD(); tauPerp=rdD();
            K::WedgeResult w = K::vertexWedge2D(tP, tN, sigma, xs, Xv, Lref, tolIn, tauSeg, tauPerp);
            std::printf("W %d %d %a %a %a %a\n", w.inWedge, w.corner, w.sideSign,
                        w.turnSin, w.aPrev, w.aNext);
        } else if (tag[0] == 'V') {
            double Xv[2], xs[2], sgn, tolLen;
            Xv[0]=rdD(); Xv[1]=rdD(); xs[0]=rdD(); xs[1]=rdD(); sgn=rdD(); tolLen=rdD();
            double gap, n[2], B[4];
            int st = K::vertexEval2D(Xv, xs, sgn, gap, n, tolLen);
            K::bOperatorVertex2D(n, B);
            std::printf("V %d %a %a %a %a %a %a %a\n", st, gap, n[0], n[1],
                        B[0], B[1], B[2], B[3]);
        } else if (tag[0] == 'R') {
            double t[2], ref[2], tauPerp;
            t[0]=rdD(); t[1]=rdD(); ref[0]=rdD(); ref[1]=rdD(); tauPerp=rdD();
            std::printf("R %a\n", K::sigmaFromRef2D(t, ref, tauPerp));
        } else {
            std::fprintf(stderr, "t1a_cpp_check: unknown record tag '%s'\n", tag);
            return 2;
        }
    }
    return 0;
}

// ================================================================ self-check utils
// One corner = the polyline Xprev -> Xv -> Xnext with orientation sigma. `which` is
// 0 = prev segment, 1 = next segment, 2 = vertex, -1 = nobody; the owner is resolved by
// the ordered rule documented in LadrunoContact2DKernel.h (see ownerOf below).
struct Corner {
    double Xprev[2], Xv[2], Xnext[2], tP[2], tN[2], sigma, Lref;
};

static Corner makeCorner(const double Xp[2], const double Xv[2], const double Xn[2], double sigma) {
    Corner c;
    for (int d = 0; d < 2; d++) { c.Xprev[d]=Xp[d]; c.Xv[d]=Xv[d]; c.Xnext[d]=Xn[d]; }
    for (int d = 0; d < 2; d++) { c.tP[d]=Xv[d]-Xp[d]; c.tN[d]=Xn[d]-Xv[d]; }
    c.sigma = sigma;
    double lp = K::norm2(c.tP), ln = K::norm2(c.tN);
    c.Lref = (lp < ln) ? lp : ln;
    return c;
}

struct Owned {
    int nOwners, which, stP, stN;
    double gap, n[2];
    double sideSign;
};

// The ORDERED OWNERSHIP RULE documented in LadrunoContact2DKernel.h (segment precedence,
// then the UNSLACKED wedge). A raw XOR of "segment in-bounds" against the slacked
// inWedge is NOT floating-point exact -- it double-counts and holes inside a +/-ulp(1)
// band at the seam (oracle family CC). This rule is total and unique with no tolerance
// coupling. nOwners counts the CANDIDATES (>=1 always, in the band it can be 2).
static Owned ownerOf(const Corner &c, const double xs[2], double committedSide) {
    Owned o;
    double xiP, gP, nP[2], xiN, gN, nN[2];
    o.stP = K::projectSegment2D(c.Xprev, c.Xv,    xs, c.sigma, c.Lref, xiP, gP, nP);
    o.stN = K::projectSegment2D(c.Xv,    c.Xnext, xs, c.sigma, c.Lref, xiN, gN, nN);
    K::WedgeResult w = K::vertexWedge2D(c.tP, c.tN, c.sigma, xs, c.Xv, c.Lref);
    o.sideSign = w.sideSign;
    bool vtxClaim = (w.aPrev > 0.0 && w.aNext < 0.0);       // UNSLACKED fallback claim
    o.nOwners = (o.stP == K::SEG2D_IN_BOUNDS ? 1 : 0)
              + (o.stN == K::SEG2D_IN_BOUNDS ? 1 : 0)
              + (vtxClaim ? 1 : 0);
    o.which = -1; o.gap = 0.0; o.n[0] = o.n[1] = 0.0;
    if (o.stP == K::SEG2D_IN_BOUNDS) {
        o.which = 0; o.gap = gP; o.n[0]=nP[0]; o.n[1]=nP[1];
    } else if (o.stN == K::SEG2D_IN_BOUNDS) {
        o.which = 1; o.gap = gN; o.n[0]=nN[0]; o.n[1]=nN[1];
    } else if (vtxClaim) {
        double side = (committedSide != 0.0) ? committedSide : w.sideSign;
        if (K::vertexEval2D(c.Xv, xs, side, o.gap, o.n, 1e-8*c.Lref) == K::VTX2D_OK) o.which = 2;
    }
    return o;
}

// offset sweep at signed distance off (off > 0 = exterior, off < 0 = into the material):
// along the previous segment, around the vertex arc, along the next segment. Returns the
// point at path parameter u in [0,1] together with the exact arclength step used.
struct Sweep {
    const Corner *c;
    double off, nP[2], nN[2], a0, a1, L1, Larc, L2, Ltot;
};

static Sweep makeSweep(const Corner &c, double off, double frac) {
    Sweep s;
    s.c = &c; s.off = off;
    double lp = K::norm2(c.tP), ln = K::norm2(c.tN);
    s.nP[0] = c.sigma*(-c.tP[1])/lp;  s.nP[1] = c.sigma*(c.tP[0])/lp;
    s.nN[0] = c.sigma*(-c.tN[1])/ln;  s.nN[1] = c.sigma*(c.tN[0])/ln;
    double sgn = (off >= 0.0) ? 1.0 : -1.0;
    s.a0 = std::atan2(sgn*s.nP[1], sgn*s.nP[0]);
    s.a1 = std::atan2(sgn*s.nN[1], sgn*s.nN[0]);
    double d = s.a1 - s.a0;                       // take the SHORT way round
    while (d >  PI) d -= 2.0*PI;
    while (d < -PI) d += 2.0*PI;
    s.a1 = s.a0 + d;
    double R = std::fabs(off);
    s.L1 = frac*lp; s.L2 = frac*ln; s.Larc = R*std::fabs(d);
    s.Ltot = s.L1 + s.Larc + s.L2;
    return s;
}

// phase: 0 = on the prev segment's offset line, 1 = on the arc, 2 = on the next segment's
static void sweepPoint(const Sweep &s, double sArc, double xs[2], int *phase) {
    const Corner &c = *s.c;
    double lp = K::norm2(c.tP), ln = K::norm2(c.tN);
    if (sArc <= s.L1) {                                   // walk toward the vertex
        double back = s.L1 - sArc;
        for (int d = 0; d < 2; d++) xs[d] = c.Xv[d] - back*c.tP[d]/lp + s.off*s.nP[d];
        *phase = 0;
    } else if (sArc <= s.L1 + s.Larc) {
        double R = std::fabs(s.off);
        double u = (s.Larc > 0.0) ? (sArc - s.L1)/s.Larc : 0.0;
        double a = s.a0 + u*(s.a1 - s.a0);
        xs[0] = c.Xv[0] + R*std::cos(a);
        xs[1] = c.Xv[1] + R*std::sin(a);
        *phase = 1;
    } else {
        double fwd = sArc - s.L1 - s.Larc;
        for (int d = 0; d < 2; d++) xs[d] = c.Xv[d] + fwd*c.tN[d]/ln + s.off*s.nN[d];
        *phase = 2;
    }
}

// ==================================================================== self-check
static void selfCheck() {
    std::printf("ADR-85 T1a C++ check -- LadrunoContact2DKernel vs the numpy oracle\n");

    // ---- T1: flat segment closed form ----
    {
        double X0[2]={0,0}, X1[2]={2,0}, xs[2]={0.5,-0.3}, xi, gap, n[2];
        int st = K::projectSegment2D(X0, X1, xs, 1.0, 2.0, xi, gap, n);
        check("flat segment status IN_BOUNDS", st == K::SEG2D_IN_BOUNDS);
        check("xi == 0.25", std::fabs(xi - 0.25) < 1e-15);
        check("n == (0,1)", std::fabs(n[0]) < 1e-15 && std::fabs(n[1]-1.0) < 1e-15);
        check("gap == -0.3 (penetrating)", std::fabs(gap + 0.3) < 1e-15 && gap < 0.0);
    }
    // ---- T2: out-of-bounds REFUSAL, never a clamp ----
    {
        double X0[2]={0,0}, X1[2]={2,0}, xi, gap, n[2];
        double lo[2]={-1.0,0.5}, hi[2]={3.0,0.5};
        check("xi < 0 -> SEG2D_OUT_LOW",
              K::projectSegment2D(X0,X1,lo,1.0,2.0,xi,gap,n) == K::SEG2D_OUT_LOW);
        check("OUT_LOW xi NOT clamped to [0,1]", xi < 0.0);
        check("xi > 1 -> SEG2D_OUT_HIGH",
              K::projectSegment2D(X0,X1,hi,1.0,2.0,xi,gap,n) == K::SEG2D_OUT_HIGH);
        check("OUT_HIGH xi NOT clamped to [0,1]", xi > 1.0);
    }
    // ---- T3: zero-length + bad-Lref refusal, outputs zeroed ----
    {
        double X0[2]={1,1}, X1[2]={1,1}, xs[2]={2,2}, xi=9, gap=9, n[2]={9,9};
        check("collapsed segment -> DEGENERATE",
              K::projectSegment2D(X0,X1,xs,1.0,1.0,xi,gap,n) == K::SEG2D_DEGENERATE);
        check("DEGENERATE zeroes the outputs", xi==0.0 && gap==0.0 && n[0]==0.0 && n[1]==0.0);
        double A[2]={0,0}, B[2]={1,0};
        check("Lref <= 0 -> DEGENERATE (contract violation, not a fallback)",
              K::projectSegment2D(A,B,xs,1.0,0.0,xi,gap,n) == K::SEG2D_DEGENERATE);
    }
    // ---- T4: sigma-flip antisymmetry is EXACT ----
    {
        double X0[2]={0.3,-1.2}, X1[2]={2.7,0.9}, xs[2]={1.1,0.4};
        double xiA, gA, nA[2], xiB, gB, nB[2];
        K::projectSegment2D(X0,X1,xs, 1.0,1.0,xiA,gA,nA);
        K::projectSegment2D(X0,X1,xs,-1.0,1.0,xiB,gB,nB);
        check("sigma flip: n(-s) == -n(s) EXACTLY", nB[0]==-nA[0] && nB[1]==-nA[1]);
        check("sigma flip: gap(-s) == -gap(s) EXACTLY", gB == -gA);
        check("sigma flip: xi unchanged EXACTLY", xiA == xiB);
    }
    // ---- T5: CONVEX corner -- both segments refuse, the wedge owns, gap is POSITIVE ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,-1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);        // material in the third quadrant
        double xs[2]={0.3,0.4};                        // inside the exterior wedge
        double xi,g,n[2];
        int sP = K::projectSegment2D(c.Xprev,c.Xv,xs,c.sigma,c.Lref,xi,g,n);
        int sN = K::projectSegment2D(c.Xv,c.Xnext,xs,c.sigma,c.Lref,xi,g,n);
        check("convex wedge: prev segment refuses OUT_HIGH", sP == K::SEG2D_OUT_HIGH);
        check("convex wedge: next segment refuses OUT_LOW",  sN == K::SEG2D_OUT_LOW);
        K::WedgeResult w = K::vertexWedge2D(c.tP,c.tN,c.sigma,xs,c.Xv,c.Lref);
        check("convex wedge: vertexWedge2D accepts", w.inWedge == 1);
        check("convex wedge: corner == C2D_CONVEX", w.corner == K::C2D_CONVEX);
        check("convex wedge: sideSign == -corner", w.sideSign == -(double)w.corner);
        double gv, nv[2];
        int st = K::vertexEval2D(c.Xv, xs, w.sideSign, gv, nv, 1e-8*c.Lref);
        double d = std::sqrt(0.3*0.3+0.4*0.4);
        char buf[80]; std::snprintf(buf,sizeof buf,"gap=%.6f d=%.6f", gv, d);
        check("convex wedge: gap == +||xs-Xv|| > 0 (NEVER penetrating)",
              st == K::VTX2D_OK && std::fabs(gv - d) < 1e-15 && gv > 0.0, buf);
    }
    // ---- T5b: out of bounds on BOTH segments but at the FAR ends is NOT the wedge ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,-1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);
        double xs[2]={-3.0,-3.0};
        double xi,g,n[2];
        int sP = K::projectSegment2D(c.Xprev,c.Xv,xs,c.sigma,c.Lref,xi,g,n);
        int sN = K::projectSegment2D(c.Xv,c.Xnext,xs,c.sigma,c.Lref,xi,g,n);
        K::WedgeResult w = K::vertexWedge2D(c.tP,c.tN,c.sigma,xs,c.Xv,c.Lref);
        check("far-side double refusal is NOT a wedge (patch boundary, owned by nobody)",
              sP == K::SEG2D_OUT_LOW && sN == K::SEG2D_OUT_HIGH && w.inWedge == 0);
    }
    // ---- T6: CONCAVE corner -- both segments refuse and the point IS penetrating ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);        // free space = second quadrant
        double xs[2]={0.3,-0.4};                       // behind the corner, INSIDE material
        double xi,g,n[2];
        int sP = K::projectSegment2D(c.Xprev,c.Xv,xs,c.sigma,c.Lref,xi,g,n);
        int sN = K::projectSegment2D(c.Xv,c.Xnext,xs,c.sigma,c.Lref,xi,g,n);
        check("concave wedge: prev segment refuses OUT_HIGH", sP == K::SEG2D_OUT_HIGH);
        check("concave wedge: next segment refuses OUT_LOW",  sN == K::SEG2D_OUT_LOW);
        K::WedgeResult w = K::vertexWedge2D(c.tP,c.tN,c.sigma,xs,c.Xv,c.Lref);
        check("concave wedge: vertexWedge2D accepts", w.inWedge == 1);
        check("concave wedge: corner == C2D_CONCAVE", w.corner == K::C2D_CONCAVE);
        check("concave wedge: sideSign == -1 (interior)", w.sideSign == -1.0);
        double gv, nv[2];
        int st = K::vertexEval2D(c.Xv, xs, w.sideSign, gv, nv, 1e-8*c.Lref);
        double d = std::sqrt(0.3*0.3+0.4*0.4);
        check("concave wedge: gap == -||xs-Xv|| < 0 (PENETRATING -- the real hole)",
              st == K::VTX2D_OK && std::fabs(gv + d) < 1e-15 && gv < 0.0);
        // the ejection direction must point back toward the free space (the vertex)
        check("concave wedge: n points from the slave BACK to the vertex",
              nv[0] < 0.0 && nv[1] > 0.0);
    }
    // ---- T7: slide around a CONVEX corner at constant distance ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,-1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);
        double R = 0.25;
        Sweep sw = makeSweep(c, R, 0.8);
        const int NS = 2000;
        double ds = sw.Ltot/(NS-1);
        double maxGapErr = 0.0, maxAng = 0.0, prevN[2] = {0,0};
        // handoff RANK (not the `which` code): prev segment 0, vertex 1, next segment 2
        const int rank[3] = { 0, 2, 1 };
        int badOwner = 0, seenVertex = 0, order = 0, ordBad = 0;
        for (int i = 0; i < NS; i++) {
            double xs[2]; int ph;
            sweepPoint(sw, i*ds, xs, &ph);
            Owned o = ownerOf(c, xs, +1.0);
            if (o.which < 0 || o.nOwners < 1) badOwner++;
            double e = std::fabs(o.gap - R); if (e > maxGapErr) maxGapErr = e;
            if (o.which == 2) seenVertex++;
            if (o.which >= 0) {
                int rk = rank[o.which];
                if (rk > order) order = rk;
                else if (rk < order) ordBad++;
            }
            if (i > 0) {
                double cs = o.n[0]*prevN[0] + o.n[1]*prevN[1];
                if (cs >  1.0) cs =  1.0;
                if (cs < -1.0) cs = -1.0;
                double a = std::acos(cs); if (a > maxAng) maxAng = a;
            }
            prevN[0]=o.n[0]; prevN[1]=o.n[1];
        }
        char b1[96];
        std::snprintf(b1,sizeof b1,"maxGapErr=%.3e maxAngStep=%.3e bound=%.3e",
                      maxGapErr, maxAng, ds/R);
        check("convex slide: the ordered rule resolves exactly one owner everywhere", badOwner == 0);
        check("convex slide: handoff order segment -> vertex -> segment",
              seenVertex > 0 && ordBad == 0);
        check("convex slide: gap constant (no jump) + normal rotation <= ds/R",
              maxGapErr < 1e-14 && maxAng <= ds/R*(1.0+1e-9), b1);
    }
    // ---- T8: slide around a CONCAVE corner INSIDE the material (the force-carrying case) ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);
        double R = 0.25;
        Sweep sw = makeSweep(c, -R, 0.8);
        const int NS = 2000;
        double ds = sw.Ltot/(NS-1);
        double maxGapErr = 0.0, maxAng = 0.0, prevN[2] = {0,0};
        int badOwner = 0, nVertex = 0, allNeg = 1;
        for (int i = 0; i < NS; i++) {
            double xs[2]; int ph;
            sweepPoint(sw, i*ds, xs, &ph);
            Owned o = ownerOf(c, xs, -1.0);
            if (o.which < 0 || o.nOwners < 1) badOwner++;
            double e = std::fabs(o.gap + R); if (e > maxGapErr) maxGapErr = e;
            if (o.which == 2) nVertex++;
            if (!(o.gap < 0.0)) allNeg = 0;
            if (i > 0) {
                double cs = o.n[0]*prevN[0] + o.n[1]*prevN[1];
                if (cs >  1.0) cs =  1.0;
                if (cs < -1.0) cs = -1.0;
                double a = std::acos(cs); if (a > maxAng) maxAng = a;
            }
            prevN[0]=o.n[0]; prevN[1]=o.n[1];
        }
        char b2[112];
        std::snprintf(b2,sizeof b2,"maxGapErr=%.3e maxAngStep=%.3e unownedIfNoConcaveRule=%d/%d",
                      maxGapErr, maxAng, nVertex, NS);
        check("concave slide: the ordered rule resolves exactly one owner everywhere", badOwner == 0);
        check("concave slide: gap == -R everywhere, ALL penetrating",
              maxGapErr < 1e-14 && allNeg == 1);
        check("concave slide: normal rotation <= ds/R (continuous through the wedge)",
              maxAng <= ds/R*(1.0+1e-9), b2);
        check("concave slide: the wedge arc is NON-EMPTY (policy 'concave never activates'"
              " would leave it unowned at penetration R)", nVertex > 0);
    }
    // ---- T8b: the seam band -- the naive XOR is NOT exact, the ordered rule is ----
    {
        double Xp[2]={-1,0}, Xv[2]={0,0}, Xn[2]={0,-1};
        Corner c = makeCorner(Xp, Xv, Xn, 1.0);
        int xorBad = 0, ruleBad = 0, n = 0;
        const double tol = K::TOL_IN_BOUNDS_DEFAULT;
        for (int i = -400; i <= 400; i++) {
            // slave placed at aPrev = tol + i*1e-18 (tP = (1,0) so aPrev is just x)
            double xs[2] = { tol + i*1e-18, -0.5 };
            double xiP, gP, nP[2], xiN, gN, nN[2];
            int sP = K::projectSegment2D(c.Xprev, c.Xv,    xs, 1.0, c.Lref, xiP, gP, nP);
            int sN = K::projectSegment2D(c.Xv,    c.Xnext, xs, 1.0, c.Lref, xiN, gN, nN);
            K::WedgeResult w = K::vertexWedge2D(c.tP, c.tN, 1.0, xs, c.Xv, c.Lref);
            int nx = (sP == K::SEG2D_IN_BOUNDS ? 1:0) + (sN == K::SEG2D_IN_BOUNDS ? 1:0)
                   + (w.inWedge ? 1:0);
            if (nx != 1) xorBad++;
            Owned o = ownerOf(c, xs, +1.0);
            if (o.which < 0 || o.nOwners < 1) ruleBad++;
            n++;
        }
        char b[96]; std::snprintf(b,sizeof b,"naive XOR wrong at %d/%d seam stations", xorBad, n);
        check("seam band: the naive slacked XOR is NOT floating-point exact "
              "(documented, not assumed)", xorBad > 0, b);
        check("seam band: the ORDERED rule is total and unique across the whole band",
              ruleBad == 0);
    }
    // ---- T9: B-operators are the EXACT first variation of the gap ----
    {
        double X0[2]={0.2,-0.4}, X1[2]={2.3,0.7}, xs[2]={1.0,0.05};
        double xi,g0,n[2],N[2],B[6];
        K::projectSegment2D(X0,X1,xs,1.0,1.0,xi,g0,n);
        K::shape2D(xi,N);
        K::bOperatorSegment2D(N[0],N[1],n,B);
        double *q[6] = { &xs[0],&xs[1],&X0[0],&X0[1],&X1[0],&X1[1] };
        double h = 1e-7, maxe = 0.0;
        for (int j = 0; j < 6; j++) {
            double sav = *q[j], gp, gm, xit, nt[2];
            *q[j] = sav + h; K::projectSegment2D(X0,X1,xs,1.0,1.0,xit,gp,nt);
            *q[j] = sav - h; K::projectSegment2D(X0,X1,xs,1.0,1.0,xit,gm,nt);
            *q[j] = sav;
            double fd = (gp-gm)/(2*h), e = std::fabs(fd - B[j]);
            if (e > maxe) maxe = e;
        }
        char b3[64]; std::snprintf(b3,sizeof b3,"max err=%.2e", maxe);
        check("segment B == d(gap)/du (FD)", maxe < 1e-7, b3);

        double Xv[2]={0.4,0.9}, xp[2]={1.3,-0.2}, gv, nv[2], BV[4];
        K::vertexEval2D(Xv,xp,-1.0,gv,nv,1e-9);
        K::bOperatorVertex2D(nv,BV);
        double *qv[4] = { &xp[0],&xp[1],&Xv[0],&Xv[1] };
        double maxev = 0.0;
        for (int j = 0; j < 4; j++) {
            double sav = *qv[j], gp, gm, nt[2];
            *qv[j] = sav + h; K::vertexEval2D(Xv,xp,-1.0,gp,nt,1e-9);
            *qv[j] = sav - h; K::vertexEval2D(Xv,xp,-1.0,gm,nt,1e-9);
            *qv[j] = sav;
            double fd = (gp-gm)/(2*h), e = std::fabs(fd - BV[j]);
            if (e > maxev) maxev = e;
        }
        char b4[64]; std::snprintf(b4,sizeof b4,"max err=%.2e", maxev);
        check("vertex B == d(gap)/du (FD)", maxev < 1e-7, b4);
    }
    // ---- T10: sigmaFromRef2D orients, and REFUSES a tangent datum ----
    {
        double t[2]={2.0,0.0}, up[2]={0.0,1.0}, dn[2]={0.0,-3.0}, tang[2]={5.0,0.0};
        check("sigmaFromRef2D(+y) == +1", K::sigmaFromRef2D(t,up) == 1.0);
        check("sigmaFromRef2D(-y) == -1", K::sigmaFromRef2D(t,dn) == -1.0);
        check("tangent datum -> 0 (REFUSE, never guess)", K::sigmaFromRef2D(t,tang) == 0.0);
    }
    // ---- T11: vertexEval2D coincidence refusal ----
    {
        double Xv[2]={1,1}, xs[2]={1.0+1e-12,1.0}, g, n[2];
        check("slave on the vertex -> VTX2D_DEGENERATE",
              K::vertexEval2D(Xv,xs,1.0,g,n,1e-9) == K::VTX2D_DEGENERATE);
        check("DEGENERATE zeroes the outputs", g==0.0 && n[0]==0.0 && n[1]==0.0);
    }

    std::printf("\n================================================\n%s  (T1a C++ check)\n",
                fails ? "FAILURE(S)" : "ALL PASS");
}

int main(int argc, char **argv) {
    if (argc > 1 && std::strcmp(argv[1], "--trials") == 0) return runTrials();
    selfCheck();
    return fails ? 1 : 0;
}
