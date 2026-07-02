// Contact-review fix (2026-07) build-free self-check for LadrunoContactProjection.h +
// LadrunoContactBucketSort.h — the offset-coordinate / scale robustness gate.
//
//   g++ -std=c++11 -I ../../SRC/domain/contact proto_projection_offset.cpp -o proj_offset_check && ./proj_offset_check
//
// Sibling of proto_reemit_selfcheck.cpp. Guards three review findings:
//   1. BLOCKER — project()'s absolute residual tolerance (1e-12 on a length² residual)
//      was unreachable for models with coordinates far from the origin (a mm-unit
//      building at x~5e4 mm with h~500 mm facets has a residual noise floor ~3e-10):
//      every projection returned -1 and contact silently vanished. Fixed by the
//      scale-free parametric-step escape (|dξ|+|dη| < tolP). The sweeps below FAIL
//      200/200 on the pre-fix header.
//   2. LOW — the detK degeneracy guard compared length⁴ against length² (falsely
//      refused healthy micro-scale facets). Now a pure angle criterion.
//   3. MED — the bucket-grid cap arithmetic used 32-bit `long` (LLP64 Windows):
//      a diverging clipPct=0 re-emit feed could wrap the cell product and request
//      a huge grid allocation. Now per-axis pre-clamped + long long.

#include "LadrunoContactProjection.h"
#include "LadrunoContactBucketSort.h"
#include <cstdio>
#include <cmath>

namespace LP = LadrunoContactProjection;

static int fails = 0;
static void check(bool ok, const char *name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}

// deterministic LCG (reproducible across platforms; no <random>)
static unsigned long long rngState = 12345;
static double urand() {
    rngState = rngState * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rngState >> 11) / 9007199254740992.0;
}

// Build a facet of edge size h, rotated by (thZ=0.37, thX=0.21), centered at (L,L,L).
// warp > 0 adds out-of-plane node perturbation of warp*h (kills exact-binary luck).
// nps=4 quad or nps=3 tri. Returns the plane unit normal (of the UNWARPED facet).
static void makeFacet(int nps, double h, double L, double warp,
                      double X[4][3], double nrm[3]) {
    const double thZ = 0.37, thX = 0.21;
    double h2 = h / 2.0;
    double loc[4][3];
    if (nps == 4) {
        double q[4][3] = {{-h2,-h2,0},{h2,-h2,0},{h2,h2,0},{-h2,h2,0}};
        for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) loc[i][d] = q[i][d];
    } else {
        double t[3][3] = {{0,0,0},{h,0,0},{0,h,0}};
        for (int i = 0; i < 3; i++) for (int d = 0; d < 3; d++) loc[i][d] = t[i][d];
    }
    for (int i = 0; i < nps; i++) loc[i][2] += (urand() - 0.5) * 2.0 * warp * h;
    for (int i = 0; i < nps; i++) {
        double x = loc[i][0], y = loc[i][1], z = loc[i][2];
        double x1 = x*std::cos(thZ) - y*std::sin(thZ);
        double y1 = x*std::sin(thZ) + y*std::cos(thZ);
        double y2 = y1*std::cos(thX) - z*std::sin(thX);
        double z2 = y1*std::sin(thX) + z*std::cos(thX);
        X[i][0] = x1 + L; X[i][1] = y2 + L; X[i][2] = z2 + L;
    }
    // unwarped plane normal = Rz,Rx applied to +z
    nrm[0] = 0.0;
    nrm[1] = -std::sin(thX);
    nrm[2] =  std::cos(thX);
}

// --- sweep: warped facet at scale h / offset L, random slaves over the interior.
//     Every projection must converge (st==0) with a tangentially-orthogonal
//     residual: |R_alpha| / (|d||g_alpha|) < 1e-6 (the scale-free convergence
//     quality measure — holds for both the residual and the step escape).
static int sweepFail(int nps, double h, double L, int trials) {
    int bad = 0;
    for (int t = 0; t < trials; t++) {
        double X[4][3], nrm[3];
        makeFacet(nps, h, L, 0.05, X, nrm);
        double c[3] = {0,0,0};
        for (int i = 0; i < nps; i++) for (int d = 0; d < 3; d++) c[d] += X[i][d] / nps;
        // slave: centroid + 0.1h along the plane normal + interior in-plane jitter
        double e1[3], e2[3];
        for (int d = 0; d < 3; d++) { e1[d] = X[1][d]-X[0][d]; e2[d] = X[nps-1][d]-X[0][d]; }
        double j1 = (urand()-0.5)*0.3, j2 = (urand()-0.5)*0.3;
        double xs[3];
        for (int d = 0; d < 3; d++) xs[d] = c[d] + 0.1*h*nrm[d] + j1*e1[d]*0.5 + j2*e2[d]*0.5;
        double xi, eta;
        int st = LP::project(nps, X, xs, xi, eta, 1e-12, 10);
        if (st != 0) { bad++; continue; }
        // convergence quality: residual orthogonality, scale-free
        double N[4], dNxi[4], dNeta[4], xbar[3], g1[3], g2[3];
        LP::shape(nps, xi, eta, N, dNxi, dNeta);
        LP::interp(nps, N, X, xbar);
        LP::tangents(nps, dNxi, dNeta, X, g1, g2);
        double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
        double dn = LP::norm3(d) + 1e-300;
        if (std::fabs(LP::dot3(d,g1))/(dn*LP::norm3(g1)) > 1e-6 ||
            std::fabs(LP::dot3(d,g2))/(dn*LP::norm3(g2)) > 1e-6) bad++;
    }
    return bad;
}

// --- flat-facet analytic accuracy: for a FLAT (warp=0) facet the exact footpoint is
//     the plane projection of the slave; project()'s footpoint must match to 1e-7*h
//     at EVERY offset (this is the "the fix does not degrade accuracy" gate).
static int flatAccFail(int nps, double h, double L, int trials) {
    int bad = 0;
    for (int t = 0; t < trials; t++) {
        double X[4][3], nrm[3];
        makeFacet(nps, h, L, 0.0, X, nrm);
        double c[3] = {0,0,0};
        for (int i = 0; i < nps; i++) for (int d = 0; d < 3; d++) c[d] += X[i][d] / nps;
        double e1[3], e2[3];
        for (int d = 0; d < 3; d++) { e1[d] = X[1][d]-X[0][d]; e2[d] = X[nps-1][d]-X[0][d]; }
        double j1 = (urand()-0.5)*0.3, j2 = (urand()-0.5)*0.3;
        double foot[3], xs[3];
        for (int d = 0; d < 3; d++) {
            foot[d] = c[d] + j1*e1[d]*0.5 + j2*e2[d]*0.5;   // exact in-plane footpoint
            xs[d]   = foot[d] + 0.1*h*nrm[d];               // slave lifted off the plane
        }
        double xi, eta;
        if (LP::project(nps, X, xs, xi, eta, 1e-12, 10) != 0) { bad++; continue; }
        double N[4], dNxi[4], dNeta[4], xbar[3];
        LP::shape(nps, xi, eta, N, dNxi, dNeta);
        LP::interp(nps, N, X, xbar);
        double err = 0.0;
        for (int d = 0; d < 3; d++) err += (xbar[d]-foot[d])*(xbar[d]-foot[d]);
        if (std::sqrt(err) > 1e-7 * h) bad++;
    }
    return bad;
}

int main()
{
    // ---------------- 1. the BLOCKER sweeps (pre-fix: 200/200 FAIL at offset) ----------------
    check(sweepFail(4, 1.0,   0.0, 200) == 0, "quad h=1 @ origin (pre-fix baseline holds)");
    check(sweepFail(4, 1.0,   1e3, 200) == 0, "quad h=1 @ L=1e3");
    check(sweepFail(4, 1.0,   1e5, 200) == 0, "quad h=1 @ L=1e5           (pre-fix: 189/200 dead)");
    check(sweepFail(4, 1.0,   1e6, 200) == 0, "quad h=1 @ L=1e6           (pre-fix: 200/200 dead)");
    check(sweepFail(4, 500.0, 5e3, 200) == 0, "quad h=500 @ L=5e3 (mm bld) (pre-fix: 200/200 dead)");
    check(sweepFail(4, 500.0, 5e4, 200) == 0, "quad h=500 @ L=5e4 (mm bld) (pre-fix: 200/200 dead)");
    check(sweepFail(4, 0.5,   30.0, 200) == 0, "quad h=0.5 @ L=30 (m bld)");
    check(sweepFail(3, 500.0, 5e4, 200) == 0, "tri  h=500 @ L=5e4 (mm bld)");
    check(sweepFail(3, 1.0,   1e6, 200) == 0, "tri  h=1 @ L=1e6");

    // ---------------- 2. footpoint accuracy is scale/offset-independent ----------------
    check(flatAccFail(4, 1.0,   0.0, 50) == 0, "flat quad footpoint 1e-7*h @ origin");
    check(flatAccFail(4, 500.0, 5e4, 50) == 0, "flat quad footpoint 1e-7*h @ mm-building coords");
    check(flatAccFail(3, 500.0, 5e4, 50) == 0, "flat tri  footpoint 1e-7*h @ mm-building coords");
    check(flatAccFail(4, 1.0,   1e6, 50) == 0, "flat quad footpoint 1e-7*h @ L=1e6");

    // ---------------- 3. out-of-bounds classification preserved at offset ----------------
    {
        double X[4][3], nrm[3];
        rngState = 999;
        makeFacet(4, 500.0, 5e4, 0.0, X, nrm);
        double c[3] = {0,0,0};
        for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) c[d] += 0.25*X[i][d];
        double e1[3]; for (int d = 0; d < 3; d++) e1[d] = X[1][d]-X[0][d];
        double xs[3]; for (int d = 0; d < 3; d++) xs[d] = c[d] + 1.2*e1[d] + 50.0*nrm[d];
        double xi, eta;
        check(LP::project(4, X, xs, xi, eta, 1e-12, 10) == 1,
              "slave past the edge -> converged OUT-OF-BOUNDS (st=1) at offset");
    }

    // ---------------- 4. detK guard: angle criterion, scale-free ----------------
    // NOTE (documented behavior, pre-existing and unchanged): at micro scales the
    // ABSOLUTE residual test passes at the initial guess (R ~ |d||g| < 1e-12), so a
    // micro-facet — healthy OR collinear — "converges" to the centre footpoint
    // before the detK guard is ever evaluated. Gap on a flat facet is footpoint-
    // independent, so contact stays functional (just parametrically loose). The
    // guard is exercised only when the residual is live: an OFF-AXIS slave whose
    // in-plane offset drives R above tolR. That is what we test, at both scales.
    {
        double Xc[4][3] = {{0,0,0},{1,0,0},{2,0,0},{3,0,0}};   // collinear "quad"
        double xs[3] = {1.7, 1.0, 0.0};                        // off-axis: R0 ~ 0.3 > tolR
        double xi, eta;
        check(LP::project(4, Xc, xs, xi, eta, 1e-12, 10) == -1,
              "collinear quad + live residual refused (h~1)");
        // same geometry scaled to h~1e5 (large scale: detK cancellation noise is where
        // the OLD length²-floor sat BELOW the noise; the angle form stays exact-zero here)
        double Xb[4][3]; for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) Xb[i][d] = Xc[i][d]*1e5;
        double xsb[3] = {1.7e5, 1e5, 0.0};
        check(LP::project(4, Xb, xsb, xi, eta, 1e-12, 10) == -1,
              "collinear quad + live residual refused (h~1e5)");
        // healthy skewed facet at the SAME large scale must NOT be refused
        int bad = sweepFail(4, 1e5, 0.0, 50);
        check(bad == 0, "healthy skewed facet h=1e5 projects (guard is angle-only)");
    }

    // ---------------- 5. evalSegment end-to-end at mm-building coords ----------------
    {
        double X[4][3], nrm[3];
        rngState = 4242;
        makeFacet(4, 500.0, 5e4, 0.0, X, nrm);
        double c[3] = {0,0,0};
        for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) c[d] += 0.25*X[i][d];
        double pen = 0.5;                                   // 0.5 mm penetration
        double xs[3], orientDir[3];
        for (int d = 0; d < 3; d++) { xs[d] = c[d] - pen*nrm[d]; orientDir[d] = nrm[d]; }
        double gap, n[3], N[4];
        bool act = LP::evalSegment(4, X, xs, orientDir, gap, n, N);
        check(act, "evalSegment ACTIVE at mm-building coords (pre-fix: dead pair)");
        check(act && std::fabs(gap + pen) < 1e-6 * 500.0,
              "evalSegment gap == -0.5mm to 1e-6*h at mm-building coords");
    }

    // ---------------- 6. bucket grid survives a diverging clipPct=0 feed ----------------
    {
        // 8 unit segments + ONE segment blown to 3e8 in x AND y. 3e8 targets the exact
        // pre-fix 32-bit wrap window: per-axis cell counts ~2e8 each fit an int, but
        // the (long)NX*NY*NZ product (4e16) wraps LLP64's 32-bit long — the cap loop
        // could exit early and grid_.assign() request a ruinous allocation. (A 1e12
        // divergence instead lands in the benign per-axis double->int saturation
        // accident.) Per-axis clamp + long long keep the grid at <= cap cells.
        const int nSeg = 9, nps = 4;
        double segs[nSeg * nps * 3];
        for (int s = 0; s < 8; s++) {
            double x0 = (double)(s % 4), y0 = (double)(s / 4);
            double q[4][2] = {{x0,y0},{x0+1,y0},{x0+1,y0+1},{x0,y0+1}};
            for (int i = 0; i < 4; i++) {
                segs[(s*nps + i)*3 + 0] = q[i][0];
                segs[(s*nps + i)*3 + 1] = q[i][1];
                segs[(s*nps + i)*3 + 2] = 0.0;
            }
        }
        for (int i = 0; i < 4; i++) {                        // the runaway segment
            segs[(8*nps + i)*3 + 0] = 3e8 + (double)i;
            segs[(8*nps + i)*3 + 1] = 3e8 + (double)i;
            segs[(8*nps + i)*3 + 2] = 0.0;
        }
        LadrunoContactBucketSort::Grid bs(nSeg, nps, segs, 1.0, /*clipPct=*/0.0);
        // crash-gated by construction: the real assertion is surviving the ctor above —
        // the pre-fix header terminates with bad_alloc there (verified on MinGW LLP64)
        check(true, "bucket build with 2-axis 3e8 divergence + clipPct=0 did not blow up");
        // near-pair discovery still works on the sane part of the grid
        double slave[3] = {0.5, 0.5, 0.05};
        int cand[nSeg];
        int nc = bs.candidates(slave, cand);
        bool found = false;
        for (int k = 0; k < nc; k++) if (cand[k] == 0) found = true;
        check(found, "bucket candidates still find the true near segment");
    }

    std::printf("\n%s (%d failure%s)\n", fails == 0 ? "ALL PASS" : "FAILURES",
                fails, fails == 1 ? "" : "s");
    return fails == 0 ? 0 : 1;
}
