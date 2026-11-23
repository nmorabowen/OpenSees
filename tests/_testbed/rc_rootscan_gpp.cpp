// Standalone g++ gate for Ladruno WP-247 in LadrunoRCKernel.h: the BeamFiber view's root-scan fallback
// (driveBeamFiber, after a stalled Newton). Header-only kernel, no OpenSees. Self-asserting: non-zero exit and
// no PASS line on any failure.
//   g++ -O2 -std=c++17 -I SRC/material/nD rc_rootscan_gpp.cpp -o rc_rs && ./rc_rs
//
// The fibre is a lightly tied core fibre of the program's lab set (Sezen & Moehle 2002, 18 x 18 in, #3 ties at
// 12 in: f'c 21.1 MPa, rho_t 0.0017, f_yt 469 MPa, so rho_t f_yt = 0.80 MPa < f_cr = 1.52 MPa), with Model C's
// envelopes (Popovics to the Mander core peak, linear to 0.2 fp, the Collins-Mitchell tension curve), MCFT beta,
// the `sum` crack check, ties in y only. Proportional ramps of (eps_xx, gamma_xy) with commits.
//
//  R1 fold step:    a 60-step ramp to (-1.22e-4, -2.1e-4). Newton converges with no scan on steps 1-59; on the
//                   crack-onset step 60 it stalls and the scan finds the balanced state (rootScan = true): the
//                   pre-WP-247 kernel refused this trial.
//  R2 reference:    every step's converged stress and eps_yy against nested bracketing bisection from the same
//                   committed state (no Jacobian, no Newton, root nearest the committed strain).
//  R3 invariance:   gamma_xz with ties in z (the scan then runs along z) gives the same shear and swapped laterals.
//  R4 long ramp:    the same fibre to (-3e-4, -3e-3) in 120 steps (through cracking and tie yield): every step
//                   converges and matches the bisection reference.
//  R5 lateral box:  a cover fibre (no ties, f'c envelope): flexural tension, crushing past the peak, unloading.
//                   The pre-WP-247 kernel's last step ran to eps_yy ~ 3e4 (flat sigma_yy once fully cracked in y);
//                   in the box [-0.01, 0.15] every step converges on the bisection reference's root.
#include "LadrunoRCKernel.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace ladruno_rc_kernel;

static int nfail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nfail++; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static const double FC = 21.098, EC = 4700.0 * 4.593256;
// envelopes: the Mander core (fcc, ecc, e20) and the cover (f'c, 0.002, e20), e20 from Coleman-Spacone with Lp = 0.25 h
struct Env { double fp, ep, e20; };
static const Env CORE = {23.12, 0.00296, 0.02472}, COVER = {21.098, 0.002, 0.01580};
static const double RHO = 0.0017, FYT = 468.8, ES = 200000.0, SMX = 165.1, SMT = 304.8;

static double popovics(const Env& v, double c)
{
    if (c <= v.ep) { const double x = c / v.ep, r = EC / (EC - v.fp / v.ep); return v.fp * x * r / (r - 1 + pow(x, r)); }
    return fmax(v.fp - 0.8 * v.fp * (c - v.ep) / (v.e20 - v.ep), 0.2 * v.fp);
}

static void fill(Params& P, bool tiesInZ, bool cover = false)
{
    const Env& V = cover ? COVER : CORE;
    P.E = EC; P.nu = 0.0; P.Kc = 2.0 / 3.0; P.betaFloor = 0.05; P.cdf = 0.0; P.eta = 0.0;
    P.betaOn = true; P.lublinerTCReduced = false; P.tangentMode = 0;
    P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4;
    P.aggSize = 19.0; P.crackStrain = 0.0; P.crackSpacing = 0.0; P.lch = 0.0; P.betaSrMin = 0.01; P.sqrtFc = 0.0;
    P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
    P.implex = false; P.implexAlpha = 1.0; P.implexControl = false; P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
    P.tensStiffMode = 0; P.tensStiffC = 200.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
    P.autoReg = false; P.lchRef = 1.0; P.betaC = 170.0; P.crackedNu = -1.0;
    P.addonPlane = 1; P.crackGate = 1;
    P.tieRhoY = (cover || tiesInZ) ? 0.0 : RHO; P.tieRhoZ = (!cover && tiesInZ) ? RHO : 0.0; P.tieFy = FYT; P.tieEs = ES; P.tieB = 0.0;
    P.crackCheck = 1; P.ccSmx = SMX; P.ccSmt = SMT;
    // compression: Popovics to the peak, linear to 0.2 fp at e20, flat; damage for the Palermo-Vecchio offset
    std::vector<double> ce, cs, cd;
    ce.push_back(0.0);
    for (int k = 1; k <= 12; ++k) ce.push_back(V.ep * k / 12.0);
    ce.push_back(V.e20); ce.push_back(3.0 * V.e20);
    double dmax = 0.0;
    for (double c : ce) {
        const double f = popovics(V, c);
        cs.push_back(f);
        double d = 0.0;
        if (c > 0.0) {
            const double x = c / V.ep, p = fmin(V.ep * (0.166 * x * x + 0.132 * x), 0.9 * c);
            d = fmin(fmax(1.0 - f / (EC * (c - p)), 0.0), 0.99);
        }
        dmax = fmax(dmax, d); cd.push_back(dmax);
    }
    // tension: linear to e_cr, then f_cr / (1 + sqrt(500 (x - e_cr))); damage = secant to the origin
    const double fcr = 0.33 * sqrt(FC), ecr = fcr / EC;
    std::vector<double> te = {0.0, ecr}, ts = {0.0, fcr}, td = {0.0, 0.0};
    double tmax = 0.0;
    for (int k = 0; k < 14; ++k) {
        const double x = ecr + 1e-6 * pow(0.02 / 1e-6, k / 13.0), s = fcr / (1.0 + sqrt(500.0 * (x - ecr)));
        te.push_back(x); ts.push_back(s); tmax = fmax(tmax, fmin(1.0 - s / (EC * x), 0.995)); td.push_back(tmax);
    }
    buildBackbone(P.hc, EC, ce.data(), cs.data(), cd.data(), (int)ce.size());
    buildBackbone(P.ht, EC, te.data(), ts.data(), td.data(), (int)te.size());
    setupParams(P);
}

struct Fib {
    RCHist h; double ep[2];
    Fib() { histZero(h); ep[0] = ep[1] = 0.0; }
};

// lateral residual (tied normal, free normal) at a full strain from a committed fibre, gamma_yz = 0
static void resid(const Params& P, const Fib& f, double e[6], double r[2], double s[6], int tied)
{
    double D[6][6]; RCHist o; double tep[2] = {f.ep[0], f.ep[1]};
    returnMap3D(P, e, f.h, s, D, o, true, 0, false, 1.0);
    double Et, epn;
    const double fs = tieUpdate(e[1 + tied], f.ep[tied], P.tieEs, P.tieFy, P.tieB, Et, epn);
    applyCrackCheck(P, e, tep, s, D);
    r[0] = s[1 + tied] + (tied ? P.tieRhoZ : P.tieRhoY) * fs;
    r[1] = s[2 - tied];
}

template <class G> static bool rootNear(G g, double x0, double& x)
{
    double off[41];
    for (int k = 0; k < 41; ++k) off[k] = 1e-8 * pow(10.0, 7.0 * k / 40.0);
    double xs[83]; int n = 0;
    for (int k = 40; k >= 0; --k) xs[n++] = x0 - off[k];
    xs[n++] = x0;
    for (int k = 0; k < 41; ++k) xs[n++] = x0 + off[k];
    double best = 1e300, lo = 0, hi = 0; bool found = false;
    double gp = g(xs[0]);
    for (int i = 1; i < n; ++i) {
        const double gi = g(xs[i]);
        if (gp <= 0.0 && gi >= 0.0) {
            const double d = fabs(0.5 * (xs[i - 1] + xs[i]) - x0);
            if (d < best) { best = d; lo = xs[i - 1]; hi = xs[i]; found = true; }
        }
        gp = gi;
    }
    if (!found) return false;
    for (int it = 0; it < 90; ++it) { const double m = 0.5 * (lo + hi); if (g(m) <= 0.0) lo = m; else hi = m; }
    x = 0.5 * (lo + hi);
    return true;
}

struct Step { int st; bool scan; double sig[6], eps[6]; double refSig[6], refLat; bool refOk; };

// ramp: shear in gamma_xy (ties in y) or gamma_xz (ties in z); returns per-step results
struct Leg { double exx, g; int n; };

// piecewise-linear path of (eps_xx, gamma) legs from the origin, with commits
static std::vector<Step> path(const std::vector<Leg>& legs, bool xz, bool cover)
{
    Params P; fill(P, xz, cover); Fib f;
    const int gi = xz ? 5 : 3, tied = xz ? 1 : 0;
    std::vector<double> PX, PG;
    double x0 = 0.0, g0 = 0.0;
    for (const Leg& l : legs) {
        for (int j = 1; j <= l.n; ++j) { PX.push_back(x0 + (l.exx - x0) * j / l.n); PG.push_back(g0 + (l.g - g0) * j / l.n); }
        x0 = l.exx; g0 = l.g;
    }
    std::vector<Step> out;
    for (size_t k = 0; k < PX.size(); ++k) {
        const double exx = PX[k], g = PG[k];
        const Fib f0 = f;
        Step S; S.scan = false;
        for (int a = 0; a < 6; ++a) S.eps[a] = 0.0;
        S.eps[0] = exx; S.eps[gi] = g;
        double D[6][6], tieEp[2], tieSig[2]; RCHist o; int it = 0; bool cc = false;
        S.st = driveBeamFiber(P, S.eps, f0.h, f0.ep, S.sig, D, o, tieEp, tieSig, &it, false, 1.0, &cc, &S.scan);
        // reference: outer bisection on the tied normal, inner on the free normal, same committed state
        auto inner = [&](double xt, double& xf, double r[2], double s[6]) {
            auto gz = [&](double z) {
                double e[6] = {exx, 0, 0, 0, 0, 0}, rr[2], ss[6];
                e[gi] = g; e[1 + tied] = xt; e[2 - tied] = z; resid(P, f0, e, rr, ss, tied); return rr[1];
            };
            const bool ok = rootNear(gz, f0.h.strain[2 - tied], xf);
            double e[6] = {exx, 0, 0, 0, 0, 0};
            e[gi] = g; e[1 + tied] = xt; e[2 - tied] = xf; resid(P, f0, e, r, s, tied);
            return ok;
        };
        auto gy = [&](double y) { double z, r[2], s[6]; inner(y, z, r, s); return r[0]; };
        double xt = 0, xf = 0, r[2];
        S.refOk = rootNear(gy, f0.h.strain[1 + tied], xt) && inner(xt, xf, r, S.refSig);
        S.refLat = xt;
        if (S.st == STATUS_OK) { f.h = o; f.ep[0] = tieEp[0]; f.ep[1] = tieEp[1]; }
        out.push_back(S);
        if (S.st != STATUS_OK) break;
    }
    return out;
}

static std::vector<Step> ramp(double exx_end, double g_end, int nst, bool xz)
{
    return path({{exx_end, g_end, nst}}, xz, false);
}

static double relDiff(const Step& S, int gi)
{
    const double sc = fabs(S.sig[0]) + fabs(S.sig[gi]) + 1e-6;
    return fmax(fabs(S.sig[0] - S.refSig[0]), fabs(S.sig[gi] - S.refSig[gi])) / sc;
}

int main()
{
    // ---------------- R1 + R2: the fold step ----------------
    {
        const std::vector<Step> R = ramp(-1.22e-4, -2.1e-4, 60, false);
        int nok = 0, nscan = 0, firstScan = -1; double worst = 0.0, worstLat = 0.0;
        for (size_t k = 0; k < R.size(); ++k) {
            if (R[k].st == STATUS_OK) nok++;
            if (R[k].scan) { nscan++; if (firstScan < 0) firstScan = (int)k + 1; }
            if (R[k].st == STATUS_OK && R[k].refOk) {
                worst = fmax(worst, relDiff(R[k], 3));
                worstLat = fmax(worstLat, fabs(R[k].eps[1] - R[k].refLat));
            }
        }
        printf("R1 fold step: %d/60 converged, root scan on %d step(s), first at step %d\n", nok, nscan, firstScan);
        CHECK(nok == 60, "R1 all 60 steps converge (%d)", nok);
        CHECK(nscan == 1 && firstScan == 60, "R1 the scan runs on the crack-onset step only (n %d, first %d)", nscan, firstScan);
        printf("R2 reference: max rel stress diff %.2e, max |d eps_yy| %.2e\n", worst, worstLat);
        CHECK(worst < 1e-9 && worstLat < 1e-12, "R2 (%.2e, %.2e)", worst, worstLat);
    }
    // ---------------- R3: invariance ----------------
    {
        const std::vector<Step> A = ramp(-1.22e-4, -2.1e-4, 60, false), B = ramp(-1.22e-4, -2.1e-4, 60, true);
        double worst = 0.0; int nscanB = 0;
        const size_t n = A.size() < B.size() ? A.size() : B.size();
        for (size_t k = 0; k < n; ++k) {
            worst = fmax(worst, fabs(A[k].sig[3] - B[k].sig[5]) + fabs(A[k].sig[0] - B[k].sig[0]));
            worst = fmax(worst, 1e3 * (fabs(A[k].eps[1] - B[k].eps[2]) + fabs(A[k].eps[2] - B[k].eps[1])));
            if (B[k].scan) nscanB++;
        }
        printf("R3 invariance (gamma_xy, ties y vs gamma_xz, ties z): %zu steps, max diff %.2e, z scans %d\n", n, worst, nscanB);
        CHECK(n == 60 && worst < 1e-9 && nscanB == 1, "R3 (n %zu, diff %.2e, scans %d)", n, worst, nscanB);
    }
    // ---------------- R4: long ramp ----------------
    {
        const std::vector<Step> R = ramp(-3.0e-4, -3.0e-3, 120, false);
        int nok = 0, nscan = 0, nref = 0; double worst = 0.0;
        for (const Step& S : R) {
            if (S.st == STATUS_OK) nok++;
            if (S.scan) nscan++;
            if (S.st == STATUS_OK && S.refOk) { nref++; worst = fmax(worst, relDiff(S, 3)); }
        }
        printf("R4 long ramp: %d/120 converged, %d scans, %d referenced, max rel stress diff %.2e, end tau %.3f MPa\n",
               nok, nscan, nref, worst, R.back().sig[3]);
        CHECK(nok == 120 && nref == 120 && worst < 1e-6, "R4 (%d, %d, %.2e)", nok, nref, worst);
    }
    // ---------------- R5: the lateral box ----------------
    {
        // a cover fibre (no ties): flexural tension, crushing past the peak, then unloading. Without the box
        // the last step's lateral Newton ran to eps_yy ~ 3e4 with a small residual (fully cracked in y).
        const std::vector<Step> R = path({{6.4e-4, -1.2e-5, 30}, {-2.14e-3, 1.2e-5, 60}, {-6.96e-4, -1.4e-6, 30}}, false, true);
        int nok = 0, nref = 0; double worst = 0.0, latMax = 0.0, dLat = 0.0;
        for (const Step& S : R) {
            if (S.st == STATUS_OK) nok++;
            latMax = fmax(latMax, fmax(fabs(S.eps[1]), fabs(S.eps[2])));
            if (S.st == STATUS_OK && S.refOk) {
                nref++; worst = fmax(worst, relDiff(S, 3)); dLat = fmax(dLat, fabs(S.eps[1] - S.refLat));
            }
        }
        printf("R5 lateral box (cover fibre, tension-crush-unload): %d/120 converged, %d referenced, max |lateral| %.2e, "
               "max rel stress diff %.2e, max |d eps_yy| %.2e\n", nok, nref, latMax, worst, dLat);
        CHECK(nok == 120 && nref == 120 && latMax < 0.01 && worst < 1e-6 && dLat < 1e-9,
              "R5 (%d, %d, %.2e, %.2e, %.2e)", nok, nref, latMax, worst, dLat);
    }
    if (nfail == 0) printf("GPP ROOTSCAN GATE: PASS\n");
    else printf("GPP ROOTSCAN GATE: FAIL (%d)\n", nfail);
    return nfail == 0 ? 0 : 1;
}
