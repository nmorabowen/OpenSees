// Standalone g++ gate for Ladruno WP-243 in LadrunoRCKernel.h: the BeamFiber view driver
// (driveBeamFiber), the tie law (tieUpdate) and the elastic condensed tangent
// (beamFiberInitialTangent). Header-only kernel, no OpenSees. Self-asserting: non-zero exit and
// no PASS line on any failure.
//   g++ -O2 -std=c++17 -I SRC/material/nD rc_beamfiber_gpp.cpp -o rc_bf && ./rc_bf
//
// Gates (each against a reference that does not reuse the driver's algorithm)
//  B1 elastic, no ties:  closed form. sigma_xx = E eps_xx, tau = G gamma, eps_yy = eps_zz =
//                        -nu eps_xx, retained tangent diag(E, G, G).
//  B2 elastic, ties:     closed form. rhoY = rhoZ = rho: eps_lat = -lam eps_xx / (2(lam+mu) + rho Es),
//                        E* = lam + 2 mu - 2 lam^2 / (2(lam+mu) + rho Es); initial tangent equals it.
//  B3 two-leg:           a cracking shear path (axial compression + gamma_xy, ties in y only). Each
//                        step is re-solved from the same committed state by nested bracketing
//                        bisection (no Jacobian, no Newton, root nearest the committed strain, the
//                        rcbf-0 oracle's rule); the stresses agree.
//  B4 condensed FD:      central FD of the retained map vs the condensed tangent, at an elastic, a
//                        cracked and a tie-yielded state. The kernel's own 3D tangent is only
//                        approximate in softening, so the gate is: condensation does not degrade it
//                        (eC <= 1e-4 + 2 e3D, e3D = FD error of the 6x6 at the same strain).
//  B5 invariance:        gamma_xy with rhoY vs gamma_xz with rhoZ: identical shear, swapped laterals.
//  B6 crack gate:        an unconfined fibre in axial compression to 0.9 fc: the damage gate (the
//                        view's default) never latches cracked-nu; the strain gate does (break-once).
//  B7 ties:              bilinear return: yield at fy, Et = b Es past it, elastic unloading keeps
//                        the plastic strain; in a cracked fibre the converged tie stress is fy.
//  B8 crack check:       Model C `sum` cap vs the spec formula (06 section 2.6, 2D membrane
//                        definitions written out here) at random x-y states.
//  B9 cap binds:         TS on, shear ramp under axial compression: every step converges, the cap
//                        binds with q = cap at the converged state, the shear ends lower; the FD
//                        gradient of the cap does not depend on its step.
//  B10 invariance:       gamma_xy vs gamma_xz with the crack check on.
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace ladruno_rc_kernel;
static const double PI_ = 3.14159265358979323846;

static int nfail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nfail++; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static void fill(Params& P)
{
    P.E = 30000.0; P.nu = 0.2; P.Kc = 0.667;
    P.cdf = 0.0; P.eta = 0.0; P.betaFloor = 0.1;
    P.betaOn = true; P.lublinerTCReduced = false; P.tangentMode = 0;
    P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4;
    P.aggSize = 16.0; P.crackStrain = 0.0; P.crackSpacing = 100.0; P.lch = 0.0; P.betaSrMin = 0.01;
    P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
    P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
    P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
    P.tensStiffMode = 1; P.tensStiffC = 200.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
    P.autoReg = false; P.lchRef = 1.0;
    P.addonPlane = 1; P.crackGate = 1;                 // what the BeamFiber view selects
    P.tieRhoY = 0.0; P.tieRhoZ = 0.0; P.tieFy = 420.0; P.tieEs = 200000.0; P.tieB = 0.0;
    double Ce[4] = {0.0, 0.0007, 0.0020, 0.0100};
    double Cs[4] = {0.0, 24.0,   30.0,   5.0};
    double Cd[4] = {0,0,0,0};
    double Te[3] = {0.0, 0.0001, 0.0010};
    double Ts[3] = {0.0, 3.0,    0.5};
    double Td[3] = {0,0,0};
    buildBackbone(P.hc, P.E, Ce, Cs, Cd, 4);
    buildBackbone(P.ht, P.E, Te, Ts, Td, 3);
    setupParams(P);
}

static double rel(double a, double b) { return fabs(a - b) / (fabs(a) + fabs(b) + 1e-12); }

struct Fib {                       // a fibre: committed history + tie plastic strains
    RCHist h; double ep[2];
    Fib() { histZero(h); ep[0] = ep[1] = 0.0; }
};

// one trial at retained (exx, gxy, gxz) from the committed fibre; optionally commit
static int trial(const Params& P, Fib& f, double exx, double gxy, double gxz,
                 double sig[6], double D[6][6], double eps[6], double ts[2], bool commit)
{
    for (int a = 0; a < 6; ++a) eps[a] = 0.0;
    eps[0] = exx; eps[3] = gxy; eps[5] = gxz;
    RCHist o; double ep[2];
    int st = driveBeamFiber(P, eps, f.h, f.ep, sig, D, o, ep, ts);
    if (commit && st == STATUS_OK) { f.h = o; f.ep[0] = ep[0]; f.ep[1] = ep[1]; }
    return st;
}

// lateral residual at a full strain from a committed fibre (y-normal, z-normal), gamma_yz = 0
static void resid2(const Params& P, const Fib& f, double e[6], double r[2], double s[6])
{
    double D[6][6]; RCHist o;
    returnMap3D(P, e, f.h, s, D, o);
    double Et, epn, sy = 0.0, sz = 0.0;
    if (P.tieRhoY > 0.0) sy = tieUpdate(e[1], f.ep[0], P.tieEs, P.tieFy, P.tieB, Et, epn);
    if (P.tieRhoZ > 0.0) sz = tieUpdate(e[2], f.ep[1], P.tieEs, P.tieFy, P.tieB, Et, epn);
    r[0] = s[1] + P.tieRhoY * sy;
    r[1] = s[2] + P.tieRhoZ * sz;
}

// bracket g(x) = 0 (g increasing) nearest x0 by scanning geometric offsets, then bisect
template <class G> static bool rootNear(G g, double x0, double& x)
{
    double off[41]; for (int k = 0; k < 41; ++k) off[k] = 1e-8 * pow(10.0, 7.0 * k / 40.0);
    double best = 1e300, lo = 0, hi = 0; bool found = false;
    double xs[83]; int n = 0;
    for (int k = 40; k >= 0; --k) xs[n++] = x0 - off[k];
    xs[n++] = x0;
    for (int k = 0; k < 41; ++k) xs[n++] = x0 + off[k];
    double gprev = g(xs[0]);
    for (int i = 1; i < n; ++i) {
        double gi = g(xs[i]);
        if (gprev <= 0.0 && gi >= 0.0) {
            double d = fabs(0.5 * (xs[i - 1] + xs[i]) - x0);
            if (d < best) { best = d; lo = xs[i - 1]; hi = xs[i]; found = true; }
        }
        gprev = gi;
    }
    if (!found) return false;
    for (int it = 0; it < 80; ++it) { double m = 0.5 * (lo + hi); if (g(m) <= 0.0) lo = m; else hi = m; }
    x = 0.5 * (lo + hi);
    return true;
}

// FD of the condensed retained map; returns max rel error of D_RR vs FD
static double fdCondensed(const Params& P, Fib& f, const double eR[3], double& e3D)
{
    static const int R[3] = {0, 3, 5};
    double sig[6], D[6][6], eps[6], ts[2];
    trial(P, f, eR[0], eR[1], eR[2], sig, D, eps, ts, false);
    double num = 0.0, den = 0.0;
    for (int j = 0; j < 3; ++j) {
        double h = 1e-7, ep[3] = {eR[0], eR[1], eR[2]}, em[3] = {eR[0], eR[1], eR[2]};
        ep[j] += h; em[j] -= h;
        double sp[6], sm[6], Dx[6][6], ex[6], tx[2];
        trial(P, f, ep[0], ep[1], ep[2], sp, Dx, ex, tx, false);
        trial(P, f, em[0], em[1], em[2], sm, Dx, ex, tx, false);
        for (int i = 0; i < 3; ++i) {
            double fd = (sp[R[i]] - sm[R[i]]) / (2 * h);
            num = fmax(num, fabs(fd - D[R[i]][R[j]])); den = fmax(den, fabs(fd));
        }
    }
    // the kernel's own 3D tangent at the converged full strain
    double s0[6], D3[6][6]; RCHist o;
    returnMap3D(P, eps, f.h, s0, D3, o);
    double n3 = 0.0, d3 = 0.0;
    for (int j = 0; j < 6; ++j) {
        double h = 1e-7, ep[6], em[6];
        for (int a = 0; a < 6; ++a) { ep[a] = eps[a]; em[a] = eps[a]; }
        ep[j] += h; em[j] -= h;
        double sp[6], sm[6], Dx[6][6];
        returnMap3D(P, ep, f.h, sp, Dx, o); returnMap3D(P, em, f.h, sm, Dx, o);
        for (int i = 0; i < 6; ++i) {
            double fd = (sp[i] - sm[i]) / (2 * h);
            n3 = fmax(n3, fabs(fd - D3[i][j])); d3 = fmax(d3, fabs(fd));
        }
    }
    e3D = n3 / (d3 + 1e-30);
    return num / (den + 1e-30);
}

int main()
{
    const double G = 30000.0 / (2.0 * 1.2);
    const double lam = 30000.0 * 0.2 / (1.2 * 0.6), mu = G;

    // ---------------- B1 elastic, no ties ----------------
    {
        Params P; fill(P); Fib f;
        double sig[6], D[6][6], eps[6], ts[2];
        int st = trial(P, f, -2.0e-5, 1.0e-5, -0.5e-5, sig, D, eps, ts, false);
        CHECK(st == STATUS_OK, "B1 status %d", st);
        double e = 0.0;
        e = fmax(e, rel(sig[0], 30000.0 * -2.0e-5));
        e = fmax(e, rel(sig[3], G * 1.0e-5));
        e = fmax(e, rel(sig[5], G * -0.5e-5));
        e = fmax(e, rel(eps[1], 0.2 * 2.0e-5));
        e = fmax(e, rel(eps[2], 0.2 * 2.0e-5));
        e = fmax(e, rel(D[0][0], 30000.0)); e = fmax(e, rel(D[3][3], G)); e = fmax(e, rel(D[5][5], G));
        double off = fabs(D[0][3]) + fabs(D[0][5]) + fabs(D[3][5]);
        double Di[6][6]; beamFiberInitialTangent(P, Di);
        e = fmax(e, rel(Di[0][0], 30000.0));
        printf("B1 elastic, no ties: max rel err %.2e, off-diag %.2e\n", e, off);
        CHECK(e < 1e-9 && off < 1e-6, "B1");
    }

    // ---------------- B2 elastic, ties ----------------
    {
        Params P; fill(P); P.tieRhoY = P.tieRhoZ = 0.01; Fib f;
        const double k = 2.0 * (lam + mu) + 0.01 * 200000.0;
        const double exx = -2.0e-5;
        const double elat = -lam * exx / k, Es = lam + 2 * mu - 2 * lam * lam / k;
        double sig[6], D[6][6], eps[6], ts[2];
        trial(P, f, exx, 0.0, 0.0, sig, D, eps, ts, false);
        double Di[6][6]; beamFiberInitialTangent(P, Di);
        double e = fmax(fmax(rel(eps[1], elat), rel(eps[2], elat)), fmax(rel(sig[0], Es * exx), rel(D[0][0], Es)));
        e = fmax(e, rel(Di[0][0], Es));
        e = fmax(e, rel(ts[0], 200000.0 * elat));
        printf("B2 elastic, ties: E* = %.1f (E = 30000), max rel err %.2e\n", Es, e);
        CHECK(e < 1e-9, "B2");
    }

    // ---------------- B3 two-leg ----------------
    {
        Params P; fill(P); P.tieRhoY = 0.004; Fib f;
        double worst = 0.0; int nbad = 0, nst = 0;
        for (int k = 1; k <= 60; ++k) {
            double t = k / 60.0;
            double exx = -3.0e-4 * t, gxy = 6.0e-3 * t;
            double sig[6], D[6][6], eps[6], ts[2];
            Fib f0 = f;
            int st = trial(P, f, exx, gxy, 0.0, sig, D, eps, ts, true);
            if (st != STATUS_OK) { nbad++; continue; }
            // leg 2: outer bisection on eps_yy (y residual), inner on eps_zz (z residual)
            auto inner = [&](double eyy, double& ezz, double r[2], double s[6]) {
                auto gz = [&](double z) { double e[6] = {exx, eyy, z, gxy, 0.0, 0.0}, rr[2], ss[6]; resid2(P, f0, e, rr, ss); return rr[1]; };
                bool ok = rootNear(gz, f0.h.strain[2], ezz);
                double e[6] = {exx, eyy, ezz, gxy, 0.0, 0.0}; resid2(P, f0, e, r, s);
                return ok;
            };
            auto gy = [&](double y) { double z, r[2], s[6]; inner(y, z, r, s); return r[0]; };
            double eyy, ezz, r[2], s[6];
            bool ok = rootNear(gy, f0.h.strain[1], eyy) && inner(eyy, ezz, r, s);
            if (!ok) { nbad++; continue; }
            nst++;
            double sc = fabs(sig[0]) + fabs(sig[3]) + 1e-6;
            worst = fmax(worst, fmax(fabs(s[0] - sig[0]), fabs(s[3] - sig[3])) / sc);
        }
        printf("B3 two-leg (Newton vs nested bisection): %d steps, max rel stress diff %.2e, failures %d\n", nst, worst, nbad);
        CHECK(nbad == 0 && nst == 60 && worst < 1e-6, "B3");
        CHECK(f.h.xt > backbonePeakAbscissa(P.ht), "B3 path must crack (xt %.3e)", f.h.xt);
    }

    // ---------------- B4 condensed FD ----------------
    {
        const char* name[3] = {"elastic", "cracked", "tie-yielded"};
        for (int c = 0; c < 3; ++c) {
            Params P; fill(P); P.tieRhoY = P.tieRhoZ = (c == 2) ? 0.0015 : 0.005; P.tieB = (c == 2) ? 0.01 : 0.0;
            Fib f;
            double gEnd = (c == 0) ? 2.0e-5 : (c == 1 ? 1.5e-3 : 8.0e-3);
            int n = 40;
            double sig[6], D[6][6], eps[6], ts[2];
            for (int k = 1; k <= n; ++k)
                trial(P, f, -1.0e-4 * k / n, gEnd * k / n, 0.3 * gEnd * k / n, sig, D, eps, ts, true);
            double eR[3] = {-1.0e-4 * 1.02, gEnd * 1.02, 0.3 * gEnd * 1.02};
            double e3D, eC = fdCondensed(P, f, eR, e3D);
            trial(P, f, eR[0], eR[1], eR[2], sig, D, eps, ts, false);
            printf("B4 condensed FD (%s): eC %.2e, e3D %.2e, tie stress %.1f\n", name[c], eC, e3D, ts[0]);
            CHECK(eC <= 1e-4 + 2.0 * e3D, "B4 %s", name[c]);
            if (c == 0) CHECK(eC < 1e-5, "B4 elastic must be exact");
            if (c == 2) CHECK(fabs(ts[0]) > P.tieFy, "B4 tie-yielded state must be past yield (%.1f)", ts[0]);
        }
    }

    // ---------------- B5 invariance ----------------
    {
        Params P; fill(P); P.tieRhoY = P.tieRhoZ = 0.004;
        Fib fa, fb; double worst = 0.0;
        for (int k = 1; k <= 50; ++k) {
            double g = 5.0e-3 * k / 50.0, exx = -2.0e-4 * k / 50.0;
            double sa[6], sb[6], Da[6][6], Db[6][6], ea[6], eb[6], ta[2], tb[2];
            trial(P, fa, exx, g, 0.0, sa, Da, ea, ta, true);
            trial(P, fb, exx, 0.0, g, sb, Db, eb, tb, true);
            double sc = fabs(sa[3]) + 1e-9;
            worst = fmax(worst, fabs(sa[3] - sb[5]) / sc);
            worst = fmax(worst, fabs(sa[0] - sb[0]) / (fabs(sa[0]) + 1e-9));
            worst = fmax(worst, fabs(ea[1] - eb[2]) / (fabs(ea[1]) + 1e-12));
            worst = fmax(worst, fabs(Da[3][3] - Db[5][5]) / (fabs(Da[3][3]) + 1e-9));
        }
        // the lateral balance is met to 1e-10 of the stress scale, so tiny lateral strains carry
        // ~1e-8 relative noise; 1e-6 still separates this from the 10 % membrane break (WP-241 G2)
        printf("B5 invariance (xy vs xz): max rel diff %.2e\n", worst);
        CHECK(worst < 1e-6, "B5");
    }

    // ---------------- B6 crack gate ----------------
    {
        for (int gate = 1; gate >= 0; --gate) {
            Params P; fill(P); P.crackGate = gate; P.crackedNu = 0.0;
            Fib f; double sig[6], D[6][6], eps[6], ts[2], sLast = 0.0; int st = 0;
            for (int k = 1; k <= 60 && st == STATUS_OK; ++k) {
                st = trial(P, f, -0.0012 * k / 60.0, 0.0, 0.0, sig, D, eps, ts, true);
                if (st == STATUS_OK) sLast = sig[0];
            }
            // strain gate: Poisson expansion reaches the crack strain (~0.5 fc here), the latch and
            // the TS floor open, sigma_yy is pinned above zero and no lateral balance exists. The
            // view must REFUSE there (status 2), never report a stress.
            printf("B6 %s gate, unconfined compression: last converged %.1f, cracked-nu %g, status %d\n",
                   gate ? "damage" : "strain", -sLast, f.h.nuCracked, st);
            if (gate == 1) CHECK(f.h.nuCracked == 0.0 && st == STATUS_OK && -sLast > 24.0, "B6 damage gate must reach 0.9 fc unlatched");
            else           CHECK(st == STATUS_NO_CONVERGE && -sLast < 20.0, "B6 strain gate should refuse early (break-once)");
        }
    }

    // ---------------- B7 ties ----------------
    {
        double Et, ep, s;
        const double Es = 200000.0, fy = 420.0, b = 0.02, ey = fy / Es;
        s = tieUpdate(0.5 * ey, 0.0, Es, fy, b, Et, ep);
        bool ok = rel(s, 0.5 * fy) < 1e-12 && Et == Es && ep == 0.0;
        s = tieUpdate(3.0 * ey, 0.0, Es, fy, b, Et, ep);
        ok = ok && rel(s, fy + b * Es * 2.0 * ey) < 1e-12 && rel(Et, b * Es) < 1e-12;
        double ep1 = ep, s2 = tieUpdate(3.0 * ey - 0.5 * ey, ep1, Es, fy, b, Et, ep);
        ok = ok && rel(s2, fy + b * Es * 2.0 * ey - 0.5 * fy) < 1e-12 && ep == ep1 && Et == Es;
        printf("B7 tie law: %s\n", ok ? "bilinear return OK" : "WRONG");
        CHECK(ok, "B7 tie law");
        Params P; fill(P); P.tieRhoY = 0.001; Fib f;
        double sig[6], D[6][6], eps[6], ts[2];
        for (int k = 1; k <= 60; ++k) trial(P, f, 0.0, 0.012 * k / 60.0, 0.0, sig, D, eps, ts, true);
        printf("B7 cracked fibre, rhoY fy = %.2f: tie stress %.2f, eps_yy %.2e\n", 0.001 * 420.0, ts[0], eps[1]);
        CHECK(rel(ts[0], 420.0) < 1e-9, "B7 the tie must sit at fy in a cracked fibre with rho fy < f_cr");
    }

    // ---------------- B8 crack-check cap vs the spec formula ----------------
    // Model C `sum` (spec 06 §2.6), written here from the 2D membrane definitions: phi from atan2,
    // th = |pi/2 - |phi||, s_th, w, v_ci,max, f1 cap = v_ci tan(th) + rho (fyt - f_st). States
    // in the x-y plane with eps_zz the most compressive principal, so n1 lies in x-y.
    {
        Params P; fill(P); P.tieRhoY = 0.003; P.crackCheck = 1; P.ccSmx = 150.0; P.ccSmt = 100.0;
        srand(7); double worst = 0.0; int n = 0;
        for (int t = 0; t < 200; ++t) {
            double ex = -6e-4 + 1.2e-3 * rand() / RAND_MAX, ey = 3e-3 * rand() / RAND_MAX;
            double g = -4e-3 + 8e-3 * rand() / RAND_MAX;
            double e[6] = { ex, ey, -5e-3, g, 0.0, 0.0 };
            double c0 = 0.5 * (ex + ey), R = sqrt(0.25 * (ex - ey) * (ex - ey) + 0.25 * g * g);
            double e1 = c0 + R;
            if (!(e1 > 1e-6)) continue;
            double phi = 0.5 * atan2(g, ex - ey);
            double th = fabs(0.5 * PI_ - fabs(phi));
            double sth = 1.0 / (fabs(sin(th)) / P.ccSmx + fabs(cos(th)) / P.ccSmt + 1e-12);
            double w = e1 * sth;
            double vci = 0.18 * sqrt(30.0) / (0.31 + 24.0 * w / (P.aggSize + 16.0));
            double thc = fmin(fmax(th, 1e-3), 0.5 * PI_ - 1e-3);
            double fs = fmin(fmax(P.tieEs * ey, -P.tieFy), P.tieFy);
            double ref = vci * tan(thc) + P.tieRhoY * fmax(P.tieFy - fs, 0.0);
            double P1[6], inj[6], meas[6], e1k;
            double ep0[2] = { 0.0, 0.0 };
            double cap = crackCheckCap(P, e, ep0, P1, inj, meas, &e1k);
            worst = fmax(worst, rel(cap, ref)); n++;
        }
        printf("B8 crack-check cap vs spec 2.6 (x-y plane): %d states, max rel diff %.2e\n", n, worst);
        CHECK(n >= 100, "B8 too few tensile states (%d)", n);
        CHECK(worst < 1e-10, "B8 cap differs from the spec formula (%.3e)", worst);
    }

    // ---------------- B9 crack check binds on a shear ramp ----------------
    // TS on, ties both ways, axial compression: every step converges, the cap binds on most of the
    // ramp with q = cap at the converged state, and the shear ends below the unchecked fibre's.
    // The FD gradient of the cap (h = 1e-9) agrees with a coarser step (h = 1e-7).
    {
        double tEnd[2] = { 0.0, 0.0 }; int nbind = 0, nconv = 0; double worstPin = 0.0, worstG = 0.0;
        for (int cc = 0; cc < 2; ++cc) {
            Params P; fill(P); P.tieRhoY = P.tieRhoZ = 0.002; P.tensStiffMode = 1;
            P.crackCheck = cc; P.ccSmx = 150.0; P.ccSmt = 100.0;
            Fib f; double sig[6], D[6][6], eps[6], ts[2];
            for (int k = 1; k <= 40; ++k) {
                double e6[6] = { -3e-4, f.h.strain[1], f.h.strain[2], 4e-3 * k / 40.0, f.h.strain[4], 0.0 };
                RCHist o; double ep[2]; int it = 0; bool b = false;
                int st = driveBeamFiber(P, e6, f.h, f.ep, sig, D, o, ep, ts, &it, false, 1.0, &b);
                if (st != STATUS_OK) break;
                if (cc == 1) {
                    nconv++;
                    if (b) {
                        nbind++;
                        double P1[6], inj[6], meas[6], e1k;
                        double cap = crackCheckCap(P, e6, f.ep, P1, inj, meas, &e1k);
                        double q = 0.0; for (int a = 0; a < 6; ++a) q += meas[a] * sig[a];
                        worstPin = fmax(worstPin, fabs(q - cap) / fmax(1.0, cap));
                        for (int kk = 0; kk < 6; ++kk) {
                            double gg[2];
                            for (int hi = 0; hi < 2; ++hi) {
                                double h = hi ? 1e-7 : 1e-9, ep_[6], em_[6], d1[6], d2[6], d3[6];
                                for (int a = 0; a < 6; ++a) { ep_[a] = e6[a]; em_[a] = e6[a]; }
                                ep_[kk] += h; em_[kk] -= h;
                                gg[hi] = (crackCheckCap(P, ep_, f.ep, d1, d2, d3, 0) - crackCheckCap(P, em_, f.ep, d1, d2, d3, 0)) / (2 * h);
                            }
                            worstG = fmax(worstG, fabs(gg[0] - gg[1]) / (fabs(gg[1]) + 1e-3 * P.E));
                        }
                    }
                }
                f.h = o; f.ep[0] = ep[0]; f.ep[1] = ep[1];
            }
            tEnd[cc] = sig[3];
        }
        printf("B9 crack check on a shear ramp: converged %d/40, binding %d, |q - cap| %.2e, tau_end %.3f vs %.3f unchecked, cap-gradient step check %.2e\n",
               nconv, nbind, worstPin, tEnd[1], tEnd[0], worstG);
        CHECK(nconv == 40, "B9 a step was refused (%d/40)", nconv);
        CHECK(nbind >= 20, "B9 the cap binds on too few steps (%d)", nbind);
        CHECK(worstPin < 1e-9, "B9 q != cap at a binding converged state (%.3e)", worstPin);
        CHECK(tEnd[1] < 0.95 * tEnd[0], "B9 the check does not lower the shear (%.3f vs %.3f)", tEnd[1], tEnd[0]);
        CHECK(worstG < 1e-4, "B9 cap gradient depends on the FD step (%.3e)", worstG);
    }

    // ---------------- B10 crack-check invariance ----------------
    {
        double worst = 0.0;
        Params P; fill(P); P.tieRhoY = P.tieRhoZ = 0.002; P.crackCheck = 1; P.ccSmx = 150.0; P.ccSmt = 100.0;
        Fib fa, fb; double sa[6], sb[6], D[6][6], ea[6], eb[6], ta[2], tb[2];
        for (int k = 1; k <= 30; ++k) {
            double g = 3e-3 * k / 30.0;
            int s1 = trial(P, fa, -3e-4, g, 0.0, sa, D, ea, ta, true);
            int s2 = trial(P, fb, -3e-4, 0.0, g, sb, D, eb, tb, true);
            CHECK(s1 == STATUS_OK && s2 == STATUS_OK, "B10 refused at step %d", k);
            worst = fmax(worst, rel(sa[3], sb[5])); worst = fmax(worst, rel(sa[0], sb[0]));
            worst = fmax(worst, rel(ea[1], eb[2])); worst = fmax(worst, rel(ea[2], eb[1]));
        }
        printf("B10 crack-check invariance (xy vs xz): max rel diff %.2e\n", worst);
        CHECK(worst < 1e-6, "B10");
    }

    if (nfail == 0) printf("GPP BEAMFIBER GATE: PASS\n");
    else            printf("GPP BEAMFIBER GATE: %d FAILURE(S)\n", nfail);
    return nfail == 0 ? 0 : 1;
}
