// Standalone g++ gate for Ladruno WP-245 (program #1131, rcbf-3): the shear map of
// LadrunoNDFiberSection3d (SRC/material/section/LadrunoFiberShearKernel.h), driven with the
// header-only RC kernel's BeamFiber view as the nonlinear fibre. No OpenSees. Self-asserting:
// non-zero exit and no PASS line on any failure.
//   g++ -O2 -std=c++17 -I SRC/material/section -I SRC/material/nD rc_section_gpp.cpp -o rc_sec
//
// Section: Model C's layout (06 §1.1), b x h rectangle, n equal strips through the depth (y), fibre
// at each strip centroid. References never reuse the kernel's maps.
//
//  S1 elastic closed forms (fibre D = diag(E, G, G)), n strips, midpoint sums written out:
//       K00 = E A, K11 = E b h^3/12 (1 - 1/n^2),
//       parabolic modelC     K33 = G A (1 + 1/(2 n^2))          (-> G A, k = 1)
//       parabolic conjugate  K33 = G A 2.25 (1 - 2 mu2 + mu4),  mu2 = 1/3 - 1/(3n^2),
//                            mu4 = 1/5 - 2/(3 n^2) + 7/(15 n^4)  (-> 6/5 G A)
//       uniform alpha        K33 = alpha G A
//  S2 convergence: conjugate K33/(G A) -> 6/5 as n grows (n = 40, 400).
//  S3 symmetry: random SPD fibre D. conjugate K is symmetric; modelC K is not (by design).
//  S4 uniform = stock: the unrolled stock NDFiberSection3d formulas (copied here) vs the kernel,
//     alpha = 0.83, random D, sigma.
//  S5 chain rule: a smooth coupled toy fibre law with an exact Jacobian; central FD of s(e) vs K
//     for both resultants (tests that a feeds the strain and b the force).
//  S6 Model C assembly: RC kernel fibres, cracking shear path under axial compression with
//     commits; section (N, M, V) vs sum over strips of each strip fibre driven DIRECTLY at
//     (eps0 - y kappa, phi(y) gamma) with V = sum tau A (06 §1.4).
//  S7 invariance: square section, parabolic in y and z; a gamma_y path vs the same gamma_z path.
//  S8 band check: a fibre at the face is inside; one past it is refused.
#include "LadrunoFiberShearKernel.h"
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace ladruno_rc_kernel;
namespace sh = ladruno_shear;

static int nfail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nfail++; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static double rel(double a, double b) { return fabs(a - b) / (fabs(a) + fabs(b) + 1e-300); }

struct Strip { double y, z, A; };

static std::vector<Strip> rect(double b, double h, int n)
{
    std::vector<Strip> s(n);
    const double dy = h / n;
    for (int i = 0; i < n; ++i) { s[i].y = -0.5 * h + dy * (i + 0.5); s[i].z = 0.0; s[i].A = b * dy; }
    return s;
}

// section with a constant fibre law: s, K at deformation e (sigma = D a e for linear laws)
static void sectionLinear(const sh::Spec& S, const std::vector<Strip>& st, const double D[3][3],
                          double s[6], double K[6][6], const double* e)
{
    for (int i = 0; i < 6; ++i) { s[i] = 0.0; for (int j = 0; j < 6; ++j) K[i][j] = 0.0; }
    for (const Strip& f : st) {
        double a[3][6], b[3][6], eps[3] = {0, 0, 0}, sig[3] = {0, 0, 0};
        sh::maps(S, f.y, f.z, a, b);
        if (e) {
            sh::strain(a, e, eps);
            for (int r = 0; r < 3; ++r) { sig[r] = 0.0; for (int c = 0; c < 3; ++c) sig[r] += D[r][c] * eps[c]; }
        }
        sh::accumulate(f.A, a, b, sig, D, s, K);
    }
}

static double urand() { return (double)rand() / RAND_MAX * 2.0 - 1.0; }

static void fillRC(Params& P)
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
    P.addonPlane = 1; P.crackGate = 1;
    P.tieRhoY = 0.002; P.tieRhoZ = 0.002; P.tieFy = 420.0; P.tieEs = 200000.0; P.tieB = 0.0;
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

struct Fib { RCHist h; double ep[2]; Fib() { histZero(h); ep[0] = ep[1] = 0.0; } };

// RC fibre at retained (exx, gxy, gxz); sig3/D3 in the order-3 BeamFiber layout
static int fibre(const Params& P, Fib& f, const double eps3[3], double sig3[3], double D3[3][3], bool commit)
{
    double e6[6] = {eps3[0], 0, 0, eps3[1], 0, eps3[2]}, s6[6], D6[6][6], ep[2], ts[2];
    RCHist o;
    int st = driveBeamFiber(P, e6, f.h, f.ep, s6, D6, o, ep, ts);
    const int R[3] = {0, 3, 5};
    for (int r = 0; r < 3; ++r) { sig3[r] = s6[R[r]]; for (int c = 0; c < 3; ++c) D3[r][c] = D6[R[r]][R[c]]; }
    if (commit && st == STATUS_OK) { f.h = o; f.ep[0] = ep[0]; f.ep[1] = ep[1]; }
    return st;
}

// a nonlinear section of RC fibres
static int sectionRC(const Params& P, const sh::Spec& S, const std::vector<Strip>& st, std::vector<Fib>& F,
                     const double e[6], double s[6], double K[6][6], bool commit)
{
    for (int i = 0; i < 6; ++i) { s[i] = 0.0; for (int j = 0; j < 6; ++j) K[i][j] = 0.0; }
    for (size_t k = 0; k < st.size(); ++k) {
        double a[3][6], b[3][6], eps[3], sig[3], D[3][3];
        sh::maps(S, st[k].y, st[k].z, a, b);
        sh::strain(a, e, eps);
        if (fibre(P, F[k], eps, sig, D, commit) != STATUS_OK) return -1;
        sh::accumulate(st[k].A, a, b, sig, D, s, K);
    }
    return 0;
}

int main()
{
    const double E = 30000.0, nu = 0.2, G = E / (2.0 * (1.0 + nu));
    const double b = 457.2, h = 457.2;
    const double Del[3][3] = {{E, 0, 0}, {0, G, 0}, {0, 0, G}};

    // ---- S1 / S2 elastic closed forms
    {
        double worst = 0.0;
        double k33m = 0, k33c = 0;
        for (int n : {40, 400}) {
            std::vector<Strip> st = rect(b, h, n);
            const double A = b * h;
            sh::Spec U; U.profile = sh::UNIFORM; U.alpha = 1.0;
            sh::Spec Pm; Pm.profile = sh::PARABOLIC; Pm.resultant = sh::MODELC; Pm.hy = h; Pm.hz = b;
            sh::Spec Pc = Pm; Pc.resultant = sh::CONJUGATE;
            sh::Spec Ua; Ua.profile = sh::UNIFORM; Ua.alpha = 0.83;
            double s[6], K[6][6];
            const double mu2 = 1.0 / 3.0 - 1.0 / (3.0 * n * n);
            const double mu4 = 0.2 - 2.0 / (3.0 * n * n) + 7.0 / (15.0 * pow((double)n, 4));
            sectionLinear(U, st, Del, s, K, 0);
            worst = fmax(worst, rel(K[0][0], E * A));
            worst = fmax(worst, rel(K[1][1], E * b * h * h * h / 12.0 * (1.0 - 1.0 / (n * (double)n))));
            worst = fmax(worst, rel(K[3][3], G * A));
            sectionLinear(Ua, st, Del, s, K, 0);
            worst = fmax(worst, rel(K[3][3], 0.83 * G * A));
            sectionLinear(Pm, st, Del, s, K, 0);
            worst = fmax(worst, rel(K[3][3], G * A * (1.0 + 1.0 / (2.0 * n * n))));
            if (n == 40) k33m = K[3][3] / (G * A);
            worst = fmax(worst, rel(K[0][0], E * A));
            sectionLinear(Pc, st, Del, s, K, 0);
            worst = fmax(worst, rel(K[3][3], G * A * 2.25 * (1.0 - 2.0 * mu2 + mu4)));
            double kc = K[3][3] / (G * A);
            if (n == 40) k33c = kc;
            CHECK(fabs(kc - 1.2) < (n == 40 ? 1e-3 : 1e-5), "S2 conjugate k at n=%d: %.8f vs 6/5", n, kc);
        }
        CHECK(worst < 1e-12, "S1 closed forms: %.2e", worst);
        printf("S1 elastic closed forms: max rel diff %.2e; n=40: modelC K_VV/GA = %.6f (1 + 1/3200), conjugate %.6f\n",
               worst, k33m, k33c);
        printf("S2 conjugate k -> 6/5: n=40 %.6f\n", k33c);
    }

    // ---- S3 symmetry with a random SPD fibre D
    {
        srand(7);
        double L[3][3], D[3][3];
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) L[r][c] = (c <= r) ? urand() + (r == c ? 2.0 : 0.0) : 0.0;
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) { D[r][c] = 0.0; for (int k = 0; k < 3; ++k) D[r][c] += L[r][k] * L[c][k]; }
        std::vector<Strip> st = rect(b, h, 40);
        for (Strip& f : st) f.z = 0.3 * f.y - 20.0 * urand();
        sh::Spec Pc; Pc.profile = sh::PARABOLIC; Pc.resultant = sh::CONJUGATE; Pc.hy = h; Pc.hz = 2 * b;
        sh::Spec Pm = Pc; Pm.resultant = sh::MODELC;
        double s[6], Kc[6][6], Km[6][6], asymC = 0, asymM = 0, nrm = 0;
        sectionLinear(Pc, st, D, s, Kc, 0);
        sectionLinear(Pm, st, D, s, Km, 0);
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            asymC = fmax(asymC, fabs(Kc[i][j] - Kc[j][i])); asymM = fmax(asymM, fabs(Km[i][j] - Km[j][i]));
            nrm = fmax(nrm, fabs(Kc[i][j]));
        }
        CHECK(asymC / nrm < 1e-14, "S3 conjugate asymmetry %.2e", asymC / nrm);
        CHECK(asymM / nrm > 1e-6, "S3 modelC should be unsymmetric, got %.2e", asymM / nrm);
        printf("S3 symmetry: conjugate %.2e, modelC %.2e (unsymmetric by design)\n", asymC / nrm, asymM / nrm);
    }

    // ---- S4 uniform = stock NDFiberSection3d (unrolled formulas copied from the stock source)
    {
        srand(11);
        const double alpha = 0.83, ra = sqrt(alpha);
        double worst = 0.0;
        std::vector<Strip> st = rect(b, h, 40);
        for (Strip& f : st) f.z = 50.0 * urand();
        double Ks[6][6] = {{0}}, Ss[6] = {0}, Kk[6][6] = {{0}}, Sk[6] = {0};
        sh::Spec U; U.profile = sh::UNIFORM; U.alpha = alpha;
        for (const Strip& f : st) {
            double D[3][3], sig[3];
            for (int r = 0; r < 3; ++r) { sig[r] = urand(); for (int c = 0; c < 3; ++c) D[r][c] = urand(); }
            double a[3][6], bb[3][6];
            sh::maps(U, f.y, f.z, a, bb);
            sh::accumulate(f.A, a, bb, sig, D, Sk, Kk);
            const double y = f.y, z = f.z, A = f.A, y2 = y*y, z2 = z*z, yz = y*z;
            double d00 = D[0][0]*A, d01 = D[0][1]*A, d02 = D[0][2]*A, d10 = D[1][0]*A, d11 = D[1][1]*A,
                   d12 = D[1][2]*A, d20 = D[2][0]*A, d21 = D[2][1]*A, d22 = D[2][2]*A, tmp;
            Ks[0][0] += d00; Ks[1][1] += y2*d00; Ks[2][2] += z2*d00;
            tmp = -y*d00; Ks[0][1] += tmp; Ks[1][0] += tmp;
            tmp = z*d00;  Ks[0][2] += tmp; Ks[2][0] += tmp;
            tmp = -yz*d00; Ks[1][2] += tmp; Ks[2][1] += tmp;
            Ks[3][3] += alpha*d11; Ks[3][4] += alpha*d12; Ks[4][3] += alpha*d21; Ks[4][4] += alpha*d22;
            Ks[5][5] += z2*d11 - yz*(d12+d21) + y2*d22;
            tmp = -z*d01 + y*d02; Ks[0][5] += tmp; Ks[1][5] -= y*tmp; Ks[2][5] += z*tmp;
            tmp = -z*d10 + y*d20; Ks[5][0] += tmp; Ks[5][1] -= y*tmp; Ks[5][2] += z*tmp;
            d01 *= ra; d02 *= ra; d10 *= ra; d11 *= ra; d12 *= ra; d20 *= ra; d21 *= ra; d22 *= ra;
            Ks[0][3] += d01; Ks[0][4] += d02; Ks[1][3] -= y*d01; Ks[1][4] -= y*d02; Ks[2][3] += z*d01; Ks[2][4] += z*d02;
            Ks[3][0] += d10; Ks[4][0] += d20; Ks[3][1] -= y*d10; Ks[4][1] -= y*d20; Ks[3][2] += z*d10; Ks[4][2] += z*d20;
            double yy2 = y*d22, zz2 = -z*d11;
            Ks[5][3] += zz2 + y*d21; Ks[5][4] += -z*d12 + yy2; Ks[3][5] += zz2 + y*d12; Ks[4][5] += -z*d21 + yy2;
            double s0 = sig[0]*A, s1 = sig[1]*A, s2 = sig[2]*A;
            Ss[0] += s0; Ss[1] += -y*s0; Ss[2] += z*s0; Ss[3] += ra*s1; Ss[4] += ra*s2; Ss[5] += -z*s1 + y*s2;
        }
        double nk = 0, ns = 0;
        for (int i = 0; i < 6; ++i) { ns = fmax(ns, fabs(Ss[i])); for (int j = 0; j < 6; ++j) nk = fmax(nk, fabs(Ks[i][j])); }
        for (int i = 0; i < 6; ++i) {
            worst = fmax(worst, fabs(Sk[i] - Ss[i]) / ns);
            for (int j = 0; j < 6; ++j) worst = fmax(worst, fabs(Kk[i][j] - Ks[i][j]) / nk);
        }
        CHECK(worst < 1e-13, "S4 uniform vs stock %.2e", worst);
        printf("S4 uniform = stock NDFiberSection3d (alpha 0.83): max rel diff %.2e\n", worst);
    }

    // ---- S5 chain rule with a smooth coupled toy law
    {
        // sigma0 = E e0 + c e0 g1^2 + d g2^2 ; sigma1 = G g1 (1 + p e0) + q g2 ; sigma2 = G g2 + r e0 g1
        const double c = 3e6, d = 2e5, p = -40.0, q = 900.0, r = 5e5;
        auto law = [&](const double e3[3], double s3[3], double D3[3][3]) {
            const double e0 = e3[0], g1 = e3[1], g2 = e3[2];
            s3[0] = E * e0 + c * e0 * g1 * g1 + d * g2 * g2;
            s3[1] = G * g1 * (1.0 + p * e0) + q * g2;
            s3[2] = G * g2 + r * e0 * g1;
            D3[0][0] = E + c * g1 * g1; D3[0][1] = 2 * c * e0 * g1; D3[0][2] = 2 * d * g2;
            D3[1][0] = G * g1 * p;      D3[1][1] = G * (1.0 + p * e0); D3[1][2] = q;
            D3[2][0] = r * g1;          D3[2][1] = r * e0;             D3[2][2] = G;
        };
        std::vector<Strip> st = rect(b, h, 40);
        for (size_t k = 0; k < st.size(); ++k) st[k].z = 60.0 * sin(0.7 * k);
        auto sec = [&](const sh::Spec& S, const double e[6], double s[6], double K[6][6]) {
            for (int i = 0; i < 6; ++i) { s[i] = 0.0; for (int j = 0; j < 6; ++j) K[i][j] = 0.0; }
            for (const Strip& f : st) {
                double a[3][6], bb[3][6], eps[3], sig[3], D[3][3];
                sh::maps(S, f.y, f.z, a, bb); sh::strain(a, e, eps); law(eps, sig, D);
                sh::accumulate(f.A, a, bb, sig, D, s, K);
            }
        };
        const double e0[6] = {-4e-4, 2e-6, -1.5e-6, 3e-3, -2e-3, 4e-6};
        const double step[6] = {1e-7, 1e-9, 1e-9, 1e-6, 1e-6, 1e-9};
        double worstAll = 0.0;
        for (int res = 0; res < 2; ++res) {
            sh::Spec S; S.profile = sh::PARABOLIC; S.resultant = res; S.hy = h; S.hz = 200.0;
            double s[6], K[6][6], sp[6], sm[6], Kd[6][6], worst = 0.0, nrm = 0.0;
            sec(S, e0, s, K);
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) nrm = fmax(nrm, fabs(K[i][j]) * step[j]);
            for (int j = 0; j < 6; ++j) {
                double ep[6], em[6];
                for (int i = 0; i < 6; ++i) { ep[i] = e0[i]; em[i] = e0[i]; }
                ep[j] += step[j]; em[j] -= step[j];
                sec(S, ep, sp, Kd); sec(S, em, sm, Kd);
                for (int i = 0; i < 6; ++i) worst = fmax(worst, fabs((sp[i] - sm[i]) / 2.0 - K[i][j] * step[j]));
            }
            worst /= nrm;
            worstAll = fmax(worstAll, worst);
            CHECK(worst < 1e-7, "S5 FD (%s) %.2e", res ? "conjugate" : "modelC", worst);
        }
        printf("S5 chain rule, central FD vs K (both resultants): %.2e\n", worstAll);
    }

    // ---- S6 Model C assembly with RC kernel fibres on a cracking shear path
    {
        Params P; fillRC(P);
        const int n = 40;
        std::vector<Strip> st = rect(b, h, n);
        sh::Spec S; S.profile = sh::PARABOLIC; S.resultant = sh::MODELC; S.hy = h; S.hz = b;
        std::vector<Fib> F(n), Fr(n);
        const double eps0 = -3e-4, kap = 1.0e-6;
        double worst = 0.0, Vpeak = 0.0, Vend = 0.0;
        int ok = 0;
        for (int k = 1; k <= 40; ++k) {
            const double gam = 4e-3 * k / 40.0;
            double e[6] = {eps0, kap, 0.0, gam, 0.0, 0.0}, s[6], K[6][6];
            if (sectionRC(P, S, st, F, e, s, K, true) != 0) break;
            // reference: Model C's per-strip formulas (06 §1.2, §1.4), strips driven directly
            double N = 0, M = 0, V = 0;
            for (int i = 0; i < n; ++i) {
                const double yi = st[i].y, u = 2.0 * yi / h, g = 1.5 * (1.0 - u * u);
                double eps3[3] = {eps0 - yi * kap, gam * g, 0.0}, sig[3], D[3][3];
                if (fibre(P, Fr[i], eps3, sig, D, true) != STATUS_OK) { worst = 1.0; break; }
                N += sig[0] * st[i].A; M += -sig[0] * st[i].A * yi; V += sig[1] * st[i].A;
            }
            worst = fmax(worst, fmax(rel(s[0], N), fmax(rel(s[1], M), rel(s[3], V))));
            Vpeak = fmax(Vpeak, s[3]); Vend = s[3];
            ++ok;
        }
        CHECK(ok == 40, "S6 converged %d/40", ok);
        CHECK(worst < 1e-12, "S6 assembly vs Model C formulas %.2e", worst);
        printf("S6 Model C assembly (RC fibres, shear ramp, commits): %d/40 steps, max rel diff %.2e, V peak %.1f kN, end %.1f kN\n",
               ok, worst, Vpeak / 1e3, Vend / 1e3);
    }

    // ---- S7 invariance: square section, gamma_y path vs gamma_z path
    {
        Params P; fillRC(P);
        const int m = 10;
        std::vector<Strip> st;
        for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j) {
            Strip f; f.y = -0.5 * h + h / m * (i + 0.5); f.z = -0.5 * h + h / m * (j + 0.5); f.A = h * h / (m * m);
            st.push_back(f);
        }
        sh::Spec S; S.profile = sh::PARABOLIC; S.resultant = sh::MODELC; S.hy = h; S.hz = h;
        std::vector<Fib> Fy(st.size()), Fz(st.size());
        double worst = 0.0;
        for (int k = 1; k <= 30; ++k) {
            const double gam = 3e-3 * k / 30.0;
            double ey[6] = {-3e-4, 0, 0, gam, 0, 0}, ez[6] = {-3e-4, 0, 0, 0, gam, 0}, sy[6], sz[6], Ky[6][6], Kz[6][6];
            if (sectionRC(P, S, st, Fy, ey, sy, Ky, true) != 0 || sectionRC(P, S, st, Fz, ez, sz, Kz, true) != 0) { worst = 1.0; break; }
            worst = fmax(worst, fmax(rel(sy[3], sz[4]), rel(sy[0], sz[0])));
            worst = fmax(worst, rel(Ky[3][3], Kz[4][4]));
        }
        CHECK(worst < 1e-8, "S7 invariance %.2e", worst);
        printf("S7 invariance (square, Vy path vs Vz path): max rel diff %.2e\n", worst);
    }

    // ---- S8 band check
    {
        sh::Spec S; S.profile = sh::PARABOLIC; S.hy = 400.0; S.hz = 300.0;
        bool okFace = sh::inside(S, 200.0, -150.0), okPast = !sh::inside(S, 200.1, 0.0) && !sh::inside(S, 0.0, -150.1);
        sh::Spec U; U.profile = sh::UNIFORM;
        CHECK(okFace && okPast && sh::inside(U, 1e9, 1e9), "S8 band check");
        printf("S8 band check: face inside %s, past the face refused %s\n", okFace ? "yes" : "NO", okPast ? "yes" : "NO");
    }

    if (nfail == 0) printf("GPP SECTION GATE: PASS\n");
    else            printf("GPP SECTION GATE: %d FAILURE(S)\n", nfail);
    return nfail == 0 ? 0 : 1;
}
