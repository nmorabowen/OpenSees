// WP-263: dual-number tangent oracle for the LadrunoRCConcrete return map (LadrunoRCKernel.h).
//
// Header-only, OpenSees-free. For every step of a set of strain PATHS (committed history threaded
// step to step, as the material does), the trial state is evaluated three ways:
//
//   * mode 0  -- returnMap3D, tangentMode 0: the shipped algorithmic tangent (the secant W_B C0
//                plus the beta / tension-stiffening / interlock cross-terms);
//   * mode 3  -- returnMap3D, tangentMode 3 (-dualTangent): one Dual<6> pass of returnMap3DT,
//                d(sigma)/d(eps_eng) of the executed map;
//   * FD      -- a central finite difference of the map itself (returnMap3D, do_tangent = false,
//                the same committed state, the same do_implex), at seven steps h.
//
// and prints, per step:
//   FD    min over h of ||C_dual - C_fd(h)||_F / ||C_dual||_F      (the oracle)
//   PLAT  how many of the seven h are within 10x of that minimum   (a step plateau, not one lucky h)
//   OMIT  ||C_dual - C_mode0||_F / ||C_dual||_F                     (what mode 0 leaves out)
//   SYM   ||sym(C_dual) - C_fd||_F / ||C_dual||_F                   (mutant: must be flagged where the
//                                                                    damage rate makes C unsymmetric)
//   BITS  1 if mode 3's stress and EVERY committed-history field equal mode 0's bit for bit
//         (-1 on the plane-stress view, whose sigma_33 Newton iterates differ by tangent)
//   LOAD  1 if a damage threshold advanced in the step; DT/DC the damage; CRK/BSR the interlock state
//
// Two plane-stress cases condense sigma_33 = 0 with the material's Newton (tolerance tightened to
// 1e-13): their condensed 3x3 tangent is checked against an FD of the CONDENSED map, which re-solves
// eps_33 at every perturbed point. This is the implicit-function-theorem step of the view.
//
// Build with -DRC_DUAL_DOUBLE_ONLY to get the bit-identity digest instead: it touches only the
// pre-WP-263 double API (modes 0/1/2), runs every path and prints an FNV-1a hash of the stress,
// the tangent and every RCHist field per (case, mode). Compiled once against the current header and
// once against the pre-WP header, the two outputs must be byte-identical.
//
//   g++ -O2 -std=c++17 -I SRC/material/nD tests/_testbed/rc_dual_gpp.cpp -o rc_dual
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <chrono>
#include <functional>
#include <string>
#include <vector>

using namespace ladruno_rc_kernel;

// ----------------------------------------------------------------- material ---
static void base(Params& P)
{
    P.E = 30000.0; P.nu = 0.2; P.Kc = 2.0 / 3.0;
    P.cdf = 0.0; P.eta = 0.0; P.betaFloor = 0.1;
    P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
    P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4;
    P.aggSize = 16.0; P.crackStrain = 0.0; P.crackSpacing = 100.0; P.lch = 0.0; P.betaSrMin = 0.01;
    P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
    P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
    P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
    P.tensStiffMode = 0; P.tensStiffC = 200.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
    P.autoReg = false; P.lchRef = 1.0;
    P.betaC = 170.0; P.crackedNu = -1.0; P.addonPlane = 0; P.crackGate = 0;
    P.tieRhoY = 0.0; P.tieRhoZ = 0.0; P.tieFy = 0.0; P.tieEs = 0.0; P.tieB = 0.0;
    P.crackCheck = 0; P.ccSmx = 0.0; P.ccSmt = 0.0; P.lateralNoRoot = 0;
    // the shell / brick battery backbones (MPa), WITH damage, so dt_bar and dc_bar mix a damage
    // part (pt.d, pc.d) and a plastic part (dt_plastic, dc_plastic)
    double Ce[4] = {0.0, 0.0007, 0.0020, 0.0100};
    double Cs[4] = {0.0, 24.0,   30.0,   5.0};
    double Cd[4] = {0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0};
    double Te[3] = {0.0, 0.0001, 0.0010};
    double Ts[3] = {0.0, 3.0,    0.5};
    double Td[3] = {0.0, 0.0,    1.0 - 0.5 / 5.0};
    buildBackbone(P.hc, P.E, Ce, Cs, Cd, 4);
    buildBackbone(P.ht, P.E, Te, Ts, Td, 3);
}

struct Case {
    std::string name, cls;
    Params P;
    int n;
    std::function<void(int, double*)> path;   // step s = 1..n -> eps6 (3D) or (exx, eyy, gxy) (plane stress)
    bool implex;
    bool planeStress;
};

static void set6(double* e, double a, double b, double c, double d, double f, double g)
{ e[0] = a; e[1] = b; e[2] = c; e[3] = d; e[4] = f; e[5] = g; }

static std::vector<Case> makeCases()
{
    std::vector<Case> C;
    auto add = [&](const char* nm, const char* cls, std::function<void(Params&)> cfg, int n,
                   std::function<void(int, double*)> path, bool implex = false, bool ps = false) {
        Case c; c.name = nm; c.cls = cls; base(c.P); cfg(c.P); setupParams(c.P);
        c.n = n; c.path = path; c.implex = implex; c.planeStress = ps;
        C.push_back(c);
    };
    auto none = [](Params&) {};
    // tension past the peak, principal axes fixed (the gxy COLUMN still rotates the projectors)
    auto crackU = [](int s, double* e) { double t = 1.2e-4 + 0.5e-4 * (s - 1); set6(e, t, -0.15*t, -0.25*t, 0, 0, 0); };
    // tension + shear: rotating principal axes, cracked and softening
    auto crackS = [](int s, double* e) { double t = 1.0e-4 + 0.5e-4 * (s - 1); set6(e, t, -0.1*t, -0.2*t, 1.2*t, 0, 0); };
    // compression past the peak, with a little shear
    auto crush = [](int s, double* e) { double t = 0.5e-3 + 0.3e-3 * (s - 1); set6(e, -t, 0.15*t, 0.25*t, 0.1*t, 0, 0); };
    // tension-compression: the MCFT strut (beta on)
    auto tc = [](int s, double* e) { double t = 0.15 * s; set6(e, -1.5e-3*t, 1.0e-3*t, 0.0, 0.3e-3*t, 0, 0); };
    auto ts = [](int s, double* e) { double t = 1.5e-4 + 1.5e-4 * (s - 1); set6(e, t, -0.2*t, -0.22*t, 0.05*t, 0, 0); };
    auto tsR = [](int s, double* e) { double t = 1.5e-4 + 1.5e-4 * (s - 1); set6(e, t, -0.2*t, -0.22*t, 0.8*t, 0, 0); };
    // crack along x first, then shear at fixed opening (interlock capture, then the cap)
    auto ilk = [](int s, double* e) {
        if (s <= 4) { double t = 0.75e-4 * s; set6(e, 4.0 * t, -0.2*t, -0.2*t, 0.02*t, 0, 0); }
        else        { double g = 0.5e-3 * (s - 4); set6(e, 3.0e-4, -0.15e-4, -0.15e-4, 0.6e-5 + g, 0, 0); }
    };
    // the same, then the shear reverses (cyclic friction-slip: unload along G, re-cap)
    auto ilkC = [](int s, double* e) {
        if (s <= 4)       { double t = 0.75e-4 * s; set6(e, 4.0 * t, -0.2*t, -0.2*t, 0.02*t, 0, 0); }
        else if (s <= 10) { double g = 0.5e-3 * (s - 4); set6(e, 3.0e-4, -0.15e-4, -0.15e-4, 0.6e-5 + g, 0, 0); }
        else              { double g = 3.0e-3 - 0.6e-3 * (s - 10); set6(e, 3.0e-4, -0.15e-4, -0.15e-4, 0.6e-5 + g, 0, 0); }
    };
    // equal in-plane principal strains (eps1FromMembrane's degenerate branch, equal principal stresses)
    auto degT = [](int s, double* e) { double t = 1.0e-4 * s; set6(e, t, t, -0.1*t, 0, 0, 0); };
    auto degC = [](int s, double* e) { double t = 1.0e-3 * s; set6(e, -t, -t, 0.45*t, 0, 0, 0); };
    // 3D: shear in gxy AND gxz (the spatial add-on plane)
    auto sp = [](int s, double* e) { double t = 1.0e-4 + 1.0e-4 * (s - 1); set6(e, t, -0.2*t, -0.25*t, 0.5*t, 0.05*t, 0.4*t); };
    // plane stress (exx, eyy, gxy)
    auto psS = [](int s, double* e) { double t = 1.0e-4 + 0.5e-4 * (s - 1); e[0] = t; e[1] = -0.2*t; e[2] = 1.0*t; };
    auto psI = [](int s, double* e) {
        if (s <= 4)       { double t = 0.75e-4 * s; e[0] = 4.0*t; e[1] = -0.2*t; e[2] = 0.02*t; }
        else if (s <= 9)  { double g = 0.6e-3 * (s - 4); e[0] = 3.0e-4; e[1] = -0.15e-4; e[2] = 0.6e-5 + g; }
        else              { double g = 3.0e-3 - 0.7e-3 * (s - 9); e[0] = 3.0e-4; e[1] = -0.15e-4; e[2] = 0.6e-5 + g; }
    };

    add("elastic", "elastic", none, 3, [](int s, double* e) { set6(e, 2e-6*s, -1e-6*s, 0.5e-6*s, 3e-6*s, 1e-6*s, -2e-6*s); });
    add("crack_uniaxial", "cracking", none, 12, crackU);
    add("crack_shear", "cracking", none, 12, crackS);
    add("crackednu_beta", "cracking", [](Params& P) { P.crackedNu = 0.0; P.betaOn = true; }, 12, crackS);
    add("crush", "crushing", none, 12, crush);
    add("beta_tc", "beta", [](Params& P) { P.betaOn = true; P.lublinerTCReduced = true; }, 10, tc);
    add("tensstiff", "tension-stiffened", [](Params& P) { P.tensStiffMode = 1; }, 12, ts);
    add("tensstiff_rot", "tension-stiffened", [](Params& P) { P.tensStiffMode = 1; P.betaOn = true; }, 12, tsR);
    add("interlock_capped", "interlock-capped", [](Params& P) { P.interlockOn = true; }, 10, ilk);
    add("interlock_cyclic", "interlock-cyclic", [](Params& P) { P.interlockOn = true; P.interlockCyclic = true; }, 16, ilkC);
    add("xcrack_dsfm", "interlock-cyclic", [](Params& P) { P.interlockOn = true; P.interlockCyclic = true;
        P.xcrackOn = true; P.shearRetMode = 2; }, 16, ilkC);
    add("degenerate_tension", "degenerate", [](Params& P) { P.betaOn = true; }, 8, degT);
    add("degenerate_compression", "degenerate", [](Params& P) { P.betaOn = true; }, 7, degC);
    add("implex_crack_shear", "implex", [](Params& P) { P.implex = true; }, 12, crackS, true);
    add("implex_ts_beta", "implex", [](Params& P) { P.implex = true; P.betaOn = true; P.tensStiffMode = 1; }, 12, tsR, true);
    add("spatial_ts_beta", "spatial", [](Params& P) { P.addonPlane = 1; P.crackGate = 1; P.betaOn = true;
        P.tensStiffMode = 1; }, 12, sp);
    add("ps_crack_shear", "plane-stress", [](Params& P) { P.betaOn = true; P.tensStiffMode = 1; }, 12, psS, false, true);
    add("ps_interlock_cyclic", "plane-stress", [](Params& P) { P.interlockOn = true; P.interlockCyclic = true;
        P.betaOn = true; }, 14, psI, false, true);
    return C;
}

// --------------------------------------------------- plane stress (sigma_33 = 0) ---
// The material's guarded Newton on eps_33 (LadrunoRCConcrete::setTrialStrain), tolerance tightened
// to 1e-13 so the FD of the condensed map is not Newton noise, then the same Schur condensation.
static bool psMap(const Params& P, const double e3[3], const RCHist& in, double sig6[6], double D[6][6],
                  RCHist& out, bool dox, double tf)
{
    double eps[6] = { e3[0], e3[1], in.strain[2], e3[2], 0.0, 0.0 };
    bool ok = false;
    double prevAbs = 1.0e300;
    for (int it = 0; it < 80; ++it) {
        returnMap3D(P, eps, in, sig6, D, out, true, 0, dox, tf);
        double smag = 0.0; for (int i = 0; i < 6; i++) smag += sig6[i] * sig6[i];
        smag = sqrt(smag);
        const double tol = 1.0e-13 * (smag > 1.0 ? smag : 1.0);
        const double a33 = fabs(sig6[2]);
        if (a33 <= tol) { ok = true; break; }
        const double d22 = D[2][2];
        const double step = (fabs(d22) < 1.0e-14 * P.E) ? (sig6[2] > 0.0 ? 1.0 : -1.0) * 1.0e-6 : sig6[2] / d22;
        eps[2] -= ((a33 > prevAbs) ? 0.5 : 1.0) * step;
        prevAbs = a33;
    }
    const double d22 = D[2][2];
    double col[6], row[6];
    for (int i = 0; i < 6; i++) { col[i] = D[i][2]; row[i] = D[2][i]; }
    for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) D[I][J] -= col[I] * row[J] / d22;
    return ok;
}

static const int PSI[3] = { 0, 1, 3 };

// ------------------------------------------------------------ bit-identity digest ---
static uint64_t fnv(uint64_t h, const void* p, size_t n)
{
    const unsigned char* b = (const unsigned char*)p;
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}
static uint64_t hashHist(uint64_t h, const RCHist& o)
{
    h = fnv(h, o.stress_eff, sizeof(o.stress_eff)); h = fnv(h, o.strain, sizeof(o.strain));
    const double f[] = { o.xt, o.xc, o.dt_bar, o.dc_bar, o.beta, o.eps1, o.cracked, o.crackC, o.crackS,
                         o.wmax, o.betaSr, o.tauCr, o.gammaCr, o.cracked2, o.slipCum, o.xt_old, o.xc_old,
                         o.eps1_old, o.nuCracked, o.tsScale };
    return fnv(h, f, sizeof(f));
}

// one trial + commit for a case, in the given tangent mode (the double API only)
static void stepCase(const Case& c, const Params& P, int s, RCHist& h, double sig[6], double D[6][6], RCHist& o)
{
    double e[6];
    c.path(s, e);
    if (c.planeStress) psMap(P, e, h, sig, D, o, c.implex, 1.0);
    else               returnMap3D(P, e, h, sig, D, o, true, 0, c.implex, 1.0);
}
static void commitCase(const Case& c, const Params& P, int s, RCHist& h, RCHist& o)
{
    if (c.implex) {                       // LadrunoRCConcrete::commitState under IMPL-EX
        double e[6], sig[6], D[6][6];
        c.path(s, e);
        RCHist oi;
        returnMap3D(P, e, h, sig, D, oi, true, 0, false, 1.0);
        oi.xt_old = h.xt; oi.xc_old = h.xc; oi.eps1_old = h.eps1;
        o = oi;
    }
    h = o;
}

#ifdef RC_DUAL_DOUBLE_ONLY
static void digest(const std::vector<Case>& C)
{
    uint64_t total = 14695981039346656037ULL;
    for (const Case& c : C) {
        for (int mode = 0; mode <= 2; ++mode) {
            Params P = c.P; P.tangentMode = mode;
            RCHist h; histZero(h);
            uint64_t hh = 14695981039346656037ULL;
            for (int s = 1; s <= c.n; ++s) {
                double sig[6], D[6][6]; RCHist o;
                stepCase(c, P, s, h, sig, D, o);
                double res = 0.0; RCHist od; double sd[6], Dd[6][6];
                double e[6]; c.path(s, e);
                if (!c.planeStress) returnMap3D(P, e, h, sd, Dd, od, false, &res, c.implex, 1.0);
                hh = fnv(hh, sig, sizeof(sig)); hh = fnv(hh, D, sizeof(D)); hh = hashHist(hh, o);
                hh = fnv(hh, &res, sizeof(res));
                commitCase(c, P, s, h, o);
                hh = hashHist(hh, h);
            }
            printf("DIGEST %s mode %d %016llx\n", c.name.c_str(), mode, (unsigned long long)hh);
            total = fnv(total, &hh, sizeof(hh));
        }
    }
    // the BeamFiber view (driveBeamFiber) on the spatial case, ties on
    {
        Params P = C[0].P; P.addonPlane = 1; P.crackGate = 1; P.tensStiffMode = 1; P.betaOn = true;
        P.tieRhoY = 0.004; P.tieRhoZ = 0.004; P.tieFy = 420.0; P.tieEs = 200000.0;
        RCHist h; histZero(h);
        double ep[2] = { 0.0, 0.0 };
        uint64_t hh = 14695981039346656037ULL;
        for (int s = 1; s <= 14; ++s) {
            double eps[6] = { 1.5e-4 * s, 0, 0, 0.8e-4 * s, 0, 0.3e-4 * s }, sig[6], D[6][6], ept[2], tsg[2];
            RCHist o; int it = 0;
            const int st = driveBeamFiber(P, eps, h, ep, sig, D, o, ept, tsg, &it);
            hh = fnv(hh, &st, sizeof(st)); hh = fnv(hh, &it, sizeof(it));
            hh = fnv(hh, eps, sizeof(eps)); hh = fnv(hh, sig, sizeof(sig)); hh = fnv(hh, D, sizeof(D));
            hh = hashHist(hh, o);
            if (st == STATUS_OK) { h = o; ep[0] = ept[0]; ep[1] = ept[1]; }
        }
        printf("DIGEST beamfiber mode 0 %016llx\n", (unsigned long long)hh);
        total = fnv(total, &hh, sizeof(hh));
    }
    printf("DIGEST TOTAL %016llx\n", (unsigned long long)total);
}
#else
static bool sameHist(const RCHist& a, const RCHist& b) { return hashHist(14695981039346656037ULL, a) == hashHist(14695981039346656037ULL, b); }
// ------------------------------------------------------------------ the oracle ---
static double frob(const double* a, int n) { double s = 0.0; for (int i = 0; i < n; ++i) s += a[i] * a[i]; return sqrt(s); }

static const int NH = 9;
static const double HS[NH] = { 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10 };

static void oracle(const std::vector<Case>& C)
{
    for (const Case& c : C) {
        printf("SCENARIO %s %s\n", c.name.c_str(), c.cls.c_str());
        Params P0 = c.P; P0.tangentMode = 0;
        Params P3 = c.P; P3.tangentMode = 3;
        RCHist h; histZero(h);
        const int nc = c.planeStress ? 3 : 6;
        for (int s = 1; s <= c.n; ++s) {
            double e[6]; c.path(s, e);
            double s0[6], D0[6][6], s3[6], D3[6][6];
            RCHist o0, o3;
            int bits, nbase = 0;
            if (c.planeStress) {
                nbase += !psMap(P0, e, h, s0, D0, o0, c.implex, 1.0);
                nbase += !psMap(P3, e, h, s3, D3, o3, c.implex, 1.0);
                bits = -1;
            } else {
                returnMap3D(P0, e, h, s0, D0, o0, true, 0, c.implex, 1.0);
                returnMap3D(P3, e, h, s3, D3, o3, true, 0, c.implex, 1.0);
                bits = (std::memcmp(s0, s3, sizeof(s0)) == 0 && sameHist(o0, o3)) ? 1 : 0;
            }
            // the tangents as nc x nc (3D: all six; plane stress: the condensed {0,1,3} block)
            double Cd[36], Cm[36];
            for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) {
                const int A = c.planeStress ? PSI[a] : a, B = c.planeStress ? PSI[b] : b;
                Cd[a*nc + b] = D3[A][B]; Cm[a*nc + b] = D0[A][B];
            }
            const double nd = frob(Cd, nc*nc);
            double best = 1e300, bestFd[36] = {0}, errs[NH];
            int nfail = 0;
            for (int k = 0; k < NH; ++k) {
                double Cf[36];
                for (int b = 0; b < nc; ++b) {
                    double ep[6], em[6], sp[6], sm[6], Dd[6][6]; RCHist od;
                    for (int a = 0; a < 6; ++a) { ep[a] = e[a]; em[a] = e[a]; }
                    ep[b] += HS[k]; em[b] -= HS[k];
                    if (c.planeStress) { nfail += !psMap(P3, ep, h, sp, Dd, od, c.implex, 1.0);
                                         nfail += !psMap(P3, em, h, sm, Dd, od, c.implex, 1.0); }
                    else { returnMap3D(P0, ep, h, sp, Dd, od, false, 0, c.implex, 1.0);
                           returnMap3D(P0, em, h, sm, Dd, od, false, 0, c.implex, 1.0); }
                    for (int a = 0; a < nc; ++a) {
                        const int A = c.planeStress ? PSI[a] : a;
                        Cf[a*nc + b] = (sp[A] - sm[A]) / (2.0 * HS[k]);
                    }
                }
                double d[36]; for (int i = 0; i < nc*nc; ++i) d[i] = Cd[i] - Cf[i];
                errs[k] = frob(d, nc*nc) / nd;
                if (errs[k] < best) { best = errs[k]; for (int i = 0; i < nc*nc; ++i) bestFd[i] = Cf[i]; }
            }
            int plat = 0; for (int k = 0; k < NH; ++k) if (errs[k] <= 10.0 * best) ++plat;
            double dm[36], ds[36];
            for (int a = 0; a < nc; ++a) for (int b = 0; b < nc; ++b) {
                dm[a*nc + b] = Cd[a*nc + b] - Cm[a*nc + b];
                ds[a*nc + b] = 0.5 * (Cd[a*nc + b] + Cd[b*nc + a]) - bestFd[a*nc + b];
            }
            const double omit = frob(dm, nc*nc) / nd, sym = frob(ds, nc*nc) / nd;
            const int load = (o3.xt > h.xt || o3.xc > h.xc) ? 1 : 0;
            printf("STEP %d LOAD %d DT %.4f DC %.4f CRK %d BSR %.4f FD %.3e PLAT %d OMIT %.3e SYM %.3e BITS %d NFAIL %d\n",
                   s, load, o3.dt_bar, o3.dc_bar, (int)o3.cracked, o3.betaSr, best, plat, omit, sym, bits, nbase + nfail);
            printf("FDH");
            for (int k = 0; k < NH; ++k) printf(" %.2e", errs[k]);
            printf("\n");
            commitCase(c, P3, s, h, o3);
        }
    }
}

// per-call cost of the tangent, mode 0 vs mode 3 (g++ -O2, no -march=native)
static void timing(const std::vector<Case>& C)
{
    const Case* c = 0;
    for (const Case& x : C) if (x.name == "tensstiff_rot") c = &x;
    Params P0 = c->P; P0.tangentMode = 0;
    Params P3 = c->P; P3.tangentMode = 3;
    RCHist h; histZero(h);
    for (int s = 1; s <= 6; ++s) { double sig[6], D[6][6]; RCHist o; stepCase(*c, P0, s, h, sig, D, o); h = o; }
    double e[6]; c->path(7, e);
    const int n = 20000;
    double sink = 0.0;
    for (int pass = 0; pass < 2; ++pass) {
        const Params& P = pass ? P3 : P0;
        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < n; ++i) {
            double sig[6], D[6][6]; RCHist o;
            e[0] += 1e-14;
            returnMap3D(P, e, h, sig, D, o, true);
            sink += D[0][0];
        }
        auto t1 = std::chrono::steady_clock::now();
        const double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / n;
        printf("TIMING mode %d %.3f us/call\n", pass ? 3 : 0, us);
    }
    printf("TIMING sink %.3e\n", sink);
}
#endif

int main()
{
    std::vector<Case> C = makeCases();
#ifdef RC_DUAL_DOUBLE_ONLY
    digest(C);
#else
    oracle(C);
    timing(C);
#endif
    return 0;
}
