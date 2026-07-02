// Contact-review fix P3 (2026-07) build-free self-check — friction-cap semantics +
// mortar inverse-isomap scale robustness.
//
//   g++ -std=c++11 -I ../../SRC/domain/contact proto_friction_validation.cpp -o fric_val_check && ./fric_val_check
//
// Guards two review findings:
//   1. HIGH-2 — a τmax-only friction config (μ=c=0, τmax>0, in contact) has a ZERO
//      cone radius (min(0,τmax)=0); the old return-map cap≤0 branch returned the RAW
//      ELASTIC traction (an unbounded bond) with a consistent stick tangent. Fixed:
//      cap≤0 in contact ⇒ FREE SLIP (zero traction, slip absorbs the motion, zero
//      tangential stiffness). The command surface additionally REFUSES the config.
//   2. Gate follow-up — LadrunoMortarKernel::inverseIsomap2D converged on an ABSOLUTE
//      tolR=1e-13 with a LENGTH-unit residual (aux-plane UV ~ facet size h): GPs were
//      silently dead for h ≳ 2000 (378/400 at h=5e4). Fixed with the PR-1 project()
//      parametric-step escape.

#include "LadrunoFrictionKernel.h"
#include "LadrunoMortarKernel.h"
#include <cstdio>
#include <cmath>

namespace LF = LadrunoFrictionKernel;
namespace LM = LadrunoMortarKernel;
namespace LP = LadrunoContactProjection;

static int fails = 0;
static void check(bool ok, const char *name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) fails++;
}

int main()
{
    const double n[3] = {0.0, 0.0, 1.0};

    // ---------------- 1. τmax-only = FREE SLIP (pre-fix: unbounded elastic bond) ----------------
    {
        double gTeff[3] = {0.3, -0.1, 0.0}, gpT[3] = {0.0, 0.0, 0.0};
        double tFric[3], gpTtrial[3];
        bool slip = LF::frictionReturnMap(gTeff, gpT, /*N=*/50.0, /*kt=*/1e5, /*mu=*/0.0,
                                          tFric, gpTtrial, /*cohesion=*/0.0, /*tauMax=*/1e6);
        check(slip, "tauMax-only return map reports SLIP");
        check(std::fabs(tFric[0]) + std::fabs(tFric[1]) + std::fabs(tFric[2]) == 0.0,
              "tauMax-only traction is EXACTLY zero (pre-fix: -kt*gTeff, unbounded bond)");
        check(gpTtrial[0] == gTeff[0] && gpTtrial[1] == gTeff[1] && gpTtrial[2] == gTeff[2],
              "tauMax-only slip absorbs the whole tangential motion");
        double Kss[3][3];
        LF::frictionTangentBlock(gTeff, gpT, n, 50.0, 1e6, 1e5, 0.0, true, Kss, 0.0, 1e6);
        double s = 0.0;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s += std::fabs(Kss[i][j]);
        check(s == 0.0, "tauMax-only tangent block is EXACTLY zero (pre-fix: kt*Pt stick)");
    }

    // ---------------- 2. cap>0 branches BYTE-UNCHANGED (stick + Coulomb slip + Tresca) ----------
    {
        double gpT[3] = {0.0, 0.0, 0.0}, tFric[3], gpTtrial[3];
        // stick: ||tT*|| = kt*0.001 = 100 <= cap = mu*N = 0.3*1000 = 300
        double gS[3] = {0.001, 0.0, 0.0};
        bool slip = LF::frictionReturnMap(gS, gpT, 1000.0, 1e5, 0.3, tFric, gpTtrial);
        check(!slip && std::fabs(tFric[0] + 100.0) < 1e-12, "Coulomb stick unchanged (t=-kt*gT)");
        // slip: ||tT*|| = 1e4 > cap = 300 -> t = -300
        double gL[3] = {0.1, 0.0, 0.0};
        slip = LF::frictionReturnMap(gL, gpT, 1000.0, 1e5, 0.3, tFric, gpTtrial);
        check(slip && std::fabs(tFric[0] + 300.0) < 1e-12, "Coulomb slip capped at mu*N");
        // Tresca (mu=0, c>0): cap = c = 200
        slip = LF::frictionReturnMap(gL, gpT, 1000.0, 1e5, 0.0, tFric, gpTtrial, 200.0, 0.0);
        check(slip && std::fabs(tFric[0] + 200.0) < 1e-12, "Tresca (cohesion) cap unchanged");
        // cohesion capped by tauMax: cap = min(0*N + 500, 150) = 150
        slip = LF::frictionReturnMap(gL, gpT, 1000.0, 1e5, 0.0, tFric, gpTtrial, 500.0, 150.0);
        check(slip && std::fabs(tFric[0] + 150.0) < 1e-12, "cohesion+tauMax min() unchanged");
    }

    // ---------------- 3. inverseIsomap2D converges at LARGE facet scales ----------------
    {
        // generic (non-parallelogram) quad in aux-plane UV at scale h; forward-map a known
        // interior (xi,eta), then recover it. Pre-fix: 1e-13 absolute residual unreachable
        // for h >~ 2000 (noise floor ~eps*h) -> conv=false -> the caller SKIPs the GP.
        const double scales[4] = {1.0, 500.0, 5e3, 5e4};
        const char *names[4] = {
            "isomap h=1 (baseline)", "isomap h=500",
            "isomap h=5e3 (pre-fix: dead GPs begin)", "isomap h=5e4 (pre-fix: 378/400 dead)"};
        for (int c = 0; c < 4; c++) {
            double h = scales[c];
            double UV[4][2] = {{0.0, 0.0}, {h, 0.05*h}, {1.15*h, 0.95*h}, {-0.08*h, 1.05*h}};
            int bad = 0;
            for (int t = 0; t < 20; t++) {
                double xi0 = -0.9 + 1.8 * (t % 5) / 4.0, eta0 = -0.9 + 1.8 * (t / 5) / 3.0;
                double N4[4], dxi[4], deta[4];
                LP::shape(4, xi0, eta0, N4, dxi, deta);
                double q[2] = {0.0, 0.0};
                for (int i = 0; i < 4; i++) { q[0] += N4[i]*UV[i][0]; q[1] += N4[i]*UV[i][1]; }
                double xi, eta;
                if (!LM::inverseIsomap2D(4, UV, q, xi, eta) ||
                    std::fabs(xi - xi0) > 1e-8 || std::fabs(eta - eta0) > 1e-8) bad++;
            }
            check(bad == 0, names[c]);
        }
    }

    std::printf("\n%s (%d failure%s)\n", fails == 0 ? "ALL PASS" : "FAILURES",
                fails, fails == 1 ? "" : "s");
    return fails == 0 ? 0 : 1;
}
