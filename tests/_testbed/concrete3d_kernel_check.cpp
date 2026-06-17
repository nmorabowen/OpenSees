// concrete3d_kernel_check.cpp — standalone g++ self-check of LadrunoConcrete3DKernel.h.
// Proves the C++ surface + hardening building blocks satisfy the same identities the numpy oracle
// (concrete3d_ref.py) gates, with NO OpenSees build. This is the committed byte-harness the review
// asked for; the FULL oracle-numeric-dump diff at the ADR-§5 cross-platform tolerance floors
// (1e-7/1e-8/1e-6) is the P1 build-PR deliverable.
//
// Build + run from the repo root:
//   g++ -std=c++17 -O2 -I. tests/_testbed/concrete3d_kernel_check.cpp -o c3dchk && ./c3dchk
// Exit 0 = all checks pass.

#include "SRC/material/nD/LadrunoConcrete3DKernel.h"
#include <cmath>
#include <cstdio>
#include <initializer_list>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Ladruno::Concrete3D;

static int fails = 0;
static void check(bool ok, const char* name) {
    if (!ok) { std::printf("  FAIL: %s\n", name); ++fails; }
    else      std::printf("  ok  : %s\n", name);
}

int main() {
    // e in (0.5,1] across realistic concrete; surface identities must hold for ALL e.
    for (double e : {0.5001, 0.52, 0.6, 0.8, 1.0}) {
        Params mp; mp.fc = 30.0; mp.ft = 3.0; mp.e = e; mp.m0 = m0Of(mp.fc, mp.ft, e);
        // G1: uniaxial compression (-fc,0,0) on f=0, independent of e/m0
        double sc[6] = {-30.0, 0, 0, 0, 0, 0};
        check(std::fabs(yieldF(sc, mp)) < 1e-12, "uniaxial compression on f=0");
        // G2: uniaxial tension (+ft,0,0) on f=0 (fixes m0)
        double st[6] = {3.0, 0, 0, 0, 0, 0};
        check(std::fabs(yieldF(st, mp)) < 1e-12, "uniaxial tension on f=0");
        // Lode endpoints: r(0)=1/e, r(pi/3)=1
        check(std::fabs(lodeR(0.0, e) * e - 1.0) < 1e-12, "lodeR(0)=1/e");
        check(std::fabs(lodeR(M_PI / 3.0, e) - 1.0) < 1e-12, "lodeR(pi/3)=1");
    }
    // Eq.18 reduces to the failure surface at qh=1 (cap term vanishes) — check a general state.
    {
        Params mp; mp.fc = 30; mp.ft = 3; mp.e = 0.5229; mp.m0 = m0Of(mp.fc, mp.ft, mp.e);
        double s[6] = {-20, -8, -3, 2, -1, 0.5};
        double f_hard = yieldF(s, mp, 1.0, 1.0);
        // independent failure-surface eval (3/2 rho^2/fc^2 + m0[rho r/(sqrt6 fc)+sigV/fc] - 1)
        double xi, rho, th; invariants(s, xi, rho, th);
        double r = lodeR(th, mp.e), sigV_fc = xi / (SQRT3 * mp.fc);
        double f_fail = 1.5 * rho * rho / (mp.fc * mp.fc)
                      + mp.m0 * (rho * r / (SQRT6 * mp.fc) + sigV_fc) - 1.0;
        check(std::fabs(f_hard - f_fail) < 1e-12, "yieldF(qh=1) == failure surface Eq.21");
    }
    // Hardening building blocks (Eq.30-36)
    check(std::fabs(qh1Of(0.0, 0.3, 0.5) - 0.3) < 1e-12, "qh1(0)=qh0");
    check(std::fabs(qh1Of(1.0, 0.3, 0.5) - 1.0) < 1e-12, "qh1(1)=1");
    check(qh2Of(0.5, 0.5) == 1.0 && qh2Of(1.5, 0.5) > 1.0, "qh2: 1 then 1+Hp(kp-1)");
    check(ductilityXh(-30, 30, 0.08, 0.003, 2.0, 1e-6)
        > ductilityXh(6, 30, 0.08, 0.003, 2.0, 1e-6), "ductility: compression more ductile");
    // ductility C0 continuity at Rh=0 (sigV = -fc/3 => Rh=0)
    {
        double sigV0 = -30.0 / 3.0;                       // Rh = -sigV/fc - 1/3 = 0
        double lo = ductilityXh(sigV0 + 1e-7, 30, 0.08, 0.003, 2.0, 1e-6);
        double hi = ductilityXh(sigV0 - 1e-7, 30, 0.08, 0.003, 2.0, 1e-6);
        check(std::fabs(lo - hi) < 1e-6, "ductility C0-continuous at Rh=0");
    }
    std::printf(fails ? "\nKERNEL CHECK: %d FAIL\n" : "\nKERNEL CHECK: ALL PASS\n", fails);
    return fails ? 1 : 0;
}
