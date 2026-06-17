// concrete3d_kernel_check.cpp — standalone g++ self-check of LadrunoConcrete3DKernel.h.
// Two parts, NO OpenSees build:
//   (A) IDENTITIES — the C++ surface + hardening building blocks satisfy the same algebraic
//       identities the numpy oracle (concrete3d_ref.py) gates (self-contained, no fixture).
//   (B) ORACLE-NUMERIC-DUMP DIFF (ADR §5 deliverable) — the C++ return map (perfect-plastic +
//       hardening) AND analytic consistent tangent are diffed against the numpy oracle's dumped
//       reference numbers (concrete3d_oracle_fixture.txt, from gen_concrete3d_fixture.py) at the
//       cross-platform tolerance floors. The perfect-plastic oracle map has an ANALYTIC inner
//       Jacobian => machine-precision match; the hardening oracle map uses a NUMERICAL Jacobian
//       (handoff §6.3) so its reference is itself ~1e-8 and the floors there are relaxed to match.
//
// Build + run from the repo root:
//   g++ -std=c++17 -O2 -I. tests/_testbed/concrete3d_kernel_check.cpp -o c3dchk && ./c3dchk
//   (optional fixture path as argv[1]; defaults to tests/_testbed/concrete3d_oracle_fixture.txt)
// Exit 0 = all checks pass.

#include "SRC/material/nD/LadrunoConcrete3DKernel.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
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

// ---- (A) algebraic identities -------------------------------------------------------------
static void run_identities() {
    std::printf("[A] surface + hardening identities\n");
    for (double e : {0.5001, 0.52, 0.6, 0.8, 1.0}) {
        Params mp; mp.fc = 30.0; mp.ft = 3.0; mp.e = e; mp.m0 = m0Of(mp.fc, mp.ft, e);
        double sc[6] = {-30.0, 0, 0, 0, 0, 0};
        check(std::fabs(yieldF(sc, mp)) < 1e-12, "uniaxial compression on f=0");
        double st[6] = {3.0, 0, 0, 0, 0, 0};
        check(std::fabs(yieldF(st, mp)) < 1e-12, "uniaxial tension on f=0");
        check(std::fabs(lodeR(0.0, e) * e - 1.0) < 1e-12, "lodeR(0)=1/e");
        check(std::fabs(lodeR(M_PI / 3.0, e) - 1.0) < 1e-12, "lodeR(pi/3)=1");
    }
    {
        Params mp; mp.fc = 30; mp.ft = 3; mp.e = 0.5229; mp.m0 = m0Of(mp.fc, mp.ft, mp.e);
        double s[6] = {-20, -8, -3, 2, -1, 0.5};
        double f_hard = yieldF(s, mp, 1.0, 1.0);
        double xi, rho, th; invariants(s, xi, rho, th);
        double r = lodeR(th, mp.e), sigV_fc = xi / (SQRT3 * mp.fc);
        double f_fail = 1.5 * rho * rho / (mp.fc * mp.fc)
                      + mp.m0 * (rho * r / (SQRT6 * mp.fc) + sigV_fc) - 1.0;
        check(std::fabs(f_hard - f_fail) < 1e-12, "yieldF(qh=1) == failure surface Eq.21");
    }
    check(std::fabs(qh1Of(0.0, 0.3, 0.5) - 0.3) < 1e-12, "qh1(0)=qh0");
    check(std::fabs(qh1Of(1.0, 0.3, 0.5) - 1.0) < 1e-12, "qh1(1)=1");
    check(qh2Of(0.5, 0.5) == 1.0 && qh2Of(1.5, 0.5) > 1.0, "qh2: 1 then 1+Hp(kp-1)");
    check(ductilityXh(-30, 30, 0.08, 0.003, 2.0, 1e-6)
        > ductilityXh(6, 30, 0.08, 0.003, 2.0, 1e-6), "ductility: compression more ductile");
    {
        double sigV0 = -30.0 / 3.0;
        double lo = ductilityXh(sigV0 + 1e-7, 30, 0.08, 0.003, 2.0, 1e-6);
        double hi = ductilityXh(sigV0 - 1e-7, 30, 0.08, 0.003, 2.0, 1e-6);
        check(std::fabs(lo - hi) < 1e-6, "ductility C0-continuous at Rh=0");
    }
    // eccentricityFromKupfer recovers the canonical concrete value
    {
        double e = eccentricityFromKupfer(30.0, 3.0, 1.16);
        check(e > 0.50 && e < 0.55, "eccentricityFromKupfer ~0.52 (Kupfer fcc/fc=1.16)");
    }
}

// ---- (B) oracle-numeric-dump diff ---------------------------------------------------------
static Params makeParams(const double pb[12]) {
    Params mp;
    mp.E = pb[0]; mp.nu = pb[1]; mp.fc = pb[2]; mp.ft = pb[3];
    mp.e = pb[4]; mp.Df = pb[5]; mp.qh0 = pb[6]; mp.Hp = pb[7];
    mp.Ah = pb[8]; mp.Bh = pb[9]; mp.Ch = pb[10]; mp.Dh = pb[11];
    mp.m0 = m0Of(mp.fc, mp.ft, mp.e);
    return mp;
}

static void run_oracle_dump(const char* path) {
    std::printf("[B] oracle-numeric-dump diff (%s)\n", path);
    std::ifstream fh(path);
    if (!fh) { std::printf("  FAIL: cannot open fixture %s\n", path); ++fails; return; }
    std::string tok;
    double worst_sig = 0, worst_kp = 0, worst_tan = 0;

    fh >> tok; int npath; fh >> npath;        // NPATH
    for (int p = 0; p < npath; ++p) {
        fh >> tok; std::string label; fh >> label;   // PATH <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        int hard, nsteps; fh >> hard >> nsteps;
        Params mp = makeParams(pb);
        double sig_n[6] = {0,0,0,0,0,0}, kp_n = 0.0, maxs = 0, maxk = 0;
        for (int s = 0; s < nsteps; ++s) {
            double deps[6], sigO[6], kpO;
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 6; ++i) fh >> sigO[i];
            fh >> kpO;
            double sigC[6], kpC, Dt[6][6];
            returnMapTensor(mp, sig_n, deps, kp_n, hard != 0, sigC, kpC, Dt, false);
            for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
            maxk = std::fmax(maxk, std::fabs(kpC - kpO));
            for (int i = 0; i < 6; ++i) sig_n[i] = sigC[i];
            kp_n = kpC;
        }
        worst_sig = std::fmax(worst_sig, maxs); worst_kp = std::fmax(worst_kp, maxk);
        const double tolS = hard ? 1.0e-6 : 1.0e-9;   // hardening oracle ref ~1e-8 (numerical Jac)
        const double tolK = hard ? 1.0e-7 : 1.0e-10;
        bool ok = maxs < tolS && maxk < tolK;
        if (!ok) ++fails;
        std::printf("  %-24s sig_err=%.2e kp_err=%.2e  %s\n", label.c_str(), maxs, maxk, ok ? "ok" : "FAIL");
    }

    fh >> tok; int ntan; fh >> ntan;          // NTAN
    for (int t = 0; t < ntan; ++t) {
        fh >> tok; std::string label; fh >> label;   // TAN <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        int hard; fh >> hard;
        Params mp = makeParams(pb);
        double sig_n[6], kp_n, deps[6], CO[36];
        for (int i = 0; i < 6; ++i) fh >> sig_n[i];
        fh >> kp_n;
        for (int i = 0; i < 6; ++i) fh >> deps[i];
        for (int i = 0; i < 36; ++i) fh >> CO[i];
        double sigC[6], kpC, Dt[6][6];
        returnMapTensor(mp, sig_n, deps, kp_n, hard != 0, sigC, kpC, Dt, true);
        double nd = 0, nn = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
            double d = Dt[A][B] - CO[A*6+B]; nd += d * d; nn += CO[A*6+B] * CO[A*6+B];
        }
        double rel = std::sqrt(nd / nn);
        worst_tan = std::fmax(worst_tan, rel);
        bool ok = rel < 1.0e-6;
        if (!ok) ++fails;
        std::printf("  %-24s tan_rel_err=%.3e  %s\n", label.c_str(), rel, ok ? "ok" : "FAIL");
    }
    std::printf("  WORST  sig=%.2e (pp<1e-9/hard<1e-6)  kp=%.2e (pp<1e-10/hard<1e-7)  tan=%.3e (<1e-6)\n",
                worst_sig, worst_kp, worst_tan);
}

int main(int argc, char** argv) {
    run_identities();
    const char* fix = argc > 1 ? argv[1] : "tests/_testbed/concrete3d_oracle_fixture.txt";
    run_oracle_dump(fix);
    std::printf(fails ? "\nKERNEL CHECK: %d FAIL\n" : "\nKERNEL CHECK: ALL PASS\n", fails);
    return fails ? 1 : 0;
}
