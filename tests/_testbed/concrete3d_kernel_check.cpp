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
        // Tolerances. Perfect-plastic: both maps use an ANALYTIC inner Jacobian => machine-precision
        // (~1e-13). Hardening: per single step the C++ (analytic 4x4 Jac) and the oracle (numerical
        // 4x4 Jac) converge to the SAME root to the residual tol (~1e-11); over a committed MULTI-
        // STEP path the hardening stiffness amplifies those per-step differences to ~1e-8 (the
        // observed 3.9e-8). NOT apex contamination — the driven paths are fully converged + radial
        // (verified). Floors set just above that accumulation bound.
        const double tolS = hard ? 1.0e-6 : 1.0e-9;
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
        double nd = 0, nn = 0, cmax = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
            double d = Dt[A][B] - CO[A*6+B]; nd += d * d; nn += CO[A*6+B] * CO[A*6+B];
            cmax = std::fmax(cmax, std::fabs(CO[A*6+B]));
        }
        double rel = std::sqrt(nd / nn);
        // PER-ENTRY relative check on SIGNIFICANT entries (relative-Frobenius hides large errors on
        // the small antisymmetric off-diagonal terms — and the whole point of this kernel is its
        // NON-symmetric tangent, PR #249 review HIGH-3).
        double perEntry = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B)
            if (std::fabs(CO[A*6+B]) > 0.01 * cmax)
                perEntry = std::fmax(perEntry, std::fabs(Dt[A][B] - CO[A*6+B]) / std::fabs(CO[A*6+B]));
        // ASYMMETRY-norm match: the non-symmetric part is the unsymmetric-solver justification —
        // verify the C++ reproduces the oracle's asymmetry, not just the symmetric bulk.
        double asC = 0, asO = 0, na = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
            asC += (Dt[A][B]-Dt[B][A])*(Dt[A][B]-Dt[B][A]);
            double oab = CO[A*6+B]-CO[B*6+A]; asO += oab*oab; na += CO[A*6+B]*CO[A*6+B];
        }
        double asymC = std::sqrt(asC/na), asymO = std::sqrt(asO/na);
        double asymErr = std::fabs(asymC - asymO);
        worst_tan = std::fmax(worst_tan, rel);
        const double tolTan = hard ? 2.0e-6 : 1.0e-6;   // hardening oracle tangent = FD-of-FD, looser
        bool ok = rel < tolTan && perEntry < 5.0e-3 && asymErr < 1.0e-3;
        if (!ok) ++fails;
        std::printf("  %-24s tan_rel=%.3e perEntry=%.2e asym C/O=%.3f/%.3f  %s\n",
                    label.c_str(), rel, perEntry, asymC, asymO, ok ? "ok" : "FAIL");
    }
    // ---- (B3) P2 DAMAGE: committed damage state + deps -> NOMINAL stress (oracle damaged_step_tensor)
    //      AND (B4) the P3b ANALYTIC dual-projector DAMAGED consistent tangent. B3 pins the damage
    //      STRESS to the oracle (cross-platform). B4 is SELF-CONTAINED: the analytic damaged tangent
    //      (returnMap doTangent=true) is checked against a NUMERICAL central-difference of the SAME
    //      C++ damaged nominal stress (the operator the global Newton consumes) — exactly as the
    //      oracle run_p2e_gate checks analytic vs numerical. Since B3 pins the stress to the oracle,
    //      analytic==numerical here certifies the analytic tangent without any cross-platform FD
    //      noise. The 4 committed states are the smooth PE1/PE2/PE3/PE4 regimes (tension / confined
    //      compression / shear-nonassoc / reversal) — off the sigma_lat=0 Macaulay kink. ----
    fh >> tok; int ndmg; fh >> ndmg;          // NDMG
    double worst_dmg = 0, worst_dtan = 0;
    for (int d = 0; d < ndmg; ++d) {
        fh >> tok; std::string label; fh >> label;   // DMG <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        Params mp = makeParams(pb);
        fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As;
        State in, out;
        for (int i = 0; i < 6; ++i) fh >> in.eps[i];
        for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
        fh >> in.kp;
        fh >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2;
        double deps[6], sigO[6];
        for (int i = 0; i < 6; ++i) fh >> deps[i];
        for (int i = 0; i < 6; ++i) fh >> sigO[i];
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
        double sigC[6], sigEff[6], Da[6][6];
        returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, /*hardening=*/true);  // analytic tangent
        double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
        worst_dmg = std::fmax(worst_dmg, maxs);
        // (B3) damage = effective return (~1e-8 hardening floor) THEN the smooth analytic damage
        // recompose; no new amplification => same hardening floor as the return-map paths.
        bool oks = maxs < 1.0e-6;

        // (B4) NUMERICAL damaged tangent: central diff of returnMap's nominal stress over the strain
        // increment, committed `in` held fixed (mirrors oracle damaged_consistent_tangent, rel_step 1e-6).
        double Cn[6][6];
        const double base = mp.fc / mp.E;
        for (int j = 0; j < 6; ++j) {
            const double dd = 1.0e-6 * (std::fabs(deps[j]) + base);
            double sp[6], sm[6], se[6], junk[6][6]; State o2;
            double strp[6], strm[6];
            for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
            strp[j] += dd; strm[j] -= dd;
            returnMap(mp, strp, in, o2, sp, se, junk, false, -1.0, true);
            returnMap(mp, strm, in, o2, sm, se, junk, false, -1.0, true);
            for (int i = 0; i < 6; ++i) Cn[i][j] = (sp[i] - sm[i]) / (2.0 * dd);
        }
        double nd = 0, nn = 0, cmax = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) {
            double df = Da[A][B] - Cn[A][B]; nd += df * df; nn += Cn[A][B] * Cn[A][B];
            cmax = std::fmax(cmax, std::fabs(Cn[A][B]));
        }
        double rel = std::sqrt(nd / nn);
        double perEntry = 0;
        for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B)
            if (std::fabs(Cn[A][B]) > 0.01 * cmax)
                perEntry = std::fmax(perEntry, std::fabs(Da[A][B] - Cn[A][B]) / std::fabs(Cn[A][B]));
        worst_dtan = std::fmax(worst_dtan, rel);
        // analytic vs numerical of the SAME C++ stress: floored by the central-diff truncation + the
        // single-step hardening return-map convergence noise (the perturbed step re-solves the inner
        // Newton). Both ~1e-5; floors set just above.
        bool okt = rel < 5.0e-5 && perEntry < 5.0e-3;
        if (!oks || !okt) ++fails;
        std::printf("  %-24s nom_sig_err=%.2e %s   tan_rel=%.2e perEntry=%.2e %s\n",
                    label.c_str(), maxs, oks ? "ok" : "FAIL", rel, perEntry, okt ? "ok" : "FAIL");
    }

    // ---- (B5) P3 Tier-2 IMPL-EX: committed IMPLEX state + step -> the REPORTED explicit stress.
    //      Pins returnMap(implex=true) reported sigma to the oracle damaged_step_implex (incl. a
    //      dt-JUMP case that exercises the extrapolation-ratio clamp r<=implexRmax). ----
    fh >> tok; int nimplex; fh >> nimplex;     // NIMPLEX
    double worst_implex = 0;
    for (int d = 0; d < nimplex; ++d) {
        fh >> tok; std::string label; fh >> label;   // IMPLEX <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        Params mp = makeParams(pb);
        fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As >> mp.implexRmax;
        mp.implex = true;
        State in, out;
        for (int i = 0; i < 6; ++i) fh >> in.eps[i];
        for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
        fh >> in.kp;
        fh >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2;
        fh >> in.wt >> in.wc >> in.dwt >> in.dwc >> in.dt_n;
        for (int i = 0; i < 6; ++i) fh >> in.depl[i];
        double dt; fh >> dt;
        double deps[6], sigO[6];
        for (int i = 0; i < 6; ++i) fh >> deps[i];
        for (int i = 0; i < 6; ++i) fh >> sigO[i];
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
        double sigC[6], sigEff[6], Da[6][6];
        returnMap(mp, strain, in, out, sigC, sigEff, Da, true, dt, /*hardening=*/true);
        double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
        worst_implex = std::fmax(worst_implex, maxs);
        bool ok = maxs < 1.0e-6;
        if (!ok) ++fails;
        std::printf("  %-24s implex_sig_err=%.2e %s\n", label.c_str(), maxs, ok ? "ok" : "FAIL");
    }

    std::printf("  WORST  sig=%.2e (pp<1e-9/hard<1e-6)  kp=%.2e (pp<1e-10/hard<1e-7)  tan=%.3e (<1e-6)"
                "  dmg=%.2e (<1e-6)  dmgtan=%.2e (<5e-5)  implex=%.2e (<1e-6)\n",
                worst_sig, worst_kp, worst_tan, worst_dmg, worst_dtan, worst_implex);
}

// ---- (C) robustness / honesty regressions (PR #249 adversarial-review fixes) ----------------
static void run_robustness() {
    std::printf("[C] robustness / honesty regressions\n");
    Params mp; mp.E = 30000; mp.nu = 0.2; mp.fc = 30; mp.ft = 3;
    mp.e = eccentricityFromKupfer(30, 3, 1.16); mp.m0 = m0Of(30, 3, mp.e);
    mp.qh0 = 0.3; mp.Hp = 0.5; mp.Ah = 0.08; mp.Bh = 0.003; mp.Ch = 2.0; mp.Dh = 1e-6;

    // C1 — ADMISSIBILITY/HONESTY fuzz: NO converged plastic hardening return may have kp<kp_n,
    // dlam<0, or a non-finite/off-surface stress. (Pre-fix: ~2% of returns lied with kp<0 / a
    // sign-flipped tension apex for deep-compression trials.)
    unsigned s = 12345u; auto rnd = [&]() { s = s*1664525u + 1013904223u; return (s>>8)*(1.0/16777216.0); };
    int bad = 0, conv = 0, total = 0;
    for (double Df : {0.3, 0.85, 1.0}) {
        mp.Df = Df;
        for (int t = 0; t < 60000; ++t) {
            double w[3]; for (int i = 0; i < 3; ++i) w[i] = (rnd()*2.0 - 1.0) * 80.0;
            double kp_n = rnd() * 1.2;
            PrincipalResult pr = returnMapHardening(w, mp, kp_n);
            ++total;
            if (!pr.plastic) continue;
            bool fin = std::isfinite(pr.sp[0]) && std::isfinite(pr.sp[1]) && std::isfinite(pr.sp[2]) && std::isfinite(pr.kp);
            if (pr.converged) {
                ++conv;
                // a converged plastic return MUST be admissible + on-surface
                if (!fin || pr.kp < kp_n - 1e-9 || pr.dlam < -1e-9 || std::fabs(pr.f_after) > 1e-6*(mp.fc+1.0)) ++bad;
            }
        }
    }
    check(bad == 0, "no converged hardening return is inadmissible (kp>=kp_n, dlam>=0, on-surface, finite)");
    std::printf("       (fuzzed %d trials, %d converged-plastic, %d inadmissible)\n", total, conv, bad);

    // C2 — deep-COMPRESSION near-apex hardening trial must NOT teleport to the tension vertex and
    // claim convergence (the headline apex bug). It must report converged=false + leave kp uncorrupted.
    {
        mp.Df = 1.0;
        double w[3] = {-46.81, -50.41, 2.81}; double kp_n = 0.0702;
        PrincipalResult pr = returnMapHardening(w, mp, kp_n);
        bool teleported = pr.converged && pr.sp[0] > 0 && pr.sp[1] > 0 && pr.sp[2] > 0;  // all-tension vertex
        check(!teleported, "deep-compression trial not falsely converged at the tension apex");
        check(!pr.converged || pr.kp >= kp_n - 1e-9, "apex/overshoot does not commit kp<kp_n");
    }

    // C3 — NaN guard: degenerate params (ft=0 => m0=Inf => apex xi blows up) must NOT leak NaN;
    // returnMapTensor falls back to the finite elastic predictor with status!=0.
    {
        Params bd = mp; bd.ft = 0.0; bd.m0 = m0Of(bd.fc, bd.ft, bd.e);  // m0 = +Inf
        double sig_n[6] = {0,0,0,0,0,0}, deps[6] = {-2e-3, 5e-4, 5e-4, 0, 0, 0};
        double sigN[6], kpN, Dt[6][6];
        int st = returnMapTensor(bd, sig_n, deps, 0.0, true, sigN, kpN, Dt, true);
        bool fin = true; for (int i = 0; i < 6; ++i) { if (!std::isfinite(sigN[i])) fin = false; for (int j=0;j<6;++j) if(!std::isfinite(Dt[i][j])) fin=false; }
        check(fin && st == 2, "degenerate-param trial returns finite fallback + status 2 (no NaN leak)");
    }
}

int main(int argc, char** argv) {
    run_identities();
    const char* fix = argc > 1 ? argv[1] : "tests/_testbed/concrete3d_oracle_fixture.txt";
    run_oracle_dump(fix);
    run_robustness();
    std::printf(fails ? "\nKERNEL CHECK: %d FAIL\n" : "\nKERNEL CHECK: ALL PASS\n", fails);
    return fails ? 1 : 0;
}
