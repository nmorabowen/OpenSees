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

// Zone-A wall-time knob (WP concrete3d-damage-drive CI fix), shared with the LADRUNO_FUZZ_TRIALS
// env var the Python oracle's test_flow_potential_gate reads (tests/test_ladrunoConcrete3D_material.py):
// unset => the FULL, thorough trial count (a human running this self-check directly); Zone-A CI
// exports a small value so C1/C4/C7 below stay a fast SMOKE of the same fuzz every push instead of
// paying their full cost. Every fuzz here uses a fixed LCG seed, so a smaller N is a deterministic
// subsample (reproducible either way), never a source of flakiness.
static long fuzzTrials(long defaultN) {
    const char* e = std::getenv("LADRUNO_FUZZ_TRIALS");
    if (!e || !*e) return defaultN;
    char* end = nullptr;
    const long v = std::strtol(e, &end, 10);
    if (end == e || v <= 0) return defaultN;
    return v;
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
    double worst_dmg = 0, worst_dtan = 0, worst_t3 = 0;
    for (int d = 0; d < ndmg; ++d) {
        fh >> tok; std::string label; fh >> label;   // DMG <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        Params mp = makeParams(pb);
        fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As >> mp.ctTemper;   // P2h ctTemper mode (0/1/2)
        fh >> mp.tensionLaw >> mp.epsFc;                           // tension law (0 exp / 1 bilinear) + direct eps_fc
        State in, out;
        for (int i = 0; i < 6; ++i) fh >> in.eps[i];
        for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
        fh >> in.kp;
        fh >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2;
        fh >> in.sigtMax >> in.sigcMax;   // P2g monotone-drive history (8-field line)
        double deps[6], sigO[6];
        for (int i = 0; i < 6; ++i) fh >> deps[i];
        for (int i = 0; i < 6; ++i) fh >> sigO[i];
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
        double sigC[6], sigEff[6], Da[6][6];
        returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, /*hardening=*/true);  // analytic tangent
        double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
        worst_dmg = std::fmax(worst_dmg, maxs);
        // (B7) TIER-3 explicit: with do_tangent=false the kernel skips the tangent assembly entirely; the
        // committed NOMINAL stress + state MUST be byte-identical to the do_tangent=true path (ADR §398:
        // "Tier-3 explicit committed stress == Tier-1 implicit"). Certifies the material is exact under an
        // explicit solver (which never factorizes the indefinite softening tangent).
        { State o3; double s3[6], se3[6], junk3[6][6];
          returnMap(mp, strain, in, o3, s3, se3, junk3, false, -1.0, true);
          double t3 = 0;
          for (int i = 0; i < 6; ++i) {
              t3 = std::fmax(t3, std::fabs(s3[i] - sigC[i]));
              t3 = std::fmax(t3, std::fabs(o3.sig[i] - out.sig[i]));
              t3 = std::fmax(t3, std::fabs(o3.sigEff[i] - out.sigEff[i]));
          }
          t3 = std::fmax(t3, std::fabs(o3.kp - out.kp));
          t3 = std::fmax(t3, std::fabs(o3.et_max - out.et_max));
          t3 = std::fmax(t3, std::fabs(o3.kdt1 - out.kdt1) + std::fabs(o3.kdc1 - out.kdc1));
          worst_t3 = std::fmax(worst_t3, t3);
        }
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
        bool ok3 = worst_t3 < 1.0e-13;                 // Tier-3 do_tangent=false == do_tangent=true (byte)
        if (!oks || !okt || !ok3) ++fails;
        std::printf("  %-24s nom_sig_err=%.2e %s   tan_rel=%.2e perEntry=%.2e %s   t3=%.1e %s\n",
                    label.c_str(), maxs, oks ? "ok" : "FAIL", rel, perEntry, okt ? "ok" : "FAIL",
                    worst_t3, ok3 ? "ok" : "FAIL");
    }

    // ---- (B8) P5 CONFINED-FIBER (§4.6): a committed confined state + one axial increment -> the kernel
    //      driveConfinedFiber NOMINAL axial stress + mobilized hoop pressure + converged lateral strain,
    //      byte-matched to the oracle confined_step. cfib_free (hoopK=0) certifies the reduce-to-free
    //      condensation; cfib_confined the headline Mander-by-mechanism path; cfib_yield the hoop-yield cap.
    fh >> tok; int ncfib; fh >> ncfib;        // NCFIB
    double worst_cfib = 0, worst_cfib_el = 0, worst_cfib_p = 0, worst_cfib_t = 0;
    for (int d = 0; d < ncfib; ++d) {
        fh >> tok; std::string label; fh >> label;   // CFIB <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        Params mp = makeParams(pb);
        fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As >> mp.ctTemper;
        State in, out;
        for (int i = 0; i < 6; ++i) fh >> in.eps[i];
        for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
        fh >> in.kp;
        fh >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2;
        fh >> in.sigtMax >> in.sigcMax;
        double e11_next, hoopK, hoopFy;
        fh >> e11_next >> hoopK >> hoopFy;
        double sigO[6]; for (int i = 0; i < 6; ++i) fh >> sigO[i];
        double elO, pconfO; fh >> elO >> pconfO;
        // impose axial = e11_next on the retained comps {0,3,5}; the lateral block {1,2,4} starts from the
        // committed guess (in.eps) and driveConfinedFiber Newton-solves it against the hoop residual.
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i];
        strain[0] = e11_next; strain[3] = 0.0; strain[5] = 0.0;
        double sigC[6], sigEff[6], Da[6][6];
        int stat = driveConfinedFiber(mp, strain, in, out, sigC, sigEff, Da, true, hoopK, hoopFy, -1.0);
        double serr  = std::fabs(sigC[0] - sigO[0]);
        double elerr = std::fabs(strain[1] - elO);
        double perr  = std::fabs(hoopStress(strain[1], hoopK, hoopFy) - pconfO);
        worst_cfib    = std::fmax(worst_cfib, serr);
        worst_cfib_el = std::fmax(worst_cfib_el, elerr);
        worst_cfib_p  = std::fmax(worst_cfib_p, perr);
        // CONDENSED TANGENT (the novel P5b deliverable): C[0][0]=d(sig_axial)/d(eps_axial) via re-solving
        // the WHOLE confined update (inner lateral Newton + damage) at the perturbed axial strain — the
        // exact operator the fiber Newton consumes. SELF-CONTAINED (no oracle tangent dump), like the B4
        // damaged-tangent check: a numerical central-diff of the (oracle-pinned) confined nominal stress.
        double dC = 1.0e-6 * (std::fabs(e11_next) + mp.fc / mp.E);
        double sp[6], sm[6], seD[6], DaD[6][6]; State op, om;
        double strp[6], strm[6];
        for (int i = 0; i < 6; ++i) { strp[i] = in.eps[i]; strm[i] = in.eps[i]; }
        strp[0] = e11_next + dC; strp[3] = strp[5] = 0.0;
        strm[0] = e11_next - dC; strm[3] = strm[5] = 0.0;
        driveConfinedFiber(mp, strp, in, op, sp, seD, DaD, false, hoopK, hoopFy, -1.0);
        driveConfinedFiber(mp, strm, in, om, sm, seD, DaD, false, hoopK, hoopFy, -1.0);
        double Cfd = (sp[0] - sm[0]) / (2.0 * dC);
        double terr = std::fabs(Da[0][0] - Cfd) / (std::fabs(Cfd) + 1.0);
        worst_cfib_t = std::fmax(worst_cfib_t, terr);
        bool ok = (stat == 0) && serr < 1.0e-6 && elerr < 1.0e-8 && perr < 1.0e-6 && terr < 5.0e-3;
        if (!ok) ++fails;
        std::printf("  %-24s nom_sig_err=%.2e  el_err=%.2e  p_err=%.2e  Ctan_fd=%.2e  %s\n",
                    label.c_str(), serr, elerr, perr, terr, ok ? "ok" : "FAIL");
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
        fh >> in.sigtMax >> in.sigcMax;   // P2g monotone-drive history (8-field line)
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

    // ---- (B6) P3 Duvaut-Lions -eta: committed state + step -> the RELAXED nominal stress (viscous,
    //      eta>0/dt>0) AND the eta->0/dt<=0 inviscid fallback. Pins returnMap(eta,dt) to the oracle
    //      damaged_step_tensor(...,dt), verifies the dt<=0 byte-identity to the inviscid path, AND
    //      confirms the viscous result genuinely DIFFERS from the inviscid one (non-tautology). ----
    fh >> tok; int neta; fh >> neta;          // NETA
    double worst_eta = 0, worst_eta_inv = 0, max_eta_gap = 0;
    for (int d = 0; d < neta; ++d) {
        fh >> tok; std::string label; fh >> label;   // ETA <label>
        double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
        Params mp = makeParams(pb);
        fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As;
        double eta, dt; fh >> eta >> dt;
        mp.eta = eta;
        State in, out;
        for (int i = 0; i < 6; ++i) fh >> in.eps[i];
        for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
        fh >> in.kp;
        fh >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2;
        fh >> in.sigtMax >> in.sigcMax;   // P2g monotone-drive history (8-field line)
        double deps[6], sigVisc[6], sigInv[6];
        for (int i = 0; i < 6; ++i) fh >> deps[i];
        for (int i = 0; i < 6; ++i) fh >> sigVisc[i];
        for (int i = 0; i < 6; ++i) fh >> sigInv[i];
        double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
        double sigC[6], sigEff[6], Da[6][6];
        returnMap(mp, strain, in, out, sigC, sigEff, Da, true, dt, /*hardening=*/true);     // viscous (eta>0, dt>0)
        double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigVisc[i]));
        worst_eta = std::fmax(worst_eta, maxs);
        State o2; double sigI[6], se2[6], junk[6][6];
        returnMap(mp, strain, in, o2, sigI, se2, junk, false, -1.0, /*hardening=*/true);    // dt<=0 => inviscid fallback
        double maxi = 0, gap = 0;
        for (int i = 0; i < 6; ++i) {
            maxi = std::fmax(maxi, std::fabs(sigI[i] - sigInv[i]));
            gap = std::fmax(gap, std::fabs(sigC[i] - sigI[i]));
        }
        worst_eta_inv = std::fmax(worst_eta_inv, maxi);
        max_eta_gap = std::fmax(max_eta_gap, gap);
        bool ok = maxs < 1.0e-6 && maxi < 1.0e-6;
        if (!ok) ++fails;
        std::printf("  %-24s eta_sig_err=%.2e %s  inv_fallback_err=%.2e  visc-inv gap=%.2e\n",
                    label.c_str(), maxs, ok ? "ok" : "FAIL", maxi, gap);
    }
    // non-tautology: at least one viscous case must differ materially from its inviscid fallback
    bool eta_nontrivial = (neta == 0) || (max_eta_gap > 1.0e-3);
    if (!eta_nontrivial) ++fails;
    if (neta > 0)
        std::printf("  eta nontrivial (max viscous-inviscid gap >1e-3): %s (gap=%.2e)\n",
                    eta_nontrivial ? "ok" : "FAIL", max_eta_gap);

    // ---- (B9) Gc-as-energy calibration: C++ compressionEnergyDensity == oracle compression_energy_density ----
    if (fh >> tok && tok == "NGCT") {
        int ngct; fh >> ngct;
        for (int d = 0; d < ngct; ++d) {
            fh >> tok; std::string label; fh >> label;   // GCT <label>
            double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
            Params mp = makeParams(pb);
            double efc, gO, pkO; fh >> mp.Gf >> mp.As >> efc >> gO >> pkO;
            double pkC = 0.0;
            const double gC = compressionEnergyDensity(mp, efc, 20000, &pkC);
            const double rel = std::fabs(gC - gO) / std::fabs(gO), relp = std::fabs(pkC - pkO) / pkO;
            const bool ok = rel < 1.0e-6 && relp < 1.0e-7;
            if (!ok) ++fails;
            std::printf("  %-24s g_rel=%.2e peak_rel=%.2e  %s\n", label.c_str(), rel, relp, ok ? "ok" : "FAIL");
        }
    }

    // ---- (B10) FULL CDPM2 PLASTIC POTENTIAL (B1): cdpm2-flow paths (incl. a sub-incremented one) + tangents ----
    if (fh >> tok && tok == "NPATHF") {
        int n; fh >> n;
        for (int p = 0; p < n; ++p) {
            fh >> tok; std::string label; fh >> label;
            double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
            Params mp = makeParams(pb);
            int nsteps; fh >> mp.flowPotential >> mp.maxSubIncr >> nsteps;
            double sig_n[6] = {0,0,0,0,0,0}, kp_n = 0.0, maxs = 0, maxk = 0;
            for (int s = 0; s < nsteps; ++s) {
                double deps[6], sigO[6], kpO;
                for (int i = 0; i < 6; ++i) fh >> deps[i];
                for (int i = 0; i < 6; ++i) fh >> sigO[i];
                fh >> kpO;
                double sigC[6], kpC, Dt[6][6];
                returnMapTensor(mp, sig_n, deps, kp_n, true, sigC, kpC, Dt, false);
                for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]) * 30.0 / mp.fc);   // MPa-equivalent (SI cases)
                maxk = std::fmax(maxk, std::fabs(kpC - kpO) / (1.0 + std::fabs(kpO)));
                for (int i = 0; i < 6; ++i) sig_n[i] = sigC[i];
                kp_n = kpC;
            }
            worst_sig = std::fmax(worst_sig, maxs); worst_kp = std::fmax(worst_kp, maxk);
            const bool ok = maxs < 1.0e-6 && maxk < 1.0e-7;
            if (!ok) ++fails;
            std::printf("  %-24s sig_err=%.2e kp_err=%.2e  %s\n", label.c_str(), maxs, maxk, ok ? "ok" : "FAIL");
        }
    }
    if (fh >> tok && tok == "NTANF") {
        int n; fh >> n;
        for (int t = 0; t < n; ++t) {
            fh >> tok; std::string label; fh >> label;
            double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
            Params mp = makeParams(pb);
            fh >> mp.flowPotential;
            double sig_n[6], kp_n, deps[6], CO[36];
            for (int i = 0; i < 6; ++i) fh >> sig_n[i];
            fh >> kp_n;
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 36; ++i) fh >> CO[i];
            double sigC[6], kpC, Dt[6][6];
            returnMapTensor(mp, sig_n, deps, kp_n, true, sigC, kpC, Dt, true);
            double nd = 0, nn = 0, cmax = 0, perEntry = 0;
            for (int A = 0; A < 36; ++A) { const double d = Dt[A/6][A%6] - CO[A]; nd += d*d; nn += CO[A]*CO[A]; cmax = std::fmax(cmax, std::fabs(CO[A])); }
            for (int A = 0; A < 36; ++A) if (std::fabs(CO[A]) > 0.01 * cmax)
                perEntry = std::fmax(perEntry, std::fabs(Dt[A/6][A%6] - CO[A]) / std::fabs(CO[A]));
            const double rel = std::sqrt(nd / nn);
            worst_tan = std::fmax(worst_tan, rel);
            const bool ok = rel < 2.0e-6 && perEntry < 5.0e-3;
            if (!ok) ++fails;
            std::printf("  %-24s tan_rel=%.3e perEntry=%.2e  %s\n", label.c_str(), rel, perEntry, ok ? "ok" : "FAIL");
        }
    }

    // ---- (B11) CDPM2 COMPRESSIVE DAMAGE DRIVE (B2): nominal stress vs oracle + analytic vs numerical tangent ----
    if (fh >> tok && tok == "NDMGC") {
        int n; fh >> n;
        for (int d = 0; d < n; ++d) {
            fh >> tok; std::string label; fh >> label;
            double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
            Params mp = makeParams(pb);
            fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As >> mp.ctTemper >> mp.tensionLaw >> mp.epsFc
               >> mp.flowPotential >> mp.maxSubIncr >> mp.compDrive;
            State in, out;
            for (int i = 0; i < 6; ++i) fh >> in.eps[i];
            for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
            fh >> in.kp >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2
               >> in.sigtMax >> in.sigcMax >> in.eqc >> in.etPrev;
            double deps[6], sigO[6];
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 6; ++i) fh >> sigO[i];
            double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
            double sigC[6], sigEff[6], Da[6][6];
            returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, true);
            double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
            double Cn[6][6]; const double base = mp.fc / mp.E;
            for (int j = 0; j < 6; ++j) {
                const double dd = 1.0e-6 * (std::fabs(deps[j]) + base);
                double sp[6], sm[6], se[6], junk[6][6]; State o2; double strp[6], strm[6];
                for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
                strp[j] += dd; strm[j] -= dd;
                returnMap(mp, strp, in, o2, sp, se, junk, false, -1.0, true);
                returnMap(mp, strm, in, o2, sm, se, junk, false, -1.0, true);
                for (int i = 0; i < 6; ++i) Cn[i][j] = (sp[i] - sm[i]) / (2.0 * dd);
            }
            double nd = 0, nn = 0;
            for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) { const double df = Da[A][B] - Cn[A][B]; nd += df * df; nn += Cn[A][B] * Cn[A][B]; }
            const double rel = std::sqrt(nd / nn);
            worst_dmg = std::fmax(worst_dmg, maxs); worst_dtan = std::fmax(worst_dtan, rel);
            const bool ok = maxs < 1.0e-6 && rel < 5.0e-5;
            if (!ok) ++fails;
            std::printf("  %-24s nom_sig_err=%.2e  tan_rel=%.2e  wc=%.3f  %s\n", label.c_str(), maxs, rel, out.wc, ok ? "ok" : "FAIL");
        }
    }

    // ---- (B12) PV20 TENSION->COMPRESSION TEMPER (-tcTemper): nominal stress vs oracle + analytic vs FD tangent ----
    if (fh >> tok && tok == "NDMGT") {
        int n; fh >> n;
        for (int d = 0; d < n; ++d) {
            fh >> tok; std::string label; fh >> label;
            double pb[12]; for (int i = 0; i < 12; ++i) fh >> pb[i];
            Params mp = makeParams(pb);
            fh >> mp.Gf >> mp.Gc >> mp.lch >> mp.As >> mp.ctTemper >> mp.tensionLaw >> mp.epsFc
               >> mp.flowPotential >> mp.maxSubIncr >> mp.compDrive >> mp.tcTemper;
            State in, out;
            for (int i = 0; i < 6; ++i) fh >> in.eps[i];
            for (int i = 0; i < 6; ++i) fh >> in.sigEff[i];
            fh >> in.kp >> in.et_max >> in.kdt1 >> in.kdt2 >> in.kdc >> in.kdc1 >> in.kdc2
               >> in.sigtMax >> in.sigcMax >> in.eqc >> in.etPrev;
            double deps[6], sigO[6];
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 6; ++i) fh >> sigO[i];
            double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
            double sigC[6], sigEff[6], Da[6][6];
            returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, true);
            double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
            double Cn[6][6]; const double base = mp.fc / mp.E;
            for (int j = 0; j < 6; ++j) {
                const double dd = 1.0e-6 * (std::fabs(deps[j]) + base);
                double sp[6], sm[6], se[6], junk[6][6]; State o2; double strp[6], strm[6];
                for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
                strp[j] += dd; strm[j] -= dd;
                returnMap(mp, strp, in, o2, sp, se, junk, false, -1.0, true);
                returnMap(mp, strm, in, o2, sm, se, junk, false, -1.0, true);
                for (int i = 0; i < 6; ++i) Cn[i][j] = (sp[i] - sm[i]) / (2.0 * dd);
            }
            double nd = 0, nn = 0;
            for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) { const double df = Da[A][B] - Cn[A][B]; nd += df * df; nn += Cn[A][B] * Cn[A][B]; }
            const double rel = std::sqrt(nd / nn);
            worst_dmg = std::fmax(worst_dmg, maxs); worst_dtan = std::fmax(worst_dtan, rel);
            const bool ok = maxs < 1.0e-6 && rel < 5.0e-5;
            if (!ok) ++fails;
            std::printf("  %-24s nom_sig_err=%.2e  tan_rel=%.2e  wc=%.3f  %s\n", label.c_str(), maxs, rel, out.wc, ok ? "ok" : "FAIL");
        }
    }

    std::printf("  WORST  sig=%.2e (pp<1e-9/hard<1e-6)  kp=%.2e (pp<1e-10/hard<1e-7)  tan=%.3e (<1e-6)"
                "  dmg=%.2e (<1e-6)  dmgtan=%.2e (<5e-5)  implex=%.2e (<1e-6)  eta=%.2e (<1e-6)  eta_inv=%.2e (<1e-6)\n",
                worst_sig, worst_kp, worst_tan, worst_dmg, worst_dtan, worst_implex, worst_eta, worst_eta_inv);
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
    const long N1 = fuzzTrials(60000);   // LADRUNO_FUZZ_TRIALS knob; full run = 180000 trials (3 Df values)
    for (double Df : {0.3, 0.85, 1.0}) {
        mp.Df = Df;
        for (long t = 0; t < N1; ++t) {
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

    // C4 — VERTEX return (WP concrete3d-oracle-diagnosis). (a) NO converged return may sign-flip a
    // compressive-mean trial onto an all-tensile state (the PR #249 gate let 7617/165641 of these through:
    // they are 'admissible'); (b) a hydrostatic-compression trial on the closed cap (kp_n<1) must YIELD at
    // the compressive vertex (pre-fix: elastic fallback, status 2 on every increment); (c) the vertex
    // tangent is purely volumetric (zero deviatoric stiffness).
    {
        unsigned s2 = 777u; auto rn = [&]() { s2 = s2*1664525u + 1013904223u; return (s2>>8)*(1.0/16777216.0); };
        int flips = 0;
        const long N4 = fuzzTrials(30000);   // LADRUNO_FUZZ_TRIALS knob; full run = 90000 trials (3 Df values)
        for (double Df : {0.3, 0.85, 1.0}) {
            mp.Df = Df;
            for (long t = 0; t < N4; ++t) {
                double w[3]; for (int i = 0; i < 3; ++i) w[i] = (rn()*2.0 - 1.0) * 80.0;
                PrincipalResult pr = returnMapHardening(w, mp, rn() * 1.2);
                if (pr.plastic && pr.converged && (w[0]+w[1]+w[2]) < 0.0 && pr.sp[0] > 0 && pr.sp[1] > 0 && pr.sp[2] > 0) ++flips;
            }
        }
        check(flips == 0, "no converged return sign-flips a compressive-mean trial to the tension vertex");
        mp.Df = 1.0;
        double w[3] = {-150.0, -150.0, -150.0};
        PrincipalResult pr = returnMapHardening(w, mp, 0.05);
        check(pr.converged && pr.apex && pr.sp[0] < 0.0 && pr.sp[0] > -150.0 && pr.kp > 0.05
              && std::fabs(pr.f_after) < 1e-9, "hydrostatic compression yields at the closed-cap vertex");
        double sig_n[6] = {0,0,0,0,0,0}, deps[6] = {-3e-3, -3e-3, -3e-3, 1e-6, 0, 0};
        double sigN[6], kpN, Dt[6][6];
        int st = returnMapTensor(mp, sig_n, deps, 0.05, true, sigN, kpN, Dt, true);
        const double devmax = std::fmax(std::fabs(Dt[3][3]), std::fabs(Dt[0][0] - Dt[1][0]));
        check(st == 0 && devmax < 1e-9 * mp.E, "vertex tangent is purely volumetric (zero deviatoric stiffness)");
    }

    // C5 — Gc = PHYSICAL compressive energy (WP concrete3d-oracle-diagnosis): with eps_fc from the calibration
    // table, single-point uniaxial compression dissipates Gc (post-peak, per unit area over lch) within 5% for
    // three (Gc, lch) pairs; a direct epsFc equal to the legacy Gc/(fc lch) is byte-identical to the legacy path.
    {
        Params q = mp; q.Df = 1.0; q.As = 2.0; q.Gf = 0.1; q.Gc = 1.0; q.lch = 1.0;
        double efc[GC_TABLE_N], g[GC_TABLE_N];
        calibrateEpsFcTable(q, efc, g);
        const double pairs[3][2] = {{30.0, 100.0}, {10.0, 50.0}, {40.0, 50.0}};
        double worst = 0.0; int bad = 0;
        for (const auto& pr : pairs) {
            int st = 0; const double e = epsFcFromGc(efc, g, pr[0], pr[1], &st);
            const double rel = std::fabs(compressionEnergyDensity(q, e) * pr[1] - pr[0]) / pr[0];
            worst = std::fmax(worst, rel); if (st != 0) ++bad;
        }
        check(bad == 0 && worst < 0.05, "Gc calibration: uniaxial compression dissipates Gc within 5% (3 Gc/lch pairs)");
        std::printf("       (worst |dissipated/Gc - 1| = %.2e)\n", worst);
        {   // the same calibration with the B1 CDPM2 potential + sub-incrementation (the wrapper default)
            Params qc = q; qc.flowPotential = 1; qc.maxSubIncr = 10; qc.Df = 0.85; qc.compDrive = 1;   // B1 + B2
            double efc2[GC_TABLE_N], g2[GC_TABLE_N];
            calibrateEpsFcTable(qc, efc2, g2);
            double w2 = 0.0; int b2 = 0;
            for (const auto& pr : pairs) {
                int st = 0; const double e = epsFcFromGc(efc2, g2, pr[0], pr[1], &st);
                w2 = std::fmax(w2, std::fabs(compressionEnergyDensity(qc, e) * pr[1] - pr[0]) / pr[0]); if (st != 0) ++b2;
            }
            check(b2 == 0 && w2 < 0.01, "Gc calibration with the CDPM2 potential + CDPM2 drive: dissipates Gc within 1%");
            std::printf("       (worst |dissipated/Gc - 1| = %.2e)\n", w2);
        }
        Params a = q; a.Gc = 30.0; a.lch = 100.0; a.epsFc = 0.0;          // legacy Gc/(fc lch)
        Params b = a; b.epsFc = a.Gc / (a.fc * a.lch);                    // the same value, given directly
        State sa, sb; bool same = true;
        for (int k = 1; k < 80; ++k) {
            double stn[6] = {-5e-5 * k, 1.2e-5 * k, 1.2e-5 * k, 0, 0, 0};
            State oa, ob; double s1[6], s2[6], e1[6], e2[6], D1[6][6], D2[6][6];
            returnMap(a, stn, sa, oa, s1, e1, D1, true, -1.0, true);
            returnMap(b, stn, sb, ob, s2, e2, D2, true, -1.0, true);
            for (int i = 0; i < 6; ++i) if (s1[i] != s2[i]) same = false;
            sa = oa; sb = ob;
        }
        check(same, "-epsFc = Gc/(fc lch) reproduces the legacy compression path byte-for-byte");
    }

    // C6 — full crack opening under the BILINEAR law (omega_t reaches 1 EXACTLY): a single GP driven in
    // uniaxial tension (axial prescribed, lateral NOMINAL stress = 0 by a plain, UNGUARDED 2x2 Newton on the
    // material's own damaged tangent — what the global solver does) must stay finite to 5e-3 (~3x the
    // separation strain), converge every step, and dissipate Gf (+ the small pre-peak plastic work). Pre-fix
    // the tensile-direction tangent was exactly singular at omega_t = 1 => NaN (the single-brick run printed
    // thousands of return-map warnings with runaway lateral strains). SI units, fork defaults.
    {
        Params t; t.E = 30e9; t.nu = 0.2; t.fc = 30e6; t.ft = 3e6; t.e = eccentricityFromKupfer(30e6, 3e6, 1.16);
        t.m0 = m0Of(t.fc, t.ft, t.e); t.Df = 1.0; t.As = 2.0; t.qh0 = 0.3; t.Hp = 0.5;
        t.Gf = 120.0; t.Gc = 30e3; t.lch = 0.1; t.tensionLaw = 1;
        State in; double ex = 0, ey = 0, ez = 0, W = 0, sp = 0, ep = 0; const int N = 1000; const double de = 5e-3 / N;
        double Dc[6][6]; elasticC(t, Dc); bool finite = true; int unconv = 0;
        for (int n = 1; n <= N && finite; ++n) {
            ez += de;
            { const double a = Dc[0][0], b = Dc[0][1], c = Dc[1][0], d = Dc[1][1], det = a * d - b * c;
              ex += (d * (-Dc[0][2] * de) - b * (-Dc[1][2] * de)) / det;
              ey += (-c * (-Dc[0][2] * de) + a * (-Dc[1][2] * de)) / det; }
            State out; double s[6], se[6], D[6][6]; bool conv = false;
            for (int it = 0; it < 50; ++it) {
                double strain[6] = {ex, ey, ez, 0, 0, 0};
                returnMap(t, strain, in, out, s, se, D, true, -1.0, true);
                if (std::fabs(s[0]) + std::fabs(s[1]) < 1e-3) { conv = true; break; }
                const double a = D[0][0], b = D[0][1], c = D[1][0], d = D[1][1], det = a * d - b * c;
                ex -= (d * s[0] - b * s[1]) / det; ey -= (-c * s[0] + a * s[1]) / det;
            }
            if (!conv) ++unconv;
            if (!std::isfinite(ex) || !std::isfinite(s[2])) finite = false;
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Dc[i][j] = D[i][j];
            W += 0.5 * (s[2] + sp) * (ez - ep); sp = s[2]; ep = ez; in = out;
        }
        const double G = W * t.lch;
        check(finite && unconv == 0 && std::fabs(sp) < 1.0e-4 * t.ft && G > 120.0 && G < 125.0,   // B2: 1e-6 residual (OMEGA_MAX)
              "bilinear full crack opening: finite, converged, dissipates Gf (residual tangent stiffness)");
        std::printf("       (W*lch = %.2f N/m vs Gf = 120, unconverged steps %d, sigma_end = %.2e Pa)\n", G, unconv, sp);
    }

    // C7 — FULL CDPM2 POTENTIAL robustness (B1): 120k random tensor increments from the stress-free state with
    // random kappa_p in [0,1.5], Df in {0.85, 1.0, 0.6}, Hp in {0.01, 0.5} — generic, deep NEAR-AXIS compression
    // (the closed-cap regime), uniaxial-ish compression with shear, near-apex tension — with sub-incrementation:
    // 0 failed returns, 0 sign flips, 0 non-finite / kappa-decreasing states (without sub-incrementation ~30% of
    // the plastic trials fail the direct return; the legacy flow fails ~43%).
    {
        Params q = mp; q.flowPotential = 1; q.maxSubIncr = 10;
        unsigned s3 = 2024u; auto rn = [&]() { s3 = s3 * 1664525u + 1013904223u; return (s3 >> 8) * (1.0 / 16777216.0); };
        long plastic = 0, fail = 0, flips = 0, bad = 0;
        const long N7 = fuzzTrials(120000);   // LADRUNO_FUZZ_TRIALS knob (Zone-A CI: WP concrete3d-damage-drive)
        for (long t = 0; t < N7; ++t) {
            const int mode = t % 4;
            q.Df = (t % 3 == 0) ? 0.85 : ((t % 3 == 1) ? 1.0 : 0.6);
            q.Hp = (t % 2) ? 0.01 : 0.5;
            double d[6];
            if (mode == 0) { for (int i = 0; i < 6; ++i) d[i] = (rn() * 2 - 1) * 3e-3; for (int i = 3; i < 6; ++i) d[i] *= 0.5; }
            else if (mode == 1) { const double v = -rn() * 8e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 2e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 1e-4; }
            else if (mode == 2) { d[0] = -rn() * 5e-3; d[1] = (rn() * 2 - 1) * 1.5e-3; d[2] = (rn() * 2 - 1) * 1.5e-3; d[3] = (rn() * 2 - 1) * 5e-4; d[4] = d[5] = 0; }
            else { const double v = rn() * 1e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 3e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 2e-4; }
            const double kp_n = rn() * 1.5;
            double sig0[6] = {0, 0, 0, 0, 0, 0}, sn[6], kn, Dt[6][6], tr[6];
            const int st = returnMapTensor(q, sig0, d, kp_n, true, sn, kn, Dt, true);
            elasticPredTensor(sig0, d, q, tr);
            if (yieldF(tr, q, qh1Of(kp_n, q.qh0, q.Hp), qh2Of(kp_n, q.Hp)) > 1e-11 * q.fc) ++plastic;
            if (st != 0) { ++fail; continue; }
            bool fin = std::isfinite(kn);
            for (int i = 0; i < 6; ++i) { fin = fin && std::isfinite(sn[i]); for (int j = 0; j < 6; ++j) fin = fin && std::isfinite(Dt[i][j]); }
            if (!fin || kn < kp_n - 1e-12) ++bad;
            double B[3][3], ws[3], Vs[3][3]; voigtToMat(sn, B); eig3sym(B, ws, Vs);
            if (tr[0] + tr[1] + tr[2] < 0 && ws[0] > 0 && ws[1] > 0 && ws[2] > 0) ++flips;
        }
        // plastic > N7/2 is the SAME ratio as the original hardcoded 60000/120000 (measured ~72% plastic
        // at full N); parametrized so the check keeps its selectivity at any LADRUNO_FUZZ_TRIALS value
        // instead of becoming unsatisfiable (or vacuous) at a reduced N.
        check(fail == 0 && flips == 0 && bad == 0 && plastic > N7 / 2,
              "CDPM2 potential fuzz: 0 failed returns, 0 sign flips, 0 non-finite");
        std::printf("       (N=%ld, plastic %ld, failed %ld, flips %ld, bad %ld)\n", N7, plastic, fail, flips, bad);
    }

    // C8 — B1 + B2 end-to-end vs the paper: Imran & Pantazopoulou triaxial compression (Grassl 2013 Fig.11
    // parameters: E=30 GPa nu=0.2 fc=47.4 ft=4.74 MPa As=15 eps_fc=1e-4, CDPM2 defaults Hp=0.01 Df=0.85),
    // hydrostatic ramp to -s3 then axial strain with the lateral NOMINAL stress held at -s3: the peak axial stress
    // equals the paper's own CDPM2 curve (vector-extracted) to 1 % at s3 = 0 / 8.6 / 43 MPa (-47.3 / -88.7 / -194.9).
    {
        Params q; q.E = 30000; q.nu = 0.2; q.fc = 47.4; q.ft = 4.74; q.As = 15.0; q.epsFc = 1e-4; q.Hp = 0.01; q.Df = 0.85;
        q.qh0 = 0.3; q.e = eccentricityFromKupfer(q.fc, q.ft, 1.16); q.m0 = m0Of(q.fc, q.ft, q.e);
        q.flowPotential = 1; q.maxSubIncr = 10; q.compDrive = 1; q.tensionLaw = 1; q.Gf = 0.158; q.Gc = 1.0; q.lch = 100.0;
        const double s3s[3] = {0.0, 8.6, 43.0}, ref[3] = {-47.3, -88.7, -194.9};
        double worst = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double s3 = s3s[c]; State in; double peak = 0.0;
            const double K = q.E / (3.0 * (1.0 - 2.0 * q.nu)); double eh = 0.0;
            for (int i = 1; i <= 60 && s3 > 0.0; ++i) {
                const double tgt = -s3 * i / 60.0; double lo = eh - 20.0 * s3 / (3.0 * K) / 60.0, hi = eh;
                State o; double sg[6], se[6], D[6][6];
                for (int it = 0; it < 100; ++it) { const double mid = 0.5 * (lo + hi); double e[6] = {mid, mid, mid, 0, 0, 0};
                    returnMap(q, e, in, o, sg, se, D, false, -1.0, true); if (sg[0] > tgt) hi = mid; else lo = mid; }
                eh = 0.5 * (lo + hi); double e[6] = {eh, eh, eh, 0, 0, 0}; returnMap(q, e, in, o, sg, se, D, false, -1.0, true); in = o;
            }
            double lat = in.eps[1]; const double e0 = in.eps[0];
            for (int k = 1; k <= 460; ++k) {
                const double ex = e0 - 1e-4 * k; State o; double sg[6], se[6], D[6][6];
                auto g = [&](double l) { double e[6] = {ex, l, l, 0, 0, 0}; returnMap(q, e, in, o, sg, se, D, false, -1.0, true); return sg[1] + s3; };
                double h = 5e-6, a = lat - h, b = lat + h + 3e-5, fa = g(a), fb = g(b);
                for (int n = 0; fa * fb > 0 && n < 80; ++n) { h *= 1.5; a = lat - h; b = lat + h; fa = g(a); fb = g(b); }
                double m = 0.5 * (a + b);
                for (int n = 0; n < 100; ++n) { m = 0.5 * (a + b); const double fm = g(m); if (std::fabs(fm) < 1e-9 * q.fc) break; if (fa * fm <= 0) { b = m; fb = fm; } else { a = m; fa = fm; } }
                g(m); lat = m; in = o;
                if (sg[0] < peak) peak = sg[0];
            }
            worst = std::fmax(worst, std::fabs(peak / ref[c] - 1.0));
            std::printf("       (IP s3=%5.1f MPa: peak %.2f MPa vs paper CDPM2 %.1f)\n", s3, peak, ref[c]);
        }
        check(worst < 0.01, "Imran-Pantazopoulou triaxial peaks == the paper's CDPM2 curves to 1% (B1+B2)");
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

// ---- (D) apeGmsh #2621 F2: the sub-step chain-rule tangent + the analytic damaged tangent ---------------------
static double relF6(const double A[6][6], const double R[6][6]) {
    double n = 0, d = 0;
    for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { n += (A[i][j] - R[i][j]) * (A[i][j] - R[i][j]); d += R[i][j] * R[i][j]; }
    return d > 0 ? std::sqrt(n / d) : std::sqrt(n);
}
static double relV6(const double a[6], const double r[6]) {
    double n = 0, d = 0;
    for (int i = 0; i < 6; ++i) { n += (a[i] - r[i]) * (a[i] - r[i]); d += r[i] * r[i]; }
    return d > 0 ? std::sqrt(n / d) : std::sqrt(n);
}
// the failure-driven piece schedule of returnMapTensor, re-derived independently (the oracle side of D1)
static int failureSchedule(const Params& q, const double s0[6], const double d[6], double k0, double fr[64]) {
    double sn[6], kn, Dt[6][6];
    if (returnMapTensor1(q, s0, d, k0, true, sn, kn, Dt, false) == 0) return 0;
    double s[6], k = k0, done = 0.0, frac = 0.5; int n = 0;
    for (int i = 0; i < 6; ++i) s[i] = s0[i];
    while (done < 1.0 && n < 64) {
        const double f = (frac < 1.0 - done) ? frac : 1.0 - done; double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] * f;
        if (returnMapTensor1(q, s, dd, k, true, sn, kn, Dt, false) == 0) { for (int i = 0; i < 6; ++i) s[i] = sn[i]; k = kn; done += f; fr[n++] = f; frac = (2 * frac < 1) ? 2 * frac : 1; }
        else { frac *= 0.5; if (frac < std::ldexp(1.0, -q.maxSubIncr)) return -1; }
    }
    return n;
}
// central FD (h in {1e-8, 1e-9, 1e-10}) of the GIVEN-schedule map, effective stress and kp. Returns the best error of
// (D, dk); `smooth` = false when two step sizes disagree by > 1e-6 (the stencil crosses a yield/apex switch).
static void fdSchedule(const Params& q, const double s0[6], const double d[6], double k0, const double* fr, int n,
                       const double D[6][6], const double dk[6], double& eD, double& eK, bool& smooth) {
    eD = eK = 1e300; smooth = true; double prev[6][6]; bool havePrev = false;
    for (double h : {1e-8, 1e-9, 1e-10}) {
        double F[6][6], FK[6]; bool ok = true;
        for (int j = 0; j < 6 && ok; ++j) {
            double dp[6], dm[6], sp[6], sm[6], kp, km, J[6][6]; for (int i = 0; i < 6; ++i) dp[i] = dm[i] = d[i]; dp[j] += h; dm[j] -= h;
            ok = returnMapTensorSchedule(q, s0, dp, k0, fr, n, sp, kp, J, false) == 0 && returnMapTensorSchedule(q, s0, dm, k0, fr, n, sm, km, J, false) == 0;
            for (int i = 0; i < 6; ++i) F[i][j] = (sp[i] - sm[i]) / (2 * h);
            FK[j] = (kp - km) / (2 * h);
        }
        if (!ok) { smooth = false; continue; }
        if (havePrev && relF6(F, prev) > 1.0e-6) smooth = false;
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) prev[i][j] = F[i][j];
        havePrev = true;
        eD = std::fmin(eD, relF6(D, F)); eK = std::fmin(eK, relV6(dk, FK));
    }
}

static void run_f2_chain() {
    std::printf("[D] F2 (apeGmsh #2621): sub-step chain-rule tangent + analytic damaged tangent\n");
    Params q; q.E = 30000; q.nu = 0.2; q.fc = 30; q.ft = 3; q.e = eccentricityFromKupfer(30, 3, 1.16); q.m0 = m0Of(30, 3, q.e);
    q.qh0 = 0.3; q.Hp = 0.01; q.Df = 0.85; q.Ah = 0.08; q.Bh = 0.003; q.Ch = 2.0; q.Dh = 1e-6;
    q.flowPotential = 1; q.maxSubIncr = 10;
    // D1 — failure-driven 2- and 4-piece returns: returnMapTensor's Dtan6 and d kp/dε == the central FD of the map with
    // that schedule frozen, to 1e-7 (the pre-F2 last-piece tangent was O(1) off here). Same trials pushed through FORCED
    // {1/2,1/2} and {1/4 x4} schedules (returnMapTensorSchedule).
    unsigned s3 = 2621u; auto rn = [&]() { s3 = s3 * 1664525u + 1013904223u; return (s3 >> 8) * (1.0 / 16777216.0); };
    int n2 = 0, n4 = 0, nF = 0, skipped = 0, tested = 0; double worstD = 0, worstK = 0, worstF = 0;
    double mutBK = 1e300, mutLast = 1e300;                 // D2: how far the two mutants sit from the FD
    double kerBK = 1e300, kerLast = 1e300;                 // D2: how far the KERNEL's Dtan6 sits from each mutant
    int nMut = 0;
    struct Case { double s0[6], d[6], k0; };
    Case c4[40]; int nc4 = 0;
    for (long t = 0; t < 400000 && (n2 < 40 || n4 < 40); ++t) {
        const int mode = t % 4; double d[6];
        if (mode == 0) { for (int i = 0; i < 6; ++i) d[i] = (rn() * 2 - 1) * 3e-3; for (int i = 3; i < 6; ++i) d[i] *= 0.5; }
        else if (mode == 1) { const double v = -rn() * 8e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 2e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 1e-4; }
        else if (mode == 2) { d[0] = -rn() * 5e-3; d[1] = (rn() * 2 - 1) * 1.5e-3; d[2] = (rn() * 2 - 1) * 1.5e-3; d[3] = (rn() * 2 - 1) * 5e-4; d[4] = d[5] = 0; }
        else { const double v = rn() * 1e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 3e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 2e-4; }
        const double k0 = rn() * 1.5;
        double s0[6] = {0, 0, 0, 0, 0, 0};
        if (t % 5 == 1) { const double p = -(rn()) * 10.0; s0[0] = s0[1] = s0[2] = p; s0[3] = (rn() * 2 - 1) * 3.0; }   // pre-stressed
        double sn[6], kn, Dt[6][6], dk[6]; int sub = 0;
        if (returnMapTensor(q, s0, d, k0, true, sn, kn, Dt, true, &sub, dk) != 0 || sub == 0) continue;
        if (!((sub == 2 && n2 < 40) || (sub == 4 && n4 < 40))) continue;
        double fr[64]; const int nf = failureSchedule(q, s0, d, k0, fr);
        if (nf != sub) { ++skipped; continue; }
        double eD, eK; bool smooth;
        fdSchedule(q, s0, d, k0, fr, nf, Dt, dk, eD, eK, smooth);
        if (!smooth) { ++skipped; continue; }
        ++tested; if (sub == 2) ++n2; else { ++n4; if (nc4 < 40) { Case& c = c4[nc4++]; for (int i = 0; i < 6; ++i) { c.s0[i] = s0[i]; c.d[i] = d[i]; } c.k0 = k0; } }
        worstD = std::fmax(worstD, eD); worstK = std::fmax(worstK, eK);
        // forced schedules through the public schedule entry point
        const double h2[2] = {0.5, 0.5}, h4[4] = {0.25, 0.25, 0.25, 0.25};
        for (int f = 0; f < 2; ++f) {
            const double* fr2 = f ? h4 : h2; const int nn = f ? 4 : 2;
            double sF[6], kF, DF[6][6], dkF[6];
            if (returnMapTensorSchedule(q, s0, d, k0, fr2, nn, sF, kF, DF, true, dkF) != 0) continue;
            double e1, e2; bool sm; fdSchedule(q, s0, d, k0, fr2, nn, DF, dkF, e1, e2, sm);
            if (!sm) continue;
            ++nF; worstF = std::fmax(worstF, std::fmax(e1, e2));
        }
        // D2 — the kernel's Dtan6 is neither of two wrong builds, rebuilt here from the kernel's own pieces: the chain
        // without b (x) K, and the last piece's tangent alone (pre-F2). A kernel carrying either mutation sits at
        // distance 0 from it and fails; the FD distances show the 1e-7 D1 tolerance would catch it as well.
        {
            double Cel[6][6]; elasticC(q, Cel);
            double J[6][6] = {{0}}, Kr[6] = {0}, Jn[6][6] = {{0}}, s[6], k = k0, sp[6], kp2, Dl[6][6];
            for (int i = 0; i < 6; ++i) s[i] = s0[i];
            for (int p = 0; p < nf; ++p) {
                double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] * fr[p];
                PieceJac pj; returnMapTensor1(q, s, dd, k, true, sp, kp2, Dl, p == nf - 1, &pj);
                double T[6][6], Jm[6][6], Km[6];
                for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) T[i][j] = J[i][j] + fr[p] * Cel[i][j];
                for (int j = 0; j < 6; ++j) { double kv = pj.e * Kr[j]; for (int i = 0; i < 6; ++i) kv += pj.a[i] * T[i][j]; Km[j] = kv;
                    for (int i = 0; i < 6; ++i) { double v = 0.0; for (int m = 0; m < 6; ++m) v += pj.A[i][m] * T[m][j]; Jm[i][j] = v; } }   // no b (x) K
                for (int i = 0; i < 6; ++i) { Kr[i] = Km[i]; for (int j = 0; j < 6; ++j) J[i][j] = Jm[i][j]; }
                for (int i = 0; i < 6; ++i) s[i] = sp[i]; k = kp2;
            }
            for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Jn[i][j] = J[i][j];
            double e1, e2, e3, e4; bool sm1, sm2; double dz[6] = {0, 0, 0, 0, 0, 0};
            fdSchedule(q, s0, d, k0, fr, nf, Jn, dz, e1, e2, sm1);
            fdSchedule(q, s0, d, k0, fr, nf, Dl, dz, e3, e4, sm2);
            if (sm1 && sm2) { ++nMut; mutBK = std::fmin(mutBK, e1); mutLast = std::fmin(mutLast, e3);
                              kerBK = std::fmin(kerBK, relF6(Dt, Jn)); kerLast = std::fmin(kerLast, relF6(Dt, Dl)); }
        }
    }
    check(n2 >= 20 && n4 >= 20 && skipped <= tested && worstD < 1.0e-7 && worstK < 1.0e-7,
          "D1 failure-driven 2-/4-piece returns: chain-rule Dtan6 and d kp/de == frozen-schedule FD to 1e-7");
    std::printf("       (2-piece %d, 4-piece %d, skipped (kink in the stencil) %d; worst Dtan6 %.2e, worst dkp %.2e)\n", n2, n4, skipped, worstD, worstK);
    check(nF >= 40 && worstF < 1.0e-7, "D1 forced {1/2,1/2} and {1/4 x4} schedules: chain-rule tangent == FD to 1e-7");
    std::printf("       (forced cases %d, worst %.2e)\n", nF, worstF);
    check(nMut >= 20 && kerBK > 1.0e-5 && kerLast > 1.0e-3 && mutBK > 1.0e-5 && mutLast > 1.0e-3,
          "D2 the kernel's Dtan6 differs from a chain without b(x)K (> 1e-5) and from the last-piece tangent (> 1e-3)");
    std::printf("       (cases %d; kernel vs no-b(x)K chain >= %.2e, vs last piece >= %.2e; those mutants miss the FD by >= %.2e / %.2e)\n",
                nMut, kerBK, kerLast, mutBK, mutLast);

    // D3 — N = k +/- 1e-12 (a rounded piece count): pieces of 1/N plus the remainder piece; stress and tangent are
    // continuous in N at N = 4 (a tiny remainder piece adds nothing).
    {
        double worstS = 0, worstT = 0; int n = 0;
        for (int c = 0; c < nc4; ++c) {
            const Case& cs = c4[c];
            const double f4[4] = {0.25, 0.25, 0.25, 0.25};
            double sA[6], kA, DA[6][6], dkA[6];
            if (returnMapTensorSchedule(q, cs.s0, cs.d, cs.k0, f4, 4, sA, kA, DA, true, dkA) != 0) continue;
            bool all = true; double es = 0, et = 0;
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double N = 4.0 + sgn * 1.0e-12; const int nfull = (int)std::floor(N);
                double fr[6]; int m = 0; for (int i = 0; i < nfull; ++i) fr[m++] = 1.0 / N;
                const double rem = 1.0 - nfull / N; if (rem > 0.0) fr[m++] = rem;
                double sB[6], kB, DB[6][6], dkB[6];
                if (returnMapTensorSchedule(q, cs.s0, cs.d, cs.k0, fr, m, sB, kB, DB, true, dkB) != 0) { all = false; break; }
                for (int i = 0; i < 6; ++i) es = std::fmax(es, std::fabs(sB[i] - sA[i]) / q.fc);
                et = std::fmax(et, relF6(DB, DA));
            }
            if (!all) continue;
            ++n; worstS = std::fmax(worstS, es); worstT = std::fmax(worstT, et);
        }
        check(n >= 10 && worstS < 1.0e-9 && worstT < 1.0e-7, "D3 N = 4 +/- 1e-12 pieces: stress and chain-rule tangent continuous in N");
        std::printf("       (cases %d, worst |dsigma|/fc %.2e, worst tangent %.2e)\n", n, worstS, worstT);
    }

    // D4 — the call counter: damagedTangent runs NO return map (every mode whose old tangent micro-FD'd the map:
    // the beta_c gradient under compressive loading, ctTemper proj, tcTemper proj), and returnMap runs exactly one.
    {
        Params p = q; p.tensionLaw = 1; p.compDrive = 1; p.tcTemper = 2; p.Gf = 0.1; p.Gc = 5.0; p.lch = 50.0;
        p.epsFc = 1e-3; p.As = 2.0;
        unsigned long long inDT = 0, perMap = 0; int nmaps = 0, nDT = 0;
        for (int ct = 0; ct <= 2; ct += 2) {
            p.ctTemper = ct;
            State in; double sg[6], se[6], D[6][6];
            // crack in tension, then compress with shear (tensile + compressive histories, sub-incremented steps)
            for (int k = 1; k <= 40; ++k) {
                double e[6] = {k <= 15 ? 1.2e-5 * k : 1.8e-4 - 1.5e-4 * (k - 15), -2e-5 * k, 0.0, 4e-5 * k, 0, 0};
                State o;
                KernelCallCounters before = kernelCallCounters();
                returnMap(p, e, in, o, sg, se, D, true, -1.0, true);
                { const unsigned long long dm = kernelCallCounters().returnMapTensor - before.returnMapTensor; if (dm > perMap) perMap = dm; }
                ++nmaps;
                double deps[6], sn[6], kn, Ce[6][6], dk[6];
                for (int i = 0; i < 6; ++i) deps[i] = e[i] - in.eps[i];
                if (returnMapTensor(p, in.sigEff, deps, in.kp, true, sn, kn, Ce, true, nullptr, dk) == 0) {
                    KernelCallCounters b2 = kernelCallCounters();
                    double Dd[6][6]; damagedTangent(p, in, sn, e, kn, Ce, dk, Dd);
                    inDT += (kernelCallCounters().returnMapTensor - b2.returnMapTensor) + (kernelCallCounters().pieceReturns - b2.pieceReturns);
                    ++nDT;
                }
                in = o;
            }
        }
        check(nDT > 0 && inDT == 0 && perMap == 1, "D4 call counter: damagedTangent runs 0 return maps; returnMap runs exactly 1");
        std::printf("       (damagedTangent calls %d: return maps inside %llu; returnMap calls %d: max maps per call %llu)\n", nDT, inDT, nmaps, perMap);
    }

    // D5 — the closed-form proj-weight gradients == a central FD of the closed-form weights themselves.
    {
        Params p = q; p.ctTemper = 2; p.tcTemper = 2;
        unsigned s5 = 77u; auto r5 = [&]() { s5 = s5 * 1664525u + 1013904223u; return (s5 >> 8) * (1.0 / 16777216.0); };
        double worstC = 0, worstT = 0; int nC = 0, nT = 0;
        for (int t = 0; t < 400; ++t) {
            double sb[6], dp[6];
            sb[0] = (r5() * 2 - 1) * 20.0; sb[1] = p.ft * (0.002 + 0.045 * r5()); sb[2] = (r5() * 2 - 1) * 5.0;   // one principal in the tcPhi band
            sb[3] = (r5() * 2 - 1) * 0.02; sb[4] = (r5() * 2 - 1) * 0.02; sb[5] = (r5() * 2 - 1) * 0.5;
            for (int i = 0; i < 6; ++i) dp[i] = (r5() * 2 - 1) * 1e-4;
            double A[3][3], w[3], V[3][3]; voigtToMat(sb, A); eig3sym(A, w, V);
            double gap = 1e300; for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) gap = std::fmin(gap, std::fabs(w[a] - w[b]));
            if (gap < 0.5) continue;                       // keep the FD off near-repeated eigenvalues (1/gap eigenvector sensitivity)
            double gs[6], gd[6], hs[6], hd[6];
            compressiveDamageWeightGrad(p, dp, w, V, gs, gd);
            tensileDamageWeightGrad(p, dp, w, V, hs, hd);
            double fsC[6], fdC[6], fsT[6], fdT[6]; bool kink = false;
            for (int k = 0; k < 6; ++k) {
                const double hsg = 1e-6 * p.fc, hdp = 1e-6 * 1e-4;
                double a[6], b[6]; for (int i = 0; i < 6; ++i) { a[i] = sb[i]; b[i] = sb[i]; } a[k] += hsg; b[k] -= hsg;
                double Aa[3][3], wa[3], Va[3][3], Ab[3][3], wb[3], Vb[3][3]; voigtToMat(a, Aa); eig3sym(Aa, wa, Va); voigtToMat(b, Ab); eig3sym(Ab, wb, Vb);
                for (int c = 0; c < 3; ++c) { const double lo = TC_DEAD * p.ft, hi = TC_BAND * p.ft;
                    const double fl = 1e-6 * p.ft; if ((wa[c] - lo) * (wb[c] - lo) <= 0 || (wa[c] - hi) * (wb[c] - hi) <= 0 || wa[c] * wb[c] <= 0 || (wa[c] - fl) * (wb[c] - fl) <= 0) kink = true; }
                fsC[k] = (compressiveDamageWeight(p, dp, wa, Va) - compressiveDamageWeight(p, dp, wb, Vb)) / (2 * hsg);
                fsT[k] = (tensileDamageWeight(p, 0.0, dp, wa, Va) - tensileDamageWeight(p, 0.0, dp, wb, Vb)) / (2 * hsg);
                double c1[6], c2[6]; for (int i = 0; i < 6; ++i) { c1[i] = dp[i]; c2[i] = dp[i]; } c1[k] += hdp; c2[k] -= hdp;
                fdC[k] = (compressiveDamageWeight(p, c1, w, V) - compressiveDamageWeight(p, c2, w, V)) / (2 * hdp);
                fdT[k] = (tensileDamageWeight(p, 0.0, c1, w, V) - tensileDamageWeight(p, 0.0, c2, w, V)) / (2 * hdp);
            }
            if (kink) continue;
            double vC[12], vCf[12], vT[12], vTf[12];
            for (int k = 0; k < 6; ++k) { vC[k] = gs[k] * p.fc; vCf[k] = fsC[k] * p.fc; vC[6 + k] = gd[k] * 1e-4; vCf[6 + k] = fdC[k] * 1e-4;
                                         vT[k] = hs[k] * p.fc; vTf[k] = fsT[k] * p.fc; vT[6 + k] = hd[k] * 1e-4; vTf[6 + k] = fdT[k] * 1e-4; }
            double nC2 = 0, dC2 = 0, nT2 = 0, dT2 = 0;
            for (int k = 0; k < 12; ++k) { nC2 += (vC[k] - vCf[k]) * (vC[k] - vCf[k]); dC2 += vCf[k] * vCf[k]; nT2 += (vT[k] - vTf[k]) * (vT[k] - vTf[k]); dT2 += vTf[k] * vTf[k]; }
            if (dC2 > 0) { ++nC; worstC = std::fmax(worstC, std::sqrt(nC2 / dC2)); }
            if (dT2 > 0) { ++nT; worstT = std::fmax(worstT, std::sqrt(nT2 / dT2)); }
        }
        check(nC >= 50 && nT >= 50 && worstC < 1.0e-6 && worstT < 1.0e-6,
              "D5 closed-form tcTemper/ctTemper proj-weight gradients == FD of the weights to 1e-6");
        std::printf("       (tc cases %d worst %.2e; ct cases %d worst %.2e)\n", nC, worstC, nT, worstT);
    }
}

int main(int argc, char** argv) {
    run_identities();
    const char* fix = argc > 1 ? argv[1] : "tests/_testbed/concrete3d_oracle_fixture.txt";
    run_oracle_dump(fix);
    run_robustness();
    run_f2_chain();
    std::printf(fails ? "\nKERNEL CHECK: %d FAIL\n" : "\nKERNEL CHECK: ALL PASS\n", fails);
    return fails ? 1 : 0;
}
