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
        fh >> in.sigtMax >> in.sigcMax >> in.wt >> in.wc >> mp.omegaDead;   // P2g monotone-drive history (8-field line)
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
        Params mpf = mp; mpf.subIncrForceN = returnMapPieces(mp, in, strain);   // both FD legs use the CENTRAL n (the map is discontinuous at n boundaries)
        for (int j = 0; j < 6; ++j) {
            const double dd = 1.0e-6 * (std::fabs(deps[j]) + base);
            double sp[6], sm[6], se[6], junk[6][6]; State o2;
            double strp[6], strm[6];
            for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
            strp[j] += dd; strm[j] -= dd;
            returnMap(mpf, strp, in, o2, sp, se, junk, false, -1.0, true);
            returnMap(mpf, strm, in, o2, sm, se, junk, false, -1.0, true);
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
        fh >> in.sigtMax >> in.sigcMax >> in.wt >> in.wc >> mp.omegaDead;   // + committed omega_t/omega_c + dead threshold
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
        fh >> in.sigtMax >> in.sigcMax >> in.wt >> in.wc >> mp.omegaDead;   // P2g monotone-drive history (8-field line)
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
            int nsteps; fh >> mp.flowPotential >> mp.maxSubIncr >> mp.subIncrMode >> nsteps;   // subIncrMode: 0 deterministic (default), 1 adaptive
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
               >> in.sigtMax >> in.sigcMax >> in.eqc >> in.etPrev >> in.wt >> in.wc >> mp.omegaDead;
            double deps[6], sigO[6];
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 6; ++i) fh >> sigO[i];
            double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
            double sigC[6], sigEff[6], Da[6][6];
            returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, true);
            double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
            const double base = mp.fc / mp.E;
            Params mpf = mp; mpf.subIncrForceN = returnMapPieces(mp, in, strain);   // both FD legs use the CENTRAL n (the map is discontinuous at n boundaries)
            auto fdAt = [&](double fac, double Cn[6][6]) {
                for (int j = 0; j < 6; ++j) {
                    const double dd = fac * (std::fabs(deps[j]) + base);
                    double sp[6], sm[6], se[6], junk[6][6]; State o2; double strp[6], strm[6];
                    for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
                    strp[j] += dd; strm[j] -= dd;
                    returnMap(mpf, strp, in, o2, sp, se, junk, false, -1.0, true);
                    returnMap(mpf, strm, in, o2, sm, se, junk, false, -1.0, true);
                    for (int i = 0; i < 6; ++i) Cn[i][j] = (sp[i] - sm[i]) / (2.0 * dd);
                }
            };
            auto relTo = [](const double A[6][6], const double R[6][6]) {
                double nd = 0, nn = 0;
                for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) { const double df = A[a][b] - R[a][b]; nd += df * df; nn += R[a][b] * R[a][b]; }
                return std::sqrt(nd / nn);
            };
            // The FD step: 1e-6 of the strain scale, except on the omega_c ONSET step (dmgc_uniax_onset, where omega_c first
            // leaves 0). There the nominal stress carries the kappa_p Newton's round-off (~1e-10 in kappa_p, ~2.5e-9 MPa through
            // omega_c, measured on a 1e-10 strain grid), so at 1e-6 the FD is noise-limited: the sweep 1e-4 / 3e-5 / 1e-5 /
            // 3e-6 / 1e-6 reads 3.04e-5 / 2.99e-5 / 3.09e-5 / 4.48e-5 / 1.89e-5, UCRT and libmmt alike. It is not a kink: no
            // flag flips inside the stencil (the onset is >= 1.3e-5 strain away, 1e4 steps). From 1e-4 to 1e-5 the FD has
            // converged (Richardson at 1e-4: 3.08e-5), so the onset step is checked against the 1e-4 FD, guarded by the 1e-5 FD
            // agreeing with it (the plateau), with the same 5e-5 bound. That plateau first read 3.04e-5, a real residual: the
            // kappa row of the hardening return was 1.6e-3 off on this uniaxial state (principalJacobian's micro-FD of the
            // R4 Lode factor on the compressive meridian, D12) and entered through d omega_c. With the closed form the onset
            // step reads 2.3e-6 (sweep 1e-4 / 3e-5 / 1e-5 / 3e-6 / 1e-6: 2.3e-6 / 6.0e-7 / 1.0e-6 / 1.7e-5 / 1.5e-5).
            const bool onset = (label == "dmgc_uniax_onset");
            double Cn[6][6]; fdAt(onset ? 1.0e-4 : 1.0e-6, Cn);
            double plateau = 0.0;
            if (onset) { double C5[6][6]; fdAt(1.0e-5, C5); plateau = relTo(C5, Cn); }
            const double rel = relTo(Da, Cn);
            worst_dmg = std::fmax(worst_dmg, maxs); worst_dtan = std::fmax(worst_dtan, rel);
            // tangent tolerance 5e-5 for EVERY record (direct returns and sub-incremented chains alike): the reported tangent of an
            // n-piece chain is the chain's own algorithmic tangent (accumulated forward incl. the kappa history, C8), and the FD
            // reference pins n on both legs (mpf), so a chain is as consistent as a direct return. (It was loosened to 5e-2 for
            // maxSubIncr > 0 while the last piece's tangent was reported: 6/41/121 % off at n = 2/4/64 in compression.)
            const double tanTol = 5.0e-5;
            const bool ok = maxs < 1.0e-6 && rel < tanTol && plateau < 1.0e-5;
            if (!ok) ++fails;
            if (onset) std::printf("  %-24s FD step 1e-4; plateau |FD(1e-5) - FD(1e-4)| / |FD(1e-4)| = %.2e (< 1e-5)\n", label.c_str(), plateau);
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
               >> in.sigtMax >> in.sigcMax >> in.eqc >> in.etPrev >> in.wt >> in.wc >> mp.omegaDead;
            double deps[6], sigO[6];
            for (int i = 0; i < 6; ++i) fh >> deps[i];
            for (int i = 0; i < 6; ++i) fh >> sigO[i];
            double strain[6]; for (int i = 0; i < 6; ++i) strain[i] = in.eps[i] + deps[i];
            double sigC[6], sigEff[6], Da[6][6];
            returnMap(mp, strain, in, out, sigC, sigEff, Da, true, -1.0, true);
            double maxs = 0; for (int i = 0; i < 6; ++i) maxs = std::fmax(maxs, std::fabs(sigC[i] - sigO[i]));
            double Cn[6][6]; const double base = mp.fc / mp.E;
    Params mpf = mp; mpf.subIncrForceN = returnMapPieces(mp, in, strain);   // both FD legs use the CENTRAL n (the map is discontinuous at n boundaries)
            for (int j = 0; j < 6; ++j) {
                const double dd = 1.0e-6 * (std::fabs(deps[j]) + base);
                double sp[6], sm[6], se[6], junk[6][6]; State o2; double strp[6], strm[6];
                for (int i = 0; i < 6; ++i) { strp[i] = strain[i]; strm[i] = strain[i]; }
                strp[j] += dd; strm[j] -= dd;
                returnMap(mpf, strp, in, o2, sp, se, junk, false, -1.0, true);
                returnMap(mpf, strm, in, o2, sm, se, junk, false, -1.0, true);
                for (int i = 0; i < 6; ++i) Cn[i][j] = (sp[i] - sm[i]) / (2.0 * dd);
            }
            double nd = 0, nn = 0;
            for (int A = 0; A < 6; ++A) for (int B = 0; B < 6; ++B) { const double df = Da[A][B] - Cn[A][B]; nd += df * df; nn += Cn[A][B] * Cn[A][B]; }
            const double rel = std::sqrt(nd / nn);
            worst_dmg = std::fmax(worst_dmg, maxs); worst_dtan = std::fmax(worst_dtan, rel);
            // tangent tolerance 5e-5 for EVERY record (direct returns and sub-incremented chains alike): the reported tangent of an
            // n-piece chain is the chain's own algorithmic tangent (accumulated forward incl. the kappa history, C8), and the FD
            // reference pins n on both legs (mpf), so a chain is as consistent as a direct return. (It was loosened to 5e-2 for
            // maxSubIncr > 0 while the last piece's tangent was reported: 6/41/121 % off at n = 2/4/64 in compression.)
        const double tanTol = 5.0e-5;
            const bool ok = maxs < 1.0e-6 && rel < tanTol;
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

    // C8 — CHAIN TANGENT of the deterministic sub-incrementation (WP concrete3d-hang-diagnosis, review M1). The reported
    // tangent of an n-piece chain is the chain's own algorithmic tangent, accumulated forward through the pieces
    // (G_p = (Dt C^-1) G_(p-1) + Sk (x) g_(p-1) + Dt/n, g_p = R G_(p-1) + Kk g_(p-1) + R C/n; kappa history included),
    // NOT the last piece's (which the review measured 6 % / 41 % / 121 % off at n = 2 / 4 / 64 in compression). Reference:
    // central differences of the SAME chain with n PINNED on both legs (Params::subIncrForceN; the map is discontinuous
    // where n changes). Virgin state, the reviewer's compression direction, and the tensile direction; every state whose
    // ladder did not escalate (subInfo == n). Bound 1e-6 (measured 1e-8..1e-10, the floor is the FD truncation and the
    // scalar Lode-angle central differences inside the principal Jacobian).
    {
        Params q = mp; q.Gf = 0.1; q.Gc = 5.0; q.Df = 0.85; q.As = 2.0; q.qh0 = 0.3; q.Hp = 0.01; q.lch = 50; q.lch_ref = 50;
        q.tensionLaw = 1; q.flowPotential = 1; q.compDrive = 1; q.tcTemper = 2; q.maxSubIncr = 10; q.epsFc = 5.0 / (30 * 50);
        const double dirs[2][6] = {{-1.0, 0.25, 0.25, 0.2, 0, 0}, {1.0, -0.2, -0.2, 0, 0, 0}};
        const int targets[6] = {1, 2, 4, 9, 26, 64};
        double worst = 0.0; int nchecked = 0; bool all = true;
        for (int dc = 0; dc < 2; ++dc) {
            for (int ti = 0; ti < 6; ++ti) {
                const int target = targets[ti];
                double z[6] = {0, 0, 0, 0, 0, 0}, e[6] = {0, 0, 0, 0, 0, 0}; bool found = false;
                for (double sc = 2e-4; sc < 5e-2 && !found; sc *= 1.004) {
                    double d[6]; for (int a = 0; a < 6; ++a) d[a] = sc * dirs[dc][a];
                    if (detPieces(q, z, d, 0.0) != target) continue;
                    double sg[6], kk, DD[6][6]; int inf = -9;
                    returnMapTensor(q, z, d, 0.0, true, sg, kk, DD, false, &inf);
                    if (inf == (target == 1 ? 0 : target)) { for (int a = 0; a < 6; ++a) e[a] = d[a]; found = true; }
                }
                if (!found) continue;                       // (tension: n = 1..4 always escalate the ladder near first cracking)
                double sg[6], kk, D[6][6]; int inf;
                returnMapTensor(q, z, e, 0.0, true, sg, kk, D, true, &inf);
                Params qf = q; qf.subIncrForceN = target;
                double nd = 0, nn = 0; bool straddle = false;
                for (int j = 0; j < 6; ++j) {
                    const double h = (target == 1) ? 1e-8 : 1e-10; double ep[6], em[6];
                    for (int a = 0; a < 6; ++a) { ep[a] = e[a]; em[a] = e[a]; }
                    ep[j] += h; em[j] -= h;
                    double s1[6], s2[6], k1, k2, Dd[6][6]; int i1 = -9, i2 = -9;
                    returnMapTensor(qf, z, ep, 0.0, true, s1, k1, Dd, false, &i1); returnMapTensor(qf, z, em, 0.0, true, s2, k2, Dd, false, &i2);
                    if (i1 != inf || i2 != inf) straddle = true;   // a leg escalated the ladder: the map jumps within 1e-10 here
                    for (int a = 0; a < 6; ++a) { const double f = (s1[a] - s2[a]) / (2 * h); nd += (D[a][j] - f) * (D[a][j] - f); nn += f * f; }
                }
                const double rel = std::sqrt(nd / nn);
                worst = std::fmax(worst, rel); ++nchecked; if (!(rel < 1.0e-6)) all = false;
                std::printf("       %s n=%2d (subInfo %2d): ||T - D_fd||/||D_fd|| = %.2e%s\n", dc == 0 ? "compression" : "tension    ", target, inf, rel, straddle ? " (straddle)" : "");
            }
        }
        check(nchecked >= 12 && all, "chain tangent == FD of the same n-piece chain, n = 1..64, compression + tension (rel < 1e-6)");
        std::printf("       (checked %d states, worst ||T - D_fd||/||D_fd|| = %.2e)\n", nchecked, worst);
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

// ---- (D) tension return (r4 rescue of #877; review 3 MAJOR 1 / probe6 / probe5d) ---------------
// The wrapper defaults of review 3 (fc 30, ft 3, E 30000, Gf 0.1, lch 50, bilinear, cdpm2 flow + drive, tcTemper proj,
// maxSubIncr 10, deterministic). Before the r4 port the DIRECT return failed on 645 / 904 / 1953 of 2000 points in a
// 2e-6-wide window at lateral -0.2, x = 1.5 / 3 / 6 eps_cr, and one window at lateral +0.3, x = 3 returned two
// different admissible roots (vertex kp 23 1984 times, regular kp 113 16 times).
static const int D4_MAX_JUMPS = 20;   // measured 14 on the r4 kernel (11 kappa_p-up, 3 kappa_p-down; headroom for platform round-off)
static Params tensionReturnParams() {
    Params q; q.E = 30000; q.nu = 0.2; q.fc = 30; q.ft = 3; q.e = eccentricityFromKupfer(30, 3, 1.16);
    q.m0 = m0Of(q.fc, q.ft, q.e); q.Df = 0.85; q.As = 2.0; q.qh0 = 0.3; q.Hp = 0.01; q.Gf = 0.1; q.Gc = 5.0;
    q.lch = 50; q.lch_ref = 50; q.tensionLaw = 1; q.flowPotential = 1; q.compDrive = 1; q.tcTemper = 2;
    q.maxSubIncr = 10; q.epsFc = 5.0 / (30 * 50); return q;
}
static void run_tension_return() {
    std::printf("[D] tension return (review 3 of #877, MAJOR 1)\n");
    const double z[6] = {0, 0, 0, 0, 0, 0};
    for (int lode = 0; lode <= 1; ++lode) {
        Params q = tensionReturnParams(); q.lodeFactor = lode;
        const double ecr = q.ft / q.E;
        const double lats[3] = {-0.2, 0.0, 0.3}, xs[3] = {1.5, 3.0, 6.0};
        int fails = 0, flips = 0, multi = 0; double worstJump = 0.0;
        for (double lat : lats) for (double x : xs) {
            const int N = 400; int prevInf = -99; double kmin = 1e300, kmax = -1e300, prev[6]; bool hp = false;
            for (int i = 0; i < N; ++i) {
                const double sc = x * ecr * (1.0 + 2e-6 * (double)i / N);
                const double d[6] = {sc, lat * sc, lat * sc, 0, 0, 0};
                double s[6], k, D[6][6];
                if (returnMapTensor1(q, z, d, 0.0, true, s, k, D, false) != 0) ++fails;
                else { kmin = std::fmin(kmin, k); kmax = std::fmax(kmax, k); }
                double s2[6], k2, D2[6][6]; int inf = -9;
                returnMapTensor(q, z, d, 0.0, true, s2, k2, D2, false, &inf);
                if (prevInf != -99 && inf != prevInf) ++flips; prevInf = inf;
                if (hp) { double j = 0; for (int a = 0; a < 6; ++a) j += (s2[a] - prev[a]) * (s2[a] - prev[a]); worstJump = std::fmax(worstJump, std::sqrt(j)); }
                for (int a = 0; a < 6; ++a) prev[a] = s2[a]; hp = true;
            }
            if (kmax - kmin > 1e-3 * (1.0 + kmax)) ++multi;   // two roots in one 2e-6 window
        }
        std::printf("       lodeFactor %d: direct-return failures %d / 3600, ladder flips %d, windows with two roots %d, max |dsig_eff| %.2e MPa\n",
                    lode, fails, flips, multi, worstJump);
        check(fails == 0 && flips == 0 && multi == 0 && worstJump < 1e-6,
              lode == 0 ? "tension return: no direct-return failure, no ladder flip, one root per 2e-6 window (lodeFactor 0)"
                        : "tension return: no direct-return failure, no ladder flip, one root per 2e-6 window (lodeFactor 1)");
    }
    // D2 — FULL-RESIDUAL HONESTY: every converged regular hardening return solves the 4-equation system (R1, R2 flow, R3
    // yield, R4 hardening), not only f = 0 (db02ccc70's nested fallback was accepted on f alone after 5 unchecked polish
    // steps). Random tension-dominated first-crack trials, virgin and preloaded; the vertex returns are checked on f.
    {
        Params q = tensionReturnParams();
        const double ecr = q.ft / q.E, z6[6] = {0, 0, 0, 0, 0, 0};
        unsigned long long sd = 777; auto rnd = [&]() { sd = sd * 6364136223846793005ULL + 1442695040888963407ULL; return ((sd >> 11) * (1.0 / 9007199254740992.0)); };
        int n = 0, conv = 0, bad = 0, unconv = 0; double worst = 0.0;
        for (int t = 0; t < 6000; ++t) {
            const double kpn = (t % 2) ? 0.0 : 5.0 * rnd(), x = 1.2 + 30 * rnd();
            const double d[6] = {x * ecr, (-0.3 + 0.6 * rnd()) * x * ecr, (-0.3 + 0.6 * rnd()) * x * ecr,
                                 0.3 * (rnd() - 0.5) * x * ecr, 0.3 * (rnd() - 0.5) * x * ecr, 0.3 * (rnd() - 0.5) * x * ecr};
            double sig_tr[6]; elasticPredTensor(z6, d, q, sig_tr);
            double A[3][3], w[3], V[3][3]; voigtToMat(sig_tr, A); eig3sym(A, w, V);
            PrincipalResult pr = returnMapHardening(w, q, kpn);
            if (!pr.plastic) continue; ++n;
            if (!pr.converged) { ++unconv; continue; } ++conv;
            bool ok = pr.dlam >= -1e-12 && pr.kp >= kpn - 1e-12 && std::fabs(pr.f_after) < 3.1e-6;
            if (!pr.apex) {
                double sv[6] = {w[0], w[1], w[2], 0, 0, 0}, xi_tr, rho_tr, th_tr; invariants(sv, xi_tr, rho_tr, th_tr);
                double Rr[4];
                cdpm2HardeningResidual(pr.xi, pr.rho, pr.dlam, pr.kp, xi_tr, rho_tr, kpn, lodeR(th_tr, q.e),
                                       hardeningLodeFactor(rho_tr, th_tr, q), q, Rr, nullptr, q.flowPotential != 1);
                const double rel = std::fmax(std::fmax(std::fabs(Rr[0]) / q.fc, std::fabs(Rr[1]) / q.fc),
                                             std::fmax(std::fabs(Rr[2]), std::fabs(Rr[3]) / (1.0 + pr.kp)));
                worst = std::fmax(worst, rel);
                if (!(rel < NESTED_RES_TOL)) ok = false;
            }
            if (!ok) ++bad;
        }
        std::printf("       %d plastic trials: %d converged, %d honest failures, %d dishonest; worst scaled residual %.2e\n", n, conv, unconv, bad, worst);
        check(bad == 0 && unconv == 0, "every converged tension return solves the full residual (R1, R2, R3, R4), none fails");
    }
    // D3 -- the nested fallback's acceptance gate checks the FULL residual: a nested candidate is accepted as returned,
    // and the same candidate with dlam moved 1 % (still admissible, R1/R2/R4 violated) is REJECTED when not polished;
    // polishing it back to the root makes it acceptable again. (A gate on f alone would accept the perturbed one: the
    // perturbation leaves (xi, rho, kp) and so f untouched.)
    {
        Params q = tensionReturnParams();
        const double ecr = q.ft / q.E, z6[6] = {0, 0, 0, 0, 0, 0}, sc = 3.0 * ecr;
        const double d[6] = {sc, -0.2 * sc, -0.2 * sc, 0, 0, 0};
        double sig_tr[6]; elasticPredTensor(z6, d, q, sig_tr);
        double A[3][3], w[3], V[3][3]; voigtToMat(sig_tr, A); eig3sym(A, w, V);
        double sv[6] = {w[0], w[1], w[2], 0, 0, 0}, xi_tr, rho_tr, th_tr; invariants(sv, xi_tr, rho_tr, th_tr);
        const double r = lodeR(th_tr, q.e), c2 = hardeningLodeFactor(rho_tr, th_tr, q);
        double x[4];
        const bool got = returnMapNested(xi_tr, rho_tr, r, c2, q, 0.0, x[0], x[1], x[2], x[3]);
        double y[4] = {x[0], x[1], x[2] * 1.01, x[3]}, z[4] = {x[0], x[1], x[2] * 1.01, x[3]};
        const double f_y = yfInvHard(y[0], y[1], r, y[3], q);
        const bool accRaw = got && polishHardeningRoot(x, xi_tr, rho_tr, 0.0, r, c2, q, 1.0e-11, 0);
        const bool accBad = polishHardeningRoot(y, xi_tr, rho_tr, 0.0, r, c2, q, 1.0e-11, 0);
        const bool accFix = polishHardeningRoot(z, xi_tr, rho_tr, 0.0, r, c2, q, 1.0e-11, 8);
        std::printf("       nested candidate kp %.4g: accepted %d; dlam +1%% (|f| = %.1e) accepted unpolished %d, polished %d\n",
                    x[3], (int)accRaw, std::fabs(f_y), (int)accBad, (int)accFix);
        check(got && accRaw && !accBad && accFix, "nested fallback gate: full residual, not f alone (perturbed-flow candidate rejected)");
    }
    // D4 -- MAGNITUDE SWEEP (r4 review of #877, item 2): 40 random tension rays, virgin, x = 1..41 eps_cr in steps of
    // 0.01 eps_cr (lateral -0.4..0.4, shear +-0.2 of the axial strain). The first r4 cut ran the plain Newton first and
    // the nested solve as the rescue, and the two converged to DIFFERENT genuine roots on neighbouring trials (the
    // review: 351 neighbour jumps > 1 MPa, up to 102 MPa, flipping back and forth). With the fixed root order (vertex,
    // nested with the smallest-kappa_p root, then the Newtons) what is left are ONE-SIDED jumps where the least-hardening
    // branch ends (a fold of the backward-Euler system) and the vertex/regular switch of lodeFactor 0. Gate: no failure,
    // no ray jumps > 1 MPa more than twice (no flip-flopping: the principled check), and the total stays near the measured
    // count (14, worst 23 MPa; limit 20). The plain-Newton-first order fails it (110 jumps, 10 on one ray, 99.7 MPa).
    {
        Params q = tensionReturnParams();
        const double ecr = q.ft / q.E, z6[6] = {0, 0, 0, 0, 0, 0};
        unsigned long long sd = 4242; auto U = [&]() { sd = sd * 2862933555777941757ULL + 3037000493ULL; return (double)(sd >> 11) / 9007199254740992.0; };
        int total = 0, worstRay = 0, fails = 0; double worstJump = 0.0;
        for (int dir = 0; dir < 40; ++dir) {
            const double d[6] = {1.0, -0.4 + 0.8 * U(), -0.4 + 0.8 * U(), 0.4 * (U() - 0.5), 0.4 * (U() - 0.5), 0.4 * (U() - 0.5)};
            double prev[6]; bool hp = false; int nray = 0; const int N = 4000;
            for (int i = 0; i <= N; ++i) {
                const double x = 1.0 + 40.0 * i / N;
                double de[6]; for (int a = 0; a < 6; ++a) de[a] = d[a] * x * ecr;
                double sg[6], k, D[6][6];
                if (returnMapTensor1(q, z6, de, 0.0, true, sg, k, D, false) != 0) { ++fails; hp = false; continue; }
                if (hp) {
                    double j = 0; for (int a = 0; a < 6; ++a) j += (sg[a] - prev[a]) * (sg[a] - prev[a]); j = std::sqrt(j);
                    if (j > 1.0) { ++nray; ++total; worstJump = std::fmax(worstJump, j); }
                }
                for (int a = 0; a < 6; ++a) prev[a] = sg[a]; hp = true;
            }
            worstRay = nray > worstRay ? nray : worstRay;
        }
        std::printf("       40 rays x 4001 points: direct-return failures %d, neighbour jumps > 1 MPa %d (worst %.1f MPa), most on one ray %d\n",
                    fails, total, worstJump, worstRay);
        check(fails == 0 && worstRay <= 2 && total <= D4_MAX_JUMPS,
              "tension magnitude sweep: no failure, at most 2 one-sided jumps per ray, total <= measured");
    }
}

// ---- (E) apeGmsh #2621 F2: the sub-step chain-rule tangent + the analytic damaged tangent ---------------------
// Ported from the F2 branch (written on d14fb6958 as D1-D9). F2's chain checks D1-D5 are E1-E5 here, because #877's
// tension-return checks own D1-D4; F2's damaged-tangent checks keep their names D6-D9 (the names apeGmsh #2621 uses),
// and D10-D11 cover the two ladruno-only chain gaps (the deterministic mode's adaptive rescue, the tension-dead cut).
// On #877 the default sub-incrementation is DETERMINISTIC (Params::subIncrMode = 0, n equal pieces from the trial
// overshoot, see returnMapTensorDet); the failure-driven halving path F2 was written against is the ADAPTIVE mode
// (subIncrMode = 1, and the deterministic mode's rescue). E1-E3 run the adaptive mode (-subIncr adaptive), E1d the
// deterministic one; E4 and D6-D9 run the shipped default (deterministic) through returnMap.
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
// the failure-driven piece schedule of the ADAPTIVE returnMapTensor (subIncrMode = 1), re-derived independently (the
// oracle side of E1): #877's rule, doubling back to at most 2x the size that just succeeded, maxSubAttempts budget
static int failureSchedule(const Params& q, const double s0[6], const double d[6], double k0, double fr[64]) {
    double sn[6], kn, Dt[6][6];
    if (returnMapTensor1(q, s0, d, k0, true, sn, kn, Dt, false) == 0) return 0;
    const int maxAttempts = (q.maxSubAttempts > 0) ? q.maxSubAttempts : 64;
    double s[6], k = k0, done = 0.0, frac = 0.5; int n = 0, attempts = 0;
    for (int i = 0; i < 6; ++i) s[i] = s0[i];
    while (done < 1.0 && n < 64) {
        if (attempts >= maxAttempts) return -1;
        ++attempts;
        const double f = (frac < 1.0 - done) ? frac : 1.0 - done; double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] * f;
        if (returnMapTensor1(q, s, dd, k, true, sn, kn, Dt, false) == 0) { for (int i = 0; i < 6; ++i) s[i] = sn[i]; k = kn; done += f; fr[n++] = f; frac = (2 * f < 1) ? 2 * f : 1; }
        else { frac *= 0.5; if (frac < std::ldexp(1.0, -q.maxSubIncr)) return -1; }
    }
    return done < 1.0 ? -1 : n;
}
// the schedule the CURRENT mode integrated at a centre that reported subInfo = sub >= 2: adaptive -> the failure-driven
// fractions; deterministic -> sub equal pieces (the ladder level that succeeded; a rescue reports the adaptive count,
// which the caller tells apart by re-deriving the adaptive schedule)
static int modeSchedule(const Params& q, const double s0[6], const double d[6], double k0, int sub, double fr[64]) {
    if (q.subIncrMode == 1) return failureSchedule(q, s0, d, k0, fr);
    if (sub < 2 || sub > 64) return -1;
    for (int i = 0; i < sub; ++i) fr[i] = 1.0 / sub;
    return sub;
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
    std::printf("[E] F2 (apeGmsh #2621): sub-step chain-rule tangent + analytic damaged tangent\n");
    Params q; q.E = 30000; q.nu = 0.2; q.fc = 30; q.ft = 3; q.e = eccentricityFromKupfer(30, 3, 1.16); q.m0 = m0Of(30, 3, q.e);
    q.qh0 = 0.3; q.Hp = 0.01; q.Df = 0.85; q.Ah = 0.08; q.Bh = 0.003; q.Ch = 2.0; q.Dh = 1e-6;
    q.flowPotential = 1; q.maxSubIncr = 10;
    q.subIncrMode = 1;                                     // E1-E3: the ADAPTIVE (failure-driven) path
    // E1 — failure-driven 2- and 4-piece returns: returnMapTensor's Dtan6 and d kp/dε == the central FD of the map with
    // that schedule frozen, to 1e-7 (the pre-F2 last-piece tangent was O(1) off here). Same trials pushed through FORCED
    // {1/2,1/2} and {1/4 x4} schedules (returnMapTensorSchedule).
    // On #877 the cdpm2-flow direct return rarely fails at F2's increment sizes (0 of 40000 trials at scale 1; ~10 % at
    // scale 20-30), while the legacy flow still fails ~3 %: the trials alternate legacy flow at F2's scale and cdpm2 flow
    // at 8x it, so both flows reach the adaptive sub-stepped path. (At 20x the cdpm2 trials land at kappa_p ~ 1e4-8e4 with
    // a trial overshoot f_tr ~ 1e4, where a piece of the frozen schedule sits on the edge of its Newton's convergence: the
    // centre returns on that schedule but a stencil leg 1e-8, 1e-10 or even 1e-13 away fails (rc 2), so the FD reference is
    // unusable and the sample is skipped -- 82 skipped vs 80 tested under UCRT, 76 under libmmt. At 8x: 3 skipped, the
    // same 3 under both, 37 of the 80 tested on the cdpm2 flow.)
    unsigned s3 = 2621u; auto rn = [&]() { s3 = s3 * 1664525u + 1013904223u; return (s3 >> 8) * (1.0 / 16777216.0); };
    int n2 = 0, n4 = 0, nF = 0, skipped = 0, tested = 0; double worstD = 0, worstK = 0, worstF = 0;
    double mutBK = 1e300, mutLast = 1e300;                 // E2: how far the two mutants sit from the FD
    double kerBK = 1e300, kerLast = 1e300;                 // E2: how far the KERNEL's Dtan6 sits from each mutant
    int nMut = 0;
    struct Case { double s0[6], d[6], k0; int fl; };
    Case c4[40]; int nc4 = 0;
    Params qLeg = q; qLeg.flowPotential = 0;
    const Params qCd = q;
    for (long t = 0; t < 400000 && (n2 < 40 || n4 < 40); ++t) {
        const bool cdT = ((t >> 2) & 1) != 0;
        const Params& q = cdT ? qCd : qLeg;               // (shadows the base parameters for this trial)
        const double scT = cdT ? 8.0 : 1.0;
        const int mode = t % 4; double d[6];
        if (mode == 0) { for (int i = 0; i < 6; ++i) d[i] = (rn() * 2 - 1) * 3e-3; for (int i = 3; i < 6; ++i) d[i] *= 0.5; }
        else if (mode == 1) { const double v = -rn() * 8e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 2e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 1e-4; }
        else if (mode == 2) { d[0] = -rn() * 5e-3; d[1] = (rn() * 2 - 1) * 1.5e-3; d[2] = (rn() * 2 - 1) * 1.5e-3; d[3] = (rn() * 2 - 1) * 5e-4; d[4] = d[5] = 0; }
        else { const double v = rn() * 1e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 3e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 2e-4; }
        for (int i = 0; i < 6; ++i) d[i] *= scT;
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
        ++tested; if (sub == 2) ++n2; else { ++n4; if (nc4 < 40) { Case& c = c4[nc4++]; for (int i = 0; i < 6; ++i) { c.s0[i] = s0[i]; c.d[i] = d[i]; } c.k0 = k0; c.fl = q.flowPotential; } }
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
        // E2 — the kernel's Dtan6 is neither of two wrong builds, rebuilt here from the kernel's own pieces: the chain
        // without b (x) K, and the last piece's tangent alone (pre-F2). A kernel carrying either mutation sits at
        // distance 0 from it and fails; the FD distances show the 1e-7 E1 tolerance would catch it as well.
        {
            // the chain WITHOUT the Sk (x) g coupling (G_p = A_p (G_(p-1) + f C), A = Dt C^-1) and the last piece's Dt
            double Cel[6][6], Cinv[6][6]; elasticC(q, Cel); invert6(Cel, Cinv);
            double J[6][6] = {{0}}, Jn[6][6] = {{0}}, s[6], k = k0, sp[6], kp2, Dl[6][6];
            for (int i = 0; i < 6; ++i) s[i] = s0[i];
            for (int p = 0; p < nf; ++p) {
                double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] * fr[p];
                PieceSens ps; returnMapTensor1(q, s, dd, k, true, sp, kp2, Dl, true, &ps);
                double A[6][6], T[6][6], Jm[6][6];
                for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { double v = 0.0; for (int m = 0; m < 6; ++m) v += Dl[i][m] * Cinv[m][j]; A[i][j] = v; }
                for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) T[i][j] = J[i][j] + fr[p] * Cel[i][j];
                for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { double v = 0.0; for (int m = 0; m < 6; ++m) v += A[i][m] * T[m][j]; Jm[i][j] = v; }   // no Sk (x) g
                for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) J[i][j] = Jm[i][j];
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
    // n2/n4 >= 20 is the coverage floor; skipped <= tested catches a generator that mostly misses (the 20x cdpm2 scale
    // did: 82 skipped vs 80 tested under UCRT, see above; at 8x 3 vs 80 under UCRT and libmmt).
    check(n2 >= 20 && n4 >= 20 && skipped <= tested && worstD < 1.0e-7 && worstK < 1.0e-7,
          "E1 failure-driven 2-/4-piece returns: chain-rule Dtan6 and d kp/de == frozen-schedule FD to 1e-7");
    std::printf("       (2-piece %d, 4-piece %d, skipped (kink in the stencil) %d; worst Dtan6 %.2e, worst dkp %.2e)\n", n2, n4, skipped, worstD, worstK);
    check(nF >= 40 && worstF < 1.0e-7, "E1 forced {1/2,1/2} and {1/4 x4} schedules: chain-rule tangent == FD to 1e-7");
    std::printf("       (forced cases %d, worst %.2e)\n", nF, worstF);
    check(nMut >= 20 && kerBK > 1.0e-5 && kerLast > 1.0e-3 && mutBK > 1.0e-5 && mutLast > 1.0e-3,
          "E2 the kernel's Dtan6 differs from a chain without b(x)K (> 1e-5) and from the last-piece tangent (> 1e-3)");
    std::printf("       (cases %d; kernel vs no-b(x)K chain >= %.2e, vs last piece >= %.2e; those mutants miss the FD by >= %.2e / %.2e)\n",
                nMut, kerBK, kerLast, mutBK, mutLast);

    // E1d — the DETERMINISTIC mode (the shipped default): Dtan6 (C8 checks it too) AND the kappa row d kp/de that the
    // damaged tangent consumes == the FD of the same L-piece chain (L = the reported subInfo), to 1e-7, n = 2 .. 64.
    {
        Params qd = q; qd.subIncrMode = 0;
        unsigned sd = 4271u; auto rd = [&]() { sd = sd * 1664525u + 1013904223u; return (sd >> 8) * (1.0 / 16777216.0); };
        int nd = 0, skd = 0, nBig = 0; double wD = 0, wK = 0;
        for (long t = 0; t < 200000 && nd < 80; ++t) {
            const int mode = t % 4; double d[6];
            if (mode == 0) { for (int i = 0; i < 6; ++i) d[i] = (rd() * 2 - 1) * 3e-3; for (int i = 3; i < 6; ++i) d[i] *= 0.5; }
            else if (mode == 1) { const double v = -rd() * 8e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rd() * 2 - 1) * 2e-4; for (int i = 3; i < 6; ++i) d[i] = (rd() * 2 - 1) * 1e-4; }
            else if (mode == 2) { d[0] = -rd() * 5e-3; d[1] = (rd() * 2 - 1) * 1.5e-3; d[2] = (rd() * 2 - 1) * 1.5e-3; d[3] = (rd() * 2 - 1) * 5e-4; d[4] = d[5] = 0; }
            else { const double v = rd() * 1e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rd() * 2 - 1) * 3e-4; for (int i = 3; i < 6; ++i) d[i] = (rd() * 2 - 1) * 2e-4; }
            const double k0 = rd() * 1.5;
            double s0[6] = {0, 0, 0, 0, 0, 0};
            if (t % 5 == 1) { const double p = -(rd()) * 10.0; s0[0] = s0[1] = s0[2] = p; s0[3] = (rd() * 2 - 1) * 3.0; }
            double sn[6], kn, Dt[6][6], dk[6]; int sub = 0;
            if (returnMapTensor(qd, s0, d, k0, true, sn, kn, Dt, true, &sub, dk) != 0 || sub < 2) continue;
            double fr[64]; const int nf = modeSchedule(qd, s0, d, k0, sub, fr);
            // a rescue (adaptive) state: its pieces are not equal; skip it here (E1 covers the adaptive path)
            double sR[6], kR, JR[6][6];
            if (nf < 2 || returnMapTensorSchedule(qd, s0, d, k0, fr, nf, sR, kR, JR, false) != 0) { ++skd; continue; }
            bool same = (kR == kn); for (int i = 0; i < 6; ++i) same = same && (sR[i] == sn[i]);
            if (!same) { ++skd; continue; }
            double eD, eK; bool smooth;
            fdSchedule(qd, s0, d, k0, fr, nf, Dt, dk, eD, eK, smooth);
            if (!smooth) { ++skd; continue; }
            ++nd; if (sub >= 8) ++nBig; wD = std::fmax(wD, eD); wK = std::fmax(wK, eK);
        }
        // bound 1e-6, C8's for the same chains (the FD floor grows with L: measured 4.6e-8 / 1.1e-7 worst at L = 64)
        check(nd >= 40 && nBig >= 5 && skd <= nd && wD < 1.0e-6 && wK < 1.0e-6,
              "E1d deterministic L-piece chains: Dtan6 and d kp/de == the FD of the same L-piece chain to 1e-6");
        std::printf("       (cases %d (L >= 8: %d), skipped (rescue / kink in the stencil) %d; worst Dtan6 %.2e, worst dkp %.2e)\n", nd, nBig, skd, wD, wK);
    }

    // E3 — N = k +/- 1e-12 (a rounded piece count): pieces of 1/N plus the remainder piece; stress and tangent are
    // continuous in N at N = 4 (a tiny remainder piece adds nothing).
    {
        double worstS = 0, worstT = 0; int n = 0;
        for (int c = 0; c < nc4; ++c) {
            const Case& cs = c4[c];
            Params qc = q; qc.flowPotential = cs.fl;
            const double f4[4] = {0.25, 0.25, 0.25, 0.25};
            double sA[6], kA, DA[6][6], dkA[6];
            if (returnMapTensorSchedule(qc, cs.s0, cs.d, cs.k0, f4, 4, sA, kA, DA, true, dkA) != 0) continue;
            bool all = true; double es = 0, et = 0;
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double N = 4.0 + sgn * 1.0e-12; const int nfull = (int)std::floor(N);
                double fr[6]; int m = 0; for (int i = 0; i < nfull; ++i) fr[m++] = 1.0 / N;
                const double rem = 1.0 - nfull / N; if (rem > 0.0) fr[m++] = rem;
                double sB[6], kB, DB[6][6], dkB[6];
                if (returnMapTensorSchedule(qc, cs.s0, cs.d, cs.k0, fr, m, sB, kB, DB, true, dkB) != 0) { all = false; break; }
                for (int i = 0; i < 6; ++i) es = std::fmax(es, std::fabs(sB[i] - sA[i]) / q.fc);
                et = std::fmax(et, relF6(DB, DA));
            }
            if (!all) continue;
            ++n; worstS = std::fmax(worstS, es); worstT = std::fmax(worstT, et);
        }
        check(n >= 10 && worstS < 1.0e-9 && worstT < 1.0e-7, "E3 N = 4 +/- 1e-12 pieces: stress and chain-rule tangent continuous in N");
        std::printf("       (cases %d, worst |dsigma|/fc %.2e, worst tangent %.2e)\n", n, worstS, worstT);
    }

    // E4 — the call counter: damagedTangent runs NO return map (every mode whose old tangent micro-FD'd the map:
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
        check(nDT > 0 && inDT == 0 && perMap == 1, "E4 call counter: damagedTangent runs 0 return maps; returnMap runs exactly 1");
        std::printf("       (damagedTangent calls %d: return maps inside %llu; returnMap calls %d: max maps per call %llu)\n", nDT, inDT, nmaps, perMap);
    }

    // E5 — the closed-form proj-weight gradients == a central FD of the closed-form weights themselves.
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
              "E5 closed-form tcTemper/ctTemper proj-weight gradients == FD of the weights to 1e-6");
        std::printf("       (tc cases %d worst %.2e; ct cases %d worst %.2e)\n", nC, worstC, nT, worstT);
    }

}

// Units (apeGmsh #2621 F2 review r2): D6 and D7 run in three unit systems. s scales every stress-dimensioned parameter
// (E, fc, ft; Gf and Gc as s*L), L every length (lch; Gf and Gc): MPa-mm (s = 1, L = 1: the original cases), kPa-m
// (s = 1e3, L = 1e-3) and Pa-m (s = 1e6, L = 1e-3, the fork's SI default). Strain, omega and the relative tangent
// error are dimensionless, so a unit-free tangent gives the same answer in all three up to round-off. (The STRESS
// path is not exactly unit-free: its pre-existing absolute guards move a converged kp in the 5th digit between
// units on D6a, and an isolated direct-return failure (below) lands on a different step in each unit system.)
static Params scaledBase(double s) {
    Params q; q.E = 30000 * s; q.nu = 0.2; q.fc = 30 * s; q.ft = 3 * s; q.e = eccentricityFromKupfer(30 * s, 3 * s, 1.16);
    q.m0 = m0Of(30 * s, 3 * s, q.e);
    q.qh0 = 0.3; q.Hp = 0.01; q.Df = 0.85; q.Ah = 0.08; q.Bh = 0.003; q.Ch = 2.0; q.Dh = 1e-6;
    return q;
}
static std::string label(const char* tag, const char* name) { return tag[0] ? std::string(name) + " [" + tag + "]" : std::string(name); }

// An ISOLATED direct-return failure (D6, D8): the centre's direct return fails and it is sub-stepped, while its FD
// stencil returns on another schedule, so the stencil's FD reads a different branch of a discontinuous map. The test
// that applies there: the effective chain-rule tangent and d kp/de at the centre (what damagedTangent consumes) ==
// the FD of the map with the centre's own failure-driven schedule frozen (failureSchedule / fdSchedule, E1's 1e-7).
static bool isolatedScheduleOk(const Params& p, const double e[6], const State& in, double& eD, double& eK) {
    double deps[6]; for (int i = 0; i < 6; ++i) deps[i] = e[i] - in.eps[i];
    double sb[6], kpn, Ce[6][6], dk[6]; int sub = 0;
    returnMapTensor(p, in.sigEff, deps, in.kp, true, sb, kpn, Ce, true, &sub, dk);
    double fr[64]; const int n = modeSchedule(p, in.sigEff, deps, in.kp, sub, fr);
    eD = eK = 1e300;
    if (n < 2 || sub != n) return false;
    bool smooth;
    fdSchedule(p, in.sigEff, deps, in.kp, fr, n, Ce, dk, eD, eK, smooth);
    return smooth && eD < 1.0e-7 && eK < 1.0e-7;
}

static void run_f2_cap(double s, double L, const char* tag) {
    std::printf("[D6] F2 (apeGmsh #2621): damaged tangent at the omega cap%s%s\n", tag[0] ? ", units " : "", tag);
    Params q = scaledBase(s); q.flowPotential = 1; q.maxSubIncr = 10;
    // D6 — the omega cap: an omega AT OMEGA_MAX is flat (d omega = 0) in damagedTangent. Before F2 a bilinear omega_t
    // (and a CDPM2-drive omega_c) at the cap fell through to the legacy exponential-law IFT, an O(1) wrong tangent on
    // a fully open crack. Drive the material INTO the cap and require the full damaged tangent (returnMap's Dtan6) ==
    // the central FD of the full damaged map (returnMap's nominal stress) to 1e-6 on EVERY converged step whose omega
    // sits at the cap. The one thing excluded is an h whose stencil leaves the cap: an h counts only when all 12
    // perturbed maps converge with that omega still == OMEGA_MAX (at the cap onset the map is one-sided). A cap step
    // with no such h is an onset step; it is printed and FAILS the check, so every cap step is in fact tested.
    // (A kernel with both F2 cap guards reverted misses here by O(1): the omega_t paths by ~2e1, omega_c by ~0.9.)
    // An h also needs every stencil map on the CENTRE's piece schedule (State::subInfo), or its FD differentiates a
    // different map. The one case this leaves is an ISOLATED direct-return failure (a pre-existing defect of the
    // direct return, filed separately): the centre's direct return fails and it is sub-stepped, while every stencil
    // neighbour at +/-1e-8..1e-10 returns directly. The map is discontinuous there and the neighbours' FD reads the
    // direct branch (1.7e-3 off on D6a in Pa). Such a step is counted, printed, and tested against the right
    // reference: the effective tangent and d kp/de that the damaged tangent consumes == the FD of the map with the
    // centre's schedule frozen (E1's 1e-7). At most 1 isolated step per unit system (0 in MPa and kPa, 1 in Pa).
    {
        // the central FD of the full damaged map at strain e from the committed state `in`; best error over the
        // h whose stencil stays at the cap of omega_t (which == 0) or omega_c (which == 1) on the centre's schedule.
        // 0 = tested (err); 1 = onset (some stencil map leaves the cap or fails); 2 = isolated schedule (every stencil
        // map converged at the cap, but on a piece schedule other than the centre's for every h).
        auto fdAtCap = [](const Params& p, const double e[6], const State& in, const double D[6][6], int which, int sub0,
                          double& err) {
            err = 1e300; bool any = false, onlySched = true;
            for (double h : {1e-8, 1e-9, 1e-10}) {
                double F[6][6]; bool stay = true, sched = true;
                for (int j = 0; j < 6 && stay; ++j) {
                    double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                    State op, om; double sp[6], sm[6], se[6], J[6][6], wtp, wcp, wtm, wcm;
                    const int a = returnMap(p, ep, in, op, sp, se, J, false, -1.0, true, &wtp, &wcp);
                    const int b = returnMap(p, em, in, om, sm, se, J, false, -1.0, true, &wtm, &wcm);
                    const double wp = which ? wcp : wtp, wm = which ? wcm : wtm;
                    if (a != 0 || b != 0 || wp != OMEGA_MAX || wm != OMEGA_MAX) stay = false;
                    if (op.subInfo != sub0 || om.subInfo != sub0) sched = false;
                    for (int i = 0; i < 6; ++i) F[i][j] = (sp[i] - sm[i]) / (2 * h);
                }
                if (!stay) { onlySched = false; continue; }
                if (!sched) continue;
                any = true; err = std::fmin(err, relF6(D, F));
            }
            return any ? 0 : (onlySched ? 2 : 1);
        };
        int nIso = 0; bool isoOk = true;
        auto isoStep = [&](const Params& p, const double e[6], const State& in, const State& o, const char* which, int k) {
            double eD, eK; const bool ok = isolatedScheduleOk(p, e, in, eD, eK);
            ++nIso; isoOk = isoOk && ok;
            std::printf("       (%s step %d: isolated direct-return failure, centre on %d pieces, the stencil on another schedule; "
                        "frozen-schedule FD: Ceff %.2e, dkp %.2e%s)\n", which, k, o.subInfo, eD, eK, ok ? "" : " FAIL");
        };
        Params p = q; p.tensionLaw = 1; p.compDrive = 1; p.tcTemper = 2; p.Gf = 0.1 * s * L; p.Gc = 5.0 * s * L; p.lch = 50.0 * L;
        p.epsFc = 1e-3; p.As = 2.0;
        // D6a — bilinear omega_t: two biaxial/shear tension paths far past wf (200 steps each, 179 + 176 at the cap;
        // path 0 enters it at step 22, where the guard-reverted kernel misses the FD by 2.6)
        int nT = 0, onT = 0; double worstT = 0;
        for (int path = 0; path < 2; ++path) {
            State in;
            for (int k = 1; k <= 200; ++k) {
                const double t = k / 200.0; double e[6] = {0, 0, 0, 0, 0, 0};
                if (path == 0) { e[0] = 2e-2 * t; e[1] = 1.5e-2 * t; e[2] = 0.5e-2 * t; }
                else           { e[0] = 2e-2 * t; e[1] = 1e-2 * t; e[3] = 2e-3 * t; }
                State o; double sg[6], se[6], D[6][6], wt, wc;
                const int st = returnMap(p, e, in, o, sg, se, D, true, -1.0, true, &wt, &wc);
                if (st == 0 && wt == OMEGA_MAX) {
                    double err; const int r = fdAtCap(p, e, in, D, 0, o.subInfo, err);
                    if (r == 0) { ++nT; worstT = std::fmax(worstT, err); } else if (r == 2) isoStep(p, e, in, o, "D6a", k); else ++onT;
                }
                in = o;
            }
        }
        check(nT >= 300 && onT == 0 && worstT < 1.0e-6,
              label(tag, "D6a bilinear omega_t at OMEGA_MAX: damaged tangent == FD of the full damaged map to 1e-6 (every cap step)").c_str());
        std::printf("       (cap steps checked %d, onset steps (stencil leaves the cap) %d; worst %.2e)\n", nT, onT, worstT);
        // D6b — CDPM2-drive omega_c: compression with lateral tension and shear (distinct principals). A brittle
        // compressive law (epsFc = 1e-4) enters the cap at step 78 (3.9% axial strain; 123 cap steps). With the
        // default epsFc = 1e-3 the cap needs ~30% strain. omega_t sits at its cap on these steps as well, so D6b
        // fails on either revert: the omega_c guard alone by 0.88, the omega_t guard alone by 3e-2 (D6a reads
        // the same with the omega_c guard reverted, so it isolates the omega_t guard).
        Params pc = p; pc.epsFc = 1e-4;
        // #877's dead-point treatment freezes a point once omega_c >= omegaDead (0.998 < OMEGA_MAX): a frozen point never
        // reaches damagedTangent, so the live omega_c-at-the-cap branch this case tests needs omegaDead = 1 (off), the
        // regime F2 measured on d14fb6958 (which had no dead points)
        pc.omegaDead = 1.0;
        int nC = 0, onC = 0; double worstC = 0;
        {
            State in;
            for (int k = 1; k <= 200; ++k) {
                const double t = k / 200.0, M = 1e-1;
                double e[6] = {-M * t, M * t, 0.8 * M * t, 0.2 * M * t, 0, 0};
                State o; double sg[6], se[6], D[6][6], wt, wc;
                const int st = returnMap(pc, e, in, o, sg, se, D, true, -1.0, true, &wt, &wc);
                if (st == 0 && wc == OMEGA_MAX) {
                    double err; const int r = fdAtCap(pc, e, in, D, 1, o.subInfo, err);
                    if (r == 0) { ++nC; worstC = std::fmax(worstC, err); } else if (r == 2) isoStep(pc, e, in, o, "D6b", k); else ++onC;
                }
                in = o;
            }
        }
        check(nC >= 100 && onC == 0 && worstC < 1.0e-6,
              label(tag, "D6b CDPM2-drive omega_c at OMEGA_MAX: damaged tangent == FD of the full damaged map to 1e-6 (every cap step)").c_str());
        std::printf("       (cap steps checked %d, onset steps (stencil leaves the cap) %d; worst %.2e)\n", nC, onC, worstC);
        check(nIso <= 1 && isoOk,
              label(tag, "D6 isolated direct-return failures (<= 1): chain-rule tangent == frozen-schedule FD to 1e-7").c_str());
        std::printf("       (isolated steps %d)\n", nIso);
    }
}

// D7 — the tension-apex return (apeGmsh #2621 F2 onset audit): on the hydrostatic axis the damage drivers (eps_tilde,
// x_s, alpha_c, eps_tilde_c) are cones, and the per-component central micro-FD of their gradient averaged the tensile-
// and compressive-meridian slopes into a spurious hydrostatic gradient (O(1): 1.00 on 6b7ebf3a8 in MPa, kPa and Pa). The apex return keeps
// sig_bar on the axis, so the true gradient is the one ALONG the axis (dscalarDsigAxis). Hydrostatic tension from the
// virgin state, legacy + cdpm2 flow x ctTemper none/proj x tcTemper none/proj, every damage-loading apex step: the
// damaged tangent == the central FD (best of h 1e-8/1e-9/1e-10, every stencil point converged and still on the apex).
static void run_f2_apex(double s, double L, const char* tag) {
    std::printf("[D7] F2 (apeGmsh #2621): damaged tangent at the tension-apex return%s%s\n", tag[0] ? ", units " : "", tag);
    Params q = scaledBase(s);
    q.tensionLaw = 1; q.compDrive = 1; q.Gf = 0.1 * s * L; q.Gc = 5.0 * s * L; q.lch = 50.0 * L; q.epsFc = 1e-3; q.As = 2.0;
    auto rhoOf = [](const double s[6]) { double xi, rho, th; invariants(s, xi, rho, th); return rho; };
    int n = 0, skipped = 0; double worst = 0;
    for (int fl = 0; fl <= 1; ++fl) for (int ct = 0; ct <= 2; ct += 2) for (int tc = 0; tc <= 2; tc += 2) {
        Params p = q; p.flowPotential = fl; p.maxSubIncr = fl ? 10 : 0; p.ctTemper = ct; p.tcTemper = tc;
        State in;
        for (int k = 1; k <= 30; ++k) {
            const double a = 6e-4 * k / 30.0;
            double e[6] = {a, a, a, 0, 0, 0};
            State o; double sg[6], se[6], D[6][6];
            const int st = returnMap(p, e, in, o, sg, se, D, true, -1.0, true);
            if (st == 0 && rhoOf(se) <= 1e-12 * s && se[0] > 0 && o.wt > 0 && o.wt < OMEGA_MAX && o.et_max > in.et_max) {
                double best = 1e300; bool anyOk = false;
                for (double h : {1e-8, 1e-9, 1e-10}) {
                    double F[6][6]; bool ok = true;
                    for (int j = 0; j < 6 && ok; ++j) {
                        double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                        State t; double a1[6], b1[6], c1[6], d1[6], J[6][6];
                        ok = returnMap(p, ep, in, t, a1, b1, J, false, -1.0, true) == 0 && returnMap(p, em, in, t, c1, d1, J, false, -1.0, true) == 0
                             && rhoOf(b1) <= 1e-10 * s && rhoOf(d1) <= 1e-10 * s;
                        for (int i = 0; i < 6; ++i) F[i][j] = (a1[i] - c1[i]) / (2 * h);
                    }
                    if (!ok) continue;
                    anyOk = true; best = std::fmin(best, relF6(D, F));
                }
                if (anyOk) { ++n; worst = std::fmax(worst, best); } else ++skipped;
            }
            in = o;
        }
    }
    // skipped = an apex state for which no h keeps every stencil point converged and on the apex: an untested state,
    // so none is allowed (216 tested in each unit system).
    check(n >= 50 && skipped == 0 && worst < 1.0e-6,
          label(tag, "D7 tension-apex return: damaged tangent == FD of the damaged map to 1e-6, no apex state skipped").c_str());
    std::printf("       (apex damage-loading states %d, skipped %d, worst %.2e)\n", n, skipped, worst);
}

// D8 — d kp/de under Duvaut-Lions (apeGmsh #2621 F2 review r2). The relaxed kp is (1 - beta) kp_n + beta kp, so the
// d kp/de that beta_c's qh2(kp) consumes is beta times the inviscid one (returnMap: dkp *= beta). It matters only where
// all three hold: kp >= 1 (qh2 is flat below), beta_c active (0 < beta_c < 1 under compressive loading) and
// -eta > 0. Triaxial compression with lateral expansion (cdpm2 flow, Hp = 0.5 so qh2 moves, eta = dt => beta = 1/2):
// every step that starts at kp >= 1 (qh2 on its linear branch over the whole step) with 0 < beta_c < 1 and kp
// advancing, damaged tangent == the FD of the damaged map at the same dt, every stencil map on the centre's piece
// schedule, none skipped. (The first-yield step, kp 0 -> 1.6, starts below 1; it is also an isolated direct-return
// failure, the pre-existing defect filed separately.) Deleting the dkp *= beta line misses here by 2.2e-3. Run in the
// three unit systems since #2639 made the driver-gradient step relative (with the absolute step the MPa residual was
// 1.8e-7 and Pa was not run).
static void run_f2_viscous_kp(double s, double L, const char* tag) {
    std::printf("[D8] F2 (apeGmsh #2621): d kp/de under Duvaut-Lions (kp >= 1, beta_c active, eta > 0)%s%s\n",
                tag[0] ? ", units " : "", tag);
    Params p = scaledBase(s); p.Hp = 0.5;
    p.tensionLaw = 1; p.compDrive = 1; p.Gf = 0.1 * s * L; p.Gc = 5.0 * s * L; p.lch = 50.0 * L; p.epsFc = 1e-3; p.As = 2.0;
    p.flowPotential = 1; p.maxSubIncr = 10; p.eta = 1.0;
    const double dt = 1.0;
    State in; int n = 0, skipped = 0; double worst = 0;
    for (int k = 1; k <= 120; ++k) {
        const double t = k / 120.0, M = 3e-2;
        double e[6] = {-M * t, 0.5 * M * t, 0.5 * M * t, 0.1 * M * t, 0, 0};
        State o; double sg[6], se[6], D[6][6];
        const int st = returnMap(p, e, in, o, sg, se, D, true, dt, true);
        double A[3][3], w[3], V[3][3]; voigtToMat(se, A); eig3sym(A, w, V);
        const double bc = betaC(w, o.kp, p);
        if (st == 0 && in.kp >= 1.0 && bc > 0.0 && bc < 1.0 && o.kp > in.kp) {
            double best = 1e300;
            for (double h : {1e-8, 1e-9, 1e-10}) {
                double F[6][6]; bool ok = true, sched = true;
                for (int j = 0; j < 6 && ok; ++j) {
                    double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                    State a, c; double a1[6], b1[6], c1[6], d1[6], J[6][6];
                    ok = returnMap(p, ep, in, a, a1, b1, J, false, dt, true) == 0 && returnMap(p, em, in, c, c1, d1, J, false, dt, true) == 0;
                    if (a.subInfo != o.subInfo || c.subInfo != o.subInfo) sched = false;
                    for (int i = 0; i < 6; ++i) F[i][j] = (a1[i] - c1[i]) / (2 * h);
                }
                if (ok && sched) best = std::fmin(best, relF6(D, F));
            }
            if (best < 1e300) { ++n; worst = std::fmax(worst, best); } else ++skipped;
        }
        in = o;
    }
    check(n >= 100 && skipped == 0 && worst < 1.0e-6,
          label(tag, "D8 kp >= 1, beta_c active, eta > 0: damaged tangent == FD of the damaged map to 1e-6 (d kp/de scaled by beta)").c_str());
    std::printf("       (states %d, skipped %d, worst %.2e)\n", n, skipped, worst);
}


// D8b — the proj damage weights under Duvaut-Lions (apeGmsh #2621, item 4 of the F2 port). Before the port the
// ctTemper-proj and tcTemper-proj weight gradients were composite micro-FDs through the INVISCID return map (rmForFD),
// while the stress they differentiate is the RELAXED one (eta > 0): the same defect as beta_c's (D8). They are now closed
// forms chained through the relaxed Ceff and I - C^-1 Ceff, so they follow the relaxation by construction. Two cdpm2-flow
// paths with ctTemper proj and tcTemper proj, eta = dt (beta = 1/2): (0) tension with shear, (1) compression with lateral
// expansion and shear. Every step that advances an omega already > 0: the damaged tangent == the FD of the damaged map at
// the same dt (every stencil map on the centre's piece schedule) to 1e-6, none skipped, and at least 20 of the states
// have a non-zero proj-weight gradient (the term under test), in MPa, kPa and Pa.
static void run_proj_viscous(double s, double L, const char* tag) {
    std::printf("[D8b] apeGmsh #2621: the proj damage weights under Duvaut-Lions (eta > 0)%s%s\n", tag[0] ? ", units " : "", tag);
    Params q = scaledBase(s);
    q.tensionLaw = 1; q.compDrive = 1; q.Gf = 0.1 * s * L; q.Gc = 5.0 * s * L; q.lch = 50.0 * L; q.epsFc = 1e-3; q.As = 2.0;
    q.flowPotential = 1; q.maxSubIncr = 10; q.eta = 1.0; q.ctTemper = 2; q.tcTemper = 2;
    const double dt = 1.0;
    int n = 0, skipped = 0, nAct = 0; double worst = 0;
    for (int path = 0; path < 2; ++path) {
        State in;
        for (int k = 1; k <= 120; ++k) {
            const double t = k / 120.0; double e[6] = {0, 0, 0, 0, 0, 0};
            if (path == 0) { e[0] = 1.5e-3 * t; e[1] = -3e-4 * t; e[3] = 6e-4 * t; e[5] = 2e-4 * t; }
            else { const double M = 3e-2; e[0] = -M * t; e[1] = 0.5 * M * t; e[2] = 0.4 * M * t; e[3] = 0.1 * M * t; }
            State o; double sg[6], se[6], D[6][6];
            const int st = returnMap(q, e, in, o, sg, se, D, true, dt, true);
            if (st == 0 && in.wt < q.omegaDead && ((o.wc > in.wc && in.wc > 0.0) || (o.wt > in.wt && in.wt > 0.0))) {
                double best = 1e300;
                for (double h : {1e-8, 1e-9, 1e-10}) {
                    double F[6][6]; bool ok = true;
                    for (int j = 0; j < 6 && ok; ++j) {
                        double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                        State a, c; double a1[6], b1[6], c1[6], d1[6], J[6][6];
                        ok = returnMap(q, ep, in, a, a1, b1, J, false, dt, true) == 0 && returnMap(q, em, in, c, c1, d1, J, false, dt, true) == 0
                             && a.subInfo == o.subInfo && c.subInfo == o.subInfo;
                        for (int i = 0; i < 6; ++i) F[i][j] = (a1[i] - c1[i]) / (2 * h);
                    }
                    if (ok) best = std::fmin(best, relF6(D, F));
                }
                if (best < 1e300) { ++n; worst = std::fmax(worst, best); } else ++skipped;
                // is a proj-weight gradient active on this state?
                double ep1[6], ep0[6], dp[6]; plasticStrain6(se, e, q, ep1); plasticStrain6(in.sigEff, in.eps, q, ep0);
                for (int i = 0; i < 6; ++i) dp[i] = ep1[i] - ep0[i];
                double A[3][3], w[3], V[3][3]; voigtToMat(se, A); eig3sym(A, w, V);
                double gs[6], gd[6], hs[6], hd[6], g2 = 0;
                compressiveDamageWeightGrad(q, dp, w, V, gs, gd); tensileDamageWeightGrad(q, dp, w, V, hs, hd);
                for (int i = 0; i < 6; ++i) g2 += gs[i] * gs[i] + gd[i] * gd[i] + hs[i] * hs[i] + hd[i] * hd[i];
                if (g2 > 0.0) ++nAct;
            }
            in = o;
        }
    }
    check(n >= 40 && skipped == 0 && nAct >= 20 && worst < 1.0e-6,
          label(tag, "D8b ctTemper/tcTemper proj, eta > 0: damaged tangent == FD of the damaged map to 1e-6, none skipped").c_str());
    std::printf("       (damage-advancing states %d, skipped %d, with an active proj-weight gradient %d; worst %.2e)\n", n, skipped, nAct, worst);
}

// D9 — the driver-gradient step (apeGmsh #2639). damagedTangent differentiates the damage drivers (eps_tilde, x_s,
// alpha_c, eps_tilde_c) with respect to sig_bar by a central micro-FD. Its step was an ABSOLUTE 1e-6 in stress units:
// far below the stress scale on far-compressive states and in Pa, where round-off dominates. It is now DRIVER_FD_REL *
// S, S = max(ft, max|w_a|). Three cdpm2-flow paths, every step that ADVANCES an omega that is already > 0 (the onset
// step is a kink of omega and is not differentiable): (0) far compressive, oedometric compression with shear, |sig_bar|
// up to 2.5e3 MPa; (1) compression with lateral expansion and shear; (2) tension with shear. The damaged tangent ==
// the central FD of the damaged map (best of h 1e-8/1e-9/1e-10, every stencil map converged on the centre's piece
// schedule) to 1e-6, none skipped. On 6b7ebf3a8 (the absolute step) path 0 misses in every unit system (MPa 7.5e-5,
// kPa 6.6e-2, Pa 8.1e1) and paths 1 and 2 miss in Pa (5.3e-4, 6.3e-4); with the relative step all read <= 2.2e-7.
static void run_f2_driver_step(double s, double L, const char* tag) {
    std::printf("[D9] F2 (apeGmsh #2639): damaged tangent with the relative driver-gradient step%s%s\n",
                tag[0] ? ", units " : "", tag);
    Params q = scaledBase(s);
    q.tensionLaw = 1; q.compDrive = 1; q.Gf = 0.1 * s * L; q.Gc = 5.0 * s * L; q.lch = 50.0 * L; q.epsFc = 1e-3; q.As = 2.0;
    q.flowPotential = 1; q.maxSubIncr = 10;
    const char* names[3] = {"far compressive", "compression + lateral expansion", "tension + shear"};
    const int minN[3] = {40, 20, 100};
    for (int path = 0; path < 3; ++path) {
        State in; int n = 0, skipped = 0; double worst = 0, sMax = 0;
        for (int k = 1; k <= 120; ++k) {
            const double t = k / 120.0; double e[6] = {0, 0, 0, 0, 0, 0};
            if (path == 0)      { const double M = 1.5e-1; e[0] = -M * t; e[3] = 0.04 * M * t; e[4] = 0.02 * M * t; }
            else if (path == 1) { const double M = 3e-2; e[0] = -M * t; e[1] = 0.5 * M * t; e[2] = 0.4 * M * t; e[3] = 0.1 * M * t; }
            else                { e[0] = 1.5e-3 * t; e[1] = -3e-4 * t; e[3] = 6e-4 * t; e[5] = 2e-4 * t; }
            State o; double sg[6], se[6], D[6][6];
            const int st = returnMap(q, e, in, o, sg, se, D, true, -1.0, true);
            if (st == 0 && ((o.wc > in.wc && in.wc > 0.0) || (o.wt > in.wt && in.wt > 0.0))) {
                double best = 1e300;
                for (double h : {1e-8, 1e-9, 1e-10}) {
                    double F[6][6]; bool ok = true;
                    for (int j = 0; j < 6 && ok; ++j) {
                        double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                        State a, c; double a1[6], b1[6], c1[6], d1[6], J[6][6];
                        ok = returnMap(q, ep, in, a, a1, b1, J, false, -1.0, true) == 0 && returnMap(q, em, in, c, c1, d1, J, false, -1.0, true) == 0
                             && a.subInfo == o.subInfo && c.subInfo == o.subInfo;
                        for (int i = 0; i < 6; ++i) F[i][j] = (a1[i] - c1[i]) / (2 * h);
                    }
                    if (ok) best = std::fmin(best, relF6(D, F));
                }
                if (best < 1e300) { ++n; worst = std::fmax(worst, best); } else ++skipped;
                for (int i = 0; i < 3; ++i) sMax = std::fmax(sMax, std::fabs(se[i]) / s);
            }
            in = o;
        }
        check(n >= minN[path] && skipped == 0 && worst < 1.0e-6,
              label(tag, (std::string("D9 ") + names[path] + ": damaged tangent == FD of the damaged map to 1e-6, none skipped").c_str()).c_str());
        std::printf("       (damage-advancing states %d, skipped %d, worst %.2e, max |sig_bar| %.3g MPa)\n", n, skipped, worst, sMax);
    }
}

// D10 — the deterministic mode's ADAPTIVE RESCUE (ladder exhausted: n, 2n and 4n equal pieces all fail; apeGmsh #2621,
// ladruno-only gap 1). The rescue integrates the failure-driven halving schedule (returnMapTensorAdaptive), and its
// tangent must be the chain rule over THAT schedule. Before the port it reported the last piece's one-piece tangent
// (the gap review's rescuechain probe on the corbel record dumps: p50 0.25-0.36, p90 0.82-0.84 against the chain).
// subIncrMaxPieces = 1 keeps the ladder at {2, 4} pieces so the rescue is reachable with small trials; every case
// is checked to be a rescue (direct, {1/2 x 2} and {1/4 x 4} all fail) and its schedule is re-derived independently
// (failureSchedule) and replayed bit-exactly.
//  (a) MEASURABLE cases: Dtan6 and d kp/de == the central FD of the map with that schedule frozen (fdSchedule) to 1e-6
//      (E1d's and C8's bound for multi-piece chains: the FD floor grows with the piece count; measured worst 2.4e-7 at
//      7 pieces in kPa, <= 2.3e-8 elsewhere).
//  (b) UNMEASURABLE cases: on the cdpm2 flow most rescues sit where a piece of the frozen schedule converges at the
//      centre but not at +/-h for any h in {1e-8, 1e-9, 1e-10} (the direct-return convergence boundary; the map has no
//      two-sided neighbourhood there). Those are verified against the chain rule rebuilt independently over the same
//      schedule from the pieces' own sensitivities (returnMapTensor1 + PieceSens; the recurrence written out here),
//      to 1e-12 relative. The piece sensitivities themselves are E1's and C8's subject.
static void run_rescue_chain(double s, const char* tag) {
    std::printf("[D10] apeGmsh #2621: the deterministic mode's adaptive rescue reports the chain over its own schedule%s%s\n",
                tag[0] ? ", units " : "", tag);
    Params qb = scaledBase(s); qb.maxSubIncr = 10; qb.subIncrMode = 0; qb.subIncrMaxPieces = 1; qb.subIncrRescue = true;
    unsigned sr = 2621001u; auto rn = [&]() { sr = sr * 1664525u + 1013904223u; return (sr >> 8) * (1.0 / 16777216.0); };
    int n = 0, nU = 0, skipped = 0, notRescue = 0, nLeg = 0, nCd = 0; double worstD = 0, worstK = 0, worstU = 0, lastMin = 1e300;
    for (long t = 0; t < 400000 && (n < 40 || nU < 20); ++t) {
        Params q = qb; q.flowPotential = ((t >> 2) & 1) ? 1 : 0;
        const double scT = q.flowPotential ? 25.0 : 1.5;
        const int mode = t % 4; double d[6];
        if (mode == 0) { for (int i = 0; i < 6; ++i) d[i] = (rn() * 2 - 1) * 3e-3; for (int i = 3; i < 6; ++i) d[i] *= 0.5; }
        else if (mode == 1) { const double v = -rn() * 8e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 2e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 1e-4; }
        else if (mode == 2) { d[0] = -rn() * 5e-3; d[1] = (rn() * 2 - 1) * 1.5e-3; d[2] = (rn() * 2 - 1) * 1.5e-3; d[3] = (rn() * 2 - 1) * 5e-4; d[4] = d[5] = 0; }
        else { const double v = rn() * 1e-3; for (int i = 0; i < 3; ++i) d[i] = v + (rn() * 2 - 1) * 3e-4; for (int i = 3; i < 6; ++i) d[i] = (rn() * 2 - 1) * 2e-4; }
        for (int i = 0; i < 6; ++i) d[i] *= scT;
        const double k0 = rn() * 1.5;
        double s0[6] = {0, 0, 0, 0, 0, 0};
        if (t % 5 == 1) { const double p = -(rn()) * 10.0 * s; s0[0] = s0[1] = s0[2] = p; s0[3] = (rn() * 2 - 1) * 3.0 * s; }
        double sn[6], kn, Dt[6][6], dk[6]; int sub = 0;
        if (returnMapTensor(q, s0, d, k0, true, sn, kn, Dt, true, &sub, dk) != 0 || sub < 2) continue;
        // a rescue: the direct return and both ladder levels fail
        { double a[6], b, c[6][6]; if (returnMapTensor1(q, s0, d, k0, true, a, b, c, false) == 0) continue; }
        const double h2[2] = {0.5, 0.5}, h4[4] = {0.25, 0.25, 0.25, 0.25};
        { double a[6], b, c[6][6];
          if (returnMapTensorSchedule(q, s0, d, k0, h2, 2, a, b, c, false) == 0 || returnMapTensorSchedule(q, s0, d, k0, h4, 4, a, b, c, false) == 0) { ++notRescue; continue; } }
        double fr[64]; const int nf = failureSchedule(q, s0, d, k0, fr);
        if (nf != sub) { ++skipped; continue; }
        // the chain rebuilt independently over the schedule (and the last piece's tangent, the pre-port report)
        double G[6][6] = {{0}}, g[6] = {0}, Dl[6][6], C0m[6][6], Cinv[6][6]; elasticC(q, C0m); invert6(C0m, Cinv);
        bool same = true;
        { double sp[6], k = k0, sq[6], kq; for (int i = 0; i < 6; ++i) sp[i] = s0[i];
          for (int p = 0; p < nf && same; ++p) {
              double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] * fr[p];
              PieceSens ps; if (returnMapTensor1(q, sp, dd, k, true, sq, kq, Dl, true, &ps) != 0) { same = false; break; }
              double A[6][6], RC[6], Gn[6][6], gn[6];
              for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { double v = 0; for (int m = 0; m < 6; ++m) v += Dl[i][m] * Cinv[m][j]; A[i][j] = v; }
              for (int j = 0; j < 6; ++j) { double v = 0; for (int m = 0; m < 6; ++m) v += ps.R[m] * C0m[m][j]; RC[j] = v; }
              for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { double v = fr[p] * Dl[i][j] + ps.Sk[i] * g[j]; for (int m = 0; m < 6; ++m) v += A[i][m] * G[m][j]; Gn[i][j] = v; }
              for (int j = 0; j < 6; ++j) { double v = ps.Kk * g[j] + fr[p] * RC[j]; for (int m = 0; m < 6; ++m) v += ps.R[m] * G[m][j]; gn[j] = v; }
              for (int i = 0; i < 6; ++i) { g[i] = gn[i]; sp[i] = sq[i]; for (int j = 0; j < 6; ++j) G[i][j] = Gn[i][j]; }
              k = kq;
          }
          same = same && k == kn; for (int i = 0; i < 6; ++i) same = same && sp[i] == sn[i]; }
        if (!same) { ++skipped; continue; }
        double eD, eK; bool smooth;
        fdSchedule(q, s0, d, k0, fr, nf, Dt, dk, eD, eK, smooth);
        if (!smooth || eD >= 1e300) {
            if (nU >= 20) continue;
            ++nU; worstU = std::fmax(worstU, std::fmax(relF6(Dt, G), relV6(dk, g)));
            continue;
        }
        if (n >= 40) continue;
        ++n; worstD = std::fmax(worstD, eD); worstK = std::fmax(worstK, eK);
        if (q.flowPotential) ++nCd; else ++nLeg;
        { double e1, e2; bool sm; const double dz[6] = {0, 0, 0, 0, 0, 0};
          fdSchedule(q, s0, d, k0, fr, nf, Dl, dz, e1, e2, sm); if (sm) lastMin = std::fmin(lastMin, e1); }
    }
    check(n >= 20 && worstD < 1.0e-6 && worstK < 1.0e-6,
          label(tag, "D10a adaptive rescue: Dtan6 and d kp/de == the FD of the map on its own frozen schedule to 1e-6").c_str());
    std::printf("       (measurable rescue cases %d (legacy %d, cdpm2 %d); worst Dtan6 %.2e, worst dkp %.2e; the last-piece tangent misses the"
                " same FD by >= %.2e)\n", n, nLeg, nCd, worstD, worstK, lastMin);
    check(nU >= 10 && worstU < 1.0e-12,
          label(tag, "D10b adaptive rescue with no two-sided FD: Dtan6 and d kp/de == the chain rebuilt over its schedule to 1e-12").c_str());
    std::printf("       (FD-unmeasurable rescue cases %d, worst %.2e; schedule not re-derived %d, ladder-solvable trials passed over %d)\n",
                nU, worstU, skipped, notRescue);
}

// D11 — the tension-DEAD cut with a sub-stepped compressive remainder (apeGmsh #2621, ladruno-only gap 2). A cracked
// point (committed omega_t >= omegaDead) carries the tensile part of the trial elastically and returns the compressive
// remainder Minus(sig_tr) with a ZERO increment; when that return takes n >= 2 deterministic pieces, the effective
// tangent must be the derivative of the map with respect to its START (wrtStart), not d/d(deps) C^-1 of the
// zero-increment chain (exact for one piece only; the gap review's cutchain probe read 0.21-2.06 off at n >= 2).
// Synthetic cracked points (omega_t = 0.999, kappa_p 0 / 0.5 / 2) compressed along three axes with shear, inviscid and
// under Duvaut-Lions (eta = dt), with the ctTemper / tcTemper proj weights off and on: the effective tangent and d kp/de (effectiveReturn) == the FD of the effective map with
// n pinned (returnMapPieces + subIncrForceN) to 1e-7, and the full damaged tangent (returnMap) == the FD of the full
// damaged map with n pinned to 1e-6, every stencil map converged on the centre's piece count; at least 10 states with
// n >= 2, in MPa, kPa and Pa.
static void run_cut_chain(double s, double L, const char* tag) {
    std::printf("[D11] apeGmsh #2621: the tension-dead cut with a sub-stepped compressive remainder%s%s\n", tag[0] ? ", units " : "", tag);
    Params pb = scaledBase(s); pb.flowPotential = 1; pb.maxSubIncr = 10; pb.tensionLaw = 1; pb.compDrive = 1;
    pb.Gf = 0.1 * s * L; pb.Gc = 5.0 * s * L; pb.lch = 50.0 * L; pb.epsFc = 1e-3; pb.As = 2.0;
    int nMulti = 0, nOne = 0, skipped = 0; double wE = 0, wK = 0, wD = 0;
    for (int pj = 0; pj <= 2; pj += 2) for (int visc = 0; visc <= 1; ++visc) for (double kp0 : {0.0, 0.5, 2.0}) {
        Params p = pb; p.eta = visc ? 1.0 : 0.0; const double dt = visc ? 1.0 : -1.0; p.ctTemper = pj; p.tcTemper = pj;
        State in; in.wt = 0.999; in.kp = kp0; in.et_max = 1e-2;
        for (double c : {-1e-3, -3e-3, -6e-3, -1.2e-2, -2.4e-2}) for (int dir = 0; dir < 3; ++dir) {
            double e[6]; for (int i = 0; i < 6; ++i) e[i] = in.eps[i]; e[dir] += c; e[3] += 0.3 * c; e[5] += 0.1 * c;
            State cb; const State* ie = nullptr; bool cut = false; double sE[6], kE, Ce[6][6], dk[6]; int sub = 0;
            if (effectiveReturn(p, e, in, true, dt, true, cb, ie, cut, sE, kE, Ce, &sub, dk) != 0 || !cut) { ++skipped; continue; }
            State o; double sg[6], se[6], D[6][6];
            if (returnMap(p, e, in, o, sg, se, D, true, dt, true) != 0) { ++skipped; continue; }
            const int n0 = returnMapPieces(p, in, e); Params q = p; q.subIncrForceN = n0;
            double bE = 1e300, bK = 1e300, bD = 1e300;
            for (double h : {1e-8, 1e-9, 1e-10}) {
                double FE[6][6], FK[6], FD[6][6]; bool ok = true;
                for (int j = 0; j < 6 && ok; ++j) {
                    double ep[6], em[6]; for (int i = 0; i < 6; ++i) ep[i] = em[i] = e[i]; ep[j] += h; em[j] -= h;
                    State c1, c2; const State* i1; const State* i2; bool t1, t2; double a[6], b[6], k1, k2, J[6][6]; int s1 = 0, s2 = 0;
                    ok = effectiveReturn(q, ep, in, true, dt, false, c1, i1, t1, a, k1, J, &s1) == 0
                      && effectiveReturn(q, em, in, true, dt, false, c2, i2, t2, b, k2, J, &s2) == 0 && s1 == sub && s2 == sub && t1 && t2;
                    State op, om; double gp[6], gm[6], x[6];
                    ok = ok && returnMap(q, ep, in, op, gp, x, J, false, dt, true) == 0 && returnMap(q, em, in, om, gm, x, J, false, dt, true) == 0
                         && op.subInfo == o.subInfo && om.subInfo == o.subInfo;
                    for (int i = 0; i < 6; ++i) { FE[i][j] = (a[i] - b[i]) / (2 * h); FD[i][j] = (gp[i] - gm[i]) / (2 * h); }
                    FK[j] = (k1 - k2) / (2 * h);
                }
                if (!ok) continue;
                bE = std::fmin(bE, relF6(Ce, FE)); bK = std::fmin(bK, relV6(dk, FK)); bD = std::fmin(bD, relF6(D, FD));
            }
            if (bE == 1e300) { ++skipped; continue; }
            if (n0 >= 2) ++nMulti; else ++nOne;
            wE = std::fmax(wE, bE); wK = std::fmax(wK, bK); wD = std::fmax(wD, bD);
        }
    }
    check(nMulti >= 10 && skipped == 0 && wE < 1.0e-7 && wK < 1.0e-7 && wD < 1.0e-6,
          label(tag, "D11 tension-dead cut: effective tangent and d kp/de == n-pinned FD to 1e-7, damaged tangent to 1e-6").c_str());
    std::printf("       (cut states n >= 2: %d, n = 1: %d, skipped %d; worst Ceff %.2e, dkp %.2e, damaged %.2e)\n", nMulti, nOne, skipped, wE, wK, wD);
}

// D12 â€” the kappa row of the hardening return ON THE COMPRESSIVE MERIDIAN (apeGmsh #2621 F2 onset audit). With two
// equal principal trial stresses (every uniaxial and axisymmetric compression state) cos 3theta sits at its extremum,
// where theta = acos(cos 3theta)/3 has a square-root kink. principalJacobian took d cos2/dw (R4's Lode factor) by a
// central micro-FD there, which read the round-off in cos 3theta as a slope (2.8e-4 for an exact 0 along the axial
// principal): the stress rows were right (closed-form dr/dw), the kappa row was not (1.6e-3 on a uniaxial direct
// return, the residual of dmgc_uniax_onset).
// Reference: along the two principal directions that keep the meridian (the distinct principal; the equal pair together)
// theta stays EXACTLY pi/3, so dR3/dw = dR4/dw = 0 and the trial enters R1, R2 only through xi_tr and rho_tr. The kappa
// derivative is then the IFT of the residual cdpm2HardeningResidual (its own analytic Jacobian, separate code from
// principalJacobian's) with that exact right-hand side. Every piece of schedules {1} and {1/2, 1/2}, both flows, kappa_n
// below and above 1, MPa / kPa / Pa: principalJacobian's d kappa/dw along both directions == the reference, 2-vector, to
// 1e-7. (A central FD of the map cannot certify this at 1e-7 on the meridian: the same sqrt(eps) theta round-off leaves a
// random ~2.4e-8 relative offset in the map's own cos2, so an FD of kappa there reads ~1e-6 at h = 1e-5; it is printed.)
static void run_meridian_kappa(double s, const char* tag) {
    std::printf("[D12] apeGmsh #2621: the kappa row on the compressive meridian%s%s\n", tag[0] ? ", units " : "", tag);
    int n = 0, nLo = 0, nHi = 0, n1 = 0, n2 = 0; double worst = 0.0, worstMapFD = 0.0;
    for (int flow = 0; flow < 2; ++flow) for (double k0 : {0.3, 0.95, 1.05, 1.5}) for (double a : {-1.0e-3, -2.0e-3, -4.0e-3})
    for (double lat : {0.0, 0.2, -0.1}) for (double pre : {0.0, -0.3}) for (int sch = 0; sch < 2; ++sch) {
        Params q = scaledBase(s); q.flowPotential = flow; q.maxSubIncr = 10;
        const double s0[6] = {pre * q.fc, 0.05 * pre * q.fc, 0.05 * pre * q.fc, 0, 0, 0};
        const double d[6] = {a, -lat * a, -lat * a, 0, 0, 0};
        const int nf = sch ? 2 : 1;
        double sp[6]; for (int i = 0; i < 6; ++i) sp[i] = s0[i]; double kpp = k0;
        for (int p = 0; p < nf; ++p) {
            double dd[6]; for (int i = 0; i < 6; ++i) dd[i] = d[i] / nf;
            double str[6]; elasticPredTensor(sp, dd, q, str);
            double Am[3][3], w[3], V[3][3]; voigtToMat(str, Am); eig3sym(Am, w, V);
            const PrincipalResult pr = returnMapHardening(w, q, kpp);
            if (pr.converged && pr.plastic && !pr.apex && pr.kp > kpp) {
                double Dp[3][3], hs[7]; principalJacobian(w, pr, q, true, Dp, hs);
                // the distinct principal (the axial one) and the equal pair
                int i0 = 0; { double best = -1.0; for (int i = 0; i < 3; ++i) { const double g = std::fabs(w[i] - w[(i + 1) % 3]) + std::fabs(w[i] - w[(i + 2) % 3]); if (g > best) { best = g; i0 = i; } } }
                double dv[2][3] = {{0, 0, 0}, {1, 1, 1}}; dv[0][i0] = 1.0; dv[1][i0] = 0.0;
                double sv[6] = {w[0], w[1], w[2], 0, 0, 0}, xitr, rhotr, thtr; invariants(sv, xitr, rhotr, thtr);
                double R[4], J[4][4];
                cdpm2HardeningResidual(pr.xi, pr.rho, pr.dlam, pr.kp, xitr, rhotr, kpp, lodeR(thtr, q.e), hardeningLodeFactor(rhotr, thtr, q),
                                       q, R, J, q.flowPotential != 1);
                const double sV = (w[0] + w[1] + w[2]) / 3.0, st[3] = {w[0] - sV, w[1] - sV, w[2] - sV};
                double an[2], ref[2];
                for (int r = 0; r < 2; ++r) {
                    an[r] = 0.0; double dxi = 0.0, drho = 0.0;
                    for (int b = 0; b < 3; ++b) { an[r] += hs[3 + b] * dv[r][b]; dxi += dv[r][b] / SQRT3; drho += st[b] / rhotr * dv[r][b]; }
                    double M[4][5]; const double rhs[4] = {dxi, drho, 0.0, 0.0};   // -dR/dw: R1 = xi - xi_tr .., R2 = rho - rho_tr ..
                    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = J[i][j]; M[i][4] = rhs[i]; }
                    for (int c = 0; c < 4; ++c) {
                        int piv = c; for (int qq = c + 1; qq < 4; ++qq) if (std::fabs(M[qq][c]) > std::fabs(M[piv][c])) piv = qq;
                        for (int j = 0; j < 5; ++j) { const double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
                        for (int qq = 0; qq < 4; ++qq) if (qq != c) { const double f = M[qq][c] / M[c][c]; for (int j = c; j < 5; ++j) M[qq][j] -= f * M[c][j]; }
                    }
                    ref[r] = M[3][4] / M[3][3];
                }
                const double e0 = an[0] - ref[0], e1 = an[1] - ref[1];
                const double err = std::sqrt(e0 * e0 + e1 * e1) / std::sqrt(ref[0] * ref[0] + ref[1] * ref[1]);
                ++n; worst = std::fmax(worst, err);
                if (kpp < 1.0) ++nLo; else ++nHi;
                if (sch) ++n2; else ++n1;
                // informational: the central FD of the map's kappa in w (h = 1e-5 fc), the meridian's round-off floor
                double fdv[2];
                for (int r = 0; r < 2; ++r) { const double h = 1.0e-5 * q.fc; double wp[3], wm[3];
                    for (int b = 0; b < 3; ++b) { wp[b] = w[b] + h * dv[r][b]; wm[b] = w[b] - h * dv[r][b]; }
                    fdv[r] = (returnMapHardening(wp, q, kpp).kp - returnMapHardening(wm, q, kpp).kp) / (2.0 * h); }
                worstMapFD = std::fmax(worstMapFD, std::sqrt((an[0] - fdv[0]) * (an[0] - fdv[0]) + (an[1] - fdv[1]) * (an[1] - fdv[1]))
                                                   / std::sqrt(fdv[0] * fdv[0] + fdv[1] * fdv[1]));
            }
            double Dt[6][6], sn[6], kn; returnMapTensor1(q, sp, dd, kpp, true, sn, kn, Dt, false);
            for (int i = 0; i < 6; ++i) sp[i] = sn[i]; kpp = kn;
        }
    }
    check(n >= 100 && nLo >= 30 && nHi >= 30 && n1 >= 30 && n2 >= 30 && worst < 1.0e-7,
          label(tag, "D12 kappa row on the compressive meridian (pieces of n = 1 and 2, kappa_n < 1 and >= 1): d kappa/dw == exact-Lode IFT to 1e-7").c_str());
    std::printf("       (pieces compared %d, kappa_n < 1: %d, >= 1: %d, of n = 1: %d, of n = 2: %d; worst %.2e; the map's own FD at h = 1e-5 fc reads %.2e)\n",
                n, nLo, nHi, n1, n2, worst, worstMapFD);
}

int main(int argc, char** argv) {
    run_identities();
    const char* fix = argc > 1 ? argv[1] : "tests/_testbed/concrete3d_oracle_fixture.txt";
    run_oracle_dump(fix);
    run_robustness();
    run_tension_return();
    run_f2_chain();
    run_f2_cap(1.0, 1.0, "");
    run_f2_apex(1.0, 1.0, "");
    run_f2_cap(1.0e3, 1.0e-3, "kPa-m");      // units sweep (F2 review r2): the same cases in kPa and in Pa (SI)
    run_f2_apex(1.0e3, 1.0e-3, "kPa-m");
    run_f2_cap(1.0e6, 1.0e-3, "Pa-m");
    run_f2_apex(1.0e6, 1.0e-3, "Pa-m");
    run_f2_viscous_kp(1.0, 1.0, "");
    run_f2_viscous_kp(1.0e3, 1.0e-3, "kPa-m");
    run_f2_viscous_kp(1.0e6, 1.0e-3, "Pa-m");
    run_proj_viscous(1.0, 1.0, "");
    run_proj_viscous(1.0e3, 1.0e-3, "kPa-m");
    run_proj_viscous(1.0e6, 1.0e-3, "Pa-m");
    run_f2_driver_step(1.0, 1.0, "");
    run_f2_driver_step(1.0e3, 1.0e-3, "kPa-m");
    run_f2_driver_step(1.0e6, 1.0e-3, "Pa-m");
    run_rescue_chain(1.0, "");
    run_rescue_chain(1.0e3, "kPa-m");
    run_rescue_chain(1.0e6, "Pa-m");
    run_cut_chain(1.0, 1.0, "");
    run_cut_chain(1.0e3, 1.0e-3, "kPa-m");
    run_cut_chain(1.0e6, 1.0e-3, "Pa-m");
    run_meridian_kappa(1.0, "");
    run_meridian_kappa(1.0e3, "kPa-m");
    run_meridian_kappa(1.0e6, "Pa-m");
    std::printf(fails ? "\nKERNEL CHECK: %d FAIL\n" : "\nKERNEL CHECK: ALL PASS\n", fails);
    return fails ? 1 : 0;
}
