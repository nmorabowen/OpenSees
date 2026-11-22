// Standalone g++ gate for Ladruno WP-241 in LadrunoRCKernel.h: the spatial add-on plane
// (-addonPlane spatial), the damage crack gate (-crackGate damage) and the interlock tangent
// spanning all six columns. Header-only kernel, no OpenSees. Self-asserting: non-zero exit and
// no PASS line on any failure.
//   g++ -O2 -std=c++17 -I SRC/material/nD rc_spatial_gpp.cpp -o rc_sp && ./rc_sp
//
// Gates
//  G1 reduce:     an x-y strain path (eps_zz the most compressive principal, no transverse
//                 shear) gives the SAME stress and tangent with spatial and membrane
//                 (beta + tension stiffening on), to 1e-12 relative.
//  G2 invariance: spatial, pure shear gamma_xy vs gamma_xz (same history) -> tau_xy == tau_xz
//                 and sigma_yy == sigma_zz. Membrane must FAIL this (break-once evidence).
//  G3 dual:       central FD of eps1Spatial vs its P1 at random multiplicity-1 states.
//  G4 pinning:    with the TS floor binding, meas . Dtan == d sigma_ts/d e1 * P1 (analytic) and
//                 the FD of sigma_ts(e1(eps)), at random 3D states.
//  G5 interlock:  membrane -interlock (2a capped, and -cyclic capped), betaSrMin = 0: the
//                 crack-plane shear is pinned, so m_sigma . Dtan must vanish in ALL six columns
//                 (FD confirms). The pre-WP-241 tangent left column 2 (eps_zz) at its raw
//                 value; the gate prints that raw value to show it is not vacuous.
//  G6 crack gate: uniaxial compression with Poisson laterals: the strain gate latches
//                 cracked-nu (the WP-241 probe defect); the damage gate never does and adds no
//                 TS floor. Uniaxial tension past the peak: the damage gate latches.
//  G7 TS onset:   compression + shear, damage gate: the floor opens with no stress jump (scaled
//                 by the onset stress, RCHist::tsScale); forcing the full floor jumps (break-once).
//                 Uniaxial tension: onset scale 1, stress identical to the strain gate.
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace ladruno_rc_kernel;

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

// drive a strain history eps(k) = k/n * target, committing each step
static void drive(const Params& P, const double tgt[6], int n, RCHist& h, double sig[6], double D[6][6])
{
    for (int k = 1; k <= n; ++k) {
        double e[6]; for (int a = 0; a < 6; ++a) e[a] = tgt[a] * k / n;
        RCHist o; returnMap3D(P, e, h, sig, D, o); h = o;
    }
}

int main()
{
    // ---------------- G1 reduce ----------------
    {
        Params Pm, Ps; fill(Pm); fill(Ps); Ps.addonPlane = 1;
        RCHist hm, hs; histZero(hm); histZero(hs);
        double worst = 0.0;
        for (int k = 1; k <= 300; ++k) {
            double t = k / 300.0;
            double e[6] = { 0.003 * t, -0.0004 * t, -0.002 * t, 0.002 * t, 0.0, 0.0 };
            double sm[6], ss[6], Dm[6][6], Ds[6][6]; RCHist om, os;
            returnMap3D(Pm, e, hm, sm, Dm, om); returnMap3D(Ps, e, hs, ss, Ds, os);
            hm = om; hs = os;
            double sc = 0.0, dc = 0.0;
            for (int a = 0; a < 6; ++a) { sc = fmax(sc, fabs(sm[a])); for (int b = 0; b < 6; ++b) dc = fmax(dc, fabs(Dm[a][b])); }
            for (int a = 0; a < 6; ++a) {
                worst = fmax(worst, fabs(sm[a] - ss[a]) / (sc + 1e-30));
                for (int b = 0; b < 6; ++b) worst = fmax(worst, fabs(Dm[a][b] - Ds[a][b]) / (dc + 1e-30));
            }
        }
        printf("G1 reduce (spatial == membrane on x-y path): worst rel=%.2e\n", worst);
        CHECK(worst < 1e-12, "G1 spatial differs from membrane on an x-y path (%.3e)", worst);
    }

    // ---------------- G2 invariance ----------------
    {
        double worst[2] = {0, 0};
        for (int mode = 0; mode < 2; ++mode) {
            Params P; fill(P); P.addonPlane = mode;
            RCHist ha, hb; histZero(ha); histZero(hb);
            for (int k = 1; k <= 300; ++k) {
                double g = 0.004 * k / 300.0;
                double ea[6] = { 0, 0, 0, g, 0, 0 }, eb[6] = { 0, 0, 0, 0, 0, g };
                double sa[6], sb[6], D[6][6]; RCHist oa, ob;
                returnMap3D(P, ea, ha, sa, D, oa); returnMap3D(P, eb, hb, sb, D, ob);
                ha = oa; hb = ob;
                // x-y shear <-> x-z shear: swap y and z
                worst[mode] = fmax(worst[mode], rel(sa[3], sb[5]));
                worst[mode] = fmax(worst[mode], rel(sa[1], sb[2]));
                worst[mode] = fmax(worst[mode], rel(sa[0], sb[0]));
            }
        }
        printf("G2 invariance gxy vs gxz: spatial worst rel=%.2e, membrane worst rel=%.2e\n", worst[1], worst[0]);
        CHECK(worst[1] < 1e-10, "G2 spatial is not invariant (%.3e)", worst[1]);
        CHECK(worst[0] > 1e-3, "G2 membrane should be asymmetric (break-once), got %.3e", worst[0]);
    }

    // ---------------- G3 dual ----------------
    {
        srand(12345);
        double worst = 0.0; int ntested = 0;
        for (int t = 0; t < 200; ++t) {
            double e[6]; for (int a = 0; a < 6; ++a) e[a] = 0.004 * (2.0 * rand() / RAND_MAX - 1.0);
            double P1[6], inj[6], meas[6]; int m;
            double e1 = eps1Spatial(e, P1, inj, meas, &m);
            if (m != 1 || e1 <= 1e-6) continue;
            ntested++;
            for (int b = 0; b < 6; ++b) {
                double h = 1e-7, ep[6], em[6], d1[6], d2[6], d3[6]; int mm;
                for (int a = 0; a < 6; ++a) { ep[a] = e[a]; em[a] = e[a]; }
                ep[b] += h; em[b] -= h;
                double fd = (eps1Spatial(ep, d1, d2, d3, &mm) - eps1Spatial(em, d1, d2, d3, &mm)) / (2 * h);
                worst = fmax(worst, fabs(fd - P1[b]) / (fabs(P1[b]) + 1e-3));
            }
            double sc = 0.0; for (int a = 0; a < 6; ++a) sc += meas[a] * inj[a];
            CHECK(fabs(sc - 1.0) < 1e-12, "G3 meas.inj != 1 (%.3e)", sc);
        }
        printf("G3 dual (FD of eps1Spatial vs P1): states=%d worst rel=%.2e\n", ntested, worst);
        CHECK(ntested > 50, "G3 too few states");
        CHECK(worst < 1e-6, "G3 P1 is not d e1/d eps (%.3e)", worst);
        // multiplicity-2 lateral state (eps_yy = eps_zz > eps_xx): meas.inj = 1, P1 = mean
        double e[6] = { -0.002, 0.0004, 0.0004, 0, 0, 0 }, P1[6], inj[6], meas[6]; int m;
        eps1Spatial(e, P1, inj, meas, &m);
        double sc = 0.0; for (int a = 0; a < 6; ++a) sc += meas[a] * inj[a];
        printf("G3 lateral double eigenvalue: m=%d P1=(%.3f,%.3f,%.3f) meas.inj=%.3f\n", m, P1[0], P1[1], P1[2], sc);
        CHECK(m == 2 && fabs(P1[1] - 0.5) < 1e-12 && fabs(P1[2] - 0.5) < 1e-12 && fabs(sc - 1.0) < 1e-12,
              "G3 double-eigenvalue projector wrong");
    }

    // ---------------- G4 pinning ----------------
    {
        Params P; fill(P); P.addonPlane = 1; P.betaOn = false;
        srand(777);
        double worstA = 0.0, worstF = 0.0; int nbind = 0;
        for (int t = 0; t < 400 && nbind < 60; ++t) {
            double e[6];
            for (int a = 0; a < 6; ++a) e[a] = 0.003 * (2.0 * rand() / RAND_MAX - 1.0);
            e[0] = fabs(e[0]) + 0.001;                       // a tensile state
            RCHist h; histZero(h);
            double s[6], D[6][6]; RCHist o;
            returnMap3D(P, e, h, s, D, o);
            double P1[6], inj[6], meas[6]; int m;
            double e1 = eps1Spatial(e, P1, inj, meas, &m);
            if (m != 1 || e1 < P.crackStrain) continue;
            double ds; double sts = tensStiff(e1, 1, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, &ds);
            double q = 0.0; for (int a = 0; a < 6; ++a) q += meas[a] * s[a];
            if (fabs(q - sts) > 1e-9 * (fabs(sts) + 1.0)) continue;   // floor not binding here
            nbind++;
            for (int c = 0; c < 6; ++c) {
                double an = 0.0; for (int a = 0; a < 6; ++a) an += meas[a] * D[a][c];
                worstA = fmax(worstA, fabs(an - ds * P1[c]) / (fabs(ds) + 1e-9));
                double h = 1e-8, ep[6], em[6], d1[6], d2[6], d3[6]; int mm;
                for (int a = 0; a < 6; ++a) { ep[a] = e[a]; em[a] = e[a]; }
                ep[c] += h; em[c] -= h;
                double sp_ = tensStiff(eps1Spatial(ep, d1, d2, d3, &mm), 1, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, 0);
                double sm_ = tensStiff(eps1Spatial(em, d1, d2, d3, &mm), 1, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, 0);
                worstF = fmax(worstF, fabs(an - (sp_ - sm_) / (2 * h)) / (fabs(ds) + 1e-9));
            }
        }
        printf("G4 TS pinning (spatial): binding states=%d analytic worst=%.2e FD worst=%.2e\n", nbind, worstA, worstF);
        CHECK(nbind >= 20, "G4 too few binding states (%d)", nbind);
        CHECK(worstA < 1e-9, "G4 pinned row != dsig*P1 (%.3e)", worstA);
        CHECK(worstF < 1e-5, "G4 pinned row != FD (%.3e)", worstF);
    }

    // ---------------- G5 interlock all columns ----------------
    for (int cyc = 0; cyc < 2; ++cyc) {
        Params P; fill(P); P.betaOn = false; P.tensStiffMode = 0;
        P.interlockOn = true; P.interlockCyclic = (cyc == 1); P.betaSrMin = 0.0;
        // crack in x-y under tension + shear, then hold the opening and slide until capped
        RCHist h; histZero(h);
        double s[6], D[6][6];
        double tgt1[6] = { 0.002, 0.0, -0.0004, 0.0, 0.0, 0.0 };
        drive(P, tgt1, 50, h, s, D);
        int capped = 0; double worstA = 0.0, worstF = 0.0, raw2 = 0.0;
        for (int k = 1; k <= 200; ++k) {
            double e[6] = { 0.002, 0.0, -0.0004, 0.006 * k / 200.0, 0.0, 0.0 };
            RCHist o; returnMap3D(P, e, h, s, D, o);
            double c = o.crackC, sn = o.crackS, cs = c * sn, c2 = c * c, s2 = sn * sn;
            double tau_cur = cs * (s[1] - s[0]) + (c2 - s2) * s[3];
            double vmax = fabs(o.betaSr) ;
            bool isCap = (o.cracked > 0.5) && ((cyc == 0) ? (fabs(o.betaSr) < 0.999) : (fabs(fabs(o.betaSr) - 1.0) < 1e-12));
            (void)vmax;
            if (isCap && k > 20) {
                capped++;
                double msig[6] = { -cs, cs, 0.0, c2 - s2, 0.0, 0.0 };
                for (int col = 0; col < 6; ++col) {
                    double an = 0.0; for (int a = 0; a < 6; ++a) an += msig[a] * D[a][col];
                    worstA = fmax(worstA, fabs(an));
                    // FD of the pinned crack shear (committed history h, trial only)
                    double hh = 1e-8, ep[6], em[6], sp_[6], sm_[6], Dd[6][6]; RCHist oo;
                    for (int a = 0; a < 6; ++a) { ep[a] = e[a]; em[a] = e[a]; }
                    ep[col] += hh; em[col] -= hh;
                    returnMap3D(P, ep, h, sp_, Dd, oo, false);
                    double tp = cs * (sp_[1] - sp_[0]) + (c2 - s2) * sp_[3];
                    returnMap3D(P, em, h, sm_, Dd, oo, false);
                    double tm = cs * (sm_[1] - sm_[0]) + (c2 - s2) * sm_[3];
                    if (col == 2) worstF = fmax(worstF, fabs((tp - tm) / (2 * hh)));
                }
                // raw (pre-WP-241) column-2 sensitivity of the smeared crack shear: interlock off
                Params Q = P; Q.interlockOn = false; Q.interlockCyclic = false;
                double sq[6], Dq[6][6]; RCHist oq; returnMap3D(Q, e, h, sq, Dq, oq);
                double r2 = 0.0;
                for (int a = 0; a < 6; ++a) r2 += msig[a] * Dq[a][2];
                raw2 = fmax(raw2, fabs(r2));
            }
            (void)tau_cur;
            h = o;
        }
        printf("G5 interlock (%s, capped steps=%d): |m_sigma.Dtan| worst over 6 cols=%.2e, FD d tau/d eps_zz=%.2e, raw col-2=%.2e\n",
               cyc ? "cyclic" : "2a", capped, worstA, worstF, raw2);
        CHECK(capped >= 10, "G5 too few capped steps (%d)", capped);
        CHECK(worstA < 1e-6 * P.E, "G5 crack-shear tangent not removed in all columns (%.3e)", worstA);
        CHECK(worstF < 1e-6 * P.E, "G5 FD says the pinned shear depends on eps_zz (%.3e)", worstF);
        CHECK(raw2 > 1e-3 * P.E, "G5 raw column-2 sensitivity is ~0: the gate would be vacuous (%.3e)", raw2);
    }

    // ---------------- G6 crack gate ----------------
    {
        int latchS = -1, latchD = -1, tsD = 0;
        for (int gate = 0; gate < 2; ++gate) {
            for (int sp = 0; sp < 2; ++sp) {
                Params P; fill(P); P.crackedNu = 0.0; P.crackGate = gate; P.addonPlane = sp;
                RCHist h; histZero(h);
                for (int k = 1; k <= 200; ++k) {
                    double ex = -0.0018 * k / 200.0;      // up to ~0.9 fc, laterals by Poisson
                    double e[6] = { ex, -P.nu * ex, -P.nu * ex, 0, 0, 0 };
                    double s[6], D[6][6]; RCHist o; returnMap3D(P, e, h, s, D, o);
                    if (o.nuCracked > 0.5) { if (gate == 0 && latchS < 0) latchS = k; if (gate == 1 && latchD < 0) latchD = k; }
                    if (gate == 1 && (fabs(s[1]) > 1e-9 || fabs(s[2]) > 1e-9)) tsD++;   // lateral stress = floor
                    h = o;
                }
            }
        }
        printf("G6 compression: strain gate latches at step %d (of 200); damage gate latch step %d, lateral-stress steps %d\n",
               latchS, latchD, tsD);
        CHECK(latchS > 0, "G6 strain gate did not reproduce the Poisson latch (the probe defect)");
        CHECK(latchD < 0, "G6 damage gate latched in plain compression");
        CHECK(tsD == 0, "G6 damage gate injected a lateral TS stress in compression");
        // tension past the peak: the damage gate latches (one step after the crack)
        Params P; fill(P); P.crackedNu = 0.0; P.crackGate = 1; P.addonPlane = 1;
        RCHist h; histZero(h); int latchT = -1;
        for (int k = 1; k <= 100; ++k) {
            double e[6] = { 0.0006 * k / 100.0, 0, 0, 0, 0, 0 };
            double s[6], D[6][6]; RCHist o; returnMap3D(P, e, h, s, D, o);
            if (o.nuCracked > 0.5 && latchT < 0) latchT = k;
            h = o;
        }
        printf("G6 tension: damage gate latches at step %d (peak strain at step ~%d)\n", latchT,
               (int)ceil(backbonePeakAbscissa(P.ht) / 0.0006 * 100.0));
        CHECK(latchT > 0, "G6 damage gate never latched in tension");
    }

    // ---------------- G7 damage-gated TS onset is continuous ----------------
    // Compression + shear (spatial): the damage gate opens below ft, where the bare principal
    // stress is under sigma_ts(e1). A fine eps_yy scan through the onset, from one committed
    // state: the largest step-to-step change of sigma_yy with TS on must be no larger than with
    // TS off (plus 5 %). Break-once: forcing tsScale = 1 (the full floor at onset, the pre-fix
    // rule) must jump by more than 10x that.
    {
        Params P; fill(P); P.crackGate = 1; P.addonPlane = 1;
        RCHist h0; histZero(h0);
        { double e[6] = { -3e-4, 0, 0, 2e-4, 0, 0 }, s[6], D[6][6]; RCHist o; returnMap3D(P, e, h0, s, D, o); h0 = o; }
        double jump[3] = { 0, 0, 0 }; int opened = 0;
        for (int mode = 0; mode < 3; ++mode) {             // 0 TS on, 1 TS off, 2 TS on forced s = 1
            Params Q = P; if (mode == 1) Q.tensStiffMode = 0;
            RCHist h = h0;                                 // every point is a trial from h0
            if (mode == 2) h.tsScale = 1.0;
            double prev = 0.0;
            for (int i = 0; i <= 240; ++i) {
                double e[6] = { -3e-4, 4.0e-5 + i * 1.0e-7, 6.152e-5, 4e-4, 0, 0 }, s[6], D[6][6]; RCHist o;
                returnMap3D(Q, e, h, s, D, o);
                if (i > 0) jump[mode] = fmax(jump[mode], fabs(s[1] - prev));
                if (mode == 0 && o.tsScale >= 0.0) opened++;
                prev = s[1];
            }
        }
        printf("G7 TS onset (compression+shear, damage gate): max |d sigma_yy| per 1e-7 step: TS on %.3e, TS off %.3e, forced full floor %.3e (onset crossings %d)\n",
               jump[0], jump[1], jump[2], opened);
        CHECK(opened > 0, "G7 the scan never crossed the damage-gated TS onset (vacuous)");
        CHECK(jump[0] <= 1.05 * jump[1] + 1e-9, "G7 TS onset jumps (%.3e vs %.3e TS off)", jump[0], jump[1]);
        CHECK(jump[2] > 10.0 * jump[1], "G7 forced full floor does not jump: the gate would be vacuous (%.3e)", jump[2]);
        // uniaxial tension: onset at the peak, where q = ft >= sigma_ts -> s = 1, and the damage
        // gate's stress equals the strain gate's (same crack strain ft/E here) once open
        double worst = 0.0; double sc = -2.0;
        for (int sp = 0; sp < 2; ++sp) {
            Params A = P, B = P; A.addonPlane = B.addonPlane = sp; A.crackGate = 1; B.crackGate = 0;
            RCHist ha, hb; histZero(ha); histZero(hb);
            for (int k = 1; k <= 100; ++k) {
                double e[6] = { 0.0008 * k / 100.0, 0, 0, 0, 0, 0 }, sa[6], sb[6], D[6][6]; RCHist oa, ob;
                returnMap3D(A, e, ha, sa, D, oa); returnMap3D(B, e, hb, sb, D, ob);
                if (oa.tsScale >= 0.0) { sc = oa.tsScale; for (int a = 0; a < 6; ++a) worst = fmax(worst, rel(sa[a], sb[a])); }
                ha = oa; hb = ob;
            }
        }
        printf("G7 uniaxial tension: onset scale %.6f, damage vs strain gate stress after onset %.2e\n", sc, worst);
        CHECK(fabs(sc - 1.0) < 1e-12, "G7 uniaxial-tension onset scale != 1 (%.6f)", sc);
        CHECK(worst < 1e-12, "G7 damage gate differs from strain gate in uniaxial tension (%.3e)", worst);
    }

    if (nfail == 0) printf("GPP SPATIAL GATE: PASS\n");
    else printf("GPP SPATIAL GATE: FAIL (%d)\n", nfail);
    return nfail == 0 ? 0 : 1;
}
