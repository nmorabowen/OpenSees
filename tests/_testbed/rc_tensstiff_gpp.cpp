// Standalone g++ gate for Phase-3 tension stiffening in LadrunoRCKernel.h.
// Mirrors tests/_testbed/rc_shell_ref.py run_T1 + an FD tangent check on the
// floored direction. No OpenSees: header-only kernel only.
//   g++ -std=c++17 -I SRC/material/nD rc_tensstiff_gpp.cpp -o rc_ts && ./rc_ts
#include "LadrunoRCKernel.h"
#include <cstdio>
#include <cmath>
using namespace ladruno_rc_kernel;

static void fill(Params& P, int tsMode)
{
    P.E = 30000.0; P.nu = 0.2; P.Kc = 0.667;
    P.cdf = 0.0; P.eta = 0.0; P.betaFloor = 0.1;
    P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
    P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4;
    P.aggSize = 16.0; P.crackStrain = 0.0; P.crackSpacing = 0.0; P.lch = 0.0; P.betaSrMin = 0.01;
    P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
    P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
    P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
    P.tensStiffMode = tsMode; P.tensStiffC = 500.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
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

int main()
{
    Params Poff, Pon; fill(Poff, 0); fill(Pon, 1);
    RCHist hOff, hOn; histZero(hOff); histZero(hOn);
    int nfloor = 0, nbad = 0;
    int n = 400;
    for (int k = 1; k <= n; ++k) {
        double exx = 0.004 * k / n;
        double eps[6] = { exx, 0,0,0,0,0 };
        double sOff[6], sOn[6], D[6][6]; RCHist oOff, oOn;
        returnMap3D(Poff, eps, hOff, sOff, D, oOff);
        returnMap3D(Pon,  eps, hOn,  sOn,  D, oOn);
        hOff = oOff; hOn = oOn;
        double e1 = oOn.eps1, ds;
        double sts = tensStiff(e1, 1, Pon.ftPeak, Pon.tensStiffC, Pon.tensStiffAlpha, &ds);
        double ecr = Pon.crackStrain;
        if (exx < ecr) {                          // (a) pre-crack untouched
            if (fabs(sOn[0] - sOff[0]) > 1e-9 * (fabs(sOff[0]) + 1.0)) { nbad++; }
        } else {
            if (sOn[0] + 1e-9 < sOff[0]) nbad++;  // never below bare
            if (sts > sOff[0] + 1e-9) {           // (b) floor binds -> closed form
                nfloor++;
                if (fabs(sOn[0] - sts) > 1e-6 * (fabs(sts) + 1.0)) {
                    nbad++;
                    if (nbad < 5) printf("  FLOOR MISS exx=%.5f sOn=%.5f sts=%.5f\n", exx, sOn[0], sts);
                }
            }
        }
    }
    printf("floor-binding steps (closed form): %d ; mismatches: %d\n", nfloor, nbad);

    // FD tangent at a floored state (uniaxial tension, p1 = x). The baseline material uses a
    // W_B secant tangent that intentionally omits d(dt_bar)/deps, so FD!=analytic in the
    // softening regime for BOTH on and off. Two meaningful checks: (1) the TS-PINNED direction
    // D[0][0] = dsig0/deps0 must match FD tightly (sig0 == sigma_ts(e1=exx) is smooth, the pin
    // removes the bare-damage dependence); (2) TS must not DEGRADE tangent quality vs baseline
    // (worst-rel-err on ~ off).
    auto worstFD = [](Params& P, double eps[6], double& d00) {
        RCHist h; histZero(h);
        double s0[6], D[6][6]; RCHist o;
        returnMap3D(P, eps, h, s0, D, o, true);
        d00 = D[0][0];
        double worst = 0.0;
        for (int b = 0; b < 6; ++b) {
            double hp = 1e-7;
            double ep[6]; for (int a=0;a<6;a++) ep[a]=eps[a]; ep[b]+=hp;
            double sp[6], Dd[6][6]; RCHist od;
            returnMap3D(P, ep, h, sp, Dd, od, false);
            for (int a = 0; a < 6; ++a) {
                double fd=(sp[a]-s0[a])/hp, an=D[a][b];
                double rel=fabs(fd-an)/(fabs(an)+fabs(fd)+1.0);
                if (rel>worst) worst=rel;
            }
        }
        return worst;
    };
    {
        double eps[6] = { 0.0025, 0,0,0,0,0 };   // post-crack, floor binding
        double d00on, d00off;
        Params Pon2; fill(Pon2, 1); double won = worstFD(Pon2, eps, d00on);
        Params Poff2; fill(Poff2, 0); double woff = worstFD(Poff2, eps, d00off);
        // analytic ts_dsig at e1 = exx
        double ds; tensStiff(0.0025, 1, Pon2.ftPeak, Pon2.tensStiffC, Pon2.tensStiffAlpha, &ds);
        // FD of the pinned direction sig0 = sigma_ts(exx)
        double s0a[6], Dx[6][6]; RCHist ha; histZero(ha); RCHist oa;
        returnMap3D(Pon2, eps, ha, s0a, Dx, oa, false);
        double ep[6]={0.0025+1e-7,0,0,0,0,0}; double s0b[6]; RCHist ob;
        returnMap3D(Pon2, ep, ha, s0b, Dx, ob, false);
        double fd00 = (s0b[0]-s0a[0])/1e-7;
        printf("pinned D[0][0]: analytic=%.4f  FD=%.4f  ts_dsig=%.4f\n", d00on, fd00, ds);
        printf("worst FD rel-err: on=%.3e  off=%.3e (baseline secant in softening)\n", won, woff);
        if (fabs(d00on - fd00) > 1e-3*(fabs(fd00)+1.0)) { printf("PINNED TANGENT FAIL\n"); return 1; }
        if (won > woff + 1e-6) { printf("TS DEGRADES TANGENT FAIL\n"); return 1; }
    }

    // Equibiaxial (degenerate) floor: drive exx==eyy post-crack; BOTH in-plane normals
    // must be pinned to sigma_ts(e1=exx) (the degen-branch fix; the old 0.5-coeff reached
    // only halfway).
    {
        Params P; fill(P, 1); RCHist h; histZero(h);
        int nbe = 0; double maxerr = 0.0;
        for (int k = 1; k <= 300; ++k) {
            double e = 0.004 * k / 300;
            double eps[6] = { e, e, 0,0,0,0 };
            double s[6], D[6][6]; RCHist o;
            returnMap3D(P, eps, h, s, D, o); h = o;
            if (e < P.crackStrain) continue;
            double ds; double sts = tensStiff(e, 1, P.ftPeak, P.tensStiffC, P.tensStiffAlpha, &ds);
            // floor binds when sigma_ts exceeds the (equal) bare in-plane normals
            if (sts > s[0] - 1e-9 && sts > 0.05) {       // crude "binding" proxy at equibiaxial
                double e0 = fabs(s[0]-sts)/(fabs(sts)+1.0);
                double e1e = fabs(s[1]-sts)/(fabs(sts)+1.0);
                if (e0 < 1e-6 && e1e < 1e-6) nbe++;
                if (e0 > maxerr) maxerr = e0; if (e1e > maxerr) maxerr = e1e;
            }
        }
        printf("equibiaxial: pinned-both steps=%d  worst rel-err=%.2e\n", nbe, maxerr);
        if (nbe < 50) { printf("EQUIBIAXIAL FLOOR FAIL (degen under-delivers?)\n"); return 1; }
    }

    // Degenerate (equibiaxial) consistent-tangent check. FD is ill-defined at the degenerate
    // point (any single-component perturbation breaks exx==eyy and leaves the degenerate branch),
    // so instead ANALYTICALLY cross-check the TS contribution: D_on - D_bare must equal the
    // degen-branch formula  ts_inj_a * (ts_dsig*P1eps[c] - row[c]),  ts_inj=(1,1,0),
    // P1eps=(0.5,0.5,0,0,0,0), row[c] = 0.5*(D_bare[0][c]+D_bare[1][c]) (= ts_meas . D_bare).
    // This catches a wrong degen P1eps / ts_inj / ts_meas that the FD uniaxial leg cannot.
    {
        double e = 0.0025;                            // post-crack, floor binding at equibiaxial
        double eps[6] = { e, e, 0,0,0,0 };
        Params Pon; fill(Pon, 1);  RCHist hon;  histZero(hon);
        Params Poff; fill(Poff, 0); RCHist hoff; histZero(hoff);
        double son[6], Don[6][6], soff[6], Dbare[6][6]; RCHist oon, ooff;
        returnMap3D(Pon,  eps, hon,  son,  Don,   oon, true);
        returnMap3D(Poff, eps, hoff, soff, Dbare, ooff, true);
        double ds; tensStiff(e, 1, Pon.ftPeak, Pon.tensStiffC, Pon.tensStiffAlpha, &ds);
        double P1eps[6] = { 0.5, 0.5, 0, 0, 0, 0 };
        const int vidx[3] = { 0, 1, 3 };
        double tinj[3] = { 1.0, 1.0, 0.0 };
        double worst = 0.0;
        for (int a = 0; a < 3; ++a) {
            int r = vidx[a];
            for (int c = 0; c < 6; ++c) {
                double row = 0.5 * (Dbare[0][c] + Dbare[1][c]);     // ts_meas . D_bare
                double expected = Dbare[r][c] + tinj[a] * (ds * P1eps[c] - row);
                double rel = fabs(Don[r][c] - expected) / (fabs(expected) + fabs(Don[r][c]) + 1.0);
                if (rel > worst) worst = rel;
            }
        }
        printf("equibiaxial degen tangent analytic cross-check: worst rel-err=%.2e\n", worst);
        if (worst > 1e-9) { printf("DEGEN TANGENT FAIL\n"); return 1; }
    }

    if (nbad == 0 && nfloor > 10) { printf("GPP TENSSTIFF GATE: PASS\n"); return 0; }
    printf("GPP TENSSTIFF GATE: FAIL\n");
    return 1;
}
