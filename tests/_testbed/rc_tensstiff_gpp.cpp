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

    if (nbad == 0 && nfloor > 10) { printf("GPP TENSSTIFF GATE: PASS\n"); return 0; }
    printf("GPP TENSSTIFF GATE: FAIL\n");
    return 1;
}
