// WP-265 (D-IMPLEX): the IMPL-EX tangent-identity gate, kernel level, OpenSees-free.
//
// Claim under test. An IMPL-EX material ships, for the open step, a tangent that
// must be the EXACT derivative of the stress map it actually evaluates: the
// extrapolated map, with the extrapolation "clock" latched (every extrapolated
// quantity is a function of the COMMITTED history only, so it is a constant of the
// step). If anything in the extrapolated map reads the trial strain -- a live
// damage, an unlatched step ratio -- the shipped tangent is no longer the
// derivative of the map, and Newton loses the property IMPL-EX exists to buy.
//
// Covered here: the two IMPL-EX materials whose kernels are templated on the
// scalar today (WP-262).
//
//   LadrunoJ2        ladruno_j2_kernel::returnMapDamagedT  (Lemaitre layer, 3D)
//   LadrunoUniaxialJ2 ladruno_uniaxial_j2::returnMap        (Lemaitre layer, 1D)
//
// Both extrapolate the damage only:  D~ = D_n + dD_n  (LadrunoJ2.cpp:366,
// LadrunoUniaxialJ2Kernel.h "dScale = C.D + CdD"), with dD_n the committed
// IMPLICIT increment (LadrunoJ2.cpp:463, LadrunoUniaxialJ2.cpp:304). Neither reads
// LadrunoImplexClock: the step ratio is 1 by construction (the STEP clock). The
// history threading below repeats those commit lines exactly.
//
// Per step, for each material, the check prints:
//
//   TAN    max |K_shipped - d(sigma_implex)/d(eps)| / max|K_shipped|   (the gate;
//          d/d(eps) by a Dual seeded on the strain: exact, no step size)
//   FD     the same for a central FD of the double IMPL-EX map (reference only)
//   LIVED  MUTANT 1 "live damage": the dual derivative of the map degraded by the
//          IMPLICIT damage D(eps) instead of the frozen D~ (the material forgot to
//          freeze D) against the shipped IMPL-EX tangent. Must be flagged wherever
//          damage grows on the step.
//          (Its mirror, "the fully-implicit tangent shipped under IMPL-EX", is
//          the same number: by the WP-262 oracle the implicit tangent IS the
//          derivative of the implicit map, so it is not printed twice.)
//   CLOCK  MUTANT 2 "unlatched clock": D~ = D_n + f(eps) dD_n with a strain clock
//          f = |eps - eps_n| / |eps_n - eps_{n-1}| read live at the trial (the
//          value of f is 1 on these uniform paths, so the mutant has the SAME
//          stress; only its derivative differs). Must be flagged wherever
//          dD_n > 0.
//   BITS   1 if the Dual IMPL-EX run's value part equals the double run bit for bit
//   DOUT   1 if the damage the IMPL-EX call returns for commit equals the
//          fully-implicit call's, bit for bit (IMPL-EX must commit a true
//          backward-Euler history)
//   GROW   1 if the implicit damage grows on the step (D > D_n)
//   EXTRAP D~ - D_n (> 0: the extrapolation is active on the step)
//   DAMAGE the committed (implicit) damage after the step
//
// Build: g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/uniaxial
//            tests/ladruno_implex_tangent_gate_check.cpp -o limx

#include <LadrunoJ2Kernel.h>
#include <LadrunoUniaxialJ2Kernel.h>
#include <cmath>
#include <cstdio>
#include <cstring>

typedef ladruno_dual::Dual<6> D6;
typedef ladruno_dual::Dual<1> D1;

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

static double clampD(double d, double Dc) { if (d >= Dc) d = Dc; if (d < 0.0) d = 0.0; return d; }

// ===========================================================================
//  LadrunoJ2 (3D)
// ===========================================================================
namespace j2 {
using namespace ladruno_j2_kernel;

static void seedStrain(const double eps[6], D6 epsD[6])
{
  // the kernel takes TENSOR components; seed the ENGINEERING directions
  for (int i = 0; i < 6; i++) { epsD[i] = D6(eps[i]); epsD[i].d[i] = (i < 3) ? 1.0 : 0.5; }
}

static void run(const char* name, const Params& p, const double* path, int nstep)
{
  printf("SCENARIO j2_%s\n", name);
  double epsP_n[6] = {0,0,0,0,0,0}, ebarP_n = 0.0, D_n = 0.0, dD_n = 0.0;
  double alpha_n[MAXBACK][6];
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = 0.0;
  double eps_n[6] = {0,0,0,0,0,0}, deps_n = 0.0;   // committed strain and |its last increment|

  for (int s = 0; s < nstep; s++) {
    const double* eps = path + 6*s;
    const double Dtil = D_n + dD_n;                  // LadrunoJ2.cpp:366 (clock latched: f = 1)

    // ---- production double path, IMPL-EX ----
    double sig[6], K[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG, D;
    returnMapDamaged(p, eps, epsP_n, ebarP_n, alpha_n, D_n, sig, K, epsP, ebarP, alpha, dG, D, 0, Dtil);

    // ---- production double path, fully implicit (DOUT) ----
    double sigI[6], KI[6][6], eI[6], bI, aI[MAXBACK][6], gI, DI;
    returnMapDamaged(p, eps, epsP_n, ebarP_n, alpha_n, D_n, sigI, KI, eI, bI, aI, gI, DI, 0, -1.0);

    // ---- the IMPL-EX map in Dual<6> ----
    D6 epsD[6], epsPnD[6], ebarPnD(ebarP_n), DnD(D_n), alphaNd[MAXBACK][6];
    seedStrain(eps, epsD);
    for (int i = 0; i < 6; i++) epsPnD[i] = D6(epsP_n[i]);
    for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alphaNd[k][i] = D6(alpha_n[k][i]);
    D6 sD[6], KD[6][6], epsPD[6], ebarPD, alphaD[MAXBACK][6], dGD, DD;
    returnMapDamagedT<D6>(p, epsD, epsPnD, ebarPnD, alphaNd, DnD,
                          sD, KD, epsPD, ebarPD, alphaD, dGD, DD, 0, Dtil);

    // ---- MUTANT 1: the map degraded by the live implicit damage ----
    D6 sL[6], KL[6][6], eL[6], bL, aL[MAXBACK][6], gL, DL;
    returnMapDamagedT<D6>(p, epsD, epsPnD, ebarPnD, alphaNd, DnD,
                          sL, KL, eL, bL, aL, gL, DL, 0, -1.0);

    // ---- MUTANT 2: unlatched strain clock, D~(eps) = D_n + f(eps) dD_n ----
    //      sigma = (1 - D~(eps)) sigma~(eps), sigma~ from the undamaged effective map
    D6 sE[6], KE[6][6], eE[6], bE, aE[MAXBACK][6], gE;
    returnMapT<D6>(p, epsD, epsPnD, ebarPnD, alphaNd, sE, KE, eE, bE, aE, gE);
    D6 dn2(0.0);
    for (int i = 0; i < 6; i++) {
      const double w = (i < 3) ? 1.0 : 2.0;          // engineering strain magnitude
      D6 di = w*(epsD[i] - eps_n[i]);
      dn2 += di*di;
    }
    D6 f = (deps_n > 0.0 && value(dn2) > 0.0) ? D6(sqrt(dn2)/deps_n) : D6(1.0);
    D6 DtilLive = D_n + f*dD_n;
    if (value(DtilLive) >= p.dmg.Dc) DtilLive = p.dmg.Dc;
    if (value(DtilLive) < 0.0) DtilLive = 0.0;
    D6 sC[6];
    for (int i = 0; i < 6; i++) sC[i] = (1.0 - DtilLive)*sE[i];

    // ---- compare ----
    double scale = 0.0, eTan = 0.0, eLive = 0.0, eClock = 0.0;
    for (int I = 0; I < 6; I++)
      for (int J = 0; J < 6; J++) {
        scale = std::fmax(scale, std::fabs(K[I][J]));
        eTan   = std::fmax(eTan,   std::fabs(K[I][J]  - sD[I].d[J]));
        eLive  = std::fmax(eLive,  std::fabs(K[I][J]  - sL[I].d[J]));
        eClock = std::fmax(eClock, std::fabs(K[I][J]  - sC[I].d[J]));
      }

    double eFD = 0.0;
    const double h = 1.0e-7;
    for (int J = 0; J < 6; J++) {
      double ep[6], em[6];
      for (int i = 0; i < 6; i++) { ep[i] = eps[i]; em[i] = eps[i]; }
      const double dt = (J < 3) ? h : 0.5*h;
      ep[J] += dt; em[J] -= dt;
      double sp[6], sm[6], Kx[6][6], ex[6], bx, ax[MAXBACK][6], gx, Dx;
      returnMapDamaged(p, ep, epsP_n, ebarP_n, alpha_n, D_n, sp, Kx, ex, bx, ax, gx, Dx, 0, Dtil);
      returnMapDamaged(p, em, epsP_n, ebarP_n, alpha_n, D_n, sm, Kx, ex, bx, ax, gx, Dx, 0, Dtil);
      for (int I = 0; I < 6; I++) eFD = std::fmax(eFD, std::fabs(K[I][J] - (sp[I] - sm[I])/(2.0*h)));
    }

    bool bits = sameBits(DD.v, D) && sameBits(dGD.v, dG) && sameBits(ebarPD.v, ebarP);
    for (int i = 0; i < 6; i++) bits = bits && sameBits(sD[i].v, sig[i]) && sameBits(epsPD[i].v, epsP[i]);
    const bool dout = sameBits(D, DI);

    printf("STEP %d TAN %.3e FD %.3e LIVED %.3e CLOCK %.3e BITS %d DOUT %d GROW %d "
           "EXTRAP %.6e DAMAGE %.6e\n",
           s, eTan/scale, eFD/scale, eLive/scale, eClock/scale, bits ? 1 : 0,
           dout ? 1 : 0, (D > D_n) ? 1 : 0, clampD(Dtil, p.dmg.Dc) - D_n, D);

    // ---- commit (LadrunoJ2::commitState, LadrunoJ2.cpp:459-466) ----
    double dn = 0.0;
    for (int i = 0; i < 6; i++) { const double w = (i < 3) ? 1.0 : 2.0; dn += (w*(eps[i]-eps_n[i]))*(w*(eps[i]-eps_n[i])); }
    deps_n = sqrt(dn);
    for (int i = 0; i < 6; i++) { epsP_n[i] = epsP[i]; eps_n[i] = eps[i]; }
    ebarP_n = ebarP;
    dD_n = D - D_n;
    D_n = D;
    for (int k = 0; k < p.nBack; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha[k][i];
  }
}
} // namespace j2

// ===========================================================================
//  LadrunoUniaxialJ2 (1D)
// ===========================================================================
namespace uj2 {
using namespace ladruno_uniaxial_j2;
typedef ParamsT<double> Pd;

static StateT<D1> lift(const StateT<double>& c)
{
  StateT<D1> r;
  r.plasticStrain = c.plasticStrain; r.ebarp = c.ebarp; r.dGamma = c.dGamma; r.wp = c.wp; r.D = c.D;
  for (int k = 0; k < MAXBACK; k++) r.X[k] = c.X[k];
  return r;
}

static void run(const char* name, const Pd& pi, const double* path, int nstep)
{
  printf("SCENARIO uj2_%s\n", name);
  Pd pimp = pi; pimp.dmgImplex = true;     // the shipped IMPL-EX map
  Pd plive = pi; plive.dmgImplex = false;  // fully implicit (mutant 1, DOUT)

  StateT<double> C;
  C.plasticStrain = 0.0; C.ebarp = 0.0; C.dGamma = 0.0; C.wp = 0.0; C.D = 0.0;
  for (int k = 0; k < MAXBACK; k++) C.X[k] = 0.0;
  double CdD = 0.0, eps_n = 0.0, deps_n = 0.0;

  for (int s = 0; s < nstep; s++) {
    const double eps = path[s];

    StateT<double> Tr; double sig, K;
    returnMap(pimp, eps, C, CdD, Tr, sig, K);
    StateT<double> TrI; double sigI, KI;                     // fully implicit (DOUT)
    returnMap(plive, eps, C, CdD, TrI, sigI, KI);

    const StateT<D1> Cd = lift(C);
    D1 epsD = D1::variable(eps, 0);

    StateT<D1> TdD; D1 sD, KD;
    returnMap(pimp, epsD, Cd, D1(CdD), TdD, sD, KD);

    StateT<D1> TdL; D1 sL, KL;                               // mutant 1 (live D)
    returnMap(plive, epsD, Cd, D1(CdD), TdL, sL, KL);

    // mutant 2: the kernel's own D~ = C.D + CdD with a LIVE strain clock in CdD
    D1 f = (deps_n > 0.0 && eps != eps_n) ? D1(fabs(epsD - eps_n)/deps_n) : D1(1.0);
    StateT<D1> TdC; D1 sC, KC;
    returnMap(pimp, epsD, Cd, f*CdD, TdC, sC, KC);

    const double scale = std::fabs(K);
    const double eTan = std::fabs(K - sD.d[0]), eLive = std::fabs(K - sL.d[0]),
                 eClock = std::fabs(K - sC.d[0]);

    const double h = 1.0e-8;
    StateT<double> Tx; double sp, sm, Kx;
    returnMap(pimp, eps + h, C, CdD, Tx, sp, Kx);
    returnMap(pimp, eps - h, C, CdD, Tx, sm, Kx);
    const double eFD = std::fabs(K - (sp - sm)/(2.0*h));

    const bool bits = sameBits(sD.v, sig) && sameBits(TdD.D.v, Tr.D) && sameBits(TdD.ebarp.v, Tr.ebarp)
                   && sameBits(TdD.plasticStrain.v, Tr.plasticStrain);
    const bool dout = sameBits(Tr.D, TrI.D);
    const double Dtil = clampD(C.D + CdD, pi.dmgDc);

    printf("STEP %d TAN %.3e FD %.3e LIVED %.3e CLOCK %.3e BITS %d DOUT %d GROW %d "
           "EXTRAP %.6e DAMAGE %.6e\n",
           s, eTan/scale, eFD/scale, eLive/scale, eClock/scale, bits ? 1 : 0,
           dout ? 1 : 0, (Tr.D > C.D) ? 1 : 0, Dtil - C.D, Tr.D);

    // ---- commit (LadrunoUniaxialJ2::commitState, LadrunoUniaxialJ2.cpp:296-306) ----
    deps_n = fabs(eps - eps_n);
    eps_n = eps;
    CdD = Tr.D - C.D;
    C = Tr;
  }
}
} // namespace uj2

int main()
{
  // ---- LadrunoJ2: the Lemaitre parameter sets of tests/ladrunoj2_dual_check.cpp ----
  const double K = 1500.0, G = 700.0;
  ladruno_j2_kernel::Params pd{K, G, 5.0, 3.0, 25.0, 8.0, 2, {500.0, 250.0}, {60.0, 12.0}};
  pd.dmg.on = true; pd.dmg.r = 3.0; pd.dmg.s = 1.0; pd.dmg.pD = 0.0; pd.dmg.Dc = 0.99;
  ladruno_j2_kernel::Params pd2 = pd; pd2.dmg.s = 2.0; pd2.dmg.r = 0.5;   // fast damage, D -> 0.9

  // uniform increments (the STEP clock's exact case), then unload and reload
  const double base[6] = {0.5, -0.2, -0.3, 0.4, 0.1, 0.2};
  double MIX[9*6];
  const double mmul[9] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.05, 0.06, 0.07};
  for (int s = 0; s < 9; s++)
    for (int i = 0; i < 6; i++) MIX[6*s + i] = mmul[s]*base[i]*(i >= 3 ? 0.5 : 1.0);
  double UNI[8*6];
  const double umul[8] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.04, 0.05, 0.06};
  for (int s = 0; s < 8; s++) {
    for (int i = 0; i < 6; i++) UNI[6*s + i] = 0.0;
    UNI[6*s] = umul[s];
  }
  j2::run("lemaitre_mixed", pd, MIX, 9);
  j2::run("lemaitre_fast_uniaxial", pd2, UNI, 8);

  // ---- LadrunoUniaxialJ2 ----
  uj2::Pd pu;
  pu.E = 2.0e5; pu.sig0 = 250.0; pu.Qinf = 100.0; pu.bIso = 20.0; pu.Hiso = 1000.0;
  pu.nBack = 2; pu.C[0] = 2.0e4; pu.gam[0] = 200.0; pu.C[1] = 3.0e3; pu.gam[1] = 20.0;
  for (int k = 2; k < ladruno_uniaxial_j2::MAXBACK; k++) { pu.C[k] = 0.0; pu.gam[k] = 0.0; }
  pu.dmgOn = true; pu.dmgR = 2.0; pu.dmgS = 1.0; pu.dmgPD = 0.0; pu.dmgDc = 0.99;
  pu.dmgImplex = true; pu.lchScale = 1.0;
  uj2::Pd pu2 = pu; pu2.dmgR = 0.3; pu2.dmgS = 2.0;
  const double UP[10] = {0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.035, 0.04, 0.05, 0.06};
  const double CYC[10] = {0.01, 0.02, 0.03, 0.02, 0.0, -0.01, -0.02, -0.03, -0.02, 0.0};
  uj2::run("lemaitre_tension", pu, UP, 10);
  uj2::run("lemaitre_cyclic", pu2, CYC, 10);
  return 0;
}
