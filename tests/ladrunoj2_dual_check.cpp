// WP-262: dual-number tangent oracle for the LadrunoJ2 return-map kernel.
//
// Includes ONLY SRC/material/nD/LadrunoJ2Kernel.h (pure, OpenSees-free). For
// every step of a set of strain PATHS (threading the committed history), it runs
// the SAME templated return map twice:
//
//   * returnMap()                      -- double, the production path: stress,
//                                         state and the ANALYTIC tangent Dtan;
//   * returnMapT<Dual<6>>()            -- the strain seeded along the six
//                                         ENGINEERING strain directions, so
//                                         d(stress_I)/d(eps_eng,J) is the exact
//                                         algorithmic tangent, with no step size.
//
// and prints, per step:
//   TAN   max |Dtan - Ddual| / max |Dtan|      (the oracle: round-off level)
//   FD    the same for a central FD of the double map (what the old gates used)
//   SYM   the same against the SYMMETRIZED Dtan (mutant: must be flagged where
//         the AF cross-term makes the tangent unsymmetric)
//   BITS  1 if the Dual run's VALUE part equals the double run bit for bit
//         (stress, ebarP, dGamma, epsP, every backstress), else 0
//   PLASTIC 1 if the step yielded
//
// Build:  g++ -O2 -std=c++17 -I SRC/material/nD ladrunoj2_dual_check.cpp -o lj2d

#include <LadrunoJ2Kernel.h>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace ladruno_j2_kernel;
typedef ladruno_dual::Dual<6> D6;

struct Scenario {
  const char* name;
  Params p;
  double dir[6];
  const double* mult;
  int nstep;
  const double* path;   // nstep*6 or null
};

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

static void run(const Scenario& sc)
{
  printf("SCENARIO %s\n", sc.name);

  double epsP_n[6] = {0,0,0,0,0,0};
  double ebarP_n = 0.0;
  double alpha_n[MAXBACK][6];
  for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = 0.0;

  for (int s = 0; s < sc.nstep; s++) {
    double eps[6];
    for (int i = 0; i < 6; i++)
      eps[i] = sc.path ? sc.path[6*s + i] : sc.mult[s] * sc.dir[i];

    // ---- production double path ----
    double stress[6], Dtan[6][6], epsP[6], ebarP, alpha[MAXBACK][6], dG;
    returnMap(sc.p, eps, epsP_n, ebarP_n, alpha_n, stress, Dtan, epsP, ebarP, alpha, dG);

    // ---- the same map in Dual<6>: seed the ENGINEERING strain directions.
    //      The kernel takes TENSOR components; tensor shear = engineering / 2.
    D6 epsD[6], epsPnD[6], ebarPnD(ebarP_n), alphaNd[MAXBACK][6];
    for (int i = 0; i < 6; i++) {
      epsD[i] = D6(eps[i]);
      epsD[i].d[i] = (i < 3) ? 1.0 : 0.5;
      epsPnD[i] = D6(epsP_n[i]);
    }
    for (int k = 0; k < MAXBACK; k++) for (int i = 0; i < 6; i++) alphaNd[k][i] = D6(alpha_n[k][i]);

    D6 sD[6], DtD[6][6], epsPD[6], ebarPD, alphaD[MAXBACK][6], dGD;
    returnMapT<D6>(sc.p, epsD, epsPnD, ebarPnD, alphaNd, sD, DtD, epsPD, ebarPD, alphaD, dGD);

    // ---- oracle: analytic Dtan vs the dual derivative ----
    // MUTANT: the same comparison against the SYMMETRIZED Dtan. Under
    // non-proportional loading with Armstrong-Frederick backstress the
    // algorithmic tangent is unsymmetric (the betaMpN Mperp (x) n term), so a
    // working oracle must flag the symmetrized tangent there (break-the-gate leg).
    double scale = 0.0, errDual = 0.0, errSym = 0.0;
    for (int I = 0; I < 6; I++)
      for (int J = 0; J < 6; J++) {
        scale = std::fmax(scale, std::fabs(Dtan[I][J]));
        errDual = std::fmax(errDual, std::fabs(Dtan[I][J] - sD[I].d[J]));
        errSym = std::fmax(errSym, std::fabs(0.5*(Dtan[I][J] + Dtan[J][I]) - sD[I].d[J]));
      }

    // ---- reference: central FD of the double map from the same committed state ----
    double errFD = 0.0;
    const double h = 1.0e-7;
    for (int J = 0; J < 6; J++) {
      double ep[6], em[6];
      for (int i = 0; i < 6; i++) { ep[i] = eps[i]; em[i] = eps[i]; }
      const double dt = (J < 3) ? h : 0.5*h;        // engineering step h -> tensor comp
      ep[J] += dt; em[J] -= dt;
      double sp[6], sm[6], Dx[6][6], ex[6], bx, ax[MAXBACK][6], gx;
      returnMap(sc.p, ep, epsP_n, ebarP_n, alpha_n, sp, Dx, ex, bx, ax, gx);
      returnMap(sc.p, em, epsP_n, ebarP_n, alpha_n, sm, Dx, ex, bx, ax, gx);
      for (int I = 0; I < 6; I++)
        errFD = std::fmax(errFD, std::fabs(Dtan[I][J] - (sp[I] - sm[I])/(2.0*h)));
    }

    // ---- value part of the Dual run == double run, bit for bit ----
    bool bits = sameBits(ebarPD.v, ebarP) && sameBits(dGD.v, dG);
    for (int i = 0; i < 6; i++)
      bits = bits && sameBits(sD[i].v, stress[i]) && sameBits(epsPD[i].v, epsP[i]);
    for (int k = 0; k < sc.p.nBack; k++)
      for (int i = 0; i < 6; i++) bits = bits && sameBits(alphaD[k][i].v, alpha[k][i]);

    printf("STEP %d PLASTIC %d TAN %.3e FD %.3e SYM %.3e BITS %d SCALE %.6e\n",
           s, dG > 0.0 ? 1 : 0, errDual/scale, errFD/scale, errSym/scale, bits ? 1 : 0, scale);

    // commit
    for (int i = 0; i < 6; i++) epsP_n[i] = epsP[i];
    ebarP_n = ebarP;
    for (int k = 0; k < sc.p.nBack; k++) for (int i = 0; i < 6; i++) alpha_n[k][i] = alpha[k][i];
  }
}

// The scenarios of tests/ladrunoj2_kernel_check.cpp, plus a strong-Voce one.
static const double MONO[6] = {0.002, 0.005, 0.01, 0.02, 0.03, 0.05};
static const double CYC[5]  = {0.03, 0.0, -0.03, 0.0, 0.03};
static const double NONPROP[8*6] = {
  0.010,0,0, 0,0,0,
  0.020,0,0, 0,0,0,
  0.030,0,0, 0,0,0,
  0.040,0,0, 0,0,0,
  0.040,0,0, 0.010,0,0,
  0.040,0,0, 0.020,0,0,
  0.040,0,0, 0.030,0,0,
  0.040,0,0, 0.040,0,0,
};

int main()
{
  const double K = 1500.0, G = 700.0;
  Scenario scen[7];
  scen[0] = {"elastic_iso",
             Params{K, G, 5.0, 4.0, 30.0, 10.0, 0, {0}, {0}},
             {1,0,0, 0,0,0}, MONO, 6, 0};
  scen[1] = {"single_af_cyclic",
             Params{K, G, 5.0, 0.0, 0.0, 0.0, 1, {400.0}, {80.0}},
             {1,0,0, 0,0,0}, CYC, 5, 0};
  scen[2] = {"chaboche3_3d",
             Params{K, G, 8.0, 3.0, 20.0, 5.0, 3, {500.0,300.0,150.0}, {100.0,50.0,10.0}},
             {0.5,-0.2,-0.3, 0.4, 0.1, 0.2}, MONO, 6, 0};
  scen[3] = {"kin_linear_cyclic",
             Params{K, G, 5.0, 0.0, 0.0, 0.0, 1, {120.0}, {0.0}},
             {1,0,0, 0,0,0}, CYC, 5, 0};
  scen[4] = {"linkin_cf",
             Params{K, G, 5.0, 0.0, 0.0, 30.0, 2, {200.0,100.0}, {0.0,0.0}},
             {0.5,-0.2,-0.3, 0.4, 0.1, 0.2}, MONO, 6, 0};
  scen[5] = {"chaboche_nonprop",
             Params{K, G, 6.0, 0.0, 0.0, 0.0, 2, {500.0,250.0}, {60.0,12.0}},
             {0,0,0, 0,0,0}, 0, 8, NONPROP};
  // 7) strong Voce saturation (b = 300) + AF, 3D shear-inclusive: several local
  //    Newton iterations per step, the case where a derivative carried through
  //    the iterations lags its value (the IFT step in returnMapT closes it).
  scen[6] = {"voce_strong_3d",
             Params{K, G, 5.0, 20.0, 300.0, 2.0, 1, {300.0}, {40.0}},
             {0.5,-0.2,-0.3, 0.4, 0.1, 0.2}, MONO, 6, 0};

  for (int s = 0; s < 7; s++) run(scen[s]);
  return 0;
}
