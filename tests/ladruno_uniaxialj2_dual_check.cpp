// WP-262: the extracted LadrunoUniaxialJ2 kernel (SRC/material/uniaxial/
// LadrunoUniaxialJ2Kernel.h) checked three ways, OpenSees-free:
//
//  BITS   the kernel in double vs a VERBATIM transcription of the pre-WP-262
//         LadrunoUniaxialJ2::setTrialStrain (legacyMap below): stress, tangent
//         and every trial-state field must agree bit for bit on every step.
//  TAN    the analytic consistent tangent vs a Dual<1> seeded on the strain.
//  SENS   the DDM history machinery: per parameter (E, sig0, Qinf, b, Hiso), the
//         stress sensitivity propagated step by step exactly as the material does
//         it (conditional stress sensitivity, then commitSensitivity with the
//         total strain gradient) vs a central FD of the WHOLE path re-run at
//         theta +/- h. The strain path depends on theta (eps_n = eps0_n +
//         (theta - theta0) c_n) so the strain-gradient leg of commitSensitivity
//         is exercised.
//
//  NOHIST the same DDM with the history sensitivities dropped every step (the
//         break-the-gate leg: it must miss, or the history storage is untested).
//
// Build: g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/uniaxial
//            tests/ladruno_uniaxialj2_dual_check.cpp -o luj2d

#include <LadrunoUniaxialJ2Kernel.h>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace ladruno_uniaxial_j2;
typedef ladruno_dual::Dual<1> D1;
typedef ParamsT<double> Pd;

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// --------------------------------------------------------------------------
// Verbatim transcription of LadrunoUniaxialJ2::setTrialStrain as it stood on
// origin/ladruno @ 757c6845d (members -> struct fields; opserr removed; the
// DBL_EPSILON early-out kept by the caller, not here). DO NOT "tidy" it: it is
// the bit-identity reference.
// --------------------------------------------------------------------------
struct Legacy {
  double E, sig0, Qinf, bIso, Hiso; int nBack; double Ckin[MAXBACK], gKin[MAXBACK];
  bool dmgOn; double dmgR, dmgS, dmgPD, dmgDc; bool dmgImplex; double lchScale;
  double CplasticStrain, Cebarp, CX[MAXBACK], CdGamma, CD, CdD, Cwp;
  double TplasticStrain, Tebarp, TX[MAXBACK], TdGamma, Twp, TD;
  double Tstrain, Tstress, Ttangent;
  double yieldStress(double p) const { return Ladruno::yieldStressVoceLinear(p, sig0, Qinf, bIso, Hiso); }
  double yieldSlope(double p) const { return Ladruno::yieldSlopeVoceLinear(p, sig0, Qinf, bIso, Hiso); }

  void legacyMap(double strain)
  {
    Tstrain = strain;
    double sigTr = E*(Tstrain - CplasticStrain);
    double Xtot = 0.0;
    for (int k = 0; k < nBack; k++) Xtot += CX[k];
    double eta = sigTr - Xtot;
    double absEta = fabs(eta);
    double syn = yieldStress(Cebarp);
    double f = absEta - syn;
    double stressScale = fmax(absEta, syn);
    double tolY = 1.0e-12 * fmax(stressScale, 1.0e-300);
    if (f <= tolY) {
      TplasticStrain = CplasticStrain;
      Tebarp = Cebarp;
      TdGamma = 0.0;
      Twp = Cwp;
      for (int k = 0; k < nBack; k++) TX[k] = CX[k];
      Tstress = sigTr;
      Ttangent = E;
      TD = CD;
      if (dmgOn) {
        double dScale = dmgImplex ? (CD + CdD) : CD;
        if (dScale > dmgDc) dScale = dmgDc;
        if (dScale < 0.0)   dScale = 0.0;
        if (dScale > 0.0) { double omd = 1.0 - dScale; Tstress *= omd; Ttangent *= omd; }
      }
      return;
    }
    int s = (eta >= 0.0) ? 1 : -1;
    double h0 = yieldSlope(Cebarp);
    for (int k = 0; k < nBack; k++) h0 += Ckin[k];
    double dp = f / (E + h0);
    if (dp < 0.0) dp = 0.0;
    const int maxIter = 50;
    int iter = 0;
    double R = 1.0;
    for (; iter < maxIter; iter++) {
      double R0 = absEta - E*dp - yieldStress(Cebarp + dp);
      double dR = -E - yieldSlope(Cebarp + dp);
      for (int k = 0; k < nBack; k++) {
        double ak = Ckin[k] - gKin[k]*s*CX[k];
        double den = 1.0 + gKin[k]*dp;
        R0 -= dp*ak/den;
        dR -= ak/(den*den);
      }
      R = R0;
      if (fabs(R) <= tolY) break;
      if (dR == 0.0) break;
      dp -= R/dR;
      if (dp < 0.0) dp = 0.0;
    }
    double pbar = Cebarp + dp;
    Tebarp = pbar;
    TdGamma = dp;
    TplasticStrain = CplasticStrain + dp*s;
    for (int k = 0; k < nBack; k++)
      TX[k] = (CX[k] + Ckin[k]*dp*s) / (1.0 + gKin[k]*dp);
    Tstress = sigTr - E*dp*s;
    Twp = Cwp + Tstress*(dp*s);
    double h = yieldSlope(pbar);
    for (int k = 0; k < nBack; k++) {
      double ak = Ckin[k] - gKin[k]*s*CX[k];
      double den = 1.0 + gKin[k]*dp;
      h += ak/(den*den);
    }
    Ttangent = E*h/(E + h);
    TD = CD;
    if (dmgOn) {
      const double sigEff  = Tstress;
      const double EalgEff = Ttangent;
      Ladruno::DamageParams dpar;
      dpar.on = dmgOn; dpar.r = dmgR; dpar.s = dmgS; dpar.pD = dmgPD; dpar.Dc = dmgDc;
      const double Y = (E > 0.0) ? sigEff*sigEff/(2.0*E) : 0.0;
      const double dDraw = lchScale * Ladruno::damageIncrement(Y, Cebarp, pbar, dpar);
      double Dnew = CD + dDraw;
      bool capped = false;
      if (Dnew >= dmgDc) { Dnew = dmgDc; capped = true; }
      if (Dnew < 0.0) Dnew = 0.0;
      TD = Dnew;
      double dScale = dmgImplex ? (CD + CdD) : Dnew;
      if (dScale > dmgDc) dScale = dmgDc;
      if (dScale < 0.0)   dScale = 0.0;
      const double omd = 1.0 - dScale;
      double dDdeps = 0.0;
      if (!dmgImplex && !capped && dDraw > 0.0) {
        double ddD_dY, ddD_dp;
        Ladruno::damageIncrementPartials(Y, Cebarp, pbar, dpar, ddD_dY, ddD_dp);
        const double dpbar_deps = s*E/(E + h);
        const double dY_deps    = (sigEff/E)*EalgEff;
        dDdeps = lchScale * (ddD_dp*dpbar_deps + ddD_dY*dY_deps);
      }
      Tstress  = omd*sigEff;
      Ttangent = omd*EalgEff - sigEff*dDdeps;
    }
  }
  void commit()
  {
    CplasticStrain = TplasticStrain; Cebarp = Tebarp; CdGamma = TdGamma; Cwp = Twp;
    CdD = TD - CD; CD = TD;
    for (int k = 0; k < nBack; k++) CX[k] = TX[k];
  }
};

// ---- a material case: constants + a strain path ----------------------------
struct Case {
  const char* name;
  Pd p;
  int nstep;
  double eps[16];      // eps0_n
  double c[16];        // d eps_n / d theta (strain-path coupling for SENS)
};

static Legacy makeLegacy(const Pd& p)
{
  Legacy L;
  std::memset(&L, 0, sizeof(L));
  L.E = p.E; L.sig0 = p.sig0; L.Qinf = p.Qinf; L.bIso = p.bIso; L.Hiso = p.Hiso; L.nBack = p.nBack;
  for (int k = 0; k < MAXBACK; k++) { L.Ckin[k] = p.C[k]; L.gKin[k] = p.gam[k]; }
  L.dmgOn = p.dmgOn; L.dmgR = p.dmgR; L.dmgS = p.dmgS; L.dmgPD = p.dmgPD; L.dmgDc = p.dmgDc;
  L.dmgImplex = p.dmgImplex; L.lchScale = p.lchScale;
  return L;
}

// BITS + TAN over the path
static void bitsAndTangent(const Case& cs)
{
  Legacy L = makeLegacy(cs.p);
  StateT<double> C;
  std::memset(&C, 0, sizeof(C));
  double CdD = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];
    L.legacyMap(eps);

    StateT<double> Tr; double stress, tangent;
    returnMap(cs.p, eps, C, CdD, Tr, stress, tangent);

    bool bits = sameBits(stress, L.Tstress) && sameBits(tangent, L.Ttangent)
             && sameBits(Tr.plasticStrain, L.TplasticStrain) && sameBits(Tr.ebarp, L.Tebarp)
             && sameBits(Tr.dGamma, L.TdGamma) && sameBits(Tr.wp, L.Twp) && sameBits(Tr.D, L.TD);
    for (int k = 0; k < cs.p.nBack; k++) bits = bits && sameBits(Tr.X[k], L.TX[k]);

    // TAN: Dual<1> on the strain, constants double
    StateT<D1> Cd, Td; D1 CdDd(CdD), sd, td;
    Cd.plasticStrain = C.plasticStrain; Cd.ebarp = C.ebarp; Cd.dGamma = C.dGamma; Cd.wp = C.wp; Cd.D = C.D;
    for (int k = 0; k < MAXBACK; k++) Cd.X[k] = C.X[k];
    D1 e = D1::variable(eps, 0);
    returnMap(cs.p, e, Cd, CdDd, Td, sd, td);
    const double tanErr = std::fabs(sd.d[0] - tangent) / std::fmax(std::fabs(cs.p.E), 1.0);

    printf("BITS %s %d %d PLASTIC %d DAMAGE %.6e TAN %.3e\n", cs.name, s, bits ? 1 : 0,
           Tr.dGamma > 0.0 ? 1 : 0, Tr.D, tanErr);

    // commit (kernel copy and legacy in lock step)
    CdD = Tr.D - C.D;
    C = Tr;
    L.commit();
  }
}

// one full path at constants p (double), strain eps0 + dtheta*c; returns stresses
static void pathStresses(const Case& cs, const Pd& p, double theta, double* out)
{
  StateT<double> C; std::memset(&C, 0, sizeof(C));
  double CdD = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    StateT<double> Tr; double stress, tangent;
    returnMap(p, cs.eps[s] + theta*cs.c[s], C, CdD, Tr, stress, tangent);
    out[s] = stress;
    CdD = Tr.D - C.D; C = Tr;
  }
}

// SENS: dual DDM propagation vs FD of the whole path, for parameter `which`
//   1 E, 2 sig0, 3 Qinf, 4 bIso, 5 Hiso
static void sensitivity(const Case& cs, int which, const char* pname, bool keepHistory)
{
  // ---- dual DDM, exactly the material's sequence ----
  ParamsT<D1> pd;
  pd.E = cs.p.E; pd.sig0 = cs.p.sig0; pd.Qinf = cs.p.Qinf; pd.bIso = cs.p.bIso; pd.Hiso = cs.p.Hiso;
  pd.nBack = cs.p.nBack;
  for (int k = 0; k < MAXBACK; k++) { pd.C[k] = cs.p.C[k]; pd.gam[k] = cs.p.gam[k]; }
  pd.dmgOn = cs.p.dmgOn; pd.dmgR = cs.p.dmgR; pd.dmgS = cs.p.dmgS; pd.dmgPD = cs.p.dmgPD;
  pd.dmgDc = cs.p.dmgDc; pd.dmgImplex = cs.p.dmgImplex; pd.lchScale = cs.p.lchScale;
  D1* seed = (which == 1) ? &pd.E : (which == 2) ? &pd.sig0 : (which == 3) ? &pd.Qinf
           : (which == 4) ? &pd.bIso : &pd.Hiso;
  seed->d[0] = 1.0;
  const double theta0 = seed->v;

  StateT<D1> C; double CdD = 0.0; D1 CdDd(0.0);
  for (int k = 0; k < MAXBACK; k++) C.X[k] = 0.0;
  C.plasticStrain = 0.0; C.ebarp = 0.0; C.dGamma = 0.0; C.wp = 0.0; C.D = 0.0;
  double ddm[16];
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];          // eps_n(theta) = eps0_n + (theta - theta0) c_n
    // (a) conditional stress sensitivity: strain held fixed
    StateT<D1> T1; D1 s1, t1;
    returnMap(pd, D1(eps), C, CdDd, T1, s1, t1);
    // (b) the structure would now give d eps/d theta = c_n; total stress sensitivity
    const double tangent = t1.v;
    ddm[s] = s1.d[0] + tangent*cs.c[s];
    // (c) commitSensitivity: re-run with the total strain gradient, keep the history
    D1 e(eps); e.d[0] = cs.c[s];
    StateT<D1> T2; D1 s2, t2;
    returnMap(pd, e, C, CdDd, T2, s2, t2);
    CdDd = T2.D - C.D;
    C = T2;
    if (!keepHistory) {                    // mutant: forget d(history)/d(theta)
      C.plasticStrain.d[0] = 0.0; C.ebarp.d[0] = 0.0; C.dGamma.d[0] = 0.0;
      C.wp.d[0] = 0.0; C.D.d[0] = 0.0; CdDd.d[0] = 0.0;
      for (int k = 0; k < MAXBACK; k++) C.X[k].d[0] = 0.0;
    }
    (void)CdD;
  }

  // ---- FD of the whole path ----
  double sp[16], sm[16];
  const double h = 1.0e-6 * std::fmax(std::fabs(theta0), 1.0);
  Pd pp = cs.p, pm = cs.p;
  double* tp = (which == 1) ? &pp.E : (which == 2) ? &pp.sig0 : (which == 3) ? &pp.Qinf
             : (which == 4) ? &pp.bIso : &pp.Hiso;
  double* tm = (which == 1) ? &pm.E : (which == 2) ? &pm.sig0 : (which == 3) ? &pm.Qinf
             : (which == 4) ? &pm.bIso : &pm.Hiso;
  *tp += h; *tm -= h;
  pathStresses(cs, pp, h, sp);
  pathStresses(cs, pm, -h, sm);

  double worst = 0.0, scale = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double fd = (sp[s] - sm[s]) / (2.0*h);
    scale = std::fmax(scale, std::fabs(fd));
    worst = std::fmax(worst, std::fabs(ddm[s] - fd));
  }
  printf("%s %s %s %.3e %.6e\n", keepHistory ? "SENS" : "NOHIST", cs.name, pname,
         worst / std::fmax(scale, 1e-30), scale);
}

static Pd base(double E, double s0, double Qi, double b, double Hi, int nb,
               double C0, double g0, double C1, double g1)
{
  Pd p;
  std::memset(&p, 0, sizeof(p));
  p.E = E; p.sig0 = s0; p.Qinf = Qi; p.bIso = b; p.Hiso = Hi; p.nBack = nb;
  p.C[0] = C0; p.gam[0] = g0; p.C[1] = C1; p.gam[1] = g1;
  p.dmgOn = false; p.dmgR = 0.0; p.dmgS = 1.0; p.dmgPD = 0.0; p.dmgDc = 1.0;
  p.dmgImplex = false; p.lchScale = 1.0;
  return p;
}

int main()
{
  // strain paths (cyclic, into plasticity both ways); c = coupling for SENS
  const double CYC[12]  = {0.001, 0.003, 0.006, 0.004, 0.0, -0.004, -0.008, -0.005, 0.0, 0.005, 0.009, 0.012};
  const double CPL[12]  = {1e-6, 2e-6, -1e-6, 3e-6, 1e-6, -2e-6, 2e-6, 1e-6, -1e-6, 2e-6, 3e-6, 1e-6};

  Case cases[6];
  int n = 0;
  auto add = [&](const char* name, const Pd& p) {
    Case& c = cases[n++];
    c.name = name; c.p = p; c.nstep = 12;
    for (int i = 0; i < 12; i++) { c.eps[i] = CYC[i]; c.c[i] = CPL[i]; }
  };
  add("iso_voce",    base(200000.0, 250.0, 120.0, 25.0, 800.0, 0, 0, 0, 0, 0));
  add("chaboche2",   base(200000.0, 250.0, 60.0, 15.0, 300.0, 2, 30000.0, 300.0, 5000.0, 30.0));
  add("kin_only",    base(200000.0, 0.0, 0.0, 0.0, 0.0, 1, 20000.0, 0.0, 0, 0));   // sig0 = 0
  Pd pdm = base(200000.0, 250.0, 60.0, 15.0, 300.0, 1, 20000.0, 200.0, 0, 0);
  pdm.dmgOn = true; pdm.dmgR = 0.05; pdm.dmgS = 1.0; pdm.dmgPD = 0.002; pdm.dmgDc = 0.99;
  add("lemaitre", pdm);
  Pd pdx = pdm; pdx.dmgImplex = true;
  add("lemaitre_implex", pdx);
  Pd pdr = pdm; pdr.lchScale = 1.7; pdr.dmgS = 2.0; pdr.dmgR = 0.12;
  add("lemaitre_reg_s2", pdr);

  for (int i = 0; i < n; i++) bitsAndTangent(cases[i]);

  const char* pn[6] = {"", "E", "sig0", "Qinf", "b", "Hiso"};
  for (int i = 0; i < n; i++)
    for (int w = 1; w <= 5; w++) {
      // kin_only has Qinf = b = Hiso = sig0 = 0: still differentiable (sig0 seeds the
      // yield surface), b is inert when Qinf = 0 -- skip it there.
      if (std::strcmp(cases[i].name, "kin_only") == 0 && w == 4) continue;
      sensitivity(cases[i], w, pn[w], true);
      sensitivity(cases[i], w, pn[w], false);
    }
  return 0;
}
