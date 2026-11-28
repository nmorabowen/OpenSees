// WP-267: the extracted LadrunoCohesiveHinge kernel (SRC/material/uniaxial/
// LadrunoCohesiveHingeKernel.h) checked four ways, OpenSees-free:
//
//  BITS   the kernel in double vs a VERBATIM transcription of the pre-WP-267
//         LadrunoCohesiveHinge (computeDerived, getPMMultiplier, envelope,
//         setTrialStrain, the kappaMax part of commitState): the calibration,
//         the P-M multiplier, stress, tangent and trial kappaMax must agree bit
//         for bit on every step.
//  TAN    the tangent the material reports vs a Dual<1> seeded on the jump.
//  SENS   the DDM history machinery: per parameter, the moment sensitivity
//         propagated step by step exactly as the material does it (conditional
//         stress sensitivity, then commitSensitivity with the total jump
//         gradient) vs a central FD of the WHOLE path re-run at theta +/- h,
//         where the perturbed run re-derives the constants the way
//         updateParameter does. The jump path depends on theta
//         (k_n = k0_n + (theta - theta0) c_n) so the strain-gradient leg runs.
//  NOHIST the same DDM with the history sensitivities dropped every step (the
//         break-the-gate leg: it must miss, or the history storage is untested).
//
// Build: g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/uniaxial
//            tests/ladruno_cohesivehinge_dual_check.cpp -o lchd

#include <LadrunoCohesiveHingeKernel.h>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace ladruno_cohesive_hinge;
typedef ladruno_dual::Dual<1> D1;
typedef ParamsT<double> Pd;

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// --------------------------------------------------------------------------
// Verbatim transcription of LadrunoCohesiveHinge as it stood on origin/ladruno
// @ f183b9cd1 (members -> struct fields; opserr removed; the IMPL-EX clock
// replaced by its factor, passed in). DO NOT "tidy" it: it is the bit-identity
// reference.
// --------------------------------------------------------------------------
struct Legacy {
  enum { EXP = 0, LINEAR = 1 };
  double Mc, Gf, Kpen; bool userKpen; int shape; double penaltyRatio;
  double P0, alphaPM, P0t, Pbal, muBal; bool implex; double P_axial;
  double kappa0, aExp, kappaf;
  double CkappaMax, CkappaMax_old;
  int implexCommittedSteps;
  double Tstrain, Tstress, Ttangent, TkappaMax;

  double getPMMultiplier(void) const
  {
    if (P0 <= 0.0) return 1.0;
    // 4-parameter C^1 piecewise-parabolic RC balance-point law (-pmBal)
    if (Pbal > 0.0 && P0 > Pbal && P0t > 0.0) {
      double f = 0.0;
      if (P_axial >= 0.0 && P_axial <= P0t) {
        f = 1.0 - P_axial / P0t;
      } else if (P_axial < 0.0 && P_axial >= -Pbal) {
        double u = (P_axial + Pbal) / Pbal;
        f = 1.0 + (muBal - 1.0) * (1.0 - u * u);
      } else if (P_axial < -Pbal && P_axial >= -P0) {
        double denom = P0 - Pbal;
        double u = (P_axial + Pbal) / denom;
        f = muBal * (1.0 - u * u);
      }
      if (f < 0.05) f = 0.05;
      else if (f > 3.0) f = 3.0;
      return f;
    }
    double n = P_axial / P0;
    double f = 1.0 - alphaPM * n - (1.0 - alphaPM) * n * n;
    if (f < 0.05) f = 0.05;
    else if (f > 2.5) f = 2.5;
    return f;
  }

  void computeDerived(void)
  {
    if (Mc <= 0.0 || Gf <= 0.0) { kappa0 = aExp = kappaf = 0.0; return; }

    const double KpenFloor = Mc*Mc / (2.0*Gf);   // pre-peak energy == Gf at the floor

    if (!userKpen || Kpen <= 0.0) {
      double r = (penaltyRatio > 1.0) ? penaltyRatio : 1000.0;
      Kpen = r * KpenFloor;
    }

    if (Kpen <= KpenFloor) {
      double rSafe = (penaltyRatio > 1.0) ? penaltyRatio : 1000.0;
      double KpenSafe = rSafe * KpenFloor;
      Kpen = KpenSafe;
    }

    kappa0 = Mc / Kpen;                       // jump at activation
    double Esoft = Gf - Mc*Mc/(2.0*Kpen);     // energy left for the softening branch (>0 by floor)
    aExp   = Mc / Esoft;                       // exp: area Mc/aExp == Esoft
    kappaf = 2.0*Esoft / Mc;                   // lin: area Mc*kappaf/2 == Esoft
  }

  void envelope(double km, double& Menv, double& dMenv) const
  {
    if (km <= kappa0) {                        // pre-peak penalty (hinge closed)
      Menv  = Kpen * km;
      dMenv = Kpen;
    }
    else if (shape == LINEAR) {
      double d = km - kappa0;
      if (d >= kappaf) {                       // fully broken: zero residual moment
        Menv = 0.0; dMenv = 0.0;
      } else {
        Menv  = Mc * (1.0 - d/kappaf);
        dMenv = -Mc/kappaf;
      }
    } else {                                   // EXP
      double e = exp(-aExp * (km - kappa0));
      Menv  =  Mc * e;
      dMenv = -aExp * Mc * e;
    }
    if (P0 > 0.0) {
      double f = getPMMultiplier();
      Menv  *= f;
      dMenv *= f;
    }
  }

  // clockFactor: what implexClock.factor(dtime_n_commit) returned
  int setTrialStrain(double strain, double clockFactor)
  {
    Tstrain = strain;
    double akappa = fabs(strain);

    bool loading = (akappa >= CkappaMax);
    TkappaMax = loading ? akappa : CkappaMax;

    if (implex) {
      double time_factor = (implexCommittedSteps == 0) ? 0.0 : clockFactor;
      double dKappaMax_n = CkappaMax - CkappaMax_old;
      double kappaMax_tilde = CkappaMax + time_factor * dKappaMax_n;
      if (kappaMax_tilde < CkappaMax) kappaMax_tilde = CkappaMax;

      double Menv_tilde = 0.0, dMenv_tilde = 0.0;
      envelope(kappaMax_tilde, Menv_tilde, dMenv_tilde);
      double fpm = (P0 > 0.0) ? getPMMultiplier() : 1.0;
      double KpenEff = Kpen * fpm;
      double Ksec_tilde = (kappaMax_tilde > 0.0) ? (Menv_tilde / kappaMax_tilde) : KpenEff;
      double Kfloor = 1.0e-6 * KpenEff;
      if (Ksec_tilde < Kfloor) Ksec_tilde = Kfloor;

      Tstress  = Ksec_tilde * strain;
      Ttangent = Ksec_tilde;
      return 0;
    }

    double Menv, dMenv;
    envelope(TkappaMax, Menv, dMenv);

    double KpenEff = (P0 > 0.0) ? (Kpen * getPMMultiplier()) : Kpen;
    double Ksec = (TkappaMax > 0.0) ? (Menv / TkappaMax) : KpenEff;

    Tstress = Ksec * strain;                   // sign carried by `strain`

    if (loading && akappa > 0.0)
      Ttangent = dMenv;                        // consistent envelope slope (may be < 0)
    else
      Ttangent = Ksec;                         // unload/reload secant (>= 0)
    return 0;
  }

  void commit(void)
  {
    CkappaMax_old = CkappaMax;
    CkappaMax     = TkappaMax;
    implexCommittedSteps++;
  }
};

// --------------------------------------------------------------------------
struct Case {
  const char* name;
  Pd p;                 // the user inputs (Kpen < 0 => default); computeDerived'd in main
  bool implex;
  double clock;         // IMPL-EX step ratio (STEP = 1; a PSEUDO ratio otherwise)
  int nstep;
  double eps[16];       // jump path
  double c[16];         // d jump / d theta coupling (per unit theta / theta0)
};

static double timeFactor(const Case& cs, int committedSteps)
{
  return (committedSteps == 0) ? 0.0 : cs.clock;
}

// BITS + TAN
static void bitsAndTangent(const Case& cs)
{
  Legacy L;
  L.Mc = cs.p.Mc; L.Gf = cs.p.Gf; L.Kpen = cs.p.Kpen; L.userKpen = cs.p.userKpen;
  L.shape = cs.p.shape; L.penaltyRatio = cs.p.penaltyRatio;
  L.P0 = cs.p.P0; L.alphaPM = cs.p.alphaPM; L.P0t = cs.p.P0t; L.Pbal = cs.p.Pbal; L.muBal = cs.p.muBal;
  L.implex = cs.implex; L.P_axial = cs.p.P_axial;
  L.CkappaMax = L.CkappaMax_old = 0.0; L.implexCommittedSteps = 0;
  L.computeDerived();

  Pd p = cs.p;
  computeDerived(p);
  bool calib = sameBits(p.Kpen, L.Kpen) && sameBits(p.kappa0, L.kappa0) &&
               sameBits(p.aExp, L.aExp) && sameBits(p.kappaf, L.kappaf) &&
               sameBits(pmMultiplier(p), L.getPMMultiplier());
  printf("CALIB %s %d\n", cs.name, calib ? 1 : 0);

  StateT<double> C; C.kappaMax = 0.0; C.kappaMaxOld = 0.0;
  int steps = 0;
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];
    const double tf = timeFactor(cs, steps);
    double km, stress, tangent;
    trialState(p, cs.implex, tf, eps, C, km, stress, tangent);
    L.setTrialStrain(eps, cs.clock);
    const bool bits = sameBits(stress, L.Tstress) && sameBits(tangent, L.Ttangent) &&
                      sameBits(km, L.TkappaMax);

    // TAN: the reported tangent vs the dual jump derivative (P = double, T = Dual)
    D1 e = D1::variable(eps, 0);
    StateT<D1> Cd; Cd.kappaMax = C.kappaMax; Cd.kappaMaxOld = C.kappaMaxOld;
    D1 kmd, sd, td;
    trialState(p, cs.implex, tf, e, Cd, kmd, sd, td);
    const double scale = std::fmax(std::fabs(p.Kpen), 1e-30);
    const double tanErr = std::fabs(sd.d[0] - tangent) / scale;
    const double dmg = (km > 0.0) ? 1.0 - (std::fabs(stress) / std::fmax(std::fabs(eps), 1e-300)) /
                                          (p.Kpen * pmMultiplier(p)) : 0.0;
    printf("BITS %s %d %d SOFT %d DAMAGE %.6e TAN %.3e\n", cs.name, s, bits ? 1 : 0,
           km > p.kappa0 ? 1 : 0, dmg, tanErr);

    C.kappaMaxOld = C.kappaMax; C.kappaMax = km;      // commitState
    L.commit();
    steps++;
  }
}

static double* field(Pd& p, int which)
{
  switch (which) {
    case 1: return &p.Mc;   case 2: return &p.Gf;      case 3: return &p.Kpen;
    case 4: return &p.P0;   case 5: return &p.alphaPM; case 6: return &p.P0t;
    case 7: return &p.Pbal; default: return &p.muBal;
  }
}
static D1* fieldD(ParamsT<D1>& p, int which)
{
  switch (which) {
    case 1: return &p.Mc;   case 2: return &p.Gf;      case 3: return &p.Kpen;
    case 4: return &p.P0;   case 5: return &p.alphaPM; case 6: return &p.P0t;
    case 7: return &p.Pbal; default: return &p.muBal;
  }
}

// emulates updateParameter(which, value) on the stored (derived) constants
static Pd updated(const Pd& stored, int which, double v)
{
  Pd q = stored;
  *field(q, which) = v;
  if (which == 3) q.userKpen = (v > 0.0);
  if (which <= 3) computeDerived(q);
  return q;
}

static void pathStresses(const Case& cs, const Pd& p, double dtheta, double scaleC, double* out)
{
  StateT<double> C; C.kappaMax = 0.0; C.kappaMaxOld = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    double km, stress, tangent;
    trialState(p, cs.implex, timeFactor(cs, s), cs.eps[s] + dtheta*cs.c[s]*scaleC, C, km, stress, tangent);
    out[s] = stress;
    C.kappaMaxOld = C.kappaMax; C.kappaMax = km;
  }
}

// SENS / NOHIST for parameter `which`
static void sensitivity(const Case& cs, int which, const char* pname, bool keepHistory)
{
  Pd stored = cs.p;
  computeDerived(stored);                       // the class's state after construction

  // ---- dual DDM, exactly the material's sequence (seededParams + passes) ----
  ParamsT<D1> pd;
  pd.Mc = stored.Mc; pd.Gf = stored.Gf; pd.Kpen = stored.Kpen; pd.userKpen = stored.userKpen;
  pd.shape = stored.shape; pd.penaltyRatio = stored.penaltyRatio;
  pd.P0 = stored.P0; pd.alphaPM = stored.alphaPM; pd.P0t = stored.P0t; pd.Pbal = stored.Pbal;
  pd.muBal = stored.muBal; pd.P_axial = stored.P_axial;
  pd.kappa0 = stored.kappa0; pd.aExp = stored.aExp; pd.kappaf = stored.kappaf;
  fieldD(pd, which)->d[0] = 1.0;
  if (which == 3) pd.userKpen = true;
  computeDerived(pd);
  const double theta0 = *field(stored, which);
  const double scaleC = 1.0 / std::fmax(std::fabs(theta0), 1.0);

  StateT<D1> C; C.kappaMax = 0.0; C.kappaMaxOld = 0.0;
  double ddm[16];
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];
    const double tf = timeFactor(cs, s);
    const double cn = cs.c[s]*scaleC;
    // (a) conditional stress sensitivity: jump held fixed
    D1 k1, s1, t1;
    trialState(pd, cs.implex, tf, D1(eps), C, k1, s1, t1);
    // (b) the structure would give d jump/d theta = c_n: total sensitivity
    ddm[s] = s1.d[0] + t1.v*cn;
    // (c) commitSensitivity: re-run with the total jump gradient, keep the history
    D1 e(eps); e.d[0] = cn;
    D1 k2, s2, t2;
    trialState(pd, cs.implex, tf, e, C, k2, s2, t2);
    C.kappaMaxOld = C.kappaMax; C.kappaMax = k2;
    if (!keepHistory) { C.kappaMax.d[0] = 0.0; C.kappaMaxOld.d[0] = 0.0; }
  }

  // ---- FD of the whole path (updateParameter semantics) ----
  double sp[16], sm[16];
  const double h = 1.0e-6 * std::fmax(std::fabs(theta0), 1.0);
  pathStresses(cs, updated(stored, which, theta0 + h), h, scaleC, sp);
  pathStresses(cs, updated(stored, which, theta0 - h), -h, scaleC, sm);

  double worst = 0.0, scale = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double fd = (sp[s] - sm[s]) / (2.0*h);
    scale = std::fmax(scale, std::fabs(fd));
    worst = std::fmax(worst, std::fabs(ddm[s] - fd));
  }
  printf("%s %s %s %.3e %.6e\n", keepHistory ? "SENS" : "NOHIST", cs.name, pname,
         worst / std::fmax(scale, 1e-30), scale);
}

static Pd base(double Mc, double Gf, double Kpen, int shape)
{
  Pd p;
  std::memset(&p, 0, sizeof(p));
  p.Mc = Mc; p.Gf = Gf; p.Kpen = Kpen; p.userKpen = (Kpen > 0.0); p.shape = shape;
  p.penaltyRatio = 1000.0;
  p.P0 = 0.0; p.alphaPM = 0.5; p.P0t = 0.0; p.Pbal = 0.0; p.muBal = 1.0; p.P_axial = 0.0;
  return p;
}

int main()
{
  // Mc = 100, Gf = 0.5: floor 1e4, default Kpen 1e7, kappa0 = 1e-5, exp a ~ 200,
  // linear kappaf ~ 0.01. Cyclic: pre-peak, softening, unload, reverse, reload
  // past the old peak, fully broken (linear) at the end.
  const double CYC[12] = {2e-6, 8e-6, 0.002, 0.001, 0.0, -0.003, -0.001, 0.004,
                          0.006, 0.002, -0.007, 0.012};
  const double CPL[12] = {1e-8, 2e-8, -1e-5, 3e-5, 1e-5, -2e-5, 2e-5, 1e-5,
                          -1e-5, 2e-5, 3e-5, 1e-5};
  Case cases[10];
  int n = 0;
  auto add = [&](const char* name, const Pd& p, bool implex, double clock) {
    Case& c = cases[n++];
    c.name = name; c.p = p; c.implex = implex; c.clock = clock; c.nstep = 12;
    for (int i = 0; i < 12; i++) { c.eps[i] = CYC[i]; c.c[i] = CPL[i]; }
  };
  add("exp_default", base(100.0, 0.5, -1.0, EXP), false, 1.0);
  add("lin_user",    base(100.0, 0.5, 2.0e6, LINEAR), false, 1.0);
  add("lin_clamped", base(100.0, 0.5, 5.0e3, LINEAR), false, 1.0);  // below the floor 1e4
  Pd ppm = base(100.0, 0.5, -1.0, EXP); ppm.P0 = 1000.0; ppm.alphaPM = 0.4; ppm.P_axial = -300.0;
  add("exp_pm", ppm, false, 1.0);
  Pd pb1 = base(100.0, 0.5, 3.0e6, LINEAR);
  pb1.P0 = 1000.0; pb1.P0t = 200.0; pb1.Pbal = 400.0; pb1.muBal = 1.4; pb1.P_axial = -150.0;
  add("lin_pmbal_mid", pb1, false, 1.0);
  Pd pb2 = pb1; pb2.shape = EXP; pb2.P_axial = -600.0;                 // P < -Pbal branch
  add("exp_pmbal_hi", pb2, false, 1.0);
  add("exp_implex", base(100.0, 0.5, -1.0, EXP), true, 1.0);
  Pd pix = ppm; pix.shape = LINEAR;
  add("lin_pm_implex_pseudo", pix, true, 1.5);                         // a PSEUDO step ratio

  for (int i = 0; i < n; i++) bitsAndTangent(cases[i]);

  const char* pn[9] = {"", "Mc", "Gf", "Kpen", "P0", "alphaPM", "P0t", "Pbal", "muBal"};
  for (int i = 0; i < n; i++) {
    const Pd& p = cases[i].p;
    const bool pmBal = p.Pbal > 0.0;
    for (int w = 1; w <= 8; w++) {
      if (w >= 4 && p.P0 <= 0.0) continue;            // no P-M: P-M constants are inert
      if (w == 5 && pmBal) continue;                  // alphaPM is inert under -pmBal
      if (w >= 6 && !pmBal) continue;
      sensitivity(cases[i], w, pn[w], true);
      sensitivity(cases[i], w, pn[w], false);
    }
  }
  return 0;
}
