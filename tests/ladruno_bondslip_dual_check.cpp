// WP-267: the extracted LadrunoBondSlip kernel (SRC/material/uniaxial/
// LadrunoBondSlipKernel.h) checked four ways, OpenSees-free:
//
//  BITS   the kernel in double vs a VERBATIM transcription of the pre-WP-267
//         LadrunoBondSlip (setDerived, envelope, envelopeTangent,
//         setTrialStrain, commitState): the derived constants, stress, tangent
//         and trial peak slip must agree bit for bit on every step.
//  TAN    the reported tangent vs a Dual<1> seeded on the slip, per branch.
//         It must agree on the envelope and on elastic unload/reload. On the
//         plateau, the residual branch and while sliding at the cap the class
//         deliberately reports k0 (or -1e-12 k0), not the consistent 0: the
//         check prints the branch so the test can tell the two apart.
//  SENS   the DDM history machinery: per parameter, the bond-stress
//         sensitivity propagated step by step exactly as the material does it
//         (a conditional pass, then commitSensitivity with the total slip
//         gradient, whose stress derivative is the total one) vs a central FD
//         of the WHOLE path re-run at theta +/- h with updateParameter's
//         re-derivation. The slip path depends on theta
//         (s_n = s0_n + (theta - theta0) c_n).
//  NOHIST the same DDM with the history sensitivities (Cslip, Cstress,
//         CslipMaxAbs) dropped every step: it must miss.
//
// Build: g++ -O2 -std=c++17 -I SRC/material/nD -I SRC/material/uniaxial
//            tests/ladruno_bondslip_dual_check.cpp -o lbsd

#include <LadrunoBondSlipKernel.h>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace ladruno_bond_slip;
typedef ladruno_dual::Dual<1> D1;
typedef ParamsT<double> Pd;

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// --------------------------------------------------------------------------
// Verbatim transcription of LadrunoBondSlip as it stood on origin/ladruno @
// f183b9cd1 (members -> struct fields). DO NOT "tidy" it: it is the
// bit-identity reference.
// --------------------------------------------------------------------------
struct Legacy {
  double tau_max, s1, s2, s3, tau_f, alpha, Gf, s0;
  double k0;
  double Cslip, Cstress, CslipMaxAbs;
  double Tslip, Tstress, Ttangent, TslipMaxAbs;

  void setDerived(void)
  {
    if (s1 <= 0.0) s1 = 1.0e-12;
    if (s0 <= 0.0) s0 = 0.1 * s1;
    if (s0 > s1)   s0 = s1;
    double tau_s0 = tau_max * pow(s0 / s1, alpha);
    k0 = (s0 > 0.0) ? tau_s0 / s0 : tau_max / s1;
    if (Gf > 0.0 && tau_max > tau_f) {
      s3 = s2 + 2.0 * Gf / (tau_max - tau_f);
    }
    if (s2 < s1) s2 = s1;
    if (s3 < s2) s3 = s2;
  }

  double envelope(double a) const
  {
    if (a <= s0)  return k0 * a;
    if (a <= s1)  return tau_max * pow(a / s1, alpha);
    if (a <= s2)  return tau_max;
    if (a <= s3)  return tau_max - (tau_max - tau_f) * (a - s2) / (s3 - s2);
    return tau_f;
  }

  double envelopeTangent(double a) const
  {
    if (a <= s0)  return k0;
    if (a <= s1)  return tau_max * alpha / s1 * pow(a / s1, alpha - 1.0);
    if (a <= s2)  return 0.0;
    if (a <= s3)  return -(tau_max - tau_f) / (s3 - s2);
    return 0.0;
  }

  int setTrialStrain(double strain)
  {
    Tslip = strain;
    double a = fabs(Tslip);
    double sgn = (Tslip >= 0.0) ? 1.0 : -1.0;

    if (a >= CslipMaxAbs) {
      TslipMaxAbs = a;
      Tstress  = sgn * envelope(a);
      Ttangent = envelopeTangent(a);
      if (Ttangent == 0.0) Ttangent = (a <= s2) ? k0 : -1.0e-12 * k0;
    }
    else {
      TslipMaxAbs = CslipMaxAbs;
      const double tauCap = envelope(CslipMaxAbs);
      Tstress = Cstress + k0 * (Tslip - Cslip);
      if (Tstress > tauCap)       Tstress = tauCap;
      else if (Tstress < -tauCap) Tstress = -tauCap;
      Ttangent = k0;
    }
    return 0;
  }

  void commit(void) { Cslip = Tslip; Cstress = Tstress; CslipMaxAbs = TslipMaxAbs; }
};

// --------------------------------------------------------------------------
struct Case {
  const char* name;
  Pd p;                 // user inputs; setDerived'd where used
  int nstep;
  double eps[16];       // slip path
  double c[16];         // d slip / d theta coupling (per unit theta / theta0)
};

static Pd derived(const Pd& in) { Pd p = in; setDerived(p); return p; }

// BITS + TAN
static void bitsAndTangent(const Case& cs)
{
  Legacy L;
  L.tau_max = cs.p.tau_max; L.s1 = cs.p.s1; L.s2 = cs.p.s2; L.s3 = cs.p.s3;
  L.tau_f = cs.p.tau_f; L.alpha = cs.p.alpha; L.Gf = cs.p.Gf; L.s0 = cs.p.s0; L.k0 = 0.0;
  L.Cslip = L.Cstress = L.CslipMaxAbs = 0.0;
  L.setDerived();

  const Pd p = derived(cs.p);
  const bool calib = sameBits(p.s1, L.s1) && sameBits(p.s2, L.s2) && sameBits(p.s3, L.s3) &&
                     sameBits(p.s0, L.s0) && sameBits(p.k0, L.k0);
  printf("CALIB %s %d\n", cs.name, calib ? 1 : 0);

  StateT<double> C; C.slip = 0.0; C.stress = 0.0; C.slipMaxAbs = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];
    double smax, stress, tangent;
    const int br = trialState(p, eps, C, smax, stress, tangent);
    L.setTrialStrain(eps);
    const bool bits = sameBits(stress, L.Tstress) && sameBits(tangent, L.Ttangent) &&
                      sameBits(smax, L.TslipMaxAbs);

    D1 e = D1::variable(eps, 0);
    StateT<D1> Cd; Cd.slip = C.slip; Cd.stress = C.stress; Cd.slipMaxAbs = C.slipMaxAbs;
    D1 smd, sd, td;
    trialState(p, e, Cd, smd, sd, td);
    const double tanErr = std::fabs(sd.d[0] - tangent) / p.k0;
    printf("BITS %s %d %d BRANCH %d TAN %.3e\n", cs.name, s, bits ? 1 : 0, br, tanErr);

    C.slip = eps; C.stress = stress; C.slipMaxAbs = smax;    // commitState
    L.commit();
  }
}

static double* field(Pd& p, int which)
{
  switch (which) {
    case 1: return &p.tau_max; case 2: return &p.s1;    case 3: return &p.s2;
    case 4: return &p.s3;      case 5: return &p.tau_f; case 6: return &p.alpha;
    case 7: return &p.Gf;      default: return &p.s0;
  }
}
static D1* fieldD(ParamsT<D1>& p, int which)
{
  switch (which) {
    case 1: return &p.tau_max; case 2: return &p.s1;    case 3: return &p.s2;
    case 4: return &p.s3;      case 5: return &p.tau_f; case 6: return &p.alpha;
    case 7: return &p.Gf;      default: return &p.s0;
  }
}

static void pathStresses(const Case& cs, const Pd& p, double dtheta, double scaleC, double* out)
{
  StateT<double> C; C.slip = 0.0; C.stress = 0.0; C.slipMaxAbs = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s] + dtheta*cs.c[s]*scaleC;
    double smax, stress, tangent;
    trialState(p, eps, C, smax, stress, tangent);
    out[s] = stress;
    C.slip = eps; C.stress = stress; C.slipMaxAbs = smax;
  }
}

static void sensitivity(const Case& cs, int which, const char* pname, bool keepHistory)
{
  const Pd stored = derived(cs.p);              // the class's state after construction

  // ---- dual DDM, exactly the material's sequence (seededParams + passes) ----
  ParamsT<D1> pd;
  pd.tau_max = stored.tau_max; pd.s1 = stored.s1; pd.s2 = stored.s2; pd.s3 = stored.s3;
  pd.tau_f = stored.tau_f; pd.alpha = stored.alpha; pd.Gf = stored.Gf; pd.s0 = stored.s0;
  pd.k0 = stored.k0;
  fieldD(pd, which)->d[0] = 1.0;
  setDerived(pd);
  Pd chk = stored;
  const double theta0 = *field(chk, which);
  const double scaleC = 1.0 / std::fmax(std::fabs(theta0), 1.0);

  StateT<D1> C; C.slip = 0.0; C.stress = 0.0; C.slipMaxAbs = 0.0;
  double ddm[16];
  for (int s = 0; s < cs.nstep; s++) {
    const double eps = cs.eps[s];
    const double cn = cs.c[s]*scaleC;
    // (a) conditional pass (getStressSensitivity): slip held fixed
    D1 m1, s1, t1;
    trialState(pd, D1(eps), C, m1, s1, t1);
    // (b) commitSensitivity: the total slip gradient; its stress derivative is
    //     the total d tau/d theta (= conditional + consistent tangent * c_n)
    D1 e(eps); e.d[0] = cn;
    D1 m2, s2, t2;
    trialState(pd, e, C, m2, s2, t2);
    ddm[s] = s2.d[0];
    C.slip = e; C.stress = s2; C.slipMaxAbs = m2;
    if (!keepHistory) { C.slip.d[0] = 0.0; C.stress.d[0] = 0.0; C.slipMaxAbs.d[0] = 0.0; }
    (void)s1;
  }

  // ---- FD of the whole path (updateParameter: set, then setDerived) ----
  double sp[16], sm[16];
  const double h = 1.0e-6 * std::fmax(std::fabs(theta0), 1.0);
  Pd pp = stored, pm = stored;
  *field(pp, which) = theta0 + h; setDerived(pp);
  *field(pm, which) = theta0 - h; setDerived(pm);
  pathStresses(cs, pp, h, scaleC, sp);
  pathStresses(cs, pm, -h, scaleC, sm);

  double worst = 0.0, scale = 0.0;
  for (int s = 0; s < cs.nstep; s++) {
    const double fd = (sp[s] - sm[s]) / (2.0*h);
    scale = std::fmax(scale, std::fabs(fd));
    worst = std::fmax(worst, std::fabs(ddm[s] - fd));
  }
  printf("%s %s %s %.3e %.6e\n", keepHistory ? "SENS" : "NOHIST", cs.name, pname,
         worst / std::fmax(scale, 1e-30), scale);
}

static Pd base(double tau_max, double s1, double s2, double s3, double tau_f,
               double alpha, double Gf, double s0)
{
  Pd p;
  std::memset(&p, 0, sizeof(p));
  p.tau_max = tau_max; p.s1 = s1; p.s2 = s2; p.s3 = s3; p.tau_f = tau_f;
  p.alpha = alpha; p.Gf = Gf; p.s0 = s0; p.k0 = 0.0;
  return p;
}

int main()
{
  // slip path (mm): linear segment, power law, plateau, unload, slide, reload,
  // softening, reversal past the old peak, residual, unload
  const double CYC[14] = {0.04, 0.6, 1.5, 1.45, 1.0, 0.9, 1.2, 2.5,
                          2.9, -1.0, -3.5, -5.0, -4.0, 4.5};
  const double CPL[14] = {1e-4, -2e-4, 1e-4, 3e-4, -1e-4, 2e-4, -2e-4, 1e-4,
                          3e-4, -1e-4, 2e-4, -3e-4, 1e-4, 2e-4};
  Case cases[4];
  int n = 0;
  auto add = [&](const char* name, const Pd& p) {
    Case& c = cases[n++];
    c.name = name; c.p = p; c.nstep = 14;
    for (int i = 0; i < 14; i++) { c.eps[i] = CYC[i]; c.c[i] = CPL[i]; }
  };
  add("mc2010",         base(10.0, 1.0, 2.0, 4.0, 3.0, 0.4, 0.0, 0.0));   // s0 default 0.1 s1
  add("gf_override",    base(10.0, 1.0, 2.0, 4.0, 3.0, 0.4, 8.0, 0.0));   // s3 = s2 + 2Gf/(tmax-tf)
  add("short_plateau",  base(12.0, 1.0, 1.2, 3.0, 2.0, 0.3, 0.0, 0.2));   // user s0
  add("steep_alpha",    base(8.0, 0.8, 1.6, 3.5, 1.5, 0.7, 6.0, 0.05));

  for (int i = 0; i < n; i++) bitsAndTangent(cases[i]);

  const char* pn[9] = {"", "tau_max", "s1", "s2", "s3", "tau_f", "alpha", "Gf", "s0"};
  for (int i = 0; i < n; i++) {
    const bool gf = cases[i].p.Gf > 0.0;
    for (int w = 1; w <= 8; w++) {
      if (w == 4 && gf) continue;        // Gf overrides s3: inert
      if (w == 7 && !gf) continue;       // Gf = 0: inert (and a kink at 0)
      sensitivity(cases[i], w, pn[w], true);
      sensitivity(cases[i], w, pn[w], false);
    }
  }
  return 0;
}
