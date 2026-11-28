/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno: 1D bond-slip tau-s UniaxialMaterial. See LadrunoBondSlip.h and
// Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md (D4).
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoBondSlip.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoBondSlip tag tau_max s1 s2 s3 tau_f alpha
//                    <-Gf Gf> <-s0 s0>
// ===========================================================================
void* OPS_LadrunoBondSlip(void)
{
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING: insufficient args\n"
           << "Want: uniaxialMaterial LadrunoBondSlip tag? tau_max? s1? s2? s3? "
           << "tau_f? alpha? <-Gf Gf?> <-s0 s0?>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoBondSlip tag\n";
    return 0;
  }

  double d[6];                 // tau_max, s1, s2, s3, tau_f, alpha
  numData = 6;
  if (OPS_GetDoubleInput(&numData, d) < 0) {
    opserr << "WARNING invalid LadrunoBondSlip backbone (tau_max s1 s2 s3 tau_f alpha)\n";
    return 0;
  }

  double Gf = 0.0;             // 0 => use s3 as given; >0 => override s3
  double s0 = 0.0;             // 0 => default 0.1*s1

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-Gf") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &Gf) < 0) {
        opserr << "WARNING LadrunoBondSlip: -Gf wants a value\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-s0") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &s0) < 0) {
        opserr << "WARNING LadrunoBondSlip: -s0 wants a value\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoBondSlip: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* m = new LadrunoBondSlip(tag, d[0], d[1], d[2], d[3],
                                            d[4], d[5], Gf, s0);
  if (m == 0) {
    opserr << "WARNING could not create LadrunoBondSlip\n";
    return 0;
  }
  return m;
}

// ===========================================================================
//  construction
// ===========================================================================
LadrunoBondSlip::LadrunoBondSlip(int tag, double tau_max_, double s1_, double s2_,
                                 double s3_, double tau_f_, double alpha_,
                                 double Gf_, double s0_)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoBondSlip),
    tau_max(tau_max_), s1(s1_), s2(s2_), s3(s3_), tau_f(tau_f_), alpha(alpha_),
    Gf(Gf_), s0(s0_), k0(0.0), s0User(s0_ > 0.0),
    Cslip(0.0), Cstress(0.0), CslipMaxAbs(0.0), Cwork(0.0),
    Tslip(0.0), Tstress(0.0), Ttangent(0.0), TslipMaxAbs(0.0), Twork(0.0),
    parameterID(0)
{
  setDerived();
  Ttangent = k0;
}

LadrunoBondSlip::LadrunoBondSlip()
  : UniaxialMaterial(0, MAT_TAG_LadrunoBondSlip),
    tau_max(0.0), s1(1.0), s2(1.0), s3(1.0), tau_f(0.0), alpha(1.0),
    Gf(0.0), s0(0.0), k0(0.0), s0User(false),
    Cslip(0.0), Cstress(0.0), CslipMaxAbs(0.0), Cwork(0.0),
    Tslip(0.0), Tstress(0.0), Ttangent(0.0), TslipMaxAbs(0.0), Twork(0.0),
    parameterID(0)
{
}

LadrunoBondSlip::~LadrunoBondSlip() {}

// fill s0 default, k0, and apply the Gf override of s3
// WP-267: LadrunoBondSlipKernel.h::setDerived (double instantiation)
void LadrunoBondSlip::setDerived(void)
{
  // WP-267: a defaulted s0 is re-derived from the current s1 (at construction s0
  // is still 0 here, so this changes no value; it keeps s0 = 0.1*s1 tracking s1 on
  // updateParameter and in the DDM pass)
  if (!s0User) s0 = 0.0;
  ladruno_bond_slip::ParamsT<double> p;
  fillKernelParams(p);
  ladruno_bond_slip::setDerived(p);
  s1 = p.s1; s2 = p.s2; s3 = p.s3; s0 = p.s0; k0 = p.k0;
}

void LadrunoBondSlip::fillKernelParams(ladruno_bond_slip::ParamsT<double>& p) const
{
  p.tau_max = tau_max; p.s1 = s1; p.s2 = s2; p.s3 = s3;
  p.tau_f = tau_f; p.alpha = alpha; p.Gf = Gf; p.s0 = s0;
  p.k0 = k0;
}

// ===========================================================================
//  backbone (a = |slip| >= 0): LadrunoBondSlipKernel.h
// ===========================================================================
double LadrunoBondSlip::envelope(double a) const
{
  ladruno_bond_slip::ParamsT<double> p;
  fillKernelParams(p);
  return ladruno_bond_slip::envelope(p, a);
}

double LadrunoBondSlip::envelopeTangent(double a) const
{
  ladruno_bond_slip::ParamsT<double> p;
  fillKernelParams(p);
  return ladruno_bond_slip::envelopeTangent(p, a);
}

// ===========================================================================
//  state update: backbone on loading past the peak |slip|; inside it, elastic
//  (k0) from the committed point, with |tau| limited to the envelope at the peak
//  (friction slider; since envelope(a) <= k0*a, the path reaches +-peak already
//  on the cap, so it rejoins the backbone without a jump). WP-267: the update
//  lives in LadrunoBondSlipKernel.h::trialState; this is its double
//  instantiation, which performs exactly the operations of the former code.
// ===========================================================================
int LadrunoBondSlip::setTrialStrain(double strain, double)
{
  Tslip = strain;
  ladruno_bond_slip::ParamsT<double> p;
  fillKernelParams(p);
  ladruno_bond_slip::StateT<double> C;
  C.slip = Cslip; C.stress = Cstress; C.slipMaxAbs = CslipMaxAbs;
  ladruno_bond_slip::trialState(p, Tslip, C, TslipMaxAbs, Tstress, Ttangent);
  // cumulative bond work per unit interface area, W = ∫ tau ds, accumulated
  // from the committed point by the trapezoidal rule (ADR 20 §10.2b — the
  // PHYSICAL bond energy the embedded element nets against the artificial
  // penalty energy). Path-dependent: committed in commitState.
  Twork = Cwork + 0.5 * (Cstress + Tstress) * (Tslip - Cslip);
  return 0;
}

// ===========================================================================
//  commit cycle
// ===========================================================================
int LadrunoBondSlip::commitState(void)
{
  Cslip = Tslip;
  Cstress = Tstress;
  CslipMaxAbs = TslipMaxAbs;
  Cwork = Twork;
  return 0;
}

int LadrunoBondSlip::revertToLastCommit(void)
{
  Tslip = Cslip;
  Tstress = Cstress;
  TslipMaxAbs = CslipMaxAbs;
  Twork = Cwork;
  Ttangent = k0;
  return 0;
}

int LadrunoBondSlip::revertToStart(void)
{
  SHV.clear();                       // WP-267: history sensitivities
  Cslip = Cstress = CslipMaxAbs = Cwork = 0.0;
  Tslip = Tstress = TslipMaxAbs = Twork = 0.0;
  Ttangent = k0;
  return 0;
}

UniaxialMaterial* LadrunoBondSlip::getCopy(void)
{
  LadrunoBondSlip* c = new LadrunoBondSlip(this->getTag(), tau_max, s1, s2, s3,
                                           tau_f, alpha, Gf, s0);
  c->Cslip = Cslip; c->Cstress = Cstress; c->CslipMaxAbs = CslipMaxAbs;
  c->Tslip = Tslip; c->Tstress = Tstress; c->TslipMaxAbs = TslipMaxAbs;
  c->Ttangent = Ttangent;
  c->Cwork = Cwork; c->Twork = Twork;
  c->parameterID = parameterID;      // WP-267: DDM
  c->s0User = s0User;                // the ctor saw the derived s0 (> 0)
  c->SHV = SHV;
  return c;
}

// ===========================================================================
//  serialization:  tag + 8 params + 4 committed state (incl. Cwork) = 13 doubles
// ===========================================================================
int LadrunoBondSlip::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(13);
  data(0)  = this->getTag();
  data(1)  = tau_max; data(2) = s1; data(3) = s2; data(4) = s3;
  data(5)  = tau_f;   data(6) = alpha; data(7) = Gf; data(8) = s0;
  data(9)  = Cslip;   data(10) = Cstress; data(11) = CslipMaxAbs;
  data(12) = Cwork;
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoBondSlip::sendSelf - failed to send\n";
    return -1;
  }
  return 0;
}

int LadrunoBondSlip::recvSelf(int commitTag, Channel& theChannel,
                              FEM_ObjectBroker& theBroker)
{
  static Vector data(13);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoBondSlip::recvSelf - failed to recv\n";
    return -1;
  }
  this->setTag((int)data(0));
  tau_max = data(1); s1 = data(2); s2 = data(3); s3 = data(4);
  tau_f = data(5); alpha = data(6); Gf = data(7); s0 = data(8);
  s0User = true;   // WP-267: the sent s0 is the derived one; keep it (not serialized)
  setDerived();
  Cslip = data(9); Cstress = data(10); CslipMaxAbs = data(11);
  Cwork = data(12);
  Tslip = Cslip; Tstress = Cstress; TslipMaxAbs = CslipMaxAbs;
  Twork = Cwork;
  Ttangent = k0;
  return 0;
}

void LadrunoBondSlip::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoBondSlip, tag: " << this->getTag() << "\n";
  s << "  tau_max: " << tau_max << "  s1: " << s1 << "  s2: " << s2
    << "  s3: " << s3 << "  tau_f: " << tau_f << "  alpha: " << alpha << "\n";
  s << "  Gf: " << Gf << "  s0: " << s0 << "  k0: " << k0 << "\n";
}

// ===========================================================================
//  responses: stress/strain/tangent (slip alias) + bondStress / slip
// ===========================================================================
Response* LadrunoBondSlip::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "bondStress") == 0)
    return new MaterialResponse(this, 1, 0.0);
  if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "slip") == 0)
    return new MaterialResponse(this, 2, 0.0);
  if (strcmp(argv[0], "tangent") == 0)
    return new MaterialResponse(this, 3, 0.0);
  // cumulative bond work per unit interface area, W = ∫ tau ds (ADR 20 §10.2b;
  // the physical bond energy, the standard UniaxialMaterial "energy" channel).
  if (strcmp(argv[0], "energy") == 0 || strcmp(argv[0], "dissipatedEnergy") == 0)
    return new MaterialResponse(this, 4, 0.0);
  return this->UniaxialMaterial::setResponse(argv, argc, s);
}

int LadrunoBondSlip::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
  case 1: matInfo.setDouble(Tstress);  return 0;
  case 2: matInfo.setDouble(Tslip);    return 0;
  case 3: matInfo.setDouble(Ttangent); return 0;
  case 4: matInfo.setDouble(Twork);    return 0;
  default: return -1;
  }
}

// ===========================================================================
//  WP-267: parameters and DDM sensitivities by dual numbers
//  (the WP-262 LadrunoUniaxialJ2 recipe)
//
//  Parameters: 1 tau_max, 2 s1, 3 s2, 4 s3, 5 tau_f, 6 alpha, 7 Gf, 8 s0.
//  updateParameter re-runs setDerived on the stored constants: a defaulted s0
//  (none given) is re-derived as 0.1*s1, as at construction, and a Gf > 0 keeps
//  overriding s3 (s3 is then inert).  After recvSelf a defaulted s0 is frozen
//  at its sent value (the flag is not serialized). The DDM pass does the same in Dual<1>, so it is the derivative of
//  what updateParameter + the run would do.
//
//  OpenSees DDM contract (computeSensitivities BEFORE commit):
//  getStressSensitivity(g) is d tau/d theta at the CURRENT trial slip held
//  fixed, given the committed state and its stored sensitivities;
//  commitSensitivity(d slip/d theta, g, n) advances them with the converged
//  slip derivative. The history is (Cslip, Cstress, CslipMaxAbs): the slider
//  reads all three. The DDM equals the derivative of the discrete response only
//  where the reported tangent is the consistent one: it is k0, not 0, on the
//  plateau and while sliding at the cap (kept so the pivot stays finite).
// ===========================================================================
int LadrunoBondSlip::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  if (strcmp(argv[0], "tau_max") == 0 || strcmp(argv[0], "tauMax") == 0) { param.setValue(tau_max); return param.addObject(1, this); }
  if (strcmp(argv[0], "s1") == 0)    { param.setValue(s1);    return param.addObject(2, this); }
  if (strcmp(argv[0], "s2") == 0)    { param.setValue(s2);    return param.addObject(3, this); }
  if (strcmp(argv[0], "s3") == 0)    { param.setValue(s3);    return param.addObject(4, this); }
  if (strcmp(argv[0], "tau_f") == 0 || strcmp(argv[0], "tauF") == 0) { param.setValue(tau_f); return param.addObject(5, this); }
  if (strcmp(argv[0], "alpha") == 0) { param.setValue(alpha); return param.addObject(6, this); }
  if (strcmp(argv[0], "Gf") == 0)    { param.setValue(Gf);    return param.addObject(7, this); }
  if (strcmp(argv[0], "s0") == 0)    { param.setValue(s0);    return param.addObject(8, this); }
  return -1;
}

int LadrunoBondSlip::updateParameter(int paramID, Information& info)
{
  switch (paramID) {
    case 1: tau_max = info.theDouble; break;
    case 2: s1      = info.theDouble; break;
    case 3: s2      = info.theDouble; break;
    case 4: s3      = info.theDouble; break;
    case 5: tau_f   = info.theDouble; break;
    case 6: alpha   = info.theDouble; break;
    case 7: Gf      = info.theDouble; break;
    case 8: s0      = info.theDouble; s0User = (s0 > 0.0); break;
    default: return -1;
  }
  setDerived();
  return 0;
}

int LadrunoBondSlip::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

void LadrunoBondSlip::seededParams(ladruno_bond_slip::ParamsT<ladruno_dual::Dual<1> >& p) const
{
  p.tau_max = tau_max; p.s1 = s1; p.s2 = s2; p.s3 = s3;
  p.tau_f = tau_f; p.alpha = alpha; p.Gf = Gf; p.s0 = s0;
  p.k0 = k0;
  // a defaulted s0 follows s1 (as in setDerived); seeding s0 itself makes it
  // user-set, as updateParameter(8) does
  if (!s0User && parameterID != 8) p.s0 = 0.0;
  switch (parameterID) {
    case 1: p.tau_max.d[0] = 1.0; break;
    case 2: p.s1.d[0]      = 1.0; break;
    case 3: p.s2.d[0]      = 1.0; break;
    case 4: p.s3.d[0]      = 1.0; break;
    case 5: p.tau_f.d[0]   = 1.0; break;
    case 6: p.alpha.d[0]   = 1.0; break;
    case 7: p.Gf.d[0]      = 1.0; break;
    case 8: p.s0.d[0]      = 1.0; break;
    default: break;                  // not this material's parameter
  }
  ladruno_bond_slip::setDerived(p);  // k0 and a Gf-driven s3 carry their derivatives
}

void LadrunoBondSlip::sensitivityPass(int gradIndex, double strainGradient,
                                      double& dStress, double& dTangent, double* dTrialState)
{
  typedef ladruno_dual::Dual<1> D1;
  ladruno_bond_slip::ParamsT<D1> p;
  seededParams(p);

  ladruno_bond_slip::StateT<D1> C;
  C.slip = Cslip; C.stress = Cstress; C.slipMaxAbs = CslipMaxAbs;
  const size_t base = (size_t)gradIndex * SHV_N;
  if (gradIndex >= 0 && SHV.size() >= base + SHV_N) {
    C.slip.d[0]       = SHV[base];
    C.stress.d[0]     = SHV[base + 1];
    C.slipMaxAbs.d[0] = SHV[base + 2];
  }

  D1 slip(Tslip);
  slip.d[0] = strainGradient;
  D1 smax, stress, tangent;
  ladruno_bond_slip::trialState(p, slip, C, smax, stress, tangent);

  dStress  = stress.d[0];
  dTangent = tangent.d[0];
  if (dTrialState) {                 // commitState: Cslip, Cstress, CslipMaxAbs = trial
    dTrialState[0] = slip.d[0];
    dTrialState[1] = stress.d[0];
    dTrialState[2] = smax.d[0];
  }
}

double LadrunoBondSlip::getStressSensitivity(int gradIndex, bool)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dS;
}

double LadrunoBondSlip::getTangentSensitivity(int gradIndex)
{
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, 0.0, dS, dT, 0);
  return dT;
}

double LadrunoBondSlip::getInitialTangentSensitivity(int)
{
  ladruno_bond_slip::ParamsT<ladruno_dual::Dual<1> > p;
  seededParams(p);
  return p.k0.d[0];                  // getInitialTangent() == k0
}

int LadrunoBondSlip::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
  if (gradIndex < 0 || numGrads <= gradIndex) return -1;
  if (SHV.size() < (size_t)numGrads * SHV_N) SHV.resize((size_t)numGrads * SHV_N, 0.0);
  double dTr[SHV_N];
  double dS = 0.0, dT = 0.0;
  this->sensitivityPass(gradIndex, strainGradient, dS, dT, dTr);
  for (int i = 0; i < SHV_N; i++) SHV[(size_t)gradIndex * SHV_N + i] = dTr[i];
  return 0;
}
