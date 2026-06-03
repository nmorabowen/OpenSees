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
    Gf(Gf_), s0(s0_), k0(0.0),
    Cslip(0.0), Cstress(0.0), CslipMaxAbs(0.0), Cwork(0.0),
    Tslip(0.0), Tstress(0.0), Ttangent(0.0), TslipMaxAbs(0.0), Twork(0.0)
{
  setDerived();
  Ttangent = k0;
}

LadrunoBondSlip::LadrunoBondSlip()
  : UniaxialMaterial(0, MAT_TAG_LadrunoBondSlip),
    tau_max(0.0), s1(1.0), s2(1.0), s3(1.0), tau_f(0.0), alpha(1.0),
    Gf(0.0), s0(0.0), k0(0.0),
    Cslip(0.0), Cstress(0.0), CslipMaxAbs(0.0), Cwork(0.0),
    Tslip(0.0), Tstress(0.0), Ttangent(0.0), TslipMaxAbs(0.0), Twork(0.0)
{
}

LadrunoBondSlip::~LadrunoBondSlip() {}

// fill s0 default, k0, and apply the Gf override of s3
void LadrunoBondSlip::setDerived(void)
{
  if (s1 <= 0.0) s1 = 1.0e-12;
  if (s0 <= 0.0) s0 = 0.1 * s1;
  if (s0 > s1)   s0 = s1;
  // initial stiffness = secant to the power-law value at s0 (D4.1 linear segment)
  double tau_s0 = tau_max * pow(s0 / s1, alpha);
  k0 = (s0 > 0.0) ? tau_s0 / s0 : tau_max / s1;
  // D4.3: Gf overrides s3 so the softening triangle dissipates Gf per unit area
  if (Gf > 0.0 && tau_max > tau_f) {
    s3 = s2 + 2.0 * Gf / (tau_max - tau_f);
  }
  if (s2 < s1) s2 = s1;
  if (s3 < s2) s3 = s2;
}

// ===========================================================================
//  backbone (a = |slip| >= 0)
// ===========================================================================
double LadrunoBondSlip::envelope(double a) const
{
  if (a <= s0)  return k0 * a;
  if (a <= s1)  return tau_max * pow(a / s1, alpha);
  if (a <= s2)  return tau_max;
  if (a <= s3)  return tau_max - (tau_max - tau_f) * (a - s2) / (s3 - s2);
  return tau_f;
}

double LadrunoBondSlip::envelopeTangent(double a) const
{
  if (a <= s0)  return k0;
  if (a <= s1)  return tau_max * alpha / s1 * pow(a / s1, alpha - 1.0);
  if (a <= s2)  return 0.0;
  if (a <= s3)  return -(tau_max - tau_f) / (s3 - s2);
  return 0.0;
}

// ===========================================================================
//  state update: backbone on loading past the peak |slip|, elastic (k0) inside
// ===========================================================================
int LadrunoBondSlip::setTrialStrain(double strain, double)
{
  Tslip = strain;
  double a = fabs(Tslip);
  double sgn = (Tslip >= 0.0) ? 1.0 : -1.0;

  if (a >= CslipMaxAbs) {              // virgin loading on the envelope
    TslipMaxAbs = a;
    Tstress  = sgn * envelope(a);
    Ttangent = envelopeTangent(a);    // (>=0 ascending, <0 softening, 0 plateau)
    if (Ttangent == 0.0) Ttangent = (a <= s2) ? k0 : -1.0e-12 * k0;  // avoid a hard 0 pivot
  }
  else {                              // unload/reload elastically toward/from the peak
    TslipMaxAbs = CslipMaxAbs;
    double sgnMax = (CslipMaxAbs == 0.0) ? sgn
                    : ((Cslip >= 0.0 || (Cslip == 0.0 && Tslip >= 0.0)) ? 1.0 : -1.0);
    // unload line through the committed-peak envelope point with slope k0
    double tauPeak = sgnMax * envelope(CslipMaxAbs);
    Tstress  = tauPeak - k0 * (sgnMax * CslipMaxAbs - Tslip);
    Ttangent = k0;
  }
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
