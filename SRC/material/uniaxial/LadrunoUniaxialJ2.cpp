/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: uniaxial combined isotropic (Voce+linear) + Chaboche Armstrong-
// Frederick kinematic von Mises (J2) UniaxialMaterial. See LadrunoUniaxialJ2.h
// and Ladruno_implementation/13_ladruno_uniaxial_j2_adr.md for the derivation.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoUniaxialJ2.h>
#include <LadrunoHardening.h>     // shared sig_y(pbar) -- oracle contract w/ LadrunoJ2
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

// ===========================================================================
//  OPS parser
//   uniaxialMaterial LadrunoUniaxialJ2 tag E -iso voce s0 Qinf b Hiso
//                                              -kin N C1 g1 C2 g2 ...
// ===========================================================================
void* OPS_LadrunoUniaxialJ2(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: uniaxialMaterial LadrunoUniaxialJ2 tag? E? "
           << "-iso voce s0? Qinf? b? Hiso? -kin N? C1? g1? ...\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoUniaxialJ2 tag\n";
    return 0;
  }

  double E;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &E) < 0) {
    opserr << "WARNING invalid LadrunoUniaxialJ2 E\n";
    return 0;
  }

  // defaults: perfectly plastic, no kinematic hardening
  double s0 = 0.0, Qinf = 0.0, bIso = 0.0, Hiso = 0.0;
  int nBack = 0;
  double C[LadrunoUniaxialJ2::MAXBACK], gam[LadrunoUniaxialJ2::MAXBACK];
  for (int i = 0; i < LadrunoUniaxialJ2::MAXBACK; i++) { C[i] = 0.0; gam[i] = 0.0; }

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();

    if (strcmp(flag, "-iso") == 0) {
      const char* law = OPS_GetString();
      if (strcmp(law, "voce") == 0 || strcmp(law, "Voce") == 0) {
        double d[4];
        numData = 4;
        if (OPS_GetDoubleInput(&numData, d) < 0) {
          opserr << "WARNING LadrunoUniaxialJ2: -iso voce wants s0 Qinf b Hiso\n";
          return 0;
        }
        s0 = d[0]; Qinf = d[1]; bIso = d[2]; Hiso = d[3];
      } else {
        opserr << "WARNING LadrunoUniaxialJ2: unknown -iso law '" << law
               << "' (only 'voce' in v1)\n";
        return 0;
      }
    }
    else if (strcmp(flag, "-kin") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nBack) < 0) {
        opserr << "WARNING LadrunoUniaxialJ2: -kin wants N\n";
        return 0;
      }
      if (nBack < 0 || nBack > LadrunoUniaxialJ2::MAXBACK) {
        opserr << "WARNING LadrunoUniaxialJ2: -kin N out of range [0,"
               << LadrunoUniaxialJ2::MAXBACK << "]\n";
        return 0;
      }
      for (int k = 0; k < nBack; k++) {
        double cg[2];
        numData = 2;
        if (OPS_GetDoubleInput(&numData, cg) < 0) {
          opserr << "WARNING LadrunoUniaxialJ2: -kin term " << k+1 << " wants C gamma\n";
          return 0;
        }
        C[k] = cg[0]; gam[k] = cg[1];
      }
    }
    else {
      opserr << "WARNING LadrunoUniaxialJ2: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  UniaxialMaterial* mat =
    new LadrunoUniaxialJ2(tag, E, s0, Qinf, bIso, Hiso, nBack, C, gam);
  if (mat == 0) {
    opserr << "WARNING LadrunoUniaxialJ2: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoUniaxialJ2::LadrunoUniaxialJ2(int tag, double e,
                                     double s0, double Qi, double b, double Hi,
                                     int nb, const double* C, const double* gam)
  : UniaxialMaterial(tag, MAT_TAG_LadrunoUniaxialJ2),
    E(e), sig0(s0), Qinf(Qi), bIso(b), Hiso(Hi), nBack(nb), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) {
    Ckin[k] = (k < nb) ? C[k]   : 0.0;
    gKin[k] = (k < nb) ? gam[k] : 0.0;
  }
  this->revertToStart();
}

LadrunoUniaxialJ2::LadrunoUniaxialJ2()
  : UniaxialMaterial(0, MAT_TAG_LadrunoUniaxialJ2),
    E(0.0), sig0(0.0), Qinf(0.0), bIso(0.0), Hiso(0.0), nBack(0), parameterID(0)
{
  for (int k = 0; k < MAXBACK; k++) { Ckin[k] = 0.0; gKin[k] = 0.0; }
  this->revertToStart();
}

LadrunoUniaxialJ2::~LadrunoUniaxialJ2() {}

// ===========================================================================
//  hardening (shared Ladruno law -- byte-identical with 3D LadrunoJ2)
// ===========================================================================
double LadrunoUniaxialJ2::yieldStress(double p) const
{
  return Ladruno::yieldStressVoceLinear(p, sig0, Qinf, bIso, Hiso);
}

double LadrunoUniaxialJ2::yieldSlope(double p) const
{
  return Ladruno::yieldSlopeVoceLinear(p, sig0, Qinf, bIso, Hiso);
}

// ===========================================================================
//  the return map: fixed direction s = sign(sigma_tr - X_n) + scalar Newton on
//  dp (= accumulated plastic strain increment). No M(dGamma) tensor solve.
// ===========================================================================
int LadrunoUniaxialJ2::setTrialStrain(double strain, double)
{
  if (fabs(Tstrain - strain) < DBL_EPSILON)
    return 0;

  Tstrain = strain;

  // trial (elastic predictor): freeze plastic state
  double sigTr = E*(Tstrain - CplasticStrain);
  double Xtot = 0.0;
  for (int k = 0; k < nBack; k++) Xtot += CX[k];
  double eta = sigTr - Xtot;                 // relative stress
  double absEta = fabs(eta);
  double syn = yieldStress(Cebarp);          // sig_y at committed pbar
  double f = absEta - syn;

  // unit-aware tolerance (never collapses to a fixed absolute when sig0 == 0,
  // e.g. kinematic-only); scaled to the operative stress magnitude
  double stressScale = fmax(absEta, syn);
  double tolY = 1.0e-12 * fmax(stressScale, 1.0e-300);

  if (f <= tolY) {
    // ---- elastic step: state frozen ----
    TplasticStrain = CplasticStrain;
    Tebarp = Cebarp;
    TdGamma = 0.0;
    Twp = Cwp;
    for (int k = 0; k < nBack; k++) TX[k] = CX[k];
    Tstress = sigTr;
    Ttangent = E;
    return 0;
  }

  // ---- plastic corrector ----
  int s = (eta >= 0.0) ? 1 : -1;

  // initial plastic modulus (for the Newton seed): iso slope + sum C_k
  double h0 = yieldSlope(Cebarp);
  for (int k = 0; k < nBack; k++) h0 += Ckin[k];
  double dp = f / (E + h0);                  // perfectly-plastic-style seed
  if (dp < 0.0) dp = 0.0;

  // scalar Newton:  R(dp) = |eta| - E dp
  //   - sum_k dp (C_k - gamma_k s X_k,n)/(1 + gamma_k dp) - sig_y(pbar_n + dp) = 0
  //   a_k = C_k - gamma_k s X_k,n   (constant in dp)
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
    if (dR == 0.0) {
      opserr << "WARNING LadrunoUniaxialJ2: singular local Jacobian (tag "
             << this->getTag() << ")\n";
      break;
    }
    dp -= R/dR;
    if (dp < 0.0) dp = 0.0;
  }
  if (iter >= maxIter)
    opserr << "WARNING LadrunoUniaxialJ2: local Newton did not converge (tag "
           << this->getTag() << ", |R|=" << fabs(R) << ")\n";

  // ---- state update ----
  double pbar = Cebarp + dp;
  Tebarp = pbar;
  TdGamma = dp;
  TplasticStrain = CplasticStrain + dp*s;
  for (int k = 0; k < nBack; k++)
    TX[k] = (CX[k] + Ckin[k]*dp*s) / (1.0 + gKin[k]*dp);

  // ---- stress ----
  Tstress = sigTr - E*dp*s;

  // ---- plastic work increment (fracture/LCF hook): dWp = sigma * deps^p ----
  Twp = Cwp + Tstress*(dp*s);

  // ---- consistent tangent: E_alg = E h/(E+h),  h = sig_y' + sum a_k/(1+g_k dp)^2
  double h = yieldSlope(pbar);
  for (int k = 0; k < nBack; k++) {
    double ak = Ckin[k] - gKin[k]*s*CX[k];
    double den = 1.0 + gKin[k]*dp;
    h += ak/(den*den);
  }
  Ttangent = E*h/(E + h);

  return 0;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoUniaxialJ2::commitState(void)
{
  CplasticStrain = TplasticStrain;
  Cebarp = Tebarp;
  CdGamma = TdGamma;
  Cwp = Twp;
  for (int k = 0; k < nBack; k++) CX[k] = TX[k];
  return 0;
}

int LadrunoUniaxialJ2::revertToLastCommit(void)
{
  TplasticStrain = CplasticStrain;
  Tebarp = Cebarp;
  TdGamma = CdGamma;
  Twp = Cwp;
  for (int k = 0; k < nBack; k++) TX[k] = CX[k];
  return 0;
}

int LadrunoUniaxialJ2::revertToStart(void)
{
  CplasticStrain = 0.0; Cebarp = 0.0; CdGamma = 0.0; Cwp = 0.0;
  TplasticStrain = 0.0; Tebarp = 0.0; TdGamma = 0.0; Twp = 0.0;
  for (int k = 0; k < MAXBACK; k++) { CX[k] = 0.0; TX[k] = 0.0; }
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = E;
  return 0;
}

// ===========================================================================
//  copy
// ===========================================================================
UniaxialMaterial* LadrunoUniaxialJ2::getCopy(void)
{
  LadrunoUniaxialJ2* theCopy =
    new LadrunoUniaxialJ2(this->getTag(), E, sig0, Qinf, bIso, Hiso,
                          nBack, Ckin, gKin);

  theCopy->CplasticStrain = CplasticStrain;
  theCopy->Cebarp = Cebarp;
  theCopy->CdGamma = CdGamma;
  theCopy->Cwp = Cwp;
  theCopy->TplasticStrain = TplasticStrain;
  theCopy->Tebarp = Tebarp;
  theCopy->TdGamma = TdGamma;
  theCopy->Twp = Twp;
  for (int k = 0; k < MAXBACK; k++) { theCopy->CX[k] = CX[k]; theCopy->TX[k] = TX[k]; }
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  theCopy->parameterID = parameterID;

  return theCopy;
}

// ===========================================================================
//  parallel
// ===========================================================================
int LadrunoUniaxialJ2::sendSelf(int cTag, Channel& theChannel)
{
  // tag + 5 params + nBack + 2*MAXBACK (C,gamma) + committed history
  //   (CplasticStrain, Cebarp, CdGamma, Cwp + MAXBACK CX) + Tstrain
  static Vector data(1 + 5 + 1 + 2*MAXBACK + 4 + MAXBACK + 1);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = E;
  data(c++) = sig0;
  data(c++) = Qinf;
  data(c++) = bIso;
  data(c++) = Hiso;
  data(c++) = nBack;
  for (int k = 0; k < MAXBACK; k++) data(c++) = Ckin[k];
  for (int k = 0; k < MAXBACK; k++) data(c++) = gKin[k];
  data(c++) = CplasticStrain;
  data(c++) = Cebarp;
  data(c++) = CdGamma;
  data(c++) = Cwp;
  for (int k = 0; k < MAXBACK; k++) data(c++) = CX[k];
  data(c++) = Tstrain;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoUniaxialJ2::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoUniaxialJ2::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(1 + 5 + 1 + 2*MAXBACK + 4 + MAXBACK + 1);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "LadrunoUniaxialJ2::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  E    = data(c++);
  sig0 = data(c++);
  Qinf = data(c++);
  bIso = data(c++);
  Hiso = data(c++);
  nBack = (int)data(c++);
  for (int k = 0; k < MAXBACK; k++) Ckin[k] = data(c++);
  for (int k = 0; k < MAXBACK; k++) gKin[k] = data(c++);
  CplasticStrain = data(c++);
  Cebarp = data(c++);
  CdGamma = data(c++);
  Cwp = data(c++);
  for (int k = 0; k < MAXBACK; k++) CX[k] = data(c++);
  Tstrain = data(c++);

  // sync trial to committed
  this->revertToLastCommit();
  Tstress = E*(Tstrain - CplasticStrain);   // recovered elastic predictor stress
  Ttangent = E;
  return 0;
}

// ===========================================================================
//  parameters / print
// ===========================================================================
int LadrunoUniaxialJ2::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1) return -1;
  // Scalar elastic / isotropic parameters (FE sensitivity + updateMaterialStage).
  // NOTE: the dispatch is intentionally a flat strcmp ladder so future layers
  // (e.g. a local-buckling slenderness key, ADR 13 sec 6.2) can be added here
  // without reshaping it -- setTrialStrain stays scalar.
  if (strcmp(argv[0], "E") == 0)                                  { param.setValue(E);    return param.addObject(1, this); }
  if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "sig0") == 0 ||
      strcmp(argv[0], "fy") == 0 || strcmp(argv[0], "Fy") == 0)   { param.setValue(sig0); return param.addObject(2, this); }
  if (strcmp(argv[0], "Qinf") == 0)                               { param.setValue(Qinf); return param.addObject(3, this); }
  if (strcmp(argv[0], "b") == 0)                                  { param.setValue(bIso); return param.addObject(4, this); }
  if (strcmp(argv[0], "Hiso") == 0 || strcmp(argv[0], "H_iso") == 0) { param.setValue(Hiso); return param.addObject(5, this); }
  return -1;
}

int LadrunoUniaxialJ2::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 1: E    = info.theDouble; return 0;
    case 2: sig0 = info.theDouble; return 0;
    case 3: Qinf = info.theDouble; return 0;
    case 4: bIso = info.theDouble; return 0;
    case 5: Hiso = info.theDouble; return 0;
    default: return -1;
  }
}

int LadrunoUniaxialJ2::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}

void LadrunoUniaxialJ2::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LadrunoUniaxialJ2\", ";
    s << "\"E\": " << E << ", ";
    s << "\"sig0\": " << sig0 << ", \"Qinf\": " << Qinf << ", \"b\": " << bIso
      << ", \"Hiso\": " << Hiso << ", ";
    s << "\"nBack\": " << nBack << "}";
    return;
  }
  s << "LadrunoUniaxialJ2 (uniaxial combined iso + Chaboche AF kinematic J2)" << endln;
  s << "  tag    : " << this->getTag() << endln;
  s << "  E      : " << E << endln;
  s << "  iso    : sig0=" << sig0 << " Qinf=" << Qinf << " b=" << bIso << " Hiso=" << Hiso << endln;
  s << "  nBack  : " << nBack << endln;
  for (int k = 0; k < nBack; k++)
    s << "    term " << k+1 << ": C=" << Ckin[k] << " gamma=" << gKin[k] << endln;
}

// ===========================================================================
//  recordable responses (stress/strain/tangent + internal state hooks)
// ===========================================================================
Response* LadrunoUniaxialJ2::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return UniaxialMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "backStress") == 0 || strcmp(a, "backstress") == 0 || strcmp(a, "alpha") == 0)
    return new MaterialResponse(this, 100, 0.0);
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 101, 0.0);
  if (strcmp(a, "equivalentPlasticStrain") == 0 || strcmp(a, "plasticStrainEq") == 0 ||
      strcmp(a, "ebarP") == 0)
    return new MaterialResponse(this, 102, 0.0);
  if (strcmp(a, "plasticWork") == 0 || strcmp(a, "plasticEnergy") == 0)
    return new MaterialResponse(this, 103, 0.0);

  return UniaxialMaterial::setResponse(argv, argc, s);
}

int LadrunoUniaxialJ2::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 100: {                              // total backstress X = sum_k X_k
      double tot = 0.0;
      for (int k = 0; k < nBack; k++) tot += TX[k];
      matInfo.theDouble = tot;
      return 0;
    }
    case 101: matInfo.theDouble = TplasticStrain; return 0;
    case 102: matInfo.theDouble = Tebarp;         return 0;
    case 103: matInfo.theDouble = Twp;            return 0;
    default:
      return UniaxialMaterial::getResponse(responseID, matInfo);
  }
}
