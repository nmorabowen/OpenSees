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

// Ladruno RC plastic-damage nDMaterial — see LadrunoRCConcrete.h. classTag 33015.

#include "LadrunoRCConcrete.h"
#include <elementAPI.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace ladruno_rc_kernel;

// ===========================================================================
//  parser:  nDMaterial LadrunoRCConcrete $tag $E $nu
//             -Ce {..} -Cs {..} <-Cd {..}>  (compression backbone: strain/nom-stress/damage)
//             -Te {..} -Ts {..} <-Td {..}>  (tension backbone)
//             <-Kc $Kc> <-beta> <-betaFloor $f> <-lublinerReduced>
//             <-tensStiff {vc|cm}> <-tensStiffC $c> <-tensStiffAlpha $a>  (Phase 3)
//             <-rho $rho> <-secant | -numericalTangent>
// ===========================================================================
void* OPS_LadrunoRCConcrete(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "nDMaterial LadrunoRCConcrete error: too few args.\n"
           << "  want: tag E nu -Ce {..} -Cs {..} <-Cd {..}> -Te {..} -Ts {..} <-Td {..}>"
           << " <-Kc Kc> <-beta> <-betaFloor f> <-lublinerReduced> <-rho rho> <-secant|-numericalTangent>\n";
    return 0;
  }

  int numData = 1;
  int tag = 0;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "nDMaterial LadrunoRCConcrete error: invalid tag.\n";
    return 0;
  }
  double E = 0.0, nu = 0.0;
  if (OPS_GetDoubleInput(&numData, &E) < 0)  { opserr << "LadrunoRCConcrete: invalid E.\n";  return 0; }
  if (OPS_GetDoubleInput(&numData, &nu) < 0) { opserr << "LadrunoRCConcrete: invalid nu.\n"; return 0; }

  std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
  double Kc = 2.0 / 3.0, betaFloor = 0.1, rho = 0.0;
  bool betaOn = false, lubRed = false;
  int  tanMode = 0;
  // Phase 2a/2b: aggregate-interlock shear retention (default off)
  bool   interlockOn = false, interlockCyclic = false;
  double aggSize = 16.0, crackStrain = 0.0, crackSpacing = 0.0, lch = 0.0, betaSrMin = 0.01;
  int    shearRetMode = 0;       // 0 mcft (default) | 1 const | 2 dsfm | 3 rots
  double shearRetFactor = 0.4;   // const-mode retention mu in (0,1]
  // Phase 2b.2b: second orthogonal crack (X-cracking) + slip-driven interlock wear
  bool   xcrackOn = false;
  double degKappa = 0.5, degSlipRef = 0.01, degMin = 0.1;
  // Phase 4: IMPL-EX (implicit-explicit) robustness for cyclic softening
  bool   implexOn = false, implexCtrl = false;
  double implexAlpha = 1.0, implexErrTol = 0.05, implexTimeRedLim = 0.01;
  // Phase 3: tension stiffening (default off => baseline-identical)
  int    tensStiffMode = 0;       // 0 off | 1 vc (Bentz) | 2 cm (Collins-Mitchell)
  double tensStiffC = 500.0, tensStiffAlpha = 1.0;

  auto readList = [](std::vector<double>& v) {
    v.clear();
    int nd = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      double item;
      int oldn = OPS_GetNumRemainingInputArgs();
      if (OPS_GetDoubleInput(&nd, &item) < 0) {
        int newn = OPS_GetNumRemainingInputArgs();
        if (newn < oldn) OPS_ResetCurrentInputArg(-1);   // un-consume the flag/string
        break;
      }
      v.push_back(item);
    }
  };

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if      (strcmp(opt, "-Ce") == 0) readList(Ce);
    else if (strcmp(opt, "-Cs") == 0) readList(Cs);
    else if (strcmp(opt, "-Cd") == 0) readList(Cd);
    else if (strcmp(opt, "-Te") == 0) readList(Te);
    else if (strcmp(opt, "-Ts") == 0) readList(Ts);
    else if (strcmp(opt, "-Td") == 0) readList(Td);
    else if (strcmp(opt, "-Kc") == 0)        { int nd = 1; OPS_GetDoubleInput(&nd, &Kc); }
    else if (strcmp(opt, "-betaFloor") == 0) { int nd = 1; OPS_GetDoubleInput(&nd, &betaFloor); }
    else if (strcmp(opt, "-rho") == 0)       { int nd = 1; OPS_GetDoubleInput(&nd, &rho); }
    else if (strcmp(opt, "-beta") == 0)              betaOn = true;
    else if (strcmp(opt, "-lublinerReduced") == 0)   lubRed = true;
    else if (strcmp(opt, "-secant") == 0)            tanMode = 1;
    else if (strcmp(opt, "-numericalTangent") == 0)  tanMode = 2;
    else if (strcmp(opt, "-interlock") == 0)         interlockOn = true;
    else if (strcmp(opt, "-cyclic") == 0)            interlockCyclic = true;
    else if (strcmp(opt, "-agg") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &aggSize) < 0)      { opserr << "LadrunoRCConcrete: -agg needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-crackStrain") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackStrain) < 0)  { opserr << "LadrunoRCConcrete: -crackStrain needs a value.\n";  return 0; } }
    else if (strcmp(opt, "-crackSpacing") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackSpacing) < 0) { opserr << "LadrunoRCConcrete: -crackSpacing needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-lch") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &lch) < 0)          { opserr << "LadrunoRCConcrete: -lch needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-betaSrMin") == 0)    { int nd = 1; if (OPS_GetDoubleInput(&nd, &betaSrMin) < 0)    { opserr << "LadrunoRCConcrete: -betaSrMin needs a value.\n";    return 0; } }
    else if (strcmp(opt, "-xcrack") == 0)            xcrackOn = true;
    else if (strcmp(opt, "-degKappa") == 0)     { int nd = 1; if (OPS_GetDoubleInput(&nd, &degKappa) < 0)     { opserr << "LadrunoRCConcrete: -degKappa needs a value.\n";     return 0; } }
    else if (strcmp(opt, "-degSlipRef") == 0)   { int nd = 1; if (OPS_GetDoubleInput(&nd, &degSlipRef) < 0)   { opserr << "LadrunoRCConcrete: -degSlipRef needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-degMin") == 0)       { int nd = 1; if (OPS_GetDoubleInput(&nd, &degMin) < 0)       { opserr << "LadrunoRCConcrete: -degMin needs a value.\n";       return 0; } }
    else if (strcmp(opt, "-implex") == 0)            implexOn = true;
    else if (strcmp(opt, "-implexAlpha") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &implexAlpha) < 0)   { opserr << "LadrunoRCConcrete: -implexAlpha needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-implexControl") == 0) {
      implexCtrl = true; int nd = 1;
      if (OPS_GetDoubleInput(&nd, &implexErrTol) < 0 || OPS_GetDoubleInput(&nd, &implexTimeRedLim) < 0) {
        opserr << "LadrunoRCConcrete: -implexControl needs $errTol $timeReductionLimit.\n"; return 0; }
    }
    // Phase 2b.2c: -shearRetention {mcft|const|dsfm|rots} crack-shear retention curve.
    else if (strcmp(opt, "-shearRetention") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "mcft")  == 0) shearRetMode = 0;
      else if (m && strcmp(m, "const") == 0) shearRetMode = 1;
      else if (m && strcmp(m, "dsfm")  == 0) shearRetMode = 2;
      else if (m && strcmp(m, "rots")  == 0) shearRetMode = 3;
      else { opserr << "LadrunoRCConcrete: -shearRetention needs {mcft|const|dsfm|rots}.\n"; return 0; }
    }
    else if (strcmp(opt, "-shearRetFactor") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &shearRetFactor) < 0) { opserr << "LadrunoRCConcrete: -shearRetFactor needs a value.\n"; return 0; } }
    // Phase 3: -tensStiff {vc|cm} tension-stiffening mode (+ -tensStiffC c, -tensStiffAlpha a).
    else if (strcmp(opt, "-tensStiff") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "vc") == 0) tensStiffMode = 1;
      else if (m && strcmp(m, "cm") == 0) tensStiffMode = 2;
      else { opserr << "LadrunoRCConcrete: -tensStiff needs {vc|cm}.\n"; return 0; }
    }
    else if (strcmp(opt, "-tensStiffC") == 0)     { int nd = 1; if (OPS_GetDoubleInput(&nd, &tensStiffC) < 0)     { opserr << "LadrunoRCConcrete: -tensStiffC needs a value.\n";     return 0; } }
    else if (strcmp(opt, "-tensStiffAlpha") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &tensStiffAlpha) < 0) { opserr << "LadrunoRCConcrete: -tensStiffAlpha needs a value.\n"; return 0; } }
    // unknown tokens are ignored (forward-compat)
  }

  // -cyclic implies -interlock; -xcrack implies -cyclic (the X-crack/wear law lives inside
  // the cyclic friction-slip block, which lives inside the interlock block).
  if (xcrackOn) interlockCyclic = true;
  // a non-default shear-retention curve is an interlock feature; const/dsfm parametrize the
  // CYCLIC slip stiffness, so they imply -interlock -cyclic. rots only needs -interlock
  // (it disables the fixed-crack shear). mcft (0) implies nothing (it is the default).
  if (shearRetMode == 1 || shearRetMode == 2) interlockCyclic = true;
  if (shearRetMode != 0) interlockOn = true;
  if (interlockCyclic && !interlockOn) interlockOn = true;

  if (Ce.size() < 2 || Cs.size() != Ce.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Ce and -Cs of equal length (>=2).\n";
    return 0;
  }
  if (Te.size() < 2 || Ts.size() != Te.size()) {
    opserr << "nDMaterial LadrunoRCConcrete error: need -Te and -Ts of equal length (>=2).\n";
    return 0;
  }

  Params P;
  P.E = E; P.nu = nu; P.Kc = Kc; P.betaFloor = betaFloor;
  P.cdf = 0.0; P.eta = 0.0;
  P.betaOn = betaOn; P.lublinerTCReduced = lubRed; P.tangentMode = tanMode;
  P.interlockOn = interlockOn; P.interlockCyclic = interlockCyclic;
  P.shearRetMode = shearRetMode; P.shearRetFactor = shearRetFactor;
  P.aggSize = aggSize; P.crackStrain = crackStrain; P.crackSpacing = crackSpacing;
  P.lch = lch; P.betaSrMin = betaSrMin; P.sqrtFc = 0.0;
  P.xcrackOn = xcrackOn; P.degKappa = degKappa; P.degSlipRef = degSlipRef; P.degMin = degMin;
  P.implex = implexOn; P.implexAlpha = implexAlpha; P.implexControl = implexCtrl;
  P.implexErrTol = implexErrTol; P.implexTimeRedLim = implexTimeRedLim;
  P.tensStiffMode = tensStiffMode; P.tensStiffC = tensStiffC;
  P.tensStiffAlpha = tensStiffAlpha; P.ftPeak = 0.0;   // ftPeak set by setupParams

  // build backbones via the faithful ASDConcrete3D HardeningLaw c-tor + adjust()
  // (elastic-consistent q). -Cd/-Td are optional -> pad with zeros to match length.
  std::vector<double> Cdf(Ce.size(), 0.0), Tdf(Te.size(), 0.0);
  for (size_t i = 0; i < Cd.size() && i < Cdf.size(); ++i) Cdf[i] = Cd[i];
  for (size_t i = 0; i < Td.size() && i < Tdf.size(); ++i) Tdf[i] = Td[i];
  buildBackbone(P.hc, E, Ce.data(), Cs.data(), Cdf.data(), (int)Ce.size());
  buildBackbone(P.ht, E, Te.data(), Ts.data(), Tdf.data(), (int)Te.size());
  setupParams(P);

  return new LadrunoRCConcrete(tag, P, rho, LadrunoRCConcrete::DIM_3D);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoRCConcrete::LadrunoRCConcrete()
  : NDMaterial(0, ND_TAG_LadrunoRCConcrete),
    rho(0.0), dim(DIM_3D), ncomp(6), condense(false), cEps33(0.0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  // safe defaults until recvSelf populates P
  P.E = 1.0; P.nu = 0.0; P.Kc = 2.0/3.0; P.fcft_ratio = 5.0; P.betaFloor = 0.1;
  P.cdf = 0.0; P.eta = 0.0; P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
  P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4; P.aggSize = 16.0;
  P.crackStrain = 0.0; P.crackSpacing = 0.0; P.lch = 0.0; P.betaSrMin = 0.01; P.sqrtFc = 0.0;
  P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
  P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
  P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
  P.tensStiffMode = 0; P.tensStiffC = 500.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
  P.ht.n = 0; P.hc.n = 0;
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::LadrunoRCConcrete(int tag, const Params& P_, double rho_, int dimMode)
  : NDMaterial(tag, ND_TAG_LadrunoRCConcrete),
    P(P_), rho(rho_), dim(dimMode), ncomp(6), condense(false), cEps33(0.0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->setupDim();
  this->revertToStart();
}

LadrunoRCConcrete::~LadrunoRCConcrete() {}

// ===========================================================================
//  dimensional view setup
// ===========================================================================
void LadrunoRCConcrete::setupDim(void)
{
  switch (dim) {
    case DIM_PSTRESS:    ncomp = 3; { int m[3]={0,1,3};     for(int a=0;a<3;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_PLATEFIBER: ncomp = 5; { int m[5]={0,1,3,4,5}; for(int a=0;a<5;a++) vmap[a]=m[a]; } condense=true;  break;
    case DIM_3D:
    default:             ncomp = 6; { int m[6]={0,1,2,3,4,5}; for(int a=0;a<6;a++) vmap[a]=m[a]; } condense=false; break;
  }
  stressOut.resize(ncomp);
  strainOut.resize(ncomp);
  tangentOut.resize(ncomp, ncomp);
}

void LadrunoRCConcrete::condenseTangent(void)
{
  double d22 = Dtan[2][2];
  if (fabs(d22) < 1.0e-300) return;
  double col[6], row[6];
  for (int i = 0; i < 6; i++) { col[i] = Dtan[i][2]; row[i] = Dtan[2][i]; }
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++)
      Dtan[I][J] -= col[I]*row[J]/d22;
}

// ===========================================================================
//  strain interface  (ENGINEERING shear — no tensor conversion)
// ===========================================================================
extern double ops_Dt;   // current (pseudo-)time increment, from OPS_Globals

void LadrunoRCConcrete::integrate(bool do_implex, double tfac)
{
  status = returnMap3D(P, strain6, histN, stress6, Dtan, histTr, true, 0, do_implex, tfac);
}

double LadrunoRCConcrete::implexTimeFactor(void) const
{
  // Extrapolation factor tf = (dt_n / dt_{n-1}) * alpha. In DYNAMIC analysis with a
  // smooth time step this is ~alpha. In STATIC analysis the "time" is the load factor,
  // whose increment is erratic (and resets across loadConst), so the raw ratio can be
  // garbage (huge/negative) and would blow up the threshold extrapolation. Guard it:
  // fall back to alpha whenever the ratio is non-finite/non-positive, and clamp the
  // acceleration to a safe band so a one-step time spike cannot detonate the damage.
  if (!P.implex || !commitDone) return P.implexAlpha;
  if (dtime_n_commit <= 0.0 || dtime_n <= 0.0) return P.implexAlpha;
  double tf = dtime_n / dtime_n_commit * P.implexAlpha;
  if (!(tf > 0.0)) return P.implexAlpha;        // NaN / non-positive guard
  const double tfMax = 2.0 * (P.implexAlpha > 0.0 ? P.implexAlpha : 1.0);
  if (tf > tfMax) tf = tfMax;
  return tf;
}

int LadrunoRCConcrete::setTrialStrain(const Vector& e)
{
  if (P.implex) {
    dtime_n = ops_Dt;
    if (!commitDone) { dtime_0 = dtime_n; dtime_n_commit = dtime_n; }
  }
  const bool   dox = P.implex;
  const double tf  = this->implexTimeFactor();

  double eps33 = strain6[2];
  for (int i = 0; i < 6; i++) strain6[i] = 0.0;
  if (condense) strain6[2] = eps33;
  for (int a = 0; a < ncomp; a++) strain6[vmap[a]] = e(a);   // engineering shear, no factor

  if (!condense) { this->integrate(dox, tf); return (status == STATUS_NO_CONVERGE) ? -1 : 0; }

  // enforce sigma_33 = 0: guarded Newton on eps_33 (= strain6[2]); dSNPO sec 9.4
  const int maxIt = 25;
  double prevAbs = 1.0e300;
  bool converged = false;
  for (int it = 0; it < maxIt; it++) {
    this->integrate(dox, tf);
    double smag = 0.0; for (int i = 0; i < 6; i++) smag += stress6[i]*stress6[i];
    smag = sqrt(smag);
    double tol22 = 1.0e-10 * (smag > 1.0 ? smag : 1.0);
    double abs33 = fabs(stress6[2]);
    if (abs33 <= tol22) { converged = true; break; }
    double d22 = Dtan[2][2];
    double step;
    if (fabs(d22) < 1.0e-14 * P.E) step = (stress6[2] > 0.0 ? 1.0 : -1.0) * 1.0e-6;  // bisection nudge
    else                          step = stress6[2] / d22;
    double lambda = (abs33 > prevAbs) ? 0.5 : 1.0;   // damp if residual grew
    strain6[2] -= lambda * step;
    prevAbs = abs33;
  }
  if (!converged) {
    status = STATUS_NO_CONVERGE;
    opserr << "WARNING LadrunoRCConcrete: sigma_33 condensation did not converge (tag "
           << this->getTag() << ", |s33|=" << fabs(stress6[2]) << ")\n";
  }
  this->condenseTangent();
  return (status == STATUS_NO_CONVERGE) ? -1 : 0;
}

int LadrunoRCConcrete::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v)
{
  Vector ne(ncomp);
  for (int a = 0; a < ncomp; a++) ne(a) = strain6[vmap[a]] + v(a);
  return this->setTrialStrain(ne);
}
int LadrunoRCConcrete::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoRCConcrete::getStress(void)
{
  for (int a = 0; a < ncomp; a++) stressOut(a) = stress6[vmap[a]];   // true stress
  return stressOut;
}
const Vector& LadrunoRCConcrete::getStrain(void)
{
  for (int a = 0; a < ncomp; a++) strainOut(a) = strain6[vmap[a]];   // engineering shear
  return strainOut;
}
const Matrix& LadrunoRCConcrete::getTangent(void)
{
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = Dtan[vmap[a]][vmap[b]];
  return tangentOut;
}
const Matrix& LadrunoRCConcrete::getInitialTangent(void)
{
  double C0[6][6];
  elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < ncomp; a++)
    for (int b = 0; b < ncomp; b++) tangentOut(a, b) = C0[vmap[a]][vmap[b]];
  return tangentOut;
}

const char* LadrunoRCConcrete::getType(void) const
{
  switch (dim) {
    case DIM_PSTRESS:    return "PlaneStress";
    case DIM_PLATEFIBER: return "PlateFiber";
    default:             return "ThreeDimensional";
  }
}
int LadrunoRCConcrete::getOrder(void) const { return ncomp; }

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRCConcrete::commitState(void)
{
  if (P.implex) {
    // IMPL-EX commit: the trial currently holds the EXPLICIT response (extrapolated,
    // frozen damage). Re-integrate IMPLICITLY at the converged strain to advance the
    // TRUE thresholds for next step's extrapolation and measure the implex error.
    double dt_ex = histTr.dt_bar, dc_ex = histTr.dc_bar;
    this->integrate(false, 1.0);                       // implicit pass -> histTr = true state
    implexError = fmax(fabs(dt_ex - histTr.dt_bar), fabs(dc_ex - histTr.dc_bar));
    // roll n -> n-1: the previously committed thresholds become the 'old' generation
    histTr.xt_old = histN.xt; histTr.xc_old = histN.xc; histTr.eps1_old = histN.eps1;
    dtime_n_commit = dtime_n;
    // advisory error control: warn if the implex extrapolation error is large (the user
    // should reduce the time/step size). Automatic dt-reduction is deferred.
    if (P.implexControl && implexError > P.implexErrTol)
      opserr << "WARNING LadrunoRCConcrete (tag " << this->getTag() << "): IMPL-EX error "
             << implexError << " > tol " << P.implexErrTol << " — reduce the step size.\n";
  }
  histN = histTr;
  cEps33 = strain6[2];
  commitDone = true;
  return 0;
}
int LadrunoRCConcrete::revertToLastCommit(void)
{
  histTr = histN;
  if (condense) strain6[2] = cEps33;
  return 0;
}
int LadrunoRCConcrete::revertToStart(void)
{
  histZero(histN);
  histZero(histTr);
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  cEps33 = 0.0;
  implexError = 0.0; dtime_n = dtime_n_commit = dtime_0 = 0.0; commitDone = false;
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  return 0;
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoRCConcrete::getCopy(void)
{
  return new LadrunoRCConcrete(this->getTag(), P, rho, dim);
}
NDMaterial* LadrunoRCConcrete::getCopy(const char* type)
{
  int d = -1;
  if      (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) d = DIM_3D;
  else if (strcmp(type, "PlateFiber") == 0)                                  d = DIM_PLATEFIBER;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0) d = DIM_PSTRESS;
  if (d < 0) return NDMaterial::getCopy(type);
  return new LadrunoRCConcrete(this->getTag(), P, rho, d);
}

// ===========================================================================
//  parallel  (serialize params + backbones + committed history)
// ===========================================================================
static const int RC_SCHEMA_VERSION = 4;    // bump when the wire layout changes (hard-checked in recvSelf); v2 = +IMPL-EX; v3 = +shearRetFactor; v4 = +tension stiffening
static const int RC_NSCALAR = 1 /*schemaVersion*/ + 3 /*tag,dim,rho*/
                            + 10 /*E,nu,Kc,fcft,betaFloor,cdf,eta,betaOn,lubRed,tanMode*/
                            + 9 /*interlockOn,shearRetMode,shearRetFactor,aggSize,crackStrain,crackSpacing,lch,betaSrMin,sqrtFc*/
                            + 1 /*interlockCyclic*/
                            + 4 /*xcrackOn,degKappa,degSlipRef,degMin*/
                            + 5 /*implex,implexAlpha,implexControl,implexErrTol,implexTimeRedLim*/
                            + 4 /*tensStiffMode,tensStiffC,tensStiffAlpha,ftPeak*/
                            + 5 /*dtime_n,dtime_n_commit,dtime_0,commitDone,implexError*/;
static const int RC_BACK = 1 + 3*MAXPTS;   // n + x[]+y[]+q[]
static const int RC_HIST = 6 + 6 + 6        // stress_eff, strain, (xt,xc,dt_bar,dc_bar,beta,eps1)
                         + 5                // + (cracked,crackC,crackS,wmax,betaSr)
                         + 2                // + Phase-2b (tauCr,gammaCr)
                         + 2                // + Phase-2b.2b (cracked2,slipCum)
                         + 3;               // + Phase-4 IMPL-EX (xt_old,xc_old,eps1_old)
static const int RC_DATA = RC_NSCALAR + 2*RC_BACK + RC_HIST + 1 /*cEps33*/;

int LadrunoRCConcrete::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(RC_DATA);
  int c = 0;
  data(c++) = RC_SCHEMA_VERSION;
  data(c++) = this->getTag();
  data(c++) = dim;
  data(c++) = rho;
  data(c++) = P.E; data(c++) = P.nu; data(c++) = P.Kc; data(c++) = P.fcft_ratio;
  data(c++) = P.betaFloor; data(c++) = P.cdf; data(c++) = P.eta;
  data(c++) = P.betaOn ? 1.0 : 0.0;
  data(c++) = P.lublinerTCReduced ? 1.0 : 0.0;
  data(c++) = P.tangentMode;
  data(c++) = P.interlockOn ? 1.0 : 0.0;
  data(c++) = P.interlockCyclic ? 1.0 : 0.0;
  data(c++) = P.shearRetMode; data(c++) = P.shearRetFactor;
  data(c++) = P.aggSize; data(c++) = P.crackStrain; data(c++) = P.crackSpacing;
  data(c++) = P.lch; data(c++) = P.betaSrMin; data(c++) = P.sqrtFc;
  data(c++) = P.xcrackOn ? 1.0 : 0.0;
  data(c++) = P.degKappa; data(c++) = P.degSlipRef; data(c++) = P.degMin;
  data(c++) = P.implex ? 1.0 : 0.0;
  data(c++) = P.implexAlpha;
  data(c++) = P.implexControl ? 1.0 : 0.0;
  data(c++) = P.implexErrTol; data(c++) = P.implexTimeRedLim;
  data(c++) = P.tensStiffMode; data(c++) = P.tensStiffC;
  data(c++) = P.tensStiffAlpha; data(c++) = P.ftPeak;
  data(c++) = dtime_n; data(c++) = dtime_n_commit; data(c++) = dtime_0;
  data(c++) = commitDone ? 1.0 : 0.0; data(c++) = implexError;
  data(c++) = P.ht.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.ht.q[i];
  data(c++) = P.hc.n;
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.x[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.y[i];
  for (int i = 0; i < MAXPTS; i++) data(c++) = P.hc.q[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.stress_eff[i];
  for (int i = 0; i < 6; i++) data(c++) = histN.strain[i];
  data(c++) = histN.xt; data(c++) = histN.xc;
  data(c++) = histN.dt_bar; data(c++) = histN.dc_bar;
  data(c++) = histN.beta; data(c++) = histN.eps1;
  data(c++) = histN.cracked; data(c++) = histN.crackC; data(c++) = histN.crackS;
  data(c++) = histN.wmax; data(c++) = histN.betaSr;
  data(c++) = histN.tauCr; data(c++) = histN.gammaCr;
  data(c++) = histN.cracked2; data(c++) = histN.slipCum;
  data(c++) = histN.xt_old; data(c++) = histN.xc_old; data(c++) = histN.eps1_old;
  data(c++) = cEps33;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoRCConcrete::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(RC_DATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCConcrete::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  int ver = (int)data(c++);
  if (ver != RC_SCHEMA_VERSION) {
    opserr << "LadrunoRCConcrete::recvSelf - unsupported schema version " << ver
           << " (expected " << RC_SCHEMA_VERSION << ")\n";
    return -1;
  }
  this->setTag((int)data(c++));
  dim = (int)data(c++);
  rho = data(c++);
  P.E = data(c++); P.nu = data(c++); P.Kc = data(c++); P.fcft_ratio = data(c++);
  P.betaFloor = data(c++); P.cdf = data(c++); P.eta = data(c++);
  P.betaOn = (data(c++) != 0.0);
  P.lublinerTCReduced = (data(c++) != 0.0);
  P.tangentMode = (int)data(c++);
  P.interlockOn = (data(c++) != 0.0);
  P.interlockCyclic = (data(c++) != 0.0);
  P.shearRetMode = (int)data(c++); P.shearRetFactor = data(c++);
  P.aggSize = data(c++); P.crackStrain = data(c++); P.crackSpacing = data(c++);
  P.lch = data(c++); P.betaSrMin = data(c++); P.sqrtFc = data(c++);
  P.xcrackOn = (data(c++) != 0.0);
  P.degKappa = data(c++); P.degSlipRef = data(c++); P.degMin = data(c++);
  P.implex = (data(c++) != 0.0);
  P.implexAlpha = data(c++);
  P.implexControl = (data(c++) != 0.0);
  P.implexErrTol = data(c++); P.implexTimeRedLim = data(c++);
  P.tensStiffMode = (int)data(c++); P.tensStiffC = data(c++);
  P.tensStiffAlpha = data(c++); P.ftPeak = data(c++);
  dtime_n = data(c++); dtime_n_commit = data(c++); dtime_0 = data(c++);
  commitDone = (data(c++) != 0.0); implexError = data(c++);
  P.ht.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.ht.q[i] = data(c++);
  P.hc.n = (int)data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.x[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.y[i] = data(c++);
  for (int i = 0; i < MAXPTS; i++) P.hc.q[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.stress_eff[i] = data(c++);
  for (int i = 0; i < 6; i++) histN.strain[i] = data(c++);
  histN.xt = data(c++); histN.xc = data(c++);
  histN.dt_bar = data(c++); histN.dc_bar = data(c++);
  histN.beta = data(c++); histN.eps1 = data(c++);
  histN.cracked = data(c++); histN.crackC = data(c++); histN.crackS = data(c++);
  histN.wmax = data(c++); histN.betaSr = data(c++);
  histN.tauCr = data(c++); histN.gammaCr = data(c++);
  histN.cracked2 = data(c++); histN.slipCum = data(c++);
  histN.xt_old = data(c++); histN.xc_old = data(c++); histN.eps1_old = data(c++);
  cEps33 = data(c++);

  this->setupDim();
  histTr = histN;
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  if (condense) strain6[2] = cEps33;
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  return 0;
}

void LadrunoRCConcrete::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoRCConcrete (RC plastic-damage + MCFT compression softening)" << endln;
  s << "  tag   : " << this->getTag() << endln;
  s << "  E, nu : " << P.E << ", " << P.nu << endln;
  s << "  Kc    : " << P.Kc << "   fc/ft ratio: " << P.fcft_ratio << endln;
  s << "  beta  : " << (P.betaOn ? "ON" : "off") << "  floor=" << P.betaFloor
    << "  lublinerReduced=" << (P.lublinerTCReduced ? 1 : 0) << endln;
  s << "  interlock: " << (P.interlockOn ? "ON" : "off")
    << (P.interlockCyclic ? " (cyclic)" : "") << (P.xcrackOn ? " (xcrack)" : "")
    << "  a_g=" << P.aggSize
    << "  eps_cr=" << P.crackStrain << "  s_theta=" << P.crackSpacing
    << "  betaSrMin=" << P.betaSrMin;
  if (P.xcrackOn) s << "  degKappa=" << P.degKappa << " degSlipRef=" << P.degSlipRef << " degMin=" << P.degMin;
  if (P.interlockOn) {
    const char* srn = (P.shearRetMode == 1) ? "const" : (P.shearRetMode == 2) ? "dsfm"
                    : (P.shearRetMode == 3) ? "rots" : "mcft";
    s << "  shearRetention=" << srn;
    if (P.shearRetMode == 1) s << "(mu=" << P.shearRetFactor << ")";
  }
  s << endln;
  s << "  implex: " << (P.implex ? "ON" : "off");
  if (P.implex) s << "  alpha=" << P.implexAlpha
                  << (P.implexControl ? " (control)" : "") << "  lastError=" << implexError;
  s << endln;
  s << "  tensStiff: " << (P.tensStiffMode == 1 ? "vc" : P.tensStiffMode == 2 ? "cm" : "off");
  if (P.tensStiffMode == 1) s << "  c=" << P.tensStiffC;
  if (P.tensStiffMode == 2) s << "  alpha=" << P.tensStiffAlpha;
  if (P.tensStiffMode != 0) s << "  ft=" << P.ftPeak;
  s << endln;
  s << "  view  : " << this->getType() << " (order " << ncomp << ")" << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoRCConcrete::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];
  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "damage") == 0 || strcmp(a, "Damage") == 0)
    return new MaterialResponse(this, 4, Vector(2));   // (dt_bar, dc_bar)
  if (strcmp(a, "beta") == 0 || strcmp(a, "compressionSoftening") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "eps1") == 0 || strcmp(a, "transverseTensileStrain") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "betaShear") == 0 || strcmp(a, "shearRetention") == 0)
    return new MaterialResponse(this, 7, Vector(1));
  if (strcmp(a, "crackState") == 0 || strcmp(a, "crackNormal") == 0)
    return new MaterialResponse(this, 8, Vector(4));   // (cracked, crackC, crackS, wmax)
  if (strcmp(a, "crackShear") == 0)
    return new MaterialResponse(this, 9, Vector(2));   // (tauCr, gammaCr) Phase-2b
  if (strcmp(a, "xcrackState") == 0)
    return new MaterialResponse(this, 10, Vector(2));  // (cracked2, slipCum) Phase-2b.2b
  if (strcmp(a, "implexError") == 0 || strcmp(a, "ImplexError") == 0)
    return new MaterialResponse(this, 11, Vector(1));  // IMPL-EX |dt_ex - dt_im| Phase-4
  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoRCConcrete::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1: if (matInfo.theVector) *(matInfo.theVector) = this->getStress();  return 0;
    case 2: if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();  return 0;
    case 3: if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent(); return 0;
    case 4: if (matInfo.theVector) { Vector& v = *(matInfo.theVector); v(0) = histTr.dt_bar; v(1) = histTr.dc_bar; } return 0;
    case 5: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.beta;   return 0;
    case 6: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.eps1;   return 0;
    case 7: if (matInfo.theVector) (*(matInfo.theVector))(0) = histTr.betaSr; return 0;
    case 8: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.cracked; v(1) = histTr.crackC; v(2) = histTr.crackS; v(3) = histTr.wmax; } return 0;
    case 9: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.tauCr; v(1) = histTr.gammaCr; } return 0;
    case 10: if (matInfo.theVector) { Vector& v = *(matInfo.theVector);
              v(0) = histTr.cracked2; v(1) = histTr.slipCum; } return 0;
    case 11: if (matInfo.theVector) (*(matInfo.theVector))(0) = implexError; return 0;
    default: return -1;
  }
}
