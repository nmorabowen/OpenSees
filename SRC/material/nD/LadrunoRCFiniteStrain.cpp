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

// Ladruno RC plastic-damage FINITE-STRAIN view — see LadrunoRCFiniteStrain.h. classTag 33018.

#include "LadrunoRCFiniteStrain.h"
#include "LogStrainKernel.h"
#include <elementAPI.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <OPS_Globals.h>   // ops_Dt (IMPL-EX) + ops_TheActiveElement (Phase 3b lch latch)
#include <Element.h>       // Element::getCharacteristicLength()
#include <string.h>
#include <math.h>
#include <vector>

using namespace ladruno_rc_kernel;

// ===========================================================================
//  parser:  nDMaterial LadrunoRCFiniteStrain $tag $E $nu  (same grammar as
//  LadrunoRCConcrete; consumed by element LadrunoBrick -geom finite). All RC
//  physics flags (-beta/-interlock/-cyclic/-xcrack/-shearRetention/-tensStiff/
//  -autoRegularization/-implex) ride the SHARED kernel unchanged.
// ===========================================================================
void* OPS_LadrunoRCFiniteStrain(void)
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: too few args.\n"
           << "  want: tag E nu -Ce {..} -Cs {..} <-Cd {..}> -Te {..} -Ts {..} <-Td {..}>"
           << " <-Kc Kc> <-beta> <-betaFloor f> <-lublinerReduced> <-rho rho> <flags...>\n";
    return 0;
  }

  int numData = 1;
  int tag = 0;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: invalid tag.\n";
    return 0;
  }
  double E = 0.0, nu = 0.0;
  if (OPS_GetDoubleInput(&numData, &E) < 0)  { opserr << "LadrunoRCFiniteStrain: invalid E.\n";  return 0; }
  if (OPS_GetDoubleInput(&numData, &nu) < 0) { opserr << "LadrunoRCFiniteStrain: invalid nu.\n"; return 0; }

  std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
  double Kc = 2.0 / 3.0, betaFloor = 0.1, rho = 0.0;
  bool betaOn = false, lubRed = false;
  int  tanMode = 0;
  bool   interlockOn = false, interlockCyclic = false;
  double aggSize = 16.0, crackStrain = 0.0, crackSpacing = 0.0, lch = 0.0, betaSrMin = 0.01;
  int    shearRetMode = 0;
  double shearRetFactor = 0.4;
  bool   xcrackOn = false;
  double degKappa = 0.5, degSlipRef = 0.01, degMin = 0.1;
  bool   implexOn = false, implexCtrl = false;
  double implexAlpha = 1.0, implexErrTol = 0.05, implexTimeRedLim = 0.01;
  int    tensStiffMode = 0;
  double tensStiffC = 500.0, tensStiffAlpha = 1.0;
  bool   autoReg = false;
  double lchRef = 1.0;

  auto readList = [](std::vector<double>& v) {
    v.clear();
    int nd = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      double item;
      int oldn = OPS_GetNumRemainingInputArgs();
      if (OPS_GetDoubleInput(&nd, &item) < 0) {
        int newn = OPS_GetNumRemainingInputArgs();
        if (newn < oldn) OPS_ResetCurrentInputArg(-1);
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
    else if (strcmp(opt, "-agg") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &aggSize) < 0)      { opserr << "LadrunoRCFiniteStrain: -agg needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-crackStrain") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackStrain) < 0)  { opserr << "LadrunoRCFiniteStrain: -crackStrain needs a value.\n";  return 0; } }
    else if (strcmp(opt, "-crackSpacing") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &crackSpacing) < 0) { opserr << "LadrunoRCFiniteStrain: -crackSpacing needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-lch") == 0)          { int nd = 1; if (OPS_GetDoubleInput(&nd, &lch) < 0)          { opserr << "LadrunoRCFiniteStrain: -lch needs a value.\n";          return 0; } }
    else if (strcmp(opt, "-betaSrMin") == 0)    { int nd = 1; if (OPS_GetDoubleInput(&nd, &betaSrMin) < 0)    { opserr << "LadrunoRCFiniteStrain: -betaSrMin needs a value.\n";    return 0; } }
    else if (strcmp(opt, "-xcrack") == 0)            xcrackOn = true;
    else if (strcmp(opt, "-degKappa") == 0)     { int nd = 1; if (OPS_GetDoubleInput(&nd, &degKappa) < 0)     { opserr << "LadrunoRCFiniteStrain: -degKappa needs a value.\n";     return 0; } }
    else if (strcmp(opt, "-degSlipRef") == 0)   { int nd = 1; if (OPS_GetDoubleInput(&nd, &degSlipRef) < 0)   { opserr << "LadrunoRCFiniteStrain: -degSlipRef needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-degMin") == 0)       { int nd = 1; if (OPS_GetDoubleInput(&nd, &degMin) < 0)       { opserr << "LadrunoRCFiniteStrain: -degMin needs a value.\n";       return 0; } }
    else if (strcmp(opt, "-implex") == 0)            implexOn = true;
    else if (strcmp(opt, "-implexAlpha") == 0)  { int nd = 1; if (OPS_GetDoubleInput(&nd, &implexAlpha) < 0)   { opserr << "LadrunoRCFiniteStrain: -implexAlpha needs a value.\n";   return 0; } }
    else if (strcmp(opt, "-implexControl") == 0) {
      implexCtrl = true; int nd = 1;
      if (OPS_GetDoubleInput(&nd, &implexErrTol) < 0 || OPS_GetDoubleInput(&nd, &implexTimeRedLim) < 0) {
        opserr << "LadrunoRCFiniteStrain: -implexControl needs $errTol $timeReductionLimit.\n"; return 0; }
    }
    else if (strcmp(opt, "-shearRetention") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "mcft")  == 0) shearRetMode = 0;
      else if (m && strcmp(m, "const") == 0) shearRetMode = 1;
      else if (m && strcmp(m, "dsfm")  == 0) shearRetMode = 2;
      else if (m && strcmp(m, "rots")  == 0) shearRetMode = 3;
      else { opserr << "LadrunoRCFiniteStrain: -shearRetention needs {mcft|const|dsfm|rots}.\n"; return 0; }
    }
    else if (strcmp(opt, "-shearRetFactor") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &shearRetFactor) < 0) { opserr << "LadrunoRCFiniteStrain: -shearRetFactor needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-tensStiff") == 0) {
      const char* m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
      if      (m && strcmp(m, "vc") == 0) tensStiffMode = 1;
      else if (m && strcmp(m, "cm") == 0) tensStiffMode = 2;
      else { opserr << "LadrunoRCFiniteStrain: -tensStiff needs {vc|cm}.\n"; return 0; }
    }
    else if (strcmp(opt, "-tensStiffC") == 0)     { int nd = 1; if (OPS_GetDoubleInput(&nd, &tensStiffC) < 0)     { opserr << "LadrunoRCFiniteStrain: -tensStiffC needs a value.\n";     return 0; } }
    else if (strcmp(opt, "-tensStiffAlpha") == 0) { int nd = 1; if (OPS_GetDoubleInput(&nd, &tensStiffAlpha) < 0) { opserr << "LadrunoRCFiniteStrain: -tensStiffAlpha needs a value.\n"; return 0; } }
    else if (strcmp(opt, "-autoRegularization") == 0) {
      autoReg = true; int nd = 1;
      if (OPS_GetDoubleInput(&nd, &lchRef) < 0) { opserr << "LadrunoRCFiniteStrain: -autoRegularization needs $lch_ref.\n"; return 0; }
    }
    // unknown tokens ignored (forward-compat)
  }

  if (xcrackOn) interlockCyclic = true;
  if (shearRetMode == 1 || shearRetMode == 2) interlockCyclic = true;
  if (shearRetMode != 0) interlockOn = true;
  if (interlockCyclic && !interlockOn) interlockOn = true;

  if (tensStiffMode == 1 && tensStiffC <= 0.0) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: -tensStiffC must be > 0.\n";
    return 0;
  }
  if (tensStiffMode == 2 && tensStiffC != 500.0)
    opserr << "WARNING nDMaterial LadrunoRCFiniteStrain: -tensStiffC is ignored in cm mode "
              "(cm uses the fixed Collins-Mitchell 500); use -tensStiff vc for a tunable c.\n";
  if (autoReg && lchRef <= 0.0) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: -autoRegularization $lch_ref must be > 0.\n";
    return 0;
  }
  if (Ce.size() < 2 || Cs.size() != Ce.size()) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: need -Ce and -Cs of equal length (>=2).\n";
    return 0;
  }
  if (Te.size() < 2 || Ts.size() != Te.size()) {
    opserr << "nDMaterial LadrunoRCFiniteStrain error: need -Te and -Ts of equal length (>=2).\n";
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
  P.tensStiffAlpha = tensStiffAlpha; P.ftPeak = 0.0;
  P.autoReg = autoReg; P.lchRef = lchRef;

  std::vector<double> Cdf(Ce.size(), 0.0), Tdf(Te.size(), 0.0);
  for (size_t i = 0; i < Cd.size() && i < Cdf.size(); ++i) Cdf[i] = Cd[i];
  for (size_t i = 0; i < Td.size() && i < Tdf.size(); ++i) Tdf[i] = Td[i];
  buildBackbone(P.hc, E, Ce.data(), Cs.data(), Cdf.data(), (int)Ce.size());
  buildBackbone(P.ht, E, Te.data(), Ts.data(), Tdf.data(), (int)Te.size());
  setupParams(P);

  return new LadrunoRCFiniteStrain(tag, P, rho);
}

// ===========================================================================
//  ctors / dtor
// ===========================================================================
LadrunoRCFiniteStrain::LadrunoRCFiniteStrain()
  : FiniteStrainNDMaterial(0, ND_TAG_LadrunoRCFiniteStrain),
    rho(0.0), Jdet(1.0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    regularizationDone(false), regLch(0.0), regWarned(false),
    sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), tangentOut(6, 6)
{
  P.E = 1.0; P.nu = 0.0; P.Kc = 2.0/3.0; P.fcft_ratio = 5.0; P.betaFloor = 0.1;
  P.cdf = 0.0; P.eta = 0.0; P.betaOn = false; P.lublinerTCReduced = false; P.tangentMode = 0;
  P.interlockOn = false; P.interlockCyclic = false; P.shearRetMode = 0; P.shearRetFactor = 0.4; P.aggSize = 16.0;
  P.crackStrain = 0.0; P.crackSpacing = 0.0; P.lch = 0.0; P.betaSrMin = 0.01; P.sqrtFc = 0.0;
  P.xcrackOn = false; P.degKappa = 0.5; P.degSlipRef = 0.01; P.degMin = 0.1;
  P.implex = false; P.implexAlpha = 1.0; P.implexControl = false;
  P.implexErrTol = 0.05; P.implexTimeRedLim = 0.01;
  P.tensStiffMode = 0; P.tensStiffC = 500.0; P.tensStiffAlpha = 1.0; P.ftPeak = 0.0;
  P.autoReg = false; P.lchRef = 1.0;
  P.ht.n = 0; P.hc.n = 0;
  this->revertToStart();
}

LadrunoRCFiniteStrain::LadrunoRCFiniteStrain(int tag, const Params& P_, double rho_)
  : FiniteStrainNDMaterial(tag, ND_TAG_LadrunoRCFiniteStrain),
    P(P_), rho(rho_), Jdet(1.0), status(STATUS_OK),
    implexError(0.0), dtime_n(0.0), dtime_n_commit(0.0), dtime_0(0.0), commitDone(false),
    regularizationDone(false), regLch(0.0), regWarned(false),
    sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), tangentOut(6, 6)
{
  this->revertToStart();
}

LadrunoRCFiniteStrain::~LadrunoRCFiniteStrain() {}

// ===========================================================================
//  IMPL-EX time factor + crack-band regularization latch (mirror LadrunoRCConcrete)
// ===========================================================================
double LadrunoRCFiniteStrain::implexTimeFactor(void) const
{
  if (!P.implex || !commitDone) return P.implexAlpha;
  if (dtime_n_commit <= 0.0 || dtime_n <= 0.0) return P.implexAlpha;
  double tf = dtime_n / dtime_n_commit * P.implexAlpha;
  if (!(tf > 0.0)) return P.implexAlpha;
  const double tfMax = 2.0 * (P.implexAlpha > 0.0 ? P.implexAlpha : 1.0);
  if (tf > tfMax) tf = tfMax;
  return tf;
}

int LadrunoRCFiniteStrain::regularizeIfNeeded(void)
{
  if (!P.autoReg || regularizationDone) return 0;
  double lch = (P.lch > 0.0) ? P.lch
             : (ops_TheActiveElement ? ops_TheActiveElement->getCharacteristicLength() : 0.0);
  if (!(lch > 0.0)) {
    if (!regWarned) {
      opserr << "FATAL LadrunoRCFiniteStrain (tag " << this->getTag()
             << "): -autoRegularization requested but no characteristic length is available "
                "(no active element and no -lch). Refusing to run with a non-objective softening "
                "law — supply -lch $l or use the material inside an element.\n";
      regWarned = true;
    }
    regLch = 0.0;
    return -1;
  }
  regLch = lch;
  ladruno_rc_kernel::regularize(P.ht, lch, P.lchRef, P.E);
  ladruno_rc_kernel::regularize(P.hc, lch, P.lchRef, P.E);
  regularizationDone = true;
  return 0;
}

void LadrunoRCFiniteStrain::integrate(bool do_implex, double tfac)
{
  status = returnMap3D(P, strain6, histN, stress6, Dtan, histTr, true, 0, do_implex, tfac);
}

// ===========================================================================
//  the finite-strain seam: F -> Hencky -> kernel -> Cauchy + spatial tangent
// ===========================================================================
int LadrunoRCFiniteStrain::setTrialF(const Matrix& F)
{
  using namespace logstrain_kernel;

  if (this->regularizeIfNeeded() < 0) { status = STATUS_NO_CONVERGE; return -1; }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Ftrial9[3*i+j] = F(i, j);

  // total left Cauchy–Green B = F Fᵀ (the elastic stretch — no plastic strain in the
  // RC spine) and J = det F. Reject a non-positive Jacobian up front (an inverted /
  // degenerate element under a large step): the element treats <0 as "cut the step".
  Jdet = mat3_det(Ftrial9);
  if (Jdet <= 0.0) {
    opserr << "LadrunoRCFiniteStrain::setTrialF - non-positive det F (" << Jdet
           << "); element inverted/degenerate\n";
    return -1;
  }
  double Ft[9]; mat3_transpose(Ftrial9, Ft);
  mat3_mul(Ftrial9, Ft, Btotal9);
  // symmetrise against round-off (B = F Fᵀ is analytically symmetric)
  Btotal9[1] = Btotal9[3] = 0.5*(Btotal9[1]+Btotal9[3]);
  Btotal9[2] = Btotal9[6] = 0.5*(Btotal9[2]+Btotal9[6]);
  Btotal9[5] = Btotal9[7] = 0.5*(Btotal9[5]+Btotal9[7]);

  // εᵉ = ½ ln B (engineering-shear Voigt {00,11,22,01,12,02} — matches the RC kernel)
  hencky_voigt(Btotal9, strain6);

  // IMPL-EX bookkeeping (same as the small-strain material)
  if (P.implex) {
    dtime_n = ops_Dt;
    if (!commitDone) { dtime_0 = dtime_n; dtime_n_commit = dtime_n; }
  }
  const bool   dox = P.implex;
  const double tf  = this->implexTimeFactor();

  this->integrate(dox, tf);              // stress6 = Kirchhoff τ, Dtan = ∂τ/∂εᵉ
  if (status == STATUS_NO_CONVERGE) return -1;

  // Cauchy σ = τ/J and material spatial tangent c = (1/2J)[D:L:B] at B
  double D6[36];
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) D6[6*I+J] = Dtan[I][J];
  double sig6[6], c6[36];
  assemble_material(Btotal9, D6, stress6, Jdet, sig6, c6);
  for (int a = 0; a < 6; a++) { sigmaCauchy(a) = sig6[a]; henckyStrain(a) = strain6[a]; }
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) aTangent(I, J) = c6[6*I+J];
  return 0;
}

const Vector& LadrunoRCFiniteStrain::getStress(void)  { return sigmaCauchy; }
const Vector& LadrunoRCFiniteStrain::getStrain(void)  { return henckyStrain; }
const Matrix& LadrunoRCFiniteStrain::getTangent(void) { return aTangent; }

int LadrunoRCFiniteStrain::getSpatialTangentTensor(double c[3][3][3][3])
{
  double D6[36];
  for (int I = 0; I < 6; I++) for (int J = 0; J < 6; J++) D6[6*I+J] = Dtan[I][J];
  logstrain_kernel::spatial_tangent_full(Btotal9, D6, Jdet, c);
  return 0;
}

const Matrix& LadrunoRCFiniteStrain::getInitialTangent(void)
{
  double C0[6][6];
  elasticTangent(P.E, P.nu, C0);          // spatial elastic modulus at F = I (B = I, J = 1)
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) tangentOut(a, b) = C0[a][b];
  return tangentOut;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoRCFiniteStrain::commitState(void)
{
  if (P.implex) {
    double dt_ex = histTr.dt_bar, dc_ex = histTr.dc_bar;
    this->integrate(false, 1.0);                       // implicit pass -> histTr = true state
    implexError = fmax(fabs(dt_ex - histTr.dt_bar), fabs(dc_ex - histTr.dc_bar));
    histTr.xt_old = histN.xt; histTr.xc_old = histN.xc; histTr.eps1_old = histN.eps1;
    dtime_n_commit = dtime_n;
    if (P.implexControl && implexError > P.implexErrTol)
      opserr << "WARNING LadrunoRCFiniteStrain (tag " << this->getTag() << "): IMPL-EX error "
             << implexError << " > tol " << P.implexErrTol << " — reduce the step size.\n";
  }
  histN = histTr;
  commitDone = true;
  return 0;
}
int LadrunoRCFiniteStrain::revertToLastCommit(void)
{
  histTr = histN;
  return 0;
}
int LadrunoRCFiniteStrain::revertToStart(void)
{
  histZero(histN);
  histZero(histTr);
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; sigmaCauchy(i) = 0.0; henckyStrain(i) = 0.0; }
  Jdet = 1.0;
  for (int i = 0; i < 9; i++) { Ftrial9[i] = (i % 4 == 0) ? 1.0 : 0.0; Btotal9[i] = (i % 4 == 0) ? 1.0 : 0.0; }
  implexError = 0.0; dtime_n = dtime_n_commit = dtime_0 = 0.0; commitDone = false;
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) { Dtan[a][b] = C0[a][b]; aTangent(a, b) = C0[a][b]; }
  return 0;
}

// ===========================================================================
//  copies — propagate the regularization latch (mirror LadrunoRCConcrete)
// ===========================================================================
NDMaterial* LadrunoRCFiniteStrain::getCopy(void)
{
  LadrunoRCFiniteStrain* c = new LadrunoRCFiniteStrain(this->getTag(), P, rho);
  c->regularizationDone = regularizationDone; c->regLch = regLch;
  return c;
}
NDMaterial* LadrunoRCFiniteStrain::getCopy(const char* type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    return this->getCopy();
  return NDMaterial::getCopy(type);
}

// ===========================================================================
//  parallel  (serialize params + backbones + committed history; no dim/cEps33,
//  the finite view is 3D-only and condensation-free — its own schema)
// ===========================================================================
static const int RCF_SCHEMA_VERSION = 1;
static const int RCF_NSCALAR = 1 /*schemaVersion*/ + 2 /*tag,rho*/
                            + 10 /*E,nu,Kc,fcft,betaFloor,cdf,eta,betaOn,lubRed,tanMode*/
                            + 9 /*interlockOn,shearRetMode,shearRetFactor,aggSize,crackStrain,crackSpacing,lch,betaSrMin,sqrtFc*/
                            + 1 /*interlockCyclic*/
                            + 4 /*xcrackOn,degKappa,degSlipRef,degMin*/
                            + 5 /*implex,implexAlpha,implexControl,implexErrTol,implexTimeRedLim*/
                            + 4 /*tensStiffMode,tensStiffC,tensStiffAlpha,ftPeak*/
                            + 4 /*autoReg,lchRef,regularizationDone,regLch*/
                            + 5 /*dtime_n,dtime_n_commit,dtime_0,commitDone,implexError*/;
static const int RCF_BACK = 1 + 3*MAXPTS;
static const int RCF_HIST = 6 + 6 + 6 + 5 + 2 + 2 + 3;
static const int RCF_DATA = RCF_NSCALAR + 2*RCF_BACK + RCF_HIST;

int LadrunoRCFiniteStrain::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(RCF_DATA);
  int c = 0;
  data(c++) = RCF_SCHEMA_VERSION;
  data(c++) = this->getTag();
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
  data(c++) = P.autoReg ? 1.0 : 0.0; data(c++) = P.lchRef;
  data(c++) = regularizationDone ? 1.0 : 0.0; data(c++) = regLch;
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

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCFiniteStrain::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoRCFiniteStrain::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(RCF_DATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoRCFiniteStrain::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  int ver = (int)data(c++);
  if (ver != RCF_SCHEMA_VERSION) {
    opserr << "LadrunoRCFiniteStrain::recvSelf - unsupported schema version " << ver
           << " (expected " << RCF_SCHEMA_VERSION << ")\n";
    return -1;
  }
  this->setTag((int)data(c++));
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
  P.autoReg = (data(c++) != 0.0); P.lchRef = data(c++);
  regularizationDone = (data(c++) != 0.0); regLch = data(c++);
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

  histTr = histN;
  for (int i = 0; i < 6; i++) { strain6[i] = 0.0; stress6[i] = 0.0; }
  double C0[6][6]; elasticTangent(P.E, P.nu, C0);
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan[a][b] = C0[a][b];
  return 0;
}

void LadrunoRCFiniteStrain::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoRCFiniteStrain (RC plastic-damage, finite-strain Hencky view)" << endln;
  s << "  tag   : " << this->getTag() << endln;
  s << "  E, nu : " << P.E << ", " << P.nu << endln;
  s << "  Kc    : " << P.Kc << "   fc/ft ratio: " << P.fcft_ratio << endln;
  s << "  beta  : " << (P.betaOn ? "ON" : "off") << "  floor=" << P.betaFloor
    << "  lublinerReduced=" << (P.lublinerTCReduced ? 1 : 0) << endln;
  s << "  interlock: " << (P.interlockOn ? "ON" : "off")
    << (P.interlockCyclic ? " (cyclic)" : "") << (P.xcrackOn ? " (xcrack)" : "") << endln;
  s << "  implex: " << (P.implex ? "ON" : "off");
  if (P.implex) s << "  alpha=" << P.implexAlpha << "  lastError=" << implexError;
  s << endln;
  s << "  autoRegularization: " << (P.autoReg ? "ON" : "off");
  if (P.autoReg) s << "  lch_ref=" << P.lchRef << (regularizationDone ? "  lch=" : "  (pending)");
  if (P.autoReg && regularizationDone) s << regLch;
  s << endln;
  s << "  view  : ThreeDimensional finite (order 6, Cauchy/J + spatial tangent)" << endln;
}

// ===========================================================================
//  recordable responses (same set as LadrunoRCConcrete; getStress=Cauchy, getStrain=Hencky)
// ===========================================================================
Response* LadrunoRCFiniteStrain::setResponse(const char** argv, int argc, OPS_Stream& s)
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
    return new MaterialResponse(this, 4, Vector(2));
  if (strcmp(a, "beta") == 0 || strcmp(a, "compressionSoftening") == 0)
    return new MaterialResponse(this, 5, Vector(1));
  if (strcmp(a, "eps1") == 0 || strcmp(a, "transverseTensileStrain") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "betaShear") == 0 || strcmp(a, "shearRetention") == 0)
    return new MaterialResponse(this, 7, Vector(1));
  if (strcmp(a, "crackState") == 0 || strcmp(a, "crackNormal") == 0)
    return new MaterialResponse(this, 8, Vector(4));
  if (strcmp(a, "crackShear") == 0)
    return new MaterialResponse(this, 9, Vector(2));
  if (strcmp(a, "xcrackState") == 0)
    return new MaterialResponse(this, 10, Vector(2));
  if (strcmp(a, "implexError") == 0 || strcmp(a, "ImplexError") == 0)
    return new MaterialResponse(this, 11, Vector(1));
  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoRCFiniteStrain::getResponse(int responseID, Information& matInfo)
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
