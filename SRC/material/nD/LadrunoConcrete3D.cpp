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

// Ladruno: CDPM2-grade solid-concrete plastic-damage nDMaterial (3D). See LadrunoConcrete3D.h and
// Ladruno_implementation/31_ladruno_concrete3d_adr.md. classTag 33017.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoConcrete3D.h>
#include <LadrunoConcrete3DKernel.h>   // the header-only OpenSees-free CDPM2 kernel
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <Element.h>                   // ops_TheActiveElement->getCharacteristicLength()
#include <OPS_Globals.h>
#include <string.h>
#include <math.h>
#include <elementAPI.h>

using Ladruno::Concrete3D::Params;
using Ladruno::Concrete3D::State;

// ===========================================================================
//  OPS parser
//   nDMaterial LadrunoConcrete3D tag? E? nu? fc? ft? Gf? Gc?
//       <-e e? | -kupfer ratio?> <-Df Df?> <-As As?> <-rho rho?>
//       <-hardening qh0? Hp?> <-ductility Ah? Bh? Ch? Dh?>
//       <-lch lch?> <-autoRegularization>
//  (fc, ft are POSITIVE magnitudes; compression is negative in the model.)
// ===========================================================================
void* OPS_LadrunoConcrete3D(void)
{
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING: insufficient args\n";
    opserr << "Want: nDMaterial LadrunoConcrete3D tag? E? nu? fc? ft? Gf? Gc? "
           << "<-e e? | -kupfer ratio?> <-Df Df?> <-As As?> <-rho rho?> "
           << "<-hardening qh0? Hp?> <-ductility Ah? Bh? Ch? Dh?> <-lch lch?> <-autoRegularization>\n";
    return 0;
  }

  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid LadrunoConcrete3D tag\n";
    return 0;
  }

  double d6[6];
  numData = 6;
  if (OPS_GetDoubleInput(&numData, d6) < 0) {
    opserr << "WARNING invalid LadrunoConcrete3D E nu fc ft Gf Gc\n";
    return 0;
  }
  double E = d6[0], nu = d6[1], fc = d6[2], ft = d6[3], Gf = d6[4], Gc = d6[5];

  // defaults (oracle make_material / ADR 6): e from Kupfer 1.16; CDPM2 published ductility table
  double ecc = -1.0;            // <0 => derive from Kupfer ratio below
  double kupfer = 1.16;
  bool haveE = false;
  double Df = 1.0, As = 2.0, rho = 0.0;
  double qh0 = 0.3, Hp = 0.5, Ah = 0.08, Bh = 0.003, Ch = 2.0, Dh = 1.0e-6;
  double lch = 1.0;
  bool autoReg = false;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* flag = OPS_GetString();
    if (strcmp(flag, "-e") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &ecc) < 0) { opserr << "WARNING LadrunoConcrete3D: -e wants e\n"; return 0; }
      haveE = true;
    }
    else if (strcmp(flag, "-kupfer") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &kupfer) < 0) { opserr << "WARNING LadrunoConcrete3D: -kupfer wants fcc/fc\n"; return 0; }
    }
    else if (strcmp(flag, "-Df") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &Df) < 0) { opserr << "WARNING LadrunoConcrete3D: -Df wants Df\n"; return 0; }
    }
    else if (strcmp(flag, "-As") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &As) < 0) { opserr << "WARNING LadrunoConcrete3D: -As wants As\n"; return 0; }
    }
    else if (strcmp(flag, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) { opserr << "WARNING LadrunoConcrete3D: -rho wants rho\n"; return 0; }
    }
    else if (strcmp(flag, "-hardening") == 0) {
      double h[2]; numData = 2;
      if (OPS_GetDoubleInput(&numData, h) < 0) { opserr << "WARNING LadrunoConcrete3D: -hardening wants qh0 Hp\n"; return 0; }
      qh0 = h[0]; Hp = h[1];
    }
    else if (strcmp(flag, "-ductility") == 0) {
      double h[4]; numData = 4;
      if (OPS_GetDoubleInput(&numData, h) < 0) { opserr << "WARNING LadrunoConcrete3D: -ductility wants Ah Bh Ch Dh\n"; return 0; }
      Ah = h[0]; Bh = h[1]; Ch = h[2]; Dh = h[3];
    }
    else if (strcmp(flag, "-lch") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &lch) < 0 || lch <= 0.0) { opserr << "WARNING LadrunoConcrete3D: -lch wants lch > 0\n"; return 0; }
    }
    else if (strcmp(flag, "-autoRegularization") == 0) {
      autoReg = true;
    }
    else {
      opserr << "WARNING LadrunoConcrete3D: unknown option '" << flag << "'\n";
      return 0;
    }
  }

  // ---- foot-gun guards (the kernel ASSUMES valid params; reject at the boundary) ----
  if (E <= 0.0)                      { opserr << "WARNING LadrunoConcrete3D: E must be > 0\n"; return 0; }
  if (nu < 0.0 || nu >= 0.5)         { opserr << "WARNING LadrunoConcrete3D: nu must be in [0,0.5)\n"; return 0; }
  if (fc <= 0.0 || ft <= 0.0)        { opserr << "WARNING LadrunoConcrete3D: fc, ft must be > 0 (POSITIVE magnitudes)\n"; return 0; }
  if (ft >= fc)                      { opserr << "WARNING LadrunoConcrete3D: need ft < fc (m0 = 3(fc^2-ft^2)/(fc ft) e/(e+1) > 0)\n"; return 0; }
  if (Gf <= 0.0 || Gc <= 0.0)        { opserr << "WARNING LadrunoConcrete3D: Gf, Gc must be > 0\n"; return 0; }
  if (As < 1.0)                      { opserr << "WARNING LadrunoConcrete3D: As must be >= 1 (ductility amplitude)\n"; return 0; }
  if (kupfer <= 1.0)                 { opserr << "WARNING LadrunoConcrete3D: -kupfer fcc/fc must be > 1\n"; return 0; }

  if (!haveE)
    ecc = Ladruno::Concrete3D::eccentricityFromKupfer(fc, ft, kupfer);   // (0.5,1] by construction
  if (!(ecc > 0.5 && ecc <= 1.0)) {
    opserr << "WARNING LadrunoConcrete3D: eccentricity e=" << ecc << " must be in (0.5, 1]\n";
    return 0;
  }

  // The CDPM2 consistent tangent is NON-SYMMETRIC (non-associated flow, and ~2% even associated
  // from the semi-implicit theta-freeze) — a symmetric solver will under-converge or diverge.
  opserr << "LadrunoConcrete3D (tag " << tag << "): the consistent tangent is NON-SYMMETRIC; "
         << "use an unsymmetric solver (e.g. `system FullGeneral` or `system UmfPack`).\n";

  NDMaterial* mat = new LadrunoConcrete3D(tag, E, nu, fc, ft, Gf, Gc, ecc, Df, As,
                                          qh0, Hp, Ah, Bh, Ch, Dh, rho, lch, autoReg);
  if (mat == 0) {
    opserr << "WARNING LadrunoConcrete3D: failed to allocate material\n";
    return 0;
  }
  return mat;
}

// ===========================================================================
//  constructors / destructor
// ===========================================================================
LadrunoConcrete3D::LadrunoConcrete3D()
  : NDMaterial(0, ND_TAG_LadrunoConcrete3D),
    E(0.0), nu(0.0), fc(0.0), ft(0.0), Gf(0.0), Gc(0.0), ecc(0.0), m0(0.0),
    Df(0.0), As(2.0), qh0(0.3), Hp(0.5), Ah(0.08), Bh(0.003), Ch(2.0), Dh(1.0e-6),
    rho(0.0), lchFixed(1.0), autoReg(false),
    kp_n(0.0), etmax_n(0.0), kdt1_n(0.0), kdt2_n(0.0), kdc_n(0.0), kdc1_n(0.0), kdc2_n(0.0),
    kp_t(0.0), etmax_t(0.0), kdt1_t(0.0), kdt2_t(0.0), kdc_t(0.0), kdc1_t(0.0), kdc2_t(0.0),
    omegaT(0.0), omegaC(0.0), lastStatus(0),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->revertToStart();
}

LadrunoConcrete3D::LadrunoConcrete3D(int tag, double E_, double nu_, double fc_, double ft_,
                                     double Gf_, double Gc_, double e_, double Df_, double As_,
                                     double qh0_, double Hp_, double Ah_, double Bh_, double Ch_, double Dh_,
                                     double rho_, double lch_, bool autoReg_)
  : NDMaterial(tag, ND_TAG_LadrunoConcrete3D),
    E(E_), nu(nu_), fc(fc_), ft(ft_), Gf(Gf_), Gc(Gc_), ecc(e_),
    m0(Ladruno::Concrete3D::m0Of(fc_, ft_, e_)),
    Df(Df_), As(As_), qh0(qh0_), Hp(Hp_), Ah(Ah_), Bh(Bh_), Ch(Ch_), Dh(Dh_),
    rho(rho_), lchFixed(lch_), autoReg(autoReg_),
    kp_n(0.0), etmax_n(0.0), kdt1_n(0.0), kdt2_n(0.0), kdc_n(0.0), kdc1_n(0.0), kdc2_n(0.0),
    kp_t(0.0), etmax_t(0.0), kdt1_t(0.0), kdt2_t(0.0), kdc_t(0.0), kdc1_t(0.0), kdc2_t(0.0),
    omegaT(0.0), omegaC(0.0), lastStatus(0),
    stressOut(6), strainOut(6), tangentOut(6, 6)
{
  this->revertToStart();
}

LadrunoConcrete3D::~LadrunoConcrete3D() {}

// ===========================================================================
//  parameter pack + the kernel call
// ===========================================================================
void LadrunoConcrete3D::integrate(bool doTangent)
{
  Params p;
  p.E = E; p.nu = nu; p.rho = rho;
  p.fc = fc; p.ft = ft;
  p.e = ecc; p.m0 = m0;
  p.Gf = Gf; p.Gc = Gc;
  p.Df = Df; p.As = As;
  p.qh0 = qh0; p.Hp = Hp; p.Ah = Ah; p.Bh = Bh; p.Ch = Ch; p.Dh = Dh;
  p.eta = 0.0; p.implex = false;                 // Tier-1 implicit (P3 robustness tiers later)

  // crack-band: lch from the active element (mesh-objective) when -autoRegularization is on,
  // else the fixed -lch the input was calibrated for.
  double lch = lchFixed;
  if (autoReg && ops_TheActiveElement != 0) {
    double l = ops_TheActiveElement->getCharacteristicLength();
    if (l > 0.0) lch = l;
  }
  p.lch = lch; p.lch_ref = lch;

  State in;
  for (int i = 0; i < 6; i++) { in.eps[i] = eps_n[i]; in.sig[i] = sig_n[i]; in.sigEff[i] = sigEff_n[i]; }
  in.kp = kp_n; in.et_max = etmax_n;
  in.kdt1 = kdt1_n; in.kdt2 = kdt2_n;
  in.kdc = kdc_n; in.kdc1 = kdc1_n; in.kdc2 = kdc2_n;

  State out;
  double sigEffImpl[6];
  lastStatus = Ladruno::Concrete3D::returnMap(p, strain6, in, out, stress6, sigEffImpl, Dtan6,
                                              doTangent, -1.0, /*hardening=*/true, &omegaT, &omegaC);

  for (int i = 0; i < 6; i++) sigEff6[i] = out.sigEff[i];
  kp_t = out.kp; etmax_t = out.et_max;
  kdt1_t = out.kdt1; kdt2_t = out.kdt2;
  kdc_t = out.kdc; kdc1_t = out.kdc1; kdc2_t = out.kdc2;

  if (lastStatus != 0)
    opserr << "WARNING LadrunoConcrete3D: return map did not converge (tag "
           << this->getTag() << ") — step-cut\n";
}

// ===========================================================================
//  strain interface (element passes ENGINEERING shear; kernel wants TENSOR shear)
// ===========================================================================
int LadrunoConcrete3D::setTrialStrain(const Vector& e)
{
  for (int i = 0; i < 6; i++) strain6[i] = (i < 3) ? e(i) : 0.5 * e(i);  // gamma -> eps_tensor
  this->integrate(true);
  return 0;
}

int LadrunoConcrete3D::setTrialStrain(const Vector& v, const Vector&) { return this->setTrialStrain(v); }

int LadrunoConcrete3D::setTrialStrainIncr(const Vector& v)
{
  Vector ne(6);
  for (int i = 0; i < 6; i++) {
    double cur = (i < 3) ? strain6[i] : 2.0 * strain6[i];   // tensor -> engineering
    ne(i) = cur + v(i);
  }
  return this->setTrialStrain(ne);
}

int LadrunoConcrete3D::setTrialStrainIncr(const Vector& v, const Vector&) { return this->setTrialStrainIncr(v); }

// ===========================================================================
//  responses
// ===========================================================================
const Vector& LadrunoConcrete3D::getStress(void)
{
  for (int i = 0; i < 6; i++) stressOut(i) = stress6[i];   // true tensor stress == engineering stress
  return stressOut;
}

const Vector& LadrunoConcrete3D::getStrain(void)
{
  for (int i = 0; i < 6; i++) strainOut(i) = (i < 3) ? strain6[i] : 2.0 * strain6[i];  // tensor -> gamma
  return strainOut;
}

const Matrix& LadrunoConcrete3D::getTangent(void)
{
  // kernel tangent is in the TENSOR convention (dsig_ij = 2G deps_ij); the element wants
  // d(sigma)/d(engineering strain), so the shear COLUMNS (deps_eng = 2 deps_tensor) are halved.
  for (int a = 0; a < 6; a++)
    for (int b = 0; b < 6; b++)
      tangentOut(a, b) = (b < 3) ? Dtan6[a][b] : 0.5 * Dtan6[a][b];
  return tangentOut;
}

const Matrix& LadrunoConcrete3D::getInitialTangent(void)
{
  Params p;
  p.E = E; p.nu = nu;
  double Cel[6][6];
  Ladruno::Concrete3D::elasticC(p, Cel);
  for (int a = 0; a < 6; a++)
    for (int b = 0; b < 6; b++)
      tangentOut(a, b) = (b < 3) ? Cel[a][b] : 0.5 * Cel[a][b];
  return tangentOut;
}

// ===========================================================================
//  state cycle
// ===========================================================================
int LadrunoConcrete3D::commitState(void)
{
  for (int i = 0; i < 6; i++) { eps_n[i] = strain6[i]; sig_n[i] = stress6[i]; sigEff_n[i] = sigEff6[i]; }
  kp_n = kp_t; etmax_n = etmax_t;
  kdt1_n = kdt1_t; kdt2_n = kdt2_t;
  kdc_n = kdc_t; kdc1_n = kdc1_t; kdc2_n = kdc2_t;
  return 0;
}

int LadrunoConcrete3D::revertToLastCommit(void)
{
  // restore the trial buffers to the committed state (a re-issued setTrialStrain overwrites them)
  for (int i = 0; i < 6; i++) { strain6[i] = eps_n[i]; stress6[i] = sig_n[i]; sigEff6[i] = sigEff_n[i]; }
  kp_t = kp_n; etmax_t = etmax_n;
  kdt1_t = kdt1_n; kdt2_t = kdt2_n;
  kdc_t = kdc_n; kdc1_t = kdc1_n; kdc2_t = kdc2_n;
  omegaT = 0.0; omegaC = 0.0;
  this->getInitialTangent();
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) Dtan6[a][b] = 0.0;
  Params p; p.E = E; p.nu = nu;
  Ladruno::Concrete3D::elasticC(p, Dtan6);
  return 0;
}

int LadrunoConcrete3D::revertToStart(void)
{
  for (int i = 0; i < 6; i++) {
    eps_n[i] = 0.0; sig_n[i] = 0.0; sigEff_n[i] = 0.0;
    strain6[i] = 0.0; stress6[i] = 0.0; sigEff6[i] = 0.0;
  }
  kp_n = 0.0; etmax_n = 0.0; kdt1_n = 0.0; kdt2_n = 0.0; kdc_n = 0.0; kdc1_n = 0.0; kdc2_n = 0.0;
  kp_t = 0.0; etmax_t = 0.0; kdt1_t = 0.0; kdt2_t = 0.0; kdc_t = 0.0; kdc1_t = 0.0; kdc2_t = 0.0;
  omegaT = 0.0; omegaC = 0.0; lastStatus = 0;
  Params p; p.E = E; p.nu = nu;
  Ladruno::Concrete3D::elasticC(p, Dtan6);
  return 0;
}

// ===========================================================================
//  copies
// ===========================================================================
NDMaterial* LadrunoConcrete3D::getCopy(void)
{
  return new LadrunoConcrete3D(this->getTag(), E, nu, fc, ft, Gf, Gc, ecc, Df, As,
                               qh0, Hp, Ah, Bh, Ch, Dh, rho, lchFixed, autoReg);
}

NDMaterial* LadrunoConcrete3D::getCopy(const char* type)
{
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    return this->getCopy();
  // v1 is 3D-only (the finite-strain view is via LogStrain wrapping this 3D material; the
  // PlaneStrain/AxiSymmetric/PlateFiber reduced views are deferred to Phase 2).
  opserr << "LadrunoConcrete3D::getCopy - type '" << type
         << "' not supported in v1 (3D only; reduced views are Phase 2)\n";
  return 0;
}

// ===========================================================================
//  parallel / serialization (flat Vector — the kernel state is all fixed-size scalars)
// ===========================================================================
static const int LC3D_NDATA = 1 + 18 + 1 + 25;   // tag + 18 params + autoReg + 25 committed state

int LadrunoConcrete3D::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(LC3D_NDATA);
  int c = 0;
  data(c++) = this->getTag();
  data(c++) = E; data(c++) = nu; data(c++) = fc; data(c++) = ft; data(c++) = Gf; data(c++) = Gc;
  data(c++) = ecc; data(c++) = m0; data(c++) = Df; data(c++) = As;
  data(c++) = qh0; data(c++) = Hp; data(c++) = Ah; data(c++) = Bh; data(c++) = Ch; data(c++) = Dh;
  data(c++) = rho; data(c++) = lchFixed;
  data(c++) = autoReg ? 1.0 : 0.0;
  for (int i = 0; i < 6; i++) data(c++) = eps_n[i];
  for (int i = 0; i < 6; i++) data(c++) = sig_n[i];
  for (int i = 0; i < 6; i++) data(c++) = sigEff_n[i];
  data(c++) = kp_n; data(c++) = etmax_n;
  data(c++) = kdt1_n; data(c++) = kdt2_n;
  data(c++) = kdc_n; data(c++) = kdc1_n; data(c++) = kdc2_n;

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoConcrete3D::sendSelf - failed to send vector\n";
    return -1;
  }
  return 0;
}

int LadrunoConcrete3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
  static Vector data(LC3D_NDATA);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "LadrunoConcrete3D::recvSelf - failed to recv vector\n";
    return -1;
  }
  int c = 0;
  this->setTag((int)data(c++));
  E = data(c++); nu = data(c++); fc = data(c++); ft = data(c++); Gf = data(c++); Gc = data(c++);
  ecc = data(c++); m0 = data(c++); Df = data(c++); As = data(c++);
  qh0 = data(c++); Hp = data(c++); Ah = data(c++); Bh = data(c++); Ch = data(c++); Dh = data(c++);
  rho = data(c++); lchFixed = data(c++);
  autoReg = (data(c++) != 0.0);
  for (int i = 0; i < 6; i++) eps_n[i] = data(c++);
  for (int i = 0; i < 6; i++) sig_n[i] = data(c++);
  for (int i = 0; i < 6; i++) sigEff_n[i] = data(c++);
  kp_n = data(c++); etmax_n = data(c++);
  kdt1_n = data(c++); kdt2_n = data(c++);
  kdc_n = data(c++); kdc1_n = data(c++); kdc2_n = data(c++);

  this->revertToLastCommit();   // sync trial buffers + tangent to the received committed state
  return 0;
}

// ===========================================================================
//  print
// ===========================================================================
void LadrunoConcrete3D::Print(OPS_Stream& s, int)
{
  s << endln;
  s << "LadrunoConcrete3D (CDPM2-grade solid concrete, 3D)" << endln;
  s << "  tag      : " << this->getTag() << endln;
  s << "  E, nu    : " << E << ", " << nu << endln;
  s << "  fc, ft   : " << fc << ", " << ft << "   (e=" << ecc << ", m0=" << m0 << ")" << endln;
  s << "  Gf, Gc   : " << Gf << ", " << Gc << "   (lch=" << lchFixed
    << (autoReg ? ", auto-regularized" : "") << ")" << endln;
  s << "  Df, As   : " << Df << ", " << As << endln;
  s << "  harden   : qh0=" << qh0 << " Hp=" << Hp
    << "   ductility Ah=" << Ah << " Bh=" << Bh << " Ch=" << Ch << " Dh=" << Dh << endln;
  s << "  damage   : omega_t=" << omegaT << " omega_c=" << omegaC << "   kappa_p=" << kp_t << endln;
  s << "  rho      : " << rho << endln;
}

// ===========================================================================
//  recordable responses
// ===========================================================================
Response* LadrunoConcrete3D::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return NDMaterial::setResponse(argv, argc, s);
  const char* a = argv[0];

  if (strcmp(a, "stress") == 0 || strcmp(a, "stresses") == 0)
    return new MaterialResponse(this, 1, this->getStress());
  if (strcmp(a, "strain") == 0 || strcmp(a, "strains") == 0)
    return new MaterialResponse(this, 2, this->getStrain());
  if (strcmp(a, "tangent") == 0 || strcmp(a, "Tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  if (strcmp(a, "effectiveStress") == 0 || strcmp(a, "sigEff") == 0)
    return new MaterialResponse(this, 4, Vector(6));
  if (strcmp(a, "damage") == 0 || strcmp(a, "omega") == 0)
    return new MaterialResponse(this, 5, Vector(2));
  if (strcmp(a, "kappaP") == 0 || strcmp(a, "kappa_p") == 0 || strcmp(a, "kappaPlastic") == 0)
    return new MaterialResponse(this, 6, Vector(1));
  if (strcmp(a, "plasticStrain") == 0 || strcmp(a, "plasticStrains") == 0)
    return new MaterialResponse(this, 7, Vector(6));

  return NDMaterial::setResponse(argv, argc, s);
}

int LadrunoConcrete3D::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
    case 1:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStress();
      return 0;
    case 2:
      if (matInfo.theVector) *(matInfo.theVector) = this->getStrain();
      return 0;
    case 3:
      if (matInfo.theMatrix) *(matInfo.theMatrix) = this->getTangent();
      return 0;
    case 4:                                   // effective (undamaged) stress, 6 tensor comps
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        for (int i = 0; i < 6; i++) v(i) = sigEff6[i];
      }
      return 0;
    case 5:                                   // dual damage [omega_t, omega_c]
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        v(0) = omegaT; v(1) = omegaC;
      }
      return 0;
    case 6:                                   // hardening variable kappa_p
      if (matInfo.theVector) (*(matInfo.theVector))(0) = kp_t;
      return 0;
    case 7:                                   // plastic strain eps - C^-1 sigEff (engineering shear)
      if (matInfo.theVector) {
        Vector& v = *(matInfo.theVector);
        double epl[6];
        Params p; p.E = E; p.nu = nu;
        Ladruno::Concrete3D::plasticStrain6(sigEff6, strain6, p, epl);
        for (int i = 0; i < 6; i++) v(i) = (i < 3) ? epl[i] : 2.0 * epl[i];
      }
      return 0;
    default:
      return -1;
  }
}
