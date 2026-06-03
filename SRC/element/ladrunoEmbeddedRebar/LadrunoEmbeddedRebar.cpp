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

// Ladruno: embedded-reinforcement coupling element. See LadrunoEmbeddedRebar.h
// and Ladruno_implementation/20_ladruno_embedded_reinforcement_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedRebar.h>
#include <Domain.h>
#include <Node.h>
#include <UniaxialMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>

// ===========================================================================
//  construction
// ===========================================================================
LadrunoEmbeddedRebar::LadrunoEmbeddedRebar(int tag, int ndm_, int rebarNode,
        const ID& hostNodes, const Vector& shape, const Vector& dir_,
        double kt_, double bondScale_, UniaxialMaterial* bm, double kAxialPerfect_,
        int hostEleTag_, bool ktAuto_, double ktAlpha_)
  : Element(tag, ELE_TAG_LadrunoEmbeddedRebar),
    ndm(ndm_), nHost(hostNodes.Size()), nDOF((1 + hostNodes.Size()) * ndm_),
    connectedNodes(1 + hostNodes.Size()), Nshape(shape), dir(ndm_),
    kt(kt_), bondScale(bondScale_), kAxialPerfect(kAxialPerfect_),
    bondMat(0),
    hostEleTag(hostEleTag_), ktAuto(ktAuto_), ktAlpha(ktAlpha_), ktResolved(false),
    theNodes(0), K(0), P(0), M0(0), bondEnergyResp(0)
{
  connectedNodes(0) = rebarNode;
  for (int i = 0; i < nHost; i++)
    connectedNodes(1 + i) = hostNodes(i);

  // normalize the bar axis
  double nrm = 0.0;
  for (int k = 0; k < ndm; k++) nrm += dir_(k) * dir_(k);
  nrm = (nrm > 0.0) ? sqrt(nrm) : 1.0;
  for (int k = 0; k < ndm; k++) dir(k) = dir_(k) / nrm;

  if (bm != 0)
    bondMat = bm->getCopy();

  theNodes = new Node*[1 + nHost];
  for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
  K  = new Matrix(nDOF, nDOF);
  P  = new Vector(nDOF);
  M0 = new Matrix(nDOF, nDOF);   // stays zero
}

LadrunoEmbeddedRebar::LadrunoEmbeddedRebar()
  : Element(0, ELE_TAG_LadrunoEmbeddedRebar),
    ndm(0), nHost(0), nDOF(0), connectedNodes(0), Nshape(), dir(),
    kt(0.0), bondScale(1.0), kAxialPerfect(0.0),
    bondMat(0),
    hostEleTag(-1), ktAuto(false), ktAlpha(0.0), ktResolved(false),
    theNodes(0), K(0), P(0), M0(0), bondEnergyResp(0)
{
}

LadrunoEmbeddedRebar::~LadrunoEmbeddedRebar()
{
  if (bondMat != 0) delete bondMat;
  if (theNodes != 0) delete[] theNodes;
  if (K != 0) delete K;
  if (P != 0) delete P;
  if (M0 != 0) delete M0;
  if (bondEnergyResp != 0) delete bondEnergyResp;
}

// ===========================================================================
//  domain
// ===========================================================================
int LadrunoEmbeddedRebar::getNumExternalNodes(void) const { return 1 + nHost; }
const ID& LadrunoEmbeddedRebar::getExternalNodes(void) { return connectedNodes; }
Node** LadrunoEmbeddedRebar::getNodePtrs(void) { return theNodes; }
int LadrunoEmbeddedRebar::getNumDOF(void) { return nDOF; }

void LadrunoEmbeddedRebar::setDomain(Domain* theDomain)
{
  if (theDomain == 0) {
    for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
    return;
  }
  for (int i = 0; i < 1 + nHost; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "LadrunoEmbeddedRebar::setDomain - node " << connectedNodes(i)
             << " not found\n";
      return;
    }
    if (theNodes[i]->getNumberDOF() < ndm) {
      opserr << "LadrunoEmbeddedRebar::setDomain - node " << connectedNodes(i)
             << " has fewer than ndm=" << ndm << " DOFs\n";
      return;
    }
  }
  this->DomainComponent::setDomain(theDomain);
}

// ADR 20 §10.2a — resolve the auto transverse penalty from the host element's
// own initial stiffness:  kt = ktAlpha * max_i |K_host(i,i)|  (~ ktAlpha*E*lch).
// Done lazily (first assembly) so the host has already setDomain'd — element
// setDomain order is not guaranteed, so we must NOT touch the host's K here.
void LadrunoEmbeddedRebar::resolveAutoKt(void)
{
  if (!ktAuto || ktResolved) return;
  Domain* theDomain = this->getDomain();
  if (theDomain == 0) return;                 // not in a domain yet; retry later
  if (hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedRebar " << this->getTag()
           << ": -kt auto needs the -host form; keeping kt=" << kt << "\n";
    ktResolved = true;
    return;
  }
  Element* host = theDomain->getElement(hostEleTag);
  if (host == 0) {
    opserr << "WARNING LadrunoEmbeddedRebar " << this->getTag()
           << ": -kt auto host element " << hostEleTag
           << " not found; keeping kt=" << kt << "\n";
    ktResolved = true;
    return;
  }
  const Matrix& Kh = host->getInitialStiff();
  int nh = Kh.noRows();
  double scale = 0.0;                          // max abs diagonal ~ E_host*lch
  for (int i = 0; i < nh; i++) {
    double d = fabs(Kh(i, i));
    if (d > scale) scale = d;
  }
  if (scale <= 0.0) {
    opserr << "WARNING LadrunoEmbeddedRebar " << this->getTag()
           << ": -kt auto host " << hostEleTag
           << " gave a zero stiffness scale; keeping kt=" << kt << "\n";
    ktResolved = true;
    return;
  }
  kt = ktAlpha * scale;
  ktResolved = true;
}

// ===========================================================================
//  state
// ===========================================================================
int LadrunoEmbeddedRebar::commitState(void)
{
  int ok = this->Element::commitState();
  if (bondMat != 0) ok += bondMat->commitState();
  return ok;
}

int LadrunoEmbeddedRebar::revertToLastCommit(void)
{
  return (bondMat != 0) ? bondMat->revertToLastCommit() : 0;
}

int LadrunoEmbeddedRebar::revertToStart(void)
{
  return (bondMat != 0) ? bondMat->revertToStart() : 0;
}

// gap g, axial slip s, transverse gt, and the axial force/tangent
void LadrunoEmbeddedRebar::formBandTraction(Vector& g, double& s, Vector& gt,
                                            double& Faxial, double& kAxial)
{
  g.resize(ndm); g.Zero();
  const Vector& uR = theNodes[0]->getTrialDisp();
  for (int k = 0; k < ndm; k++) g(k) = uR(k);
  for (int i = 0; i < nHost; i++) {
    const Vector& uH = theNodes[1 + i]->getTrialDisp();
    for (int k = 0; k < ndm; k++) g(k) -= Nshape(i) * uH(k);
  }
  s = 0.0;
  for (int k = 0; k < ndm; k++) s += g(k) * dir(k);
  gt.resize(ndm);
  for (int k = 0; k < ndm; k++) gt(k) = g(k) - s * dir(k);

  if (bondMat != 0) {
    bondMat->setTrialStrain(s);
    Faxial = bondMat->getStress()  * bondScale;
    kAxial = bondMat->getTangent() * bondScale;
  } else {
    Faxial = kAxialPerfect * s;
    kAxial = kAxialPerfect;
  }
}

int LadrunoEmbeddedRebar::update(void)
{
  this->resolveAutoKt();
  Vector g(ndm), gt(ndm);
  double s, Faxial, kAxial;
  this->formBandTraction(g, s, gt, Faxial, kAxial);
  return 0;
}

// ===========================================================================
//  matrices / forces:  P = B^T t,  K = B^T D B,  coeff c = [1, -N_1, ..., -N_M]
// ===========================================================================
const Vector& LadrunoEmbeddedRebar::getResistingForce(void)
{
  this->resolveAutoKt();
  P->Zero();
  Vector g(ndm), gt(ndm);
  double s, Faxial, kAxial;
  this->formBandTraction(g, s, gt, Faxial, kAxial);

  // traction t = Faxial*dir + kt*gt
  Vector t(ndm);
  for (int k = 0; k < ndm; k++) t(k) = Faxial * dir(k) + kt * gt(k);

  // block 0 (rebar): +t ; block (1+i): -N_i t
  for (int k = 0; k < ndm; k++) (*P)(k) = t(k);
  for (int i = 0; i < nHost; i++) {
    double c = -Nshape(i);
    for (int k = 0; k < ndm; k++) (*P)((1 + i) * ndm + k) = c * t(k);
  }
  return *P;
}

const Vector& LadrunoEmbeddedRebar::getResistingForceIncInertia(void)
{
  return this->getResistingForce();   // no element mass / damping
}

const Matrix& LadrunoEmbeddedRebar::getTangentStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  Vector g(ndm), gt(ndm);
  double s, Faxial, kAxial;
  this->formBandTraction(g, s, gt, Faxial, kAxial);

  // D = kAxial dir(x)dir + kt (I - dir(x)dir)
  Matrix D(ndm, ndm);
  for (int a = 0; a < ndm; a++)
    for (int b = 0; b < ndm; b++) {
      double dd = dir(a) * dir(b);
      D(a, b) = kAxial * dd + kt * ((a == b ? 1.0 : 0.0) - dd);
    }

  // c = [1, -N_1, ..., -N_M]; K block(p,q) = c_p c_q D
  Vector c(1 + nHost);
  c(0) = 1.0;
  for (int i = 0; i < nHost; i++) c(1 + i) = -Nshape(i);

  for (int p = 0; p < 1 + nHost; p++)
    for (int q = 0; q < 1 + nHost; q++) {
      double cc = c(p) * c(q);
      if (cc == 0.0) continue;
      for (int a = 0; a < ndm; a++)
        for (int b = 0; b < ndm; b++)
          (*K)(p * ndm + a, q * ndm + b) += cc * D(a, b);
    }
  return *K;
}

const Matrix& LadrunoEmbeddedRebar::getInitialStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  double kAxial0 = (bondMat != 0) ? bondMat->getInitialTangent() * bondScale
                                  : kAxialPerfect;
  Matrix D(ndm, ndm);
  for (int a = 0; a < ndm; a++)
    for (int b = 0; b < ndm; b++) {
      double dd = dir(a) * dir(b);
      D(a, b) = kAxial0 * dd + kt * ((a == b ? 1.0 : 0.0) - dd);
    }
  Vector c(1 + nHost);
  c(0) = 1.0;
  for (int i = 0; i < nHost; i++) c(1 + i) = -Nshape(i);
  for (int p = 0; p < 1 + nHost; p++)
    for (int q = 0; q < 1 + nHost; q++) {
      double cc = c(p) * c(q);
      if (cc == 0.0) continue;
      for (int a = 0; a < ndm; a++)
        for (int b = 0; b < ndm; b++)
          (*K)(p * ndm + a, q * ndm + b) += cc * D(a, b);
    }
  return *K;
}

const Matrix& LadrunoEmbeddedRebar::getMass(void)
{
  M0->Zero();
  return *M0;
}

// ===========================================================================
//  serialization
// ===========================================================================
int LadrunoEmbeddedRebar::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  // header: tag, ndm, nHost, kt, bondScale, kAxialPerfect, hasBond, bondClassTag,
  //         bondDbTag, hostEleTag, ktAuto, ktAlpha
  static Vector hdr(12);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nHost;
  hdr(3) = kt;
  hdr(4) = bondScale;
  hdr(5) = kAxialPerfect;
  hdr(6) = (bondMat != 0) ? 1.0 : 0.0;
  int bondClassTag = 0, bondDbTag = 0;
  if (bondMat != 0) {
    bondClassTag = bondMat->getClassTag();
    bondDbTag = bondMat->getDbTag();
    if (bondDbTag == 0) {
      bondDbTag = theChannel.getDbTag();
      bondMat->setDbTag(bondDbTag);
    }
  }
  hdr(7) = bondClassTag;
  hdr(8) = bondDbTag;
  hdr(9) = hostEleTag;
  hdr(10) = ktAuto ? 1.0 : 0.0;
  hdr(11) = ktAlpha;
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoEmbeddedRebar::sendSelf - header failed\n";
    return -1;
  }
  // payload: connectedNodes (1+nHost ints) via an ID, then Nshape + dir via a Vector
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedRebar::sendSelf - ID failed\n";
    return -1;
  }
  Vector payload(nHost + ndm);
  for (int i = 0; i < nHost; i++) payload(i) = Nshape(i);
  for (int k = 0; k < ndm; k++) payload(nHost + k) = dir(k);
  if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedRebar::sendSelf - payload failed\n";
    return -1;
  }
  if (bondMat != 0 && bondMat->sendSelf(commitTag, theChannel) < 0) {
    opserr << "LadrunoEmbeddedRebar::sendSelf - bondMat failed\n";
    return -1;
  }
  return 0;
}

int LadrunoEmbeddedRebar::recvSelf(int commitTag, Channel& theChannel,
                                   FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  static Vector hdr(12);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoEmbeddedRebar::recvSelf - header failed\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nHost = (int)hdr(2);
  kt = hdr(3);
  bondScale = hdr(4);
  kAxialPerfect = hdr(5);
  bool hasBond = (hdr(6) != 0.0);
  int bondClassTag = (int)hdr(7);
  int bondDbTag = (int)hdr(8);
  hostEleTag = (int)hdr(9);
  ktAuto = (hdr(10) != 0.0);
  ktAlpha = hdr(11);
  ktResolved = false;   // re-resolve against the host on first assembly

  nDOF = (1 + nHost) * ndm;
  connectedNodes.resize(1 + nHost);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedRebar::recvSelf - ID failed\n";
    return -1;
  }
  Vector payload(nHost + ndm);
  if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedRebar::recvSelf - payload failed\n";
    return -1;
  }
  Nshape.resize(nHost);
  for (int i = 0; i < nHost; i++) Nshape(i) = payload(i);
  dir.resize(ndm);
  for (int k = 0; k < ndm; k++) dir(k) = payload(nHost + k);

  if (bondEnergyResp != 0) { delete bondEnergyResp; bondEnergyResp = 0; }
  if (bondMat != 0) { delete bondMat; bondMat = 0; }
  if (hasBond) {
    bondMat = theBroker.getNewUniaxialMaterial(bondClassTag);
    if (bondMat == 0) {
      opserr << "LadrunoEmbeddedRebar::recvSelf - broker no mat " << bondClassTag << "\n";
      return -1;
    }
    bondMat->setDbTag(bondDbTag);
    if (bondMat->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "LadrunoEmbeddedRebar::recvSelf - bondMat recv failed\n";
      return -1;
    }
  }

  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node*[1 + nHost];
  for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
  if (K != 0) delete K;  K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;  P = new Vector(nDOF);
  if (M0 != 0) delete M0; M0 = new Matrix(nDOF, nDOF);
  return 0;
}

void LadrunoEmbeddedRebar::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoEmbeddedRebar, tag: " << this->getTag() << "\n";
  s << "  rebar node: " << connectedNodes(0) << "  host nodes: ";
  for (int i = 0; i < nHost; i++) s << connectedNodes(1 + i) << " ";
  s << "\n  mode: " << (bondMat != 0 ? "bond-slip" : "perfect(penalty)")
    << "  kt: " << kt << "  bondScale: " << bondScale << "\n";
}

// ===========================================================================
//  responses: slip / bondStress / force
// ===========================================================================
Response* LadrunoEmbeddedRebar::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  if (strcmp(argv[0], "slip") == 0)
    return new ElementResponse(this, 1, 0.0);
  if (strcmp(argv[0], "bondStress") == 0)
    return new ElementResponse(this, 2, 0.0);
  if (strcmp(argv[0], "axialForce") == 0)
    return new ElementResponse(this, 3, 0.0);
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "localForce") == 0)
    return new ElementResponse(this, 4, Vector(ndm));
  // ADR 20 §10.2b diagnostics — the ARTIFICIAL penalty energy currently stored
  // (transverse 1/2 kt|gt|^2, plus the perfect-bond axial 1/2 k s^2 when there
  // is no bond law), the transverse constraint violation |gt|, and the resolved
  // transverse penalty (useful with -kt auto). These let the user net the
  // artificial penalty energy out of a global EnergyBalance and audit the tie.
  if (strcmp(argv[0], "penaltyEnergy") == 0)
    return new ElementResponse(this, 5, 0.0);
  if (strcmp(argv[0], "constraintViolation") == 0)
    return new ElementResponse(this, 6, 0.0);
  if (strcmp(argv[0], "kt") == 0 || strcmp(argv[0], "penalty") == 0)
    return new ElementResponse(this, 7, 0.0);
  // physical bond energy = perimeter*L_trib * (cumulative material bond work),
  // single-sourced from the bond law's own "energy" response (bond-law-agnostic;
  // 0 for perfect bond, where the axial energy is artificial and lives in
  // penaltyEnergy instead). ADR 20 §10.2b.
  if (strcmp(argv[0], "bondEnergy") == 0 || strcmp(argv[0], "bondDissipation") == 0) {
    if (bondEnergyResp != 0) { delete bondEnergyResp; bondEnergyResp = 0; }
    if (bondMat != 0) {
      const char* a[1] = {"energy"};
      bondEnergyResp = bondMat->setResponse(a, 1, s);
    }
    return new ElementResponse(this, 8, 0.0);
  }
  return 0;
}

int LadrunoEmbeddedRebar::getResponse(int responseID, Information& eleInfo)
{
  this->resolveAutoKt();
  Vector g(ndm), gt(ndm);
  double s, Faxial, kAxial;
  this->formBandTraction(g, s, gt, Faxial, kAxial);
  switch (responseID) {
  case 1: return eleInfo.setDouble(s);
  case 2: return eleInfo.setDouble((bondMat != 0) ? bondMat->getStress()
                                                  : kAxialPerfect * s);
  case 3: return eleInfo.setDouble(Faxial);
  case 4: {
    Vector t(ndm);
    for (int k = 0; k < ndm; k++) t(k) = Faxial * dir(k) + kt * gt(k);
    return eleInfo.setVector(t);
  }
  case 5: {
    // artificial penalty energy currently stored (recoverable). The transverse
    // hold is always artificial; the axial slot is artificial ONLY for perfect
    // bond (with a bond law the axial tau-s energy is physical, owned by bondMat).
    double gt2 = 0.0;
    for (int k = 0; k < ndm; k++) gt2 += gt(k) * gt(k);
    double E = 0.5 * kt * gt2;
    if (bondMat == 0) E += 0.5 * kAxialPerfect * s * s;
    return eleInfo.setDouble(E);
  }
  case 6: {
    double gt2 = 0.0;
    for (int k = 0; k < ndm; k++) gt2 += gt(k) * gt(k);
    return eleInfo.setDouble(sqrt(gt2));
  }
  case 7: return eleInfo.setDouble(kt);
  case 8: {
    double w = 0.0;
    if (bondMat != 0 && bondEnergyResp != 0) {
      if (bondEnergyResp->getResponse() == 0)
        w = bondEnergyResp->getInformation().theDouble;
    }
    return eleInfo.setDouble(bondScale * w);
  }
  default: return -1;
  }
}
