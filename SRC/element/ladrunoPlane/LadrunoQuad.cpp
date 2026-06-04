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

// LadrunoQuad — unified 4-node bilinear plane continuum element (Ladruno fork).
// See LadrunoQuad.h and Ladruno_implementation/25_ladruno_plane_elements_adr.md
//
// The std formulation mirrors upstream FourNodeQuad (same 2x2 Gauss rule and
// bilinear shape functions) so it reduces to it to ~1e-9. The bbar formulation
// is mean-dilatation B-bar with the 2D dilatation factor 1/2 (NOT the brick's
// 1/3) — PlaneStrain only.

#include <LadrunoQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Response.h>
#include <DummyStream.h>
#include <elementAPI.h>
#include <math.h>

// Tier-A damage-scaled hourglass: floor on the stabilization multiplier so a
// fully-damaged element keeps a little hourglass control (mirrors LadrunoBrick).
static const double HG_DAMAGE_FLOOR = 0.01;

// 2x2 dyadic product helper (SSP base-vector second-moment tensor).
static inline void dyad2(const Vector &a, const Vector &b, Matrix &out)
{
  out(0, 0) = a(0) * b(0); out(0, 1) = a(0) * b(1);
  out(1, 0) = a(1) * b(0); out(1, 1) = a(1) * b(1);
}

double LadrunoQuad::matrixData[64];
Matrix LadrunoQuad::K(matrixData, 8, 8);
Vector LadrunoQuad::P(8);
double LadrunoQuad::shp[3][4];
double LadrunoQuad::shpBar[2][4];
double LadrunoQuad::pts[4][2];
double LadrunoQuad::wts[4];

LadrunoQuad::LadrunoQuad(int tag, int nd1, int nd2, int nd3, int nd4,
                         NDMaterial &m, const char *type, double t,
                         Formulation form, double r, double b1, double b2,
                         double p)
 :Element(tag, ELE_TAG_LadrunoQuad),
  theMaterial(0), connectedExternalNodes(4),
  Q(8), pressureLoad(8), thickness(t), pressure(p), rho(r),
  formulation(form), planeType(1),
  Mmem(3, 8), Kstab(8, 8), J0(0.0), J1(0.0), J2(0.0), damageResponse(0), Ki(0)
{
  pts[0][0] = -0.5773502691896258; pts[0][1] = -0.5773502691896258;
  pts[1][0] =  0.5773502691896258; pts[1][1] = -0.5773502691896258;
  pts[2][0] =  0.5773502691896258; pts[2][1] =  0.5773502691896258;
  pts[3][0] = -0.5773502691896258; pts[3][1] =  0.5773502691896258;
  wts[0] = wts[1] = wts[2] = wts[3] = 1.0;

  if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)
    planeType = 1;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0)
    planeType = 2;
  else {
    opserr << "LadrunoQuad::LadrunoQuad -- improper material type: " << type << "\n";
    exit(-1);
  }

  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b1;
  b[1] = b2;

  theMaterial = new NDMaterial *[4];
  for (int i = 0; i < 4; i++) {
    theMaterial[i] = m.getCopy(type);
    if (theMaterial[i] == 0) {
      opserr << "LadrunoQuad::LadrunoQuad -- failed to get a copy of material model\n";
      exit(-1);
    }
  }

  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;
  connectedExternalNodes(3) = nd4;

  for (int i = 0; i < 4; i++)
    theNodes[i] = 0;
}

LadrunoQuad::LadrunoQuad()
 :Element(0, ELE_TAG_LadrunoQuad),
  theMaterial(0), connectedExternalNodes(4),
  Q(8), pressureLoad(8), thickness(0.0), pressure(0.0), rho(0.0),
  formulation(Formulation::STD), planeType(1),
  Mmem(3, 8), Kstab(8, 8), J0(0.0), J1(0.0), J2(0.0), damageResponse(0), Ki(0)
{
  pts[0][0] = -0.5773502691896258; pts[0][1] = -0.5773502691896258;
  pts[1][0] =  0.5773502691896258; pts[1][1] = -0.5773502691896258;
  pts[2][0] =  0.5773502691896258; pts[2][1] =  0.5773502691896258;
  pts[3][0] = -0.5773502691896258; pts[3][1] =  0.5773502691896258;
  wts[0] = wts[1] = wts[2] = wts[3] = 1.0;
  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b[1] = 0.0;
  for (int i = 0; i < 4; i++)
    theNodes[i] = 0;
}

LadrunoQuad::~LadrunoQuad()
{
  if (theMaterial) {
    for (int i = 0; i < 4; i++)
      if (theMaterial[i]) delete theMaterial[i];
    delete[] theMaterial;
  }
  if (damageResponse) delete damageResponse;
  if (Ki) delete Ki;
}

const char *LadrunoQuad::typeString(void) const
{
  return (planeType == 2) ? "PlaneStress" : "PlaneStrain";
}

int LadrunoQuad::getNumExternalNodes(void) const { return 4; }
const ID &LadrunoQuad::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoQuad::getNodePtrs(void) { return theNodes; }
int LadrunoQuad::getNumDOF(void) { return 8; }

void LadrunoQuad::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    for (int i = 0; i < 4; i++) theNodes[i] = 0;
    return;
  }

  for (int i = 0; i < 4; i++)
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));

  for (int i = 0; i < 4; i++)
    if (theNodes[i] == 0) return;

  for (int i = 0; i < 4; i++)
    if (theNodes[i]->getNumberDOF() != 2) {
      opserr << "LadrunoQuad::setDomain -- node " << connectedExternalNodes(i)
             << " must have 2 DOFs\n";
      return;
    }

  // SSP: build the (geometry-fixed) stabilization once; rebuilt on the receive
  // side here too, so sendSelf carries nothing extra.
  if (formulation == Formulation::SSP)
    this->computeSSP();

  // Tier-A: cache the centroid material's "damage" probe so the constant elastic
  // Kstab can degrade under softening (ssp only). Materials with no "damage"
  // channel return 0 -> damageScale() falls back to 1.0 (elastic/J2 unchanged).
  if (damageResponse) { delete damageResponse; damageResponse = 0; }
  if (formulation == Formulation::SSP && theMaterial[0] != 0) {
    DummyStream dmgStream;
    const char *dmgArgv[1] = {"damage"};
    damageResponse = theMaterial[0]->setResponse(dmgArgv, 1, dmgStream);
  }

  this->DomainComponent::setDomain(theDomain);
  this->setPressureLoadAtNodes();
}

double LadrunoQuad::getCharacteristicLength(void)
{
  // crack-band size = sqrt(area); area = integral of detJ over the element.
  double A = 0.0;
  for (int i = 0; i < 4; i++)
    A += this->shapeFunction(pts[i][0], pts[i][1]) * wts[i];
  if (A <= 0.0)
    return this->Element::getCharacteristicLength();
  return sqrt(A);
}

double LadrunoQuad::damageScale(void)
{
  if (damageResponse == 0)
    return 1.0;
  if (damageResponse->getResponse() < 0)
    return 1.0;
  Information &info = damageResponse->getInformation();
  const Vector *d = info.theVector;
  if (d == 0)
    return 1.0;
  double dmax = 0.0;
  for (int i = 0; i < d->Size(); i++)
    if ((*d)(i) > dmax) dmax = (*d)(i);   // max(d_tension, d_compression)
  double s = 1.0 - dmax;
  if (s < HG_DAMAGE_FLOOR) s = HG_DAMAGE_FLOOR;
  if (s > 1.0) s = 1.0;
  return s;
}

// Port of UWelements/SSPquad::GetStab — single-point membrane operator Mmem
// (3x8) at the centroid + statically-condensed hourglass stabilization Kstab
// (8x8) from the INITIAL material tangent. Geometry-fixed (small strain) so it
// is computed once in setDomain.
void LadrunoQuad::computeSSP(void)
{
  // node coordinate matrix (2x4)
  static Matrix nodeCrd(2, 4);
  for (int i = 0; i < 4; i++) {
    const Vector &X = theNodes[i]->getCrds();
    nodeCrd(0, i) = X(0);
    nodeCrd(1, i) = X(1);
  }

  // jacobian terms (SSPquad convention)
  J0 = ((nodeCrd(0,1)-nodeCrd(0,3))*(nodeCrd(1,2)-nodeCrd(1,0))+(nodeCrd(0,2)-nodeCrd(0,0))*(nodeCrd(1,3)-nodeCrd(1,1)))/8.0;
  J1 = ((nodeCrd(0,1)-nodeCrd(0,0))*(nodeCrd(1,2)-nodeCrd(1,3))+(nodeCrd(0,2)-nodeCrd(0,3))*(nodeCrd(1,0)-nodeCrd(1,1)))/24.0;
  J2 = ((nodeCrd(0,0)-nodeCrd(0,3))*(nodeCrd(1,2)-nodeCrd(1,1))+(nodeCrd(0,2)-nodeCrd(0,1))*(nodeCrd(1,3)-nodeCrd(1,0)))/24.0;

  // shape-function derivatives (local crd) at center
  static Matrix dNloc(4, 2);
  dNloc(0,0) = -0.25; dNloc(1,0) =  0.25; dNloc(2,0) =  0.25; dNloc(3,0) = -0.25;
  dNloc(0,1) = -0.25; dNloc(1,1) = -0.25; dNloc(2,1) =  0.25; dNloc(3,1) =  0.25;

  static Matrix Jmat(2, 2), Jinv(2, 2), dN(4, 2);
  Jmat = nodeCrd * dNloc;     // 2x2
  Jmat.Invert(Jinv);
  dN = dNloc * Jinv;          // 4x2 global gradients

  // hourglass stabilization vector gamma = 0.25*(h - (h.x)bx - (h.y)by)
  double hx = nodeCrd(0,0) - nodeCrd(0,1) + nodeCrd(0,2) - nodeCrd(0,3);
  double hy = nodeCrd(1,0) - nodeCrd(1,1) + nodeCrd(1,2) - nodeCrd(1,3);
  double gamma[4];
  gamma[0] = 0.25*( 1.0 - hx*dN(0,0) - hy*dN(0,1));
  gamma[1] = 0.25*(-1.0 - hx*dN(1,0) - hy*dN(1,1));
  gamma[2] = 0.25*( 1.0 - hx*dN(2,0) - hy*dN(2,1));
  gamma[3] = 0.25*(-1.0 - hx*dN(3,0) - hy*dN(3,1));

  Mmem.Zero();
  static Matrix Mben(2, 8);
  Mben.Zero();
  for (int i = 0; i < 4; i++) {
    Mmem(0, 2*i)   = dN(i,0);
    Mmem(1, 2*i+1) = dN(i,1);
    Mmem(2, 2*i)   = dN(i,1);
    Mmem(2, 2*i+1) = dN(i,0);
    Mben(0, 2*i)   = gamma[i];
    Mben(1, 2*i+1) = gamma[i];
  }

  // normalized base vectors
  static Vector g1(2), g2(2);
  g1(0) = Jmat(0,0); g1(1) = Jmat(1,0);
  g2(0) = Jmat(0,1); g2(1) = Jmat(1,1);
  g1.Normalize();
  g2.Normalize();

  // second moment of area tensor I = (4/3) t J0 (g1 (x) g1 + g2 (x) g2)
  static Matrix I(2, 2), tmp(2, 2);
  dyad2(g1, g1, I);
  dyad2(g2, g2, tmp);
  I += tmp;
  I *= (4.0 / 3.0) * thickness * J0;

  double Hss = (I(0,0)*Jinv(1,0)*Jinv(1,0) + I(0,1)*Jinv(0,0)*Jinv(1,0) + I(1,1)*Jinv(0,0)*Jinv(0,0))*0.25;
  double Htt = (I(0,0)*Jinv(1,1)*Jinv(1,1) + I(0,1)*Jinv(0,1)*Jinv(1,1) + I(1,1)*Jinv(0,1)*Jinv(0,1))*0.25;
  double Hst = (I(0,0)*Jinv(1,1)*Jinv(1,0) + I(0,1)*(Jinv(1,0)*Jinv(0,1) + Jinv(1,1)*Jinv(0,0)) + I(1,1)*Jinv(0,1)*Jinv(0,0))*0.25;

  const Matrix &C = theMaterial[0]->getInitialTangent();
  static Matrix FCF(2, 2);
  FCF(0,0) = (C(0,0) - (C(0,1) + C(1,0)) + C(1,1))*Hss;
  FCF(0,1) = (C(0,1) - (C(0,0) + C(1,1)) + C(1,0))*Hst;
  FCF(1,0) = (C(1,0) - (C(0,0) + C(1,1)) + C(0,1))*Hst;
  FCF(1,1) = (C(1,1) - (C(0,1) + C(1,0)) + C(0,0))*Htt;

  Kstab.Zero();
  Kstab.addMatrixTripleProduct(1.0, Mben, FCF, 1.0);
}

int LadrunoQuad::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0)
    opserr << "LadrunoQuad::commitState() - failed in base class\n";
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->commitState();
  return retVal;
}

int LadrunoQuad::revertToLastCommit(void)
{
  int retVal = 0;
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->revertToLastCommit();
  return retVal;
}

int LadrunoQuad::revertToStart(void)
{
  int retVal = 0;
  for (int i = 0; i < 4; i++)
    retVal += theMaterial[i]->revertToStart();
  return retVal;
}

void LadrunoQuad::computeShapeBar(void)
{
  for (int a = 0; a < 4; a++) {
    shpBar[0][a] = 0.0;
    shpBar[1][a] = 0.0;
  }
  double A = 0.0;
  for (int i = 0; i < 4; i++) {
    double dv = this->shapeFunction(pts[i][0], pts[i][1]) * wts[i];
    A += dv;
    for (int a = 0; a < 4; a++) {
      shpBar[0][a] += shp[0][a] * dv;
      shpBar[1][a] += shp[1][a] * dv;
    }
  }
  for (int a = 0; a < 4; a++) {
    shpBar[0][a] /= A;
    shpBar[1][a] /= A;
  }
}

void LadrunoQuad::formB(Matrix &B)
{
  B.Zero();
  bool bbar = (formulation == Formulation::BBAR);
  for (int a = 0; a < 4; a++) {
    double bx = shp[0][a];
    double by = shp[1][a];
    if (bbar) {
      // mean-dilatation B-bar: distribute the (mean - local) dilatation
      // equally over the two in-plane normal strains (2D factor 1/2).
      double cx = 0.5 * (shpBar[0][a] - bx);
      double cy = 0.5 * (shpBar[1][a] - by);
      B(0, 2 * a)     = bx + cx;  B(0, 2 * a + 1) = cy;
      B(1, 2 * a)     = cx;       B(1, 2 * a + 1) = by + cy;
    } else {
      B(0, 2 * a)     = bx;       B(0, 2 * a + 1) = 0.0;
      B(1, 2 * a)     = 0.0;      B(1, 2 * a + 1) = by;
    }
    B(2, 2 * a)     = by;
    B(2, 2 * a + 1) = bx;
  }
}

int LadrunoQuad::update(void)
{
  static Vector u(8);
  for (int a = 0; a < 4; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u(2 * a)     = d(0);
    u(2 * a + 1) = d(1);
  }

  static Vector eps(3);

  if (formulation == Formulation::SSP) {
    // single material evaluation at the centroid: strain = Mmem * u
    eps.addMatrixVector(0.0, Mmem, u, 1.0);
    return theMaterial[0]->setTrialStrain(eps);
  }

  if (formulation == Formulation::BBAR)
    this->computeShapeBar();

  static Matrix B(3, 8);
  int ret = 0;
  for (int i = 0; i < 4; i++) {
    this->shapeFunction(pts[i][0], pts[i][1]);
    this->formB(B);
    eps.addMatrixVector(0.0, B, u, 1.0);
    ret += theMaterial[i]->setTrialStrain(eps);
  }
  return ret;
}

const Matrix &LadrunoQuad::getTangentStiff(void)
{
  K.Zero();

  if (formulation == Formulation::SSP) {
    // K = s*Kstab + 4*J0*t * Mmem^T C Mmem  (s = Tier-A damage scale)
    const Matrix &C = theMaterial[0]->getTangent();
    K.addMatrix(0.0, Kstab, this->damageScale());
    K.addMatrixTripleProduct(1.0, Mmem, C, 4.0 * J0 * thickness);
    return K;
  }

  if (formulation == Formulation::BBAR)
    this->computeShapeBar();

  static Matrix B(3, 8);
  for (int i = 0; i < 4; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getTangent();
    this->formB(B);
    K.addMatrixTripleProduct(1.0, B, D, dvol);
  }
  return K;
}

const Matrix &LadrunoQuad::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  K.Zero();

  if (formulation == Formulation::SSP) {
    const Matrix &C = theMaterial[0]->getInitialTangent();
    K.addMatrix(0.0, Kstab, 1.0);   // initial state: undamaged
    K.addMatrixTripleProduct(1.0, Mmem, C, 4.0 * J0 * thickness);
    Ki = new Matrix(K);
    return *Ki;
  }

  if (formulation == Formulation::BBAR)
    this->computeShapeBar();

  static Matrix B(3, 8);
  for (int i = 0; i < 4; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getInitialTangent();
    this->formB(B);
    K.addMatrixTripleProduct(1.0, B, D, dvol);
  }
  Ki = new Matrix(K);
  return K;
}

const Matrix &LadrunoQuad::getMass(void)
{
  K.Zero();

  if (formulation == Formulation::SSP) {
    double density = (rho == 0.0) ? theMaterial[0]->getRho() : rho;
    if (density == 0.0)
      return K;
    const double xi[4]  = {-1.0,  1.0, 1.0, -1.0};
    const double eta[4] = {-1.0, -1.0, 1.0,  1.0};
    for (int i = 0; i < 4; i++) {
      double m = density * thickness * (J0 + J1 * xi[i] + J2 * eta[i]);
      K(2 * i,     2 * i)     += m;
      K(2 * i + 1, 2 * i + 1) += m;
    }
    return K;
  }

  static double rhoi[4];
  double sum = 0.0;
  for (int i = 0; i < 4; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }
  if (sum == 0.0)
    return K;

  for (int i = 0; i < 4; i++) {
    double rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);
    rhodvol *= (rhoi[i] * thickness * wts[i]);
    for (int a = 0, ia = 0; a < 4; a++, ia += 2) {
      double Nrho = shp[2][a] * rhodvol;
      K(ia, ia)         += Nrho;
      K(ia + 1, ia + 1) += Nrho;
    }
  }
  return K;
}

void LadrunoQuad::zeroLoad(void)
{
  Q.Zero();
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
}

int LadrunoQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    return 0;
  }
  opserr << "LadrunoQuad::addLoad - load type unknown for ele " << this->getTag() << "\n";
  return -1;
}

int LadrunoQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  static double rhoi[4];
  double sum = 0.0;
  for (int i = 0; i < 4; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }
  if (sum == 0.0)
    return 0;

  static double ra[8];
  for (int a = 0; a < 4; a++) {
    const Vector &Raccel = theNodes[a]->getRV(accel);
    if (Raccel.Size() != 2) {
      opserr << "LadrunoQuad::addInertiaLoadToUnbalance - matrix/vector sizes incompatible\n";
      return -1;
    }
    ra[2 * a]     = Raccel(0);
    ra[2 * a + 1] = Raccel(1);
  }

  this->getMass();
  for (int i = 0; i < 8; i++)
    Q(i) += -K(i, i) * ra[i];
  return 0;
}

const Vector &LadrunoQuad::getResistingForce(void)
{
  P.Zero();
  static Vector sigma(3);

  if (formulation == Formulation::SSP) {
    // fint = s*Kstab*d + 4*t*J0*Mmem^T*sigma - body - Q
    static Vector d(8);
    for (int a = 0; a < 4; a++) {
      const Vector &di = theNodes[a]->getTrialDisp();
      d(2 * a)     = di(0);
      d(2 * a + 1) = di(1);
    }
    P.addMatrixVector(0.0, Kstab, d, this->damageScale());
    sigma = theMaterial[0]->getStress();
    P.addMatrixTransposeVector(1.0, Mmem, sigma, 4.0 * thickness * J0);

    const double xi[4]  = {-1.0,  1.0, 1.0, -1.0};
    const double eta[4] = {-1.0, -1.0, 1.0,  1.0};
    double bx = (applyLoad == 0) ? b[0] : appliedB[0];
    double by = (applyLoad == 0) ? b[1] : appliedB[1];
    for (int i = 0; i < 4; i++) {
      double w = thickness * (J0 + J1 * xi[i] + J2 * eta[i]);
      P(2 * i)     -= bx * w;
      P(2 * i + 1) -= by * w;
    }
    if (pressure != 0.0)
      P.addVector(1.0, pressureLoad, -1.0);
    P.addVector(1.0, Q, -1.0);
    return P;
  }

  if (formulation == Formulation::BBAR)
    this->computeShapeBar();

  static Matrix B(3, 8);
  for (int i = 0; i < 4; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    this->formB(B);
    sigma = theMaterial[i]->getStress();
    P.addMatrixTransposeVector(1.0, B, sigma, dvol);

    for (int a = 0, ia = 0; a < 4; a++, ia += 2) {
      if (applyLoad == 0) {
        P(ia)     -= dvol * shp[2][a] * b[0];
        P(ia + 1) -= dvol * shp[2][a] * b[1];
      } else {
        P(ia)     -= dvol * shp[2][a] * appliedB[0];
        P(ia + 1) -= dvol * shp[2][a] * appliedB[1];
      }
    }
  }

  if (pressure != 0.0)
    P.addVector(1.0, pressureLoad, -1.0);

  P.addVector(1.0, Q, -1.0);
  return P;
}

const Vector &LadrunoQuad::getResistingForceIncInertia(void)
{
  static double rhoi[4];
  double sum = 0.0;
  for (int i = 0; i < 4; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }

  if (sum == 0.0) {
    this->getResistingForce();
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
    return P;
  }

  static double a[8];
  for (int n = 0; n < 4; n++) {
    const Vector &accel = theNodes[n]->getTrialAccel();
    a[2 * n]     = accel(0);
    a[2 * n + 1] = accel(1);
  }

  this->getResistingForce();
  this->getMass();
  for (int i = 0; i < 8; i++)
    P(i) += K(i, i) * a[i];

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();
  return P;
}

int LadrunoQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(12);
  data(0)  = this->getTag();
  data(1)  = thickness;
  data(2)  = b[0];
  data(3)  = b[1];
  data(4)  = pressure;
  data(5)  = rho;
  data(6)  = static_cast<int>(formulation);
  data(7)  = planeType;
  data(8)  = alphaM;
  data(9)  = betaK;
  data(10) = betaK0;
  data(11) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoQuad::sendSelf() - failed to send Vector\n";
    return res;
  }

  static ID idData(12);
  for (int i = 0; i < 4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    int matDbTag = theMaterial[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i + 4) = matDbTag;
  }
  idData(8)  = connectedExternalNodes(0);
  idData(9)  = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoQuad::sendSelf() - failed to send ID\n";
    return res;
  }

  for (int i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LadrunoQuad::sendSelf() - failed to send material\n";
      return res;
    }
  }
  return res;
}

int LadrunoQuad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(12);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  thickness   = data(1);
  b[0]        = data(2);
  b[1]        = data(3);
  pressure    = data(4);
  rho         = data(5);
  formulation = static_cast<Formulation>((int)data(6));
  planeType   = (int)data(7);
  alphaM      = data(8);
  betaK       = data(9);
  betaK0      = data(10);
  betaKc      = data(11);

  static ID idData(12);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoQuad::recvSelf() - failed to receive ID\n";
    return res;
  }
  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);

  // the cached damage probe points into theMaterial[0], which is about to be
  // re-brokered; drop it so it can't dangle (setDomain rebuilds it).
  if (damageResponse) { delete damageResponse; damageResponse = 0; }

  if (theMaterial == 0) {
    theMaterial = new NDMaterial *[4];
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
        opserr << "LadrunoQuad::recvSelf() - broker could not create NDMaterial " << matClassTag << "\n";
        return -1;
      }
      theMaterial[i]->setDbTag(idData(i + 4));
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) return res;
    }
  } else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      if (theMaterial[i]->getClassTag() != matClassTag) {
        delete theMaterial[i];
        theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial[i] == 0) return -1;
      }
      theMaterial[i]->setDbTag(idData(i + 4));
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) return res;
    }
  }
  return res;
}

void LadrunoQuad::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoQuad, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tformulation:  " << static_cast<int>(formulation) << "  type: " << typeString() << "\n";
    s << "\tthickness:  " << thickness << "  rho: " << rho << "\n";
    theMaterial[0]->Print(s, flag);
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoQuad\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << ", " << connectedExternalNodes(2) << ", " << connectedExternalNodes(3) << "], ";
    s << "\"thickness\": " << thickness << ", ";
    s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int LadrunoQuad::displaySelf(Renderer &theViewer, int displayMode, float fact,
                             const char **modes, int numMode)
{
  static Vector v1(3), v2(3), v3(3), v4(3);
  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  theNodes[2]->getDisplayCrds(v3, fact, displayMode);
  theNodes[3]->getDisplayCrds(v4, fact, displayMode);

  static Matrix coords(4, 3);
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v2(i);
    coords(2, i) = v3(i);
    coords(3, i) = v4(i);
  }
  static Vector values(4);
  const bool sp = this->isSinglePoint();
  if (displayMode < 4 && displayMode > 0)
    for (int i = 0; i < 4; i++)
      values(i) = theMaterial[sp ? 0 : i]->getStress()(displayMode - 1);
  else
    for (int i = 0; i < 4; i++) values(i) = 0.0;

  return theViewer.drawPolygon(coords, values, this->getTag());
}

Response *LadrunoQuad::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoQuad");
  output.attr("eleTag", this->getTag());

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
    theResponse = new ElementResponse(this, 1, P);
  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    // SSP evaluates the material only at the centroid (slot 0).
    int idx = this->isSinglePoint() ? 0 : (pointNum - 1);
    if (pointNum > 0 && pointNum <= 4)
      theResponse = theMaterial[idx]->setResponse(&argv[2], argc - 2, output);
  } else if (strcmp(argv[0], "stresses") == 0 || strcmp(argv[0], "stress") == 0) {
    theResponse = new ElementResponse(this, 3, Vector(12));
  } else if (strcmp(argv[0], "strains") == 0 || strcmp(argv[0], "strain") == 0) {
    theResponse = new ElementResponse(this, 4, Vector(12));
  } else if (strcmp(argv[0], "charLength") == 0 || strcmp(argv[0], "characteristicLength") == 0) {
    theResponse = new ElementResponse(this, 5, 0.0);
  }

  output.endTag();
  return theResponse;
}

int LadrunoQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  if (responseID == 3 || responseID == 4) {
    static Vector v(12);
    const bool sp = this->isSinglePoint();
    int cnt = 0;
    for (int i = 0; i < 4; i++) {
      int idx = sp ? 0 : i;   // SSP: mirror the single centroid material
      const Vector &s = (responseID == 3) ? theMaterial[idx]->getStress()
                                          : theMaterial[idx]->getStrain();
      v(cnt) = s(0); v(cnt + 1) = s(1); v(cnt + 2) = s(2);
      cnt += 3;
    }
    return eleInfo.setVector(v);
  }

  if (responseID == 5)
    return eleInfo.setDouble(this->getCharacteristicLength());

  return -1;
}

int LadrunoQuad::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "pressure") == 0)
    return param.addObject(2, this);

  if (strstr(argv[0], "material") != 0) {
    if (argc < 3) return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4)
      return theMaterial[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
    return -1;
  }

  for (int i = 0; i < 4; i++) {
    int matRes = theMaterial[i]->setParameter(argv, argc, param);
    if (matRes != -1) res = matRes;
  }
  return res;
}

int LadrunoQuad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
    case 2:
      pressure = info.theDouble;
      this->setPressureLoadAtNodes();
      return 0;
    default:
      return -1;
  }
}

double LadrunoQuad::shapeFunction(double xi, double eta)
{
  const Vector &nd1 = theNodes[0]->getCrds();
  const Vector &nd2 = theNodes[1]->getCrds();
  const Vector &nd3 = theNodes[2]->getCrds();
  const Vector &nd4 = theNodes[3]->getCrds();

  double oneMinuseta = 1.0 - eta;
  double onePluseta  = 1.0 + eta;
  double oneMinusxi  = 1.0 - xi;
  double onePlusxi   = 1.0 + xi;

  shp[2][0] = 0.25 * oneMinusxi * oneMinuseta;
  shp[2][1] = 0.25 * onePlusxi  * oneMinuseta;
  shp[2][2] = 0.25 * onePlusxi  * onePluseta;
  shp[2][3] = 0.25 * oneMinusxi * onePluseta;

  double J[2][2];
  J[0][0] = 0.25 * (-nd1(0) * oneMinuseta + nd2(0) * oneMinuseta + nd3(0) * onePluseta - nd4(0) * onePluseta);
  J[0][1] = 0.25 * (-nd1(0) * oneMinusxi  - nd2(0) * onePlusxi   + nd3(0) * onePlusxi  + nd4(0) * oneMinusxi);
  J[1][0] = 0.25 * (-nd1(1) * oneMinuseta + nd2(1) * oneMinuseta + nd3(1) * onePluseta - nd4(1) * onePluseta);
  J[1][1] = 0.25 * (-nd1(1) * oneMinusxi  - nd2(1) * onePlusxi   + nd3(1) * onePlusxi  + nd4(1) * oneMinusxi);

  double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  double oneOverdetJ = 1.0 / detJ;
  double L[2][2];
  L[0][0] =  J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] =  J[0][0] * oneOverdetJ;

  double L00 = 0.25 * L[0][0];
  double L10 = 0.25 * L[1][0];
  double L01 = 0.25 * L[0][1];
  double L11 = 0.25 * L[1][1];

  double L00oneMinuseta = L00 * oneMinuseta;
  double L00onePluseta  = L00 * onePluseta;
  double L01oneMinusxi  = L01 * oneMinusxi;
  double L01onePlusxi   = L01 * onePlusxi;
  double L10oneMinuseta = L10 * oneMinuseta;
  double L10onePluseta  = L10 * onePluseta;
  double L11oneMinusxi  = L11 * oneMinusxi;
  double L11onePlusxi   = L11 * onePlusxi;

  shp[0][0] = -L00oneMinuseta - L01oneMinusxi;
  shp[0][1] =  L00oneMinuseta - L01onePlusxi;
  shp[0][2] =  L00onePluseta  + L01onePlusxi;
  shp[0][3] = -L00onePluseta  + L01oneMinusxi;

  shp[1][0] = -L10oneMinuseta - L11oneMinusxi;
  shp[1][1] =  L10oneMinuseta - L11onePlusxi;
  shp[1][2] =  L10onePluseta  + L11onePlusxi;
  shp[1][3] = -L10onePluseta  + L11oneMinusxi;

  return detJ;
}

void LadrunoQuad::setPressureLoadAtNodes(void)
{
  pressureLoad.Zero();
  if (pressure == 0.0)
    return;

  const Vector &n1 = theNodes[0]->getCrds();
  const Vector &n2 = theNodes[1]->getCrds();
  const Vector &n3 = theNodes[2]->getCrds();
  const Vector &n4 = theNodes[3]->getCrds();

  double x1 = n1(0), y1 = n1(1);
  double x2 = n2(0), y2 = n2(1);
  double x3 = n3(0), y3 = n3(1);
  double x4 = n4(0), y4 = n4(1);

  double dx12 = x2 - x1, dy12 = y2 - y1;
  double dx23 = x3 - x2, dy23 = y3 - y2;
  double dx34 = x4 - x3, dy34 = y4 - y3;
  double dx41 = x1 - x4, dy41 = y1 - y4;

  double po2 = pressure / 2.0;
  pressureLoad(0) += po2 * dy12;  pressureLoad(2) += po2 * dy12;
  pressureLoad(1) += po2 * -dx12; pressureLoad(3) += po2 * -dx12;
  pressureLoad(2) += po2 * dy23;  pressureLoad(4) += po2 * dy23;
  pressureLoad(3) += po2 * -dx23; pressureLoad(5) += po2 * -dx23;
  pressureLoad(4) += po2 * dy34;  pressureLoad(6) += po2 * dy34;
  pressureLoad(5) += po2 * -dx34; pressureLoad(7) += po2 * -dx34;
  pressureLoad(6) += po2 * dy41;  pressureLoad(0) += po2 * dy41;
  pressureLoad(7) += po2 * -dx41; pressureLoad(1) += po2 * -dx41;
}
