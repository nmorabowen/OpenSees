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

// LadrunoCST — 3-node constant-strain triangle (Ladruno fork). std mirrors
// upstream Tri31 (single centroid Gauss point, linear shape functions) so it
// reduces to it to ~1e-9. See LadrunoCST.h and ADR 25.

#include <LadrunoCST.h>
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
#include <elementAPI.h>

double LadrunoCST::matrixData[36];
Matrix LadrunoCST::K(matrixData, 6, 6);
Vector LadrunoCST::P(6);
double LadrunoCST::shp[3][3];
double LadrunoCST::pts[1][2];
double LadrunoCST::wts[1];

LadrunoCST::LadrunoCST(int tag, int nd1, int nd2, int nd3,
                       NDMaterial &m, const char *type, double t,
                       double r, double b1, double b2, double p)
 :Element(tag, ELE_TAG_LadrunoCST),
  theMaterial(0), connectedExternalNodes(3),
  Q(6), pressureLoad(6), thickness(t), pressure(p), rho(r),
  planeType(1), Ki(0)
{
  pts[0][0] = 0.333333333333333;
  pts[0][1] = 0.333333333333333;
  wts[0] = 0.5;

  if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)
    planeType = 1;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0)
    planeType = 2;
  else {
    opserr << "LadrunoCST::LadrunoCST -- improper material type: " << type << "\n";
    exit(-1);
  }

  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b1;
  b[1] = b2;

  theMaterial = new NDMaterial *[numgp];
  for (int i = 0; i < numgp; i++) {
    theMaterial[i] = m.getCopy(type);
    if (theMaterial[i] == 0) {
      opserr << "LadrunoCST::LadrunoCST -- failed to get a copy of material model\n";
      exit(-1);
    }
  }

  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;

  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoCST::LadrunoCST()
 :Element(0, ELE_TAG_LadrunoCST),
  theMaterial(0), connectedExternalNodes(3),
  Q(6), pressureLoad(6), thickness(0.0), pressure(0.0), rho(0.0),
  planeType(1), Ki(0)
{
  pts[0][0] = 0.333333333333333;
  pts[0][1] = 0.333333333333333;
  wts[0] = 0.5;
  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b[1] = 0.0;
  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoCST::~LadrunoCST()
{
  if (theMaterial) {
    for (int i = 0; i < numgp; i++)
      if (theMaterial[i]) delete theMaterial[i];
    delete[] theMaterial;
  }
  if (Ki) delete Ki;
}

const char *LadrunoCST::typeString(void) const
{
  return (planeType == 2) ? "PlaneStress" : "PlaneStrain";
}

int LadrunoCST::getNumExternalNodes(void) const { return numnodes; }
const ID &LadrunoCST::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoCST::getNodePtrs(void) { return theNodes; }
int LadrunoCST::getNumDOF(void) { return 2 * numnodes; }

void LadrunoCST::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    for (int i = 0; i < numnodes; i++) theNodes[i] = 0;
    return;
  }
  for (int i = 0; i < numnodes; i++)
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
  for (int i = 0; i < numnodes; i++)
    if (theNodes[i] == 0) return;
  for (int i = 0; i < numnodes; i++)
    if (theNodes[i]->getNumberDOF() != 2) {
      opserr << "LadrunoCST::setDomain -- node " << connectedExternalNodes(i)
             << " must have 2 DOFs\n";
      return;
    }
  this->DomainComponent::setDomain(theDomain);
  this->setPressureLoadAtNodes();
}

int LadrunoCST::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0)
    opserr << "LadrunoCST::commitState() - failed in base class\n";
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->commitState();
  return retVal;
}

int LadrunoCST::revertToLastCommit(void)
{
  int retVal = 0;
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->revertToLastCommit();
  return retVal;
}

int LadrunoCST::revertToStart(void)
{
  int retVal = 0;
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->revertToStart();
  return retVal;
}

int LadrunoCST::update(void)
{
  static double u[2][3];
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u[0][a] = d(0);
    u[1][a] = d(1);
  }

  static Vector eps(3);
  int ret = 0;
  for (int i = 0; i < numgp; i++) {
    this->shapeFunction(pts[i][0], pts[i][1]);
    eps.Zero();
    for (int beta = 0; beta < numnodes; beta++) {
      eps(0) += shp[0][beta] * u[0][beta];
      eps(1) += shp[1][beta] * u[1][beta];
      eps(2) += shp[0][beta] * u[1][beta] + shp[1][beta] * u[0][beta];
    }
    ret += theMaterial[i]->setTrialStrain(eps);
  }
  return ret;
}

const Matrix &LadrunoCST::getTangentStiff(void)
{
  K.Zero();
  double DB[3][2];
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getTangent();
    double D00 = D(0,0), D01 = D(0,1), D02 = D(0,2);
    double D10 = D(1,0), D11 = D(1,1), D12 = D(1,2);
    double D20 = D(2,0), D21 = D(2,1), D22 = D(2,2);
    for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {
      for (int beta = 0, ib = 0; beta < numnodes; beta++, ib += 2) {
        DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
        DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
        DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
        DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
        DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
        DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
        K(ia,   ib)     += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
        K(ia,   ib + 1) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
        K(ia+1, ib)     += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
        K(ia+1, ib + 1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }
  return K;
}

const Matrix &LadrunoCST::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;
  K.Zero();
  double DB[3][2];
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getInitialTangent();
    double D00 = D(0,0), D01 = D(0,1), D02 = D(0,2);
    double D10 = D(1,0), D11 = D(1,1), D12 = D(1,2);
    double D20 = D(2,0), D21 = D(2,1), D22 = D(2,2);
    for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {
      for (int beta = 0, ib = 0; beta < numnodes; beta++, ib += 2) {
        DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
        DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
        DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
        DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
        DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
        DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
        K(ia,   ib)     += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
        K(ia,   ib + 1) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
        K(ia+1, ib)     += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
        K(ia+1, ib + 1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }
  Ki = new Matrix(K);
  return K;
}

const Matrix &LadrunoCST::getMass(void)
{
  K.Zero();
  double r = (rho == 0.0) ? theMaterial[0]->getRho() : rho;
  if (r == 0.0)
    return K;

  // lumped mass: total mass / 3 to each node
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * r * thickness * wts[i];
    for (int a = 0, ia = 0; a < numnodes; a++, ia += 2) {
      double Nrho = shp[2][a] * dvol;
      K(ia, ia)         += Nrho;
      K(ia + 1, ia + 1) += Nrho;
    }
  }
  return K;
}

void LadrunoCST::zeroLoad(void)
{
  Q.Zero();
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
}

int LadrunoCST::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    return 0;
  }
  opserr << "LadrunoCST::addLoad - load type unknown for ele " << this->getTag() << "\n";
  return -1;
}

int LadrunoCST::addInertiaLoadToUnbalance(const Vector &accel)
{
  double r = (rho == 0.0) ? theMaterial[0]->getRho() : rho;
  if (r == 0.0)
    return 0;

  static double ra[6];
  for (int a = 0; a < numnodes; a++) {
    const Vector &Raccel = theNodes[a]->getRV(accel);
    if (Raccel.Size() != 2) {
      opserr << "LadrunoCST::addInertiaLoadToUnbalance - incompatible sizes\n";
      return -1;
    }
    ra[2 * a]     = Raccel(0);
    ra[2 * a + 1] = Raccel(1);
  }
  this->getMass();
  for (int i = 0; i < 6; i++)
    Q(i) += -K(i, i) * ra[i];
  return 0;
}

const Vector &LadrunoCST::getResistingForce(void)
{
  P.Zero();
  static Vector sigma(3);
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    sigma = theMaterial[i]->getStress();
    for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {
      P(ia)     += dvol * (shp[0][alpha] * sigma(0) + shp[1][alpha] * sigma(2));
      P(ia + 1) += dvol * (shp[1][alpha] * sigma(1) + shp[0][alpha] * sigma(2));
      if (applyLoad == 0) {
        P(ia)     -= dvol * shp[2][alpha] * b[0];
        P(ia + 1) -= dvol * shp[2][alpha] * b[1];
      } else {
        P(ia)     -= dvol * shp[2][alpha] * appliedB[0];
        P(ia + 1) -= dvol * shp[2][alpha] * appliedB[1];
      }
    }
  }
  if (pressure != 0.0)
    P.addVector(1.0, pressureLoad, -1.0);
  P.addVector(1.0, Q, -1.0);
  return P;
}

const Vector &LadrunoCST::getResistingForceIncInertia(void)
{
  double r = (rho == 0.0) ? theMaterial[0]->getRho() : rho;

  if (r == 0.0) {
    this->getResistingForce();
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
    return P;
  }

  static double a[6];
  for (int n = 0; n < numnodes; n++) {
    const Vector &accel = theNodes[n]->getTrialAccel();
    a[2 * n]     = accel(0);
    a[2 * n + 1] = accel(1);
  }
  this->getResistingForce();
  this->getMass();
  for (int i = 0; i < 6; i++)
    P(i) += K(i, i) * a[i];
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();
  return P;
}

int LadrunoCST::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(11);
  data(0)  = this->getTag();
  data(1)  = thickness;
  data(2)  = b[0];
  data(3)  = b[1];
  data(4)  = pressure;
  data(5)  = rho;
  data(6)  = planeType;
  data(7)  = alphaM;
  data(8)  = betaK;
  data(9)  = betaK0;
  data(10) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoCST::sendSelf - send Vector failed\n"; return res; }

  static ID idData(5);
  idData(0) = theMaterial[0]->getClassTag();
  int matDbTag = theMaterial[0]->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0) theMaterial[0]->setDbTag(matDbTag);
  }
  idData(1) = matDbTag;
  idData(2) = connectedExternalNodes(0);
  idData(3) = connectedExternalNodes(1);
  idData(4) = connectedExternalNodes(2);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoCST::sendSelf - send ID failed\n"; return res; }

  res += theMaterial[0]->sendSelf(commitTag, theChannel);
  if (res < 0) { opserr << "WARNING LadrunoCST::sendSelf - send material failed\n"; return res; }
  return res;
}

int LadrunoCST::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(11);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoCST::recvSelf - recv Vector failed\n"; return res; }

  this->setTag((int)data(0));
  thickness = data(1);
  b[0]      = data(2);
  b[1]      = data(3);
  pressure  = data(4);
  rho       = data(5);
  planeType = (int)data(6);
  alphaM    = data(7);
  betaK     = data(8);
  betaK0    = data(9);
  betaKc    = data(10);

  static ID idData(5);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoCST::recvSelf - recv ID failed\n"; return res; }
  connectedExternalNodes(0) = idData(2);
  connectedExternalNodes(1) = idData(3);
  connectedExternalNodes(2) = idData(4);

  int matClassTag = idData(0);
  if (theMaterial == 0) {
    theMaterial = new NDMaterial *[numgp];
    theMaterial[0] = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial[0] == 0) {
      opserr << "LadrunoCST::recvSelf - broker could not create NDMaterial " << matClassTag << "\n";
      return -1;
    }
  } else if (theMaterial[0]->getClassTag() != matClassTag) {
    delete theMaterial[0];
    theMaterial[0] = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial[0] == 0) return -1;
  }
  theMaterial[0]->setDbTag(idData(1));
  res += theMaterial[0]->recvSelf(commitTag, theChannel, theBroker);
  return res;
}

void LadrunoCST::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoCST, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\ttype: " << typeString() << "  thickness:  " << thickness << "  rho: " << rho << "\n";
    theMaterial[0]->Print(s, flag);
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoCST\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << ", " << connectedExternalNodes(2) << "], ";
    s << "\"thickness\": " << thickness << ", ";
    s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int LadrunoCST::displaySelf(Renderer &theViewer, int displayMode, float fact,
                            const char **modes, int numMode)
{
  static Vector v1(3), v2(3), v3(3);
  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);
  theNodes[2]->getDisplayCrds(v3, fact, displayMode);

  static Matrix coords(3, 3);
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v2(i);
    coords(2, i) = v3(i);
  }
  static Vector values(3);
  if (displayMode < 4 && displayMode > 0)
    for (int i = 0; i < numnodes; i++)
      values(i) = theMaterial[0]->getStress()(displayMode - 1);
  else
    for (int i = 0; i < numnodes; i++) values(i) = 0.0;

  return theViewer.drawPolygon(coords, values, this->getTag());
}

Response *LadrunoCST::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoCST");
  output.attr("eleTag", this->getTag());

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
    theResponse = new ElementResponse(this, 1, P);
  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum == 1)
      theResponse = theMaterial[0]->setResponse(&argv[2], argc - 2, output);
  } else if (strcmp(argv[0], "stresses") == 0 || strcmp(argv[0], "stress") == 0) {
    theResponse = new ElementResponse(this, 3, Vector(3));
  } else if (strcmp(argv[0], "strains") == 0 || strcmp(argv[0], "strain") == 0) {
    theResponse = new ElementResponse(this, 4, Vector(3));
  }

  output.endTag();
  return theResponse;
}

int LadrunoCST::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  if (responseID == 3)
    return eleInfo.setVector(theMaterial[0]->getStress());
  if (responseID == 4)
    return eleInfo.setVector(theMaterial[0]->getStrain());
  return -1;
}

int LadrunoCST::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  if (argc < 1)
    return -1;
  if (strcmp(argv[0], "pressure") == 0)
    return param.addObject(2, this);
  if (strstr(argv[0], "material") != 0) {
    if (argc < 3) return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum == 1)
      return theMaterial[0]->setParameter(&argv[2], argc - 2, param);
    return -1;
  }
  for (int i = 0; i < numgp; i++) {
    int matRes = theMaterial[i]->setParameter(argv, argc, param);
    if (matRes != -1) res = matRes;
  }
  return res;
}

int LadrunoCST::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1: return -1;
    case 2:
      pressure = info.theDouble;
      this->setPressureLoadAtNodes();
      return 0;
    default: return -1;
  }
}

double LadrunoCST::shapeFunction(double xi, double eta)
{
  const Vector &nd1 = theNodes[0]->getCrds();
  const Vector &nd2 = theNodes[1]->getCrds();
  const Vector &nd3 = theNodes[2]->getCrds();

  shp[2][0] = xi;
  shp[2][1] = eta;
  shp[2][2] = 1.0 - xi - eta;

  double J[2][2];
  J[0][0] = nd1(0) - nd3(0);
  J[0][1] = nd2(0) - nd3(0);
  J[1][0] = nd1(1) - nd3(1);
  J[1][1] = nd2(1) - nd3(1);

  double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  double oneOverdetJ = 1.0 / detJ;
  double L[2][2];
  L[0][0] =  J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] =  J[0][0] * oneOverdetJ;

  shp[0][0] = L[0][0];
  shp[0][1] = L[0][1];
  shp[0][2] = -(L[0][0] + L[0][1]);

  shp[1][0] = L[1][0];
  shp[1][1] = L[1][1];
  shp[1][2] = -(L[1][0] + L[1][1]);

  return detJ;
}

void LadrunoCST::setPressureLoadAtNodes(void)
{
  pressureLoad.Zero();
  if (pressure == 0.0)
    return;

  const Vector &n1 = theNodes[0]->getCrds();
  const Vector &n2 = theNodes[1]->getCrds();
  const Vector &n3 = theNodes[2]->getCrds();

  double x1 = n1(0), y1 = n1(1);
  double x2 = n2(0), y2 = n2(1);
  double x3 = n3(0), y3 = n3(1);

  double dx12 = x2 - x1, dy12 = y2 - y1;
  double dx23 = x3 - x2, dy23 = y3 - y2;
  double dx31 = x1 - x3, dy31 = y1 - y3;

  double po2 = pressure / 2.0;
  // edge 1-2
  pressureLoad(0) += po2 * dy12;  pressureLoad(2) += po2 * dy12;
  pressureLoad(1) += po2 * -dx12; pressureLoad(3) += po2 * -dx12;
  // edge 2-3
  pressureLoad(2) += po2 * dy23;  pressureLoad(4) += po2 * dy23;
  pressureLoad(3) += po2 * -dx23; pressureLoad(5) += po2 * -dx23;
  // edge 3-1
  pressureLoad(4) += po2 * dy31;  pressureLoad(0) += po2 * dy31;
  pressureLoad(5) += po2 * -dx31; pressureLoad(1) += po2 * -dx31;
}
