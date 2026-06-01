/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LadrunoBrick — unified 8-node hexahedral solid element (Ladruno fork).
//
// v1 implements the small-strain std + bbar formulations (uri/eas reserved).
// std/bbar kernels are adapted from Ed Love's Brick / BbarBrick so that
//   -formulation std  reproduces upstream Brick     to ~1e-12
//   -formulation bbar reproduces upstream bbarBrick  to ~1e-12
// which is our correctness anchor.  See Ladruno_implementation/09_ladruno_brick.md.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <LadrunoBrick.h>
#include <SolidTransformation.h>          // Ladruno — geometry-method seam (2/3)
#include <SolidTransformationLinear.h>    // Ladruno — v1 identity method
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//static data
double  LadrunoBrick::xl[3][8];

Matrix  LadrunoBrick::stiff(24, 24);
Vector  LadrunoBrick::resid(24);
Matrix  LadrunoBrick::mass(24, 24);

//quadrature data (full 2x2x2)
const double  LadrunoBrick::root3 = sqrt(3.0);
const double  LadrunoBrick::one_over_root3 = 1.0 / root3;
const double  LadrunoBrick::sg[] = { -one_over_root3, one_over_root3 };
const double  LadrunoBrick::wg[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

//file-scope B workspace (computeB returns a reference to this, as in Brick)
static Matrix B(6, 3);

// The 4 raw hourglass base vectors for the 8-node hex (Flanagan-Belytschko
// 1981): the nodal values of the bilinear/trilinear modes {xy, yz, zx, xyz}
// the constant centroid gradient cannot represent. Node order matches Brick's
// (nodes 0-3 = z- face, 4-7 = z+ face), natural coords:
//   0(-,-,-) 1(+,-,-) 2(+,+,-) 3(-,+,-) 4(-,-,+) 5(+,-,+) 6(+,+,+) 7(-,+,+)
const double LadrunoBrick::hg0[4][8] = {
  { 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0 },  // xy
  { 1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0 },  // yz
  { 1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  1.0, -1.0 },  // zx
  {-1.0,  1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0 }   // xyz
};

const char *
LadrunoBrick::formulationName(Formulation f)
{
  switch (f) {
  case Formulation::STD:  return "std";
  case Formulation::BBAR: return "bbar";
  case Formulation::URI:  return "uri";
  case Formulation::EAS:  return "eas";
  }
  return "std";
}

//null constructor (broker)
LadrunoBrick::LadrunoBrick()
  :Element(0, ELE_TAG_LadrunoBrick),
   connectedExternalNodes(8),
   formulation(Formulation::STD),
   hourglassType(Hourglass::PHYSICAL), hourglassCoeff(0.0),
   applyLoad(0), load(0), Ki(0), massType(0),
   theGeom(new SolidTransformationLinear())   // Ladruno — v1 identity geometry
{
  B.Zero();

  for (int i = 0; i < 8; i++) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
    theDamping[i] = 0;
  }

  b[0] = 0.0; b[1] = 0.0; b[2] = 0.0;
}

//full constructor
LadrunoBrick::LadrunoBrick(int tag,
                           int node1, int node2, int node3, int node4,
                           int node5, int node6, int node7, int node8,
                           NDMaterial &theMaterial,
                           Formulation form,
                           double b1, double b2, double b3,
                           int matype,
                           Hourglass hgType, double hgCoeff,
                           Damping *damping)
  :Element(tag, ELE_TAG_LadrunoBrick),
   connectedExternalNodes(8),
   formulation(form),
   hourglassType(hgType), hourglassCoeff(hgCoeff),
   applyLoad(0), load(0), Ki(0), massType(matype),
   theGeom(new SolidTransformationLinear())   // Ladruno — v1 identity geometry
{
  B.Zero();

  connectedExternalNodes(0) = node1;
  connectedExternalNodes(1) = node2;
  connectedExternalNodes(2) = node3;
  connectedExternalNodes(3) = node4;
  connectedExternalNodes(4) = node5;
  connectedExternalNodes(5) = node6;
  connectedExternalNodes(6) = node7;
  connectedExternalNodes(7) = node8;

  for (int i = 0; i < 8; i++) {
    materialPointers[i] = theMaterial.getCopy("ThreeDimensional");
    if (materialPointers[i] == 0) {
      opserr << "LadrunoBrick::constructor - failed to get a material of type: ThreeDimensional\n";
    }
    nodePointers[i] = 0;
  }

  b[0] = b1; b[1] = b2; b[2] = b3;

  if (damping) {
    for (int i = 0; i < 8; i++) {
      theDamping[i] = (*damping).getCopy();
      if (!theDamping[i])
        opserr << "LadrunoBrick::LadrunoBrick -- failed to get copy of damping\n";
    }
  } else {
    for (int i = 0; i < 8; i++) theDamping[i] = 0;
  }
}

//destructor
LadrunoBrick::~LadrunoBrick()
{
  for (int i = 0; i < 8; i++) {
    if (materialPointers[i]) delete materialPointers[i];
  }

  if (load != 0) delete load;
  if (Ki != 0) delete Ki;

  for (int i = 0; i < 8; i++) {
    if (theDamping[i]) {
      delete theDamping[i];
      theDamping[i] = 0;
    }
  }

  if (theGeom) delete theGeom;   // Ladruno — geometry-method layer
}

//set domain
void  LadrunoBrick::setDomain(Domain *theDomain)
{
  for (int i = 0; i < 8; i++)
    nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

  for (int i = 0; i < 8; i++) {
    if (theDamping[i] && theDamping[i]->setDomain(theDomain, 6)) {
      opserr << "LadrunoBrick::setDomain -- Error initializing damping\n";
      return;
    }
  }

  this->DomainComponent::setDomain(theDomain);
}

int
LadrunoBrick::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping) {
    for (int i = 0; i < 8; i++) {
      if (theDamping[i]) delete theDamping[i];
      theDamping[i] = (*damping).getCopy();
      if (!theDamping[i]) {
        opserr << "LadrunoBrick::setDamping -- failed to get copy of damping\n";
        return -1;
      }
      if (theDamping[i]->setDomain(theDomain, 6)) {
        opserr << "LadrunoBrick::setDamping -- Error initializing damping\n";
        return -2;
      }
    }
  }
  return 0;
}

int  LadrunoBrick::getNumExternalNodes(void) const { return 8; }
const ID &  LadrunoBrick::getExternalNodes(void) { return connectedExternalNodes; }
Node **  LadrunoBrick::getNodePtrs(void) { return nodePointers; }
int  LadrunoBrick::getNumDOF(void) { return 24; }

//commit state
int  LadrunoBrick::commitState(void)
{
  int success = 0;

  if ((success = this->Element::commitState()) != 0)
    opserr << "LadrunoBrick::commitState () - failed in base class";

  for (int i = 0; i < 8; i++)
    success += materialPointers[i]->commitState();

  for (int i = 0; i < 8; i++)
    if (theDamping[i]) success += theDamping[i]->commitState();

  return success;
}

int  LadrunoBrick::revertToLastCommit(void)
{
  int success = 0;
  for (int i = 0; i < 8; i++)
    success += materialPointers[i]->revertToLastCommit();
  for (int i = 0; i < 8; i++)
    if (theDamping[i]) success += theDamping[i]->revertToLastCommit();
  return success;
}

int  LadrunoBrick::revertToStart(void)
{
  int success = 0;
  for (int i = 0; i < 8; i++)
    success += materialPointers[i]->revertToStart();
  for (int i = 0; i < 8; i++)
    if (theDamping[i]) success += theDamping[i]->revertToStart();
  return success;
}

//print
void  LadrunoBrick::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "LadrunoBrick (8-node hexahedron), formulation: " << formulationName(formulation) << "\n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Nodes: " << connectedExternalNodes;
    s << "Material Information : \n ";
    materialPointers[0]->Print(s, flag);
    s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << endln;
    s << "Resisting Force (no inertia): " << this->getResistingForce();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoBrick\", ";
    s << "\"formulation\": \"" << formulationName(formulation) << "\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
    for (int i = 1; i < 7; i++)
      s << connectedExternalNodes(i) << ", ";
    s << connectedExternalNodes(7) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
  }
}

//tangent stiffness
const Matrix &  LadrunoBrick::getTangentStiff(void)
{
  formResidAndTangent(1);
  return stiff;
}

//mass
const Matrix &  LadrunoBrick::getMass(void)
{
  formInertiaTerms(1);
  return mass;
}

//----------------------------------------------------------------------
// bbar helper: element-averaged shape-function gradients (mean dilatation).
// shpBar[p][q] = (1/V) * sum_g dvol[g] * Shape[p][q][g]   (matches BbarBrick)
//----------------------------------------------------------------------
void
LadrunoBrick::computeShapeBar(double shpBar[4][8], double Shape[4][8][8],
                              const double dvol[8], double volume)
{
  for (int p = 0; p < 4; p++)
    for (int q = 0; q < 8; q++)
      shpBar[p][q] = 0.0;

  for (int g = 0; g < 8; g++)
    for (int p = 0; p < 4; p++)
      for (int q = 0; q < 8; q++)
        shpBar[p][q] += dvol[g] * Shape[p][q][g];

  for (int p = 0; p < 4; p++)
    for (int q = 0; q < 8; q++)
      shpBar[p][q] /= volume;
}

//initial stiffness
const Matrix &  LadrunoBrick::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  if (formulation == Formulation::URI) {
    formUri(1, true);
    Ki = new Matrix(stiff);
    return *Ki;
  }

  static const int ndf = 3;
  static const int nstress = 6;
  static const int numberNodes = 8;
  static const int numberGauss = 8;
  static const int nShape = 4;

  double xsj;
  double volume = 0.0;
  static double dvol[numberGauss];
  static double gaussPoint[3];
  static double shp[nShape][numberNodes];
  static double Shape[nShape][numberNodes][numberGauss];
  static double shpBar[nShape][numberNodes];
  static Matrix stiffJK(ndf, ndf);
  static Matrix dd(nstress, nstress);
  static Matrix BJ(nstress, ndf);
  static Matrix BJtran(ndf, nstress);
  static Matrix BK(nstress, ndf);
  static Matrix BJtranD(ndf, nstress);

  const bool useBbar = (formulation == Formulation::BBAR);

  stiff.Zero();
  computeBasis();

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        gaussPoint[0] = sg[i]; gaussPoint[1] = sg[j]; gaussPoint[2] = sg[k];
        shp3d(gaussPoint, xsj, shp, xl);
        for (int p = 0; p < nShape; p++)
          for (int q = 0; q < numberNodes; q++)
            Shape[p][q][count] = shp[p][q];
        dvol[count] = wg[count] * xsj;
        volume += dvol[count];
        count++;
      }
    }
  }

  if (useBbar)
    computeShapeBar(shpBar, Shape, dvol, volume);

  for (int i = 0; i < numberGauss; i++) {
    for (int p = 0; p < nShape; p++)
      for (int q = 0; q < numberNodes; q++)
        shp[p][q] = Shape[p][q][i];

    dd = materialPointers[i]->getInitialTangent();
    if (theDamping[i]) dd *= theDamping[i]->getStiffnessMultiplier();
    dd *= dvol[i];

    int jj = 0;
    for (int j = 0; j < numberNodes; j++) {
      BJ = useBbar ? computeBbar(j, shp, shpBar) : computeB(j, shp);

      for (int p = 0; p < ndf; p++)
        for (int q = 0; q < nstress; q++)
          BJtran(p, q) = BJ(q, p);

      BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);

      int kk = 0;
      for (int k = 0; k < numberNodes; k++) {
        BK = useBbar ? computeBbar(k, shp, shpBar) : computeB(k, shp);
        stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);
        for (int p = 0; p < ndf; p++)
          for (int q = 0; q < ndf; q++)
            stiff(jj + p, kk + q) += stiffJK(p, q);
        kk += ndf;
      }
      jj += ndf;
    }
  }

  // seam 3: globalize the initial stiffness (zero stress => no K_geo).
  // identity for -geom linear.  // Ladruno
  static Vector zeroF(24);
  theGeom->globalizeStiff(stiff, zeroF, stiff);

  Ki = new Matrix(stiff);
  return stiff;
}

//----------------------------------------------------------------------
void  LadrunoBrick::zeroLoad(void)
{
  if (load != 0) load->Zero();
  applyLoad = 0;
  appliedB[0] = 0.0; appliedB[1] = 0.0; appliedB[2] = 0.0;
}

int
LadrunoBrick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * b[0];
    appliedB[1] += loadFactor * b[1];
    appliedB[2] += loadFactor * b[2];
    return 0;
  } else if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    appliedB[2] += loadFactor * data(2) * b[2];
    return 0;
  } else {
    opserr << "LadrunoBrick::addLoad() - ele with tag: " << this->getTag()
           << " does not deal with load type: " << type << "\n";
    return -1;
  }
}

int
LadrunoBrick::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberNodes = 8;
  static const int numberGauss = 8;
  static const int ndf = 3;

  int haveRho = 0;
  for (int i = 0; i < numberGauss; i++)
    if (materialPointers[i]->getRho() != 0.0) haveRho = 1;
  if (haveRho == 0) return 0;

  formInertiaTerms(1);

  int count = 0;
  for (int i = 0; i < numberNodes; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j = 0; j < ndf; j++)
      resid(count++) = Raccel(j);
  }

  if (load == 0) load = new Vector(numberNodes * ndf);
  load->addMatrixVector(1.0, mass, resid, -1.0);
  return 0;
}

//residual
const Vector &  LadrunoBrick::getResistingForce(void)
{
  formResidAndTangent(0);
  if (load != 0) resid -= *load;
  return resid;
}

const Vector &  LadrunoBrick::getResistingForceIncInertia(void)
{
  static Vector res(24);

  formResidAndTangent(0);
  formInertiaTerms(0);

  res = resid;

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  if (load != 0) res -= *load;

  return res;
}

//inertia terms (mass is formulation-independent: full 2x2x2 of rho*N^T N)
void   LadrunoBrick::formInertiaTerms(int tangFlag)
{
  static const int ndm = 3;
  static const int ndf = 3;
  static const int numberNodes = 8;
  static const int numberGauss = 8;
  static const int nShape = 4;
  static const int massIndex = nShape - 1;

  double xsj;
  double dvol[numberGauss];
  static double shp[nShape][numberNodes];
  static double Shape[nShape][numberNodes][numberGauss];
  static double gaussPoint[ndm];
  static Vector momentum(ndf);

  double temp, rho, massJK;

  mass.Zero();
  computeBasis();

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        gaussPoint[0] = sg[i]; gaussPoint[1] = sg[j]; gaussPoint[2] = sg[k];
        shp3d(gaussPoint, xsj, shp, xl);
        for (int p = 0; p < nShape; p++)
          for (int q = 0; q < numberNodes; q++)
            Shape[p][q][count] = shp[p][q];
        dvol[count] = wg[count] * xsj;
        count++;
      }
    }
  }

  for (int i = 0; i < numberGauss; i++) {
    for (int p = 0; p < nShape; p++)
      for (int q = 0; q < numberNodes; q++)
        shp[p][q] = Shape[p][q][i];

    momentum.Zero();
    for (int j = 0; j < numberNodes; j++)
      momentum.addVector(1.0, nodePointers[j]->getTrialAccel(), shp[massIndex][j]);

    rho = materialPointers[i]->getRho();
    momentum *= rho;

    int jj = 0;
    for (int j = 0; j < numberNodes; j++) {
      temp = shp[massIndex][j] * dvol[i];

      for (int p = 0; p < ndf; p++)
        resid(jj + p) += (temp * momentum(p));

      if (tangFlag == 1) {
        temp *= rho;
        int kk = 0;
        for (int k = 0; k < numberNodes; k++) {
          if (massType == 0) {
            massJK = temp * shp[massIndex][k];
            for (int p = 0; p < ndf; p++)
              mass(jj + p, kk + p) += massJK;
          } else {
            if (j == k) {
              massJK = temp;
              for (int p = 0; p < ndf; p++)
                mass(jj + p, kk + p) += massJK;
            }
          }
          kk += ndf;
        }
      }
      jj += ndf;
    }
  }
}

//update — seam 1 (kinematics ledger): set the material trial strain
int
LadrunoBrick::update(void)
{
  // uri: single centroid integration point; set its strain on all 8 material
  // pointers so stress/strain output reports the (uniform) centroid value.
  if (formulation == Formulation::URI) {
    computeBasis();
    const Vector &uCore = this->computeLocalDisp();   // seam 0+2 (identity for linear)
    double gp[3] = {0.0, 0.0, 0.0};
    double detJ;
    static double shpC[4][8];
    shp3d(gp, detJ, shpC, xl);

    static Vector strainC(6);
    static Matrix Bc(6, 3);
    static Vector ulj(3);
    strainC.Zero();
    for (int J = 0; J < 8; J++) {
      Bc = computeB(J, shpC);
      ulj(0) = uCore(3 * J); ulj(1) = uCore(3 * J + 1); ulj(2) = uCore(3 * J + 2);
      strainC.addMatrixVector(1.0, Bc, ulj, 1.0);
    }
    for (int i = 0; i < 8; i++)
      materialPointers[i]->setTrialStrain(strainC);
    return 0;
  }

  static const int ndf = 3;
  static const int nstress = 6;
  static const int numberNodes = 8;
  static const int numberGauss = 8;
  static const int nShape = 4;

  double xsj;
  double volume = 0.0;
  static double dvol[numberGauss];
  static double gaussPoint[3];
  static Vector strain(nstress);
  static double shp[nShape][numberNodes];
  static double Shape[nShape][numberNodes][numberGauss];
  static double shpBar[nShape][numberNodes];
  static Matrix BJ(nstress, ndf);

  const bool useBbar = (formulation == Formulation::BBAR);

  computeBasis();

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        gaussPoint[0] = sg[i]; gaussPoint[1] = sg[j]; gaussPoint[2] = sg[k];
        shp3d(gaussPoint, xsj, shp, xl);
        for (int p = 0; p < nShape; p++)
          for (int q = 0; q < numberNodes; q++)
            Shape[p][q][count] = shp[p][q];
        dvol[count] = wg[count] * xsj;
        volume += dvol[count];
        count++;
      }
    }
  }

  if (useBbar)
    computeShapeBar(shpBar, Shape, dvol, volume);

  const Vector &uCore = this->computeLocalDisp();   // seam 0+2 (identity for linear)
  static Vector ulj(3);

  for (int i = 0; i < numberGauss; i++) {
    for (int p = 0; p < nShape; p++)
      for (int q = 0; q < numberNodes; q++)
        shp[p][q] = Shape[p][q][i];

    strain.Zero();

    for (int j = 0; j < numberNodes; j++) {
      BJ = useBbar ? computeBbar(j, shp, shpBar) : computeB(j, shp);
      ulj(0) = uCore(3 * j); ulj(1) = uCore(3 * j + 1); ulj(2) = uCore(3 * j + 2);
      strain.addMatrixVector(1.0, BJ, ulj, 1.0);
    }

    materialPointers[i]->setTrialStrain(strain);
  }

  return 0;
}

//form residual and tangent
void  LadrunoBrick::formResidAndTangent(int tang_flag)
{
  if (formulation == Formulation::URI) {
    formUri(tang_flag, false);
    return;
  }

  static const int ndf = 3;
  static const int nstress = 6;
  static const int numberNodes = 8;
  static const int numberGauss = 8;
  static const int nShape = 4;

  double xsj;
  double volume = 0.0;
  static double dvol[numberGauss];
  static double gaussPoint[3];
  static double shp[nShape][numberNodes];
  static double Shape[nShape][numberNodes][numberGauss];
  static double shpBar[nShape][numberNodes];
  static Vector residJ(ndf);
  static Matrix stiffJK(ndf, ndf);
  static Vector stress(nstress);
  static Vector dampingStress(nstress);
  static Matrix dd(nstress, nstress);
  static Matrix BJ(nstress, ndf);
  static Matrix BJtran(ndf, nstress);
  static Matrix BK(nstress, ndf);
  static Matrix BJtranD(ndf, nstress);

  const bool useBbar = (formulation == Formulation::BBAR);

  stiff.Zero();
  resid.Zero();
  computeBasis();

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        gaussPoint[0] = sg[i]; gaussPoint[1] = sg[j]; gaussPoint[2] = sg[k];
        shp3d(gaussPoint, xsj, shp, xl);
        for (int p = 0; p < nShape; p++)
          for (int q = 0; q < numberNodes; q++)
            Shape[p][q][count] = shp[p][q];
        dvol[count] = wg[count] * xsj;
        volume += dvol[count];
        count++;
      }
    }
  }

  if (useBbar)
    computeShapeBar(shpBar, Shape, dvol, volume);

  for (int i = 0; i < numberGauss; i++) {
    for (int p = 0; p < nShape; p++)
      for (int q = 0; q < numberNodes; q++)
        shp[p][q] = Shape[p][q][i];

    stress = materialPointers[i]->getStress();

    if (theDamping[i]) {
      theDamping[i]->update(stress);
      dampingStress = theDamping[i]->getDampingForce();
      dampingStress *= dvol[i];
    }

    stress *= dvol[i];

    if (tang_flag == 1) {
      dd = materialPointers[i]->getTangent();
      if (theDamping[i]) dd *= theDamping[i]->getStiffnessMultiplier();
      dd *= dvol[i];
    }

    int jj = 0;
    for (int j = 0; j < numberNodes; j++) {
      BJ = useBbar ? computeBbar(j, shp, shpBar) : computeB(j, shp);

      for (int p = 0; p < ndf; p++)
        for (int q = 0; q < nstress; q++)
          BJtran(p, q) = BJ(q, p);

      //residual: B^T * sigma
      residJ.addMatrixVector(0.0, BJtran, stress, 1.0);

      if (theDamping[i]) {
        static Vector dForce(ndf);
        dForce.addMatrixVector(0.0, BJtran, dampingStress, 1.0);
        residJ += dForce;
      }

      for (int p = 0; p < ndf; p++) {
        resid(jj + p) += residJ(p);
        if (applyLoad == 0)
          resid(jj + p) -= dvol[i] * b[p] * shp[3][j];
        else
          resid(jj + p) -= dvol[i] * appliedB[p] * shp[3][j];
      }

      if (tang_flag == 1) {
        BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);
        int kk = 0;
        for (int k = 0; k < numberNodes; k++) {
          BK = useBbar ? computeBbar(k, shp, shpBar) : computeB(k, shp);
          stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);
          for (int p = 0; p < ndf; p++)
            for (int q = 0; q < ndf; q++)
              stiff(jj + p, kk + q) += stiffJK(p, q);
          kk += ndf;
        }
      }
      jj += ndf;
    }
  }

  // seam 3: globalize the core-frame f/K back to global DOFs, adding K_geo.
  // identity for -geom linear (fGlobal=fCore, kGlobal=kCore, K_geo=0).  // Ladruno
  if (tang_flag == 1)
    theGeom->globalizeStiff(stiff, resid, stiff);
  theGeom->globalizeForce(resid, resid);
}

//compute local nodal coordinates
void   LadrunoBrick::computeBasis(void)
{
  for (int i = 0; i < 8; i++) {
    const Vector &coorI = nodePointers[i]->getCrds();
    xl[0][i] = coorI(0);
    xl[1][i] = coorI(1);
    xl[2][i] = coorI(2);
  }
}

//----------------------------------------------------------------------
// Seam 0+2 (geometry method): refresh theGeom from the current geometry and
// return the localized (core-frame) 24-dof trial displacement.
//   * theGeom->update() is the one-shot per-evaluation refresh (linear: no-op;
//     corot/finite: build the frame / deformation gradient from ref+cur coords).
//   * theGeom->localizeDisp() strips the rigid rotation (linear: identity, so
//     uCore == uGlobal and the strain below is bit-for-bit the direct kernel).
//----------------------------------------------------------------------
const Vector &
LadrunoBrick::computeLocalDisp(void)
{
  static Matrix refCrds(8, 3), curCrds(8, 3);
  static Vector uGlobal(24), uCore(24);

  for (int i = 0; i < 8; i++) {
    const Vector &X = nodePointers[i]->getCrds();
    const Vector &u = nodePointers[i]->getTrialDisp();
    for (int d = 0; d < 3; d++) {
      refCrds(i, d) = X(d);
      curCrds(i, d) = X(d) + u(d);
      uGlobal(3 * i + d) = u(d);
    }
  }

  theGeom->update(8, refCrds, curCrds);
  theGeom->localizeDisp(uGlobal, uCore);
  return uCore;
}

//compute the standard small-strain B at a node (std / uri / shear rows of bbar)
//
//               | N,1      0     0    |
//   B       =   |   0     N,2    0    |
//               |   0      0     N,3  |   (6x3)
//               | N,2     N,1     0   |
//               |   0     N,3    N,2  |
//               | N,3      0     N,1  |
const Matrix &
LadrunoBrick::computeB(int node, const double shp[4][8])
{
  B(0, 0) = shp[0][node];
  B(1, 1) = shp[1][node];
  B(2, 2) = shp[2][node];

  B(3, 0) = shp[1][node];
  B(3, 1) = shp[0][node];

  B(4, 1) = shp[2][node];
  B(4, 2) = shp[1][node];

  B(5, 0) = shp[2][node];
  B(5, 2) = shp[0][node];

  return B;
}

//compute the mean-dilatation B-bar at a node (bbar formulation).
//Bbar = Bdev + Bvol, with the volumetric part using element-mean gradients.
const Matrix &
LadrunoBrick::computeBbar(int node, const double shp[4][8], const double shpBar[4][8])
{
  static Matrix Bbar(6, 3);
  static double Bdev[3][3];
  static double BbarVol[3][3];
  static const double one3 = 1.0 / 3.0;

  Bbar.Zero();

  //deviatoric
  Bdev[0][0] = 2.0 * shp[0][node];
  Bdev[0][1] =      -shp[1][node];
  Bdev[0][2] =      -shp[2][node];

  Bdev[1][0] =      -shp[0][node];
  Bdev[1][1] = 2.0 * shp[1][node];
  Bdev[1][2] =      -shp[2][node];

  Bdev[2][0] =      -shp[0][node];
  Bdev[2][1] =      -shp[1][node];
  Bdev[2][2] = 2.0 * shp[2][node];

  //volumetric (mean dilatation)
  BbarVol[0][0] = shpBar[0][node];
  BbarVol[0][1] = shpBar[1][node];
  BbarVol[0][2] = shpBar[2][node];

  BbarVol[1][0] = shpBar[0][node];
  BbarVol[1][1] = shpBar[1][node];
  BbarVol[1][2] = shpBar[2][node];

  BbarVol[2][0] = shpBar[0][node];
  BbarVol[2][1] = shpBar[1][node];
  BbarVol[2][2] = shpBar[2][node];

  //extensional terms
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Bbar(i, j) = one3 * (Bdev[i][j] + BbarVol[i][j]);

  //shear terms
  Bbar(3, 0) = shp[1][node];
  Bbar(3, 1) = shp[0][node];

  Bbar(4, 1) = shp[2][node];
  Bbar(4, 2) = shp[1][node];

  Bbar(5, 0) = shp[2][node];
  Bbar(5, 2) = shp[0][node];

  return Bbar;
}

//----------------------------------------------------------------------
// Flanagan-Belytschko corrected hourglass vectors. gamma is orthogonal to the
// linear displacement field, so hourglass control is consistent (adds exactly
// zero force/stiffness to any rigid-body or constant-strain state -> the patch
// test passes for any kappa). Uses the current nodal coords xl (set by
// computeBasis) and the centroid gradients bC[i][I] = dN_I/dx_i.
//----------------------------------------------------------------------
void
LadrunoBrick::computeGammaHourglass(const double bC[3][8], double gamma[4][8])
{
  for (int a = 0; a < 4; a++) {
    double hx[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++)
      for (int J = 0; J < 8; J++)
        hx[i] += hg0[a][J] * xl[i][J];

    for (int I = 0; I < 8; I++) {
      double corr = 0.0;
      for (int i = 0; i < 3; i++)
        corr += bC[i][I] * hx[i];
      gamma[a][I] = hg0[a][I] - corr;
    }
  }
}

//----------------------------------------------------------------------
// uri: 1-point (centroid) reduced integration + Flanagan-Belytschko stiffness
// hourglass control. Assembles stiff (tang_flag==1) and, when
// !useInitialTangent, resid. The single integration point captures the 6
// constant-strain modes; hourglass control restores stiffness to the 12
// otherwise-zero-energy modes.
//----------------------------------------------------------------------
void
LadrunoBrick::formUri(int tang_flag, bool useInitialTangent)
{
  static const int ndf = 3;
  static const int nstress = 6;
  static const int numberNodes = 8;

  computeBasis();

  // centroid shape functions & gradients (the single integration point)
  double gp[3] = {0.0, 0.0, 0.0};
  double detJ;
  static double shpC[4][8];
  shp3d(gp, detJ, shpC, xl);
  const double vol = 8.0 * detJ;   // 1-point Gauss weight (2^3) x |J| at centroid

  double bC[3][8];                 // centroid gradients dN_I/dx_i
  for (int i = 0; i < 3; i++)
    for (int I = 0; I < numberNodes; I++)
      bC[i][I] = shpC[i][I];

  stiff.Zero();
  if (!useInitialTangent)
    resid.Zero();

  // material tangent at the centroid (drives K and the hourglass modulus)
  static Matrix dd(nstress, nstress);
  dd = useInitialTangent ? materialPointers[0]->getInitialTangent()
                         : materialPointers[0]->getTangent();

  // hourglass stiffness coefficient (FB "stiffness" form):
  //   kappa = scale * G * vol * sum_iI bC_iI^2 ,  G ~ dd(3,3) (shear row)
  double bb = 0.0;
  for (int i = 0; i < 3; i++)
    for (int I = 0; I < numberNodes; I++)
      bb += bC[i][I] * bC[i][I];
  const double Ghg = dd(3, 3);
  const double scale = (hourglassCoeff > 0.0) ? hourglassCoeff : 0.05;
  const double kappa = scale * Ghg * vol * bb;

  double gamma[4][8];
  computeGammaHourglass(bC, gamma);

  // --- 1-point constant-strain contribution ---
  static Matrix ddVol(nstress, nstress);
  if (tang_flag == 1) {
    ddVol = dd;
    ddVol *= vol;
  }

  static Vector stress0(nstress);
  if (!useInitialTangent) {
    stress0 = materialPointers[0]->getStress();
    stress0 *= vol;
  }

  static Matrix BJ(nstress, ndf), BK(nstress, ndf);
  static Matrix BJtran(ndf, nstress), BJtranD(ndf, nstress);
  static Matrix stiffJK(ndf, ndf);
  static Vector residJ(ndf);

  int jj = 0;
  for (int j = 0; j < numberNodes; j++) {
    BJ = computeB(j, shpC);
    for (int p = 0; p < ndf; p++)
      for (int q = 0; q < nstress; q++)
        BJtran(p, q) = BJ(q, p);

    if (!useInitialTangent) {
      residJ.addMatrixVector(0.0, BJtran, stress0, 1.0);
      for (int p = 0; p < ndf; p++) {
        resid(jj + p) += residJ(p);
        if (applyLoad == 0)
          resid(jj + p) -= vol * b[p] * shpC[3][j];      // b = member body force; Nc = shpC[3][j]
        else
          resid(jj + p) -= vol * appliedB[p] * shpC[3][j];
      }
    }

    if (tang_flag == 1) {
      BJtranD.addMatrixProduct(0.0, BJtran, ddVol, 1.0);
      int kk = 0;
      for (int k = 0; k < numberNodes; k++) {
        BK = computeB(k, shpC);
        stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);
        for (int p = 0; p < ndf; p++)
          for (int q = 0; q < ndf; q++)
            stiff(jj + p, kk + q) += stiffJK(p, q);
        kk += ndf;
      }
    }
    jj += ndf;
  }

  // --- Flanagan-Belytschko stiffness hourglass control ---
  if (!useInitialTangent) {
    // generalized hourglass strains q_ai = sum_J gamma_aJ u_iJ
    double q[4][3];
    for (int a = 0; a < 4; a++)
      for (int i = 0; i < 3; i++)
        q[a][i] = 0.0;
    for (int J = 0; J < numberNodes; J++) {
      const Vector &ul = nodePointers[J]->getTrialDisp();
      for (int a = 0; a < 4; a++)
        for (int i = 0; i < 3; i++)
          q[a][i] += gamma[a][J] * ul(i);
    }
    // hourglass nodal forces f_iI = kappa * sum_a gamma_aI q_ai
    for (int I = 0; I < numberNodes; I++)
      for (int i = 0; i < 3; i++) {
        double f = 0.0;
        for (int a = 0; a < 4; a++)
          f += gamma[a][I] * q[a][i];
        resid(3 * I + i) += kappa * f;
      }
  }

  if (tang_flag == 1) {
    // hourglass stiffness K_(iI)(jJ) = delta_ij * kappa * sum_a gamma_aI gamma_aJ
    for (int I = 0; I < numberNodes; I++)
      for (int J = 0; J < numberNodes; J++) {
        double g = 0.0;
        for (int a = 0; a < 4; a++)
          g += gamma[a][I] * gamma[a][J];
        const double kg = kappa * g;
        for (int i = 0; i < 3; i++)
          stiff(3 * I + i, 3 * J + i) += kg;
      }
  }

  // seam 3: globalize the core-frame f/K back to global DOFs, adding K_geo.
  // identity for -geom linear. resid is assembled only when !useInitialTangent;
  // the initial-stiffness call uses zero internal force (zero stress => no K_geo).
  static Vector zeroF(24);                                                  // Ladruno
  const Vector &fCore = useInitialTangent ? zeroF : resid;
  if (tang_flag == 1)
    theGeom->globalizeStiff(stiff, fCore, stiff);
  if (!useInitialTangent)
    theGeom->globalizeForce(resid, resid);
}

//----------------------------------------------------------------------

Matrix  LadrunoBrick::transpose(int dim1, int dim2, const Matrix &M)
{
  Matrix Mtran(dim2, dim1);
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim2; j++)
      Mtran(j, i) = M(i, j);
  return Mtran;
}

//----------------------------------------------------------------------
// sendSelf / recvSelf
//----------------------------------------------------------------------
int  LadrunoBrick::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();
  int matDbTag;

  static ID idData(29);

  idData(24) = this->getTag();
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0)
    idData(25) = 1;
  else
    idData(25) = 0;

  for (int i = 0; i < 8; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i + 8) = matDbTag;
  }

  for (int i = 0; i < 8; i++)
    idData(16 + i) = connectedExternalNodes(i);

  idData(26) = 0;
  idData(27) = 0;
  if (theDamping[0]) {
    idData(26) = theDamping[0]->getClassTag();
    int dbTag = theDamping[0]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        for (int i = 0; i < 8; i++)
          theDamping[i]->setDbTag(dbTag);
    }
    idData(27) = dbTag;
  }

  // formulation + massType + hourglass type + geometry method packed into one
  // int slot. geometry method id is x1000; Linear=0 => byte-identical to a
  // pre-seam stream.  // Ladruno
  idData(28) = static_cast<int>(formulation)
             + 10 * massType
             + 100 * static_cast<int>(hourglassType)
             + 1000 * theGeom->getMethodID();

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoBrick::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector dData(8);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];
  dData(7) = hourglassCoeff;

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoBrick::sendSelf() - failed to send double data\n";
    return -1;
  }

  for (int i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LadrunoBrick::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  if (theDamping[0]) {
    for (int i = 0; i < 8; i++) {
      res += theDamping[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "LadrunoBrick::sendSelf -- could not send Damping\n";
        return res;
      }
    }
  }

  return res;
}

int  LadrunoBrick::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static ID idData(29);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoBrick::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(24));

  static Vector dData(8);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoBrick::recvSelf() - failed to recv double data\n";
    return -1;
  }
  alphaM = dData(0);
  betaK  = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  b[0]   = dData(4);
  b[1]   = dData(5);
  b[2]   = dData(6);
  hourglassCoeff = dData(7);

  for (int i = 0; i < 8; i++)
    connectedExternalNodes(i) = idData(16 + i);

  int packed = idData(28);
  formulation   = static_cast<Formulation>(packed % 10);
  massType      = (packed / 10) % 10;
  hourglassType = static_cast<Hourglass>((packed / 100) % 10);

  // rebuild the geometry-method layer from its serialized method id.  // Ladruno
  int geomID = (packed / 1000) % 10;
  if (theGeom) delete theGeom;
  theGeom = SolidTransformation::create(geomID);
  if (theGeom == 0)
    theGeom = new SolidTransformationLinear();   // safe fallback (unknown id)

  if (materialPointers[0] == 0) {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 8);
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
        opserr << "LadrunoBrick::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
        return -1;
      }
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoBrick::recvSelf() - material " << i << " failed to recv itself\n";
        return res;
      }
    }
  } else {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + 8);
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (materialPointers[i] == 0) {
          opserr << "LadrunoBrick::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
          return -1;
        }
        materialPointers[i]->setDbTag(matDbTag);
      }
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoBrick::recvSelf() - material " << i << " failed to recv itself\n";
        return res;
      }
    }
  }

  int dmpTag = (int)idData(26);
  if (dmpTag) {
    for (int i = 0; i < 8; i++) {
      if (theDamping[i] == 0 || theDamping[i]->getClassTag() != dmpTag) {
        if (theDamping[i]) delete theDamping[i];
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "LadrunoBrick::recvSelf -- could not get a Damping\n";
          return -1;
        }
      }
      theDamping[i]->setDbTag((int)idData(27));
      res += theDamping[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LadrunoBrick::recvSelf -- could not receive Damping\n";
        return res;
      }
    }
  } else {
    for (int i = 0; i < 8; i++) {
      if (theDamping[i]) {
        delete theDamping[i];
        theDamping[i] = 0;
      }
    }
  }

  return res;
}

//----------------------------------------------------------------------
int
LadrunoBrick::displaySelf(Renderer &theViewer, int displayMode, float fact,
                          const char **modes, int numMode)
{
  static Matrix coords(8, 3);
  static Vector values(8);
  static Vector vN(3);

  for (int n = 0; n < 8; n++) {
    nodePointers[n]->getDisplayCrds(vN, fact, displayMode);
    for (int i = 0; i < 3; i++)
      coords(n, i) = vN(i);
  }

  if (displayMode < 3 && displayMode > 0) {
    int index = displayMode - 1;
    for (int n = 0; n < 8; n++) {
      const Vector &stressN = materialPointers[n]->getStress();
      values(n) = stressN(index);
    }
  } else if (displayMode < 0) {
    for (int n = 0; n < 8; n++)
      values(n) = 0.0;
  }

  return theViewer.drawCube(coords, values, this->getTag());
}

//----------------------------------------------------------------------
Response *
LadrunoBrick::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoBrick");
  output.attr("eleTag", this->getTag());
  for (int i = 1; i <= 8; i++) {
    sprintf(outputData, "node%d", i);
    output.attr(outputData, nodePointers[i - 1]->getTag());
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
    for (int i = 1; i <= 8; i++) {
      sprintf(outputData, "P1_%d", i); output.tag("ResponseType", outputData);
      sprintf(outputData, "P2_%d", i); output.tag("ResponseType", outputData);
      sprintf(outputData, "P3_%d", i); output.tag("ResponseType", outputData);
    }
    theResponse = new ElementResponse(this, 1, resid);

  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8) {
      output.tag("GaussPoint");
      output.attr("number", pointNum);
      theResponse = materialPointers[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
      output.endTag();
    }

  } else if (strcmp(argv[0], "stresses") == 0) {
    for (int i = 0; i < 8; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma33");
      output.tag("ResponseType", "sigma12");
      output.tag("ResponseType", "sigma23");
      output.tag("ResponseType", "sigma13");
      output.endTag();
      output.endTag();
    }
    theResponse = new ElementResponse(this, 3, Vector(48));

  } else if (strcmp(argv[0], "strains") == 0) {
    for (int i = 0; i < 8; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      output.tag("ResponseType", "eps11");
      output.tag("ResponseType", "eps22");
      output.tag("ResponseType", "eps33");
      output.tag("ResponseType", "eps12");
      output.tag("ResponseType", "eps23");
      output.tag("ResponseType", "eps13");
      output.endTag();
      output.endTag();
    }
    theResponse = new ElementResponse(this, 4, Vector(48));

  } else if (strcmp(argv[0], "stress3D6") == 0) {
    output.tag("GaussPoint");
    output.attr("number", 1);
    output.tag("NdMaterialOutput");
    output.attr("classType", materialPointers[0]->getClassTag());
    output.attr("tag", materialPointers[0]->getTag());
    output.tag("ResponseType", "sigma11");
    output.tag("ResponseType", "sigma22");
    output.tag("ResponseType", "sigma33");
    output.tag("ResponseType", "sigma12");
    output.tag("ResponseType", "sigma23");
    output.tag("ResponseType", "sigma13");
    output.endTag();
    output.endTag();
    theResponse = new ElementResponse(this, 6, Vector(6));

  } else if (strcmp(argv[0], "strain3D6") == 0) {
    output.tag("GaussPoint");
    output.attr("number", 1);
    output.tag("NdMaterialOutput");
    output.attr("classType", materialPointers[0]->getClassTag());
    output.attr("tag", materialPointers[0]->getTag());
    output.tag("ResponseType", "eps11");
    output.tag("ResponseType", "eps22");
    output.tag("ResponseType", "eps33");
    output.tag("ResponseType", "eps12");
    output.tag("ResponseType", "eps23");
    output.tag("ResponseType", "eps13");
    output.endTag();
    output.endTag();
    theResponse = new ElementResponse(this, 7, Vector(6));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}

int
LadrunoBrick::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());

  else if (responseID == 3) {
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      const Vector &sigma = materialPointers[i]->getStress();
      for (int j = 0; j < 6; j++) stresses(cnt++) = sigma(j);
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 4) {
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      const Vector &eps = materialPointers[i]->getStrain();
      for (int j = 0; j < 6; j++) stresses(cnt++) = eps(j);
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 6) {
    Vector tmpStress(6);
    for (int i = 0; i < 8; i++) {
      const Vector &sigma = materialPointers[i]->getStress();
      for (int j = 0; j < 6; j++) tmpStress(j) += sigma(j) * 0.125;
    }
    return eleInfo.setVector(tmpStress);

  } else if (responseID == 7) {
    Vector tmpStrain(6);
    for (int i = 0; i < 8; i++) {
      const Vector &eps = materialPointers[i]->getStrain();
      for (int j = 0; j < 6; j++) tmpStrain(j) += eps(j);
    }
    tmpStrain /= 8.0;
    return eleInfo.setVector(tmpStrain);
  }

  return -1;
}

int
LadrunoBrick::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  if (argc < 1) return -1;

  // damping  (// Ladruno: upstream Brick looped i<4 over an 8-entry array — bug.
  //           LadrunoBrick iterates all 8.)
  if (strstr(argv[0], "damp") != 0) {
    if (argc < 2) return -1;
    for (int i = 0; i < 8; i++) {
      if (theDamping[i]) {
        int dmpRes = theDamping[i]->setParameter(argv, argc, param);
        if (dmpRes != -1) res = dmpRes;
      }
    }
    return res;
  }

  // specific material point
  if (strstr(argv[0], "material") != 0) {
    if (argc < 3) return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8)
      return materialPointers[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  // all material points
  for (int i = 0; i < 8; i++) {
    int matRes = materialPointers[i]->setParameter(argv, argc, param);
    if (matRes != -1) res = matRes;
  }
  return res;
}

int
LadrunoBrick::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;

  if (parameterID == res)
    return -1;

  for (int i = 0; i < 8; i++)
    matRes = materialPointers[i]->updateParameter(parameterID, info);

  if (matRes != -1) res = matRes;
  return res;
}
