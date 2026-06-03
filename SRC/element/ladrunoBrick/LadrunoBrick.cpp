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

// LadrunoBrick — unified 8-node hexahedral solid element (Ladruno fork).
//
// Implements the small-strain std + bbar + uri + ssp formulations. True
// Simo-Rifai enhanced assumed strain ('-formulation eas') is reserved (ADR 19).
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
#include <FiniteStrainNDMaterial.h>       // Ladruno — v3 finite: setTrialF(F) seam
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Response.h>        // Ladruno — cached material "damage" query (Tier-A Kstab)
#include <Information.h>
#include <DummyStream.h>

// Ladruno — Tier-A damage-scaled hourglass stabilization: residual fraction of
// the elastic Kstab kept in a fully-damaged element so it never goes completely
// unstabilized (avoids spurious hourglassing in the cracked band). 1% per the
// 1-5% range in Ladruno_implementation/11_brick_asdconcrete_integration.md §3;
// tune against the "hourglassEnergy" report (§5).
static const double HG_DAMAGE_FLOOR = 0.01;

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

// Node natural coordinates, Brick order (nodes 0-3 = zeta- face, 4-7 = zeta+):
//   0(-,-,-) 1(+,-,-) 2(+,+,-) 3(-,+,-) 4(-,-,+) 5(+,-,+) 6(+,+,+) 7(-,+,+)
const double LadrunoBrick::natCoord[8][3] = {
  {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
  {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
};

// The 4 hourglass base vectors, in Belytschko's mode order (8.7.25):
//   h_alpha = [ eta*zeta, zeta*xi, xi*eta, xi*eta*zeta ]  (alpha = 0..3)
// i.e. nodal values of the bilinear/trilinear modes the constant centroid
// gradient cannot represent. The order matters for the assumed-strain shear
// subsets g^{12}/g^{13}/g^{23} (physical); it is immaterial for the stiffness
// control (which sums all four).
const double LadrunoBrick::hg0[4][8] = {
  { 1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0 },  // h1 = eta*zeta  (yz)
  { 1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  1.0, -1.0 },  // h2 = zeta*xi   (zx)
  { 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0 },  // h3 = xi*eta    (xy)
  {-1.0,  1.0, -1.0,  1.0,  1.0, -1.0,  1.0, -1.0 }   // h4 = xi*eta*zeta
};

const char *
LadrunoBrick::formulationName(Formulation f)
{
  switch (f) {
  case Formulation::STD:  return "std";
  case Formulation::BBAR: return "bbar";
  case Formulation::URI:  return "uri";
  case Formulation::SSP:  return "ssp";
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
   theGeom(new SolidTransformationLinear()),  // Ladruno — v1 identity geometry
   damageResponse(0),                          // Ladruno — Tier-A Kstab (built in setDomain)
   sspBnot(0), sspKstab(0), sspVol(0.0),      // Ladruno — ssp (built in setDomain)
   alpha(9), alphaCommit(9), easJ0inv(3, 3), easJ0det(0.0)   // Ladruno — eas state
{
  B.Zero();
  alpha.Zero();
  alphaCommit.Zero();

  for (int i = 0; i < 8; i++) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
    theDamping[i] = 0;
  }

  b[0] = 0.0; b[1] = 0.0; b[2] = 0.0;

  // Ladruno — viscous-hourglass dissipation accumulator
  hgDissipated = 0.0;
  hgPrevValid  = false;
  for (int i = 0; i < 24; i++) uPrevCommit[i] = 0.0;
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
                           Damping *damping,
                           int geomMethodID)
  :Element(tag, ELE_TAG_LadrunoBrick),
   connectedExternalNodes(8),
   formulation(form),
   hourglassType(hgType), hourglassCoeff(hgCoeff),
   applyLoad(0), load(0), Ki(0), massType(matype),
   theGeom(0),                                // Ladruno — set below from geomMethodID
   damageResponse(0),                          // Ladruno — Tier-A Kstab (built in setDomain)
   sspBnot(0), sspKstab(0), sspVol(0.0),      // Ladruno — ssp (built in setDomain)
   alpha(9), alphaCommit(9), easJ0inv(3, 3), easJ0det(0.0)   // Ladruno — eas state
{
  alpha.Zero();
  alphaCommit.Zero();

  // Ladruno — geometry-method layer: linear (default, identity) / finite (UL).
  theGeom = SolidTransformation::create(geomMethodID);
  if (theGeom == 0)
    theGeom = new SolidTransformationLinear();   // safe fallback (unknown id)

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

  // Ladruno — viscous-hourglass dissipation accumulator
  hgDissipated = 0.0;
  hgPrevValid  = false;
  for (int i = 0; i < 24; i++) uPrevCommit[i] = 0.0;

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

  if (sspBnot)  delete sspBnot;  // Ladruno — ssp
  if (sspKstab) delete sspKstab;

  if (damageResponse) delete damageResponse;   // Ladruno — Tier-A Kstab
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

  // ssp: condense the single-point stabilization once, from geometry + the
  // initial material tangent (SSP-style). Deterministic, so the receive side
  // rebuilds it here after recvSelf — nothing extra travels in sendSelf.  // Ladruno
  if (formulation == Formulation::SSP) {
    bool haveNodes = true;
    for (int i = 0; i < 8; i++)
      if (nodePointers[i] == 0) haveNodes = false;
    if (haveNodes) buildSSP();
  }

  // eas (true Simo-Rifai): cache the centroid Jacobian inverse / det once for the
  // (j0/j) sym[J0^-T E_i] enhanced-mode map. Geometry-only, so the receive side
  // rebuilds it here (nothing extra in sendSelf besides the committed alpha).
  // v1 is small strain ⇒ the reference geometry is fixed and J0 never changes.  // Ladruno
  if (formulation == Formulation::EAS) {
    bool haveNodes = true;
    for (int i = 0; i < 8; i++)
      if (nodePointers[i] == 0) haveNodes = false;
    if (haveNodes) buildEAStrue();
  }

  // Ladruno — Tier-A: cache the centroid material's "damage" query for the
  // STIFFNESS-stabilized single-point formulations (ssp, uri+stiffness), whose
  // constant elastic Kstab must degrade with damage under softening. Rebuilt
  // here so the receive side (recvSelf -> setDomain) re-caches against the freshly
  // brokered material. Materials that report no "damage" channel return 0 here,
  // and damageScale() falls back to the elastic Kstab. uri+viscous and the
  // multi-point forms (std/bbar/uri+physical) need no scaling, so we skip them.
  if (damageResponse) { delete damageResponse; damageResponse = 0; }
  const bool stiffnessStabilized =
      (formulation == Formulation::SSP) ||
      (formulation == Formulation::URI && hourglassType == Hourglass::STIFFNESS);
  if (stiffnessStabilized && materialPointers[0] != 0) {
    DummyStream dmgStream;
    const char *dmgArgv[1] = {"damage"};
    damageResponse = materialPointers[0]->setResponse(dmgArgv, 1, dmgStream);
  }

  this->DomainComponent::setDomain(theDomain);
}

//----------------------------------------------------------------------
// Tier-A damage-scaled hourglass stabilization (softening support). Returns the
// multiplier s in [floor,1] that degrades the constant elastic Kstab with the
// material's current damage; 1.0 when the material reports no "damage" channel
// (so std elastic/J2 behave exactly as before). See the header note + §3 of
// Ladruno_implementation/11_brick_asdconcrete_integration.md.  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick::damageScale(void)
{
  if (damageResponse == 0)
    return 1.0;                       // material has no "damage" channel -> elastic Kstab
  if (damageResponse->getResponse() < 0)
    return 1.0;
  Information &info = damageResponse->getInformation();
  const Vector *d = info.theVector;
  if (d == 0)
    return 1.0;
  double dmax = 0.0;
  for (int i = 0; i < d->Size(); i++)
    if ((*d)(i) > dmax) dmax = (*d)(i);   // ASDConcrete3D: max(d_tension, d_compression)
  double s = 1.0 - dmax;
  if (s < HG_DAMAGE_FLOOR) s = HG_DAMAGE_FLOOR;
  if (s > 1.0) s = 1.0;
  return s;
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

  // Ladruno — eas: commit the converged enhanced parameters (DB/parallel sends
  // happen post-commit, so alphaCommit is the field sendSelf ships).
  if (formulation == Formulation::EAS)
    alphaCommit = alpha;

  // Ladruno — accumulate viscous-hourglass dissipation. The FB viscous force
  // stores no energy, so we integrate the work done against it over committed
  // steps. The first commit (or the first after a parallel recv) only seeds the
  // baseline so no spurious increment is booked.
  if (formulation == Formulation::URI && hourglassType == Hourglass::VISCOUS) {
    if (hgPrevValid)
      hgDissipated += this->viscousHourglassIncrement();
    else
      hgPrevValid = true;
    for (int J = 0; J < 8; J++) {
      const Vector &ul = nodePointers[J]->getTrialDisp();
      for (int i = 0; i < 3; i++)
        uPrevCommit[3 * J + i] = ul(i);
    }
  }

  return success;
}

int  LadrunoBrick::revertToLastCommit(void)
{
  int success = 0;
  for (int i = 0; i < 8; i++)
    success += materialPointers[i]->revertToLastCommit();
  for (int i = 0; i < 8; i++)
    if (theDamping[i]) success += theDamping[i]->revertToLastCommit();
  if (formulation == Formulation::EAS)   // Ladruno — eas: restore enhanced params
    alpha = alphaCommit;
  return success;
}

int  LadrunoBrick::revertToStart(void)
{
  int success = 0;
  for (int i = 0; i < 8; i++)
    success += materialPointers[i]->revertToStart();
  for (int i = 0; i < 8; i++)
    if (theDamping[i]) success += theDamping[i]->revertToStart();

  // Ladruno — reset the viscous-hourglass dissipation accumulator
  hgDissipated = 0.0;
  hgPrevValid  = false;
  for (int i = 0; i < 24; i++) uPrevCommit[i] = 0.0;

  // Ladruno — eas: reset enhanced parameters to the unstrained state
  alpha.Zero();
  alphaCommit.Zero();

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
    if (hourglassType == Hourglass::PHYSICAL)
      formPhysical(1, true);
    else
      formUri(1, true);
    Ki = new Matrix(stiff);
    return *Ki;
  }

  if (formulation == Formulation::SSP) {
    formSSP(1, true);
    Ki = new Matrix(stiff);
    return *Ki;
  }

  if (formulation == Formulation::EAS) {
    formEAStrue(1, true);   // condensed elastic K* at alpha=0 (no inner Newton)
    Ki = new Matrix(stiff);
    return *Ki;
  }

  // -geom finite (incl. F-bar) reaches the small-strain assembly below ON
  // PURPOSE. The initial stiffness is the tangent at the REFERENCE state
  // (u = 0 ⇒ F = I ⇒ σ = 0, no geometric term, material modulus = initial
  // tangent), so for std the small-strain Bᵀ D_init B with reference gradients
  // IS that finite tangent. We deliberately do NOT assemble the F-bar centroid
  // coupling here: getInitialStiff feeds Rayleigh βK0, -initial / Krylov seeding
  // and eigen analysis, which want a SYMMETRIC, well-conditioned seed — not the
  // generally-UNSYMMETRIC F-bar tangent (the bbar branch below gives the
  // mean-dilatation B-bar seed). The consistent (F-bar) tangent that drives the
  // converged solution comes from getTangentStiff → formResidAndTangentFinite.  // Ladruno

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
  // Refresh the geometry frame to the REFERENCE configuration (cur == ref ⇒
  // R = I) so the initial stiffness is the un-rotated core tangent regardless
  // of when it is queried — never a stale current-frame R (sweep #7).  // Ladruno
  {
    static Matrix refC(8, 3), curC(8, 3);
    for (int i = 0; i < 8; i++) {
      const Vector &X = nodePointers[i]->getCrds();
      for (int d = 0; d < 3; d++) { refC(i, d) = X(d); curC(i, d) = X(d); }
    }
    theGeom->update(8, refC, curC);
  }
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

// True for the centroid-only formulations (ssp, uri+stiffness, uri+viscous):
// the constitutive update runs once on material slot 0; slots 1-7 are output
// mirrors. std/bbar/uri-physical and -geom finite use all 8 GPs.  // Ladruno
bool
LadrunoBrick::isSinglePoint(void) const
{
  if (this->isFinite()) return false;
  if (formulation == Formulation::SSP) return true;
  if (formulation == Formulation::URI && hourglassType != Hourglass::PHYSICAL) return true;
  return false;
}

//update — seam 1 (kinematics ledger): set the material trial strain
int
LadrunoBrick::update(void)
{
  if (this->isFinite())          // -geom finite (v3): F-driven, updated-Lagrangian
    return this->updateFinite();

  // eas (true Simo-Rifai): set the 8 material trial states here (the standard
  // update() contract). formEAStrue(0,...) solves the enhanced parameters alpha by
  // an inner Newton and strains each GP to B*u + M*alpha, so commitState captures
  // the correct converged state under ANY algorithm — including a single-pass
  // 'Linear' step, where no force/tangent form runs at the final u. (The residual
  // it also assembles is harmless; getResistingForce recomputes it.)  // Ladruno
  if (formulation == Formulation::EAS) {
    formEAStrue(0, false);
    return 0;
  }

  if (formulation == Formulation::SSP) {
    if (sspBnot == 0) buildSSP();   // safety (normally built in setDomain)
    computeBasis();
    const Vector &uCore = this->computeLocalDisp();   // identity for linear
    static Vector strainE(6);
    strainE.addMatrixVector(0.0, *sspBnot, uCore, 1.0);   // strain = Bnot * u
    // Single integration point: evaluate ONLY material slot 0 (the centroid).
    // formSSP reads only slot 0; slots 1-7 are output mirrors (getResponse /
    // setResponse map per-GP queries to slot 0 via isSinglePoint()) so we avoid
    // 7 redundant return maps per element — material cost matters for e.g.
    // ASDConcrete3D. See LEDGER/11_brick_asdconcrete_integration §4.  // Ladruno
    materialPointers[0]->setTrialStrain(strainE);
    return 0;
  }

  if (formulation == Formulation::URI) {
    computeBasis();
    const Vector &uCore = this->computeLocalDisp();   // seam 0+2 (identity for linear)

    if (hourglassType == Hourglass::PHYSICAL) {
      // physical: assumed-strain at the full 2x2x2 rule, one material per GP.
      double gp0[3] = {0.0, 0.0, 0.0};
      double detJ0;
      static double shpC[4][8];
      shp3d(gp0, detJ0, shpC, xl);
      double bC[3][8];
      for (int i = 0; i < 3; i++)
        for (int I = 0; I < 8; I++)
          bC[i][I] = shpC[i][I];
      double gamma[4][8];
      computeGammaHourglass(bC, gamma, 0.125);

      static double Bbar[8][6][3];
      static Vector strainG(6);
      int gpIdx = 0;
      for (int gi = 0; gi < 2; gi++)
        for (int gj = 0; gj < 2; gj++)
          for (int gk = 0; gk < 2; gk++) {
            double gp[3] = {sg[gi], sg[gj], sg[gk]};
            assumedStrainB(gp, gamma, bC, Bbar);
            strainG.Zero();
            for (int J = 0; J < 8; J++) {        // uCore = localized (core-frame) disp
              for (int r = 0; r < 6; r++)
                for (int c = 0; c < 3; c++)
                  strainG(r) += Bbar[J][r][c] * uCore(3 * J + c);
            }
            materialPointers[gpIdx++]->setTrialStrain(strainG);
          }
      return 0;
    }

    // uri (perturbation hourglass): single centroid integration point. Evaluate
    // ONLY material slot 0; formUri reads only slot 0 and slots 1-7 are output
    // mirrors (isSinglePoint()) — avoids 7 redundant return maps.  // Ladruno
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
    materialPointers[0]->setTrialStrain(strainC);
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
  if (this->isFinite()) {        // -geom finite (v3): updated-Lagrangian assembly
    this->formResidAndTangentFinite(tang_flag);
    return;
  }

  if (formulation == Formulation::URI) {
    if (hourglassType == Hourglass::PHYSICAL)
      formPhysical(tang_flag, false);
    else
      formUri(tang_flag, false);
    return;
  }

  if (formulation == Formulation::SSP) {
    formSSP(tang_flag, false);
    return;
  }

  if (formulation == Formulation::EAS) {
    formEAStrue(tang_flag, false);
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
  static Vector bodyForce(numberNodes * ndf);   // external (global-frame) body load
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
  bodyForce.Zero();
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
        // external body/self-weight load is a GLOBAL-frame dead load: accumulate it
        // OUT of the core force so globalizeForce (corot) rotates only the internal
        // force. Added back unrotated, in the global frame, after globalize.  // Ladruno (COROT-1)
        if (applyLoad == 0)
          bodyForce(jj + p) -= dvol[i] * b[p] * shp[3][j];
        else
          bodyForce(jj + p) -= dvol[i] * appliedB[p] * shp[3][j];
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
  resid += bodyForce;            // global-frame body load, added AFTER globalize  // Ladruno (COROT-1)
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

//======================================================================
// -geom finite (v3): updated-Lagrangian finite-strain path.
//
// The material is a FiniteStrainNDMaterial driven by setTrialF(F): it returns
// the Cauchy stress σ (getStress) and the spatial CONSTITUTIVE modulus c
// (getTangent). Per the LOCKED seam-3 contract the element owns the geometric
// stiffness, so the element forms the FULL consistent tangent on the current
// configuration:
//     f = ∫ Bᵀ σ dv ,   K = ∫ Bᵀ c B dv + ∫ Gᵀ Σ G dv ,   dv = J dV
// with spatial gradients ∂Nₐ/∂xⱼ = Σ_k (∂Nₐ/∂X_k)(F⁻¹)_kj and the geometric
// (initial-stress) term Σ built from σ.
//======================================================================
bool
LadrunoBrick::isFinite(void) const
{
  return theGeom != 0 &&
         theGeom->getStrainMeasure() ==
           SolidTransformation::StrainMeasure::DeformationGradient;
}

// F (row-major [9]) from reference shape gradients shpRef[i][a]=∂Nₐ/∂Xᵢ and the
// nodal trial displacements: F_iJ = δ_iJ + Σ_a u_a[i] ∂N_a/∂X_J. Returns det F.
double
LadrunoBrick::deformationGradient(const double shpRef[4][8], double F[9])
{
  for (int i = 0; i < 9; i++) F[i] = 0.0;
  F[0] = F[4] = F[8] = 1.0;
  for (int a = 0; a < 8; a++) {
    const Vector &ua = nodePointers[a]->getTrialDisp();
    for (int i = 0; i < 3; i++)
      for (int J = 0; J < 3; J++)
        F[3 * i + J] += ua(i) * shpRef[J][a];
  }
  return F[0]*(F[4]*F[8]-F[5]*F[7]) - F[1]*(F[3]*F[8]-F[5]*F[6])
       + F[2]*(F[3]*F[7]-F[4]*F[6]);
}

// Inverse of a row-major 3×3 via adjugate/det. Returns det F; fills Finv (row-
// major), zero-filled if |det| underflows to 0. Shared by the finite/F-bar GP
// and centroid spatial-gradient pushes (one source of truth for the 9-term
// adjugate, which is exactly the kind of formula that drifts under copy-paste).
static double
invert3x3(const double F[9], double Finv[9])
{
  double det = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[1]*(F[3]*F[8]-F[5]*F[6])
             + F[2]*(F[3]*F[7]-F[4]*F[6]);
  double id = (det != 0.0) ? 1.0 / det : 0.0;
  Finv[0] =  (F[4]*F[8]-F[5]*F[7])*id;  Finv[1] = -(F[1]*F[8]-F[2]*F[7])*id;
  Finv[2] =  (F[1]*F[5]-F[2]*F[4])*id;  Finv[3] = -(F[3]*F[8]-F[5]*F[6])*id;
  Finv[4] =  (F[0]*F[8]-F[2]*F[6])*id;  Finv[5] = -(F[0]*F[5]-F[2]*F[3])*id;
  Finv[6] =  (F[3]*F[7]-F[4]*F[6])*id;  Finv[7] = -(F[0]*F[7]-F[1]*F[6])*id;
  Finv[8] =  (F[0]*F[4]-F[1]*F[3])*id;
  return det;
}

// F-bar centroid data (bbar + finite). Assumes computeBasis() has set xl.
// Returns J0 = det F0, F0 = deformation gradient at the element centroid
// (ξ=0), per dSNPO eq 15.5. If G0 != 0, also fills the centroid spatial gradient
// operator G0[k][b] = ∂N_b/∂x_k|_centroid (from F0⁻¹) for the eq 15.10 coupling.
double
LadrunoBrick::centroidFbar(double (*G0)[8])
{
  static double shp0[4][8];
  double detJ0c;
  double pt0[3] = { 0.0, 0.0, 0.0 };
  shp3d(pt0, detJ0c, shp0, xl);            // shp0[0..2][a] = ∂Nₐ/∂X at the centroid

  double F0[9];
  double J0 = deformationGradient(shp0, F0);

  if (G0 != 0) {
    double Fi[9];
    invert3x3(F0, Fi);                       // F0⁻¹ (row-major)
    for (int b = 0; b < 8; b++)
      for (int kk = 0; kk < 3; kk++) {
        double s = 0.0;
        for (int m = 0; m < 3; m++)
          s += shp0[m][b] * Fi[3 * m + kk];      // ∂N_b/∂x_kk at the centroid
        G0[kk][b] = s;
      }
  }
  return J0;
}

// update — finite: per GP compute F and drive the material via setTrialF(F).
// bbar + finite = F-bar: drive the material with Fbar = (J0/J)^(1/3) F (dSNPO
// eq 15.5) so every GP shares the centroid dilatation J0 (volumetric-locking
// cure); std finite uses F directly.
int
LadrunoBrick::updateFinite(void)
{
  computeBasis();                          // xl = reference nodal coordinates

  static double shp[4][8];
  static Matrix Fm(3, 3);
  double detJ0;

  const bool useFbar = (formulation == Formulation::BBAR);
  double J0 = 1.0;
  if (useFbar) {
    J0 = this->centroidFbar();             // det F0 at the centroid
    if (J0 <= 0.0) {
      opserr << "LadrunoBrick::updateFinite - non-positive centroid det F0 = "
             << J0 << " (element " << this->getTag() << ", F-bar)\n";
      return -1;
    }
  }

  int gp = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        double pt[3] = { sg[i], sg[j], sg[k] };
        shp3d(pt, detJ0, shp, xl);         // shp[0..2][a] = ∂Nₐ/∂X
        double F[9];
        double J = deformationGradient(shp, F);
        if (useFbar) {
          if (J <= 0.0) {
            opserr << "LadrunoBrick::updateFinite - non-positive det F = " << J
                   << " at GP " << gp << " (element " << this->getTag()
                   << ", F-bar)\n";
            return -1;
          }
          double s = pow(J0 / J, 1.0 / 3.0);     // Fbar = (J0/J)^(1/3) F
          for (int n = 0; n < 9; n++) F[n] *= s;
        }
        for (int r = 0; r < 3; r++)
          for (int c = 0; c < 3; c++)
            Fm(r, c) = F[3 * r + c];

        FiniteStrainNDMaterial *fsm =
          dynamic_cast<FiniteStrainNDMaterial *>(materialPointers[gp]);
        if (fsm == 0) {
          opserr << "LadrunoBrick::updateFinite - material at GP " << gp
                 << " is not a FiniteStrainNDMaterial (element " << this->getTag()
                 << ")\n";
          return -1;
        }
        if (fsm->setTrialF(Fm) < 0) {
          opserr << "LadrunoBrick::updateFinite - setTrialF failed at GP " << gp
                 << " (element " << this->getTag() << ", det F<=0?)\n";
          return -1;
        }
        gp++;
      }
  return 0;
}

// formResidAndTangentFinite — updated-Lagrangian assembly.
//
//   f_{a,i} = ∫ σ_ij ∂Nₐ/∂x_j dv                         (current config, dv=J dV)
//   K_{(a,i)(b,k)} = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv
//
// with the FULL 4th-order spatial consistent tangent (dSNPO eq.14.99)
//   a_ijkl = c_ijkl − σ_il δ_jk ,
// where c_ijkl = (1/2J)[D:L:B]_ijkl is the material's getTangent() (returned as
// the minor-symmetric 6×6 we expand via ladrunoVI) and the −σ_il δ_jk term is
// the element-owned geometric/initial-stress contribution (the LOCKED seam-3
// contract). The geometric term is NON-minor-symmetric, so it must enter the
// full-gradient contraction above — it is NOT the conventional ∫ Gᵀ Σ G alone,
// and using c without it gives a non-symmetric, wrong tangent.
void
LadrunoBrick::formResidAndTangentFinite(int tang_flag)
{
  static const int numberNodes = 8;

  stiff.Zero();
  resid.Zero();
  computeBasis();                          // reference coordinates

  static double shpRef[4][8];              // [0..2]=∂Nₐ/∂X, [3]=Nₐ
  static double g[3][8];                   // spatial gradients ∂Nₐ/∂x_j
  static Vector stress(6);

  // F-bar (bbar + finite): the eq 15.10 tangent gains a (generally UNSYMMETRIC)
  // coupling to the element-centroid gradient operator G0. The residual is
  // unchanged — F-bar enters it only through σ̄ (set in updateFinite, eq 15.9 /
  // Remark 15.2). G0 is element-wide; compute it once.  // Ladruno
  const bool useFbar = (formulation == Formulation::BBAR);
  static double G0[3][8];
  if (useFbar && tang_flag == 1)
    this->centroidFbar(G0);

  int gp = 0;
  for (int gi = 0; gi < 2; gi++)
   for (int gj = 0; gj < 2; gj++)
    for (int gk = 0; gk < 2; gk++) {
      double pt[3] = { sg[gi], sg[gj], sg[gk] };
      double detJ0;
      shp3d(pt, detJ0, shpRef, xl);

      double F[9];
      double J = deformationGradient(shpRef, F);

      double Fi[9];
      invert3x3(F, Fi);                           // F⁻¹ (row-major)

      for (int a = 0; a < numberNodes; a++)
        for (int j = 0; j < 3; j++) {
          double s = 0.0;
          for (int k = 0; k < 3; k++)
            s += shpRef[k][a] * Fi[3 * k + j];      // ∂Nₐ/∂x_j
          g[j][a] = s;
        }

      double dv = J * detJ0;                         // current volume (wg = 1)

      stress = materialPointers[gp]->getStress();    // Cauchy σ (set via setTrialF)
      double sig[3][3];
      sig[0][0]=stress(0); sig[1][1]=stress(1); sig[2][2]=stress(2);
      sig[0][1]=sig[1][0]=stress(3);
      sig[1][2]=sig[2][1]=stress(4);
      sig[2][0]=sig[0][2]=stress(5);

      // residual: f_{a,i} = ∫ σ_ij ∂Nₐ/∂x_j dv  (− body force). Body force uses
      // the REFERENCE volume detJ0 (dead load per reference volume), consistent
      // with the reference-config mass in formInertiaTerms.
      for (int a = 0; a < numberNodes; a++)
        for (int i = 0; i < 3; i++) {
          double fi = 0.0;
          for (int j = 0; j < 3; j++)
            fi += sig[i][j] * g[j][a];
          resid(3 * a + i) += fi * dv;
          if (applyLoad == 0)
            resid(3 * a + i) -= detJ0 * b[i] * shpRef[3][a];
          else
            resid(3 * a + i) -= detJ0 * appliedB[i] * shpRef[3][a];
        }

      if (tang_flag == 1) {
        // FULL 4th-order spatial constitutive modulus c_ijkl (the material's 6×6
        // getTangent is lossy in (k,l); the consistent tangent needs the full
        // tensor — see FiniteStrainNDMaterial.h). Then a_ijkl = c_ijkl − σ_il δ_jk
        // (the −σδ is the element-owned geometric term), contracted with the full
        // nodal gradients.  // Ladruno
        FiniteStrainNDMaterial *fsm =
          dynamic_cast<FiniteStrainNDMaterial *>(materialPointers[gp]);
        static double cmat[3][3][3][3];
        if (fsm == 0 || fsm->getSpatialTangentTensor(cmat) != 0) {
          opserr << "LadrunoBrick::formResidAndTangentFinite - material at GP " << gp
                 << " did not provide a full spatial tangent (getSpatialTangentTensor); "
                    "element " << this->getTag() << "\n";
          return;
        }
        double a4[3][3][3][3];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
              for (int l = 0; l < 3; l++)
                a4[i][j][k][l] = cmat[i][j][k][l]
                               - sig[i][l] * (j == k ? 1.0 : 0.0);

        for (int a = 0; a < numberNodes; a++)
          for (int bn = 0; bn < numberNodes; bn++)
            for (int i = 0; i < 3; i++)
              for (int k = 0; k < 3; k++) {
                double s = 0.0;
                for (int j = 0; j < 3; j++)
                  for (int l = 0; l < 3; l++)
                    s += g[j][a] * a4[i][j][k][l] * g[l][bn];
                stiff(3 * a + i, 3 * bn + k) += s * dv;
              }

        // F-bar additional (generally UNSYMMETRIC) stiffness, dSNPO eq 15.10:
        //   K_{(a,i)(b,k)} += ∫ (Σ_j g[j][a] q_ij)(G0[k][b] − g[k][b]) dv,
        // with q the matrix form of the eq 15.11 tensor evaluated at F=F̄,
        //   q_ij = (1/3) a_ijpp − (2/3) σ̄_ij ,
        // using the SAME a4 = c̄ − σ̄δ modulus as the standard term above (NOT the
        // material part c̄ alone). (G0 − g) is the centroid-minus-GP gradient: the
        // term vanishes when G0 = g, collapsing F-bar to the plain-F tangent.  // Ladruno
        if (useFbar) {
          double M[3][3];
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
              M[i][j] = (a4[i][j][0][0] + a4[i][j][1][1] + a4[i][j][2][2]) / 3.0
                      - (2.0 / 3.0) * sig[i][j];

          static double Lfac[8][3];        // Lfac[a][i] = Σ_j g[j][a] q_ij
          for (int a = 0; a < numberNodes; a++)
            for (int i = 0; i < 3; i++) {
              double s = 0.0;
              for (int j = 0; j < 3; j++) s += g[j][a] * M[i][j];
              Lfac[a][i] = s;
            }

          for (int a = 0; a < numberNodes; a++)
            for (int i = 0; i < 3; i++)
              for (int bn = 0; bn < numberNodes; bn++)
                for (int kk = 0; kk < 3; kk++)
                  stiff(3 * a + i, 3 * bn + kk)
                    += Lfac[a][i] * (G0[kk][bn] - g[kk][bn]) * dv;
        }
      }
      gp++;
    }
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
LadrunoBrick::computeGammaHourglass(const double bC[3][8], double gamma[4][8], double beta)
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
      gamma[a][I] = beta * (hg0[a][I] - corr);
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

  static Vector bodyForce(numberNodes * ndf);   // external (global-frame) body load
  stiff.Zero();
  if (!useInitialTangent) {
    resid.Zero();
    bodyForce.Zero();
  }

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
  double kappa = scale * Ghg * vol * bb;

  // Ladruno — Tier-A: for a softening material (a "damage" channel is present)
  // base the STIFFNESS hourglass stiffness on the ELASTIC shear and degrade it
  // explicitly by a floored damage scale, so a fully-cracked element keeps
  // residual hourglass control (κ → floor·κ_elastic, not 0). Using the current
  // secant shear (Ghg) here would let κ → 0 at full damage and defeat the floor.
  // The VISCOUS flavour adds damping not stiffness, so it is NOT degraded (full
  // hourglass control in the cracked band is desirable there). Non-softening
  // materials (no damage channel) keep the existing current-tangent κ unchanged.
  if (!useInitialTangent && hourglassType == Hourglass::STIFFNESS && damageResponse != 0) {
    const double G0 = materialPointers[0]->getInitialTangent()(3, 3);
    kappa = scale * G0 * vol * bb * this->damageScale();
  }

  double gamma[4][8];
  computeGammaHourglass(bC, gamma, 1.0);   // stiffness: beta absorbed into kappa

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
        // global-frame body load kept out of the corotated core force  // Ladruno (COROT-1)
        if (applyLoad == 0)
          bodyForce(jj + p) -= vol * b[p] * shpC[3][j];      // b = member body force; Nc = shpC[3][j]
        else
          bodyForce(jj + p) -= vol * appliedB[p] * shpC[3][j];
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

  if (hourglassType == Hourglass::VISCOUS) {
    // --- Flanagan-Belytschko VISCOUS hourglass control (rate form, explicit) ---
    // A damping force resisting the hourglass *velocity* (Belytschko 8.7.10).
    // It adds NO stiffness (velocity-dependent), so the tangent keeps the 12
    // hourglass modes at zero energy — fine for explicit dynamics (they are
    // damped, and dt_cr depends only on the populated max eigenvalue). Vanishes
    // in statics (v=0). Coefficient c_visc = eps * rho * c_d * V^(2/3), with the
    // dilatational wave speed c_d = sqrt((lambda+2mu)/rho) ~ sqrt(dd(0,0)/rho).
    if (!useInitialTangent) {
      double rho = materialPointers[0]->getRho();
      if (rho > 0.0) {
        const double cd = sqrt(dd(0, 0) / rho);
        const double epsV = (hourglassCoeff > 0.0) ? hourglassCoeff : 0.1;
        const double cVisc = epsV * rho * cd * pow(vol, 2.0 / 3.0);
        // generalized hourglass velocities qdot_ai = sum_J gamma_aJ v_iJ
        double qd[4][3];
        for (int a = 0; a < 4; a++)
          for (int i = 0; i < 3; i++)
            qd[a][i] = 0.0;
        for (int J = 0; J < numberNodes; J++) {
          const Vector &vl = nodePointers[J]->getTrialVel();
          for (int a = 0; a < 4; a++)
            for (int i = 0; i < 3; i++)
              qd[a][i] += gamma[a][J] * vl(i);
        }
        // hourglass damping force f_iI = c_visc * sum_a gamma_aI qdot_ai
        for (int I = 0; I < numberNodes; I++)
          for (int i = 0; i < 3; i++) {
            double f = 0.0;
            for (int a = 0; a < 4; a++)
              f += gamma[a][I] * qd[a][i];
            resid(3 * I + i) += cVisc * f;
          }
      }
    }
  } else {
    // --- Flanagan-Belytschko STIFFNESS hourglass control (default) ---
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
  }

  // seam 3: globalize the core-frame f/K back to global DOFs, adding K_geo.
  // identity for -geom linear. resid is assembled only when !useInitialTangent;
  // the initial-stiffness call uses zero internal force (zero stress => no K_geo).
  static Vector zeroF(24);                                                  // Ladruno
  const Vector &fCore = useInitialTangent ? zeroF : resid;
  if (tang_flag == 1)
    theGeom->globalizeStiff(stiff, fCore, stiff);
  if (!useInitialTangent) {
    theGeom->globalizeForce(resid, resid);
    resid += bodyForce;          // global-frame body load, added AFTER globalize  // Ladruno (COROT-1)
  }
}

//----------------------------------------------------------------------
// Belytschko-Bindeman isochoric assumed-strain B-bar at a Gauss point (eq
// 8.7.26). Returns |J|. Bbar[I] is 6x3 in Voigt {xx,yy,zz,xy,yz,zx}. The
// normal-strain rows carry the deviatoric-projected hourglass dilatation
// (2/3,-1/3 split -> isochoric, no volumetric locking); the shear rows use the
// mode-subset hourglass vectors g^{12}/g^{13}/g^{23} (selective reduction ->
// no shear locking). Constant parts equal the compatible centroid B, so the
// assumed field reproduces constant strain exactly (patch test holds).
//----------------------------------------------------------------------
double
LadrunoBrick::assumedStrainB(const double gp[3], const double gamma[4][8],
                             const double bC[3][8], double Bbar[8][6][3])
{
  const double xi = gp[0], eta = gp[1], zeta = gp[2];

  // analytic dN_I/dxi_k and the Jacobian A_ik = dx_i/dxi_k
  double dNdxi[3][8];
  for (int I = 0; I < 8; I++) {
    const double x0 = natCoord[I][0], y0 = natCoord[I][1], z0 = natCoord[I][2];
    dNdxi[0][I] = 0.125 * x0 * (1.0 + y0 * eta) * (1.0 + z0 * zeta);
    dNdxi[1][I] = 0.125 * y0 * (1.0 + x0 * xi) * (1.0 + z0 * zeta);
    dNdxi[2][I] = 0.125 * z0 * (1.0 + x0 * xi) * (1.0 + y0 * eta);
  }
  double A[3][3];
  for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++) {
      double s = 0.0;
      for (int I = 0; I < 8; I++) s += dNdxi[k][I] * xl[i][I];
      A[i][k] = s;
    }
  const double det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                   - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                   + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
  const double idet = 1.0 / det;
  double Ainv[3][3];   // Ainv[k][i] = dxi_k/dx_i
  Ainv[0][0] = (A[1][1]*A[2][2]-A[1][2]*A[2][1])*idet;
  Ainv[0][1] = (A[0][2]*A[2][1]-A[0][1]*A[2][2])*idet;
  Ainv[0][2] = (A[0][1]*A[1][2]-A[0][2]*A[1][1])*idet;
  Ainv[1][0] = (A[1][2]*A[2][0]-A[1][0]*A[2][2])*idet;
  Ainv[1][1] = (A[0][0]*A[2][2]-A[0][2]*A[2][0])*idet;
  Ainv[1][2] = (A[0][2]*A[1][0]-A[0][0]*A[1][2])*idet;
  Ainv[2][0] = (A[1][0]*A[2][1]-A[1][1]*A[2][0])*idet;
  Ainv[2][1] = (A[0][1]*A[2][0]-A[0][0]*A[2][1])*idet;
  Ainv[2][2] = (A[0][0]*A[1][1]-A[0][1]*A[1][0])*idet;

  // hourglass-mode spatial derivatives h_{alpha,i} = dh_alpha/dxi_k * dxi_k/dx_i
  const double dh[4][3] = {
    {0.0,      zeta,     eta},        // h1 = eta*zeta
    {zeta,     0.0,      xi},         // h2 = zeta*xi
    {eta,      xi,       0.0},        // h3 = xi*eta
    {eta*zeta, xi*zeta,  xi*eta}      // h4 = xi*eta*zeta
  };
  double hi[4][3];
  for (int a = 0; a < 4; a++)
    for (int i = 0; i < 3; i++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += dh[a][k] * Ainv[k][i];
      hi[a][i] = s;
    }

  for (int I = 0; I < 8; I++) {
    const double bx = bC[0][I], by = bC[1][I], bz = bC[2][I];
    double g1234[3], g12[3], g13[3], g23[3];
    for (int i = 0; i < 3; i++) {
      const double a0 = hi[0][i] * gamma[0][I];   // mode 1 (yz)
      const double a1 = hi[1][i] * gamma[1][I];   // mode 2 (zx)
      const double a2 = hi[2][i] * gamma[2][I];   // mode 3 (xy)
      const double a3 = hi[3][i] * gamma[3][I];   // mode 4 (xyz)
      g1234[i] = a0 + a1 + a2 + a3;
      g12[i]   = a0 + a1;
      g13[i]   = a0 + a2;
      g23[i]   = a1 + a2;
    }
    const double gx = g1234[0], gy = g1234[1], gz = g1234[2];
    double (*B)[3] = Bbar[I];
    // normal strains: FULL compatible (no deviatoric projection). The pointwise-
    // isochoric projection of eq 8.7.26 strips volumetric bending stiffness and
    // over-softens (worse with nu); SSPbrick-style accuracy needs the full
    // bending normal strain. Volumetric locking is a separate axis (use bbar).
    B[0][0] = bx + gx;  B[0][1] = 0.0;       B[0][2] = 0.0;
    B[1][0] = 0.0;      B[1][1] = by + gy;   B[1][2] = 0.0;
    B[2][0] = 0.0;      B[2][1] = 0.0;       B[2][2] = bz + gz;
    // shears (Voigt xy,yz,zx) with mode-subset hourglass vectors (eq 8.7.26) —
    // the selective reduction that relieves shear locking
    B[3][0] = by;                 B[3][1] = bx + g23[0];         B[3][2] = g23[0];         // xy  <- g^{23}
    B[4][0] = g12[1];             B[4][1] = bz + g12[2];         B[4][2] = by;             // yz  <- g^{12}
    B[5][0] = bz + g13[2];        B[5][1] = 0.0;                 B[5][2] = bx + g13[0];    // zx  <- g^{13}
  }
  return det;
}

//----------------------------------------------------------------------
// physical (uri + assumed-strain) — full 2x2x2 integration of an assumed strain
// (8.7.20): FULL compatible normal strains + Belytschko-Bindeman mode-subset
// transverse shear (8.7.26). Rank-sufficient (no perturbation hourglass),
// tuning-free, and a verified SHEAR-locking cure (matches SSPbrick to 3 digits
// and converges in bending at compressible/moderate nu). NOTE: it does NOT cure
// VOLUMETRIC locking — near-incompressible (nu -> 0.5) it locks; use -formulation
// bbar there. (The pointwise-isochoric dev-projection of eq 8.7.26 would relieve
// volumetric locking but, combined with the reduced shear, over-softens bending
// for all nu; no single static projection is correct across nu with this shear
// field -- that needs the coupled SSP/ASQBI operator, future work.)
//----------------------------------------------------------------------
void
LadrunoBrick::formPhysical(int tang_flag, bool useInitialTangent)
{
  static const int ndf = 3, nstress = 6, numberNodes = 8;

  computeBasis();

  // centroid gradients + FB hourglass vectors (beta = 1/8: gamma is the true
  // hourglass-mode coefficient so the assumed strain matches the compatible one)
  double gp0[3] = {0.0, 0.0, 0.0};
  double detJ0;
  static double shpC[4][8];
  shp3d(gp0, detJ0, shpC, xl);
  double bC[3][8];
  for (int i = 0; i < 3; i++)
    for (int I = 0; I < numberNodes; I++)
      bC[i][I] = shpC[i][I];
  double gamma[4][8];
  computeGammaHourglass(bC, gamma, 0.125);

  stiff.Zero();
  if (!useInitialTangent)
    resid.Zero();

  static double Bbar[8][6][3];
  static Matrix dd(nstress, nstress), BJ(nstress, ndf), BJtran(ndf, nstress);
  static Matrix BK(nstress, ndf), BJtranD(ndf, nstress), stiffJK(ndf, ndf);
  static Vector stress(nstress), residJ(ndf);

  int gpIdx = 0;
  for (int gi = 0; gi < 2; gi++) {
    for (int gj = 0; gj < 2; gj++) {
      for (int gk = 0; gk < 2; gk++) {
        double gp[3] = {sg[gi], sg[gj], sg[gk]};
        double detJ = assumedStrainB(gp, gamma, bC, Bbar);
        double dvol = wg[gpIdx] * detJ;

        double N[8];
        for (int I = 0; I < numberNodes; I++)
          N[I] = 0.125 * (1.0 + natCoord[I][0]*gp[0])
                       * (1.0 + natCoord[I][1]*gp[1])
                       * (1.0 + natCoord[I][2]*gp[2]);

        if (tang_flag == 1) {
          dd = useInitialTangent ? materialPointers[gpIdx]->getInitialTangent()
                                 : materialPointers[gpIdx]->getTangent();
          dd *= dvol;
        }
        if (!useInitialTangent) {
          stress = materialPointers[gpIdx]->getStress();
          stress *= dvol;
        }

        int jj = 0;
        for (int j = 0; j < numberNodes; j++) {
          for (int r = 0; r < nstress; r++)
            for (int c = 0; c < ndf; c++)
              BJ(r, c) = Bbar[j][r][c];
          for (int p = 0; p < ndf; p++)
            for (int q = 0; q < nstress; q++)
              BJtran(p, q) = BJ(q, p);

          if (!useInitialTangent) {
            residJ.addMatrixVector(0.0, BJtran, stress, 1.0);
            for (int p = 0; p < ndf; p++) {
              resid(jj + p) += residJ(p);
              if (applyLoad == 0)
                resid(jj + p) -= dvol * b[p] * N[j];
              else
                resid(jj + p) -= dvol * appliedB[p] * N[j];
            }
          }

          if (tang_flag == 1) {
            BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);
            int kk = 0;
            for (int k = 0; k < numberNodes; k++) {
              for (int r = 0; r < nstress; r++)
                for (int c = 0; c < ndf; c++)
                  BK(r, c) = Bbar[k][r][c];
              stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);
              for (int p = 0; p < ndf; p++)
                for (int q = 0; q < ndf; q++)
                  stiff(jj + p, kk + q) += stiffJK(p, q);
              kk += ndf;
            }
          }
          jj += ndf;
        }
        gpIdx++;
      }
    }
  }
}

//----------------------------------------------------------------------
// ssp (Stabilized Single-Point, SSPbrick port) — cross product helper.
static Vector
ladrunoCross(const Vector &v1, const Vector &v2)
{
  Vector r(3);
  r(0) = v1(1) * v2(2) - v1(2) * v2(1);
  r(1) = v1(2) * v2(0) - v1(0) * v2(2);
  r(2) = v1(0) * v2(1) - v1(1) * v2(0);
  return r;
}

//----------------------------------------------------------------------
// buildSSP — compute the constant mean-dilatation B (sspBnot), element volume
// (sspVol) and the statically-condensed enhanced-strain stabilization
// (sspKstab). Ported verbatim from UWelements/SSPbrick::GetStab + the nodal
// pattern setup in SSPbrick::setDomain. The 9 enhanced-strain modes are
// condensed using the INITIAL material tangent, so the result is constant — it
// adapts to C(0) (curing both shear and volumetric locking across all nu) but
// carries no per-step internal state. Called once from setDomain.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick::buildSSP(void)
{
  // nodal coordinate matrix (3 x 8)
  Matrix mNodeCrd(3, 8);
  for (int n = 0; n < 8; n++) {
    const Vector &crd = nodePointers[n]->getCrds();
    mNodeCrd(0, n) = crd(0);
    mNodeCrd(1, n) = crd(1);
    mNodeCrd(2, n) = crd(2);
  }

  // nodal pattern vectors (Brick node order; SSPbrick scaling 0.125)
  Vector xi(8), et(8), ze(8), hut(8), hus(8), hst(8), hstu(8);
  xi(0)=-0.125; xi(1)= 0.125; xi(2)= 0.125; xi(3)=-0.125; xi(4)=-0.125; xi(5)= 0.125; xi(6)= 0.125; xi(7)=-0.125;
  et(0)=-0.125; et(1)=-0.125; et(2)= 0.125; et(3)= 0.125; et(4)=-0.125; et(5)=-0.125; et(6)= 0.125; et(7)= 0.125;
  ze(0)=-0.125; ze(1)=-0.125; ze(2)=-0.125; ze(3)=-0.125; ze(4)= 0.125; ze(5)= 0.125; ze(6)= 0.125; ze(7)= 0.125;
  hst(0)= 0.125; hst(1)=-0.125; hst(2)= 0.125; hst(3)=-0.125; hst(4)= 0.125; hst(5)=-0.125; hst(6)= 0.125; hst(7)=-0.125;
  hut(0)= 0.125; hut(1)= 0.125; hut(2)=-0.125; hut(3)=-0.125; hut(4)=-0.125; hut(5)=-0.125; hut(6)= 0.125; hut(7)= 0.125;
  hus(0)= 0.125; hus(1)=-0.125; hus(2)=-0.125; hus(3)= 0.125; hus(4)=-0.125; hus(5)= 0.125; hus(6)= 0.125; hus(7)=-0.125;
  hstu(0)=-0.125; hstu(1)= 0.125; hstu(2)=-0.125; hstu(3)= 0.125; hstu(4)= 0.125; hstu(5)=-0.125; hstu(6)= 0.125; hstu(7)=-0.125;

  // shape-function derivatives (local crd) at the element center
  Matrix dNloc(8, 3);
  dNloc(0,0)=-0.125; dNloc(1,0)= 0.125; dNloc(2,0)= 0.125; dNloc(3,0)=-0.125;
  dNloc(4,0)=-0.125; dNloc(5,0)= 0.125; dNloc(6,0)= 0.125; dNloc(7,0)=-0.125;
  dNloc(0,1)=-0.125; dNloc(1,1)=-0.125; dNloc(2,1)= 0.125; dNloc(3,1)= 0.125;
  dNloc(4,1)=-0.125; dNloc(5,1)=-0.125; dNloc(6,1)= 0.125; dNloc(7,1)= 0.125;
  dNloc(0,2)=-0.125; dNloc(1,2)=-0.125; dNloc(2,2)=-0.125; dNloc(3,2)=-0.125;
  dNloc(4,2)= 0.125; dNloc(5,2)= 0.125; dNloc(6,2)= 0.125; dNloc(7,2)= 0.125;

  // jacobian (center) and its inverse
  Matrix Jmat(3, 3), Jinv(3, 3);
  Jmat = mNodeCrd * dNloc;
  Jmat.Invert(Jinv);

  // nodal coordinate vectors
  Vector x(8), y(8), z(8);
  for (int n = 0; n < 8; n++) { x(n) = mNodeCrd(0, n); y(n) = mNodeCrd(1, n); z(n) = mNodeCrd(2, n); }

  // jacobian-determinant coefficient terms
  double a1=x^xi, a2=x^et, a3=x^ze, a4=x^hut, a5=x^hus, a6=x^hst, a7=x^hstu;
  double b1=y^xi, b2=y^et, b3=y^ze, b4=y^hut, b5=y^hus, b6=y^hst, b7=y^hstu;
  double c1=z^xi, c2=z^et, c3=z^ze, c4=z^hut, c5=z^hus, c6=z^hst, c7=z^hstu;

  Vector e1(3), e2(3), e3(3), e4(3), e5(3), e6(3), e7(3);
  e1(0)=a1; e1(1)=b1; e1(2)=c1;
  e2(0)=a2; e2(1)=b2; e2(2)=c2;
  e3(0)=a3; e3(1)=b3; e3(2)=c3;
  e4(0)=a4; e4(1)=b4; e4(2)=c4;
  e5(0)=a5; e5(1)=b5; e5(2)=c5;
  e6(0)=a6; e6(1)=b6; e6(2)=c6;
  e7(0)=a7; e7(1)=b7; e7(2)=c7;

  double J[20];
  J[0] = e1^(ladrunoCross(e2,e3));
  J[1] = 0.0;   // (SSPbrick zeroes J[1..3] after computing them)
  J[2] = 0.0;
  J[3] = 0.0;
  J[4] = (e7^(ladrunoCross(e2,e3))) + (e4^(ladrunoCross(e5,e2))) + (e4^(ladrunoCross(e3,e6)));
  J[5] = (e1^(ladrunoCross(e7,e3))) + (e4^(ladrunoCross(e5,e1))) + (e3^(ladrunoCross(e5,e6)));
  J[6] = (e1^(ladrunoCross(e2,e7))) + (e4^(ladrunoCross(e1,e6))) + (e2^(ladrunoCross(e5,e6)));
  J[7] = -1.0*e1^(ladrunoCross(e5,e6));
  J[8] = -1.0*e4^(ladrunoCross(e2,e6));
  J[9] = -1.0*e4^(ladrunoCross(e5,e3));
  J[10] = e2^(ladrunoCross(e4,e7));
  J[11] = -1.0*e3^(ladrunoCross(e4,e7));
  J[12] = e3^(ladrunoCross(e5,e7));
  J[13] = -1.0*e1^(ladrunoCross(e5,e7));
  J[14] = e1^(ladrunoCross(e6,e7));
  J[15] = -1.0*e2^(ladrunoCross(e6,e7));
  J[16] = 2.0*e4^(ladrunoCross(e5,e6));
  J[17] = e7^(ladrunoCross(e5,e6));
  J[18] = e4^(ladrunoCross(e7,e6));
  J[19] = e4^(ladrunoCross(e5,e7));

  double J0789  = 8.0*(J[0]/3.0 + J[7]/5.0 + J[8]/9.0 + J[9]/9.0);
  double J0879  = 8.0*(J[0]/3.0 + J[8]/5.0 + J[7]/9.0 + J[9]/9.0);
  double J0978  = 8.0*(J[0]/3.0 + J[9]/5.0 + J[7]/9.0 + J[8]/9.0);
  double J417   = 8.0*(J[4]/9.0 + J[17]/27.0);
  double J518   = 8.0*(J[5]/9.0 + J[18]/27.0);
  double J619   = 8.0*(J[6]/9.0 + J[19]/27.0);
  double J11215 = 8.0*(J[1]/9.0 + J[12]/15.0 + J[15]/27.0);
  double J11512 = 8.0*(J[1]/9.0 + J[15]/15.0 + J[12]/27.0);
  double J21114 = 8.0*(J[2]/9.0 + J[11]/15.0 + J[14]/27.0);
  double J21411 = 8.0*(J[2]/9.0 + J[14]/15.0 + J[11]/27.0);
  double J31013 = 8.0*(J[3]/9.0 + J[10]/15.0 + J[13]/27.0);
  double J31310 = 8.0*(J[3]/9.0 + J[13]/15.0 + J[10]/27.0);
  double J789   = 8.0*(J[0]/9.0 + J[7]/15.0 + J[8]/15.0 + J[9]/27.0);
  double J897   = 8.0*(J[0]/9.0 + J[8]/15.0 + J[9]/15.0 + J[7]/27.0);
  double J798   = 8.0*(J[0]/9.0 + J[7]/15.0 + J[9]/15.0 + J[8]/27.0);
  double J174   = 8.0*(J[4]/27.0 + 64.0*J[17]/45.0);
  double J185   = 8.0*(J[5]/27.0 + 64.0*J[18]/45.0);
  double J196   = 8.0*(J[6]/27.0 + 64.0*J[19]/45.0);
  double J16    = 8.0*J[16]/27.0;

  // element volume
  double mVol = 8.0*(J[0] + (J[7] + J[8] + J[9])/3.0);

  // kinematic (mean-dilatation, B-bar) shape gradients
  Vector bx(8), by(8), bz(8);
  bx = (8.0*((b2*c3-c2*b3)*xi + (b3*c1-c3*b1)*et + (b1*c2-c1*b2)*ze) + (8.0/3.0)*((b6*c5-c6*b5)*xi + (b4*c6-c4*b6)*et + (b5*c4-c5*b4)*ze
       + (b5*c1-c5*b1 + b2*c4-c2*b4)*hst + (b6*c2-c6*b2 + b3*c5-c3*b5)*hut + (b1*c6-c1*b6 + b4*c3-c4*b3)*hus))/mVol;
  by = (8.0*((c2*a3-a2*c3)*xi + (c3*a1-a3*c1)*et + (c1*a2-a1*c2)*ze) + (8.0/3.0)*((c6*a5-a6*c5)*xi + (c4*a6-a4*c6)*et + (c5*a4-a5*c4)*ze
       + (c5*a1-a5*c1 + c2*a4-a2*c4)*hst + (c6*a2-a6*c2 + c3*a5-a3*c5)*hut + (c1*a6-a1*c6 + c4*a3-a4*c3)*hus))/mVol;
  bz = (8.0*((a2*b3-b2*a3)*xi + (a3*b1-b3*a1)*et + (a1*b2-b1*a2)*ze) + (8.0/3.0)*((a6*b5-b6*a5)*xi + (a4*b6-b4*a6)*et + (a5*b4-b5*a4)*ze
       + (a5*b1-b5*a1 + a2*b4-b2*a4)*hst + (a6*b2-b6*a2 + a3*b5-b3*a5)*hut + (a1*b6-b1*a6 + a4*b3-b4*a3)*hus))/mVol;

  Matrix dNmod(8, 3);
  for (int i = 0; i < 8; i++) { dNmod(i,0) = bx(i); dNmod(i,1) = by(i); dNmod(i,2) = bz(i); }

  // hourglass transformation matrix G = I - dNmod * X, and gamma vectors
  Matrix I8(8, 8);
  for (int i = 0; i < 8; i++) I8(i, i) = 1.0;
  Matrix G(8, 8);
  G = I8 - dNmod * mNodeCrd;
  Vector gst  = G * hst;
  Vector gut  = G * hut;
  Vector gus  = G * hus;
  Vector gstu = G * hstu;

  // mean-dilatation B (membrane modes) and Mben (hourglass mapping)
  Matrix Bnot(6, 24), Mben(12, 24);
  Bnot.Zero();
  Mben.Zero();
  for (int i = 0; i < 8; i++) {
    Bnot(0,3*i)   = dNmod(i,0);
    Bnot(1,3*i+1) = dNmod(i,1);
    Bnot(2,3*i+2) = dNmod(i,2);
    Bnot(3,3*i)   = dNmod(i,1);
    Bnot(3,3*i+1) = dNmod(i,0);
    Bnot(4,3*i+1) = dNmod(i,2);
    Bnot(4,3*i+2) = dNmod(i,1);
    Bnot(5,3*i)   = dNmod(i,2);
    Bnot(5,3*i+2) = dNmod(i,0);

    Mben(0,3*i)    = gst(i);
    Mben(1,3*i+1)  = gst(i);
    Mben(2,3*i+2)  = gst(i);
    Mben(3,3*i)    = gut(i);
    Mben(4,3*i+1)  = gut(i);
    Mben(5,3*i+2)  = gut(i);
    Mben(6,3*i)    = gus(i);
    Mben(7,3*i+1)  = gus(i);
    Mben(8,3*i+2)  = gus(i);
    Mben(9,3*i)    = gstu(i);
    Mben(10,3*i+1) = gstu(i);
    Mben(11,3*i+2) = gstu(i);
  }

  // FCF terms (12x12)
  double HstXX = J0879*Jinv(0,0)*Jinv(0,0) + J0789*Jinv(1,0)*Jinv(1,0) + J619*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0));
  double HstXY = J0879*Jinv(0,0)*Jinv(0,1) + J0789*Jinv(1,0)*Jinv(1,1) + J619*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1));
  double HstXZ = J0879*Jinv(0,0)*Jinv(0,2) + J0789*Jinv(1,0)*Jinv(1,2) + J619*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2));
  double HstYY = J0879*Jinv(0,1)*Jinv(0,1) + J0789*Jinv(1,1)*Jinv(1,1) + J619*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1));
  double HstYZ = J0879*Jinv(0,1)*Jinv(0,2) + J0789*Jinv(1,1)*Jinv(1,2) + J619*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2));
  double HstZZ = J0879*Jinv(0,2)*Jinv(0,2) + J0789*Jinv(1,2)*Jinv(1,2) + J619*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2));

  double HutXX = J0978*Jinv(1,0)*Jinv(1,0) + J0879*Jinv(2,0)*Jinv(2,0) + J417*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
  double HutXY = J0978*Jinv(1,0)*Jinv(1,1) + J0879*Jinv(2,0)*Jinv(2,1) + J417*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
  double HutXZ = J0978*Jinv(1,0)*Jinv(1,2) + J0879*Jinv(2,0)*Jinv(2,2) + J417*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
  double HutYY = J0978*Jinv(1,1)*Jinv(1,1) + J0879*Jinv(2,1)*Jinv(2,1) + J417*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
  double HutYZ = J0978*Jinv(1,1)*Jinv(1,2) + J0879*Jinv(2,1)*Jinv(2,2) + J417*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
  double HutZZ = J0978*Jinv(1,2)*Jinv(1,2) + J0879*Jinv(2,2)*Jinv(2,2) + J417*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

  double HusXX = J0978*Jinv(0,0)*Jinv(0,0) + J0789*Jinv(2,0)*Jinv(2,0) + J518*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0));
  double HusXY = J0978*Jinv(0,0)*Jinv(0,1) + J0789*Jinv(2,0)*Jinv(2,1) + J518*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1));
  double HusXZ = J0978*Jinv(0,0)*Jinv(0,2) + J0789*Jinv(2,0)*Jinv(2,2) + J518*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2));
  double HusYY = J0978*Jinv(0,1)*Jinv(0,1) + J0789*Jinv(2,1)*Jinv(2,1) + J518*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1));
  double HusYZ = J0978*Jinv(0,1)*Jinv(0,2) + J0789*Jinv(2,1)*Jinv(2,2) + J518*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2));
  double HusZZ = J0978*Jinv(0,2)*Jinv(0,2) + J0789*Jinv(2,2)*Jinv(2,2) + J518*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2));

  double HstuXX = J897*Jinv(0,0)*Jinv(0,0) + J798*Jinv(1,0)*Jinv(1,0) + J789*Jinv(2,0)*Jinv(2,0) + J185*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0))
                  + J196*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0)) + J174*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
  double HstuXY = J897*Jinv(0,0)*Jinv(0,1) + J798*Jinv(1,0)*Jinv(1,1) + J789*Jinv(2,0)*Jinv(2,1) + J185*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1))
                  + J196*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1)) + J174*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
  double HstuXZ = J897*Jinv(0,0)*Jinv(0,2) + J798*Jinv(1,0)*Jinv(1,2) + J789*Jinv(2,0)*Jinv(2,2) + J185*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2))
                  + J196*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2)) + J174*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
  double HstuYY = J897*Jinv(0,1)*Jinv(0,1) + J798*Jinv(1,1)*Jinv(1,1) + J789*Jinv(2,1)*Jinv(2,1) + J185*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1))
                  + J196*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1)) + J174*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
  double HstuYZ = J897*Jinv(0,1)*Jinv(0,2) + J798*Jinv(1,1)*Jinv(1,2) + J789*Jinv(2,1)*Jinv(2,2) + J185*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2))
                  + J196*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2)) + J174*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
  double HstuZZ = J897*Jinv(0,2)*Jinv(0,2) + J798*Jinv(1,2)*Jinv(1,2) + J789*Jinv(2,2)*Jinv(2,2) + J185*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2))
                  + J196*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2)) + J174*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

  double IttXX = J0879*Jinv(0,0)*Jinv(2,0) + J417*Jinv(0,0)*Jinv(1,0) + J518*Jinv(1,0)*Jinv(1,0) + J619*Jinv(1,0)*Jinv(2,0);
  double IttXY = J0879*Jinv(0,0)*Jinv(2,1) + J417*Jinv(0,0)*Jinv(1,1) + J518*Jinv(1,0)*Jinv(1,1) + J619*Jinv(1,0)*Jinv(2,1);
  double IttXZ = J0879*Jinv(0,0)*Jinv(2,2) + J417*Jinv(0,0)*Jinv(1,2) + J518*Jinv(1,0)*Jinv(1,2) + J619*Jinv(1,0)*Jinv(2,2);
  double IttYX = J0879*Jinv(0,1)*Jinv(2,0) + J417*Jinv(0,1)*Jinv(1,0) + J518*Jinv(1,1)*Jinv(1,0) + J619*Jinv(1,1)*Jinv(2,0);
  double IttYY = J0879*Jinv(0,1)*Jinv(2,1) + J417*Jinv(0,1)*Jinv(1,1) + J518*Jinv(1,1)*Jinv(1,1) + J619*Jinv(1,1)*Jinv(2,1);
  double IttYZ = J0879*Jinv(0,1)*Jinv(2,2) + J417*Jinv(0,1)*Jinv(1,2) + J518*Jinv(1,1)*Jinv(1,2) + J619*Jinv(1,1)*Jinv(2,2);
  double IttZX = J0879*Jinv(0,2)*Jinv(2,0) + J417*Jinv(0,2)*Jinv(1,0) + J518*Jinv(1,2)*Jinv(1,0) + J619*Jinv(1,2)*Jinv(2,0);
  double IttZY = J0879*Jinv(0,2)*Jinv(2,1) + J417*Jinv(0,2)*Jinv(1,1) + J518*Jinv(1,2)*Jinv(1,1) + J619*Jinv(1,2)*Jinv(2,1);
  double IttZZ = J0879*Jinv(0,2)*Jinv(2,2) + J417*Jinv(0,2)*Jinv(1,2) + J518*Jinv(1,2)*Jinv(1,2) + J619*Jinv(1,2)*Jinv(2,2);

  double IssXX = J0789*Jinv(1,0)*Jinv(2,0) + J417*Jinv(0,0)*Jinv(0,0) + J518*Jinv(1,0)*Jinv(0,0) + J619*Jinv(0,0)*Jinv(2,0);
  double IssXY = J0789*Jinv(1,0)*Jinv(2,1) + J417*Jinv(0,0)*Jinv(0,1) + J518*Jinv(1,0)*Jinv(0,1) + J619*Jinv(0,0)*Jinv(2,1);
  double IssXZ = J0789*Jinv(1,0)*Jinv(2,2) + J417*Jinv(0,0)*Jinv(0,2) + J518*Jinv(1,0)*Jinv(0,2) + J619*Jinv(0,0)*Jinv(2,2);
  double IssYX = J0789*Jinv(1,1)*Jinv(2,0) + J417*Jinv(0,1)*Jinv(0,0) + J518*Jinv(1,1)*Jinv(0,0) + J619*Jinv(0,1)*Jinv(2,0);
  double IssYY = J0789*Jinv(1,1)*Jinv(2,1) + J417*Jinv(0,1)*Jinv(0,1) + J518*Jinv(1,1)*Jinv(0,1) + J619*Jinv(0,1)*Jinv(2,1);
  double IssYZ = J0789*Jinv(1,1)*Jinv(2,2) + J417*Jinv(0,1)*Jinv(0,2) + J518*Jinv(1,1)*Jinv(0,2) + J619*Jinv(0,1)*Jinv(2,2);
  double IssZX = J0789*Jinv(1,2)*Jinv(2,0) + J417*Jinv(0,2)*Jinv(0,0) + J518*Jinv(1,2)*Jinv(0,0) + J619*Jinv(0,2)*Jinv(2,0);
  double IssZY = J0789*Jinv(1,2)*Jinv(2,1) + J417*Jinv(0,2)*Jinv(0,1) + J518*Jinv(1,2)*Jinv(0,1) + J619*Jinv(0,2)*Jinv(2,1);
  double IssZZ = J0789*Jinv(1,2)*Jinv(2,2) + J417*Jinv(0,2)*Jinv(0,2) + J518*Jinv(1,2)*Jinv(0,2) + J619*Jinv(0,2)*Jinv(2,2);

  double IuuXX = J0978*Jinv(1,0)*Jinv(0,0) + J417*Jinv(2,0)*Jinv(0,0) + J518*Jinv(1,0)*Jinv(2,0) + J619*Jinv(2,0)*Jinv(2,0);
  double IuuXY = J0978*Jinv(1,0)*Jinv(0,1) + J417*Jinv(2,0)*Jinv(0,1) + J518*Jinv(1,0)*Jinv(2,1) + J619*Jinv(2,0)*Jinv(2,1);
  double IuuXZ = J0978*Jinv(1,0)*Jinv(0,2) + J417*Jinv(2,0)*Jinv(0,2) + J518*Jinv(1,0)*Jinv(2,2) + J619*Jinv(2,0)*Jinv(2,2);
  double IuuYX = J0978*Jinv(1,1)*Jinv(0,0) + J417*Jinv(2,1)*Jinv(0,0) + J518*Jinv(1,1)*Jinv(2,0) + J619*Jinv(2,1)*Jinv(2,0);
  double IuuYY = J0978*Jinv(1,1)*Jinv(0,1) + J417*Jinv(2,1)*Jinv(0,1) + J518*Jinv(1,1)*Jinv(2,1) + J619*Jinv(2,1)*Jinv(2,1);
  double IuuYZ = J0978*Jinv(1,1)*Jinv(0,2) + J417*Jinv(2,1)*Jinv(0,2) + J518*Jinv(1,1)*Jinv(2,2) + J619*Jinv(2,1)*Jinv(2,2);
  double IuuZX = J0978*Jinv(1,2)*Jinv(0,0) + J417*Jinv(2,2)*Jinv(0,0) + J518*Jinv(1,2)*Jinv(2,0) + J619*Jinv(2,2)*Jinv(2,0);
  double IuuZY = J0978*Jinv(1,2)*Jinv(0,1) + J417*Jinv(2,2)*Jinv(0,1) + J518*Jinv(1,2)*Jinv(2,1) + J619*Jinv(2,2)*Jinv(2,1);
  double IuuZZ = J0978*Jinv(1,2)*Jinv(0,2) + J417*Jinv(2,2)*Jinv(0,2) + J518*Jinv(1,2)*Jinv(2,2) + J619*Jinv(2,2)*Jinv(2,2);

  double IstXX = J31013*Jinv(0,0)*Jinv(0,0) + J31310*Jinv(1,0)*Jinv(1,0) + J21411*Jinv(2,0)*Jinv(1,0) + J11512*Jinv(2,0)*Jinv(0,0) + J16*(Jinv(0,0)*Jinv(1,0) + Jinv(1,0)*Jinv(0,0));
  double IstXY = J31013*Jinv(0,0)*Jinv(0,1) + J31310*Jinv(1,0)*Jinv(1,1) + J21411*Jinv(2,0)*Jinv(1,1) + J11512*Jinv(2,0)*Jinv(0,1) + J16*(Jinv(0,0)*Jinv(1,1) + Jinv(1,0)*Jinv(0,1));
  double IstXZ = J31013*Jinv(0,0)*Jinv(0,2) + J31310*Jinv(1,0)*Jinv(1,2) + J21411*Jinv(2,0)*Jinv(1,2) + J11512*Jinv(2,0)*Jinv(0,2) + J16*(Jinv(0,0)*Jinv(1,2) + Jinv(1,0)*Jinv(0,2));
  double IstYX = J31013*Jinv(0,1)*Jinv(0,0) + J31310*Jinv(1,1)*Jinv(1,0) + J21411*Jinv(2,1)*Jinv(1,0) + J11512*Jinv(2,1)*Jinv(0,0) + J16*(Jinv(0,1)*Jinv(1,0) + Jinv(1,1)*Jinv(0,0));
  double IstYY = J31013*Jinv(0,1)*Jinv(0,1) + J31310*Jinv(1,1)*Jinv(1,1) + J21411*Jinv(2,1)*Jinv(1,1) + J11512*Jinv(2,1)*Jinv(0,1) + J16*(Jinv(0,1)*Jinv(1,1) + Jinv(1,1)*Jinv(0,1));
  double IstYZ = J31013*Jinv(0,1)*Jinv(0,2) + J31310*Jinv(1,1)*Jinv(1,2) + J21411*Jinv(2,1)*Jinv(1,2) + J11512*Jinv(2,1)*Jinv(0,2) + J16*(Jinv(0,1)*Jinv(1,2) + Jinv(1,1)*Jinv(0,2));
  double IstZX = J31013*Jinv(0,2)*Jinv(0,0) + J31310*Jinv(1,2)*Jinv(1,0) + J21411*Jinv(2,2)*Jinv(1,0) + J11512*Jinv(2,2)*Jinv(0,0) + J16*(Jinv(0,2)*Jinv(1,0) + Jinv(1,2)*Jinv(0,0));
  double IstZY = J31013*Jinv(0,2)*Jinv(0,1) + J31310*Jinv(1,2)*Jinv(1,1) + J21411*Jinv(2,2)*Jinv(1,1) + J11512*Jinv(2,2)*Jinv(0,1) + J16*(Jinv(0,2)*Jinv(1,1) + Jinv(1,2)*Jinv(0,1));
  double IstZZ = J31013*Jinv(0,2)*Jinv(0,2) + J31310*Jinv(1,2)*Jinv(1,2) + J21411*Jinv(2,2)*Jinv(1,2) + J11512*Jinv(2,2)*Jinv(0,2) + J16*(Jinv(0,2)*Jinv(1,2) + Jinv(1,2)*Jinv(0,2));

  double IutXX = J21114*Jinv(0,0)*Jinv(1,0) + J31013*Jinv(0,0)*Jinv(2,0) + J11215*Jinv(1,0)*Jinv(1,0) + J11512*Jinv(2,0)*Jinv(2,0) + J16*(Jinv(1,0)*Jinv(2,0) + Jinv(2,0)*Jinv(1,0));
  double IutXY = J21114*Jinv(0,0)*Jinv(1,1) + J31013*Jinv(0,0)*Jinv(2,1) + J11215*Jinv(1,0)*Jinv(1,1) + J11512*Jinv(2,0)*Jinv(2,1) + J16*(Jinv(1,0)*Jinv(2,1) + Jinv(2,0)*Jinv(1,1));
  double IutXZ = J21114*Jinv(0,0)*Jinv(1,2) + J31013*Jinv(0,0)*Jinv(2,2) + J11215*Jinv(1,0)*Jinv(1,2) + J11512*Jinv(2,0)*Jinv(2,2) + J16*(Jinv(1,0)*Jinv(2,2) + Jinv(2,0)*Jinv(1,2));
  double IutYX = J21114*Jinv(0,1)*Jinv(1,0) + J31013*Jinv(0,1)*Jinv(2,0) + J11215*Jinv(1,1)*Jinv(1,0) + J11512*Jinv(2,1)*Jinv(2,0) + J16*(Jinv(1,1)*Jinv(2,0) + Jinv(2,1)*Jinv(1,0));
  double IutYY = J21114*Jinv(0,1)*Jinv(1,1) + J31013*Jinv(0,1)*Jinv(2,1) + J11215*Jinv(1,1)*Jinv(1,1) + J11512*Jinv(2,1)*Jinv(2,1) + J16*(Jinv(1,1)*Jinv(2,1) + Jinv(2,1)*Jinv(1,1));
  double IutYZ = J21114*Jinv(0,1)*Jinv(1,2) + J31013*Jinv(0,1)*Jinv(2,2) + J11215*Jinv(1,1)*Jinv(1,2) + J11512*Jinv(2,1)*Jinv(2,2) + J16*(Jinv(1,1)*Jinv(2,2) + Jinv(2,1)*Jinv(1,2));
  double IutZX = J21114*Jinv(0,2)*Jinv(1,0) + J31013*Jinv(0,2)*Jinv(2,0) + J11215*Jinv(1,2)*Jinv(1,0) + J11512*Jinv(2,2)*Jinv(2,0) + J16*(Jinv(1,2)*Jinv(2,0) + Jinv(2,2)*Jinv(1,0));
  double IutZY = J21114*Jinv(0,2)*Jinv(1,1) + J31013*Jinv(0,2)*Jinv(2,1) + J11215*Jinv(1,2)*Jinv(1,1) + J11512*Jinv(2,2)*Jinv(2,1) + J16*(Jinv(1,2)*Jinv(2,1) + Jinv(2,2)*Jinv(1,1));
  double IutZZ = J21114*Jinv(0,2)*Jinv(1,2) + J31013*Jinv(0,2)*Jinv(2,2) + J11215*Jinv(1,2)*Jinv(1,2) + J11512*Jinv(2,2)*Jinv(2,2) + J16*(Jinv(1,2)*Jinv(2,2) + Jinv(2,2)*Jinv(1,2));

  double IusXX = J21114*Jinv(0,0)*Jinv(0,0) + J11215*Jinv(1,0)*Jinv(0,0) + J31310*Jinv(1,0)*Jinv(2,0) + J21411*Jinv(2,0)*Jinv(2,0) + J16*(Jinv(0,0)*Jinv(2,0) + Jinv(2,0)*Jinv(0,0));
  double IusXY = J21114*Jinv(0,0)*Jinv(0,1) + J11215*Jinv(1,0)*Jinv(0,1) + J31310*Jinv(1,0)*Jinv(2,1) + J21411*Jinv(2,0)*Jinv(2,1) + J16*(Jinv(0,0)*Jinv(2,1) + Jinv(2,0)*Jinv(0,1));
  double IusXZ = J21114*Jinv(0,0)*Jinv(0,2) + J11215*Jinv(1,0)*Jinv(0,2) + J31310*Jinv(1,0)*Jinv(2,2) + J21411*Jinv(2,0)*Jinv(2,2) + J16*(Jinv(0,0)*Jinv(2,2) + Jinv(2,0)*Jinv(0,2));
  double IusYX = J21114*Jinv(0,1)*Jinv(0,0) + J11215*Jinv(1,1)*Jinv(0,0) + J31310*Jinv(1,1)*Jinv(2,0) + J21411*Jinv(2,1)*Jinv(2,0) + J16*(Jinv(0,1)*Jinv(2,0) + Jinv(2,1)*Jinv(0,0));
  double IusYY = J21114*Jinv(0,1)*Jinv(0,1) + J11215*Jinv(1,1)*Jinv(0,1) + J31310*Jinv(1,1)*Jinv(2,1) + J21411*Jinv(2,1)*Jinv(2,1) + J16*(Jinv(0,1)*Jinv(2,1) + Jinv(2,1)*Jinv(0,1));
  double IusYZ = J21114*Jinv(0,1)*Jinv(0,2) + J11215*Jinv(1,1)*Jinv(0,2) + J31310*Jinv(1,1)*Jinv(2,2) + J21411*Jinv(2,1)*Jinv(2,2) + J16*(Jinv(0,1)*Jinv(2,2) + Jinv(2,1)*Jinv(0,2));
  double IusZX = J21114*Jinv(0,2)*Jinv(0,0) + J11215*Jinv(1,2)*Jinv(0,0) + J31310*Jinv(1,2)*Jinv(2,0) + J21411*Jinv(2,2)*Jinv(2,0) + J16*(Jinv(0,2)*Jinv(2,0) + Jinv(2,2)*Jinv(0,0));
  double IusZY = J21114*Jinv(0,2)*Jinv(0,1) + J11215*Jinv(1,2)*Jinv(0,1) + J31310*Jinv(1,2)*Jinv(2,1) + J21411*Jinv(2,2)*Jinv(2,1) + J16*(Jinv(0,2)*Jinv(2,1) + Jinv(2,2)*Jinv(0,1));
  double IusZZ = J21114*Jinv(0,2)*Jinv(0,2) + J11215*Jinv(1,2)*Jinv(0,2) + J31310*Jinv(1,2)*Jinv(2,2) + J21411*Jinv(2,2)*Jinv(2,2) + J16*(Jinv(0,2)*Jinv(2,2) + Jinv(2,2)*Jinv(0,2));

  // constitutive constants from the INITIAL tangent (=> constant stabilization)
  const Matrix &CmatI = materialPointers[0]->getInitialTangent();
  double C1 = CmatI(0,0);
  double C2 = CmatI(0,1);
  double C3 = CmatI(3,3);
  double C4 = C2 + C3;

  Matrix FCF(12, 12);
  // block11
  FCF(0,0) = C1*HstXX + C3*(HstYY + HstZZ);
  FCF(0,1) = C4*HstXY;  FCF(0,2) = C4*HstXZ;
  FCF(1,0) = FCF(0,1);  FCF(1,1) = C1*HstYY + C3*(HstXX + HstZZ);  FCF(1,2) = C4*HstYZ;
  FCF(2,0) = FCF(0,2);  FCF(2,1) = FCF(1,2);  FCF(2,2) = C1*HstZZ + C3*(HstYY + HstXX);
  // block12
  FCF(0,3) = C1*IttXX + C3*(IttYY + IttZZ);  FCF(0,4) = C2*IttXY + C3*IttYX;  FCF(0,5) = C2*IttXZ + C3*IttZX;
  FCF(1,3) = C2*IttYX + C3*IttXY;  FCF(1,4) = C1*IttYY + C3*(IttXX + IttZZ);  FCF(1,5) = C2*IttYZ + C3*IttZY;
  FCF(2,3) = C2*IttZX + C3*IttXZ;  FCF(2,4) = C2*IttZY + C3*IttYZ;  FCF(2,5) = C1*IttZZ + C3*(IttYY + IttXX);
  // block13
  FCF(0,6) = C1*IssXX + C3*(IssYY + IssZZ);  FCF(0,7) = C2*IssXY + C3*IssYX;  FCF(0,8) = C2*IssXZ + C3*IssZX;
  FCF(1,6) = C2*IssYX + C3*IssXY;  FCF(1,7) = C1*IssYY + C3*(IssXX + IssZZ);  FCF(1,8) = C2*IssYZ + C3*IssZY;
  FCF(2,6) = C2*IssZX + C3*IssXZ;  FCF(2,7) = C2*IssZY + C3*IssYZ;  FCF(2,8) = C1*IssZZ + C3*(IssYY + IssXX);
  // block14
  FCF(0,9)  = C3*(IstYY + IstZZ);  FCF(0,10) = C3*IstXY;  FCF(0,11) = C3*IstXZ;
  FCF(1,9)  = C3*IstYX;  FCF(1,10) = C3*(IstXX + IstZZ);  FCF(1,11) = C3*IstYZ;
  FCF(2,9)  = C3*IstZX;  FCF(2,10) = C3*IstZY;  FCF(2,11) = C3*(IstYY + IstXX);
  // block21
  FCF(3,0) = C1*IttXX + C3*(IttYY + IttZZ);  FCF(3,1) = C2*IttYX + C3*IttXY;  FCF(3,2) = C2*IttZX + C3*IttXZ;
  FCF(4,0) = C2*IttXY + C3*IttYX;  FCF(4,1) = C1*IttYY + C3*(IttXX + IttZZ);  FCF(4,2) = C2*IttZY + C3*IttYZ;
  FCF(5,0) = C2*IttXZ + C3*IttZX;  FCF(5,1) = C2*IttYZ + C3*IttZY;  FCF(5,2) = C1*IttZZ + C3*(IttYY + IttXX);
  // block22
  FCF(3,3) = C1*HutXX + C3*(HutYY + HutZZ);  FCF(3,4) = C4*HutXY;  FCF(3,5) = C4*HutXZ;
  FCF(4,3) = FCF(3,4);  FCF(4,4) = C1*HutYY + C3*(HutXX + HutZZ);  FCF(4,5) = C4*HutYZ;
  FCF(5,3) = FCF(3,5);  FCF(5,4) = FCF(4,5);  FCF(5,5) = C1*HutZZ + C3*(HutYY + HutXX);
  // block23
  FCF(3,6) = C1*IuuXX + C3*(IuuYY + IuuZZ);  FCF(3,7) = C2*IuuXY + C3*IuuYX;  FCF(3,8) = C2*IuuXZ + C3*IuuZX;
  FCF(4,6) = C2*IuuYX + C3*IuuXY;  FCF(4,7) = C1*IuuYY + C3*(IuuXX + IuuZZ);  FCF(4,8) = C2*IuuYZ + C3*IuuZY;
  FCF(5,6) = C2*IuuZX + C3*IuuXZ;  FCF(5,7) = C2*IuuZY + C3*IuuYZ;  FCF(5,8) = C1*IuuZZ + C3*(IuuYY + IuuXX);
  // block24
  FCF(3,9)  = C3*(IutYY + IutZZ);  FCF(3,10) = C3*IutXY;  FCF(3,11) = C3*IutXZ;
  FCF(4,9)  = C3*IutYX;  FCF(4,10) = C3*(IutXX + IutZZ);  FCF(4,11) = C3*IutYZ;
  FCF(5,9)  = C3*IutZX;  FCF(5,10) = C3*IutZY;  FCF(5,11) = C3*(IutYY + IutXX);
  // block31
  FCF(6,0) = C1*IssXX + C3*(IssYY + IssZZ);  FCF(6,1) = C2*IssYX + C3*IssXY;  FCF(6,2) = C2*IssZX + C3*IssXZ;
  FCF(7,0) = C2*IssXY + C3*IssYX;  FCF(7,1) = C1*IssYY + C3*(IssXX + IssZZ);  FCF(7,2) = C2*IssZY + C3*IssYZ;
  FCF(8,0) = C2*IssXZ + C3*IssZX;  FCF(8,1) = C2*IssYZ + C3*IssZY;  FCF(8,2) = C1*IssZZ + C3*(IssYY + IssXX);
  // block32
  FCF(6,3) = C1*IuuXX + C3*(IuuYY + IuuZZ);  FCF(6,4) = C2*IuuYX + C3*IuuXY;  FCF(6,5) = C2*IuuZX + C3*IuuXZ;
  FCF(7,3) = C2*IuuXY + C3*IuuYX;  FCF(7,4) = C1*IuuYY + C3*(IuuXX + IuuZZ);  FCF(7,5) = C2*IuuZY + C3*IuuYZ;
  FCF(8,3) = C2*IuuXZ + C3*IuuZX;  FCF(8,4) = C2*IuuYZ + C3*IuuZY;  FCF(8,5) = C1*IuuZZ + C3*(IuuYY + IuuXX);
  // block33
  FCF(6,6) = C1*HusXX + C3*(HusYY + HusZZ);  FCF(6,7) = C4*HusXY;  FCF(6,8) = C4*HusXZ;
  FCF(7,6) = FCF(6,7);  FCF(7,7) = C1*HusYY + C3*(HusXX + HusZZ);  FCF(7,8) = C4*HusYZ;
  FCF(8,6) = FCF(6,8);  FCF(8,7) = FCF(7,8);  FCF(8,8) = C1*HusZZ + C3*(HusYY + HusXX);
  // block34
  FCF(6,9)  = C3*(IusYY + IusZZ);  FCF(6,10) = C3*IusXY;  FCF(6,11) = C3*IusXZ;
  FCF(7,9)  = C3*IusYX;  FCF(7,10) = C3*(IusXX + IusZZ);  FCF(7,11) = C3*IusYZ;
  FCF(8,9)  = C3*IusZX;  FCF(8,10) = C3*IusZY;  FCF(8,11) = C3*(IusYY + IusXX);
  // block41
  FCF(9,0)  = C3*(IstYY + IstZZ);  FCF(9,1)  = C3*IstYX;  FCF(9,2)  = C3*IstZX;
  FCF(10,0) = C3*IstXY;  FCF(10,1) = C3*(IstXX + IstZZ);  FCF(10,2) = C3*IstZY;
  FCF(11,0) = C3*IstXZ;  FCF(11,1) = C3*IstYZ;  FCF(11,2) = C3*(IstYY + IstXX);
  // block42
  FCF(9,3)  = C3*(IutYY + IutZZ);  FCF(9,4)  = C3*IutYX;  FCF(9,5)  = C3*IutZX;
  FCF(10,3) = C3*IutXY;  FCF(10,4) = C3*(IutXX + IutZZ);  FCF(10,5) = C3*IutZY;
  FCF(11,3) = C3*IutXZ;  FCF(11,4) = C3*IutYZ;  FCF(11,5) = C3*(IutYY + IutXX);
  // block43
  FCF(9,6)  = C3*(IusYY + IusZZ);  FCF(9,7)  = C3*IusYX;  FCF(9,8)  = C3*IusZX;
  FCF(10,6) = C3*IusXY;  FCF(10,7) = C3*(IusXX + IusZZ);  FCF(10,8) = C3*IusZY;
  FCF(11,6) = C3*IusXZ;  FCF(11,7) = C3*IusYZ;  FCF(11,8) = C3*(IusYY + IusXX);
  // block44
  FCF(9,9)   = C3*(HstuYY + HstuZZ);  FCF(9,10)  = C3*HstuXY;  FCF(9,11)  = C3*HstuXZ;
  FCF(10,9)  = FCF(9,10);  FCF(10,10) = C3*(HstuXX + HstuZZ);  FCF(10,11) = C3*HstuYZ;
  FCF(11,9)  = FCF(9,11);  FCF(11,10) = FCF(10,11);  FCF(11,11) = C3*(HstuYY + HstuXX);

  // enhanced-strain constitutive coefficients
  double CssXX = C1*Jinv(0,0)*Jinv(0,0) + C3*(Jinv(0,1)*Jinv(0,1) + Jinv(0,2)*Jinv(0,2));
  double CssXY = C4*Jinv(0,0)*Jinv(0,1);
  double CssXZ = C4*Jinv(0,0)*Jinv(0,2);
  double CssYY = C1*Jinv(0,1)*Jinv(0,1) + C3*(Jinv(0,0)*Jinv(0,0) + Jinv(0,2)*Jinv(0,2));
  double CssYZ = C4*Jinv(0,1)*Jinv(0,2);
  double CssZZ = C1*Jinv(0,2)*Jinv(0,2) + C3*(Jinv(0,0)*Jinv(0,0) + Jinv(0,1)*Jinv(0,1));

  double CttXX = C1*Jinv(1,0)*Jinv(1,0) + C3*(Jinv(1,1)*Jinv(1,1) + Jinv(1,2)*Jinv(1,2));
  double CttXY = C4*Jinv(1,0)*Jinv(1,1);
  double CttXZ = C4*Jinv(1,0)*Jinv(1,2);
  double CttYY = C1*Jinv(1,1)*Jinv(1,1) + C3*(Jinv(1,0)*Jinv(1,0) + Jinv(1,2)*Jinv(1,2));
  double CttYZ = C4*Jinv(1,1)*Jinv(1,2);
  double CttZZ = C1*Jinv(1,2)*Jinv(1,2) + C3*(Jinv(1,0)*Jinv(1,0) + Jinv(1,1)*Jinv(1,1));

  double CuuXX = C1*Jinv(2,0)*Jinv(2,0) + C3*(Jinv(2,1)*Jinv(2,1) + Jinv(2,2)*Jinv(2,2));
  double CuuXY = C4*Jinv(2,0)*Jinv(2,1);
  double CuuXZ = C4*Jinv(2,0)*Jinv(2,2);
  double CuuYY = C1*Jinv(2,1)*Jinv(2,1) + C3*(Jinv(2,0)*Jinv(2,0) + Jinv(2,2)*Jinv(2,2));
  double CuuYZ = C4*Jinv(2,1)*Jinv(2,2);
  double CuuZZ = C1*Jinv(2,2)*Jinv(2,2) + C3*(Jinv(2,0)*Jinv(2,0) + Jinv(2,1)*Jinv(2,1));

  double CstXX = C1*Jinv(0,0)*Jinv(1,0) + C3*(Jinv(0,1)*Jinv(1,1) + Jinv(0,2)*Jinv(1,2));
  double CstXY = C2*Jinv(0,0)*Jinv(1,1) + C3*Jinv(0,1)*Jinv(1,0);
  double CstXZ = C2*Jinv(0,0)*Jinv(1,2) + C3*Jinv(0,2)*Jinv(1,0);
  double CstYX = C2*Jinv(0,1)*Jinv(1,0) + C3*Jinv(0,0)*Jinv(1,1);
  double CstYY = C1*Jinv(0,1)*Jinv(1,1) + C3*(Jinv(0,0)*Jinv(1,0) + Jinv(0,2)*Jinv(1,2));
  double CstYZ = C2*Jinv(0,1)*Jinv(1,2) + C3*Jinv(0,2)*Jinv(1,1);
  double CstZX = C2*Jinv(0,2)*Jinv(1,0) + C3*Jinv(0,0)*Jinv(1,2);
  double CstZY = C2*Jinv(0,2)*Jinv(1,1) + C3*Jinv(0,1)*Jinv(1,2);
  double CstZZ = C1*Jinv(0,2)*Jinv(1,2) + C3*(Jinv(0,0)*Jinv(1,0) + Jinv(0,1)*Jinv(1,1));

  double CsuXX = C1*Jinv(0,0)*Jinv(2,0) + C3*(Jinv(0,1)*Jinv(2,1) + Jinv(0,2)*Jinv(2,2));
  double CsuXY = C2*Jinv(0,0)*Jinv(2,1) + C3*Jinv(0,1)*Jinv(2,0);
  double CsuXZ = C2*Jinv(0,0)*Jinv(2,2) + C3*Jinv(0,2)*Jinv(2,0);
  double CsuYX = C2*Jinv(0,1)*Jinv(2,0) + C3*Jinv(0,0)*Jinv(2,1);
  double CsuYY = C1*Jinv(0,1)*Jinv(2,1) + C3*(Jinv(0,0)*Jinv(2,0) + Jinv(0,2)*Jinv(2,2));
  double CsuYZ = C2*Jinv(0,1)*Jinv(2,2) + C3*Jinv(0,2)*Jinv(2,1);
  double CsuZX = C2*Jinv(0,2)*Jinv(2,0) + C3*Jinv(0,0)*Jinv(2,2);
  double CsuZY = C2*Jinv(0,2)*Jinv(2,1) + C3*Jinv(0,1)*Jinv(2,2);
  double CsuZZ = C1*Jinv(0,2)*Jinv(2,2) + C3*(Jinv(0,0)*Jinv(2,0) + Jinv(0,1)*Jinv(2,1));

  double CtuXX = C1*Jinv(1,0)*Jinv(2,0) + C3*(Jinv(1,1)*Jinv(2,1) + Jinv(1,2)*Jinv(2,2));
  double CtuXY = C2*Jinv(1,0)*Jinv(2,1) + C3*Jinv(1,1)*Jinv(2,0);
  double CtuXZ = C2*Jinv(1,0)*Jinv(2,2) + C3*Jinv(1,2)*Jinv(2,0);
  double CtuYX = C2*Jinv(1,1)*Jinv(2,0) + C3*Jinv(1,0)*Jinv(2,1);
  double CtuYY = C1*Jinv(1,1)*Jinv(2,1) + C3*(Jinv(1,0)*Jinv(2,0) + Jinv(1,2)*Jinv(2,2));
  double CtuYZ = C2*Jinv(1,1)*Jinv(2,2) + C3*Jinv(1,2)*Jinv(2,1);
  double CtuZX = C2*Jinv(1,2)*Jinv(2,0) + C3*Jinv(1,0)*Jinv(2,2);
  double CtuZY = C2*Jinv(1,2)*Jinv(2,1) + C3*Jinv(1,1)*Jinv(2,2);
  double CtuZZ = C1*Jinv(1,2)*Jinv(2,2) + C3*(Jinv(1,0)*Jinv(2,0) + Jinv(1,1)*Jinv(2,1));

  Matrix FeCFe(9, 9), FeCFeInv(9, 9);
  // block11
  FeCFe(0,0) = CssXX*J0789;  FeCFe(0,1) = CssXY*J0789;  FeCFe(0,2) = CssXZ*J0789;
  FeCFe(1,0) = FeCFe(0,1);   FeCFe(1,1) = CssYY*J0789;  FeCFe(1,2) = CssYZ*J0789;
  FeCFe(2,0) = FeCFe(0,2);   FeCFe(2,1) = FeCFe(1,2);   FeCFe(2,2) = CssZZ*J0789;
  // block12
  FeCFe(0,3) = CstXX*J619;  FeCFe(0,4) = CstXY*J619;  FeCFe(0,5) = CstXZ*J619;
  FeCFe(1,3) = CstYX*J619;  FeCFe(1,4) = CstYY*J619;  FeCFe(1,5) = CstYZ*J619;
  FeCFe(2,3) = CstZX*J619;  FeCFe(2,4) = CstZY*J619;  FeCFe(2,5) = CstZZ*J619;
  // block13
  FeCFe(0,6) = CsuXX*J518;  FeCFe(0,7) = CsuXY*J518;  FeCFe(0,8) = CsuXZ*J518;
  FeCFe(1,6) = CsuYX*J518;  FeCFe(1,7) = CsuYY*J518;  FeCFe(1,8) = CsuYZ*J518;
  FeCFe(2,6) = CsuZX*J518;  FeCFe(2,7) = CsuZY*J518;  FeCFe(2,8) = CsuZZ*J518;
  // block21
  FeCFe(3,0) = CstXX*J619;  FeCFe(3,1) = CstYX*J619;  FeCFe(3,2) = CstZX*J619;
  FeCFe(4,0) = CstXY*J619;  FeCFe(4,1) = CstYY*J619;  FeCFe(4,2) = CstZY*J619;
  FeCFe(5,0) = CstXZ*J619;  FeCFe(5,1) = CstYZ*J619;  FeCFe(5,2) = CstZZ*J619;
  // block22
  FeCFe(3,3) = CttXX*J0879;  FeCFe(3,4) = CttXY*J0879;  FeCFe(3,5) = CttXZ*J0879;
  FeCFe(4,3) = FeCFe(3,4);   FeCFe(4,4) = CttYY*J0879;  FeCFe(4,5) = CttYZ*J0879;
  FeCFe(5,3) = FeCFe(3,5);   FeCFe(5,4) = FeCFe(4,5);   FeCFe(5,5) = CttZZ*J0879;
  // block23
  FeCFe(3,6) = CtuXX*J417;  FeCFe(3,7) = CtuXY*J417;  FeCFe(3,8) = CtuXZ*J417;
  FeCFe(4,6) = CtuYX*J417;  FeCFe(4,7) = CtuYY*J417;  FeCFe(4,8) = CtuYZ*J417;
  FeCFe(5,6) = CtuZX*J417;  FeCFe(5,7) = CtuZY*J417;  FeCFe(5,8) = CtuZZ*J417;
  // block31
  FeCFe(6,0) = CsuXX*J518;  FeCFe(6,1) = CsuYX*J518;  FeCFe(6,2) = CsuZX*J518;
  FeCFe(7,0) = CsuXY*J518;  FeCFe(7,1) = CsuYY*J518;  FeCFe(7,2) = CsuZY*J518;
  FeCFe(8,0) = CsuXZ*J518;  FeCFe(8,1) = CsuYZ*J518;  FeCFe(8,2) = CsuZZ*J518;
  // block32
  FeCFe(6,3) = CtuXX*J417;  FeCFe(6,4) = CtuYX*J417;  FeCFe(6,5) = CtuZX*J417;
  FeCFe(7,3) = CtuXY*J417;  FeCFe(7,4) = CtuYY*J417;  FeCFe(7,5) = CtuZY*J417;
  FeCFe(8,3) = CtuXZ*J417;  FeCFe(8,4) = CtuYZ*J417;  FeCFe(8,5) = CtuZZ*J417;
  // block33
  FeCFe(6,6) = CuuXX*J0978;  FeCFe(6,7) = CuuXY*J0978;  FeCFe(6,8) = CuuXZ*J0978;
  FeCFe(7,6) = FeCFe(6,7);   FeCFe(7,7) = CuuYY*J0978;  FeCFe(7,8) = CuuYZ*J0978;
  FeCFe(8,6) = FeCFe(6,8);   FeCFe(8,7) = FeCFe(7,8);   FeCFe(8,8) = CuuZZ*J0978;

  FeCFe.Invert(FeCFeInv);

  Matrix FeCFhg(9, 12);
  FeCFhg.Zero();
  // block11
  FeCFhg(0,0) = CstXX*J0789 + CssXX*J619;  FeCFhg(0,1) = CstXY*J0789 + CssXY*J619;  FeCFhg(0,2) = CstXZ*J0789 + CssXZ*J619;
  FeCFhg(1,0) = CstYX*J0789 + CssXY*J619;  FeCFhg(1,1) = CstYY*J0789 + CssYY*J619;  FeCFhg(1,2) = CstYZ*J0789 + CssYZ*J619;
  FeCFhg(2,0) = CstZX*J0789 + CssXZ*J619;  FeCFhg(2,1) = CstZY*J0789 + CssYZ*J619;  FeCFhg(2,2) = CstZZ*J0789 + CssZZ*J619;
  // block12
  FeCFhg(0,3) = CsuXX*J619 + CstXX*J518;  FeCFhg(0,4) = CsuXY*J619 + CstXY*J518;  FeCFhg(0,5) = CsuXZ*J619 + CstXZ*J518;
  FeCFhg(1,3) = CsuYX*J619 + CstYX*J518;  FeCFhg(1,4) = CsuYY*J619 + CstYY*J518;  FeCFhg(1,5) = CsuYZ*J619 + CstYZ*J518;
  FeCFhg(2,3) = CsuZX*J619 + CstZX*J518;  FeCFhg(2,4) = CsuZY*J619 + CstZY*J518;  FeCFhg(2,5) = CsuZZ*J619 + CstZZ*J518;
  // block13
  FeCFhg(0,6) = CsuXX*J0789 + CssXX*J518;  FeCFhg(0,7) = CsuXY*J0789 + CssXY*J518;  FeCFhg(0,8) = CsuXZ*J0789 + CssXZ*J518;
  FeCFhg(1,6) = CsuYX*J0789 + CssXY*J518;  FeCFhg(1,7) = CsuYY*J0789 + CssYY*J518;  FeCFhg(1,8) = CsuYZ*J0789 + CssYZ*J518;
  FeCFhg(2,6) = CsuZX*J0789 + CssXZ*J518;  FeCFhg(2,7) = CsuZY*J0789 + CssYZ*J518;  FeCFhg(2,8) = CsuZZ*J0789 + CssZZ*J518;
  // block21
  FeCFhg(3,0) = CstXX*J0879 + CttXX*J619;  FeCFhg(3,1) = CstYX*J0879 + CttXY*J619;  FeCFhg(3,2) = CstZX*J0879 + CttXZ*J619;
  FeCFhg(4,0) = CstXY*J0879 + CttXY*J619;  FeCFhg(4,1) = CstYY*J0879 + CttYY*J619;  FeCFhg(4,2) = CstZY*J0879 + CttYZ*J619;
  FeCFhg(5,0) = CstXZ*J0879 + CttXZ*J619;  FeCFhg(5,1) = CstYZ*J0879 + CttYZ*J619;  FeCFhg(5,2) = CstZZ*J0879 + CttZZ*J619;
  // block22
  FeCFhg(3,3) = CtuXX*J0879 + CttXX*J417;  FeCFhg(3,4) = CtuXY*J0879 + CttXY*J417;  FeCFhg(3,5) = CtuXZ*J0879 + CttXZ*J417;
  FeCFhg(4,3) = CtuYX*J0879 + CttXY*J417;  FeCFhg(4,4) = CtuYY*J0879 + CttYY*J417;  FeCFhg(4,5) = CtuYZ*J0879 + CttYZ*J417;
  FeCFhg(5,3) = CtuZX*J0879 + CttXZ*J417;  FeCFhg(5,4) = CtuZY*J0879 + CttYZ*J417;  FeCFhg(5,5) = CtuZZ*J0879 + CttZZ*J417;
  // block23
  FeCFhg(3,6) = CtuXX*J619 + CstXX*J417;  FeCFhg(3,7) = CtuXY*J619 + CstYX*J417;  FeCFhg(3,8) = CtuXZ*J619 + CstZX*J417;
  FeCFhg(4,6) = CtuYX*J619 + CstXY*J417;  FeCFhg(4,7) = CtuYY*J619 + CstYY*J417;  FeCFhg(4,8) = CtuYZ*J619 + CstZY*J417;
  FeCFhg(5,6) = CtuZX*J619 + CstXZ*J417;  FeCFhg(5,7) = CtuZY*J619 + CstYZ*J417;  FeCFhg(5,8) = CtuZZ*J619 + CstZZ*J417;
  // block31
  FeCFhg(6,0) = CtuXX*J518 + CsuXX*J417;  FeCFhg(6,1) = CtuYX*J518 + CsuYX*J417;  FeCFhg(6,2) = CtuZX*J518 + CsuZX*J417;
  FeCFhg(7,0) = CtuXY*J518 + CsuXY*J417;  FeCFhg(7,1) = CtuYY*J518 + CsuYY*J417;  FeCFhg(7,2) = CtuZY*J518 + CsuZY*J417;
  FeCFhg(8,0) = CtuXZ*J518 + CsuXZ*J417;  FeCFhg(8,1) = CtuYZ*J518 + CsuYZ*J417;  FeCFhg(8,2) = CtuZZ*J518 + CsuZZ*J417;
  // block32
  FeCFhg(6,3) = CtuXX*J0978 + CuuXX*J417;  FeCFhg(6,4) = CtuYX*J0978 + CuuXY*J417;  FeCFhg(6,5) = CtuZX*J0978 + CuuXZ*J417;
  FeCFhg(7,3) = CtuXY*J0978 + CuuXY*J417;  FeCFhg(7,4) = CtuYY*J0978 + CuuYY*J417;  FeCFhg(7,5) = CtuZY*J0978 + CuuYZ*J417;
  FeCFhg(8,3) = CtuXZ*J0978 + CuuXZ*J417;  FeCFhg(8,4) = CtuYZ*J0978 + CuuYZ*J417;  FeCFhg(8,5) = CtuZZ*J0978 + CuuZZ*J417;
  // block33
  FeCFhg(6,6) = CsuXX*J0978 + CuuXX*J518;  FeCFhg(6,7) = CsuYX*J0978 + CuuXY*J518;  FeCFhg(6,8) = CsuZX*J0978 + CuuXZ*J518;
  FeCFhg(7,6) = CsuXY*J0978 + CuuXY*J518;  FeCFhg(7,7) = CsuYY*J0978 + CuuYY*J518;  FeCFhg(7,8) = CsuZY*J0978 + CuuYZ*J518;
  FeCFhg(8,6) = CsuXZ*J0978 + CuuXZ*J518;  FeCFhg(8,7) = CsuYZ*J0978 + CuuYZ*J518;  FeCFhg(8,8) = CsuZZ*J0978 + CuuZZ*J518;

  Matrix KuT(12, 9);
  KuT = transpose(9, 12, FeCFhg);

  Matrix interior(12, 12);
  interior = FCF - KuT * FeCFeInv * FeCFhg;

  Matrix Kstab(24, 24);
  Kstab.Zero();
  Kstab.addMatrixTripleProduct(1.0, Mben, interior, 1.0);

  // persist the constant operators
  if (sspBnot == 0)  sspBnot  = new Matrix(6, 24);
  if (sspKstab == 0) sspKstab = new Matrix(24, 24);
  *sspBnot  = Bnot;
  *sspKstab = Kstab;
  sspVol    = mVol;
}

//----------------------------------------------------------------------
// formSSP — assemble the ssp tangent / residual from the constant operators:
//   K = Kstab + V * Bnot^T C Bnot          (C = current centroid tangent)
//   f = Kstab * u + V * Bnot^T * sigma - f_body
// The centroid material (materialPointers[0]) drives the constitutive update;
// strain was set in update() as Bnot*u. Body force uses the 2x2x2 N integral
// (consistent with the std/bbar paths).  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick::formSSP(int tang_flag, bool useInitialTangent)
{
  if (sspKstab == 0) buildSSP();   // safety (normally built in setDomain)

  computeBasis();
  const Vector &uCore = this->computeLocalDisp();   // identity for linear

  stiff.Zero();
  resid.Zero();
  static Vector bodyForce(24);   // external (global-frame) body load  // Ladruno (COROT-1)
  bodyForce.Zero();

  // Ladruno — Tier-A: degrade the constant elastic stabilization with damage
  // (1.0 for non-softening materials). The membrane terms below are NOT scaled —
  // they already carry the true (degraded) material C and stress; only the
  // artificial enhanced-strain Kstab is softened so a cracked element can localize.
  // The initial-tangent assembly keeps the full elastic Kstab (s=1) so the
  // predictor stiffness stays well-conditioned.
  const double sHg = useInitialTangent ? 1.0 : this->damageScale();

  if (tang_flag == 1) {
    const Matrix &C = useInitialTangent ? materialPointers[0]->getInitialTangent()
                                        : materialPointers[0]->getTangent();
    stiff = *sspKstab;
    stiff *= sHg;
    stiff.addMatrixTripleProduct(1.0, *sspBnot, C, sspVol);
  }

  if (!useInitialTangent) {
    // stabilization + membrane internal force
    resid.addMatrixVector(0.0, *sspKstab, uCore, sHg);
    const Vector &stress = materialPointers[0]->getStress();
    resid.addMatrixTransposeVector(1.0, *sspBnot, stress, sspVol);

    // body force: 2x2x2 Gauss integral of b * N_I (matches std/bbar)
    double xsj;
    static double shpBF[4][8];
    int cnt = 0;
    for (int gi = 0; gi < 2; gi++)
      for (int gj = 0; gj < 2; gj++)
        for (int gk = 0; gk < 2; gk++) {
          double gp[3] = {sg[gi], sg[gj], sg[gk]};
          shp3d(gp, xsj, shpBF, xl);
          double dvolBF = wg[cnt] * xsj;
          for (int j = 0; j < 8; j++)
            for (int p = 0; p < 3; p++) {
              // global-frame body load kept out of the corotated core force  // Ladruno (COROT-1)
              if (applyLoad == 0)
                bodyForce(3 * j + p) -= dvolBF * b[p] * shpBF[3][j];
              else
                bodyForce(3 * j + p) -= dvolBF * appliedB[p] * shpBF[3][j];
            }
          cnt++;
        }
  }

  // seam 3: globalize core-frame f/K to global DOFs (identity for -geom linear)
  static Vector zeroF(24);
  const Vector &fCore = useInitialTangent ? zeroF : resid;
  if (tang_flag == 1)
    theGeom->globalizeStiff(stiff, fCore, stiff);
  if (!useInitialTangent) {
    theGeom->globalizeForce(resid, resid);
    resid += bodyForce;          // global-frame body load, added AFTER globalize  // Ladruno (COROT-1)
  }
}

//----------------------------------------------------------------------
// eas (true Simo-Rifai) — buildEAStrue: cache the centroid Jacobian inverse and
// determinant for the (j0/j) sym[J0^-T E_i] enhanced-mode map. Geometry-only and
// (small strain) constant, so computed once in setDomain and rebuilt on the
// receive side. Mirrors the Jacobian-at-center step of buildSSP.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick::buildEAStrue(void)
{
  Matrix mNodeCrd(3, 8);
  for (int n = 0; n < 8; n++) {
    const Vector &crd = nodePointers[n]->getCrds();
    mNodeCrd(0, n) = crd(0);
    mNodeCrd(1, n) = crd(1);
    mNodeCrd(2, n) = crd(2);
  }

  // dN/dxi at the element center (Brick node order; trilinear ±0.125 pattern,
  // identical to buildSSP so easJ0inv matches the shp3d geometry used by computeB)
  Matrix dNloc(8, 3);
  dNloc(0,0)=-0.125; dNloc(1,0)= 0.125; dNloc(2,0)= 0.125; dNloc(3,0)=-0.125;
  dNloc(4,0)=-0.125; dNloc(5,0)= 0.125; dNloc(6,0)= 0.125; dNloc(7,0)=-0.125;
  dNloc(0,1)=-0.125; dNloc(1,1)=-0.125; dNloc(2,1)= 0.125; dNloc(3,1)= 0.125;
  dNloc(4,1)=-0.125; dNloc(5,1)=-0.125; dNloc(6,1)= 0.125; dNloc(7,1)= 0.125;
  dNloc(0,2)=-0.125; dNloc(1,2)=-0.125; dNloc(2,2)=-0.125; dNloc(3,2)=-0.125;
  dNloc(4,2)= 0.125; dNloc(5,2)= 0.125; dNloc(6,2)= 0.125; dNloc(7,2)= 0.125;

  Matrix J0(3, 3);
  J0 = mNodeCrd * dNloc;       // J0(i,j) = dx_i/dxi_j at the centroid
  J0.Invert(easJ0inv);         // easJ0inv(b,k) = dxi_b/dx_k
  easJ0det = J0(0,0) * (J0(1,1)*J0(2,2) - J0(1,2)*J0(2,1))
           - J0(0,1) * (J0(1,0)*J0(2,2) - J0(1,2)*J0(2,0))
           + J0(0,2) * (J0(1,0)*J0(2,1) - J0(1,1)*J0(2,0));
}

//----------------------------------------------------------------------
// eas — computeMenh: the 6x9 enhanced-strain operator M at a Gauss point
// (Voigt {xx,yy,zz,xy,yz,zx}, engineering shear). E9 = 3 natural-direction
// incompatible bubbles x 3 dofs (Wilson lineage = the 3-D lift of
// EnhancedQuad::computeBenhanced). For bubble b the physical mode gradient is
//   g_k = J0inv(b,k) * (xi_b / j)         (ADR 19 eq E.8, the one-sided J0^-T map)
// and column (3b+i) is the symmetric gradient of the displacement e_i * N_b.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick::computeMenh(const double gp[3], double jdet, Matrix &M)
{
  M.Zero();
  for (int b = 0; b < 3; b++) {            // bubble = natural direction xi/eta/zeta
    double scale = gp[b] / jdet;           // (j0 folded into dvol, as in EnhancedQuad)
    double g[3];
    for (int k = 0; k < 3; k++)
      g[k] = easJ0inv(b, k) * scale;
    for (int i = 0; i < 3; i++) {          // dof = displacement direction of the bubble
      int col = 3 * b + i;
      M(0, col) = (i == 0) ? g[0] : 0.0;   // eps_xx = du_x/dx
      M(1, col) = (i == 1) ? g[1] : 0.0;   // eps_yy
      M(2, col) = (i == 2) ? g[2] : 0.0;   // eps_zz
      M(3, col) = ((i == 0) ? g[1] : 0.0) + ((i == 1) ? g[0] : 0.0);  // gamma_xy
      M(4, col) = ((i == 1) ? g[2] : 0.0) + ((i == 2) ? g[1] : 0.0);  // gamma_yz
      M(5, col) = ((i == 2) ? g[0] : 0.0) + ((i == 0) ? g[2] : 0.0);  // gamma_zx
    }
  }
}

//----------------------------------------------------------------------
// eas — formEAStrue: true Simo-Rifai assembly. Full 2x2x2 integration; at each GP
// the total strain is eps = B*u + M*alpha. An inner Newton solves the 9 enhanced
// parameters alpha so the internal condition int M^T sigma = 0 holds (d fixed),
// then the parameters are statically condensed:
//   K* = Kdd - Kda Kaa^-1 Kad,   f* = int B^T sigma   (the inner Newton drives
// h = int M^T sigma -> 0 so the f* condensation term vanishes — the EnhancedQuad
// contract). Routed through the seam-3 globalize (identity for -geom linear). The
// useInitialTangent path condenses the INITIAL elastic tangent at alpha=0 (no
// Newton), mirroring EnhancedQuad::getInitialStiff.  // Ladruno
//----------------------------------------------------------------------
void
LadrunoBrick::formEAStrue(int tang_flag, bool useInitialTangent)
{
  static const int    maxIters = 12;
  // RELATIVE inner-Newton tolerance: converge when ||int M^T sigma|| drops to
  // tolRel of its first-iteration value (+ a tiny absolute floor for the case where
  // alpha already satisfies it, e.g. a uniform-strain state where the first residual
  // is ~0). An ABSOLUTE tol fails here — ||int M^T sigma|| has units of force×length²,
  // so its converged (machine-relative) magnitude scales with E and the mesh size; a
  // fixed 1e-10 is unreachable for steel-scale problems (spurious non-convergence).  // Ladruno
  static const double tolRel = 1.0e-8;
  static const double tolAbs = 1.0e-12;

  if (easJ0det == 0.0) buildEAStrue();       // safety (normally cached in setDomain)

  computeBasis();
  const Vector &uCore = this->computeLocalDisp();   // identity for -geom linear

  stiff.Zero();
  resid.Zero();

  // precompute shape functions, |J|, dvol and natural coords at the 8 GPs
  static double Shape[4][8][8];
  static double dvol[8];
  static double jdetGP[8];
  static double gpNat[8][3];
  {
    double xsj;
    static double shpTmp[4][8];
    int cnt = 0;
    for (int gi = 0; gi < 2; gi++)
      for (int gj = 0; gj < 2; gj++)
        for (int gk = 0; gk < 2; gk++) {
          gpNat[cnt][0] = sg[gi]; gpNat[cnt][1] = sg[gj]; gpNat[cnt][2] = sg[gk];
          double gp[3] = {sg[gi], sg[gj], sg[gk]};
          shp3d(gp, xsj, shpTmp, xl);
          for (int p = 0; p < 4; p++)
            for (int q = 0; q < 8; q++)
              Shape[p][q][cnt] = shpTmp[p][q];
          jdetGP[cnt] = xsj;
          dvol[cnt]   = wg[cnt] * xsj;
          cnt++;
        }
  }

  static Matrix B(6, 24), M(6, 9);
  static Matrix dd(6, 6), DB(6, 24), DM(6, 9);
  static Matrix Kaa(9, 9), Kda(24, 9), Kad(9, 24), KaaInvKad(9, 24);
  static Vector residE(9), dalpha(9), strain(6), stress(6);
  static double shp[4][8];

  // helper macro-ish: build the compatible 6x24 B at the current saved GP `g`
  // (filled inline below to avoid a member; computeB returns the 6x3 node block)

  // -------- getInitialStiff: condensed initial elastic tangent at alpha=0 --------
  if (useInitialTangent) {
    Kaa.Zero(); Kda.Zero(); Kad.Zero();
    for (int g = 0; g < 8; g++) {
      for (int p = 0; p < 4; p++) for (int q = 0; q < 8; q++) shp[p][q] = Shape[p][q][g];
      B.Zero();
      for (int node = 0; node < 8; node++) {
        const Matrix &Bn = computeB(node, shp);
        for (int r = 0; r < 6; r++) for (int c = 0; c < 3; c++) B(r, 3*node + c) = Bn(r, c);
      }
      computeMenh(gpNat[g], jdetGP[g], M);
      dd = materialPointers[g]->getInitialTangent();
      dd *= dvol[g];
      DB.addMatrixProduct(0.0, dd, B, 1.0);
      DM.addMatrixProduct(0.0, dd, M, 1.0);
      stiff.addMatrixTransposeProduct(1.0, B, DB, 1.0);   // Kdd += B^T dd B
      Kda.addMatrixTransposeProduct(1.0, B, DM, 1.0);     //     += B^T dd M
      Kad.addMatrixTransposeProduct(1.0, M, DB, 1.0);     //     += M^T dd B
      Kaa.addMatrixTransposeProduct(1.0, M, DM, 1.0);     //     += M^T dd M
    }
    Kaa.Solve(Kad, KaaInvKad);
    stiff.addMatrixProduct(1.0, Kda, KaaInvKad, -1.0);    // K* = Kdd - Kda Kaa^-1 Kad
    static Vector zeroF(24);
    theGeom->globalizeStiff(stiff, zeroF, stiff);
    return;
  }

  // -------- inner Newton: solve alpha s.t. int M^T sigma = 0 (d fixed) --------
  int count = 0;
  double r0 = -1.0;          // first-iteration ||residE|| (sets the relative scale)
  while (true) {
    residE.Zero();
    Kaa.Zero();
    for (int g = 0; g < 8; g++) {
      for (int p = 0; p < 4; p++) for (int q = 0; q < 8; q++) shp[p][q] = Shape[p][q][g];
      B.Zero();
      for (int node = 0; node < 8; node++) {
        const Matrix &Bn = computeB(node, shp);
        for (int r = 0; r < 6; r++) for (int c = 0; c < 3; c++) B(r, 3*node + c) = Bn(r, c);
      }
      computeMenh(gpNat[g], jdetGP[g], M);
      strain.addMatrixVector(0.0, B, uCore, 1.0);     // eps = B*u
      strain.addMatrixVector(1.0, M, alpha, 1.0);     //     + M*alpha
      materialPointers[g]->setTrialStrain(strain);
      stress = materialPointers[g]->getStress();
      stress *= dvol[g];
      dd = materialPointers[g]->getTangent();
      dd *= dvol[g];
      residE.addMatrixTransposeVector(1.0, M, stress, -1.0);   // residE += -(M^T sigma)
      DM.addMatrixProduct(0.0, dd, M, 1.0);
      Kaa.addMatrixTransposeProduct(1.0, M, DM, 1.0);          // Kaa += M^T dd M
    }
    double r = residE.Norm();
    if (count == 0) r0 = r;                          // scale from the first residual
    // converged (after >=1 update): relative to r0, with a tiny absolute floor for
    // the alpha-already-satisfies (uniform-strain) case where r0 itself is ~0.
    if (count >= 1 && r <= tolRel * r0 + tolAbs)
      break;
    if (count >= maxIters) {
      opserr << "LadrunoBrick::formEAStrue - element " << this->getTag()
             << ": enhanced-strain Newton did not converge in " << maxIters
             << " iters (||r||=" << r << ", r0=" << r0 << ")\n";
      break;
    }
    dalpha.Zero();
    Kaa.Solve(residE, dalpha);     // dalpha = -Kaa^-1 h  (residE = -h)
    alpha += dalpha;
    count++;
  }

  // -------- final assembly at converged alpha (Kaa preserved from the loop) ------
  static Vector bodyForce(24);
  bodyForce.Zero();
  Kda.Zero(); Kad.Zero();
  for (int g = 0; g < 8; g++) {
    for (int p = 0; p < 4; p++) for (int q = 0; q < 8; q++) shp[p][q] = Shape[p][q][g];
    B.Zero();
    for (int node = 0; node < 8; node++) {
      const Matrix &Bn = computeB(node, shp);
      for (int r = 0; r < 6; r++) for (int c = 0; c < 3; c++) B(r, 3*node + c) = Bn(r, c);
    }
    computeMenh(gpNat[g], jdetGP[g], M);
    stress = materialPointers[g]->getStress();
    stress *= dvol[g];
    resid.addMatrixTransposeVector(1.0, B, stress, 1.0);   // f_int += B^T sigma
    // body / self-weight (2x2x2 N integral), kept OUT of the corotated core force
    for (int node = 0; node < 8; node++)
      for (int pp = 0; pp < 3; pp++) {
        if (applyLoad == 0) bodyForce(3*node + pp) -= dvol[g] * b[pp]       * shp[3][node];
        else                bodyForce(3*node + pp) -= dvol[g] * appliedB[pp] * shp[3][node];
      }
    if (tang_flag == 1) {
      dd = materialPointers[g]->getTangent();
      dd *= dvol[g];
      DB.addMatrixProduct(0.0, dd, B, 1.0);
      DM.addMatrixProduct(0.0, dd, M, 1.0);
      stiff.addMatrixTransposeProduct(1.0, B, DB, 1.0);   // Kdd
      Kda.addMatrixTransposeProduct(1.0, B, DM, 1.0);
      Kad.addMatrixTransposeProduct(1.0, M, DB, 1.0);
    }
  }

  if (tang_flag == 1) {
    Kaa.Solve(Kad, KaaInvKad);                            // Kaa from the inner loop
    stiff.addMatrixProduct(1.0, Kda, KaaInvKad, -1.0);    // K* = Kdd - Kda Kaa^-1 Kad
  }

  // seam 3: globalize core-frame f/K to global DOFs (identity for -geom linear)
  if (tang_flag == 1)
    theGeom->globalizeStiff(stiff, resid, stiff);
  theGeom->globalizeForce(resid, resid);
  resid += bodyForce;
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

  static Vector dData(9);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];
  dData(7) = hourglassCoeff;
  dData(8) = hgDissipated;   // Ladruno — viscous-hourglass dissipation accumulator

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoBrick::sendSelf() - failed to send double data\n";
    return -1;
  }

  // Ladruno — eas: ship the committed enhanced parameters as an EXTRA guarded
  // Vector(9), only for -formulation eas. std/bbar/uri/ssp send nothing here, so
  // their stream stays byte-identical. (idData(25) is the Rayleigh-damping flag,
  // NOT a version slot — the formulation ordinal in idData(28) gates the recv.)
  if (formulation == Formulation::EAS) {
    if (theChannel.sendVector(dataTag, commitTag, alphaCommit) < 0) {
      opserr << "LadrunoBrick::sendSelf() - failed to send eas alpha\n";
      return -1;
    }
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

  // Ladruno — Tier-A: drop any cached "damage" query before the materials below
  // may be deleted/rebrokered. It holds a raw pointer to materialPointers[0], so a
  // stale one would dangle on a material-class change; setDomain rebuilds it
  // against the live material on the receive side.  // Ladruno
  if (damageResponse) { delete damageResponse; damageResponse = 0; }

  static Vector dData(9);
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

  // Ladruno — restore the dissipation total; the work baseline reseeds on the
  // next commit (hgPrevValid=false), so no spurious increment is booked across
  // the migration. uPrevCommit was zeroed by the null constructor.
  hgDissipated = dData(8);
  hgPrevValid  = false;

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

  // Ladruno — eas: receive the committed enhanced parameters (sent right after
  // dData, before the materials — match that order here). Gated on the formulation
  // ordinal just decoded; non-eas streams shipped nothing extra. alpha = alphaCommit
  // = received (no trial state survives a migration; the next form pass re-solves).
  if (formulation == Formulation::EAS) {
    if (theChannel.recvVector(dataTag, commitTag, alphaCommit) < 0) {
      opserr << "LadrunoBrick::recvSelf() - failed to recv eas alpha\n";
      return -1;
    }
    alpha = alphaCommit;
  }

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
    const bool sp = this->isSinglePoint();   // single-point: mirror slot 0
    for (int n = 0; n < 8; n++) {
      const Vector &stressN = materialPointers[sp ? 0 : n]->getStress();
      values(n) = stressN(index);
    }
  } else if (displayMode < 0) {
    for (int n = 0; n < 8; n++)
      values(n) = 0.0;
  }

  return theViewer.drawCube(coords, values, this->getTag());
}

//----------------------------------------------------------------------
// Work done against the Flanagan-Belytschko VISCOUS hourglass damping force
// over the step now being committed (uri + Hourglass::VISCOUS only). The
// damping force stores no energy, so integrating its work is the only way to
// report what the viscous hourglass control dissipated. Mirrors the force
// assembled in formResidAndTangent's viscous branch:
//   f_iI = c_visc·Σ_a γ_aI·q̇_ai ,  c_visc = epsV·ρ·c_d·V^(2/3), c_d=√(dd(0,0)/ρ)
// and returns ΔE = Σ_iI f_iI·Δu_iI = c_visc·Σ_a,i q̇_ai·Δq_ai, with q̇ from the
// trial velocity and Δq = Σ_J γ_aJ·(u_now−uPrevCommit)_iJ. Clamped ≥ 0
// (dissipation is non-negative; drops sub-step numerical sign noise). Returns 0
// when the force itself is 0 (no mass).  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick::viscousHourglassIncrement(void)
{
  static const int numberNodes = 8;

  double rho = materialPointers[0]->getRho();
  if (rho <= 0.0)
    return 0.0;

  computeBasis();
  double gp[3] = {0.0, 0.0, 0.0};
  double detJ;
  static double shpC[4][8];
  shp3d(gp, detJ, shpC, xl);
  const double vol = 8.0 * detJ;

  double bC[3][8];
  for (int i = 0; i < 3; i++)
    for (int I = 0; I < numberNodes; I++)
      bC[i][I] = shpC[i][I];

  const Matrix &dd = materialPointers[0]->getTangent();
  const double cd = sqrt(dd(0, 0) / rho);
  const double epsV = (hourglassCoeff > 0.0) ? hourglassCoeff : 0.1;
  const double cVisc = epsV * rho * cd * pow(vol, 2.0 / 3.0);

  double gamma[4][8];
  computeGammaHourglass(bC, gamma, 1.0);

  // generalized hourglass velocities q̇_ai and displacement increments Δq_ai
  double qd[4][3], dq[4][3];
  for (int a = 0; a < 4; a++)
    for (int i = 0; i < 3; i++) {
      qd[a][i] = 0.0;
      dq[a][i] = 0.0;
    }
  for (int J = 0; J < numberNodes; J++) {
    const Vector &vl = nodePointers[J]->getTrialVel();
    const Vector &ul = nodePointers[J]->getTrialDisp();
    for (int a = 0; a < 4; a++)
      for (int i = 0; i < 3; i++) {
        qd[a][i] += gamma[a][J] * vl(i);
        dq[a][i] += gamma[a][J] * (ul(i) - uPrevCommit[3 * J + i]);
      }
  }

  double dE = 0.0;
  for (int a = 0; a < 4; a++)
    for (int i = 0; i < 3; i++)
      dE += qd[a][i] * dq[a][i];
  dE *= cVisc;

  return (dE > 0.0) ? dE : 0.0;
}

//----------------------------------------------------------------------
// Reference-config element volume: V = Σ_g wg·detJ over the 8 (2x2x2) Gauss
// points, with detJ from shp3d at the reference nodal coordinates (computeBasis
// sets xl = reference coords). Formulation-independent and available before
// buildSSP, so getCharacteristicLength works for every -formulation and for
// -geom finite (crack-band regularization wants the reference/material size).  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick::computeVolume(void)
{
  computeBasis();   // xl = reference nodal coordinates

  double xsj;
  static double shp[4][8];
  double gp[3];
  double V = 0.0;

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        gp[0] = sg[i]; gp[1] = sg[j]; gp[2] = sg[k];
        shp3d(gp, xsj, shp, xl);
        V += wg[count] * xsj;          // wg[count] = 1 for the 2x2x2 rule
        count++;
      }
    }
  }
  return V;
}

//----------------------------------------------------------------------
// Characteristic length for crack-band / regularized-softening materials.
// lch = cbrt(V), the edge of a cube of equal volume — recovers the true edge
// for an undistorted hex and stays geometry-true for distorted ones. The base
// default (min inter-node distance) under-sizes the band on a skewed/graded hex
// and over-softens. The cbrt-of-volume convention matches BezierTet10 (cbrt(6V),
// the factor being the reciprocal reference-element volume: tet = 1/6, hex = 1)
// and BezierTri6 (sqrt(2A)). Degenerate (V<=0) falls back to the base.  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick::getCharacteristicLength(void)
{
  double V = this->computeVolume();
  if (V <= 0.0)
    return Element::getCharacteristicLength();   // degenerate: fall back
  return cbrt(V);
}

//----------------------------------------------------------------------
// getInterpolationWeights — trilinear 8-node hex shape weights at the natural
// coordinate xi = (ξ,η,ζ). Same formula as the body-force / mass interpolation
// used throughout this element (N_I = 0.125·∏(1+ξ_I·ξ_k), natCoord = nodal
// corner signs). Lets LadrunoEmbeddedRebar embed a rebar node in this host
// without re-supplying the weights (ADR 20 §9). N is resized to 8.  // Ladruno
//----------------------------------------------------------------------
int
LadrunoBrick::getInterpolationWeights(const Vector &xi, Vector &N)
{
  if (xi.Size() < 3) {
    opserr << "LadrunoBrick::getInterpolationWeights - xi needs 3 natural "
              "coords (xi,eta,zeta)\n";
    return -1;
  }
  if (N.Size() != numberNodes)
    N.resize(numberNodes);
  for (int I = 0; I < numberNodes; I++)
    N(I) = 0.125 * (1.0 + natCoord[I][0]*xi(0))
                 * (1.0 + natCoord[I][1]*xi(1))
                 * (1.0 + natCoord[I][2]*xi(2));
  return 0;
}

//----------------------------------------------------------------------
// Recoverable elastic hourglass / stabilization energy at the current trial
// state. uri 'stiffness': ½κ·Σ q_aι² (the FB perturbation energy; q from the
// trial displacement, mirrors formUri exactly). ssp: ½·u_core·Kstab·u_core.
// 0 for std/bbar (full integration) and physical (assumed strain folded into
// the strain energy). For uri 'viscous' (dissipative, not stored) this returns
// the CUMULATIVE committed dissipation hgDissipated instead.  // Ladruno
//----------------------------------------------------------------------
double
LadrunoBrick::hourglassEnergy(void)
{
  static const int numberNodes = 8;

  // uri 'viscous' dissipates rather than stores: report the running total.
  if (formulation == Formulation::URI && hourglassType == Hourglass::VISCOUS)
    return hgDissipated;

  if (formulation == Formulation::SSP) {
    if (sspKstab == 0) buildSSP();
    const Vector &uCore = this->computeLocalDisp();   // identity for linear
    static Vector Ku(24);
    Ku.addMatrixVector(0.0, *sspKstab, uCore, 1.0);
    double e = 0.0;
    for (int i = 0; i < 24; i++) e += uCore(i) * Ku(i);
    // Ladruno — report the ACTUALLY stored (damage-degraded) stabilization energy,
    // consistent with the Kstab formSSP assembles; the instrument for tuning the
    // floor against internal energy (§5). 1.0 for non-softening materials.
    return 0.5 * e * this->damageScale();
  }

  if (formulation == Formulation::URI && hourglassType == Hourglass::STIFFNESS) {
    computeBasis();
    double gp[3] = {0.0, 0.0, 0.0};
    double detJ;
    static double shpC[4][8];
    shp3d(gp, detJ, shpC, xl);
    const double vol = 8.0 * detJ;

    double bC[3][8];
    for (int i = 0; i < 3; i++)
      for (int I = 0; I < numberNodes; I++)
        bC[i][I] = shpC[i][I];

    double bb = 0.0;
    for (int i = 0; i < 3; i++)
      for (int I = 0; I < numberNodes; I++)
        bb += bC[i][I] * bC[i][I];
    const double Ghg = materialPointers[0]->getTangent()(3, 3);
    const double scale = (hourglassCoeff > 0.0) ? hourglassCoeff : 0.05;
    double kappa = scale * Ghg * vol * bb;
    // Ladruno — mirror formUri's Tier-A κ exactly (elastic base × floored damage
    // scale for a softening material) so the reported stored stabilization energy
    // matches what the element actually assembles. Single application of the scale.
    // (already inside the URI+STIFFNESS branch, so only the damage guard remains.)
    if (damageResponse != 0) {
      const double G0 = materialPointers[0]->getInitialTangent()(3, 3);
      kappa = scale * G0 * vol * bb * this->damageScale();
    }

    double gamma[4][8];
    computeGammaHourglass(bC, gamma, 1.0);

    // generalized hourglass strains q_aι = Σ_J γ_aJ u_iJ (trial displacement)
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
    double e = 0.0;
    for (int a = 0; a < 4; a++)
      for (int i = 0; i < 3; i++)
        e += q[a][i] * q[a][i];
    return 0.5 * kappa * e;   // kappa already carries the Tier-A damage scale
  }

  return 0.0;   // std / bbar / physical (uri-viscous handled above)
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

  } else if (strcmp(argv[0], "stiff") == 0 || strcmp(argv[0], "stiffness") == 0) {
    // full 24x24 tangent — used by the FD-tangent verification (and handy generally).
    theResponse = new ElementResponse(this, 2, Matrix(24, 24));

  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8) {
      // single-point formulations: every GP query resolves to the live slot 0
      int slot = this->isSinglePoint() ? 0 : pointNum - 1;
      output.tag("GaussPoint");
      output.attr("number", pointNum);
      theResponse = materialPointers[slot]->setResponse(&argv[2], argc - 2, output);
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

  } else if (strcmp(argv[0], "hourglassEnergy") == 0 ||
             strcmp(argv[0], "hgEnergy") == 0 ||
             strcmp(argv[0], "hourglassWork") == 0 ||
             strcmp(argv[0], "hourglassDissipation") == 0 ||
             strcmp(argv[0], "hgDissipation") == 0) {
    // For uri 'stiffness'/ssp this is the STORED stabilization energy; for uri
    // 'viscous' it is the CUMULATIVE DISSIPATED energy (see hourglassEnergy()).
    output.tag("ResponseType", "Ehg");
    theResponse = new ElementResponse(this, 8, Vector(1));

  } else if (strcmp(argv[0], "charLength") == 0 ||
             strcmp(argv[0], "characteristicLength") == 0) {
    // Ladruno — the element size handed to crack-band materials (= cbrt(V)).
    output.tag("ResponseType", "lch");
    theResponse = new ElementResponse(this, 9, Vector(1));
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
    // single-point formulations: only slot 0 is live → mirror it to all 8 GPs
    const bool sp = this->isSinglePoint();
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      const Vector &sigma = materialPointers[sp ? 0 : i]->getStress();
      for (int j = 0; j < 6; j++) stresses(cnt++) = sigma(j);
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 4) {
    const bool sp = this->isSinglePoint();
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      const Vector &eps = materialPointers[sp ? 0 : i]->getStrain();
      for (int j = 0; j < 6; j++) stresses(cnt++) = eps(j);
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 6) {
    Vector tmpStress(6);
    if (this->isSinglePoint()) {           // uniform: average == centroid (slot 0)
      tmpStress = materialPointers[0]->getStress();
    } else {
      for (int i = 0; i < 8; i++) {
        const Vector &sigma = materialPointers[i]->getStress();
        for (int j = 0; j < 6; j++) tmpStress(j) += sigma(j) * 0.125;
      }
    }
    return eleInfo.setVector(tmpStress);

  } else if (responseID == 7) {
    Vector tmpStrain(6);
    if (this->isSinglePoint()) {
      tmpStrain = materialPointers[0]->getStrain();
    } else {
      for (int i = 0; i < 8; i++) {
        const Vector &eps = materialPointers[i]->getStrain();
        for (int j = 0; j < 6; j++) tmpStrain(j) += eps(j);
      }
      tmpStrain /= 8.0;
    }
    return eleInfo.setVector(tmpStrain);

  } else if (responseID == 8) {
    static Vector ehg(1);
    ehg(0) = this->hourglassEnergy();
    return eleInfo.setVector(ehg);

  } else if (responseID == 9) {
    static Vector lch(1);
    lch(0) = this->getCharacteristicLength();
    return eleInfo.setVector(lch);
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
    if (pointNum > 0 && pointNum <= 8) {
      int slot = this->isSinglePoint() ? 0 : pointNum - 1;   // single-point: live slot 0
      return materialPointers[slot]->setParameter(&argv[2], argc - 2, param);
    } else
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
