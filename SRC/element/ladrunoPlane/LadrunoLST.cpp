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

// LadrunoLST — 6-node linear-strain triangle (T6, Ladruno fork, ADR 70 P3).
// std -geom linear mirrors upstream SixNodeTri (same quadratic shape functions,
// same 3-point interior rule) so it reduces to it to ~1e-9; -geom finite (std
// only) routes through the shared LadrunoFiniteStrain2DKernel. NO bbar/F-bar
// lane — constant element-mean dilatation is RANK-DEFICIENT on the T6 (two
// conformal zero-energy modes; refuted at P3, see LadrunoLST.h and ADR 70 §9).
// See Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md

#include <LadrunoLST.h>
#include <LadrunoFiniteStrain2DKernel.h>   // Ladruno (ADR 70): shared 2D finite-strain kernel
#include <FiniteStrainND2DMaterial.h>      // Ladruno (ADR 70): setTrialF(F) seam
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
#include <LadrunoResponseTokens.h>   // Ladruno — shared recorder-token aliases
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <math.h>

double LadrunoLST::matrixData[(2 * 6) * (2 * 6)];
Matrix LadrunoLST::K(matrixData, 12, 12);
Vector LadrunoLST::P(12);
double LadrunoLST::shp[3][6];
double LadrunoLST::pts[3][2];
double LadrunoLST::wts[3];

// 3-point interior rule (barycentric (2/3,1/6),(1/6,2/3),(1/6,1/6), w = 1/6) —
// integrates the quadratic strain energy exactly and matches upstream SixNodeTri
// (the reduce-to gate); the alternative midside rule is less accurate.
static void lst_init_rule(double pts[3][2], double wts[3])
{
  pts[0][0] = 0.666666666666666667;
  pts[0][1] = 0.166666666666666667;
  pts[1][0] = 0.166666666666666667;
  pts[1][1] = 0.666666666666666667;
  pts[2][0] = 0.166666666666666667;
  pts[2][1] = 0.166666666666666667;
  wts[0] = 0.166666666666666667;
  wts[1] = 0.166666666666666667;
  wts[2] = 0.166666666666666667;
}

LadrunoLST::LadrunoLST(int tag, int nd1, int nd2, int nd3,
                       int nd4, int nd5, int nd6,
                       NDMaterial &m, const char *type, double t,
                       Formulation form, Geom g,
                       double r, double b1, double b2, double p,
                       double b1bv, double b2bv)
 :Element(tag, ELE_TAG_LadrunoLST),
  theMaterial(0), connectedExternalNodes(6),
  Q(12), pressureLoad(12), thickness(t), pressure(p), rho(r),
  formulation(form), geom(g), bulkVisc_b1(b1bv), bulkVisc_b2(b2bv),
  planeType(1), Ki(0)
{
  lst_init_rule(pts, wts);

  if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)
    planeType = 1;
  else if (strcmp(type, "PlaneStress") == 0 || strcmp(type, "PlaneStress2D") == 0)
    planeType = 2;
  else {
    opserr << "LadrunoLST::LadrunoLST -- improper material type: " << type << "\n";
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
      opserr << "LadrunoLST::LadrunoLST -- failed to get a copy of material model\n";
      exit(-1);
    }
  }

  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;
  connectedExternalNodes(3) = nd4;
  connectedExternalNodes(4) = nd5;
  connectedExternalNodes(5) = nd6;

  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoLST::LadrunoLST()
 :Element(0, ELE_TAG_LadrunoLST),
  theMaterial(0), connectedExternalNodes(6),
  Q(12), pressureLoad(12), thickness(0.0), pressure(0.0), rho(0.0),
  formulation(Formulation::STD), geom(Geom::LINEAR),
  bulkVisc_b1(0.0), bulkVisc_b2(0.0), planeType(1), Ki(0)
{
  lst_init_rule(pts, wts);
  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b[1] = 0.0;
  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoLST::~LadrunoLST()
{
  if (theMaterial) {
    for (int i = 0; i < numgp; i++)
      if (theMaterial[i]) delete theMaterial[i];
    delete[] theMaterial;
  }
  if (Ki) delete Ki;
}

const char *LadrunoLST::typeString(void) const
{
  return (planeType == 2) ? "PlaneStress" : "PlaneStrain";
}

int LadrunoLST::getNumExternalNodes(void) const { return numnodes; }
const ID &LadrunoLST::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoLST::getNodePtrs(void) { return theNodes; }
int LadrunoLST::getNumDOF(void) { return 2 * numnodes; }

void LadrunoLST::setDomain(Domain *theDomain)
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
      opserr << "LadrunoLST::setDomain -- node " << connectedExternalNodes(i)
             << " must have 2 DOFs\n";
      return;
    }
  this->DomainComponent::setDomain(theDomain);

  // Ladruno (ADR 70): the factory already refuses -geom finite with a
  // non-FiniteStrainND2DMaterial; this re-check covers broker-built (recvSelf)
  // elements. updateFinite's dynamic_cast is what actually keeps a misuse from
  // casting UB — this just fails LOUD at setup instead of on the first update.
  if (this->isFinite()) {
    for (int i = 0; i < numgp; i++)
      if (dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[i]) == 0) {
        opserr << "LadrunoLST::setDomain -- WARNING element " << this->getTag()
               << " uses -geom finite but its material is not a "
                  "FiniteStrainND2DMaterial (e.g. LogStrain2D); the analysis will fail\n";
        break;
      }
  }

  this->setPressureLoadAtNodes();
}

double LadrunoLST::getCharacteristicLength(void)
{
  // crack-band size for a triangle = sqrt(2*area); area = Σ detJ·w over the
  // 3-point rule (Σw = ½). Matches the LadrunoCST / BezierTri6 convention.
  double A = 0.0;
  for (int i = 0; i < numgp; i++)
    A += this->shapeFunction(pts[i][0], pts[i][1]) * wts[i];
  if (A <= 0.0)
    return this->Element::getCharacteristicLength();
  return sqrt(2.0 * A);
}

int LadrunoLST::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0)
    opserr << "LadrunoLST::commitState() - failed in base class\n";
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->commitState();
  return retVal;
}

int LadrunoLST::revertToLastCommit(void)
{
  int retVal = 0;
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->revertToLastCommit();
  return retVal;
}

int LadrunoLST::revertToStart(void)
{
  int retVal = 0;
  for (int i = 0; i < numgp; i++)
    retVal += theMaterial[i]->revertToStart();
  return retVal;
}

// standard displacement B (3x12). No B-bar variant: constant element-mean
// dilatation is rank-deficient on the T6 (conformal zero-energy modes) — see
// the header note; the volumetric cure is ADR 70 P4, not an element-local
// average.
void LadrunoLST::formB(Matrix &B)
{
  B.Zero();
  for (int a = 0; a < numnodes; a++) {
    double bx = shp[0][a];
    double by = shp[1][a];
    B(0, 2 * a)     = bx;
    B(1, 2 * a + 1) = by;
    B(2, 2 * a)     = by;
    B(2, 2 * a + 1) = bx;
  }
}

int LadrunoLST::update(void)
{
  if (this->isFinite())            // -geom finite (ADR 70): F-driven updated-Lagrangian
    return this->updateFinite();

  static Vector u(12);
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u(2 * a)     = d(0);
    u(2 * a + 1) = d(1);
  }

  static Vector eps(3);
  static Matrix B(3, 12);
  int ret = 0;
  for (int i = 0; i < numgp; i++) {
    this->shapeFunction(pts[i][0], pts[i][1]);
    this->formB(B);
    eps.addMatrixVector(0.0, B, u, 1.0);
    ret += theMaterial[i]->setTrialStrain(eps);
  }
  return ret;
}

const Matrix &LadrunoLST::getTangentStiff(void)
{
  if (this->isFinite()) {          // -geom finite (ADR 70): consistent spatial tangent
    this->formFinite(1);
    return K;
  }

  K.Zero();

  static Matrix B(3, 12);
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getTangent();
    this->formB(B);
    K.addMatrixTripleProduct(1.0, B, D, dvol);
  }
  return K;
}

const Matrix &LadrunoLST::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  // Ladruno (ADR 70): under -geom finite this deliberately stays the SYMMETRIC
  // small-strain reference tangent Bᵀ D₀ B dV — at F = I it equals the finite
  // consistent tangent (σ = 0, no F-bar lane on the LST), and it is a
  // well-posed symmetric seed for Rayleigh βK0 / -initial / eigen use.
  // Mirrors LadrunoCST.
  K.Zero();

  static Matrix B(3, 12);
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    const Matrix &D = theMaterial[i]->getInitialTangent();
    this->formB(B);
    K.addMatrixTripleProduct(1.0, B, D, dvol);
  }
  Ki = new Matrix(K);
  return *Ki;   // return the cached copy, not the shared static scratch K
}

const Matrix &LadrunoLST::getMass(void)
{
  // Ladruno (ADR-77 G2 ext): per-instance mass cache -- LadrunoMassCache.h.
  // Signature = rho override + numgp material rhos + thickness; coords
  // guarded inside.
  double mcSig[2 + numgp];
  mcSig[0] = rho;
  mcSig[1] = thickness;
  for (int i = 0; i < numgp; i++)
    mcSig[2 + i] = theMaterial[i]->getRho();
  if (const Matrix *Mc = massCache.lookup(mcSig, 2 + numgp, theNodes, numnodes, 2))
    return *Mc;

  K.Zero();

  static double rhoi[numgp];
  double sum = 0.0;
  for (int i = 0; i < numgp; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }
  if (sum == 0.0) {
    massCache.fill(K, mcSig, 2 + numgp, theNodes, numnodes, 2);   // Ladruno (ADR-77 G2 ext)
    return K;
  }

  // Ladruno (ADR 70 P3, deliberate divergence from SixNodeTri): HRZ lumping.
  // Upstream's plain N-lumping (Σ N_a ρ dV) gives EXACTLY ZERO corner masses
  // on the T6 with the 3-interior-point rule, and NEGATIVE corner masses on
  // distorted elements — unusable for explicit dynamics (M⁻¹ diagonal) and
  // the same failure the fork's H20 hit (ADR-72: row-sum corners = −M/8 →
  // HRZ-only there too). HRZ: diagonal of the consistent mass (∫ρ N_a² dV,
  // strictly positive) rescaled to preserve the total mass exactly.
  double d[numnodes];
  for (int a = 0; a < numnodes; a++) d[a] = 0.0;
  double total = 0.0;
  for (int i = 0; i < numgp; i++) {
    double rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);
    rhodvol *= (rhoi[i] * thickness * wts[i]);
    total += rhodvol;
    for (int a = 0; a < numnodes; a++)
      d[a] += shp[2][a] * shp[2][a] * rhodvol;
  }
  double dsum = 0.0;
  for (int a = 0; a < numnodes; a++) dsum += d[a];
  if (dsum <= 0.0) {
    massCache.fill(K, mcSig, 2 + numgp, theNodes, numnodes, 2);   // Ladruno (ADR-77 G2 ext)
    return K;
  }
  double scale = total / dsum;
  for (int a = 0, ia = 0; a < numnodes; a++, ia += 2) {
    double Nrho = scale * d[a];
    K(ia, ia)         = Nrho;
    K(ia + 1, ia + 1) = Nrho;
  }
  massCache.fill(K, mcSig, 2 + numgp, theNodes, numnodes, 2);   // Ladruno (ADR-77 G2 ext)
  return K;
}

void LadrunoLST::zeroLoad(void)
{
  Q.Zero();
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
}

int LadrunoLST::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    return 0;
  }
  opserr << "LadrunoLST::addLoad - load type unknown for ele " << this->getTag() << "\n";
  return -1;
}

int LadrunoLST::addInertiaLoadToUnbalance(const Vector &accel)
{
  static double rhoi[numgp];
  double sum = 0.0;
  for (int i = 0; i < numgp; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }
  if (sum == 0.0)
    return 0;

  static double ra[12];
  for (int a = 0; a < numnodes; a++) {
    const Vector &Raccel = theNodes[a]->getRV(accel);
    if (Raccel.Size() != 2) {
      opserr << "LadrunoLST::addInertiaLoadToUnbalance - incompatible sizes\n";
      return -1;
    }
    ra[2 * a]     = Raccel(0);
    ra[2 * a + 1] = Raccel(1);
  }
// Ladruno (ADR-77 G2 ext): consume the RETURNED matrix. The old idiom
  // called getMass() for its side effect of filling the class-static K and
  // then read K(i,i) directly -- with the per-instance cache a HIT returns
  // *Mi without touching K (which still holds the last TANGENT), so the
  // side-effect contract is dead. Caught by
  // test_dynamic_rayleigh_preserves_inertia[quad/lst].
  const Matrix &Mq = this->getMass();
  for (int i = 0; i < 2 * numnodes; i++)
    Q(i) += -Mq(i, i) * ra[i];
  return 0;
}

const Vector &LadrunoLST::getResistingForce(void)
{
  if (this->isFinite()) {          // -geom finite (ADR 70): spatial internal force
    this->formFinite(0);           // fills P with ∫ σ g dv − body force
    if (pressure != 0.0)
      P.addVector(1.0, pressureLoad, -1.0);
    P.addVector(1.0, Q, -1.0);
    return P;
  }

  P.Zero();
  static Vector sigma(3);

  // Ladruno (W2-E1, ADR-52): explicit bulk viscosity (2D). L_e is element-
  // constant -> hoisted. The bv coeffs are stripped at parse under -geom finite
  // (the finite branch returns early above), so this block only runs
  // small-strain. Off (b1=b2=0) skips the block -> default path bit-identical.
  const bool bvActive = (bulkVisc_b1 > 0.0 || bulkVisc_b2 > 0.0);
  const double bvLe = bvActive ? this->getCharacteristicLength() : 0.0;

  static Matrix B(3, 12);
  for (int i = 0; i < numgp; i++) {
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]) * thickness * wts[i];
    this->formB(B);
    sigma = theMaterial[i]->getStress();

    // Ladruno (W2-E1): viscous volumetric stress s=c_bulk*edotV added to the
    // normal comps (xx,yy) of the RESISTING FORCE only -- reported stress queries
    // the material (sigma is a local copy). c_d from the initial elastic tangent;
    // edotV = trace(D); dissipative (s*edotV >= 0).
    if (bvActive) {
      double edotV = 0.0;                       // volumetric strain rate, trace(D)
      for (int a = 0; a < numnodes; a++) {
        const Vector &va = theNodes[a]->getTrialVel();
        edotV += shp[0][a] * va(0) + shp[1][a] * va(1);
      }
      if (edotV != 0.0) {
        double rhoBV = theMaterial[i]->getRho();
        if (rhoBV > 0.0) {
          const Matrix &ddBV = theMaterial[i]->getInitialTangent();
          double cd = sqrt(fabs(ddBV(0, 0)) / rhoBV);       // elastic dilatational speed
          double cbulk = bulkVisc_b1 * rhoBV * cd * bvLe;   // linear: both signs
          if (edotV < 0.0)                                  // quadratic: compression only
            cbulk += bulkVisc_b2 * bulkVisc_b2 * rhoBV * bvLe * bvLe * (-edotV);
          double svisc = cbulk * edotV;
          sigma(0) += svisc;                                // xx
          sigma(1) += svisc;                                // yy
        }
      }
    }

    P.addMatrixTransposeVector(1.0, B, sigma, dvol);

    for (int a = 0, ia = 0; a < numnodes; a++, ia += 2) {
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

const Vector &LadrunoLST::getResistingForceIncInertia(void)
{
  static double rhoi[numgp];
  double sum = 0.0;
  for (int i = 0; i < numgp; i++) {
    rhoi[i] = (rho == 0.0) ? theMaterial[i]->getRho() : rho;
    sum += rhoi[i];
  }

  // SNAPSHOT the shared static P into a local BEFORE getRayleighDampingForces():
  // under -geom finite with stiffness-proportional Rayleigh (betaK) that call
  // re-enters getTangentStiff() -> formFinite(1) -> P.Zero()+refill, silently
  // destroying the f_int − Q + M·a just accumulated (adversarial gate, ADR 70
  // P4a review; LadrunoBrick donor pattern; LEDGER_quirks
  // "getResistingForceIncInertia MUST snapshot").  // Ladruno
  static Vector res(12);

  if (sum == 0.0) {
    this->getResistingForce();
    res = P;
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res += this->getRayleighDampingForces();
    P = res;
    return P;
  }

  static double a[12];
  for (int n = 0; n < numnodes; n++) {
    const Vector &accel = theNodes[n]->getTrialAccel();
    a[2 * n]     = accel(0);
    a[2 * n + 1] = accel(1);
  }
  this->getResistingForce();
  // Ladruno (ADR-77 G2 ext): consume the RETURNED matrix. The old idiom
  // called getMass() for its side effect of filling the class-static K and
  // then read K(i,i) directly -- with the per-instance cache a HIT returns
  // *Mi without touching K (which still holds the last TANGENT), so the
  // side-effect contract is dead. Caught by
  // test_dynamic_rayleigh_preserves_inertia[quad/lst].
  const Matrix &Mq = this->getMass();
  for (int i = 0; i < 2 * numnodes; i++)
    P(i) += Mq(i, i) * a[i];
  res = P;
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();
  P = res;
  return P;
}

// ======================================================================= //
//  -geom finite (ADR 70 P3): updated-Lagrangian large-strain assembly.    //
//  Mechanics are the shared LadrunoFiniteStrain2DKernel; this element     //
//  supplies the T6 reference gradients ∂Nₐ/∂X per GP (shapeFunction       //
//  differentiates w.r.t. the REFERENCE nodal coords). std ONLY — the      //
//  centroid-sampled F-bar was REFUTED on the T6 at P3 (conformal          //
//  zero-energy modes; see the header note / ADR 70 §9).                   //
// ======================================================================= //

// update — finite: per GP compute F, guard det F via the material
// (LogStrain2D::setTrialF rejects det F ≤ 0), and drive the material by the
// total F. Returns < 0 so the solver cuts the step.
int LadrunoLST::updateFinite(void)
{
  using namespace ladruno_fs2d;

  static double u[12];
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u[2 * a]     = d(0);
    u[2 * a + 1] = d(1);
  }

  static Matrix Fm(2, 2);
  for (int i = 0; i < numgp; i++) {
    this->shapeFunction(pts[i][0], pts[i][1]);       // shp = ∂Nₐ/∂X at GP i
    double gradRef[12], F[4];
    for (int a = 0; a < numnodes; a++) {
      gradRef[a] = shp[0][a];
      gradRef[numnodes + a] = shp[1][a];
    }
    deformationGradient2D(gradRef, u, numnodes, F);

    Fm(0, 0) = F[0]; Fm(0, 1) = F[1];
    Fm(1, 0) = F[2]; Fm(1, 1) = F[3];

    // dynamic_cast (not static) so a -geom finite element built on a non-finite
    // material fails gracefully here instead of casting UB (setDomain only warns).
    FiniteStrainND2DMaterial *fsm =
      dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[i]);
    if (fsm == 0) {
      opserr << "LadrunoLST::updateFinite - material at GP " << i
             << " is not a FiniteStrainND2DMaterial (element " << this->getTag()
             << ", -geom finite)\n";
      return -1;
    }
    // setTrialF guards det F ≤ 0 internally too (LogStrain2D), so the std path
    // (no explicit J guard above) still cuts the step on inversion.
    if (fsm->setTrialF(Fm) < 0) {
      opserr << "LadrunoLST::updateFinite - setTrialF failed at GP " << i
             << " (element " << this->getTag() << ", det F<=0?)\n";
      return -1;
    }
  }
  return 0;
}

// formFinite — fills the shared scratch P (residual) always and, when tang_flag
// == 1, K (consistent tangent). Recomputes F per GP for the spatial gradients /
// current volume, and consumes the material's Cauchy σ and full spatial modulus
// c set by updateFinite's setTrialF. No F-bar terms — see the block comment.
void LadrunoLST::formFinite(int tang_flag)
{
  using namespace ladruno_fs2d;

  double Kloc[144];
  double Ploc[12];
  for (int n = 0; n < 144; n++) Kloc[n] = 0.0;
  for (int n = 0; n < 12; n++)  Ploc[n] = 0.0;

  static double u[12];
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u[2 * a]     = d(0);
    u[2 * a + 1] = d(1);
  }

  const double bx = (applyLoad == 0) ? b[0] : appliedB[0];
  const double by = (applyLoad == 0) ? b[1] : appliedB[1];

  static Vector sigma(3);
  for (int i = 0; i < numgp; i++) {
    double detJ0 = this->shapeFunction(pts[i][0], pts[i][1]);   // reference detJ
    double gradRef[12], F[4], Finv[4], g[12];
    for (int a = 0; a < numnodes; a++) {
      gradRef[a] = shp[0][a];
      gradRef[numnodes + a] = shp[1][a];
    }
    double J = deformationGradient2D(gradRef, u, numnodes, F);
    inv2(F, Finv);
    spatialGradients2D(gradRef, Finv, numnodes, g);              // ∂Nₐ/∂x_j

    double dv = J * detJ0 * thickness * wts[i];                  // current volume

    sigma = theMaterial[i]->getStress();                        // Cauchy σ {11,22,12}
    double sig[2][2];
    sig[0][0] = sigma(0); sig[1][1] = sigma(1);
    sig[0][1] = sig[1][0] = sigma(2);

    // internal force f_{a,i} = ∫ σ_ij g_j^a dv
    addInternalForce2D(sig, g, numnodes, dv, Ploc);

    // body force (reference config, dead load per reference volume — matches
    // the reference-config lumped mass in getMass)
    double bw = detJ0 * thickness * wts[i];
    for (int a = 0, ia = 0; a < numnodes; a++, ia += 2) {
      Ploc[ia]     -= bw * shp[2][a] * bx;
      Ploc[ia + 1] -= bw * shp[2][a] * by;
    }

    if (tang_flag == 1) {
      FiniteStrainND2DMaterial *fsm =
        dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[i]);
      double c2[2][2][2][2];
      if (fsm == 0 || fsm->getSpatialTangentTensor2D(c2) != 0) {
        opserr << "LadrunoLST::formFinite - material at GP " << i
               << " is not a FiniteStrainND2DMaterial / gave no full 2D spatial "
                  "tangent (getSpatialTangentTensor2D); element " << this->getTag() << "\n";
        // K and P are SHARED static scratch (ADR 70 review): bailing out without
        // clearing them would hand the caller a stale sibling element's state.
        K.Zero();
        P.Zero();
        return;
      }
      double a4[2][2][2][2];
      spatialTangent2D(c2, sig, a4);                            // a = c − σ_il δ_jk
      addTangent2D(a4, g, numnodes, dv, Kloc, 12);
    }
  }

  // publish into the shared scratch. Kloc is ROW-major (kernel convention); the
  // OpenSees Matrix is column-major, so copy element-wise (mirrors the family
  // convention even though the std tangent happens to be symmetric).
  P.Zero();
  for (int r = 0; r < 12; r++) P(r) = Ploc[r];
  if (tang_flag == 1) {
    K.Zero();
    for (int r = 0; r < 12; r++)
      for (int cc = 0; cc < 12; cc++)
        K(r, cc) = Kloc[r * 12 + cc];
  }
}

int LadrunoLST::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(15);
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
  data(12) = bulkVisc_b1;   // Ladruno (W2-E1)
  data(13) = bulkVisc_b2;
  data(14) = static_cast<int>(geom);   // Ladruno (ADR 70)

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoLST::sendSelf - send Vector failed\n"; return res; }

  static ID idData(2 * numgp + numnodes);   // 3 mat class + 3 mat db + 6 nodes
  for (int i = 0; i < numgp; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    int matDbTag = theMaterial[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i + numgp) = matDbTag;
  }
  for (int i = 0; i < numnodes; i++)
    idData(2 * numgp + i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoLST::sendSelf - send ID failed\n"; return res; }

  for (int i = 0; i < numgp; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) { opserr << "WARNING LadrunoLST::sendSelf - send material failed\n"; return res; }
  }
  return res;
}

int LadrunoLST::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(15);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoLST::recvSelf - recv Vector failed\n"; return res; }

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
  bulkVisc_b1 = data(12);   // Ladruno (W2-E1)
  bulkVisc_b2 = data(13);
  geom        = static_cast<Geom>((int)data(14));   // Ladruno (ADR 70)

  static ID idData(2 * numgp + numnodes);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoLST::recvSelf - recv ID failed\n"; return res; }
  for (int i = 0; i < numnodes; i++)
    connectedExternalNodes(i) = idData(2 * numgp + i);

  if (theMaterial == 0) {
    theMaterial = new NDMaterial *[numgp];
    for (int i = 0; i < numgp; i++) {
      int matClassTag = idData(i);
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
        opserr << "LadrunoLST::recvSelf - broker could not create NDMaterial " << matClassTag << "\n";
        return -1;
      }
      theMaterial[i]->setDbTag(idData(i + numgp));
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) return res;
    }
  } else {
    for (int i = 0; i < numgp; i++) {
      int matClassTag = idData(i);
      if (theMaterial[i]->getClassTag() != matClassTag) {
        delete theMaterial[i];
        theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial[i] == 0) return -1;
      }
      theMaterial[i]->setDbTag(idData(i + numgp));
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) return res;
    }
  }
  // Ladruno (ADR-77 review wave): defensive, keeps the cache lifecycle uniform
  // across the G2 family -- nothing sig-exempt that recvSelf rewrites may
  // survive a re-receive into a live element.
  massCache.invalidate();
  return res;
}

void LadrunoLST::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoLST, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tformulation:  " << static_cast<int>(formulation) << "  type: " << typeString()
      << "  geom: " << (this->isFinite() ? "finite" : "linear") << "\n";
    s << "\tthickness:  " << thickness << "  rho: " << rho << "\n";
    if (bulkVisc_b1 > 0.0 || bulkVisc_b2 > 0.0)   // Ladruno (W2-E1)
      s << "\tbulk viscosity: b1=" << bulkVisc_b1 << " b2=" << bulkVisc_b2 << "\n";
    theMaterial[0]->Print(s, flag);
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoLST\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << ", " << connectedExternalNodes(2) << ", " << connectedExternalNodes(3)
      << ", " << connectedExternalNodes(4) << ", " << connectedExternalNodes(5) << "], ";
    s << "\"thickness\": " << thickness << ", ";
    s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int LadrunoLST::displaySelf(Renderer &theViewer, int displayMode, float fact,
                            const char **modes, int numMode)
{
  // draw the corner triangle (midside curvature not rendered)
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
    for (int i = 0; i < 3; i++)
      values(i) = theMaterial[i]->getStress()(displayMode - 1);
  else
    for (int i = 0; i < 3; i++) values(i) = 0.0;

  return theViewer.drawPolygon(coords, values, this->getTag());
}

Response *LadrunoLST::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  if (argc < 1) return 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoLST");
  output.attr("eleTag", this->getTag());

  if (LadrunoResp::is(argv[0], "force")) {
    theResponse = new ElementResponse(this, 1, P);
  } else if (LadrunoResp::is(argv[0], "material")) {
    if (argc > 1) {                        // guard before reading argv[1]
      int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= numgp)
        theResponse = theMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
    }
  } else if (LadrunoResp::is(argv[0], "stress")) {
    theResponse = new ElementResponse(this, 3, Vector(3 * numgp));
  } else if (LadrunoResp::is(argv[0], "stressPlaneStrain")) {
    // plane-strain stress incl. out-of-plane sigma_zz (NaN when the material
    // doesn't expose it); full GaussPoint/NdMaterialOutput tags so XML-driven
    // recorders (Ladruno/MPCO) get real component names instead of C1..C12
    for (int i = 0; i < numgp; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", pts[i][0]);
      output.attr("neta", pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma12");
      output.tag("ResponseType", "sigma33");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, 21, Vector(4 * numgp));
  } else if (LadrunoResp::is(argv[0], "strain")) {
    theResponse = new ElementResponse(this, 4, Vector(3 * numgp));
  } else if (LadrunoResp::is(argv[0], "charLength")) {
    theResponse = new ElementResponse(this, 5, 0.0);
  } else if (LadrunoResp::is(argv[0], "stiff")) {
    theResponse = new ElementResponse(this, 6, Matrix(P.Size(), P.Size()));
  } else if (LadrunoResp::is(argv[0], "stiffInitial")) {
    theResponse = new ElementResponse(this, 7, Matrix(P.Size(), P.Size()));
  }

  output.endTag();

  // Ladruno — base vocabulary (globalForce, dampingForce, dynamicForce,
  // inertialForce); Element::setResponse opens its own ElementOutput tag, so
  // this MUST come after endTag().
  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, output);
  return theResponse;
}

int LadrunoLST::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  if (responseID == 3 || responseID == 4) {
    static Vector v(3 * numgp);
    int cnt = 0;
    for (int i = 0; i < numgp; i++) {
      const Vector &s = (responseID == 3) ? theMaterial[i]->getStress()
                                          : theMaterial[i]->getStrain();
      v(cnt) = s(0); v(cnt + 1) = s(1); v(cnt + 2) = s(2);
      cnt += 3;
    }
    return eleInfo.setVector(v);
  }

  if (responseID == 21) {
    // 4-component plane-strain stress [sxx, syy, sxy, szz] per integration point
    static Vector v4(4 * numgp);
    int cnt = 0;
    for (int i = 0; i < numgp; i++) {
      const Vector &s = theMaterial[i]->getStress();
      v4(cnt) = s(0); v4(cnt + 1) = s(1); v4(cnt + 2) = s(2);
      v4(cnt + 3) = theMaterial[i]->getStressZZ();
      cnt += 4;
    }
    return eleInfo.setVector(v4);
  }

  if (responseID == 5)
    return eleInfo.setDouble(this->getCharacteristicLength());

  if (responseID == 6)
    return eleInfo.setMatrix(this->getTangentStiff());

  if (responseID == 7)
    return eleInfo.setMatrix(this->getInitialStiff());

  return this->Element::getResponse(responseID, eleInfo);
}

int LadrunoLST::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "pressure") == 0)
    return param.addObject(2, this);

  if (strstr(argv[0], "material") != 0) {
    if (argc < 3) return -1;
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= numgp)
      return theMaterial[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
    return -1;
  }

  for (int i = 0; i < numgp; i++) {
    int matRes = theMaterial[i]->setParameter(argv, argc, param);
    if (matRes != -1) res = matRes;
  }
  return res;
}

int LadrunoLST::updateParameter(int parameterID, Information &info)
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

// quadratic T6 shape functions and REFERENCE gradients — a faithful mirror of
// upstream SixNodeTri::shapeFunction (the reduce-to gate): corners at
// L1=s, L2=t, L3=1-s-t; midsides 4=(1-2), 5=(2-3), 6=(3-1).
double LadrunoLST::shapeFunction(double s, double t)
{
  const Vector &nd1 = theNodes[0]->getCrds();
  const Vector &nd2 = theNodes[1]->getCrds();
  const Vector &nd3 = theNodes[2]->getCrds();
  const Vector &nd4 = theNodes[3]->getCrds();
  const Vector &nd5 = theNodes[4]->getCrds();
  const Vector &nd6 = theNodes[5]->getCrds();

  shp[2][0] = s * (2 * s - 1);
  shp[2][1] = t * (2 * t - 1);
  shp[2][2] = (1 - s - t) * (1 - 2 * s - 2 * t);
  shp[2][3] = 4 * s * t;
  shp[2][4] = 4 * t * (1 - s - t);
  shp[2][5] = 4 * s * (1 - s - t);

  // natural derivatives (Na1 = dN/ds, Na2 = dN/dt)
  double N11 = 4 * s - 1;
  double N12 = 0;
  double N21 = 0;
  double N22 = 4 * t - 1;
  double N31 = -3 + 4 * s + 4 * t;
  double N32 = -3 + 4 * t + 4 * s;
  double N41 = 4 * t;
  double N42 = 4 * s;
  double N51 = -4 * t;
  double N52 = 4 - 4 * s - 8 * t;
  double N61 = 4 - 4 * t - 8 * s;
  double N62 = -4 * s;

  double J[2][2];
  J[0][0] = nd1(0) * N11 + nd2(0) * N21 + nd3(0) * N31 + nd4(0) * N41 +
            nd5(0) * N51 + nd6(0) * N61;
  J[0][1] = nd1(0) * N12 + nd2(0) * N22 + nd3(0) * N32 + nd4(0) * N42 +
            nd5(0) * N52 + nd6(0) * N62;
  J[1][0] = nd1(1) * N11 + nd2(1) * N21 + nd3(1) * N31 + nd4(1) * N41 +
            nd5(1) * N51 + nd6(1) * N61;
  J[1][1] = nd1(1) * N12 + nd2(1) * N22 + nd3(1) * N32 + nd4(1) * N42 +
            nd5(1) * N52 + nd6(1) * N62;

  double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  double oneOverdetJ = 1.0 / detJ;
  double L[2][2];
  L[0][0] =  J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] =  J[0][0] * oneOverdetJ;

  double L00 = L[0][0], L10 = L[1][0], L01 = L[0][1], L11 = L[1][1];

  shp[0][0] = L00 * N11 + L01 * N12;
  shp[0][1] = L00 * N21 + L01 * N22;
  shp[0][2] = L00 * N31 + L01 * N32;
  shp[0][3] = L00 * N41 + L01 * N42;
  shp[0][4] = L00 * N51 + L01 * N52;
  shp[0][5] = L00 * N61 + L01 * N62;

  shp[1][0] = L10 * N11 + L11 * N12;
  shp[1][1] = L10 * N21 + L11 * N22;
  shp[1][2] = L10 * N31 + L11 * N32;
  shp[1][3] = L10 * N41 + L11 * N42;
  shp[1][4] = L10 * N51 + L11 * N52;
  shp[1][5] = L10 * N61 + L11 * N62;

  return detJ;
}

// consistent edge pressure for the quadratic edge: corner 1/3, midside 2/3
// halves (upstream SixNodeTri distribution). NOTE (Ladruno): upstream's side-61
// segment uses dx61 = x4-x6 — a copy-paste typo (node 4 is the 1-2 midside, not
// an endpoint of edge 3-1); the correct closing segment is node6 -> node1.
// Written correctly here; documented as a deliberate divergence.
void LadrunoLST::setPressureLoadAtNodes(void)
{
  pressureLoad.Zero();
  if (pressure == 0.0)
    return;

  const Vector &node1 = theNodes[0]->getCrds();
  const Vector &node2 = theNodes[1]->getCrds();
  const Vector &node3 = theNodes[2]->getCrds();
  const Vector &node4 = theNodes[3]->getCrds();
  const Vector &node5 = theNodes[4]->getCrds();
  const Vector &node6 = theNodes[5]->getCrds();

  double x1 = node1(0), y1 = node1(1);
  double x2 = node2(0), y2 = node2(1);
  double x3 = node3(0), y3 = node3(1);
  double x4 = node4(0), y4 = node4(1);
  double x5 = node5(0), y5 = node5(1);
  double x6 = node6(0), y6 = node6(1);

  double dx14 = x4 - x1, dy14 = y4 - y1;
  double dx42 = x2 - x4, dy42 = y2 - y4;
  double dx25 = x5 - x2, dy25 = y5 - y2;
  double dx53 = x3 - x5, dy53 = y3 - y5;
  double dx36 = x6 - x3, dy36 = y6 - y3;
  double dx61 = x1 - x6, dy61 = y1 - y6;   // Ladruno: fixes upstream x4-x6 typo

  double fac1 = 0.3333333333333333;
  double fac2 = 0.6666666666666667;

  // edge 1-4-2 (corner-midside-corner), then 2-5-3, then 3-6-1
  pressureLoad(0)  += pressure * fac1 *  dy14;
  pressureLoad(6)  += pressure * fac2 *  dy14;
  pressureLoad(1)  += pressure * fac1 * -dx14;
  pressureLoad(7)  += pressure * fac2 * -dx14;

  pressureLoad(6)  += pressure * fac2 *  dy42;
  pressureLoad(2)  += pressure * fac1 *  dy42;
  pressureLoad(7)  += pressure * fac2 * -dx42;
  pressureLoad(3)  += pressure * fac1 * -dx42;

  pressureLoad(2)  += pressure * fac1 *  dy25;
  pressureLoad(8)  += pressure * fac2 *  dy25;
  pressureLoad(3)  += pressure * fac1 * -dx25;
  pressureLoad(9)  += pressure * fac2 * -dx25;

  pressureLoad(8)  += pressure * fac2 *  dy53;
  pressureLoad(4)  += pressure * fac1 *  dy53;
  pressureLoad(9)  += pressure * fac2 * -dx53;
  pressureLoad(5)  += pressure * fac1 * -dx53;

  pressureLoad(4)  += pressure * fac1 *  dy36;
  pressureLoad(10) += pressure * fac2 *  dy36;
  pressureLoad(5)  += pressure * fac1 * -dx36;
  pressureLoad(11) += pressure * fac2 * -dx36;

  pressureLoad(10) += pressure * fac2 *  dy61;
  pressureLoad(0)  += pressure * fac1 *  dy61;
  pressureLoad(11) += pressure * fac2 * -dx61;
  pressureLoad(1)  += pressure * fac1 * -dx61;
}
