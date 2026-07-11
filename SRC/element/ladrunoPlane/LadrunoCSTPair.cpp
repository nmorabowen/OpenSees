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

// LadrunoCSTPair — disjoint 2-triangle F-bar-Patch macro-element (ADR 70 P4a,
// dSNPO §15.1.9). See LadrunoCSTPair.h for the formulation summary and
// tests/cstpair_reference.py for the numpy oracle this assembly must match.

#include <LadrunoCSTPair.h>
#include <LadrunoFiniteStrain2DKernel.h>   // shared 2D finite-strain kernel (ADR 70)
#include <FiniteStrainND2DMaterial.h>      // setTrialF(F) seam
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <math.h>

double LadrunoCSTPair::matrixData[LadrunoCSTPair::ndf * LadrunoCSTPair::ndf];
Matrix LadrunoCSTPair::K(matrixData, LadrunoCSTPair::ndf, LadrunoCSTPair::ndf);
Vector LadrunoCSTPair::P(LadrunoCSTPair::ndf);

// triangle connectivity within the macro: T1=(n1,n2,n3), T2=(n1,n3,n4)
static const int TRI_NODES[2][3] = { {0, 1, 2}, {0, 2, 3} };

LadrunoCSTPair::LadrunoCSTPair(int tag, int nd1, int nd2, int nd3, int nd4,
                               NDMaterial &m, double t,
                               double r, double b1, double b2)
 :Element(tag, ELE_TAG_LadrunoCSTPair),
  theMaterial(0), connectedExternalNodes(numnodes),
  Q(ndf), thickness(t), rho(r), Ki(0)
{
  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b1;
  b[1] = b2;

  theMaterial = new NDMaterial *[numtri];
  for (int i = 0; i < numtri; i++) {
    theMaterial[i] = m.getCopy("PlaneStrain");   // finite path is PlaneStrain-only
    if (theMaterial[i] == 0) {
      opserr << "LadrunoCSTPair::LadrunoCSTPair -- failed to get a copy of material model\n";
      exit(-1);
    }
  }

  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = nd3;
  connectedExternalNodes(3) = nd4;

  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoCSTPair::LadrunoCSTPair()
 :Element(0, ELE_TAG_LadrunoCSTPair),
  theMaterial(0), connectedExternalNodes(numnodes),
  Q(ndf), thickness(0.0), rho(0.0), Ki(0)
{
  applyLoad = 0;
  appliedB[0] = appliedB[1] = 0.0;
  b[0] = b[1] = 0.0;
  for (int i = 0; i < numnodes; i++)
    theNodes[i] = 0;
}

LadrunoCSTPair::~LadrunoCSTPair()
{
  if (theMaterial) {
    for (int i = 0; i < numtri; i++)
      if (theMaterial[i]) delete theMaterial[i];
    delete[] theMaterial;
  }
  if (Ki) delete Ki;
}

int LadrunoCSTPair::getNumExternalNodes(void) const { return numnodes; }
const ID &LadrunoCSTPair::getExternalNodes(void) { return connectedExternalNodes; }
Node **LadrunoCSTPair::getNodePtrs(void) { return theNodes; }
int LadrunoCSTPair::getNumDOF(void) { return ndf; }

void LadrunoCSTPair::setDomain(Domain *theDomain)
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
      opserr << "LadrunoCSTPair::setDomain -- node " << connectedExternalNodes(i)
             << " must have 2 DOFs\n";
      return;
    }
  this->DomainComponent::setDomain(theDomain);

  // the factory already refuses a non-FiniteStrainND2DMaterial; this re-check
  // covers broker-built (recvSelf) elements — fail LOUD at setup, and
  // updatePair's dynamic_cast is the UB backstop.
  for (int i = 0; i < numtri; i++)
    if (dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[i]) == 0) {
      opserr << "LadrunoCSTPair::setDomain -- WARNING element " << this->getTag()
             << " material is not a FiniteStrainND2DMaterial (e.g. LogStrain2D); "
                "the analysis will fail\n";
      break;
    }

  // both triangles must be CCW / non-degenerate in the reference config
  double gradRef[numtri][2 * numnodes], area[numtri];
  if (this->refGradients(gradRef, area) <= 0.0)
    opserr << "LadrunoCSTPair::setDomain -- WARNING element " << this->getTag()
           << " has a non-positive-area triangle (node order n1..n4 must be a "
              "CCW quad split along the n1-n3 diagonal)\n";
}

double LadrunoCSTPair::getCharacteristicLength(void)
{
  // crack-band size: the constituents are triangles, so use the mean-triangle
  // convention sqrt(2*A_mean) (LadrunoCST uses sqrt(2*A) per triangle).
  double gradRef[numtri][2 * numnodes], area[numtri];
  if (this->refGradients(gradRef, area) <= 0.0)
    return this->Element::getCharacteristicLength();
  return sqrt(area[0] + area[1]);   // = sqrt(2 * (A1+A2)/2)
}

int LadrunoCSTPair::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0)
    opserr << "LadrunoCSTPair::commitState() - failed in base class\n";
  for (int i = 0; i < numtri; i++)
    retVal += theMaterial[i]->commitState();
  return retVal;
}

int LadrunoCSTPair::revertToLastCommit(void)
{
  int retVal = 0;
  for (int i = 0; i < numtri; i++)
    retVal += theMaterial[i]->revertToLastCommit();
  return retVal;
}

int LadrunoCSTPair::revertToStart(void)
{
  int retVal = 0;
  for (int i = 0; i < numtri; i++)
    retVal += theMaterial[i]->revertToStart();
  return retVal;
}

// per-triangle constant reference gradients in 4-node PADDED layout
// (gradRef[t][k*4+a] = dN_a/dX_k, zero on the node the triangle doesn't touch)
// + areas. Returns min(2*A_t) so a caller can reject degenerate geometry.
double LadrunoCSTPair::refGradients(double gradRef[numtri][2 * numnodes],
                                    double area[numtri]) const
{
  double minDet = 1.0e300;
  for (int t = 0; t < numtri; t++) {
    const Vector &X1 = theNodes[TRI_NODES[t][0]]->getCrds();
    const Vector &X2 = theNodes[TRI_NODES[t][1]]->getCrds();
    const Vector &X3 = theNodes[TRI_NODES[t][2]]->getCrds();

    // detJ = 2*A (CCW positive); dN_a/dX from the standard T3 closed form
    double detJ = (X2(0) - X1(0)) * (X3(1) - X1(1))
                - (X3(0) - X1(0)) * (X2(1) - X1(1));
    area[t] = 0.5 * detJ;
    if (detJ < minDet) minDet = detJ;
    if (detJ == 0.0) continue;                 // leave rows zero; caller rejects

    double id = 1.0 / detJ;
    const Vector *Xn[3] = { &X1, &X2, &X3 };
    for (int k = 0; k < 2 * numnodes; k++) gradRef[t][k] = 0.0;
    for (int a = 0; a < 3; a++) {
      int jn = (a + 1) % 3, kn = (a + 2) % 3;
      double gx = ((*Xn[jn])(1) - (*Xn[kn])(1)) * id;   // dN_a/dX
      double gy = ((*Xn[kn])(0) - (*Xn[jn])(0)) * id;   // dN_a/dY
      int node = TRI_NODES[t][a];
      gradRef[t][0 * numnodes + node] = gx;
      gradRef[t][1 * numnodes + node] = gy;
    }
  }
  return minDet;
}

int LadrunoCSTPair::update(void)
{
  return this->updatePair();
}

// update — compute both triangles' F (constant each), the patch dilatation
// J̄ = (J1 V1 + J2 V2)/(V1+V2), and drive each material at F̄_e = (J̄/J_e)^{1/2} F_e.
// Returns < 0 (step cut) on det F_e ≤ 0 — checked HERE because the F-bar scale
// needs both J's before any setTrialF, and fbarScale2D on a negative ratio is NaN.
int LadrunoCSTPair::updatePair(void)
{
  using namespace ladruno_fs2d;

  double gradRef[numtri][2 * numnodes], area[numtri];
  if (this->refGradients(gradRef, area) <= 0.0) {
    opserr << "LadrunoCSTPair::updatePair - degenerate reference triangle (element "
           << this->getTag() << ")\n";
    return -1;
  }

  static double u[ndf];
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u[2 * a]     = d(0);
    u[2 * a + 1] = d(1);
  }

  double F[numtri][4], J[numtri];
  for (int t = 0; t < numtri; t++) {
    J[t] = deformationGradient2D(gradRef[t], u, numnodes, F[t]);
    if (J[t] <= 0.0) {
      opserr << "LadrunoCSTPair::updatePair - det F <= 0 in triangle " << t + 1
             << " (element " << this->getTag() << "), cutting step\n";
      return -1;
    }
  }

  const double Vp = area[0] + area[1];
  const double Jbar = (J[0] * area[0] + J[1] * area[1]) / Vp;   // eq 15.36

  static Matrix Fm(2, 2);
  for (int t = 0; t < numtri; t++) {
    const double s = fbarScale2D(Jbar, J[t]);                   // (J̄/J)^{1/2}
    Fm(0, 0) = s * F[t][0]; Fm(0, 1) = s * F[t][1];
    Fm(1, 0) = s * F[t][2]; Fm(1, 1) = s * F[t][3];

    FiniteStrainND2DMaterial *fsm =
      dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[t]);
    if (fsm == 0) {
      opserr << "LadrunoCSTPair::updatePair - material is not a "
                "FiniteStrainND2DMaterial (element " << this->getTag() << ")\n";
      return -1;
    }
    if (fsm->setTrialF(Fm) < 0) {
      opserr << "LadrunoCSTPair::updatePair - setTrialF failed (element "
             << this->getTag() << ", triangle " << t + 1 << ")\n";
      return -1;
    }
  }
  return 0;
}

// formPair — fills the shared scratch P (residual) and, when tang_flag == 1,
// K (exact consistent tangent). Per triangle e (dSNPO Remark 15.2: residual =
// standard form at σ̄ on the UNBARRED spatial config):
//   f_e = ∫ σ̄ g_e dv_e,  dv_e = J_e A_e t
//   K  += ∫ g_eᵀ (c̄ − σ̄δ) g_e dv_e  +  ∫ p_e (ḡ − g_e) dv_e     [kernel calls]
// with ḡ = Σ_s (v_s/v_p) g_s the patch volume-weighted row — this ONE
// substitute-row call reproduces both the (v_e/v_p − 1) own block of eq 15.37
// and the (v_s/v_p) cross block of eq 15.38. All rows live in the 4-node
// padded layout, so the kernel needs no changes. GENERALLY UNSYMMETRIC.
void LadrunoCSTPair::formPair(int tang_flag)
{
  using namespace ladruno_fs2d;

  double Kloc[ndf * ndf];
  double Ploc[ndf];
  for (int n = 0; n < ndf * ndf; n++) Kloc[n] = 0.0;
  for (int n = 0; n < ndf; n++)       Ploc[n] = 0.0;

  double gradRef[numtri][2 * numnodes], area[numtri];
  if (this->refGradients(gradRef, area) <= 0.0) {
    K.Zero(); P.Zero();
    return;
  }

  static double u[ndf];
  for (int a = 0; a < numnodes; a++) {
    const Vector &d = theNodes[a]->getTrialDisp();
    u[2 * a]     = d(0);
    u[2 * a + 1] = d(1);
  }

  // pass 1: per-triangle spatial rows + deformed volumes (needed for ḡ)
  double F[4], Finv[4];
  double g[numtri][2 * numnodes], dv[numtri];
  for (int t = 0; t < numtri; t++) {
    double J = deformationGradient2D(gradRef[t], u, numnodes, F);
    inv2(F, Finv);
    spatialGradients2D(gradRef[t], Finv, numnodes, g[t]);   // padded rows
    dv[t] = J * area[t] * thickness;
  }
  const double vp = dv[0] + dv[1];
  double gbar[2 * numnodes];
  for (int k = 0; k < 2 * numnodes; k++)
    gbar[k] = (dv[0] * g[0][k] + dv[1] * g[1][k]) / vp;

  const double bx = (applyLoad == 0) ? b[0] : appliedB[0];
  const double by = (applyLoad == 0) ? b[1] : appliedB[1];

  // pass 2: material state (σ̄ set by updatePair), residual + tangent
  for (int t = 0; t < numtri; t++) {
    static Vector sigma(3);
    sigma = theMaterial[t]->getStress();          // Cauchy σ̄ {11,22,12}
    double sig[2][2];
    sig[0][0] = sigma(0); sig[1][1] = sigma(1);
    sig[0][1] = sig[1][0] = sigma(2);

    addInternalForce2D(sig, g[t], numnodes, dv[t], Ploc);

    // body force: reference-config dead load, centroid N = 1/3 per tri node
    double bw = area[t] * thickness / 3.0;
    for (int a = 0; a < 3; a++) {
      int ia = 2 * TRI_NODES[t][a];
      Ploc[ia]     -= bw * bx;
      Ploc[ia + 1] -= bw * by;
    }

    if (tang_flag == 1) {
      FiniteStrainND2DMaterial *fsm =
        dynamic_cast<FiniteStrainND2DMaterial *>(theMaterial[t]);
      double c2[2][2][2][2];
      if (fsm == 0 || fsm->getSpatialTangentTensor2D(c2) != 0) {
        opserr << "LadrunoCSTPair::formPair - material gave no full 2D spatial "
                  "tangent (element " << this->getTag() << ")\n";
        // K and P are SHARED static scratch: never hand back a sibling's state
        K.Zero();
        P.Zero();
        return;
      }
      double a4[2][2][2][2];
      spatialTangent2D(c2, sig, a4);                          // a = c − σ_il δ_jk
      addTangent2D(a4, g[t], numnodes, dv[t], Kloc, ndf);
      // patch coupling: substitute row ḡ in place of the centroid g₀
      addFbarCoupling2D(a4, sig, g[t], gbar, numnodes, dv[t], Kloc, ndf);
    }
  }

  // publish. Kloc is ROW-major (kernel convention); the OpenSees Matrix is
  // column-major — copy ELEMENT-WISE. A blind memcpy would TRANSPOSE the
  // genuinely-unsymmetric patch tangent (the P1 quad lesson, ADR 70).
  P.Zero();
  for (int r = 0; r < ndf; r++) P(r) = Ploc[r];
  if (tang_flag == 1) {
    K.Zero();
    for (int r = 0; r < ndf; r++)
      for (int c = 0; c < ndf; c++)
        K(r, c) = Kloc[r * ndf + c];
  }
}

const Matrix &LadrunoCSTPair::getTangentStiff(void)
{
  this->formPair(1);
  return K;
}

const Matrix &LadrunoCSTPair::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;

  // deliberately the SYMMETRIC small-strain reference tangent Bᵀ D₀ B dV per
  // triangle (family convention: a well-posed symmetric seed for Rayleigh βK0 /
  // -initial / eigen). NOT the consistent patch tangent — that one is
  // unsymmetric whenever σ ≠ 0 (the coupling is stress-proportional and DOES
  // vanish at F = I, unlike the quad's centroid coupling; still, keep the seed
  // symmetric and reference-config like the rest of the family).
  double gradRef[numtri][2 * numnodes], area[numtri];
  K.Zero();
  if (this->refGradients(gradRef, area) > 0.0) {
    for (int t = 0; t < numtri; t++) {
      const Matrix &D = theMaterial[t]->getInitialTangent();
      double dvol = area[t] * thickness;
      const double *gx = &gradRef[t][0];
      const double *gy = &gradRef[t][numnodes];
      for (int a = 0; a < numnodes; a++) {
        for (int bb = 0; bb < numnodes; bb++) {
          // B-matrix contraction with engineering-shear D (3x3 Voigt)
          double DB00 = D(0,0)*gx[bb] + D(0,2)*gy[bb];
          double DB10 = D(1,0)*gx[bb] + D(1,2)*gy[bb];
          double DB20 = D(2,0)*gx[bb] + D(2,2)*gy[bb];
          double DB01 = D(0,1)*gy[bb] + D(0,2)*gx[bb];
          double DB11 = D(1,1)*gy[bb] + D(1,2)*gx[bb];
          double DB21 = D(2,1)*gy[bb] + D(2,2)*gx[bb];
          K(2*a,   2*bb)   += dvol * (gx[a]*DB00 + gy[a]*DB20);
          K(2*a,   2*bb+1) += dvol * (gx[a]*DB01 + gy[a]*DB21);
          K(2*a+1, 2*bb)   += dvol * (gy[a]*DB10 + gx[a]*DB20);
          K(2*a+1, 2*bb+1) += dvol * (gy[a]*DB11 + gx[a]*DB21);
        }
      }
    }
  }
  Ki = new Matrix(K);
  return *Ki;   // cached copy, not the shared static scratch
}

const Matrix &LadrunoCSTPair::getMass(void)
{
  K.Zero();
  double gradRef[numtri][2 * numnodes], area[numtri];
  if (this->refGradients(gradRef, area) <= 0.0)
    return K;

  for (int t = 0; t < numtri; t++) {
    double r = (rho == 0.0) ? theMaterial[t]->getRho() : rho;
    if (r == 0.0) continue;
    // lumped: each triangle's mass /3 to its own nodes (LadrunoCST convention)
    double mNode = r * area[t] * thickness / 3.0;
    for (int a = 0; a < 3; a++) {
      int ia = 2 * TRI_NODES[t][a];
      K(ia, ia)         += mNode;
      K(ia + 1, ia + 1) += mNode;
    }
  }
  return K;
}

void LadrunoCSTPair::zeroLoad(void)
{
  Q.Zero();
  applyLoad = 0;
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
}

int LadrunoCSTPair::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    return 0;
  }
  opserr << "LadrunoCSTPair::addLoad - load type unknown for ele " << this->getTag() << "\n";
  return -1;
}

int LadrunoCSTPair::addInertiaLoadToUnbalance(const Vector &accel)
{
  bool haveRho = (rho != 0.0);
  for (int t = 0; t < numtri && !haveRho; t++)
    haveRho = (theMaterial[t]->getRho() != 0.0);
  if (!haveRho)
    return 0;

  static double ra[ndf];
  for (int a = 0; a < numnodes; a++) {
    const Vector &Raccel = theNodes[a]->getRV(accel);
    if (Raccel.Size() != 2) {
      opserr << "LadrunoCSTPair::addInertiaLoadToUnbalance - incompatible sizes\n";
      return -1;
    }
    ra[2 * a]     = Raccel(0);
    ra[2 * a + 1] = Raccel(1);
  }
  this->getMass();
  for (int i = 0; i < ndf; i++)
    Q(i) += -K(i, i) * ra[i];
  return 0;
}

const Vector &LadrunoCSTPair::getResistingForce(void)
{
  this->formPair(0);           // fills P with Σ_e ∫ σ̄ g_e dv_e − body force
  P.addVector(1.0, Q, -1.0);
  return P;
}

const Vector &LadrunoCSTPair::getResistingForceIncInertia(void)
{
  bool haveRho = (rho != 0.0);
  for (int t = 0; t < numtri && !haveRho; t++)
    haveRho = (theMaterial[t]->getRho() != 0.0);

  if (!haveRho) {
    this->getResistingForce();
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();
    return P;
  }

  static double a[ndf];
  for (int n = 0; n < numnodes; n++) {
    const Vector &accel = theNodes[n]->getTrialAccel();
    a[2 * n]     = accel(0);
    a[2 * n + 1] = accel(1);
  }
  this->getResistingForce();
  // P is the shared static scratch getMass() also writes through K — snapshot
  // pattern not needed here because getMass writes K, not P; but Rayleigh
  // MUST come after the inertia add (getRayleighDampingForces reuses P? no —
  // it returns its own vector). Keep the LadrunoCST ordering exactly.
  this->getMass();
  for (int i = 0; i < ndf; i++)
    P(i) += K(i, i) * a[i];
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();
  return P;
}

int LadrunoCSTPair::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(9);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = rho;
  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoCSTPair::sendSelf - send Vector failed\n"; return res; }

  static ID idData(2 * numtri + numnodes);
  for (int t = 0; t < numtri; t++) {
    idData(t) = theMaterial[t]->getClassTag();
    int matDbTag = theMaterial[t]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0) theMaterial[t]->setDbTag(matDbTag);
    }
    idData(numtri + t) = matDbTag;
  }
  for (int i = 0; i < numnodes; i++)
    idData(2 * numtri + i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoCSTPair::sendSelf - send ID failed\n"; return res; }

  for (int t = 0; t < numtri; t++) {
    res += theMaterial[t]->sendSelf(commitTag, theChannel);
    if (res < 0) { opserr << "WARNING LadrunoCSTPair::sendSelf - send material failed\n"; return res; }
  }
  return res;
}

int LadrunoCSTPair::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) { opserr << "WARNING LadrunoCSTPair::recvSelf - recv Vector failed\n"; return res; }

  this->setTag((int)data(0));
  thickness = data(1);
  b[0]      = data(2);
  b[1]      = data(3);
  rho       = data(4);
  alphaM    = data(5);
  betaK     = data(6);
  betaK0    = data(7);
  betaKc    = data(8);

  static ID idData(2 * numtri + numnodes);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) { opserr << "WARNING LadrunoCSTPair::recvSelf - recv ID failed\n"; return res; }
  for (int i = 0; i < numnodes; i++)
    connectedExternalNodes(i) = idData(2 * numtri + i);

  if (theMaterial == 0) {
    theMaterial = new NDMaterial *[numtri];
    for (int t = 0; t < numtri; t++) theMaterial[t] = 0;
  }
  for (int t = 0; t < numtri; t++) {
    int matClassTag = idData(t);
    if (theMaterial[t] == 0 || theMaterial[t]->getClassTag() != matClassTag) {
      if (theMaterial[t]) delete theMaterial[t];
      theMaterial[t] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[t] == 0) {
        opserr << "LadrunoCSTPair::recvSelf - broker could not create NDMaterial "
               << matClassTag << "\n";
        return -1;
      }
    }
    theMaterial[t]->setDbTag(idData(numtri + t));
    res += theMaterial[t]->recvSelf(commitTag, theChannel, theBroker);
  }
  return res;
}

void LadrunoCSTPair::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLadrunoCSTPair (F-bar-Patch, dSNPO 15.1.9), element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\ttype: PlaneStrain  thickness:  " << thickness << "  rho: " << rho
      << "  geom: finite (patch)\n";
    theMaterial[0]->Print(s, flag);
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoCSTPair\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << ", " << connectedExternalNodes(2) << ", " << connectedExternalNodes(3) << "], ";
    s << "\"thickness\": " << thickness << ", ";
    s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int LadrunoCSTPair::displaySelf(Renderer &theViewer, int displayMode, float fact,
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
  for (int i = 0; i < numnodes; i++) values(i) = 0.0;
  return theViewer.drawPolygon(coords, values, this->getTag());
}

Response *LadrunoCSTPair::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoCSTPair");
  output.attr("eleTag", this->getTag());

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
    theResponse = new ElementResponse(this, 1, P);
  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum >= 1 && pointNum <= numtri)
      theResponse = theMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
  } else if (strcmp(argv[0], "stresses") == 0 || strcmp(argv[0], "stress") == 0) {
    theResponse = new ElementResponse(this, 3, Vector(3 * numtri));
  } else if (strcmp(argv[0], "Jbar") == 0 || strcmp(argv[0], "jbar") == 0) {
    // the shared patch dilatation det(F-bar) = v_patch/V_patch — the quantity
    // the whole element exists to control; prime diagnostic for pressure checks
    theResponse = new ElementResponse(this, 6, 0.0);
  } else if (strcmp(argv[0], "charLength") == 0 || strcmp(argv[0], "characteristicLength") == 0) {
    theResponse = new ElementResponse(this, 5, 0.0);
  }

  output.endTag();
  return theResponse;
}

int LadrunoCSTPair::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  if (responseID == 3) {
    static Vector s6(3 * numtri);
    for (int t = 0; t < numtri; t++) {
      const Vector &s = theMaterial[t]->getStress();
      s6(3 * t) = s(0); s6(3 * t + 1) = s(1); s6(3 * t + 2) = s(2);
    }
    return eleInfo.setVector(s6);
  }
  if (responseID == 5)
    return eleInfo.setDouble(this->getCharacteristicLength());
  if (responseID == 6) {
    using namespace ladruno_fs2d;
    double gradRef[numtri][2 * numnodes], area[numtri];
    if (this->refGradients(gradRef, area) <= 0.0)
      return eleInfo.setDouble(0.0);
    static double u[ndf];
    for (int a = 0; a < numnodes; a++) {
      const Vector &d = theNodes[a]->getTrialDisp();
      u[2 * a] = d(0); u[2 * a + 1] = d(1);
    }
    double F[4];
    double J0 = deformationGradient2D(gradRef[0], u, numnodes, F);
    double J1 = deformationGradient2D(gradRef[1], u, numnodes, F);
    return eleInfo.setDouble((J0 * area[0] + J1 * area[1]) / (area[0] + area[1]));
  }
  return -1;
}
