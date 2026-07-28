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

// LadrunoSolidShell — 8-node ANS/EAS solid-shell element (Ladruno fork, ADR 66).
// See LadrunoSolidShell.h for the formulation summary and the file-level docs.
//
// Formulation references:
//   Dvorkin & Bathe (1984)  — ANS/MITC transverse-shear tying
//   Betsch & Stein (1995)   — ANS thickness-strain (E33) corner tying
//   Simo & Rifai (1990)     — EAS; the state-dependent inner-Newton +
//                             condensation machinery is adapted from
//                             LadrunoBrick::formEAStrue (ADR 66 §3)
//   Hauptmann & Schweizerhof (1998), Klinkel-Gruttmann-Wagner (1999) —
//                             the assembled ANS+EAS solid-shell.  // Ladruno

#include "LadrunoSolidShell.h"

#include <classTags.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <LadrunoResponseTokens.h>   // Ladruno — shared recorder-token aliases
#include <Information.h>
#include <Parameter.h>
#include <ElementalLoad.h>
#include <OPS_Globals.h>

#include <string.h>

// -------- static workspace --------
Matrix LadrunoSolidShell::stiff(24, 24);
Vector LadrunoSolidShell::resid(24);
Matrix LadrunoSolidShell::mass(24, 24);

// out-of-line definitions (the in-class const ints are odr-used)
const int LadrunoSolidShell::NEAS;
const int LadrunoSolidShell::NZ_MIN;
const int LadrunoSolidShell::NZ_MAX;

// node natural coordinates (brick order: 1-4 bottom zeta=-1 CCW, 5-8 top)
const double LadrunoSolidShell::natCoord[8][3] = {
  {-1.0, -1.0, -1.0}, { 1.0, -1.0, -1.0}, { 1.0,  1.0, -1.0}, {-1.0,  1.0, -1.0},
  {-1.0, -1.0,  1.0}, { 1.0, -1.0,  1.0}, { 1.0,  1.0,  1.0}, {-1.0,  1.0,  1.0}
};

// Voigt pair table {xx,yy,zz,xy,yz,zx} (0-based tensor indices)
static const int voigtPair[6][2] = {
  {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}
};

//----------------------------------------------------------------------
// constructors / destructor
//----------------------------------------------------------------------

LadrunoSolidShell::LadrunoSolidShell()
  : Element(0, ELE_TAG_LadrunoSolidShell),
    connectedExternalNodes(8),
    materialPointers(0),
    formulation(Formulation::ANS), quadz(QuadZ::GAUSS),
    nz(0), numGP(0), massType(0),
    alpha(NEAS), alphaCommit(NEAS),
    j0det(0.0), easBuilt(false),
    load(0), Ki(0)
{
  alpha.Zero();
  alphaCommit.Zero();
  for (int i = 0; i < 8; i++) {
    nodePointers[i] = 0;
    for (int k = 0; k < 3; k++) X[i][k] = 0.0;
  }
  for (int r = 0; r < 6; r++) t0col3[r] = 0.0;
}

LadrunoSolidShell::LadrunoSolidShell(int tag,
                                     int node1, int node2, int node3, int node4,
                                     int node5, int node6, int node7, int node8,
                                     NDMaterial &theMaterial,
                                     int nzPts, QuadZ zQuad, Formulation form,
                                     int matype)
  : Element(tag, ELE_TAG_LadrunoSolidShell),
    connectedExternalNodes(8),
    materialPointers(0),
    formulation(form), quadz(zQuad),
    nz(nzPts), numGP(4 * nzPts), massType(matype),
    alpha(NEAS), alphaCommit(NEAS),
    j0det(0.0), easBuilt(false),
    load(0), Ki(0)
{
  alpha.Zero();
  alphaCommit.Zero();

  // The OPS factory validates nz against the (rule, nz) tables; this guards a
  // DIRECT C++ construction (developer path) — nz outside [NZ_MIN, NZ_MAX]
  // makes numGP overflow the fixed 71-slot sendSelf ID. Warn loudly; the
  // element is still built (materialPointers is sized to numGP) but must not be
  // serialized.  // Ladruno
  if (nzPts < NZ_MIN || nzPts > NZ_MAX)
    opserr << "WARNING LadrunoSolidShell " << tag << ": nz=" << nzPts
           << " is out of range [" << NZ_MIN << ", " << NZ_MAX
           << "] — construct via the element() factory (this element is not "
              "serialization-safe)\n";

  connectedExternalNodes(0) = node1;
  connectedExternalNodes(1) = node2;
  connectedExternalNodes(2) = node3;
  connectedExternalNodes(3) = node4;
  connectedExternalNodes(4) = node5;
  connectedExternalNodes(5) = node6;
  connectedExternalNodes(6) = node7;
  connectedExternalNodes(7) = node8;

  materialPointers = new NDMaterial *[numGP];
  for (int i = 0; i < numGP; i++) {
    materialPointers[i] = theMaterial.getCopy("ThreeDimensional");
    if (materialPointers[i] == 0)
      opserr << "LadrunoSolidShell::constructor - failed to get a "
                "ThreeDimensional material copy\n";
  }

  for (int i = 0; i < 8; i++) {
    nodePointers[i] = 0;
    for (int k = 0; k < 3; k++) X[i][k] = 0.0;
  }
  for (int r = 0; r < 6; r++) t0col3[r] = 0.0;
}

LadrunoSolidShell::~LadrunoSolidShell()
{
  if (materialPointers != 0) {
    for (int i = 0; i < numGP; i++)
      if (materialPointers[i] != 0)
        delete materialPointers[i];
    delete[] materialPointers;
  }
  if (load != 0) delete load;
  if (Ki != 0) delete Ki;
}

//----------------------------------------------------------------------
// domain
//----------------------------------------------------------------------

void LadrunoSolidShell::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    for (int i = 0; i < 8; i++) nodePointers[i] = 0;
    this->DomainComponent::setDomain(0);
    return;
  }

  bool ok = true;
  for (int i = 0; i < 8; i++) {
    nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));
    if (nodePointers[i] == 0) {
      opserr << "LadrunoSolidShell::setDomain - node "
             << connectedExternalNodes(i) << " does not exist\n";
      ok = false;
      continue;
    }
    if (nodePointers[i]->getNumberDOF() != 3) {
      opserr << "LadrunoSolidShell::setDomain - node "
             << connectedExternalNodes(i) << " has "
             << nodePointers[i]->getNumberDOF() << " DOF (needs 3)\n";
      ok = false;
    }
  }

  if (ok) {
    for (int i = 0; i < 8; i++) {
      const Vector &crd = nodePointers[i]->getCrds();
      for (int k = 0; k < 3; k++) X[i][k] = crd(k);
    }
    // one-time geometry sanity: a non-positive center Jacobian means the node
    // ordering is mis-wound (top/bottom faces swapped) or the element is
    // degenerate/concave. Warn ONCE here (cheap, out of the hot loop) so a
    // silently-wrong-stiffness run is caught up front — for BOTH formulations
    // (buildEAS below only runs for ANS). Per-GP detJ can still go negative on
    // a pathologically concave hex; that stays user error in v1.  // Ladruno
    double N[8], dN[8][3], J[3][3];
    const double xi0[3] = {0.0, 0.0, 0.0};
    shapeFcn(xi0, N, dN);
    const double j0 = jacobian(X, dN, J);
    if (j0 <= 0.0)
      opserr << "WARNING LadrunoSolidShell::setDomain - element "
             << this->getTag() << ": non-positive center Jacobian (" << j0
             << "); check node ordering (1-4 bottom face CCW, 5-8 top)\n";

    // center enhanced-mode map (Simo-Rifai). Geometry-only, so the receive
    // side rebuilds it here after recvSelf — only alphaCommit travels.
    easBuilt = false;
    if (formulation == Formulation::ANS)
      buildEAS();
  }

  this->DomainComponent::setDomain(theDomain);
}

//----------------------------------------------------------------------
// topology
//----------------------------------------------------------------------

int LadrunoSolidShell::getNumExternalNodes(void) const { return 8; }

const ID &LadrunoSolidShell::getExternalNodes(void) { return connectedExternalNodes; }

Node **LadrunoSolidShell::getNodePtrs(void) { return nodePointers; }

int LadrunoSolidShell::getNumDOF(void) { return 24; }

//----------------------------------------------------------------------
// state cycle
//----------------------------------------------------------------------

int LadrunoSolidShell::commitState(void)
{
  int success = 0;
  if ((success = this->Element::commitState()) != 0)
    opserr << "LadrunoSolidShell::commitState - failed in base class\n";
  for (int i = 0; i < numGP; i++)
    success += materialPointers[i]->commitState();
  alphaCommit = alpha;
  return success;
}

int LadrunoSolidShell::revertToLastCommit(void)
{
  int success = 0;
  for (int i = 0; i < numGP; i++)
    success += materialPointers[i]->revertToLastCommit();
  alpha = alphaCommit;
  return success;
}

int LadrunoSolidShell::revertToStart(void)
{
  int success = 0;
  for (int i = 0; i < numGP; i++)
    success += materialPointers[i]->revertToStart();
  alpha.Zero();
  alphaCommit.Zero();
  return success;
}

// Set the per-GP trial states (the standard update() contract). formANS(0,...)
// solves the enhanced parameter alpha by an inner Newton and strains each GP
// to T*(Bcov_ANS*u) + G*alpha, so commitState captures the correct converged
// state under ANY algorithm — including a single-pass 'Linear' step where no
// force/tangent form runs at the final u (the LadrunoBrick::formEAStrue
// contract). The residual it also assembles is harmless.  // Ladruno
int LadrunoSolidShell::update(void)
{
  formANS(0, false);
  return 0;
}

//----------------------------------------------------------------------
// quadrature (through-thickness)
//----------------------------------------------------------------------

int LadrunoSolidShell::zRuleStatic(QuadZ q, int n, double *pts, double *wts)
{
  if (q == QuadZ::GAUSS) {
    switch (n) {
    case 2: {
      static const double p = 0.5773502691896257;
      pts[0] = -p; pts[1] = p; wts[0] = wts[1] = 1.0;
      return 0;
    }
    case 3: {
      static const double p = 0.7745966692414834;
      pts[0] = -p; pts[1] = 0.0; pts[2] = p;
      wts[0] = wts[2] = 0.5555555555555556; wts[1] = 0.8888888888888889;
      return 0;
    }
    case 4: {
      static const double p1 = 0.3399810435848563, p2 = 0.8611363115940526;
      static const double w1 = 0.6521451548625461, w2 = 0.3478548451374538;
      pts[0] = -p2; pts[1] = -p1; pts[2] = p1; pts[3] = p2;
      wts[0] = wts[3] = w2; wts[1] = wts[2] = w1;
      return 0;
    }
    case 5: {
      static const double p1 = 0.5384693101056831, p2 = 0.9061798459386640;
      static const double w0 = 0.5688888888888889;
      static const double w1 = 0.4786286704993665, w2 = 0.2369268850561891;
      pts[0] = -p2; pts[1] = -p1; pts[2] = 0.0; pts[3] = p1; pts[4] = p2;
      wts[0] = wts[4] = w2; wts[1] = wts[3] = w1; wts[2] = w0;
      return 0;
    }
    default:
      return -1;   // gauss supports nz in [2,5]
    }
  }

  // Lobatto (points on the faces; n>=3 — n=2 is the trapezoid rule, too weak
  // for the zeta^2 bending stiffness integrand: rejected at parse time too)
  switch (n) {
  case 3:
    pts[0] = -1.0; pts[1] = 0.0; pts[2] = 1.0;
    wts[0] = wts[2] = 1.0 / 3.0; wts[1] = 4.0 / 3.0;
    return 0;
  case 4: {
    static const double p = 0.4472135954999579;   // 1/sqrt(5)
    pts[0] = -1.0; pts[1] = -p; pts[2] = p; pts[3] = 1.0;
    wts[0] = wts[3] = 1.0 / 6.0; wts[1] = wts[2] = 5.0 / 6.0;
    return 0;
  }
  case 5: {
    static const double p = 0.6546536707079771;   // sqrt(3/7)
    pts[0] = -1.0; pts[1] = -p; pts[2] = 0.0; pts[3] = p; pts[4] = 1.0;
    wts[0] = wts[4] = 0.1; wts[1] = wts[3] = 49.0 / 90.0; wts[2] = 32.0 / 45.0;
    return 0;
  }
  case 6: {
    static const double p1 = 0.2852315164806451, p2 = 0.7650553239294647;
    static const double w1 = 0.5548583770354864, w2 = 0.3784749562978470;
    pts[0] = -1.0; pts[1] = -p2; pts[2] = -p1; pts[3] = p1; pts[4] = p2; pts[5] = 1.0;
    wts[0] = wts[5] = 1.0 / 15.0; wts[1] = wts[4] = w2; wts[2] = wts[3] = w1;
    return 0;
  }
  case 7: {
    static const double p1 = 0.4688487934707142, p2 = 0.8302238962785670;
    static const double w0 = 0.4876190476190476;
    static const double w1 = 0.4317453812098626, w2 = 0.2768260473615659;
    pts[0] = -1.0; pts[1] = -p2; pts[2] = -p1; pts[3] = 0.0;
    pts[4] = p1; pts[5] = p2; pts[6] = 1.0;
    wts[0] = wts[6] = 1.0 / 21.0; wts[1] = wts[5] = w2;
    wts[2] = wts[4] = w1; wts[3] = w0;
    return 0;
  }
  default:
    return -1;   // lobatto supports nz in [3,7]
  }
}

int LadrunoSolidShell::zRule(double *pts, double *wts) const
{
  int res = zRuleStatic(quadz, nz, pts, wts);
  if (res != 0)
    opserr << "LadrunoSolidShell::zRule - unsupported nz=" << nz
           << " for the selected -quadz rule (gauss: 2-5, lobatto: 3-7)\n";
  return res;
}

//----------------------------------------------------------------------
// geometry helpers
//----------------------------------------------------------------------

void LadrunoSolidShell::shapeFcn(const double xi[3], double N[8], double dN[8][3])
{
  for (int a = 0; a < 8; a++) {
    const double xa = natCoord[a][0], ya = natCoord[a][1], za = natCoord[a][2];
    const double fx = 1.0 + xa * xi[0];
    const double fy = 1.0 + ya * xi[1];
    const double fz = 1.0 + za * xi[2];
    N[a]     = 0.125 * fx * fy * fz;
    dN[a][0] = 0.125 * xa * fy * fz;
    dN[a][1] = 0.125 * fx * ya * fz;
    dN[a][2] = 0.125 * fx * fy * za;
  }
}

// covariant base J[i][j] = dx_i/dxi_j = sum_a dN_a/dxi_j * X_a[i]; returns det
double LadrunoSolidShell::jacobian(const double X8[8][3], const double dN[8][3],
                                   double J[3][3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int a = 0; a < 8; a++)
        s += dN[a][j] * X8[a][i];
      J[i][j] = s;
    }
  return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
       - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
       + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

void LadrunoSolidShell::invert3(const double J[3][3], double det, double Jinv[3][3])
{
  const double inv = 1.0 / det;
  Jinv[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv;
  Jinv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * inv;
  Jinv[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv;
  Jinv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * inv;
  Jinv[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv;
  Jinv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * inv;
  Jinv[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv;
  Jinv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * inv;
  Jinv[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv;
}

// 6x6 covariant -> Cartesian strain map, engineering Voigt both sides:
//   eps_ab = E_ij g^i_a g^j_b   with  g^i_a = Jinv[i][a]  (rows of J^-1).
// Row I = target pair (a,b), col Jc = source pair (i,j):
//   a==b:  T = q_ia q_ja
//   a!=b:  T = q_ia q_jb + q_ja q_ib
// (the source engineering-shear factor is already absorbed: source shear cols
// carry gamma_ij = 2E_ij, so their coefficient is q_ia q_ja / the symmetric sum)
void LadrunoSolidShell::buildT(const double Jinv[3][3], double T[6][6])
{
  for (int I = 0; I < 6; I++) {
    const int a = voigtPair[I][0], b = voigtPair[I][1];
    for (int Jc = 0; Jc < 6; Jc++) {
      const int i = voigtPair[Jc][0], j = voigtPair[Jc][1];
      if (a == b)
        T[I][Jc] = Jinv[i][a] * Jinv[j][a];
      else
        T[I][Jc] = Jinv[i][a] * Jinv[j][b] + Jinv[j][a] * Jinv[i][b];
    }
  }
}

// one covariant B row (engineering shear) at natural point xi:
//   normal row (i==j):  out[3a+k] = J[k][i] * dN_a/dxi_i
//   shear  row (i!=j):  out[3a+k] = J[k][i] * dN_a/dxi_j + J[k][j] * dN_a/dxi_i
void LadrunoSolidShell::covRow(int row, const double xi[3], double out[24]) const
{
  double N[8], dN[8][3], J[3][3];
  shapeFcn(xi, N, dN);
  (void)jacobian(X, dN, J);

  const int i = voigtPair[row][0], j = voigtPair[row][1];
  for (int a = 0; a < 8; a++)
    for (int k = 0; k < 3; k++) {
      if (i == j)
        out[3 * a + k] = J[k][i] * dN[a][i];
      else
        out[3 * a + k] = J[k][i] * dN[a][j] + J[k][j] * dN[a][i];
    }
}

// full 6x24 covariant B (+ J, detJ) at xi — no ANS substitution
void LadrunoSolidShell::covB(const double xi[3], double Bcov[6][24],
                             double J[3][3], double &detJ) const
{
  double N[8], dN[8][3];
  shapeFcn(xi, N, dN);
  detJ = jacobian(X, dN, J);

  for (int row = 0; row < 6; row++) {
    const int i = voigtPair[row][0], j = voigtPair[row][1];
    for (int a = 0; a < 8; a++)
      for (int k = 0; k < 3; k++) {
        if (i == j)
          Bcov[row][3 * a + k] = J[k][i] * dN[a][i];
        else
          Bcov[row][3 * a + k] = J[k][i] * dN[a][j] + J[k][j] * dN[a][i];
      }
  }
}

// Cartesian 6x24 B at a GP: covariant B, ANS-substituted rows (ANS only),
// then transformed with T built from the GP's Jinv.
void LadrunoSolidShell::cartB(const double xi[3], Matrix &B, double &detJ) const
{
  double Bcov[6][24], J[3][3];
  covB(xi, Bcov, J, detJ);

  if (formulation == Formulation::ANS) {
    double row[24], acc[24];

    // (2) E33 — Betsch-Stein bilinear corner tying at the midsurface zeta=0
    static const double cor[4][2] = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
    for (int c = 0; c < 24; c++) acc[c] = 0.0;
    for (int A = 0; A < 4; A++) {
      const double xit[3] = {cor[A][0], cor[A][1], 0.0};
      covRow(2, xit, row);
      const double w = 0.25 * (1.0 + cor[A][0] * xi[0]) * (1.0 + cor[A][1] * xi[1]);
      for (int c = 0; c < 24; c++) acc[c] += w * row[c];
    }
    for (int c = 0; c < 24; c++) Bcov[2][c] = acc[c];

    // (4) gamma_23 — Dvorkin-Bathe tying at (-+1, 0, zeta), linear in xi
    {
      const double xtA[3] = {-1.0, 0.0, xi[2]};
      const double xtB[3] = { 1.0, 0.0, xi[2]};
      double rowB[24];
      covRow(4, xtA, row);
      covRow(4, xtB, rowB);
      const double wA = 0.5 * (1.0 - xi[0]), wB = 0.5 * (1.0 + xi[0]);
      for (int c = 0; c < 24; c++) Bcov[4][c] = wA * row[c] + wB * rowB[c];
    }

    // (5) gamma_13 — Dvorkin-Bathe tying at (0, -+1, zeta), linear in eta
    {
      const double xtA[3] = {0.0, -1.0, xi[2]};
      const double xtB[3] = {0.0,  1.0, xi[2]};
      double rowB[24];
      covRow(5, xtA, row);
      covRow(5, xtB, rowB);
      const double wA = 0.5 * (1.0 - xi[1]), wB = 0.5 * (1.0 + xi[1]);
      for (int c = 0; c < 24; c++) Bcov[5][c] = wA * row[c] + wB * rowB[c];
    }
  }

  double Jinv[3][3], T[6][6];
  invert3(J, detJ, Jinv);
  buildT(Jinv, T);

  for (int I = 0; I < 6; I++)
    for (int c = 0; c < 24; c++) {
      double s = 0.0;
      for (int r = 0; r < 6; r++)
        s += T[I][r] * Bcov[r][c];
      B(I, c) = s;
    }
}

// cache the center map for the EAS mode: t0col3 = T0(:,3), j0det = detJ(0).
// The enhanced Cartesian strain is eps_enh(xi) = (j0/j) * zeta * t0col3 * alpha
// (Simo-Rifai center transform => int eps_enh dV = j0*t0col3 * SUM w*zeta = 0
// exactly for the symmetric product rule: the patch-test orthogonality).
void LadrunoSolidShell::buildEAS(void)
{
  const double xi0[3] = {0.0, 0.0, 0.0};
  double N[8], dN[8][3], J[3][3], Jinv[3][3], T[6][6];
  shapeFcn(xi0, N, dN);
  j0det = jacobian(X, dN, J);
  if (j0det <= 0.0) {
    opserr << "LadrunoSolidShell::buildEAS - element " << this->getTag()
           << ": non-positive center Jacobian (" << j0det
           << ") — check node ordering (1-4 bottom CCW, 5-8 top)\n";
    return;
  }
  invert3(J, j0det, Jinv);
  buildT(Jinv, T);
  for (int r = 0; r < 6; r++) t0col3[r] = T[r][2];
  easBuilt = true;
}

//----------------------------------------------------------------------
// formANS — the single assembly routine (adapted LadrunoBrick::formEAStrue).
// Solves the enhanced parameter by an inner Newton against the LIVE material
// stress, sets all trial states, assembles resid (always) and stiff
// (tang_flag==1) with static condensation of the enhanced mode.
//----------------------------------------------------------------------
void LadrunoSolidShell::formANS(int tang_flag, bool useInitialTangent)
{
  static const int    maxIters = 12;
  // RELATIVE inner-Newton tolerance (the LadrunoBrick::formEAStrue lesson):
  // ||int G^T sigma|| scales with E and mesh size, so an absolute tol is
  // unreachable for steel-scale problems. Converge on tolRel * first residual
  // + a tiny absolute floor for the alpha-already-satisfies case.  // Ladruno
  static const double tolRel = 1.0e-8;
  static const double tolAbs = 1.0e-12;

  const bool useEAS = (formulation == Formulation::ANS);
  if (useEAS && !easBuilt)
    buildEAS();

  stiff.Zero();
  resid.Zero();

  // trial displacement (global frame == core frame: v1 is -geom linear)
  static Vector u(24);
  for (int a = 0; a < 8; a++) {
    const Vector &d = nodePointers[a]->getTrialDisp();
    for (int k = 0; k < 3; k++) u(3 * a + k) = d(k);
  }

  // quadrature: 2x2 in-plane Gauss (weight 1) x nz through thickness
  double zp[NZ_MAX], zw[NZ_MAX];
  if (zRule(zp, zw) != 0) return;
  static const double g2 = 0.5773502691896257;
  static const double ip[4][2] = {{-g2, -g2}, {g2, -g2}, {g2, g2}, {-g2, g2}};

  static Matrix B(6, 24);
  static Matrix dd(6, 6), DB(6, 24);
  static Vector strain(6), stressV(6), Gv(6), DG(6);
  static Matrix Kda(24, NEAS), Kad(NEAS, 24);

  // -------- getInitialStiff: condensed initial tangent at alpha = 0 --------
  if (useInitialTangent) {
    double Kaa = 0.0;
    Kda.Zero(); Kad.Zero();
    int gp = 0;
    for (int l = 0; l < nz; l++)
      for (int q = 0; q < 4; q++, gp++) {
        const double xi[3] = {ip[q][0], ip[q][1], zp[l]};
        double detJ;
        cartB(xi, B, detJ);
        const double dvol = zw[l] * detJ;
        dd = materialPointers[gp]->getInitialTangent();
        dd *= dvol;
        DB.addMatrixProduct(0.0, dd, B, 1.0);
        stiff.addMatrixTransposeProduct(1.0, B, DB, 1.0);      // Kdd
        if (useEAS) {
          const double gs = zp[l] * j0det / detJ;
          for (int r = 0; r < 6; r++) Gv(r) = gs * t0col3[r];
          DG.addMatrixVector(0.0, dd, Gv, 1.0);
          for (int c = 0; c < 24; c++) {
            double sad = 0.0, sda = 0.0;
            for (int r = 0; r < 6; r++) {
              sad += Gv(r) * DB(r, c);      // Kad += G^T dd B
              sda += B(r, c) * DG(r);       // Kda += B^T dd G
            }
            Kad(0, c) += sad;
            Kda(c, 0) += sda;
          }
          double s = 0.0;
          for (int r = 0; r < 6; r++) s += Gv(r) * DG(r);
          Kaa += s;                          // Kaa += G^T dd G
        }
      }
    if (useEAS) {
      if (Kaa > 0.0) {
        const double inv = 1.0 / Kaa;
        for (int i = 0; i < 24; i++)
          for (int j = 0; j < 24; j++)
            stiff(i, j) -= Kda(i, 0) * inv * Kad(0, j);
      } else {
        opserr << "LadrunoSolidShell::formANS - element " << this->getTag()
               << ": non-positive initial Kaa (" << Kaa
               << "); skipping the EAS condensation\n";
      }
    }
    return;
  }

  // -------- inner Newton: solve alpha s.t. int G^T sigma dV = 0 (u fixed) ----
  // Each pass sets ALL trial strains with the current alpha, so on break the
  // material states are consistent with the final alpha (the update contract).
  double Kaa = 0.0;
  int count = 0;
  double r0 = -1.0;
  while (true) {
    double h = 0.0;      // int G^T sigma dV
    Kaa = 0.0;
    int gp = 0;
    for (int l = 0; l < nz; l++)
      for (int q = 0; q < 4; q++, gp++) {
        const double xi[3] = {ip[q][0], ip[q][1], zp[l]};
        double detJ;
        cartB(xi, B, detJ);
        const double dvol = zw[l] * detJ;
        strain.addMatrixVector(0.0, B, u, 1.0);                // eps = B*u
        double gs = 0.0;
        if (useEAS) {
          gs = zp[l] * j0det / detJ;
          for (int r = 0; r < 6; r++) {
            Gv(r) = gs * t0col3[r];
            strain(r) += Gv(r) * alpha(0);                     //     + G*alpha
          }
        }
        materialPointers[gp]->setTrialStrain(strain);
        if (!useEAS)
          continue;
        stressV = materialPointers[gp]->getStress();
        dd = materialPointers[gp]->getTangent();
        dd *= dvol;
        for (int r = 0; r < 6; r++) h += Gv(r) * stressV(r) * dvol;
        DG.addMatrixVector(0.0, dd, Gv, 1.0);
        double s = 0.0;
        for (int r = 0; r < 6; r++) s += Gv(r) * DG(r);
        Kaa += s;
      }

    if (!useEAS)
      break;                                    // STD: trial states set; done

    const double r = fabs(h);
    if (count == 0) r0 = r;
    if (count >= 1 && r <= tolRel * r0 + tolAbs)
      break;
    if (count >= maxIters) {
      opserr << "LadrunoSolidShell::formANS - element " << this->getTag()
             << ": enhanced-strain Newton did not converge in " << maxIters
             << " iters (|h|=" << r << ", r0=" << r0 << ")\n";
      break;
    }
    if (Kaa <= 0.0) {
      // an indefinite enhanced block (deep softening) — freeze alpha for this
      // pass rather than stepping through a singular solve; the G6 stability
      // gate (ADR 66, P5.2) owns the systematic treatment.  // Ladruno
      opserr << "LadrunoSolidShell::formANS - element " << this->getTag()
             << ": non-positive Kaa (" << Kaa << "); freezing alpha this pass\n";
      break;
    }
    alpha(0) -= h / Kaa;                        // dalpha = -Kaa^-1 h
    count++;
  }

  // -------- final assembly at converged alpha (Kaa preserved) --------
  Kda.Zero(); Kad.Zero();
  int gp = 0;
  for (int l = 0; l < nz; l++)
    for (int q = 0; q < 4; q++, gp++) {
      const double xi[3] = {ip[q][0], ip[q][1], zp[l]};
      double detJ;
      cartB(xi, B, detJ);
      const double dvol = zw[l] * detJ;
      stressV = materialPointers[gp]->getStress();
      stressV *= dvol;
      resid.addMatrixTransposeVector(1.0, B, stressV, 1.0);    // f_int += B^T sigma
      if (tang_flag == 1) {
        dd = materialPointers[gp]->getTangent();
        dd *= dvol;
        DB.addMatrixProduct(0.0, dd, B, 1.0);
        stiff.addMatrixTransposeProduct(1.0, B, DB, 1.0);      // Kdd
        if (useEAS) {
          const double gs = zp[l] * j0det / detJ;
          for (int r = 0; r < 6; r++) Gv(r) = gs * t0col3[r];
          DG.addMatrixVector(0.0, dd, Gv, 1.0);
          for (int c = 0; c < 24; c++) {
            double sad = 0.0, sda = 0.0;
            for (int r = 0; r < 6; r++) {
              sad += Gv(r) * DB(r, c);
              sda += B(r, c) * DG(r);
            }
            Kad(0, c) += sad;
            Kda(c, 0) += sda;
          }
        }
      }
    }

  if (tang_flag == 1 && useEAS) {
    if (Kaa > 0.0) {
      const double inv = 1.0 / Kaa;
      for (int i = 0; i < 24; i++)
        for (int j = 0; j < 24; j++)
          stiff(i, j) -= Kda(i, 0) * inv * Kad(0, j);          // K* = Kdd - Kda Kaa^-1 Kad
    } else {
      opserr << "LadrunoSolidShell::formANS - element " << this->getTag()
             << ": non-positive Kaa (" << Kaa
             << ") at assembly; skipping the EAS condensation\n";
    }
  }
}

//----------------------------------------------------------------------
// matrices / vectors
//----------------------------------------------------------------------

const Matrix &LadrunoSolidShell::getTangentStiff(void)
{
  formANS(1, false);
  return stiff;
}

const Matrix &LadrunoSolidShell::getInitialStiff(void)
{
  if (Ki != 0)
    return *Ki;
  formANS(1, true);
  Ki = new Matrix(stiff);
  return *Ki;
}

void LadrunoSolidShell::formInertiaTerms(int tangFlag)
{
  // zero the shared static mass FIRST, so a zRule failure (bad nz) returns a
  // clean zero matrix rather than a stale one from another element's form.
  if (tangFlag == 1)
    mass.Zero();

  double zp[NZ_MAX], zw[NZ_MAX];
  if (zRule(zp, zw) != 0) return;
  static const double g2 = 0.5773502691896257;
  static const double ip[4][2] = {{-g2, -g2}, {g2, -g2}, {g2, g2}, {-g2, g2}};

  double N[8], dN[8][3], J[3][3];
  int gp = 0;
  for (int l = 0; l < nz; l++)
    for (int q = 0; q < 4; q++, gp++) {
      const double rho = materialPointers[gp]->getRho();
      if (rho == 0.0)
        continue;
      const double xi[3] = {ip[q][0], ip[q][1], zp[l]};
      shapeFcn(xi, N, dN);
      const double detJ = jacobian(X, dN, J);
      const double dm = rho * zw[l] * detJ;

      if (tangFlag == 1) {
        if (massType == 1) {
          // -lumped: row-sum diagonal m_a = rho * int N_a dV (sum_b N_b = 1
          // absorbed) == HRZ on the trilinear shape; always positive. The
          // explicit-path mass (ADR 66 G9): system Diagonal then carries the
          // FULL element mass instead of the 8/27 the raw consistent diagonal
          // holds, and the '-lump hrz'/'rowsum' pencil matches the runtime.
          // The LadrunoBrick donor pattern (LadrunoBrick.cpp:767-773).  // Ladruno
          for (int a = 0; a < 8; a++) {
            const double m = dm * N[a];
            for (int k = 0; k < 3; k++)
              mass(3 * a + k, 3 * a + k) += m;
          }
        } else {
          for (int a = 0; a < 8; a++)
            for (int b = 0; b < 8; b++) {
              const double m = dm * N[a] * N[b];
              for (int k = 0; k < 3; k++)
                mass(3 * a + k, 3 * b + k) += m;
            }
        }
      } else if (massType == 1) {
        // lumped inertial force: m_a * accel_a (nodal accel, matching the
        // diagonal mass above)  // Ladruno
        for (int a = 0; a < 8; a++) {
          const Vector &aa = nodePointers[a]->getTrialAccel();
          for (int k = 0; k < 3; k++)
            resid(3 * a + k) += dm * N[a] * aa(k);
        }
      } else {
        // inertial force: rho * N_a * (sum_b N_b * accel_b)
        double acc[3] = {0.0, 0.0, 0.0};
        for (int b = 0; b < 8; b++) {
          const Vector &ab = nodePointers[b]->getTrialAccel();
          for (int k = 0; k < 3; k++) acc[k] += N[b] * ab(k);
        }
        for (int a = 0; a < 8; a++)
          for (int k = 0; k < 3; k++)
            resid(3 * a + k) += dm * N[a] * acc[k];
      }
    }
}

const Matrix &LadrunoSolidShell::getMass(void)
{
  // Ladruno (ADR-77 G2 ext): per-instance mass cache -- LadrunoMassCache.h;
  // this element is the brick pattern verbatim (formInertiaTerms(1) into a
  // CLASS-STATIC `mass`), so it gets the brick treatment. Signature = the
  // numGP (= 4*nz) material rhos; coords guarded inside; massType fixed at
  // construction. Bounded stack signature: numGP > 64 (nz > 16) falls back to
  // the uncached path rather than allocating -- far beyond any real nz.
  double mcSig[64];
  const bool mcOk = (numGP <= 64);
  if (mcOk) {
    for (int i = 0; i < numGP; i++)
      mcSig[i] = materialPointers[i]->getRho();
    if (const Matrix *Mc = massCache.lookup(mcSig, numGP, nodePointers, 8, 3))
      return *Mc;
  }

  formInertiaTerms(1);

  if (mcOk)
    massCache.fill(mass, mcSig, numGP, nodePointers, 8, 3);
  return mass;
}

void LadrunoSolidShell::zeroLoad(void)
{
  if (load != 0)
    load->Zero();
}

int LadrunoSolidShell::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "LadrunoSolidShell::addLoad - element " << this->getTag()
         << " does not support elemental loads in v1 (surface pressure is the "
            "ADR 66 P5.3 deliverable); use nodal loads\n";
  return -1;
}

int LadrunoSolidShell::addInertiaLoadToUnbalance(const Vector &accel)
{
  bool haveRho = false;
  for (int i = 0; i < numGP; i++)
    if (materialPointers[i]->getRho() != 0.0) {
      haveRho = true;
      break;
    }
  if (!haveRho)
    return 0;

  // Ladruno (ADR-77 G2 ext): route through getMass() so the cached M is used
  // (this function reads only the mass matrix; ra is overwritten below).
  const Matrix &M = this->getMass();

  static Vector ra(24);
  for (int i = 0; i < 8; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int k = 0; k < 3; k++) ra(3 * i + k) = Raccel(k);
  }

  if (load == 0)
    load = new Vector(24);
  load->addMatrixVector(1.0, M, ra, -1.0);
  return 0;
}

const Vector &LadrunoSolidShell::getResistingForce(void)
{
  formANS(0, false);
  if (load != 0)
    resid -= *load;
  return resid;
}

const Vector &LadrunoSolidShell::getResistingForceIncInertia(void)
{
  formANS(0, false);        // fills the shared static resid with f_int
  formInertiaTerms(0);      // adds M*accel into resid

  // SNAPSHOT resid into a local BEFORE getRayleighDampingForces(): with
  // stiffness-proportional Rayleigh (betaK/betaK0) that call re-enters
  // getTangentStiff()/getInitialStiff() -> formANS() -> resid.Zero(), which
  // would silently destroy the f_int + M*accel we just accumulated (the
  // returned unbalance would drop the inertia term and Newton would diverge /
  // converge to a wrong dynamics solution). This is the LadrunoBrick donor
  // pattern (LadrunoBrick.cpp:688-700). The dynamic-Rayleigh regression gate
  // (test_dynamic_rayleigh_preserves_inertia) fails without this snapshot.  // Ladruno
  static Vector res(24);
  res = resid;

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  if (load != 0)
    res -= *load;
  return res;
}

//----------------------------------------------------------------------
// characteristic length — IN-PLANE projected sqrt(midsurface area) (ADR 66 D6)
//----------------------------------------------------------------------
double LadrunoSolidShell::getCharacteristicLength(void)
{
  static const double g2 = 0.5773502691896257;
  static const double ip[4][2] = {{-g2, -g2}, {g2, -g2}, {g2, g2}, {-g2, g2}};

  double A = 0.0;
  double N[8], dN[8][3], J[3][3];
  for (int q = 0; q < 4; q++) {
    const double xi[3] = {ip[q][0], ip[q][1], 0.0};
    shapeFcn(xi, N, dN);
    (void)jacobian(X, dN, J);
    // |g1 x g2| at the midsurface (columns 0,1 of J)
    const double cx = J[1][0] * J[2][1] - J[2][0] * J[1][1];
    const double cy = J[2][0] * J[0][1] - J[0][0] * J[2][1];
    const double cz = J[0][0] * J[1][1] - J[1][0] * J[0][1];
    A += sqrt(cx * cx + cy * cy + cz * cz);   // in-plane weights = 1
  }
  if (A <= 0.0)
    return this->Element::getCharacteristicLength();   // degenerate fallback
  return sqrt(A);
}

//----------------------------------------------------------------------
// parallel / database
//----------------------------------------------------------------------

int LadrunoSolidShell::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  const int dataTag = this->getDbTag();

  // fixed-size ID (numGP <= 4*NZ_MAX = 28): header(7) + nodes(8) + 2*28
  static ID idData(71);
  idData.Zero();
  idData(0) = this->getTag();
  idData(1) = nz;
  idData(2) = static_cast<int>(quadz);
  idData(3) = static_cast<int>(formulation);
  idData(4) = (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ? 1 : 0;
  idData(5) = numGP;
  idData(6) = massType;   // was the spare schema slot; 0 == consistent, so
                          // pre-(-lumped) streams read back byte-compatibly  // Ladruno

  for (int i = 0; i < 8; i++)
    idData(7 + i) = connectedExternalNodes(i);

  for (int i = 0; i < numGP; i++) {
    idData(15 + i) = materialPointers[i]->getClassTag();
    int matDbTag = materialPointers[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(15 + 28 + i) = matDbTag;
  }

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "LadrunoSolidShell::sendSelf - " << this->getTag()
           << " failed to send ID\n";
    return res;
  }

  static Vector dData(4 + NEAS);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  for (int i = 0; i < NEAS; i++)
    dData(4 + i) = alphaCommit(i);

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoSolidShell::sendSelf - failed to send double data\n";
    return -1;
  }

  for (int i = 0; i < numGP; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "LadrunoSolidShell::sendSelf - " << this->getTag()
             << " failed to send material " << i << "\n";
      return res;
    }
  }
  return res;
}

int LadrunoSolidShell::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
  int res = 0;
  const int dataTag = this->getDbTag();

  static ID idData(71);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "LadrunoSolidShell::recvSelf - failed to recv ID\n";
    return res;
  }

  this->setTag(idData(0));
  const int nzNew    = idData(1);
  const int numGPNew = idData(5);

  // Validate stream-derived control fields BEFORE using numGPNew as an
  // allocation size and index bound into the fixed 71-slot ID (the material
  // classTag/dbTag blocks live at 15..15+numGP and 43..43+numGP). A corrupt or
  // version-skewed record with numGP outside [4*NZ_MIN, 4*NZ_MAX] or a bad enum
  // would otherwise overflow the ID or mis-cast. Refuse the record.  // Ladruno
  if (nzNew < NZ_MIN || nzNew > NZ_MAX || numGPNew != 4 * nzNew ||
      idData(2) < 0 || idData(2) > 1 || idData(3) < 0 || idData(3) > 1 ||
      idData(6) < 0 || idData(6) > 1) {
    opserr << "LadrunoSolidShell::recvSelf - element " << idData(0)
           << ": inconsistent stream (nz=" << nzNew << ", numGP=" << numGPNew
           << ", quadz=" << idData(2) << ", formulation=" << idData(3)
           << ", massType=" << idData(6) << ")\n";
    return -1;
  }
  quadz              = static_cast<QuadZ>(idData(2));
  formulation        = static_cast<Formulation>(idData(3));
  massType           = idData(6);

  for (int i = 0; i < 8; i++)
    connectedExternalNodes(i) = idData(7 + i);

  static Vector dData(4 + NEAS);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "LadrunoSolidShell::recvSelf - failed to recv double data\n";
    return -1;
  }
  alphaM = dData(0);
  betaK  = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  for (int i = 0; i < NEAS; i++)
    alphaCommit(i) = dData(4 + i);
  alpha = alphaCommit;

  // (re)allocate the material array if the GP count changed
  if (materialPointers != 0 && numGP != numGPNew) {
    for (int i = 0; i < numGP; i++)
      if (materialPointers[i] != 0)
        delete materialPointers[i];
    delete[] materialPointers;
    materialPointers = 0;
  }
  nz    = nzNew;
  numGP = numGPNew;
  if (materialPointers == 0) {
    materialPointers = new NDMaterial *[numGP];
    for (int i = 0; i < numGP; i++)
      materialPointers[i] = 0;
  }

  for (int i = 0; i < numGP; i++) {
    const int matClassTag = idData(15 + i);
    const int matDbTag    = idData(15 + 28 + i);
    if (materialPointers[i] == 0 ||
        materialPointers[i]->getClassTag() != matClassTag) {
      if (materialPointers[i] != 0)
        delete materialPointers[i];
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
        opserr << "LadrunoSolidShell::recvSelf - broker could not create "
                  "NDMaterial classTag " << matClassTag << "\n";
        return -1;
      }
    }
    materialPointers[i]->setDbTag(matDbTag);
    res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "LadrunoSolidShell::recvSelf - material " << i
             << " failed to recv itself\n";
      return res;
    }
  }

  easBuilt = false;   // center map rebuilt in setDomain (geometry-only)
  // drop any cached initial stiffness: nz/formulation/material class may all
  // differ from the pre-recv element, so a stale Ki would be wrong.  // Ladruno
  if (Ki != 0) { delete Ki; Ki = 0; }
  // same for the mass cache: massType/quadz are sig-exempt (construction-fixed
  // in normal flows) but a restore into a live element CAN flip them with rho
  // and coords unchanged -- a guard hit would then serve the pre-recv mass
  // (ADR-77 review wave).  // Ladruno
  massCache.invalidate();
  return res;
}

//----------------------------------------------------------------------
// output
//----------------------------------------------------------------------

void LadrunoSolidShell::Print(OPS_Stream &s, int flag)
{
  // a broker-constructed element (before recvSelf) has numGP==0 and
  // materialPointers==0; guard the material dereference.  // Ladruno
  const bool haveMat = (materialPointers != 0 && numGP > 0 &&
                        materialPointers[0] != 0);

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"LadrunoSolidShell\", ";
    s << "\"nodes\": [";
    for (int i = 0; i < 8; i++) {
      s << connectedExternalNodes(i);
      if (i < 7) s << ", ";
    }
    s << "], ";
    s << "\"nz\": " << nz << ", ";
    s << "\"quadz\": \"" << (quadz == QuadZ::GAUSS ? "gauss" : "lobatto") << "\", ";
    s << "\"formulation\": \"" << (formulation == Formulation::ANS ? "ans" : "std") << "\", ";
    s << "\"mass\": \"" << (massType == 1 ? "lumped" : "consistent") << "\", ";
    s << "\"material\": \"";
    if (haveMat) s << materialPointers[0]->getTag();
    s << "\"}";
    return;
  }

  s << "LadrunoSolidShell — 8-node ANS/EAS solid-shell (ADR 66)\n";
  s << "  tag: " << this->getTag() << "\n";
  s << "  formulation: " << (formulation == Formulation::ANS ? "ans" : "std")
    << "   nz: " << nz << " ("
    << (quadz == QuadZ::GAUSS ? "gauss" : "lobatto") << ")"
    << "   mass: " << (massType == 1 ? "lumped (row-sum/HRZ)" : "consistent")
    << "\n";
  s << "  connected nodes: " << connectedExternalNodes;
  if (haveMat)
    s << "  material: " << materialPointers[0]->getTag() << "\n";
}

Response *LadrunoSolidShell::setResponse(const char **argv, int argc,
                                         OPS_Stream &output)
{
  Response *theResponse = 0;

  if (argc < 1) return 0;

  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoSolidShell");
  output.attr("eleTag", this->getTag());
  for (int i = 0; i < 8; i++) {
    char nodeTag[16];
    sprintf(nodeTag, "node%d", i + 1);
    output.attr(nodeTag, connectedExternalNodes(i));
  }

  if (LadrunoResp::is(argv[0], "material")) {
    if (argc >= 2) {
      const int pointNum = atoi(argv[1]);
      if (pointNum > 0 && pointNum <= numGP) {
        output.tag("GaussPoint");
        output.attr("number", pointNum);
        theResponse = materialPointers[pointNum - 1]->setResponse(
            &argv[2], argc - 2, output);
        output.endTag();
      }
    }
  }
  else if (LadrunoResp::is(argv[0], "stress")) {
    theResponse = new ElementResponse(this, 1, Vector(6 * numGP));
  }
  else if (LadrunoResp::is(argv[0], "strain")) {
    theResponse = new ElementResponse(this, 2, Vector(6 * numGP));
  }
  else if (LadrunoResp::is(argv[0], "alpha")) {
    // enhanced-parameter diagnostic (the G6 stability gate reads this)
    theResponse = new ElementResponse(this, 3, Vector(NEAS));
  }
  else if (LadrunoResp::is(argv[0], "stiff")) {
    // condensed tangent stiffness, row-major 24x24 (diagnostic — lets a gate
    // compare against the initial stiffness and check symmetry)  // Ladruno
    theResponse = new ElementResponse(this, 4, Vector(24 * 24));
  }
  else if (LadrunoResp::is(argv[0], "stiffInitial")) {
    theResponse = new ElementResponse(this, 5, Vector(24 * 24));
  }
  else if (LadrunoResp::is(argv[0], "force")) {
    // Ladruno — the element had NO force response at all; the family contract
    // is force / stress / strain / stiff / stiffInitial / charLength.
    theResponse = new ElementResponse(this, 6, Vector(24));
  }
  else if (LadrunoResp::is(argv[0], "charLength")) {
    // Ladruno — the IN-PLANE projected size handed to crack-band materials.
    output.tag("ResponseType", "lch");
    theResponse = new ElementResponse(this, 7, Vector(1));
  }

  output.endTag();

  // Ladruno — base vocabulary (globalForce, dampingForce, dynamicForce,
  // inertialForce); Element::setResponse opens its own ElementOutput tag, so
  // this MUST come after endTag().
  if (theResponse == 0)
    return this->Element::setResponse(argv, argc, output);
  return theResponse;
}

int LadrunoSolidShell::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {
    Vector data(6 * numGP);
    for (int i = 0; i < numGP; i++) {
      const Vector &sig = materialPointers[i]->getStress();
      for (int r = 0; r < 6; r++)
        data(6 * i + r) = sig(r);
    }
    return eleInfo.setVector(data);
  }
  if (responseID == 2) {
    Vector data(6 * numGP);
    for (int i = 0; i < numGP; i++) {
      const Vector &eps = materialPointers[i]->getStrain();
      for (int r = 0; r < 6; r++)
        data(6 * i + r) = eps(r);
    }
    return eleInfo.setVector(data);
  }
  if (responseID == 3)
    return eleInfo.setVector(alpha);
  if (responseID == 4 || responseID == 5) {
    const Matrix &K = (responseID == 4) ? this->getTangentStiff()
                                        : this->getInitialStiff();
    Vector data(24 * 24);
    for (int i = 0; i < 24; i++)
      for (int j = 0; j < 24; j++)
        data(24 * i + j) = K(i, j);
    return eleInfo.setVector(data);
  }
  if (responseID == 6)
    return eleInfo.setVector(this->getResistingForce());
  if (responseID == 7) {
    static Vector lch(1);
    lch(0) = this->getCharacteristicLength();
    return eleInfo.setVector(lch);
  }

  return this->Element::getResponse(responseID, eleInfo);
}
