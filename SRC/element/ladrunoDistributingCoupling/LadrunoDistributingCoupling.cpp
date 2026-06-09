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

// Ladruno: RBE3 / distributing-coupling element (ADR 28). See
// LadrunoDistributingCoupling.h and
// Ladruno_implementation/28_ladruno_distributing_coupling_rbe3_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoDistributingCoupling.h>
#include <LadrunoEmbeddedKernel.h>   // shared coupling-kernel (translation block + bipenalty math)
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>

// ---------------------------------------------------------------------------
// cyclic Jacobi eigensolver for a 3x3 symmetric matrix (Numerical Recipes form).
// a[][] is destroyed; d[] = eigenvalues, v[][] columns = eigenvectors (v[j][k] =
// component j of eigenvector k). Robust for the small SPD/PSD I_c here.
// ---------------------------------------------------------------------------
static void jacobi3(double a[3][3], double d[3], double v[3][3])
{
  for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) v[i][j] = (i == j) ? 1.0 : 0.0; }
  double b[3], z[3];
  for (int i = 0; i < 3; i++) { d[i] = a[i][i]; b[i] = d[i]; z[i] = 0.0; }
#define JROT(M,i,j,k,l) { double g_=M[i][j]; double h_=M[k][l]; \
                          M[i][j]=g_-s*(h_+g_*tau); M[k][l]=h_+s*(g_-h_*tau); }
  for (int sweep = 0; sweep < 50; sweep++) {
    double sm = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (sm == 0.0) { return; }
    double tresh = (sweep < 3) ? (0.2 * sm / 9.0) : 0.0;
    for (int p = 0; p < 2; p++) {
      for (int q = p + 1; q < 3; q++) {
        double g = 100.0 * fabs(a[p][q]);
        if (sweep > 3 && (fabs(d[p]) + g) == fabs(d[p]) && (fabs(d[q]) + g) == fabs(d[q])) {
          a[p][q] = 0.0;
        } else if (fabs(a[p][q]) > tresh) {
          double h = d[q] - d[p], t;
          if ((fabs(h) + g) == fabs(h)) t = a[p][q] / h;
          else {
            double theta = 0.5 * h / a[p][q];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) t = -t;
          }
          double c = 1.0 / sqrt(1.0 + t * t), s = t * c, tau = s / (1.0 + c);
          double hh = t * a[p][q];
          z[p] -= hh; z[q] += hh; d[p] -= hh; d[q] += hh; a[p][q] = 0.0;
          for (int j = 0;     j < p; j++) JROT(a, j, p, j, q)
          for (int j = p + 1; j < q; j++) JROT(a, p, j, j, q)
          for (int j = q + 1; j < 3; j++) JROT(a, p, j, q, j)
          for (int j = 0;     j < 3; j++) JROT(v, j, p, j, q)
        }
      }
    }
    for (int i = 0; i < 3; i++) { b[i] += z[i]; d[i] = b[i]; z[i] = 0.0; }
  }
#undef JROT
}

// ===========================================================================
//  construction
// ===========================================================================
LadrunoDistributingCoupling::LadrunoDistributingCoupling(int tag, int ndm_, int refNode,
        const ID& indepNodes, const Vector& weights_, double kt,
        double kr, bool krUser_, int enforce_,
        bool bipenalty_, int bpMode_, double bpDt_, double bpBeta_,
        double kAlpha_, int hostEleTag_, bool ktAuto_, bool initGapCapture_)
  : Element(tag, ELE_TAG_LadrunoDistributingCoupling),
    ndm(ndm_), nrot((ndm_ == 3) ? 3 : 1), nIndep(indepNodes.Size()),
    connectedNodes(1 + indepNodes.Size()), weights(weights_),
    Kt(kt), Kr(kr), krUser(krUser_), ktAuto(ktAuto_), kAlpha(kAlpha_),
    hostEleTag(hostEleTag_), ktResolved(false),
    enforce(enforce_), lambda(ndm_), lambda_r((ndm_ == 3) ? 3 : 1),
    bipenalty(bipenalty_), bpMode(bpMode_), bpDt(bpDt_), bpBeta(bpBeta_),
    mPenalty(0.0), iPenalty(0.0), bpResolved(false),
    valid(false), W(0.0), ell2(0.0), dRef(ndm_), rvec(),
    Icplus((ndm_ == 3) ? 3 : 1, (ndm_ == 3) ? 3 : 1),
    Pproj((ndm_ == 3) ? 3 : 1, (ndm_ == 3) ? 3 : 1), nKept(0),
    nDOF(0), nodeNdf(1 + indepNodes.Size()), dofOffset(1 + indepNodes.Size()),
    B(0), initGapCapture(initGapCapture_), g0Computed(false),
    g0(ndm_), gr0((ndm_ == 3) ? 3 : 1),
    theNodes(0), K(0), P(0), M0(0), C0(0), dampF(0)
{
  connectedNodes(0) = refNode;
  for (int i = 0; i < nIndep; i++) connectedNodes(1 + i) = indepNodes(i);
  lambda.Zero();
  lambda_r.Zero();
  g0.Zero();
  gr0.Zero();
  theNodes = new Node*[1 + nIndep];
  for (int i = 0; i < 1 + nIndep; i++) theNodes[i] = 0;
}

LadrunoDistributingCoupling::LadrunoDistributingCoupling()
  : Element(0, ELE_TAG_LadrunoDistributingCoupling),
    ndm(0), nrot(0), nIndep(0), connectedNodes(), weights(),
    Kt(0.0), Kr(0.0), krUser(false), ktAuto(false), kAlpha(0.0),
    hostEleTag(-1), ktResolved(false),
    enforce(0), lambda(), lambda_r(),
    bipenalty(false), bpMode(0), bpDt(0.0), bpBeta(0.0),
    mPenalty(0.0), iPenalty(0.0), bpResolved(false),
    valid(false), W(0.0), ell2(0.0), dRef(), rvec(),
    Icplus(), Pproj(), nKept(0),
    nDOF(0), nodeNdf(), dofOffset(),
    B(0), initGapCapture(true), g0Computed(false), g0(), gr0(),
    theNodes(0), K(0), P(0), M0(0), C0(0), dampF(0)
{
}

LadrunoDistributingCoupling::~LadrunoDistributingCoupling()
{
  if (theNodes != 0) delete[] theNodes;
  if (B != 0) delete B;
  if (K != 0) delete K;
  if (P != 0) delete P;
  if (M0 != 0) delete M0;
  if (C0 != 0) delete C0;
  if (dampF != 0) delete dampF;
}

// ===========================================================================
//  domain
// ===========================================================================
int LadrunoDistributingCoupling::getNumExternalNodes(void) const { return 1 + nIndep; }
const ID& LadrunoDistributingCoupling::getExternalNodes(void) { return connectedNodes; }
Node** LadrunoDistributingCoupling::getNodePtrs(void) { return theNodes; }
int LadrunoDistributingCoupling::getNumDOF(void) { return nDOF; }

void LadrunoDistributingCoupling::allocate(void)
{
  if (K != 0) delete K;       K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;       P = new Vector(nDOF);
  if (M0 != 0) delete M0;     M0 = new Matrix(nDOF, nDOF);
  if (C0 != 0) delete C0;     C0 = new Matrix(nDOF, nDOF);   // always zero (getDamp)
  if (dampF != 0) delete dampF; dampF = new Vector(nDOF);    // always zero
}

void LadrunoDistributingCoupling::setDomain(Domain* theDomain)
{
  if (theDomain == 0) {
    if (theNodes != 0)
      for (int i = 0; i < 1 + nIndep; i++) theNodes[i] = 0;
    return;
  }
  // resolve nodes + lay out the per-node DOF offsets. The reference may carry more
  // DOFs (rotations) than the independents; the reference needs ndf >= ndm+nrot
  // (refuse — RBE3 without rotation is a different element), independents need ndf >= ndm.
  nodeNdf.resize(1 + nIndep);
  dofOffset.resize(1 + nIndep);
  int pos = 0;
  for (int i = 0; i < 1 + nIndep; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "LadrunoDistributingCoupling::setDomain - node " << connectedNodes(i)
             << " not found\n";
      return;
    }
    int ndf = theNodes[i]->getNumberDOF();
    int need = (i == 0) ? (ndm + nrot) : ndm;
    if (ndf < need) {
      opserr << "LadrunoDistributingCoupling::setDomain - "
             << (i == 0 ? "reference" : "independent") << " node " << connectedNodes(i)
             << " has " << ndf << " DOFs, needs >= " << need
             << (i == 0 ? " (translations + rotations)\n" : " (translations)\n");
      return;
    }
    nodeNdf(i) = ndf;
    dofOffset(i) = pos;
    pos += ndf;
  }
  nDOF = pos;

  this->resolveGeometry();   // x_c, r_i, dRef, W, ell2, I_c eigen -> Icplus/Pproj, nKept

  // Allocate B/K/P/M0 even if the geometry is ill-posed (W<=0): a zeroed B makes the
  // element inert (contributes nothing) instead of leaving null matrices that the
  // assembler would dereference. The WARNING was already printed in resolveGeometry.
  if (B != 0) delete B;
  B = new Matrix(ndm + nrot, nDOF);
  B->Zero();
  this->allocate();
  if (valid) {
    this->buildB();
    this->captureInitialGap();   // born stress-free at the (possibly deformed) current state
  }
  this->DomainComponent::setDomain(theDomain);
}

// resolve the RBE3 geometry from the node coordinates (known at setDomain): the
// weighted centroid x_c, the relative positions r_i, the transport lever dRef =
// x_c - x_R, the rotation scale ell2 = Sum w_i|r_i|^2 / W, and the spectral
// pseudo-inverse I_c+ / range projector P of the weighted position-inertia I_c
// (dropping collinear/degenerate rotation axes, ADR 28 sec 3).
void LadrunoDistributingCoupling::resolveGeometry(void)
{
  valid = false;
  W = 0.0;
  double wmax = 0.0;
  for (int i = 0; i < nIndep; i++) {
    if (weights(i) < 0.0) {
      opserr << "LadrunoDistributingCoupling " << this->getTag()
             << ": negative weight w(" << i << ")=" << weights(i) << "\n";
      return;
    }
    W += weights(i);
    if (weights(i) > wmax) wmax = weights(i);
  }
  if (W <= 1.0e-12 * (wmax > 0.0 ? wmax : 1.0)) {
    opserr << "LadrunoDistributingCoupling " << this->getTag()
           << ": total weight W=" << W << " is non-positive; refusing\n";
    return;
  }

  // weighted centroid x_c and the transport lever dRef = x_c - x_R
  Vector xc(ndm); xc.Zero();
  for (int i = 0; i < nIndep; i++) {
    const Vector& xi = theNodes[1 + i]->getCrds();
    for (int k = 0; k < ndm; k++) xc(k) += weights(i) * xi(k);
  }
  for (int k = 0; k < ndm; k++) xc(k) /= W;
  const Vector& xR = theNodes[0]->getCrds();
  dRef.resize(ndm);
  for (int k = 0; k < ndm; k++) dRef(k) = xc(k) - xR(k);

  // r_i = x_i - x_c, and the rotation scale ell2 = Sum w_i|r_i|^2 / W
  rvec.resize(nIndep, ndm);
  double Lc2 = 0.0, swr2 = 0.0;
  for (int i = 0; i < nIndep; i++) {
    const Vector& xi = theNodes[1 + i]->getCrds();
    double r2 = 0.0;
    for (int k = 0; k < ndm; k++) {
      double rk = xi(k) - xc(k);
      rvec(i, k) = rk;
      r2 += rk * rk;
    }
    swr2 += weights(i) * r2;
    if (r2 > Lc2) Lc2 = r2;
  }
  ell2 = swr2 / W;

  // weighted position-inertia I_c and its spectral pseudo-inverse / projector.
  Icplus.resize(nrot, nrot); Icplus.Zero();
  Pproj.resize(nrot, nrot);  Pproj.Zero();
  nKept = 0;

  if (ndm == 2) {
    // I_c is the scalar polar inertia Sum w_i|r_i|^2 (= ell2*W). One drilling axis.
    double Ic = swr2;
    double floor = 1.0e-8 * Ic;
    double absFloor = 1.0e-10 * W * Lc2;
    if (absFloor > floor) floor = absFloor;
    if (Ic > floor && Ic > 1.0e-300) {
      Icplus(0, 0) = 1.0 / Ic;
      Pproj(0, 0) = 1.0;
      nKept = 1;
    } else {
      opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
             << ": independent set has no in-plane spread; the reference drilling "
             << "rotation is left unconstrained\n";
    }
    valid = true;
    return;
  }

  // 3D: build I_c (3x3) = Sum w_i[(r.r)I - r(x)r], eigendecompose, drop degenerate axes.
  double Ic[3][3];
  for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) Ic[a][b] = 0.0;
  for (int i = 0; i < nIndep; i++) {
    double rx = rvec(i, 0), ry = rvec(i, 1), rz = rvec(i, 2);
    double r2 = rx * rx + ry * ry + rz * rz;
    double w = weights(i);
    Ic[0][0] += w * (r2 - rx * rx); Ic[0][1] += w * (-rx * ry); Ic[0][2] += w * (-rx * rz);
    Ic[1][1] += w * (r2 - ry * ry); Ic[1][2] += w * (-ry * rz);
    Ic[2][2] += w * (r2 - rz * rz);
  }
  Ic[1][0] = Ic[0][1]; Ic[2][0] = Ic[0][2]; Ic[2][1] = Ic[1][2];

  double eval[3], evec[3][3];
  jacobi3(Ic, eval, evec);
  double lmax = eval[0];
  for (int k = 1; k < 3; k++) if (eval[k] > lmax) lmax = eval[k];
  double floor = 1.0e-8 * lmax;
  double absFloor = 1.0e-10 * W * Lc2;
  if (absFloor > floor) floor = absFloor;

  for (int k = 0; k < 3; k++) {
    if (eval[k] > floor && eval[k] > 1.0e-300) {
      nKept++;
      double inv = 1.0 / eval[k];
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
          Icplus(a, b) += inv * evec[a][k] * evec[b][k];
          Pproj(a, b)  +=       evec[a][k] * evec[b][k];
        }
    } else {
      opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
             << ": reference rotation about axis (" << evec[0][k] << ", " << evec[1][k]
             << ", " << evec[2][k] << ") is unconstrained (degenerate independent set)\n";
    }
  }
  valid = true;
}

// cross operator C_i (nrot x ndm): (C_i u)_r = (r_i x u)_r.
double LadrunoDistributingCoupling::crossOp(int i, int r, int j) const
{
  if (ndm == 3) {
    double rx = rvec(i, 0), ry = rvec(i, 1), rz = rvec(i, 2);
    switch (r) {
    case 0: return (j == 1) ? -rz : (j == 2) ?  ry : 0.0;   // (r x u)_x = ry uz - rz uy
    case 1: return (j == 0) ?  rz : (j == 2) ? -rx : 0.0;   // (r x u)_y = rz ux - rx uz
    case 2: return (j == 0) ? -ry : (j == 1) ?  rx : 0.0;   // (r x u)_z = rx uy - ry ux
    default: return 0.0;
    }
  }
  // 2D drilling: (r x u)_z = rx uy - ry ux
  double rx = rvec(i, 0), ry = rvec(i, 1);
  return (j == 0) ? -ry : (j == 1) ? rx : 0.0;
}

// transport operator T (ndm x nrot): (T theta)_k = (theta x dRef)_k, dRef = x_c - x_R.
double LadrunoDistributingCoupling::transOp(int k, int r) const
{
  if (ndm == 3) {
    double dx = dRef(0), dy = dRef(1), dz = dRef(2);
    switch (k) {                                 // (theta x d)_k, d/dtheta_r
    case 0: return (r == 1) ?  dz : (r == 2) ? -dy : 0.0;   // (t x d)_x = ty dz - tz dy
    case 1: return (r == 0) ? -dz : (r == 2) ?  dx : 0.0;   // (t x d)_y = tz dx - tx dz
    case 2: return (r == 0) ?  dy : (r == 1) ? -dx : 0.0;   // (t x d)_z = tx dy - ty dx
    default: return 0.0;
    }
  }
  // 2D: theta (scalar) x d = theta * (-dy, dx)
  return (k == 0) ? -dRef(1) : dRef(0);
}

// build the constant gap operator B (ndm+nrot x nDOF): g = B u (then minus g0/gr0).
//   rows 0..ndm-1   = g_t:  u_R (+I) + T theta_R - Sum (w_i/W) u_i
//   rows ndm..+nrot = g_r:  P theta_R - I_c+ Sum w_i (r_i x u_i)
void LadrunoDistributingCoupling::buildB(void)
{
  B->Zero();
  int oR = dofOffset(0);
  for (int k = 0; k < ndm; k++) {
    (*B)(k, oR + k) += 1.0;                                    // u_R
    for (int r = 0; r < nrot; r++)
      (*B)(k, oR + ndm + r) += this->transOp(k, r);           // T theta_R (transport)
    for (int i = 0; i < nIndep; i++)
      (*B)(k, dofOffset(1 + i) + k) += -(weights(i) / W);     // -(w_i/W) u_i
  }
  for (int r = 0; r < nrot; r++) {
    int rowR = ndm + r;
    for (int r2 = 0; r2 < nrot; r2++)
      (*B)(rowR, oR + ndm + r2) += Pproj(r, r2);              // P theta_R
    for (int i = 0; i < nIndep; i++)
      for (int j = 0; j < ndm; j++) {
        double v = 0.0;                                        // (I_c+ C_i)(r,j)
        for (int s = 0; s < nrot; s++) v += Icplus(r, s) * this->crossOp(i, s, j);
        (*B)(rowR, dofOffset(1 + i) + j) += -weights(i) * v;  // -w_i I_c+ C_i u_i
      }
  }
}

// -k auto (needs a representative -host) + derive K_r = K_t * ell2 (unless -kr given). Lazy.
void LadrunoDistributingCoupling::resolveAutoKt(void)
{
  if (ktResolved) return;
  if (ktAuto) {
    Domain* theDomain = this->getDomain();
    Element* host = (theDomain != 0 && hostEleTag >= 0) ? theDomain->getElement(hostEleTag) : 0;
    if (host == 0) {
      opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
             << ": -k auto needs a valid -host element; keeping K_t=" << Kt << "\n";
    } else {
      double scale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
      if (scale > 0.0) Kt = kAlpha * scale;
      else
        opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
               << ": -k auto host gave a zero stiffness scale; keeping K_t=" << Kt << "\n";
    }
  }
  if (!krUser) Kr = Kt * ell2;   // derived rotational penalty (ADR 28 sec 5)
  ktResolved = true;
}

void LadrunoDistributingCoupling::resolveBipenalty(void)
{
  if (!bipenalty || bpResolved) return;
  this->resolveAutoKt();
  iPenalty = 0.0;
  double kEff = this->effectiveCouplingStiffness();
  if (kEff <= 0.0) {
    opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
           << ": -bipenalty got a zero coupling stiffness; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  if (bpMode == 0) {
    if (bpDt <= 0.0) {
      opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
             << ": -dtcr needs a positive step; m_p=0\n";
      mPenalty = 0.0; bpResolved = true; return;
    }
    mPenalty = LadrunoEmbedded::massPenaltyDtcr(kEff, bpDt);
    if (Kr > 0.0) iPenalty = LadrunoEmbedded::massPenaltyDtcr(Kr, bpDt);
    bpResolved = true;
    return;
  }
  // -wcap beta: omega_host from a representative host (needs -host).
  Domain* theDomain = this->getDomain();
  if (theDomain == 0) return;
  Element* host = (hostEleTag >= 0) ? theDomain->getElement(hostEleTag) : 0;
  if (host == 0) {
    opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
           << ": -wcap needs a valid -host for omega_host; use -dtcr; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  double kScale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
  double mScale = LadrunoEmbedded::maxAbsDiagonal(host->getMass());
  if (mScale <= 0.0 || kScale <= 0.0) {
    opserr << "WARNING LadrunoDistributingCoupling " << this->getTag()
           << ": -wcap host has no mass (or no stiffness); use -dtcr; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  double omega2_host = kScale / mScale;
  if (bpBeta <= 0.0) bpBeta = 2.0;
  mPenalty = LadrunoEmbedded::massPenaltyWcap(kEff, bpBeta, omega2_host);
  if (Kr > 0.0) iPenalty = LadrunoEmbedded::massPenaltyWcap(Kr, bpBeta, omega2_host);
  bpResolved = true;
}

int LadrunoDistributingCoupling::setRayleighDampingFactors(double, double, double, double)
{
  return 0;   // a pure penalty coupling carries no physical Rayleigh damping
}

// D ≡ 0. Overridden (with getRayleighDampingForces) so the base Element's lazy
// damping-matrix slot is never needed: the no-op setRayleighDampingFactors above
// leaves the base index = −1, so base getDamp would dereference theMatrices[-1] the
// first time a transient integrator forms the C-tangent (FE_Element::addCtoTang) or
// the residual damping force (addD_Force) — a hard crash. Returning a zeroed,
// element-owned matrix bypasses that path. (allocate() builds C0 sized nDOF.)
const Matrix& LadrunoDistributingCoupling::getDamp(void)
{
  C0->Zero();
  return *C0;
}

const Vector& LadrunoDistributingCoupling::getRayleighDampingForces(void)
{
  dampF->Zero();
  return *dampF;
}

double LadrunoDistributingCoupling::getExplicitCriticalTimeStep(void)
{
  if (!bipenalty) return -1.0;
  this->resolveBipenalty();
  double kEff = this->effectiveCouplingStiffness();
  double dt = LadrunoEmbedded::criticalTimeStep(mPenalty, kEff);   // translational class
  if (kEff > 0.0 && dt <= 0.0) return -1.0;
  if (Kr > 0.0) {                                                  // rotational class
    double dtR = LadrunoEmbedded::criticalTimeStep(iPenalty, Kr);
    if (dtR <= 0.0) return -1.0;
    if (dt <= 0.0 || dtR < dt) dt = dtR;
  }
  return dt;
}

// ===========================================================================
//  state
// ===========================================================================
int LadrunoDistributingCoupling::commitState(void)
{
  if (enforce == 1) {                       // augmented-Lagrangian Uzawa update
    this->resolveAutoKt();
    Vector g(ndm + nrot);
    this->computeGap(g);
    for (int k = 0; k < ndm; k++)  lambda(k)   += Kt * g(k);
    for (int r = 0; r < nrot; r++) lambda_r(r) += Kr * g(ndm + r);
  }
  return this->Element::commitState();
}

int LadrunoDistributingCoupling::revertToLastCommit(void) { return 0; }

int LadrunoDistributingCoupling::revertToStart(void)
{
  lambda.Zero();
  lambda_r.Zero();
  return 0;
}

int LadrunoDistributingCoupling::update(void)
{
  this->resolveAutoKt();
  return 0;
}

// full gap g = B u (then minus the captured offsets g0/gr0). g is sized ndm+nrot.
void LadrunoDistributingCoupling::computeGap(Vector& g)
{
  int nGap = ndm + nrot;
  if (g.Size() != nGap) g.resize(nGap);
  Vector uFull(nDOF); uFull.Zero();
  for (int p = 0; p < 1 + nIndep; p++) {
    const Vector& u = theNodes[p]->getTrialDisp();
    int ndf = nodeNdf(p);
    for (int d = 0; d < ndf; d++) uFull(dofOffset(p) + d) = u(d);
  }
  g.Zero();
  for (int row = 0; row < nGap; row++) {
    double s = 0.0;
    for (int c = 0; c < nDOF; c++) s += (*B)(row, c) * uFull(c);
    g(row) = s;
  }
  if (g0Computed) {
    for (int k = 0; k < ndm; k++)  g(k)       -= g0(k);
    for (int r = 0; r < nrot; r++) g(ndm + r) -= gr0(r);
  }
}

void LadrunoDistributingCoupling::captureInitialGap(void)
{
  if (!initGapCapture || g0Computed) return;
  Vector g(ndm + nrot);
  this->computeGap(g);                 // absolute gap (g0Computed still false)
  for (int k = 0; k < ndm; k++)  g0(k)  = g(k);
  for (int r = 0; r < nrot; r++) gr0(r) = g(ndm + r);
  g0Computed = true;
}

// ===========================================================================
//  matrices / forces
// ===========================================================================
const Vector& LadrunoDistributingCoupling::getResistingForce(void)
{
  this->resolveAutoKt();
  P->Zero();
  int nGap = ndm + nrot;
  Vector g(nGap);
  this->computeGap(g);
  Vector t(nGap);
  for (int k = 0; k < ndm; k++) {
    t(k) = Kt * g(k);
    if (enforce == 1) t(k) += lambda(k);
  }
  for (int r = 0; r < nrot; r++) {
    t(ndm + r) = Kr * g(ndm + r);
    if (enforce == 1) t(ndm + r) += lambda_r(r);
  }
  // P = B^T t
  for (int row = 0; row < nGap; row++) {
    double tr = t(row);
    if (tr == 0.0) continue;
    for (int c = 0; c < nDOF; c++) {
      double b = (*B)(row, c);
      if (b != 0.0) (*P)(c) += b * tr;
    }
  }
  return *P;
}

const Vector& LadrunoDistributingCoupling::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  if (bipenalty) {
    this->resolveBipenalty();
    if (theNodes[0] != 0) {
      const Vector& a = theNodes[0]->getTrialAccel();
      if (mPenalty != 0.0)
        for (int k = 0; k < ndm; k++) (*P)(dofOffset(0) + k) += mPenalty * a(k);
      if (iPenalty != 0.0)
        for (int r = 0; r < nrot; r++) (*P)(dofOffset(0) + ndm + r) += iPenalty * a(ndm + r);
    }
  }
  return *P;
}

const Matrix& LadrunoDistributingCoupling::getTangentStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  int nGap = ndm + nrot;
  for (int row = 0; row < nGap; row++) {
    double dval = (row < ndm) ? Kt : Kr;
    if (dval == 0.0) continue;
    for (int c1 = 0; c1 < nDOF; c1++) {
      double b1 = (*B)(row, c1);
      if (b1 == 0.0) continue;
      double db = dval * b1;
      for (int c2 = 0; c2 < nDOF; c2++) {
        double b2 = (*B)(row, c2);
        if (b2 != 0.0) (*K)(c1, c2) += db * b2;
      }
    }
  }
  return *K;
}

const Matrix& LadrunoDistributingCoupling::getInitialStiff(void)
{
  return this->getTangentStiff();   // B, K_t, K_r all state-independent
}

const Matrix& LadrunoDistributingCoupling::getMass(void)
{
  M0->Zero();
  if (bipenalty) {
    this->resolveBipenalty();
    for (int k = 0; k < ndm; k++)  (*M0)(dofOffset(0) + k, dofOffset(0) + k) = mPenalty;
    for (int r = 0; r < nrot; r++) (*M0)(dofOffset(0) + ndm + r, dofOffset(0) + ndm + r) = iPenalty;
  }
  return *M0;
}

// ===========================================================================
//  serialization (recompute the geometry/B on recv from coords at setDomain)
// ===========================================================================
int LadrunoDistributingCoupling::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  static Vector hdr(17);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nIndep;
  hdr(3) = Kt;
  hdr(4) = Kr;
  hdr(5) = krUser ? 1.0 : 0.0;
  hdr(6) = ktAuto ? 1.0 : 0.0;
  hdr(7) = kAlpha;
  hdr(8) = hostEleTag;
  hdr(9) = enforce;
  hdr(10) = bipenalty ? 1.0 : 0.0;
  hdr(11) = bpMode;
  hdr(12) = bpDt;
  hdr(13) = bpBeta;
  hdr(14) = initGapCapture ? 1.0 : 0.0;
  hdr(15) = g0Computed ? 1.0 : 0.0;
  hdr(16) = 1.0;                          // version
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoDistributingCoupling::sendSelf - header failed\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoDistributingCoupling::sendSelf - ID failed\n";
    return -1;
  }
  // payload: weights(N) + lambda(ndm) + lambda_r(nrot) + g0(ndm) + gr0(nrot)
  Vector payload(nIndep + 2 * ndm + 2 * nrot);
  int off = 0;
  for (int i = 0; i < nIndep; i++) payload(off + i) = weights(i);
  off += nIndep;
  for (int k = 0; k < ndm; k++)  payload(off + k) = (lambda.Size() == ndm) ? lambda(k) : 0.0;
  off += ndm;
  for (int r = 0; r < nrot; r++) payload(off + r) = (lambda_r.Size() == nrot) ? lambda_r(r) : 0.0;
  off += nrot;
  for (int k = 0; k < ndm; k++)  payload(off + k) = (g0.Size() == ndm) ? g0(k) : 0.0;
  off += ndm;
  for (int r = 0; r < nrot; r++) payload(off + r) = (gr0.Size() == nrot) ? gr0(r) : 0.0;
  off += nrot;
  if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoDistributingCoupling::sendSelf - payload failed\n";
    return -1;
  }
  return 0;
}

int LadrunoDistributingCoupling::recvSelf(int commitTag, Channel& theChannel,
                                          FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  static Vector hdr(17);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoDistributingCoupling::recvSelf - header failed\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nIndep = (int)hdr(2);
  Kt = hdr(3);
  Kr = hdr(4);
  krUser = (hdr(5) != 0.0);
  ktAuto = (hdr(6) != 0.0);
  kAlpha = hdr(7);
  hostEleTag = (int)hdr(8);
  enforce = (int)hdr(9);
  bipenalty = (hdr(10) != 0.0);
  bpMode = (int)hdr(11);
  bpDt = hdr(12);
  bpBeta = hdr(13);
  initGapCapture = (hdr(14) != 0.0);
  g0Computed = (hdr(15) != 0.0);
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  bpResolved = false;
  mPenalty = 0.0;
  iPenalty = 0.0;
  nDOF = 0;
  valid = false;

  connectedNodes.resize(1 + nIndep);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoDistributingCoupling::recvSelf - ID failed\n";
    return -1;
  }
  Vector payload(nIndep + 2 * ndm + 2 * nrot);
  if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoDistributingCoupling::recvSelf - payload failed\n";
    return -1;
  }
  int off = 0;
  weights.resize(nIndep);
  for (int i = 0; i < nIndep; i++) weights(i) = payload(off + i);
  off += nIndep;
  lambda.resize(ndm);
  for (int k = 0; k < ndm; k++) lambda(k) = payload(off + k);
  off += ndm;
  lambda_r.resize(nrot);
  for (int r = 0; r < nrot; r++) lambda_r(r) = payload(off + r);
  off += nrot;
  g0.resize(ndm);
  for (int k = 0; k < ndm; k++) g0(k) = payload(off + k);
  off += ndm;
  gr0.resize(nrot);
  for (int r = 0; r < nrot; r++) gr0(r) = payload(off + r);
  off += nrot;

  dRef.resize(ndm);
  Icplus.resize(nrot, nrot);
  Pproj.resize(nrot, nrot);
  nodeNdf.resize(1 + nIndep);
  dofOffset.resize(1 + nIndep);
  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node*[1 + nIndep];
  for (int i = 0; i < 1 + nIndep; i++) theNodes[i] = 0;
  // geometry (x_c, r_i, I_c+, P, ell2) and B are rebuilt in setDomain from coords.
  return 0;
}

void LadrunoDistributingCoupling::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoDistributingCoupling (RBE3), tag: " << this->getTag() << "\n";
  s << "  reference node: " << connectedNodes(0) << "  independent nodes: ";
  for (int i = 0; i < nIndep; i++) s << connectedNodes(1 + i) << " ";
  s << "\n  K_t: " << Kt << "  K_r: " << Kr << "  enforce: "
    << (enforce == 1 ? "al" : "penalty");
  s << "  resolvable-rotation-axes: " << nKept << "/" << nrot;
  s << (initGapCapture ? (g0Computed ? "  +initGap(captured)" : "  +initGap")
                       : "  +absolute(no initGap)");
  s << "\n";
}

// ===========================================================================
//  responses
// ===========================================================================
Response* LadrunoDistributingCoupling::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  int nGap = ndm + nrot;
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "couplingForce") == 0)
    return new ElementResponse(this, 1, Vector(nGap));
  if (strcmp(argv[0], "gap") == 0)
    return new ElementResponse(this, 2, Vector(nGap));
  if (strcmp(argv[0], "kt") == 0 || strcmp(argv[0], "k") == 0)
    return new ElementResponse(this, 3, 0.0);
  if (strcmp(argv[0], "kr") == 0)
    return new ElementResponse(this, 4, 0.0);
  if (strcmp(argv[0], "lambda") == 0 || strcmp(argv[0], "augLambda") == 0)
    return new ElementResponse(this, 5, Vector(ndm));
  if (strcmp(argv[0], "lambdaR") == 0 || strcmp(argv[0], "augLambdaR") == 0)
    return new ElementResponse(this, 6, Vector(nrot));
  if (strcmp(argv[0], "mpenalty") == 0 || strcmp(argv[0], "massPenalty") == 0)
    return new ElementResponse(this, 7, 0.0);
  if (strcmp(argv[0], "dtcr") == 0 || strcmp(argv[0], "dtCritical") == 0)
    return new ElementResponse(this, 8, 0.0);
  if (strcmp(argv[0], "nKept") == 0 || strcmp(argv[0], "rotationAxes") == 0)
    return new ElementResponse(this, 9, 0.0);
  return 0;
}

int LadrunoDistributingCoupling::getResponse(int responseID, Information& eleInfo)
{
  this->resolveAutoKt();
  int nGap = ndm + nrot;
  switch (responseID) {
  case 1: {
    Vector g(nGap); this->computeGap(g);
    Vector t(nGap);
    for (int k = 0; k < ndm; k++) {
      t(k) = Kt * g(k);
      if (enforce == 1) t(k) += lambda(k);
    }
    for (int r = 0; r < nrot; r++) {
      t(ndm + r) = Kr * g(ndm + r);
      if (enforce == 1) t(ndm + r) += lambda_r(r);
    }
    return eleInfo.setVector(t);
  }
  case 2: { Vector g(nGap); this->computeGap(g); return eleInfo.setVector(g); }
  case 3: return eleInfo.setDouble(Kt);
  case 4: return eleInfo.setDouble(Kr);
  case 5: return eleInfo.setVector(lambda);
  case 6: return eleInfo.setVector(lambda_r);
  case 7: { this->resolveBipenalty(); return eleInfo.setDouble(mPenalty); }
  case 8: { double dt = this->getExplicitCriticalTimeStep(); return eleInfo.setDouble(dt > 0.0 ? dt : 0.0); }
  case 9: return eleInfo.setDouble((double)nKept);
  default: return -1;
  }
}
