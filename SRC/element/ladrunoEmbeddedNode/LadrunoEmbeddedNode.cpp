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

// Ladruno: general node-to-host coupling element (ADR 23, Phase 1 = U). See
// LadrunoEmbeddedNode.h and Ladruno_implementation/23_ladruno_embedded_node_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoEmbeddedNode.h>
#include <LadrunoEmbeddedKernel.h>   // ADR 23 §D1 — shared coupling-kernel
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>        // ADR 23 Phase 2b (D9) — interface materials
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>

// ===========================================================================
//  construction
// ===========================================================================
LadrunoEmbeddedNode::LadrunoEmbeddedNode(int tag, int ndm_, int cNode,
        const ID& hostNodes, const Vector& shape, double ku,
        int hostEleTag_, bool ktAuto_, double ktAlpha_, int enforce_,
        bool bipenalty_, int bpMode_, double bpDt_, double bpBeta_,
        bool pressure_, double kp_,
        bool rot_, double kr_, bool krAuto_, double krAlpha_, const Matrix* gradN_,
        bool matMode_, const Vector* normalDir_, const Vector* orientDir_,
        UniaxialMaterial* matN_, UniaxialMaterial* matT1_, UniaxialMaterial* matT2_)
  : Element(tag, ELE_TAG_LadrunoEmbeddedNode),
    ndm(ndm_), nHost(hostNodes.Size()),
    connectedNodes(1 + hostNodes.Size()), Nshape(shape),
    Ku(ku), Kp(kp_), upActive(false), lambda_p(0.0),
    rflag(rot_ ? 1 : 0), nrot((ndm_ == 3) ? 3 : 1), rActive(false),
    gradN(), Kr(kr_), krAuto(krAuto_), krAlpha(krAlpha_), krResolved(false),
    lambda_r((ndm_ == 3) ? 3 : 1),
    matMode(matMode_ ? 1 : 0), normalDir(), orientDir(), haveOrient(false), frame(),
    enforce(enforce_), lambda(ndm_),
    hostEleTag(hostEleTag_), ktAuto(ktAuto_), ktAlpha(ktAlpha_), ktResolved(false),
    bipenalty(bipenalty_), bpMode(bpMode_), bpDt(bpDt_), bpBeta(bpBeta_),
    mPenalty(0.0), iPenalty(0.0), bpResolved(false),
    nDOF(0), nodeNdf(1 + hostNodes.Size()), dofOffset(1 + hostNodes.Size()),
    pflag(pressure_ ? 1 : 0),
    theNodes(0), K(0), P(0), M0(0)
{
  connectedNodes(0) = cNode;
  for (int i = 0; i < nHost; i++)
    connectedNodes(1 + i) = hostNodes(i);
  lambda.Zero();
  lambda_r.Zero();
  // host cartesian gradients ∂N_i/∂x_j (nHost × ndm) for the UR tie (-dNdx / -gradXi).
  if (rflag == 1 && gradN_ != 0)
    gradN = *gradN_;

  // ADR 23 Phase 2b (D9) — material-driven interface: clone each supplied uniaxial
  // (the element OWNS its copies), record the local frame, and build it.
  for (int d = 0; d < 3; d++) { matDir[d] = 0; hasMat[d] = false; }
  if (matMode == 1) {
    if (normalDir_ != 0) normalDir = *normalDir_;
    if (orientDir_ != 0) { orientDir = *orientDir_; haveOrient = true; }
    UniaxialMaterial* src[3] = { matN_, matT1_, matT2_ };
    for (int d = 0; d < ndm; d++)
      if (src[d] != 0) { matDir[d] = src[d]->getCopy(); hasMat[d] = (matDir[d] != 0); }
    this->buildFrame();
  }

  theNodes = new Node*[1 + nHost];
  for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
  // K/P/M0 are allocated in setDomain once the per-node ndf (⇒ nDOF) is known.
}

LadrunoEmbeddedNode::LadrunoEmbeddedNode()
  : Element(0, ELE_TAG_LadrunoEmbeddedNode),
    ndm(0), nHost(0), connectedNodes(0), Nshape(),
    Ku(0.0), Kp(0.0), upActive(false), lambda_p(0.0),
    rflag(0), nrot(0), rActive(false),
    gradN(), Kr(0.0), krAuto(false), krAlpha(0.0), krResolved(false), lambda_r(),
    matMode(0), normalDir(), orientDir(), haveOrient(false), frame(),
    enforce(0), lambda(),
    hostEleTag(-1), ktAuto(false), ktAlpha(0.0), ktResolved(false),
    bipenalty(false), bpMode(0), bpDt(0.0), bpBeta(0.0),
    mPenalty(0.0), iPenalty(0.0), bpResolved(false),
    nDOF(0), nodeNdf(), dofOffset(), pflag(0),
    theNodes(0), K(0), P(0), M0(0)
{
  for (int d = 0; d < 3; d++) { matDir[d] = 0; hasMat[d] = false; }
}

LadrunoEmbeddedNode::~LadrunoEmbeddedNode()
{
  if (theNodes != 0) delete[] theNodes;
  if (K != 0) delete K;
  if (P != 0) delete P;
  if (M0 != 0) delete M0;
  for (int d = 0; d < 3; d++) if (matDir[d] != 0) delete matDir[d];   // D9 — owned copies
}

// ===========================================================================
//  domain
// ===========================================================================
int LadrunoEmbeddedNode::getNumExternalNodes(void) const { return 1 + nHost; }
const ID& LadrunoEmbeddedNode::getExternalNodes(void) { return connectedNodes; }
Node** LadrunoEmbeddedNode::getNodePtrs(void) { return theNodes; }
int LadrunoEmbeddedNode::getNumDOF(void) { return nDOF; }

void LadrunoEmbeddedNode::allocate(void)
{
  if (K != 0) delete K;   K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;   P = new Vector(nDOF);
  if (M0 != 0) delete M0; M0 = new Matrix(nDOF, nDOF);
}

void LadrunoEmbeddedNode::setDomain(Domain* theDomain)
{
  if (theDomain == 0) {
    for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
    return;
  }
  // resolve the nodes and lay out the per-node DOF offsets. The constrained / host
  // nodes may carry more DOFs than ndm (u-p / u-r); Phase 1 couples the first ndm
  // (translational) DOFs of each (ADR 23 D8 — nDOF = Σ ndf_i).
  nodeNdf.resize(1 + nHost);
  dofOffset.resize(1 + nHost);
  int pos = 0;
  for (int i = 0; i < 1 + nHost; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "LadrunoEmbeddedNode::setDomain - node " << connectedNodes(i)
             << " not found\n";
      return;
    }
    int ndf = theNodes[i]->getNumberDOF();
    if (ndf < ndm) {
      opserr << "LadrunoEmbeddedNode::setDomain - node " << connectedNodes(i)
             << " has " << ndf << " DOFs, fewer than ndm=" << ndm << "\n";
      return;
    }
    nodeNdf(i) = ndf;
    dofOffset(i) = pos;
    pos += ndf;
  }
  nDOF = pos;

  // ADR 23 Phase 1b — resolve the pressure tie (UP). Activate only if -pressure was
  // requested AND every coupled node carries a pressure DOF at index ndm (ndf >= ndm+1,
  // the u-p convention); otherwise degrade to U-only with a warning.
  upActive = false;
  if (pflag == 1) {
    bool allUP = true;
    for (int i = 0; i < 1 + nHost; i++)
      if (nodeNdf(i) < ndm + 1) { allUP = false; break; }
    if (allUP) {
      upActive = true;
    } else {
      opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
             << ": -pressure requested but not all nodes are u-p (need ndf >= ndm+1="
             << ndm + 1 << "); coupling translations (U) only\n";
    }
  }

  // ADR 23 Phase 2 — resolve the rotation tie (UR). Activate only if -rot was set,
  // the constrained node carries the rotation DOFs (ndf ≥ ndm+nrot: 3D needs 6, 2D
  // needs 3), AND the host gradients gradN are present (nHost × ndm). Else degrade
  // to no-UR with a warning. (Host nodes need only translations, already checked.)
  rActive = false;
  if (rflag == 1) {
    if (nodeNdf(0) < ndm + nrot) {
      opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
             << ": -rot requested but the constrained node has " << nodeNdf(0)
             << " DOFs (need ndf >= ndm+nrot=" << ndm + nrot
             << " for rotations); coupling translations only\n";
    } else if (gradN.noRows() != nHost || gradN.noCols() != ndm) {
      opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
             << ": -rot requested but host gradients (∂N/∂x) are missing/mis-sized "
             << "(need " << nHost << "x" << ndm << ", got " << gradN.noRows() << "x"
             << gradN.noCols() << "); supply -dNdx or -gradXi; coupling translations only\n";
    } else {
      rActive = true;
    }
  }

  this->allocate();
  this->DomainComponent::setDomain(theDomain);
}

// ADR 23 D3 / M2 — resolve the auto translational penalty from the host element's
// own initial stiffness:  K_u = ktAlpha * max_i|K_host(i,i)|  (the rebar -kt auto
// mechanism: getInitialStiff max-abs-diagonal, NO lch). Lazy (first assembly).
void LadrunoEmbeddedNode::resolveAutoKt(void)
{
  if (!ktAuto || ktResolved) return;
  Domain* theDomain = this->getDomain();
  if (theDomain == 0) return;
  if (hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -k auto needs the -host form; keeping K_u=" << Ku << "\n";
    ktResolved = true;
    return;
  }
  Element* host = theDomain->getElement(hostEleTag);
  if (host == 0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -k auto host element " << hostEleTag
           << " not found; keeping K_u=" << Ku << "\n";
    ktResolved = true;
    return;
  }
  double scale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
  if (scale <= 0.0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -k auto host " << hostEleTag
           << " gave a zero stiffness scale; keeping K_u=" << Ku << "\n";
    ktResolved = true;
    return;
  }
  Ku = ktAlpha * scale;
  ktResolved = true;
}

// ADR 23 D3 — resolve the auto rotational penalty (after K_u). The rotation gap is
// a ROTATION (rad), so it needs a moment/rotation scale: K_r = krAlpha·K_u·lch²
// (lch from the host's getCharacteristicLength — the first such call in the embedded
// family, M2/HON-1). Dimensionally required so a unit rotation gap costs energy
// comparable to a unit translation gap. With krAuto off, K_r is the -kr numeric value.
void LadrunoEmbeddedNode::resolveKr(void)
{
  if (!rActive || krResolved) return;
  if (!krAuto) { krResolved = true; return; }   // numeric K_r already set
  this->resolveAutoKt();                         // resolve K_u first (may be auto)
  Domain* theDomain = this->getDomain();
  Element* host = (theDomain != 0 && hostEleTag >= 0)
                    ? theDomain->getElement(hostEleTag) : 0;
  if (host == 0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -kr auto needs the -host form for lch; keeping K_r=" << Kr << "\n";
    krResolved = true;
    return;
  }
  double lch = host->getCharacteristicLength();
  if (lch <= 0.0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -kr auto host " << hostEleTag << " gave lch=" << lch
           << "; keeping K_r=" << Kr << "\n";
    krResolved = true;
    return;
  }
  Kr = krAlpha * Ku * lch * lch;
  krResolved = true;
}

// host-rotation operator: ∂θ_host_r/∂u_i,j at the embedded point, from the host
// cartesian gradients g_i = ∇N_i. The continuum rotation θ = ½ curl(u) = skew(∇u):
//   3D: θ = ½ Σ_i (∇N_i × u_i)  ⇒  the block is ½·skew(∇N_i).
//   2D: θ_z = ½ Σ_i (∂N_i/∂x·u_y − ∂N_i/∂y·u_x)  (the single drilling rotation).
double LadrunoEmbeddedNode::rotOper(int r, int i, int j) const
{
  const double half = 0.5;
  if (ndm == 3) {
    double gx = gradN(i, 0), gy = gradN(i, 1), gz = gradN(i, 2);
    switch (r) {
    case 0: return (j == 1) ? -half * gz : (j == 2) ?  half * gy : 0.0;  // θ_x
    case 1: return (j == 0) ?  half * gz : (j == 2) ? -half * gx : 0.0;  // θ_y
    case 2: return (j == 0) ? -half * gy : (j == 1) ?  half * gx : 0.0;  // θ_z
    default: return 0.0;
    }
  }
  // 2D drilling
  double gx = gradN(i, 0), gy = gradN(i, 1);
  return (j == 0) ? -half * gy : (j == 1) ? half * gx : 0.0;
}

// g_r = θ_c − θ_host(ξ): the constrained node's rotation DOFs (indices ndm..ndm+nrot)
// minus the host continuum rotation interpolated at the embedded point (size nrot).
void LadrunoEmbeddedNode::computeGapR(Vector& gr)
{
  if (gr.Size() != nrot) gr.resize(nrot);
  const Vector& uc = theNodes[0]->getTrialDisp();
  for (int r = 0; r < nrot; r++) {
    double th = 0.0;
    for (int i = 0; i < nHost; i++) {
      const Vector& uh = theNodes[1 + i]->getTrialDisp();
      for (int j = 0; j < ndm; j++)
        th += this->rotOper(r, i, j) * uh(j);
    }
    gr(r) = uc(ndm + r) - th;
  }
}

// ADR 23 Phase 2b (D9) — build the orthonormal local frame from the supplied normal
// (e_0) and the optional -orient hint (else an auto reference); columns of `frame` are
// e_0 (normal), e_1, e_2 (tangents). 2D: e_1 = (−n_y, n_x). 3D: e_1 = Gram-Schmidt of
// the orient/auto reference against n, e_2 = n × e_1.
void LadrunoEmbeddedNode::buildFrame(void)
{
  frame.resize(ndm, ndm);
  frame.Zero();
  Vector n(ndm);
  double nn = 0.0;
  for (int k = 0; k < ndm; k++) {
    n(k) = (normalDir.Size() == ndm) ? normalDir(k) : 0.0;
    nn += n(k) * n(k);
  }
  if (nn <= 1.0e-20) { n.Zero(); n(0) = 1.0; nn = 1.0; }   // degenerate ⇒ global x
  nn = sqrt(nn);
  for (int k = 0; k < ndm; k++) { n(k) /= nn; frame(k, 0) = n(k); }

  if (ndm == 2) {
    frame(0, 1) = -n(1); frame(1, 1) = n(0);   // the one in-plane tangent
    return;
  }
  // 3D — first tangent from -orient (or an auto axis least aligned with n)
  Vector t1(3);
  if (haveOrient && orientDir.Size() == 3) {
    for (int k = 0; k < 3; k++) t1(k) = orientDir(k);
  } else {
    double ax = fabs(n(0)), ay = fabs(n(1)), az = fabs(n(2));
    t1.Zero();
    if (ax <= ay && ax <= az) t1(0) = 1.0; else if (ay <= az) t1(1) = 1.0; else t1(2) = 1.0;
  }
  double d = t1(0) * n(0) + t1(1) * n(1) + t1(2) * n(2);   // Gram-Schmidt ⟂ n
  for (int k = 0; k < 3; k++) t1(k) -= d * n(k);
  double t1n = sqrt(t1(0) * t1(0) + t1(1) * t1(1) + t1(2) * t1(2));
  if (t1n <= 1.0e-12) {                                    // orient ∥ n ⇒ fall back
    Vector ref(3); ref.Zero(); ref((fabs(n(0)) < 0.9) ? 0 : 1) = 1.0;
    d = ref(0) * n(0) + ref(1) * n(1) + ref(2) * n(2);
    for (int k = 0; k < 3; k++) t1(k) = ref(k) - d * n(k);
    t1n = sqrt(t1(0) * t1(0) + t1(1) * t1(1) + t1(2) * t1(2));
  }
  for (int k = 0; k < 3; k++) t1(k) /= t1n;
  Vector t2(3);                                            // t2 = n × t1
  t2(0) = n(1) * t1(2) - n(2) * t1(1);
  t2(1) = n(2) * t1(0) - n(0) * t1(2);
  t2(2) = n(0) * t1(1) - n(1) * t1(0);
  for (int k = 0; k < 3; k++) { frame(k, 1) = t1(k); frame(k, 2) = t2(k); }
}

// ADR 23 Phase 2b (D9) — translational traction t and tangent D for the gap g, either
// the isotropic penalty (matMode 0: t=K_u·g, D=K_u·I) or the material-frame interface
// (matMode 1): per local direction d, g_d=g·e_d, (t_d,k_d) from mat_d (FORCE units) or
// the penalty K_u fallback; t=Σ t_d e_d, D=Σ k_d e_d⊗e_d. (Excludes the AL λ — added by
// the caller.) Drives the materials' setTrialStrain (their state advances on commit).
void LadrunoEmbeddedNode::formTransTraction(const Vector& g, Vector& t, Matrix& D)
{
  if (t.Size() != ndm) t.resize(ndm);
  if (D.noRows() != ndm || D.noCols() != ndm) D.resize(ndm, ndm);
  t.Zero();
  D.Zero();
  if (matMode == 0) {
    for (int k = 0; k < ndm; k++) { t(k) = Ku * g(k); D(k, k) = Ku; }
    return;
  }
  for (int dd = 0; dd < ndm; dd++) {
    double gd = 0.0;
    for (int k = 0; k < ndm; k++) gd += g(k) * frame(k, dd);
    double td, kd;
    if (matDir[dd] != 0) {
      matDir[dd]->setTrialStrain(gd);
      td = matDir[dd]->getStress();
      kd = matDir[dd]->getTangent();
    } else {
      td = Ku * gd; kd = Ku;
    }
    for (int k = 0; k < ndm; k++) t(k) += td * frame(k, dd);
    for (int a = 0; a < ndm; a++)
      for (int b = 0; b < ndm; b++) D(a, b) += kd * frame(a, dd) * frame(b, dd);
  }
}

// k_eff for the bipenalty bound. Translational class: K_u, or (D9 material mode) the
// max over the active directions' INITIAL tangents (M6/ES-4 — a safe upper bound ONLY
// for non-stiffening materials; a stiffening re-contact tangent breaks the closed-form
// dt_cr, so stiffening materials on a bipenalty direction are unsupported — use -dtcr
// from the closed-contact stiffness, or -enforce al).
double LadrunoEmbeddedNode::effectiveCouplingStiffness(void)
{
  double k = fabs(Ku);
  if (matMode == 1)
    for (int d = 0; d < ndm; d++)
      if (matDir[d] != 0) {
        double ki = fabs(matDir[d]->getInitialTangent());
        if (ki > k) k = ki;
      }
  return LadrunoEmbedded::effectiveCouplingStiffness(k, 0.0);
}

// ADR 23 D5 / ADR 20 §10.6 — resolve the mass penalties (lazily, after K_u / K_r).
// PER-DOF-CLASS (stiffness, inertia) pairs (M1/ES-1): a translational-only m_p
// CANNOT bound the rotation mode (different units), so the rotation class gets its
// OWN inertia I_p via the SAME budget formula but with K_r. Both lumped on the
// constrained (slave) node only (getMass) ⇒ global mass stays diagonal. Budgets:
//   bpMode 0 (-dtcr dt): m_p = k_eff·(dt/2)²,   I_p = K_r·(dt/2)²    ⇒ dt_r = dt_u = dt.
//   bpMode 1 (-wcap β):  m_p = k_eff/(β·ω_host)², I_p = K_r/(β·ω_host)² ⇒ dt_r = dt_u.
// (lch² cancels — I_p and K_r both carry it — so the rotation mode self-bounds.)
void LadrunoEmbeddedNode::resolveBipenalty(void)
{
  if (!bipenalty || bpResolved) return;
  this->resolveAutoKt();
  if (rActive) this->resolveKr();
  iPenalty = 0.0;
  double kEff = this->effectiveCouplingStiffness();
  if (kEff <= 0.0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -bipenalty got a zero coupling stiffness; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  if (bpMode == 0) {
    if (bpDt <= 0.0) {
      opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
             << ": -dtcr needs a positive step; m_p=0\n";
      mPenalty = 0.0; bpResolved = true; return;
    }
    mPenalty = LadrunoEmbedded::massPenaltyDtcr(kEff, bpDt);
    if (rActive && Kr > 0.0)
      iPenalty = LadrunoEmbedded::massPenaltyDtcr(Kr, bpDt);
    bpResolved = true;
    return;
  }
  // bpMode 1 — -wcap β: ω_host from the host (needs -host).
  Domain* theDomain = this->getDomain();
  if (theDomain == 0) return;
  if (hostEleTag < 0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -wcap needs the -host form for ω_host; use -dtcr; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  Element* host = theDomain->getElement(hostEleTag);
  if (host == 0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -wcap host element " << hostEleTag << " not found; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  double kScale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
  double mScale = LadrunoEmbedded::maxAbsDiagonal(host->getMass());
  if (mScale <= 0.0 || kScale <= 0.0) {
    opserr << "WARNING LadrunoEmbeddedNode " << this->getTag()
           << ": -wcap host " << hostEleTag << " has no mass (or no stiffness); "
           << "use -dtcr instead; m_p=0\n";
    mPenalty = 0.0; bpResolved = true; return;
  }
  double omega2_host = kScale / mScale;
  if (bpBeta <= 0.0) bpBeta = 2.0;
  mPenalty = LadrunoEmbedded::massPenaltyWcap(kEff, bpBeta, omega2_host);
  if (rActive && Kr > 0.0)
    iPenalty = LadrunoEmbedded::massPenaltyWcap(Kr, bpBeta, omega2_host);
  bpResolved = true;
}

int LadrunoEmbeddedNode::setRayleighDampingFactors(double, double, double, double)
{
  return 0;   // a pure penalty coupling carries no physical Rayleigh damping
}

double LadrunoEmbeddedNode::getExplicitCriticalTimeStep(void)
{
  if (!bipenalty) return -1.0;
  this->resolveBipenalty();
  // min over active DOF classes of 2√(m_class/k_class) (ADR 23 D5 / M1·ES-1): the
  // translational (m_p, k_eff) and, when UR is on, the rotational (I_p, K_r) pair.
  double dt = LadrunoEmbedded::criticalTimeStep(mPenalty, this->effectiveCouplingStiffness());
  if (rActive) {
    double dtR = LadrunoEmbedded::criticalTimeStep(iPenalty, Kr);
    if (dtR > 0.0 && (dt <= 0.0 || dtR < dt)) dt = dtR;
  }
  return dt;
}

// ===========================================================================
//  state
// ===========================================================================
int LadrunoEmbeddedNode::commitState(void)
{
  // Augmented-Lagrangian per-step Uzawa update (ADR 23 D4): at the converged state,
  // accumulate the PENALTY traction into the multiplier.
  if (enforce == 1) {
    this->resolveAutoKt();
    Vector g(ndm);
    this->computeGap(g);
    if (matMode == 0) {
      // ISOTROPIC tie ⇒ NO directional re-projection of λ (no preferred axis; ES-3).
      for (int k = 0; k < ndm; k++) lambda(k) += Ku * g(k);
    } else {
      // D9 — Uzawa over the PENALTY directions only (material directions carry physical
      // force, not λ); then M4/D9-3 re-project λ off every material direction so an AL
      // step can never leak the multiplier onto a material axis.
      for (int dd = 0; dd < ndm; dd++) {
        if (hasMat[dd]) continue;
        double gd = 0.0;
        for (int k = 0; k < ndm; k++) gd += g(k) * frame(k, dd);
        for (int k = 0; k < ndm; k++) lambda(k) += Ku * gd * frame(k, dd);
      }
      for (int dd = 0; dd < ndm; dd++) {
        if (!hasMat[dd]) continue;
        double la = 0.0;
        for (int k = 0; k < ndm; k++) la += lambda(k) * frame(k, dd);
        for (int k = 0; k < ndm; k++) lambda(k) -= la * frame(k, dd);
      }
    }
    if (upActive) lambda_p += Kp * this->computeGapP();
    if (rActive) {
      this->resolveKr();
      Vector gr(nrot);
      this->computeGapR(gr);
      for (int r = 0; r < nrot; r++) lambda_r(r) += Kr * gr(r);
    }
  }
  int ok = this->Element::commitState();
  for (int d = 0; d < 3; d++) if (matDir[d] != 0) ok += matDir[d]->commitState();   // D9
  return ok;
}

int LadrunoEmbeddedNode::revertToLastCommit(void)
{
  int ok = 0;
  for (int d = 0; d < 3; d++) if (matDir[d] != 0) ok += matDir[d]->revertToLastCommit();
  return ok;
}

int LadrunoEmbeddedNode::revertToStart(void)
{
  lambda.Zero();   // AL multipliers accumulate only via commitState's Uzawa step
  lambda_p = 0.0;
  lambda_r.Zero();
  int ok = 0;
  for (int d = 0; d < 3; d++) if (matDir[d] != 0) ok += matDir[d]->revertToStart();
  return ok;
}

void LadrunoEmbeddedNode::computeGap(Vector& g)
{
  // g = u_c − Σ N_i u_host,i (translational, the first ndm DOFs of each node).
  LadrunoEmbedded::computeGap(theNodes, Nshape, nHost, ndm, g);
}

double LadrunoEmbeddedNode::computeGapP(void)
{
  // g_p = p_c − Σ N_i p_host,i  (pressure DOF = index ndm, the u-p convention).
  double gp = theNodes[0]->getTrialDisp()(ndm);
  for (int i = 0; i < nHost; i++)
    gp -= Nshape(i) * theNodes[1 + i]->getTrialDisp()(ndm);
  return gp;
}

int LadrunoEmbeddedNode::update(void)
{
  this->resolveAutoKt();
  if (rActive) this->resolveKr();
  return 0;
}

// ===========================================================================
//  matrices / forces — compact ndm-packed coupling SCATTERED to the full layout
// ===========================================================================
const Vector& LadrunoEmbeddedNode::getResistingForce(void)
{
  this->resolveAutoKt();
  P->Zero();
  Vector g(ndm);
  this->computeGap(g);

  // translational traction: isotropic K_u·g, or (D9) the material-frame interface
  // t = Σ_d f_d(g·e_d) e_d. (+ λ for AL — constant within an inner solve.)
  Vector t(ndm);
  Matrix Dscratch(ndm, ndm);
  this->formTransTraction(g, t, Dscratch);
  if (enforce == 1)
    for (int k = 0; k < ndm; k++) t(k) += lambda(k);

  // compact P = Bᵀt over the (1+M)·ndm translational space (shared kernel), then
  // scatter each node's ndm-block to its global DOF offset.
  int nc = (1 + nHost) * ndm;
  Vector Pc(nc);
  LadrunoEmbedded::assembleResistingForce(Pc, Nshape, t, ndm);
  for (int p = 0; p < 1 + nHost; p++)
    for (int k = 0; k < ndm; k++)
      (*P)(dofOffset(p) + k) = Pc(p * ndm + k);

  // ADR 23 Phase 1b — pressure tie (UP): a scalar coupling on the pressure DOF
  // (index ndm). t_p = K_p·g_p (+ λ_p for AL); same Bᵀt structure with ndm=1,
  // scattered to each node's pressure slot dofOffset(p)+ndm.
  if (upActive) {
    double tp = Kp * this->computeGapP();
    if (enforce == 1) tp += lambda_p;
    Vector tpv(1); tpv(0) = tp;
    Vector Pp(1 + nHost);
    LadrunoEmbedded::assembleResistingForce(Pp, Nshape, tpv, 1);
    for (int p = 0; p < 1 + nHost; p++)
      (*P)(dofOffset(p) + ndm) = Pp(p);
  }

  // ADR 23 Phase 2 — rotation tie (UR): t_r = K_r·g_r (+ λ_r for AL). P_r = B_rᵀ t_r
  // with B_r: cNode rotations → +I, host translations → −∂θ_host/∂u (rotOper). The
  // cNode rotation slots (ndm..ndm+nrot−1) are disjoint from U/UP (UR ⊥ UP by parse
  // guard), so += writes a zeroed slot; host translation slots get the rotation
  // contribution ON TOP of the U force.
  if (rActive) {
    this->resolveKr();
    Vector gr(nrot);
    this->computeGapR(gr);
    Vector tr(nrot);
    for (int r = 0; r < nrot; r++) {
      tr(r) = Kr * gr(r);
      if (enforce == 1) tr(r) += lambda_r(r);
    }
    for (int r = 0; r < nrot; r++)
      (*P)(dofOffset(0) + ndm + r) += tr(r);          // cNode rotation DOFs
    for (int i = 0; i < nHost; i++)
      for (int j = 0; j < ndm; j++) {
        double f = 0.0;
        for (int r = 0; r < nrot; r++) f += this->rotOper(r, i, j) * tr(r);
        (*P)(dofOffset(1 + i) + j) += -f;             // host translation DOFs
      }
  }
  return *P;
}

const Vector& LadrunoEmbeddedNode::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  // ADR 23 D5 — bipenalty inertia on the cNode only: m_p·a on the translational DOFs,
  // and (UR) I_p·α on the rotation DOFs (the per-DOF-class inertia, M1/ES-1).
  if (bipenalty) {
    this->resolveBipenalty();
    if (theNodes[0] != 0) {
      const Vector& aC = theNodes[0]->getTrialAccel();
      if (mPenalty != 0.0)
        for (int k = 0; k < ndm; k++) (*P)(dofOffset(0) + k) += mPenalty * aC(k);
      if (rActive && iPenalty != 0.0)
        for (int r = 0; r < nrot; r++)
          (*P)(dofOffset(0) + ndm + r) += iPenalty * aC(ndm + r);
    }
  }
  return *P;
}

const Matrix& LadrunoEmbeddedNode::getTangentStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  // translational tangent D_u: isotropic K_u·I, or (D9) the material-frame
  // D = Σ_d k_d e_d⊗e_d (state-dependent ⇒ needs the current gap).
  Vector g(ndm);
  this->computeGap(g);
  Vector tscratch(ndm);
  Matrix Du(ndm, ndm);
  this->formTransTraction(g, tscratch, Du);
  // compact K = Bᵀ D_u B over the (1+M)·ndm space, scattered to the full layout.
  int nc = (1 + nHost) * ndm;
  Matrix Kc(nc, nc);
  Kc.Zero();
  LadrunoEmbedded::assembleTangent(Kc, Nshape, Du, ndm);
  for (int p = 0; p < 1 + nHost; p++)
    for (int q = 0; q < 1 + nHost; q++)
      for (int a = 0; a < ndm; a++)
        for (int b = 0; b < ndm; b++)
          (*K)(dofOffset(p) + a, dofOffset(q) + b) += Kc(p * ndm + a, q * ndm + b);

  // ADR 23 Phase 1b — pressure block (UP): scalar D_p = K_p, assembled with ndm=1
  // and scattered to the pressure DOF slots (dofOffset(p)+ndm). No overlap with the
  // translational blocks (a,b < ndm).
  if (upActive) {
    Matrix Dp(1, 1); Dp(0, 0) = Kp;
    Matrix Kp_c(1 + nHost, 1 + nHost); Kp_c.Zero();
    LadrunoEmbedded::assembleTangent(Kp_c, Nshape, Dp, 1);
    for (int p = 0; p < 1 + nHost; p++)
      for (int q = 0; q < 1 + nHost; q++)
        (*K)(dofOffset(p) + ndm, dofOffset(q) + ndm) += Kp_c(p, q);
  }

  // ADR 23 Phase 2 — rotation block (UR): K_r = B_rᵀ D_r B_r, D_r = K_r·I_{nrot}.
  // cNode-rot/cNode-rot = K_r·I (diagonal); cNode-rot/host = −K_r·rotOper; host/host
  // = K_r·Σ_r rotOper·rotOper. Rotation/translation slots are disjoint (UR ⊥ UP), so
  // these ADD onto the (already-zeroed or U-populated) full layout. State-independent
  // (K_r, rotOper constant) ⇒ getInitialStiff (=getTangentStiff) is exact for UR too.
  if (rActive) {
    this->resolveKr();
    for (int r = 0; r < nrot; r++) {
      int cr = dofOffset(0) + ndm + r;
      (*K)(cr, cr) += Kr;                                       // cNode-rot diagonal
      for (int i = 0; i < nHost; i++)
        for (int j = 0; j < ndm; j++) {
          int hj = dofOffset(1 + i) + j;
          double v = -Kr * this->rotOper(r, i, j);             // cNode-rot × host
          (*K)(cr, hj) += v;
          (*K)(hj, cr) += v;
        }
    }
    for (int i = 0; i < nHost; i++)                             // host × host
      for (int j = 0; j < ndm; j++) {
        int hj = dofOffset(1 + i) + j;
        for (int i2 = 0; i2 < nHost; i2++)
          for (int j2 = 0; j2 < ndm; j2++) {
            double s = 0.0;
            for (int r = 0; r < nrot; r++)
              s += this->rotOper(r, i, j) * this->rotOper(r, i2, j2);
            (*K)(hj, dofOffset(1 + i2) + j2) += Kr * s;
          }
      }
  }
  return *K;
}

const Matrix& LadrunoEmbeddedNode::getInitialStiff(void)
{
  // U coupling is state-independent ⇒ identical to the tangent.
  return this->getTangentStiff();
}

const Matrix& LadrunoEmbeddedNode::getMass(void)
{
  M0->Zero();
  // ADR 23 D5 — bipenalty: lump m_p on the cNode's first ndm (translational) DOFs
  // (cNode is node 0, offset 0). Host blocks stay zero ⇒ M diagonal, DiagonalSOE-safe.
  if (bipenalty) {
    this->resolveBipenalty();
    for (int k = 0; k < ndm; k++) (*M0)(dofOffset(0) + k, dofOffset(0) + k) = mPenalty;
    // (UR) rotational inertia I_p on the cNode rotation DOFs (M1/ES-1 — bounds the
    // rotation mode, which the translational m_p cannot). Lumped on the slave only.
    if (rActive)
      for (int r = 0; r < nrot; r++)
        (*M0)(dofOffset(0) + ndm + r, dofOffset(0) + ndm + r) = iPenalty;
  }
  return *M0;
}

// ===========================================================================
//  serialization
// ===========================================================================
int LadrunoEmbeddedNode::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  static Vector hdr(26);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nHost;
  hdr(3) = Ku;
  hdr(4) = hostEleTag;
  hdr(5) = ktAuto ? 1.0 : 0.0;
  hdr(6) = ktAlpha;
  hdr(7) = enforce;
  hdr(8) = bipenalty ? 1.0 : 0.0;
  hdr(9) = bpMode;
  hdr(10) = bpDt;
  hdr(11) = bpBeta;
  hdr(12) = pflag;
  hdr(13) = Kp;
  hdr(14) = rflag;                    // ADR 23 Phase 2 — UR
  hdr(15) = Kr;
  hdr(16) = krAuto ? 1.0 : 0.0;
  hdr(17) = krAlpha;
  hdr(18) = matMode;                  // ADR 23 Phase 2b — D9 material interface
  hdr(19) = haveOrient ? 1.0 : 0.0;
  for (int d = 0; d < 3; d++) {       // per-direction material classTag + dbTag
    int ct = 0, dt = 0;
    if (matDir[d] != 0) {
      ct = matDir[d]->getClassTag();
      dt = matDir[d]->getDbTag();
      if (dt == 0) { dt = theChannel.getDbTag(); matDir[d]->setDbTag(dt); }
    }
    hdr(20 + d) = ct;
    hdr(23 + d) = dt;
  }
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - header failed\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - ID failed\n";
    return -1;
  }
  // payload: Nshape(nHost) + lambda(ndm) + lambda_p(1) [+ gradN(nHost·ndm)+lambda_r(nrot)
  // if UR] [+ normalDir(ndm)+orientDir(ndm) if D9].
  int extraR = (rflag == 1) ? (nHost * ndm + nrot) : 0;
  int extraM = (matMode == 1) ? (2 * ndm) : 0;
  Vector payload(nHost + ndm + 1 + extraR + extraM);
  for (int i = 0; i < nHost; i++) payload(i) = Nshape(i);
  for (int k = 0; k < ndm; k++)
    payload(nHost + k) = (lambda.Size() == ndm) ? lambda(k) : 0.0;
  payload(nHost + ndm) = lambda_p;
  int off = nHost + ndm + 1;
  if (rflag == 1) {
    for (int i = 0; i < nHost; i++)
      for (int j = 0; j < ndm; j++)
        payload(off + i * ndm + j) =
          (gradN.noRows() == nHost && gradN.noCols() == ndm) ? gradN(i, j) : 0.0;
    off += nHost * ndm;
    for (int r = 0; r < nrot; r++)
      payload(off + r) = (lambda_r.Size() == nrot) ? lambda_r(r) : 0.0;
    off += nrot;
  }
  if (matMode == 1) {
    for (int k = 0; k < ndm; k++)
      payload(off + k) = (normalDir.Size() == ndm) ? normalDir(k) : 0.0;
    for (int k = 0; k < ndm; k++)
      payload(off + ndm + k) = (haveOrient && orientDir.Size() == ndm) ? orientDir(k) : 0.0;
  }
  if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - payload failed\n";
    return -1;
  }
  // D9 — each interface material serializes itself after the element payload.
  for (int d = 0; d < 3; d++)
    if (matDir[d] != 0 && matDir[d]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "LadrunoEmbeddedNode::sendSelf - interface material " << d << " failed\n";
      return -1;
    }
  return 0;
}

int LadrunoEmbeddedNode::recvSelf(int commitTag, Channel& theChannel,
                                  FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  static Vector hdr(26);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoEmbeddedNode::recvSelf - header failed\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nHost = (int)hdr(2);
  Ku = hdr(3);
  hostEleTag = (int)hdr(4);
  ktAuto = (hdr(5) != 0.0);
  ktAlpha = hdr(6);
  enforce = (int)hdr(7);
  bipenalty = (hdr(8) != 0.0);
  bpMode = (int)hdr(9);
  bpDt = hdr(10);
  bpBeta = hdr(11);
  pflag = (int)hdr(12);
  Kp = hdr(13);
  rflag = (int)hdr(14);                  // ADR 23 Phase 2 — UR
  Kr = hdr(15);
  krAuto = (hdr(16) != 0.0);
  krAlpha = hdr(17);
  matMode = (int)hdr(18);                // ADR 23 Phase 2b — D9
  haveOrient = (hdr(19) != 0.0);
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  bpResolved = false;
  krResolved = false;
  mPenalty = 0.0;
  iPenalty = 0.0;
  upActive = false;   // re-resolved in setDomain (needs the nodes' ndf)
  rActive = false;    // re-resolved in setDomain
  nDOF = 0;           // recomputed in setDomain

  connectedNodes.resize(1 + nHost);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedNode::recvSelf - ID failed\n";
    return -1;
  }
  int extraR = (rflag == 1) ? (nHost * ndm + nrot) : 0;
  int extraM = (matMode == 1) ? (2 * ndm) : 0;
  Vector payload(nHost + ndm + 1 + extraR + extraM);
  if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedNode::recvSelf - payload failed\n";
    return -1;
  }
  Nshape.resize(nHost);
  for (int i = 0; i < nHost; i++) Nshape(i) = payload(i);
  lambda.resize(ndm);
  for (int k = 0; k < ndm; k++) lambda(k) = payload(nHost + k);
  lambda_p = payload(nHost + ndm);
  lambda_r.resize(nrot); lambda_r.Zero();
  int off = nHost + ndm + 1;
  if (rflag == 1) {
    gradN.resize(nHost, ndm);
    for (int i = 0; i < nHost; i++)
      for (int j = 0; j < ndm; j++)
        gradN(i, j) = payload(off + i * ndm + j);
    off += nHost * ndm;
    for (int r = 0; r < nrot; r++) lambda_r(r) = payload(off + r);
    off += nrot;
  }
  // D9 — restore the local frame inputs + reconstruct each interface material via the
  // broker, then rebuild the orthonormal frame.
  for (int d = 0; d < 3; d++) { if (matDir[d] != 0) delete matDir[d]; matDir[d] = 0; hasMat[d] = false; }
  if (matMode == 1) {
    normalDir.resize(ndm);
    for (int k = 0; k < ndm; k++) normalDir(k) = payload(off + k);
    orientDir.resize(ndm);
    for (int k = 0; k < ndm; k++) orientDir(k) = payload(off + ndm + k);
    for (int d = 0; d < 3; d++) {
      int ct = (int)hdr(20 + d), dt = (int)hdr(23 + d);
      if (ct != 0) {
        matDir[d] = theBroker.getNewUniaxialMaterial(ct);
        if (matDir[d] == 0) {
          opserr << "LadrunoEmbeddedNode::recvSelf - broker has no material classTag "
                 << ct << " (dir " << d << ")\n";
          return -1;
        }
        matDir[d]->setDbTag(dt);
        if (matDir[d]->recvSelf(commitTag, theChannel, theBroker) < 0) {
          opserr << "LadrunoEmbeddedNode::recvSelf - interface material " << d
                 << " recv failed\n";
          return -1;
        }
        hasMat[d] = true;
      }
    }
    this->buildFrame();
  }

  nodeNdf.resize(1 + nHost);
  dofOffset.resize(1 + nHost);
  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node*[1 + nHost];
  for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
  // K/P/M0 (re)allocated in setDomain once nDOF is known.
  return 0;
}

void LadrunoEmbeddedNode::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoEmbeddedNode, tag: " << this->getTag() << "\n";
  s << "  constrained node: " << connectedNodes(0) << "  host nodes: ";
  for (int i = 0; i < nHost; i++) s << connectedNodes(1 + i) << " ";
  s << "\n  K_u: " << Ku << "  enforce: " << (enforce == 1 ? "al" : "penalty");
  if (pflag == 1) s << "  +pressure(K_p=" << Kp << (upActive ? ", active" : ", inactive") << ")";
  if (rflag == 1) s << "  +rotation(K_r=" << Kr << (rActive ? ", active" : ", inactive") << ")";
  if (matMode == 1) {
    s << "  +material-interface(";
    for (int d = 0; d < ndm; d++)
      s << (matDir[d] != 0 ? "mat" : "penalty") << (d < ndm - 1 ? "/" : "");
    s << ")";
  }
  s << "\n";
}

// ===========================================================================
//  responses
// ===========================================================================
Response* LadrunoEmbeddedNode::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "localForce") == 0)
    return new ElementResponse(this, 1, Vector(ndm));
  if (strcmp(argv[0], "gap") == 0)
    return new ElementResponse(this, 2, Vector(ndm));
  if (strcmp(argv[0], "kt") == 0 || strcmp(argv[0], "penalty") == 0 ||
      strcmp(argv[0], "k") == 0)
    return new ElementResponse(this, 3, 0.0);
  if (strcmp(argv[0], "penaltyEnergy") == 0)
    return new ElementResponse(this, 4, 0.0);
  if (strcmp(argv[0], "constraintViolation") == 0)
    return new ElementResponse(this, 5, 0.0);
  if (strcmp(argv[0], "augLambda") == 0 || strcmp(argv[0], "lambda") == 0)
    return new ElementResponse(this, 6, Vector(ndm));
  if (strcmp(argv[0], "mpenalty") == 0 || strcmp(argv[0], "massPenalty") == 0)
    return new ElementResponse(this, 7, 0.0);
  if (strcmp(argv[0], "dtcr") == 0 || strcmp(argv[0], "dtCritical") == 0)
    return new ElementResponse(this, 8, 0.0);
  // ADR 23 Phase 1b — pressure-tie (UP) diagnostics.
  if (strcmp(argv[0], "pgap") == 0 || strcmp(argv[0], "pressureGap") == 0)
    return new ElementResponse(this, 9, 0.0);
  if (strcmp(argv[0], "pforce") == 0 || strcmp(argv[0], "pressureForce") == 0)
    return new ElementResponse(this, 10, 0.0);
  if (strcmp(argv[0], "plambda") == 0 || strcmp(argv[0], "augLambdaP") == 0)
    return new ElementResponse(this, 11, 0.0);
  // ADR 23 Phase 2 — rotation-tie (UR) diagnostics (size nrot: 3 in 3D, 1 in 2D).
  if (strcmp(argv[0], "rgap") == 0 || strcmp(argv[0], "rotationGap") == 0)
    return new ElementResponse(this, 12, Vector(nrot > 0 ? nrot : 1));
  if (strcmp(argv[0], "rforce") == 0 || strcmp(argv[0], "moment") == 0)
    return new ElementResponse(this, 13, Vector(nrot > 0 ? nrot : 1));
  if (strcmp(argv[0], "rlambda") == 0 || strcmp(argv[0], "augLambdaR") == 0)
    return new ElementResponse(this, 14, Vector(nrot > 0 ? nrot : 1));
  if (strcmp(argv[0], "kr") == 0 || strcmp(argv[0], "rotPenalty") == 0)
    return new ElementResponse(this, 15, 0.0);
  // ADR 23 Phase 2b — D9 material interface diagnostics (local-frame components).
  if (strcmp(argv[0], "localGap") == 0 || strcmp(argv[0], "frameGap") == 0)
    return new ElementResponse(this, 16, Vector(ndm));
  if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "frameForce") == 0 ||
      strcmp(argv[0], "interfaceForce") == 0)
    return new ElementResponse(this, 17, Vector(ndm));
  if (strcmp(argv[0], "normal") == 0)
    return new ElementResponse(this, 18, Vector(ndm));
  return 0;
}

int LadrunoEmbeddedNode::getResponse(int responseID, Information& eleInfo)
{
  this->resolveAutoKt();
  Vector g(ndm);
  this->computeGap(g);
  switch (responseID) {
  case 1: {
    Vector t(ndm);
    Matrix Dscr(ndm, ndm);
    this->formTransTraction(g, t, Dscr);   // isotropic or D9 material-frame
    if (enforce == 1)
      for (int k = 0; k < ndm; k++) t(k) += lambda(k);
    return eleInfo.setVector(t);
  }
  case 2: return eleInfo.setVector(g);
  case 3: return eleInfo.setDouble(Ku);
  case 4: {
    double g2 = 0.0;
    for (int k = 0; k < ndm; k++) g2 += g(k) * g(k);
    return eleInfo.setDouble(0.5 * Ku * g2);
  }
  case 5: {
    double g2 = 0.0;
    for (int k = 0; k < ndm; k++) g2 += g(k) * g(k);
    return eleInfo.setDouble(sqrt(g2));
  }
  case 6: return eleInfo.setVector(lambda);
  case 7: { this->resolveBipenalty(); return eleInfo.setDouble(mPenalty); }
  case 8: {
    double dt = this->getExplicitCriticalTimeStep();
    return eleInfo.setDouble(dt > 0.0 ? dt : 0.0);
  }
  case 9:  return eleInfo.setDouble(upActive ? this->computeGapP() : 0.0);
  case 10: {
    if (!upActive) return eleInfo.setDouble(0.0);
    double tp = Kp * this->computeGapP();
    if (enforce == 1) tp += lambda_p;
    return eleInfo.setDouble(tp);
  }
  case 11: return eleInfo.setDouble(lambda_p);
  case 12: {   // rotation gap g_r = θ_c − θ_host (size nrot)
    Vector gr(nrot);
    if (rActive) this->computeGapR(gr); else gr.Zero();
    return eleInfo.setVector(gr);
  }
  case 13: {   // rotation traction (moment) t_r = K_r·g_r (+ λ_r)
    Vector tr(nrot); tr.Zero();
    if (rActive) {
      this->resolveKr();
      Vector gr(nrot); this->computeGapR(gr);
      for (int r = 0; r < nrot; r++) {
        tr(r) = Kr * gr(r);
        if (enforce == 1) tr(r) += lambda_r(r);
      }
    }
    return eleInfo.setVector(tr);
  }
  case 14: return eleInfo.setVector(lambda_r);
  case 15: { if (rActive) this->resolveKr(); return eleInfo.setDouble(Kr); }
  case 16: {   // D9 local-frame gap components g_d = g·e_d (global g if no frame)
    Vector gl(ndm); gl.Zero();
    if (matMode == 1)
      for (int d = 0; d < ndm; d++)
        for (int k = 0; k < ndm; k++) gl(d) += g(k) * frame(k, d);
    else
      gl = g;
    return eleInfo.setVector(gl);
  }
  case 17: {   // D9 local-frame force components t_d (per direction)
    Vector tl(ndm); tl.Zero();
    if (matMode == 1) {
      for (int d = 0; d < ndm; d++) {
        double gd = 0.0;
        for (int k = 0; k < ndm; k++) gd += g(k) * frame(k, d);
        if (matDir[d] != 0) { matDir[d]->setTrialStrain(gd); tl(d) = matDir[d]->getStress(); }
        else                  tl(d) = Ku * gd;
      }
    } else {
      for (int k = 0; k < ndm; k++) tl(k) = Ku * g(k);
    }
    return eleInfo.setVector(tl);
  }
  case 18: {   // D9 frame normal e_0 (global components)
    Vector nrm(ndm); nrm.Zero();
    if (matMode == 1) for (int k = 0; k < ndm; k++) nrm(k) = frame(k, 0);
    return eleInfo.setVector(nrm);
  }
  default: return -1;
  }
}
