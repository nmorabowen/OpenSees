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
        bool pressure_, double kp_)
  : Element(tag, ELE_TAG_LadrunoEmbeddedNode),
    ndm(ndm_), nHost(hostNodes.Size()),
    connectedNodes(1 + hostNodes.Size()), Nshape(shape),
    Ku(ku), Kp(kp_), upActive(false), lambda_p(0.0),
    enforce(enforce_), lambda(ndm_),
    hostEleTag(hostEleTag_), ktAuto(ktAuto_), ktAlpha(ktAlpha_), ktResolved(false),
    bipenalty(bipenalty_), bpMode(bpMode_), bpDt(bpDt_), bpBeta(bpBeta_),
    mPenalty(0.0), bpResolved(false),
    nDOF(0), nodeNdf(1 + hostNodes.Size()), dofOffset(1 + hostNodes.Size()),
    pflag(pressure_ ? 1 : 0),
    theNodes(0), K(0), P(0), M0(0)
{
  connectedNodes(0) = cNode;
  for (int i = 0; i < nHost; i++)
    connectedNodes(1 + i) = hostNodes(i);
  lambda.Zero();

  theNodes = new Node*[1 + nHost];
  for (int i = 0; i < 1 + nHost; i++) theNodes[i] = 0;
  // K/P/M0 are allocated in setDomain once the per-node ndf (⇒ nDOF) is known.
}

LadrunoEmbeddedNode::LadrunoEmbeddedNode()
  : Element(0, ELE_TAG_LadrunoEmbeddedNode),
    ndm(0), nHost(0), connectedNodes(0), Nshape(),
    Ku(0.0), Kp(0.0), upActive(false), lambda_p(0.0),
    enforce(0), lambda(),
    hostEleTag(-1), ktAuto(false), ktAlpha(0.0), ktResolved(false),
    bipenalty(false), bpMode(0), bpDt(0.0), bpBeta(0.0),
    mPenalty(0.0), bpResolved(false),
    nDOF(0), nodeNdf(), dofOffset(), pflag(0),
    theNodes(0), K(0), P(0), M0(0)
{
}

LadrunoEmbeddedNode::~LadrunoEmbeddedNode()
{
  if (theNodes != 0) delete[] theNodes;
  if (K != 0) delete K;
  if (P != 0) delete P;
  if (M0 != 0) delete M0;
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

// k_eff for the bipenalty bound. Phase 1 (U only) has a single translational spring,
// so k_eff = K_u (the (stiffness,inertia) per-DOF-class generalization of M1/ES-1
// arrives with UR/UP).
double LadrunoEmbeddedNode::effectiveCouplingStiffness(void)
{
  return LadrunoEmbedded::effectiveCouplingStiffness(fabs(Ku), 0.0);
}

// ADR 23 D5 / ADR 20 §10.6 — resolve the mass penalty m_p (lazily, after K_u).
// Lumped on the cNode's first ndm (translational) diagonal entries (getMass), so
// the global mass stays diagonal (DiagonalSOE-safe). Budgets:
//   bpMode 0 (-dtcr dt): m_p = k_eff·(dt/2)².
//   bpMode 1 (-wcap β):  m_p = k_eff/(β·ω_host)², ω_host² ≈ ‖K_host‖_∞/‖M_host‖_∞.
void LadrunoEmbeddedNode::resolveBipenalty(void)
{
  if (!bipenalty || bpResolved) return;
  this->resolveAutoKt();
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
  double kEff = this->effectiveCouplingStiffness();
  return LadrunoEmbedded::criticalTimeStep(mPenalty, kEff);
}

// ===========================================================================
//  state
// ===========================================================================
int LadrunoEmbeddedNode::commitState(void)
{
  // Augmented-Lagrangian per-step Uzawa update (ADR 23 D4): at the converged state,
  // accumulate the penalty traction into the multiplier. ISOTROPIC tie ⇒ NO
  // directional re-projection of λ (no preferred axis; ES-3 / M4).
  if (enforce == 1) {
    this->resolveAutoKt();
    Vector g(ndm);
    this->computeGap(g);
    for (int k = 0; k < ndm; k++) lambda(k) += Ku * g(k);
    if (upActive) lambda_p += Kp * this->computeGapP();
  }
  return this->Element::commitState();
}

int LadrunoEmbeddedNode::revertToLastCommit(void) { return 0; }

int LadrunoEmbeddedNode::revertToStart(void)
{
  lambda.Zero();   // AL multipliers accumulate only via commitState's Uzawa step
  lambda_p = 0.0;
  return 0;
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

  // isotropic traction t = K_u·g (+ λ for AL — constant within an inner solve).
  Vector t(ndm);
  for (int k = 0; k < ndm; k++) t(k) = Ku * g(k);
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
  return *P;
}

const Vector& LadrunoEmbeddedNode::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  // ADR 23 D5 — bipenalty inertia m_p·a on the cNode translational DOFs only.
  if (bipenalty) {
    this->resolveBipenalty();
    if (mPenalty != 0.0 && theNodes[0] != 0) {
      const Vector& aC = theNodes[0]->getTrialAccel();
      for (int k = 0; k < ndm; k++) (*P)(dofOffset(0) + k) += mPenalty * aC(k);
    }
  }
  return *P;
}

const Matrix& LadrunoEmbeddedNode::getTangentStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  // isotropic D_u = K_u·I_{ndm}
  Matrix Du(ndm, ndm);
  Du.Zero();
  for (int a = 0; a < ndm; a++) Du(a, a) = Ku;
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
  }
  return *M0;
}

// ===========================================================================
//  serialization
// ===========================================================================
int LadrunoEmbeddedNode::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  static Vector hdr(14);
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
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - header failed\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - ID failed\n";
    return -1;
  }
  Vector payload(nHost + ndm + 1);
  for (int i = 0; i < nHost; i++) payload(i) = Nshape(i);
  for (int k = 0; k < ndm; k++)
    payload(nHost + k) = (lambda.Size() == ndm) ? lambda(k) : 0.0;
  payload(nHost + ndm) = lambda_p;
  if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedNode::sendSelf - payload failed\n";
    return -1;
  }
  return 0;
}

int LadrunoEmbeddedNode::recvSelf(int commitTag, Channel& theChannel,
                                  FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  static Vector hdr(14);
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
  ktResolved = false;
  bpResolved = false;
  mPenalty = 0.0;
  upActive = false;   // re-resolved in setDomain (needs the nodes' ndf)
  nDOF = 0;           // recomputed in setDomain

  connectedNodes.resize(1 + nHost);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoEmbeddedNode::recvSelf - ID failed\n";
    return -1;
  }
  Vector payload(nHost + ndm + 1);
  if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
    opserr << "LadrunoEmbeddedNode::recvSelf - payload failed\n";
    return -1;
  }
  Nshape.resize(nHost);
  for (int i = 0; i < nHost; i++) Nshape(i) = payload(i);
  lambda.resize(ndm);
  for (int k = 0; k < ndm; k++) lambda(k) = payload(nHost + k);
  lambda_p = payload(nHost + ndm);

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
    for (int k = 0; k < ndm; k++) t(k) = Ku * g(k);
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
  default: return -1;
  }
}
