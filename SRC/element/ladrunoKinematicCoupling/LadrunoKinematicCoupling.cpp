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

// Ladruno: RBE2 / kinematic-coupling element (ADR 29). See
// LadrunoKinematicCoupling.h and
// Ladruno_implementation/29_ladruno_kinematic_coupling_rbe2_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoKinematicCoupling.h>
#include <LadrunoEmbeddedKernel.h>   // scalar bipenalty math (maxAbsDiagonal, massPenalty*, criticalTimeStep)
#include <classTags.h>
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
LadrunoKinematicCoupling::LadrunoKinematicCoupling(int tag, int ndm_, int refNode,
        const ID& slaveNodes, const ID& dofSel_, double kt, double kr, bool krUser_,
        int enforce_, bool bipenalty_, int bpMode_, double bpDt_, double bpBeta_,
        double kAlpha_, int hostEleTag_, bool ktAuto_, bool initGapCapture_)
  : Element(tag, ELE_TAG_LadrunoKinematicCoupling),
    ndm(ndm_), nrot((ndm_ == 3) ? 3 : 1), nSlave(slaveNodes.Size()),
    connectedNodes(1 + slaveNodes.Size()), dofSel(dofSel_),
    Kt(kt), Kr(kr), krUser(krUser_), ktAuto(ktAuto_), kAlpha(kAlpha_),
    hostEleTag(hostEleTag_), ktResolved(false), ell2(0.0),
    enforce(enforce_), lambdaAL(),
    bipenalty(bipenalty_), bpMode(bpMode_), bpDt(bpDt_), bpBeta(bpBeta_),
    bpResolved(false), hasRefRot(false),
    valid(false), dvec(), nGap(0), gapNode(), gapDof(), gapIsRot(),
    nDOF(0), nodeNdf(1 + slaveNodes.Size()), dofOffset(1 + slaveNodes.Size()),
    B(0), initGapCapture(initGapCapture_), g0Computed(false), g0(),
    theNodes(0), K(0), P(0), M0(0), C0(0), dampF(0)
{
  connectedNodes(0) = refNode;
  for (int i = 0; i < nSlave; i++) connectedNodes(1 + i) = slaveNodes(i);
  theNodes = new Node*[1 + nSlave];
  for (int i = 0; i < 1 + nSlave; i++) theNodes[i] = 0;
}

LadrunoKinematicCoupling::LadrunoKinematicCoupling()
  : Element(0, ELE_TAG_LadrunoKinematicCoupling),
    ndm(0), nrot(0), nSlave(0), connectedNodes(), dofSel(),
    Kt(0.0), Kr(0.0), krUser(false), ktAuto(false), kAlpha(0.0),
    hostEleTag(-1), ktResolved(false), ell2(0.0),
    enforce(0), lambdaAL(),
    bipenalty(false), bpMode(0), bpDt(0.0), bpBeta(0.0),
    bpResolved(false), hasRefRot(false),
    valid(false), dvec(), nGap(0), gapNode(), gapDof(), gapIsRot(),
    nDOF(0), nodeNdf(), dofOffset(),
    B(0), initGapCapture(true), g0Computed(false), g0(),
    theNodes(0), K(0), P(0), M0(0), C0(0), dampF(0)
{
}

LadrunoKinematicCoupling::~LadrunoKinematicCoupling()
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
int LadrunoKinematicCoupling::getNumExternalNodes(void) const { return 1 + nSlave; }
const ID& LadrunoKinematicCoupling::getExternalNodes(void) { return connectedNodes; }
Node** LadrunoKinematicCoupling::getNodePtrs(void) { return theNodes; }
int LadrunoKinematicCoupling::getNumDOF(void) { return nDOF; }

void LadrunoKinematicCoupling::allocate(void)
{
  if (K != 0) delete K;       K = new Matrix(nDOF, nDOF);
  if (P != 0) delete P;       P = new Vector(nDOF);
  if (M0 != 0) delete M0;     M0 = new Matrix(nDOF, nDOF);
  if (C0 != 0) delete C0;     C0 = new Matrix(nDOF, nDOF);   // always zero (getDamp)
  if (dampF != 0) delete dampF; dampF = new Vector(nDOF);    // always zero
}

void LadrunoKinematicCoupling::setDomain(Domain* theDomain)
{
  if (theDomain == 0) {
    if (theNodes != 0)
      for (int i = 0; i < 1 + nSlave; i++) theNodes[i] = 0;
    return;
  }
  // resolve nodes + lay out per-node DOF offsets. R needs ndf >= ndm (translations);
  // slaves need ndf >= ndm. Rotation features need ndf >= ndm+nrot on R (hasRefRot).
  nodeNdf.resize(1 + nSlave);
  dofOffset.resize(1 + nSlave);
  int pos = 0;
  bool nodeError = false;
  for (int i = 0; i < 1 + nSlave; i++) {
    theNodes[i] = theDomain->getNode(connectedNodes(i));
    if (theNodes[i] == 0) {
      opserr << "LadrunoKinematicCoupling::setDomain - node " << connectedNodes(i)
             << " not found\n";
      nodeError = true;
      break;
    }
    int ndf = theNodes[i]->getNumberDOF();
    if (ndf < ndm) {
      opserr << "LadrunoKinematicCoupling::setDomain - "
             << (i == 0 ? "reference" : "slave") << " node " << connectedNodes(i)
             << " has " << ndf << " DOFs, needs >= " << ndm << " (translations)\n";
      nodeError = true;
      break;
    }
    nodeNdf(i) = ndf;
    dofOffset(i) = pos;
    pos += ndf;
  }
  if (nodeError) { valid = false; return; }
  nDOF = pos;

  // degeneracy / wiring refusals (ADR 29 §3, review D6): self-tie + duplicate slave.
  // Set refused (NOT early-return) so we still fall through to allocate inert B/K/P —
  // returning before allocate would leave null matrices the assembler dereferences.
  bool refused = false;
  for (int i = 0; i < nSlave && !refused; i++) {
    if (connectedNodes(1 + i) == connectedNodes(0)) {
      opserr << "LadrunoKinematicCoupling " << this->getTag()
             << ": slave node " << connectedNodes(1 + i)
             << " equals the reference node — cannot tie a node to itself; refusing\n";
      refused = true; break;
    }
    for (int j = i + 1; j < nSlave; j++)
      if (connectedNodes(1 + j) == connectedNodes(1 + i)) {
        opserr << "LadrunoKinematicCoupling " << this->getTag()
               << ": slave node " << connectedNodes(1 + i)
               << " is listed more than once; refusing (over-stiffens / over-counts)\n";
        refused = true; break;
      }
  }

  if (refused) { valid = false; nGap = 0; }
  else this->resolveGeometry();   // d_i, ragged layout (nGap, gap* arrays), ell2, hasRefRot

  // Allocate B/K/P/M0 even if ill-posed: a zeroed inert B contributes nothing instead
  // of leaving null matrices the assembler would dereference (WARNING already printed).
  if (B != 0) delete B;
  B = new Matrix((nGap > 0 ? nGap : 1), (nDOF > 0 ? nDOF : 1));
  B->Zero();
  this->allocate();
  // init lambdaAL / g0 — check whether recvSelf already restored them (size == nGap)
  // BEFORE the resize destroys the size-0 signal of a fresh (non-recv) construct.
  bool lamRestored = (nGap > 0 && lambdaAL.Size() == nGap);
  bool g0Restored = (nGap > 0 && g0.Size() == nGap && g0Computed);
  lambdaAL.resize(nGap > 0 ? nGap : 1);
  if (!lamRestored || !valid) lambdaAL.Zero();
  g0.resize(nGap > 0 ? nGap : 1);
  if (!g0Restored) g0.Zero();
  if (valid) {
    this->buildB();
    this->captureInitialGap();   // born stress-free at the (possibly deformed) current state
  }
  this->DomainComponent::setDomain(theDomain);
}

// resolve d_i = x_i − x_R, the ragged per-slave tied-DOF layout (from -dof ∩ each
// slave's ndf and R's rotation availability), and the floored rotation length scale
// ℓ² (ADR 29 §3, §6). REFUSES inertly (valid=false) on an empty effective selection
// or a rotation explicitly requested through a reference that carries none.
void LadrunoKinematicCoupling::resolveGeometry(void)
{
  valid = false;
  nGap = 0;                  // inert until the layout is built (early returns leave it 0)

  // R carries usable rotations iff its ndf covers all nrot rotation slots.
  hasRefRot = (nodeNdf(0) >= ndm + nrot);

  // d_i = x_i − x_R, and the length scale.
  dvec.resize(nSlave, ndm);
  const Vector& xR = theNodes[0]->getCrds();
  double sumD2 = 0.0, maxD2 = 0.0;
  int nTransScale = 0;
  for (int i = 0; i < nSlave; i++) {
    const Vector& xi = theNodes[1 + i]->getCrds();
    double d2 = 0.0;
    for (int k = 0; k < ndm; k++) {
      double dk = xi(k) - xR(k);
      dvec(i, k) = dk;
      d2 += dk * dk;
    }
    sumD2 += d2;
    if (d2 > maxD2) maxD2 = d2;
    nTransScale++;
  }

  // does -dof explicitly request a rotation component? (range was validated at parse)
  bool reqRot = false;
  bool useDefault = (dofSel.Size() == 0);
  if (!useDefault)
    for (int s = 0; s < dofSel.Size(); s++)
      if (dofSel(s) > ndm) reqRot = true;

  if (reqRot && !hasRefRot) {
    opserr << "LadrunoKinematicCoupling " << this->getTag()
           << ": -dof requested a slave rotation but the reference node carries no "
           << "rotation DOFs — moment transfer impossible; refusing (give R "
           << (ndm + nrot) << " DOFs or drop the rotation components)\n";
    return;
  }
  if (!hasRefRot)
    opserr << "LadrunoKinematicCoupling " << this->getTag()
           << ": reference node has no rotation DOFs — tying translations only "
           << "(multi-slave equalDOF, no moment transfer / transport)\n";

  // Build the ragged layout. A row is one scalar constraint (slave, component).
  // Pass 1: count. Pass 2: fill. Component c (1-based): 1..ndm = translation
  // (node DOF c−1), ndm+1..ndm+nrot = rotation (node DOF c−1).
  int count = 0;
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      nGap = count;
      gapNode.resize(nGap > 0 ? nGap : 1);
      gapDof.resize(nGap > 0 ? nGap : 1);
      gapIsRot.resize(nGap > 0 ? nGap : 1);
      count = 0;
    }
    for (int i = 0; i < nSlave; i++) {
      int sndf = nodeNdf(1 + i);
      // which components for this slave: explicit list or default-all-available.
      int cmax = ndm + (hasRefRot ? nrot : 0);
      for (int c = 1; c <= cmax; c++) {
        bool wanted;
        if (useDefault) wanted = true;
        else {
          wanted = false;
          for (int s = 0; s < dofSel.Size(); s++)
            if (dofSel(s) == c) { wanted = true; break; }
        }
        if (!wanted) continue;
        bool isRot = (c > ndm);
        // slave must possess this DOF (node DOF index c−1).
        if ((c - 1) >= sndf) {
          if (pass == 0 && !useDefault)
            opserr << "LadrunoKinematicCoupling " << this->getTag() << ": slave node "
                   << connectedNodes(1 + i) << " lacks DOF component " << c
                   << " — skipped\n";
          continue;
        }
        if (isRot && !hasRefRot) continue;   // no θ_R to tie against
        if (pass == 1) {
          gapNode(count) = 1 + i;
          gapDof(count) = c - 1;
          gapIsRot(count) = isRot ? 1 : 0;
        }
        count++;
      }
    }
  }

  if (nGap == 0) {
    opserr << "LadrunoKinematicCoupling " << this->getTag()
           << ": no DOFs tied after intersecting -dof with the slaves' DOFs; refusing\n";
    return;
  }

  // floored rotation length scale (ADR 29 §6 D1): never let ℓ²=0 zero a tied rotation.
  double meanD2 = (nTransScale > 0) ? (sumD2 / nTransScale) : 0.0;
  ell2 = meanD2;
  if (ell2 < maxD2) ell2 = maxD2;            // floor toward the largest lever
  if (ell2 <= 0.0) {
    // all slaves coincident with R: if a rotation is tied, use a unit length so the
    // direct rotation tie still has stiffness; warn.
    bool anyRot = false;
    for (int r = 0; r < nGap; r++) if (gapIsRot(r)) { anyRot = true; break; }
    if (anyRot) {
      ell2 = 1.0;
      opserr << "LadrunoKinematicCoupling " << this->getTag()
             << ": all slaves coincident with the reference; using unit length for the "
             << "rotational penalty scale (supply -kr to set it explicitly)\n";
    }
  }
  valid = true;
}

// transport operator: ∂(θ × d_i)_k / ∂θ_r. The g_i^t B-block is the NEGATIVE of this.
double LadrunoKinematicCoupling::transOp(int i, int k, int r) const
{
  if (ndm == 3) {
    double dx = dvec(i, 0), dy = dvec(i, 1), dz = dvec(i, 2);
    switch (k) {                                  // (θ × d)_k, d/dθ_r
    case 0: return (r == 1) ?  dz : (r == 2) ? -dy : 0.0;   // (θ×d)_x = θy dz − θz dy
    case 1: return (r == 0) ? -dz : (r == 2) ?  dx : 0.0;   // (θ×d)_y = θz dx − θx dz
    case 2: return (r == 0) ?  dy : (r == 1) ? -dx : 0.0;   // (θ×d)_z = θx dy − θy dx
    default: return 0.0;
    }
  }
  // 2D drilling: θ (scalar) × d = θ (−d_y, d_x). ∂(θ×d)_k/∂θ.
  return (k == 0) ? -dvec(i, 1) : dvec(i, 0);
}

// build the constant gap operator B (nGap × nDOF): each row is one scalar constraint.
//   translation row (slave i, node DOF j=c−1): +u_i[j] − u_R[j] − Σ_r transOp(i,j,r) θ_R[r]
//     (the B-block ∂g/∂θ_R = +[d_i]_× = −transOp — the SIGN FLIP vs RBE3, ADR 29 §2.4)
//   rotation row    (slave i, node DOF j=c−1): +θ_i[j] − θ_R[j]
void LadrunoKinematicCoupling::buildB(void)
{
  B->Zero();
  int oR = dofOffset(0);
  for (int row = 0; row < nGap; row++) {
    int i = gapNode(row) - 1;                  // slave index 0..N-1
    int oS = dofOffset(gapNode(row));
    int j = gapDof(row);                       // node-local DOF index (component − 1)
    (*B)(row, oS + j) += 1.0;                   // dependent (slave) DOF
    (*B)(row, oR + j) += -1.0;                  // reference DOF (same component index)
    if (gapIsRot(row) == 0 && hasRefRot) {
      for (int r = 0; r < nrot; r++)
        (*B)(row, oR + ndm + r) += -this->transOp(i, j, r);   // +[d_i]_×  (= −transOp)
    }
  }
}

// -k auto (needs a representative -host) + derive K_r = K_t·ℓ² (unless -kr given). Lazy.
void LadrunoKinematicCoupling::resolveAutoKt(void)
{
  if (ktResolved) return;
  if (ktAuto) {
    Domain* theDomain = this->getDomain();
    Element* host = (theDomain != 0 && hostEleTag >= 0) ? theDomain->getElement(hostEleTag) : 0;
    if (host == 0) {
      opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
             << ": -k auto needs a valid -host element; keeping K_t=" << Kt << "\n";
    } else {
      double scale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
      if (scale > 0.0) Kt = kAlpha * scale;
      else
        opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
               << ": -k auto host gave a zero stiffness scale; keeping K_t=" << Kt << "\n";
    }
  }
  if (!krUser) Kr = Kt * ell2;   // derived rotational penalty (ADR 29 §6; ℓ² floored)
  ktResolved = true;
}

// bipenalty (ADR 29 §6): scan every tied DOF of R AND every slave; where the node is
// actually massless on that DOF, lump a penalty mass sized from the Gershgorin row-sum
// of the assembled penalty tangent K (≥ λ_max ⇒ conservative, covers the rotation
// coupling without an eigensolve). Builds the diagonal *M0. Default OFF.
void LadrunoKinematicCoupling::resolveBipenalty(void)
{
  if (!bipenalty || bpResolved) return;
  this->resolveAutoKt();
  M0->Zero();
  if (!valid) { bpResolved = true; return; }

  // assemble the penalty tangent once (for Gershgorin row-sums).
  this->getTangentStiff();   // fills *K

  // ω_host² for -wcap (needs -host).
  double omega2_host = 0.0;
  if (bpMode == 1) {
    Domain* theDomain = this->getDomain();
    Element* host = (theDomain != 0 && hostEleTag >= 0) ? theDomain->getElement(hostEleTag) : 0;
    if (host == 0) {
      opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
             << ": -wcap needs a valid -host for omega_host; use -dtcr; bipenalty off\n";
      bpResolved = true; return;
    }
    double kScale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
    double mScale = LadrunoEmbedded::maxAbsDiagonal(host->getMass());
    if (mScale <= 0.0 || kScale <= 0.0) {
      opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
             << ": -wcap host has no mass (or no stiffness); use -dtcr; bipenalty off\n";
      bpResolved = true; return;
    }
    omega2_host = kScale / mScale;
    if (bpBeta <= 0.0) bpBeta = 2.0;
  } else if (bpDt <= 0.0) {
    opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
           << ": -dtcr needs a positive step; bipenalty off\n";
    bpResolved = true; return;
  }

  // per element-DOF: if the owning node is massless on that DOF, lump a penalty mass.
  for (int p = 0; p < 1 + nSlave; p++) {
    const Matrix& mN = theNodes[p]->getMass();
    int ndf = nodeNdf(p);
    bool haveNodalMass = (mN.noRows() == ndf && mN.noCols() == ndf);
    for (int d = 0; d < ndf; d++) {
      int e = dofOffset(p) + d;
      double kDof = 0.0;                          // Gershgorin row-sum of K at this DOF
      for (int c = 0; c < nDOF; c++) kDof += fabs((*K)(e, c));
      if (kDof <= 0.0) continue;                  // DOF not coupled (untied) — leave it
      double nodeMassDd = haveNodalMass ? mN(d, d) : 0.0;
      if (nodeMassDd > 0.0) continue;             // node already carries mass here — don't double-count
      double mp = (bpMode == 0) ? LadrunoEmbedded::massPenaltyDtcr(kDof, bpDt)
                                : LadrunoEmbedded::massPenaltyWcap(kDof, bpBeta, omega2_host);
      (*M0)(e, e) = mp;
    }
  }
  bpResolved = true;
}

int LadrunoKinematicCoupling::setRayleighDampingFactors(double, double, double, double)
{
  return 0;   // a pure penalty coupling carries no physical Rayleigh damping
}

// D ≡ 0. Overridden (with getRayleighDampingForces) so the base Element's lazy
// damping-matrix slot is never needed — see the header note + ADR 29 §6 (transient
// C-tangent index-landmine). allocate() builds C0/dampF sized nDOF.
const Matrix& LadrunoKinematicCoupling::getDamp(void)
{
  C0->Zero();
  return *C0;
}

const Vector& LadrunoKinematicCoupling::getRayleighDampingForces(void)
{
  dampF->Zero();
  return *dampF;
}

double LadrunoKinematicCoupling::getExplicitCriticalTimeStep(void)
{
  if (!bipenalty) return -1.0;
  this->resolveBipenalty();
  if (!valid) return -1.0;
  // min over lumped (massless) DOFs of 2√(m_p/k_dof). Each m_p was sized so this equals
  // the -dtcr budget; we take the true min in case -wcap produced a tighter bound.
  double dt = -1.0;
  for (int p = 0; p < 1 + nSlave; p++) {
    int ndf = nodeNdf(p);
    for (int d = 0; d < ndf; d++) {
      int e = dofOffset(p) + d;
      double mp = (*M0)(e, e);
      if (mp <= 0.0) continue;
      double kDof = 0.0;
      for (int c = 0; c < nDOF; c++) kDof += fabs((*K)(e, c));
      double dtd = LadrunoEmbedded::criticalTimeStep(mp, kDof);
      if (dtd > 0.0 && (dt < 0.0 || dtd < dt)) dt = dtd;
    }
  }
  return dt;
}

// ===========================================================================
//  state
// ===========================================================================
int LadrunoKinematicCoupling::commitState(void)
{
  if (enforce == 1 && valid) {                // augmented-Lagrangian Uzawa update
    this->resolveAutoKt();
    Vector g(nGap);
    this->computeGap(g);
    for (int row = 0; row < nGap; row++)
      lambdaAL(row) += this->rowPenalty(row) * g(row);
  }
  return this->Element::commitState();
}

int LadrunoKinematicCoupling::revertToLastCommit(void) { return 0; }

int LadrunoKinematicCoupling::revertToStart(void)
{
  lambdaAL.Zero();
  return 0;
}

int LadrunoKinematicCoupling::update(void)
{
  this->resolveAutoKt();
  return 0;
}

// full gap g = B u (then minus the captured offsets g0). g is sized nGap.
void LadrunoKinematicCoupling::computeGap(Vector& g)
{
  if (g.Size() != nGap) g.resize(nGap);
  g.Zero();
  if (!valid) return;
  Vector uFull(nDOF); uFull.Zero();
  for (int p = 0; p < 1 + nSlave; p++) {
    const Vector& u = theNodes[p]->getTrialDisp();
    int ndf = nodeNdf(p);
    for (int d = 0; d < ndf; d++) uFull(dofOffset(p) + d) = u(d);
  }
  for (int row = 0; row < nGap; row++) {
    double s = 0.0;
    for (int c = 0; c < nDOF; c++) s += (*B)(row, c) * uFull(c);
    g(row) = s;
  }
  if (g0Computed)
    for (int row = 0; row < nGap; row++) g(row) -= g0(row);
}

void LadrunoKinematicCoupling::captureInitialGap(void)
{
  if (!initGapCapture || g0Computed) return;
  Vector g(nGap);
  this->computeGap(g);                 // absolute gap (g0Computed still false)
  for (int row = 0; row < nGap; row++) g0(row) = g(row);
  g0Computed = true;
}

// ===========================================================================
//  matrices / forces
// ===========================================================================
const Vector& LadrunoKinematicCoupling::getResistingForce(void)
{
  this->resolveAutoKt();
  P->Zero();
  if (!valid) return *P;
  Vector g(nGap);
  this->computeGap(g);
  Vector t(nGap);
  for (int row = 0; row < nGap; row++) {
    t(row) = this->rowPenalty(row) * g(row);
    if (enforce == 1) t(row) += lambdaAL(row);
  }
  // P = Bᵀ t
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

const Vector& LadrunoKinematicCoupling::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  if (bipenalty && valid) {
    this->resolveBipenalty();
    for (int p = 0; p < 1 + nSlave; p++) {
      const Vector& a = theNodes[p]->getTrialAccel();
      int ndf = nodeNdf(p);
      for (int d = 0; d < ndf; d++) {
        int e = dofOffset(p) + d;
        double mp = (*M0)(e, e);
        if (mp != 0.0) (*P)(e) += mp * a(d);
      }
    }
  }
  return *P;
}

const Matrix& LadrunoKinematicCoupling::getTangentStiff(void)
{
  this->resolveAutoKt();
  K->Zero();
  if (!valid) return *K;
  for (int row = 0; row < nGap; row++) {
    double dval = this->rowPenalty(row);
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

const Matrix& LadrunoKinematicCoupling::getInitialStiff(void)
{
  return this->getTangentStiff();   // B, K_t, K_r all state-independent
}

const Matrix& LadrunoKinematicCoupling::getMass(void)
{
  // resolveBipenalty fills (and internally zeroes) *M0 ONCE (bpResolved guard); do NOT
  // pre-zero here or a 2nd getMass would wipe the persisted lumped masses. Non-bipenalty:
  // M0 is the fresh zero matrix from allocate().
  if (bipenalty) this->resolveBipenalty();
  else M0->Zero();
  return *M0;
}

// ===========================================================================
//  serialization (rebuild geometry/layout/B on recv from coords at setDomain)
// ===========================================================================
int LadrunoKinematicCoupling::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  static Vector hdr(20);
  hdr(0) = this->getTag();
  hdr(1) = ndm;
  hdr(2) = nSlave;
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
  hdr(16) = nGap;
  hdr(17) = dofSel.Size();
  hdr(18) = ell2;
  hdr(19) = 1.0;                          // version
  if (theChannel.sendVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoKinematicCoupling::sendSelf - header failed\n";
    return -1;
  }
  if (theChannel.sendID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoKinematicCoupling::sendSelf - nodes ID failed\n";
    return -1;
  }
  if (dofSel.Size() > 0 && theChannel.sendID(dbTag, commitTag, dofSel) < 0) {
    opserr << "LadrunoKinematicCoupling::sendSelf - dofSel ID failed\n";
    return -1;
  }
  // payload: lambdaAL(nGap) + g0(nGap). nGap is recomputed deterministically in
  // setDomain on the recv side; the ordering matches, so these align by index.
  if (nGap > 0) {
    Vector payload(2 * nGap);
    for (int r = 0; r < nGap; r++)
      payload(r) = (lambdaAL.Size() == nGap) ? lambdaAL(r) : 0.0;
    for (int r = 0; r < nGap; r++)
      payload(nGap + r) = (g0.Size() == nGap) ? g0(r) : 0.0;
    if (theChannel.sendVector(dbTag, commitTag, payload) < 0) {
      opserr << "LadrunoKinematicCoupling::sendSelf - payload failed\n";
      return -1;
    }
  }
  return 0;
}

int LadrunoKinematicCoupling::recvSelf(int commitTag, Channel& theChannel,
                                       FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  static Vector hdr(20);
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoKinematicCoupling::recvSelf - header failed\n";
    return -1;
  }
  this->setTag((int)hdr(0));
  ndm = (int)hdr(1);
  nSlave = (int)hdr(2);
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
  nGap = (int)hdr(16);
  int nDofSel = (int)hdr(17);
  ell2 = hdr(18);
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  bpResolved = false;
  nDOF = 0;
  valid = false;

  connectedNodes.resize(1 + nSlave);
  if (theChannel.recvID(dbTag, commitTag, connectedNodes) < 0) {
    opserr << "LadrunoKinematicCoupling::recvSelf - nodes ID failed\n";
    return -1;
  }
  if (nDofSel > 0) {
    dofSel.resize(nDofSel);
    if (theChannel.recvID(dbTag, commitTag, dofSel) < 0) {
      opserr << "LadrunoKinematicCoupling::recvSelf - dofSel ID failed\n";
      return -1;
    }
  } else {
    dofSel = ID();
  }
  lambdaAL.resize(nGap > 0 ? nGap : 1); lambdaAL.Zero();
  g0.resize(nGap > 0 ? nGap : 1);       g0.Zero();
  if (nGap > 0) {
    Vector payload(2 * nGap);
    if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
      opserr << "LadrunoKinematicCoupling::recvSelf - payload failed\n";
      return -1;
    }
    for (int r = 0; r < nGap; r++) lambdaAL(r) = payload(r);
    for (int r = 0; r < nGap; r++) g0(r) = payload(nGap + r);
  }

  nodeNdf.resize(1 + nSlave);
  dofOffset.resize(1 + nSlave);
  if (theNodes != 0) delete[] theNodes;
  theNodes = new Node*[1 + nSlave];
  for (int i = 0; i < 1 + nSlave; i++) theNodes[i] = 0;
  // d_i, ragged layout, ℓ², B are rebuilt in setDomain from coords. g0Computed is kept
  // (restored above) so setDomain does NOT recapture — the absolute tie is preserved.
  return 0;
}

void LadrunoKinematicCoupling::Print(OPS_Stream& s, int flag)
{
  s << "LadrunoKinematicCoupling (RBE2), tag: " << this->getTag() << "\n";
  s << "  reference node: " << connectedNodes(0) << "  slave nodes: ";
  for (int i = 0; i < nSlave; i++) s << connectedNodes(1 + i) << " ";
  s << "\n  K_t: " << Kt << "  K_r: " << Kr << "  enforce: "
    << (enforce == 1 ? "al" : "penalty");
  s << "  tied DOFs (nGap): " << nGap << (hasRefRot ? "  +moment-transfer" : "  (translation-only)");
  s << (initGapCapture ? (g0Computed ? "  +initGap(captured)" : "  +initGap")
                       : "  +absolute(no initGap)");
  s << "\n";
}

// ===========================================================================
//  responses
// ===========================================================================
Response* LadrunoKinematicCoupling::setResponse(const char** argv, int argc, OPS_Stream& s)
{
  if (argc < 1) return 0;
  int n = (nGap > 0) ? nGap : 1;
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "couplingForce") == 0)
    return new ElementResponse(this, 1, Vector(n));
  if (strcmp(argv[0], "gap") == 0)
    return new ElementResponse(this, 2, Vector(n));
  if (strcmp(argv[0], "kt") == 0 || strcmp(argv[0], "k") == 0)
    return new ElementResponse(this, 3, 0.0);
  if (strcmp(argv[0], "kr") == 0)
    return new ElementResponse(this, 4, 0.0);
  if (strcmp(argv[0], "lambda") == 0 || strcmp(argv[0], "augLambda") == 0)
    return new ElementResponse(this, 5, Vector(n));
  if (strcmp(argv[0], "dtcr") == 0 || strcmp(argv[0], "dtCritical") == 0)
    return new ElementResponse(this, 6, 0.0);
  if (strcmp(argv[0], "nGap") == 0 || strcmp(argv[0], "tiedDOFs") == 0)
    return new ElementResponse(this, 7, 0.0);
  return 0;
}

int LadrunoKinematicCoupling::getResponse(int responseID, Information& eleInfo)
{
  this->resolveAutoKt();
  switch (responseID) {
  case 1: {
    Vector g(nGap); this->computeGap(g);
    Vector t(nGap);
    for (int row = 0; row < nGap; row++) {
      t(row) = this->rowPenalty(row) * g(row);
      if (enforce == 1) t(row) += lambdaAL(row);
    }
    return eleInfo.setVector(t);
  }
  case 2: { Vector g(nGap); this->computeGap(g); return eleInfo.setVector(g); }
  case 3: return eleInfo.setDouble(Kt);
  case 4: return eleInfo.setDouble(Kr);
  case 5: return eleInfo.setVector(lambdaAL);
  case 6: { double dt = this->getExplicitCriticalTimeStep(); return eleInfo.setDouble(dt > 0.0 ? dt : 0.0); }
  case 7: return eleInfo.setDouble((double)nGap);
  default: return -1;
  }
}
