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
#include <LadrunoResponseTokens.h>   // Ladruno — shared recorder-token aliases
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <EquiSolnAlgo.h>        // Ladruno (WP-101 r1): -alUpdate iter guard
#include <StaticIntegrator.h>    // Ladruno (WP-101 r1)
#include <TransientIntegrator.h> // Ladruno (WP-101 r1)
#include <math.h>
#include <string.h>

// ===========================================================================
//  construction
// ===========================================================================
LadrunoKinematicCoupling::LadrunoKinematicCoupling(int tag, int ndm_, int refNode,
        const ID& slaveNodes, const ID& dofSel_, double kt, double kr, bool krUser_,
        int enforce_, bool bipenalty_, int bpMode_, double bpDt_, double bpBeta_,
        double kAlpha_, int hostEleTag_, bool ktAuto_, bool initGapCapture_,
        int alUpdate_)
  : Element(tag, ELE_TAG_LadrunoKinematicCoupling),
    ndm(ndm_), nrot((ndm_ == 3) ? 3 : 1), nSlave(slaveNodes.Size()),
    connectedNodes(1 + slaveNodes.Size()), dofSel(dofSel_),
    Kt(kt), Kr(kr), krUser(krUser_), ktAuto(ktAuto_), kAlpha(kAlpha_),
    hostEleTag(hostEleTag_), ktResolved(false), ktHostMissWarned(false), ell2(0.0),
    enforce(enforce_), alUpdate(alUpdate_), alGuardWarned(false),
    alWarnedTransient(false), lambdaAL(), lambdaCommitted(), alSkipUpdate(false),
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
    hostEleTag(-1), ktResolved(false), ktHostMissWarned(false), ell2(0.0),
    enforce(0), alUpdate(0), alGuardWarned(false),
    alWarnedTransient(false), lambdaAL(), lambdaCommitted(), alSkipUpdate(false),
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
  // Ladruno (WP-101): lambdaCommitted mirrors lambdaAL — it is the revertToLastCommit
  // target, so it must exist and match in size before the first update()/commitState().
  // Same first-call-detection-before-resize rule as lambdaAL (the resize destroys the
  // size-0 fresh-vs-recv signal): a recv-restored pair is kept as received.
  bool lamCommRestored = (nGap > 0 && lambdaCommitted.Size() == nGap && lamRestored);
  lambdaCommitted.resize(nGap > 0 ? nGap : 1);
  if (!lamCommRestored || !valid) lambdaCommitted = lambdaAL;
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

  // Ladruno (WP-101 r1): resolve the -host ONCE, up front, and DO NOT latch if it is not in
  // the domain yet. Domain::addElement() calls update() -> resolveAutoKt() at DECLARATION
  // time; a deck that declares the coupling before its host element used to latch on that
  // failed lookup, so `-k auto` silently stayed at the 1e12 default and the conditioning
  // warning below could never fire. Derive K_r provisionally and retry on the next call.
  Element* host = 0;
  if (hostEleTag >= 0) {
    Domain* theDomain = this->getDomain();
    if (theDomain != 0) host = theDomain->getElement(hostEleTag);
    if (host == 0) {
      if (!krUser) Kr = Kt * ell2;           // provisional, recomputed once the host lands
      // Warn once -- but only when an analysis exists, i.e. we are past deck construction and
      // the host is genuinely missing rather than merely not declared YET.
      EquiSolnAlgo** algoPtr = OPS_GetAlgorithm();
      if (!ktHostMissWarned && algoPtr != 0 && *algoPtr != 0) {
        ktHostMissWarned = true;
        opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
               << ": -host element " << hostEleTag << " is not in the domain"
               << (ktAuto ? "; -k auto cannot resolve and K_t stays "
                          : "; K_t stays ") << Kt
               << " (the conditioning check is skipped too)\n";
      }
      return;                                 // NO latch -- retry next call
    }
  }

  if (ktAuto) {
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

  // Ladruno (WP-101): conditioning guard on a NUMERIC -k when a -host is named. The
  // rigidity error of a penalty tie falls as c/K_t while the condition number of the
  // assembled system rises linearly with K_t, so there is a band, not a "bigger is
  // better" (guide §3.1). Measured on the TIMs strip (E = 45 MPa soil, B = 1.5 m
  // footing, 101 583 DOF): K_t = 1e12 drove Pardiso onto perturbed pivots and then to
  // failure (SuperLU failed its first factorisation) while K_t = 5e9 carried the same
  // leg to a clean plateau at a rigidity error of 6.7e-7. 1e6x the host's diagonal
  // stiffness is the loud end of the recommended 1e2..1e4x band.
  //
  // This lives here and NOT in the parser on purpose: the check needs the host element's
  // assembled initial stiffness, and at parse time the host may not exist yet (element
  // ordering in the deck is free) nor have had setDomain() called — calling
  // getInitialStiff() there is a null-node dereference waiting to happen. resolveAutoKt()
  // is the first point where the value is genuinely available, is ktResolved-guarded (so
  // the warning fires exactly once), and still runs before the first factorisation.
  if (!ktAuto && host != 0 && Kt > 0.0) {
    {
      double scale = LadrunoEmbedded::maxAbsDiagonal(host->getInitialStiff());
      if (scale > 0.0 && Kt > 1.0e6 * scale)
        opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
               << ": -k " << Kt << " is " << (Kt / scale)
               << "x the -host element's stiffness scale (" << scale
               << "). Above ~1e6x the penalty block dominates the global matrix and the "
               << "solve goes ill-conditioned (perturbed pivots / failed factorisation) "
               << "while the rigidity error stops improving. Recommended band: 1e2..1e4x "
               << "the host diagonal; use -enforce al to tighten the tie at a MODERATE "
               << "K_t instead. See LadrunoKinematicCoupling_guide.md section 3.1\n";
    }
  }
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
  if (enforce == 1 && valid) {
    this->resolveAutoKt();
    // Ladruno (WP-101): the LEGACY per-commit Uzawa step (-alUpdate commit). Under the
    // default -alUpdate iter the recursion has already advanced once per Newton iteration
    // inside update(), so adding another D·g here would double-count the last iterate.
    if (alUpdate == 0) {
      Vector g(nGap);
      this->computeGap(g);
      for (int row = 0; row < nGap; row++)
        lambdaAL(row) += this->rowPenalty(row) * g(row);
    }
    // snapshot: the multipliers this step converged with, and the revert target.
    if (lambdaCommitted.Size() != lambdaAL.Size())
      lambdaCommitted.resize(lambdaAL.Size());
    lambdaCommitted = lambdaAL;
  }
  return this->Element::commitState();
}

// Ladruno (WP-101): a FAILED / retried step must not inherit the multipliers accumulated
// on the discarded trial state. With the per-iteration recursion lambdaAL moves inside the
// step, so the (previously empty) revert MUST roll it back to the last committed value.
int LadrunoKinematicCoupling::revertToLastCommit(void)
{
  if (enforce == 1 && lambdaCommitted.Size() == lambdaAL.Size())
    lambdaAL = lambdaCommitted;
  alSkipUpdate = true;      // Domain::revertToLastCommit ends with update() — see the header
  return 0;
}

int LadrunoKinematicCoupling::revertToStart(void)
{
  lambdaAL.Zero();
  lambdaCommitted.Zero();
  alSkipUpdate = true;      // Domain::revertToStart ends with update() too
  return 0;
}

// Called by Domain::update() once per equilibrium iteration, on the CURRENT trial
// displacements. Ladruno (WP-101): this is where the augmented-Lagrangian recursion
// advances under the default -alUpdate iter — Uzawa nested in Newton.
//
//   λ_{k+1} = λ_k + D g(u_k),   r = f − S u_k − Bᵀ(λ_{k+1} + D g(u_k)),   T = S + BᵀDB
//
// The tangent deliberately stays the penalty operator (λ is frozen w.r.t. u), so with a
// linear structure one Newton solve IS the Uzawa inner solve and the multiplier error
// contracts by (I + D B S⁻¹Bᵀ)⁻¹ per iteration. Stationarity forces D g = 0, i.e. the
// converged state satisfies the constraint EXACTLY (to the algorithm's tolerance) rather
// than to the penalty's O(1/K) floor — that is the whole point of the change.
//
// Caveat (documented, not defended against): a line-search / trial-and-discard algorithm
// calls update() on states it then throws away, so λ rides those iterates too. The fixed
// point is unchanged; only the path is. revertToLastCommit restores the step's start value.
int LadrunoKinematicCoupling::update(void)
{
  this->resolveAutoKt();
  this->warnAlUnderTransient();
  if (alSkipUpdate) {       // consume the post-revert / post-recv call; do NOT advance on it
    alSkipUpdate = false;
    return 0;
  }
  if (enforce == 1 && alUpdate == 1 && valid) {
    if (this->refuseIterCadence()) return -1;
    Vector g(nGap);
    this->computeGap(g);
    for (int row = 0; row < nGap; row++)
      lambdaAL(row) += this->rowPenalty(row) * g(row);
  }
  return 0;
}

// Ladruno (WP-101 r1): gate for the opt-in per-iteration cadence. `iter` makes the residual
// path-dependent (see the header note), which is only survivable under FULL NEWTON +
// LoadControl. Everything else was MEASURED to fail or stagnate, so refuse loudly instead of
// returning a silently wrong answer. Returns true when the caller must abort.
//
// Deliberately NOT latched on a pass: a deck can swap algorithm/integrator between analyze()
// calls, and re-reading two pointers per iteration is free. Says nothing when it cannot tell
// (no analysis configured yet -- Domain::addElement calls update() at declaration time, long
// before `analysis Static` exists).
bool LadrunoKinematicCoupling::refuseIterCadence(void)
{
  EquiSolnAlgo** algoPtr = OPS_GetAlgorithm();
  EquiSolnAlgo* algo = (algoPtr != 0) ? *algoPtr : 0;
  if (algo == 0) return false;                       // no opinion yet

  StaticIntegrator** siPtr = OPS_GetStaticIntegrator();
  StaticIntegrator* si = (siPtr != 0) ? *siPtr : 0;
  TransientIntegrator** tiPtr = OPS_GetTransientIntegrator();
  TransientIntegrator* ti = (tiPtr != 0) ? *tiPtr : 0;

  bool algoOK = (algo->getClassTag() == EquiALGORITHM_TAGS_NewtonRaphson);
  bool integOK = (si != 0 && si->getClassTag() == INTEGRATOR_TAGS_LoadControl && ti == 0);
  if (algoOK && integOK) return false;

  if (!alGuardWarned) {
    alGuardWarned = true;
    opserr << "LadrunoKinematicCoupling " << this->getTag()
           << ": -enforce al -alUpdate iter is REFUSED here. The per-iteration Uzawa update "
           << "makes the residual path-dependent (the tie force carries lambda_k + 2*D*g while "
           << "the tangent linearises one D*g), so only FULL NEWTON under LoadControl is safe. "
           << "This analysis has algorithm classTag " << algo->getClassTag() << " and "
           << (ti != 0 ? "a TRANSIENT integrator" : (si != 0 ? "a static integrator" : "NO integrator"));
    if (ti == 0 && si != 0) opserr << " classTag " << si->getClassTag();
    opserr << ". Measured: DisplacementControl fails 5/5 steps with EVERY algorithm; "
           << "KrylovNewton / BFGS / Broyden fail or diverge under LoadControl; ModifiedNewton "
           << "fails 10/10 on a nonlinear host. Use the DEFAULT -alUpdate commit and, to close "
           << "the constraint WITHIN a step, wrap the step in the ADR-41 held-load augmentation "
           << "sweep (ladrunoBeginAugment; integrator LoadControl 0.0; analyze 1 repeatedly "
           << "until eleResponse <tag> constraintViolation is small; ladrunoEndAugment) -- that "
           << "is a proper OUTER Uzawa loop and was measured to close this gate with every "
           << "algorithm and integrator. See LadrunoKinematicCoupling_guide.md section 4.2\n";
  }
  return true;
}

// Ladruno (WP-101 r2): is this integrator class an EXPLICIT one? The r1 cut tested
// `TransientIntegrator != 0`, which also caught IMPLICIT transient integrators -- and
// -enforce al under Newmark + Newton works perfectly (measured 0/10 failed steps, max gap
// 3.4e-21), so that warning was actively wrong. Only the explicit family has the problem,
// because only there is there no equilibrium iteration AND no mass source (see the warning
// text). Enumerated from SRC/classTags.h rather than guessed; unknown transient tags are
// treated as implicit (say nothing) so a new integrator never inherits a false warning.
static bool ladrunoIsExplicitIntegratorTag(int tag)
{
  switch (tag) {
  case INTEGRATOR_TAGS_CentralDifference:                 // 5
  case INTEGRATOR_TAGS_CentralDifferenceAlternative:      // 17
  case INTEGRATOR_TAGS_CentralDifferenceNoDamping:        // 18
  case INTEGRATOR_TAGS_ExplicitDifference:                // 55
  case INTEGRATOR_TAGS_ExplicitBathe:                     // 33000 (+ its collapsed aliases)
  case INTEGRATOR_TAGS_ExplicitDifferenceStatic:          // 33001
  case INTEGRATOR_TAGS_ExplicitBatheLNVD:                 // 33002
  case INTEGRATOR_TAGS_CentralDifferenceLadruno:          // 33003
  case INTEGRATOR_TAGS_CentralDifferenceSMS:              // 33007
  case INTEGRATOR_TAGS_CentralDifferenceSMSConsistent:    // 33008
  case INTEGRATOR_TAGS_ExplicitBatheSMS:                  // 33009
  case INTEGRATOR_TAGS_ExplicitBatheSMSConsistent:        // 33010
  case INTEGRATOR_TAGS_ExplicitBatheLNVDSMS:              // 33011
  case INTEGRATOR_TAGS_ExplicitBatheLNVDSMSConsistent:    // 33012
    return true;
  default:
    return false;
  }
}

// Ladruno (WP-101 r1, narrowed in r2): one-time note that -enforce al cannot work under an
// EXPLICIT integrator. A PARSE-TIME refusal is impossible -- the integrator is unknown when
// the element is declared -- so the check lives at the first update(), where the active
// integrator is finally visible. Warn, do not refuse: the combination is pre-existing.
//
// IMPLICIT transient (Newmark, HHT, GeneralizedAlpha, ...) is NOT warned about: it has
// equilibrium iterations, so the commit-cadence Uzawa behaves exactly as it does in statics.
void LadrunoKinematicCoupling::warnAlUnderTransient(void)
{
  if (alWarnedTransient || enforce != 1) return;
  TransientIntegrator** tiPtr = OPS_GetTransientIntegrator();
  StaticIntegrator** siPtr = OPS_GetStaticIntegrator();
  if (tiPtr == 0 || *tiPtr == 0) return;
  if (siPtr != 0 && *siPtr != 0) return;             // a static analysis is the active one
  if (!ladrunoIsExplicitIntegratorTag((*tiPtr)->getClassTag())) return;   // implicit: fine
  alWarnedTransient = true;
  opserr << "WARNING LadrunoKinematicCoupling " << this->getTag()
         << ": -enforce al has no effect under an EXPLICIT integrator and leaves the tie "
         << "without a mass source. The Uzawa update needs equilibrium iterations to converge "
         << "against, and -bipenalty is DROPPED when -enforce al is given, so a massless tied "
         << "DOF gets no penalty mass -- measured under CentralDifferenceLadruno: "
         << "`-enforce penalty -bipenalty -dtcr <dt>` runs, while `-enforce al -bipenalty` "
         << "fails at step 0. Use -enforce penalty with -bipenalty for explicit runs. "
         << "(Implicit transient -- Newmark, HHT, GeneralizedAlpha -- is fine and is not "
         << "warned about.) See LadrunoKinematicCoupling_guide.md section 4.5\n";
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
  static Vector hdr(21);        // Ladruno (WP-101): +alUpdate at slot 20
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
  hdr(19) = 2.0;                          // version (2 = +alUpdate, +lambdaCommitted)
  hdr(20) = alUpdate;                     // Ladruno (WP-101)
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
  // payload: lambdaAL(nGap) + g0(nGap) + lambdaCommitted(nGap). nGap is recomputed
  // deterministically in setDomain on the recv side; the ordering matches, so these align
  // by index. Ladruno (WP-101): lambdaCommitted rides along or a partition that receives
  // mid-step would revert to a zero multiplier on the next failed step.
  if (nGap > 0) {
    Vector payload(3 * nGap);
    for (int r = 0; r < nGap; r++)
      payload(r) = (lambdaAL.Size() == nGap) ? lambdaAL(r) : 0.0;
    for (int r = 0; r < nGap; r++)
      payload(nGap + r) = (g0.Size() == nGap) ? g0(r) : 0.0;
    for (int r = 0; r < nGap; r++)
      payload(2 * nGap + r) = (lambdaCommitted.Size() == nGap) ? lambdaCommitted(r)
                            : ((lambdaAL.Size() == nGap) ? lambdaAL(r) : 0.0);
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
  static Vector hdr(21);        // Ladruno (WP-101): must match sendSelf's size
  if (theChannel.recvVector(dbTag, commitTag, hdr) < 0) {
    opserr << "LadrunoKinematicCoupling::recvSelf - header failed\n";
    return -1;
  }
  // Ladruno (WP-101 r1): the version field means something now -- refuse a payload written
  // by a NEWER layout rather than mis-reading it.
  if (hdr(19) > 2.0) {
    opserr << "LadrunoKinematicCoupling::recvSelf - payload version " << hdr(19)
           << " is newer than this build understands (2); refusing\n";
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
  alUpdate = (int)hdr(20);                 // Ladruno (WP-101)
  nrot = (ndm == 3) ? 3 : 1;
  ktResolved = false;
  ktHostMissWarned = false;
  bpResolved = false;
  alGuardWarned = false;
  alWarnedTransient = false;
  // Ladruno (WP-101 r1): ARM the latch, do not clear it. Domain::recv() calls theEle->update()
  // immediately after recvSelf, which under `iter` advanced lambda on the just-restored state
  // (measured: a database save/restore moved lambda by 6.6e-9). Same one-shot contract as the
  // post-revert call.
  alSkipUpdate = true;
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
  lambdaCommitted.resize(nGap > 0 ? nGap : 1); lambdaCommitted.Zero();
  g0.resize(nGap > 0 ? nGap : 1);       g0.Zero();
  if (nGap > 0) {
    Vector payload(3 * nGap);
    if (theChannel.recvVector(dbTag, commitTag, payload) < 0) {
      opserr << "LadrunoKinematicCoupling::recvSelf - payload failed\n";
      return -1;
    }
    for (int r = 0; r < nGap; r++) lambdaAL(r) = payload(r);
    for (int r = 0; r < nGap; r++) g0(r) = payload(nGap + r);
    for (int r = 0; r < nGap; r++) lambdaCommitted(r) = payload(2 * nGap + r);
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
  if (enforce == 1)                          // Ladruno (WP-101)
    s << " (alUpdate: " << (alUpdate == 1 ? "iter" : "commit") << ")";
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
  if (LadrunoResp::is(argv[0], "force"))
    return new ElementResponse(this, 1, Vector(n));
  if (LadrunoResp::is(argv[0], "gap"))
    return new ElementResponse(this, 2, Vector(n));
  if (LadrunoResp::is(argv[0], "penalty"))      // kt / k / penalty
    return new ElementResponse(this, 3, 0.0);
  if (LadrunoResp::is(argv[0], "rotPenalty"))   // kr
    return new ElementResponse(this, 4, 0.0);
  if (LadrunoResp::is(argv[0], "lambda"))
    return new ElementResponse(this, 5, Vector(n));
  if (LadrunoResp::is(argv[0], "dtCritical"))
    return new ElementResponse(this, 6, 0.0);
  if (LadrunoResp::is(argv[0], "tiedDOFs"))
    return new ElementResponse(this, 7, 0.0);
  // Ladruno — the tie-family audit pair, present on LadrunoEmbeddedNode/Rebar
  // but missing here: the ARTIFICIAL penalty energy currently stored (net it
  // out of a global EnergyBalance) and the constraint violation ||g||.
  // NOTE: no massPenalty response — this element lumps a per-DOF DIAGONAL
  // bipenalty mass (*M0), not the single scalar the rest of the family reports.
  if (LadrunoResp::is(argv[0], "penaltyEnergy"))
    return new ElementResponse(this, 8, 0.0);
  if (LadrunoResp::is(argv[0], "constraintViolation"))
    return new ElementResponse(this, 9, 0.0);
  return this->Element::setResponse(argv, argc, s);
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
  case 8: {   // artificial penalty energy 1/2 sum_row k_row g_row^2
    Vector g(nGap); this->computeGap(g);
    double E = 0.0;
    for (int row = 0; row < nGap; row++)
      E += this->rowPenalty(row) * g(row) * g(row);
    return eleInfo.setDouble(0.5 * E);
  }
  case 9: {   // constraint violation ||g||
    Vector g(nGap); this->computeGap(g);
    double g2 = 0.0;
    for (int row = 0; row < nGap; row++) g2 += g(row) * g(row);
    return eleInfo.setDouble(sqrt(g2));
  }
  default: return this->Element::getResponse(responseID, eleInfo);
  }
}
