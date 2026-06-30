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

// ADR-39 P1a — LadrunoContactHandler implementation. See header for design.
// Replicates the PlainHandler DOF_Group/FE loop (mirrors the proven in-fork
// LadrunoProjectionHandler) and appends the contact FE adapter(s).

#include "LadrunoContactHandler.h"
#include "LadrunoContactFE.h"

#include <LadrunoContactDomain.h>    // Ladruno: ADR-39 (adapter count from the engine)
#include <LadrunoContactSurface.h>   // Ladruno: ADR-39 P2a (slave node-set)
#include <LadrunoContactKernel.h>    // Ladruno: ADR-39 P2b-2b (reference normal for -kn auto)
#include <LadrunoContactProjection.h> // Ladruno: ADR-41 A2 (normalOriented projection geometry)
#include <LadrunoContactBucketSort.h>// Ladruno: ADR-39 P2.5 (broad-phase pairing)
#include <LadrunoContactReemit.h>   // Ladruno: ADR-60 (finite-sliding re-emit band/metric)
#include <LadrunoEdgeKernel.h>       // Ladruno: ADR-57 E2 (edge-edge routing + closest point)
#include <Matrix.h>                  // Ladruno: ADR-39 P2b-2b (master getInitialStiff)
#include <Domain.h>
#include <AnalysisModel.h>
#include <Integrator.h>
#include <Node.h>
#include <Vector.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <FE_Element.h>
#include <Subdomain.h>
#include <ID.h>
#include <classTags.h>
#include <elementAPI.h>
#include <map>
#include <set>
#include <vector>
#include <cmath>

// ----------------------------------------------------------------------------
// command factory: constraints LadrunoContact
// ----------------------------------------------------------------------------
void *OPS_LadrunoContactHandler(void)
{
    return new LadrunoContactHandler();
}

// ----------------------------------------------------------------------------
LadrunoContactHandler::LadrunoContactHandler()
  : ConstraintHandler(HANDLER_TAG_LadrunoContactHandler)
{
}

LadrunoContactHandler::~LadrunoContactHandler()
{
}

// P2b-2b: auto-size the penalty kn for ONE (slave, master-segment) pair from the
// OWNING solid element's initial stiffness, sourced GENERICALLY through base-Element
// virtuals (getExternalNodes / getNumDOF / getInitialStiff) so any 3-DOF/node solid
// works with NO element-type coupling and NO vanilla edit. Reduction (validated in
// proto_p2b2b_autokn.py): kn = f_si * mean_over_seg_nodes( nᵀ K_block_node n ), where
// K_block_node is the 3x3 diagonal block at that node's DOFs and n is the reference
// segment normal. This is the element-stiffness form of the LS-DYNA 26.14a penalty
// f·K·A²/V (K_diag ~ E·L ~ E·A²/V), the A²/V geometry absorbed exactly into the
// assembled matrix. f_si = 0.10 (LS-DYNA SLSFAC default). Returns <= 0 on failure
// (no owning solid found / non-3-DOF element / ambiguous reference normal) so the
// caller skips the pair with a warning. The SOFT Courant floor (26.15) is P2b-2c.
static double
ladrunoResolveAutoKn(Domain *theDomain, Node **segNodes, int nps,
                     const double orientDir[3])
{
    const double F_SI = 0.10;
    // reference (undeformed) segment node coords
    double Xref[4][3];
    for (int k = 0; k < nps; k++) {
        const Vector &Xk = segNodes[k]->getCrds();
        for (int d = 0; d < 3; d++) Xref[k][d] = Xk(d);
    }
    // reference outward normal at the segment parametric center, oriented by orientDir
    // (winding-immune; refuses an ambiguous in-plane reference direction).
    double xiC  = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double etaC = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double n[3];
    if (!LadrunoContactProjection::normalOriented(nps, xiC, etaC, Xref, orientDir, n))
        return -1.0;
    int segTags[4];
    for (int k = 0; k < nps; k++) segTags[k] = segNodes[k]->getTag();
    // owning solid = the first Domain element that (a) carries 3 translational DOFs
    // per node and (b) contains ALL the segment's nodes. Brute force over the Domain
    // elements (fine for the gate meshes; the bucket-sort broad phase is P2.5).
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0 || e->getNumDOF() != 3 * nn) continue;   // not a 3-DOF/node solid
        int loc[4]; bool allIn = true;
        for (int k = 0; k < nps; k++) {
            loc[k] = en.getLocation(segTags[k]);
            if (loc[k] < 0) { allIn = false; break; }
        }
        if (!allIn) continue;
        const Matrix &K = e->getInitialStiff();
        if (K.noRows() != 3 * nn || K.noCols() != 3 * nn) return -1.0;
        double acc = 0.0;
        for (int k = 0; k < nps; k++) {
            int c = 3 * loc[k];
            double Kn[3];
            for (int i = 0; i < 3; i++)
                Kn[i] = K(c + i, c + 0) * n[0] + K(c + i, c + 1) * n[1] + K(c + i, c + 2) * n[2];
            acc += n[0] * Kn[0] + n[1] * Kn[1] + n[2] * Kn[2];
        }
        double kn = F_SI * acc / nps;
        return (kn > 0.0) ? kn : -1.0;
    }
    return -1.0;   // no owning solid element found for this segment
}

// B1 (P4) — build the ASSEMBLED translational nodal-mass cache the SOFT=1 penalty needs. The explicit
// integrator inverts the assembled global diagonal M = Σ_elements diag(M_e) + the nodal `mass`, but
// Node::getMass() holds ONLY the nodal `mass` (element-density mass never reaches the node). So here we
// reconstruct, per node, m[d] = nodal mass(d) + Σ_e diag(M_e) at that node's translational DOFs (the
// SAME diagonal `system Diagonal` inverts) and stash it on the engine for the stateless adapter to
// read. One pass over the elements; called only when a SOFT contact exists (masses are constant
// between domain changes, so rebuild-per-handle is correct). Translation-first DOF order (OpenSees
// nodes order [u | θ]) ⇒ a node's first 3 element DOFs are translational.
static void
ladrunoBuildNodalMass(Domain *theDomain, LadrunoContactDomain *cd)
{
    cd->clearNodalMass();
    // seed every node with its nodal `mass` (the Node::getMass diagonal, translational DOFs)
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0) {
        const Matrix &M = nodPtr->getMass();
        int nr = M.noRows();
        double m[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < 3 && d < nr; d++) m[d] = M(d, d);
        cd->setNodalMass(nodPtr->getTag(), m);
    }
    // add each element's diagonal mass contribution to its nodes
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0) continue;
        int ndof = e->getNumDOF();
        if (ndof <= 0 || ndof % nn != 0) continue;          // non-uniform DOF layout ⇒ cannot map
        int dpn = ndof / nn;                                // DOFs per node (3 solid, 6 beam/shell)
        const Matrix &Me = e->getMass();
        if (Me.noRows() != ndof || Me.noCols() != ndof) continue;   // no/odd mass matrix ⇒ skip
        for (int k = 0; k < nn; k++) {
            double m[3];
            if (!cd->getNodalMass(en(k), m)) { m[0] = m[1] = m[2] = 0.0; }
            for (int d = 0; d < 3 && d < dpn; d++) m[d] += Me(k * dpn + d, k * dpn + d);
            cd->setNodalMass(en(k), m);
        }
    }
}

// ============================================================ ADR-57 E2 edge-edge routing geometry
// Handler-local helpers porting the proto_e1_router oracle (16/16): the facet-normal align_fn regime
// split + the proximity-gated NTS ≻ edge-edge ≻ face-mortar precedence. Built from REFERENCE coords
// at handle() (the §7 scoping — large sliding needs re-emission, the inherited subsystem fence).

// Newell unit normal of a planar convex tri/quad facet (orientation-independent; align_fn = |nS·nM|).
static void ladrunoFacetNormalNewell(int npts, const double X[4][3], double n[3])
{
    n[0] = n[1] = n[2] = 0.0;
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nn > 1e-300) { n[0] /= nn; n[1] /= nn; n[2] /= nn; }
}

// does point x project in-bounds onto the convex facet X (npts) with |gap| ≤ band? gap by ref.
static bool ladrunoPointInFacet(const double x[3], int npts, const double X[4][3],
                                const double nf[3], double band, double &gap)
{
    double c[3] = {0, 0, 0};
    for (int i = 0; i < npts; i++) for (int d = 0; d < 3; d++) c[d] += X[i][d] / npts;
    gap = 0.0;
    for (int d = 0; d < 3; d++) gap += (x[d] - c[d]) * nf[d];
    double proj[3];
    for (int d = 0; d < 3; d++) proj[d] = x[d] - gap * nf[d];
    int sgn = 0;                                   // consistent cross-product sign vs each edge (convex)
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        double e[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        double w[3] = {proj[0]-a[0], proj[1]-a[1], proj[2]-a[2]};
        double cr[3] = {e[1]*w[2]-e[2]*w[1], e[2]*w[0]-e[0]*w[2], e[0]*w[1]-e[1]*w[0]};
        double dp = cr[0]*nf[0] + cr[1]*nf[1] + cr[2]*nf[2];
        int s = (dp > 1e-300) ? 1 : ((dp < -1e-300) ? -1 : 0);
        if (s != 0) { if (sgn == 0) sgn = s; else if (s != sgn) return false; }
    }
    return std::fabs(gap) <= band;
}

// min distance between two facets: vertices-projected-in-other + all boundary edge-edge closest
// points (the proto_e1 min_facet_dist). Reference coords ⇒ the proximity gate.
static double ladrunoMinFacetDist(int npS, const double S[4][3], const double nS[3],
                                  int npM, const double M[4][3], const double nM[3])
{
    double dmin = HUGE_VAL, gap;
    for (int i = 0; i < npS; i++)
        if (ladrunoPointInFacet(S[i], npM, M, nM, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int i = 0; i < npM; i++)
        if (ladrunoPointInFacet(M[i], npS, S, nS, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int ie = 0; ie < npS; ie++)
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            double g = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            if (g < dmin) dmin = g;
        }
    return dmin;
}

// d_band for a master facet: explicit -edgeBand override, else a fraction of the mean edge length.
static double ladrunoEdgeBand(int npM, const double M[4][3], double override_)
{
    if (override_ > 0.0) return override_;
    double el = 0.0;
    for (int k = 0; k < npM; k++) {
        const double *a = M[k]; const double *b = M[(k + 1) % npM];
        double dx = b[0]-a[0], dy = b[1]-a[1], dz = b[2]-a[2];
        el += std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    return 0.25 * (el / npM);
}

// The SINGLE SOURCE OF TRUTH for "edge-edge OWNS this (slave facet, master facet) pair" — the
// proto_e1 routing decision: in proximity, NON-FACING (align_fn < τ_face), NOT owned by an NTS
// slave-vertex-on-face poke, AND ≥1 boundary-edge pair is a well-conditioned margin-interior in-band
// crossing. Called by BOTH the mortar injection loop (which STANDS DOWN when this is true — the ADR §2
// "face-mortar vacates only if edge-edge succeeds", so the mortar clip and the edge-edge force never
// superpose in the oblique band where the clip is non-degenerate) AND the edge-edge injection loop
// (which injects only when this is true). One decision, two consumers ⇒ no double-count (the gate's
// EE-1 fix). Reference coords (§7 scoping).
static bool ladrunoEdgeEdgeOwns(int npS, const double S[4][3], int npM, const double M[4][3],
                                double dBand)
{
    const double TAU_FACE = std::cos(50.0 * 3.14159265358979323846 / 180.0);
    double nS[3], nM[3];
    ladrunoFacetNormalNewell(npS, S, nS);
    ladrunoFacetNormalNewell(npM, M, nM);
    if (ladrunoMinFacetDist(npS, S, nS, npM, M, nM) > dBand) return false;   // separated ⇒ NONE
    double alignFn = std::fabs(nS[0]*nM[0] + nS[1]*nM[1] + nS[2]*nM[2]);
    if (alignFn >= TAU_FACE) return false;                                   // facing ⇒ FACE-MORTAR owns
    double gap;
    for (int k = 0; k < npS; k++)                                           // slave vertex poke ⇒ NTS owns
        if (ladrunoPointInFacet(S[k], npM, M, nM, dBand, gap)) return false;
    for (int ie = 0; ie < npS; ie++)                                        // ≥1 in-band edge crossing
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            if (std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) <= dBand) return true;
        }
    return false;
}

int
LadrunoContactHandler::handle(const ID *nodesLast)
{
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();
    if (theDomain == 0 || theModel == 0 || theIntegrator == 0) {
        opserr << "WARNING LadrunoContactHandler::handle() - setLinks() not called\n";
        return -1;
    }

    // P1a: MP constraints (equalDOF / rigidLink / rigidDiaphragm) are NOT enforced
    // by the contact handler yet (Plain-style numbering). Warn rather than silently
    // ignore. Delegation to a base handler lands in P1b.
    MP_ConstraintIter &theMPcheck = theDomain->getMPs();
    if (theMPcheck() != 0)
        opserr << "WARNING LadrunoContactHandler::handle() - MP constraints (e.g. equalDOF) "
                  "are present but NOT enforced by the P1a contact handler; use a model with "
                  "SP (fix) constraints only, or 'constraints Transformation' for those models.\n";

    // collect SPs (domain + load-pattern), keyed by node. Like PlainHandler this is a
    // Plain-style handler: it REMOVES a constrained DOF but cannot enforce a NON-homogeneous
    // (imposed-displacement) SP value — that needs the Transformation/Penalty/Lagrange family,
    // which this handler can't be (it must inject the contact FE adapters in handle()). Warn
    // the same way PlainHandler does so a displacement-driven contact model fails LOUD, not
    // silently-zero; point the user to the DisplacementControl integrator (works with this
    // handler — it drives a DOF via the load factor, no imposed SP needed). (B3 follow-up.)
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0) {
        if (theSP->isHomogeneous() == false)
            opserr << "WARNING LadrunoContactHandler::handle() - non-homogeneous SP (imposed "
                      "displacement) at node " << theSP->getNodeTag() << " is NOT enforced by "
                      "the contact handler (Plain-style; homogeneous constraint assumed). Use the "
                      "DisplacementControl integrator for displacement-driven contact.\n";
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));
    }

    // DOF_Groups: free = -2, single-point-constrained = -1.
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    int numDOF = 0, countDOF = 0, count3 = 0;
    while ((nodPtr = theNod()) != 0) {
        DOF_Group *dofPtr = new DOF_Group(numDOF++, nodPtr);
        if (dofPtr == 0) {
            opserr << "WARNING LadrunoContactHandler::handle() - out of memory (DOF_Group)\n";
            return -4;
        }
        const ID &id = dofPtr->getID();
        for (int j = 0; j < id.Size(); j++) {
            dofPtr->setID(j, -2);
            countDOF++;
        }
        int nodeID = nodPtr->getTag();
        std::multimap<int, SP_Constraint *>::iterator first = allSPs.lower_bound(nodeID);
        std::multimap<int, SP_Constraint *>::iterator last = allSPs.upper_bound(nodeID);
        for (std::multimap<int, SP_Constraint *>::iterator it = first; it != last; it++) {
            int dof = it->second->getDOF_Number();
            const ID &cid = dofPtr->getID();
            if (dof >= 0 && dof < cid.Size() && cid(dof) == -2) {
                dofPtr->setID(dof, -1);
                countDOF--;
            } else {
                opserr << "WARNING LadrunoContactHandler::handle() - invalid/duplicate SP at DOF "
                       << dof << " node " << nodeID << endln;
            }
        }
        nodPtr->setDOF_GroupPtr(dofPtr);
        theModel->addDOF_Group(dofPtr);
    }
    theModel->setNumEqn(countDOF);

    // nodesLast (subdomain boundary): mark -3 (as PlainHandler)
    if (nodesLast != 0) {
        for (int i = 0; i < nodesLast->Size(); i++) {
            Node *np = theDomain->getNode((*nodesLast)(i));
            if (np != 0) {
                DOF_Group *dp = np->getDOF_GroupPtr();
                const ID &id = dp->getID();
                for (int j = 0; j < id.Size(); j++)
                    if (id(j) == -2) { dp->setID(j, -3); count3++; }
            }
        }
    }

    // FE_Elements for every Domain element (as PlainHandler). Track numFe so the
    // contact adapter tags continue ABOVE the structural FE tags (addFE_Element
    // silently drops a duplicate tag).
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    int numFe = 0;
    // ADR-60 R6 (D7 serial-only refusal): detect a partitioned host — this Domain is itself a
    // worker Subdomain (dynamic_cast) OR it owns Subdomain elements (the SP coordinator). Under
    // either, the finite-sliding migration trigger would re-sort uncoordinated across ranks, so
    // re-emit reverts to the shipped frozen-config NTS broad phase below (a one-time warning).
    // Sequential ⇒ false ⇒ every re-emit branch is inert ⇒ byte-identical to the shipped path.
    bool hostPartitioned = (dynamic_cast<Subdomain *>(theDomain) != 0);
    while ((elePtr = theEle()) != 0) {
        if (elePtr->isSubdomain() == true) {
            hostPartitioned = true;
            Subdomain *theSub = (Subdomain *)elePtr;
            if (theSub->doesIndependentAnalysis() == false) {
                FE_Element *fePtr = new FE_Element(numFe++, elePtr);
                if (fePtr == 0) return -5;
                theModel->addFE_Element(fePtr);
                theSub->setFE_ElementPtr(fePtr);
            }
        } else {
            FE_Element *fePtr = new FE_Element(numFe++, elePtr);
            if (fePtr == 0) return -5;
            theModel->addFE_Element(fePtr);
        }
    }

    // --- inject the contact FE adapter(s) ---
    // None if no contact engine is attached (pure Plain -> byte-identical to stock).
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();

    // B1 (P4): if ANY NTS contact / rigid plane requests the SOFT=1 penalty, pre-compute the
    // assembled translational nodal-mass cache (nodal `mass` + element-density mass) the adapters
    // read to size k_soft = SOFSCL·4·m_eff/dt². Only when soft is in play ⇒ no cost (and byte-
    // identical) otherwise. Mortar SOFT is out of scope (refused at the command surface).
    if (cd != 0) {
        bool anySoft = false;
        for (int c = 0; c < cd->getNumContacts() && !anySoft; c++)
            if (cd->getContact(c).softScale > 0.0) anySoft = true;
        for (int p = 0; p < cd->getNumRigidPlanes() && !anySoft; p++)
            if (cd->getRigidPlane(p).softScale > 0.0) anySoft = true;
        // B2 (P5): a SOFT=2 mortar contact also needs the assembled nodal-mass cache (the segment-
        // based m_eff). The cache is over ALL nodes (one element pass), shared with the NTS SOFT=1 lane.
        // ADR-57 E5: a -edgeSoft edge-edge contact needs the SAME cache (the edge gap-mode m_eff).
        for (int c = 0; c < cd->getNumMortarContacts() && !anySoft; c++)
            if (cd->getMortarContact(c).softScale > 0.0 ||
                cd->getMortarContact(c).edgeSoftScale > 0.0) anySoft = true;
        if (anySoft)
            ladrunoBuildNodalMass(theDomain, cd);
    }

    // P2b: faceted node-to-segment penalty contact. For each contact definition
    // (MASTER_SEGMENTS surface vs SLAVE_NODES surface) build ONE bound adapter per
    // (slave node, master segment) pair — the gate-sanctioned brute-force pairing
    // (bucket-sort broad phase = P2.5). A degenerate/non-projecting segment yields
    // zero force inside the adapter (so a topology-only contact stays graph-neutral
    // under a connectivity-independent numberer — the P1b regression relies on this).
    if (cd != 0) {
        // P3: rebuild the live friction-state key-set this handle() so the engine can
        // GC slots from a previous (re-meshed/re-paired) analysis. Mark every built
        // frictional pair; sweep at the end.
        cd->frictionGCBegin();
        cd->clearReemit();   // ADR-60: drop last epoch's re-emit anchors (re-registered per opted-in contact below)

        // P3 cross-contact shared-slave warning: a slave node in two contacts' slave
        // sets gets two adapters → double traction (a modeling error the engine can't
        // auto-resolve). Cheap multiset scan, mirrors the ndf!=3 warning.
        {
            std::multimap<int, int> slaveSeen;   // slaveNodeTag -> contactTag
            for (int c = 0; c < cd->getNumContacts(); c++) {
                const LadrunoContactDomain::Contact &ct = cd->getContact(c);
                LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
                if (ss == 0 || ss->getKind() != LadrunoContactSurface::SLAVE_NODES)
                    continue;
                const ID &st = ss->getNodeTags();
                for (int i = 0; i < st.Size(); i++) {
                    int tag = st(i);
                    if (slaveSeen.find(tag) != slaveSeen.end())
                        opserr << "WARNING LadrunoContactHandler::handle() - slave node "
                               << tag << " appears in more than one contact's slave set; "
                                  "its contact tractions will be ADDED (check the model)\n";
                    slaveSeen.insert(std::make_pair(tag, ct.tag));
                }
            }
        }

        for (int c = 0; c < cd->getNumContacts(); c++) {
            const LadrunoContactDomain::Contact &ct = cd->getContact(c);
            LadrunoContactSurface *ms = cd->getSurface(ct.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": master/slave surface not defined; skipped\n";
                continue;
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_NODES) {
                opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": need a MASTER_SEGMENTS master + SLAVE_NODES slave; skipped\n";
                continue;
            }
            int nps = ms->getNodesPerSeg();
            if (nps < 3 || nps > 4) {
                opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": nodesPerSeg " << nps << " unsupported (need 3 or 4); skipped\n";
                continue;
            }
            if (!ct.knAuto && ct.kn <= 0.0) {
                // A SEGMENT contact needs a positive penalty (P2b-1 requires -kn).
                // kn == 0 (e.g. `contact ... -outward` with the kn omitted) is inert;
                // warn rather than silently build dead adapters. (Gate PARSE-2/H1.)
                // `-kn auto` (P2b-2b) carries a 0 placeholder and is resolved per pair
                // below, so it bypasses this guard.
                opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": segment contact needs kn > 0 (got " << ct.kn << "); skipped\n";
                continue;
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSeg = mTags.Size() / nps;

            // ADR-60 R6 (serial-only): one-time warning when -reemit is requested under a
            // partitioned host. The deformed feed + trigger are disabled below (reemitActive=false
            // ⇒ the shipped frozen-config NTS path), so this only surfaces the silent downgrade.
            if (ct.enableReemit && hostPartitioned) {
                static bool warnedReemitPartition = false;
                if (!warnedReemitPartition) {
                    opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": -reemit is serial-only but the host is partitioned (Subdomain); "
                              "finite-sliding re-emit DISABLED (reverting to the frozen-config NTS "
                              "broad phase) to avoid an uncoordinated cross-rank re-sort.\n";
                    warnedReemitPartition = true;
                }
            }
            // ADR-60 R1 (BLOCKER-MEMBERSHIP): fingerprint this contact's master mTags ORDERING. If
            // it changed since the last handle (element removal / re-mesh reordered the surface),
            // the segment-ordinal friction key would alias onto a different physical segment — drop
            // this contact's friction slots (they re-engage fresh, traction-continuous per D4) and
            // refuse to ARM the trigger this step (named error). Only evaluated when -reemit is on
            // ⇒ OFF leaves theReemitFp untouched ⇒ byte-identical.
            bool membershipChanged = false;
            if (ct.enableReemit) {
                std::vector<int> mt((size_t)mTags.Size());
                for (int i = 0; i < mTags.Size(); i++) mt[i] = mTags(i);
                unsigned long long fp = LadrunoContactReemit::membershipFingerprint(
                    mt.empty() ? 0 : &mt[0], (int)mt.size(), nps);
                if (cd->reemitMembershipChanged(ct.tag, fp)) {
                    membershipChanged = true;
                    opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": master surface node ordering changed since the last handle; the "
                              "segment-keyed re-emit friction history is invalid. Dropping this "
                              "contact's friction state (re-engages fresh) and refusing to arm "
                              "re-emit this step.\n";
                    cd->dropFrictionForContact(ct.tag);
                }
            }
            // The deformed broad-phase feed is active when -reemit is on AND the host is serial
            // (R6). The migration TRIGGER additionally needs a STABLE master membership this
            // handle (R1) — under a membership change we keep the deformed feed (no pass-through)
            // but do not arm an automatic re-handle off the just-changed mesh.
            bool reemitActive = ct.enableReemit && !hostPartitioned;
            bool reemitArmed  = reemitActive && !membershipChanged;

            // P2.5 BROAD PHASE: bucket-sort the master segments (reference coords)
            // once, then per slave query only the 27-neighbour candidate segments
            // instead of all nSeg (brute force). The kept set is a SUPERSET of every
            // near pair, so the contact result is identical to brute force; a huge
            // -cell (=> 1 bucket) reproduces brute force exactly (the equivalence
            // gate). Missing master nodes (malformed model) are filled with the
            // segment's first valid coord — superset-preserving, and the narrow loop
            // below re-fetches + skips them anyway.
            // segCoords is filled with REFERENCE coords first (the shipped P0 path). When re-emit is
            // on it is shifted to the DEFORMED config IN PLACE below — but the R2 search band is taken
            // from these REFERENCE coords (deformation-invariant; it must not drift vs the grid the
            // next re-emit builds — gate D6 + the referenceBand() contract).
            std::vector<double> segCoords((size_t)nSeg * nps * 3, 0.0);
            for (int seg = 0; seg < nSeg; seg++) {
                double *S = &segCoords[(size_t)seg * nps * 3];
                bool haveFirst = false; double first[3] = {0.0, 0.0, 0.0};
                for (int k = 0; k < nps; k++) {
                    Node *mn = theDomain->getNode(mTags(seg * nps + k));
                    if (mn != 0) {
                        const Vector &X = mn->getCrds();
                        S[k*3+0] = X(0); S[k*3+1] = X(1); S[k*3+2] = X(2);
                        if (!haveFirst) { first[0]=X(0); first[1]=X(1); first[2]=X(2); haveFirst = true; }
                    } else {
                        S[k*3+0] = S[k*3+1] = S[k*3+2] = HUGE_VAL;   // mark missing
                    }
                }
                for (int k = 0; k < nps; k++)   // backfill missing with first valid
                    if (S[k*3+0] == HUGE_VAL)
                        for (int d = 0; d < 3; d++) S[k*3+d] = first[d];
            }

            // ADR-60 R2: the deformation-invariant re-emit search band = reference median segment
            // diagonal × cellFrac, from the REFERENCE segCoords ABOVE (before the deformed shift).
            double reemitBand = reemitActive
                ? LadrunoContactReemit::referenceBand(nSeg, nps, segCoords.data(), ct.cellFrac) : 0.0;

            // ADR-60 P1: when re-emit is on, shift segCoords to the DEFORMED (committed) config IN
            // PLACE so the rebuilt candidate set tracks the slid geometry. OFF ⇒ segCoords stays the
            // reference value (the verbatim shipped path) ⇒ byte-identical. A missing-node slot keeps
            // its backfilled reference value (the narrow loop skips that pair anyway).
            if (reemitActive) {
                for (int seg = 0; seg < nSeg; seg++)
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        if (mn == 0) continue;
                        const Vector &u = mn->getDisp();
                        double *S = &segCoords[((size_t)seg * nps + k) * 3];
                        S[0] += u(0); S[1] += u(1); S[2] += u(2);
                    }
            }

            // ADR-60 P1 (BLOCKER-CLIP): the runaway percentile clip targets a STATIC reference
            // outlier; on the deformed feed a legitimately-diverging node would collapse the grid and
            // silently drop real pairs. Disable it (clipPct=0) when re-emit is on — the
            // NX*NY*NZ≤min(nSeg,5000) cap still bounds memory. OFF ⇒ the shipped 1.0 clip (identical).
            // ADR-60 R8: Grid::runawayGuardFired() is intentionally NOT auto-warned. As implemented the
            // guard fires whenever the 1/99-percentile clamp moves a bound — i.e. for ANY mesh with >100
            // segment-centroids (the tails are clipped by design) — so it is not a clean "a node ran away"
            // signal and surfacing it would spam normal large models. The safety Finding-3 cared about (the
            // grid never collapsing during instability) is provided by clipPct=0 here, not by a warning;
            // the accessor stays available for debugging/tests.
            LadrunoContactBucketSort::Grid grid(nSeg, nps, segCoords.data(), ct.cellFrac,
                                                reemitActive ? 0.0 : 1.0);
            std::vector<int> cand(nSeg);

            // ADR-60: register this contact for finite-sliding re-emit (opt-in via -reemit) with the
            // R2 deformation-invariant reference band. Gated on reemitArmed: serial host (R6) AND a
            // stable master membership this handle (R1). OFF ⇒ nothing registered ⇒ byte-identical.
            if (reemitArmed)
                cd->setReemitContact(ct.tag, reemitBand, ct.resortFrac, ct.resortEvery);

            for (int si = 0; si < sTags.Size(); si++) {
                Node *sn = theDomain->getNode(sTags(si));
                if (sn == 0 || sn->getNumberDOF() != 3) {
                    if (sn != 0)
                        opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                               << " slave node " << sTags(si) << " ndf=" << sn->getNumberDOF()
                               << " != 3; skipped (P2b is 3D translational)\n";
                    continue;
                }
                const Vector &Xs0 = sn->getCrds();
                double slavePt[3] = { Xs0(0), Xs0(1), Xs0(2) };
                if (reemitActive) {
                    // ADR-60 P1: query the grid at the slave's DEFORMED (committed) position so a slid
                    // slave finds the segment it is ACTUALLY over (the no-pass-through fix); the SAME
                    // deformed point is the re-emit anchor (drift resets each re-emit ⇒ sane cadence).
                    // OFF ⇒ slavePt stays Xs0 ⇒ byte-identical.
                    const Vector &us = sn->getDisp();
                    slavePt[0] = Xs0(0) + us(0); slavePt[1] = Xs0(1) + us(1); slavePt[2] = Xs0(2) + us(2);
                    // R1: register the anchor only when the trigger is armed (stable membership);
                    // on a membership-change step we keep the deformed query but do not arm a re-sort.
                    if (reemitArmed)
                        cd->addReemitAnchor(ct.tag, sTags(si), slavePt);
                }
                int nCand = grid.candidates(slavePt, cand.data());
                for (int ci = 0; ci < nCand; ci++) {
                    int seg = cand[ci];
                    Node *segNodes[4]; bool ok = true;
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        segNodes[k] = mn;
                    }
                    if (!ok) continue;   // missing/incompatible master node -> skip pair
                    // ADR-60 P1 (BLOCKER-SELF-CONTACT): finite-sliding re-emit can pair a slave with a
                    // segment that CONTAINS it (a folded declared surface). Skip a candidate whose node
                    // set includes the slave node — a self-penetration guard. Re-emit-only ⇒ identical.
                    if (reemitActive) {
                        bool selfSeg = false;
                        for (int k = 0; k < nps; k++)
                            if (mTags(seg * nps + k) == sTags(si)) { selfSeg = true; break; }
                        if (selfSeg) continue;
                    }
                    // orientation direction toward the slave's allowed half-space:
                    // explicit -outward if given, else auto = (slave ref) − (segment
                    // ref centroid). The derived normal is flipped to n·orientDir>0
                    // (winding-immune). Use -outward for just-penetrated starts.
                    double orientDir[3];
                    if (ct.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = ct.outward[d];
                    } else {
                        double cen[3] = {0.0, 0.0, 0.0};
                        for (int k = 0; k < nps; k++) {
                            const Vector &Xk = segNodes[k]->getCrds();
                            for (int d = 0; d < 3; d++) cen[d] += Xk(d);
                        }
                        const Vector &Xs = sn->getCrds();
                        for (int d = 0; d < 3; d++) orientDir[d] = Xs(d) - cen[d] / nps;
                    }
                    // P2b-2b: resolve `-kn auto` per (slave, segment) pair from the
                    // owning solid element's stiffness; a fixed `-kn $val` rides ct.kn.
                    double knUse = ct.kn;
                    if (ct.knAuto) {
                        knUse = ladrunoResolveAutoKn(theDomain, segNodes, nps, orientDir);
                        if (knUse <= 0.0) {
                            opserr << "WARNING LadrunoContactHandler::handle() - contact "
                                   << ct.tag << " slave node " << sTags(si) << " segment "
                                   << seg << ": -kn auto could not size a penalty (no owning "
                                      "solid element / non-3-DOF element / ambiguous normal); "
                                      "pair skipped\n";
                            continue;
                        }
                    }
                    // P3: thread kt/mu + the engine ptr + the friction-state key
                    // (contactTag, GLOBAL segment ordinal `seg`). mu<=0 ⇒ the adapter
                    // short-circuits friction (byte-identical to frictionless P2b).
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sn, segNodes, nps, knUse, orientDir,
                                             ct.kt, ct.mu, theDomain, ct.tag, seg,
                                             ct.consistentTan, ct.muc,       // D2 viscous (0 ⇒ off)
                                             ct.consistentNormal,            // B3 ∂n/∂u geom tangent
                                             ct.softScale);                  // B1 SOFT=1 (0 ⇒ off)
                    if (fe == 0) return -5;
                    theModel->addFE_Element(fe);
                    if (ct.mu > 0.0)
                        cd->frictionGCMark(ct.tag, sTags(si), seg);  // live this handle()
                }
            }
        }

        cd->frictionGCEnd();   // prune friction slots no live pair referenced
    }

    // --- ADR-41 C2.1: MORTAR pairing -> ONE clipped-GP penalty adapter per (slave facet,
    //     candidate master facet). Mirrors the NTS loop but FACET-vs-FACET: broad-phase
    //     the master facets, query per slave-facet CENTROID, and let the kernel's overlap
    //     clip reject non-overlapping candidates (status<0 ⇒ the adapter is inert that
    //     eval). Frictionless penalty (λ_N Uzawa = C2.2; friction = C3). No path state. ---
    if (cd != 0) {
        cd->mortarNormalGCBegin();   // C2.2: rebuild the live λ_N node-set this handle()
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": master/slave surface not defined; skipped\n";
                continue;
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) {
                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": need a MASTER_SEGMENTS master + SLAVE_SEGMENTS slave; skipped\n";
                continue;
            }
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) {
                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": nodesPerSeg must be 3 or 4 (slave " << npsS << ", master " << npsM
                       << "); skipped\n";
                continue;
            }
            // penalty: epsN if given, else the kn slot; auto-size per pair if requested.
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            if (!epsAuto && epsFixed <= 0.0) {
                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": needs a penalty (-epsN val|auto or kn > 0); skipped\n";
                continue;
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;

            // C2.1 broad phase = BRUTE FORCE (every master facet a candidate per slave
            // facet); the kernel's overlap clip rejects non-overlaps (status<0 ⇒ inert).
            // The NTS bucket-sort superset contract is proven for POINT slaves (slave
            // NODES) only — a CENTROID query for an EXTENDED slave facet silently drops
            // overlaps when the slave facet is coarser than the master bucket (gate MAJOR,
            // C2.1 review: ~44% area lost on a coarse-slave/fine-master pairing). A
            // slave-aware mortar broad phase (corner-span query / max-diag cell) is a
            // follow-up ("mortar P2.5"); correctness-first brute force here, like NTS P2b-1.
            for (int sf = 0; sf < nSegS; sf++) {
                Node *sNodes[4]; bool ok = true; double scen[3] = {0.0, 0.0, 0.0};
                double Sx[4][3];                                  // ref coords (ADR-57 E2 stand-down)
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() != 3) {
                        if (sn != 0)
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave node " << sTags(sf * npsS + k) << " ndf="
                                   << sn->getNumberDOF() << " != 3; facet skipped (mortar is 3D)\n";
                        ok = false; break;
                    }
                    sNodes[k] = sn;
                    const Vector &Xk = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = Xk(d); scen[d] += Xk(d); }
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) scen[d] /= npsS;     // slave-facet centroid

                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mNodes[4]; ok = true; double mcen[3] = {0.0, 0.0, 0.0};
                    double Mx[4][3];                              // ref coords (ADR-57 E2 stand-down)
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        mNodes[k] = mn;
                        const Vector &Xk = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = Xk(d); mcen[d] += Xk(d); }
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) mcen[d] /= npsM;
                    // ADR-57 E2 — FACE-MORTAR STANDS DOWN when the edge-edge lane owns this pair (the
                    // gate EE-1 fix): in the oblique band (50°–90°) the mortar clip is NON-degenerate,
                    // so without this the clip pressure and the injected edge-edge force would
                    // SUPERPOSE on the same region. The SAME routing decision the edge-edge loop uses
                    // (one source of truth) ⇒ exactly one lane owns the pair. mc.edgeEdge off ⇒ the
                    // shipped mortar behavior, unchanged.
                    if (mc.edgeEdge &&
                        ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx, ladrunoEdgeBand(npsM, Mx, mc.edgeBand)))
                        continue;
                    // orientation toward the slave's allowed half-space: explicit -outward,
                    // else (slave facet centroid − master facet centroid) at the ref config.
                    double orientDir[3];
                    if (mc.hasOutward)
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    else
                        for (int d = 0; d < 3; d++) orientDir[d] = scen[d] - mcen[d];

                    double epsUse = epsFixed;
                    if (epsAuto) {
                        epsUse = ladrunoResolveAutoKn(theDomain, mNodes, npsM, orientDir);
                        if (epsUse <= 0.0) {
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave facet " << sf << " master facet " << seg
                                   << ": auto penalty could not be sized; pair skipped\n";
                            continue;
                        }
                    }
                    // C3.1 friction: epsT auto ⇒ size from the normal penalty (epsUse); else the
                    // given value. mu/cohesion/tauMax ≤0 ⇒ the adapter short-circuits (frictionless).
                    double epsTuse = mc.epsTAuto ? epsUse : mc.epsT;
                    // gate MINOR-1: friction REQUESTED but no tangential penalty given ⇒ the return
                    // map would stick at ZERO force (silently inert). Default epsT to the normal
                    // penalty and warn ONCE rather than ship dead friction.
                    bool wantFric = (mc.mu > 0.0 || mc.cohesion > 0.0 || mc.tauMax > 0.0);
                    // C3.3: -consistanttan opts into the NON-SYMMETRIC consistent Coulomb friction
                    // tangent (Csl), which REQUIRES a non-symmetric solver (system FullGeneral /
                    // UmfPack / BandGeneral); a symmetric SOE silently drops the lower triangle and
                    // corrupts it. The default (symmetric) tangent is correct on any solver. Warn once.
                    if (wantFric && mc.consistentTan) {
                        static bool warnedCsl = false;
                        if (!warnedCsl) {
                            warnedCsl = true;
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": -consistanttan (non-symmetric Coulomb friction "
                                      "tangent) needs a non-symmetric solver (system FullGeneral/"
                                      "UmfPack/BandGeneral); symmetric solvers will silently corrupt it.\n";
                        }
                    }
                    if (wantFric && epsTuse <= 0.0) {
                        epsTuse = epsUse;
                        static bool warnedEpsT = false;
                        if (!warnedEpsT) {
                            warnedEpsT = true;
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": friction requested without -epsT; defaulting the "
                                      "tangential penalty to the normal penalty (-epsT auto). Set "
                                      "-epsT explicitly to control it.\n";
                        }
                    }
                    // C4 — a mesh-tie pair (mc.isTie): epsUse is epsTie; friction params are 0
                    // (refused with -tie upstream). The adapter runs the full-3-vec r→0 bond.
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sNodes, npsS, mNodes, npsM, epsUse,
                                             orientDir, mc.tag, sf, theDomain,
                                             mc.mu, epsTuse, mc.cohesion, mc.tauMax, mc.consistentTan,
                                             mc.isTie, mc.muc,    // D2.2 mortar viscous (0 ⇒ off)
                                             mc.softScale);       // B2 SOFT=2 segment-based penalty (0 ⇒ off)
                    if (fe == 0) return -5;
                    theModel->addFE_Element(fe);
                    // C2.2: this pair's slave nodes have a live λ_N slot this handle().
                    for (int k = 0; k < npsS; k++)
                        cd->mortarNormalGCMark(mc.tag, sTags(sf * npsS + k));
                }
            }
        }
        cd->mortarNormalGCEnd();   // prune λ_N slots no live mortar pair referenced
    }

    // --- ADR-57 E2: perpendicular edge-edge routing + adapter injection (the cos_t→0 fallback) ---
    // For each -mortar contact with -edgeedge, run the E1 routing detector over the SAME brute-force
    // (slave facet, master facet) enumeration the mortar lane uses (NOT the bucket sort — that is
    // NTS point-vs-segment only) and, on an EDGE_EDGE route, inject one EDGE_EDGE adapter per kept
    // boundary-edge pair. The mortar loop STANDS DOWN on a pair edge-edge owns (the same
    // ladrunoEdgeEdgeOwns decision), so the mortar clip and the edge-edge force never superpose —
    // even in the oblique band where the clip is non-degenerate (the gate EE-1 fix). Reference-config
    // geometry (§7 scoping; large sliding ⇒ re-emit). -edgeedge absent ⇒ zero adapters ⇒ byte-
    // identical; -edgeedge present but no pair routes ⇒ also zero adapters.
    if (cd != 0) {
        cd->edgeGCBegin();
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            if (!mc.edgeEdge) continue;                 // opt-in off ⇒ skip (byte-identical)
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) continue;           // already warned by the mortar loop
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) continue;
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) continue;
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;
            // de-dup: a facet edge is shared by two facets, so key each injected adapter by the
            // ORDERED (slave edge nodes, master edge nodes) 4-tuple and inject it once (the §7 step-3
            // set, mirroring the bucket-sort stamp idiom — a NEW handler set, not a reuse).
            std::set<std::vector<int> > injectedEdges;
            for (int sf = 0; sf < nSegS; sf++) {
                Node *sN[4]; bool ok = true; double Scen[3] = {0, 0, 0}, Sx[4][3];
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() != 3) { ok = false; break; }
                    sN[k] = sn;
                    const Vector &X = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = X(d); Scen[d] += X(d); }
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) Scen[d] /= npsS;
                double nS[3]; ladrunoFacetNormalNewell(npsS, Sx, nS);
                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mN[4]; ok = true; double Mcen[3] = {0, 0, 0}, Mx[4][3];
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        mN[k] = mn;
                        const Vector &X = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = X(d); Mcen[d] += X(d); }
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) Mcen[d] /= npsM;

                    // does edge-edge OWN this pair? (the single-source-of-truth routing — the SAME
                    // decision the mortar stand-down uses, so exactly one lane owns the pair). d_band
                    // from -edgeBand, else the facet edge length.
                    double dBand = ladrunoEdgeBand(npsM, Mx, mc.edgeBand);
                    if (!ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx, dBand)) continue;

                    // orientation toward the slave's allowed half-space (the §2 first-capture sign
                    // reference): explicit -outward, else the SLAVE FACET NORMAL (always non-degenerate,
                    // unlike a centroid difference that can vanish near-coincident centroids — the gate
                    // E2-2 fix; ADR §2 "orient once from {n̂_s, n̂_m}") flipped from master toward slave.
                    double orientDir[3];
                    if (mc.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    } else {
                        double cdv[3] = {Scen[0]-Mcen[0], Scen[1]-Mcen[1], Scen[2]-Mcen[2]};
                        double dc = nS[0]*cdv[0] + nS[1]*cdv[1] + nS[2]*cdv[2];
                        double sgn = (dc >= 0.0) ? 1.0 : -1.0;
                        for (int d = 0; d < 3; d++) orientDir[d] = sgn * nS[d];
                    }
                    // edge penalty: explicit -edgeKn, else auto per master facet (if -edgeKn auto or the
                    // mortar penalty is auto), else the resolved fixed mortar penalty (the ADR default).
                    double edgeKnUse = mc.edgeKn;
                    if (edgeKnUse <= 0.0) {
                        if (mc.edgeKnAuto || epsAuto)
                            edgeKnUse = ladrunoResolveAutoKn(theDomain, mN, npsM, orientDir);
                        else
                            edgeKnUse = epsFixed;
                    }
                    if (edgeKnUse <= 0.0) {
                        opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                               << " slave facet " << sf << " master facet " << seg
                               << ": -edgeedge penalty could not be sized; pair skipped\n";
                        continue;
                    }
                    // E3 friction: edgeKt (tangential penalty). Friction requested but no -edgeKt ⇒
                    // the return map would stick at ZERO force (silently inert); default it to the
                    // edge normal penalty and warn ONCE (mirrors the mortar -epsT auto-default).
                    double edgeKtUse = mc.edgeKt;
                    bool edgeFric = (mc.edgeMu > 0.0 || mc.edgeCohesion > 0.0 || mc.edgeTauMax > 0.0);
                    if (edgeFric && edgeKtUse <= 0.0) {
                        edgeKtUse = edgeKnUse;
                        static bool warnedEdgeKt = false;
                        if (!warnedEdgeKt) {
                            warnedEdgeKt = true;
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                                   << ": edge-edge friction requested without -edgeKt; defaulting the "
                                      "tangential penalty to the edge normal penalty. Set -edgeKt to control it.\n";
                        }
                    }

                    // enumerate the (≤4)×(≤4) boundary-edge pairs; inject the margin-interior, well-
                    // conditioned, in-band crossings (the E0 kernel decides interior + non-parallel).
                    for (int ie = 0; ie < npsS; ie++) {
                        Node *sa = sN[ie], *sb = sN[(ie + 1) % npsS];
                        for (int je = 0; je < npsM; je++) {
                            Node *ma = mN[je], *mb = mN[(je + 1) % npsM];
                            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                                Sx[ie], Sx[(ie + 1) % npsS], Mx[je], Mx[(je + 1) % npsM]);
                            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
                            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
                            double g = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                            if (g > dBand) continue;
                            int sat = sa->getTag(), sbt = sb->getTag();
                            int mat = ma->getTag(), mbt = mb->getTag();
                            std::vector<int> key(4);
                            key[0] = (sat < sbt) ? sat : sbt; key[1] = (sat < sbt) ? sbt : sat;
                            key[2] = (mat < mbt) ? mat : mbt; key[3] = (mat < mbt) ? mbt : mat;
                            if (injectedEdges.find(key) != injectedEdges.end()) continue;
                            injectedEdges.insert(key);
                            LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sa, sb, ma, mb,
                                                       edgeKnUse, orientDir, mc.tag, theDomain,
                                                       mc.edgeMu, edgeKtUse, mc.edgeCohesion,
                                                       mc.edgeTauMax, mc.edgeConsistentTan,
                                                       mc.edgeSoftScale,    // E5 explicit SOFT (0 ⇒ off)
                                                       mc.edgeAlm);         // E6 one-scalar ALM (off ⇒ penalty)
                            if (fe == 0) return -5;
                            theModel->addFE_Element(fe);
                            cd->edgeGCMark(mc.tag, sat, sbt, mat, mbt);  // live this handle()
                        }
                    }
                }
            }
        }
        cd->edgeGCEnd();   // prune edge-edge slots no live pair referenced
    }

    // P2a: rigid analytical-plane contacts -> ONE bound adapter per slave node
    // (connectivity = that node). ndm derived per-node from its coordinates.
    if (cd != 0) {
        for (int p = 0; p < cd->getNumRigidPlanes(); p++) {
            const LadrunoContactDomain::RigidPlane &rp = cd->getRigidPlane(p);
            LadrunoContactSurface *surf = cd->getSurface(rp.slaveSurfTag);
            if (surf == 0) continue;
            const ID &nodeTags = surf->getNodeTags();
            for (int k = 0; k < nodeTags.Size(); k++) {
                Node *sn = theDomain->getNode(nodeTags(k));
                if (sn == 0) {
                    opserr << "WARNING LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " not in domain; skipped\n";
                    continue;
                }
                int nd = sn->getCrds().Size();
                if (sn->getNumberDOF() < nd) {
                    // The adapter couples the slave's first nd translational DOFs;
                    // a node with fewer DOFs than coordinates (ndf < ndm) would let
                    // base setID() leave trailing myID slots at 0 (contact assembled
                    // onto equation 0) and rigidPlaneGap() read past getTrialDisp().
                    // Skip with a warning rather than silently corrupt the system.
                    opserr << "WARNING LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " has ndf=" << sn->getNumberDOF()
                           << " < ndm=" << nd << "; skipped\n";
                    continue;
                }
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn,
                                                            rp.muc,         // D2 viscous (0 ⇒ off)
                                                            rp.softScale,   // B1 SOFT=1 (0 ⇒ off)
                                                            theDomain);     // B1 engine mass cache
                if (fe == 0) return -5;
                theModel->addFE_Element(fe);
            }
        }
    }

    return count3;
}

void
LadrunoContactHandler::clearAll(void)
{
    // reset node DOF_Group pointers (as PlainHandler); the AnalysisModel owns and
    // deletes the DOF_Groups + FE_Elements (incl. our contact adapters).
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return;
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0)
        nodPtr->setDOF_GroupPtr(0);
}

int
LadrunoContactHandler::sendSelf(int, Channel &)
{
    // P1a single-process stub (no state to serialize); parallel/database in a
    // later phase rides LadrunoContactDomain::sendSelf, not the handler.
    return 0;
}

int
LadrunoContactHandler::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}
