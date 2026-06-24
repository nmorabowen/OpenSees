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

    // collect SPs (domain + load-pattern), keyed by node
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0)
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));

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
    while ((elePtr = theEle()) != 0) {
        if (elePtr->isSubdomain() == true) {
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

            // P2.5 BROAD PHASE: bucket-sort the master segments (reference coords)
            // once, then per slave query only the 27-neighbour candidate segments
            // instead of all nSeg (brute force). The kept set is a SUPERSET of every
            // near pair, so the contact result is identical to brute force; a huge
            // -cell (=> 1 bucket) reproduces brute force exactly (the equivalence
            // gate). Missing master nodes (malformed model) are filled with the
            // segment's first valid coord — superset-preserving, and the narrow loop
            // below re-fetches + skips them anyway.
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
            LadrunoContactBucketSort::Grid grid(nSeg, nps, segCoords.data(), ct.cellFrac, 1.0);
            std::vector<int> cand(nSeg);

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
                                             ct.consistentNormal);           // B3 ∂n/∂u geom tangent
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
                    for (int d = 0; d < 3; d++) scen[d] += Xk(d);
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) scen[d] /= npsS;     // slave-facet centroid

                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mNodes[4]; ok = true; double mcen[3] = {0.0, 0.0, 0.0};
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        mNodes[k] = mn;
                        const Vector &Xk = mn->getCrds();
                        for (int d = 0; d < 3; d++) mcen[d] += Xk(d);
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) mcen[d] /= npsM;
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
                                             mc.isTie, mc.muc);   // D2.2 mortar viscous (0 ⇒ off)
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
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn, rp.muc);   // D2 viscous (0 ⇒ off)
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
