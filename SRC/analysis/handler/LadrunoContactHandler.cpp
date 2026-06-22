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
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSeg = mTags.Size() / nps;
            for (int si = 0; si < sTags.Size(); si++) {
                Node *sn = theDomain->getNode(sTags(si));
                if (sn == 0 || sn->getNumberDOF() != 3) {
                    if (sn != 0)
                        opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                               << " slave node " << sTags(si) << " ndf=" << sn->getNumberDOF()
                               << " != 3; skipped (P2b is 3D translational)\n";
                    continue;
                }
                for (int seg = 0; seg < nSeg; seg++) {
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
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sn, segNodes, nps, ct.kn, orientDir);
                    if (fe == 0) return -5;
                    theModel->addFE_Element(fe);
                }
            }
        }
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
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn);
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
