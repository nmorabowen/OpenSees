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

// ADR-30 P1 — LadrunoProjectionHandler implementation. See the header.

#include <LadrunoProjectionHandler.h>
#include <LadrunoConstraintProjector.h>
#include <LadrunoProjectionConsumer.h>

#include <stdlib.h>
#include <string.h>
#include <map>
#include <utility>

#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Integrator.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <classTags.h>

// ---- command parser -------------------------------------------------------
void *OPS_LadrunoProjectionHandler(void)
{
    // constraints LadrunoProjection <-verbose> <-projectICs> <-icTol $tol>
    bool verbose = false;
    bool projectICs = false;
    double icTol = 1.0e-8;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0)
            verbose = true;
        else if (strcmp(arg, "-projectICs") == 0)
            projectICs = true;
        else if (strcmp(arg, "-icTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int n = 1;
                OPS_GetDoubleInput(&n, &icTol);
            }
        } else {
            opserr << "WARNING LadrunoProjection - unknown option " << arg << " (ignored)\n";
        }
    }
    return new LadrunoProjectionHandler(verbose, projectICs, icTol);
}

LadrunoProjectionHandler::LadrunoProjectionHandler()
  : ConstraintHandler(HANDLER_TAG_LadrunoProjectionHandler),
    vtxNode(), vtxDof(), theGroups(), theProjector(0),
    verbose(false), projectICs(false), icTol(1.0e-8)
{
}

LadrunoProjectionHandler::LadrunoProjectionHandler(bool v, bool p, double t)
  : ConstraintHandler(HANDLER_TAG_LadrunoProjectionHandler),
    vtxNode(), vtxDof(), theGroups(), theProjector(0),
    verbose(v), projectICs(p), icTol(t)
{
}

LadrunoProjectionHandler::~LadrunoProjectionHandler()
{
    if (theProjector != 0) delete theProjector;
}

// ===========================================================================
// handle(): Plain-style assembly WITHOUT MP elimination + group construction
// ===========================================================================
int
LadrunoProjectionHandler::handle(const ID *nodesLast)
{
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();
    if (theDomain == 0 || theModel == 0 || theIntegrator == 0) {
        opserr << "WARNING LadrunoProjectionHandler::handle() - setLinks() not called\n";
        return -1;
    }

    // collect homogeneous SPs (Plain-style); warn on non-homogeneous (not v1)
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0) {
        if (theSP->isHomogeneous() == false)
            opserr << "WARNING LadrunoProjectionHandler::handle() - non-homogeneous SP "
                      "at node " << theSP->getNodeTag()
                   << " not enforced in v1 (homogeneous assumed)\n";
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));
    }

    // DOF_Groups: free=-2, homogeneous-SP=-1. KEY: MP slave DOFs are NOT set to -4
    // (they keep their own equation + diagonal mass; the projector enforces the tie).
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    int numDOF = 0, countDOF = 0, count3 = 0;
    while ((nodPtr = theNod()) != 0) {
        DOF_Group *dofPtr = new DOF_Group(numDOF++, nodPtr);
        if (dofPtr == 0) {
            opserr << "WARNING LadrunoProjectionHandler::handle() - out of memory (DOF_Group)\n";
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
            if (cid(dof) == -2) {
                dofPtr->setID(dof, -1);
                countDOF--;
            } else {
                opserr << "WARNING LadrunoProjectionHandler::handle() - multiple SPs at DOF "
                       << dof << " node " << nodeID << endln;
            }
        }
        nodPtr->setDOF_GroupPtr(dofPtr);
        theModel->addDOF_Group(dofPtr);
    }
    theModel->setNumEqn(countDOF);

    // nodesLast (subdomain boundary): mark -3 (as PlainHandler)
    if (nodesLast != 0)
        for (int i = 0; i < nodesLast->Size(); i++) {
            Node *np = theDomain->getNode((*nodesLast)(i));
            if (np != 0) {
                DOF_Group *dp = np->getDOF_GroupPtr();
                const ID &id = dp->getID();
                for (int j = 0; j < id.Size(); j++)
                    if (id(j) == -2) { dp->setID(j, -3); count3++; }
            }
        }

    // FE_Elements for every element (as PlainHandler)
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

    // build the constraint groups + run the topological diagnostics
    if (this->buildGroups() < 0)
        return -6;
    if (this->consistentMassGuard() < 0)
        return -7;

    if (verbose)
        opserr << "LadrunoProjectionHandler: " << (int)theGroups.size()
               << " constraint group(s), " << countDOF << " equations.\n";

    return count3;
}

// ---- union-find over (node,dof) vertices + classify masters/slaves --------
int
LadrunoProjectionHandler::buildGroups(void)
{
    Domain *theDomain = this->getDomainPtr();
    theGroups.clear();
    vtxNode.clear();
    vtxDof.clear();

    std::map<std::pair<int, int>, int> vtxOf;   // (node,dof) -> vertex id
    std::vector<int> parent;                     // union-find

    // local lambdas replaced by helper indices (C++ in this tree is pre-lambda-safe)
    // find with path compression
    // (implemented inline below)

    // gather edges; one edge per (constrained dof i, retained dof j) with C(i,j)!=0
    struct Edge { int slaveVtx, masterVtx; double coeff; double delta; };
    std::vector<Edge> edges;
    std::vector<int> slaveMPcount;   // per vertex, #MPs it is a slave in
    std::vector<char> isMaster, isSlave;

    MP_ConstraintIter &theMPs = theDomain->getMPs();
    MP_Constraint *mpPtr;
    while ((mpPtr = theMPs()) != 0) {
        int cNode = mpPtr->getNodeConstrained();
        int rNode = mpPtr->getNodeRetained();
        const ID &cDofs = mpPtr->getConstrainedDOFs();
        const ID &rDofs = mpPtr->getRetainedDOFs();
        const Matrix &Ccr = mpPtr->getConstraint();
        const Vector &Uc0 = mpPtr->getConstrainedDOFsInitialDisplacement();
        const Vector &Ur0 = mpPtr->getRetainedDOFsInitialDisplacement();

        for (int i = 0; i < cDofs.Size(); i++) {
            std::pair<int, int> cp(cNode, cDofs(i));
            int cv;
            std::map<std::pair<int,int>,int>::iterator itc = vtxOf.find(cp);
            if (itc == vtxOf.end()) {
                cv = (int)vtxNode.size();
                vtxOf[cp] = cv; vtxNode.push_back(cNode); vtxDof.push_back(cDofs(i));
                parent.push_back(cv); slaveMPcount.push_back(0);
                isMaster.push_back(0); isSlave.push_back(0);
            } else cv = itc->second;
            isSlave[cv] = 1;
            slaveMPcount[cv] += 1;

            double delta = (i < Uc0.Size()) ? Uc0(i) : 0.0;
            for (int j = 0; j < rDofs.Size(); j++) {
                double c = Ccr(i, j);
                if (c == 0.0) continue;
                std::pair<int, int> rp(rNode, rDofs(j));
                int rv;
                std::map<std::pair<int,int>,int>::iterator itr = vtxOf.find(rp);
                if (itr == vtxOf.end()) {
                    rv = (int)vtxNode.size();
                    vtxOf[rp] = rv; vtxNode.push_back(rNode); vtxDof.push_back(rDofs(j));
                    parent.push_back(rv); slaveMPcount.push_back(0);
                    isMaster.push_back(0); isSlave.push_back(0);
                } else rv = itr->second;
                isMaster[rv] = 1;
                if (j < Ur0.Size()) delta -= c * Ur0(j);
                Edge e; e.slaveVtx = cv; e.masterVtx = rv; e.coeff = c; e.delta = 0.0;
                edges.push_back(e);
            }
            // stash this slave-row's delta on the FIRST edge of the row (the row is
            // rebuilt per slave below; we re-derive delta there). Carry via a parallel map.
            // (delta is finalized in the per-slave assembly using stored row data.)
        }
    }

    int nv = (int)vtxNode.size();
    if (nv == 0)
        return 0;   // no MPs -> nothing to project (still a valid analysis)

    // union-find connect
    for (int i = 0; i < nv; i++) parent[i] = i;
    // re-add unions from edges
    // find
    std::vector<int> &par = parent;
    // path-compressed find via loop
    // union each edge
    for (int e = 0; e < (int)edges.size(); e++) {
        int a = edges[e].slaveVtx, b = edges[e].masterVtx;
        while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
        while (par[b] != b) { par[b] = par[par[b]]; b = par[b]; }
        if (a != b) par[a] = b;
    }
    // diagnostics: chain (vtx both master & slave) and double-slave
    for (int v = 0; v < nv; v++) {
        if (isMaster[v] && isSlave[v]) {
            opserr << "LadrunoProjectionHandler - CHAIN: DOF (node " << vtxNode[v]
                   << ", dof " << vtxDof[v] << ") is both a retained master and a "
                      "constrained slave. v1 refuses MP chains; compose them or use "
                      "constraints Transformation.\n";
            return -1;
        }
        if (slaveMPcount[v] > 1) {
            opserr << "LadrunoProjectionHandler - DOUBLE CONSTRAINT: DOF (node "
                   << vtxNode[v] << ", dof " << vtxDof[v] << ") is constrained by more "
                      "than one MP. Remove the redundant tie.\n";
            return -1;
        }
    }

    // group by root
    std::map<int, int> rootToGroup;
    for (int v = 0; v < nv; v++) {
        int r = v;
        while (par[r] != r) r = par[r];
        if (rootToGroup.find(r) == rootToGroup.end()) {
            rootToGroup[r] = (int)theGroups.size();
            theGroups.push_back(RawGroup());
        }
    }
    // assign masters
    for (int v = 0; v < nv; v++) {
        if (!isMaster[v]) continue;
        int r = v; while (par[r] != r) r = par[r];
        theGroups[rootToGroup[r]].masterVtx.push_back(v);
    }
    // assign slaves with their masters/coeffs/delta (re-walk MPs for the row data)
    MP_ConstraintIter &theMPs2 = theDomain->getMPs();
    while ((mpPtr = theMPs2()) != 0) {
        int cNode = mpPtr->getNodeConstrained();
        int rNode = mpPtr->getNodeRetained();
        const ID &cDofs = mpPtr->getConstrainedDOFs();
        const ID &rDofs = mpPtr->getRetainedDOFs();
        const Matrix &Ccr = mpPtr->getConstraint();
        const Vector &Uc0 = mpPtr->getConstrainedDOFsInitialDisplacement();
        const Vector &Ur0 = mpPtr->getRetainedDOFsInitialDisplacement();
        for (int i = 0; i < cDofs.Size(); i++) {
            int cv = vtxOf[std::make_pair(cNode, cDofs(i))];
            int r = cv; while (par[r] != r) r = par[r];
            RawGroup &g = theGroups[rootToGroup[r]];
            SlaveRec s;
            s.vtx = cv;
            s.delta = (i < Uc0.Size()) ? Uc0(i) : 0.0;
            for (int j = 0; j < rDofs.Size(); j++) {
                double c = Ccr(i, j);
                if (c == 0.0) continue;
                int rv = vtxOf[std::make_pair(rNode, rDofs(j))];
                s.masterVtx.push_back(rv);
                s.coeff.push_back(c);
                if (j < Ur0.Size()) s.delta -= c * Ur0(j);
            }
            g.slaves.push_back(s);
        }
    }
    return 0;
}

// ---- refuse consistent (off-diagonal) element mass on any tied DOF (R5) ----
int
LadrunoProjectionHandler::consistentMassGuard(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theGroups.empty()) return 0;

    // tied (node,dof) set
    std::map<std::pair<int,int>, char> tied;
    for (int g = 0; g < (int)theGroups.size(); g++) {
        for (int m = 0; m < (int)theGroups[g].masterVtx.size(); m++) {
            int v = theGroups[g].masterVtx[m];
            tied[std::make_pair(vtxNode[v], vtxDof[v])] = 1;
        }
        for (int s = 0; s < (int)theGroups[g].slaves.size(); s++) {
            int v = theGroups[g].slaves[s].vtx;
            tied[std::make_pair(vtxNode[v], vtxDof[v])] = 1;
        }
    }
    std::map<int, char> tiedNode;
    for (std::map<std::pair<int,int>,char>::iterator it = tied.begin(); it != tied.end(); ++it)
        tiedNode[it->first.first] = 1;

    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    while ((elePtr = theEle()) != 0) {
        const ID &enodes = elePtr->getExternalNodes();
        bool incident = false;
        for (int k = 0; k < enodes.Size(); k++)
            if (tiedNode.find(enodes(k)) != tiedNode.end()) { incident = true; break; }
        if (!incident) continue;

        const Matrix &M = elePtr->getMass();
        int n = M.noRows();
        if (n == 0 || M.noCols() != n) continue;

        // local index -> (node,dof): nodes concatenated, each contributing its numDOF
        std::vector<int> idxNode(n), idxDof(n);
        int pos = 0;
        for (int k = 0; k < enodes.Size() && pos < n; k++) {
            Node *np = theDomain->getNode(enodes(k));
            int nd = (np != 0) ? np->getNumberDOF() : 0;
            for (int d = 0; d < nd && pos < n; d++) { idxNode[pos] = enodes(k); idxDof[pos] = d; pos++; }
        }
        if (pos != n) continue;   // could not map (unexpected ndf); skip conservatively

        double maxDiag = 0.0;
        for (int p = 0; p < n; p++) if (fabs(M(p,p)) > maxDiag) maxDiag = fabs(M(p,p));
        if (maxDiag == 0.0) continue;
        double tol = 1.0e-8 * maxDiag;

        for (int p = 0; p < n; p++)
            for (int q = 0; q < n; q++) {
                if (p == q) continue;
                if (fabs(M(p,q)) <= tol) continue;
                bool pTied = tied.find(std::make_pair(idxNode[p], idxDof[p])) != tied.end();
                bool qTied = tied.find(std::make_pair(idxNode[q], idxDof[q])) != tied.end();
                if (pTied || qTied) {
                    opserr << "LadrunoProjectionHandler - element " << elePtr->getTag()
                           << " assembles CONSISTENT (off-diagonal) mass coupling a tied "
                              "DOF (node " << idxNode[p] << " dof " << idxDof[p] << " <-> node "
                           << idxNode[q] << " dof " << idxDof[q] << "). The projection handler "
                              "requires lumped mass on tied DOFs; recreate the element with "
                              "-cMass 0, or use constraints Transformation with a banded SOE.\n";
                    return -1;
                }
            }
    }
    return 0;
}

// ===========================================================================
// doneNumberingDOF(): freeze equation indices, build the projector, push it
// ===========================================================================
int
LadrunoProjectionHandler::doneNumberingDOF(void)
{
    // base behavior: propagate the new equation numbers to every FE_Element
    // (setID). PlainHandler relies on this via the un-overridden base; we override
    // doneNumberingDOF, so we must invoke it explicitly or elements never assemble.
    this->ConstraintHandler::doneNumberingDOF();

    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return -1;

    if (theProjector != 0) { delete theProjector; theProjector = 0; }
    theProjector = new LadrunoConstraintProjector();
    theProjector->setICPolicy(projectICs, icTol);

    for (int gi = 0; gi < (int)theGroups.size(); gi++) {
        RawGroup &g = theGroups[gi];

        // surviving master columns (drop SP-fixed masters, eqn < 0)
        std::map<int, int> colOf;     // master vtx -> column index
        ID retEqnTmp(0);
        std::vector<int> survMasterVtx;
        for (int m = 0; m < (int)g.masterVtx.size(); m++) {
            int v = g.masterVtx[m];
            Node *np = theDomain->getNode(vtxNode[v]);
            int e = (np != 0) ? np->getDOF_GroupPtr()->getID()(vtxDof[v]) : -1;
            if (e >= 0) {
                colOf[v] = (int)survMasterVtx.size();
                survMasterVtx.push_back(v);
            }
        }
        int nRet = (int)survMasterVtx.size();

        // classify slave rows: surviving (>=1 master col left) vs orphaned (all fixed)
        std::vector<int> rowSlaveEqn;          // eqn per kept slave row
        std::vector<std::vector<int> > rowCols;
        std::vector<std::vector<double> > rowCoef;
        std::vector<double> rowDelta;
        std::vector<int> fixedEqns;
        for (int s = 0; s < (int)g.slaves.size(); s++) {
            SlaveRec &sr = g.slaves[s];
            Node *np = theDomain->getNode(vtxNode[sr.vtx]);
            int ec = (np != 0) ? np->getDOF_GroupPtr()->getID()(vtxDof[sr.vtx]) : -1;
            if (ec < 0) {
                opserr << "LadrunoProjectionHandler - SP-on-SLAVE: constrained DOF (node "
                       << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") is also SP-fixed. "
                          "Remove the fix or the tie.\n";
                return -2;
            }
            std::vector<int> cols; std::vector<double> coef;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                std::map<int,int>::iterator it = colOf.find(sr.masterVtx[k]);
                if (it != colOf.end()) { cols.push_back(it->second); coef.push_back(sr.coeff[k]); }
            }
            if (cols.empty()) {
                fixedEqns.push_back(ec);     // every master SP-fixed -> a_c = 0
            } else {
                rowSlaveEqn.push_back(ec);
                rowCols.push_back(cols);
                rowCoef.push_back(coef);
                rowDelta.push_back(sr.delta);
            }
        }

        int nCon = (int)rowSlaveEqn.size();
        // a group with no free columns AND no kept rows is fully fixed -> only fixedEqns
        if (nRet == 0 && nCon == 0) {
            if (!fixedEqns.empty()) {
                ID fe((int)fixedEqns.size());
                for (int i = 0; i < (int)fixedEqns.size(); i++) fe(i) = fixedEqns[i];
                Matrix L0(0, 0); Vector d0(0); ID a0(0);
                theProjector->addGroup(a0, 0, L0, d0, fe);
            }
            continue;
        }
        if (nRet == 0) {
            // kept slave rows but no surviving master: impossible unless masters fixed;
            // those rows were routed to fixedEqns above, so this is a logic guard.
            opserr << "LadrunoProjectionHandler - group " << gi
                   << " has constrained rows but no free retained DOF.\n";
            return -3;
        }

        // assemble allEqn (retained-first), L, delta, fixedEqn
        int nrows = nRet + nCon;
        ID allEqn(nrows);
        Matrix L(nrows, nRet);
        L.Zero();
        Vector delta(nCon > 0 ? nCon : 1);
        for (int p = 0; p < nRet; p++) {
            int v = survMasterVtx[p];
            allEqn(p) = theDomain->getNode(vtxNode[v])->getDOF_GroupPtr()->getID()(vtxDof[v]);
            L(p, p) = 1.0;
        }
        for (int j = 0; j < nCon; j++) {
            allEqn(nRet + j) = rowSlaveEqn[j];
            for (int k = 0; k < (int)rowCols[j].size(); k++)
                L(nRet + j, rowCols[j][k]) = rowCoef[j][k];
            delta(j) = rowDelta[j];
        }
        ID fe(fixedEqns.empty() ? 0 : (int)fixedEqns.size());
        for (int i = 0; i < (int)fixedEqns.size(); i++) fe(i) = fixedEqns[i];

        theProjector->addGroup(allEqn, nRet, L, delta, fe);
    }

    // Push the projector to the integrator via the abstract interface. ALWAYS re-bind
    // the consumer after a rebuild — including the zero-group case — because we just
    // deleted the previous projector at the top of this function; skipping the push on
    // the zero-group path (e.g. all MPs removed mid-run) would leave the integrator's
    // non-owning pointer dangling -> use-after-free in its domainChanged()/newStep()/
    // update() (Gate-B blocker). A zero-group projector is a valid no-op; passing 0
    // makes the integrator skip projection cleanly (every deref is `if (theProjector)`).
    Integrator *integ = this->getIntegratorPtr();
    LadrunoProjectionConsumer *consumer = dynamic_cast<LadrunoProjectionConsumer *>(integ);
    if (consumer == 0) {
        // Only fatal if there are constraints to enforce; with zero groups a
        // non-projection-aware integrator (Transformation/Penalty path) is legitimate.
        if (theProjector->numGroups() > 0) {
            opserr << "LadrunoProjectionHandler::doneNumberingDOF() - this handler requires "
                      "an explicit projection-aware integrator (e.g. CentralDifferenceLadruno). "
                      "Use constraints Transformation/Penalty for implicit analyses.\n";
            return -4;
        }
    } else {
        consumer->setConstraintProjector(theProjector->numGroups() > 0 ? theProjector : 0);
    }
    return 0;
}

void
LadrunoProjectionHandler::clearAll(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain != 0) {
        NodeIter &theNod = theDomain->getNodes();
        Node *nodPtr;
        while ((nodPtr = theNod()) != 0)
            nodPtr->setDOF_GroupPtr(0);
    }
    theGroups.clear();
    vtxNode.clear();
    vtxDof.clear();
    // doneNumberingDOF() now re-binds the integrator UNCONDITIONALLY (incl. a 0 push on
    // the zero-group path), so the integrator never retains a freed projector after a
    // rebuild. We keep ownership and free in the dtor / next doneNumberingDOF.
}

int
LadrunoProjectionHandler::sendSelf(int cTag, Channel &theChannel)
{
    static ID data(2);
    data(0) = projectICs ? 1 : 0;
    data(1) = verbose ? 1 : 0;
    theChannel.sendID(0, cTag, data);
    static Vector ddata(1);
    ddata(0) = icTol;
    theChannel.sendVector(0, cTag, ddata);
    return 0;
}

int
LadrunoProjectionHandler::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static ID data(2);
    theChannel.recvID(0, cTag, data);
    projectICs = (data(0) == 1);
    verbose = (data(1) == 1);
    static Vector ddata(1);
    theChannel.recvVector(0, cTag, ddata);
    icTol = ddata(0);
    return 0;
}
