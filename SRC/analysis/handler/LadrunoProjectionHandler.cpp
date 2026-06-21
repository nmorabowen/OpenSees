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
#include <math.h>     // Ladruno: fabs (libstdc++/gcc needs the explicit include; MSVC pulled it transitively)
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
#include <EQ_Constraint.h>
#include <EQ_ConstraintIter.h>
#include <Integrator.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <classTags.h>

// Ladruno (ADR-30 P6): a frozen small-rotation Ccr (rigidLink-beam / rigidDiaphragm lever arm)
// loses validity once the master rotation accumulates past ~0.1 rad. Warn (once per node) above.
#define LADRUNO_CCR_STALE_RAD 0.1

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

    // P4b/P4c: a fresh handle() rebuilds the prescribed-motion + derived-slave records.
    prescribedDOFs.clear();
    prescribedKey.clear();
    derivedSlaves.clear();
    derivedKey.clear();

    // Collect every SP (Plain-style). Homogeneous SPs (`fix`) are excluded from the
    // equation set below. Non-homogeneous SPs / imposedMotion are ALSO excluded from the
    // equation set (eqn = -1), but their prescribed DISPLACEMENT is imposed on the node each
    // step by applyLoad() (P4b); ImposedMotionSP supplies vel/accel itself. They are recorded
    // in the node loop below (where the Node* and the -1 mark are both in hand).
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0)
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));

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
                // P4b: a non-homogeneous SP (or imposedMotion) carries a prescribed motion.
                // Record it so applyLoad() can impose the displacement on the node each step.
                if (it->second->isHomogeneous() == false) {
                    PrescribedDOF p;
                    p.nodeTag = nodeID;        // resolved fresh in applyLoad() (no stored Node*)
                    p.dof     = dof;
                    p.sp      = it->second;
                    prescribedDOFs.push_back(p);
                    prescribedKey.insert(std::make_pair(nodeID, dof));
                }
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
    // P4c: a slave driven ONLY by prescribed master(s) is fully determined by the prescribed
    // motion -> exclude it from the equation set (eqn=-1) and impose it kinematically in
    // applyLoad(). Refuses a slave tied to BOTH a free and a prescribed master (mixed). This
    // runs after buildGroups (needs the groups) and re-sizes the equation count below.
    if (this->classifyDerivedSlaves(countDOF) < 0)
        return -8;
    theModel->setNumEqn(countDOF);
    // Zero free equations (every DOF fixed / prescribed / derived-from-prescribed) leaves the
    // explicit integrator nothing to solve and a 0-size SOE (a segfault downstream). Refuse with
    // a clean, named error (ADR §5.2 "zero free DOFs -> clean error, not segfault"). count3 are
    // subdomain-boundary DOFs (numbered later), so they still count as free here.
    if (countDOF == 0 && count3 == 0) {
        opserr << "LadrunoProjectionHandler::handle() - no free equations: every DOF is fixed, "
                  "prescribed (SP/imposedMotion), or driven by a prescribed master. There is "
                  "nothing for the explicit integrator to solve. Add a free DOF, or this model is "
                  "pure kinematic playback (not supported by the projection handler).\n";
        return -9;
    }
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
    rotMonitors.clear();   // P6: rebuilt below from the rotational lever-arm couplings

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
                // P6: a nonzero Ccr entry coupling a master ROTATIONAL DOF into a slave
                // TRANSLATIONAL DOF is a frozen lever arm (rigidLink-beam / rigidDiaphragm)
                // that goes stale under finite rotation -> monitor it.
                this->flagRotMonitor(rNode, rDofs(j), cNode, cDofs(i));
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

    // EQ_Constraint (P3): one constrained DOF tied to a coefficient VECTOR of retained
    // DOFs (u_c = sum_k C[k] u_{r_k}), possibly across several retained nodes — i.e. the
    // multi-master general-C row the projector already handles. Mirror the MP edge build.
    EQ_ConstraintIter &theEQs = theDomain->getEQs();
    EQ_Constraint *eqPtr;
    while ((eqPtr = theEQs()) != 0) {
        int cNode = eqPtr->getNodeConstrained();
        int cDof  = eqPtr->getConstrainedDOFs();
        const ID &rNodes = eqPtr->getNodeRetained();
        const ID &rDofs  = eqPtr->getRetainedDOFs();
        const Vector &C  = eqPtr->getConstraint();

        std::pair<int, int> cp(cNode, cDof);
        int cv;
        std::map<std::pair<int,int>,int>::iterator itc = vtxOf.find(cp);
        if (itc == vtxOf.end()) {
            cv = (int)vtxNode.size();
            vtxOf[cp] = cv; vtxNode.push_back(cNode); vtxDof.push_back(cDof);
            parent.push_back(cv); slaveMPcount.push_back(0);
            isMaster.push_back(0); isSlave.push_back(0);
        } else cv = itc->second;
        isSlave[cv] = 1;
        slaveMPcount[cv] += 1;

        for (int k = 0; k < rDofs.Size(); k++) {
            double c = C(k);
            if (c == 0.0) continue;
            std::pair<int, int> rp(rNodes(k), rDofs(k));
            int rv;
            std::map<std::pair<int,int>,int>::iterator itr = vtxOf.find(rp);
            if (itr == vtxOf.end()) {
                rv = (int)vtxNode.size();
                vtxOf[rp] = rv; vtxNode.push_back(rNodes(k)); vtxDof.push_back(rDofs(k));
                parent.push_back(rv); slaveMPcount.push_back(0);
                isMaster.push_back(0); isSlave.push_back(0);
            } else rv = itr->second;
            isMaster[rv] = 1;
            Edge e; e.slaveVtx = cv; e.masterVtx = rv; e.coeff = c; e.delta = 0.0;
            edges.push_back(e);
        }
    }

    int nv = (int)vtxNode.size();
    if (nv == 0)
        return 0;   // no MPs/EQs -> nothing to project (still a valid analysis)

    // Partition-locality guard (ADR §4): every DOF a constraint references must resolve
    // to a node in THIS domain. A node tag absent from the domain = a partition-crossing
    // MP/EQ under OpenSeesSP/MP (v1 is partition-interior only). Refuse with a named error
    // rather than silently mis-assemble — doneNumberingDOF() would otherwise treat the
    // missing master as SP-fixed and drop its column (a silent wrong answer).
    for (int v = 0; v < nv; v++) {
        if (theDomain->getNode(vtxNode[v]) == 0) {
            opserr << "LadrunoProjectionHandler - constraint references node " << vtxNode[v]
                   << " not present in this domain. v1 requires partition-LOCAL constraint "
                      "groups (no partition-crossing MP/EQ); use constraints Transformation "
                      "for distributed analysis.\n";
            return -1;
        }
    }

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

    // EQ re-walk: populate the slave rows (masters/coeffs/delta) for each EQ_Constraint.
    EQ_ConstraintIter &theEQs2 = theDomain->getEQs();
    while ((eqPtr = theEQs2()) != 0) {
        int cNode = eqPtr->getNodeConstrained();
        int cDof  = eqPtr->getConstrainedDOFs();
        const ID &rNodes = eqPtr->getNodeRetained();
        const ID &rDofs  = eqPtr->getRetainedDOFs();
        const Vector &C  = eqPtr->getConstraint();
        double Uc0 = eqPtr->getConstrainedDOFsInitialDisplacement();
        const Vector &Ur0 = eqPtr->getRetainedDOFsInitialDisplacement();

        int cv = vtxOf[std::make_pair(cNode, cDof)];
        int r = cv; while (par[r] != r) r = par[r];
        RawGroup &g = theGroups[rootToGroup[r]];
        SlaveRec s;
        s.vtx = cv;
        s.delta = Uc0;
        for (int k = 0; k < rDofs.Size(); k++) {
            double c = C(k);
            if (c == 0.0) continue;
            int rv = vtxOf[std::make_pair(rNodes(k), rDofs(k))];
            s.masterVtx.push_back(rv);
            s.coeff.push_back(c);
            if (k < Ur0.Size()) s.delta -= c * Ur0(k);
        }
        g.slaves.push_back(s);
    }
    return 0;
}

// ---- P4c: classify slaves driven by prescribed masters; exclude + record them ----------
// A slave tied ONLY to prescribed master(s) (no free retained master) is fully determined by
// the prescribed motion. Exclude it from the equation set (eqn=-1, decrement countDOF) and
// record it for kinematic imposition in applyLoad() (u/v/a set from the masters, exact — no
// projection drift). A slave tied to BOTH a free and a prescribed master is REFUSED (mixed:
// the displacement tie to the prescribed master cannot be held by the free-DOF projection).
// Masters are classified by their DOF mark at this point (free = -2/-3 not-yet-numbered; SP =
// -1, prescribed iff in prescribedKey, else homogeneous-fixed). A homogeneous-fixed master
// contributes 0 (its -Ur0 already folded into the slave's delta), so it does not make a slave
// "free". Runs after buildGroups(), before equation numbering.
int
LadrunoProjectionHandler::classifyDerivedSlaves(int &countDOF)
{
    Domain *theDomain = this->getDomainPtr();
    for (int gi = 0; gi < (int)theGroups.size(); gi++) {
        RawGroup &g = theGroups[gi];
        for (int s = 0; s < (int)g.slaves.size(); s++) {
            SlaveRec &sr = g.slaves[s];
            Node *snp = theDomain->getNode(vtxNode[sr.vtx]);
            int ec = (snp != 0) ? snp->getDOF_GroupPtr()->getID()(vtxDof[sr.vtx]) : -1;
            // At this point (post-buildGroups, pre-numbering) a normal MP/EQ slave DOF is marked
            // -2 (free, awaiting an equation); -1 = the slave is itself SP-fixed/prescribed (the
            // SP-on-slave / prescribed-slave overconstraint, refused in doneNumberingDOF); -3 =
            // subdomain boundary. Only a -2 slave is a candidate for prescribed-driven exclusion.
            if (ec != -2)
                continue;

            bool hasFree = false, hasPres = false;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                int mv = sr.masterVtx[k];
                Node *mnp = theDomain->getNode(vtxNode[mv]);
                int me = (mnp != 0) ? mnp->getDOF_GroupPtr()->getID()(vtxDof[mv]) : -1;
                bool mPres = prescribedKey.find(std::make_pair(vtxNode[mv], vtxDof[mv]))
                             != prescribedKey.end();
                if (mPres) hasPres = true;
                else if (me != -1) hasFree = true;   // -2/-3 = a free DOF that will be numbered
                // me==-1 && !mPres -> homogeneous fix: contributes 0, not "free"
            }

            if (hasFree && hasPres) {
                opserr << "LadrunoProjectionHandler - MIXED MASTER: constrained DOF (node "
                       << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") is tied to BOTH a "
                          "free DOF and a prescribed-motion DOF. Driving one slave from a mix of "
                          "free and prescribed masters is not supported (the prescribed "
                          "displacement tie cannot be held by the free-DOF acceleration "
                          "projection). Split the constraint, or use constraints Transformation.\n";
                return -1;
            }
            if (!hasPres)
                continue;          // free-only or all-homogeneous-fixed -> normal projection path

            // derived-prescribed slave: exclude from the equation set + record for applyLoad()
            snp->getDOF_GroupPtr()->setID(vtxDof[sr.vtx], -1);
            countDOF--;
            DerivedSlave d;
            d.nodeTag = vtxNode[sr.vtx];
            d.dof     = vtxDof[sr.vtx];
            d.delta   = sr.delta;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                int mv = sr.masterVtx[k];
                if (prescribedKey.find(std::make_pair(vtxNode[mv], vtxDof[mv])) == prescribedKey.end())
                    continue;      // skip homogeneous-fixed masters (their -Ur0 is in delta)
                d.masterNode.push_back(vtxNode[mv]);
                d.masterDof.push_back(vtxDof[mv]);
                d.coeff.push_back(sr.coeff[k]);
            }
            derivedSlaves.push_back(d);
            derivedKey.insert(std::make_pair(d.nodeTag, d.dof));

            if (verbose)
                opserr << "  LadrunoProjection: DOF (node " << d.nodeTag << " dof " << d.dof
                       << ") driven by " << (int)d.masterNode.size()
                       << " prescribed master(s) -> kinematically imposed.\n";
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
            // P4c: a prescribed-motion master (eqn<0, in prescribedKey) is DROPPED here, like a
            // homogeneous fix — it drives only derived-prescribed slaves, which were excluded by
            // classifyDerivedSlaves() and are imposed kinematically in applyLoad(). (A mix of free
            // and prescribed masters on one slave was already refused there.) A homogeneous-fixed
            // master (eqn<0, not prescribed) is likewise dropped (its slave -> a=0 via fixedEqns).
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
                // P4c: a derived-prescribed slave (driven purely by prescribed masters) was
                // excluded by classifyDerivedSlaves() and is imposed kinematically in
                // applyLoad() — skip it here (it owns no equation / projector row).
                if (derivedKey.find(std::make_pair(vtxNode[sr.vtx], vtxDof[sr.vtx]))
                    != derivedKey.end())
                    continue;
                if (prescribedKey.find(std::make_pair(vtxNode[sr.vtx], vtxDof[sr.vtx]))
                    != prescribedKey.end())
                    opserr << "LadrunoProjectionHandler - PRESCRIBED SLAVE: constrained DOF "
                              "(node " << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") "
                              "carries a prescribed motion (non-homogeneous SP / imposedMotion) "
                              "AND is tied by an MP/EQ. A DOF cannot be both prescribed and "
                              "constraint-driven (overconstraint). Remove the tie or the "
                              "prescribed motion.\n";
                else
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

        if (verbose)
            opserr << "  LadrunoProjection group " << gi << ": " << nRet
                   << " retained DOF(s), " << nCon << " constrained, "
                   << fe.Size() << " SP-fixed-master slave(s) zeroed.\n";
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
                      "an explicit projection-aware integrator (e.g. CentralDifferenceLadruno, "
                      "ExplicitBathe, ExplicitBatheLNVD, or their SMS variants). "
                      "Use constraints Transformation/Penalty for implicit analyses.\n";
            return -4;
        }
    } else {
        consumer->setConstraintProjector(theProjector->numGroups() > 0 ? theProjector : 0);
    }
    return 0;
}

// ===========================================================================
// applyLoad(): impose prescribed-motion displacement on the recorded SP DOFs (P4b)
// ===========================================================================
// Called from AnalysisModel::updateDomain(time,dt) AFTER theDomain->applyLoad(time) — so an
// ImposedMotionSP has already set the node's trial vel/accel from its GroundMotion, and we set
// the displacement here (the disp is "the responsibility of the constraint handler", per
// ImposedMotionSP). A plain non-homogeneous SP supplies only the displacement; its vel/accel
// stay 0 (the same limitation constraints Transformation carries — not a regression). The DOF
// is SP-excluded (eqn = -1), so the integrator's setResponse() scatter skips it and the imposed
// value survives (DOF_Group::setNodeDisp guards on loc>=0). Mirrors
// TransformationDOF_Group::enforceSPs(doMP=1).

// P6: register a frozen lever-arm coupling for staleness monitoring. Called from the transport
// loop with a master (node,dof) -> slave (node,dof) whose Ccr entry is nonzero. The lever arm
// that goes stale under finite rotation is the MASTER-ROTATION -> SLAVE-TRANSLATION cross term
// (rigidLink-beam / rigidDiaphragm). A rotation->rotation tie (equalDOF on rz) or a pure
// translation tie (rigidLink-bar) carries no lever arm and is ignored. Deduped by master node.
void
LadrunoProjectionHandler::flagRotMonitor(int masterNode, int masterDof,
                                         int slaveNode, int slaveDof)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0) return;
    Node *mp = theDomain->getNode(masterNode);
    Node *sp = theDomain->getNode(slaveNode);
    if (mp == 0 || sp == 0) return;
    int ndmM = mp->getCrds().Size();
    int ndfM = mp->getNumberDOF();
    int ndmS = sp->getCrds().Size();
    // OpenSees convention: 2D node (ndm==2) -> dof 2 is rz; 3D 6-dof node -> dofs 3,4,5 are rotations;
    // translational DOFs are 0..(ndm-1).
    bool masterRot  = (ndmM == 2 && masterDof >= 2) || (ndmM == 3 && ndfM >= 6 && masterDof >= 3);
    bool slaveTrans = (slaveDof >= 0 && slaveDof < ndmS);
    if (!masterRot || !slaveTrans) return;
    for (int i = 0; i < (int)rotMonitors.size(); i++)
        if (rotMonitors[i].nodeTag == masterNode) return;   // already monitored
    RotMonitor m;
    m.nodeTag = masterNode;
    m.warned = false;
    if (ndmM == 2) { m.rotDof.push_back(2); }
    else           { m.rotDof.push_back(3); m.rotDof.push_back(4); m.rotDof.push_back(5); }
    const Vector &u = mp->getTrialDisp();   // at handle() time trial == committed == initial config
    for (int k = 0; k < (int)m.rotDof.size(); k++) {
        int d = m.rotDof[k];
        m.theta0.push_back((d < u.Size()) ? u(d) : 0.0);
    }
    rotMonitors.push_back(m);
}

int
LadrunoProjectionHandler::applyLoad(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return 0;
    // Resolve the node fresh each step (store the tag, not a Node* — matches the handler's
    // vtxNode convention and is robust to a node removed without a re-handle).
    for (int i = 0; i < (int)prescribedDOFs.size(); i++) {
        PrescribedDOF &p = prescribedDOFs[i];
        if (p.sp == 0)
            continue;
        Node *np = theDomain->getNode(p.nodeTag);
        if (np == 0)
            continue;
        np->setTrialDisp(p.sp->getValue() + p.sp->getInitialValue(), p.dof);
    }

    // P4c: impose the DERIVED-prescribed slaves (driven purely by prescribed masters). Runs
    // AFTER the masters' own disp is set above, so u_master is current; v/a of an imposedMotion
    // master were set by ImposedMotionSP in Domain::applyLoad (a plain SP master -> v=a=0). The
    // slave is eqn-excluded, so setResponse() skips it and these imposed values survive.
    //   u_c = sum_k C_k u_{m_k} + delta ;  v_c = sum_k C_k v_{m_k} ;  a_c = sum_k C_k a_{m_k}
    for (int i = 0; i < (int)derivedSlaves.size(); i++) {
        DerivedSlave &d = derivedSlaves[i];
        Node *snp = theDomain->getNode(d.nodeTag);
        if (snp == 0)
            continue;
        double uc = d.delta, vc = 0.0, ac = 0.0;
        for (int k = 0; k < (int)d.masterNode.size(); k++) {
            Node *mnp = theDomain->getNode(d.masterNode[k]);
            if (mnp == 0)
                continue;
            int md = d.masterDof[k];
            const Vector &mu = mnp->getTrialDisp();
            const Vector &mv = mnp->getTrialVel();
            const Vector &ma = mnp->getTrialAccel();
            double c = d.coeff[k];
            if (md >= 0 && md < mu.Size()) uc += c * mu(md);
            if (md >= 0 && md < mv.Size()) vc += c * mv(md);
            if (md >= 0 && md < ma.Size()) ac += c * ma(md);
        }
        // one bounds check guards all three setters consistently (d.dof is the recorded
        // constrained dof, so this holds in practice; defensive against a malformed record).
        int snd = snp->getNumberDOF();
        if (d.dof < 0 || d.dof >= snd)
            continue;
        snp->setTrialDisp(uc, d.dof);                       // single-dof disp setter
        // vel / accel have only whole-Vector setters: read-modify-write the slave's component
        Vector vv = snp->getTrialVel();  vv(d.dof) = vc;  snp->setTrialVel(vv);
        Vector av = snp->getTrialAccel(); av(d.dof) = ac; snp->setTrialAccel(av);
    }

    // P6: frozen-Ccr staleness warning. Watch each flagged master's rotation drift since the
    // constraints were built; warn ONCE per node once it exceeds the small-rotation validity of
    // the baked-in lever arm. (Behavior is unchanged — this is a diagnostic only.)
    for (int i = 0; i < (int)rotMonitors.size(); i++) {
        RotMonitor &m = rotMonitors[i];
        if (m.warned)
            continue;
        Node *np = theDomain->getNode(m.nodeTag);
        if (np == 0)
            continue;
        const Vector &u = np->getTrialDisp();
        double drift = 0.0;
        for (int k = 0; k < (int)m.rotDof.size(); k++) {
            int dd = m.rotDof[k];
            double th = (dd < u.Size()) ? u(dd) : 0.0;
            double a = fabs(th - m.theta0[k]);
            if (a > drift) drift = a;
        }
        if (drift > LADRUNO_CCR_STALE_RAD) {
            opserr << "LadrunoProjectionHandler - NOTE master node " << m.nodeTag
                   << " has rotated " << drift << " rad since the constraints were built; the "
                      "rigidLink -beam / rigidDiaphragm lever arm (frozen small-rotation Ccr) is "
                      "now stale and the tie will drift. For finite rotation use the RBE2 element "
                      "(LadrunoKinematicCoupling) or constraints Transformation.\n";
            m.warned = true;
        }
    }
    return 0;
}

double
LadrunoProjectionHandler::getTieForce(int nodeTag, int dof) const
{
    if (theProjector == 0)
        return 0.0;
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return 0.0;
    Node *np = theDomain->getNode(nodeTag);
    if (np == 0)
        return 0.0;
    DOF_Group *dg = np->getDOF_GroupPtr();
    if (dg == 0)
        return 0.0;
    const ID &id = dg->getID();
    if (dof < 0 || dof >= id.Size())
        return 0.0;
    int eqn = id(dof);
    if (eqn < 0)
        return 0.0;          // SP-fixed / unnumbered DOF carries no projected tie force
    return theProjector->tieForceAtEqn(eqn);
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
    prescribedDOFs.clear();   // P4b: drop the prescribed-motion records (rebuilt at next handle())
    prescribedKey.clear();
    derivedSlaves.clear();    // P4c: drop the derived-prescribed-slave records
    derivedKey.clear();
    rotMonitors.clear();      // P6: drop the frozen-Ccr staleness monitors
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
