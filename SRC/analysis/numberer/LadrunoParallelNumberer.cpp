/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Ladruno — LadrunoParallelNumberer (ADR-74). N2/T0: the O(V) engine.
//
// The algorithm is a faithful transcription of ParallelNumberer::numberDOF(int)
// — same gather, same vertex-creation order, same getFreeTag sequence, same
// append order, same RCM input, same scatter — with every O(V²) linear scan
// replaced by an O(1) lookup:
//   pass 1 dedup:  vertexRefs.getLocation(ref)        -> RefIndex (dense/hash)
//   pass 2 remap:  theSubdomainMap.getLocation(tag)   -> per-subgraph hash map
//   plain branch:  theOrderedRefs->getLocation(tag)   -> seen[] flags
//   -4 fixup:      full MP scan per constrained group -> constrainedNode index
// Bit-identity for the RCM verb is a THEOREM of "lookups only" (identical
// merged graph => identical RCM => identical numbers) and is ENFORCED by gate
// G1 (tests/test_adr74_numberer_1.py). The plain branch intentionally fixes
// the upstream tag-0 quirk (vertex 0 falsely "already ordered" against a
// zero-filled ID, numbered last) => G1b gates it by bijection + same-physics,
// NOT bit-identity. Measured baseline this replaces: dc.numberDOF ~ N^2.01
// (39.9 / 578.9 / 2477.5 s over 0.25/1/2 M hex, np8; 16.35 h at 19.18 M np240).
//
// Theory + gates: Ladruno_implementation/74_ladruno_parallel_numberer_adr.md

#include <LadrunoParallelNumberer.h>
#include <classTags.h>
#include <OPS_Globals.h>

#include <AnalysisModel.h>
#include <Domain.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Node.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <GraphNumberer.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <profiler/ProfilerMacros.h>

#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace {

// O(1) ref -> merge-position lookup (replaces vertexRefs.getLocation).
// Dense direct-address primary — gmsh/apeGmsh node tags are contiguous 1..N —
// with an unordered_map fallback for sparse/hand-authored tag spaces (the
// kappa guard; G1c-gated). All sizing math in 64-bit (ADR-74 risk item).
class RefIndex {
 public:
  explicit RefIndex(int expectedVertices)
      : expected_(expectedVertices > 0 ? expectedVertices : 1) {}

  // ref MUST be >= 0 — caller guards (Lagrange node-less DOF groups return -1
  // and stock code silently FUSES them; we error loudly instead).
  int lookup(int ref) const {
    if (useHash_) {
      std::unordered_map<int, int>::const_iterator it = hash_.find(ref);
      return (it == hash_.end()) ? -1 : it->second;
    }
    if (static_cast<std::size_t>(ref) < dense_.size())
      return dense_[ref];
    return -1;
  }

  void insert(int ref, int loc) {
    ++count_;
    if (useHash_) {
      hash_[ref] = loc;
      return;
    }
    if (static_cast<std::size_t>(ref) >= dense_.size()) {
      // kappa guard: if the tag space is far sparser than the vertex count,
      // a dense table wastes memory — migrate to the hash (one-time O(n)).
      const std::int64_t need = static_cast<std::int64_t>(ref) + 1;
      const std::int64_t budget =
          KAPPA * (static_cast<std::int64_t>(expected_) > count_
                       ? static_cast<std::int64_t>(expected_)
                       : count_);
      if (need > budget && need > 1024) {
        // one-line observable so G1c can assert the fallback actually engaged
        // (an oracle that cannot fail is not an oracle — N0 discipline)
        opserr << "LadrunoParallelNumberer: RefIndex -> hash fallback "
               << "(sparse tag space, ref=" << ref << ")\n";
        useHash_ = true;
        hash_.reserve(static_cast<std::size_t>(expected_) * 2u);
        for (std::size_t t = 0; t < dense_.size(); ++t)
          if (dense_[t] >= 0)
            hash_[static_cast<int>(t)] = dense_[t];
        dense_.clear();
        hash_[ref] = loc;
        return;
      }
      std::size_t newSize = dense_.empty() ? 1024 : dense_.size();
      while (static_cast<std::int64_t>(newSize) < need)
        newSize *= 2;
      dense_.resize(newSize, -1);
    }
    dense_[ref] = loc;
  }

  bool usingHash() const { return useHash_; }

 private:
  static const std::int64_t KAPPA = 4;
  std::vector<std::int32_t> dense_;
  std::unordered_map<int, int> hash_;
  bool useHash_ = false;
  std::int64_t count_ = 0;
  int expected_;
};

}  // namespace

LadrunoParallelNumberer::LadrunoParallelNumberer(GraphNumberer &theGraphNumberer)
  : ParallelNumberer(NUMBERER_TAG_LadrunoParallelNumberer, theGraphNumberer)
{

}

LadrunoParallelNumberer::LadrunoParallelNumberer()
  : ParallelNumberer(NUMBERER_TAG_LadrunoParallelNumberer)
{

}

LadrunoParallelNumberer::~LadrunoParallelNumberer()
{
  // base dtor owns theNumberer + the channel array
}

// The T0 merge: faithful to ParallelNumberer::mergeSubGraph in every mutation
// (vertex creation order, getFreeTag sequence, vertexTags/vertexRefs appends,
// theSubdomainMap layout [subTags | mergedTags], addEdge call sequence) —
// only the SEARCHES are replaced. Returns <0 on a negative ref (Lagrange
// node-less group) instead of silently fusing them like stock.
static int
ladrunoMergeSubGraph(Graph &theGraph, Graph &theSubGraph,
                     ID &vertexTags, ID &vertexRefs, ID &theSubdomainMap,
                     RefIndex &refIndex, int &numVertexInOut)
{
  Vertex *subVertexPtr;
  VertexIter &theSubGraphIter1 = theSubGraph.getVertices();
  int count = 0;
  int numVertex = numVertexInOut;
  int numVertexSub = theSubGraph.getNumVertex();

  // pass 2's subTag -> mergedTag lookup, built during pass 1 (replaces the
  // theSubdomainMap.getLocation linear scans — the pass-2 quadratic, which
  // DOMINATES at low np: cost ~ V^2/P).
  std::unordered_map<int, int> subToMerged;
  subToMerged.reserve(static_cast<std::size_t>(numVertexSub) * 2u);

  while ((subVertexPtr = theSubGraphIter1()) != 0) {
    int vertexTagSub = subVertexPtr->getTag();
    int vertexTagRef = subVertexPtr->getRef();

    if (vertexTagRef < 0) {
      opserr << "ERROR LadrunoParallelNumberer - DOF group with no node "
             << "(ref=" << vertexTagRef << ", e.g. a Lagrange multiplier "
             << "group) cannot be merged by node tag; stock code silently "
             << "fuses these across ranks (wrong numbering). Use the "
             << "Transformation/Penalty handler under MP, or extend ADR-74.\n";
      return -2;
    }

    int loc = refIndex.lookup(vertexTagRef);

    int vertexTagMerged;
    if (loc < 0) {
      // not present: create the merged vertex exactly as stock does
      vertexTagMerged = theGraph.getFreeTag();
      vertexTags[numVertex] = vertexTagMerged;
      vertexRefs[numVertex] = vertexTagRef;
      Vertex *newVertex = new Vertex(vertexTagMerged, vertexTagRef,
                                     subVertexPtr->getWeight(),
                                     subVertexPtr->getColor());
      theGraph.addVertex(newVertex);
      refIndex.insert(vertexTagRef, numVertex);
      numVertex++;
    } else
      vertexTagMerged = vertexTags[loc];

    theSubdomainMap[count] = vertexTagSub;
    theSubdomainMap[count + numVertexSub] = vertexTagMerged;
    subToMerged[vertexTagSub] = vertexTagMerged;
    count++;
  }

  // pass 2: adjacency into the merged graph — same addEdge sequence as stock,
  // O(1) remaps instead of getLocation scans.
  VertexIter &theSubGraphIter2 = theSubGraph.getVertices();
  while ((subVertexPtr = theSubGraphIter2()) != 0) {
    int vertexTagSub = subVertexPtr->getTag();
    int vertexTagMerged = subToMerged[vertexTagSub];

    const ID &adjacency = subVertexPtr->getAdjacency();
    for (int i = 0; i < adjacency.Size(); i++) {
      int vertexTagSubAdjacent = adjacency(i);
      // find(), NOT operator[]: a dangling adjacency tag (malformed subgraph)
      // must fail loudly — operator[] would default-insert 0 and silently
      // wire a spurious edge to a real vertex (review MAJOR). Stock corrupts
      // differently-but-silently here; neither is acceptable.
      std::unordered_map<int, int>::const_iterator it =
          subToMerged.find(vertexTagSubAdjacent);
      if (it == subToMerged.end()) {
        opserr << "FATAL LadrunoParallelNumberer - subgraph adjacency "
               << "references unknown vertex " << vertexTagSubAdjacent
               << " (malformed graph) — aborting.\n";
        exit(-1);
      }
      theGraph.addEdge(vertexTagMerged, it->second);
    }
  }

  numVertexInOut = numVertex;
  return 0;
}

int
LadrunoParallelNumberer::numberDOF(int lastDOF)
{
  int result = 0;

  AnalysisModel *theModel = this->getAnalysisModelPtr();
  Domain *theDomain = 0;
  if (theModel != 0) theDomain = theModel->getDomainPtr();

  if (theModel == 0 || theDomain == 0) {
    opserr << "WARNING LadrunoParallelNumberer::numberDOF(int) -";
    opserr << " - no AnalysisModel - has setLinks() been invoked?\n";
    return -1;
  }

  if (lastDOF != -1) {
    opserr << "WARNING LadrunoParallelNumberer::numberDOF(int lastDOF):";
    opserr << " does not use the lastDOF as requested\n";
  }

  Graph &theGraph = theModel->getDOFGroupGraph();

  if (processID != 0) {
    // -------- non-root: identical to stock (send graph, apply scatter) -----
    Channel *theChannel = theChannels[0];
    int numVertex = theGraph.getNumVertex();

    { OPS_PROFILE_SCOPE("dc.n.gather");
    theGraph.sendSelf(0, *theChannel);
    }

    ID theID(2*numVertex);
    { OPS_PROFILE_SCOPE("dc.n.scatter");
    theChannel->recvID(0, 0, theID);

    for (int i = 0; i < numVertex; i++) {
      int vertexTag = theID(i);
      int startID = theID(i + numVertex);
      int dofTag = vertexTag;
      DOF_Group *dofPtr;
      dofPtr = theModel->getDOF_GroupPtr(dofTag);
      if (dofPtr == 0) {
        opserr << "WARNING LadrunoParallelNumberer::numberDOF - ";
        opserr << "DOF_Group " << dofTag << "not in AnalysisModel!\n";
        result = -4;
      } else {
        const ID &theDOFID = dofPtr->getID();
        int idSize = theDOFID.Size();
        for (int j = 0; j < idSize; j++)
          if (theDOFID(j) == -2 || theDOFID(j) == -3)
            dofPtr->setID(j, startID++);
      }
    }

    theChannel->sendID(0, 0, theID);
    }

  } else {
    // -------- P0: gather, merge O(V), order, number, scatter ---------------
    int numVertex = theGraph.getNumVertex();
    int numVertexP0 = numVertex;

    ID vertexTags(numVertex);
    ID vertexRefs(numVertex);
    Vertex *vertexPtr;
    int loc = 0;
    // 64-bit sizing hint, clamped (review MAJOR-3: int*int overflows on
    // unbalanced decks; the hint only tunes the kappa budget, never correctness)
    const std::int64_t expect64 =
        static_cast<std::int64_t>(numVertex) * (numChannels + 1);
    RefIndex refIndex(expect64 > 0x7fffffff ? 0x7fffffff
                                            : static_cast<int>(expect64));
    VertexIter &theVertices = theGraph.getVertices();
    while ((vertexPtr = theVertices()) != 0) {
      int ref = vertexPtr->getRef();
      if (ref < 0) {
        // fail-STOP (AnalysisModel::getDOFGroupGraph convention): a return
        // here would leave every non-root rank blocked in recvID forever
        // while the production caller ignores the result (review MAJOR-1).
        opserr << "FATAL LadrunoParallelNumberer - P0 DOF group with no node "
               << "(ref=" << ref << ", e.g. a Lagrange multiplier group) — "
               << "see the Lagrange note in ADR-74.\n";
        exit(-1);
      }
      vertexTags[loc] = vertexPtr->getTag();
      vertexRefs[loc] = ref;
      // NB duplicate refs on P0 would be last-wins here vs stock's first-wins
      // getLocation — unreachable in-tree (one DOF group per node; the only
      // duplicate ref is -1, rejected above). Comment per review.
      refIndex.insert(ref, loc);
      loc++;
    }

    ID **theSubdomainIDs = new ID *[numChannels];
    FEM_ObjectBroker theBroker;

    for (int j = 0; j < numChannels; j++) {
      Channel *theChannel = theChannels[j];
      Graph *theSubGraph = new Graph();

      { OPS_PROFILE_SCOPE("dc.n.gather");
      theSubGraph->recvSelf(0, *theChannel, theBroker);
      }

      theSubdomainIDs[j] = new ID(theSubGraph->getNumVertex()*2);

      { OPS_PROFILE_SCOPE("dc.n.merge");
      if (ladrunoMergeSubGraph(theGraph, *theSubGraph, vertexTags, vertexRefs,
                               *theSubdomainIDs[j], refIndex, numVertex) < 0) {
        // fail-STOP: mid-collective return = distributed deadlock + a
        // partially-merged cached graph (review MAJOR-1).
        opserr << "FATAL LadrunoParallelNumberer - subgraph merge failed "
               << "(channel " << j << ") — aborting all ranks.\n";
        exit(-1);
      }
      }

      delete theSubGraph;
    }

    // ---- ordering: RCM (or any GraphNumberer) exactly as stock; the plain
    // branch replaces the theOrderedRefs->getLocation scans with seen flags
    // (this FIXES the upstream tag-0 quirk — G1b, not bit-identity) ----------
    ID *theOrderedRefs = new ID(theGraph.getNumVertex());

    { OPS_PROFILE_SCOPE("dc.n.order");
    if (theNumberer != 0) {
      *theOrderedRefs = theNumberer->number(theGraph, lastDOF);
    } else {
      // order by subdomain: subdomain 1's vertices first, then those of 2 not
      // already ordered, ...; finally P0's own not yet ordered. Merged tags
      // are getFreeTag-contiguous [0, numVertex) => dense seen flags.
      std::vector<char> seen(static_cast<std::size_t>(theGraph.getNumVertex()) + 1u, 0);
      int nOrdered = 0;
      for (int l = 0; l < numChannels; l++) {
        const ID &theSubdomain = *theSubdomainIDs[l];
        int numVertexSubdomain = theSubdomain.Size()/2;
        for (int i = 0; i < numVertexSubdomain; i++) {
          int vertexTagMerged = theSubdomain(i + numVertexSubdomain);
          if (static_cast<std::size_t>(vertexTagMerged) < seen.size() &&
              seen[vertexTagMerged] == 0) {
            seen[vertexTagMerged] = 1;
            (*theOrderedRefs)[nOrdered++] = vertexTagMerged;
          }
        }
      }
      for (int j = 0; j < numVertexP0; j++) {
        int refTagP0 = vertexTags[j];
        if (static_cast<std::size_t>(refTagP0) < seen.size() &&
            seen[refTagP0] == 0) {
          seen[refTagP0] = 1;
          (*theOrderedRefs)[nOrdered++] = refTagP0;
        }
      }
      // Live assertion of the tag-contiguity theorem the seen[] flags rest on
      // (review MAJOR-2): a silently-skipped vertex would leave trailing zeros
      // in theOrderedRefs and corrupt vertex 0 downstream.
      if (nOrdered != theGraph.getNumVertex()) {
        opserr << "FATAL LadrunoParallelNumberer - ordered " << nOrdered
               << " of " << theGraph.getNumVertex() << " vertices (DOF-group "
               << "tags not contiguous from 0?) — aborting.\n";
        exit(-1);
      }
    }

    // cumulative equation numbers onto vertex Tmp — verbatim from stock
    int count = 0;
    for (int i = 0; i < theOrderedRefs->Size(); i++) {
      int vertexTag = (*theOrderedRefs)(i);
      Vertex *vPtr = theGraph.getVertexPtr(vertexTag);
      int numDOF = vPtr->getColor();
      vPtr->setTmp(count);
      count += numDOF;
    }
    }

    delete theOrderedRefs;

    // number own dof's — verbatim from stock
    for (int i = 0; i < numVertexP0; i++) {
      int vertexTag = vertexTags(i);
      Vertex *vPtr = theGraph.getVertexPtr(vertexTag);

      int startID = vPtr->getTmp();
      int dofTag = vertexTag;
      DOF_Group *dofPtr;
      dofPtr = theModel->getDOF_GroupPtr(dofTag);
      if (dofPtr == 0) {
        opserr << "WARNING LadrunoParallelNumberer::numberDOF - ";
        opserr << "DOF_Group (P0) " << dofTag << "not in AnalysisModel!\n";
        result = -4;
      } else {
        const ID &theDOFID = dofPtr->getID();
        int idSize = theDOFID.Size();
        for (int j = 0; j < idSize; j++)
          if (theDOFID(j) == -2 || theDOFID(j) == -3)
            dofPtr->setID(j, startID++);
      }
    }

    // scatter the merged startDOFs back — verbatim from stock
    { OPS_PROFILE_SCOPE("dc.n.scatter");
    for (int k = 0; k < numChannels; k++) {
      Channel *theChannel = theChannels[k];
      ID &theSubdomain = *theSubdomainIDs[k];
      int numVertexSubdomain = theSubdomain.Size()/2;

      for (int i = 0; i < numVertexSubdomain; i++) {
        int vertexTagMerged = theSubdomain[numVertexSubdomain + i];
        Vertex *vPtr = theGraph.getVertexPtr(vertexTagMerged);
        int startDOF = vPtr->getTmp();
        theSubdomain[i + numVertexSubdomain] = startDOF;
      }

      theChannel->sendID(0, 0, theSubdomain);
      theChannel->recvID(0, 0, theSubdomain);
      delete theSubdomainIDs[k];
    }
    delete [] theSubdomainIDs;
    }
  }

  // ---- the -4 (MP-constrained) fixup: same assignments as stock, but via a
  // constrainedNode -> MPs index instead of a full MP scan per group (the
  // "second quadratic" — O(#groups x #MPs) upstream, invisible on MP-free
  // decks, hours on tie-heavy ones). Push order preserves getMPs() iteration
  // order per node, so multi-constraint application order matches stock. ----
  {
    std::unordered_map<int, std::vector<MP_Constraint*> > mpIndex;
    MP_ConstraintIter &theMPs = theDomain->getMPs();
    MP_Constraint *mpPtr;
    while ((mpPtr = theMPs()) != 0)
      mpIndex[mpPtr->getNodeConstrained()].push_back(mpPtr);

    if (!mpIndex.empty()) {
      DOF_GrpIter &tDOFs = theModel->getDOFs();
      DOF_Group *dofPtr;
      while ((dofPtr = tDOFs()) != 0) {
        const ID &theID = dofPtr->getID();
        int have4s = 0;
        for (int i = 0; i < theID.Size(); i++)
          if (theID(i) == -4) have4s = 1;

        if (have4s == 1) {
          int nodeID = dofPtr->getNodeTag();
          std::unordered_map<int, std::vector<MP_Constraint*> >::iterator it =
              mpIndex.find(nodeID);
          if (it == mpIndex.end())
            continue;
          for (std::size_t m = 0; m < it->second.size(); ++m) {
            MP_Constraint *mp = it->second[m];
            int nodeRetained = mp->getNodeRetained();
            Node *nodeRetainedPtr = theDomain->getNode(nodeRetained);
            DOF_Group *retainedDOF = nodeRetainedPtr->getDOF_GroupPtr();
            const ID &retainedDOFIDs = retainedDOF->getID();
            const ID &constrainedDOFs = mp->getConstrainedDOFs();
            const ID &retainedDOFs = mp->getRetainedDOFs();
            for (int i = 0; i < constrainedDOFs.Size(); i++) {
              int dofC = constrainedDOFs(i);
              int dofR = retainedDOFs(i);
              int dofID = retainedDOFIDs(dofR);
              dofPtr->setID(dofC, dofID);
            }
          }
        }
      }
    }
  }

  // FE elements cache their IDs — verbatim from stock
  FE_EleIter &theEle = theModel->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEle()) != 0)
    elePtr->setID();

  theModel->clearDOFGroupGraph();

  return result;
}

int
LadrunoParallelNumberer::numberDOF(ID &lastDOFs)
{
  // Upstream ParallelNumberer::numberDOF(ID&) is non-functional (the numbering
  // call is commented out and the recv layout mismatches — see ADR-74 §Where).
  // It is reachable only from the SP/DomainDecomposition lane, which ADR-74
  // scopes out. Fail loudly rather than wrap code that never worked.
  opserr << "ERROR LadrunoParallelNumberer::numberDOF(ID&) — the lastDOFs "
         << "variant is unsupported (upstream implementation is broken; "
         << "SP lane out of scope, ADR-74). Use numberDOF(int).\n";
  return -1;
}
