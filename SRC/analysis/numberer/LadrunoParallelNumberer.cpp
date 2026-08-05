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
// N3/T1: the merge + ordering run on owned dense structures.
//
// The algorithm is a faithful transcription of ParallelNumberer::numberDOF(int)
// — same gather, same vertex-creation order, same getFreeTag sequence, same
// append order, same RCM input, same scatter — with every O(V²) linear scan
// replaced by an O(1) lookup:
//   pass 1 dedup:  vertexRefs.getLocation(ref)        -> RefIndex (dense/hash)
//   pass 2 remap:  theSubdomainMap.getLocation(tag)   -> per-subgraph hash map
//   plain branch:  theOrderedRefs->getLocation(tag)   -> seen[] flags
//   -4 fixup:      full MP scan per constrained group -> constrainedNode index
// and (T1) every std::map pointer-chase on the merged graph replaced by dense
// array access:
//   the merged graph lives on ArrayOfTaggedObjects storage owned HERE (dense
//   tags sit at position == tag => Graph::getVertexPtr is O(1), so RCM's BFS
//   reads — 1-2 per directed adjacency entry — stop being map finds), and a
//   tag -> Vertex* vector is maintained through the merge so edge insertion
//   calls Vertex::addEdge on both endpoints directly (no Graph::addEdge
//   lookups at all). The model's own map-backed DOF-group graph is left
//   UNMERGED — it only sources P0's vertices, and clearDOFGroupGraph() frees
//   it at the end exactly as before.
// Bit-identity for the RCM verb is a THEOREM of "lookups only" (identical
// merged graph => identical RCM => identical numbers) and is ENFORCED by gate
// G1 (tests/test_adr74_numberer_1.py). For T1 the theorem needs three more
// legs, each verified in source: adjacency is a sorted ID::insert set
// (insertion-order canonical), getFreeTag sequence is preserved (nextFreeTag
// == numVertexP0 after the P0 copy, exactly stock's state at merge start),
// and vertex iteration order is ascending-tag in BOTH storages for dense
// tags (RCM's start vertex + disconnected-component restarts depend on it).
// The plain branch intentionally fixes the upstream tag-0 quirk (vertex 0
// falsely "already ordered" against a zero-filled ID, numbered last) => G1b
// gates it by bijection + same-physics, NOT bit-identity. Measured baseline
// T0 replaced: dc.numberDOF ~ N^2.01 (39.9 / 578.9 / 2477.5 s over
// 0.25/1/2 M hex, np8; 16.35 h at 19.18 M np240); T1 targets the ~N^1.1 map
// residual left after T0 (16.5 s at 2 M np8: merge 7.9 + order 5.2).
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
#include <ArrayOfTaggedObjects.h>   // T1: dense storage for the merged graph
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

// T1 edge insertion: what Graph::addEdge does, minus its two getVertexPtr
// lookups — both endpoints come from the owned tag -> Vertex* vector.
// Vertex::addEdge is a sorted ID::insert (returns 0 added / 1 already there),
// with the self-edge no-op preserved because we call the SAME member function
// stock calls. Asymmetric adjacency (one side present, the other not) is the
// same corruption stock exits on — kept as a fatal check.
static inline void
ladrunoAddEdgeDirect(std::vector<Vertex*> &tagToVertex, int tag1, int tag2)
{
  int r = tagToVertex[static_cast<std::size_t>(tag1)]->addEdge(tag2);
  if (r == 0) {
    if (tagToVertex[static_cast<std::size_t>(tag2)]->addEdge(tag1) != 0) {
      opserr << "FATAL LadrunoParallelNumberer - asymmetric adjacency "
             << tag1 << " <-> " << tag2 << " (corrupt graph) — aborting.\n";
      exit(-1);
    }
  } else if (r < 0) {
    opserr << "FATAL LadrunoParallelNumberer - Vertex::addEdge(" << tag1
           << ", " << tag2 << ") failed — aborting.\n";
    exit(-1);
  }
}

// The T0 merge: faithful to ParallelNumberer::mergeSubGraph in every mutation
// (vertex creation order, getFreeTag sequence, vertexTags/vertexRefs appends,
// theSubdomainMap layout [subTags | mergedTags], addEdge call sequence) —
// only the SEARCHES are replaced. T1: the target graph is the owned
// dense-storage merged graph, new vertices are registered in tagToVertex,
// and pass-2 edges go in via ladrunoAddEdgeDirect (no Graph::addEdge
// lookups). Returns <0 on a negative ref (Lagrange node-less group) instead
// of silently fusing them like stock.
static int
ladrunoMergeSubGraph(Graph &mergedGraph, std::vector<Vertex*> &tagToVertex,
                     Graph &theSubGraph,
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
      vertexTagMerged = mergedGraph.getFreeTag();
      vertexTags[numVertex] = vertexTagMerged;
      vertexRefs[numVertex] = vertexTagRef;
      Vertex *newVertex = new Vertex(vertexTagMerged, vertexTagRef,
                                     subVertexPtr->getWeight(),
                                     subVertexPtr->getColor());
      mergedGraph.addVertex(newVertex);
      // tag contiguity is the theorem tagToVertex indexing rests on —
      // live-checked (N2 discipline: an assumption without an assertion
      // is a latent silent-wrong-answer).
      if (static_cast<std::size_t>(vertexTagMerged) != tagToVertex.size()) {
        opserr << "FATAL LadrunoParallelNumberer - merged tag "
               << vertexTagMerged << " != expected "
               << static_cast<int>(tagToVertex.size())
               << " (getFreeTag not contiguous?) — aborting.\n";
        exit(-1);
      }
      tagToVertex.push_back(newVertex);
      refIndex.insert(vertexTagRef, numVertex);
      numVertex++;
    } else
      vertexTagMerged = vertexTags[loc];

    theSubdomainMap[count] = vertexTagSub;
    theSubdomainMap[count + numVertexSub] = vertexTagMerged;
    subToMerged[vertexTagSub] = vertexTagMerged;
    count++;
  }

  // pass 2: adjacency into the merged graph — same edge-insertion sequence as
  // stock, O(1) remaps instead of getLocation scans, direct Vertex::addEdge
  // instead of Graph::addEdge's per-call vertex lookups (T1).
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
      ladrunoAddEdgeDirect(tagToVertex, vertexTagMerged, it->second);
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

    // ---- T1: the merged graph lives HERE, on dense array-backed storage,
    // with a tag -> Vertex* mirror. The Graph ctor takes ownership of the
    // storage; the Graph dtor deletes the copied vertices with it. The model's
    // map-backed group graph is only read (P0 vertex source) — it is never
    // mutated, and clearDOFGroupGraph() disposes of it at the end as before.
    const int denseHint = expect64 > 0x7fffffff ? 0x7fffffff
                          : (expect64 > 0 ? static_cast<int>(expect64) : 1024);
    Graph *mergedGraph = new Graph(*(new ArrayOfTaggedObjects(denseHint)));
    std::vector<Vertex*> tagToVertex;
    tagToVertex.reserve(static_cast<std::size_t>(denseHint));

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
      int p0Tag = vertexPtr->getTag();
      // T1 copy of the P0 vertex (tag/ref/weight/color; adjacency below once
      // all copies exist). The group graph's tags are 0..numVertexP0-1 in
      // ascending iteration order — the theorem tagToVertex indexing (and the
      // stock getFreeTag sequence, nextFreeTag == numVertexP0 after this
      // loop) rests on; live-checked.
      if (p0Tag != loc) {
        opserr << "FATAL LadrunoParallelNumberer - P0 group-graph tag "
               << p0Tag << " != iteration position " << loc
               << " (tags not contiguous-ascending?) — aborting.\n";
        exit(-1);
      }
      Vertex *copyVertex = new Vertex(p0Tag, ref, vertexPtr->getWeight(),
                                      vertexPtr->getColor());
      mergedGraph->addVertex(copyVertex);
      tagToVertex.push_back(copyVertex);
      vertexTags[loc] = p0Tag;
      vertexRefs[loc] = ref;
      // NB duplicate refs on P0 would be last-wins here vs stock's first-wins
      // getLocation — unreachable in-tree (one DOF group per node; the only
      // duplicate ref is -1, rejected above). Comment per review.
      refIndex.insert(ref, loc);
      loc++;
    }

    // T1: replicate P0's own adjacency onto the copies (sorted ID::insert is
    // insertion-order canonical, so per-edge direct insertion reproduces the
    // exact adjacency sets stock starts the merge with).
    { OPS_PROFILE_SCOPE("dc.n.merge");
    VertexIter &theVerticesAdj = theGraph.getVertices();
    while ((vertexPtr = theVerticesAdj()) != 0) {
      int p0Tag = vertexPtr->getTag();
      Vertex *copyVertex = tagToVertex[static_cast<std::size_t>(p0Tag)];
      const ID &adjacency = vertexPtr->getAdjacency();
      for (int i = 0; i < adjacency.Size(); i++)
        copyVertex->addEdge(adjacency(i));   // symmetric source => symmetric copy
    }
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
      if (ladrunoMergeSubGraph(*mergedGraph, tagToVertex, *theSubGraph,
                               vertexTags, vertexRefs,
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

    // ---- ordering: RCM (or any GraphNumberer) exactly as stock, but running
    // on the OWNED dense-storage merged graph (T1) so every BFS getVertexPtr
    // is an array index, not a map find; the plain branch replaces the
    // theOrderedRefs->getLocation scans with seen flags (this FIXES the
    // upstream tag-0 quirk — G1b, not bit-identity) --------------------------
    ID *theOrderedRefs = new ID(mergedGraph->getNumVertex());

    { OPS_PROFILE_SCOPE("dc.n.order");
    if (theNumberer != 0) {
      *theOrderedRefs = theNumberer->number(*mergedGraph, lastDOF);
    } else {
      // order by subdomain: subdomain 1's vertices first, then those of 2 not
      // already ordered, ...; finally P0's own not yet ordered. Merged tags
      // are getFreeTag-contiguous [0, numVertex) => dense seen flags.
      std::vector<char> seen(static_cast<std::size_t>(mergedGraph->getNumVertex()) + 1u, 0);
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
      if (nOrdered != mergedGraph->getNumVertex()) {
        opserr << "FATAL LadrunoParallelNumberer - ordered " << nOrdered
               << " of " << mergedGraph->getNumVertex() << " vertices (DOF-group "
               << "tags not contiguous from 0?) — aborting.\n";
        exit(-1);
      }
    }

    // cumulative equation numbers onto vertex Tmp — same assignments as stock,
    // O(1) tagToVertex reads on the owned copies (T1; the Tmp values written
    // here live on the merged-graph copies, so the loops below must read the
    // SAME copies — they do, via tagToVertex).
    int count = 0;
    for (int i = 0; i < theOrderedRefs->Size(); i++) {
      int vertexTag = (*theOrderedRefs)(i);
      Vertex *vPtr = tagToVertex[static_cast<std::size_t>(vertexTag)];
      int numDOF = vPtr->getColor();
      vPtr->setTmp(count);
      count += numDOF;
    }
    }

    delete theOrderedRefs;

    // number own dof's — same assignments as stock, tagToVertex reads (T1)
    for (int i = 0; i < numVertexP0; i++) {
      int vertexTag = vertexTags(i);
      Vertex *vPtr = tagToVertex[static_cast<std::size_t>(vertexTag)];

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
        Vertex *vPtr = tagToVertex[static_cast<std::size_t>(vertexTagMerged)];
        int startDOF = vPtr->getTmp();
        theSubdomain[i + numVertexSubdomain] = startDOF;
      }

      theChannel->sendID(0, 0, theSubdomain);
      theChannel->recvID(0, 0, theSubdomain);
      delete theSubdomainIDs[k];
    }
    delete [] theSubdomainIDs;
    }

    // T1: the owned merged graph (and every copied vertex with it) dies here;
    // tagToVertex pointers dangle past this line by construction — nothing
    // below touches them.
    delete mergedGraph;
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
