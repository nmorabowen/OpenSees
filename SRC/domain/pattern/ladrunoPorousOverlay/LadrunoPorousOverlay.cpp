/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoPorousOverlay implementation — see the class header and ADR 73 for the
// architecture. This file owns: the geometry snapshot, the CSR fluid operator
// (H, S*, L) assembled cell-wise through the ADR-71 kernel, the hand-rolled
// Jacobi-PCG SPD solver, the fixed-stress RATE-form fluid advance, element
// removal rescan, applyLoad force injection, -record CSV output, and full
// send/recvSelf serialization.

#include "LadrunoPorousOverlay.h"

// The ADR-71 Biot u-p block kernel + stateless shape providers — reused VERBATIM
// (no copies), relative include (their include dir is scoped to OPS_Element).
#include "../../../element/ladrunoUP/LadrunoUPKernel.h"
#include "../../../element/ladrunoUP/LadrunoUPShapes.h"

#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <LoadPatternIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>     // Ladruno (ADR-73 P3): getUndrainedAugmentation Kadd
#include <ID.h>
#include <Parameter.h>
#include <Information.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <set>

using namespace ladruno_up;

// ---------------------------------------------------------------------------
// small local helpers
// ---------------------------------------------------------------------------
namespace {

// drained bulk / shear / oedometric (constrained) modulus from (E, nu).
inline void moduliOf(double E, double nu, double& Kdr, double& G, double& Moed) {
  Kdr  = E / (3.0 * (1.0 - 2.0 * nu));
  G    = E / (2.0 * (1.0 + nu));
  Moed = Kdr + 4.0 * G / 3.0;
}

// is this element classTag a monolithic u-p (or u-p-U) element? Its region
// would double-count the pore fluid ⟨A-13⟩.
inline bool isUPElementTag(int ct) {
  return ct == ELE_TAG_LadrunoUP ||
         ct == ELE_TAG_FourNodeQuadUP ||
         ct == ELE_TAG_BrickUP ||
         ct == ELE_TAG_Nine_Four_Node_QuadUP ||
         ct == ELE_TAG_Twenty_Eight_Node_BrickUP ||
         ct == ELE_TAG_BBarFourNodeQuadUP ||
         ct == ELE_TAG_BBarBrickUP ||
         ct == ELE_TAG_EightNodeBrick_u_p_U ||
         ct == ELE_TAG_TwentyNodeBrick_u_p_U ||
         ct == ELE_TAG_EightNode_LDBrick_u_p ||
         ct == ELE_TAG_EightNode_Brick_u_p;
}

} // anonymous namespace

// ===========================================================================
//  Constructors / destructor
// ===========================================================================
LadrunoPorousOverlay::LadrunoPorousOverlay(
    int tag,
    const std::vector<int>& regionEleTags,
    const std::vector<int>& drainedNodeTags,
    double Kf, double rhoF, double biotAlpha, double Ks,
    double poroGlobal, double EGlobal, double nuGlobal,
    const double permGlobal[3], double thick,
    int stabMode, double stabValue,
    int fsLMode, double fsLScale,
    int staticMode, int onRemovalMode, double onRemovalKf,
    int fluidUpdate, int subcycleMode, int subcycleN, int dynSeepage,
    int pInitMode, double pInitGw, double pInitZw,
    const std::vector<int>& pInitListNodes,
    const std::vector<double>& pInitListVals,
    const double fluidBody[3],
    const std::vector<Layer>& layers,
    const char* recordFile, int recordEveryN)
 : LoadPattern(tag, PATTERN_TAG_LadrunoPorousOverlay),
   ndm_(0),
   regionEleTags_(regionEleTags), drainedNodeTags_(drainedNodeTags),
   Kf_(Kf), rhoF_(rhoF), biotAlpha_(biotAlpha), Ks_(Ks),
   poroGlobal_(poroGlobal), EGlobal_(EGlobal), nuGlobal_(nuGlobal),
   thick_(thick),
   stabMode_(stabMode), stabValue_(stabValue),
   fsLMode_(fsLMode), fsLScale_(fsLScale),
   staticMode_(staticMode), onRemovalMode_(onRemovalMode), onRemovalKf_(onRemovalKf),
   fluidUpdate_(fluidUpdate), subcycleMode_(subcycleMode), subcycleN_(subcycleN),
   dynSeepage_(dynSeepage),
   pInitMode_(pInitMode), pInitGw_(pInitGw), pInitZw_(pInitZw),
   pInitListNodes_(pInitListNodes), pInitListVals_(pInitListVals),
   layers_(layers),
   recordFile_(recordFile ? recordFile : ""), recordEveryN_(recordEveryN),
   snapshotOk_(false), uSnapshotValid_(false), pLoadedFromRestart_(false),
   nRegionNodes_(0),
   minEdge_(0.0), maxCv_(0.0),
   subcycleResolved_(subcycleMode == SC_FIXED),
   commitsSinceSync_(0), dtAccum_(0.0), recordCount_(0),
   seriesNoticeShown_(false), fsLFloorWarned_(false), hasSeries_(false),
   marchNoticeShown_(false),
   externalStepping_(false), pScale_(1.0), relChange_(0.0),
   moduliDirty_(false), driverSubcycleNoticeShown_(false),
   augSingularWarned_(false),
   dbTag2_(0)
{
  for (int i = 0; i < 3; i++) {
    permGlobal_[i] = permGlobal[i];
    fluidBody_[i]  = fluidBody[i];
  }
  if (subcycleN_ < 1) subcycleN_ = 1;
}

LadrunoPorousOverlay::LadrunoPorousOverlay()
 : LoadPattern(0, PATTERN_TAG_LadrunoPorousOverlay),
   ndm_(0),
   Kf_(0.0), rhoF_(0.0), biotAlpha_(1.0), Ks_(-1.0),
   poroGlobal_(0.0), EGlobal_(0.0), nuGlobal_(0.0), thick_(1.0),
   stabMode_(STAB_AUTO), stabValue_(0.25),
   fsLMode_(FSL_CLASSIC), fsLScale_(1.0),
   staticMode_(SM_MARCH), onRemovalMode_(RM_KEEP), onRemovalKf_(1.0),
   fluidUpdate_(FU_IMPLICIT), subcycleMode_(SC_FIXED), subcycleN_(1),
   dynSeepage_(0),
   pInitMode_(PI_ZERO), pInitGw_(0.0), pInitZw_(0.0),
   recordEveryN_(1),
   snapshotOk_(false), uSnapshotValid_(false), pLoadedFromRestart_(false),
   nRegionNodes_(0), minEdge_(0.0), maxCv_(0.0),
   subcycleResolved_(true), commitsSinceSync_(0), dtAccum_(0.0), recordCount_(0),
   seriesNoticeShown_(false), fsLFloorWarned_(false), hasSeries_(false),
   marchNoticeShown_(false),
   externalStepping_(false), pScale_(1.0), relChange_(0.0),
   moduliDirty_(false), driverSubcycleNoticeShown_(false),
   augSingularWarned_(false),
   dbTag2_(0)
{
  for (int i = 0; i < 3; i++) { permGlobal_[i] = 0.0; fluidBody_[i] = 0.0; }
}

LadrunoPorousOverlay::~LadrunoPorousOverlay()
{
  if (recordStream_.is_open())
    recordStream_.close();
}

// ===========================================================================
//  setDomain — triggers the geometry snapshot
// ===========================================================================
void LadrunoPorousOverlay::setDomain(Domain* theDomain)
{
  this->LoadPattern::setDomain(theDomain);
  if (theDomain == 0)                 // detached from the domain — nothing to do
    return;

  snapshotOk_ = false;
  if (this->buildSnapshot() < 0) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- geometry snapshot failed; the overlay is INVALID and every "
              "commit will ABORT the analysis until the region is corrected "
              "(see the errors above).\n";
    return;
  }
  snapshotOk_ = true;

  if (!recordFile_.empty())
    this->openRecord();
}

// ===========================================================================
//  Snapshot construction
// ===========================================================================
int LadrunoPorousOverlay::classifyCell(Element* ele, int& shapeKind, int& nNodes) const
{
  int nn = ele->getNumExternalNodes();
  Node** nds = ele->getNodePtrs();
  if (nds == 0 || nds[0] == 0) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag() << " -- region element "
           << ele->getTag() << " has unresolved nodes at snapshot time\n";
    return -1;
  }
  int cdim = nds[0]->getCrds().Size();
  nNodes = nn;
  if (cdim == 2 && nn == 3)      shapeKind = SHK_T3;
  else if (cdim == 2 && nn == 4) shapeKind = SHK_Q4;
  else if (cdim == 3 && nn == 8) shapeKind = SHK_H8;
  else {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- region element " << ele->getTag() << " is an unsupported cell "
              "(ndm=" << cdim << ", " << nn << " nodes). v1 supports only "
              "(2,3) T3, (2,4) Q4, (3,8) H8 equal-order solids\n";
    shapeKind = SHK_NONE;
    return -1;
  }
  return 0;
}

int LadrunoPorousOverlay::buildSnapshot(void)
{
  Domain* d = this->getDomain();
  if (d == 0) return -1;

  if (regionEleTags_.empty()) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- empty region (no elements). Give -region {..}\n";
    return -1;
  }

  // ---- region-overlap guard: other overlays claiming any of my elements -----
  {
    LoadPatternIter& lps = d->getLoadPatterns();
    LoadPattern* lp;
    while ((lp = lps()) != 0) {
      if (lp == this) continue;
      if (lp->getClassTag() != PATTERN_TAG_LadrunoPorousOverlay) continue;
      LadrunoPorousOverlay* other = (LadrunoPorousOverlay*)lp;
      for (size_t e = 0; e < regionEleTags_.size(); e++)
        if (other->ownsElement(regionEleTags_[e])) {
          opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
                 << " -- region element " << regionEleTags_[e]
                 << " is already claimed by LadrunoPorousOverlay " << lp->getTag()
                 << " (one overlay per hydraulically connected water body ⟨A-13⟩)\n";
          return -1;
        }
    }
  }

  // ---- build cell list + region-node numbering ------------------------------
  cells_.clear();
  regionNodeTags_.clear();
  nodeTagToIndex_.clear();
  regionNodePtrs_.clear();
  // Ladruno (ADR-73 P3): a rebuilt snapshot may keep the cell COUNT (e.g. a DB
  // restore onto the same instance) — the lazy size check alone would miss it.
  augBlocks_.clear();
  augState_.clear();
  augIndex_.clear();
  cells_.reserve(regionEleTags_.size());

  // per-element layer lookup (fatal on the same element in two layers)
  std::map<int,int> eleToLayer;
  for (size_t L = 0; L < layers_.size(); L++)
    for (size_t k = 0; k < layers_[L].eleTags.size(); k++) {
      int et = layers_[L].eleTags[k];
      if (eleToLayer.count(et)) {
        opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
               << " -- element " << et << " is assigned to two -layer blocks "
                  "(" << eleToLayer[et] << " and " << (int)L << ")\n";
        return -1;
      }
      eleToLayer[et] = (int)L;
    }

  // duplicate-in-region guard (robustness panel 1.F MAJOR): a repeated tag
  // would build TWO cells for one element → 2× stiffness/force/storage, all
  // silent. Fatal, not dedup — a repeat is always a user error.
  std::set<int> seenEle;
  for (size_t e = 0; e < regionEleTags_.size(); e++) {
    if (!seenEle.insert(regionEleTags_[e]).second) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- region element " << regionEleTags_[e] << " listed more "
                "than once in -region (would double-count that cell)\n";
      return -1;
    }
  }

  // -layer elements must live in the region (else silently material-mismatched)
  for (std::map<int,int>::iterator it = eleToLayer.begin();
       it != eleToLayer.end(); ++it)
    if (!seenEle.count(it->first)) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- -layer element " << it->first << " is not in -region\n";
      return -1;
    }

  ndm_ = 0;
  for (size_t e = 0; e < regionEleTags_.size(); e++) {
    int etag = regionEleTags_[e];
    Element* ele = d->getElement(etag);
    if (ele == 0) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- region element " << etag << " not found in the domain "
                "(define the solid elements before the pattern)\n";
      return -1;
    }
    if (isUPElementTag(ele->getClassTag())) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- region element " << etag << " is a monolithic u-p element "
                "(classTag " << ele->getClassTag() << "); it already carries the "
                "pore fluid — the overlay would double-count it ⟨A-13⟩\n";
      return -1;
    }

    int shk = SHK_NONE, nn = 0;
    if (this->classifyCell(ele, shk, nn) < 0)
      return -1;

    Node** nds = ele->getNodePtrs();
    int cdim = nds[0]->getCrds().Size();
    if (ndm_ == 0) ndm_ = cdim;
    else if (cdim != ndm_) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- region mixes dimensions: element " << etag << " is " << cdim
             << "D but the region started " << ndm_ << "D (one overlay = one ndm)\n";
      return -1;
    }

    Cell c;
    c.eleTag = etag;
    c.shapeKind = shk;
    c.nNodes = nn;
    c.alive = true;
    c.layerIndex = (eleToLayer.count(etag) ? eleToLayer[etag] : -1);

    // per-cell material params (layer override or global)
    const Layer* Ly = (c.layerIndex >= 0) ? &layers_[c.layerIndex] : 0;
    for (int i = 0; i < 3; i++)
      c.kbar[i] = (Ly ? Ly->perm[i] : permGlobal_[i]);
    c.poro = (Ly && Ly->poro > 0.0) ? Ly->poro : poroGlobal_;
    this->resolveCellModuli(c);              // c.E / c.nu (layer override or global)

    // per-cell poro vs Biot α (robustness panel 1.F MINOR): the global check
    // is in the parser, but a per-layer -poro slips past it. n > α makes the
    // storage 1/Q̄ = n/Kf + (α−n)/Ks physically wrong (and, with finite Ks,
    // can push S* indefinite → a downstream SPD abort). Catch it here, named.
    if (c.poro > biotAlpha_ + 1e-12) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- element " << etag << " porosity " << c.poro
             << " exceeds Biot alpha " << biotAlpha_
             << " (n must be <= alpha)\n";
      return -1;
    }

    // region-node numbering + resolved node pointers
    c.rIdx.resize(nn);
    c.nodes.resize(nn);
    const ID& conn = ele->getExternalNodes();
    for (int a = 0; a < nn; a++) {
      int nt = conn(a);
      std::map<int,int>::iterator it = nodeTagToIndex_.find(nt);
      int ridx;
      if (it == nodeTagToIndex_.end()) {
        ridx = (int)regionNodeTags_.size();
        nodeTagToIndex_[nt] = ridx;
        regionNodeTags_.push_back(nt);
        regionNodePtrs_.push_back(nds[a]);
      } else {
        ridx = it->second;
      }
      c.rIdx[a] = ridx;
      c.nodes[a] = nds[a];
    }
    cells_.push_back(c);
  }

  nRegionNodes_ = (int)regionNodeTags_.size();
  if (nRegionNodes_ == 0) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- region produced no nodes\n";
    return -1;
  }

  // ---- drained set (row/col elimination targets) ----------------------------
  isDrained_.assign(nRegionNodes_, 0);
  if (drainedNodeTags_.empty()) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- empty -drained set. A drained (p-fixed) boundary is required "
              "(>=1 region node); a sealed body is out of v1 scope\n";
    return -1;
  }
  for (size_t k = 0; k < drainedNodeTags_.size(); k++) {
    std::map<int,int>::iterator it = nodeTagToIndex_.find(drainedNodeTags_[k]);
    if (it == nodeTagToIndex_.end()) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- drained node " << drainedNodeTags_[k]
             << " is not a region node\n";
      return -1;
    }
    isDrained_[it->second] = 1;
  }

  // ---- CSR sparsity + block assembly ---------------------------------------
  this->buildSparsity();
  aH_.assign(colind_.size(), 0.0);
  aS_.assign(colind_.size(), 0.0);
  aL_.assign(colind_.size(), 0.0);
  fseep_.assign(nRegionNodes_, 0.0);

  minEdge_ = DBL_MAX;
  maxCv_   = 0.0;
  for (size_t ci = 0; ci < cells_.size(); ci++)
    if (this->assembleCell((int)ci) < 0)
      return -1;
  if (minEdge_ == DBL_MAX) minEdge_ = 0.0;

  // ---- fluid unknowns -------------------------------------------------------
  pTrial_.assign(nRegionNodes_, 0.0);
  if (pLoadedFromRestart_) {
    // remap loaded committed p / dp by node tag (restart-exact)
    pCommitted_.assign(nRegionNodes_, 0.0);
    dpCommitted_.assign(nRegionNodes_, 0.0);
    for (size_t k = 0; k < loadedPNodeTags_.size(); k++) {
      std::map<int,int>::iterator it = nodeTagToIndex_.find(loadedPNodeTags_[k]);
      if (it != nodeTagToIndex_.end()) {
        pCommitted_[it->second]  = loadedP_[k];
        dpCommitted_[it->second] = loadedDp_[k];
      }
    }
  } else {
    pCommitted_.assign(nRegionNodes_, 0.0);
    dpCommitted_.assign(nRegionNodes_, 0.0);
    this->applyPInit();
  }
  pTrial_ = pCommitted_;

  // uSnapshot_ = current committed region displacements. On a DB restore the
  // node state is NOT guaranteed restored yet at snapshot time (measured
  // non-deterministic divergence, 1.E-ii finding) — defer the ⟨A-10⟩
  // re-derivation to the first use (applyLoad/onDomainCommit), when the
  // committed disps are reliably the restore-instant values.
  if (pLoadedFromRestart_) {
    uSnapshot_.assign((size_t)nRegionNodes_ * ndm_, 0.0);
    uSnapshotValid_ = false;
  } else {
    this->readRegionDisp(uSnapshot_);
    uSnapshotValid_ = true;
  }

  commitsSinceSync_ = 0;
  dtAccum_ = 0.0;
  subcycleResolved_ = (subcycleMode_ == SC_FIXED);
  return 0;
}

void LadrunoPorousOverlay::buildSparsity(void)
{
  int N = nRegionNodes_;
  std::vector< std::vector<int> > rows(N);
  for (size_t ci = 0; ci < cells_.size(); ci++) {
    const Cell& c = cells_[ci];
    for (int a = 0; a < c.nNodes; a++)
      for (int b = 0; b < c.nNodes; b++)
        rows[c.rIdx[a]].push_back(c.rIdx[b]);
  }
  rowptr_.assign(N + 1, 0);
  colind_.clear();
  diagPtr_.assign(N, -1);
  for (int i = 0; i < N; i++) {
    std::vector<int>& r = rows[i];
    std::sort(r.begin(), r.end());
    r.erase(std::unique(r.begin(), r.end()), r.end());
    rowptr_[i + 1] = rowptr_[i] + (int)r.size();
    for (size_t k = 0; k < r.size(); k++) {
      if (r[k] == i) diagPtr_[i] = (int)colind_.size();
      colind_.push_back(r[k]);
    }
  }
}

// per-cell E/nu resolution — the ONE resolution logic (snapshot + rebuild).
void LadrunoPorousOverlay::resolveCellModuli(Cell& c) const
{
  const Layer* Ly = (c.layerIndex >= 0) ? &layers_[c.layerIndex] : 0;
  c.E  = (Ly && Ly->E  > 0.0) ? Ly->E  : EGlobal_;
  c.nu = (Ly && Ly->nu > 0.0) ? Ly->nu : nuGlobal_;
}

// per-cell geometry + the MODULI-INDEPENDENT blocks (H, Q, f_seep) via the
// ADR-71 kernel (templated on the shape struct). Stores Np/dNp/dv/hEdge so the
// moduli-dependent aS_/aL_ can be (re)assembled from retained data by
// assembleStorageAndL — called by assembleCell right after this, and again by
// rebuildModuliCaches on the parameter route (no shape re-evaluation).
template <class SH>
void LadrunoPorousOverlay::fillCell(Cell& c)
{
  const int nNp = SH::nNodes;
  const int ndm = SH::ndm;
  const int nGP = SH::nGP;

  // node coordinates
  double xy[SH::nNodes * SH::ndm];
  for (int a = 0; a < nNp; a++) {
    const Vector& X = c.nodes[a]->getCrds();
    for (int i = 0; i < ndm; i++)
      xy[a * ndm + i] = X(i);
  }

  c.hEdge = SH::elemSizeH(xy);                   // for -stab auto α (geometry)
  if (c.hEdge < minEdge_) minEdge_ = c.hEdge;

  double drive[3] = {fluidBody_[0], fluidBody_[1], fluidBody_[2]};

  // per-GP storage + dense cell blocks (H, Q, f_seep — no S*/L here)
  c.nGP = nGP;
  c.Np.assign(nGP * nNp, 0.0);
  c.dNp.assign(nGP * nNp * ndm, 0.0);
  c.dv.assign(nGP, 0.0);
  c.Qe.assign(nNp * ndm * nNp, 0.0);            // (nNu*ndm) x nNp

  std::vector<double> He(nNp * nNp, 0.0), fse(nNp, 0.0);

  double Nu[SH::nNodes], dNu[SH::nNodes * SH::ndm];
  double NpG[SH::nNodes], dNpG[SH::nNodes * SH::ndm];
  for (int g = 0; g < nGP; g++) {
    double detJ;
    SH::eval(g, xy, Nu, dNu, &detJ);
    if (detJ <= 0.0) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- region element " << c.eleTag << " has non-positive detJ ("
             << detJ << ") at GP " << g << " (check node ordering/geometry)\n";
      // signal via a poisoned block: leave zeros; caller aborts on the return
      // path below by re-checking. Use a sentinel:
      c.nGP = -1;
      return;
    }
    SH::evalP(g, xy, false, NpG, dNpG);         // equal-order p ≡ u
    double dv = SH::gpWeight[g] * detJ * (ndm == 2 ? thick_ : 1.0);
    for (int b = 0; b < nNp; b++) {
      c.Np[g * nNp + b] = NpG[b];
      for (int i = 0; i < ndm; i++)
        c.dNp[(g * nNp + b) * ndm + i] = dNpG[b * ndm + i];
    }
    c.dv[g] = dv;

    ladruno_up::addQ(dNu, NpG, nNp, nNp, ndm, biotAlpha_, dv, &c.Qe[0], nNp);
    ladruno_up::addH(dNpG, nNp, ndm, c.kbar, dv, &He[0], nNp);
    ladruno_up::addFseep(dNpG, nNp, ndm, c.kbar, rhoF_, drive, dv, &fse[0]);
  }

  // scatter the moduli-independent blocks into the shared CSR value arrays
  this->scatterMat(c, &He[0], aH_);
  this->scatterVecN(c, &fse[0], fseep_);
}

// (re)assemble the moduli-dependent aS_ (S + αstab·H̃) and aL_ (L) blocks for one
// cell from its retained per-GP data (Np/dNp/dv/hEdge). The ONE storage/L path:
// fillCell calls it at snapshot; rebuildModuliCaches calls it on the parameter
// route. Also owns the moduli-dependent maxCv_ reduction. Assumes the target
// CSR values were zeroed by the caller before the first cell.
void LadrunoPorousOverlay::assembleStorageAndL(Cell& c)
{
  if (c.nGP <= 0) return;                        // poisoned/empty cell
  const int nNp = c.nNodes;

  this->resolveCellModuli(c);                    // pick up E/nu (layer/global)
  double Kdr, G, Moed;
  moduliOf(c.E, c.nu, Kdr, G, Moed);

  double a2 = biotAlpha_ * biotAlpha_;
  double Lclassic = a2 / Kdr;
  double Loed     = a2 / Moed;
  double Lcoef;
  if (fsLMode_ == FSL_ZERO) {
    // Ladruno (ADR-73 P3 §3.1): the explicit lane — L ≡ 0, the shipped rate
    // form with aL_ = 0 (both L terms vanish; the reference rule is inert).
    // The oedometric floor below guards the fs1/iterated lanes only.
    Lcoef = 0.0;
  }
  else if (fsLMode_ == FSL_OEDOMETRIC) Lcoef = Loed;
  else if (fsLMode_ == FSL_MANUAL) {
    Lcoef = fsLScale_ * Lclassic;                 // $scale multiplies the classic L
    if (Lcoef < Loed) {                           // floor at oedometric ⟨A-2⟩
      Lcoef = Loed;
      if (!fsLFloorWarned_) {
        fsLFloorWarned_ = true;
        opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
               << " -- manual -fsL scale drove L below the oedometric floor "
                  "alpha^2/(K_dr+4G/3); clamped to the floor (printed once)\n";
      }
    }
  } else {
    Lcoef = Lclassic;                             // FSL_CLASSIC (default)
  }

  // consolidation coefficient for the -subcycle auto formula (worst = largest)
  double kmax = c.kbar[0];
  for (int i = 1; i < ndm_; i++) if (c.kbar[i] > kmax) kmax = c.kbar[i];
  double cv = kmax * Moed;
  if (cv > maxCv_) maxCv_ = cv;

  std::vector<double> Se, Le(nNp * nNp, 0.0);
  this->buildCellStorageBlock(c, Se);            // S* = S + αstab·H̃ (shared build)
  for (int g = 0; g < c.nGP; g++) {
    const double* NpG = &c.Np[g * nNp];
    double dv = c.dv[g];
    ladruno_up::addS(NpG, nNp, Lcoef, dv, &Le[0], nNp);      // L block
  }
  this->scatterMat(c, &Se[0], aS_);
  this->scatterMat(c, &Le[0], aL_);
}

// dense per-cell S* block (S + αstab·H̃) from retained per-GP data — the ONE S*
// definition, shared by the CSR scatter above and the ADR-73 P3 undrained
// augmentation (per-cell dense factor; never read back from CSR). Requires
// c.E/c.nu already resolved (assembleStorageAndL / rebuildModuliCaches did).
void LadrunoPorousOverlay::buildCellStorageBlock(const Cell& c,
                                                 std::vector<double>& Se) const
{
  const int nNp = c.nNodes;
  double Kdr, G, Moed;
  moduliOf(c.E, c.nu, Kdr, G, Moed);
  double ooq = ladruno_up::oneOverQbar(c.poro, Kf_, biotAlpha_, Ks_);

  double stabAlpha = 0.0;
  if (stabMode_ == STAB_AUTO)
    stabAlpha = ladruno_up::stabAlphaAuto(c.hEdge, Kdr, G, stabValue_);
  else if (stabMode_ == STAB_MANUAL)
    stabAlpha = stabValue_;

  Se.assign(nNp * nNp, 0.0);
  for (int g = 0; g < c.nGP; g++) {
    const double* NpG  = &c.Np[g * nNp];
    const double* dNpG = &c.dNp[(g * nNp) * ndm_];
    double dv = c.dv[g];
    ladruno_up::addS(NpG, nNp, ooq, dv, &Se[0], nNp);        // S storage
    if (stabAlpha != 0.0)
      ladruno_up::addHtilde(dNpG, nNp, ndm_, stabAlpha, dv, &Se[0], nNp); // + αstab H̃
  }
}

int LadrunoPorousOverlay::assembleCell(int ci)
{
  Cell& c = cells_[ci];
  switch (c.shapeKind) {
    case SHK_T3: this->fillCell<shapes::T3>(c); break;
    case SHK_Q4: this->fillCell<shapes::Q4>(c); break;
    case SHK_H8: this->fillCell<shapes::H8>(c); break;
    default:
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- internal: unclassified cell " << c.eleTag << "\n";
      return -1;
  }
  if (c.nGP < 0)      // detJ<=0 sentinel from fillCell
    return -1;
  this->assembleStorageAndL(c);   // moduli-dependent aS_/aL_ from stored per-GP data
  return 0;
}

// parameter-route moduli transport (ADR-73 P2 §2.4): re-resolve per-cell E/nu,
// re-derive α_stab / L, re-scatter aS_/aL_ from the retained per-GP snapshot
// data; aH_/fseep_/Qe untouched (perm/α only). maxCv_ recomputed. Honored
// LAZILY at the next fluid use (advanceTrial / SM_STEADY solve).
void LadrunoPorousOverlay::rebuildModuliCaches(void)
{
  std::fill(aS_.begin(), aS_.end(), 0.0);
  std::fill(aL_.begin(), aL_.end(), 0.0);
  maxCv_ = 0.0;
  for (size_t ci = 0; ci < cells_.size(); ci++)
    this->assembleStorageAndL(cells_[ci]);
  moduliDirty_ = false;
  // Ladruno (ADR-73 P3): S* changed (moduli → α-stab) — drop the cached
  // undrained-augmentation blocks (rebuilt lazily on the next request).
  // NOTE removal does NOT invalidate: -onRemoval drain touches only H/f_seep
  // (rebuildHFseep; storage untouched), and dead cells are refused by the
  // alive check in getUndrainedAugmentation — live-cell blocks stay valid.
  augBlocks_.clear();
  augState_.clear();
  augIndex_.clear();
}

int LadrunoPorousOverlay::csrPos(int i, int j) const
{
  int lo = rowptr_[i], hi = rowptr_[i + 1] - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    int cj = colind_[mid];
    if (cj == j) return mid;
    if (cj < j) lo = mid + 1; else hi = mid - 1;
  }
  return -1;
}

void LadrunoPorousOverlay::scatterMat(const Cell& c, const double* blk,
                                      std::vector<double>& tgt)
{
  int nNp = c.nNodes;
  for (int a = 0; a < nNp; a++)
    for (int b = 0; b < nNp; b++) {
      int k = this->csrPos(c.rIdx[a], c.rIdx[b]);
      if (k >= 0) tgt[k] += blk[a * nNp + b];
    }
}

void LadrunoPorousOverlay::scatterVecN(const Cell& c, const double* v,
                                       std::vector<double>& tgt)
{
  for (int b = 0; b < c.nNodes; b++)
    tgt[c.rIdx[b]] += v[b];
}

// rebuild H and f_seep after a drain-policy removal (dead cells' k̄ ×= kFactor).
// Storage S* and L are untouched (fluid measure never shrinks). Uses the stored
// per-GP data — no shape re-evaluation.
void LadrunoPorousOverlay::rebuildHFseep(void)
{
  std::fill(aH_.begin(), aH_.end(), 0.0);
  std::fill(fseep_.begin(), fseep_.end(), 0.0);
  double drive[3] = {fluidBody_[0], fluidBody_[1], fluidBody_[2]};
  for (size_t ci = 0; ci < cells_.size(); ci++) {
    Cell& c = cells_[ci];
    int nNp = c.nNodes;
    double kbarEff[3];
    for (int i = 0; i < 3; i++)
      kbarEff[i] = c.alive ? c.kbar[i] : c.kbar[i] * onRemovalKf_;
    std::vector<double> He(nNp * nNp, 0.0), fse(nNp, 0.0);
    for (int g = 0; g < c.nGP; g++) {
      const double* dNpG = &c.dNp[(g * nNp) * ndm_];
      ladruno_up::addH(dNpG, nNp, ndm_, kbarEff, c.dv[g], &He[0], nNp);
      ladruno_up::addFseep(dNpG, nNp, ndm_, kbarEff, rhoF_, drive, c.dv[g], &fse[0]);
    }
    this->scatterMat(c, &He[0], aH_);
    this->scatterVecN(c, &fse[0], fseep_);
  }
}

void LadrunoPorousOverlay::applyPInit(void)
{
  if (pInitMode_ == PI_ZERO)
    return;                                   // already zeroed

  if (pInitMode_ == PI_HYDRO) {
    // p = gamma_w * (z_w - z), clipped at 0 (z = last coordinate component)
    for (int i = 0; i < nRegionNodes_; i++) {
      const Vector& X = regionNodePtrs_[i]->getCrds();
      double z = X(ndm_ - 1);
      double p = pInitGw_ * (pInitZw_ - z);
      pCommitted_[i] = (p > 0.0) ? p : 0.0;
    }
    return;
  }

  if (pInitMode_ == PI_LIST) {
    for (size_t k = 0; k < pInitListNodes_.size(); k++) {
      std::map<int,int>::iterator it = nodeTagToIndex_.find(pInitListNodes_[k]);
      if (it != nodeTagToIndex_.end())
        pCommitted_[it->second] = pInitListVals_[k];
      else
        opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
               << " -- -pInit list node " << pInitListNodes_[k]
               << " is not a region node; ignored\n";
    }
    return;
  }

  if (pInitMode_ == PI_STEADY) {
    // solve H p = f_seep with the drained set fixed at 0 (steady seepage)
    std::vector<double> x(nRegionNodes_, 0.0);   // warm start 0; drained = 0
    if (!this->solveFluid(fseep_, 0.0, 0.0, 1.0, x))
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- -pInit steady solve did not converge; p left at 0\n";
    else
      pCommitted_ = x;
    return;
  }
}

void LadrunoPorousOverlay::readRegionDisp(std::vector<double>& u, bool trial) const
{
  u.assign((size_t)nRegionNodes_ * ndm_, 0.0);
  for (int i = 0; i < nRegionNodes_; i++) {
    // committed disp by default; trial disp under the driver latch (the mid-step
    // solve iterate). regionNodePtrs_ holds Node* (mutable target through a const
    // method), so the non-const getTrialDisp()/getDisp() accessors are legal here.
    Node* nd = regionNodePtrs_[i];
    const Vector& disp = trial ? nd->getTrialDisp() : nd->getDisp();
    int nc = disp.Size();
    for (int c = 0; c < ndm_; c++)
      u[(size_t)i * ndm_ + c] = (c < nc) ? disp(c) : 0.0;
  }
}

// ===========================================================================
//  Fluid solver:  A = cS·S* + cL·L + cH·H  (drained rows/cols eliminated)
// ===========================================================================
void LadrunoPorousOverlay::applyOp(double cS, double cL, double cH,
                                   const std::vector<double>& x,
                                   std::vector<double>& y) const
{
  int N = nRegionNodes_;
  for (int i = 0; i < N; i++) {
    double s = 0.0;
    int k0 = rowptr_[i], k1 = rowptr_[i + 1];
    for (int k = k0; k < k1; k++) {
      double a = cS * aS_[k] + cL * aL_[k] + cH * aH_[k];
      s += a * x[colind_[k]];
    }
    y[i] = s;
  }
}

bool LadrunoPorousOverlay::solveFluid(const std::vector<double>& rhsFull,
                                      double cS, double cL, double cH,
                                      std::vector<double>& x) const
{
  int N = nRegionNodes_;
  // xD carries the prescribed drained values (x on drained rows on input)
  std::vector<double> xD(N, 0.0);
  for (int i = 0; i < N; i++) if (isDrained_[i]) xD[i] = x[i];

  std::vector<double> AxD(N, 0.0);
  this->applyOp(cS, cL, cH, xD, AxD);

  std::vector<double> rhs(N, 0.0);
  double rnorm0 = 0.0;
  for (int i = 0; i < N; i++) {
    if (isDrained_[i]) { rhs[i] = 0.0; }
    else { rhs[i] = rhsFull[i] - AxD[i]; rnorm0 += rhs[i] * rhs[i]; }
  }
  rnorm0 = sqrt(rnorm0);
  if (rnorm0 < 1e-300) rnorm0 = 1.0;            // near-zero rhs → absolute test

  // Jacobi diagonal (drained rows → identity)
  std::vector<double> diag(N, 1.0);
  for (int i = 0; i < N; i++) {
    if (isDrained_[i]) { diag[i] = 1.0; continue; }
    int dk = diagPtr_[i];
    double d = (dk >= 0) ? (cS * aS_[dk] + cL * aL_[dk] + cH * aH_[dk]) : 0.0;
    if (d <= 0.0) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- non-positive fluid diagonal at region node " << i
             << " (d=" << d << "); the SPD assumption is violated\n";
      return false;
    }
    diag[i] = d;
  }

  // warm start (drained pinned to 0 in the free-space iterate)
  std::vector<double> p(N, 0.0);
  for (int i = 0; i < N; i++) p[i] = isDrained_[i] ? 0.0 : x[i];

  std::vector<double> Ap(N, 0.0), r(N, 0.0), z(N, 0.0), dir(N, 0.0);
  this->applyOp(cS, cL, cH, p, Ap);
  for (int i = 0; i < N; i++) {
    if (isDrained_[i]) { r[i] = 0.0; }
    else               { r[i] = rhs[i] - Ap[i]; }
  }
  const double tol = 1e-10;
  int maxit = 4 * N + 1000;
  bool converged = false;

  // warm start may already satisfy the system → skip CG (avoids a spurious
  // pAp==0 "breakdown" on a zero search direction).
  {
    double r0 = 0.0;
    for (int i = 0; i < N; i++) r0 += r[i] * r[i];
    if (sqrt(r0) <= tol * rnorm0) {
      for (int i = 0; i < N; i++) x[i] = isDrained_[i] ? xD[i] : p[i];
      return true;
    }
  }

  for (int i = 0; i < N; i++) z[i] = isDrained_[i] ? 0.0 : r[i] / diag[i];
  dir = z;
  double rz = 0.0;
  for (int i = 0; i < N; i++) rz += r[i] * z[i];

  for (int it = 0; it < maxit; it++) {
    this->applyOp(cS, cL, cH, dir, Ap);
    for (int i = 0; i < N; i++) if (isDrained_[i]) Ap[i] = 0.0;
    double pAp = 0.0;
    for (int i = 0; i < N; i++) pAp += dir[i] * Ap[i];
    if (pAp <= 0.0) {
      opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
             << " -- fluid CG breakdown (dir^T A dir = " << pAp
             << " <= 0); operator not SPD\n";
      return false;
    }
    double alpha = rz / pAp;
    double rnorm2 = 0.0;
    for (int i = 0; i < N; i++) {
      p[i] += alpha * dir[i];
      r[i] -= alpha * Ap[i];
      rnorm2 += r[i] * r[i];
    }
    if (sqrt(rnorm2) <= tol * rnorm0) { converged = true; break; }
    for (int i = 0; i < N; i++) z[i] = isDrained_[i] ? 0.0 : r[i] / diag[i];
    double rznew = 0.0;
    for (int i = 0; i < N; i++) rznew += r[i] * z[i];
    double beta = rznew / rz;
    for (int i = 0; i < N; i++) dir[i] = z[i] + beta * dir[i];
    rz = rznew;
  }

  if (!converged) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- fluid CG failed to reach rel-tol " << tol << " in " << maxit
           << " iterations (N=" << N << ")\n";
    return false;
  }

  for (int i = 0; i < N; i++) x[i] = isDrained_[i] ? xD[i] : p[i];
  return true;
}

void LadrunoPorousOverlay::assembleQtDu(std::vector<double>& QtDu) const
{
  QtDu.assign(nRegionNodes_, 0.0);
  std::vector<double> ucur;
  // Under the driver latch, advanceTrial runs MID-STEP (after solveCurrentStep,
  // before commit): read the solid's TRIAL disp so Δu = u_iterate − uSnapshot_ is
  // the current fixed-point increment. Unlatched (P1), advanceTrial runs at the
  // commit hook (post-commit, trial==committed) → identical to the old committed
  // read, so P1 is bit-untouched.
  this->readRegionDisp(ucur, externalStepping_);
  for (size_t ci = 0; ci < cells_.size(); ci++) {
    const Cell& c = cells_[ci];
    if (!c.alive) continue;                     // dead cells drop their Q coupling
    int nNp = c.nNodes;
    for (int b = 0; b < nNp; b++) {
      double g = 0.0;
      for (int a = 0; a < nNp; a++) {
        int ri = c.rIdx[a];
        for (int i = 0; i < ndm_; i++) {
          int row = a * ndm_ + i;
          double du = ucur[(size_t)ri * ndm_ + i] - uSnapshot_[(size_t)ri * ndm_ + i];
          g += c.Qe[row * nNp + b] * du;
        }
      }
      QtDu[c.rIdx[b]] += g;
    }
  }
}

// ===========================================================================
//  Fluid state machine — the RATE form (ADR §12 P0 pin 1)
// ===========================================================================
int LadrunoPorousOverlay::advanceTrial(bool firstOfStep, double dt)
{
  if (!snapshotOk_) return -1;
  if (moduliDirty_) this->rebuildModuliCaches();   // lazy parameter-route refresh
  if (dt <= 0.0) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- advanceTrial called with non-positive dt=" << dt << "\n";
    return -1;
  }
  int N = nRegionNodes_;

  // reference increment: first advance of a step → committed increment (rate
  // form = fs1); a repeated advance → last trial increment (fixed-point / P2).
  std::vector<double> dpref(N, 0.0);
  if (firstOfStep) dpref = dpCommitted_;
  else for (int i = 0; i < N; i++) dpref[i] = pTrial_[i] - pCommitted_[i];

  std::vector<double> SLp0(N, 0.0), Ldp(N, 0.0), QtDu(N, 0.0);
  this->applyOp(1.0, 1.0, 0.0, pCommitted_, SLp0);   // (S*+L) p0
  this->applyOp(0.0, 1.0, 0.0, dpref, Ldp);          // L Δp_ref
  this->assembleQtDu(QtDu);                          // Qᵀ Δu

  std::vector<double> rhs(N, 0.0);
  for (int i = 0; i < N; i++)
    rhs[i] = SLp0[i] + Ldp[i] - QtDu[i] + dt * fseep_[i];

  // warm start = last trial; drained rows carry the prescribed (committed) value
  std::vector<double> x = pTrial_;
  for (int i = 0; i < N; i++) if (isDrained_[i]) x[i] = pCommitted_[i];

  if (!this->solveFluid(rhs, 1.0, 1.0, dt, x))
    return -1;

  // relative change of this advance (toy `dcon` on the FREE system; the
  // driver's residual). Reference = pTrial_ BEFORE the solve (the previous
  // iterate / step-start p). Drained rows are excluded from BOTH sums: their
  // numerator term is 0 by construction, and including their (possibly large
  // prescribed) values in the denominator would under-report the residual
  // relative to the toy's reduced-system norm (panel 2.D math-6). NOTE the
  // first-of-step advance references dpCommitted_ (rate-form warm start), a
  // deliberate improvement over the toy fs_iter's cold Δp_ref = 0 first
  // iterate — same fixed point, typically fewer iterations (panel 2.D math-2).
  {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < N; i++) {
      if (isDrained_[i]) continue;
      double d = x[i] - pTrial_[i];
      num += d * d;
      den += x[i] * x[i];
    }
    relChange_ = sqrt(num) / (sqrt(den) + pScale_);
  }

  pTrial_ = x;
  return 0;
}

void LadrunoPorousOverlay::commitFluid(void)
{
  int N = nRegionNodes_;
  for (int i = 0; i < N; i++) dpCommitted_[i] = pTrial_[i] - pCommitted_[i];
  pCommitted_ = pTrial_;
  this->readRegionDisp(uSnapshot_);              // refresh Δu reference
}

void LadrunoPorousOverlay::revertFluid(void)
{
  // drop the trial; the reference rule resets to committed (next firstOfStep
  // advance references dpCommitted_ again).
  pTrial_ = pCommitted_;
}

// ===========================================================================
//  Domain-commit hook (ADR §4.2 seam)
// ===========================================================================
void LadrunoPorousOverlay::resolveSubcycleN(double dt)
{
  subcycleResolved_ = true;
  if (dt <= 0.0 || maxCv_ <= 0.0 || minEdge_ <= 0.0) { subcycleN_ = 1; return; }
  double val = 0.089 * minEdge_ * minEdge_ / (maxCv_ * dt);
  int N = (int)floor(val);
  if (N < 1) N = 1;
  subcycleN_ = N;
  opserr << "LadrunoPorousOverlay " << this->getTag()
         << ": -subcycle auto resolved N=" << subcycleN_
         << " (theta=0.089, h_min=" << minEdge_ << ", c_v_max=" << maxCv_
         << ", dt=" << dt << ")\n";
}

bool LadrunoPorousOverlay::rescanRemoval(void)
{
  Domain* d = this->getDomain();
  if (d == 0) return false;
  bool changed = false;
  for (size_t ci = 0; ci < cells_.size(); ci++) {
    Cell& c = cells_[ci];
    if (c.alive && d->getElement(c.eleTag) == 0) {
      c.alive = false;
      changed = true;
    }
  }
  if (changed && onRemovalMode_ == RM_DRAIN)
    this->rebuildHFseep();
  return changed;
}

int LadrunoPorousOverlay::onDomainCommit(double domainTime, double dT)
{
  if (!snapshotOk_) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- onDomainCommit on an INVALID overlay (snapshot failed)\n";
    return -1;
  }
  if (!uSnapshotValid_) {          // restore happened and no applyLoad yet:
    this->readRegionDisp(uSnapshot_);   // committed disps = restore instant
    uSnapshotValid_ = true;
  }
  if (fluidUpdate_ == FU_EXPLICIT) {
    opserr << "ERROR LadrunoPorousOverlay " << this->getTag()
           << " -- -fluidUpdate explicit is not implemented until P3b\n";
    return -1;
  }

  bool removal = this->rescanRemoval();

  // -staticMode: domain "time" is the load factor under static analyses ⟨A-3⟩.
  if (staticMode_ == SM_HOLD) {
    // freeze p; the applied +Q·p forces (from committed p) stay in effect.
    return 0;
  }
  if (staticMode_ == SM_STEADY) {
    if (moduliDirty_) this->rebuildModuliCaches();   // lazy parameter-route refresh
    std::vector<double> x = pCommitted_;         // warm start; drained held
    if (!this->solveFluid(fseep_, 0.0, 0.0, 1.0, x)) return -1;
    for (int i = 0; i < nRegionNodes_; i++)
      dpCommitted_[i] = x[i] - pCommitted_[i];
    pCommitted_ = x;
    pTrial_     = x;
    this->readRegionDisp(uSnapshot_);
    if (!recordFile_.empty()) {
      recordCount_++;
      if ((recordCount_ % recordEveryN_) == 0) this->recordRow(domainTime);
    }
    return 0;
  }

  // ---- transient march (SM_MARCH) -------------------------------------------
  // Driver latch (ADR-73 P2 §2.3): the fluid trial is ALREADY the converged
  // iterate (the driver ran advanceTrial to convergence this step); this commit
  // only SYNCS it. No advanceTrial (would double-march), no window accumulation.
  // rescanRemoval above already updated c.alive / H for the next step. This is
  // the ONLY commit that fires while latched.
  if (externalStepping_) {
    if (!driverSubcycleNoticeShown_ &&
        (subcycleN_ > 1 || subcycleMode_ == SC_AUTO)) {
      driverSubcycleNoticeShown_ = true;
      opserr << "LadrunoPorousOverlay " << this->getTag()
             << ": LadrunoStaggeredAnalyze overrides -subcycle while driving "
                "(single sync per driver step; printed once)\n";
    }
    this->commitFluid();
    commitsSinceSync_ = 0;
    dtAccum_ = 0.0;
    if (!recordFile_.empty()) {
      recordCount_++;
      if ((recordCount_ % recordEveryN_) == 0) this->recordRow(domainTime);
    }
    return 0;
  }

  if (subcycleMode_ == SC_AUTO && !subcycleResolved_)
    this->resolveSubcycleN(dT);

  commitsSinceSync_++;
  dtAccum_ += dT;
  bool doAdvance = removal || (commitsSinceSync_ >= subcycleN_);
  if (!doAdvance)
    return 0;                                    // accumulate Δu across the window

  double dt = dtAccum_;
  if (dt <= 0.0) {
    // non-positive step — cannot march; freeze p (no warning: a genuinely
    // static run is caught by the loud advisory below on the FIRST march;
    // a dt<=0 is just a degenerate step to skip).
    commitsSinceSync_ = 0;
    dtAccum_ = 0.0;
    return 0;
  }

  // ⟨A-3⟩ silent-static-march guard (panel 1.F MAJOR): a LoadPattern cannot see
  // the integrator, and a static LoadControl step hands us dT = Δλ > 0 (NOT
  // dt<=0 as an earlier comment wrongly assumed) — so SM_MARCH would silently
  // consolidate the fluid on load-factor pseudo-time during a gravity/initial-
  // stress stage. Convert silent-wrong to a loud one-time advisory on the first
  // real march; the user sets -staticMode hold|steady if that stage is static.
  if (!marchNoticeShown_) {
    marchNoticeShown_ = true;
    opserr << "LadrunoPorousOverlay " << this->getTag()
           << ": marching the pore fluid with dt=" << dt << " (domain-time "
              "increment). If this is a STATIC / load-factor stage (e.g. a "
              "LoadControl gravity stage), the fluid is consolidating on "
              "pseudo-time — set -staticMode hold or steady for that stage "
              "⟨A-3⟩. (printed once)\n";
  }

  if (this->advanceTrial(true, dt) < 0) return -1;   // fs1: advance ...
  this->commitFluid();                               // ... + commit back-to-back
  commitsSinceSync_ = 0;
  dtAccum_ = 0.0;

  if (!recordFile_.empty()) {
    recordCount_++;
    if ((recordCount_ % recordEveryN_) == 0) this->recordRow(domainTime);
  }
  return 0;
}

// ===========================================================================
//  Force injection — applyLoad (once per step at newStep, ⟨A-5⟩)
// ===========================================================================
void LadrunoPorousOverlay::applyLoad(double time)
{
  (void)time;
  if (!snapshotOk_) return;
  if (!uSnapshotValid_) {          // deferred post-restore baseline (⟨A-10⟩)
    this->readRegionDisp(uSnapshot_);
    uSnapshotValid_ = true;
  }

  if (hasSeries_ && !seriesNoticeShown_) {
    seriesNoticeShown_ = true;
    opserr << "LadrunoPorousOverlay " << this->getTag()
           << ": a TimeSeries/load-factor was assigned but is IGNORED by design "
              "(the overlay manages its own pore-pressure force amplitudes)\n";
  }

  // +Q·p, applied cell-wise into the nodal unbalance (no global Q). Under the
  // driver latch the iterate's TRIAL pressure drives the forces (+Q·pTrial_);
  // outside iteration pTrial_ == pCommitted_ (commitFluid syncs) so the unlatched
  // P1 path is bit-identical. Covers Domain::revertToLastCommit's re-applyLoad too
  // (harmless — unbalance is zeroed and re-applied at the next newStep).
  const std::vector<double>& pSrc = externalStepping_ ? pTrial_ : pCommitted_;
  for (size_t ci = 0; ci < cells_.size(); ci++) {
    const Cell& c = cells_[ci];
    if (!c.alive) continue;
    int nNp = c.nNodes;
    int nrows = nNp * ndm_;
    for (int a = 0; a < nNp; a++) {
      Node* nd = c.nodes[a];
      int ndf = nd->getNumberDOF();
      Vector load(ndf);
      load.Zero();
      for (int i = 0; i < ndm_; i++) {
        int row = a * ndm_ + i;
        double f = 0.0;
        for (int b = 0; b < nNp; b++)
          f += c.Qe[row * nNp + b] * pSrc[c.rIdx[b]];
        if (i < ndf) load(i) = f;
      }
      nd->addUnbalancedLoad(load);
    }
    (void)nrows;
  }
}

void LadrunoPorousOverlay::setTimeSeries(TimeSeries* theSeries)
{
  hasSeries_ = (theSeries != 0);
  this->LoadPattern::setTimeSeries(theSeries);
}

bool LadrunoPorousOverlay::ownsElement(int eleTag) const
{
  for (size_t e = 0; e < regionEleTags_.size(); e++)
    if (regionEleTags_[e] == eleTag) return true;
  return false;
}

// ===========================================================================
//  ADR-73 P3 — per-element undrained stiffness augmentation ΔK_e = Q_e·S_e⁻¹·Q_eᵀ
//  (the element-local undrained condensation; yields the DISCRETE undrained
//  pencil the Δt_cr advisory is pinned to — plan §3.2, E7.4). Cached lazily
//  per cell; invalidated by rebuildModuliCaches. Dead cells refuse (post-
//  removal the drained skeleton pencil is the honest bound for that cell).
// ===========================================================================
bool LadrunoPorousOverlay::getUndrainedAugmentation(int eleTag, Matrix& Kadd) const
{
  if (!snapshotOk_) return false;
  // parameter-route staleness: the SAME lazy-rebuild idiom advanceTrial uses —
  // a stale S* here would silently mis-price the pencil after a stage flip.
  if (moduliDirty_)
    const_cast<LadrunoPorousOverlay*>(this)->rebuildModuliCaches();

  // lazy (re)size — also covers a snapshot rebuild (cells_ changed size).
  if (augState_.size() != cells_.size()) {
    augState_.assign(cells_.size(), 0);
    augBlocks_.assign(cells_.size(), std::vector<double>());
    augIndex_.clear();
    for (size_t ci = 0; ci < cells_.size(); ci++)
      augIndex_[cells_[ci].eleTag] = (int)ci;
  }

  std::map<int,int>::const_iterator it = augIndex_.find(eleTag);
  if (it == augIndex_.end()) return false;       // not an owned cell
  const Cell& c = cells_[it->second];
  if (!c.alive || c.nGP <= 0) return false;      // dead/poisoned cell
  if (augState_[it->second] == 2) return false;  // singular S_e (already warned)

  const int nNp = c.nNodes;
  const int m   = nNp * ndm_;                    // node-major, first-ndm layout

  if (augState_[it->second] == 0) {
    // ---- factor the dense per-cell S* (SPD; small: nNp <= 8) ---------------
    std::vector<double> A;                       // S* row-major; becomes chol L
    this->buildCellStorageBlock(c, A);

    double maxDiag = 0.0;
    for (int i = 0; i < nNp; i++)
      if (A[i * nNp + i] > maxDiag) maxDiag = A[i * nNp + i];
    const double pivFloor = 1.0e-14 * maxDiag;   // relative pivot floor

    bool singular = (maxDiag <= 0.0);
    for (int k = 0; k < nNp && !singular; k++) {
      double d = A[k * nNp + k];
      for (int j = 0; j < k; j++) d -= A[k * nNp + j] * A[k * nNp + j];
      if (d <= pivFloor) { singular = true; break; }
      double Lkk = sqrt(d);
      A[k * nNp + k] = Lkk;
      for (int i = k + 1; i < nNp; i++) {
        double s = A[i * nNp + k];
        for (int j = 0; j < k; j++) s -= A[i * nNp + j] * A[k * nNp + j];
        A[i * nNp + k] = s / Lkk;
      }
    }
    if (singular) {
      augState_[it->second] = 2;
      if (!augSingularWarned_) {
        augSingularWarned_ = true;
        opserr << "LadrunoPorousOverlay " << this->getTag()
               << ": undrained Dt_cr augmentation SKIPPED for element "
               << eleTag << " -- cell storage S* is singular/ill-conditioned "
                  "(Kf -> infinity with -stab off?). The reported pencil is "
                  "the DRAINED one for such cells -- OPTIMISTIC; the material "
                  "bound sqrt(1+Kf/(n*M_oed)) is effectively infinite there. "
                  "(printed once)\n";
      }
      return false;
    }

    // ---- ΔK_e = Q_e·S_e⁻¹·Q_eᵀ: y_s = S_e⁻¹ q_s per Qe row, then dot -------
    // Qe is (nNu*ndm) x nNp row-major (ld = nNp): row r is q_r (length nNp).
    std::vector<double> Y((size_t)m * nNp, 0.0);
    for (int r = 0; r < m; r++) {
      double* y = &Y[(size_t)r * nNp];
      const double* q = &c.Qe[(size_t)r * nNp];
      for (int i = 0; i < nNp; i++) {            // forward: L z = q
        double s = q[i];
        for (int j = 0; j < i; j++) s -= A[i * nNp + j] * y[j];
        y[i] = s / A[i * nNp + i];
      }
      for (int i = nNp - 1; i >= 0; i--) {       // backward: Lᵀ y = z
        double s = y[i];
        for (int j = i + 1; j < nNp; j++) s -= A[j * nNp + i] * y[j];
        y[i] = s / A[i * nNp + i];
      }
    }
    std::vector<double>& blk = augBlocks_[it->second];
    blk.assign((size_t)m * m, 0.0);
    for (int r = 0; r < m; r++) {
      const double* q = &c.Qe[(size_t)r * nNp];
      for (int s = r; s < m; s++) {              // symmetric: fill both halves
        const double* y = &Y[(size_t)s * nNp];
        double v = 0.0;
        for (int b = 0; b < nNp; b++) v += q[b] * y[b];
        blk[(size_t)r * m + s] = v;
        blk[(size_t)s * m + r] = v;
      }
    }
    augState_[it->second] = 1;
  }

  const std::vector<double>& blk = augBlocks_[it->second];
  Kadd.resize(m, m);
  for (int r = 0; r < m; r++)
    for (int s = 0; s < m; s++)
      Kadd(r, s) = blk[(size_t)r * m + s];
  return true;
}

// ===========================================================================
//  ADR-73 P2 driver seam (LadrunoStaggeredAnalyze / LadrunoStaggeredDriver)
// ===========================================================================
void   LadrunoPorousOverlay::setExternalStepping(bool on) { externalStepping_ = on; }
bool   LadrunoPorousOverlay::externalStepping(void) const { return externalStepping_; }
void   LadrunoPorousOverlay::setResidualPScale(double s)  { pScale_ = s; }
double LadrunoPorousOverlay::lastAdvanceRelChange(void) const { return relChange_; }
int    LadrunoPorousOverlay::staticModeCode(void) const   { return staticMode_; }
int    LadrunoPorousOverlay::fluidUpdateCode(void) const  { return fluidUpdate_; }
int    LadrunoPorousOverlay::fsLModeCode(void) const      { return fsLMode_; }

// ===========================================================================
//  ADR-73 P4 recorder-channel seam — committed p at one region node tag.
//  Non-region tag => 0.0, no error (callers pre-validate; ADR §4.1 pin).
// ===========================================================================
double LadrunoPorousOverlay::pAtNodeTag(int nodeTag) const
{
  std::map<int,int>::const_iterator it = nodeTagToIndex_.find(nodeTag);
  if (it == nodeTagToIndex_.end()) return 0.0;
  int idx = it->second;
  if (idx < 0 || idx >= (int)pCommitted_.size()) return 0.0;
  return pCommitted_[(size_t)idx];
}

// Early fs1 sync of a pending -subcycle window before the driver takes over
// (plan §2.2 step 0) — otherwise the first driver advance would pair a
// multi-commit Δu window with the driver's single dt. Returns 1 if it fired,
// 0 if nothing was pending (or a degenerate dt), <0 on advance failure.
int LadrunoPorousOverlay::catchUpPendingWindow(void)
{
  if (!snapshotOk_)            return 0;
  if (commitsSinceSync_ <= 0)  return 0;         // nothing accumulated
  double dt = dtAccum_;
  if (dt <= 0.0) {                               // degenerate window: just reset
    commitsSinceSync_ = 0;
    dtAccum_ = 0.0;
    return 0;
  }
  if (this->advanceTrial(true, dt) < 0) return -1;   // fs1 sync (committed ref)
  this->commitFluid();
  commitsSinceSync_ = 0;
  dtAccum_ = 0.0;
  // mirror the window-end sync's record row (panel 2.D framework-5: a caught-up
  // window must not leave a gap in the -record CSV). Driver runs between steps,
  // so current == committed time here.
  if (!recordFile_.empty()) {
    recordCount_++;
    Domain* d = this->getDomain();
    if ((recordCount_ % recordEveryN_) == 0 && d != 0)
      this->recordRow(d->getCurrentTime());
  }
  return 1;
}

// ===========================================================================
//  Moduli / stage transport — the parameter route (ADR-73 P2 §2.4)
//    parameter $p loadPattern $overlay  E | nu | layerE $i | layerNu $i
//  updateMaterialStage does NOT and CANNOT reach a LoadPattern
//  (MaterialStageParameter registers only the first accepting ELEMENT and never
//  scans load patterns); the user re-sets overlay moduli here after a stage flip.
// ===========================================================================
int LadrunoPorousOverlay::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "E") == 0)
    return param.addObject(1, this);
  if (strcmp(argv[0], "nu") == 0)
    return param.addObject(2, this);

  // per-layer moduli: layerE $i / layerNu $i (1-based -layer declaration index)
  if (strcmp(argv[0], "layerE") == 0 || strcmp(argv[0], "layerNu") == 0) {
    if (argc < 2) {
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- parameter " << argv[0] << " needs a 1-based layer index\n";
      return -1;
    }
    int i = atoi(argv[1]);
    if (i < 1 || i > (int)layers_.size()) {
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- parameter " << argv[0] << " layer index " << i
             << " out of range (1.." << (int)layers_.size() << ")\n";
      return -1;
    }
    int base = (strcmp(argv[0], "layerE") == 0) ? 100 : 200;
    return param.addObject(base + i, this);      // encode E:100+i / nu:200+i
  }

  return -1;                                     // unrecognized target
}

int LadrunoPorousOverlay::updateParameter(int parameterID, Information& info)
{
  double val = info.theDouble;

  if (parameterID == 1) {                        // global E
    if (!(val > 0.0)) {                          // NaN-rejecting form (panel 2.D)
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- E must be > 0 (got " << val << "); ignored\n";
      return -1;
    }
    EGlobal_ = val;
    moduliDirty_ = true;
    return 0;
  }
  if (parameterID == 2) {                        // global nu
    if (!(val >= 0.0 && val < 0.5)) {            // NaN-rejecting form (panel 2.D)
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- nu must satisfy 0 <= nu < 0.5 (got " << val << "); ignored\n";
      return -1;
    }
    nuGlobal_ = val;
    moduliDirty_ = true;
    return 0;
  }
  if (parameterID >= 200) {                       // layer nu (id = 200 + i)
    int i = parameterID - 200;
    if (i < 1 || i > (int)layers_.size()) return -1;
    // NaN-rejecting form; STRICTLY positive: the layer-override storage encodes
    // "inherit global" as nu <= 0 (the P1 sentinel), so a layer nu of exactly
    // 0.0 would be silently re-interpreted as "unset" by resolveCellModuli —
    // reject it loudly instead of no-opting (panel 2.D robustness-7; quirks row).
    if (!(val > 0.0 && val < 0.5)) {
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- layerNu " << i << " must satisfy 0 < nu < 0.5 via the "
                "parameter route (nu = 0 collides with the inherit-global "
                "sentinel; set the global nu instead) (got " << val
             << "); ignored\n";
      return -1;
    }
    layers_[i - 1].nu = val;                      // establishes/updates the override
    moduliDirty_ = true;
    return 0;
  }
  if (parameterID >= 100) {                       // layer E (id = 100 + i)
    int i = parameterID - 100;
    if (i < 1 || i > (int)layers_.size()) return -1;
    if (val <= 0.0) {
      opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
             << " -- layerE " << i << " must be > 0 (got " << val
             << "); ignored\n";
      return -1;
    }
    layers_[i - 1].E = val;                       // establishes/updates the override
    moduliDirty_ = true;
    return 0;
  }
  return -1;
}

int LadrunoPorousOverlay::activateParameter(int parameterID)
{
  (void)parameterID;
  return 0;
}

// ===========================================================================
//  -record CSV
// ===========================================================================
void LadrunoPorousOverlay::openRecord(void)
{
  if (recordFile_.empty()) return;
  if (recordStream_.is_open()) recordStream_.close();
  recordStream_.open(recordFile_.c_str(), std::ios::out | std::ios::trunc);
  if (!recordStream_.is_open()) {
    opserr << "WARNING LadrunoPorousOverlay " << this->getTag()
           << " -- could not open -record file '" << recordFile_.c_str()
           << "'; recording disabled\n";
    return;
  }
  // full-f64 rows (max_digits10): the CSV doubles as the recorder-channel twin
  // oracle (ADR-73 P4 battery gate a) — default 6-sig-fig text cannot gate 1e-12
  recordStream_.precision(17);
  recordStream_ << "time";
  for (int i = 0; i < nRegionNodes_; i++)
    recordStream_ << "," << regionNodeTags_[i];
  recordStream_ << "\n";
  recordStream_.flush();
  recordCount_ = 0;
}

void LadrunoPorousOverlay::recordRow(double t)
{
  if (!recordStream_.is_open()) return;
  recordStream_ << t;
  for (int i = 0; i < nRegionNodes_; i++)
    recordStream_ << "," << pCommitted_[i];
  recordStream_ << "\n";
  recordStream_.flush();
}

// ===========================================================================
//  Serialization (WP1.D) — every ctor arg + region/layers/drained/options +
//  pCommitted_ AND dpCommitted_ (⟨A-10⟩: dpCommitted_ is NOT re-derivable;
//  uSnapshot_ is NOT serialized — re-derived from committed node disps).
// ===========================================================================
// 38 config/count slots + slot 38 = the payload Vector's OWN dbTag.
// FileDatastore keys stored objects by (dbTag, commitTag, size) per type: two
// Vectors on ONE dbTag whose sizes coincide silently overwrite each other
// (measured: a payload of exactly 38 doubles clobbered this scalar block on
// restore — 1.E-ii pInit-list crash). The payload therefore travels under its
// own channel dbTag, transmitted here (the upstream matDbTag idiom).
static const int OV_NDATA = 39;

int LadrunoPorousOverlay::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  int nEle = (int)regionEleTags_.size();
  int nDr  = (int)drainedNodeTags_.size();
  int nLay = (int)layers_.size();
  int nPI  = (int)pInitListNodes_.size();
  int nRN  = nRegionNodes_;
  int fLen = (int)recordFile_.size();
  int sumLayEle = 0;
  for (int L = 0; L < nLay; L++) sumLayEle += (int)layers_[L].eleTags.size();

  Vector data(OV_NDATA);
  data(0)  = this->getTag();
  data(1)  = ndm_;
  data(2)  = Kf_;
  data(3)  = rhoF_;
  data(4)  = biotAlpha_;
  data(5)  = Ks_;
  data(6)  = poroGlobal_;
  data(7)  = EGlobal_;
  data(8)  = nuGlobal_;
  data(9)  = permGlobal_[0];
  data(10) = permGlobal_[1];
  data(11) = permGlobal_[2];
  data(12) = thick_;
  data(13) = stabMode_;
  data(14) = stabValue_;
  data(15) = fsLMode_;
  data(16) = fsLScale_;
  data(17) = staticMode_;
  data(18) = onRemovalMode_;
  data(19) = onRemovalKf_;
  data(20) = fluidUpdate_;
  data(21) = subcycleMode_;
  data(22) = subcycleN_;
  data(23) = dynSeepage_;
  data(24) = pInitMode_;
  data(25) = pInitGw_;
  data(26) = pInitZw_;
  data(27) = fluidBody_[0];
  data(28) = fluidBody_[1];
  data(29) = fluidBody_[2];
  data(30) = recordEveryN_;
  data(31) = nEle;
  data(32) = nDr;
  data(33) = nLay;
  data(34) = nPI;
  data(35) = nRN;
  data(36) = fLen;
  data(37) = sumLayEle;
  if (dbTag2_ == 0)
    dbTag2_ = theChannel.getDbTag();        // payload Vector's own tag (see OV_NDATA note)
  data(38) = dbTag2_;
  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::sendSelf " << this->getTag()
           << " -- failed to send scalar data\n";
    return -1;
  }

  // ---- ID payload -----------------------------------------------------------
  int idSize = nEle + nDr + nPI + nLay + sumLayEle + nRN + fLen;
  ID idData(idSize > 0 ? idSize : 1);
  int q = 0;
  for (int i = 0; i < nEle; i++)  idData(q++) = regionEleTags_[i];
  for (int i = 0; i < nDr; i++)   idData(q++) = drainedNodeTags_[i];
  for (int i = 0; i < nPI; i++)   idData(q++) = pInitListNodes_[i];
  for (int L = 0; L < nLay; L++)  idData(q++) = (int)layers_[L].eleTags.size();
  for (int L = 0; L < nLay; L++)
    for (size_t k = 0; k < layers_[L].eleTags.size(); k++)
      idData(q++) = layers_[L].eleTags[k];
  for (int i = 0; i < nRN; i++)   idData(q++) = regionNodeTags_[i];
  for (int i = 0; i < fLen; i++)  idData(q++) = (int)(unsigned char)recordFile_[i];
  if (idSize > 0 && theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::sendSelf " << this->getTag()
           << " -- failed to send ID\n";
    return -1;
  }

  // ---- Vector payload -------------------------------------------------------
  int vSize = nPI + 6 * nLay + 2 * nRN;
  Vector vdata(vSize > 0 ? vSize : 1);
  int v = 0;
  for (int i = 0; i < nPI; i++) vdata(v++) = pInitListVals_[i];
  for (int L = 0; L < nLay; L++) {
    vdata(v++) = layers_[L].perm[0];
    vdata(v++) = layers_[L].perm[1];
    vdata(v++) = layers_[L].perm[2];
    vdata(v++) = layers_[L].poro;
    vdata(v++) = layers_[L].E;
    vdata(v++) = layers_[L].nu;
  }
  for (int i = 0; i < nRN; i++) vdata(v++) = pCommitted_[i];
  for (int i = 0; i < nRN; i++) vdata(v++) = dpCommitted_[i];
  if (vSize > 0 && theChannel.sendVector(dbTag2_, commitTag, vdata) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::sendSelf " << this->getTag()
           << " -- failed to send Vector payload\n";
    return -1;
  }
  return 0;
}

int LadrunoPorousOverlay::recvSelf(int commitTag, Channel& theChannel,
                                   FEM_ObjectBroker& theBroker)
{
  (void)theBroker;                          // no sub-objects to broker
  int dbTag = this->getDbTag();
  Vector data(OV_NDATA);
  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::recvSelf -- failed to recv scalar data\n";
    return -1;
  }
  this->setTag((int)data(0));
  ndm_          = (int)data(1);
  Kf_           = data(2);
  rhoF_         = data(3);
  biotAlpha_    = data(4);
  Ks_           = data(5);
  poroGlobal_   = data(6);
  EGlobal_      = data(7);
  nuGlobal_     = data(8);
  permGlobal_[0]= data(9);
  permGlobal_[1]= data(10);
  permGlobal_[2]= data(11);
  thick_        = data(12);
  stabMode_     = (int)data(13);
  stabValue_    = data(14);
  fsLMode_      = (int)data(15);
  fsLScale_     = data(16);
  staticMode_   = (int)data(17);
  onRemovalMode_= (int)data(18);
  onRemovalKf_  = data(19);
  fluidUpdate_  = (int)data(20);
  subcycleMode_ = (int)data(21);
  subcycleN_    = (int)data(22);
  dynSeepage_   = (int)data(23);
  pInitMode_    = (int)data(24);
  pInitGw_      = data(25);
  pInitZw_      = data(26);
  fluidBody_[0] = data(27);
  fluidBody_[1] = data(28);
  fluidBody_[2] = data(29);
  recordEveryN_ = (int)data(30);
  int nEle = (int)data(31);
  int nDr  = (int)data(32);
  int nLay = (int)data(33);
  int nPI  = (int)data(34);
  int nRN  = (int)data(35);
  int fLen = (int)data(36);
  int sumLayEle = (int)data(37);          // transmitted so the ID is exactly sized
  dbTag2_  = (int)data(38);               // payload Vector's own tag
  if (recordEveryN_ < 1) recordEveryN_ = 1;
  if (subcycleN_ < 1)    subcycleN_ = 1;

  // corruption guard (robustness panel 1.F MINOR): a negative/garbage count
  // from a damaged datastore would blow up the resize()/index below. Truncation
  // itself is caught by recvID/recvVector<0; this catches non-truncation garbage.
  if (nEle < 0 || nDr < 0 || nLay < 0 || nPI < 0 || nRN < 0 || fLen < 0 ||
      sumLayEle < 0) {
    opserr << "WARNING LadrunoPorousOverlay::recvSelf -- negative payload count "
              "(corrupt datastore); aborting restore\n";
    return -1;
  }

  // ---- ID payload (exact size known from the scalar block) ------------------
  int idSize = nEle + nDr + nPI + nLay + sumLayEle + nRN + fLen;
  ID idData(idSize > 0 ? idSize : 1);
  if (idSize > 0 && theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::recvSelf -- failed to recv ID\n";
    return -1;
  }
  if (idData.Size() < idSize) {           // received fewer than the header claims
    opserr << "WARNING LadrunoPorousOverlay::recvSelf -- ID payload short ("
           << idData.Size() << " < " << idSize << "); corrupt datastore\n";
    return -1;
  }
  int q = 0;
  regionEleTags_.resize(nEle);
  for (int i = 0; i < nEle; i++) regionEleTags_[i] = idData(q++);
  drainedNodeTags_.resize(nDr);
  for (int i = 0; i < nDr; i++)  drainedNodeTags_[i] = idData(q++);
  pInitListNodes_.resize(nPI);
  for (int i = 0; i < nPI; i++)  pInitListNodes_[i] = idData(q++);
  std::vector<int> layerNEle(nLay);
  for (int L = 0; L < nLay; L++) layerNEle[L] = idData(q++);
  layers_.clear();
  layers_.resize(nLay);
  for (int L = 0; L < nLay; L++) {
    layers_[L].eleTags.resize(layerNEle[L]);
    for (int k = 0; k < layerNEle[L]; k++)
      layers_[L].eleTags[k] = idData(q++);
  }
  loadedPNodeTags_.resize(nRN);
  for (int i = 0; i < nRN; i++)  loadedPNodeTags_[i] = idData(q++);
  recordFile_.clear();
  for (int i = 0; i < fLen; i++) recordFile_.push_back((char)idData(q++));

  // ---- Vector payload -------------------------------------------------------
  int vSize = nPI + 6 * nLay + 2 * nRN;
  Vector vdata(vSize > 0 ? vSize : 1);
  if (vSize > 0 && theChannel.recvVector(dbTag2_, commitTag, vdata) < 0) {
    opserr << "WARNING LadrunoPorousOverlay::recvSelf -- failed to recv Vector\n";
    return -1;
  }
  int v = 0;
  pInitListVals_.resize(nPI);
  for (int i = 0; i < nPI; i++) pInitListVals_[i] = vdata(v++);
  for (int L = 0; L < nLay; L++) {
    layers_[L].perm[0] = vdata(v++);
    layers_[L].perm[1] = vdata(v++);
    layers_[L].perm[2] = vdata(v++);
    layers_[L].poro    = vdata(v++);
    layers_[L].E       = vdata(v++);
    layers_[L].nu      = vdata(v++);
  }
  loadedP_.resize(nRN);
  for (int i = 0; i < nRN; i++)  loadedP_[i]  = vdata(v++);
  loadedDp_.resize(nRN);
  for (int i = 0; i < nRN; i++)  loadedDp_[i] = vdata(v++);

  nRegionNodes_ = nRN;
  pLoadedFromRestart_ = (nRN > 0);
  snapshotOk_ = false;                    // setDomain rebuilds geometry + remaps p

  // transient driver latch is NEVER serialized (a driver run is synchronous, so
  // no send can observe it true; pinned defensively). moduliDirty_ is transient
  // too — the restore rebuilds the snapshot from the config values anyway.
  externalStepping_ = false;
  pScale_ = 1.0;
  relChange_ = 0.0;
  moduliDirty_ = false;
  driverSubcycleNoticeShown_ = false;
  return 0;
}

// ===========================================================================
//  Print / getCopy
// ===========================================================================
void LadrunoPorousOverlay::Print(OPS_Stream& s, int flag)
{
  (void)flag;
  s << "LadrunoPorousOverlay " << this->getTag()
    << " (ADR 73, PATTERN 33022)\n";
  s << "  region elements: " << (int)regionEleTags_.size()
    << "   region nodes: " << nRegionNodes_
    << "   drained nodes: " << (int)drainedNodeTags_.size() << "\n";
  s << "  Kf=" << Kf_ << " rhoF=" << rhoF_ << " alpha=" << biotAlpha_
    << " Ks=" << Ks_ << "  moduli E=" << EGlobal_ << " nu=" << nuGlobal_ << "\n";
  const char* fsL = (fsLMode_ == FSL_CLASSIC) ? "classic" :
                    (fsLMode_ == FSL_OEDOMETRIC) ? "oedometric" :
                    (fsLMode_ == FSL_ZERO) ? "zero" : "manual";
  const char* sm  = (staticMode_ == SM_MARCH) ? "march" :
                    (staticMode_ == SM_HOLD) ? "hold" : "steady";
  const char* rm  = (onRemovalMode_ == RM_KEEP) ? "keep" : "drain";
  s << "  fsL=" << fsL << " staticMode=" << sm << " onRemoval=" << rm
    << " subcycleN=" << subcycleN_ << " snapshotOk=" << (snapshotOk_ ? 1 : 0) << "\n";
}

LoadPattern* LadrunoPorousOverlay::getCopy(void)
{
  LadrunoPorousOverlay* theCopy = new LadrunoPorousOverlay(
      this->getTag(), regionEleTags_, drainedNodeTags_,
      Kf_, rhoF_, biotAlpha_, Ks_, poroGlobal_, EGlobal_, nuGlobal_,
      permGlobal_, thick_, stabMode_, stabValue_, fsLMode_, fsLScale_,
      staticMode_, onRemovalMode_, onRemovalKf_, fluidUpdate_,
      subcycleMode_, subcycleN_, dynSeepage_, pInitMode_, pInitGw_, pInitZw_,
      pInitListNodes_, pInitListVals_, fluidBody_, layers_,
      recordFile_.c_str(), recordEveryN_);
  return theCopy;
}
