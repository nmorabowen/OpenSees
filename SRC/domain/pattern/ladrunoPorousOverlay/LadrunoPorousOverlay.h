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
// LadrunoPorousOverlay — a persistent pore-fluid pressure field with its own
// life-cycle, decoupled from the solid elements it couples to (ADR 73,
// PATTERN_TAG 33022). A LoadPattern subclass that owns a pore-pressure field
// OUTSIDE the Domain's DOF graph: its own p vector, its own H/S*/H̃/L blocks
// (assembled once from LadrunoUPKernel.h over a GEOMETRY SNAPSHOT of a named
// region of solid elements), its own drained-BC set, and its own small SPD
// CG solver. The skeleton never sees p as a DOF — it feels only the nodal
// force +Q·p, applied through applyLoad() via node->addUnbalancedLoad() (the
// verified H5DRM mechanism, H5DRMLoadPattern.cpp:897).
//
// Coupling = fixed-stress staggered split, solid-first (ADR §3.1). Per step the
// solid analysis runs with the overlay's frozen forces +Q·pⁿ; at Domain::commit
// the sibling hook calls onDomainCommit(), which advances the fluid system with
// the FIXED-STRESS RATE FORM (ADR §12 P0 pin 1 — governs over §3.1's drafted
// fixed-point wording):
//
//   (S* + L + Δt·H)·p_trial = (S* + L)·p₀ + L·Δp_ref − Qᵀ·Δu + Δt·f_seep
//
// with Δp_ref = the previous COMMITTED fluid increment (dpCommitted_) on the
// first advance of a step (rate form = fs1, O(Δt), single pass). A repeated
// advance within one step references the last TRIAL increment (pTrial−p₀ → the
// §3.1 fixed-point iteration; L cancels at convergence = monolithic BE) — that
// is the P2 iterated-driver path, wired here but unused by v1's plain-analyze
// flow. Never run the fixed-POINT form single-pass (O(1) storage error,
// measured 41–61 %).
//
// v1 lane: fs1 single-pass implicit. -fluidUpdate explicit is fatal-NYI (P3b).
// See Ladruno_implementation/73_ladruno_porous_overlay_adr.md and
// Ladruno_implementation/_adr73_implementation_plan.md (WP1 pins).

#ifndef LadrunoPorousOverlay_h
#define LadrunoPorousOverlay_h

#include <LoadPattern.h>
#include <Vector.h>
#include <ID.h>

#include <vector>
#include <map>
#include <string>
#include <fstream>

class Domain;
class Node;
class Element;
class Channel;
class FEM_ObjectBroker;
class OPS_Stream;

class LadrunoPorousOverlay : public LoadPattern
{
 public:
  // ---- option enumerations (also the serialized integer codes) --------------
  enum StabMode     { STAB_OFF = 0, STAB_AUTO = 1, STAB_MANUAL = 2 };
  enum FsLMode      { FSL_CLASSIC = 0, FSL_OEDOMETRIC = 1, FSL_MANUAL = 2 };
  // SM_MARCH is the default (flag absent): the overlay marches the fluid at each
  // qualifying transient commit using dT. The -staticMode flag selects HOLD
  // (freeze p, forces still applied from committed p) or STEADY (re-solve
  // H·p = f_seep) for static stages, where the domain "time" is the load factor
  // and a dT-march would be silently wrong physics (ADR §8 ⟨A-3⟩). Detection is
  // the FLAG, never the integrator.
  enum StaticMode   { SM_MARCH = 0, SM_HOLD = 1, SM_STEADY = 2 };
  enum RemovalMode  { RM_KEEP = 0, RM_DRAIN = 1 };
  enum FluidUpdate  { FU_IMPLICIT = 0, FU_EXPLICIT = 1 };  // explicit = P3b NYI
  enum SubcycleMode { SC_FIXED = 0, SC_AUTO = 1 };
  enum PInitMode    { PI_ZERO = 0, PI_STEADY = 1, PI_HYDRO = 2, PI_LIST = 3 };

  // per-layer perm / porosity / moduli override (ADR §4.1 -layer). A negative
  // poro / E / nu means "use the overlay-global value".
  struct Layer {
    double perm[3];
    double poro;              // <=0 => use global
    double E, nu;             // <=0 => use global
    std::vector<int> eleTags;
    Layer() : poro(-1.0), E(-1.0), nu(-1.0) { perm[0] = perm[1] = perm[2] = 0.0; }
  };

  // ---- constructors ---------------------------------------------------------
  LadrunoPorousOverlay(int tag,
                       const std::vector<int>& regionEleTags,
                       const std::vector<int>& drainedNodeTags,
                       double Kf, double rhoF, double biotAlpha, double Ks,
                       double poroGlobal, double EGlobal, double nuGlobal,
                       const double permGlobal[3], double thick,
                       int stabMode, double stabValue,
                       int fsLMode, double fsLScale,
                       int staticMode, int onRemovalMode, double onRemovalKf,
                       int fluidUpdate, int subcycleMode, int subcycleN,
                       int dynSeepage,
                       int pInitMode, double pInitGw, double pInitZw,
                       const std::vector<int>& pInitListNodes,
                       const std::vector<double>& pInitListVals,
                       const double fluidBody[3],
                       const std::vector<Layer>& layers,
                       const char* recordFile, int recordEveryN);
  LadrunoPorousOverlay();                    // for FEM_ObjectBroker
  ~LadrunoPorousOverlay();

  // ---- LoadPattern overrides ------------------------------------------------
  void   setDomain(Domain* theDomain);
  void   applyLoad(double time);
  void   setTimeSeries(TimeSeries* theSeries);   // ignored-by-design; one notice
  int    sendSelf(int commitTag, Channel& theChannel);
  int    recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void   Print(OPS_Stream& s, int flag = 0);
  LoadPattern* getCopy(void);

  // ---- ADR-73 seam (called by the sibling's guarded Domain::commit block) ---
  // 0 on success, negative = fatal (printed loudly). Honors -subcycle, applies
  // -staticMode semantics, rescans element removal, then advanceTrial+commitFluid.
  int    onDomainCommit(double domainTime, double dT);

  // fluid state-machine (public: the hook drives fs1 back-to-back; the P2
  // iterated driver / domain-revert path call these directly — ADR §3.1 pins).
  int    advanceTrial(bool firstOfStep, double dt);
  void   commitFluid(void);
  void   revertFluid(void);

  // region-overlap discovery ⟨A-13⟩: does this overlay claim eleTag?
  bool   ownsElement(int eleTag) const;

 private:
  struct Cell;                               // defined below (used in decls)

  // ---- geometry snapshot ----------------------------------------------------
  int    buildSnapshot(void);                // 0 ok, <0 fatal (prints)
  int    classifyCell(Element* ele, int& shapeKind, int& nNodes) const;
  int    assembleCell(int ci);               // per-cell blocks via the kernel
  template <class SH> void fillCell(Cell& c);// templated per-shape assembly
  void   buildSparsity(void);                // CSR pattern from region connectivity
  int    csrPos(int i, int j) const;         // CSR index of (i,j) or -1
  void   scatterMat(const Cell& c, const double* blk, std::vector<double>& tgt);
  void   scatterVecN(const Cell& c, const double* v, std::vector<double>& tgt);
  void   rebuildHFseep(void);                // reassemble H/f_seep (drain-removal)
  void   applyPInit(void);                   // -pInit: zero/steady/hydro/list
  void   readRegionDisp(std::vector<double>& u) const;  // committed region u
  bool   rescanRemoval(void);                // true if any cell newly died
  void   resolveSubcycleN(double dt);        // -subcycle auto formula

  // ---- fluid solver (CSR operator A = cS·S* + cL·L + cH·H) -------------------
  void   applyOp(double cS, double cL, double cH,
                 const std::vector<double>& x, std::vector<double>& y) const;
  bool   solveFluid(const std::vector<double>& rhsFull,
                    double cS, double cL, double cH,
                    std::vector<double>& x) const;   // x: warm start in / soln out
  void   assembleQtDu(std::vector<double>& QtDu) const;   // Qᵀ·(u−uSnapshot_)

  void   recordRow(double t);
  void   openRecord(void);

  // ---- one cell of the region ----------------------------------------------
  struct Cell {
    int    eleTag;
    int    shapeKind;              // SHK_T3 / SHK_Q4 / SHK_H8
    int    nNodes;                 // = nNp = nNu
    int    nGP;
    bool   alive;
    int    layerIndex;            // -1 = global
    double kbar[3];
    double poro;
    double E, nu;
    std::vector<int>    rIdx;      // local node -> region-node index
    std::vector<Node*>  nodes;     // resolved at snapshot
    std::vector<double> Np;        // nGP*nNp
    std::vector<double> dNp;       // nGP*nNp*ndm
    std::vector<double> dv;        // nGP
    std::vector<double> Qe;        // (nNu*ndm) x nNp, row-major, ld = nNp
    Cell() : eleTag(0), shapeKind(0), nNodes(0), nGP(0), alive(true),
             layerIndex(-1), poro(0.0), E(0.0), nu(0.0)
    { kbar[0] = kbar[1] = kbar[2] = 0.0; }
  };

  // ---- configuration (mirrors the ctor; serialized in full) -----------------
  int    ndm_;                               // resolved from region node coords
  std::vector<int>    regionEleTags_;
  std::vector<int>    drainedNodeTags_;
  double Kf_, rhoF_, biotAlpha_, Ks_;
  double poroGlobal_, EGlobal_, nuGlobal_;
  double permGlobal_[3];
  double thick_;
  int    stabMode_;   double stabValue_;
  int    fsLMode_;    double fsLScale_;
  int    staticMode_;
  int    onRemovalMode_;  double onRemovalKf_;
  int    fluidUpdate_;
  int    subcycleMode_;   int subcycleN_;
  int    dynSeepage_;
  int    pInitMode_;  double pInitGw_, pInitZw_;
  std::vector<int>    pInitListNodes_;
  std::vector<double> pInitListVals_;
  double fluidBody_[3];
  std::vector<Layer>  layers_;
  std::string recordFile_;   int recordEveryN_;

  // ---- snapshot / analysis state -------------------------------------------
  bool   snapshotOk_;
  bool   pLoadedFromRestart_;                // recvSelf loaded p; skip pInit
  int    nRegionNodes_;
  std::vector<int>    regionNodeTags_;
  std::vector<Node*>  regionNodePtrs_;       // resolved at snapshot (size N)
  std::map<int,int>   nodeTagToIndex_;
  std::vector<Cell>   cells_;
  double minEdge_, maxCv_;                   // -subcycle auto inputs

  // restart carriers (recvSelf → remapped in buildSnapshot by node tag)
  std::vector<int>    loadedPNodeTags_;
  std::vector<double> loadedP_, loadedDp_;

  // CSR (region-node × region-node; single shared pattern for H, S*, L)
  std::vector<int>    rowptr_;
  std::vector<int>    colind_;
  std::vector<int>    diagPtr_;              // per-row index of the diagonal entry
  std::vector<double> aH_, aS_, aL_;         // H, S* (=S+αstab·H̃), L values
  std::vector<double> fseep_;                // seepage source vector (size N)
  std::vector<char>   isDrained_;            // size N

  // fluid unknowns (size N) — internal double storage (packed to Vector on wire)
  std::vector<double> pCommitted_;
  std::vector<double> dpCommitted_;          // previous committed increment (rate ref)
  std::vector<double> pTrial_;
  std::vector<double> uSnapshot_;            // region u at last fluid advance (N*ndm)

  // subcycle / record bookkeeping
  bool   subcycleResolved_;
  int    commitsSinceSync_;
  double dtAccum_;
  int    recordCount_;
  std::ofstream recordStream_;

  // one-time notices
  bool   seriesNoticeShown_;
  bool   fsLFloorWarned_;
  bool   hasSeries_;

  // shape-kind codes
  enum { SHK_NONE = 0, SHK_T3 = 1, SHK_Q4 = 2, SHK_H8 = 3 };

  // non-copyable (owns an ofstream + snapshot); use getCopy() for a config clone
  LadrunoPorousOverlay(const LadrunoPorousOverlay&);
  LadrunoPorousOverlay& operator=(const LadrunoPorousOverlay&);
};

#endif
