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
// Fluid-update lanes (-fluidUpdate): the default fs1 single-pass IMPLICIT lane
// above, and the P3b EXPLICIT lane (toy-exact vs the E7.4 oracle):
// lumped S* (sLump_ = row-sum of the assembled S*; H̃·1 = 0 so the stab matrix
// cannot enter the lump) and a forward p-step from COMMITTED state,
//
//   p_trial = p₀ + (−Δt·H·p₀ − Qᵀ·Δu + Δt·f_seep) / sLump    (free rows;
//                                                drained rows hold committed)
//
// — no CG, no factorization in the transient march. L is INERT on this lane
// (no iteration; -fsL gets a one-time notice). SM_STEADY keeps its CG solve at
// static commits (setup/gravity-init, not a march cost). The dual-CFL advisory
// (explicitFluidDiffusionDt below) bounds the -subcycle window: N·Δt ≤ Δt_diff.
//
// SEQUENCING (ADR-73 §12 P3b entry — the hook-time pin was REFUTED, measured):
// the explicit lane advances at LOAD-APPLICATION time (applyLoad), on the
// TRIAL window Δu, and the Domain::commit hook only SYNCS. Why: the CD-family
// integrators (CentralDifferenceLadruno) advance u in newStep with the
// PREVIOUS step's solve acceleration, so a hook-time (post-commit) fluid
// advance pairs the solid with a TWO-commit-stale p — one force lag beyond the
// toy/ZPC contract, and the explicit fluid is unconditionally unstable with it
// (growth ∝ dt, no stable dt; toy replica reproduces the C++ divergence to 7
// digits with exactly one added lag step). Advancing inside applyLoad — after
// newStep set the trial u, before the residual forms — restores the toy
// pairing. The advance stays FROM COMMITTED and idempotent, so repeated
// applyLoads within one step (Noh-Bathe sub-steps, revert re-applies) simply
// recompute the same window with the latest trial Δu; the last one before
// commit wins. CAVEAT: under an IMPLICIT integrator the applyLoad-time trial
// is the predictor state, so the explicit fluid would pair with predictor Δu —
// the CD family is the supported lane (the guide owns this).
// See Ladruno_implementation/73_ladruno_porous_overlay_adr.md and
// Ladruno_implementation/_adr73_implementation_plan.md (WP1 + §3b pins).

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
class Parameter;
class Information;
class Matrix;

class LadrunoPorousOverlay : public LoadPattern
{
 public:
  // ---- option enumerations (also the serialized integer codes) --------------
  enum StabMode     { STAB_OFF = 0, STAB_AUTO = 1, STAB_MANUAL = 2 };
  // FSL_ZERO (ADR-73 P3 §3.1): L ≡ 0 — the EXPLICIT-lane setting (CD at
  // Δt ≤ 0.5× the discrete undrained pencil). Bypasses the oedometric floor
  // (the floor guards the fs1/iterated lanes); refused by the P2 driver
  // (iterating with L = 0 IS the drained split, KTJ-divergent at soil coupling).
  enum FsLMode      { FSL_CLASSIC = 0, FSL_OEDOMETRIC = 1, FSL_MANUAL = 2,
                      FSL_ZERO = 3 };
  // SM_MARCH is the default (flag absent): the overlay marches the fluid at each
  // qualifying transient commit using dT. The -staticMode flag selects HOLD
  // (freeze p, forces still applied from committed p) or STEADY (re-solve
  // H·p = f_seep) for static stages, where the domain "time" is the load factor
  // and a dT-march would be silently wrong physics (ADR §8 ⟨A-3⟩). Detection is
  // the FLAG, never the integrator.
  enum StaticMode   { SM_MARCH = 0, SM_HOLD = 1, SM_STEADY = 2 };
  enum RemovalMode  { RM_KEEP = 0, RM_DRAIN = 1 };
  // FU_EXPLICIT (ADR-73 P3b §3b.2 + §12 P3b): lumped-S* forward p-step,
  // matrix-free transient march (no CG). Advances at LOAD-APPLICATION time on
  // the trial window Δu (the hook only syncs — see the top-doc; the hook-time
  // advance was refuted, one force lag = unconditionally unstable). -fsL / L
  // are inert; SM_STEADY keeps its CG solve.
  enum FluidUpdate  { FU_IMPLICIT = 0, FU_EXPLICIT = 1 };
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

  // ---- moduli / stage transport (ADR-73 P2 §2.4, parameter route) -----------
  // Reached via `parameter $p loadPattern $overlay E|nu|layerE $i|layerNu $i`.
  int    setParameter(const char** argv, int argc, Parameter& param);
  int    updateParameter(int parameterID, Information& info);
  int    activateParameter(int parameterID);

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

  // ---- ADR-73 P3 overlay-aware Δt_cr advisory seam (CriticalTimeStep.cpp) ----
  // Per-element UNDRAINED stiffness augmentation ΔK_e = Q_e·S_e⁻¹·Q_eᵀ, the
  // element-local undrained condensation that turns the drained per-element
  // pencil into the DISCRETE undrained pencil (plan §3.2 pins; E7.4: material
  // formula √(1+K_f/(n·M_oed)) is documentation-only, ~1.85× conservative).
  // S_e = this cell's dense storage block S* (= S + αstab·H̃) — the SAME block
  // assembleStorageAndL scatters (built per-cell on demand, never read back
  // from CSR). Kadd is resized to (nNodes·ndm)² — node-major, first-ndm-DOFs
  // layout (the cell's own Qe layout, the same assumption HRZ lumping makes).
  // Returns false for unowned or dead cells, and (with a loud one-time
  // advisory naming the cell) for a singular/ill-conditioned S_e — never
  // silently optimistic. Blocks are cached lazily per cell; the cache is
  // invalidated by rebuildModuliCaches (S* depends on moduli via α-stab).
  // State-independent: valid for both useTangent and initial-stiff pencils.
  bool   getUndrainedAugmentation(int eleTag, Matrix& Kadd) const;

  // ---- ADR-73 P3b dual-CFL advisory seam (CriticalTimeStep.cpp) --------------
  // Per-fluid-advance forward-Euler diffusion bound for the EXPLICIT fluid lane:
  // Δt_diff = minEdge² / (2·ndm·c_v,max) (plan §3b.3 pin; conservative —
  // min-edge with max-c_v; ~7e3× slack vs the solid pencil at realistic k̄, so
  // it advises rather than governs). The -subcycle window must keep
  // N·Δt ≤ Δt_diff. Returns −1.0 unless this overlay is FU_EXPLICIT with a
  // completed snapshot and finite positive minEdge_/maxCv_ (zero new state:
  // both already exist for -subcycle auto).
  double explicitFluidDiffusionDt() const {
    if (fluidUpdate_ != FU_EXPLICIT || !snapshotOk_ ||
        maxCv_ <= 0.0 || minEdge_ <= 0.0)
      return -1.0;
    return minEdge_ * minEdge_ / (2.0 * (double)ndm_ * maxCv_);
  }

  // ---- ADR-73 P2 driver seam (LadrunoStaggeredAnalyze, LadrunoStaggeredDriver) --
  // The iterated fixed-stress driver latches the overlay into external-stepping
  // mode: applyLoad then injects +Q·pTrial_ (the current iterate's forces, not
  // the committed ones) and the driver's converged Domain::commit only SYNCS the
  // already-advanced fluid (NO extra advanceTrial → no double-march). The driver
  // reads lastAdvanceRelChange() (the toy `dcon` residual) after each
  // advanceTrial and transports staged moduli through the parameter route.
  void   setExternalStepping(bool on);
  bool   externalStepping(void) const;
  void   setResidualPScale(double s);        // toy `dcon` denominator offset
  double lastAdvanceRelChange(void) const;   // residual of the last advanceTrial
  int    staticModeCode(void) const;         // driven-set filter (== SM_MARCH?)
  int    fluidUpdateCode(void) const;        // driven-set filter (== FU_IMPLICIT?)
  int    fsLModeCode(void) const;            // driver refusal (== FSL_ZERO? ADR-73 P3)
  // Early fs1 sync of a pending -subcycle window before the driver takes over
  // (catch-up, plan §2.2 step 0). Returns 1 if it fired, 0 if nothing pending,
  // <0 on advance failure. Keeps the subcycle counters private to the overlay.
  int    catchUpPendingWindow(void);

  // ---- ADR-73 P4 recorder-channel seam (const reads only; ADR §4.1 pins) ----
  // The overlay p-field recorder channels (recorder ladruno -overlay / Monitor
  // -overlay) read the committed field through these trivial const accessors.
  // Zero physics, nothing mutates, nothing serializes. getRegionNodeTags() is
  // THE canonical id order for every channel; getCommittedP() is aligned to it.
  const std::vector<int>&    getRegionNodeTags() const { return regionNodeTags_; }
  const std::vector<double>& getCommittedP()     const { return pCommitted_; }
  const std::vector<int>&    getDrainedNodeTags() const { return drainedNodeTags_; }
  // Committed p at one region node tag. Returns 0.0 (and sets no error) for a
  // non-region tag: callers MUST pre-validate against getRegionNodeTags() (the
  // Monitor parser fatals at construction on a non-region tag).
  double pAtNodeTag(int nodeTag) const;
  // The lazy geometry snapshot has completed. The recorder must NOT fire the
  // snapshot from registration/writeModel (timing pin 4.3) — it guards every
  // channel build with this instead.
  bool   snapshotReady() const { return snapshotOk_; }

 private:
  struct Cell;                               // defined below (used in decls)

  // ---- geometry snapshot ----------------------------------------------------
  int    buildSnapshot(void);                // 0 ok, <0 fatal (prints)
  int    classifyCell(Element* ele, int& shapeKind, int& nNodes) const;
  int    assembleCell(int ci);               // per-cell blocks via the kernel
  template <class SH> void fillCell(Cell& c);// templated per-shape geometry+H/Q/f
  // per-cell E/nu resolution (layer override or global) — the ONE resolution
  // logic shared by the snapshot and the parameter-route rebuild.
  void   resolveCellModuli(Cell& c) const;
  // assemble the moduli-dependent aS_ (S + αstab·H̃) and aL_ (L) blocks from the
  // retained per-GP snapshot data (Np/dNp/dv). Called by fillCell at snapshot
  // AND by rebuildModuliCaches — the ONE storage/L assembly path (no drift).
  void   assembleStorageAndL(Cell& c);
  // dense per-cell S* block (S + αstab·H̃) from the retained per-GP data —
  // hoisted out of assembleStorageAndL so getUndrainedAugmentation can build
  // the SAME block for ONE cell on demand (ADR-73 P3; never read back from
  // CSR). Requires c.E/c.nu already resolved (resolveCellModuli).
  void   buildCellStorageBlock(const Cell& c, std::vector<double>& Se) const;
  // parameter-route moduli transport: re-resolve E/nu, re-derive α_stab / L, and
  // re-scatter aS_/aL_; aH_/fseep_/Qe untouched (perm/α only). Lazy at next use.
  void   rebuildModuliCaches(void);
  // ADR-73 P3b: (re)build sLump_ = per-row sums of the assembled aS_ CSR rows
  // (= lumped physical S; H̃·1 = 0 exactly so rowsum(S*) = rowsum(S_phys) —
  // the stab matrix cannot enter the explicit march, plan §3b.1 pin). Called
  // at snapshot completion and by rebuildModuliCaches, FU_EXPLICIT only.
  // Returns 0 ok; <0 (loud fatal naming the node tag) on any row-sum ≤ 0 —
  // the toy's slump.min() > 0 assert.
  int    buildStorageLump(void);
  // ADR-73 P3b (§12 P3b sequencing fix): the explicit lane's load-application-
  // time advance. Called from applyLoad under (FU_EXPLICIT && SM_MARCH &&
  // !externalStepping_). Owns the window bookkeeping (dtW = time − tLastSync_),
  // the -subcycle auto resolution, the ⟨A-3⟩ march advisory, and the
  // advancedThisStep_ selector applyLoad's injection reads. Every exit sets
  // advancedThisStep_ explicitly (no stale selector after aborted steps).
  void   explicitAdvanceAtApplyLoad(double time);
  void   buildSparsity(void);                // CSR pattern from region connectivity
  int    csrPos(int i, int j) const;         // CSR index of (i,j) or -1
  void   scatterMat(const Cell& c, const double* blk, std::vector<double>& tgt);
  void   scatterVecN(const Cell& c, const double* v, std::vector<double>& tgt);
  void   rebuildHFseep(void);                // reassemble H/f_seep (drain-removal)
  void   applyPInit(void);                   // -pInit: zero/steady/hydro/list
  // region u: committed by default. trial=true reads the solid's TRIAL disp —
  // needed by the P2 driver's coupling term mid-step (pre-commit, the current
  // solve iterate); post-commit trial==committed so the P1 hook is unaffected.
  void   readRegionDisp(std::vector<double>& u, bool trial = false) const;
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
    double hEdge;                 // element size h (for -stab auto α; geometry)
    std::vector<int>    rIdx;      // local node -> region-node index
    std::vector<Node*>  nodes;     // resolved at snapshot
    std::vector<double> Np;        // nGP*nNp
    std::vector<double> dNp;       // nGP*nNp*ndm
    std::vector<double> dv;        // nGP
    std::vector<double> Qe;        // (nNu*ndm) x nNp, row-major, ld = nNp
    Cell() : eleTag(0), shapeKind(0), nNodes(0), nGP(0), alive(true),
             layerIndex(-1), poro(0.0), E(0.0), nu(0.0), hEdge(0.0)
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
  int    dbTag2_;                            // payload Vector's own channel tag
                                             // (FileDatastore same-size clobber
                                             //  guard — see OV_NDATA note)

  // CSR (region-node × region-node; single shared pattern for H, S*, L)
  std::vector<int>    rowptr_;
  std::vector<int>    colind_;
  std::vector<int>    diagPtr_;              // per-row index of the diagonal entry
  std::vector<double> aH_, aS_, aL_;         // H, S* (=S+αstab·H̃), L values
  std::vector<double> fseep_;                // seepage source vector (size N)
  std::vector<char>   isDrained_;            // size N
  // ADR-73 P3b explicit-lane storage lump (size N; FU_EXPLICIT only, else
  // empty). DERIVED state — never serialized; rebuilt at snapshot and by
  // rebuildModuliCaches (buildStorageLump). NOT touched by removal rescan
  // (the fluid measure never shrinks — P1 pin).
  std::vector<double> sLump_;

  // fluid unknowns (size N) — internal double storage (packed to Vector on wire)
  std::vector<double> pCommitted_;
  std::vector<double> dpCommitted_;          // previous committed increment (rate ref)
  std::vector<double> pTrial_;
  std::vector<double> uSnapshot_;            // region u at last fluid advance (N*ndm)
  bool   uSnapshotValid_;                    // false after DB restore: re-derive
                                             // LAZILY at first use (⟨A-10⟩ —
                                             // restore-time node state is not
                                             // ordering-safe, measured 1.E-ii)

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
  bool   marchNoticeShown_;      // per-instance ⟨A-3⟩ first-march advisory latch

  // ---- ADR-73 P3b explicit-lane march state (transient; NOT serialized; -----
  // recvSelf → invalid; tLastSync_ re-derived at first applyLoad — the
  // uSnapshotValid_ pattern). The explicit lane advances at LOAD-APPLICATION
  // time (top-doc / ADR-73 §12 P3b): tLastSync_ anchors the window
  // dtW = time − tLastSync_; advancedThisStep_ selects the +Q·pTrial_
  // injection and tells the hook to sync; removalSyncPending_ forces the next
  // applyLoad advance after a removal commit regardless of the window.
  double tLastSync_;
  bool   tLastSyncValid_;
  bool   advancedThisStep_;
  bool   removalSyncPending_;
  // set by recvSelf ONLY: the first applyLoad anchors at time − Domain::getDT()
  // (the restored committed time — saves happen at sync commits) instead of at
  // `time`, so a restored explicit march continues bit-exactly with no
  // one-step hold. Dedicated flag — uSnapshotValid_ is re-derived by applyLoad
  // BEFORE the helper runs and cannot key this.
  bool   anchorFromRestore_;

  // ---- ADR-73 P2 driver latch (transient; NOT serialized; recvSelf → defaults) --
  bool   externalStepping_;      // driver latch (applyLoad→pTrial_, commit=sync)
  double pScale_;                // residual denominator offset (toy `dcon`)
  double relChange_;             // relative change of the last advanceTrial
  bool   moduliDirty_;           // parameter route touched E/nu → rebuild lazily
  bool   driverSubcycleNoticeShown_;  // one-time "driver overrides -subcycle"

  // ---- ADR-73 P3 undrained-augmentation cache (transient; lazily sized; ------
  // invalidated by rebuildModuliCaches — S* depends on moduli via α-stab).
  // augState_: 0 = unset, 1 = cached block valid, 2 = S_e singular (warned).
  mutable std::vector<std::vector<double> > augBlocks_;   // per-cell dense ΔK_e
  mutable std::vector<char> augState_;
  mutable std::map<int,int> augIndex_;                    // eleTag → cell index
  mutable bool   augSingularWarned_;   // one-time singular-S_e advisory latch

  // shape-kind codes
  enum { SHK_NONE = 0, SHK_T3 = 1, SHK_Q4 = 2, SHK_H8 = 3 };

  // non-copyable (owns an ofstream + snapshot); use getCopy() for a config clone
  LadrunoPorousOverlay(const LadrunoPorousOverlay&);
  LadrunoPorousOverlay& operator=(const LadrunoPorousOverlay&);
};

#endif
