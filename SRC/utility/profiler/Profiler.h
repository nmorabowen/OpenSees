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

// Ladruno Stack Profiler — Phase P1 (timing core).
//
// File: SRC/utility/profiler/Profiler.h
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: Profiler timing core. Owns the per-thread profile trees, the
//   process-global singleton, the RAII ScopedTimer, and the cold-path report
//   assembly that produces the plain structs (RunMeta / Series / MemorySnapshot
//   / ProfileNode) consumed by the HDF5 writer.
//
//   HDF5-FREE (locked, 06_profiler.md gating model): this header and Profiler.cpp
//   never include <hdf5.h>. They include ProfilerHDF5Writer.h ONLY for the plain
//   result structs (that header is itself HDF5-free). When the writer TU is
//   compiled out (HDF5 absent), coarse phase timing still works; report() then
//   becomes a no-op + warning.
//
//   TIME UNITS: uint64_t NANOSECONDS everywhere on the hot path, sampled via
//   PerfClock. Converted to the schema's i8 wall_ns/cpu_ns (and ms totals) only
//   at report assembly.
//
//   THREAD MODEL (locked, P1#9): each thread owns its own tree + thread-local
//   cursor; there are NO shared atomics on the node counters (avoids Heisenberg
//   inflation). Per-thread roots are folded into one merged tree at report time
//   (single-threaded, off the hot path).
//
//   STRUCT RECONCILIATION: RunMeta / Series / MemorySnapshot / ComponentLive /
//   ElemByType / ProfileNode are the wire structs defined in ProfilerHDF5Writer.h.
//   The LIVE mutable tree here is ProfileNodeLive (a distinct type) with an
//   internal ElemTypeRow accumulator; mergedRollup() lowers the merged live tree
//   into the POD ProfileNode the writer consumes.

#ifndef Profiler_h
#define Profiler_h

#include "PerfClock.h"
#include "ProfilerHDF5Writer.h"   // HDF5-free: RunMeta, Series, MemorySnapshot, ProfileNode, ElemByType

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <atomic>
#include <optional>

namespace ops_profiler {

// ===========================================================================
// Live tree
// ===========================================================================

// One row of the optional per-element-class breakdown (deep mode, P0#1), held
// live during accumulation. wall_ns_per_elem is DERIVED at lowering time
// (wall_ns / count) and is therefore NOT stored here.
struct ElemTypeRow {
    int64_t  count      = 0;
    uint64_t wall_ns    = 0;
    uint8_t  fb_coupled = 0;   // 1 = force-based (coupled tangent/residual), 0 = displ-based
};

// Live mutable node. Children are stored in a vector (ownership + DETERMINISTIC
// insertion order = the diff key) plus an index map keyed by string_view into
// each child's own stable name_ for O(1) find-or-create on the hot path.
// elemByType_ uses std::map so rows emit in ascending classTag order
// (deterministic dataset ordering for diff).
class ProfileNodeLive {
public:
    explicit ProfileNodeLive(const char* nm, ProfileNodeLive* parent = nullptr);

    // find-or-create a child by name; O(1) amortized; deterministic order.
    ProfileNodeLive* child(const char* nm);

    // accumulate one completed scope sample (called from ScopedTimer dtor).
    void accumulate(uint64_t wall, uint64_t cpu) noexcept;

    // deep: add to a per-classTag bucket (lazy-allocates elemByType_).
    void addElem(int classTag, uint64_t wall, uint8_t fbCoupled);

    // merge another (per-thread) tree of identical shape into this one
    // (recurses by name/path; sums counters, min/max-folds, unions elem rows).
    void mergeFrom(const ProfileNodeLive& other);

    ProfileNodeLive* parent() const noexcept { return parent_; }

    // ---- read accessors used by report lowering ----
    const std::string& name()      const noexcept { return name_; }
    uint64_t calls()               const noexcept { return calls_; }
    uint64_t wallNsTotal()         const noexcept { return wallTotal_; }
    uint64_t wallNsMin()           const noexcept { return calls_ ? wallMin_ : 0; }
    uint64_t wallNsMax()           const noexcept { return wallMax_; }
    uint64_t cpuNsTotal()          const noexcept { return cpuTotal_; }
    const std::vector<std::unique_ptr<ProfileNodeLive>>& children() const noexcept { return children_; }
    const std::map<int, ElemTypeRow>* elemByType() const noexcept { return elemByType_.get(); }

private:
    std::string      name_;             // owned, stable; index keys are views into this
    ProfileNodeLive* parent_ = nullptr;

    uint64_t calls_     = 0;
    uint64_t wallTotal_ = 0;
    uint64_t wallMin_   = UINT64_MAX;   // sentinel; first sample sets it
    uint64_t wallMax_   = 0;
    uint64_t cpuTotal_  = 0;

    std::vector<std::unique_ptr<ProfileNodeLive>>      children_;   // ownership + DETERMINISTIC order
    std::unordered_map<std::string_view, ProfileNodeLive*> index_;  // fast find-or-create

    std::unique_ptr<std::map<int, ElemTypeRow>>        elemByType_; // lazy; deep mode only
};

// ===========================================================================
// Per-thread state
// ===========================================================================

// Numeric-buffer allocation classes the P4 counters attribute bytes to. The
// values index ThreadState::liveBytesByType[] and lower into the snapshot's
// matrix_live / vector_live / id_live fields. Any out-of-range value (an
// uncategorised alloc) is folded into the aggregate liveBytes only.
enum AllocType { ALLOC_MATRIX = 0, ALLOC_VECTOR = 1, ALLOC_ID = 2, NUM_ALLOC_TYPES = 3 };

// One per worker thread. The cursor `current` is the only per-scope mutable
// state and it is thread-local, so push/pop never touches another thread's
// memory. Memory counters (P4) live here too — summed at report time, no atomics.
struct ThreadState {
    std::unique_ptr<ProfileNodeLive> root;              // this thread's tree
    ProfileNodeLive*                 current = nullptr; // nesting cursor (== root.get() at rest)

    // P4 thread-local memory counters (no global atomic; P1#9 fix).
    int64_t liveBytes = 0;          // aggregate live across all types
    int64_t peakBytes = 0;          // high-water mark of aggregate liveBytes
    int64_t cumAllocs = 0;          // cumulative alloc events (never decremented)
    int64_t liveCount = 0;          // live alloc count (alloc++ / free--)
    int64_t liveBytesByType[NUM_ALLOC_TYPES] = {0, 0, 0};  // per-type live split

    // P4 TaggedObject census: live MovableObject count keyed by classTag
    // (born++ in the MovableObject ctor, died-- in the dtor). Thread-local like
    // the byte counters; merged at report time. classTag is a plain member of
    // MovableObject so it is valid through the dtor (TaggedObject::getClassTag is
    // virtual and would be unsafe in a base ctor/dtor).
    std::map<int, int64_t> componentLive;
};

// ===========================================================================
// Singleton config / run clock
// ===========================================================================

struct ProfilerConfig {
    bool perStep     = false;   // also keep the per-step Series buffer
    int  warmupSteps = 0;       // drop first K steps from the rollup (P1#8)
};

// Run-level wall/cpu/timestamp markers, set at start()/stop(); buildMeta()
// converts ns -> ms for wall_ms_total / cpu_ms_total.
struct RunClock {
    int64_t     startWallNs = 0, endWallNs = 0;
    int64_t     startCpuNs  = 0, endCpuNs  = 0;
    std::string timestamp;      // wall-clock ISO-8601 at start()
};

// ===========================================================================
// ScopedTimer — RAII hot-path timer (~50-100 ns)
// ===========================================================================
//
// ctor: read wall(+cpu), find-or-create child under the thread-local cursor,
// push. dtor: read wall(+cpu), accumulate the delta, pop. No allocation after
// first touch, no lock, no atomic. The macro constructs this only when the
// runtime gate is enabled (see ScopedTimerOpt + ProfilerMacros.h), so the timer
// itself does not re-check the flag.
class ScopedTimer {
public:
    explicit ScopedTimer(const char* name) noexcept;
    ~ScopedTimer() noexcept;

    // deep: tag a per-element-class sub-bucket on the *current* node (P0#1).
    void addElem(int classTag, uint64_t wall, uint8_t fbCoupled) noexcept;

    // Non-copyable and non-movable: exactly one dtor must pair with one ctor, or
    // the cursor-pop invariant (ts_.current = node_->parent()) corrupts. A user-
    // declared copy ctor already suppresses implicit moves; spell it out so a
    // future maintainer can't silently re-enable a move (correctness review M2).
    ScopedTimer(const ScopedTimer&)            = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&)                 = delete;
    ScopedTimer& operator=(ScopedTimer&&)      = delete;

private:
    ThreadState&     ts_;
    ProfileNodeLive* node_   = nullptr;
    uint64_t         startW_ = 0;
    uint64_t         startC_ = 0;
};

// Gated wrapper used by OPS_PROFILE_SCOPE: constructs the real timer only on a
// non-null name (zero cost when the runtime gate is off — one predicted branch,
// no clock reads, no node lookup). std::optional gives small-buffer storage,
// no heap.
class ScopedTimerOpt {
public:
    explicit ScopedTimerOpt(const char* name) noexcept {
        if (name) t_.emplace(name);     // gated construct
    }
    // Forward per-element attribution to the underlying timer when engaged; a
    // no-op when the gate was off (so OPS_PROFILE_ELEM is always safe to call).
    void addElem(int classTag, uint64_t wall, uint8_t fbCoupled) noexcept {
        if (t_) t_->addElem(classTag, wall, fbCoupled);
    }
    // True when the deep gate was on at construction (the timer is live). Callers
    // use this to short-circuit the per-element classTag lookup so an unprofiled
    // run never touches FE_Element::getElement() (P3 hot-path guard).
    bool engaged() const noexcept { return t_.has_value(); }
private:
    std::optional<ScopedTimer> t_;
};

// Best-effort force-based (flexibility-formulation) element classification for
// the elem_by_type fb_coupled tag (P0#1 / P3). Force-based beam-columns form
// their tangent from an internal coupled state-determination, so their
// tangent-vs-residual split is an artifact the viewer flags. This is a classTag
// whitelist of the standard force-based families (values mirror SRC/classTags.h);
// custom/unknown elements default to displacement-based (0), which is correct for
// the overwhelming majority. Kept here as a heuristic rather than a virtual on
// Element to avoid a broad vanilla footprint.
inline uint8_t isForceBasedClassTag(int classTag) noexcept {
    switch (classTag) {
        case 28:  case 29:                       // NLBeamColumn2d/3d
        case 34:  case 35:                       // BeamWithHinges2d/3d
        case 73:  case 731: case 74:             // ForceBeamColumn2d/Warping2d/3d
        case 75:  case 751: case 76:             // ElasticForceBeamColumn2d/Warping2d/3d
        case 77:  case 78:                       // ForceBeamColumnCBDI2d/3d
        case 171: case 172:                      // ForceBeamColumn2d/3dThermal
        case 192: case 193:                      // GradientInelasticBeamColumn2d/3d
        case 30765: case 30766: case 30767:      // MixedBeamColumn3d/2d/Asym3d
            return 1;
        default:
            return 0;
    }
}

// RAII per-element timer (P3): times one element kernel call and folds the wall
// into a named deep scope's per-classTag bucket (elem_by_type). Entirely inert —
// NO clock read, NO bucket touch — when the deep scope is off (tmr not engaged) or
// the classTag is negative (a constraint FE_Element with no underlying Element).
// The fb_coupled tag is derived from the classTag. Non-copyable/movable: exactly
// one ctor must pair with one dtor.
class ElemScope {
public:
    ElemScope(ScopedTimerOpt& tmr, int classTag) noexcept
        : tmr_(tmr), classTag_(classTag), startW_(0), active_(false) {
        if (classTag >= 0 && tmr.engaged()) {
            startW_ = static_cast<uint64_t>(PerfClock::wall_ns());
            active_ = true;
        }
    }
    ~ElemScope() noexcept {
        if (!active_) return;
        const uint64_t endW = static_cast<uint64_t>(PerfClock::wall_ns());
        const uint64_t dW   = (endW >= startW_) ? (endW - startW_) : 0;
        tmr_.addElem(classTag_, dW, isForceBasedClassTag(classTag_));
    }
    ElemScope(const ElemScope&)            = delete;
    ElemScope& operator=(const ElemScope&) = delete;
    ElemScope(ElemScope&&)                 = delete;
    ElemScope& operator=(ElemScope&&)      = delete;
private:
    ScopedTimerOpt& tmr_;
    int      classTag_;
    uint64_t startW_;
    bool     active_;
};

// ===========================================================================
// Profiler singleton
// ===========================================================================
//
// Storage is owned by OpenSeesCommands (constructed in its ctor; registered via
// setInstance) so the report flush has a controlled call site rather than
// static-dtor order (heeds the MPCO exit(-1) incident). theProfiler() returns a
// lazily-constructed default instance when no owner has registered (e.g. a
// standalone core unit test), keeping the core testable without the interpreter.
class Profiler {
public:
    Profiler();
    ~Profiler();

    // ---- runtime gates (cheap relaxed loads; the SINGLE-BINARY model) ----
    // enabled() gates coarse phase timing; deep() additionally gates per-element
    // buckets; mem() gates the alloc counters. All three live in one shipped binary
    // and are toggled per run by the profiler command (no special build needed).
    // The deep/mem layers are extra runtime flags rather than a compile flag so the
    // user has full control without rebuilding; cost when off is one predicted
    // branch. (A build may still hard-remove everything with -DOPS_PROFILE_DISABLE.)
    bool enabled() const noexcept { return enabled_.load(std::memory_order_relaxed); }
    bool deep()    const noexcept { return deepEnabled_.load(std::memory_order_relaxed); }
    bool mem()     const noexcept { return memEnabled_.load(std::memory_order_relaxed); }
    void setEnabled(bool on) noexcept { enabled_.store(on, std::memory_order_relaxed); }
    void setDeep(bool on)    noexcept { deepEnabled_.store(on, std::memory_order_relaxed); }
    void setMem(bool on)     noexcept { memEnabled_.store(on, std::memory_order_relaxed); }

    // ---- command verbs (driven by OPS_profiler in P5) ----
    void start();   // setEnabled(true) + mark run start (wall/cpu/timestamp)
    void stop();    // setEnabled(false) + mark run end
    void reset();   // clear all per-thread trees + series (config preserved)
    ProfilerConfig& config() noexcept { return config_; }

    // ---- per-thread access (hot path) ----
    ThreadState& threadState();   // thread_local fast path; registers on first call

    // ---- report assembly (cold path, single-threaded) ----
    // Fold all per-thread live trees into one merged live root, lower it into the
    // POD wire ProfileNode, and return it (cached, owned by Profiler).
    const ProfileNode& mergedRollup();
    RunMeta            buildMeta() const;
    const Series&      series() const noexcept { return series_; }
    MemorySnapshot     buildMemorySnapshot() const;

    // ---- per-step hook (P2 calls this from the analysis driver) ----
    void recordStep(int64_t step, double t, double dt, int32_t iters);

    // ---- ownership wiring (called by OpenSeesCommands) ----
    static void      setInstance(Profiler* p) noexcept;
    static Profiler& instance();   // == theProfiler()

    Profiler(const Profiler&)            = delete;
    Profiler& operator=(const Profiler&) = delete;

private:
    // registers the calling thread's ThreadState (mutex, once per thread).
    ThreadState* registerThread();

    std::atomic<bool>            enabled_{false};
    std::atomic<bool>            deepEnabled_{false};   // per-element buckets (runtime)
    std::atomic<bool>            memEnabled_{false};    // alloc counters (runtime)
    ProfilerConfig               config_{};

    std::mutex                   registryMtx_;   // guards threads_/owned_ registration ONLY
    std::vector<std::unique_ptr<ThreadState>> owned_;  // OWNS the ThreadState objects
    std::vector<ThreadState*>    threads_;        // flat registry view (one per worker thread)

    std::unique_ptr<ProfileNodeLive> mergedLive_; // built at report time
    ProfileNode                      mergedPod_;   // lowered wire form (cached)

    Series                       series_;          // optional per-step buffer (perStep)
    RunClock                     runClock_;
    int64_t                      warmupRemaining_ = 0;
};

// Process-global accessor.
inline Profiler& theProfiler() { return Profiler::instance(); }

// ===========================================================================
// Memory-counter hooks (P4 wiring; interface only). Update the THREAD-LOCAL
// ThreadState counters — no global atomic (P1#9). Summed at report time like the
// trees. Always compiled (the OPS_PROFILE_COUNT_* macros runtime-gate the call on
// theProfiler().mem()); only -DOPS_PROFILE_DISABLE removes the macros entirely.
// ===========================================================================
#ifndef OPS_PROFILE_DISABLE
inline void countAlloc(int64_t bytes, int type) noexcept {
    ThreadState& ts = theProfiler().threadState();
    ts.liveBytes += bytes;
    ts.cumAllocs += 1;
    ++ts.liveCount;
    if (type >= 0 && type < NUM_ALLOC_TYPES) ts.liveBytesByType[type] += bytes;
    if (ts.liveBytes > ts.peakBytes) ts.peakBytes = ts.liveBytes;
}
inline void countFree(int64_t bytes, int type) noexcept {
    ThreadState& ts = theProfiler().threadState();
    ts.liveBytes -= bytes;
    --ts.liveCount;
    if (type >= 0 && type < NUM_ALLOC_TYPES) ts.liveBytesByType[type] -= bytes;
}

// TaggedObject census (P4): bump/decrement the live count for a classTag. Called
// from the MovableObject ctor/dtor via OPS_PROFILE_CENSUS_BORN/DIED (which apply
// the mem() gate). One std::map op on a profiling run; off-path cost is the gate.
inline void countComponentBorn(int classTag) noexcept {
    ThreadState& ts = theProfiler().threadState();
    ++ts.componentLive[classTag];
}
inline void countComponentDied(int classTag) noexcept {
    ThreadState& ts = theProfiler().threadState();
    --ts.componentLive[classTag];
}
#endif // !OPS_PROFILE_DISABLE

} // namespace ops_profiler

#endif // Profiler_h
