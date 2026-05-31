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

// Ladruno Stack Profiler — Phase P1 (timing core implementation).
//
// File: SRC/utility/profiler/Profiler.cpp
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: Implementation of the Profiler timing core declared in Profiler.h.
//   Owns the per-thread live profile trees, the process-global singleton, the
//   RAII ScopedTimer, the per-thread tree merge, and the cold-path report
//   assembly that lowers the merged live tree into the POD wire structs
//   (RunMeta / Series / MemorySnapshot / ProfileNode) the HDF5 writer consumes.
//
//   HDF5-FREE: this TU does NOT include <hdf5.h>. It includes Profiler.h (which
//   pulls in ProfilerHDF5Writer.h only for the plain result structs). Only
//   ProfilerHDF5Writer.cpp touches <hdf5.h>.
//
//   TIME UNITS: uint64_t NANOSECONDS on the hot path (PerfClock), converted to
//   the schema's i8 ns counters / ms totals only at report assembly.
//
//   THREAD MODEL (P1#9): each thread owns its own tree + a thread-local cursor;
//   node counters carry NO shared atomics. Per-thread roots are folded into one
//   merged tree at report time (single-threaded, off the hot path).

#include "Profiler.h"

#include <algorithm>   // std::min, std::max
#include <chrono>
#include <ctime>
#include <cstdio>

namespace ops_profiler {

// ===========================================================================
// ProfileNodeLive
// ===========================================================================
//
// children_ holds unique_ptr<ProfileNodeLive>, so each child's name_ string has
// a STABLE heap address independent of the children_ vector's own reallocations.
// index_ therefore keys on string_view(child->name_) safely: growing children_
// moves the unique_ptr objects (pointers), never the pointed-to node, so the
// string buffer the view references does not move. (Small-string-optimized
// std::string would move its buffer if the *node* moved — it does not, because
// only the owning pointer in the vector moves.)

ProfileNodeLive::ProfileNodeLive(const char* nm, ProfileNodeLive* parent)
    : name_(nm ? nm : ""), parent_(parent)
{
}

ProfileNodeLive*
ProfileNodeLive::child(const char* nm)
{
    // O(1) amortized find. Keyed by string_view over the candidate name; on a
    // hit the stored view points into the existing child's stable name_.
    std::string_view key(nm ? nm : "");
    auto it = index_.find(key);
    if (it != index_.end())
        return it->second;

    // Miss: create a new child, preserving INSERTION ORDER in children_ (the
    // diff key). Index it by a view into the child's own name_ (stable heap).
    children_.push_back(std::make_unique<ProfileNodeLive>(nm, this));
    ProfileNodeLive* raw = children_.back().get();
    index_.emplace(std::string_view(raw->name_), raw);
    return raw;
}

void
ProfileNodeLive::accumulate(uint64_t wall, uint64_t cpu) noexcept
{
    ++calls_;
    wallTotal_ += wall;
    if (wall < wallMin_) wallMin_ = wall;
    if (wall > wallMax_) wallMax_ = wall;
    cpuTotal_  += cpu;
}

void
ProfileNodeLive::addElem(int classTag, uint64_t wall, uint8_t fbCoupled)
{
    if (!elemByType_)
        elemByType_ = std::make_unique<std::map<int, ElemTypeRow>>();
    ElemTypeRow& r = (*elemByType_)[classTag];
    r.count      += 1;
    r.wall_ns    += wall;
    r.fb_coupled  = fbCoupled;   // a class is uniformly force- or displ-based
}

void
ProfileNodeLive::mergeFrom(const ProfileNodeLive& other)
{
    // Sum scalar counters; min/max-fold the extrema. Note: do NOT bump calls_
    // for the merge itself — we are summing two independent accumulations of the
    // same logical node.
    calls_     += other.calls_;
    wallTotal_ += other.wallTotal_;
    cpuTotal_  += other.cpuTotal_;
    if (other.calls_) {
        if (other.wallMin_ < wallMin_) wallMin_ = other.wallMin_;
        if (other.wallMax_ > wallMax_) wallMax_ = other.wallMax_;
    }

    // Merge per-classTag element rows.
    if (other.elemByType_) {
        if (!elemByType_)
            elemByType_ = std::make_unique<std::map<int, ElemTypeRow>>();
        for (const auto& kv : *other.elemByType_) {
            ElemTypeRow& dst = (*elemByType_)[kv.first];
            dst.count      += kv.second.count;
            dst.wall_ns    += kv.second.wall_ns;
            dst.fb_coupled  = kv.second.fb_coupled;
        }
    }

    // Recurse children by name/path; child() creates the matching node if this
    // tree did not see it (e.g. a phase only one thread executed).
    for (const auto& c : other.children_) {
        this->child(c->name_.c_str())->mergeFrom(*c);
    }
}

// ===========================================================================
// ScopedTimer
// ===========================================================================
//
// ctor: sample wall(+cpu), find-or-create the child under the thread-local
// cursor, descend the cursor into it. dtor: sample the delta, accumulate, and
// restore the cursor to the parent. No allocation after first touch, no lock,
// no atomic on the counters.

ScopedTimer::ScopedTimer(const char* name) noexcept
    : ts_(Profiler::instance().threadState())
{
    // child() does a first-touch allocation on a cache miss (push_back + map
    // insert) which can throw bad_alloc. The ctor is noexcept (the macro/optional
    // path relies on it), so degrade to an untimed scope on failure rather than
    // letting a throw escape -> std::terminate (correctness review M1). A null
    // node_ leaves the cursor untouched; the dtor and addElem null-check it.
    try {
        node_ = ts_.current->child(name);
    } catch (...) {
        node_ = nullptr;
        return;
    }
    ts_.current = node_;
    startW_ = static_cast<uint64_t>(PerfClock::wall_ns());
    startC_ = static_cast<uint64_t>(PerfClock::cpu_ns());
}

ScopedTimer::~ScopedTimer() noexcept
{
    if (node_ == nullptr)
        return;   // ctor allocation failed; nothing was pushed, nothing to pop
    const uint64_t endW = static_cast<uint64_t>(PerfClock::wall_ns());
    const uint64_t endC = static_cast<uint64_t>(PerfClock::cpu_ns());
    const uint64_t dW   = (endW >= startW_) ? (endW - startW_) : 0;
    const uint64_t dC   = (endC >= startC_) ? (endC - startC_) : 0;
    node_->accumulate(dW, dC);
    // Pop the cursor back to the parent. node_->parent() is the cursor's value
    // at ctor time (we only ever descend by one), so this restores it exactly.
    ts_.current = node_->parent();
}

void
ScopedTimer::addElem(int classTag, uint64_t wall, uint8_t fbCoupled) noexcept
{
    if (node_ != nullptr)
        node_->addElem(classTag, wall, fbCoupled);
}

// ===========================================================================
// Profiler — singleton storage
// ===========================================================================

namespace {
// Process-global registered instance (set by OpenSeesCommands via setInstance).
Profiler* g_instance = nullptr;
}

Profiler::Profiler()  = default;
Profiler::~Profiler() = default;

void
Profiler::setInstance(Profiler* p) noexcept
{
    g_instance = p;
}

Profiler&
Profiler::instance()
{
    if (g_instance)
        return *g_instance;
    // Lazy function-local fallback so the core is usable in a standalone unit
    // test without the interpreter wiring it up. Constructed on first use,
    // destroyed at process exit AFTER any thread that touched it (Meyers
    // singleton; we deliberately do NOT flush in a static dtor — heeds the
    // MPCO exit(-1) incident; flushing is the command layer's job).
    static Profiler fallback;
    return fallback;
}

// ---------------------------------------------------------------------------
// Per-thread registry
// ---------------------------------------------------------------------------

ThreadState&
Profiler::threadState()
{
    // One ThreadState per OS thread, cached thread-locally. The pointer is owned
    // by owned_ (under the registry mutex); threads_ is the flat view the report
    // path folds over.
    thread_local ThreadState* ts = nullptr;
    if (!ts)
        ts = registerThread();
    return *ts;
}

ThreadState*
Profiler::registerThread()
{
    std::lock_guard<std::mutex> lock(registryMtx_);
    auto owned = std::make_unique<ThreadState>();
    owned->root    = std::make_unique<ProfileNodeLive>("root");
    owned->current = owned->root.get();
    ThreadState* raw = owned.get();
    owned_.push_back(std::move(owned));   // owner (see header note: .cpp owns)
    threads_.push_back(raw);              // registry view used by the report path
    return raw;
}

// ---------------------------------------------------------------------------
// Command verbs
// ---------------------------------------------------------------------------

void
Profiler::start()
{
    runClock_.startWallNs = PerfClock::wall_ns();
    runClock_.startCpuNs  = PerfClock::cpu_ns();
    runClock_.endWallNs   = runClock_.startWallNs;
    runClock_.endCpuNs    = runClock_.startCpuNs;

    // ISO-8601 wall-clock timestamp (local time) at start().
    {
        std::time_t now = std::time(nullptr);
        std::tm     tmv{};
#if defined(_WIN32)
        localtime_s(&tmv, &now);
#else
        localtime_r(&now, &tmv);
#endif
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tmv);
        runClock_.timestamp = buf;
    }

    warmupRemaining_ = config_.warmupSteps;
    setEnabled(true);
}

void
Profiler::stop()
{
    setEnabled(false);
    runClock_.endWallNs = PerfClock::wall_ns();
    runClock_.endCpuNs  = PerfClock::cpu_ns();
}

void
Profiler::reset()
{
    // Clear all per-thread trees + the series buffer; preserve config_ and the
    // enabled flag. Each thread's cursor is reset to a fresh root so an in-flight
    // (but currently inactive) thread re-anchors cleanly on its next scope.
    {
        std::lock_guard<std::mutex> lock(registryMtx_);
        for (ThreadState* ts : threads_) {
            ts->root    = std::make_unique<ProfileNodeLive>("root");
            ts->current = ts->root.get();
            ts->liveBytes = 0;
            ts->peakBytes = 0;
            ts->cumAllocs = 0;
            ts->liveCount = 0;
            for (int t = 0; t < NUM_ALLOC_TYPES; ++t)
                ts->liveBytesByType[t] = 0;
            ts->componentLive.clear();
        }
    }
    mergedLive_.reset();
    mergedPod_ = ProfileNode{};

    series_ = Series{};
    runClock_ = RunClock{};
    warmupRemaining_ = config_.warmupSteps;
}

// ---------------------------------------------------------------------------
// Report assembly (cold path, single-threaded)
// ---------------------------------------------------------------------------

namespace {
// Recursively lower a merged live node into the POD wire ProfileNode the writer
// consumes. Children are emitted in children_ vector order (= insertion order =
// the stable diff key). elem_by_type rows emit in ascending classTag order
// (std::map iteration) with wall_ns_per_elem derived here (divide-by-zero
// guarded).
void lowerNode(const ProfileNodeLive& live, ProfileNode& out)
{
    out.name        = live.name();
    out.calls       = static_cast<int64_t>(live.calls());
    out.wall_ns     = static_cast<int64_t>(live.wallNsTotal());
    out.wall_ns_min = static_cast<int64_t>(live.wallNsMin());
    out.wall_ns_max = static_cast<int64_t>(live.wallNsMax());
    out.cpu_ns      = static_cast<int64_t>(live.cpuNsTotal());

    if (const std::map<int, ElemTypeRow>* rows = live.elemByType()) {
        out.elem_by_type.reserve(rows->size());
        for (const auto& kv : *rows) {                 // ascending classTag
            const ElemTypeRow& r = kv.second;
            ElemByType e;
            e.classTag         = kv.first;
            e.count            = r.count;
            e.wall_ns          = static_cast<int64_t>(r.wall_ns);
            e.wall_ns_per_elem = (r.count > 0)
                ? (static_cast<double>(r.wall_ns) / static_cast<double>(r.count))
                : 0.0;
            e.fb_coupled       = r.fb_coupled;
            out.elem_by_type.push_back(e);
        }
    }

    const auto& kids = live.children();
    out.children.reserve(kids.size());
    for (const auto& c : kids) {
        out.children.emplace_back();
        lowerNode(*c, out.children.back());
    }
}
} // namespace

const ProfileNode&
Profiler::mergedRollup()
{
    // PRECONDITION (correctness review H1): the caller must guarantee no worker
    // thread is currently inside a profiled scope. The per-thread live trees are
    // read here WITHOUT hot-path synchronization (that lock is deliberately absent
    // by design, P1#9), so a concurrent ScopedTimer mutating a tree (child() can
    // push_back/realloc children_) would race this read. registryMtx_ guards only
    // thread registration, NOT tree mutation. In practice report is called from the
    // single-threaded command layer after profiler('stop'); P5 must enforce that the
    // active analysis is idle before reporting (assert/guard added when seams land).
    if (enabled_) {
        // Reporting while still accumulating is a usage error under any parallel
        // solver; warn rather than risk a silent race. (Coarse single-thread runs
        // are safe, so this is a warning, not a hard stop.) fprintf keeps the core
        // dependency-free (no OPS_Globals/opserr) so it stays standalone-compilable.
        std::fprintf(stderr, "[profiler] mergedRollup() called while still enabled; "
                     "call profiler('stop') first to avoid racing live worker threads.\n");
    }

    // Fold all per-thread live trees into one merged live root...
    mergedLive_ = std::make_unique<ProfileNodeLive>("root");
    {
        std::lock_guard<std::mutex> lock(registryMtx_);
        for (ThreadState* ts : threads_) {
            if (ts->root)
                mergedLive_->mergeFrom(*ts->root);
        }
    }
    // ...then lower it into the POD wire form. The viewer reads the rollup root
    // node under /rollup/<root.name>; we keep the taxonomy root named "root".
    mergedPod_ = ProfileNode{};
    lowerNode(*mergedLive_, mergedPod_);
    return mergedPod_;
}

RunMeta
Profiler::buildMeta() const
{
    RunMeta m;

    // ns deltas -> ms totals (÷1e6). Guard against an unstamped clock.
    const int64_t wallNs = runClock_.endWallNs - runClock_.startWallNs;
    const int64_t cpuNs  = runClock_.endCpuNs  - runClock_.startCpuNs;
    m.wall_ms_total = (wallNs > 0) ? static_cast<double>(wallNs) / 1.0e6 : 0.0;
    m.cpu_ms_total  = (cpuNs  > 0) ? static_cast<double>(cpuNs)  / 1.0e6 : 0.0;

    m.timestamp = runClock_.timestamp;
    m.threads   = static_cast<int64_t>(threads_.size());

    // model / engine_sha / integrator / algorithm / solver / units and the
    // size normalizers (nSteps/nDOF/nElem/nNode/nnz, dt_*, oversample_ratio) are
    // populated by the P5 command layer from the live SimulationInformation /
    // SOE / census before handing meta to the writer. P1 leaves them defaulted.
    if (config_.perStep)
        m.nSteps = static_cast<int64_t>(series_.nSteps());

    return m;
}

MemorySnapshot
Profiler::buildMemorySnapshot() const
{
    // Always returns a valid struct (the writer always writes /memory). Sum the
    // per-thread live/peak byte counters (P4 fills these; zero in P1).
    MemorySnapshot snap;
    int64_t peak = 0;
    int64_t byType[NUM_ALLOC_TYPES] = {0, 0, 0};
    std::map<int, int64_t> census;   // merged per-classTag live counts (P4 census)
    {
        std::lock_guard<std::mutex> lock(
            const_cast<std::mutex&>(registryMtx_));   // read-only fold; mutex guards threads_
        for (const ThreadState* ts : threads_) {
            peak += ts->peakBytes;
            for (int t = 0; t < NUM_ALLOC_TYPES; ++t)
                byType[t] += ts->liveBytesByType[t];
            for (const auto& kv : ts->componentLive)
                census[kv.first] += kv.second;
        }
    }
    // Per-type live splits (P4 alloc counters in Matrix/Vector/ID); peak_bytes is
    // the high-water mark of the aggregate live bytes.
    snap.matrix_live = byType[ALLOC_MATRIX];
    snap.vector_live = byType[ALLOC_VECTOR];
    snap.id_live     = byType[ALLOC_ID];
    snap.peak_bytes  = peak;

    // TaggedObject census (P4): one ComponentLive row per classTag with a non-zero
    // live count. std::map iterates in ascending classTag order, giving a
    // deterministic dataset ordering. Zero (and net-negative) counts are skipped.
    for (const auto& kv : census) {
        if (kv.second > 0)
            snap.components_live.push_back(ComponentLive{
                static_cast<int32_t>(kv.first), kv.second});
    }
    return snap;
}

// ---------------------------------------------------------------------------
// Per-step hook
// ---------------------------------------------------------------------------

void
Profiler::recordStep(int64_t step, double t, double dt, int32_t iters)
{
    if (!config_.perStep)
        return;

    // Honor warmup: silently drop the first K recorded steps (P1#8) so the
    // first-touch outlier + turbo ramp don't pollute the series/rollup window.
    if (warmupRemaining_ > 0) {
        --warmupRemaining_;
        return;
    }

    series_.step.push_back(step);
    series_.t.push_back(t);
    series_.dt.push_back(dt);
    series_.iters.push_back(iters);

    // Memory snapshot per step (P4 fills the thread counters; zero in P1).
    int64_t live = 0;
    int64_t peak = 0;
    {
        std::lock_guard<std::mutex> lock(registryMtx_);
        for (const ThreadState* ts : threads_) {
            live += ts->liveBytes;
            peak += ts->peakBytes;
        }
    }
    series_.mem_live_bytes.push_back(live);
    series_.mem_peak_bytes.push_back(peak);

    // Per-phase wall_ms row: the phase column set is fixed (series_.phases) and
    // is established by the P2 driver wiring. In P1 phases is empty (nPhase()==0)
    // so we append a zero-width row — wall_ms stays consistent at
    // nSteps()*nPhase(). Once phases is populated, append one float per phase in
    // the fixed column order.
    const std::size_t nPhase = series_.phases.size();
    for (std::size_t p = 0; p < nPhase; ++p)
        series_.wall_ms.push_back(0.0f);
}

} // namespace ops_profiler
