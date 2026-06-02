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

// Ladruno Stack Profiler — Phase P1 (HDF5 writer).
//
// File: SRC/utility/profiler/ProfilerHDF5Writer.h
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: ProfilerHDF5Writer — translates the profiler's plain in-memory
//   result structs into the on-disk HDF5 contract defined (authoritatively) by
//   Ladruno_tools/profiler_viewer/profiler_schema.py. The C++ writer MUST be
//   byte-structurally identical to that reference writer: same group names,
//   nesting, attribute names, compound dtypes, and NANOSECOND integer time
//   units.
//
//   HDF5 ISOLATION (locked, 06_profiler.md gating model):
//     This header is deliberately HDF5-FREE. It declares only plain structs and
//     a thin class whose single private handle is stored as `long long` (hid_t
//     is int64_t on every supported build). ONLY ProfilerHDF5Writer.cpp pulls in
//     <hdf5.h>. The Profiler timing core therefore links and runs (coarse phase
//     timing) even when the writer translation unit is compiled out because HDF5
//     is unavailable.
//
//   STRUCT OWNERSHIP (reconciled):
//     The wire structs below (ElemByType, ComponentLive, ProfileNode, RunMeta,
//     Series, MemorySnapshot) are the SINGLE definition consumed by the writer.
//     The Profiler core (Profiler.h) includes THIS header to reuse RunMeta /
//     Series / MemorySnapshot, and builds a `ProfileNode` POD tree at report
//     time from its internal live tree. Field ORDER and NAMES in the compound
//     structs are load-bearing (they back H5Tinsert/HOFFSET) — do not reorder.

#ifndef ProfilerHDF5Writer_h
#define ProfilerHDF5Writer_h

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>

namespace ops_profiler {

// -------------------------------------------------------------------------
// Compound-dataset row structs. The field order below MUST match the numpy
// dtypes in profiler_schema.py so HOFFSET-based H5Tinsert produces a binary
// layout numpy reads back by field name.
// -------------------------------------------------------------------------

// Mirror of ELEM_BY_TYPE_DT (profiler_schema.py:50-56). One row of the optional
// per-element-class breakdown attached to deep-mode nodes (formTangent / ...).
struct ElemByType {                 // schema field   numpy dtype
    int32_t  classTag;              //  "classTag"          i4
    int64_t  count;                 //  "count"             i8
    int64_t  wall_ns;               //  "wall_ns"           i8
    double   wall_ns_per_elem;      //  "wall_ns_per_elem"  f8
    uint8_t  fb_coupled;            //  "fb_coupled"        u1  (1=force-based)
};

// Mirror of COMPONENTS_LIVE_DT (profiler_schema.py:58-61). TaggedObject census.
struct ComponentLive {
    int32_t  classTag;              //  "classTag"  i4
    int64_t  count;                 //  "count"     i8
};

// -------------------------------------------------------------------------
// Rollup-tree node (POD wire form). One node => one HDF5 group named `name`,
// carrying NODE_ATTRS_INT (profiler_schema.py:68) and an optional
// `elem_by_type` child dataset. Children recurse as subgroups in INSERTION
// ORDER — the stable diff key (P0#6); the core must emit them in a fixed
// taxonomy order so two runs align by path.
// -------------------------------------------------------------------------
struct ProfileNode {
    std::string name;                          // => HDF5 group name (schema:119)
    int64_t calls       = 0;                   //  "calls"
    int64_t wall_ns     = 0;                   //  "wall_ns"
    int64_t wall_ns_min = 0;                   //  "wall_ns_min"
    int64_t wall_ns_max = 0;                   //  "wall_ns_max"
    int64_t cpu_ns      = 0;                   //  "cpu_ns"
    std::vector<ElemByType>  elem_by_type;     // optional child dataset (schema:122-124)
    std::vector<ProfileNode> children;         // deterministic order = insertion order
};

// -------------------------------------------------------------------------
// Run header. The three groups (STR / INT / FLOAT) match RUN_ATTRS_STR /
// RUN_ATTRS_INT / RUN_ATTRS_FLOAT (profiler_schema.py:64-66) exactly, by name.
// -------------------------------------------------------------------------
struct RunMeta {
    // RUN_ATTRS_STR
    std::string model, engine_sha, integrator, algorithm, solver, units, timestamp;
    // RUN_ATTRS_INT  (written as i8 attrs; reader does not pin width)
    int64_t nSteps = 0, nDOF = 0, nElem = 0, nNode = 0, nnz = 0, threads = 0;
    // RUN_ATTRS_FLOAT  (P0#5 normalizers; ns->ms conversions done by the core)
    double  dt_min = 0, dt_max = 0, dt_cr = 0, oversample_ratio = 0;
    double  wall_ms_total = 0, cpu_ms_total = 0;
};

// -------------------------------------------------------------------------
// Per-step time history. If the run has no per-step buffer (perStep off) the
// core passes nullptr OR an empty Series and NO /series group is written
// (profiler_schema.py:99 `if series:`). All scalar vectors are length nSteps;
// wall_ms is row-major [nSteps * nPhase] float32, labeled by the `phases` attr.
// -------------------------------------------------------------------------
struct Series {
    // SERIES_SCALAR (profiler_schema.py:70-73) — all 1-D, length nSteps()
    std::vector<int64_t> step;            // "step"            i8
    std::vector<double>  t;               // "t"               f8
    std::vector<double>  dt;              // "dt"              f8
    std::vector<int32_t> iters;           // "iters"           i4  (P0#3)
    std::vector<int64_t> mem_live_bytes;  // "mem_live_bytes"  i8
    std::vector<int64_t> mem_peak_bytes;  // "mem_peak_bytes"  i8
    // wall_ms[nSteps][nPhase] row-major float32 (profiler_schema.py:105-107)
    std::vector<float>       wall_ms;     // size == nSteps()*nPhase()
    std::vector<std::string> phases;      // attr "phases" on wall_ms ds (schema:108)

    std::size_t nSteps() const { return step.size(); }
    std::size_t nPhase() const { return phases.size(); }
};

// -------------------------------------------------------------------------
// End-of-run memory snapshot (profiler_schema.py:110-115). The /memory group
// and its components_live dataset are ALWAYS written (the dataset may be
// zero-length when deep/memory mode is off).
// -------------------------------------------------------------------------
struct MemorySnapshot {
    int64_t matrix_live = 0, vector_live = 0, id_live = 0, peak_bytes = 0;  // attrs (schema:112)
    std::vector<ComponentLive> components_live;                            // dataset (schema:114-115)
};

// -------------------------------------------------------------------------
// Public writer surface. HDF5-free signature: callers never see hid_t.
// -------------------------------------------------------------------------
class ProfilerHDF5Writer {
public:
    ProfilerHDF5Writer() = default;
    ~ProfilerHDF5Writer();   // calls close() (final H5Fflush) if still open

    ProfilerHDF5Writer(const ProfilerHDF5Writer&)            = delete;
    ProfilerHDF5Writer& operator=(const ProfilerHDF5Writer&) = delete;

    // Open-or-create the file: H5Fopen(H5F_ACC_RDWR) if the path already exists
    // on disk, else H5Fcreate(H5F_ACC_EXCL). Sets schema_version="0.1.0" via
    // setdefault (only if absent). Returns false on failure — the caller keeps
    // its in-memory result tree (v1 = write-on-report) and never exit()s.
    bool open(const std::string& path);

    // Final H5Fflush + H5Fclose. Safe to call twice / when not open.
    void close();

    // Append exactly one run as /runs/<run_id>. FAILS (returns false) if that
    // group already exists -> run immutability (profiler_schema.py:86-87).
    // Writes run-header attrs, the rollup tree, optional /series, and /memory,
    // then a single H5Fflush(H5F_SCOPE_GLOBAL). `series` may be nullptr or have
    // nSteps()==0 to skip the /series group. On any internal HDF5 error it
    // closes whatever it opened and returns false (NEVER exit()).
    bool writeRun(const std::string&    run_id,
                  const ProfileNode&    rollup,
                  const RunMeta&        meta,
                  const Series*         series,        // nullptr / empty => skip /series
                  const MemorySnapshot& memory);

    bool isOpen() const { return file_ >= 0; }

private:
    // hid_t kept opaque (long long == hid_t/int64_t on all builds) so <hdf5.h>
    // stays out of this header. The .cpp does `hid_t fid = (hid_t)file_;`.
    long long file_ = -1;
};

} // namespace ops_profiler

#endif // ProfilerHDF5Writer_h
