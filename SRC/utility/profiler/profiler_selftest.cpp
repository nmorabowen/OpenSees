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

// Ladruno Stack Profiler — Phase P1 (C++<->Python parity self-test).
//
// File: SRC/utility/profiler/profiler_selftest.cpp
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: A STANDALONE main() (buildable independently of OpenSees) that
//   synthesizes the SAME "caseA" run that Ladruno_tools/profiler_viewer/make_sample.py
//   builds in Python, then drives ProfilerHDF5Writer::writeRun to emit
//   `cpp_selftest.h5`. This is the C++<->Python parity fixture: the structure,
//   group names, attribute names, compound dtypes and (deterministic) numeric
//   values must match what profiler_schema.write_run() produces, so that
//   ProfilerResults reads cpp_selftest.h5 identically to profile_sample.h5/caseA.
//
//   It does NOT link the Profiler timing core (PerfClock.cpp / Profiler.cpp): it
//   builds the POD wire structs directly and exercises ONLY the writer + hdf5.
//   That keeps the fixture minimal and independent of thread/registry machinery.
//
//   Build (see the bottom of this file / the task return value for exact cmds):
//     g++ -std=c++17 profiler_selftest.cpp ProfilerHDF5Writer.cpp \
//         -I<hdf5_inc> -L<hdf5_lib> -lhdf5 -o profiler_selftest
//   then:  ./profiler_selftest   (writes cpp_selftest.h5 in the cwd)
//   verify: python -c "from profiler_results import ProfilerResults; \
//                      pr=ProfilerResults('cpp_selftest.h5'); print(pr.runs())"

#include "ProfilerHDF5Writer.h"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace ops_profiler;

// ns per ms (mirror make_sample.py MS = 1.0e6).
static const double MS = 1.0e6;

// -------------------------------------------------------------------------
// node()  — mirror make_sample.py:25-35.
//   wall_ns     = int(wall_ms * MS)
//   wall_ns_min = int(wall_ns / max(calls,1) * 0.7)
//   wall_ns_max = int(wall_ns / max(calls,1) * 1.6)
//   cpu_ns      = int(wall_ns * cpu_frac)        (cpu_frac default 0.99)
// Integer truncation (C++ static_cast<int64_t>) matches Python int(): both
// truncate toward zero for the non-negative values used here.
// -------------------------------------------------------------------------
static ProfileNode mkNode(const std::string& name, double wall_ms, int64_t calls,
                          std::vector<ProfileNode> children = {},
                          std::vector<ElemByType> elem_by_type = {},
                          double cpu_frac = 0.99)
{
    ProfileNode n;
    n.name  = name;
    n.calls = calls;
    const int64_t wall_ns = static_cast<int64_t>(wall_ms * MS);
    const int64_t denom   = calls > 1 ? calls : 1;          // max(calls, 1)
    n.wall_ns     = wall_ns;
    n.wall_ns_min = static_cast<int64_t>(static_cast<double>(wall_ns) / denom * 0.7);
    n.wall_ns_max = static_cast<int64_t>(static_cast<double>(wall_ns) / denom * 1.6);
    n.cpu_ns      = static_cast<int64_t>(static_cast<double>(wall_ns) * cpu_frac);
    n.children     = std::move(children);
    n.elem_by_type = std::move(elem_by_type);
    return n;
}

// -------------------------------------------------------------------------
// elem_rows()  — mirror make_sample.py:38-44.
//   rows: [(classTag, count, wall_ms, fb_coupled)]
//   wns = int(wall_ms * MS); wall_ns_per_elem = wns / max(count,1)  (float)
// -------------------------------------------------------------------------
struct ElemSpec { int32_t classTag; int64_t count; double wall_ms; uint8_t fb; };

static std::vector<ElemByType> elemRows(const std::vector<ElemSpec>& rows)
{
    std::vector<ElemByType> out;
    out.reserve(rows.size());
    for (const ElemSpec& r : rows) {
        ElemByType e;
        e.classTag = r.classTag;
        e.count    = r.count;
        e.wall_ns  = static_cast<int64_t>(r.wall_ms * MS);
        const int64_t denom = r.count > 1 ? r.count : 1;     // max(count, 1)
        e.wall_ns_per_elem  = static_cast<double>(e.wall_ns) / static_cast<double>(denom);
        e.fb_coupled        = r.fb;
        out.push_back(e);
    }
    return out;
}

// tangent_elems(scale) — make_sample.py:48-53.
static std::vector<ElemByType> tangentElems(double scale)
{
    return elemRows({
        {74,   200, scale * 0.80, 1},   // forceBeamColumn (fb_coupled)
        {53,   400, scale * 0.15, 0},   // ShellMITC4
        {12, 20000, scale * 0.05, 0},   // truss
    });
}

// residual_elems(scale) — make_sample.py:56-61.
static std::vector<ElemByType> residualElems(double scale)
{
    return elemRows({
        {74,   200, scale * 0.78, 1},
        {53,   400, scale * 0.16, 0},
        {12, 20000, scale * 0.06, 0},
    });
}

// implicit_rollup(form_tangent_ms, factor_ms) — make_sample.py:64-79.
// caseA uses form_tangent_ms=3950.0, factor_ms=2250.0 (make_sample.py:128).
static ProfileNode implicitRollup(double form_tangent_ms, double factor_ms)
{
    return mkNode("step", form_tangent_ms + factor_ms + 2326.0, 1200, {
        mkNode("newStep", 50.0, 1200),
        mkNode("solveCurrentStep", form_tangent_ms + factor_ms + 2006.0, 1200, {
            mkNode("formTangent", form_tangent_ms, 3600, {}, tangentElems(form_tangent_ms)),
            mkNode("linearSolve", factor_ms + 350.0, 3600, {
                mkNode("solve.factor",   factor_ms, 3600),    // P1#7 split
                mkNode("solve.triSolve", 350.0,     3600),
            }),
            mkNode("update", 1000.0, 3600),
            mkNode("formUnbalance", 1200.0, 7200, {}, residualElems(1200.0)),
            mkNode("convTest", 60.0, 7200),
        }),
        mkNode("commit", 270.0, 1200),
    });
}

// -------------------------------------------------------------------------
// A small, DETERMINISTIC per-step Series (caseA's phase order, 4 steps).
// make_sample.py's series uses random jitter; for a C++<->Python *structural*
// parity fixture we use fixed values so the file is reproducible. nPhase=5 in
// caseA's phase order: formTangent, linearSolve, update, formUnbalance, convTest.
// -------------------------------------------------------------------------
static Series mkSeries()
{
    Series s;
    const int n = 4;
    s.phases = {"formTangent", "linearSolve", "update", "formUnbalance", "convTest"};
    for (int i = 0; i < n; ++i) {
        s.step.push_back(static_cast<int64_t>(i + 1));
        s.t.push_back((i + 1) * 0.005);
        s.dt.push_back(0.005);
        s.iters.push_back(3);
        s.mem_live_bytes.push_back(10000000 + i * 1000);
        s.mem_peak_bytes.push_back(12582912);
    }
    // wall_ms[n][nPhase] row-major, caseA phase means (make_sample.py:130-131).
    const float means[5] = {3.29f, 2.17f, 0.83f, 1.0f, 0.05f};
    for (int i = 0; i < n; ++i)
        for (int p = 0; p < 5; ++p)
            s.wall_ms.push_back(means[p]);
    return s;
}

// A series with steps but NO phase columns (nPhase()==0): the default P1 state
// until P2 wires per-phase wall times. wall_ms is logically [n,0]; the writer must
// fall back to CONTIGUOUS storage (a chunk of [n,1] over a 0-extent dim is illegal
// in HDF5). This is the C1 regression fixture (correctness review).
static Series mkSeriesEmptyPhases()
{
    Series s;
    const int n = 4;
    for (int i = 0; i < n; ++i) {
        s.step.push_back(static_cast<int64_t>(i + 1));
        s.t.push_back((i + 1) * 0.005);
        s.dt.push_back(0.005);
        s.iters.push_back(3);
        s.mem_live_bytes.push_back(10000000 + i * 1000);
        s.mem_peak_bytes.push_back(12582912);
    }
    // s.phases and s.wall_ms left empty -> nPhase()==0.
    return s;
}

// caseA memory snapshot — make_sample.py:132-133.
static MemorySnapshot mkMemory()
{
    MemorySnapshot m;
    m.matrix_live = 1422;
    m.vector_live = 880;
    m.id_live     = 300;
    m.peak_bytes  = 12582912;
    m.components_live = {
        {74,   200},
        {53,   400},
        {12, 20000},
        { 1,  8400},
    };
    return m;
}

// caseA run header — make_sample.py base_meta(...) + caseA overrides (124-127).
static RunMeta mkMeta()
{
    RunMeta m;
    m.model      = "frame_5story";
    m.engine_sha = "deadbeef";
    m.integrator = "Newmark";
    m.algorithm  = "Newton";
    m.solver     = "Mumps";
    m.units      = "N,mm,s";
    m.timestamp  = "2026-05-30T12:00:00";
    m.nSteps  = 1200;
    m.nDOF    = 50000;
    m.nElem   = 20600;
    m.nNode   = 8400;
    m.nnz     = 2400000;
    m.threads = 8;
    m.dt_min = 0.005;
    m.dt_max = 0.005;
    m.dt_cr  = 0.0;
    m.oversample_ratio = 0.0;
    m.wall_ms_total = 8421.0;
    m.cpu_ms_total  = 0.0;
    return m;
}

int main(int argc, char** argv)
{
    const std::string path = (argc > 1) ? argv[1] : "cpp_selftest.h5";

    ProfileNode    rollup = implicitRollup(3950.0, 2250.0);   // caseA
    RunMeta        meta   = mkMeta();
    Series         series = mkSeries();
    MemorySnapshot memory = mkMemory();

    ProfilerHDF5Writer w;
    if (!w.open(path)) {
        std::fprintf(stderr, "[selftest] FAILED to open '%s'\n", path.c_str());
        return 1;
    }

    if (!w.writeRun("caseA", rollup, meta, &series, memory)) {
        std::fprintf(stderr, "[selftest] writeRun('caseA') FAILED\n");
        w.close();
        return 2;
    }

    // C1 regression: a run whose series has steps but NO phases (the default P1
    // state until P2 wires phase columns). wall_ms is [n,0] -> must use contiguous
    // storage, not a chunk of [n,1] (which HDF5 rejects). Must succeed.
    Series emptyPhases = mkSeriesEmptyPhases();
    if (!w.writeRun("caseP1empty", rollup, meta, &emptyPhases, memory)) {
        std::fprintf(stderr, "[selftest] writeRun('caseP1empty') FAILED (C1 regression)\n");
        w.close();
        return 4;
    }

    // Immutability check: a second writeRun with the same id MUST fail (false),
    // and MUST NOT exit() the process (heeds the MPCO exit(-1) incident).
    if (w.writeRun("caseA", rollup, meta, &series, memory)) {
        std::fprintf(stderr, "[selftest] expected duplicate-run rejection, got success\n");
        w.close();
        return 3;
    }

    w.close();
    std::printf("[selftest] OK — wrote %s (run 'caseA'); duplicate correctly rejected\n",
                path.c_str());
    return 0;
}
