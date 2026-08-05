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

// Ladruno ADR-75 P1i — run-header population shared by the two command layers.
//
// File: SRC/utility/profiler/ProfilerRunMeta.h
//
// Why a header and not a .cpp in the profiler lib:
//
//   * it needs Domain, and the profiler core (OPS_Utility) deliberately does not
//     depend on the domain layer;
//   * it needs mkl_get_max_threads() to be authoritative about the thread cap,
//     and -D_PARDISO is applied to the two COMMAND targets (OPS_InterpPyCmds and
//     the OpenSees exe), not to OPS_Utility -- so the MKL branch can only compile
//     where this header is included, which is exactly the right place.
//
// Why it exists at all: `profiler report` and `profiler checkpoint` are
// implemented FOUR times -- twice in SRC/interpreter/OpenSeesCommands.cpp
// (Python) and twice in SRC/tcl/commands.cpp (Tcl), which are completely separate
// command ladders. That split is a banked trap on this fork (wiring one does not
// wire the other; it left `system Pardiso` Python-only for a release). The run
// header had drifted the same way: nElem/nNode were promised by a comment in
// Profiler.cpp and populated by nobody. One helper, four call sites, no drift.

#ifndef ProfilerRunMeta_h
#define ProfilerRunMeta_h

#include <Domain.h>
#include "Profiler.h"          // RunMeta (via ProfilerHDF5Writer.h)

#ifdef _PARDISO
#include <mkl_service.h>       // mkl_get_max_threads
#endif

// Fill the model-size normalizers and finalize the thread count.
//
// `theDomain` may be null (no model built yet) -- the size fields then keep the
// 0 they already carry, which is the honest answer.
static inline void
ops_profiler_fillModelMeta(ops_profiler::RunMeta &meta, Domain *theDomain)
{
    if (theDomain != 0) {
        meta.nElem = static_cast<long long>(theDomain->getNumElements());
        meta.nNode = static_cast<long long>(theDomain->getNumNodes());
    }
    // `nnz` is intentionally NOT filled: LinearSOE has no size-agnostic nnz
    // accessor (only some concrete SOEs carry the member), so populating it means
    // a virtual on an upstream base class. Left a uniform 0 rather than filled
    // for some solvers and not others -- a normalizer that exists only sometimes
    // is worse than one that never does.

#ifdef _PARDISO
    // Profiler::buildMeta() has already resolved the ENV-declared cap. MKL is
    // authoritative where it exists: it sees a programmatic mkl_set_num_threads()
    // and it reports PHYSICAL cores when nothing is declared, which is the number
    // the env fallback can only estimate. Guarded > 0 so a surprising MKL answer
    // can never make the attribute worse than what buildMeta() resolved.
    const int mklThreads = mkl_get_max_threads();
    if (mklThreads > 0)
        meta.threads = static_cast<long long>(mklThreads);
#endif
}

#endif // ProfilerRunMeta_h
