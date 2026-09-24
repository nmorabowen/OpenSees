// LADRUNO-HEADER-START
// ==========================================================================
//  Ladruno — a research fork of OpenSees
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno WP-107 red-team B2 — the MPI fence for the threaded element loop.
//
// WHY THIS IS ITS OWN TRANSLATION UNIT, and why the #ifdef cannot live in
// Domain.cpp. `_PARALLEL_INTERPRETERS` / `_PARALLEL_PROCESSING` are PUBLIC
// compile definitions on the OpenSeesSP / OpenSeesMP / OpenSeesPyMP TARGETS.
// Domain.cpp is compiled into the OPS_Domain OBJECT library exactly ONCE, with
// neither define, and that single set of objects is folded into all five
// targets — so an `#ifdef _PARALLEL_*` block written in Domain.cpp compiles to
// nothing in every build, OpenSeesMP included.
//
// ADR-78 P1 paid for that lesson already: see LadrunoContactAbort.h, whose
// guarded MPI_Abort silently vanished from OpenSeesMP for exactly this reason
// and whose fix — a ~one-function file listed in OPS_MPI_PER_TARGET_SOURCES and
// added to each executable/module directly — this file copies verbatim.
//
// WHAT IT FENCES. WP-107 threads the element state-determination loop inside a
// single shared-memory process. Under `_PARALLEL_INTERPRETERS` (OpenSeesMP /
// OpenSeesPyMP) every rank holds a PLAIN `Domain` — SRC/tcl/commands.cpp's
// `#elif _PARALLEL_INTERPRETERS -> Domain theDomain;` and
// SRC/interpreter/OpenSeesCommands.cpp's `theDomain = new Domain;` — not a
// PartitionedDomain, so `Domain::ladrunoThreadedUpdateAllowed()` answers true
// and the PartitionedDomain/Subdomain overrides never see those binaries. The
// thread count is seeded from an ENVIRONMENT VARIABLE and `mpiexec` propagates
// the environment to every rank, so one `set LADRUNO_THREADS=8` in a job script
// would silently turn an np-8 run into 64-way oversubscription, with the single
// per-process announcement lost in rank-interleaved stdout. Hybrid MPI+threads
// is ADR-75b §11 q6 and is DEFERRED; this WP must not ship it by accident.
//
// THE FENCE IS THE BUILD, NOT THE RANK COUNT. It refuses on the compile
// definition alone, before asking MPI anything: `np == 1` under OpenSeesMP is
// still the parallel code path (parallel numberer, parallel interpreter
// semantics) and is not the desktop case this WP measured. `ladruno_parallelRanks()`
// exists only to make the refusal message say how many ranks are actually
// running, and it is safe to call before MPI_Init.

#include <LadrunoThreads.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>
#endif

bool
ladruno_parallelBuild(void)
{
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    return true;
#else
    return false;
#endif
}

int
ladruno_parallelRanks(void)
{
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    // MPI_Comm_size before MPI_Init is undefined, so ask MPI_Initialized first.
    // The fence itself does not depend on the answer (see the header comment);
    // this only decorates the message.
    int inited = 0;
    if (MPI_Initialized(&inited) != MPI_SUCCESS || !inited)
        return 1;
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    return np;
#else
    return 0;
#endif
}
