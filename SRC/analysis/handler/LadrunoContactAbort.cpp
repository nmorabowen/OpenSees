// LADRUNO-HEADER-START
// ==========================================================================
//  Ladruno — a research fork of OpenSees
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-78 P1 — see LadrunoContactAbort.h for why this is its own TU.

#include "LadrunoContactAbort.h"
#include <OPS_Globals.h>

// Guard style copied from AutoConstraintHandler.cpp. NB that file is itself in
// the OPS_Analysis OBJECT library and is NOT listed per-target, so its own
// MPI_Allreduce is compiled out of every build -- a separate, pre-existing bug
// (constraints Auto sizes its penalty from rank-local stiffness under MPI).
// Flagged, deliberately not fixed here: it is unrelated to contact and deserves
// its own change rather than riding along in P1.
#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#endif // _PARALLEL_PROCESSING
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif // _PARALLEL_INTERPRETERS

int
ladrunoContactFatal()
{
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > 1) {
        // np == 1 falls through to the plain return: no peers to strand, and it
        // keeps single-rank MP runs debuggable.
        opserr << "LadrunoContactHandler: tearing down all " << np << " MPI ranks -- a "
                  "contact failure on one rank leaves every peer blocked in the next "
                  "collective, so returning alone would HANG the job (ADR-78 P1).\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    return -1;
}

int
ladrunoContactNumRanks()
{
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    return np;
#else
    return 1;
#endif
}
