// LADRUNO-HEADER-START
// ==========================================================================
//  Ladruno — a research fork of OpenSees
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-78 follow-up — see LadrunoAutoPenaltyReduce.h for why this is its own TU.
//
// The body below is the block that used to sit in AutoConstraintHandler.cpp,
// moved verbatim in behaviour: same four collectives, same order, same
// OPS_PARTITIONED quick-return, same warning strings. The only thing that
// changed is WHERE it is compiled -- which is the entire bug.

#include "LadrunoAutoPenaltyReduce.h"
#include <OPS_Globals.h>

#ifdef _PARALLEL_PROCESSING
extern bool OPS_PARTITIONED;
#include <mpi.h>
#endif // _PARALLEL_PROCESSING
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif // _PARALLEL_INTERPRETERS

int
ladrunoAutoPenaltyReduce(double &gpMin, double &gpMax, double &gpAvg, double &gpCnt)
{
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    int pid = 0;
    int np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // quick return for 1 process or for non-partitioned cases (should not happen)
    bool do_allreduce = true;
    if (np == 1) do_allreduce = false;
#if defined(_PARALLEL_PROCESSING)
    // OpenSeesSP: the peers are ActorSubdomain processes that never run this
    // handler, so an unpartitioned master must NOT enter a collective alone.
    if (pid == 0 && !OPS_PARTITIONED) do_allreduce = false;
#endif // defined(_PARALLEL_PROCESSING)
    if (do_allreduce) {
        double local_min = gpMin;
        double local_max = gpMax;
        double local_avg = gpAvg;
        double local_cnt = gpCnt;
        if (MPI_Allreduce(&local_min, &gpMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS) {
            opserr << "AutoConstraintHandler Warning: MPI_Allreduce failed to get MIN\n";
        }
        if (MPI_Allreduce(&local_max, &gpMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
            opserr << "AutoConstraintHandler Warning: MPI_Allreduce failed to get MAX\n";
        }
        if (MPI_Allreduce(&local_avg, &gpAvg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
            opserr << "AutoConstraintHandler Warning: MPI_Allreduce failed to get SUM\n";
        }
        if (MPI_Allreduce(&local_cnt, &gpCnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
            opserr << "AutoConstraintHandler Warning: MPI_Allreduce failed to get CNT\n";
        }
        return np;
    }
#else
    // Serial build: nothing to reduce. Referenced so the signature stays honest
    // and the compiler does not warn about unused parameters.
    (void)gpMin; (void)gpMax; (void)gpAvg; (void)gpCnt;
#endif // defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    return 1;
}
