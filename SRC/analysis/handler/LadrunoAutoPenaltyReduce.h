// LADRUNO-HEADER-START
// ==========================================================================
//  Ladruno — a research fork of OpenSees
// ==========================================================================
// LADRUNO-HEADER-END

// ADR-78 follow-up — the cross-rank reduction behind `constraints Auto`.
//
// Same TU-isolation reason as LadrunoContactAbort.h, and the same trap.
// AutoConstraintHandler.cpp is compiled into the OPS_Analysis OBJECT library,
// which is built ONCE with no parallel define and folded into every target, so
// the `#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)`
// block that used to hold this reduction compiled to nothing in EVERY binary,
// OpenSeesMP included. It has been dead since the handler was contributed.
//
// The consequence is silent and physical: `constraints Auto` sizes its penalty
// from the order of magnitude of the AVERAGE diagonal stiffness. With the
// reduction gone each rank averages only its OWN elements, so a partitioned
// model gives every rank a different m_global_penalty. That penalty is the
// value used for any constraint whose nodes are not in the local node map --
// precisely the interface constraints that straddle a partition -- so the two
// ranks enforce the SAME constraint with DIFFERENT stiffnesses and the
// distributed SOE sums two mismatched rows. No warning is printed anywhere.
//
// So the one parallel-sensitive step lives here, in its own TU, listed in
// OPS_MPI_PER_TARGET_SOURCES and compiled per target with that target's
// defines. AutoConstraintHandler.cpp is left with NO `#ifdef _PARALLEL_*` at
// all, which is what stops the trap from silently reopening.

#ifndef LadrunoAutoPenaltyReduce_h
#define LadrunoAutoPenaltyReduce_h

// Reduce the four global-stiffness statistics that AutoConstraintHandler's
// penalty evaluator accumulates from its rank-local elements, so every rank
// finishes with the same numbers and therefore the same global penalty:
//
//   gpMin <- MIN over ranks, gpMax <- MAX over ranks,
//   gpAvg <- SUM over ranks, gpCnt <- SUM over ranks
//
// (gpAvg/gpCnt are the running sum and count; the caller forms the mean after
// this returns, so summing both is what makes the mean global.)
//
// Returns the number of ranks the statistics were reduced over. 1 means no
// reduction happened and the arguments are untouched -- a serial build, a
// single-rank MPI job, or the _PARALLEL_PROCESSING master running an
// unpartitioned model (where the peers are actors that never reach this code,
// so entering a collective alone would deadlock).
int ladrunoAutoPenaltyReduce(double &gpMin, double &gpMax, double &gpAvg, double &gpCnt);

#endif
