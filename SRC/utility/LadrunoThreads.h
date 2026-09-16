/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno WP-107 (ADR-75b L3-1, desktop-scoped). Shared-memory thread policy
// for the fork's ONE threaded loop: the element state-determination loop in
// Domain::update() -- "loop A" in ADR-75b section 2.
//
// WHY THIS FILE IS THE ONLY KNOB (ADR-75b P-5 / stage L3-5). MKL solver threads
// x assembly threads x MPI ranks oversubscribe and make every bench lie, so the
// thread count lives in exactly one place and is set by exactly one verb
// (`ladrunoThreads` in Tcl, `ops.ladrunoThreads` in Python), seeded from the
// LADRUNO_THREADS environment variable.
//
// DEFAULT IS 1 -- ADR-40's standing anti-goal is "OpenMP-by-default". A run is
// threaded only when explicitly asked for, and at 1 thread Domain::update takes
// the byte-identical serial path (not a one-thread parallel region).
//
// THE GLOBAL OpenMP THREAD COUNT IS NEVER RAISED. The fork carries 7
// pre-existing `#pragma omp` lines in PFEM (ADR-75b section 1) that were dead
// no-ops while nothing passed /openmp. Compiling the serial targets WITH
// OpenMP would silently activate them. So nothing here calls
// omp_set_num_threads(); the count is applied ONLY through an explicit
// num_threads(n) clause on Domain::update's own region. PFEM's pragmas
// therefore keep running on the runtime default, exactly as before this WP.

#ifndef LadrunoThreads_h
#define LadrunoThreads_h

// The requested element-loop thread count (>= 1). The first call seeds from the
// LADRUNO_THREADS environment variable.
//
// WP-107 red-team S7: the seed is VALIDATED and every rejection is announced.
// It used to fail silently -- `LADRUNO_THREADS=abc` and `=0` both ran serial
// with no message at all, and `=99999` was honoured as 1024 threads on an
// 8-core box (measured 5x SLOWER than serial). A typo'd env var is exactly
// "how a bench lies", which is the anti-goal this file exists to serve. Now:
// unparsable, trailing garbage, or < 1 -> warn and fall back to 1; a value
// above the box's hardware concurrency -> warn and clamp to it.
int ladruno_getNumThreads(void);

// Set the requested element-loop thread count. Values < 1 are clamped to 1 and
// values above hardware concurrency are clamped down, both with a warning
// (WP-107 red-team S7). Returns the value actually stored.
int ladruno_setNumThreads(int n);

// True when this binary was compiled with LADRUNO_OPENMP=ON, so the verb can
// report honestly instead of silently accepting a count it cannot honour.
bool ladruno_openmpCompiledIn(void);

// What the OpenMP runtime would give us, or 1 when OpenMP is compiled out.
// Diagnostic only -- it does not gate anything. (Quoted in the S7 clamp
// warning, which is what finally gave it a caller -- red-team N1.)
int ladruno_openmpMaxThreads(void);

// The box's hardware concurrency: omp_get_num_procs() when OpenMP is compiled
// in, std::thread::hardware_concurrency() otherwise, and 0 when neither can
// answer -- in which case nothing is clamped. This is the S7 ceiling.
int ladruno_hardwareThreads(void);

// ---------------------------------------------------------------------------
// WP-107 red-team B2 -- the MPI fence. DEFINED IN LadrunoParallelBuild.cpp,
// declared here because Domain.cpp already includes this header.
//
// WHY IT IS NOT AN #ifdef IN Domain.cpp. `_PARALLEL_INTERPRETERS` /
// `_PARALLEL_PROCESSING` are PUBLIC compile definitions on the OpenSeesSP /
// OpenSeesMP / OpenSeesPyMP TARGETS. Domain.cpp lives in the OPS_Domain OBJECT
// library, which is compiled exactly ONCE with neither define and folded into
// all five targets -- so an `#ifdef _PARALLEL_*` written in Domain.cpp compiles
// to nothing in every build, OpenSeesMP included. That is not a guess: ADR-78
// P1 measured it (see SRC/analysis/handler/LadrunoContactAbort.h, whose guarded
// MPI_Abort silently vanished for exactly this reason), and the answer there --
// a one-function translation unit listed in OPS_MPI_PER_TARGET_SOURCES and
// compiled once per target -- is the answer here.
// ---------------------------------------------------------------------------

// True when THIS target was compiled with _PARALLEL_INTERPRETERS or
// _PARALLEL_PROCESSING, i.e. it is OpenSeesSP / OpenSeesMP / OpenSeesPyMP.
bool ladruno_parallelBuild(void);

// MPI_COMM_WORLD's size on a parallel build (1 when MPI is not initialised),
// and 0 on a serial build to mean "no MPI in this binary at all".
int ladruno_parallelRanks(void);

#endif
