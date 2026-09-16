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
// LADRUNO_THREADS environment variable; an absent or unparsable value gives 1.
int ladruno_getNumThreads(void);

// Set the requested element-loop thread count. Values < 1 are clamped to 1.
// Returns the value actually stored.
int ladruno_setNumThreads(int n);

// True when this binary was compiled with LADRUNO_OPENMP=ON, so the verb can
// report honestly instead of silently accepting a count it cannot honour.
bool ladruno_openmpCompiledIn(void);

// What the OpenMP runtime would give us, or 1 when OpenMP is compiled out.
// Diagnostic only -- it does not gate anything.
int ladruno_openmpMaxThreads(void);

#endif
