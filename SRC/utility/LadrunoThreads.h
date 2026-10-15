/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// ---------------------------------------------------------------------------
// Ladruno WP-192 -- the shared machinery of the threaded element loops.
//
// WP-107 threaded one loop (Domain::update) with an inline `#pragma omp`. WP-192
// threads three more (the element tangent and residual loops of the
// integrators, and the element commit loop of Domain::commit), and two of them
// live in OPS_Analysis, which is NOT compiled with /openmp (and must not be:
// see CMakeLists.txt, PFEM's dormant pragmas). So the parallel-for itself lives
// here, in OPS_Utilities, and the callers hand it a plain function pointer.
// ---------------------------------------------------------------------------

// Nonzero on a thread while it is executing a body of ladruno_parallelFor().
// It is how the class-wide scratch pools (FE_Element, TransformationFE,
// Element's Rayleigh pool, Matrix's solve work area) know to hand out the
// calling thread's private copy instead of the shared one. Zero everywhere
// else, so every serial path keeps using the original objects.
extern thread_local int ladruno_tlInRegion;
inline bool ladruno_inThreadedRegion(void) { return ladruno_tlInRegion != 0; }

// One iteration of a ladruno_parallelFor loop. Returns 0 on success.
typedef int (*LadrunoIndexBody)(int i, void *ctx);

// Run body(i, ctx) for every i in [0, n) on nThreads OpenMP threads
// (schedule(dynamic, chunk)). Every executing thread has ladruno_tlInRegion set
// for the duration. Returns the INTEGER sum of the body results (order-free,
// so identical at every thread count); *firstFail, when given, receives the
// lowest i whose body returned nonzero (n when none did). With OpenMP compiled
// out, or nThreads < 2, it runs the loop serially on the calling thread (with
// the region flag set, so the same scratch selection is exercised).
int ladruno_parallelFor(int n, int nThreads, int chunk,
                        LadrunoIndexBody body, void *ctx, int *firstFail);

// Hybrid MPI + threads (ADR-75b section 11 q6), WP-192. WP-107 refuses every
// MPI binary because an `mpiexec`-inherited LADRUNO_THREADS would oversubscribe
// the job silently. WP-192 keeps that refusal as the DEFAULT and adds an
// explicit acknowledgement: `ladrunoThreads <n> -hybrid` in the deck, or
// LADRUNO_THREADS_HYBRID=1 in the environment. Without it an MPI binary still
// runs serial and says so.
bool ladruno_hybridAllowed(void);
void ladruno_setHybridAllowed(bool on);

// The serial-path guard every WP-192 loop runs before threading: the requested
// count (> 1), OpenMP compiled in, and on an MPI binary the hybrid
// acknowledgement. Returns the thread count to use, or 1 for "run serial".
// *why receives a short reason when it returns 1 (0 = the plain default of
// 1 thread, nothing to announce).
int ladruno_loopThreads(const char **why);


// Serialises a diagnostic written to opserr from inside a threaded loop
// (opserr is ONE shared stream; WP-107 rule 4). A no-op outside a region, so
// every serial path prints exactly as before. Use as a scoped guard:
//   { LadrunoOutputGuard g; opserr << ...; }
void ladruno_outputLock(void);
void ladruno_outputUnlock(void);
struct LadrunoOutputGuard {
  bool on;
  LadrunoOutputGuard() : on(ladruno_inThreadedRegion()) { if (on) ladruno_outputLock(); }
  ~LadrunoOutputGuard() { if (on) ladruno_outputUnlock(); }
};


// Diagnostic mask of the WP-192 loops that may thread: LADRUNO_THREADS_LOOPS,
// a bitmask read once (1 = element update, 2 = tangent, 4 = residual,
// 8 = commit, 16 = the EnergyBalance element sweep (WP-196); default 31 = all). For bisecting a suspected race to one loop
// without a rebuild; production runs leave it unset.
enum { LADRUNO_LOOP_UPDATE = 1, LADRUNO_LOOP_TANGENT = 2,
       LADRUNO_LOOP_RESIDUAL = 4, LADRUNO_LOOP_COMMIT = 8,
       LADRUNO_LOOP_ENERGY = 16 };   // Ladruno WP-196
bool ladruno_loopEnabled(int which);

#endif
