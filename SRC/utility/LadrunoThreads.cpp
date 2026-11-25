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

// Ladruno WP-107 (ADR-75b L3-1, desktop-scoped). See LadrunoThreads.h for the
// policy this file implements.

#include "LadrunoThreads.h"

#include <OPS_Globals.h>    // Ladruno WP-107 red-team S7: opserr for the refusals

#include <stdlib.h>
#include <climits>          // Ladruno WP-253: INT_MAX (KMP_MAX_BLOCKTIME)
#include <cstring>          // Ladruno WP-253: memcpy (import-stub decoding)
#include <stdio.h>
#include <ctype.h>
#include <thread>
#include <mutex>
#include <atomic>    // Ladruno WP-258: the loop serial number

#ifdef _LADRUNO_OPENMP
#include <omp.h>
#  ifdef _WIN32                 // Ladruno WP-253: run-time lookup of kmp_set_blocktime
#    ifndef WIN32_LEAN_AND_MEAN
#      define WIN32_LEAN_AND_MEAN
#    endif
#    ifndef NOMINMAX
#      define NOMINMAX
#    endif
#    ifndef NOGDI
#      define NOGDI
#    endif
#    include <windows.h>
#  else
#    include <dlfcn.h>
#  endif
#endif

// -1 = "not seeded yet". Seeded once, lazily, from LADRUNO_THREADS so that a
// bench harness can set the count without editing a deck.
static int ladrunoNumThreads = -1;

// An absolute ceiling kept below the hardware clamp: a box with an enormous
// core count still does not want a per-element parallel region that wide, and
// it bounds the num_threads() clause no matter what the runtime reports.
static const int LADRUNO_THREADS_HARD_CAP = 1024;

// WP-107 red-team S7. Clamp a requested count into [1, hardware], announcing
// every correction. `who` names the path so the message is actionable
// ("LADRUNO_THREADS" vs "ladrunoThreads").
static int
ladruno_clamp(long v, const char *who)
{
    if (v < 1) {
        opserr << "WARNING " << who << ": requested " << (int)v
               << " threads, which is < 1 -- using 1 (SERIAL). The element "
               << "loop is never threaded below 2.\n";
        return 1;
    }

    // Display guard only: strtol can hand us something far outside int range,
    // and the warning must still quote what the user actually asked for.
    const int shown = (int)((v > 2000000000L) ? 2000000000L : v);

    const int hw = ladruno_hardwareThreads();
    if (hw > 0 && v > (long)hw) {
        opserr << "WARNING " << who << ": requested " << shown
               << " threads but this box reports " << hw
               << " hardware threads (OpenMP max " << ladruno_openmpMaxThreads()
               << ") -- CLAMPED to " << hw
               << ". Oversubscribing the element loop is measurably slower "
               << "than serial, and a bench that does it lies.\n";
        return hw;
    }

    if (v > (long)LADRUNO_THREADS_HARD_CAP) {
        opserr << "WARNING " << who << ": requested " << shown
               << " threads -- CLAMPED to the hard cap " << LADRUNO_THREADS_HARD_CAP
               << " (hardware concurrency could not be determined).\n";
        return LADRUNO_THREADS_HARD_CAP;
    }

    return (int)v;
}

static void
ladruno_seedFromEnv(void)
{
    ladrunoNumThreads = 1;

    const char *env = getenv("LADRUNO_THREADS");
    if (env == 0 || env[0] == '\0')
        return;

    // strtol, not atoi: atoi cannot tell "0" from garbage, and silently
    // threading on a typo is exactly the class of surprise ADR-40's anti-goal
    // exists to prevent.
    //
    // WP-107 red-team S7: the three silent-failure exits below are now loud.
    // The env path is the one a bench harness and every job script use, so a
    // typo here used to produce a serial run labelled as a threaded one.
    char *end = 0;
    long v = strtol(env, &end, 10);

    if (end == env) {
        opserr << "WARNING LADRUNO_THREADS: \"" << env
               << "\" is not an integer -- IGNORED, the element loop runs "
               << "SERIAL (1 thread).\n";
        return;
    }

    // Trailing garbage ("4x", "2.7") is a typo, not a count. Whitespace is not.
    while (*end != '\0' && isspace((unsigned char)*end))
        end++;
    if (*end != '\0') {
        opserr << "WARNING LADRUNO_THREADS: \"" << env
               << "\" has trailing text after the number -- IGNORED, the "
               << "element loop runs SERIAL (1 thread).\n";
        return;
    }

    ladrunoNumThreads = ladruno_clamp(v, "LADRUNO_THREADS");
}

int
ladruno_getNumThreads(void)
{
    if (ladrunoNumThreads < 0)
        ladruno_seedFromEnv();
    return ladrunoNumThreads;
}

int
ladruno_setNumThreads(int n)
{
    if (ladrunoNumThreads < 0)
        ladruno_seedFromEnv();
    ladrunoNumThreads = ladruno_clamp((long)n, "ladrunoThreads");
    return ladrunoNumThreads;
}

bool
ladruno_openmpCompiledIn(void)
{
#ifdef _LADRUNO_OPENMP
    return true;
#else
    return false;
#endif
}

int
ladruno_openmpMaxThreads(void)
{
#ifdef _LADRUNO_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// WP-107 red-team S7. omp_get_num_procs() rather than omp_get_max_threads():
// the latter reports OMP_NUM_THREADS if the environment set it, which is a
// POLICY, not the hardware. Clamping a deliberate `ladrunoThreads 8` down to a
// stray OMP_NUM_THREADS=2 would be a surprise; clamping it to the core count is
// not. Falls back to the C++11 query when OpenMP is compiled out, and returns 0
// ("unknown") when that also answers 0, in which case the caller clamps nothing.
int
ladruno_hardwareThreads(void)
{
#ifdef _LADRUNO_OPENMP
    const int n = omp_get_num_procs();
    if (n > 0)
        return n;
#endif
    const unsigned hc = std::thread::hardware_concurrency();
    return (hc == 0u) ? 0 : (int)hc;
}

// ---------------------------------------------------------------------------
// Ladruno WP-192 -- shared machinery of the threaded element loops (see the
// header for the contract).
// ---------------------------------------------------------------------------

thread_local int ladruno_tlInRegion = 0;
thread_local long long ladruno_tlLoopId = 0;   // Ladruno WP-258

// Ladruno WP-258: one number per ladruno_parallelFor call (see the header).
static long long ladruno_nextLoopId(void)
{
    static std::atomic<long long> serial(0);
    return ++serial;
}

// -1 = not seeded yet (from LADRUNO_THREADS_HYBRID), 0 = refused, 1 = allowed.
static int ladrunoHybrid = -1;

bool
ladruno_hybridAllowed(void)
{
    if (ladrunoHybrid < 0) {
        const char *env = getenv("LADRUNO_THREADS_HYBRID");
        ladrunoHybrid = (env != 0 && env[0] == '1' && env[1] == '\0') ? 1 : 0;
    }
    return ladrunoHybrid == 1;
}

void
ladruno_setHybridAllowed(bool on)
{
    ladrunoHybrid = on ? 1 : 0;
}

int
ladruno_loopThreads(const char **why)
{
    if (why != 0)
        *why = 0;
#ifndef _LADRUNO_OPENMP
    return 1;
#else
    const int n = ladruno_getNumThreads();
    if (n <= 1)
        return 1;
    if (ladruno_parallelBuild() && !ladruno_hybridAllowed()) {
        if (why != 0)
            *why = "this is an MPI binary and hybrid MPI+threads was not "
                   "acknowledged (`ladrunoThreads <n> -hybrid` or "
                   "LADRUNO_THREADS_HYBRID=1)";
        return 1;
    }
    return n;
#endif
}

int
ladruno_parallelFor(int n, int nThreads, int chunk,
                    LadrunoIndexBody body, void *ctx, int *firstFail)
{
    int sum = 0;
    int first = n;
    if (chunk < 1)
        chunk = 1;
    const long long loopId = ladruno_nextLoopId();   // Ladruno WP-258

#ifdef _LADRUNO_OPENMP
    if (nThreads >= 2 && n >= 2) {
        // An INTEGER reduction: order-free, identical at every thread count.
        // The lowest failing index is found through a critical section that
        // is entered only on a failure (MSVC is OpenMP 2.0: no reduction(min)).
        #pragma omp parallel num_threads(nThreads) reduction(+:sum)
        {
            ladruno_tlInRegion = 1;
            ladruno_tlLoopId = loopId;   // Ladruno WP-258
            #pragma omp for schedule(dynamic, chunk)
            for (int i = 0; i < n; i++) {
                const int r = body(i, ctx);
                sum += r;
                if (r != 0) {
                    #pragma omp critical (ladruno_wp192_firstFail)
                    {
                        if (i < first)
                            first = i;
                    }
                }
            }
            ladruno_tlInRegion = 0;
            ladruno_tlLoopId = 0;   // Ladruno WP-258
        }
        if (firstFail != 0)
            *firstFail = first;
        return sum;
    }
#endif

    (void)nThreads;
    const int saved = ladruno_tlInRegion;
    const long long savedLoop = ladruno_tlLoopId;   // Ladruno WP-258
    ladruno_tlInRegion = 1;
    ladruno_tlLoopId = loopId;
    for (int i = 0; i < n; i++) {
        const int r = body(i, ctx);
        sum += r;
        if (r != 0) {
            if (i < first)
                first = i;
        }
    }
    ladruno_tlInRegion = saved;
    ladruno_tlLoopId = savedLoop;   // Ladruno WP-258
    if (firstFail != 0)
        *firstFail = first;
    return sum;
}

// ---------------------------------------------------------------------------
// Ladruno WP-253 -- the OpenMP quiet span (see LadrunoThreads.h). Relies on the LLVM / Intel
// OpenMP runtime (libiomp5 / libomp), read in its source (kmp_ftn_entry.h, kmp_runtime.cpp,
// kmp_barrier.cpp): kmp_set_blocktime(v) converts v from KMP_BLOCKTIME's units and stores it in
// the CALLING thread's ICVs (td_icvs.blocktime) with bt_set = TRUE; the fork copies the master's
// ICVs to the team's implicit tasks; each worker sets its wait interval (th_team_bt_intervals)
// from them in the join barrier and waits with it for the next fork. kmp_get_blocktime()
// returns the effective value in the same units, so set(get()) restores the VALUE for every
// later region (WP-192's, MKL's). It does not restore bt_set: from the first span on, the
// master's blocktime counts as user-set, which switches off libiomp's oversubscription auto-zero
// (KMP_ADJUST_BLOCKTIME: blocktime 0 while threads > procs and bt_set is false). Timing only
// (LEDGER quirk WP-253): no result depends on how long a worker waits.
//
// The entry points are taken from the module that OWNS the running OpenMP runtime -- the one
// defining omp_get_max_threads as this translation unit resolved it -- never from the first DLL
// found by name. A runtime without them (libgomp, vcomp) leaves the span a no-op.
// ---------------------------------------------------------------------------
namespace {
typedef void (*LadrunoKmpSetBt)(int);
typedef int (*LadrunoKmpGetBt)(void);
struct LadrunoKmpBt {
    LadrunoKmpSetBt set = 0;
    LadrunoKmpGetBt get = 0;
};

#if defined(_LADRUNO_OPENMP) && defined(_WIN32)
// The module containing address a, or 0 (reference count untouched).
HMODULE ladruno_moduleOf(const void *a)
{
    HMODULE h = 0;
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            reinterpret_cast<LPCSTR>(a), &h))
        return 0;
    return h;
}

// &omp_get_max_threads is this module's own import stub when omp.h does not declare the function
// dllimport: an incremental-link `jmp rel32` (E9) and/or the x64 import thunk `jmp [rip+rel32]`
// (FF 25) that reads the IAT slot. Follow at most two hops while the address stays inside this
// module; anything else is returned unchanged (and then owns no kmp_*: a no-op).
const void *ladruno_followImportStub(const void *a)
{
#  if defined(_M_X64) || defined(__x86_64__)
    const HMODULE self = ladruno_moduleOf(reinterpret_cast<const void *>(&ladruno_moduleOf));
    for (int hop = 0; hop < 2 && a != 0 && self != 0 && ladruno_moduleOf(a) == self; hop++) {
        const unsigned char *c = static_cast<const unsigned char *>(a);
        int rel;
        if (c[0] == 0xE9) {
            memcpy(&rel, c + 1, sizeof rel);
            a = c + 5 + rel;
        } else if (c[0] == 0xFF && c[1] == 0x25) {
            memcpy(&rel, c + 2, sizeof rel);
            const void *slot;
            memcpy(&slot, c + 6 + rel, sizeof slot);
            a = slot;
        } else {
            break;
        }
    }
#  endif
    return a;
}
#endif

const LadrunoKmpBt &ladruno_kmpBlocktime(void)
{
    static const LadrunoKmpBt bt = [] {
        LadrunoKmpBt b;
#ifdef _LADRUNO_OPENMP
#  ifdef _WIN32
        const void *ompFn = ladruno_followImportStub(
            reinterpret_cast<const void *>(&omp_get_max_threads));
        const HMODULE owner = ladruno_moduleOf(ompFn);
        if (owner != 0) {
            b.set = reinterpret_cast<LadrunoKmpSetBt>(
                reinterpret_cast<void *>(GetProcAddress(owner, "kmp_set_blocktime")));
            b.get = reinterpret_cast<LadrunoKmpGetBt>(
                reinterpret_cast<void *>(GetProcAddress(owner, "kmp_get_blocktime")));
        }
#  else
        Dl_info info;
        if (dladdr(reinterpret_cast<void *>(&omp_get_max_threads), &info) != 0 &&
            info.dli_fname != 0) {
            // the owner is loaded already: RTLD_NOLOAD only takes a handle (dlclose drops it)
            void *owner = dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD);
            if (owner != 0) {
                b.set = reinterpret_cast<LadrunoKmpSetBt>(dlsym(owner, "kmp_set_blocktime"));
                b.get = reinterpret_cast<LadrunoKmpGetBt>(dlsym(owner, "kmp_get_blocktime"));
                dlclose(owner);
            }
        }
#  endif
        if (b.set == 0 || b.get == 0) {
            b.set = 0;
            b.get = 0;
        }
#endif
        return b;
    }();
    return bt;
}
}

int
ladruno_ompQuietBegin(void)
{
    const LadrunoKmpBt &bt = ladruno_kmpBlocktime();
    if (bt.set == 0)
        return -1;
    const int previous = bt.get();
    // KMP_BLOCKTIME=infinite: get() reports INT_MAX whatever is set, and set(INT_MAX) prints an
    // "OMP: Info #109" clamp line on every call (measured, oneAPI 2026.1). A user who asked for
    // an infinite spin keeps it: nothing is changed.
    if (previous < 0 || previous == INT_MAX)
        return -1;
    bt.set(0);
    return previous;
}

void
ladruno_ompQuietEnd(int token)
{
    const LadrunoKmpBt &bt = ladruno_kmpBlocktime();
    if (bt.set == 0 || token < 0)
        return;
    bt.set(token);
}

static std::mutex &ladruno_outputMutex(void) { static std::mutex m; return m; }
void ladruno_outputLock(void) { ladruno_outputMutex().lock(); }
void ladruno_outputUnlock(void) { ladruno_outputMutex().unlock(); }

bool
ladruno_loopEnabled(int which)
{
    static int mask = -1;
    if (mask < 0) {
        const char *env = getenv("LADRUNO_THREADS_LOOPS");
        mask = (env != 0 && env[0] != 0) ? (int)strtol(env, 0, 10) : 31;   // WP-196: + energy sweep
    }
    return (mask & which) != 0;
}
