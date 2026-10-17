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
#include <stdio.h>
#include <ctype.h>
#include <thread>
#include <mutex>

#ifdef _LADRUNO_OPENMP
#include <omp.h>
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

#ifdef _LADRUNO_OPENMP
    if (nThreads >= 2 && n >= 2) {
        // An INTEGER reduction: order-free, identical at every thread count.
        // The lowest failing index is found through a critical section that
        // is entered only on a failure (MSVC is OpenMP 2.0: no reduction(min)).
        #pragma omp parallel num_threads(nThreads) reduction(+:sum)
        {
            ladruno_tlInRegion = 1;
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
        }
        if (firstFail != 0)
            *firstFail = first;
        return sum;
    }
#endif

    (void)nThreads;
    const int saved = ladruno_tlInRegion;
    ladruno_tlInRegion = 1;
    for (int i = 0; i < n; i++) {
        const int r = body(i, ctx);
        sum += r;
        if (r != 0) {
            if (i < first)
                first = i;
        }
    }
    ladruno_tlInRegion = saved;
    if (firstFail != 0)
        *firstFail = first;
    return sum;
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
