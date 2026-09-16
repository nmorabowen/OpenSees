/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno WP-107 (ADR-75b L3-1, desktop-scoped). See LadrunoThreads.h for the
// policy this file implements.

#include "LadrunoThreads.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef _LADRUNO_OPENMP
#include <omp.h>
#endif

// -1 = "not seeded yet". Seeded once, lazily, from LADRUNO_THREADS so that a
// bench harness can set the count without editing a deck.
static int ladrunoNumThreads = -1;

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
    char *end = 0;
    long v = strtol(env, &end, 10);
    if (end == env || v < 1)
        return;
    if (v > 1024)
        v = 1024;

    ladrunoNumThreads = (int)v;
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
    ladrunoNumThreads = (n < 1) ? 1 : ((n > 1024) ? 1024 : n);
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
