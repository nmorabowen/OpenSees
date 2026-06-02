/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
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

// Ladruno Stack Profiler — Phase P1 (timing core).
//
// File: SRC/utility/profiler/PerfClock.cpp
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: Implementation of PerfClock — a tiny, self-contained, portable
//   wall+CPU clock. wall_ns() is backed by std::chrono::steady_clock; cpu_ns()
//   forks on the platform: GetProcessTimes on Windows, getrusage(RUSAGE_SELF)
//   on POSIX. All platform headers (<windows.h>) are confined to this
//   translation unit so PerfClock.h stays clean and dependency-free.
//
//   Both methods return nanoseconds as signed int64 (HDF5 contract). Neither
//   throws; on API failure cpu_ns() returns 0.

#include "PerfClock.h"

#include <chrono>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

// ---------------------------------------------------------------------------
// wall_ns(): monotonic wall-clock time in nanoseconds since an unspecified
// epoch. steady_clock never runs backwards; only differences are meaningful.
// ---------------------------------------------------------------------------
int64_t
PerfClock::wall_ns()
{
    using namespace std::chrono;
    return static_cast<int64_t>(
        duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count());
}

// ---------------------------------------------------------------------------
// cpu_ns(): process CPU time (kernel + user, this process, all threads) in
// nanoseconds. Returns 0 on API failure; never throws.
// ---------------------------------------------------------------------------
int64_t
PerfClock::cpu_ns()
{
#if defined(_WIN32)

    // GetProcessTimes reports kernel+user time as FILETIME in 100 ns ticks.
    // Creation/exit times are ignored. Multiply the 100 ns tick count by 100
    // to obtain nanoseconds. NOTE: the Windows scheduler accounts CPU time in
    // ~15.6 ms quanta, so this is COARSE — meaningful only over a whole run.
    FILETIME createTime, exitTime, kernelTime, userTime;
    if (GetProcessTimes(GetCurrentProcess(),
                        &createTime, &exitTime,
                        &kernelTime, &userTime) == 0) {
        return 0; // API failure: do not throw, report zero.
    }

    // Reassemble each FILETIME (low/high DWORD) into a 64-bit tick count.
    ULARGE_INTEGER k, u;
    k.LowPart  = kernelTime.dwLowDateTime;
    k.HighPart = kernelTime.dwHighDateTime;
    u.LowPart  = userTime.dwLowDateTime;
    u.HighPart = userTime.dwHighDateTime;

    // (kernel + user) ticks * 100 ns/tick -> nanoseconds.
    const uint64_t ticks100ns = k.QuadPart + u.QuadPart;
    return static_cast<int64_t>(ticks100ns * 100ULL);

#else

    // getrusage(RUSAGE_SELF) sums CPU time over all threads of this process.
    // ru_utime (user) + ru_stime (system), each a timeval {sec, usec}.
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) {
        return 0; // API failure: do not throw, report zero.
    }

    const int64_t userNs =
        static_cast<int64_t>(ru.ru_utime.tv_sec)  * 1000000000LL +
        static_cast<int64_t>(ru.ru_utime.tv_usec) * 1000LL;
    const int64_t sysNs =
        static_cast<int64_t>(ru.ru_stime.tv_sec)  * 1000000000LL +
        static_cast<int64_t>(ru.ru_stime.tv_usec) * 1000LL;

    return userNs + sysNs;

#endif
}
