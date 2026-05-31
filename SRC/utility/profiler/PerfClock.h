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

// Ladruno Stack Profiler — Phase P1 (timing core).
//
// File: SRC/utility/profiler/PerfClock.h
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: PerfClock — a tiny, self-contained, portable wall+CPU clock.
//   Replaces the no-op Win32 stub in SRC/utility/Timer.cpp (which returns 0.0
//   for getReal()/getCPU() on Windows). This is the timing primitive the
//   Profiler registry and the OPS_PROFILE_SCOPE macros are built on.
//
//   Design constraints (locked, 06_profiler.md):
//     - NO OpenSees dependencies. NO HDF5. NO dynamic allocation.
//     - NO global/static mutable state -> reentrant, thread-safe.
//     - <windows.h> MUST NOT leak through this header (it is included only in
//       the .cpp). Keeping the header clean lets the timing core be pulled into
//       hot paths without dragging the Win32 API into every translation unit.
//
//   Time units: NANOSECONDS, signed 64-bit (int64), matching the HDF5 contract
//   (profiler_schema.py: wall_ns / cpu_ns are i8). A monotonic int64 ns counter
//   does not overflow for ~292 years, so plain subtraction of two reads is safe.

#ifndef PerfClock_h
#define PerfClock_h

#include <cstdint>

// PerfClock: stateless namespace of static clock reads. No instances needed;
// callers sample wall_ns()/cpu_ns() at scope entry/exit and subtract.
struct PerfClock {

    // Monotonic wall-clock time in nanoseconds since an unspecified epoch.
    // Backed by std::chrono::steady_clock, which on MSVC maps to
    // QueryPerformanceCounter (sub-microsecond, never runs backwards).
    // Only DIFFERENCES of two reads are meaningful (the epoch is arbitrary).
    static int64_t wall_ns();

    // Process CPU time (kernel + user, this process, all threads) in nanoseconds.
    //   Windows: GetProcessTimes (100 ns FILETIME ticks).
    //   POSIX:   getrusage(RUSAGE_SELF) utime + stime.
    //
    // NOTE ON GRANULARITY: the Windows scheduler accounts CPU time in
    // ~15.6 ms quanta, so cpu_ns() is COARSE on Windows. It is meaningful only
    // over a whole run (the run-level wall-vs-cpu parallelism signal); it is NOT
    // reliable for a single analysis step. Per-step time series therefore record
    // WALL time only (see 06_profiler.md overhead budget).
    static int64_t cpu_ns();
};

#endif // PerfClock_h
