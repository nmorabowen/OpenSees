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

// Ladruno Stack Profiler — Phase P1 (instrumentation macros).
//
// File: SRC/utility/profiler/ProfilerMacros.h
//
// Written: Ladruno / nmora
// Created: 2026-05
//
// Description: The instrumentation surface engine code uses. SINGLE-BINARY model
//   (06_profiler.md, revised 2026-05-30): everything is compiled into one shipped
//   binary and gated by RUNTIME flags, so the user controls all three layers per
//   run without a special build. Cost when off is one predicted-not-taken branch.
//
//     OPS_PROFILE_SCOPE(name)         COARSE phase timing. Gated by enabled().
//     OPS_PROFILE_SCOPE_DEEP(name)    DEEP per-element timing. Gated by
//     OPS_PROFILE_SCOPE_DEEP_NAMED    enabled() && deep() — the extra runtime flag,
//     OPS_PROFILE_ELEM(t, ...)        NOT a compile flag.
//     OPS_PROFILE_COUNT_ALLOC(bytes)  ALLOC counters. Gated by mem(). The branch on
//     OPS_PROFILE_COUNT_FREE(bytes)   mem() costs ~1ns against a heap new[] (~100ns).
//
//   The OPTIONAL build switch -DOPS_PROFILE_DISABLE compiles EVERY macro to
//   `((void)0)` — no symbol reference, no theProfiler() call, literally zero (not
//   even a branch). Use only for a benchmark/purist build; the default binary keeps
//   the profiler in and runtime-gated.

#ifndef ProfilerMacros_h
#define ProfilerMacros_h

#include "Profiler.h"   // ScopedTimerOpt, ScopedTimer, theProfiler, countAlloc/countFree

// ---- token-paste helpers (unique per-line variable names) ----
#define OPS_PROF_CAT_(a, b) a##b
#define OPS_PROF_CAT(a, b)  OPS_PROF_CAT_(a, b)
#define OPS_PROF_UNIQ(base) OPS_PROF_CAT(base, __LINE__)

#ifndef OPS_PROFILE_DISABLE

  // Coarse: the flag is checked BEFORE the timer is constructed — ScopedTimerOpt
  // gets the name only when enabled, else nullptr, so a disabled run pays just the
  // relaxed-atomic load + an empty optional. Uniquely-named local so two scopes in
  // one function don't collide.
  #define OPS_PROFILE_SCOPE(name)                                                 \
      ::ops_profiler::ScopedTimerOpt OPS_PROF_UNIQ(_ops_prof_)(                    \
          ::ops_profiler::theProfiler().enabled() ? (name) : nullptr)

  // Deep per-element scope: gated by enabled() && deep(). Same ScopedTimerOpt
  // wrapper (no-op when off), so the deep layer is a runtime toggle, not a rebuild.
  #define OPS_PROFILE_SCOPE_DEEP(name)                                            \
      ::ops_profiler::ScopedTimerOpt OPS_PROF_UNIQ(_ops_prof_deep_)(              \
          (::ops_profiler::theProfiler().enabled() &&                             \
           ::ops_profiler::theProfiler().deep()) ? (name) : nullptr)

  // Named deep scope: yields a handle `tmr` usable with OPS_PROFILE_ELEM. Gated
  // identically; when off `tmr` holds an empty optional and addElem is a no-op.
  #define OPS_PROFILE_SCOPE_DEEP_NAMED(tmr, name)                                 \
      ::ops_profiler::ScopedTimerOpt tmr(                                         \
          (::ops_profiler::theProfiler().enabled() &&                             \
           ::ops_profiler::theProfiler().deep()) ? (name) : nullptr)

  // Per-element-class attribution onto a named deep-scope handle (P0#1).
  #define OPS_PROFILE_ELEM(tmr, classTag, wall_ns, fb)                            \
      do { (tmr).addElem((classTag), (wall_ns), (fb)); } while (0)

  // RAII per-element timer around an FE_Element kernel call (P3). Expands in the
  // caller's TU, where FE_Element/Element are complete. The `tmr.engaged()` guard
  // short-circuits so `getElement()` (a non-inline call) is touched ONLY when the
  // deep scope is live — an unprofiled run pays nothing. A constraint FE_Element
  // (getElement()==0) yields classTag -1, which ElemScope treats as inert.
  #define OPS_PROFILE_FE_ELEM_SCOPE(tmr, fe_elem)                                 \
      ::ops_profiler::ElemScope OPS_PROF_UNIQ(_ops_elem_)(                         \
          (tmr),                                                                  \
          ((tmr).engaged() && (fe_elem)->getElement())                           \
              ? (fe_elem)->getElement()->getClassTag() : -1)

  // Alloc counters (P4 wires Matrix/Vector/ID ctors/dtors). Runtime-gated on mem().
  // `type` is an ops_profiler::AllocType (ALLOC_MATRIX/ALLOC_VECTOR/ALLOC_ID) for
  // the per-type live split; an out-of-range value counts toward the aggregate only.
  #define OPS_PROFILE_COUNT_ALLOC(bytes, type)                                    \
      do { if (::ops_profiler::theProfiler().mem())                               \
               ::ops_profiler::countAlloc(static_cast<int64_t>(bytes), (type)); } while (0)
  #define OPS_PROFILE_COUNT_FREE(bytes, type)                                     \
      do { if (::ops_profiler::theProfiler().mem())                               \
               ::ops_profiler::countFree(static_cast<int64_t>(bytes), (type)); } while (0)

  // Per-step series hook (P0#3): record one row of the per-step time history
  // (step index, pseudo-time, dt, iteration count) from the analysis driver.
  // Gated on enabled() so the argument expressions (getCommitTag / getCurrentTime
  // / getNumIterations) are evaluated ONLY when a run is active; recordStep itself
  // additionally early-returns unless the run armed -perStep. Negligible vs a step.
  #define OPS_PROFILE_STEP(step, t, dt, iters)                                    \
      do { if (::ops_profiler::theProfiler().enabled())                          \
               ::ops_profiler::theProfiler().recordStep((step), (t), (dt), (iters)); } while (0)

  // TaggedObject census (P4): one live-object count per classTag, bumped in the
  // MovableObject ctor and decremented in its dtor. Runtime-gated on mem() like
  // the byte counters; off-path cost is one predicted branch on a hot path.
  #define OPS_PROFILE_CENSUS_BORN(classTag)                                       \
      do { if (::ops_profiler::theProfiler().mem())                               \
               ::ops_profiler::countComponentBorn((classTag)); } while (0)
  #define OPS_PROFILE_CENSUS_DIED(classTag)                                       \
      do { if (::ops_profiler::theProfiler().mem())                               \
               ::ops_profiler::countComponentDied((classTag)); } while (0)

#else // OPS_PROFILE_DISABLE -> the profiler is compiled fully out

  #define OPS_PROFILE_SCOPE(name)                     ((void)0)
  #define OPS_PROFILE_SCOPE_DEEP(name)                ((void)0)
  #define OPS_PROFILE_SCOPE_DEEP_NAMED(tmr, name)     ((void)0)
  #define OPS_PROFILE_ELEM(tmr, classTag, w, fb)      ((void)0)
  #define OPS_PROFILE_FE_ELEM_SCOPE(tmr, fe_elem)     ((void)0)
  #define OPS_PROFILE_COUNT_ALLOC(bytes, type)        ((void)0)
  #define OPS_PROFILE_COUNT_FREE(bytes, type)         ((void)0)
  #define OPS_PROFILE_STEP(step, t, dt, iters)        ((void)0)
  #define OPS_PROFILE_CENSUS_BORN(classTag)           ((void)0)
  #define OPS_PROFILE_CENSUS_DIED(classTag)           ((void)0)

#endif // !OPS_PROFILE_DISABLE

#endif // ProfilerMacros_h
