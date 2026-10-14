/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno WP-192 -- the threaded element tangent / residual loops ("loops B and
// C" of ADR-75b section 2) of the incremental integrators.
//
// THE ASSEMBLY STRATEGY: chunked compute, ordered serial scatter. The FE_Element
// list is walked in its serial order and cut into segments: a maximal run of
// consecutive thread-safe FEs (capped at a chunk size) is computed in parallel
// -- each FE's getTangent()/getResidual() copied into its own slot of a
// per-call buffer -- and then scattered into the SOE by ONE thread, in the
// serial order, with the very same addA/addB calls the serial loop makes. An FE
// that is not thread-safe is computed and scattered alone, in its serial place.
//
// So the SOE receives the identical sequence of (matrix, ID) contributions as
// the serial loop -- the replay order is (FE index, i, j), ADR-75b 4.2 (i) --
// and the assembled A and b are BIT-IDENTICAL to serial at every thread count.
// No graph colouring, no atomics, no per-thread copies of A; the buffer is
// bounded by the chunk (C x idSize^2 doubles, not ADR-75b's sum over all
// elements). The price is that the scatter itself stays serial.
//
// The DOF_Group loops (nodal mass / damping into A, nodal unbalance into b) are
// untouched and keep their serial position (ADR-75b 4.2 (ii)).

#ifndef LadrunoThreadedAssembly_h
#define LadrunoThreadedAssembly_h

class AnalysisModel;
class LinearSOE;
class Integrator;

enum LadrunoFELoopKind {
  LADRUNO_FE_TANGENT  = 0,   // getTangent  -> addA
  LADRUNO_FE_RESIDUAL = 1    // getResidual -> addB
};

// Runs the element loop threaded when every precondition holds and returns true
// (result then carries the loop's return code: 0, or failCode if any addA/addB
// failed, with the same warning text the serial loop prints, prefixed by
// `who`). Returns false when the caller must run its serial loop: 1 thread
// (the default), OpenMP compiled out, an MPI binary without the hybrid
// acknowledgement, the deep profiler armed, or fewer than two thread-safe FEs.
bool ladrunoThreadedFELoop(AnalysisModel *theModel, LinearSOE *theSOE,
                           Integrator *theIntegrator, int kind,
                           int failCode, const char *who, int &result);

#endif
