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

// Ladruno WP-194 -- mixed-precision MUMPS for MumpsParallelSolver
// (`system Mumps ... -precision single [-refine k]`). ADR:
// Ladruno_implementation/194_ladruno_mumps_mixed_precision_adr.md.
//
// WHAT IT IS. The distributed double-precision matrix the SOE assembled is
// copied to float and factored by SMUMPS (single precision). Each solve sends
// the double right-hand side through the float factor and returns a double
// correction. Assembly, the residual and the Newton update stay double: the
// nonlinear algorithm evaluates the exact residual every iteration, so the
// single-precision factor only changes the iteration path (an inexact Newton
// direction), as BLR and WP-188 factor reuse already do. `-refine k` adds k
// classical mixed-precision refinement steps inside the solve (residual
// b - A x in double against the assembled double matrix, correction through
// the float factor); default 0, because Newton refines.
//
// COLLECTIVE. Every method below that talks to MUMPS or reduces is collective
// over the solver's communicator, exactly as the DMUMPS calls it replaces.
//
// Compiled only when the build found single-precision MUMPS (libsmumps +
// smumps_c.h; CMake defines _LADRUNO_SMUMPS). Without it the option is refused
// at the `system` command (MumpsParallelSolver::singlePrecisionAvailable()).
#ifndef LadrunoMumpsMixed_h
#define LadrunoMumpsMixed_h

#ifdef _LADRUNO_SMUMPS

#include <mpi.h>
extern "C" {
#include <smumps_c.h>
}

struct LadrunoMumpsControls {
  int    icntl7;
  int    icntl14;
  int    icntl35;
  double cntl7;
  int    printStats;
};

class LadrunoMumpsMixed
{
 public:
  explicit LadrunoMumpsMixed(int refineSteps);
  ~LadrunoMumpsMixed();

  // JOB=-2 on the live instance, if any (collective). Next solve re-initializes.
  void terminate(void);
  // the sparsity pattern changed (SOE setSize): analyse again at the next solve
  void requestAnalysis(void) { needsAnalysis = true; }

  // Collective over comm. rowA/colA/A: this rank's local triplets, 0-based
  // (temporarily shifted to 1-based around the MUMPS calls, as the double path
  // does). B and X are read/written on rank 0 only. `factored` is the SOE flag:
  // false -> convert + factor (JOB=2) first; true -> reuse the kept float factor.
  int solve(MPI_Comm comm, int sym, int n, int nnz, int *rowA, int *colA,
	    const double *A, const double *B, double *X, bool &factored,
	    const LadrunoMumpsControls &c);

  int refineSteps(void) const { return refine; }

 private:
  int  initialize(MPI_Comm comm, int sym, const LadrunoMumpsControls &c);
  void setControls(const LadrunoMumpsControls &c);
  int  convert(int nnz, const double *A, MPI_Comm comm);
  int  report(const char *phase);
  int  floatSolve(int n);                       // JOB=3 on rhsF (rank 0)
  void residual(MPI_Comm comm, int sym, int n, int nnz, const int *rowA,
		const int *colA, const double *A, const double *B,
		const double *X);               // r = B - A X into rbuf (rank 0)

  SMUMPS_STRUC_C sid;
  bool init;
  bool needsAnalysis;
  int  refine;
  int  rank;

  float  *aF;   int aFsize;                     // float copy of the local A
  float  *rhsF; int rhsSize;                    // rank-0 float RHS / solution
  double *xw;   double *yw; double *rbuf; int wSize;  // refinement work (n)

  // `-stats` bookkeeping (wall clock, max over ranks), WP-190 labels
  double tAnalysis, tSolveSum;
  int    nSolveSum;
  bool   firstSolveAfterFactor;
};

#endif // _LADRUNO_SMUMPS
#endif
