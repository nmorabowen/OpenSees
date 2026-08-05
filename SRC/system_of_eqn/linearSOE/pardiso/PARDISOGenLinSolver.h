// File: ~/system_of_eqn/pardiso/PARDISOGenLinSolver.h
//
// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A
//
// Description: This file contains the class definition for 
// PARDISOLinSolver. It solves the Sparse General SOE by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) PARDISOGenLinSolver.h, revA"


#ifndef PARDISOGenLinSolver_h
#define PARDISOGenLinSolver_h

#include <LinearSOESolver.h>
#include <PARDISOGenLinSOE.h>

//https://software.intel.com/en-us/articles/pardiso-parameter-table

class PARDISOGenLinSolver : public LinearSOESolver
{
  public:
	  PARDISOGenLinSolver();
    ~PARDISOGenLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(PARDISOGenLinSOE &theSOE);

    // Ladruno ADR-75 P1d: `system Pardiso -stats` — dump PARDISO's own memory
    // counters once per sparsity pattern. Mirrors the shipped MUMPS `-stats`
    // (INFOG/RINFOG); without it the symmetric path's MEMORY claim cannot be
    // checked, and P1c measured memory — not time — as the binding constraint.
    void setStats(int on);

    // Ladruno ADR-75 P1e: `system Pardiso -krylov <digits>` — factorization-
    // preconditioned CGS (iparm[3]). `digits` is Intel's L: the iteration stops
    // at ||dx_i||/||dx_0|| < 10^-digits. 0 disables (the default). This is the
    // ONE reuse lever the `factored` flag cannot reach: it reuses a retained
    // L/U across a tangent that HAS changed, i.e. it pays under full Newton.
    void setKrylov(int digits);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag,
		 Channel &theChannel,
		 FEM_ObjectBroker &theBroker);
  protected:

  private:
	  PARDISOGenLinSOE *theSOE;

	  // Ladruno ADR-75 P1: factorization reuse. The stock prototype declared
	  // pt[64]/iparm as LOCALS of solve() and ran phases 11->22->33->-1 on EVERY
	  // call, so it re-did the METIS symbolic reorder + numeric factorization and
	  // then freed everything, every solve. The handle must persist across calls
	  // for reuse to be possible at all, so it lives here now.
	  void *pt[64];        // PARDISO internal handle — OPAQUE, never touch/copy
	  int   iparm[64];     // control array (was leaked via `new` each solve)
	  // Ladruno ADR-75 P1d: `mtype` is DERIVED from the SOE's matType at the
	  // symbolic phase (11 unsym / 2 SPD / -2 symmetric indefinite) and then
	  // cached here, so the storage format and the factorization mode cannot
	  // disagree — there is deliberately no setter. It must survive into the
	  // destructor for the phase -1 release, which is why it is a member.
	  int   mtype;         // 11 = real unsymmetric (matches the SOE's full CSR)
	  bool  init;          // true once phase 11 has run (=> phase -1 owed)
	  bool  needsSymbolic; // sparsity changed => redo the reorder
	  // Ladruno ADR-75 P1: the destructor MUST NOT dereference theSOE.
	  // ~LinearSOE() deletes theSolver, so by the time ~PARDISOGenLinSolver runs
	  // the SOE's derived destructor has ALREADY freed A/B/X/rowStartA/colA and
	  // the object is partially destroyed — reading it is use-after-free. Cache
	  // the order here and hand PARDISO dummy arrays for the release phase.
	  int   cachedN;       // matrix order, captured at the symbolic phase
	  int   reportStats;   // Ladruno ADR-75 P1d: -stats requested
	  bool  statsDone;     // ...and already printed for this pattern

	  // ---- Ladruno ADR-75 P1e: factorization-preconditioned CGS ------------
	  int   krylovL;       // -krylov <digits>; 0 = off (eps_CGS = 10^-krylovL)
	  int   krylovK;       // Intel's K: 1 = LU-CGS (mtype 11), 2 = LL^T-CG
	                       // (mtype 2 ONLY — K=2 is documented for symmetric
	                       // POSITIVE DEFINITE). 0 = disabled for this mtype.
	  // haveFactors/factorsCurrent are NOT the same question, and conflating
	  // them is the silent-wrong-answer in this feature:
	  //   haveFactors    - some L/U exists in the handle (a preconditioner is
	  //                    available at all).
	  //   factorsCurrent - that L/U is of exactly the A now in the SOE.
	  // A CGS *win* returns a correct solution while leaving the STALE factors
	  // in place, so a later phase-33-only call would solve the previous
	  // tangent. factorsCurrent is what forbids that shortcut.
	  bool  haveFactors;
	  bool  factorsCurrent;
	  int   cgsCalls;      // phase-23 calls made
	  int   cgsWins;       // ...that were answered by CGS (no refactorization)
	  bool  cgsAdviceDone; // cgs_error 5 ("turn this off") reported once
};

#endif

