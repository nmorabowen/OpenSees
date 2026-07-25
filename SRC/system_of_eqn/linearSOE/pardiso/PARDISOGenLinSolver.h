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
};

#endif

