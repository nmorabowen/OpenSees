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
// What: "@(#) PARDISOGenLinSolver.cpp, revA"
//
// Ladruno ADR-75 P1 (2026-07): factorization reuse. As contributed, solve() ran
// the full PARDISO cycle 11 (reorder) -> 22 (factor) -> 33 (solve) -> -1 (free)
// on EVERY call, with pt[]/iparm as function locals — so the METIS symbolic
// reorder and the numeric factorization were repeated for every single solve and
// then thrown away, and iparm was leaked (`new int[64]`, never deleted). The
// handle now persists as a member and the three phases are split:
//   * phase 11 once per sparsity pattern  (driven by setSize())
//   * phase 22 only when A changed        (driven by the SOE `factored` flag,
//                                          mirroring MumpsSolver's job=5/job=3)
//   * phase 33 every call
//   * phase -1 once, in the destructor
// This is what makes the solver pay off under tangent-reusing algorithms
// (ModifiedNewton / Initial / Krylov / IMPL-EX).


#include <PARDISOGenLinSolver.h>
#include <PARDISOGenLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>

PARDISOGenLinSolver::PARDISOGenLinSolver()
:LinearSOESolver(SOLVER_TAGS_PARDISOGenLinSolver),
 theSOE(0), mtype(11), init(false), needsSymbolic(false), cachedN(0)
{
	// Ladruno ADR-75 P1: the handle and control array are members now (see the
	// header note). PARDISO REQUIRES pt[] to be zeroed before the first call.
	for (int i = 0; i < 64; i++) {
		pt[i] = 0;
		iparm[i] = 0;
	}
}


PARDISOGenLinSolver::~PARDISOGenLinSolver()
{
	// Ladruno ADR-75 P1: release PARDISO's internal memory exactly ONCE, here,
	// instead of after every solve. Skipped when no factorization was ever set
	// up (calling phase -1 on a virgin handle is not meaningful).
	//
	// CRITICAL (found in adversarial review): do NOT touch theSOE here.
	// ~LinearSOE() deletes theSolver, so ~PARDISOGenLinSOE() has already run and
	// freed rowStartA/colA/A/B/X — dereferencing them would be use-after-free.
	// The release phase does not read the user's matrix, so cachedN + dummies
	// are sufficient (and safe).
	if (init == true) {
		int maxfct = 1, mnum = 1, phase = -1, error = 0, msglvl = 0, nrhs = 1;
		int n = cachedN;
		double ddum; int idum = 0;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, &idum, &idum, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		init = false;
	}
}


// Ladruno ADR-75 P1: decode the common MKL PARDISO error codes (the prototype
// printed a bare integer).
static void
ops_pardiso_report(const char *whatPhase, int error)
{
	opserr << "WARNING PARDISOGenLinSolver::solve() - error " << error
	       << " during " << whatPhase << ": ";
	switch (error) {
	case  -1: opserr << "input inconsistent\n"; break;
	case  -2: opserr << "not enough memory\n"; break;
	case  -3: opserr << "reordering problem\n"; break;
	case  -4: opserr << "zero pivot / singular matrix — check your model\n"; break;
	case  -5: opserr << "unclassified internal error\n"; break;
	case  -6: opserr << "reordering failed\n"; break;
	case  -7: opserr << "diagonal matrix is singular\n"; break;
	case  -8: opserr << "32-bit integer overflow\n"; break;
	case -10: opserr << "could not open the OOC file\n"; break;
	default:  opserr << "see the MKL PARDISO error table\n"; break;
	}
}


int
PARDISOGenLinSolver::solve(void)
{
	if (theSOE == 0) {
		opserr << "WARNING PARDISOLinSolver::solve(void)- ";
		opserr << " No LinearSOE object has been set\n";
		return -1;
	}

	int     n    = theSOE->size;
	int    *ia   = theSOE->rowStartA;   // 1-based CSR row pointers (see the SOE)
	int    *ja   = theSOE->colA;        // 1-based CSR column indices
	double *a    = theSOE->A;
	double *Xptr = theSOE->X;
	double *Bptr = theSOE->B;

	int maxfct = 1, mnum = 1, msglvl = 0, nrhs = 1, error = 0;
	double ddum; int idum;

	// ---- symbolic: ONCE per sparsity pattern (Ladruno ADR-75 P1) ----------
	// `init == false` is part of the guard deliberately (adversarial review):
	// if solve() is ever reached without setSize() having run, needsSymbolic is
	// still false and we would otherwise drive phase 22 on a virgin handle with
	// an all-zero iparm. Re-analyzing is the safe fallback.
	if (needsSymbolic == true || init == false) {

		if (init == true) {   // discard the previous pattern's factors first
			int phase = -1;
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja,
				&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
			init = false;
			for (int i = 0; i < 64; i++) pt[i] = 0;
		}

		for (int i = 0; i < 64; i++) iparm[i] = 0;
		iparm[0]  =  1;  /* do not use the solver defaults; the values below apply */
		iparm[1]  =  2;  /* fill-reducing reordering from METIS */
		iparm[2]  =  0;  /* reserved in current MKL — the thread count comes from
		                    MKL_NUM_THREADS / mkl_set_num_threads, NOT from here
		                    (the prototype's iparm[2]=1 did not pin 1 thread) */
		iparm[3]  =  0;  /* no iterative-direct algorithm */
		iparm[4]  =  0;  /* no user fill-in reducing permutation */
		iparm[5]  =  0;  /* write the solution into x, leave b intact */
		iparm[7]  =  2;  /* max steps of iterative refinement */
		iparm[9]  = 13;  /* perturb pivots with 1e-13 */
		iparm[10] =  1;  /* nonsymmetric permutation and MPS scaling */
		iparm[12] =  1;  /* maximum weighted matching (MKL's default for mtype 11;
		                    matters on badly conditioned unsymmetric tangents) */
		iparm[17] =  0;  /* no nnz-in-factors report (msglvl=0: never printed) */
		iparm[18] =  0;  /* no Mflops report — Intel documents -1 as INCREASING
		                    reordering time, and with msglvl=0 we never read it
		                    (the prototype paid that cost for nothing) */
		iparm[34] =  0;  /* ONE-based indexing — the SOE builds Fortran-style CSR */

		int phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ops_pardiso_report("symbolic factorization", error);
			return -1;
		}

		init = true;
		needsSymbolic = false;
		cachedN = n;                // for the destructor; see the header note
		theSOE->factored = false;   // a new pattern always owes a numeric pass
	}

	// ---- numeric: ONLY when A changed — this is the reuse win --------------
	if (theSOE->factored == false) {
		int phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ops_pardiso_report("numerical factorization", error);
			return -2;
		}
		theSOE->factored = true;
	}

	// ---- triangular solve + iterative refinement: every call ---------------
	{
		int phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		if (error != 0) {
			ops_pardiso_report("solution", error);
			return -3;
		}
	}

	return 0;
}


int
PARDISOGenLinSolver::setSize()
{
	// Ladruno ADR-75 P1: the SOE re-derived the sparsity pattern, so the stored
	// symbolic factorization is stale. Defer the actual phase-11 call to solve()
	// (setSize can be invoked before A holds anything meaningful).
	needsSymbolic = true;
	return 0;
}


int
PARDISOGenLinSolver::setLinearSOE(PARDISOGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}


int
PARDISOGenLinSolver::sendSelf(int cTAg, Channel &theChannel)
{
    // doing nothing
    return 0;
}


int
PARDISOGenLinSolver::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}
