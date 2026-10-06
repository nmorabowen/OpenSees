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
// The PARDISO phases are split so that work is reused across solves:
//   phase 11 (reordering, symbolic)  once per sparsity pattern (setSize)
//   phase 22 (numerical factorization) only when A has changed
//   phase 33 (solve, refinement)      every call
//   phase -1 (release)                once, in the destructor
// Algorithms that keep the tangent (ModifiedNewton, Initial, KrylovNewton)
// then pay one factorization per tangent instead of one per solve.
//
// Reference: Intel oneMKL Developer Reference, "pardiso" and "pardiso iparm
// Parameter".


#include <PARDISOGenLinSolver.h>
#include <PARDISOGenLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_service.h>

PARDISOGenLinSolver::PARDISOGenLinSolver()
:LinearSOESolver(SOLVER_TAGS_PARDISOGenLinSolver),
 theSOE(0), mtype(11), init(false), needsSymbolic(false), cachedN(0),
 reportStats(0),
 krylovL(0), krylovK(0), haveFactors(false), factorsCurrent(false),
 cgsCalls(0), cgsWins(0), cgsAdviceDone(false)
{
	for (int i = 0; i < 64; i++) {
		pt[i] = 0;
		iparm[i] = 0;
	}
}


PARDISOGenLinSolver::~PARDISOGenLinSolver()
{
	// Release PARDISO's memory once. theSOE is not touched here (see the
	// header); the release phase does not read the matrix.
	if (init == true) {
		int maxfct = 1, mnum = 1, phase = -1, error = 0, msglvl = 0, nrhs = 1;
		int n = cachedN;
		double ddum = 0.0; int idum = 0;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, &idum, &idum, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		init = false;
	}
}


static void
pardiso_report(const char *whatPhase, int error, int mtype)
{
	opserr << "WARNING PARDISOGenLinSolver::solve() - error " << error
	       << " during " << whatPhase << " (mtype " << mtype << "): ";
	switch (error) {
	case  -1: opserr << "input inconsistent\n"; break;
	case  -2: opserr << "not enough memory\n"; break;
	case  -3: opserr << "reordering problem\n"; break;
	case  -4:
		opserr << "zero pivot, numerical factorization or iterative "
		          "refinement problem\n";
		if (mtype == 2)
			opserr << "     -matrixType 1 assumes a positive definite matrix; "
			          "use -matrixType 2 for an indefinite\n     (softening "
			          "or buckling) tangent.\n";
		break;
	case  -5: opserr << "unclassified internal error\n"; break;
	case  -6: opserr << "reordering failed\n"; break;
	case  -7: opserr << "diagonal matrix is singular\n"; break;
	case  -8: opserr << "32-bit integer overflow problem\n"; break;
	case -10: opserr << "error opening OOC files\n"; break;
	default:  opserr << "see the oneMKL PARDISO error table\n"; break;
	}
}


// iparm[9] makes PARDISO replace pivots smaller than 10^-iparm[9]*||A||
// instead of failing, so a nearly singular matrix returns error 0 and the
// solution of a perturbed matrix. Report it once per factorization.
static void
pardiso_perturbed(const int *iparm, int mtype)
{
	if (iparm[13] > 0)
		opserr << "WARNING PARDISOGenLinSolver: PARDISO perturbed "
		       << iparm[13] << " pivot(s) during factorization (mtype "
		       << mtype << ", threshold 1e-" << iparm[9]
		       << "); the matrix is nearly singular\n";
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
	int    *ia   = theSOE->rowStartA;   // one-based CSR row pointers
	int    *ja   = theSOE->colA;        // one-based CSR column indices
	double *a    = theSOE->A;
	double *Xptr = theSOE->X;
	double *Bptr = theSOE->B;

	if (n <= 0 || ia == 0 || ja == 0 || a == 0) {
		opserr << "WARNING PARDISOGenLinSolver::solve() - the SOE has no "
		          "equations (n=" << n << "); nothing to factor\n";
		return -1;
	}

	int maxfct = 1, mnum = 1, msglvl = 0, nrhs = 1, error = 0;
	double ddum = 0.0; int idum = 0;

	// ---- reordering and symbolic factorization, once per pattern ----------
	if (needsSymbolic == true || init == false) {

		if (init == true) {   // discard the previous pattern first
			int phase = -1;
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja,
				&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
			init = false;
			for (int i = 0; i < 64; i++) pt[i] = 0;
		}

		// mtype follows the storage the SOE built, so the two cannot disagree
		const int soeMatType = theSOE->getMatType();
		mtype = (soeMatType == 1) ? 2 : (soeMatType == 2 ? -2 : 11);
		const bool symmetric = (mtype != 11);

		// iparm is rebuilt for every new pattern; any option stored in it
		// must be re-applied here.
		for (int i = 0; i < 64; i++) iparm[i] = 0;
		iparm[0]  =  1;  /* use the values below, not the solver defaults */
		iparm[1]  =  2;  /* nested dissection from METIS */
		iparm[3]  =  0;  /* no iterative-direct algorithm */
		iparm[4]  =  0;  /* no user permutation */
		iparm[5]  =  0;  /* solution in x, b unchanged */
		iparm[7]  =  2;  /* maximum iterative refinement steps */
		/* Pivoting and scaling use Intel's documented defaults per mtype:
		   unsymmetric: perturbation 1e-13, scaling and weighted matching on;
		   symmetric:   perturbation 1e-8, Bunch-Kaufman 1x1/2x2 pivoting,
		                scaling and matching off. */
		iparm[9]  = symmetric ?  8 : 13;
		iparm[10] = symmetric ?  0 :  1;
		iparm[12] = symmetric ?  0 :  1;
		if (symmetric)
			iparm[20] = 1;  /* Bunch-Kaufman pivoting for indefinite matrices */
		/* nonzeros in the factors and factorization Mflops are reported only
		   when requested with a negative value; the flop count costs extra
		   analysis time, so only with -stats. */
		iparm[17] = reportStats ? -1 : 0;
		iparm[18] = reportStats ? -1 : 0;
		iparm[34] =  0;  /* one-based indexing */

		int phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			pardiso_report("symbolic factorization", error, mtype);
			return -1;
		}

		init = true;
		needsSymbolic = false;
		cachedN = n;
		theSOE->factored = false;   // a new pattern always needs phase 22

		haveFactors = false;
		factorsCurrent = false;
		// Intel documents K=1 (CGS) for unsymmetric and K=2 (CG) for
		// symmetric positive definite matrices only; there is no mode for
		// mtype -2.
		if (krylovL > 0) {
			krylovK = (mtype == 11) ? 1 : (mtype == 2 ? 2 : 0);
			if (krylovK == 0)
				opserr << "WARNING PARDISOGenLinSolver: -krylov is not available "
				          "for -matrixType 2 (mtype -2); PARDISO documents the "
				          "preconditioned iteration for mtype 11 and 2 only. "
				          "Continuing with direct factorization.\n";
		}
	}

	// ---- preconditioned CGS with the retained factors (phase 23) ----------
	// Used when the stored factors are not a factorization of the current A:
	// either A was reassembled, or the last solve was a CGS success. Phase 23
	// is required: PARDISO refactors automatically when the iteration fails
	// only in phase 23.
	bool solvedByKrylov = false;
	bool didFactorNow = false;

	if (krylovK != 0 && haveFactors == true &&
	    (theSOE->factored == false || factorsCurrent == false)) {

		iparm[3] = 10 * krylovL + krylovK;
		int phase = 23;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		iparm[3] = 0;

		if (error != 0) {
			factorsCurrent = false;
			pardiso_report("CGS solve (phase 23)", error, mtype);
			return -2;
		}

		cgsCalls++;
		solvedByKrylov = true;
		theSOE->factored = true;

		if (iparm[19] > 0) {
			// converged; the stored factors are still those of an older A
			cgsWins++;
			factorsCurrent = false;
			if (reportStats)
				opserr << "PARDISO -krylov: CGS converged in " << iparm[19]
				       << " iteration(s); factorization reused\n";
		} else {
			// PARDISO refactored; iparm[19] = -it_cgs*10 - cgs_error
			factorsCurrent = true;
			const int cgsErr = -iparm[19] % 10;
			const int cgsIts = -iparm[19] / 10;
			if (reportStats)
				opserr << "PARDISO -krylov: CGS stopped after " << cgsIts
				       << " iteration(s) (cgs_error " << cgsErr
				       << "); PARDISO refactored  [" << cgsWins << "/"
				       << cgsCalls << " CGS successes so far]\n";
			// cgs_error 5: factorization is fast enough that CGS does not pay
			if (cgsErr == 5 && cgsAdviceDone == false) {
				cgsAdviceDone = true;
				opserr << "WARNING PARDISOGenLinSolver: PARDISO reports "
				          "cgs_error 5 (the factorization is fast enough that "
				          "CGS is not worthwhile); consider dropping -krylov. "
				          "Reported once.\n";
			}
			pardiso_perturbed(iparm, mtype);
		}
	}

	// ---- numerical factorization, only when A changed ---------------------
	else if (theSOE->factored == false) {
		int phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			pardiso_report("numerical factorization", error, mtype);
			return -2;
		}
		theSOE->factored = true;
		haveFactors = true;
		factorsCurrent = true;
		didFactorNow = true;
		pardiso_perturbed(iparm, mtype);
	}

	// ---- back substitution and iterative refinement, every call -----------
	if (solvedByKrylov == false) {
		int phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		if (error != 0) {
			pardiso_report("solution", error, mtype);
			return -3;
		}
	}

	// ---- -stats, after every numerical factorization ----------------------
	// Read after phase 33 because iparm(17) is the peak over factorization
	// and solution. Labels use Intel's one-based iparm() numbering.
	if (reportStats && didFactorNow) {
		opserr << "PARDISO stats: n=" << n << " nnz(A)=" << theSOE->nnz
		       << " matrixType=" << mtype
		       << " threads=" << mkl_get_max_threads() << "\n";
		opserr << "  factor entries iparm(18)  = " << iparm[17] << "\n";
		opserr << "  peak memory KB iparm(15)  = " << iparm[14] << "\n";
		opserr << "  perm memory KB iparm(16)  = " << iparm[15] << "\n";
		opserr << "  fact memory KB iparm(17)  = " << iparm[16] << "\n";
		opserr << "  factor Mflops  iparm(19)  = " << iparm[18] << "\n";
		opserr << "  perturbed pivots = " << iparm[13]
		       << "   refinement steps = " << iparm[6] << "\n";
	}

	return 0;
}


int
PARDISOGenLinSolver::setSize()
{
	// the pattern changed; phase 11 runs on the next solve()
	needsSymbolic = true;
	return 0;
}


int
PARDISOGenLinSolver::setLinearSOE(PARDISOGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    // state derived from the previous SOE no longer applies
    needsSymbolic = true;
    haveFactors = false;
    factorsCurrent = false;
    return 0;
}


void
PARDISOGenLinSolver::setStats(int on)
{
    reportStats = on;
}


void
PARDISOGenLinSolver::setKrylov(int digits)
{
    // iparm[3] = 10*L + K must stay non-negative
    if (digits < 0) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - negative digits ("
               << digits << "); -krylov disabled\n";
        digits = 0;
    }
    // PARDISO caps the iteration at 150 steps; a tolerance tighter than the
    // iteration can reach only forces the refactorization fallback.
    if (digits > 9) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - digits (" << digits
               << ") clamped to 9\n";
        digits = 9;
    }
    krylovL = digits;
    // K is otherwise set only at the symbolic phase; clear it when disabling
    if (krylovL == 0)
        krylovK = 0;
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
