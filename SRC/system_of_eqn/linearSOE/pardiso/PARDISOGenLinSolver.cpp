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
//
// Ladruno ADR-75 P1d (2026-07): symmetric factorization. `mtype` was hardcoded
// 11; it is now DERIVED from the SOE's matType (see PARDISOGenLinSOE.h) so the
// upper-triangle storage and the factorization mode can never disagree, and the
// pivoting/scaling iparm entries branch on it.


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
 theSOE(0), mtype(11), init(false), needsSymbolic(false), cachedN(0),
 reportStats(0), statsDone(false)
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
ops_pardiso_report(const char *whatPhase, int error, int mtype)
{
	opserr << "WARNING PARDISOGenLinSolver::solve() - error " << error
	       << " during " << whatPhase << " (mtype " << mtype << "): ";
	switch (error) {
	case  -1: opserr << "input inconsistent\n"; break;
	case  -2: opserr << "not enough memory\n"; break;
	case  -3: opserr << "reordering problem\n"; break;
	case  -4:
		opserr << "zero pivot / singular matrix — check your model\n";
		// Ladruno ADR-75 P1d: by far the likeliest cause of -4 under mtype 2.
		if (mtype == 2)
			opserr << "     NOTE: -matrixType 1 asserts the tangent is POSITIVE "
			          "DEFINITE. A softening,\n     buckling or otherwise "
			          "indefinite tangent fails here; use -matrixType 2 "
			          "(mtype -2).\n";
		break;
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

	// Ladruno ADR-75 P1d (adversarial review): a size-0 SOE (a fully constrained
	// model, or solve() reached before a successful setSize) would otherwise
	// drive phase 11 with n=0 and NULL ia/ja/a.
	if (n <= 0 || ia == 0 || ja == 0 || a == 0) {
		opserr << "WARNING PARDISOGenLinSolver::solve() - the SOE has no "
		          "equations (n=" << n << "); nothing to factor\n";
		return -1;
	}

	int maxfct = 1, mnum = 1, msglvl = 0, nrhs = 1, error = 0;
	// Initialized, not just declared: they are passed to PARDISO as the unused
	// perm/rhs dummies, and MSVC /RTCu traps on reading an uninitialized local.
	double ddum = 0.0; int idum = 0;

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

		// Ladruno ADR-75 P1d: the factorization mode follows the STORAGE the SOE
		// actually built — full CSR => 11, upper-triangle CSR => 2 / -2. Reading
		// it here (rather than accepting it through a setter) makes a
		// storage/mtype mismatch unrepresentable. mtype persists for the
		// destructor's phase -1.
		const int soeMatType = theSOE->getMatType();
		mtype = (soeMatType == 1) ? 2 : (soeMatType == 2 ? -2 : 11);
		const bool symmetric = (mtype != 11);

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
		/* Pivoting/scaling differ by mtype — these are Intel's documented
		   per-mtype DEFAULTS, not a shared default:
		     unsymmetric (11): eps 1e-13 + MPS scaling + weighted matching,
		                       which is what keeps badly conditioned unsymmetric
		                       tangents factorizable;
		     symmetric (±2):   eps 1e-8 with Bunch-Kaufman 1x1/2x2 pivoting,
		                       scaling and matching OFF.
		   CORRECTION (adversarial review): an earlier version of this comment
		   claimed MKL applies scaling/matching as an UNSYMMETRIC permutation for
		   ±2 and that enabling them is the classic "symmetric PARDISO returns
		   garbage" cause. That is NOT Intel's position — Intel explicitly
		   supports iparm[10]=1 + iparm[12]=1 for mtype -2 and RECOMMENDS it for
		   highly indefinite symmetric systems, saddle-point structures in
		   particular (`constraints Lagrange` produces exactly that). The values
		   below are still the right conservative default, but the reason is
		   "Intel's documented default", not "the alternative is broken". Since
		   iparm is hardcoded and mtype is deliberately un-settable, there is
		   currently no way to opt into Intel's own mitigation for the model
		   class that most needs it — a known gap, not a defect. */
		iparm[9]  = symmetric ?  8 : 13;
		iparm[10] = symmetric ?  0 :  1;
		iparm[12] = symmetric ?  0 :  1;
		if (symmetric)
			iparm[20] = 1;  /* Bunch-Kaufman pivoting (required for indefinite
			                   tangents: a softening/buckling structure has
			                   negative eigenvalues, so mtype -2 — NOT 2 — is the
			                   safe symmetric choice, see the -matrixType docs) */
		/* nnz-in-factors report: -1 asks PARDISO to fill iparm[17] during the
		   reorder. Off by default (msglvl=0 means it is never printed anyway);
		   -stats turns it on. Unlike iparm[18], Intel does NOT document this
		   one as slowing the reordering. */
		iparm[17] = reportStats ? -1 : 0;
		iparm[18] =  0;  /* no Mflops report — Intel documents -1 as INCREASING
		                    reordering time, and with msglvl=0 we never read it
		                    (the prototype paid that cost for nothing) */
		iparm[34] =  0;  /* ONE-based indexing — the SOE builds Fortran-style CSR */

		int phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ops_pardiso_report("symbolic factorization", error, mtype);
			return -1;
		}

		init = true;
		needsSymbolic = false;
		statsDone = false;          // new pattern => report its memory once
		cachedN = n;                // for the destructor; see the header note
		theSOE->factored = false;   // a new pattern always owes a numeric pass
	}

	// ---- numeric: ONLY when A changed — this is the reuse win --------------
	if (theSOE->factored == false) {
		int phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			ops_pardiso_report("numerical factorization", error, mtype);
			return -2;
		}
		theSOE->factored = true;

		// ---- Ladruno ADR-75 P1d: PERTURBED PIVOTS ARE NOT AN ERROR ---------
		// iparm[9] tells PARDISO to REPLACE any pivot below eps*||A|| rather
		// than fail, so a near-singular tangent comes back with error == 0 and
		// a solution to a matrix that is NOT the one we assembled. The
		// symmetric path runs eps = 1e-8 (Intel's documented recommendation for
		// mtype +-2) — five orders looser than the unsymmetric 1e-13 — so
		// choosing -matrixType 2 materially widens this window, and iparm[7]=2
		// iterative refinement then hides small perturbations.
		//
		// Silent-wrong-answer class: without this the only symptom is Newton
		// convergence quietly degrading, which gets blamed on the model. This
		// fork already treats perturbed pivots as report-worthy in the FEAST
		// `-certify` Sturm counts (ADR-43 P2) — same hazard, same answer.
		// Warned once per factorization rather than refused: at a limit point a
		// perturbed pivot may be exactly what lets a run continue, and that is
		// the author's call to make, not ours. (Found by adversarial review.)
		if (iparm[13] > 0)
			opserr << "WARNING PARDISOGenLinSolver: PARDISO perturbed "
			       << iparm[13] << " pivot(s) during factorization (mtype "
			       << mtype << ", threshold 1e-" << iparm[9]
			       << "). The solve is to a PERTURBED matrix — treat a slow or "
			          "stalling Newton here as a near-singular tangent, not a "
			          "solver hiccup.\n";
	}

	// ---- triangular solve + iterative refinement: every call ---------------
	{
		int phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		if (error != 0) {
			ops_pardiso_report("solution", error, mtype);
			return -3;
		}
	}

	// ---- Ladruno ADR-75 P1d: `-stats`, ONCE per sparsity pattern -----------
	// Deliberately AFTER phase 33, not after phase 22 (adversarial review):
	// Intel documents iparm[16] as the peak over the numerical factorization
	// *and solution* phases, so reading it before the first solve reports a
	// LOWER BOUND. The ratio between two configurations measured the same way
	// survives that error, but the absolute figure is billed as "the number
	// that decides whether a model fits" and so has to be the real one.
	//   iparm[14] peak KB during the symbolic phase
	//   iparm[15] permanent KB kept after the symbolic phase
	//   iparm[16] peak KB during numeric factorization + solve
	//   total peak = max(iparm[14], iparm[15] + iparm[16])
	// Counters are pattern-determined, so reprinting per Newton iteration would
	// be noise. Report peak AND factor nnz, never just nnz — the MUMPS BLR study
	// (P2b) found the analogous INFOG(21) barely moved while stored factors
	// shrank 21.8%.
	if (reportStats && statsDone == false) {
		statsDone = true;
		const double peakSym  = iparm[14] / 1024.0;
		const double permSym  = iparm[15] / 1024.0;
		const double peakFact = iparm[16] / 1024.0;
		const double totalMB  = (iparm[14] > iparm[15] + iparm[16])
		                        ? peakSym : (permSym + peakFact);
		opserr << "PARDISO stats: n=" << n << " nnz(A)=" << theSOE->nnz
		       << " mtype=" << mtype
		       << (mtype == 11 ? " (unsymmetric, full CSR)"
		                       : " (symmetric, upper-triangle CSR)") << "\n";
		opserr << "  peak symbolic  = " << peakSym  << " MB\n";
		opserr << "  permanent      = " << permSym  << " MB\n";
		opserr << "  peak numeric   = " << peakFact << " MB   (factor + solve)\n";
		opserr << "  TOTAL PEAK     = " << totalMB  << " MB   <- the fit/no-fit number\n";
		if (iparm[17] > 0)
			opserr << "  nnz in factors = " << iparm[17] << "\n";
		opserr << "  perturbed pivots = " << iparm[13]
		       << "   refinement steps = " << iparm[6] << "\n";
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
    // Ladruno ADR-75 P1d (adversarial review): the solver now carries per-SOE
    // state (mtype is derived from the SOE, statsDone latches per pattern), so
    // re-pointing it at a DIFFERENT SOE has to invalidate that state. Not
    // reachable through today's command paths — each factory pairs one solver
    // with one SOE for life — but the class no longer tolerates the assumption
    // being broken silently.
    needsSymbolic = true;
    statsDone = false;
    return 0;
}


// Ladruno ADR-75 P1d
void
PARDISOGenLinSolver::setStats(int on)
{
    reportStats = on;
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
