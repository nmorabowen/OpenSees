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
//
// Ladruno ADR-75 P1e (2026-07): `-krylov <digits>` — the OTHER reuse axis.
// The `factored` gate above can only skip work when A is UNCHANGED, so it pays
// nothing under full Newton, which is what most decks actually run. iparm[3]
// turns the retained L/U into a PRECONDITIONER for a tangent that has changed:
// phase 23 runs CGS (Intel's K=1, mtype 11) or CG (K=2, mtype 2 SPD only) and
// falls back to a real refactorization by itself if the iteration stalls.
// Two things make this more than a one-line iparm change:
//   * phase 23 is mandatory — the automatic direct fallback is documented for
//     phase 23 only, and the same failure under phase 33 is error = -4;
//   * a CGS WIN leaves the stored factors one tangent behind, so the phase-33
//     shortcut has to be forbidden afterwards (see `factorsCurrent`). Getting
//     that wrong is a silent wrong answer, not a crash.
// Off by default. Intel: "other values are only recommended for an advanced
// user."


#include <PARDISOGenLinSolver.h>
#include <PARDISOGenLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_service.h>              // Ladruno ADR-75 P1k: mkl_get_max_threads()
                                      // for the -stats "threads=" field. This TU
                                      // is only compiled into the build when MKL
                                      // was found (see the pardiso/CMakeLists.txt
                                      // guard) — unlike ProfilerRunMeta.h's use
                                      // of the same call, no _PARDISO ifdef is
                                      // needed here.
#include <profiler/ProfilerMacros.h>  // Ladruno ADR-75: phase-split brackets
                                      // (UmfPack parity, ADR-40 rank 8/10)

PARDISOGenLinSolver::PARDISOGenLinSolver()
:LinearSOESolver(SOLVER_TAGS_PARDISOGenLinSolver),
 theSOE(0), mtype(11), init(false), needsSymbolic(false), cachedN(0),
 reportStats(0),
 krylovL(0), krylovK(0), haveFactors(false), factorsCurrent(false),
 cgsCalls(0), cgsWins(0), cgsAdviceDone(false)
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


// ---- Ladruno ADR-75 P1d: PERTURBED PIVOTS ARE NOT AN ERROR -----------------
// iparm[9] tells PARDISO to REPLACE any pivot below eps*||A|| rather than fail,
// so a near-singular tangent comes back with error == 0 and a solution to a
// matrix that is NOT the one we assembled. The symmetric path runs eps = 1e-8
// (Intel's documented recommendation for mtype +-2) — five orders looser than
// the unsymmetric 1e-13 — so choosing -matrixType 2 materially widens this
// window, and iparm[7]=2 iterative refinement then hides small perturbations.
//
// Silent-wrong-answer class: without this the only symptom is Newton
// convergence quietly degrading, which gets blamed on the model. This fork
// already treats perturbed pivots as report-worthy in the FEAST `-certify`
// Sturm counts (ADR-43 P2) — same hazard, same answer. Warned once per
// factorization rather than refused: at a limit point a perturbed pivot may be
// exactly what lets a run continue, and that is the author's call to make, not
// ours. (Found by adversarial review.)
//
// Ladruno ADR-75 P1e: hoisted out of solve() — the CGS path (phase 23) can also
// end in a factorization, and that one needs the same check.
static void
ops_pardiso_perturbed(const int *iparm, int mtype)
{
	if (iparm[13] > 0)
		opserr << "WARNING PARDISOGenLinSolver: PARDISO perturbed "
		       << iparm[13] << " pivot(s) during factorization (mtype "
		       << mtype << ", threshold 1e-" << iparm[9]
		       << "). The solve is to a PERTURBED matrix — treat a slow or "
		          "stalling Newton here as a near-singular tangent, not a "
		          "solver hiccup.\n";
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
		/* nnz-in-factors (iparm 18 Fortran / iparm[17] here) and Mflops-of-
		   factorization (iparm 19 Fortran / iparm[18] here) reports: both are
		   IN/OUT controls — a negative value requested BEFORE the call is what
		   makes PARDISO fill the same slot with the real count afterwards.
		   Intel documents the Mflops report as costing extra analysis time, so
		   -- Ladruno ADR-75 P1k -- both are left at the safe 0 (disabled)
		   unless `-stats`/`-pardisoStats` asked for them, matching the MUMPS
		   `-stats` rule of paying that cost only when asked. */
		iparm[17] = reportStats ? -1 : 0;
		iparm[18] = reportStats ? -1 : 0;
		iparm[34] =  0;  /* ONE-based indexing — the SOE builds Fortran-style CSR */

		int phase = 11;
		// Ladruno ADR-75: same bracket names as UmfpackGenLinSolver, so a
		// profile comparing the two solvers lines the phases up column-for-
		// column. This one is the METIS reorder — once per sparsity pattern.
		{ OPS_PROFILE_SCOPE("soe.symbolic");
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		}
		if (error != 0) {
			ops_pardiso_report("symbolic factorization", error, mtype);
			return -1;
		}

		init = true;
		needsSymbolic = false;
		cachedN = n;                // for the destructor; see the header note
		theSOE->factored = false;   // a new pattern always owes a numeric pass

		// Ladruno ADR-75 P1e: a new pattern discards the factors, so there is no
		// preconditioner to hand CGS until the next phase 22 has run.
		haveFactors = false;
		factorsCurrent = false;
		// The CGS mode follows mtype, which is only known here. Intel documents
		// K=1 for nonsymmetric/structurally symmetric and K=2 for symmetric
		// POSITIVE DEFINITE — there is NO documented K for mtype -2, so the
		// symmetric-indefinite path (which is the right one for a softening or
		// buckling tangent) simply cannot use this lever.
		if (krylovL > 0) {
			krylovK = (mtype == 11) ? 1 : (mtype == 2 ? 2 : 0);
			if (krylovK == 0)
				opserr << "WARNING PARDISOGenLinSolver: -krylov is not available "
				          "for -matrixType 2 (mtype -2). Intel documents the CGS "
				          "preconditioner for\n     unsymmetric (mtype 11) and "
				          "symmetric POSITIVE DEFINITE (mtype 2) only; a "
				          "symmetric-indefinite LDL^T has no\n     documented CGS "
				          "mode. Continuing with direct factorization.\n";
		}
	}

	// ---- Ladruno ADR-75 P1e: CGS preconditioned by the RETAINED factors ----
	// The `factored` flag can only skip work when A is UNCHANGED. iparm[3] is
	// the complementary lever: it reuses the previous L/U as a preconditioner
	// for a tangent that HAS changed, which is the full-Newton case the flag
	// cannot touch. Phase 23 is mandatory here — Intel documents the automatic
	// direct fallback ("the factorization for a given A is automatically
	// recomputed in cases where the Krylov Subspace iteration failed") for
	// phase 23 ONLY; the same failure under phase 33 is just error = -4.
	//
	// Entered when the retained factors are unusable as a *direct* answer,
	// which is either of:
	//   theSOE->factored == false  - A was re-assembled since the last solve
	//   factorsCurrent   == false  - the last solve was a CGS win, so the
	//                                stored L/U belong to an OLDER A
	// The second is the subtle one: without it, a second solve() against the
	// same A (a second RHS, a recorder-driven re-solve) would take the phase-33
	// shortcut and silently answer with the previous tangent.
	bool solvedByKrylov = false;
	// Ladruno ADR-75 P1k: true only when THIS call ran phase 22 (a real numeric
	// factorization, first-time or a refactorization) — the -stats block below
	// is gated on this, not on a "have we ever printed for this pattern" latch,
	// so it fires once per factorization event, matching MUMPS `-stats` (which
	// prints inside `if (theMumpsSOE->factored == false)`, i.e. every job=5).
	bool didFactorNow = false;

	if (krylovK != 0 && haveFactors == true &&
	    (theSOE->factored == false || factorsCurrent == false)) {

		iparm[3] = 10 * krylovL + krylovK;
		int phase = 23;
		// Ladruno ADR-75: its own bracket, NOT soe.trisolve — the whole point
		// of -krylov is that this call replaces a factor+trisolve pair, and a
		// profile has to show which of the two regimes carried the run. NB a
		// CGS give-up refactorizes INSIDE this call (Intel's automatic
		// fallback), so that cost is billed here, not to soe.factor; iparm[19]
		// (the -stats win/fallback tally) says how often that happened.
		{ OPS_PROFILE_SCOPE("soe.cgs");
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		}
		iparm[3] = 0;   // leave the control array in its direct-solve state

		if (error != 0) {
			// The handle's factors are in an indeterminate state after a failed
			// phase 23 (PARDISO may have been part-way through the fallback
			// refactorization). Forbid the phase-33 shortcut on any subsequent
			// call rather than trusting them. (Adversarial review.)
			factorsCurrent = false;
			ops_pardiso_report("CGS solve (phase 23)", error, mtype);
			return -2;
		}

		cgsCalls++;
		solvedByKrylov = true;
		theSOE->factored = true;

		if (iparm[19] > 0) {
			// CGS answered it. The factors were NOT recomputed and are now one
			// tangent behind — see the header note on factorsCurrent.
			cgsWins++;
			factorsCurrent = false;
			if (reportStats)
				opserr << "PARDISO -krylov: CGS converged in " << iparm[19]
				       << " iteration(s); factorization reused\n";
		} else {
			// PARDISO gave up and refactored, so the handle now holds L/U of
			// THIS A. iparm[19] = -it_cgs*10 - cgs_error.
			factorsCurrent = true;
			const int cgsErr = -iparm[19] % 10;
			const int cgsIts = -iparm[19] / 10;
			// The running tally rides the FALLBACK line specifically: a fallback
			// is the moment the win rate matters, and it is the number that
			// decides whether -krylov is earning its place on this model.
			if (reportStats)
				opserr << "PARDISO -krylov: CGS gave up after " << cgsIts
				       << " iteration(s) (cgs_error " << cgsErr
				       << "); PARDISO refactored  [" << cgsWins << "/"
				       << cgsCalls << " CGS wins so far]\n";
			// cgs_error 5 is PARDISO telling us this flag is counterproductive
			// on this matrix — surface it even without -stats, once, because
			// the only other symptom is a run that is quietly SLOWER than the
			// default. Diagnostic, not a refusal: on a nonlinear path the
			// verdict can differ step to step.
			if (cgsErr == 5 && cgsAdviceDone == false) {
				cgsAdviceDone = true;
				opserr << "WARNING PARDISOGenLinSolver: PARDISO reports "
				          "cgs_error 5 — factorization on this matrix is fast "
				          "enough that\n     CGS preconditioning is not worth "
				          "it. Consider dropping -krylov (this is Intel's own "
				          "advice for\n     iparm[19] = -...5). Reported once; "
				          "later occurrences are silent.\n";
			}
			ops_pardiso_perturbed(iparm, mtype);
		}
	}

	// ---- numeric: ONLY when A changed — this is the reuse win --------------
	else if (theSOE->factored == false) {
		int phase = 22;
		{ OPS_PROFILE_SCOPE("soe.factor");   // Ladruno ADR-75 (UmfPack parity)
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		}
		if (error != 0) {
			ops_pardiso_report("numerical factorization", error, mtype);
			return -2;
		}
		theSOE->factored = true;
		haveFactors = true;      // Ladruno ADR-75 P1e: CGS now has a
		factorsCurrent = true;   // preconditioner, and it matches A
		didFactorNow = true;     // Ladruno ADR-75 P1k: -stats fires below
		ops_pardiso_perturbed(iparm, mtype);
	}

	// ---- triangular solve + iterative refinement: every call ---------------
	// Skipped after a phase-23 call, which already produced X.
	if (solvedByKrylov == false) {
		int phase = 33;
		{ OPS_PROFILE_SCOPE("soe.trisolve");   // Ladruno ADR-75 (UmfPack parity)
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		}
		if (error != 0) {
			ops_pardiso_report("solution", error, mtype);
			return -3;
		}
	}

	// ---- Ladruno ADR-75 P1k: `-stats`/`-pardisoStats`, EVERY numeric factor -
	// Printed once per phase-22 call (`didFactorNow`), i.e. the first
	// factorization AND every refactorization — not latched to "once per
	// sparsity pattern" the way P1d originally had it: a Newton run
	// refactorizes the SAME pattern repeatedly, and TIMs PM-01 D26 wants a
	// factorization-memory/fill number in every desktop leg's log, matching
	// what the shipped MUMPS `-stats` already does for every job=5 call (see
	// MumpsParallelSolver.cpp — printed inside `if (factored == false)`, the
	// same "every real factorization, not every solve" rule).
	//
	// Deliberately read AFTER phase 33 (not right after phase 22): Intel
	// documents iparm[16] (Fortran iparm(17)) as the peak over the numerical
	// factorization *and solution* phases, so reading it before the first
	// solve would report a LOWER BOUND.
	//
	// Labels use the Fortran 1-based iparm() numbering (iparm(N) == this
	// array's iparm[N-1]) so the printed numbers are checkable directly
	// against the MKL Developer Guide's PARDISO iparm table:
	//   iparm(15) = iparm[14]  peak memory (KB) during symbolic factorization
	//   iparm(16) = iparm[15]  permanent memory (KB) kept after phase 11
	//   iparm(17) = iparm[16]  memory (KB) for numerical factorization + solve
	//   iparm(18) = iparm[17]  nonzeros in the factors (L+U) -- reported only
	//                          because iparm[17] was set to -1 before phase 11
	//   iparm(19) = iparm[18]  Mflops of factorization -- likewise only
	//                          because iparm[18] was set to -1 before phase 11
	// Both negatives cost extra analysis time, which is exactly why they are
	// set only when `-stats` is on (see the symbolic-phase block above) --
	// byte-identical output otherwise.
	if (reportStats && didFactorNow) {
		opserr << "PARDISO stats: n=" << n << " nnz(A)=" << theSOE->nnz
		       << " matrixType=" << mtype
		       << " threads=" << mkl_get_max_threads() << "\n";
		opserr << "  factor entries iparm(18)  = " << iparm[17] << "\n";
		opserr << "  peak memory KB iparm(15)  = " << iparm[14] << "\n";
		opserr << "  perm memory KB iparm(16)  = " << iparm[15] << "\n";
		opserr << "  fact memory KB iparm(17)  = " << iparm[16] << "\n";
		opserr << "  factor Mflops  iparm(19)  = " << iparm[18] << "\n";
		// Supplementary, not part of the MUMPS-matching block above: the same
		// perturbed-pivot/refinement-step diagnostic P1d already reported.
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
    // state (mtype is derived from the SOE), so re-pointing it at a DIFFERENT
    // SOE has to invalidate that state. Not reachable through today's command
    // paths — each factory pairs one solver with one SOE for life — but the
    // class no longer tolerates the assumption being broken silently.
    needsSymbolic = true;
    // Ladruno ADR-75 P1e: the retained factors belong to the OLD SOE's matrix.
    haveFactors = false;
    factorsCurrent = false;
    return 0;
}


// Ladruno ADR-75 P1d
void
PARDISOGenLinSolver::setStats(int on)
{
    reportStats = on;
}


// Ladruno ADR-75 P1e
void
PARDISOGenLinSolver::setKrylov(int digits)
{
    // Intel's L. eps_CGS = 10^-L, and iparm[3] = 10*L + K must stay >= 0, so a
    // negative L is not merely useless — it would encode a different K.
    if (digits < 0) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - negative digits ("
               << digits << ") is meaningless; -krylov disabled\n";
        digits = 0;
    }
    // Bound L at 9. NOT an encoding limit — 10*L+K decodes unambiguously for
    // any L (K is the units digit, K <= 2), so iparm[3] = 101 would be a
    // perfectly legal L=10 (an earlier version of this comment claimed
    // otherwise; adversarial review). The real reason is behavioural: PARDISO
    // fixes a 150-iteration ceiling, so a criterion tighter than the
    // preconditioned iteration can actually reach does not buy accuracy — it
    // buys a guaranteed trip down the failure path, i.e. wasted iterations
    // followed by the factorization we were trying to avoid.
    if (digits > 9) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - digits (" << digits
               << ") is tighter than the CGS iteration can deliver and would "
                  "just force fallbacks; clamped to 9\n";
        digits = 9;
    }
    krylovL = digits;
    // Disabling must also clear K. K is otherwise only assigned in the symbolic
    // phase under `krylovL > 0`, so a stale K=1 left behind by an earlier
    // enable would keep the CGS branch live with iparm[3] = 10*0 + 1 = 1, i.e.
    // eps_CGS = 10^0 = 1 — an "accept almost anything" tolerance. Not reachable
    // through today's construct-once factories, but the class no longer relies
    // on that. (Adversarial review.)
    if (krylovL == 0)
        krylovK = 0;
    // Otherwise krylovK stays 0 until the symbolic phase derives mtype — that
    // is where the "no CGS mode for mtype -2" refusal lives.
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
