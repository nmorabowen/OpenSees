/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Written: nmb (UANDES), 05/2026
//
// Shared per-element critical-time-step machinery for the explicit Noh-Bathe
// integrators (ExplicitBathe, ExplicitBatheLNVD). Hoisted out of ExplicitBathe.cpp
// so the two integrators reference ONE definition instead of a hand-copied
// `extern` (which could silently drift out of sync). See the explicit-dynamics
// ADR, item 2 (header hoist) and item 6 (DSYGV + lumping robustness).
//
// The estimate is the governing (minimum over the local domain) critical time
// step from the per-element generalized eigenproblem  K v = lambda M v, reported
// as the CONSERVATIVE central-difference limit dt_cd = 2/omega_max. The Noh-Bathe
// scheme is stable up to ~EB_NB_STABILITY_FACTOR * dt_cd (see ExplicitBathe.h).

#ifndef CriticalTimeStep_h
#define CriticalTimeStep_h

class AnalysisModel;

// Sentinels returned by the integrators' getCriticalTimeStep() so a caller (the
// criticalTimeStep() Py/Tcl command) can tell *why* there is no value, instead
// of collapsing every "no value" case onto a single negative number (ADR item 4).
namespace CriticalTimeStep {
    const double DISABLED       = -1.0; // estimation not requested (no -cfl/-cflAbort)
    const double NOT_COMPUTED   = -2.0; // requested, but domainChanged() has not run yet
    const double NOT_APPLICABLE = -3.0; // computed, but no element produced a finite
                                        //   estimate (e.g. a pure nodal-mass model with
                                        //   no element mass+stiffness pencil)
}

// How the element mass is lumped to the diagonal M of the per-element pencil.
enum class CTSLumping {
    RowSum,  // sum each row onto the diagonal. Cheap, but non-conservative for the
             //   rotational DOFs of beams/shells (row sums there can be ~0).
    Diagonal // diagonal-of-consistent: M_ii directly. Dimensionally consistent and
             //   strictly positive for a consistent mass, so it avoids the row-sum
             //   rotational pathology (the "HRZ-style" option in the ADR).
};

// Result of one scan over the (local) domain's elements.
struct CTSResult {
    double damped_dt;      // conservative CD limit incl. element Rayleigh damping;
                           //   +inf if nothing contributed
    double undamped_dt;    // 2/omega_max; +inf if nothing contributed
    int    damped_tag;     // governing element tag (damped); -1 if none
    int    undamped_tag;   // governing element tag (undamped); -1 if none
    int    n_contributing; // elements that produced a finite, positive estimate
    int    n_scanned;      // elements iterated
};

// Scan all elements of theModel's domain and return the governing critical step.
//   useTangent : eigensolve uses getTangentStiff() (current state) instead of
//                getInitialStiff(). Lets -recompute track a stiffening/softening
//                model (contact, plasticity).
//   lumping    : how the element mass is lumped to the diagonal M.
//
// Numerics: the lumped M is diagonal and (usually) positive, so the symmetric-
// definite solver DSYGV is tried first; if M is not positive definite (zero-mass
// DOFs make the pencil indefinite) the routine falls back to the general DGGEV.
// A relative-beta threshold discards near-singular eigenpairs (spurious huge
// lambda from beta ~ 0). Under a parallel build the result is reduced (MPI_MIN)
// across ranks so every rank sees the global governing step.
//
// Mass is read and copied to the diagonal BEFORE the stiffness is fetched: many
// elements return a reference to a shared static matrix from getMass()/
// getInitialStiff()/getTangentStiff(), so holding both references at once would
// alias (the second call clobbers the first). See ADR decision D8.
CTSResult computeCriticalTimeStep(AnalysisModel *theModel,
                                  bool useTangent,
                                  CTSLumping lumping);

#endif
