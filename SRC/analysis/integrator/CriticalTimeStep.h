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
class Element;
class Matrix;
class Vector;

// How the element mass is lumped to the diagonal M of the per-element pencil.
enum class CTSLumping {
    RowSum,   // (0) sum each row onto the diagonal. Cheap, but non-conservative for
              //   the rotational DOFs of beams/shells (row sums there can be ~0).
    Diagonal, // (1) diagonal-of-consistent: M_ii directly. Strictly positive for a
              //   consistent mass, so it avoids the row-sum rotational pathology.
              //   NOTE: this is NOT mass-conserving (sum of M_ii need not equal the
              //   element mass), so it can be NON-conservative for translational DOFs.
    HRZ       // (2) Hinton-Rock-Zienkiewicz mass-conserving lump (Ladruno, ADR 35):
              //   per translational direction d, scale the consistent diagonal so it
              //   sums to the element mass; rotational DOFs take the mean of the
              //   translational scales. Positive AND mass-conserving AND rotation-aware.
              //   See LadrunoMassLumping.h (Ladruno::hrzLump).
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
    // Ladruno (ADR-73 P3): overlay-aware undrained pencil. When a
    // LadrunoPorousOverlay owns elements of the domain, their per-element
    // pencils are augmented by dK_e = Q_e S_e^-1 Q_e^T (the element-local
    // undrained condensation), so undamped_dt/damped_dt above are the
    // UNDRAINED (corrected) limits. governing_drained_dt is the governing
    // (undamped) element's skeleton-only value, kept so reports can print
    // both and the implied factor; +inf when unknown (no augmentation, or
    // the governing element lives on another MPI rank).
    bool   overlayAugmented;      // any element's pencil was augmented
    double governing_drained_dt;  // drained 2/omega_max of the governing element
    // Ladruno (ADR-73 P3b): dual-CFL advisory for the EXPLICIT fluid lane
    // (-fluidUpdate explicit). When any FU_EXPLICIT overlay reports a finite
    // forward-Euler diffusion bound Δt_diff = h_min²/(2·ndm·c_v,max),
    // explicitFluid is set and fluid_diffusion_dt carries the min over
    // overlays; the bound is then min-folded into BOTH damped_dt and
    // undamped_dt (if diffusion ever governs — absurd k̄ — the returned
    // advisory is honest; the governing TAGS stay the solid element's, the
    // report clarifies). At realistic k̄ the slack is ~7e3× (E7.4), so this
    // advises rather than governs. Bounds the -subcycle window: N·Δt ≤ Δt_diff.
    bool   explicitFluid;         // any FU_EXPLICIT overlay contributed
    double fluid_diffusion_dt;    // min Δt_diff over overlays; +inf if none
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

// ---- per-element kernel (reused by computeCriticalTimeStep AND the selective
//      mass-scaling integrator, ADR 36) ------------------------------------
//
// Lump one element's CONSISTENT mass `M` into the diagonal `mdiag` (caller-sized
// to M.noRows()) per `lumping`. For HRZ the per-DOF translation/rotation tags
// are built from the element's node layout (first ndm DOFs of each node are
// translational; the element mass is assumed ordered node-by-node — true for
// node-major element matrices). Mass-first / D8: the caller passes the M it
// already fetched; this copies it fully into mdiag before returning, so the
// caller may then fetch the stiffness without aliasing a shared static matrix.
void lumpElementMass(Element *ele, const Matrix &M, CTSLumping lumping,
                     double *mdiag);

// Largest generalized eigenvalue lambda_max (= omega_max^2) of K v = lambda M v
// for one element, with M the diagonal `mdiag` (length n). useTangent selects
// getTangentStiff() over getInitialStiff(). Returns < 0 on DOF mismatch or
// eigensolve failure. NOTE: this is the EIGENSOLVE half only — it does NOT
// account for an element's self-reported bound (getExplicitCriticalTimeStep);
// callers that need the authoritative per-element step must use
// elementCriticalDt (below), not this directly.
//   Kadd (Ladruno, ADR-73 P3): optional state-independent ADDITIVE stiffness
//   (n x n) folded into K before the eigensolve — the porous-overlay undrained
//   augmentation dK_e = Q_e S_e^-1 Q_e^T. Null (the default) preserves the
//   drained pencil bit-exactly for every pre-P3 call site.
double elementLambdaMax(Element *ele, bool useTangent,
                        const double *mdiag, int n,
                        const Matrix *Kadd = 0);

// Authoritative UNDAMPED per-element critical step, with the same precedence
// computeCriticalTimeStep applies: a positive self-reported bound
// (getExplicitCriticalTimeStep) wins over the eigensolve (some elements — e.g.
// bipenalty couplings — carry a stability limit their per-element pencil can't
// express, where elementLambdaMax would return ~0 and miss it). Otherwise
// returns 2/sqrt(lambda_max). Returns < 0 if neither is available. This is the
// accessor the selective mass-scaling integrator (CentralDifferenceSMS, ADR 36)
// must call to size its per-element scale factor.
//   Kadd (Ladruno, ADR-73 P3 seam, WIRED at P3b): same optional additive
//   stiffness as elementLambdaMax. Callers: computeCriticalTimeStep (P3) and
//   BOTH SMS builders (LadrunoMassScaling::buildMassScaling /
//   buildMassScalingConsistent, P3b §3b.4) — SMS sizing now prices the
//   UNDRAINED per-element pencil for overlay-owned elements. The closed-form
//   scale s = T² + 2Tc is unchanged and remains exact: dK_e is mass- and
//   state-independent, so λ_max still scales as 1/s under mass injection.
//   A positive self-report still wins UNCORRECTED (the scan advises loudly
//   for that combo). History: pre-P3b no caller passed Kadd and SMS +
//   LadrunoPorousOverlay was UNSUPPORTED (blanket warning, since retired) —
//   ADR-73 §12 P3/P3b entries + LEDGER_quirks row.
double elementCriticalDt(Element *ele, bool useTangent,
                         const double *mdiag, int n,
                         const Matrix *Kadd = 0);

// Total STIFFNESS-proportional Rayleigh coefficient of an element, for the
// damped stable-step estimate and the SMS damped sizing. The element damping is
// C = alphaM*M + betaK*K + betaK0*K_init + betaKc*K_commit: all three beta slots
// produce the same high-frequency stiffness-proportional damping force
// (identical at the initial state, where K == K0 == Kc; conservative under
// softening), so ALL of them shrink the explicit step and ALL must enter the
// damped closed form. Reading coefs(1) alone silently UNDER-scaled
// `rayleigh 0 0 betaKinit 0` models — the most common form in practice, chosen
// precisely to avoid the current tangent (review 2026-07-01). Each slot is
// clamped at 0 before summing (matches the old betaK<0 guard).
double stiffnessRayleighBeta(Element *ele);

// true iff every entry of v is finite. The explicit integrators' divergence
// circuit breakers MUST use this, not Vector::pNorm(0): pNorm's max-compare
// `(data > value) ? data : value` is false for NaN entries, so NaN is silently
// SKIPPED and only +/-Inf ever tripped the old `A_max != A_max` check.
bool vectorIsFinite(const Vector &v);

#endif
