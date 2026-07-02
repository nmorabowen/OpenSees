/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024
// Revision: C  (Ladruno ADR-52 W1-E2: collapsed the ExplicitBathe* family
//               6 -> 1 with orthogonal -lnvd / -sms / -consistent flags; nmb 06/2026)
//
// Description: This file contains the class definition for ExplicitBathe.
// ExplicitBathe is an algorithmic class for performing a transient analysis
// using the explicit Bathe (Noh-Bathe 2013) time integration scheme.
//
// This algorithm is similar to a central difference scheme but perturbed to
// introduce numerical damping. It is a second-order accurate explicit scheme.
// Unlike the CentralDifference class, this one only assembles the mass matrix
// on the right hand side, making it much easier to use with diagonal mass matrices.
//
// The time-step required for stability is approximately twice that of
// standard explicit central difference schemes, making it more efficient for
// many applications. The optional critical-time-step computation below reports
// the CONSERVATIVE central-difference limit dt_cd = 2/omega_max (a guaranteed-
// safe lower bound for this scheme) together with the Noh-Bathe limit
// EB_NB_STABILITY_FACTOR * dt_cd (the published ~2x advantage at the optimal p).
//
// ── Ladruno ADR-52 W1-E2: capability flags ──────────────────────────────────
// One class now hosts what used to be SIX sibling/subclass tags. The three
// capability axes are ORTHOGONAL and selected by flags on the integrator command:
//
//   -lnvd <alpha>     FLAC-style Local Non-Viscous Damping for dynamic relaxation
//                     / quasi-static solving. At every solve the assembled residual
//                     r is modified to r_i <- r_i - alpha*|r_i|*sign(v_i) (opposes
//                     local motion, vanishes as v->0, so it does NOT bias the
//                     converged static solution). alpha in [0,1), classic FLAC
//                     default 0.8. Applied symmetrically to BOTH Noh-Bathe sub-steps
//                     via the formUnbalance() override. (was ExplicitBatheLNVD.)
//
//   -sms <dtTarget>   Selective MASS SCALING (LS-DYNA DT2MS-style): at domainChanged()
//                     inflate the mass of the elements that throttle the timestep so a
//                     larger dtTarget is stable. LUMPED by default (additive nodal mass,
//                     no solve-path change — the RHS mass assembly sees it).
//                     (was ExplicitBatheSMS.)
//
//   -consistent       CONSISTENT (Olovsson) mass scaling on top of -sms: a centroidal
//                     scaling mass M_bar = beta[diag(m)-mm^T/M_e], solved matrix-free
//                     (PCG) at BOTH Noh-Bathe sub-steps via the refineAccel() hook.
//                     REQUIRES -sms and `system Diagonal`/`system MPIDiagonal`.
//                     (was ExplicitBatheSMSConsistent / *LNVDSMSConsistent.)
//
// The six historical command names + class tags are retained one release as
// deprecated aliases (each forces the corresponding fixed flag set); the flag
// combo deterministically selects the matching legacy classTag so saved-DB /
// parallel recvSelf of any retired tag still reconstructs. See
// Ladruno_implementation/52_ladruno_integrator_strengthening_adr.md (W1-E2).
//
// Recommended configuration (explicit recipe):
//   - lumped (diagonal) element mass, e.g. -lMass / "-mass" on elements
//   - system Diagonal      (trivial M^{-1}; consistent mass + sparse solver is
//                           wrong/pointless here)
//   - algorithm Linear     (the scheme needs exactly 2 solves/step)
//   - integrator ExplicitBathe 0.54
// Pair with `recorder EnergyBalance ... -file energy.txt` to watch for spurious
// energy growth (the instability signature of an over-large dt).
//
// Critical-time-step (-cfl / -tangent / -recompute / criticalTimeStep()) caveats:
//   - It is a per-element estimate using row-sum-lumped element mass. This is
//     exact/conservative for translational DOFs (bars, solids) but row-sum
//     lumping can be non-conservative for ROTATIONAL DOFs (beams/shells); treat
//     dt_cr as a guide there, not a guarantee.
//   - It ignores constraints (equalDOF/rigid diaphragms/MP) and pure nodal mass,
//     so a constrained or nodal-mass-only model may report a non-binding or
//     <=0 (criticalTimeStep() returns <=0) value.
//   - -recompute N / -tangent only UPDATE the reported dt_cr; they do NOT change
//     the analysis dt. For stiffening/contact, either pair with -cflAbort or
//     re-query criticalTimeStep() from the driver and sub-divide (see the
//     adaptive_substep example). It is O(N_elements) DGGEV per refresh -- use a
//     large N. NOTE: with -sms active, -cflAbort/-recompute are downgraded to
//     REPORT-ONLY (the un-augmented element pencil cannot see the scaling mass).
//
// Key features:
// - Second-order accurate
// - Conditionally stable (dt <= 2/omega_max central-difference reference)
// - Built-in numerical damping (controlled by parameter p, p=0.54 is a good choice)
// - Only mass matrix on RHS (no tangent matrix assembly required)
// - Optional automatic critical time step calculation
// - Optional FLAC local non-viscous damping (-lnvd) and selective mass scaling
//   (-sms [-consistent]) — orthogonal, composable
//
// Reference:
// Gunwoo Noh, Klaus-Jürgen Bathe,
// "An explicit time integration scheme for the analysis of wave propagations",
// Computers & Structures, Volume 129, 2013, Pages 178-193,
// ISSN 0045-7949,
// https://doi.org/10.1016/j.compstruc.2013.06.007.

#ifndef ExplicitBathe_h
#define ExplicitBathe_h

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>   // CTSLumping, CTSResult, computeCriticalTimeStep()
#include <LadrunoProjectionConsumer.h>   // Ladruno (ADR-30 P5): projection seam
#include <Vector.h>             // Ladruno (W1-E2): SMS injected-mass map value type
#include <map>                  // Ladruno (W1-E2): per-node injected diagonal ΔM
#include <vector>               // Ladruno (W1-E2): consistent-scaling blocks

// Published Noh-Bathe stability advantage over central difference (~2x at the
// optimal sub-step parameter). dt_cr_NB = EB_NB_STABILITY_FACTOR * (2/omega_max).
// Reference: Noh & Bathe (2013), Computers & Structures 129, 178-193.
#define EB_NB_STABILITY_FACTOR 2.0

class DOF_Group;
class FE_Element;
class Domain;                       // Ladruno (W1-E2): SMS nodal-mass injection target
class LadrunoConstraintProjector;   // Ladruno (ADR-30 P5)
namespace Ladruno { struct ConsistentBlock; }   // Ladruno (W1-E2): consistent SMS

// Ladruno (ADR-30 P5): ExplicitBathe is a LadrunoProjectionConsumer, so the
// LadrunoProjectionHandler can enforce MP/EQ constraints by M-orthogonal accel
// projection under the Noh-Bathe scheme (same contract as CentralDifferenceLadruno).
class ExplicitBathe : public TransientIntegrator, public LadrunoProjectionConsumer
{
public:
    // Broker default constructor (classTag 33000, all flags off).
    ExplicitBathe();

    // Ladruno (W1-E2): unified constructor. The plain-ExplicitBathe call (omitting the
    // trailing LNVD/SMS args) is BYTE-IDENTICAL to the historical class. The actual
    // integrator classTag is DERIVED from the {lnvd, sms, consistent} flag combo
    // (tagForFlags) so serialization/broker round-trips reconstruct the right legacy tag.
    ExplicitBathe(double p, int compute_critical_timestep_ = 0,
                  bool verbose = false, bool cflAbort = false,
                  double divergenceFactor = 0.0,
                  bool cflUseTangent = false, int cflRecomputeEvery = 0,
                  CTSLumping lumping = CTSLumping::RowSum,
                  // --- LNVD (-lnvd) ---
                  bool useLNVD = false, double alpha_flac = 0.0,
                  // --- selective mass scaling (-sms [-consistent]) ---
                  bool useSMS = false, bool useConsistent = false,
                  double dtTarget = 0.0, double maxAddedMassFrac = 0.05,
                  double pcgTol = 1.0e-10, int pcgMaxIt = 200);

    // Destructor
    ~ExplicitBathe();

    // Ladruno (W1-E2): broker factory — reconstruct a default object whose flag set is
    // decoded from a (possibly retired) integrator classTag, so getNewIntegrator() can
    // route all six legacy tags to this one class before recvSelf() fills the params.
    static ExplicitBathe *makeForBroker(int classTag);

    // Conservative (central-difference) critical time step, available after a
    // step has run with compute_critical_timestep enabled; <=0 if not computed.
    double getCriticalTimeStep(void) const override;

    // Ladruno (W1-E2, was ExplicitBatheLNVD): infinity-norm of the most recent
    // unbalanced force (dynamic-relaxation convergence indicator); <0 until first solve.
    double getUnbalanceNorm(void) const;

    // Methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    // Ladruno (W1-E2, was ExplicitBatheLNVD): assemble the residual and, when -lnvd is
    // active, add the FLAC local non-viscous damping. A pure pass-through to the base
    // IncrementalIntegrator::formUnbalance() when -lnvd is off (byte-identical).
    int formUnbalance(void);

    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Re-sync the integrator with a REVERTED Domain after a failed step. The
    // analysis calls Domain::revertToLastCommit() (node/element trial state,
    // time — including the mid-step p*dt / (1-p)*dt advances — and loads) and
    // then this. Noh-Bathe defers its cross-step (u,v,a)_t advance to commit(),
    // so newStep/update aborts leave the vectors clean; this override closes the
    // remaining windows (a commitDomain() failure after the advance, stale
    // prevKE / lastUnbalanceNorm from the aborted solve) by re-seeding
    // (u,v,a)_t from the COMMITTED node state and restoring the per-step
    // scalars. Time is deliberately NOT touched (revertToLastCommit covers it).
    int revertToLastStep(void);

    // Method to obtain current velocity
    const Vector &getVel(void);

    // Ladruno (ADR-30 P5): the LadrunoProjectionHandler pushes its (non-owning)
    // constraint projector through this seam at doneNumberingDOF().
    void setConstraintProjector(LadrunoConstraintProjector *theProjector) override;

    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:
    // Ladruno (W1-E2): master constructor (the real init). The integrator classTag is
    // passed in (derived from the flag combo by the public ctors / makeForBroker), so a
    // serialized object reports the matching legacy tag.
    ExplicitBathe(int classTag, double p, int compute_critical_timestep_,
                  bool verbose, bool cflAbort, double divergenceFactor,
                  bool cflUseTangent, int cflRecomputeEvery, CTSLumping lumping,
                  bool useLNVD, double alpha_flac,
                  bool useSMS, bool useConsistent,
                  double dtTarget, double maxAddedMassFrac,
                  double pcgTol, int pcgMaxIt);

    // Ladruno (W1-E2): map the orthogonal flag combo <-> the six legacy integrator tags.
    static int tagForFlags(bool lnvd, bool sms, bool consistent);
    static void flagsForTag(int classTag, bool &lnvd, bool &sms, bool &consistent);

    // Ladruno (W1-E2): consistent-scaling acceleration refinement, invoked after BOTH
    // Noh-Bathe sub-step diagonal solves. refinesAccel() is true only when
    // useSMS && useConsistent; otherwise the leap-frog stays byte-identical.
    virtual bool refinesAccel(void) const;
    virtual int  refineAccel(Vector &a);

private:
    // Ladruno (W1-E2, was ExplicitBatheLNVD): add the FLAC local non-viscous damping
    // to the assembled SOE residual (only reached when useLNVD && alpha_flac>0).
    void addLocalDamping(void);

    // Ladruno (W1-E2, was ExplicitBatheSMS): apply / undo the additive nodal mass and
    // (re-)build the selective mass scaling for the current mesh. No-ops when !useSMS.
    void removeScaling(void);     // undo lumped nodal injection (re-baseline)
    int  applyMassScalingSMS(void);   // build + apply for the current domainChanged

    // Time step
    double deltaT;

    // Explicit Bathe parameters
    double p;           // Sub-step parameter (0 < p < 1), controls damping
    double q0;          // Integration coefficient q0
    double q1;          // Integration coefficient q1
    double q2;          // Integration coefficient q2

    // State vectors at different time levels
    // The scheme has two sub-steps per time step:
    // Step 1: t -> t + p*dt (using acceleration A_tpdt)
    // Step 2: t + p*dt -> t + dt (using accelerations A_tpdt and A_tdt)

    Vector *U_t;        // Displacement at time t
    Vector *V_t;        // Velocity at time t
    Vector *A_t;        // Acceleration at time t

    Vector *U_tpdt;     // Displacement at time t + p*dt
    Vector *V_tpdt;     // Velocity at time t + p*dt
    Vector *V_fake;     // Temporary velocity for setting response
    Vector *A_tpdt;     // Acceleration at time t + p*dt

    Vector *U_tdt;      // Displacement at time t + dt
    Vector *V_tdt;      // Velocity at time t + dt
    Vector *A_tdt;      // Acceleration at time t + dt
    Vector *R_tdt;      // Forces at time t + dt (unused, kept for compatibility)

    // Update counter
    int updateCount;    // Counts external update() calls per step (exactly ONE —
                        //   the second Noh-Bathe solve is internal to update())

    // Integration coefficients (computed from p)
    double a0;          // = p * dt
    double a1;          // = (p * dt)^2 / 2
    double a2;          // = a0 / 2
    double a3;          // = (1 - p) * dt
    double a4;          // = ((1 - p) * dt)^2 / 2
    double a5;          // = q0 * a3
    double a6;          // = (0.5 + q1) * a3
    double a7;          // = q2 * a3

    // Critical time step computation
    int compute_critical_timestep;          // Control flag for critical dt computation
    double damped_minimum_critical_timestep;    // Minimum critical dt (damped)
    double undamped_minimum_critical_timestep;  // Minimum critical dt (undamped)
    int damped_critical_element_tag;            // Element with minimum damped dt
    int undamped_critical_element_tag;          // Element with minimum undamped dt

    // Diagnostics / run control
    bool verbose;             // gate the per-step opserr reporting (default off)
    bool cflAbort;            // abort the step if dt exceeds the Noh-Bathe limit
    double divergenceFactor;  // >0: abort if kinetic energy grows by this factor
                              //     in one step (spurious-energy circuit breaker)
    double prevKE;            // previous-step kinetic-energy proxy (0.5*v.v)
    double committedPrevKE;   // prevKE at newStep() entry (revertToLastStep; transient)
    bool firstStep;           // first newStep() seeds prevKE / checks cold start

    // critical-time-step refresh policy
    bool cflUseTangent;       // use getTangentStiff() instead of getInitialStiff()
    int cflRecomputeEvery;    // recompute dt_cr every N steps (0 = once); tracks
                              // stiffness changes in nonlinear runs
    int cflStepCount;         // step counter for periodic recompute
    int committedCflStepCount;// cflStepCount at newStep() entry (revert cadence)
    bool cflFirstComputation; // gate the detailed report to the first computation
    CTSLumping lumping;       // element-mass lumping for the dt_cr pencil (-lump)

    // Ladruno (W1-E2, -lnvd): FLAC local non-viscous damping (was ExplicitBatheLNVD)
    bool   useLNVD;           // route the residual through addLocalDamping()
    double alpha_flac;        // FLAC local damping coefficient (0 <= alpha < 1)
    double lastUnbalanceNorm; // |r|_inf at the most recent solve (-1 until first)
    double committedUnbalanceNorm; // lastUnbalanceNorm at newStep() entry (revert)

    // Ladruno (W1-E2, -sms / -consistent): selective mass scaling (was ExplicitBatheSMS /
    // ExplicitBatheSMSConsistent). Inert unless useSMS.
    bool   useSMS;            // selective mass scaling active
    bool   useConsistent;     // consistent (Olovsson) variant (requires useSMS)
    double dtTarget;          // target stable step the scaling sizes to
    double maxAddedMassFrac;  // soft cap on added-mass fraction (warn if exceeded)
    double smsEffectiveLimit; // POST-SCALING effective stable step (dtTarget capped by
                              //   still-governing excluded/self-reported elements); set
                              //   by applyMassScalingSMS, consumed by the newStep()
                              //   dt_cr report so it stops warning against the
                              //   PRE-scaling pencil (transient; NOT serialized)
    CTSLumping lumpingSMS;    // lumping used by the scaling sizing
    bool   useTangentSMS;     // size scaling from the tangent stiffness
    double pcgTol;            // consistent PCG tolerance
    int    pcgMaxIt;          // consistent PCG max iterations
    // lumped path state (additive nodal injection)
    std::map<int, Vector> injected;  // per-node-tag injected diagonal ΔM (re-baseline)
    bool   scaled;                   // injection currently applied to appliedDomain
    Domain *appliedDomain;           // the Domain the injection sits on
    // consistent path state (matrix-free PCG)
    std::vector<Ladruno::ConsistentBlock> *blocks;  // per-element M_bar blocks
    bool   warnedLimitations;        // one-time SMS limitations warning
    bool   warnedSolver;             // one-time consistent SOE-type warning
    bool   warnedPCG;                // one-time consistent PCG-failure warning

    // Ladruno (ADR-30 P5): constraint projection (non-owning; owned by the handler).
    LadrunoConstraintProjector *theProjector;  // 0 unless LadrunoProjection is active
    bool massBuilt;                            // has the projector read diag(M) yet?
};

#endif
