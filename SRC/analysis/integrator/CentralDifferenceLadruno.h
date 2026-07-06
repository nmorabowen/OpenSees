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

#ifndef CentralDifferenceLadruno_h
#define CentralDifferenceLadruno_h

// Written: nmb (UANDES), 05/2026  (Ladruno sibling-fork explicit integrator)
// Created: May 2026
//
// Description: A single, clean EXPLICIT LEAP-FROG central-difference integrator.
// It is the genuinely-missing combination in OpenSees: the explicit leap-frog
// scheme (as in ExplicitDifference) done right -- with a CORRECT first step, a
// built-in critical-time-step (dt_cr) guard, CLEAN full-step velocity output, and
// energy-balance discipline -- WITHOUT touching any upstream class (sibling-fork
// policy). See Ladruno_implementation/05_robust_central_difference.md (the ADR).
//
// Scheme (mass-only LHS, viscous damping lagged at the known half-step velocity):
//
//     a_n      = M^{-1} ( P_n - C * v_{n-1/2} - F_int(u_n) )    // diagonal M^{-1}
//     v_{n+1/2} = v_{n-1/2} + dt * a_n
//     u_{n+1}   = u_n + dt * v_{n+1/2}
//
// This is Newmark(beta=0, gamma=1/2): 2nd-order accurate, ZERO algorithmic
// dissipation (slight period SHORTENING ~ -Omega^2/24, not elongation). Good for
// wave propagation / energy conservation; for controllable high-frequency
// dissipation use ExplicitBathe instead.
//
// Two velocities, kept separate (ADR decision C5):
//   - getVel()                : the HALF-step velocity v_{n-1/2} known at solve
//                               time. This is the modal-damping hook
//                               (IncrementalIntegrator::addModalDampingForce ->
//                               residual via setB), so it MUST be the lagged
//                               velocity. EnergyBalanceRecorder's damping-work
//                               term integrates against this same velocity.
//   - node/recorder velocity  : the FULL-step v_n = 1/2 (v_{n-1/2} + v_{n+1/2})
//                               = v_{n-1/2} + 1/2 dt a_n, pushed via setResponse().
//                               Fixes the legacy half-step / crude-reconstruct
//                               output of ExplicitDifference.
//
// Correct first step (ADR decisions C3 / B1). The starter
//     a_0   = M^{-1} ( P_0 - C v_0 - F_int(u_0) )
//     v_{-1/2} = v_0 - 1/2 dt a_0       // seed the half-step velocity with THIS
// is computed on the FIRST newStep() (where dt is set and the SOE can be formed
// and factored) -- NOT in domainChanged(), which has neither a dt nor a factored
// SOE (exactly why the legacy CD classes punt with "assuming Ut-1 = Ut"). This is
// the correctness fix the legacy classes never made (acceptance test #8).
//
// Critical time step (ADR decision C4). Reuses CriticalTimeStep::
// computeCriticalTimeStep() (its own LAPACK eigensolve of K v = lambda M v, so it
// does NOT need the global SOE) and reports the central-difference limit
// dt_cr = 2/omega_max with STABILITY FACTOR 1.0 (NO Noh-Bathe 2x bonus -- the CD
// stability limit is exactly 2/omega_max). When damping is active, damping REDUCES
// the step to dt_cr = (2/omega_max)(sqrt(1+xi^2) - xi), and the DAMPED value is
// reported. getCriticalTimeStep() is overridden (base returns -1.0); the existing
// criticalTimeStep() Py/Tcl command dispatches to it.
//
// beta-K trap (ADR decision C6). Stiffness-proportional Rayleigh damping
// (xi = beta*omega/2 grows with frequency) COLLAPSES the explicit step
// ~ 2/(beta*omega_max^2). We detect this data-driven -- when the damped dt_cr is
// markedly below the undamped one -- and WARN (prefer mass-proportional alphaM,
// for which xi = alpha/2omega is benign at high frequency). We do not refuse;
// -cflAbort is the hard guard.
//
// Required recipe (explicit): lumped/diagonal element mass, `system Diagonal`
// (trivial M^{-1}), `algorithm Linear` (exactly ONE solve/step -- guarded at
// updateCount > 1), and dt < dt_cr. For the coupled / implicit-damped /
// consistent-mass central-difference case, use the existing
// `integrator NewmarkExplicit 0.5` instead -- this class is the explicit
// leap-frog scheme only.

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>   // CTSLumping, CTSResult, computeCriticalTimeStep()
#include <LadrunoProjectionConsumer.h>   // Ladruno: ADR-30 explicit constraint projection

class DOF_Group;
class FE_Element;
class Vector;
class LinearSOE;
class LadrunoConstraintProjector;

class CentralDifferenceLadruno : public TransientIntegrator,
                                 public LadrunoProjectionConsumer
{
public:
    // Constructors
    CentralDifferenceLadruno();
    CentralDifferenceLadruno(int compute_critical_timestep,
                             bool verbose, bool cflAbort,
                             double divergenceFactor,
                             bool cflUseTangent, int cflRecomputeEvery,
                             CTSLumping lumping);

    // Destructor
    ~CentralDifferenceLadruno();

    // Central-difference critical time step (2/omega_max, or the damped limit when
    // damping is active); <=0 if not computed / not applicable. ADR C4.
    double getCriticalTimeStep(void) const override;

    // Methods which define what the FE_Element and DOF_Groups add to the SOE.
    // Explicit M-only: zeroTangent(); addMtoTang();  (NO damping on the LHS).
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    // Ladruno (ADR-67 P-NEW-1): constant-mass tangent cache. The assembled
    // tangent is PURE MASS (see formEleTangent/formNodTangent above — no K, no C,
    // no dt; the modal-damping branch in the base formTangent is scaled by the
    // inherited getCFactor()==0), so it is constant between M-changing /
    // SOE-resizing events. Once formed, skip the whole assembly loop; the
    // Diagonal solver's factored-in-place state then persists exactly as under
    // `algorithm Linear -factorOnce` (measured −17.9% explicit wall,
    // md5-bit-identical). Invalidated in domainChanged() (contact re-emission,
    // element removal/birth, staged activation, SMS injection — every SOE resize),
    // setConstraintProjector() and revertToLastStep() (coupled with massBuilt so
    // the two mass-re-read latches never disagree). KNOWN EXCLUSION (documented,
    // adversarial gate item 8): updateParameter on a density changes M WITHOUT a
    // domainChanged — use -noMassCache (or fire a manual domainChange) for
    // mid-run density updates.
    int formTangent(int statFlag);

    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Re-sync the integrator with a REVERTED Domain after a failed step. The
    // analysis calls Domain::revertToLastCommit() (node/element trial state,
    // time, loads) and then this; the leap-frog advance in newStep() has already
    // mutated Ut/Vhalf/Aprev by the time any failure exit (incl. the -divergence
    // and NaN circuit breakers) fires, so without this override an adaptive
    // driver that halves dt and retries continues from poisoned state. Strategy:
    // re-seed from the COMMITTED DOF state and re-arm the first-step starter, so
    // the retry rebuilds v_{-1/2} = v_n - dt_retry/2 a_n for WHATEVER dt it uses
    // (a literal snapshot restore would hand it a half-step centered for the OLD
    // dt) — retry-after-abort == fresh restart from the committed state. NOTE
    // the restart is machine-exact for undamped models; with damping it is a
    // well-defined second-order restart (the starter evaluates C·v at the
    // full-step v_n instead of the lagged v_{n+1/2}).
    int revertToLastStep(void);

    // Modal-damping hook: returns the LAGGED half-step velocity (ADR C5).
    const Vector &getVel(void);

    // Ladruno (ADR-39 B1): the current step's Δt, for the SOFT=1 Courant-stable contact
    // penalty (k_soft = SOFSCL·4·m_eff/Δt²). The contact adapter dynamic_casts the
    // integrator it is handed in getResidual() to this (explicit-only) class and reads Δt
    // so the contact's own ω·Δt ≤ 2 ⇒ the contact never throttles dt_cr below the global
    // step. > 0 only once newStep() has set it; the cast fails under any implicit
    // integrator ⇒ SOFT is inert there (the configured kn is used ⇒ byte-identical).
    double getCurrentDeltaT(void) const { return deltaT; }

    // Ladruno (ADR-67 P-NEW-1): parser hook for -noMassCache (escape hatch —
    // restores today's every-step mass reassembly; needed for mid-run
    // updateParameter-on-density decks).
    void setMassCacheEnabled(bool on) { useMassCache = on; if (!on) massTangentValid = false; }

    // Ladruno (ADR-30): LadrunoProjectionConsumer — the handler pushes its
    // (non-owning) acceleration projector here; we project a0 in the starter and
    // a_{n+1} in update() so MP constraints are enforced without penalty/elimination.
    void setConstraintProjector(LadrunoConstraintProjector *theProjector_) override;

    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:
    // Delegating constructor for subclasses (e.g. CentralDifferenceSMS) that must
    // carry their OWN integrator classTag through the FEM_ObjectBroker round-trip.
    CentralDifferenceLadruno(int classTag,
                             int compute_critical_timestep,
                             bool verbose, bool cflAbort,
                             double divergenceFactor,
                             bool cflUseTangent, int cflRecomputeEvery,
                             CTSLumping lumping);

    // Ladruno (ADR-38): acceleration-refinement hook for CONSISTENT (Olovsson) mass
    // scaling. The diagonal `system Diagonal` solve gives a = M_lump^-1 r; the
    // consistent path needs a = M_tilde^-1 r with M_tilde = M_lump + sum_e M_bar_e
    // (non-diagonal). A subclass (CentralDifferenceSMSConsistent) overrides
    // refinesAccel()->true and refineAccel() to recover r = M_lump .* a and replace
    // a by the matrix-free PCG solve. Called at the TWO sites that consume the diagonal
    // solve: the newStep() starter and update(). The default is a no-op, so the lumped
    // path (this base + CentralDifferenceSMS) stays BYTE-IDENTICAL.
    virtual bool refinesAccel(void) const { return false; }
    virtual int  refineAccel(Vector &a)   { return 0; }

    // Ladruno (ADR-36/52 review 2026-07-01): the SMS subclasses report their
    // POST-SCALING effective stable step here (dtTarget capped by any excluded /
    // self-reported element that still governs). > 0 flags "mass scaling active"
    // to newStep(), which then (a) labels the step-1 dt_cr report a PRE-SCALING
    // estimate and (b) warns "expect INSTABILITY" only when dt exceeds THIS limit
    // — the un-augmented element pencil cannot see the injected mass, so warning
    // against it cried wolf on every correctly-scaled run.
    void setSMSEffectiveLimit(double lim) { smsEffectiveLimit = lim; }

private:
    // Time step
    double deltaT;

    // Primary leap-frog state
    Vector *Ut;        // displacement u_n (advanced to u_{n+1} during newStep)
    Vector *Vhalf;     // HALF-step velocity v_{n-1/2} (primary state; what getVel
                       //   returns at solve time)
    Vector *Aprev;     // acceleration a_n carried between steps (drives the newStep
                       //   velocity advance); reset each solve
    Vector *Vfull;     // FULL-step velocity output buffer v_n = 1/2(v_-+v_+)
    Vector *Azero;     // a persistent zero vector -> setAccel(0) before each solve
                       //   so the residual carries NO inertia term (M a on RHS)

    // Update counter (CD allows EXACTLY one solve/step; guard at > 1)
    int updateCount;

    // First-step starter gate (ADR C3 / B1)
    bool firstStep;

    // Critical-time-step computation (mirrors ExplicitBathe surface)
    int    compute_critical_timestep;     // 0 = off, 1 = (re)compute, 2 = computed
    double damped_minimum_critical_timestep;    // 2/omega_max incl. element damping
    double undamped_minimum_critical_timestep;  // 2/omega_max, undamped
    int    damped_critical_element_tag;
    int    undamped_critical_element_tag;

    // Diagnostics / run control
    bool   verbose;            // per-step dt/energy reporting (default off)
    bool   cflAbort;           // abort the step if dt exceeds the CD limit
    double divergenceFactor;   // >0: abort on runaway kinetic-energy growth
    double prevKE;             // RUNNING-MAX kinetic-energy proxy (0.5 v.v) --
                               //   the breaker baseline (a previous-step baseline
                               //   false-tripped at free-vibration KE troughs)
    double committedPrevKE;    // prevKE at newStep() entry — restored by
                               //   revertToLastStep so the breaker stays armed
                               //   with the pre-fault baseline (transient; NOT
                               //   serialized)

    // critical-time-step refresh policy
    bool   cflUseTangent;      // use getTangentStiff() instead of getInitialStiff()
    int    cflRecomputeEvery;  // refresh dt_cr every N steps (0 = once); tracks
                               //   stiffening/softening in nonlinear runs
    int    cflStepCount;       // step counter for periodic recompute
    int    committedCflStepCount; // cflStepCount at newStep() entry (revert cadence)
    bool   cflFirstComputation;// gate the detailed report to the first computation
    CTSLumping lumping;        // element-mass lumping for the dt_cr pencil (-lump)
    bool   betaKWarned;        // emit the beta-K-trap warning at most once
    double smsEffectiveLimit;  // post-scaling effective stable step pushed by an SMS
                               //   subclass (0 = no mass scaling); see the protected
                               //   setter (transient; NOT serialized)

    // Shared DOF-loop: seed Ut/Vhalf/Aprev from the committed node state (used by
    // domainChanged and revertToLastStep; Vhalf gets the full-step v_n — the
    // first-step starter re-centers it to v_{-1/2}).
    void seedFromCommitted(AnalysisModel *theModel);

    // Ladruno (ADR-30): explicit constraint projection (null unless the handler is
    // LadrunoProjectionHandler). theProjector is NON-owning (handler owns it).
    LadrunoConstraintProjector *theProjector;
    bool   massBuilt;          // projector mass read this domainChanged() yet?
    Vector *Aproj;             // projected-acceleration scratch for update()

    // Ladruno (ADR-67 P-NEW-1): constant-mass tangent cache state. ALL transient —
    // NOT serialized (a recvSelf-built integrator starts invalid and re-forms).
    bool   useMassCache;       // default true; -noMassCache restores every-step assembly
    bool   massTangentValid;   // tangent assembled since the last invalidation?
    bool   massCacheNoted;     // one-time density-parameter contract note emitted?
};

#endif
