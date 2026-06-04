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

class DOF_Group;
class FE_Element;
class Vector;
class LinearSOE;

class CentralDifferenceLadruno : public TransientIntegrator
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

    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Modal-damping hook: returns the LAGGED half-step velocity (ADR C5).
    const Vector &getVel(void);

    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:

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
    double prevKE;             // previous-step kinetic-energy proxy (0.5 v.v)

    // critical-time-step refresh policy
    bool   cflUseTangent;      // use getTangentStiff() instead of getInitialStiff()
    int    cflRecomputeEvery;  // refresh dt_cr every N steps (0 = once); tracks
                               //   stiffening/softening in nonlinear runs
    int    cflStepCount;       // step counter for periodic recompute
    bool   cflFirstComputation;// gate the detailed report to the first computation
    CTSLumping lumping;        // element-mass lumping for the dt_cr pencil (-lump)
    bool   betaKWarned;        // emit the beta-K-trap warning at most once
};

#endif
