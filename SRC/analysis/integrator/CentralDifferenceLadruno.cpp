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

// Written: nmb (UANDES), 05/2026
//
// Implementation of CentralDifferenceLadruno -- a single, clean explicit
// leap-frog central-difference integrator (classTag 33003, Ladruno private band).
// See the class header
// and Ladruno_implementation/05_robust_central_difference.md for the full ADR.

#include <CentralDifferenceLadruno.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Domain.h>        // Ladruno (ADR-30 P4): tie-force scatter to nodes
#include <Node.h>
#include <NodeIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <LadrunoConstraintProjector.h>   // Ladruno: ADR-30
#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

#define OPS_Export

// The central-difference stability limit is EXACTLY 2/omega_max -- no Noh-Bathe 2x
// bonus (ADR C4). Stability factor is therefore 1.0 and is left implicit below.

// When the damping-reduced critical step drops below this fraction of the
// undamped one, we flag the beta-K trap. Mass-proportional alphaM (xi = alpha/2w)
// barely moves the ratio at omega_max, so this fires for stiffness-proportional
// betaK (xi = beta*w/2, which collapses dt_cr ~ 2/(beta*w_max^2)) but not for a
// benign alphaM (ADR C6 / acceptance test #9).
#define CDL_BETAK_WARN_RATIO 0.90

void *OPS_CentralDifferenceLadruno(void)
{
    // Usage: integrator CentralDifferenceLadruno <-cfl> <-cflAbort> <-tangent>
    //                                            <-recompute N> <-lump rowsum|diagonal|hrz>
    //                                            <-verbose> <-divergence f>
    // No positional argument and NO -damping flag: this is ONE explicit scheme.
    // (For coupled / implicit-damped CD use `integrator NewmarkExplicit 0.5`.)
    int compute_critical_timestep = 0;
    bool verbose = false;
    bool cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    CTSLumping lumping = CTSLumping::Diagonal;   // explicit default: diagonal-of-
                                                 // consistent is robust for the
                                                 // rotational DOFs of beams/shells

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verbose = true;
        } else if (strcmp(arg, "-cflAbort") == 0) {
            cflAbort = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-cfl") == 0 ||
                   strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-divergence") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetDoubleInput(&nd, &divergenceFactor);
            }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-recompute") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetIntInput(&nd, &cflRecomputeEvery);
            }
            if (cflRecomputeEvery <= 0)
                opserr << "WARNING CentralDifferenceLadruno -recompute expects a positive "
                          "integer N (steps between dt_cr refreshes); dt_cr computed once\n";
            cflUseTangent = true;   // recomputing the initial stiffness every N steps is pointless
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING CentralDifferenceLadruno - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else {
            opserr << "WARNING CentralDifferenceLadruno - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    TransientIntegrator *theIntegrator =
        new CentralDifferenceLadruno(compute_critical_timestep, verbose, cflAbort,
                                     divergenceFactor, cflUseTangent,
                                     cflRecomputeEvery, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CentralDifferenceLadruno integrator\n";
    return theIntegrator;
}

// Default constructor (broker)
CentralDifferenceLadruno::CentralDifferenceLadruno()
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceLadruno),
      deltaT(0.0),
      Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
      updateCount(0), firstStep(true),
      compute_critical_timestep(0),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(false), cflAbort(false), divergenceFactor(0.0), prevKE(0.0),
      committedPrevKE(0.0),
      cflUseTangent(false), cflRecomputeEvery(0), cflStepCount(0),
      committedCflStepCount(0),
      cflFirstComputation(true), lumping(CTSLumping::Diagonal), betaKWarned(false),
      smsEffectiveLimit(0.0),
      theProjector(0), massBuilt(false), Aproj(0)
{
}

// Main constructor (public) — delegates to the classTag-carrying protected ctor.
CentralDifferenceLadruno::CentralDifferenceLadruno(
    int compute_critical_timestep_, bool verbose_, bool cflAbort_,
    double divergenceFactor_, bool cflUseTangent_, int cflRecomputeEvery_,
    CTSLumping lumping_)
    : CentralDifferenceLadruno(INTEGRATOR_TAGS_CentralDifferenceLadruno,
                               compute_critical_timestep_, verbose_, cflAbort_,
                               divergenceFactor_, cflUseTangent_,
                               cflRecomputeEvery_, lumping_)
{
}

// Protected constructor — subclasses (CentralDifferenceSMS) pass their own classTag.
CentralDifferenceLadruno::CentralDifferenceLadruno(
    int classTag,
    int compute_critical_timestep_, bool verbose_, bool cflAbort_,
    double divergenceFactor_, bool cflUseTangent_, int cflRecomputeEvery_,
    CTSLumping lumping_)
    : TransientIntegrator(classTag),
      deltaT(0.0),
      Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
      updateCount(0), firstStep(true),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(verbose_), cflAbort(cflAbort_), divergenceFactor(divergenceFactor_),
      prevKE(0.0), committedPrevKE(0.0), cflUseTangent(cflUseTangent_),
      cflRecomputeEvery(cflRecomputeEvery_), cflStepCount(0),
      committedCflStepCount(0),
      cflFirstComputation(true), lumping(lumping_), betaKWarned(false),
      smsEffectiveLimit(0.0),
      theProjector(0), massBuilt(false), Aproj(0)
{
}

CentralDifferenceLadruno::~CentralDifferenceLadruno()
{
    if (Ut    != 0) delete Ut;
    if (Vhalf != 0) delete Vhalf;
    if (Aprev != 0) delete Aprev;
    if (Vfull != 0) delete Vfull;
    if (Azero != 0) delete Azero;
    if (Aproj != 0) delete Aproj;
    // theProjector is owned by LadrunoProjectionHandler — do NOT delete here.
}

// Ladruno (ADR-30): the projection handler pushes its (non-owning) projector. We
// re-read the mass next solve (massBuilt=false) since a new projector implies a new
// constraint set / numbering.
void
CentralDifferenceLadruno::setConstraintProjector(LadrunoConstraintProjector *p)
{
    theProjector = p;
    massBuilt = false;
}

// Explicit: only the mass matrix on the LHS -> a trivial diagonal M^{-1} solve.
// Damping is NOT on the LHS (that would be the coupled scheme = NewmarkExplicit);
// it enters the residual at the lagged half-step velocity set on the nodes.
int CentralDifferenceLadruno::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

int CentralDifferenceLadruno::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// The reported critical step: the damped limit when damping is active, otherwise
// the undamped 2/omega_max. Both come from the shared CriticalTimeStep eigensolve.
static inline double cdl_reported(double damped, double undamped)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (undamped == inf || undamped <= 0.0)
        return -1.0;
    if (damped != inf && damped > 0.0 && damped < undamped)
        return damped;   // damping reduces the step (ADR T3) -> report damped_dt
    return undamped;
}

int CentralDifferenceLadruno::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (theModel == 0 || theLinSOE == 0) {
        opserr << "CentralDifferenceLadruno::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    if (size == 0) {
        opserr << "CentralDifferenceLadruno::domainChanged - invalid size\n";
        return -1;
    }

    // (Re)allocate state if the size changed.
    if (Ut == 0 || Ut->Size() != size) {
        if (Ut    != 0) delete Ut;
        if (Vhalf != 0) delete Vhalf;
        if (Aprev != 0) delete Aprev;
        if (Vfull != 0) delete Vfull;
        if (Azero != 0) delete Azero;
        if (Aproj != 0) delete Aproj;

        Ut    = new Vector(size);
        Vhalf = new Vector(size);
        Aprev = new Vector(size);
        Vfull = new Vector(size);
        Azero = new Vector(size);   // stays zero for the lifetime of the object
        Aproj = new Vector(size);   // Ladruno (ADR-30): projected-accel scratch

        if (Ut == 0 || Ut->Size() != size ||
            Vhalf == 0 || Vhalf->Size() != size ||
            Aprev == 0 || Aprev->Size() != size ||
            Vfull == 0 || Vfull->Size() != size ||
            Azero == 0 || Azero->Size() != size ||
            Aproj == 0 || Aproj->Size() != size) {
            opserr << "CentralDifferenceLadruno::domainChanged - out of memory\n";
            return -1;
        }
    }
    Azero->Zero();

    // Seed u_n, v_{n-1/2} (provisional = v_0) and a_n (provisional = committed a_0)
    // from the committed DOF state. The TRUE starter a_0 and v_{-1/2} = v_0 - dt/2 a_0
    // are computed on the first newStep() (ADR C3 / B1) -- NOT here, where there is
    // neither a dt nor a factored SOE.
    this->seedFromCommitted(theModel);

    firstStep = true;
    cflStepCount = 0;

    // The KE-proxy divergence baseline is stale across a mesh change (a smaller /
    // re-numbered model would compare its fresh KE against the old model's) — the
    // breaker re-arms on the first nonzero KE of the new mesh.
    prevKE = 0.0;
    committedPrevKE = 0.0;

    // Ladruno (ADR-30): the projected-mass cache is invalid after a domain change;
    // re-read diag(M) at the next starter. Then check (or project) IC compliance of
    // the committed (u0, v0) against every constraint group — Ut/Vhalf currently hold
    // them. enforceIC() needs no mass (it uses L, delta only), so it is safe here.
    massBuilt = false;
    if (theProjector != 0) {
        int nViol = theProjector->enforceIC(*Ut, *Vhalf);
        if (nViol > 0) {
            opserr << "CentralDifferenceLadruno::domainChanged - " << nViol
                   << " initial-condition constraint violation(s); aborting. Fix the "
                      "ICs or add -projectICs to the constraints command.\n";
            return -1;
        }
    } else {
        // A LadrunoProjectionHandler would have pushed a projector at
        // doneNumberingDOF(); none here means either a different handler (fine) or a
        // mid-run integrator swap that bypassed doneNumberingDOF (then the IC check is
        // silently skipped — warn once is handled by the handler path).
    }

    // dt_cr is computed here unconditionally (ADR C4): CriticalTimeStep runs its
    // OWN LAPACK eigensolve of K v = lambda M v, so it needs neither a dt nor the
    // factored global SOE. This makes the criticalTimeStep() query valid BEFORE the
    // first analyze. The compute_critical_timestep flag controls only the per-step
    // reporting / -cflAbort / -recompute behavior, not whether dt_cr exists.
    CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
    damped_minimum_critical_timestep   = r.damped_dt;
    undamped_minimum_critical_timestep = r.undamped_dt;
    damped_critical_element_tag        = r.damped_tag;
    undamped_critical_element_tag      = r.undamped_tag;
    if (compute_critical_timestep == 2)   // recompute cadence: mark fresh-needed
        compute_critical_timestep = 1;

    const double inf = std::numeric_limits<double>::infinity();
    const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                      undamped_minimum_critical_timestep);

    // beta-K trap detection (ADR C6 / test #9): damping markedly below the undamped
    // limit is the stiffness-proportional-Rayleigh signature.
    if (!betaKWarned &&
        damped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep > 0.0 &&
        damped_minimum_critical_timestep <
            CDL_BETAK_WARN_RATIO * undamped_minimum_critical_timestep) {
        opserr << "CentralDifferenceLadruno - NOTE damping is significantly reducing the "
                  "critical time step (damped dt_cr = " << damped_minimum_critical_timestep
               << " vs undamped " << undamped_minimum_critical_timestep << ").\n"
                  "  In an explicit scheme stiffness-proportional (betaK) Rayleigh damping "
                  "collapses dt_cr ~ 2/(betaK*omega_max^2); prefer mass-proportional alphaM "
                  "(xi = alphaM/2omega is benign at high frequency). Reporting the damped dt_cr.\n";
        betaKWarned = true;
    }

    if (compute_critical_timestep || verbose) {
        opserr << "CentralDifferenceLadruno: critical time step estimate"
               << (cflUseTangent ? " (tangent)" : "") << "\n"
               << "  central-difference limit (2/omega_max): "
               << undamped_minimum_critical_timestep
               << " @ element #" << undamped_critical_element_tag << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            opserr << "  damped limit (reported): " << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n";
        if (dt_cr <= 0.0)
            opserr << "  note: dt_cr not applicable (no element produced a finite, "
                      "positive estimate -- e.g. a pure nodal-mass model).\n";
    }

    return 0;
}

// Shared DOF-loop (domainChanged + revertToLastStep): seed the leap-frog state
// from the COMMITTED node state. Vhalf receives the full-step v_n; the first-step
// starter re-centers it to the half-step for whatever dt the next step uses.
void CentralDifferenceLadruno::seedFromCommitted(AnalysisModel *theModel)
{
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Ut)(loc) = disp(i);
        }
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Vhalf)(loc) = vel(i);
        }
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Aprev)(loc) = accel(i);
        }
    }
}

// Re-sync with a reverted Domain after a failed step (see the header comment).
// Contract: the analysis calls Domain::revertToLastCommit() FIRST (node/element
// trial state, currentTime, dT, loads — including everything newStep/update
// scattered via setResponse/updateDomain); here we restore ONLY the private
// leap-frog state, and we must NOT touch time (that would double-revert). The
// re-seed is valid even standalone: getCommittedDisp/Vel/Accel read committed
// values, which a failed (uncommitted) step never touched.
int CentralDifferenceLadruno::revertToLastStep(void)
{
    if (Ut == 0)                                  // domainChanged never ran
        return 0;
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)
        return 0;

    this->seedFromCommitted(theModel);

    // Re-arm the starter: the next newStep(dt_retry) rebuilds the TRUE
    // v_{-1/2} = v_n - dt_retry/2 a_n at the retry dt (this also repairs a
    // starter that failed halfway through, after mutating Vhalf).
    firstStep = true;
    updateCount = 0;
    prevKE = committedPrevKE;                     // keep the -divergence breaker armed
    cflStepCount = committedCflStepCount;         // keep the -recompute cadence
    return 0;
}

int CentralDifferenceLadruno::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (deltaT <= 0.0) {
        opserr << "CentralDifferenceLadruno::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -1;
    }

    if (Ut == 0 || Vhalf == 0 || Aprev == 0) {
        opserr << "CentralDifferenceLadruno::newStep() - domainChanged() failed or not called\n";
        return -2;
    }

    // Snapshot the per-step scalar state for revertToLastStep (the vector state
    // needs no snapshot — the committed node state IS the restart point).
    committedPrevKE       = prevKE;
    committedCflStepCount = cflStepCount;

    // Each step performs exactly ONE solve; reset here so a failed/retried step
    // cannot poison the next step's update() counter (ADR: guard at > 1).
    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    // -recompute N: refresh dt_cr from the (tangent) stiffness so it tracks
    // softening/stiffening in nonlinear runs.
    if (cflRecomputeEvery > 0 && compute_critical_timestep == 2) {
        if (++cflStepCount >= cflRecomputeEvery) {
            compute_critical_timestep = 1;
            cflStepCount = 0;
        }
    }
    if (compute_critical_timestep == 1) {
        CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
        damped_minimum_critical_timestep   = r.damped_dt;
        undamped_minimum_critical_timestep = r.undamped_dt;
        damped_critical_element_tag        = r.damped_tag;
        undamped_critical_element_tag      = r.undamped_tag;
        compute_critical_timestep = 2;

        const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (cflFirstComputation || verbose) {
            opserr << "CentralDifferenceLadruno: critical time step (2/omega_max)"
                   << (cflUseTangent ? " (tangent)" : "")
                   << (smsEffectiveLimit > 0.0 ? " [PRE-SCALING estimate]" : "")
                   << ": " << dt_cr
                   << " @ element #"
                   << ((dt_cr == damped_minimum_critical_timestep)
                           ? damped_critical_element_tag
                           : undamped_critical_element_tag)
                   << "\n";
            if (smsEffectiveLimit > 0.0) {
                // Mass scaling raised every scaled element's stable step to dtTarget;
                // the pre-scaling pencil is NOT the run's limit (warning against it
                // cried "expect INSTABILITY" on every correctly-scaled run). Compare
                // dt against the POST-SCALING effective limit instead (dtTarget
                // capped by any excluded / self-reported element that still governs).
                if (deltaT > smsEffectiveLimit)
                    opserr << "  WARNING dt = " << deltaT
                           << " exceeds the POST-SCALING effective stable step "
                           << smsEffectiveLimit
                           << " (dtTarget capped by excluded/self-reported elements)"
                              " - expect INSTABILITY.\n";
            } else if (dt_cr > 0.0 && deltaT > dt_cr)
                opserr << "  WARNING dt = " << deltaT
                       << " exceeds the central-difference stability limit " << dt_cr
                       << " - expect INSTABILITY.\n";
        }
        cflFirstComputation = false;
    }

    // Optional hard CFL guard (stability factor 1.0 -- no Noh-Bathe 2x).
    if (cflAbort) {
        const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (dt_cr > 0.0 && deltaT > dt_cr) {
            opserr << "CentralDifferenceLadruno::newStep() - ABORT: dt = " << deltaT
                   << " > central-difference stability limit " << dt_cr << "\n";
            return -2;
        }
    }

    // ---- First-step starter (ADR C3 / B1) -------------------------------------
    // Compute the TRUE initial acceleration at the committed configuration and seed
    // the back-half-step velocity v_{-1/2} = v_0 - dt/2 a_0. domainChanged() left
    // Vhalf = v_0 and Ut = u_0, so a load already applied at rest gets a consistent
    // a_0 (this is the fix the legacy CD classes never made). The SOE is formed and
    // factored HERE (it is not in domainChanged()).
    if (firstStep) {
        firstStep = false;

        // residual at (u_0, v_0, accel = 0): r = P_0 - C v_0 - F_int(u_0)
        theModel->setResponse(*Ut, *Vhalf, *Azero);
        double t0 = theModel->getCurrentDomainTime();
        if (theModel->updateDomain(t0, deltaT) < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: failed to update domain\n";
            return -3;
        }
        if (this->formTangent(CURRENT_TANGENT) < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: formTangent failed\n";
            return -3;
        }
        // Ladruno (ADR-30): read diag(M) for the projector from the assembled SOE NOW,
        // before solve() factors it in place (DiagonalDirectSolver overwrites A[i] with
        // 1/A[i] on the factor pass). Mass is constant between domain changes.
        if (theProjector != 0 && !massBuilt) {
            if (theProjector->buildMass(theLinSOE) < 0) {
                opserr << "CentralDifferenceLadruno::newStep() - starter: projector "
                          "buildMass failed (massless DOF / non-Diagonal SOE)\n";
                return -3;
            }
            massBuilt = true;
        }
        if (this->formUnbalance() < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: formUnbalance failed\n";
            return -3;
        }
        if (theLinSOE->solve() < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: SOE solve failed\n";
            return -3;
        }
        *Aprev = theLinSOE->getX();              // a_0 = M^{-1}(P_0 - C v_0 - F_int(u_0))
        if (this->refinesAccel())                // Ladruno (ADR-38): consistent solve
            this->refineAccel(*Aprev);           //   a_0 = M_tilde^-1 r (matrix-free PCG)
        if (theProjector != 0)
            theProjector->project(*Aprev);       // Ladruno: project a_0 onto the manifold
        Vhalf->addVector(1.0, *Aprev, -0.5 * deltaT);   // v_{-1/2} = v_0 - dt/2 a_0 (projected)
    }

    // ---- Leap-frog advance (ADR How) ------------------------------------------
    // v_{n+1/2} = v_{n-1/2} + dt a_n ;   u_{n+1} = u_n + dt v_{n+1/2}
    Vhalf->addVector(1.0, *Aprev, deltaT);   // Vhalf is now v_{n+1/2}
    Ut->addVector(1.0, *Vhalf, deltaT);      // Ut    is now u_{n+1}

    // Set the trial state for the (M-only) solve of a_{n+1}. Accel = 0 so the
    // residual carries NO inertia term; damping uses the lagged v_{n+1/2}.
    theModel->setResponse(*Ut, *Vhalf, *Azero);

    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0) {
        opserr << "CentralDifferenceLadruno::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

int CentralDifferenceLadruno::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1) {
        opserr << "WARNING CentralDifferenceLadruno::update() - called more than once -";
        opserr << " the central-difference scheme requires a LINEAR solution algorithm "
                  "(exactly one solve/step)\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceLadruno::update() - no AnalysisModel set\n";
        return -2;
    }
    if (Ut == 0) {
        opserr << "WARNING CentralDifferenceLadruno::update() - domainChanged() failed or not called\n";
        return -3;
    }
    if (U.Size() != Ut->Size()) {
        opserr << "WARNING CentralDifferenceLadruno::update() - Vectors of incompatible size ";
        opserr << " expecting " << Ut->Size() << " obtained " << U.Size() << endln;
        return -4;
    }

    // U is the solved acceleration a_{n+1}. Circuit breaker: catch a blown-up
    // (NaN/Inf) acceleration before it silently propagates. NOT via pNorm(0):
    // its max-compare skips NaN entries, so the old A_max != A_max test only
    // ever fired on +/-Inf (vectorIsFinite scans with std::isfinite).
    const double A_max = U.pNorm(0);   // diagnostics only (verbose print below)
    if (!vectorIsFinite(U)) {
        opserr << "CentralDifferenceLadruno::update() - ABORT: non-finite acceleration "
                  "(NaN/Inf) - the integration has diverged (dt likely too large).\n";
        return -5;
    }

    // Ladruno (ADR-30): project the solved acceleration onto the constraint manifold.
    // The projected a_{n+1} (Aproj) drives the full-step velocity, the node state, AND
    // — load-bearing for manifold invariance — next step's leap-frog via *Aprev below.
    // The NaN/Inf circuit breaker above deliberately ran on the RAW U.
    const Vector *Aused = &U;
    // Ladruno (ADR-38): consistent (Olovsson) mass scaling replaces the diagonal solve
    // a = M_lump^-1 r by a = M_tilde^-1 r (matrix-free PCG) BEFORE projection. Default
    // no-op (refinesAccel()==false) -> the lumped path is byte-identical.
    if (this->refinesAccel()) {
        *Aproj = U;
        this->refineAccel(*Aproj);
        Aused = Aproj;
    }
    if (theProjector != 0) {
        if (Aused != Aproj) *Aproj = *Aused;   // (no-refine case) stage raw U into Aproj
        theProjector->project(*Aproj);
        Aused = Aproj;
    }

    // Clean FULL-step velocity output at t_{n+1}:
    //   v_{n+1} = 1/2 (v_{n+1/2} + v_{n+3/2}) = v_{n+1/2} + dt/2 a_{n+1}
    // (the extra half-step is exactly the centered (u_{n+2}-u_n)/2dt identity).
    // Vhalf currently holds v_{n+1/2}.
    *Vfull = *Vhalf;
    Vfull->addVector(1.0, *Aused, 0.5 * deltaT);

    // Push the consistent full-step snapshot (u_{n+1}, v_{n+1}, a_{n+1}) to the nodes.
    theModel->setResponse(*Ut, *Vfull, *Aused);
    if (theModel->updateDomain() < 0) {
        opserr << "CentralDifferenceLadruno::update() - failed to update the domain\n";
        return -6;
    }

    // Optional circuit breaker: abort on runaway kinetic-energy growth. The
    // baseline is the RUNNING MAX of the KE proxy, NOT the previous step: in
    // free vibration KE passes through ~0 at every displacement extreme, and a
    // previous-step baseline left prevKE ~ eps there -- the quadratic regrowth
    // off that floor is an unbounded ratio (phase-luck false trips, review
    // 2026-07-02). For monotonic divergence the previous step IS the max, so
    // the trip behavior there is unchanged.
    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*Vfull) ^ (*Vfull));   // velocity-based KE proxy
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "CentralDifferenceLadruno::update() - ABORT: kinetic-energy proxy grew to "
                   << (ke / prevKE) << "x its running maximum (> " << divergenceFactor
                   << ") - spurious energy growth / instability.\n";
            return -7;
        }
        if (ke > prevKE) prevKE = ke;   // running max
    }

    if (verbose)
        opserr << "CentralDifferenceLadruno::update() max|a| = " << A_max << endln;

    // Carry a_{n+1} as a_n for the next step's leap-frog velocity advance. This is the
    // load-bearing manifold write (Gate-A R3): newStep() advances Vhalf/Ut with *Aprev,
    // so it MUST be the PROJECTED acceleration.
    *Aprev = *Aused;

    return 0;
}

int CentralDifferenceLadruno::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceLadruno::commit() - no AnalysisModel set\n";
        return -1;
    }

    // Ladruno (ADR-30 P4): scatter the projector's per-equation tie force
    // f = M(a_raw - a_proj) (cached during the step's update()) onto the nodes, so the
    // node-based LadrunoRecorder can record a `constraintTieForce` field. Done before
    // commitDomain() (recorders run after commit). Untied/SP-fixed DOFs -> 0.
    if (theProjector != 0 && theProjector->isMassReady()) {
        Domain *theDomain = theModel->getDomainPtr();
        if (theDomain != 0) {
            NodeIter &theNodes = theDomain->getNodes();
            Node *nodePtr;
            while ((nodePtr = theNodes()) != 0) {
                DOF_Group *dg = nodePtr->getDOF_GroupPtr();
                if (dg == 0) continue;
                const ID &id = dg->getID();
                int ndf = nodePtr->getNumberDOF();
                Vector tf(ndf);
                for (int d = 0; d < ndf && d < id.Size(); d++) {
                    int eqn = id(d);
                    tf(d) = (eqn >= 0) ? theProjector->tieForceAtEqn(eqn) : 0.0;
                }
                nodePtr->setProjectionTieForce(tf);
            }
        }
    }

    // Time was already advanced inside newStep()/update() via updateDomain(newtime, dt).
    return theModel->commitDomain();
}

// Modal-damping hook (ADR C5): returns the LAGGED half-step velocity. At every
// solve the node velocity was set to this same vector, so the residual's modal-
// damping force (addModalDampingForce -> setB) and the element Rayleigh force are
// both evaluated at v_{n-1/2}/v_{n+1/2}, consistent with the explicit scheme.
const Vector &CentralDifferenceLadruno::getVel(void)
{
    if (Vhalf == 0) {
        // called before domainChanged() allocated the state (integrator swapped in
        // mid-session and queried immediately) — dereferencing would crash.
        opserr << "WARNING CentralDifferenceLadruno::getVel() - called before "
                  "domainChanged(); returning an empty vector\n";
        static Vector emptyVel(0);
        return emptyVel;
    }
    return *Vhalf;
}

double CentralDifferenceLadruno::getCriticalTimeStep(void) const
{
    return cdl_reported(damped_minimum_critical_timestep,
                        undamped_minimum_critical_timestep);
}

int CentralDifferenceLadruno::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(7);
    data(0) = (double)compute_critical_timestep;
    data(1) = verbose ? 1.0 : 0.0;
    data(2) = cflAbort ? 1.0 : 0.0;
    data(3) = divergenceFactor;
    data(4) = cflUseTangent ? 1.0 : 0.0;
    data(5) = (double)cflRecomputeEvery;
    data(6) = (double)(int)lumping;   // 0=RowSum, 1=Diagonal, 2=HRZ (legacy: 1=Diag, 0=RowSum)

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceLadruno::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}

int CentralDifferenceLadruno::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceLadruno::recvSelf() - could not receive data\n";
        return -1;
    }
    compute_critical_timestep = (int)data(0);
    verbose          = (data(1) != 0.0);
    cflAbort         = (data(2) != 0.0);
    divergenceFactor = data(3);
    cflUseTangent    = (data(4) != 0.0);
    cflRecomputeEvery = (int)data(5);
    {   // 3-way decode (legacy binaries sent 1.0=Diagonal, 0.0=RowSum; both map back
        // correctly here). HRZ=2 from a newer peer; an OLD binary reading our HRZ would
        // see 2.0!=0 -> Diagonal, the only (documented, safe) lossy pairing.
        int lc = (int)data(6);
        lumping = (lc == 2) ? CTSLumping::HRZ
                : (lc == 0) ? CTSLumping::RowSum
                            : CTSLumping::Diagonal;
        if (lc < 0 || lc > 2)
            opserr << "WARNING CentralDifferenceLadruno::recvSelf - unknown lumping "
                   << "code " << lc << " (using Diagonal)\n";
    }

    return 0;
}

void CentralDifferenceLadruno::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    s << "CentralDifferenceLadruno - explicit leap-frog central difference\n";
    if (theModel != 0)
        s << "  currentTime: " << theModel->getCurrentDomainTime() << endln;
    s << "  recipe: lumped/diagonal mass, system Diagonal, algorithm Linear, dt < dt_cr\n";
    if (compute_critical_timestep > 0) {
        s << "  critical dt (2/omega_max): " << undamped_minimum_critical_timestep << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            s << "  critical dt (damped, reported): " << damped_minimum_critical_timestep << "\n";
    }
}
