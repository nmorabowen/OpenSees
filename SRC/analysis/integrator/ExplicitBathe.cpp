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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024

#include <ExplicitBathe.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>

#define OPS_Export

// LAPACK eigenvalue solver declarations
#ifdef _WIN32
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);
#else
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);
#endif

// OPS interface for creating the ExplicitBathe integrator
//
// Usage: integrator ExplicitBathe <p> <computeCriticalTimestep>
//                                 <-verbose> <-cfl> <-cflAbort>
//                                 <-divergence factor>
//   p (optional, default 0.54)         sub-step parameter, 0<p<1
//   computeCriticalTimestep (optional) 1 to estimate dt_cr (default 0)
//   -cfl                               alias to enable dt_cr estimation
//   -verbose                           per-step dt/energy reporting (default off)
//   -cflAbort                          stop if dt exceeds the Noh-Bathe limit
//   -divergence <f>                    stop if kinetic energy grows by factor f
//                                      in one step (spurious-energy guard)
//   -tangent                           estimate dt_cr from the current tangent
//                                      stiffness (default: initial stiffness)
//   -recompute <N>                     refresh dt_cr every N steps (implies -tangent)
//   -lump <rowsum|diagonal>            element-mass lumping for dt_cr (default
//                                      rowsum; diagonal = diagonal-of-consistent,
//                                      better for rotational DOFs)
void *OPS_ExplicitBathe(void) {
    double p = 0.54;            // default: good high-frequency dissipation
    int compute_critical_timestep = 0;
    bool verbose = false;
    bool cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    CTSLumping lumping = CTSLumping::RowSum;

    // p is the leading numeric positional. Read it with the typed getter so it
    // works under BOTH Tcl and OpenSeesPy (OPS_GetString mis-reads a numeric
    // Python argument).
    if (OPS_GetNumRemainingInputArgs() > 0) {
        int nd = 1;
        if (OPS_GetDoubleInput(&nd, &p) < 0) {
            opserr << "WARNING ExplicitBathe - could not read p (give p first, "
                      "e.g. integrator ExplicitBathe 0.54)\n";
            return 0;
        }
    }
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
                opserr << "WARNING ExplicitBathe -recompute expects a positive integer N "
                          "(steps between dt_cr refreshes); dt_cr will be computed once\n";
            cflUseTangent = true;            // recomputing K0 every N steps is pointless
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else opserr << "WARNING ExplicitBathe - unknown -lump " << m
                            << " (use rowsum|diagonal; keeping rowsum)\n";
            }
        } else {
            opserr << "WARNING ExplicitBathe - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBathe - p must be in (0,1); typical 0.50-0.95, "
                  "recommended 0.54 (got " << p << ")\n";
        return 0;
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBathe(p, compute_critical_timestep, verbose, cflAbort,
                          divergenceFactor, cflUseTangent, cflRecomputeEvery, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    return theIntegrator;
}

// Default constructor
ExplicitBathe::ExplicitBathe()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), 
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), R_tdt(0),
      updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(0),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(false), cflAbort(false), divergenceFactor(0.0),
      prevKE(0.0), firstStep(true),
      cflUseTangent(false), cflRecomputeEvery(0), cflStepCount(0),
      cflFirstComputation(true), lumping(CTSLumping::RowSum)
{}

// Main constructor with parameters
//
// Integration coefficients q0, q1, q2 are computed from p:
// q1 = (1 - 2p) / (2p(1-p))
// q2 = 0.5 - p*q1
// q0 = -q1 - q2 + 0.5
ExplicitBathe::ExplicitBathe(double _p, int compute_critical_timestep_,
                             bool verbose_, bool cflAbort_, double divergenceFactor_,
                             bool cflUseTangent_, int cflRecomputeEvery_,
                             CTSLumping lumping_)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), R_tdt(0),
      updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(verbose_), cflAbort(cflAbort_), divergenceFactor(divergenceFactor_),
      prevKE(0.0), firstStep(true),
      cflUseTangent(cflUseTangent_), cflRecomputeEvery(cflRecomputeEvery_),
      cflStepCount(0), cflFirstComputation(true), lumping(lumping_)
{
    // Calculate integration coefficients from p parameter
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    if (verbose)
        opserr << "ExplicitBathe: p = " << p
               << ", compute_critical_timestep = " << compute_critical_timestep << endln;
}

// Destructor - clean up allocated memory
ExplicitBathe::~ExplicitBathe() {
    if (U_t) delete U_t;
    if (V_t) delete V_t;
    if (A_t) delete A_t;
    if (U_tpdt) delete U_tpdt;
    if (V_tpdt) delete V_tpdt;
    if (V_fake) delete V_fake;
    if (A_tpdt) delete A_tpdt;
    if (U_tdt) delete U_tdt;
    if (V_tdt) delete V_tdt;
    if (A_tdt) delete A_tdt;
    if (R_tdt) delete R_tdt;
}

// Handle domain changes (mesh modifications, initial conditions, etc.)
int ExplicitBathe::domainChanged() {
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBathe::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    if (size == 0) {
        opserr << "ExplicitBathe::domainChanged - invalid size\n";
        return -1;
    }

    // Allocate or reallocate state vectors if size has changed
    if (!U_t || U_t->Size() != size) {
        // Delete old vectors
        if (U_t) delete U_t;
        if (V_t) delete V_t;
        if (A_t) delete A_t;
        if (U_tpdt) delete U_tpdt;
        if (V_tpdt) delete V_tpdt;
        if (V_fake) delete V_fake;
        if (A_tpdt) delete A_tpdt;
        if (U_tdt) delete U_tdt;
        if (V_tdt) delete V_tdt;
        if (A_tdt) delete A_tdt;

        // Create new vectors
        U_t = new Vector(size);
        V_t = new Vector(size);
        A_t = new Vector(size);
        U_tpdt = new Vector(size);
        V_tpdt = new Vector(size);
        V_fake = new Vector(size);
        A_tpdt = new Vector(size);
        U_tdt = new Vector(size);
        V_tdt = new Vector(size);
        A_tdt = new Vector(size);

        // Verify allocation was successful
        if (!U_t || !V_t || !A_t || !U_tpdt || !V_tpdt || !V_fake ||
            !A_tpdt || !U_tdt || !V_tdt || !A_tdt) {
            opserr << "ExplicitBathe::domainChanged - out of memory\n";
            return -1;
        }
    }

    // Initialize state vectors from committed DOF values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != nullptr) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        // Initialize displacements
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*U_t)(loc) = disp(i);
                (*U_tpdt)(loc) = disp(i);
                (*U_tdt)(loc) = disp(i);
            }
        }

        // Initialize velocities
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*V_t)(loc) = vel(i);
            }
        }

        // Initialize accelerations
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*A_t)(loc) = accel(i);
            }
        }
    }

    // Initialize critical timestep values
    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    // Reset computation flag if it was already computed; also reset the
    // periodic-recompute counter so the cadence restarts cleanly after a
    // mesh change / staged construction.
    if (compute_critical_timestep == 2) {
        compute_critical_timestep = 1;
    }
    cflStepCount = 0;

    return 0;
}

// The per-element critical-time-step eigensolve now lives in the shared
// SRC/analysis/integrator/CriticalTimeStep.{h,cpp} (pulled in via ExplicitBathe.h),
// shared with ExplicitBatheLNVD instead of a hand-copied `extern`. That version
// also adds DSYGV (with a DGGEV fallback + relative-beta threshold), the -lump
// option, and a guarded cross-rank MPI_MIN reduction. See computeCriticalTimeStep().

// Advance to a new time step
//
// This method implements the first sub-step of the Bathe scheme:
// From time t to time t + p*dt
// 
// u_{t+p*dt} = u_t + p*dt*v_t + (p*dt)^2/2 * a_t
// v_{t+p*dt} = v_t + p*dt*a_t
// 
// Then solve: M*a_{t+p*dt} = R_{t+p*dt} - C*v_{t+p*dt} - K*u_{t+p*dt}
int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;

    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
        return -1;
    }

    // Each step performs exactly two solves; reset here so a failed/retried
    // step cannot poison the next step's update() counter (matches CentralDifference).
    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();

    // Cold-start diagnostic: an explicit scheme started from rest with zero
    // committed initial acceleration produces no motion on the first step if a
    // load is already applied. Warn once so the user runs an initial
    // equilibrium step or supplies a consistent a0.
    if (firstStep) {
        firstStep = false;
        if (A_t->pNorm(0) == 0.0 && V_t->pNorm(0) == 0.0)
            opserr << "ExplicitBathe - NOTE starting from rest with zero initial "
                      "acceleration; if a load is already applied, run one "
                      "equilibrium step (or set the initial acceleration) so a0 is "
                      "consistent.\n";
    }

    // Critical-time-step estimate (per-element generalized eigenproblem).
    // damped_minimum_critical_timestep is the CONSERVATIVE central-difference
    // limit; the Noh-Bathe scheme is stable up to ~EB_NB_STABILITY_FACTOR times it.
    // With -recompute N it is refreshed every N steps from the tangent stiffness
    // so it tracks softening/stiffening in nonlinear runs.
    if (cflRecomputeEvery > 0 && compute_critical_timestep == 2) {
        if (++cflStepCount >= cflRecomputeEvery) {
            compute_critical_timestep = 1;
            cflStepCount = 0;
        }
    }
    if (compute_critical_timestep == 1) {
        // shared eigensolve (DSYGV/DGGEV, relative-beta, chosen lumping, MPI_MIN
        // reduced across ranks). Fresh CTSResult, so no manual reset needed.
        CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
        damped_minimum_critical_timestep   = r.damped_dt;
        undamped_minimum_critical_timestep = r.undamped_dt;
        damped_critical_element_tag        = r.damped_tag;
        undamped_critical_element_tag      = r.undamped_tag;
        compute_critical_timestep = 2;  // mark as computed

        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (cflFirstComputation || verbose) {
            opserr << "ExplicitBathe: critical time step estimate"
                   << (cflUseTangent ? " (tangent)" : "") << "\n"
                   << "  central-difference limit (conservative): "
                   << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n"
                   << "  Noh-Bathe limit (~" << EB_NB_STABILITY_FACTOR << "x): " << dt_nb << "\n";
            if (deltaT > dt_nb)
                opserr << "  WARNING dt = " << deltaT
                       << " exceeds the Noh-Bathe stability limit - expect INSTABILITY.\n";
            else if (deltaT > damped_minimum_critical_timestep)
                opserr << "  note: dt exceeds the conservative limit but is within the "
                          "Noh-Bathe stable range.\n";
        }
        cflFirstComputation = false;
    }

    // Optional hard CFL guard.
    if (cflAbort && damped_minimum_critical_timestep > 0.0) {
        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (deltaT > dt_nb) {
            opserr << "ExplicitBathe::newStep() - ABORT: dt = " << deltaT
                   << " > Noh-Bathe stability limit " << dt_nb << "\n";
            return -2;
        }
    }

    // Optional per-step reporting (off by default).
    if (verbose && damped_minimum_critical_timestep > 0.0) {
        const double f = deltaT / damped_minimum_critical_timestep;
        opserr << "ExplicitBathe::newStep() dt = " << deltaT
               << " dt_cd = " << damped_minimum_critical_timestep
               << " factor = " << f << (f < 1.0 ? " [OK]" : " [> CD limit]") << endln;
    }

    // Compute integration coefficients for this time step
    a0 = p * deltaT;
    a1 = std::pow(p * deltaT, 2) / 2.0;
    a2 = a0 / 2.0;
    a3 = (1.0 - p) * deltaT;
    a4 = std::pow((1.0 - p) * deltaT, 2) / 2.0;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // Predict displacement and velocity at t + p*dt
    *U_tpdt = *U_t;
    U_tpdt->addVector(1.0, *V_t, a0);      // += v_t * p*dt
    U_tpdt->addVector(1.0, *A_t, a1);      // += a_t * (p*dt)^2/2
    
    *V_fake = *V_t;
    V_fake->addVector(1.0, *A_t, a0);      // += a_t * p*dt
    
    A_tpdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + p * deltaT;

    if (theModel->updateDomain(newtime, p * deltaT) < 0) {
        opserr << "ExplicitBathe::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

// Update the response quantities
//
// This method is called twice per time step:
// 1st call: After solving at t + p*dt, updates to prepare for t + dt
// 2nd call: After solving at t + dt, finalizes the time step
int ExplicitBathe::update(const Vector &U) {
    updateCount++;
    if (updateCount > 2) {
        opserr << "WARNING ExplicitBathe::update() - called more than twice in a step.\n";
        opserr << "  ExplicitBathe is an explicit scheme and requires 'algorithm Linear' "
                  "(exactly 2 solves per step).\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitBathe::update() - no AnalysisModel set\n";
        return -2;
    }

    if (U_t == 0) {
        opserr << "WARNING ExplicitBathe::update() - domainChange() failed or not called\n";
        return -3;
    }

    if (U.Size() != A_t->Size()) {
        opserr << "WARNING ExplicitBathe::update() - Vector size mismatch\n";
        return -4;
    }

    LinearSOE *theLinSOE = this->getLinearSOE();
    
    // Store acceleration at t + p*dt
    *A_tpdt = U;

    // Update velocity at t + p*dt (corrected)
    // v_{t+p*dt} = v_t + (a_t + a_{t+p*dt}) * p*dt/2
    *V_tpdt = *V_t;
    V_tpdt->addVector(1.0, *A_t, a2);      // += a_t * p*dt/2
    V_tpdt->addVector(1.0, *A_tpdt, a2);   // += a_{t+p*dt} * p*dt/2

    // Prepare fake velocity for response setting
    *V_fake = *V_tpdt;
    V_fake->addVector(1.0, *A_tpdt, a3);   // += a_{t+p*dt} * (1-p)*dt

    // Predict displacement at t + dt
    // u_{t+dt} = u_{t+p*dt} + v_{t+p*dt} * (1-p)*dt + a_{t+p*dt} * ((1-p)*dt)^2/2
    *U_tdt = *U_tpdt;
    U_tdt->addVector(1.0, *V_tpdt, a3);    // += v_{t+p*dt} * (1-p)*dt
    U_tdt->addVector(1.0, *A_tpdt, a4);    // += a_{t+p*dt} * ((1-p)*dt)^2/2

    A_tdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + (1.0 - p) * deltaT;

    if (theModel->updateDomain(newtime, (1.0 - p) * deltaT) < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -3;
    }

    // Solve for acceleration at t + dt
    this->formUnbalance();
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();

    // Circuit breaker 1: catch a blown-up (NaN/Inf) acceleration before it
    // silently propagates through the rest of the run.
    const double A_max = A_tdt->pNorm(0);
    if (A_max != A_max || A_max == std::numeric_limits<double>::infinity()) {
        opserr << "ExplicitBathe::update() - ABORT: non-finite acceleration "
                  "(NaN/Inf) - the integration has diverged (dt likely too large).\n";
        return -5;
    }

    // Circuit breaker 2 (opt-in): abort on runaway kinetic energy growth.
    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*V_tpdt) ^ (*V_tpdt));   // velocity-based KE proxy
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "ExplicitBathe::update() - ABORT: kinetic-energy proxy grew by "
                   << (ke / prevKE) << "x in one step (> " << divergenceFactor
                   << ") - spurious energy growth / instability.\n";
            return -6;
        }
        if (ke > 0.0) prevKE = ke;
    }

    if (verbose)
        opserr << "ExplicitBathe::update() max|a| = " << A_max << endln;

    // Final velocity update at t + dt using all three accelerations
    // v_{t+dt} = v_{t+p*dt} + q0*a_t*(1-p)*dt + (0.5+q1)*a_{t+p*dt}*(1-p)*dt + q2*a_{t+dt}*(1-p)*dt
    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);       // += q0 * a_t * (1-p)*dt
    V_tdt->addVector(1.0, *A_tpdt, a6);    // += (0.5+q1) * a_{t+p*dt} * (1-p)*dt
    V_tdt->addVector(1.0, *A_tdt, a7);     // += q2 * a_{t+dt} * (1-p)*dt

    // Set final response
    theModel->setResponse(*U_tdt, *V_tdt, *A_tdt);
    
    if (theModel->updateDomain() < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

// Form element tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

// Form nodal tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// Commit the state for this time step
//
// This method is called after both sub-steps are complete.
// It updates the state vectors for the next time step.
int ExplicitBathe::commit() {
    updateCount = 0;  // Reset update counter for next step

    // Update state vectors: t+dt becomes new t
    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBathe::commit() - no AnalysisModel set\n";
        return -1;
    }

    return theModel->commitDomain();
}

// Get current velocity (for modal damping interface)
const Vector &ExplicitBathe::getVel() {
    return *V_t;
}

// Conservative (central-difference) critical time step; valid once a step has
// run with critical-time-step computation enabled. <=0 if not yet computed.
double ExplicitBathe::getCriticalTimeStep(void) const {
    if (compute_critical_timestep == 0)
        return -1.0;
    if (damped_minimum_critical_timestep == std::numeric_limits<double>::infinity())
        return -1.0;
    return damped_minimum_critical_timestep;
}

// Send object state for parallel processing
int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(1);
    data(0) = p;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

// Receive object state for parallel processing
int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);

    // Recalculate integration coefficients from received p
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    return 0;
}

// Print integrator information
void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  Damping parameter p: " << p << "\n";
    stream << "  Integration coefficients: q0 = " << q0 
           << ", q1 = " << q1 << ", q2 = " << q2 << "\n";
    if (compute_critical_timestep > 0) {
        stream << "  Critical timestep (damped): " << damped_minimum_critical_timestep << "\n";
        stream << "  Critical timestep (undamped): " << undamped_minimum_critical_timestep << "\n";
    }
}
