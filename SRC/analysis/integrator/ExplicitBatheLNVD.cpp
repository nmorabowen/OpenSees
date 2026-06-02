/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA), 12/2024
// Local-damping wired/activated + resync: nmb 05/2026
//
// See ExplicitBatheLNVD.h for the scheme and the FLAC local-damping definition.

#include <ExplicitBatheLNVD.h>
#include <ExplicitBathe.h>           // EB_NB_STABILITY_FACTOR
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <Vector.h>
#include <ID.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

// The per-element critical-time-step eigensolve is shared with ExplicitBathe and
// lives in CriticalTimeStep.{h,cpp} (pulled in via ExplicitBatheLNVD.h). This
// replaces the former hand-copied `extern`, which could silently drift from the
// definition.

// OPS interface for creating the ExplicitBatheLNVD integrator
//
// Usage: integrator ExplicitBatheLNVD <p> <alpha> <computeCriticalTimestep>
//                                     <-verbose> <-cfl> <-cflAbort>
//                                     <-divergence factor>
//   p     (optional, default 0.54)  sub-step parameter, 0<p<1
//   alpha (optional, default 0.80)  FLAC local damping, 0<=alpha<1
void *OPS_ExplicitBatheLNVD(void) {
    double p = 0.54;
    double alpha_flac = 0.8;          // classic FLAC default
    int compute_critical_timestep = 0;
    bool verbose = false, cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    CTSLumping lumping = CTSLumping::RowSum;

    // p and alpha are the leading numeric positionals; read them with the typed
    // getter (Tcl- and OpenSeesPy-safe). Flags follow.
    int navail = OPS_GetNumRemainingInputArgs();
    if (navail >= 2) {
        int nd = 2; double vals[2];
        if (OPS_GetDoubleInput(&nd, vals) < 0) {
            opserr << "WARNING ExplicitBatheLNVD - could not read p and alpha "
                      "(e.g. integrator ExplicitBatheLNVD 0.54 0.8)\n";
            return 0;
        }
        p = vals[0]; alpha_flac = vals[1];
    } else if (navail == 1) {
        int nd = 1;
        if (OPS_GetDoubleInput(&nd, &p) < 0) {
            opserr << "WARNING ExplicitBatheLNVD - could not read p\n";
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
                opserr << "WARNING ExplicitBatheLNVD -recompute expects a positive integer N "
                          "(steps between dt_cr refreshes); dt_cr will be computed once\n";
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else opserr << "WARNING ExplicitBatheLNVD - unknown -lump " << m
                            << " (use rowsum|diagonal; keeping rowsum)\n";
            }
        } else {
            opserr << "WARNING ExplicitBatheLNVD - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVD - p must be in (0,1), recommended 0.54 "
                  "(got " << p << ")\n";
        return 0;
    }
    if (alpha_flac < 0.0 || alpha_flac >= 1.0) {
        opserr << "WARNING ExplicitBatheLNVD - alpha must be in [0,1); FLAC default "
                  "0.8 (got " << alpha_flac << ")\n";
        return 0;
    }

    TransientIntegrator *theIntegrator =
        new ExplicitBatheLNVD(p, alpha_flac, compute_critical_timestep,
                              verbose, cflAbort, divergenceFactor,
                              cflUseTangent, cflRecomputeEvery, lumping);
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitBatheLNVD integrator\n";
    return theIntegrator;
}


ExplicitBatheLNVD::ExplicitBatheLNVD()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBatheLNVD),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), alpha_flac(0.0),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(0),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(false), cflAbort(false), divergenceFactor(0.0),
      prevKE(0.0), firstStep(true), lastUnbalanceNorm(-1.0),
      cflUseTangent(false), cflRecomputeEvery(0), cflStepCount(0),
      cflFirstComputation(true), lumping(CTSLumping::RowSum)
{}

ExplicitBatheLNVD::ExplicitBatheLNVD(double _p, double _alpha_flac,
                                     int compute_critical_timestep_,
                                     bool verbose_, bool cflAbort_,
                                     double divergenceFactor_,
                                     bool cflUseTangent_, int cflRecomputeEvery_,
                                     CTSLumping lumping_)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBatheLNVD),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0), alpha_flac(_alpha_flac),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(verbose_), cflAbort(cflAbort_), divergenceFactor(divergenceFactor_),
      prevKE(0.0), firstStep(true), lastUnbalanceNorm(-1.0),
      cflUseTangent(cflUseTangent_), cflRecomputeEvery(cflRecomputeEvery_),
      cflStepCount(0), cflFirstComputation(true), lumping(lumping_)
{
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    if (verbose)
        opserr << "ExplicitBatheLNVD - p = " << p << ", alpha = " << alpha_flac << endln;
}

ExplicitBatheLNVD::~ExplicitBatheLNVD() {
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
}

int ExplicitBatheLNVD::domainChanged() {
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBatheLNVD::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    if (size == 0) {
        opserr << "ExplicitBatheLNVD::domainChanged - invalid size\n";
        return -1;
    }

    if (!U_t || U_t->Size() != size) {
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

        if (!U_t || !V_t || !A_t || !U_tpdt || !V_tpdt || !V_fake ||
            !A_tpdt || !U_tdt || !V_tdt || !A_tdt) {
            opserr << "ExplicitBatheLNVD::domainChanged - out of memory\n";
            return -1;
        }
    }

    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*U_t)(loc) = disp(i);
                (*U_tpdt)(loc) = disp(i);
                (*U_tdt)(loc) = disp(i);
            }
        }
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) (*V_t)(loc) = vel(i);
        }
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) (*A_t)(loc) = accel(i);
        }
    }

    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    if (compute_critical_timestep == 2)
        compute_critical_timestep = 1;
    cflStepCount = 0;   // restart periodic-recompute cadence after a mesh change

    return 0;
}

int ExplicitBatheLNVD::newStep(double _deltaT) {
    deltaT = _deltaT;

    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBatheLNVD::newStep() - state variables not initialized\n";
        return -1;
    }

    updateCount = 0;   // exactly two solves per step; reset for robustness

    AnalysisModel *theModel = this->getAnalysisModel();

    if (firstStep) {
        firstStep = false;
        if (A_t->pNorm(0) == 0.0 && V_t->pNorm(0) == 0.0 && verbose)
            opserr << "ExplicitBatheLNVD - NOTE starting from rest with zero initial "
                      "acceleration.\n";
    }

    // critical time step (conservative central-difference limit); with
    // -recompute N it is refreshed every N steps from the tangent stiffness
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
        compute_critical_timestep = 2;

        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (cflFirstComputation || verbose) {
            opserr << "ExplicitBatheLNVD: critical time step estimate"
                   << (cflUseTangent ? " (tangent)" : "") << "\n"
                   << "  central-difference limit (conservative): "
                   << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n"
                   << "  Noh-Bathe limit (~" << EB_NB_STABILITY_FACTOR << "x): " << dt_nb << "\n";
        }
        cflFirstComputation = false;
    }

    if (cflAbort && damped_minimum_critical_timestep > 0.0) {
        const double dt_nb = EB_NB_STABILITY_FACTOR * damped_minimum_critical_timestep;
        if (deltaT > dt_nb) {
            opserr << "ExplicitBatheLNVD::newStep() - ABORT: dt = " << deltaT
                   << " > Noh-Bathe stability limit " << dt_nb << "\n";
            return -2;
        }
    }

    // integration coefficients
    a0 = p * deltaT;
    a1 = (p * deltaT) * (p * deltaT) / 2.0;
    a2 = a0 / 2.0;
    a3 = (1.0 - p) * deltaT;
    a4 = ((1.0 - p) * deltaT) * ((1.0 - p) * deltaT) / 2.0;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // predictor for sub-step 1 (t -> t + p*dt); the solution algorithm then
    // calls formUnbalance()/solve()/update() (sub-step-1 solve sees the damping)
    *U_tpdt = *U_t;
    U_tpdt->addVector(1.0, *V_t, a0);
    U_tpdt->addVector(1.0, *A_t, a1);
    *V_fake = *V_t;
    V_fake->addVector(1.0, *A_t, a0);
    A_tpdt->Zero();

    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);

    double oldtime = theModel->getCurrentDomainTime();
    if (theModel->updateDomain(oldtime + p*deltaT, p*deltaT) < 0) {
        opserr << "ExplicitBatheLNVD::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

int ExplicitBatheLNVD::update(const Vector &U) {
    updateCount++;
    if (updateCount > 2) {
        opserr << "WARNING ExplicitBatheLNVD::update() - called more than twice in a step.\n";
        opserr << "  ExplicitBatheLNVD is explicit and requires 'algorithm Linear'.\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitBatheLNVD::update() - no AnalysisModel set\n";
        return -2;
    }
    if (U_t == 0) {
        opserr << "WARNING ExplicitBatheLNVD::update() - domainChanged() not called\n";
        return -3;
    }
    if (U.Size() != A_t->Size()) {
        opserr << "WARNING ExplicitBatheLNVD::update() - vector size mismatch\n";
        return -4;
    }

    LinearSOE *theLinSOE = this->getLinearSOE();

    // store acceleration at t + p*dt
    *A_tpdt = U;

    // velocity at t + p*dt
    *V_tpdt = *V_t;
    V_tpdt->addVector(1.0, *A_t, a2);
    V_tpdt->addVector(1.0, *A_tpdt, a2);

    *V_fake = *V_tpdt;
    V_fake->addVector(1.0, *A_tpdt, a3);

    // displacement at t + dt
    *U_tdt = *U_tpdt;
    U_tdt->addVector(1.0, *V_tpdt, a3);
    U_tdt->addVector(1.0, *A_tpdt, a4);

    A_tdt->Zero();
    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    double oldtime = theModel->getCurrentDomainTime();
    if (theModel->updateDomain(oldtime + (1.0-p)*deltaT, (1.0-p)*deltaT) < 0) {
        opserr << "ExplicitBatheLNVD::update() - failed to update the domain\n";
        return -3;
    }

    // sub-step-2 solve (formUnbalance() adds the FLAC local damping)
    this->formUnbalance();
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();

    const double A_max = A_tdt->pNorm(0);
    if (A_max != A_max || A_max == std::numeric_limits<double>::infinity()) {
        opserr << "ExplicitBatheLNVD::update() - ABORT: non-finite acceleration (NaN/Inf).\n";
        return -5;
    }
    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*V_tpdt) ^ (*V_tpdt));
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "ExplicitBatheLNVD::update() - ABORT: kinetic energy grew "
                   << (ke / prevKE) << "x in one step (> " << divergenceFactor << ").\n";
            return -6;
        }
        if (ke > 0.0) prevKE = ke;
    }
    if (verbose)
        opserr << "ExplicitBatheLNVD::update() max|a| = " << A_max
               << "  |unbalance| = " << lastUnbalanceNorm << endln;

    // velocity at t + dt
    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);
    V_tdt->addVector(1.0, *A_tpdt, a6);
    V_tdt->addVector(1.0, *A_tdt, a7);

    theModel->setResponse(*U_tdt, *V_tdt, *A_tdt);
    if (theModel->updateDomain() < 0) {
        opserr << "ExplicitBatheLNVD::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

// Assemble the residual and add FLAC local non-viscous damping. Called by the
// solution algorithm for sub-step 1 and by update() for sub-step 2, so the
// damping is applied symmetrically to both Noh-Bathe sub-steps.
int ExplicitBatheLNVD::formUnbalance(void) {
    int res = this->IncrementalIntegrator::formUnbalance();
    if (res < 0) return res;

    LinearSOE *theSOE = this->getLinearSOE();
    if (theSOE != 0)
        lastUnbalanceNorm = theSOE->getB().pNorm(0);

    if (alpha_flac > 0.0)
        this->addLocalDamping();

    return res;
}

void ExplicitBatheLNVD::addLocalDamping(void) {
    LinearSOE *theSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theSOE == 0 || theModel == 0) return;

    const Vector &B = theSOE->getB();
    int n = B.Size();
    if (n == 0) return;

    // gather the global trial velocity at the instant being solved
    Vector vglob(n);
    vglob.Zero();
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        const Vector &v = dofPtr->getTrialVel();
        for (int i = 0; i < id.Size(); ++i) {
            int loc = id(i);
            if (loc >= 0 && loc < n)
                vglob(loc) = v(i);
        }
    }

    // r_i <- r_i - alpha*|r_i|*sign(v_i)  (opposes motion, vanishes at v->0)
    Vector damped(B);
    for (int i = 0; i < n; ++i) {
        double v = vglob(i);
        double sgn = (v > 0.0) ? 1.0 : ((v < 0.0) ? -1.0 : 0.0);
        damped(i) = B(i) - alpha_flac * std::abs(B(i)) * sgn;
    }
    theSOE->setB(damped);
}

int ExplicitBatheLNVD::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

int ExplicitBatheLNVD::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

int ExplicitBatheLNVD::commit() {
    updateCount = 0;
    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "ExplicitBatheLNVD::commit() - no AnalysisModel set\n";
        return -1;
    }
    return theModel->commitDomain();
}

const Vector &ExplicitBatheLNVD::getVel() {
    return *V_t;
}

double ExplicitBatheLNVD::getCriticalTimeStep(void) const {
    if (compute_critical_timestep == 0)
        return -1.0;
    if (damped_minimum_critical_timestep == std::numeric_limits<double>::infinity())
        return -1.0;
    return damped_minimum_critical_timestep;
}

double ExplicitBatheLNVD::getUnbalanceNorm(void) const {
    return lastUnbalanceNorm;
}

int ExplicitBatheLNVD::sendSelf(int cTag, Channel &theChannel) {
    Vector data(2);
    data(0) = p;
    data(1) = alpha_flac;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVD::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}

int ExplicitBatheLNVD::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(2);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVD::recvSelf() - could not receive data\n";
        return -1;
    }
    p = data(0);
    alpha_flac = data(1);
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;
    return 0;
}

void ExplicitBatheLNVD::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method with FLAC Local Non-Viscous Damping\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  p: " << p << ", q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << "\n";
    stream << "  alpha (FLAC local damping): " << alpha_flac << "\n";
}
