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

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell
//
// See CentralDifferenceExplicit.h for the scheme and the references.

#include <CentralDifferenceExplicit.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Domain.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <classTags.h>
#include <elementAPI.h>
#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

// A damped limit below this fraction of the undamped limit prints the
// stiffness-proportional damping note. Mass-proportional damping barely moves
// the ratio at omega_max, because xi = alphaM/(2 omega) is small there.
#define CDE_BETAK_WARN_RATIO 0.90

void *OPS_CentralDifferenceExplicit(void)
{
    // integrator CentralDifferenceExplicit <-cfl> <-cflAbort> <-tangent>
    //     <-recompute N> <-lump rowsum|diagonal> <-verbose> <-divergence f>
    //     <-noMassCache> <-commitSolveState>
    int compute_critical_timestep = 0;
    bool verbose = false;
    bool cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    bool massCache = true;
    bool commitSolve = false;
    CTSLumping lumping = CTSLumping::Diagonal;

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
                opserr << "WARNING CentralDifferenceExplicit -recompute expects a positive "
                          "integer N (steps between dt_cr refreshes); dt_cr computed once\n";
            cflUseTangent = true;   // recomputing the initial stiffness has no effect
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-noMassCache") == 0) {
            massCache = false;
        } else if (strcmp(arg, "-commitSolveState") == 0) {
            commitSolve = true;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else opserr << "WARNING CentralDifferenceExplicit - unknown -lump " << m
                            << " (use rowsum|diagonal; keeping diagonal)\n";
            }
        } else {
            opserr << "WARNING CentralDifferenceExplicit - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    CentralDifferenceExplicit *theIntegrator =
        new CentralDifferenceExplicit(compute_critical_timestep, verbose, cflAbort,
                                      divergenceFactor, cflUseTangent,
                                      cflRecomputeEvery, lumping);
    theIntegrator->setMassCacheEnabled(massCache);
    theIntegrator->setCommitSolveState(commitSolve);
    return theIntegrator;
}

CentralDifferenceExplicit::CentralDifferenceExplicit()
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceExplicit),
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
      useMassCache(true), massTangentValid(false), massCacheNoted(false),
      commitSolveState(false), commitSolveStateNoted(false)
{
}

CentralDifferenceExplicit::CentralDifferenceExplicit(
    int compute_critical_timestep_, bool verbose_, bool cflAbort_,
    double divergenceFactor_, bool cflUseTangent_, int cflRecomputeEvery_,
    CTSLumping lumping_)
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceExplicit),
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
      useMassCache(true), massTangentValid(false), massCacheNoted(false),
      commitSolveState(false), commitSolveStateNoted(false)
{
}

CentralDifferenceExplicit::~CentralDifferenceExplicit()
{
    if (Ut    != 0) delete Ut;
    if (Vhalf != 0) delete Vhalf;
    if (Aprev != 0) delete Aprev;
    if (Vfull != 0) delete Vfull;
    if (Azero != 0) delete Azero;
}

// The tangent is pure mass and getCFactor() is 0, so the modal damping matrix
// branch of TransientIntegrator::formTangent contributes nothing. Skipping the
// reassembly leaves the factored diagonal in place.
int CentralDifferenceExplicit::formTangent(int statFlag)
{
    if (useMassCache && massTangentValid) {
        statusFlag = statFlag;
        return 0;
    }
    int res = this->TransientIntegrator::formTangent(statFlag);
    if (res == 0 && useMassCache) {
        massTangentValid = true;
        if (!massCacheNoted) {
            massCacheNoted = true;
            AnalysisModel *theModel = this->getAnalysisModel();
            Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
            if (theDomain != 0 && theDomain->getNumParameters() > 0)
                opserr << "CentralDifferenceExplicit: the mass tangent is reused between "
                          "domain changes and the model declares parameters. An "
                          "updateParameter that changes a density during the run needs "
                          "-noMassCache.\n";
        }
    }
    return res;
}

int CentralDifferenceExplicit::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

int CentralDifferenceExplicit::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

static inline double cde_reported(double damped, double undamped)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (undamped == inf || undamped <= 0.0)
        return -1.0;
    if (damped != inf && damped > 0.0 && damped < undamped)
        return damped;
    return undamped;
}

int CentralDifferenceExplicit::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (theModel == 0 || theLinSOE == 0) {
        opserr << "CentralDifferenceExplicit::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    if (size == 0) {
        opserr << "CentralDifferenceExplicit::domainChanged - invalid size\n";
        return -1;
    }

    if (Ut == 0 || Ut->Size() != size) {
        if (Ut    != 0) delete Ut;
        if (Vhalf != 0) delete Vhalf;
        if (Aprev != 0) delete Aprev;
        if (Vfull != 0) delete Vfull;
        if (Azero != 0) delete Azero;

        Ut    = new Vector(size);
        Vhalf = new Vector(size);
        Aprev = new Vector(size);
        Vfull = new Vector(size);
        Azero = new Vector(size);
    }
    Azero->Zero();

    // Vhalf holds v_0 until the first-step starter centers it.
    this->seedFromCommitted(theModel);

    firstStep = true;
    cflStepCount = 0;

    // The kinetic energy baseline does not survive a change of mesh.
    prevKE = 0.0;
    committedPrevKE = 0.0;

    massTangentValid = false;

    // The estimate has its own eigensolve and needs neither dt nor a factored
    // system, so the value is available before the first analyze. The -cfl flag
    // controls reporting, -cflAbort and -recompute only.
    CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
    damped_minimum_critical_timestep   = r.damped_dt;
    undamped_minimum_critical_timestep = r.undamped_dt;
    damped_critical_element_tag        = r.damped_tag;
    undamped_critical_element_tag      = r.undamped_tag;
    if (compute_critical_timestep == 2)
        compute_critical_timestep = 1;

    const double inf = std::numeric_limits<double>::infinity();
    const double dt_cr = cde_reported(damped_minimum_critical_timestep,
                                      undamped_minimum_critical_timestep);

    if (!betaKWarned &&
        damped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep > 0.0 &&
        damped_minimum_critical_timestep <
            CDE_BETAK_WARN_RATIO * undamped_minimum_critical_timestep) {
        opserr << "CentralDifferenceExplicit - NOTE damping is significantly reducing the "
                  "critical time step (damped dt_cr = " << damped_minimum_critical_timestep
               << " vs undamped " << undamped_minimum_critical_timestep << ").\n"
                  "  Stiffness-proportional Rayleigh damping scales dt_cr with "
                  "2/(betaK*omega_max^2); mass-proportional alphaM is benign at high "
                  "frequency. Reporting the damped dt_cr.\n";
        betaKWarned = true;
    }

    if (compute_critical_timestep || verbose) {
        opserr << "CentralDifferenceExplicit: critical time step estimate"
               << (cflUseTangent ? " (tangent)" : "") << "\n"
               << "  central-difference limit (2/omega_max): "
               << undamped_minimum_critical_timestep
               << " @ element #" << undamped_critical_element_tag << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            opserr << "  damped limit (reported): " << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n";
        if (dt_cr <= 0.0)
            opserr << "  note: no element produced a finite, positive estimate.\n";
    }

    return 0;
}

void CentralDifferenceExplicit::seedFromCommitted(AnalysisModel *theModel)
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

int CentralDifferenceExplicit::revertToLastStep(void)
{
    if (Ut == 0)
        return 0;
    massTangentValid = false;
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)
        return 0;

    this->seedFromCommitted(theModel);

    firstStep = true;
    updateCount = 0;
    prevKE = committedPrevKE;
    cflStepCount = committedCflStepCount;
    return 0;
}

int CentralDifferenceExplicit::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (deltaT <= 0.0) {
        opserr << "CentralDifferenceExplicit::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -1;
    }

    if (Ut == 0 || Vhalf == 0 || Aprev == 0) {
        opserr << "CentralDifferenceExplicit::newStep() - domainChanged() failed or not called\n";
        return -2;
    }

    committedPrevKE       = prevKE;
    committedCflStepCount = cflStepCount;

    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

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

        const double dt_cr = cde_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (cflFirstComputation || verbose) {
            opserr << "CentralDifferenceExplicit: critical time step (2/omega_max)"
                   << (cflUseTangent ? " (tangent)" : "")
                   << ": " << dt_cr
                   << " @ element #"
                   << ((dt_cr == damped_minimum_critical_timestep)
                           ? damped_critical_element_tag
                           : undamped_critical_element_tag)
                   << "\n";
            if (dt_cr > 0.0 && deltaT > dt_cr)
                opserr << "  WARNING dt = " << deltaT
                       << " exceeds the central-difference stability limit " << dt_cr
                       << " - expect INSTABILITY.\n";
        }
        cflFirstComputation = false;
    }

    if (cflAbort) {
        const double dt_cr = cde_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (dt_cr > 0.0 && deltaT > dt_cr) {
            opserr << "CentralDifferenceExplicit::newStep() - ABORT: dt = " << deltaT
                   << " > central-difference stability limit " << dt_cr << "\n";
            return -2;
        }
    }

    // First-step starter: a_0 at the committed configuration, v_{-1/2} = v_0 - dt/2 a_0.
    if (firstStep) {
        firstStep = false;

        theModel->setResponse(*Ut, *Vhalf, *Azero);
        double t0 = theModel->getCurrentDomainTime();
        if (theModel->updateDomain(t0, deltaT) < 0) {
            opserr << "CentralDifferenceExplicit::newStep() - starter: failed to update domain\n";
            return -3;
        }
        if (this->formTangent(CURRENT_TANGENT) < 0) {
            opserr << "CentralDifferenceExplicit::newStep() - starter: formTangent failed\n";
            return -3;
        }
        if (this->formUnbalance() < 0) {
            opserr << "CentralDifferenceExplicit::newStep() - starter: formUnbalance failed\n";
            return -3;
        }
        if (theLinSOE->solve() < 0) {
            opserr << "CentralDifferenceExplicit::newStep() - starter: SOE solve failed\n";
            return -3;
        }
        *Aprev = theLinSOE->getX();
        Vhalf->addVector(1.0, *Aprev, -0.5 * deltaT);
    }

    // v_{n+1/2} = v_{n-1/2} + dt a_n ;   u_{n+1} = u_n + dt v_{n+1/2}
    Vhalf->addVector(1.0, *Aprev, deltaT);
    Ut->addVector(1.0, *Vhalf, deltaT);

    // Trial state of the mass-only solve. Acceleration 0: no inertia in the
    // residual. Damping sees the lagged v_{n+1/2}.
    theModel->setResponse(*Ut, *Vhalf, *Azero);

    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0) {
        opserr << "CentralDifferenceExplicit::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

int CentralDifferenceExplicit::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1) {
        opserr << "WARNING CentralDifferenceExplicit::update() - called more than once -";
        opserr << " the central-difference scheme requires a LINEAR solution algorithm "
                  "(exactly one solve/step)\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceExplicit::update() - no AnalysisModel set\n";
        return -2;
    }
    if (Ut == 0) {
        opserr << "WARNING CentralDifferenceExplicit::update() - domainChanged() failed or not called\n";
        return -3;
    }
    if (U.Size() != Ut->Size()) {
        opserr << "WARNING CentralDifferenceExplicit::update() - Vectors of incompatible size ";
        opserr << " expecting " << Ut->Size() << " obtained " << U.Size() << endln;
        return -4;
    }

    // U is the solved a_{n+1}. Vector::pNorm(0) cannot serve as the guard: its
    // maximum compares with `>`, which is false for NaN, so NaN entries are skipped.
    if (!vectorIsFinite(U)) {
        opserr << "CentralDifferenceExplicit::update() - ABORT: non-finite acceleration "
                  "(NaN/Inf) - the integration has diverged (dt likely too large).\n";
        return -5;
    }

    // v_{n+1} = v_{n+1/2} + dt/2 a_{n+1}, the centered (u_{n+2}-u_n)/(2 dt).
    *Vfull = *Vhalf;
    Vfull->addVector(1.0, U, 0.5 * deltaT);

    // Required under -commitSolveState as well: committed nodal output reads this.
    theModel->setResponse(*Ut, *Vfull, U);
    if (!commitSolveState) {
        if (theModel->updateDomain() < 0) {
            opserr << "CentralDifferenceExplicit::update() - failed to update the domain\n";
            return -6;
        }
    } else if (!commitSolveStateNoted) {
        commitSolveStateNoted = true;
        opserr << "CentralDifferenceExplicit: -commitSolveState - element state is committed "
                  "at the half-step velocity v_{n+1/2}. Bit-identical for rate-independent "
                  "materials; rate-dependent materials with memory commit a v_{n+1/2}-rate "
                  "state.\n";
    }

    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*Vfull) ^ (*Vfull));
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "CentralDifferenceExplicit::update() - ABORT: kinetic-energy proxy grew to "
                   << (ke / prevKE) << "x its running maximum (> " << divergenceFactor
                   << ") - spurious energy growth / instability.\n";
            return -7;
        }
        if (ke > prevKE) prevKE = ke;
    }

    if (verbose)
        opserr << "CentralDifferenceExplicit::update() max|a| = " << U.pNorm(0) << endln;

    *Aprev = U;

    return 0;
}

int CentralDifferenceExplicit::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceExplicit::commit() - no AnalysisModel set\n";
        return -1;
    }
    return theModel->commitDomain();
}

const Vector &CentralDifferenceExplicit::getVel(void)
{
    if (Vhalf == 0) {
        opserr << "WARNING CentralDifferenceExplicit::getVel() - called before "
                  "domainChanged(); returning an empty vector\n";
        static Vector emptyVel(0);
        return emptyVel;
    }
    return *Vhalf;
}

double CentralDifferenceExplicit::getCriticalTimeStep(void) const
{
    return cde_reported(damped_minimum_critical_timestep,
                        undamped_minimum_critical_timestep);
}

int CentralDifferenceExplicit::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(9);
    data(0) = (double)compute_critical_timestep;
    data(1) = verbose ? 1.0 : 0.0;
    data(2) = cflAbort ? 1.0 : 0.0;
    data(3) = divergenceFactor;
    data(4) = cflUseTangent ? 1.0 : 0.0;
    data(5) = (double)cflRecomputeEvery;
    data(6) = (double)(int)lumping;
    data(7) = useMassCache ? 1.0 : 0.0;
    data(8) = commitSolveState ? 1.0 : 0.0;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceExplicit::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}

int CentralDifferenceExplicit::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(9);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceExplicit::recvSelf() - could not receive data\n";
        return -1;
    }
    compute_critical_timestep = (int)data(0);
    verbose           = (data(1) != 0.0);
    cflAbort          = (data(2) != 0.0);
    divergenceFactor  = data(3);
    cflUseTangent     = (data(4) != 0.0);
    cflRecomputeEvery = (int)data(5);
    lumping           = ((int)data(6) == 0) ? CTSLumping::RowSum : CTSLumping::Diagonal;
    useMassCache      = (data(7) != 0.0);
    commitSolveState  = (data(8) != 0.0);

    return 0;
}

void CentralDifferenceExplicit::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    s << "CentralDifferenceExplicit - explicit leap-frog central difference\n";
    if (theModel != 0)
        s << "  currentTime: " << theModel->getCurrentDomainTime() << endln;
    s << "  recipe: lumped/diagonal mass, system Diagonal, algorithm Linear, dt < dt_cr\n";
    if (compute_critical_timestep > 0) {
        s << "  critical dt (2/omega_max): " << undamped_minimum_critical_timestep << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            s << "  critical dt (damped, reported): " << damped_minimum_critical_timestep << "\n";
    }
}
