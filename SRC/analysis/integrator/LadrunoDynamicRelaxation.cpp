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

// Ladruno: LadrunoDynamicRelaxation — quasi-static dynamic relaxation
// (TransientIntegrator, classTag 33005). Leap-frog skeleton adapted from
// CentralDifferenceLadruno; fictitious Gershgorin mass + Cundall kinetic damping.
// See Ladruno_implementation/21_ladruno_dynamic_relaxation_adr.md.

#include <LadrunoDynamicRelaxation.h>
#include <LadrunoFictitiousMass.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Element.h>
#include <classTags.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>
#include <limits>

void *OPS_LadrunoDynamicRelaxation(void)
{
    // Usage: integrator LadrunoDynamicRelaxation
    //          <-mass gershgorin|lumped $scale|unity> <-dt $dt>
    //          <-recompute $N> <-interp> <-divergence $f> <-verbose>
    int    massMode = 0;            // gershgorin
    double dtPseudo = 1.0, massScale = 1.0, divergenceFactor = 0.0;
    int    recomputeEvery = 0;
    bool   interp = false, verbose = false;
    int    dampMode = 0;            // 0 kinetic (default) | 1 viscous-critical
    double zetaTarget = 1.0;        // viscous damping ratio (1 = critical)
    bool   autoRefresh = true;      // rebuild M* at KE peaks (gershgorin)

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING LadrunoDynamicRelaxation -mass expects gershgorin|lumped|unity\n";
                return 0;
            }
            const char *m = OPS_GetString();
            if (strcmp(m, "gershgorin") == 0)      massMode = 0;
            else if (strcmp(m, "lumped") == 0) {
                massMode = 1;
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    int nd = 1;
                    if (OPS_GetDoubleInput(&nd, &massScale) < 0) massScale = 1.0;
                }
            } else if (strcmp(m, "unity") == 0)    massMode = 2;
            else opserr << "WARNING LadrunoDynamicRelaxation - unknown -mass " << m
                        << " (keeping gershgorin)\n";
        } else if (strcmp(opt, "-dt") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &dtPseudo) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -dt failed to read value\n";
                return 0;
            }
        } else if (strcmp(opt, "-recompute") == 0) {
            int nd = 1;
            if (OPS_GetIntInput(&nd, &recomputeEvery) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -recompute failed to read N\n";
                return 0;
            }
        } else if (strcmp(opt, "-damping") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING LadrunoDynamicRelaxation -damping expects kinetic|viscous\n";
                return 0;
            }
            const char *d = OPS_GetString();
            if (strcmp(d, "kinetic") == 0)      dampMode = 0;
            else if (strcmp(d, "viscous") == 0) {
                dampMode = 1;
                if (OPS_GetNumRemainingInputArgs() > 0) {   // optional trailing zeta
                    int nd = 1; double z;
                    if (OPS_GetDoubleInput(&nd, &z) == 0) zetaTarget = z;
                }
            } else opserr << "WARNING LadrunoDynamicRelaxation - unknown -damping " << d
                          << " (keeping kinetic)\n";
        } else if (strcmp(opt, "-noAutoRefresh") == 0) {
            autoRefresh = false;
        } else if (strcmp(opt, "-interp") == 0) {
            interp = true;
        } else if (strcmp(opt, "-divergence") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &divergenceFactor) < 0) divergenceFactor = 0.0;
        } else if (strcmp(opt, "-verbose") == 0) {
            verbose = true;
        } else {
            opserr << "WARNING LadrunoDynamicRelaxation - unknown option " << opt
                   << " (ignored)\n";
        }
    }

    return new LadrunoDynamicRelaxation(massMode, dtPseudo, massScale,
                                        recomputeEvery, interp, divergenceFactor,
                                        verbose, dampMode, zetaTarget, autoRefresh);
}

LadrunoDynamicRelaxation::LadrunoDynamicRelaxation(int massMode_, double dtPseudo_,
                                                   double massScale_, int recomputeEvery_,
                                                   bool interp_, double divergenceFactor_,
                                                   bool verbose_, int dampMode_,
                                                   double zetaTarget_, bool autoRefresh_)
 :TransientIntegrator(INTEGRATOR_TAGS_LadrunoDynamicRelaxation),
  massMode(massMode_), dtPseudo(dtPseudo_), massScale(massScale_),
  recomputeEvery(recomputeEvery_), interp(interp_),
  divergenceFactor(divergenceFactor_), verbose(verbose_),
  Mstar(0), Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
  firstStep(true), updateCount(0), size(0), stepCount(0),
  prevKE(0.0), kineticEnergy(0.0), residualNorm(0.0), deltaT(0.0),
  dampMode(dampMode_), zetaTarget(zetaTarget_), cVisc(0.0), autoRefresh(autoRefresh_)
{

}

LadrunoDynamicRelaxation::~LadrunoDynamicRelaxation()
{
    if (Mstar != 0) delete Mstar;
    if (Ut    != 0) delete Ut;
    if (Vhalf != 0) delete Vhalf;
    if (Aprev != 0) delete Aprev;
    if (Vfull != 0) delete Vfull;
    if (Azero != 0) delete Azero;
}

// ---- matrix-free LHS: assemble the integrator-owned fictitious M* onto the
// diagonal SOE so theSOE->solve() degenerates to a M*^{-1} apply. This is the
// seam that lets DR use an ARTIFICIAL mass instead of the getMass()-bound
// addMtoTang path (ADR §1.2 / §4.2).
int
LadrunoDynamicRelaxation::formTangent(int statFlag)
{
    LinearSOE *theSOE = this->getLinearSOE();
    if (theSOE == 0 || Mstar == 0) {
        opserr << "WARNING LadrunoDynamicRelaxation::formTangent() - no SOE or M* not built\n";
        return -1;
    }
    // poke M* onto the diagonal (zeroFirst => the solve degenerates to M*^{-1})
    Ladruno::addDiagonalToSOE(theSOE, *Mstar, size, 1.0, /*zeroFirst=*/true);
    return 0;
}

int LadrunoDynamicRelaxation::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();   // unused: formTangent() assembles M* directly
    return 0;
}

int LadrunoDynamicRelaxation::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();   // unused: formTangent() assembles M* directly
    return 0;
}

// ---- build the fictitious lumped diagonal M* (ADR decision 2 / §4.2) --------
// gershgorin: m_i = (dt^2/4) * sum_j |K_ij|  (scale-free, omega*dt ~ 2 by
//             construction) — a ONE-TIME element-stiffness probe, SOE-free.
// lumped:     row-sum of the real element mass, * massScale (needs nonzero mass).
// unity:      m_i = 1.
int
LadrunoDynamicRelaxation::buildFictitiousMass(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Mstar == 0) return -1;

    // gershgorin prefactor (mode 0): f = dt^2/4, scale-free so omega*dt ~ 2 by
    // construction. Viscous: grow M* so the DAMPED step stays stable. Critical
    // viscous damping shrinks the explicit safe step by (sqrt(1+z^2)-z); grow M*
    // by its inverse-square so omega*dt stays <= 2. (Passed as dt2Quarter; the
    // shared builder applies it. lumped/unity ignore it.)
    double f = 0.25 * dtPseudo * dtPseudo;
    if (dampMode == 1) {
        double s = sqrt(1.0 + zetaTarget * zetaTarget) - zetaTarget;
        if (s > 0.0) f /= (s * s);
    }
    // build M* (row-sum of K or real mass, or unity) + mmax*1e-8 zero-floor.
    if (Ladruno::buildGershgorinDiagonal(theModel, *Mstar, massMode, massScale, f) < 0)
        return -1;

    // viscous-critical coefficient C* = cVisc*M*. The M* rescale above grew M* by
    // 1/s^2 (s = sqrt(1+z^2)-z), so the RESCALED omega1*dt ~ 2s (not 2). Critical
    // damping for that omega is cVisc = 2*zeta*omega1 = 4*zeta*s/dtPseudo. Using the
    // unrescaled 2/dt here overshoots (dt*cVisc = 4 > 2 => explicit-damping blowup);
    // the s factor keeps dt*cVisc = 4*zeta*s < 2 (=1.66 at zeta=1). gershgorin is the
    // recommended pairing (closed-form omega); unity/lumped is approximate.
    if (dampMode == 1) {
        double s = sqrt(1.0 + zetaTarget * zetaTarget) - zetaTarget;
        cVisc = 4.0 * zetaTarget * s / dtPseudo;
    }

    return 0;
}

int
LadrunoDynamicRelaxation::domainChanged(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    size = x.Size();
    if (size == 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - invalid size\n";
        return -1;
    }

    if (Ut == 0 || Ut->Size() != size) {
        if (Ut    != 0) delete Ut;
        if (Vhalf != 0) delete Vhalf;
        if (Aprev != 0) delete Aprev;
        if (Vfull != 0) delete Vfull;
        if (Azero != 0) delete Azero;
        if (Mstar != 0) delete Mstar;
        Ut    = new Vector(size);
        Vhalf = new Vector(size);
        Aprev = new Vector(size);
        Vfull = new Vector(size);
        Azero = new Vector(size);
        Mstar = new Vector(size);
        if (Ut == 0 || Vhalf == 0 || Aprev == 0 || Vfull == 0 || Azero == 0 ||
            Mstar == 0 || Mstar->Size() != size) {
            opserr << "LadrunoDynamicRelaxation::domainChanged - out of memory\n";
            return -1;
        }
    }
    Azero->Zero();

    // seed u_n, v_{n-1/2}=v_0, a_n=a_0 from the committed DOF state (CDL pattern);
    // the true starter a_0 / v_{-1/2} are computed on the first newStep().
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Ut)(loc) = disp(i); }
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Vhalf)(loc) = vel(i); }
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Aprev)(loc) = accel(i); }
    }

    firstStep = true;
    stepCount = 0;
    prevKE = 0.0;

    // build the fictitious mass (SOE-free element probe)
    if (this->buildFictitiousMass() < 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - failed to build fictitious mass\n";
        return -1;
    }
    if (verbose) {
        double mn = (*Mstar)(0), mx = (*Mstar)(0);
        for (int i = 1; i < size; i++) { if ((*Mstar)(i) < mn) mn = (*Mstar)(i);
                                         if ((*Mstar)(i) > mx) mx = (*Mstar)(i); }
        opserr << "LadrunoDynamicRelaxation: fictitious M* (mode " << massMode
               << ") range [" << mn << ", " << mx << "], dt = " << dtPseudo << endln;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (deltaT <= 0.0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - dt = " << deltaT << " <= 0\n";
        return -1;
    }
    if (Ut == 0 || Vhalf == 0 || Aprev == 0 || Mstar == 0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - domainChanged() not called\n";
        return -2;
    }
    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    // -recompute N: refresh M* from the (softened) tangent and re-zero velocities
    // so the refresh injects no energy (ADR R6 / test DR-10).
    if (recomputeEvery > 0 && !firstStep) {
        if (++stepCount >= recomputeEvery) {
            stepCount = 0;
            this->buildFictitiousMass();
            Vhalf->Zero();
            prevKE = 0.0;
        }
    }

    // ---- first-step starter: a_0 = M*^{-1}(P_0 - F_int(u_0)), v_{-1/2}=v_0-dt/2 a_0
    if (firstStep) {
        firstStep = false;
        theModel->setResponse(*Ut, *Vhalf, *Azero);
        double t0 = theModel->getCurrentDomainTime();
        if (theModel->updateDomain(t0, deltaT) < 0) return -3;
        if (this->formTangent(CURRENT_TANGENT) < 0) return -3;
        if (this->formUnbalance() < 0) return -3;
        if (theLinSOE->solve() < 0) return -3;
        *Aprev = theLinSOE->getX();
        Vhalf->addVector(1.0, *Aprev, -0.5 * deltaT);
    }

    // ---- leap-frog advance: v_{n+1/2}=v_{n-1/2}+dt a_n ; u_{n+1}=u_n+dt v_{n+1/2}
    Vhalf->addVector(1.0, *Aprev, deltaT);
    Ut->addVector(1.0, *Vhalf, deltaT);

    // trial state for the M*-only solve of a_{n+1} (accel = 0, no inertia/damping
    // in the residual -> B = P - F_int = the true static unbalance).
    theModel->setResponse(*Ut, *Vhalf, *Azero);
    double time = theModel->getCurrentDomainTime() + deltaT;
    if (theModel->updateDomain(time, deltaT) < 0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - failed to update domain\n";
        return -3;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1) {
        opserr << "WARNING LadrunoDynamicRelaxation::update() - called >1; DR requires "
                  "`algorithm Linear` (exactly one solve/step)\n";
        return -1;
    }
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Ut == 0) return -2;
    if (U.Size() != size) {
        opserr << "WARNING LadrunoDynamicRelaxation::update() - size mismatch\n";
        return -4;
    }

    // U is the solved acceleration a_{n+1}. NaN/Inf circuit breaker.
    const double A_max = U.pNorm(0);
    if (A_max != A_max || A_max == std::numeric_limits<double>::infinity()) {
        opserr << "LadrunoDynamicRelaxation::update() - ABORT: non-finite acceleration "
                  "(dt likely too large / M* too small)\n";
        return -5;
    }

    // true static residual norm = ||f_ext - f_int||_inf = ||M* .* a||_inf, from the
    // PRE-damping solved accel U -> the genuine static unbalance (no f_v pollution).
    residualNorm = 0.0;
    for (int i = 0; i < size; i++) {
        double r = fabs((*Mstar)(i) * U(i));
        if (r > residualNorm) residualNorm = r;
    }

    // viscous-critical damping (dampMode 1): lag the damping force at the known
    // half-step velocity so it stays fully explicit (LHS untouched):
    //   a_work = a_static - cVisc * v_{n+1/2}   (f_v = -cVisc*M**Vhalf, /M* => -cVisc*Vhalf)
    Vector aWork(U);
    if (dampMode == 1)
        for (int i = 0; i < size; i++) aWork(i) -= cVisc * (*Vhalf)(i);

    // full-step velocity v_{n+1} = v_{n+1/2} + dt/2 a_work
    *Vfull = *Vhalf;
    Vfull->addVector(1.0, aWork, 0.5 * deltaT);

    // kinetic energy KE = 1/2 v^T M* v   (mass-weighted, matches EnergyBalance)
    kineticEnergy = 0.0;
    for (int i = 0; i < size; i++) {
        double vi = (*Vfull)(i);
        kineticEnergy += (*Mstar)(i) * vi * vi;
    }
    kineticEnergy *= 0.5;

    if (dampMode == 0) {
        // ---- Cundall kinetic damping: zero ALL velocities at each KE peak ----
        if (prevKE > 0.0 && kineticEnergy < prevKE) {
            Vhalf->Zero();
            Vfull->Zero();
            prevKE = 0.0;            // restart the KE sawtooth from rest
            // auto-refresh: rebuild M* from the CURRENT (stiffened) tangent at this
            // instant — velocities are zero so KE = 0 before & after => no energy
            // injected. Tracks far-branch stiffening so snap-through is stable with
            // NO manual -recompute (gershgorin only).
            if (autoRefresh && massMode == 0)
                this->buildFictitiousMass();
        } else {
            prevKE = kineticEnergy;
        }
    } else {
        prevKE = kineticEnergy;      // viscous: feed the watchdog, NEVER zero
    }

    // push the (possibly damped) full-step snapshot to the nodes
    theModel->setResponse(*Ut, *Vfull, aWork);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoDynamicRelaxation::update() - failed to update domain\n";
        return -6;
    }

    if (divergenceFactor > 0.0 && prevKE > 0.0) {
        if (kineticEnergy > divergenceFactor * prevKE) {
            opserr << "LadrunoDynamicRelaxation::update() - ABORT: KE grew "
                   << (kineticEnergy / prevKE) << "x in one step (> " << divergenceFactor << ")\n";
            return -7;
        }
    }
    if (verbose)
        opserr << "LadrunoDynamicRelaxation::update() max|a|=" << A_max
               << " RES=" << residualNorm << " KE=" << kineticEnergy << endln;

    // carry the (damped) accel so the next leap-frog half-step advance uses it
    *Aprev = aWork;
    return 0;
}

int
LadrunoDynamicRelaxation::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING LadrunoDynamicRelaxation::commit() - no AnalysisModel\n";
        return -1;
    }
    return theModel->commitDomain();
}

const Vector &
LadrunoDynamicRelaxation::getVel(void)
{
    return *Vhalf;
}

int
LadrunoDynamicRelaxation::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(10);
    data(0) = (double)massMode;
    data(1) = dtPseudo;
    data(2) = massScale;
    data(3) = (double)recomputeEvery;
    data(4) = interp ? 1.0 : 0.0;
    data(5) = divergenceFactor;
    data(6) = verbose ? 1.0 : 0.0;
    data(7) = (double)dampMode;
    data(8) = zetaTarget;
    data(9) = autoRefresh ? 1.0 : 0.0;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoDynamicRelaxation::sendSelf() - failed to send\n";
        return -1;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(10);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoDynamicRelaxation::recvSelf() - failed to receive\n";
        return -1;
    }
    massMode        = (int)data(0);
    dtPseudo        = data(1);
    massScale       = data(2);
    recomputeEvery  = (int)data(3);
    interp          = data(4) != 0.0;
    divergenceFactor = data(5);
    verbose         = data(6) != 0.0;
    dampMode        = (int)data(7);
    zetaTarget      = data(8);
    autoRefresh     = data(9) != 0.0;
    return 0;
}

void
LadrunoDynamicRelaxation::Print(OPS_Stream &s, int flag)
{
    s << "LadrunoDynamicRelaxation - quasi-static dynamic relaxation\n";
    s << "\t mass mode: " << (massMode == 0 ? "gershgorin" : massMode == 1 ? "lumped" : "unity")
      << "  dt: " << dtPseudo << "  recompute: " << recomputeEvery
      << "  autoRefresh: " << (autoRefresh ? "on" : "off") << endln;
    s << "\t damping: " << (dampMode == 0 ? "kinetic" : "viscous")
      << "  zeta: " << zetaTarget << "  cVisc: " << cVisc << endln;
    s << "\t last RES: " << residualNorm << "  KE: " << kineticEnergy << endln;
}
