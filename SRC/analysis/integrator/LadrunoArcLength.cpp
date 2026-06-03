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

// Ladruno: LadrunoArcLength — fork-authored static integrator (classTag 33004).
// Core algorithm copied from ArcLength (fmk, 07/98); the DDM sensitivity
// machinery is intentionally omitted. The only new behaviour is the optional
// Ramm desired-iteration radius adaptation, gated by -adapt (default OFF =>
// bit-identical to ArcLength). See
// Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.

#include <LadrunoArcLength.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Element.h>
#include <math.h>
#include <string.h>
#include <elementAPI.h>
#include <classTags.h>

void* OPS_LadrunoArcLength()
{
    // Usage: integrator LadrunoArcLength arcLength alpha
    //                   <-adapt Jd ellMin ellMax> <-p exp>
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING integrator LadrunoArcLength arcLength alpha "
                  "<-adapt Jd ellMin ellMax> <-p exp>\n";
        return 0;
    }

    int numdata = 1;
    double arcLength = 0.0, alpha = 1.0;
    if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
        opserr << "WARNING integrator LadrunoArcLength failed to read arcLength\n";
        return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
        opserr << "WARNING integrator LadrunoArcLength failed to read alpha\n";
        return 0;
    }

    bool   adapt = false;
    int    Jd = 5;
    double ellMin = 0.0, ellMax = 0.0, pExp = 1.0;
    bool   stabilize = false, adaptStab = false;
    double fTarget = 2.0e-4, cVisc = 0.0, dtPseudo = 1.0, massScale = 1.0;
    int    massMode = 0;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-adapt") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING integrator LadrunoArcLength -adapt expects "
                          "Jd ellMin ellMax\n";
                return 0;
            }
            int idata = 1;
            if (OPS_GetIntInput(&idata, &Jd) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -adapt failed to read Jd\n";
                return 0;
            }
            int ndata = 2;
            double mm[2] = {0.0, 0.0};
            if (OPS_GetDoubleInput(&ndata, mm) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -adapt failed to read "
                          "ellMin ellMax\n";
                return 0;
            }
            ellMin = mm[0];
            ellMax = mm[1];
            adapt = true;
        } else if (strcmp(opt, "-p") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING integrator LadrunoArcLength -p expects exp\n";
                return 0;
            }
            int idata = 1;
            if (OPS_GetDoubleInput(&idata, &pExp) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -p failed to read exp\n";
                return 0;
            }
        } else if (strcmp(opt, "-stabilize") == 0) {
            stabilize = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {   // optional trailing f fraction
                int nd = 1; double f;
                if (OPS_GetDoubleInput(&nd, &f) == 0) fTarget = f;
            }
        } else if (strcmp(opt, "-adaptStab") == 0) {
            adaptStab = true;
        } else if (strcmp(opt, "-cVisc") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &cVisc) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -cVisc failed to read value\n";
                return 0;
            }
        } else if (strcmp(opt, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING integrator LadrunoArcLength -mass expects gershgorin|lumped|unity\n";
                return 0;
            }
            const char *m = OPS_GetString();
            if (strcmp(m, "gershgorin") == 0)      massMode = 0;
            else if (strcmp(m, "lumped") == 0) {
                massMode = 1;
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    int nd = 1; if (OPS_GetDoubleInput(&nd, &massScale) < 0) massScale = 1.0;
                }
            } else if (strcmp(m, "unity") == 0)    massMode = 2;
            else opserr << "WARNING integrator LadrunoArcLength - unknown -mass " << m
                        << " (keeping gershgorin)\n";
        } else {
            opserr << "WARNING integrator LadrunoArcLength - unknown option " << opt
                   << " (ignored)\n";
        }
    }

    return new LadrunoArcLength(arcLength, alpha, adapt, Jd, ellMin, ellMax, pExp,
                                stabilize, fTarget, adaptStab, massMode, massScale,
                                cVisc, dtPseudo);
}

LadrunoArcLength::LadrunoArcLength(double arcLength, double alpha,
                                   bool adapt_, int Jd_,
                                   double ellMin_, double ellMax_, double pExp_,
                                   bool stabilize_, double fTarget_, bool adaptStab_,
                                   int massMode_, double massScale_, double cVisc_,
                                   double dtPseudo_)
 :StaticIntegrator(INTEGRATOR_TAGS_LadrunoArcLength),
  arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
  a(0.0), b(0.0), c(0.0), b24ac(0.0),
  deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), phat(0),
  deltaLambdaStep(0.0), currentLambda(0.0), signLastDeltaLambdaStep(1),
  adapt(adapt_), Jd(Jd_ > 0 ? Jd_ : 1), ellMin(ellMin_), ellMax(ellMax_),
  pExp(pExp_),
  // init Jlast = Jd so the FIRST step's factor (Jd/Jlast)^p == 1 (no jump),
  // mirroring LoadControl's numIncrLastStep = numIncr.
  numIncrLastStep(Jd_ > 0 ? Jd_ : 1),
  stabilize(stabilize_), adaptStab(adaptStab_), fTarget(fTarget_), cVisc(cVisc_),
  dtPseudo(dtPseudo_ > 0.0 ? dtPseudo_ : 1.0),
  cOverDt(cVisc_ / (dtPseudo_ > 0.0 ? dtPseudo_ : 1.0)),
  // if cVisc supplied (-cVisc), latch calibration off (deterministic); else
  // calibrate once from fTarget at the first stabilized commit.
  cCalibrated(cVisc_ != 0.0), Estrain0(0.0), dissipVisc(0.0), residualTrueNorm(0.0),
  dLambda0(sqrt(arcLength*arcLength)), massMode(massMode_), massScale(massScale_),
  nEqn(0), Mstar(0),
  // Layer-B snapshot: seed with the construction radius so a failure before the
  // first commit() reverts to a sane state (haveSnapshot latches at commit()).
  cArcLength2(arcLength*arcLength), cDeltaLambdaStep(0.0), cCurrentLambda(0.0),
  cSign(1), cDeltaUstep(0), haveSnapshot(false)
{

}

LadrunoArcLength::~LadrunoArcLength()
{
    if (deltaUhat  != 0) delete deltaUhat;
    if (deltaU     != 0) delete deltaU;
    if (deltaUstep != 0) delete deltaUstep;
    if (deltaUbar  != 0) delete deltaUbar;
    if (phat       != 0) delete phat;
    if (Mstar      != 0) delete Mstar;
    if (cDeltaUstep != 0) delete cDeltaUstep;
}

// ---- integrator-owned artificial diagonal mass (lifted from
// LadrunoDynamicRelaxation::buildFictitiousMass; NO dt^2/4 prefactor -- the scale
// lives in cOverDt). Gershgorin row-sum of K, or lumped real mass, or unity.
// NOT Element::getMass() (zero on zero-density models -- the ADR-20 BLOCKER).
int
LadrunoArcLength::buildArtificialMass(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Mstar == 0) return -1;
    Mstar->Zero();

    if (massMode == 2) {                 // unity
        for (int i = 0; i < nEqn; i++) (*Mstar)(i) = 1.0;
        return 0;
    }
    FE_EleIter &theEles = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEles()) != 0) {
        Element *e = elePtr->getElement();
        if (e == 0) continue;
        const ID &id = elePtr->getID();
        const Matrix &Ke = (massMode == 1) ? e->getMass() : e->getTangentStiff();
        int n = id.Size();
        if (Ke.noRows() != n) continue;
        for (int i = 0; i < n; i++) {
            int loc = id(i);
            if (loc < 0) continue;
            double rs = 0.0;
            for (int j = 0; j < n; j++) rs += fabs(Ke(i, j));
            (*Mstar)(loc) += rs;
        }
    }
    if (massMode == 1)
        for (int i = 0; i < nEqn; i++) (*Mstar)(i) *= massScale;

    double mmax = 0.0;
    for (int i = 0; i < nEqn; i++) if ((*Mstar)(i) > mmax) mmax = (*Mstar)(i);
    double fl = (mmax > 0.0) ? mmax * 1.0e-8 : 1.0;
    for (int i = 0; i < nEqn; i++) if (!((*Mstar)(i) > 0.0)) (*Mstar)(i) = fl;
    return 0;
}

// ---- regularize the tangent: K + cOverDt*M* on the diagonal (additive; the DR
// diagonal-poke idiom). Identity when !stabilize.
int
LadrunoArcLength::formTangent(int statFlag)
{
    this->IncrementalIntegrator::formTangent(statFlag);   // assemble real K onto the SOE
    if (stabilize && Mstar != 0) {
        LinearSOE *soe = this->getLinearSOE();
        static Matrix m1(1, 1);
        static ID id1(1);
        for (int i = 0; i < nEqn; i++) {
            m1(0, 0) = cOverDt * (*Mstar)(i);
            id1(0) = i;
            soe->addA(m1, id1);
        }
    }
    return 0;
}

// ---- inject the artificial viscous force: B <- (lambda p - f_int) - cOverDt*M**deltaUstep.
// Records the TRUE static unbalance (without f_v) for the watchdog. Identity when !stabilize.
int
LadrunoArcLength::formUnbalance(void)
{
    this->IncrementalIntegrator::formUnbalance();         // B = lambda p - f_int
    if (stabilize && Mstar != 0 && deltaUstep != 0) {
        LinearSOE *soe = this->getLinearSOE();
        const Vector &B = soe->getB();
        residualTrueNorm = B.Norm();                      // true unbalance BEFORE f_v
        Vector R(B);
        for (int i = 0; i < nEqn; i++)
            R(i) -= cOverDt * (*Mstar)(i) * (*deltaUstep)(i);
        soe->setB(R);
    }
    return 0;
}

int
LadrunoArcLength::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::newStep() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

    // base sign of the load change on what happened last step
    if (deltaLambdaStep < 0)
        signLastDeltaLambdaStep = -1;
    else
        signLastDeltaLambdaStep = +1;

    // ======== -stabilize: viscous-regularized incremental LOAD CONTROL ========
    // Mutually exclusive with the arc-length quadratic (ADR-20 decision 3b): no
    // predictor phat-solve, lambda is NOT a Newton unknown; the regularized K
    // (formTangent) lets Newton clear the limit point. -adapt drives the LOAD
    // increment dLambda0 here (not arcLength2).
    if (stabilize) {
        if (adapt) {
            int jlast = numIncrLastStep > 0 ? numIncrLastStep : 1;
            dLambda0 *= pow((double)Jd / (double)jlast, pExp);
            if (dLambda0 < ellMin) dLambda0 = ellMin;
            else if (ellMax > 0.0 && dLambda0 > ellMax) dLambda0 = ellMax;
        }
        numIncrLastStep = 0;
        deltaUstep->Zero();            // per-step increment accumulator for f_v
        deltaLambdaStep = signLastDeltaLambdaStep * dLambda0;
        currentLambda += deltaLambdaStep;
        theModel->applyLoadDomain(currentLambda);
        if (theModel->updateDomain() < 0) {
            opserr << "LadrunoArcLength::newStep() - model failed to update (stabilize)\n";
            return -1;
        }
        return 0;
    }

    // -------- Layer A: Ramm desired-iteration radius adaptation --------
    // arcLength <- clamp( arcLength * (Jd / Jlast)^p , ellMin, ellMax )
    // numIncrLastStep holds the corrector-iteration count of the LAST step.
    if (adapt) {
        int jlast = numIncrLastStep > 0 ? numIncrLastStep : 1;  // guard div-by-zero
        double factor = pow((double)Jd / (double)jlast, pExp);
        double ell = sqrt(arcLength2) * factor;
        if (ell < ellMin) ell = ellMin;
        else if (ell > ellMax) ell = ellMax;
        arcLength2 = ell * ell;
    }
    numIncrLastStep = 0;   // reset for this step's correctors (counted in update)

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
        opserr << "LadrunoArcLength::newStep(void) - failed in solver\n";
        return -1;
    }
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    // determine delta lambda(1) == dlambda
    double dLambda = sqrt(arcLength2 / ((dUhat ^ dUhat) + alpha2));
    dLambda *= signLastDeltaLambdaStep;
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoArcLength::newStep() - model failed to update for new dU\n";
        return -1;
    }

    return 0;
}

int
LadrunoArcLength::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::update() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // -stabilize: pure load-control corrector. lambda is fixed within the step; the
    // regularized K (formTangent) + viscous residual (formUnbalance) make Newton
    // converge. Accumulate the per-step increment deltaUstep that drives f_v.
    if (stabilize) {
        (*deltaU) = dU;
        (*deltaUstep) += dU;
        theModel->incrDisp(dU);
        theModel->applyLoadDomain(currentLambda);
        if (theModel->updateDomain() < 0) {
            opserr << "LadrunoArcLength::update() - model failed to update (stabilize)\n";
            return -1;
        }
        theLinSOE->setX(dU);
        numIncrLastStep++;
        return 0;
    }

    (*deltaUbar) = dU;   // have to do this as the SOE is gonna change

    // determine dUhat
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();

    // determine the coefficients of our quadratic equation
    a = alpha2 + ((*deltaUhat) ^ (*deltaUhat));
    b = alpha2 * deltaLambdaStep
        + ((*deltaUhat) ^ (*deltaUbar))
        + ((*deltaUstep) ^ (*deltaUhat));
    b *= 2.0;
    c = 2.0 * ((*deltaUstep) ^ (*deltaUbar)) + ((*deltaUbar) ^ (*deltaUbar));

    // check for a solution to quadratic
    b24ac = b * b - 4.0 * a * c;
    if (b24ac < 0) {
        opserr << "LadrunoArcLength::update() - imaginary roots due to multiple "
                  "instability directions - initial load increment was too large\n";
        opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << endln;
        return -1;
    }
    double a2 = 2.0 * a;
    if (a2 == 0.0) {
        opserr << "LadrunoArcLength::update() - zero denominator,"
                  " alpha was set to 0.0 and zero reference load\n";
        return -2;
    }

    // determine the roots of the quadratic
    double sqrtb24ac = sqrt(b24ac);
    double dlambda1 = (-b + sqrtb24ac) / a2;
    double dlambda2 = (-b - sqrtb24ac) / a2;

    double val = (*deltaUhat) ^ (*deltaUstep);
    double theta1 = ((*deltaUstep) ^ (*deltaUstep)) + ((*deltaUbar) ^ (*deltaUstep));
    theta1 += dlambda1 * val;

    // choose dLambda based on angle between incremental displacement before
    // and after this step -- want positive
    double dLambda;
    if (theta1 > 0)
        dLambda = dlambda1;
    else
        dLambda = dlambda2;

    // determine delta U(i)
    (*deltaU) = (*deltaUbar);
    deltaU->addVector(1.0, *deltaUhat, dLambda);

    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoArcLength::update() - model failed to update for new dU\n";
        return -1;
    }

    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    // Layer A: count this corrector iteration
    numIncrLastStep++;

    return 0;
}

int
LadrunoArcLength::domainChanged(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::domainChanged() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }
    int size = theModel->getNumEqn();   // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) {
        if (deltaUhat != 0) delete deltaUhat;
        deltaUhat = new Vector(size);
    }
    if (deltaUbar == 0 || deltaUbar->Size() != size) {
        if (deltaUbar != 0) delete deltaUbar;
        deltaUbar = new Vector(size);
    }
    if (deltaU == 0 || deltaU->Size() != size) {
        if (deltaU != 0) delete deltaU;
        deltaU = new Vector(size);
    }
    if (deltaUstep == 0 || deltaUstep->Size() != size) {
        if (deltaUstep != 0) delete deltaUstep;
        deltaUstep = new Vector(size);
    }
    if (phat == 0 || phat->Size() != size) {
        if (phat != 0) delete phat;
        phat = new Vector(size);
    }
    // Layer-B committed snapshot of deltaUstep (sized with the others)
    if (cDeltaUstep == 0 || cDeltaUstep->Size() != size) {
        if (cDeltaUstep != 0) delete cDeltaUstep;
        cDeltaUstep = new Vector(size);
    }
    if (deltaUhat == 0 || deltaUbar == 0 || deltaU == 0 ||
        deltaUstep == 0 || phat == 0 ||
        deltaUhat->Size() != size || phat->Size() != size) {
        opserr << "FATAL LadrunoArcLength::domainChanged() - ran out of memory "
                  "for Vectors of size " << size << endln;
        return -1;
    }

    // now we have to determine phat:
    // do this by incrementing lambda by 1, applying load and getting phat
    // from the unbalance.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance();   // NOTE: assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);

    // check there is a reference load
    int haveLoad = 0;
    for (int i = 0; i < size; i++)
        if ((*phat)(i) != 0.0) {
            haveLoad = 1;
            i = size;
        }
    if (haveLoad == 0) {
        opserr << "WARNING LadrunoArcLength::domainChanged() - zero reference load\n";
        return -1;
    }

    // -stabilize: build the artificial M* AFTER the phat probe (the formUnbalance
    // override above is a no-op while Mstar==0, so phat stays uncorrupted).
    nEqn = size;
    if (stabilize) {
        if (Mstar == 0 || Mstar->Size() != size) {
            if (Mstar != 0) delete Mstar;
            Mstar = new Vector(size);
        }
        this->buildArtificialMass();
    }

    // Seed the Layer-B snapshot with the (committed) starting state so a failure
    // before the first commit() reverts somewhere sane. haveSnapshot stays false
    // until a real commit() captures a converged step.
    cArcLength2      = arcLength2;
    cDeltaLambdaStep = deltaLambdaStep;
    cCurrentLambda   = currentLambda;
    cSign            = signLastDeltaLambdaStep;
    cDeltaUstep->Zero();

    return 0;
}

int
LadrunoArcLength::commit(void)
{
    if (stabilize && Mstar != 0 && deltaUstep != 0) {
        // accumulate viscous work W_v = sum cOverDt*M*_i*deltaUstep_i^2 (watchdog)
        double Wv = 0.0;
        for (int i = 0; i < nEqn; i++)
            Wv += cOverDt * (*Mstar)(i) * (*deltaUstep)(i) * (*deltaUstep)(i);
        dissipVisc += Wv;

        // one-shot calibration from the energy fraction f (skipped if -cVisc given)
        if (!cCalibrated) {
            Estrain0 = 0.5 * fabs((*deltaUstep) ^ (*phat)) * fabs(currentLambda);
            double Wunit = 0.0;
            for (int i = 0; i < nEqn; i++)
                Wunit += (*Mstar)(i) * (*deltaUstep)(i) * (*deltaUstep)(i) / dtPseudo;
            cVisc   = (Wunit > 0.0) ? fTarget * Estrain0 / Wunit : 0.0;
            cOverDt = cVisc / dtPseudo;
            cCalibrated = true;
        } else if (adaptStab && Estrain0 > 0.0) {
            double ratio = dissipVisc / Estrain0;
            if (ratio > 0.0) { cVisc *= (fTarget / ratio); cOverDt = cVisc / dtPseudo; }
        }
    }

    // -------- Layer B: snapshot the converged per-step state (ADR-20 §4.3) --------
    // This is the state stock ArcLength pollutes in newStep()/update() BEFORE its
    // b24ac<0 bail; revertToLastStep() restores exactly this on a failed retry.
    cArcLength2      = arcLength2;
    cDeltaLambdaStep = deltaLambdaStep;
    cCurrentLambda   = currentLambda;
    cSign            = signLastDeltaLambdaStep;
    if (cDeltaUstep != 0 && deltaUstep != 0 &&
        cDeltaUstep->Size() == deltaUstep->Size())
        *cDeltaUstep = *deltaUstep;
    haveSnapshot = true;

    return this->StaticIntegrator::commit();
}

int
LadrunoArcLength::revertToLastStep(void)
{
    // Restore the committed per-step state captured at the last commit(). Stock
    // ArcLength inherits a no-op revertToLastStep(); ours undoes the in-flight
    // mutation of newStep()/update() so a script's reduceStep()+retry starts from
    // the genuine last-converged state, not the polluted one.
    //
    // We restore ONLY the four per-step items of ADR-20 decision 4
    // {deltaLambdaStep, deltaUstep, currentLambda, signLastDeltaLambdaStep} —
    // deliberately NOT arcLength2 / dLambda0: the radius is SCRIPT-owned policy,
    // and restoring it would wipe each reduceStep() so the documented
    //   while ok!=0: reduceStep(0.5); analyze(1)
    // cut-and-retry loop could never make progress (it would re-test the same
    // 0.5x radius forever). -adapt re-derives the radius from numIncrLastStep in
    // the next newStep() regardless, so leaving arcLength2 alone is safe there too.
    //
    // No-op until the first commit() so identity (no failure observed before any
    // converged step) is unperturbed.
    if (!haveSnapshot)
        return 0;

    deltaLambdaStep         = cDeltaLambdaStep;
    currentLambda           = cCurrentLambda;
    signLastDeltaLambdaStep = cSign;
    if (deltaUstep != 0 && cDeltaUstep != 0 &&
        deltaUstep->Size() == cDeltaUstep->Size())
        *deltaUstep = *cDeltaUstep;
    return 0;
}

// ---- Layer-B scalar mutators (ADR-20 §4.3). Pure scalar mutation of the arc
// radius, no reallocation; a script drives the cut-and-retry loop:
//   ok = analyze(1); while ok!=0: ladrunoArcLength reduceStep 0.5; ok = analyze(1)
// StaticAnalysis::analyze already calls revertToLastStep() on a failed step, so
// these operate on the restored committed radius.
int
LadrunoArcLength::reduceStep(double f)
{
    if (f <= 0.0) {
        opserr << "WARNING LadrunoArcLength::reduceStep - factor must be > 0 "
                  "(got " << f << "); ignored\n";
        return -1;
    }
    double newAL2 = arcLength2 * f * f;
    // ell >= floor guard: use ellMin if the user set one, else a tiny absolute
    // floor so the radius cannot collapse to zero.
    double floor2 = (ellMin > 0.0) ? ellMin * ellMin : 1.0e-30;
    if (newAL2 < floor2) newAL2 = floor2;
    arcLength2 = newAL2;
    dLambda0  *= f;                       // keep the stabilize load-control predictor in sync
    return 0;
}

int
LadrunoArcLength::increaseStep(double f)
{
    if (f <= 0.0) {
        opserr << "WARNING LadrunoArcLength::increaseStep - factor must be > 0 "
                  "(got " << f << "); ignored\n";
        return -1;
    }
    double newAL2 = arcLength2 * f * f;
    if (ellMax > 0.0) {                   // ell <= ellMax ceiling (only if set)
        double ceil2 = ellMax * ellMax;
        if (newAL2 > ceil2) newAL2 = ceil2;
    }
    arcLength2 = newAL2;
    dLambda0  *= f;
    return 0;
}

int
LadrunoArcLength::setArcLength(double arcLength)
{
    if (arcLength <= 0.0) {
        opserr << "WARNING LadrunoArcLength::setArcLength - arcLength must be > 0 "
                  "(got " << arcLength << "); ignored\n";
        return -1;
    }
    arcLength2 = arcLength * arcLength;
    dLambda0   = arcLength;
    return 0;
}

double LadrunoArcLength::getArcLength(void) const           { return sqrt(arcLength2); }
double LadrunoArcLength::getDeltaLambdaStep(void) const     { return deltaLambdaStep; }
double LadrunoArcLength::getCurrentLambda(void) const       { return currentLambda; }
int    LadrunoArcLength::getSignLastDeltaLambdaStep(void) const { return signLastDeltaLambdaStep; }
double LadrunoArcLength::getDeltaUstepNorm(void) const
{
    return (deltaUstep != 0) ? deltaUstep->Norm() : 0.0;
}

int
LadrunoArcLength::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(27);
    data(0) = arcLength2;
    data(1) = alpha2;
    data(2) = deltaLambdaStep;
    data(3) = currentLambda;
    data(4) = signLastDeltaLambdaStep;
    data(5) = adapt ? 1.0 : 0.0;
    data(6) = (double)Jd;
    data(7) = ellMin;
    data(8) = ellMax;
    data(9) = pExp;
    data(10) = (double)numIncrLastStep;
    data(11) = stabilize ? 1.0 : 0.0;
    data(12) = adaptStab ? 1.0 : 0.0;
    data(13) = fTarget;
    data(14) = cVisc;
    data(15) = dtPseudo;
    data(16) = cOverDt;
    data(17) = cCalibrated ? 1.0 : 0.0;
    data(18) = (double)massMode;
    data(19) = massScale;
    // Layer-B committed snapshot scalars (cDeltaUstep is rebuilt by domainChanged,
    // matching stock ArcLength which never serializes its per-step vectors).
    data(20) = dLambda0;
    data(21) = cArcLength2;
    data(22) = cDeltaLambdaStep;
    data(23) = cCurrentLambda;
    data(24) = (double)cSign;
    data(25) = haveSnapshot ? 1.0 : 0.0;
    data(26) = 0.0;   // reserved

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoArcLength::sendSelf() - failed to send the data\n";
        return -1;
    }
    return 0;
}

int
LadrunoArcLength::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(27);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoArcLength::recvSelf() - failed to receive the data\n";
        return -1;
    }
    arcLength2              = data(0);
    alpha2                  = data(1);
    deltaLambdaStep         = data(2);
    currentLambda           = data(3);
    signLastDeltaLambdaStep = (int)data(4);
    adapt                   = data(5) != 0.0;
    Jd                      = (int)data(6);
    ellMin                  = data(7);
    ellMax                  = data(8);
    pExp                    = data(9);
    numIncrLastStep         = (int)data(10);
    stabilize               = data(11) != 0.0;
    adaptStab               = data(12) != 0.0;
    fTarget                 = data(13);
    cVisc                   = data(14);
    dtPseudo                = data(15);
    cOverDt                 = data(16);
    cCalibrated             = data(17) != 0.0;
    massMode                = (int)data(18);
    massScale               = data(19);
    dLambda0                = data(20);
    cArcLength2             = data(21);
    cDeltaLambdaStep        = data(22);
    cCurrentLambda          = data(23);
    cSign                   = (int)data(24);
    haveSnapshot            = data(25) != 0.0;
    return 0;
}

void
LadrunoArcLength::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double cLambda = theModel->getCurrentDomainTime();
        s << "\t LadrunoArcLength - currentLambda: " << cLambda;
        s << "  arcLength: " << sqrt(arcLength2) << "  alpha: " << sqrt(alpha2);
        if (adapt)
            s << "  [adapt Jd=" << Jd << " ell=[" << ellMin << "," << ellMax
              << "] p=" << pExp << "]";
        if (stabilize)
            s << "  [stabilize f=" << fTarget << " c=" << cVisc
              << " diss/E=" << (Estrain0 > 0.0 ? dissipVisc / Estrain0 : 0.0)
              << " trueRes=" << residualTrueNorm << "]";
        s << endln;
    } else
        s << "\t LadrunoArcLength - no associated AnalysisModel\n";
}
