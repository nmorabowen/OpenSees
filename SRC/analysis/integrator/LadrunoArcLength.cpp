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
#include <Channel.h>
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
        } else {
            opserr << "WARNING integrator LadrunoArcLength - unknown option " << opt
                   << " (ignored)\n";
        }
    }

    return new LadrunoArcLength(arcLength, alpha, adapt, Jd, ellMin, ellMax, pExp);
}

LadrunoArcLength::LadrunoArcLength(double arcLength, double alpha,
                                   bool adapt_, int Jd_,
                                   double ellMin_, double ellMax_, double pExp_)
 :StaticIntegrator(INTEGRATOR_TAGS_LadrunoArcLength),
  arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
  a(0.0), b(0.0), c(0.0), b24ac(0.0),
  deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), phat(0),
  deltaLambdaStep(0.0), currentLambda(0.0), signLastDeltaLambdaStep(1),
  adapt(adapt_), Jd(Jd_ > 0 ? Jd_ : 1), ellMin(ellMin_), ellMax(ellMax_),
  pExp(pExp_),
  // init Jlast = Jd so the FIRST step's factor (Jd/Jlast)^p == 1 (no jump),
  // mirroring LoadControl's numIncrLastStep = numIncr.
  numIncrLastStep(Jd_ > 0 ? Jd_ : 1)
{

}

LadrunoArcLength::~LadrunoArcLength()
{
    if (deltaUhat  != 0) delete deltaUhat;
    if (deltaU     != 0) delete deltaU;
    if (deltaUstep != 0) delete deltaUstep;
    if (deltaUbar  != 0) delete deltaUbar;
    if (phat       != 0) delete phat;
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

    return 0;
}

int
LadrunoArcLength::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(11);
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

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoArcLength::sendSelf() - failed to send the data\n";
        return -1;
    }
    return 0;
}

int
LadrunoArcLength::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(11);
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
        s << endln;
    } else
        s << "\t LadrunoArcLength - no associated AnalysisModel\n";
}
