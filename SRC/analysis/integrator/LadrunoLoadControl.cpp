/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno (ADR-80 S1) — N. Mora-Bowen
// LadrunoLoadControl — stock LoadControl + an optional linear displacement
// predictor. See LadrunoLoadControl.h for the mechanism this exists to fix.

#include <LadrunoLoadControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Domain.h>
#include <classTags.h>
#include <elementAPI.h>
#include <math.h>

void* OPS_LadrunoLoadControl()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING LadrunoLoadControl - insufficient arguments\n";
        opserr << "Want: integrator LadrunoLoadControl lambda? "
               << "<numIter? minLambda? maxLambda?> <-extrapolate frac?>\n";
        return 0;
    }

    double lambda;
    int numData = 1;
    if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
        opserr << "WARNING LadrunoLoadControl - failed to read double lambda\n";
        return 0;
    }

    int numIter = 1;
    double mLambda[2] = {lambda, lambda};
    double frac = 0.0;

    // The stock parser reads the optional numIter/min/max triple only when
    // MORE THAN TWO args remain. Keep that rule, but count only the POSITIONAL
    // args so that `-extrapolate f` appended after `lambda` cannot be mistaken
    // for the triple.
    int nRemaining = OPS_GetNumRemainingInputArgs();
    bool hasTriple = false;
    if (nRemaining > 2) {
        // peek: the triple is present only if the next token parses as a number
        const char* nxt = OPS_GetString();
        OPS_ResetCurrentInputArg(-1);
        if (nxt != 0 && nxt[0] != '-')
            hasTriple = true;
    }

    if (hasTriple) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &numIter) < 0) {
            opserr << "WARNING LadrunoLoadControl - failed to read int numIter\n";
            return 0;
        }
        numData = 2;
        if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
            opserr << "WARNING LadrunoLoadControl - failed to read min and max\n";
            return 0;
        }
    }

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-extrapolate") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING LadrunoLoadControl -extrapolate needs a value; "
                       << "0.0 (stock behaviour) assumed\n";
                break;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &frac) < 0) {
                opserr << "WARNING LadrunoLoadControl - failed to read "
                       << "-extrapolate frac; 0.0 assumed\n";
                frac = 0.0;
            }
        } else {
            // Degrade, do not return null (ADR-75 P1d rule): an unknown token
            // must be loud, not silently dropped the way ModifiedNewton's
            // parser used to drop options (ADR-76 R4).
            opserr << "WARNING LadrunoLoadControl - unknown option '" << opt
                   << "' ignored\n";
        }
    }

    if (frac < 0.0) {
        opserr << "WARNING LadrunoLoadControl - negative -extrapolate ("
               << frac << ") clamped to 0.0\n";
        frac = 0.0;
    }

    return new LadrunoLoadControl(lambda, numIter, mLambda[0], mLambda[1], frac);
}

LadrunoLoadControl::LadrunoLoadControl(double dLambda, int numIncr,
                                       double min, double max, double frac)
 : StaticIntegrator(INTEGRATOR_TAGS_LadrunoLoadControl),
   deltaLambda(dLambda),
   specNumIncrStep(numIncr), numIncrLastStep(numIncr),
   dLambdaMin(min), dLambdaMax(max),
   extrapolate(frac),
   dUprev(0), dUaccum(0), dLambdaPrev(dLambda),
   havePrev(false), stepOpen(false), justCutBack(false)
{
  // stock: avoid divide-by-zero on the first update()
  if (numIncr == 0) {
    opserr << "WARNING LadrunoLoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}

LadrunoLoadControl::~LadrunoLoadControl()
{
  if (dUprev  != 0) delete dUprev;
  if (dUaccum != 0) delete dUaccum;
}

// Ladruno (ADR-80 S1): the predictor's stored increment is expressed in the
// CURRENT equation numbering. Any domainChanged() renumbers, so the stored
// vector becomes meaningless and must not be reused.
void
LadrunoLoadControl::invalidatePredictor(void)
{
  if (dUprev  != 0) { delete dUprev;  dUprev  = 0; }
  if (dUaccum != 0) { delete dUaccum; dUaccum = 0; }
  havePrev    = false;
  stepOpen    = false;
  justCutBack = false;
}

int
LadrunoLoadControl::domainChanged(void)
{
  this->invalidatePredictor();
  return this->StaticIntegrator::domainChanged();
}

int
LadrunoLoadControl::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "LadrunoLoadControl::newStep() - no associated AnalysisModel\n";
        return -1;
    }

    // ---- Ladruno (ADR-80 S1): cutback detection ---------------------------
    // newStep() reaching here with a step still OPEN means the previous attempt
    // never reached commit(), i.e. the caller's driver cut back. The partial dU
    // of that failed attempt is garbage and must not seed a prediction; Abaqus
    // likewise abandons extrapolation when the increment changes sharply.
    if (stepOpen)
        justCutBack = true;
    stepOpen = true;
    if (dUaccum != 0)
        dUaccum->Zero();

    // ---- stock LoadControl: adapt dlambda by last step's iteration count ---
    double factor = specNumIncrStep/numIncrLastStep;
    deltaLambda *= factor;

    if (deltaLambda < dLambdaMin)
      deltaLambda = dLambdaMin;
    else if (deltaLambda > dLambdaMax)
      deltaLambda = dLambdaMax;

    double currentLambda = theModel->getCurrentDomainTime();
    currentLambda += deltaLambda;

    // ---- Ladruno (ADR-80 S1): the predictor -------------------------------
    // Ordering copied from DisplacementControl::newStep (:287-289):
    //   incrDisp -> applyLoadDomain -> updateDomain
    // so that when enforceSPs pre-updates the driven element layer, that
    // layer's NEIGHBOURS have already moved and its strain increment is the
    // physical one instead of delta/h.
    bool doPredict = (extrapolate > 0.0) && havePrev && !justCutBack &&
                     (dUprev != 0) && (dLambdaPrev != 0.0);

    if (doPredict) {
        double ratio = deltaLambda/dLambdaPrev;
        // Abaqus abandons extrapolation when the increment changes sharply.
        // Ours is driven by a caller's adaptive ladder, which resizes often, so
        // this guard fires more than Abaqus's would — deliberately conservative:
        // skipping the prediction only costs the stock path.
        if (ratio < 0.5 || ratio > 2.0)
            doPredict = false;
        else {
            Vector predicted(*dUprev);
            predicted *= (extrapolate * ratio);
            theModel->incrDisp(predicted);
        }
    }

    theModel->applyLoadDomain(currentLambda);

    if (doPredict) {
        if (theModel->updateDomain() < 0) {
            opserr << "LadrunoLoadControl::newStep() - model failed to update "
                   << "after extrapolation\n";
            return -1;
        }
    }

    justCutBack = false;
    numIncrLastStep = 0;

    return 0;
}

int
LadrunoLoadControl::update(const Vector &deltaU)
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    if (myModel == 0 || theSOE == 0) {
        opserr << "WARNING LadrunoLoadControl::update() - No AnalysisModel or "
               << "LinearSOE has been set\n";
        return -1;
    }

    myModel->incrDisp(deltaU);
    if (myModel->updateDomain() < 0) {
      opserr << "LadrunoLoadControl::update - model failed to update for new dU\n";
      return -1;
    }

    // Ladruno (ADR-80 S1): accumulate the step's total corrector so commit()
    // can promote it. Costs one vector add per iteration and only when the
    // predictor is enabled.
    if (extrapolate > 0.0) {
        if (dUaccum == 0)
            dUaccum = new Vector(deltaU);
        else if (dUaccum->Size() != deltaU.Size())
            { delete dUaccum; dUaccum = new Vector(deltaU); }
        else
            *dUaccum += deltaU;
    }

    // Set deltaU for the convergence test
    theSOE->setX(deltaU);

    numIncrLastStep++;

    return 0;
}

int
LadrunoLoadControl::commit(void)
{
    // Ladruno (ADR-80 S1): a step that commits is a step whose accumulated dU
    // is a legitimate basis for predicting the next one. Promote, then close
    // the step so the next newStep() does not read it as a cutback.
    if (extrapolate > 0.0 && dUaccum != 0) {
        if (dUprev == 0)
            dUprev = new Vector(*dUaccum);
        else if (dUprev->Size() != dUaccum->Size())
            { delete dUprev; dUprev = new Vector(*dUaccum); }
        else
            *dUprev = *dUaccum;
        dLambdaPrev = deltaLambda;
        havePrev = true;
    }
    stepOpen    = false;
    justCutBack = false;

    return this->StaticIntegrator::commit();
}

int
LadrunoLoadControl::setDeltaLambda(double newValue)
{
  // we set the #incr at last step = #incr so get newValue incr
  numIncrLastStep = specNumIncrStep;
  deltaLambda = newValue;
  return 0;
}

int
LadrunoLoadControl::sendSelf(int cTag, Channel &theChannel)
{
  // Ladruno (ADR-80 S1): parameters only. The predictor's stored increment is
  // tied to this process's equation numbering and is ADVISORY — losing it
  // across a channel is safe by construction (the receiver simply starts
  // without a prediction), so it is deliberately not transmitted.
  Vector data(6);
  data(0) = deltaLambda;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = dLambdaMin;
  data(4) = dLambdaMax;
  data(5) = extrapolate;
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LadrunoLoadControl::sendSelf() - failed to send the Vector\n";
      return -1;
  }
  return 0;
}

int
LadrunoLoadControl::recvSelf(int cTag, Channel &theChannel,
                             FEM_ObjectBroker &theBroker)
{
  Vector data(6);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LadrunoLoadControl::recvSelf() - failed to receive the Vector\n";
      deltaLambda = 0;
      return -1;
  }
  deltaLambda     = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  dLambdaMin      = data(3);
  dLambdaMax      = data(4);
  extrapolate     = data(5);

  // Ladruno (ADR-80 S1): whatever predictor state this object held belonged to
  // the sender's numbering. Start clean.
  this->invalidatePredictor();

  return 0;
}

void
LadrunoLoadControl::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentLambda = theModel->getCurrentDomainTime();
        s << "\t LadrunoLoadControl - currentLambda: " << currentLambda;
        s << "  deltaLambda: " << deltaLambda;
        s << "  extrapolate: " << extrapolate << endln;
    } else
        s << "\t LadrunoLoadControl - no associated AnalysisModel\n";
}
