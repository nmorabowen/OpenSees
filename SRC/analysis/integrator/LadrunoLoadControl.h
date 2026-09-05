/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno (ADR-80 S1) — N. Mora-Bowen
//
// LadrunoLoadControl: a strict superset of stock LoadControl adding an optional
// linear DISPLACEMENT PREDICTOR (`-extrapolate <frac>`), the Abaqus
// EXTRAPOLATION=LINEAR analogue.
//
// WHY. Under `constraints Transformation` a non-homogeneous `sp` has its DOF
// genuinely eliminated (setID(dof,-1)), so the only way the prescribed motion
// reaches the RHS is TransformationConstraintHandler::enforceSPs pre-updating
// the elements attached to the driven node (:518-521). Stock LoadControl has
// NOTHING in front of applyLoadDomain (LoadControl.cpp:130), so the first
// constitutive evaluation of every increment happens with the driven face
// advanced by the full increment and the interior still at the last commit —
// an overstrain of order L/h in one element layer. Harmless for an elastic law;
// for a path-dependent one the layer yields SPURIOUSLY, its consistent tangent
// collapses toward 2G*H/(3G+H), and the march degenerates into cutbacks.
// Measured x18.7 iterations and 23 cutbacks synthetically (ADR-80 G5), x28.6
// and 52 cutbacks on the Cerro Lindo fuse — for an identical converged answer.
//
// The remedy is a static predictor, not a new constraint handler. The correct
// ordering already exists in the tree: DisplacementControl::newStep does
// incrDisp -> applyLoadDomain -> updateDomain (DisplacementControl.cpp:287-289).
// This class puts that ordering in front of LoadControl.
//
// NOT a subclass of LoadControl: stock LoadControl's members are private, and
// LadrunoArcLength sets the fork precedent of deriving from StaticIntegrator
// directly. That keeps the vanilla footprint at ZERO (no header promotion).
// The stock behaviour below is replicated verbatim so that `-extrapolate 0.0`
// is bit-identical to `integrator LoadControl`.
//
// See Ladruno_implementation/80_ladruno_sp_imposition_strengthening_adr.md

#ifndef LadrunoLoadControl_h
#define LadrunoLoadControl_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

class LadrunoLoadControl : public StaticIntegrator
{
  public:
    LadrunoLoadControl(double deltaLambda, int numIncr,
                       double minLambda, double maxLambda,
                       double extrapolateFrac = 0.0,
                       bool tangentPredictor = false);

    ~LadrunoLoadControl();

    int newStep(void);
    int update(const Vector &deltaU);
    int commit(void);
    int domainChanged(void);
    int formUnbalance(void);
    int setDeltaLambda(double newDeltaLambda);

    // Ladruno (ADR-80 S1): runtime accessors for `ladrunoLoadControl`. An
    // adaptive caller MUST resize the step through setDeltaLambda() rather than
    // re-issuing `integrator LadrunoLoadControl ...`, because re-issuing
    // CONSTRUCTS A NEW OBJECT and destroys the predictor state -- measured: the
    // predictor then never fires at all. See 80c.
    double getDeltaLambda(void) const { return deltaLambda; }
    double getExtrapolate(void) const { return extrapolate; }
    bool   predictorArmed(void) const { return havePrev && !justCutBack; }
    bool   getTangentPredictor(void) const { return tangentPredictReq; }
    int    getTangentPredictCount(void) const { return numTangentPredicts; }

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    // ---- stock LoadControl state (replicated verbatim) --------------------
    double deltaLambda;                       // dlambda at step (i-1)
    double specNumIncrStep, numIncrLastStep;  // Jd & J(i-1)
    double dLambdaMin, dLambdaMax;            // min & max dlambda

    // ---- Ladruno (ADR-80 S1): predictor state ----------------------------
    // frac == 0.0 disables every line of it; the class is then stock.
    double extrapolate;      // user fraction; 1.0 == Abaqus EXTRAPOLATION=LINEAR
    Vector *dUprev;          // accumulated dU of the last COMMITTED step
    Vector *dUaccum;         // dU accumulated so far in the OPEN step
    double  dLambdaPrev;     // the dlambda that produced dUprev
    bool    havePrev;        // dUprev/dLambdaPrev are meaningful
    bool    stepOpen;        // newStep() ran without an intervening commit()
    bool    justCutBack;     // the previous attempt at this step FAILED

    // ---- Ladruno (ADR-80 P3): the TANGENT predictor -----------------------
    // The Kratos `b -= K*du_D` route (use_old_stiffness_in_first_iteration),
    // and the answer to why -extrapolate failed its acceptance gate (80c):
    // linear extrapolation REDUCES the driven layer's overstrain, it does not
    // REMOVE it, and what survives still trips the spurious yield.
    //
    // newStep() applies the domain loads -- which SETS every sp value -- but
    // does NOT let the constraint handler enforce them, so the elements stay at
    // the COMMITTED state. Iteration 1 therefore forms K at the committed state
    // and gets its prescribed-motion forcing from -K*du_D assembled directly,
    // with NO constitutive evaluation in the lagging state anywhere. The sps
    // are enforced in update(), on an interior the predictor solve has already
    // distributed.
    //
    // Unlike -extrapolate this carries NO state across steps: it fires on the
    // first increment, after a cutback, and when an adaptive caller re-issues
    // `integrator LadrunoLoadControl` every step -- the three holes in S1.
    bool    tangentPredictReq;   // what the user asked for (-tangentPredictor)
    bool    tangentPredict;      // what is ACTIVE (cleared by the fallback below)
    bool    spNotYetEnforced;    // between newStep() and the step's first update()
    bool    warnedNoSPTerm;      // the "nothing contributed" note fired once
    double  stepLambda;          // lambda applied by newStep()
    int     numTangentPredicts;  // diagnostic: increments that used the route

    void invalidatePredictor(void);   // after domainChanged / recvSelf
};

#endif
