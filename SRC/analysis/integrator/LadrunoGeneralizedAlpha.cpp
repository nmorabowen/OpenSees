/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno fork — N. Mora-Bowen
// ADR-52 W3-I2: sensitivity-carrying (DDM) subclass of GeneralizedAlpha.
//
// IMPORTANT — base-class tangent inconsistency (drives the DDM design here):
// OpenSees GeneralizedAlpha::update() sets the model acceleration to *Udotdot
// (GeneralizedAlpha.cpp), so the PRIMAL dynamic residual integrates inertia at
// the FULL step: R = F(t+alphaF*dt) - P(Ualpha) - C*Ualphadot - M*Udotdot. Its
// consistent Jacobian is ∂R/∂U = alphaF*K + alphaF*c2*C + c3*M (M-coef c3). But
// GeneralizedAlpha::formEleTangent emits alphaF*K + alphaF*c2*C + alphaM*c3*M
// (M-coef alphaM*c3) — i.e. the base TANGENT is inconsistent with its own
// residual for alphaM != 1 (a long-standing base quirk; harmless to the converged
// LINEAR solution but it degrades Newton's rate and breaks tangent-reuse DDM).
//
// Consequences for DDM (validated by the FD oracle, which failed ~2e-3 when the
// sensitivity residual/tangent were built tangent-consistent instead of
// primal-consistent):
//   * the sensitivity RESIDUAL must match the PRIMAL — M acts at Udotdot with NO
//     alphaM weighting (exactly like Newmark/LadrunoHHT). K and C genuinely act at
//     Ualpha/Ualphadot, so their alphaF weighting + the -K*(1-alphaF)*dUn term stay.
//   * the sensitivity SOLVE must use the PRIMAL Jacobian (M-coef c3), so we
//     RE-FORM the tangent with c3*M (sensTangentFlag) in computeSensitivities
//     rather than reuse the inconsistent factored alphaM*c3*M primal tangent.
// The non-sensitivity (primal) path delegates to the base unchanged ⇒ a plain
// transient run is byte-identical to vanilla GeneralizedAlpha. See
// Ladruno_implementation/LEDGER_quirks.md for the base-class note.

#include <LadrunoGeneralizedAlpha.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <elementAPI.h>


void *
OPS_LadrunoGeneralizedAlpha(void)
{
    TransientIntegrator *theIntegrator = 0;

    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 2 && argc != 4) {
        opserr << "WARNING - incorrect number of args want LadrunoGeneralizedAlpha $alphaM $alphaF <$gamma $beta>\n";
        return 0;
    }

    double dData[4];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want LadrunoGeneralizedAlpha $alphaM $alphaF <$gamma $beta>\n";
        return 0;
    }

    if (argc == 2)
        theIntegrator = new LadrunoGeneralizedAlpha(dData[0], dData[1]);
    else
        theIntegrator = new LadrunoGeneralizedAlpha(dData[0], dData[1], dData[2], dData[3]);

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating LadrunoGeneralizedAlpha integrator\n";

    return theIntegrator;
}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha()
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, 0.0, 0.0, 0.25, 0.5),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha(double _alphaM, double _alphaF)
    // same beta/gamma the base derives from (alphaM, alphaF)
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, _alphaM, _alphaF,
                       (1+_alphaM-_alphaF)*(1+_alphaM-_alphaF)*0.25, 0.5+_alphaM-_alphaF),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::LadrunoGeneralizedAlpha(double _alphaM, double _alphaF,
                                                 double _beta, double _gamma)
    : GeneralizedAlpha(INTEGRATOR_TAGS_LadrunoGeneralizedAlpha, _alphaM, _alphaF, _beta, _gamma),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS(), sensTangentFlag(0)
{

}


LadrunoGeneralizedAlpha::~LadrunoGeneralizedAlpha()
{
    if (massMatrixMultiplicator != 0)
        delete massMatrixMultiplicator;
    if (dampingMatrixMultiplicator != 0)
        delete dampingMatrixMultiplicator;
}


int
LadrunoGeneralizedAlpha::formEleResidual(FE_Element *theEle)
{
    if (sensitivityFlag == 0) {  // no sensitivity -> identical to the base path
        this->TransientIntegrator::formEleResidual(theEle);
    } else {

        theEle->zeroResidual();

        // Generalized-alpha residual is evaluated at:
        //   K, C at the alphaF state (Ualpha, Ualphadot); M at the alphaM state
        //   (Ualphadotdot). U/Udot/Udotdot evolve by the SAME Newmark relations,
        //   so the a2..a8 constants are the Newmark ones:
        //     dUdotdot/dh = c3*dU/dh + a2*dUn + a3*dVn + a4*dAn
        //     dUdot/dh    = c2*dU/dh + a6*dUn + a7*dVn + a8*dAn
        //   and the alpha-intermediate sensitivities:
        //     dUalpha/dh       = (1-alphaF)*dUn + alphaF*dU/dh
        //     dUalphadot/dh    = (1-alphaF)*dVn + alphaF*dUdot/dh
        //     dUalphadotdot/dh = (1-alphaM)*dAn + alphaM*dUdotdot/dh
        //   The dU/dh chain terms (alphaF*K, alphaF*c2*C, alphaM*c3*M) match
        //   formEleTangent exactly and stay on the LHS.
        double a2 = -c3;
        double a3 = -c2/gamma;
        double a4 = 1.0 - 1.0/(2.0*beta);
        double a6 = -c2;
        double a7 = 1.0 - gamma/beta;
        double dt = gamma/(beta*c2);
        double a8 = dt*(1.0 - gamma/(2.0*beta));

        int vectorSize = U->Size();
        Vector dUn(vectorSize);
        Vector dVn(vectorSize);
        Vector dAn(vectorSize);
        int i, loc;

        AnalysisModel *myModel = this->getAnalysisModel();
        DOF_GrpIter &theDOFs = myModel->getDOFs();
        DOF_Group *dofPtr;
        while ((dofPtr = theDOFs()) != 0) {
            const ID &id = dofPtr->getID();
            int idSize = id.Size();

            const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dUn(loc) = dispSens(i);
            }
            const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dVn(loc) = velSens(i);
            }
            const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
            for (i = 0; i < idSize; i++) {
                loc = id(i);
                if (loc >= 0)
                    dAn(loc) = accelSens(i);
            }
        }

        // mass multiplicator: the PRIMAL integrates inertia at the FULL step
        // (M*Udotdot, see the file header) — so this is the Newmark form with NO
        // alphaM weighting (dUdotdot/dh = c3*dU/dh + a2*dUn + a3*dVn + a4*dAn; the
        // c3*dU/dh chain term goes to the re-formed consistent LHS):
        //   tmp1 = a2*dUn + a3*dVn + a4*dAn
        Vector tmp1(vectorSize);
        tmp1.addVector(0.0, dUn, a2);
        tmp1.addVector(1.0, dVn, a3);
        tmp1.addVector(1.0, dAn, a4);

        // damping multiplicator: C acts at Ualphadot, whose non-(dU/dh) part is
        //   (1-alphaF)*dVn + alphaF*(a6*dUn + a7*dVn + a8*dAn)
        Vector tmp2(vectorSize);
        tmp2.addVector(0.0, dUn, alphaF*a6);
        tmp2.addVector(1.0, dVn, (1.0-alphaF) + alphaF*a7);
        tmp2.addVector(1.0, dAn, alphaF*a8);

        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(tmp1.Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(tmp2.Size());

        (*massMatrixMultiplicator) = tmp1;
        (*dampingMatrixMultiplicator) = tmp2;

        // -dPint/dh|u fixed (element evaluated at Ualpha, the current trial state)
        theEle->addResistingForceSensitivity(gradNumber);

        // -dM/dh*acc  (acc at Udotdot — the PRIMAL inertia state, NOT Ualphadotdot)
        theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

        // generalized-alpha stiffness term: -K*(1-alphaF)*dUn (the (1-alphaF) part
        // of d(Ualpha)/dh; addK_Force uses the current consistent tangent)
        theEle->addK_Force(dUn, -(1.0-alphaF));

        // -M*(alphaM-weighted mass multiplicator)
        theEle->addM_Force(*massMatrixMultiplicator, -1.0);

        // -C*(alphaF-weighted damping multiplicator)
        theEle->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh*vel  (vel at Ualphadot)
        theEle->addD_ForceSensitivity(gradNumber, *Ualphadot, -1.0);
    }

    return 0;
}


int
LadrunoGeneralizedAlpha::formNodUnbalance(DOF_Group *theDof)
{
    if (sensitivityFlag == 0) {  // no sensitivity -> identical to the base path
        this->TransientIntegrator::formNodUnbalance(theDof);
    } else {

        theDof->zeroUnbalance();

        // Safety: the global multiplicators are built in formEleResidual; a model
        // with nodal mass/damping but ZERO contributing FE_Elements never runs it,
        // so allocate-and-zero here rather than deref null (a fresh Vector(size) is
        // zero-initialized). Normal models with elements never hit this branch.
        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(U->Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(U->Size());

        // -M*(mass multiplicator, full-step / no-alphaM form set by formEleResidual)
        theDof->addM_Force(*massMatrixMultiplicator, -1.0);

        // -dM/dh*acc  (acc at Udotdot — the PRIMAL inertia state)
        theDof->addM_ForceSensitivity(*Udotdot, -1.0);

        // -C*(alphaF-weighted damping multiplicator)
        theDof->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh*vel  (vel at Ualphadot)
        theDof->addD_ForceSensitivity(*Ualphadot, -1.0);

        // random nodal loads (already formed by applyLoadSensitivity)
        theDof->addPtoUnbalance();
    }

    return 0;
}


int
LadrunoGeneralizedAlpha::formEleTangent(FE_Element *theEle)
{
    if (sensTangentFlag == 0)   // primal path -> base tangent (byte-identical)
        return this->GeneralizedAlpha::formEleTangent(theEle);

    // DDM solve path: the PRIMAL-CONSISTENT Jacobian of R = F - P(Ualpha) -
    // C*Ualphadot - M*Udotdot, i.e. alphaF*K + alphaF*c2*C + c3*M (M-coef c3, NOT
    // the base's alphaM*c3). Re-formed in computeSensitivities so the sensitivity
    // solve uses the true ∂R/∂U; the primal run never sees this branch.
    theEle->zeroTangent();
    if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(alphaF*c1);
    else
        theEle->addKtToTang(alphaF*c1);
    theEle->addCtoTang(alphaF*c2);
    theEle->addMtoTang(c3);

    return 0;
}


int
LadrunoGeneralizedAlpha::formNodTangent(DOF_Group *theDof)
{
    if (sensTangentFlag == 0)   // primal path -> base tangent (byte-identical)
        return this->GeneralizedAlpha::formNodTangent(theDof);

    // DDM solve path: nodal M/C with the primal-consistent M-coef c3.
    theDof->zeroTangent();
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(c3);

    return 0;
}


int
LadrunoGeneralizedAlpha::revertToStart()
{
    if (Ut != 0)         Ut->Zero();
    if (Utdot != 0)      Utdot->Zero();
    if (Utdotdot != 0)   Utdotdot->Zero();
    if (U != 0)          U->Zero();
    if (Udot != 0)       Udot->Zero();
    if (Udotdot != 0)    Udotdot->Zero();
    if (Ualpha != 0)     Ualpha->Zero();
    if (Ualphadot != 0)  Ualphadot->Zero();
    if (Ualphadotdot != 0) Ualphadotdot->Zero();

    return 0;
}


int
LadrunoGeneralizedAlpha::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;
    gradNumber = passedGradNumber;

    LinearSOE *theSOE = this->getLinearSOE();

    if (assemblyFlag != 0)
        theSOE->setB(independentRHS);

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();

    Node *nodePtr;
    NodeIter &theNodeIter = theDomain->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
        nodePtr->zeroUnbalancedLoad();

    // randomness in external load (incl. time series).
    // NOTE on the evaluation time: the generalized-alpha dynamic residual is
    // formed at t_n + alphaF*deltaT. With the standard -computeAtEachStep
    // (pre-commit) sensitivity path getCurrentTime() returns exactly that, so the
    // load gradient is consistent. A post-commit -computeByCommand caller would see
    // t_n + deltaT (off by (1-alphaF)*deltaT) — only affects RANDOM-LOAD
    // sensitivity (dF/dh != 0); material-parameter DDM (dF/dh == 0) is unaffected.
    LoadPattern *loadPatternPtr;
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    double time;
    while ((loadPatternPtr = thePatterns()) != 0) {
        time = theDomain->getCurrentTime();
        loadPatternPtr->applyLoadSensitivity(time);
    }

    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while ((elePtr = theEles()) != 0)
        theSOE->addB(elePtr->getResidual(this), elePtr->getID());

    DOF_Group *dofPtr;
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    while ((dofPtr = theDOFs()) != 0)
        theSOE->addB(dofPtr->getUnbalance(this), dofPtr->getID());

    sensitivityFlag = 0;

    return 0;
}


int
LadrunoGeneralizedAlpha::formIndependentSensitivityRHS()
{
    // mirrors Newmark: not used for now
    return 0;
}


int
LadrunoGeneralizedAlpha::saveSensitivity(const Vector &vNew, int gradNum, int numGrads)
{
    // U / Udot / Udotdot evolve by the SAME Newmark relations, so the propagation
    // of the response sensitivities is identical to Newmark (no alpha here).
    double a1 = c3;
    double a2 = -c3;
    double a3 = -c2/gamma;
    double a4 = 1.0 - 1.0/(2.0*beta);
    double a5 = c2;
    double a6 = -c2;
    double a7 = 1.0 - gamma/beta;
    double dt = gamma/(beta*c2);
    double a8 = dt*(1.0 - gamma/(2.0*beta));

    int vectorSize = U->Size();
    Vector dUn(vectorSize);
    Vector dVn(vectorSize);
    Vector dAn(vectorSize);
    int i, loc;

    AnalysisModel *myModel = this->getAnalysisModel();
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dUn(loc) = dispSens(i);
        }
        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dVn(loc) = velSens(i);
        }
        const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
        for (i = 0; i < idSize; i++) {
            loc = id(i);
            if (loc >= 0)
                dAn(loc) = accelSens(i);
        }
    }

    Vector vdotNew(vectorSize);
    Vector vdotdotNew(vectorSize);
    vdotdotNew.addVector(0.0, vNew, a1);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    vdotNew.addVector(0.0, vNew, a5);
    vdotNew.addVector(1.0, dUn, a6);
    vdotNew.addVector(1.0, dVn, a7);
    vdotNew.addVector(1.0, dAn, a8);

    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group *dofPtr1;
    while ((dofPtr1 = theDOFGrps()) != 0)
        dofPtr1->saveSensitivity(vNew, vdotNew, vdotdotNew, gradNum, numGrads);

    return 0;
}


int
LadrunoGeneralizedAlpha::commitSensitivity(int gradNum, int numGrads)
{
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while ((elePtr = theEles()) != 0)
        elePtr->commitSensitivity(gradNum, numGrads);

    return 0;
}


int
LadrunoGeneralizedAlpha::computeSensitivities()
{
    LinearSOE *theSOE = this->getLinearSOE();

    theSOE->zeroB();
    this->formIndependentSensitivityRHS();

    // Re-form the system tangent with the PRIMAL-CONSISTENT Jacobian (M-coef c3,
    // not the base's alphaM*c3) so every sensitivity solve below uses the true
    // ∂R/∂U. The base GeneralizedAlpha tangent that the primal Newton left
    // factored in the SOE is inconsistent for alphaM != 1 (see the file header);
    // reusing it (as Newmark's computeSensitivities does, where it IS consistent)
    // would bias the gradients (FD oracle: ~2e-3). Re-formed once here; the
    // per-parameter solves below reuse this factorization. The flag is cleared
    // immediately so nothing else sees the sensitivity tangent.
    sensTangentFlag = 1;
    this->formTangent(CURRENT_TANGENT);
    sensTangentFlag = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    while ((theParam = paramIter()) != 0)
        theParam->activate(false);

    int numGrads = theDomain->getNumParameters();
    paramIter = theDomain->getParameters();
    while ((theParam = paramIter()) != 0) {
        theParam->activate(true);

        theSOE->zeroB();
        int gradIndex = theParam->getGradIndex();
        this->formSensitivityRHS(gradIndex);
        theSOE->solve();
        this->saveSensitivity(theSOE->getX(), gradIndex, numGrads);
        this->commitSensitivity(gradIndex, numGrads);

        theParam->activate(false);
    }

    return 0;
}


const char *
LadrunoGeneralizedAlpha::getClassType(void) const
{
    return "LadrunoGeneralizedAlpha";
}


void
LadrunoGeneralizedAlpha::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "LadrunoGeneralizedAlpha (sensitivity/DDM) - currentTime: " << currentTime << endln;
        s << "  alphaM: " << alphaM << "  alphaF: " << alphaF;
        s << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "LadrunoGeneralizedAlpha - no associated AnalysisModel\n";
}
