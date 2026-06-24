/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno fork — N. Mora-Bowen
// ADR-52 W3-I2: sensitivity-carrying (DDM) subclass of HHT. See LadrunoHHT.h
// for the derivation of the alpha-weighted sensitivity residual.

#include <LadrunoHHT.h>
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
OPS_LadrunoHHT(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;

    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 3) {
        opserr << "WARNING - incorrect number of args want LadrunoHHT $alpha <$gamma $beta>\n";
        return 0;
    }

    double dData[3];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want LadrunoHHT $alpha <$gamma $beta>\n";
        return 0;
    }

    if (argc == 1)
        theIntegrator = new LadrunoHHT(dData[0]);
    else
        theIntegrator = new LadrunoHHT(dData[0], dData[2], dData[1]); // alpha, beta, gamma for ctor

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating LadrunoHHT integrator\n";

    return theIntegrator;
}


LadrunoHHT::LadrunoHHT()
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, 1.0, 0.25, 0.5),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::LadrunoHHT(double _alpha)
    // same beta/gamma defaults HHT derives from a single alpha
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, _alpha, (2-_alpha)*(2-_alpha)*0.25, 1.5-_alpha),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::LadrunoHHT(double _alpha, double _beta, double _gamma)
    : HHT(INTEGRATOR_TAGS_LadrunoHHT, _alpha, _beta, _gamma),
      sensitivityFlag(0), gradNumber(0),
      massMatrixMultiplicator(0), dampingMatrixMultiplicator(0),
      assemblyFlag(0), independentRHS()
{

}


LadrunoHHT::~LadrunoHHT()
{
    if (massMatrixMultiplicator != 0)
        delete massMatrixMultiplicator;
    if (dampingMatrixMultiplicator != 0)
        delete dampingMatrixMultiplicator;
}


int
LadrunoHHT::formEleResidual(FE_Element *theEle)
{
    if (sensitivityFlag == 0) {  // no sensitivity -> identical to the base HHT path
        this->TransientIntegrator::formEleResidual(theEle);
    } else {

        theEle->zeroResidual();

        // HHT residual is evaluated at the alpha-intermediate state:
        //   R = F(t+alpha*dt) - P(Ualpha) - C*Ualphadot - M*Udotdot
        // U / Udot / Udotdot evolve by the SAME Newmark relations, so the
        // a2..a8 constants are the Newmark ones (see Newmark::formEleResidual):
        //   dUdotdot/dh = c3*dU/dh + a2*dUn + a3*dVn + a4*dAn
        //   dUdot/dh    = c2*dU/dh + a6*dUn + a7*dVn + a8*dAn
        double a2 = -c3;
        double a3 = -c2/gamma;
        double a4 = 1.0 - 1.0/(2.0*beta);
        double a6 = -c2;
        double a7 = 1.0 - gamma/beta;
        double dt = gamma/(beta*c2);
        double a8 = dt*(1.0 - gamma/(2.0*beta));

        // gather the previous-step response sensitivities dUn/dVn/dAn
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

        // mass multiplicator: M acts at the FULL step -> Newmark form (no alpha)
        //   tmp1 = a2*dUn + a3*dVn + a4*dAn
        Vector tmp1(vectorSize);
        tmp1.addVector(0.0, dUn, a2);
        tmp1.addVector(1.0, dVn, a3);
        tmp1.addVector(1.0, dAn, a4);

        // damping multiplicator: C acts at Ualphadot, whose non-(dU/dh) part is
        //   (1-alpha)*dVn + alpha*(a6*dUn + a7*dVn + a8*dAn)
        Vector tmp2(vectorSize);
        tmp2.addVector(0.0, dUn, alpha*a6);
        tmp2.addVector(1.0, dVn, (1.0-alpha) + alpha*a7);
        tmp2.addVector(1.0, dAn, alpha*a8);

        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(tmp1.Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(tmp2.Size());

        (*massMatrixMultiplicator) = tmp1;
        (*dampingMatrixMultiplicator) = tmp2;

        // -dPint/dh|u fixed (element evaluated at Ualpha, the current trial state)
        theEle->addResistingForceSensitivity(gradNumber);

        // -dM/dh*acc  (acc at full step)
        theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

        // HHT-only term: -K*(1-alpha)*dUn  (the (1-alpha) part of d(Ualpha)/dh;
        // addK_Force uses the current consistent tangent, matching the LHS).
        theEle->addK_Force(dUn, -(1.0-alpha));

        // -M*(a2*dUn + a3*dVn + a4*dAn)
        theEle->addM_Force(*massMatrixMultiplicator, -1.0);

        // -C*(alpha-weighted damping multiplicator)
        theEle->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh*vel  (vel at Ualphadot)
        theEle->addD_ForceSensitivity(gradNumber, *Ualphadot, -1.0);
    }

    return 0;
}


int
LadrunoHHT::formNodUnbalance(DOF_Group *theDof)
{
    if (sensitivityFlag == 0) {  // no sensitivity -> identical to the base HHT path
        this->TransientIntegrator::formNodUnbalance(theDof);
    } else {

        theDof->zeroUnbalance();

        // Safety: the global multiplicators are built in formEleResidual. A model
        // with nodal mass/damping but ZERO contributing FE_Elements never runs it,
        // so allocate-and-zero here rather than deref a null pointer (a fresh
        // Vector(size) is zero-initialized ⇒ the element contribution is correctly
        // treated as zero). Normal models with elements never hit this branch.
        if (massMatrixMultiplicator == 0)
            massMatrixMultiplicator = new Vector(U->Size());
        if (dampingMatrixMultiplicator == 0)
            dampingMatrixMultiplicator = new Vector(U->Size());

        // nodal contributions carry only M and C (no element stiffness); the
        // global multiplicators were set by formEleResidual (same for all DOFs).

        // -M*(a2*dUn + a3*dVn + a4*dAn)
        theDof->addM_Force(*massMatrixMultiplicator, -1.0);

        // -dM/dh*acc
        theDof->addM_ForceSensitivity(*Udotdot, -1.0);

        // -C*(alpha-weighted damping multiplicator)
        theDof->addD_Force(*dampingMatrixMultiplicator, -1.0);

        // -dC/dh*vel  (vel at Ualphadot)
        theDof->addD_ForceSensitivity(*Ualphadot, -1.0);

        // random nodal loads (already formed by applyLoadSensitivity)
        theDof->addPtoUnbalance();
    }

    return 0;
}


int
LadrunoHHT::revertToStart()
{
    if (Ut != 0)       Ut->Zero();
    if (Utdot != 0)    Utdot->Zero();
    if (Utdotdot != 0) Utdotdot->Zero();
    if (U != 0)        U->Zero();
    if (Udot != 0)     Udot->Zero();
    if (Udotdot != 0)  Udotdot->Zero();
    if (Ualpha != 0)   Ualpha->Zero();
    if (Ualphadot != 0) Ualphadot->Zero();

    return 0;
}


int
LadrunoHHT::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;
    gradNumber = passedGradNumber;

    LinearSOE *theSOE = this->getLinearSOE();

    // optional independent part of the RHS
    if (assemblyFlag != 0)
        theSOE->setB(independentRHS);

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();

    // zero nodal unbalanced loads
    Node *nodePtr;
    NodeIter &theNodeIter = theDomain->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
        nodePtr->zeroUnbalancedLoad();

    // randomness in external load (incl. time series).
    // NOTE on the evaluation time: the HHT dynamic residual is formed at
    // t_n + alpha*deltaT. With the standard -computeAtEachStep (pre-commit)
    // sensitivity path, getCurrentTime() returns exactly that, so the load
    // gradient is consistent with the residual. A post-commit -computeByCommand
    // caller would instead see t_n + deltaT (off by (1-alpha)*deltaT) — this only
    // affects RANDOM-LOAD sensitivity (dF/dh != 0); material-parameter sensitivity
    // (dF/dh == 0, the fragility/FORM-on-material case) is unaffected either way.
    LoadPattern *loadPatternPtr;
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    double time;
    while ((loadPatternPtr = thePatterns()) != 0) {
        time = theDomain->getCurrentTime();
        loadPatternPtr->applyLoadSensitivity(time);
    }

    // randomness in element/material contributions
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while ((elePtr = theEles()) != 0)
        theSOE->addB(elePtr->getResidual(this), elePtr->getID());

    // DOF groups (MUST be last)
    DOF_Group *dofPtr;
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    while ((dofPtr = theDOFs()) != 0)
        theSOE->addB(dofPtr->getUnbalance(this), dofPtr->getID());

    sensitivityFlag = 0;

    return 0;
}


int
LadrunoHHT::formIndependentSensitivityRHS()
{
    // mirrors Newmark: not used for now
    return 0;
}


int
LadrunoHHT::saveSensitivity(const Vector &vNew, int gradNum, int numGrads)
{
    // U / Udot / Udotdot evolve by the SAME Newmark relations as Newmark, so the
    // propagation of the response sensitivities is identical (no alpha here).
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
    // vdotdotNew = vNew*a1 + dUn*a2 + dVn*a3 + dAn*a4
    vdotdotNew.addVector(0.0, vNew, a1);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    // vdotNew = vNew*a5 + dUn*a6 + dVn*a7 + dAn*a8
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
LadrunoHHT::commitSensitivity(int gradNum, int numGrads)
{
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();
    while ((elePtr = theEles()) != 0)
        elePtr->commitSensitivity(gradNum, numGrads);

    return 0;
}


int
LadrunoHHT::computeSensitivities()
{
    LinearSOE *theSOE = this->getLinearSOE();

    theSOE->zeroB();
    this->formIndependentSensitivityRHS();

    AnalysisModel *theModel = this->getAnalysisModel();
    Domain *theDomain = theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    // de-activate all parameters
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
LadrunoHHT::getClassType(void) const
{
    return "LadrunoHHT";
}


void
LadrunoHHT::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "LadrunoHHT (sensitivity/DDM) - currentTime: " << currentTime << endln;
        s << "  alpha: " << alpha;
        s << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "LadrunoHHT - no associated AnalysisModel\n";
}
