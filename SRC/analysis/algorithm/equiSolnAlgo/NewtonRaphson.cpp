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
                                                                        
// $Revision: 1.10 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/NewtonRaphson.C 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// NewtonRaphson. NewtonRaphson is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)NewtonRaphson.C, revA"

#include <NewtonRaphson.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>
#include <elementAPI.h>
#include <profiler/ProfilerMacros.h>
#include <string>


void* OPS_NewtonRaphsonAlgorithm()
{
    int formTangent = CURRENT_TANGENT;
    double iFactor = 0;
    double cFactor = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* type = OPS_GetString();
      if(strcmp(type,"-secant")==0 || strcmp(type,"-Secant")==0) {
	formTangent = CURRENT_SECANT;
	iFactor = 0;
	cFactor = 1.0;
      } else if(strcmp(type,"-initial")==0 || strcmp(type,"-Initial")==0) {
	formTangent = INITIAL_TANGENT;
	iFactor = 1.;
	cFactor = 0;
      } else if(strcmp(type,"-intialThenCurrent")==0 || strcmp(type,"-intialCurrent")==0) {
	formTangent = INITIAL_THEN_CURRENT_TANGENT;
	iFactor = 0;
	cFactor = 1.0;
      } else if(strcmp(type,"-hall")==0 || strcmp(type,"-Hall")==0) {
	formTangent = HALL_TANGENT;
	iFactor = 0.1;
	cFactor = 0.9;
	if (OPS_GetNumRemainingInputArgs() == 2) {
	  double data[2];
	  int numData = 2;
	  if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	    opserr << "WARNING invalid data reading 2 hall factors\n";
	    return 0;
	  }
	  iFactor = data[0];
	  cFactor = data[1];
	}
      }
    }

    return new NewtonRaphson(formTangent, iFactor, cFactor);

}

// Constructor
NewtonRaphson::NewtonRaphson(int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{

}

NewtonRaphson::NewtonRaphson()
	:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
	tangent(CURRENT_TANGENT), iFactor(0.), cFactor(1.)
{

}


NewtonRaphson::NewtonRaphson(ConvergenceTest &theT, int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{

}

// Destructor
NewtonRaphson::~NewtonRaphson()
{
  

}


int 
NewtonRaphson::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
    //IncrementalIntegrator *theIntegratorSens=this->getIncrementalIntegratorPtr();//Abbas
    LinearSOE  *theSOE = this->getLinearSOEptr();

    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING NewtonRaphson::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    { OPS_PROFILE_SCOPE("formUnbalance");
    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";
      return -2;
    }
    }

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "NewtonRaphson::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    int result = -1;
    numIterations = 0;

    // Ladruno ADR-76 (R1 — documentation only, no behaviour change here):
    //
    // `formTangent` below sits INSIDE the iteration loop, so `algorithm Newton
    // -initial` performs a full tangent assembly AND a full numeric
    // factorization on every iteration. `INITIAL_THEN_CURRENT_TANGENT` is the
    // only flag given special treatment (immediately below); `INITIAL_TANGENT`
    // falls through to the generic branch and is re-formed unconditionally.
    //
    // Under `StaticIntegrator` that re-formation is `zeroTangent(); addKiToTang();`
    // (StaticIntegrator.cpp:75-76) — for MANY models an invariant matrix, rebuilt
    // and re-factorized once per iteration for the whole analysis.
    //
    // READ THE CAVEAT BEFORE THE RECOMMENDATION. "Initial tangent" is a naming
    // convention here, not a contract, so the matrix is NOT always invariant:
    //   * `addKiToTang()` reaches `Element::getInitialStiff()`, which for a 3D
    //     corotational transformation is material initial stiffness on the
    //     CURRENT geometry — CorotCrdTransf3d::getInitialGlobalStiffMatrix
    //     (:1452) triple-products through the member `T`, which `update()`
    //     rebuilds at its tail (:549). LadrunoContactFE::addKiToTang likewise
    //     reads the current active set.
    //   * under a transient integrator the `c2*C` term carries `betaK*Kt`
    //     (Newmark.cpp:295-298, SRC/element/Element.cpp:222-223), so with
    //     `betaK != 0` the matrix genuinely differs between iterations.
    // In those models `Newton -initial` and `ModifiedNewton -initial` are NOT the
    // same algorithm: the former re-linearizes on the current configuration each
    // iteration, the latter freezes at step start, so iteration counts and
    // convergence can differ.
    //
    // ==> Outside those cases — i.e. when the assembled initial tangent really is
    //     state-independent — `algorithm ModifiedNewton -initial` is the cheap
    //     spelling of the same iteration: one assembly and one factorization per
    //     step, forming the tangent BEFORE the loop (see
    //     ModifiedNewton::solveCurrentStep). Add `-factorOnce` for one per
    //     analysis, but see the staleness trap in LEDGER_quirks first.
    //
    // The cost, so the number is not mistaken for a redundancy proof: on a
    // 39k-DOF SSPbrickUP/ShellMITC4 model with MKL PARDISO, `Newton -initial`
    // cost 1.504 s/step against `ModifiedNewton -initial` at 0.936 s/step (1.61x)
    // for an identical converged answer in an identical iteration count. That run
    // was TRANSIENT with betaK != 0 — one of the cases above where the extra
    // assembly is defensible. It measures what the extra work COSTS, not that it
    // was redundant; it simply bought nothing there.
    //
    // Not fixed here on purpose: the skip is only sound when invariance is known,
    // and that is an ELEMENT property, not an integrator-family one. See ADR-76
    // for the tangent-version design that makes the skip safe, and for the
    // audit of which transformations are and are not configuration-dependent.
    do {

      { OPS_PROFILE_SCOPE("formTangent");
      if (tangent == INITIAL_THEN_CURRENT_TANGENT) {
	if (numIterations == 0) {
	  SOLUTION_ALGORITHM_tangentFlag = INITIAL_TANGENT;
	  if (theIntegrator->formTangent(INITIAL_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  }
	} else {
	  SOLUTION_ALGORITHM_tangentFlag = CURRENT_TANGENT;
	  if (theIntegrator->formTangent(CURRENT_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  }
	}
      }	else {

	SOLUTION_ALGORITHM_tangentFlag = tangent;
	if (theIntegrator->formTangent(tangent, iFactor, cFactor) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	}
      }
      }
      { OPS_PROFILE_SCOPE("linearSolve");
      if (theSOE->solve() < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the LinearSysOfEqn failed in solve()\n";
	return -3;
      }
      }

      { OPS_PROFILE_SCOPE("update");
      if (theIntegrator->update(theSOE->getX()) < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";
	return -4;
      }
      }
      { OPS_PROFILE_SCOPE("formUnbalance");
      if (theIntegrator->formUnbalance() < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";
	return -2;
      }
      }

      { OPS_PROFILE_SCOPE("convTest"); result = theTest->test(); }
       numIterations++;
      this->record(numIterations);

    } while (result == -1);

    if (result == -2) {
      opserr << "NewtonRaphson::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }
// note - if positive result we are returning what the convergence test returned
    // which should be the number of iterations
    
        return result;
}


int
NewtonRaphson::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(3);
  data(0) = tangent;
  data(1) = iFactor;
  data(2) = cFactor;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
NewtonRaphson::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  static Vector data(3);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  tangent = int(data(0));
  iFactor = data(1);
  cFactor = data(2);
  return 0;
}


void
NewtonRaphson::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "NewtonRaphson" << endln;
  }
}


int
NewtonRaphson::getNumIterations(void)
{
  return numIterations;
}


