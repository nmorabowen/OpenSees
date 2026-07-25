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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ModifiedNewton.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/ModifiedNewton.C 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
// Description: This file contains the class definition for 
// ModifiedNewton. ModifiedNewton is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)ModifiedNewton.C, revA"

#include <ModifiedNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
//#include <Timer.h>
#include <elementAPI.h>
#include <profiler/ProfilerMacros.h>  // Ladruno (ADR-40 rank 8/10): phase scopes, parity with NewtonRaphson

void* OPS_ModifiedNewton()
{
  int formTangent = CURRENT_TANGENT;
  int factoronce = 0;
  double iFactor = 0;
  double cFactor = 1;

  // Ladruno ADR-76 (R4): loop over EVERY remaining option instead of reading
  // exactly one. Upstream's `if (OPS_GetNumRemainingInputArgs() > 0)` consumed a
  // single option string, so `algorithm ModifiedNewton -initial -factoronce`
  // silently honoured `-initial` and dropped the rest — yet that pair is the
  // natural spelling of "initial stiffness, assembled and factorized once",
  // which for a static analysis is the cheapest robust algorithm the framework
  // offers. `OPS_NewtonRaphsonAlgorithm()` already uses this `while` form.
  //
  // NOT behaviour-preserving in general, and deliberately so: any deck passing
  // two or more recognised options now gets the LAST one to set a given field
  // rather than the first, which changes the iteration matrix. That is the whole
  // point of the fix, but it is a silent results change for such decks — no
  // in-tree deck passes more than one option, so the blast radius here is zero.
  // Single-option decks are byte-identical.
  //
  // Three things the adversarial review of this change forced (see the ADR):
  //   * `-factorOnce` is now accepted alongside `-factoronce`/`-FactorOnce`.
  //     Linear, ExpressNewton and both commands.cpp sites all spell it
  //     `-factorOnce`, and it is what the quirks ledger writes — so the camelCase
  //     carry-over spelling was still being silently dropped by the very fix
  //     meant to stop options being silently dropped. Same for `-Initial`/`-Secant`.
  //   * `-secant`/`-initial` now reset iFactor/cFactor the way
  //     OPS_NewtonRaphsonAlgorithm() does. With an `if` only one branch could ever
  //     run, so the reset was unnecessary; a loop is exactly the context where a
  //     later option must be able to override an earlier one's factors.
  //   * a failed `-hall` factor read no longer returns null (see below).
  //
  // Left alone: `-hall`'s trailing-factor read still gates on
  // `OPS_GetNumRemainingInputArgs() == 2`, so the factors are only picked up when
  // `-hall a b` ends the command. NewtonRaphson carries the identical gate, and
  // fixing it here alone would make two parsers for the same flag disagree —
  // which is the failure mode this fork's ledger keeps recording. The unknown-token
  // warning below makes it self-diagnosing in the meantime.
  //
  // Trap that this newly makes reachable, see [[LEDGER_quirks]]: `factorOnce`
  // has NO domainChanged reset, so `-initial -factoronce` must not be combined
  // with anything that resizes the SOE mid-analysis (element removal, contact
  // re-emission, staged construction).
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* type = OPS_GetString();
    if (strcmp(type,"-secant") == 0 || strcmp(type,"-Secant") == 0) {
      formTangent = CURRENT_SECANT;
      iFactor = 0;
      cFactor = 1.0;
    } else if (strcmp(type,"-factoronce")==0 || strcmp(type,"-factorOnce")==0 ||
               strcmp(type,"-FactorOnce")==0)  {
      factoronce = 1;
    } else if (strcmp(type,"-initial") == 0 || strcmp(type,"-Initial") == 0) {
      formTangent = INITIAL_TANGENT;
      iFactor = 1.;
      cFactor = 0;
    } else if(strcmp(type,"-hall")==0 || strcmp(type,"-Hall")==0) {
      formTangent = HALL_TANGENT;
      iFactor = 0.1;
      cFactor = 0.9;
      if (OPS_GetNumRemainingInputArgs() == 2) {
        double data[2];
        int numData = 2;
        if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          // Ladruno ADR-76: KEEP `return 0` — do not degrade to the defaults.
          //
          // An earlier cut of this fix warned-and-continued here, reasoning that
          // a null is silently swallowed on the openseespy path. The premise was
          // right; the cure was worse. The classic Tcl caller has ALWAYS handled
          // null correctly (`commands.cpp:4469`, `if (theNewtonAlgo == 0) return
          // TCL_ERROR`), so degrading turned a typo'd Tcl deck that used to abort
          // loudly into one that runs to completion with the WRONG Hall factors.
          // It also made this parser disagree with OPS_NewtonRaphsonAlgorithm on
          // the same flag's error handling — the exact objection used two hunks
          // above to leave the `== 2` arity gate alone.
          //
          // The openseespy hole is now closed at its actual source: OPS_Algorithm
          // returns -1 on a null instead of reporting success
          // (`OpenSeesCommands.cpp`), which fixes EVERY algorithm factory at once
          // rather than loosening this one.
          opserr << "WARNING ModifiedNewton - invalid data reading 2 hall factors\n";
          return 0;
        }
        iFactor = data[0];
        cFactor = data[1];
      }
    } else {
      // Ladruno ADR-76 (R4, adversarial review): the point of the option loop is
      // that nothing gets silently dropped, so an unrecognised token must say so.
      // This also makes the `-hall` arity quirk above self-diagnosing: in
      // `-hall 0.2 0.8 -factoronce` the `== 2` gate does not fire, and the two
      // orphaned numerics land here and are reported instead of vanishing.
      //
      // KNOWN LIMIT under openseespy: `type` here is the literal string
      // "Invalid String Input!" whenever the token was NOT a Python string —
      // OPS_GetString() maps a non-`PyUnicode` argument to that sentinel
      // (PythonModule.cpp:246-258 -> OpenSeesCommands.cpp:1199). So for the case
      // this warning most exists for — the orphaned numerics of
      // `-hall 0.2 0.8 -factoronce` — Python users get two warnings that name
      // nothing. The warning still FIRES, which is the point; it just cannot
      // identify the token. Reading through OPS_GetStringFromAll would fix that,
      // but it must replace the OPS_GetString() at the TOP of the loop — calling
      // it here would advance the cursor and silently eat the next option. Left
      // as-is deliberately: changing the loop's read path is a wider change than
      // this warning justifies, and Tcl (where `type` is always the real token)
      // is unaffected.
      opserr << "WARNING ModifiedNewton - ignoring unrecognised option " << type << "\n";
    }
  }
  
  return new ModifiedNewton(formTangent, iFactor, cFactor,factoronce);
}

// Constructor
ModifiedNewton::ModifiedNewton(int theTangentToUse, double iFact, double cFact, int factOnce)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact), factorOnce(factOnce)
{
  
}


ModifiedNewton::ModifiedNewton(ConvergenceTest &theT, int theTangentToUse, double iFact, double cFact, int factOnce)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact), factorOnce(factOnce)
{

}

// Destructor
ModifiedNewton::~ModifiedNewton()
{

}


int 
ModifiedNewton::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel       *theAnalysisModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIncIntegratorr = this->getIncrementalIntegratorPtr();
    LinearSOE	        *theSOE = this->getLinearSOEptr();

    if ((theAnalysisModel == 0) || (theIncIntegratorr == 0) || (theSOE == 0) || (theTest == 0)){
      opserr << "WARNING ModifiedNewton::solveCurrentStep() - setLinks() has";
      opserr << " not been called - or no ConvergenceTest has been set\n";
      return -5;
    }	

    // we form the tangent
    //    Timer timer1;
    // timer1.start();

    { OPS_PROFILE_SCOPE("formUnbalance");   // Ladruno (ADR-40 rank 8/10)
    if (theIncIntegratorr->formUnbalance() < 0) {
      opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";
      return -2;
    }
    }

    SOLUTION_ALGORITHM_tangentFlag = tangent;
    if (factorOnce!=2) {
      OPS_PROFILE_SCOPE("formTangent");   // Ladruno (ADR-40 rank 8/10)
      if (theIncIntegratorr->formTangent(tangent, iFactor, cFactor) < 0){
        opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
        opserr << "the Integrator failed in formTangent()\n";
        return -1;
      }
      if (factorOnce==1) {
        factorOnce =2;
      }
    }

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "ModifiedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    numIterations = 0;
    do {
      //Timer timer2;
      //timer2.start();
	{ OPS_PROFILE_SCOPE("linearSolve");   // Ladruno (ADR-40 rank 8/10)
	if (theSOE->solve() < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";
	    return -3;
	}
	}
	//timer2.pause();
	//opserr << "TIMER::SOLVE()- " << timer2;

	{ OPS_PROFILE_SCOPE("update");   // Ladruno (ADR-40 rank 8/10)
	if (theIncIntegratorr->update(theSOE->getX()) < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in update()\n";
	    return -4;
	}
	}

	{ OPS_PROFILE_SCOPE("formUnbalance");   // Ladruno (ADR-40 rank 8/10)
	if (theIncIntegratorr->formUnbalance() < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in formUnbalance()\n";
	    return -2;
	}
	}

	this->record(numIterations++);
	{ OPS_PROFILE_SCOPE("convTest"); result = theTest->test(); }   // Ladruno (ADR-40 rank 8/10)


    } while (result == -1);

    //timer1.pause();
    //opserr << "TIMER::solveCurrentStep - " << timer1;

    if (result == -2) {
      opserr << "ModifiedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      if (factorOnce ==2) {
        factorOnce = 1;
      }
      return -3;
    }
    return result;
}

int
ModifiedNewton::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(4);
  data(0) = tangent;
  data(1) = iFactor;
  data(2) = cFactor;
  data(3) = factorOnce;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
ModifiedNewton::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  tangent = data(0);
  iFactor = data(1);
  cFactor = data(2);
  factorOnce = data(3);
  return 0;
}

void
ModifiedNewton::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
	s << "ModifiedNewton";
    }
}

int
ModifiedNewton::getNumIterations(void)
{
  return numIterations;
}
