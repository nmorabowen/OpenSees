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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/LoadControl.h,v $
                                                                        
                                                                        
#ifndef LoadControl_h
#define LoadControl_h

// File: ~/analysis/integrator/LoadControl.h
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for LoadControl.
// LoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
//
// What: "@(#) LoadControl.h, revA"

#include <StaticIntegrator.h>
#include <classTags.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class EquiSolnAlgo;
class ReliabilityDomain;

class LoadControl : public StaticIntegrator
{
  public:
    LoadControl(double deltaLambda, int numIncr, 
		double minLambda, double maxlambda, int classtag=INTEGRATOR_TAGS_LoadControl);

    ~LoadControl();

    int newStep(void);    
    int update(const Vector &deltaU);
    int setDeltaLambda(double newDeltaLambda);
    // Ladruno WP-232f: the MPM quasi-static driver (`mpm analyze` in mode A) sets
    // the step itself. setDeltaLambda alone is not enough: newStep clamps the
    // increment to [dLambdaMin, dLambdaMax], which default to dLambda. The driver
    // reads the user's (dLambda, min, max), steps with min = max = its dt, and
    // puts the user's values back at the end.
    void ladrunoGetIncrement(double &d, double &dMin, double &dMax) const   // Ladruno WP-232f
      { d = deltaLambda; dMin = dLambdaMin; dMax = dLambdaMax; }
    void ladrunoSetIncrement(double d, double dMin, double dMax)            // Ladruno WP-232f
      { deltaLambda = d; dLambdaMin = dMin; dLambdaMax = dMax; numIncrLastStep = specNumIncrStep; }

    // Public methods for Output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    int formEleResidual(FE_Element *theEle);
    
    // Adding sensitivity
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity(const Vector &v, int gradNum, int numGrads);
    int commitSensitivity(int gradNum, int numGrads);
    int computeSensitivities(void);//Abbas
    bool computeSensitivityAtEachIteration();


    ///////////////////////
    
protected:
    
  private:
    double deltaLambda;  // dlambda at step (i-1)
    
    double specNumIncrStep, numIncrLastStep; // Jd & J(i-1) 
    double dLambdaMin, dLambdaMax; // min & max values for dlambda at step (i)
    
    // Adding sensitivity
    int gradNumber;
    int sensitivityFlag;
   // EquiSolnAlgo *theAlgorithm;
    ReliabilityDomain *theDomain;
    ////////////////////
};

#endif

