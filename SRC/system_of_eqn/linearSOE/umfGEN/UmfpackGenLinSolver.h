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
// $Date: 2009-05-11 20:56:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.h, revA"

#ifndef UmfpackGenLinSolver_h
#define UmfpackGenLinSolver_h

#include <LinearSOESolver.h>
#include <vector>
#include "../../../../OTHER/UMFPACK/umfpack.h"

class UmfpackGenLinSOE;

class UmfpackGenLinSolver : public LinearSOESolver
{
  public:
    // Ladruno (ADR-40 rank 2): strategy + pivot tolerance are configurable.
    // Default strategy is UMFPACK_STRATEGY_AUTO (was hardcoded SYMMETRIC, which
    // mis-orders unsymmetric tangents); default pivotTol keeps the legacy 1.0.
    // Ladruno (ADR-40 rank 8/10): persist=true (default) reuses the Numeric LU
    // across solves while theSOE->factored stays true; -noNumericPersist sets it
    // false to restore the legacy free-every-solve behavior (bit-for-bit).
    UmfpackGenLinSolver(bool useLongIndices = false,
                        int strategy = UMFPACK_STRATEGY_AUTO,
                        double pivotTol = 1.0,
                        bool persist = true);
    ~UmfpackGenLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(UmfpackGenLinSOE &theSOE);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

  protected:

  private:
    void syncIndexBuffers(void);

    void freeNumeric(void);   // Ladruno (ADR-40 rank 8/10)

    bool useLongIndices;
    int strategy;        // Ladruno (ADR-40 rank 2)
    double pivotTol;     // Ladruno (ADR-40 rank 2)
    bool persist;        // Ladruno (ADR-40 rank 8/10): reuse Numeric across solves
    void *Symbolic;
    void *Numeric;       // Ladruno (ADR-40 rank 8/10): persisted LU factorization
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    UmfpackGenLinSOE *theSOE;
    std::vector<SuiteSparse_long> Ap64;
    std::vector<SuiteSparse_long> Ai64;
};

#endif
