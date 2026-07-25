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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-04-15 07:15:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h,v $
                                                                        
                                                                        
#ifndef MumpsParallelSolver_h
#define MumpsParallelSolver_h

// Written: fmk 
// Created: 02/06
//
// Description: This file contains the class definition for Mumps.
// A Mumps object can be constructed to solve a MumpsSOE. It obtains 
// the solution by making calls on the the Mumps library, which is based 
// on public domain software developed during the Esprit IV European project 
// PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. Since this first 
// public domain version in 1999, the developments are supported by CERFACS, 
// ENSEEIHT-IRIT, and INRIA Rhone-Alpes. 
// Main contributors to MUMPS are Patrick Amestoy, Iain Duff, Abdou Guermouche,
// Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

// Up-to-date copies of the MUMPS package can be obtained from the Web pages 
// http://www.enseeiht.fr/apo/MUMPS/ or http://graal.ens-lyon.fr/MUMPS

// What: "@(#) Mumps.h, revA"

#include <mpi.h>

#include <LinearSOESolver.h>
extern "C" {
#include <dmumps_c.h>
}

class MumpsParallelSOE;

class MumpsParallelSolver : public LinearSOESolver
{
  public:
  // Ladruno ADR-75 P2: ICNTL35/CNTL7 are the Block Low-Rank controls (default
  // OFF => byte-identical to the pre-BLR solver).
  MumpsParallelSolver(int ICNTL7 = 7, int ICNTL14 = 20,
		      int ICNTL35 = 0, double CNTL7 = 0.0,
		      int printStats = 0);

  MumpsParallelSolver(int MPI_COMM,
		      int ICNTL7,
		      int ICNTL14,
		      int ICNTL35 = 0,
		      double CNTL7 = 0.0,
		      int printStats = 0);

  virtual ~MumpsParallelSolver();
  
  int solve(void);
  int setSize(void);
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  int setLinearSOE(MumpsParallelSOE &theSOE);

  // Ladruno ADR43 P3a: run MUMPS on an arbitrary communicator (e.g. an
  // MPI_Comm_split sub-communicator) instead of the hardwired MPI_COMM_WORLD.
  // Must be called before the first solve; forces re-initialization.
  int setCommunicator(MPI_Comm comm);

 protected:

 private:

  int initializeMumps(void);
  int solveAfterInitialization(void);

  bool init;
  MumpsParallelSOE *theMumpsSOE;

  int rank;
  int np;
  int icntl14;
  int icntl7;
  // Ladruno ADR-75 P2: Block Low-Rank. ICNTL(35)=0 off / 1 BLR facto+solve /
  // 2 BLR facto, full solve / 3 BLR facto, no compression of the solve.
  // CNTL(7) is the dropping tolerance eps (0.0 => MUMPS default). BLR is an
  // APPROXIMATE factorization: it trades accuracy for factor memory/flops, so
  // it is opt-in and must never be enabled on a byte-identical/oracle lane.
  int icntl35;
  double cntl7;
  // Ladruno ADR-75 P2b: `-stats` -- dump MUMPS INFOG/RINFOG after each
  // factorization. Without this, BLR's whole purpose (factor MEMORY) is
  // invisible: a user cannot tell whether compression did anything.
  int printStats;

  MPI_Comm theComm;  // Ladruno ADR43 P3a: defaults to MPI_COMM_WORLD

  DMUMPS_STRUC_C id;

  bool needsSetSize;
};

#endif

