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

// $Revision: 1.7 $
// $Date: 2008-04-01 00:35:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp,v $

// Written: fmk 
// Created: 02/06
                                                                        
// Description: This file contains the implementation of MumpsParallelSolver

#include <MumpsParallelSolver.h>
#include <MumpsParallelSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <ID.h>
#include <iostream>

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I)  cntl[(I)-1]  /* Ladruno ADR-75 P2: same for the real controls */

#include <mpi.h>

MumpsParallelSolver::MumpsParallelSolver(int ICNTL7, int ICNTL14,
					 int ICNTL35, double CNTL7)
  :LinearSOESolver(SOLVER_TAGS_MumpsParallelSolver),
   theMumpsSOE(0), rank(0), np(0), theComm(MPI_COMM_WORLD)  // Ladruno ADR43 P3a
{
  memset(&id, 0, sizeof(id));

  icntl7 = ICNTL7;
  icntl14 = ICNTL14;
  icntl35 = ICNTL35;   // Ladruno ADR-75 P2 (BLR; 0 = off = stock behavior)
  cntl7 = CNTL7;
  init = false;
  needsSetSize = false;
}

MumpsParallelSolver::MumpsParallelSolver(int mpi_comm, int ICNTL7, int ICNTL14,
					 int ICNTL35, double CNTL7)
  :LinearSOESolver(SOLVER_TAGS_MumpsParallelSolver),
   theMumpsSOE(0), rank(0), np(0), theComm(MPI_COMM_WORLD)  // Ladruno ADR43 P3a
{
  memset(&id, 0, sizeof(id));

  icntl14 = ICNTL14;
  icntl7 = ICNTL7;
  icntl35 = ICNTL35;   // Ladruno ADR-75 P2
  cntl7 = CNTL7;
  init = false;
  needsSetSize = false;
}

// Ladruno ADR43 P3a: point MUMPS at an arbitrary communicator. Any existing
// MUMPS instance is torn down so the next solve re-initializes on the new comm.
int
MumpsParallelSolver::setCommunicator(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL) {
    opserr << "WARNING MumpsParallelSolver::setCommunicator - MPI_COMM_NULL "
	   << "(this rank is not a member of the sub-communicator); ignored\n";
    return -1;
  }

  if (init == true) {
    id.job = -2;
    dmumps_c(&id); /* Terminate instance on the old comm */
    init = false;
  }

  theComm = comm;
  needsSetSize = true;

  // A factorization from the old comm is meaningless on the new one; without
  // this, a reconfigure-after-solve with no intervening zeroA() would issue a
  // solve-only JOB=3 on a re-analyzed-but-unfactored instance.
  if (theMumpsSOE != 0)
    theMumpsSOE->factored = false;

  return 0;
}


MumpsParallelSolver::~MumpsParallelSolver()
{
  id.job=-2; 
  if (init == true)
    dmumps_c(&id); /* Terminate instance */
}

int
MumpsParallelSolver::initializeMumps()
{
  if (needsSetSize == false)	{
    return 0;
  }
  else {

    if (init == true) {
      id.job=-2; 
      dmumps_c(&id); /* Terminate instance */
      init = false;
    } 
    
    id.job = -1;
      
    id.par = 1; // host involved in calcs
    id.sym = theMumpsSOE->matType;
    
    // Ladruno ADR43 P3a: honor a sub-communicator set via setCommunicator().
    // The historical _OPENMPI "0" magic value means "MUMPS's own COMM_WORLD",
    // so it is only valid when we are actually running on the world comm.
#ifdef _OPENMPI
    //    id.comm_fortran=-987654;
    if (theComm == MPI_COMM_WORLD)
      id.comm_fortran = 0;
    else
      id.comm_fortran = MPI_Comm_c2f(theComm);
#else
    id.comm_fortran = MPI_Comm_c2f(theComm);
#endif

    id.ICNTL(5) = 0; id.ICNTL(18) = 3;

    dmumps_c(&id);

    MPI_Comm_rank(theComm, &rank);   // Ladruno ADR43 P3a
    MPI_Comm_size(theComm, &np);     // Ladruno ADR43 P3a
    
    init = true;
    
    // parallel solver; distributed i/p matrix A
    id.ICNTL(5) = 0; id.ICNTL(18) = 3;
    
    // No outputs 
    //  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0; 
    id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 3;
    //id.ICNTL(1) = 1; id.ICNTL(2) = 1; id.ICNTL(3) = 1; id.ICNTL(4) = 3;
    
    id.ICNTL(14) = icntl14;
    id.ICNTL(7) = icntl7;
    // Ladruno ADR-75 P2: Block Low-Rank. MUST be set before the ANALYSIS
    // (job=1) as well as the factorization -- MUMPS uses it when planning the
    // assembly tree. icntl35=0 leaves the stock (full-rank) path untouched.
    id.ICNTL(35) = icntl35;
    if (cntl7 != 0.0)
      id.CNTL(7) = cntl7;

    int nnz = theMumpsSOE->nnz;
    int *colA = theMumpsSOE->colA;
    int *rowA = theMumpsSOE->rowA;
    
    // increment row and col A values by 1 for mumps fortran indexing
    for (int i = 0; i < nnz; i++) {
      rowA[i]++;
      colA[i]++;
    }
    
    // analyze the matrix
    id.n = theMumpsSOE->size;
    id.nz_loc = theMumpsSOE->nnz;
    id.irn_loc = theMumpsSOE->rowA;
    id.jcn_loc = theMumpsSOE->colA;
    id.a_loc = theMumpsSOE->A;
    
    // Call the MUMPS package to analyze the system
    id.job = 1;
    dmumps_c(&id);
    
    // decrement row and col A values by 1 to return to C++ indexing
    for (int i = 0; i < nnz; i++) {
      rowA[i]--;
      colA[i]--;
    }
    
    int info = id.infog[0];
    int info2   = id.infog[1];
    if (info != 0) {	
      opserr << "WARNING MumpsParallelSolver::setSize(void)- ";
      opserr << " Error " << info << " returned in substitution dmumps()\n";
      switch(info) {
      case -2:
	opserr << "nz " << info2 << " out of range\n";
	break;
      case -5:
	opserr << " out of memory allocation error\n";
	break;
      case -6:  
	opserr << " cause: Matrix is Singular in Structure: check your model\n";
	break;
      case -7:
	opserr << " out of memory allocation error\n";
	break;
      case -8:
	opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
	break;
      case -9:
	opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
	break;
      case -10:  
	opserr << " cause: Matrix is Singular Numerically\n";
	break;
      case -13:
	opserr << " out of memory wanted " << info2 << " (if < 0 mult absolute by 1 million)\n";
	break;
      default:
	opserr << " mumps returned infog[0] and infog[1] error codes: " << info << " and " << info2;
      }
    }

    if (info < 0)
      return info;
    
    needsSetSize = false;
    
    return info;
  }
}

int
MumpsParallelSolver::solveAfterInitialization(void)
{
  int n = theMumpsSOE->size;
  int nnz = theMumpsSOE->nnz;
  int *rowA = theMumpsSOE->rowA;
  int *colA = theMumpsSOE->colA;
  double *A = theMumpsSOE->A;
  double *X = theMumpsSOE->X;
  double *B = theMumpsSOE->B;

  // parallel solver; distributed i/p matrix A
  id.ICNTL(5)=0; id.ICNTL(18)=3; 

  // No outputs 
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  //id.ICNTL(1) = 1; id.ICNTL(2) = 1; id.ICNTL(3) = 1; id.ICNTL(4) = 3;


  id.ICNTL(14)=icntl14;
  id.ICNTL(7)=icntl7;
  id.ICNTL(35)=icntl35;                 // Ladruno ADR-75 P2 (BLR)
  if (cntl7 != 0.0)
    id.CNTL(7) = cntl7;

  // increment row and col A values by 1 for mumps fortran indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]++;
    colA[i]++;
  }
  
  if (rank == 0) {
    id.n   = n; 
    for (int i=0; i<n; i++) {
      X[i] = B[i];
    }
    id.rhs = X;
  } 

  // factor the matrix
  id.nz_loc  = nnz; 
  id.irn_loc = rowA;
  id.jcn_loc = colA;
  id.a_loc   = A; 

  if (theMumpsSOE->factored == false) {

    // Call the MUMPS package to factor & solve the system
    id.job = 5;
    dmumps_c(&id);
    theMumpsSOE->factored = true;

  } else {

    // Call the MUMPS package to solve the system
    id.job = 3;
    dmumps_c(&id);
  }	

  int info = id.infog[0];
  int info2   = id.infog[1];
  if (info != 0) {	
    opserr << "WARNING MumpsParallelSolver::solve(void)- ";
    opserr << " Error " << info << " returned in substitution dmumps()\n";
    switch(info) {
    case -2:
      opserr << "nz " << info2 << " out of range\n";
      break;
    case -5:
      opserr << " out of memory allocation error\n";
      break;
    case -6:  
      opserr << " cause: Matrix is Singular in Structure: check your model\n";
      break;
    case -7:
      opserr << " out of memory allocation error\n";
      break;
    case -8:
      opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
      break;
    case -9:
      opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
      break;
    case -10:  
      opserr << " cause: Matrix is Singular Numerically\n";
      break;
    case -13:
      opserr << " out of memory wanted " << info2 << " (if < 0 mult absolute by 1 million)\n";
      break;
    default:
      opserr << " mumps returned infog[0] and infog[1] error codes: " << info << " and " << info2;
    }
    return info;
  }

  // decrement row and col A values by 1 to return to C++ indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]--;
    colA[i]--;
  }

  return 0;
}

int
MumpsParallelSolver::solve(void)
{
	int initializationResult = initializeMumps();

	if (initializationResult == 0)
		return solveAfterInitialization();
	else
		return initializationResult;
}

int
MumpsParallelSolver::setSize(void)
{
  needsSetSize = true;
  return 0;
}

int
MumpsParallelSolver::sendSelf(int cTag, Channel &theChannel)
{
  // Ladruno ADR-75 P2: icntl35 (BLR mode) rides along, else the subordinate
  // ranks would factor FULL-RANK while P0 factors BLR -- an inconsistent
  // distributed factorization. cntl7 (a double) goes in its own Vector.
  ID icntlData(3);

  icntlData(0) = icntl7;
  icntlData(1) = icntl14;
  icntlData(2) = icntl35;
  theChannel.sendID(0, cTag, icntlData);

  Vector cntlData(1);
  cntlData(0) = cntl7;
  theChannel.sendVector(0, cTag, cntlData);

  return 0;
}

int
MumpsParallelSolver::recvSelf(int ctag,
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  // Ladruno ADR-75 P2: must mirror sendSelf exactly (3 ints + 1 double).
  ID icntlData(3);

  theChannel.recvID(0, ctag, icntlData);

  icntl7 = icntlData(0);
  icntl14 = icntlData(1);
  icntl35 = icntlData(2);

  Vector cntlData(1);
  theChannel.recvVector(0, ctag, cntlData);
  cntl7 = cntlData(0);

  return 0;
}

int 
MumpsParallelSolver::setLinearSOE(MumpsParallelSOE &theSOE)
{
  theMumpsSOE = &theSOE;
  return 0;
}













