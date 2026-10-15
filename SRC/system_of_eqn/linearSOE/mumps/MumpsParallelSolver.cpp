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
#include <LadrunoMumpsMixed.h>   // Ladruno WP-194 (empty unless _LADRUNO_SMUMPS)
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <ID.h>
#include <iostream>

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I)  cntl[(I)-1]  /* Ladruno ADR-75 P2: same for the real controls */

#include <mpi.h>

MumpsParallelSolver::MumpsParallelSolver(int ICNTL7, int ICNTL14,
					 int ICNTL35, double CNTL7,
					 int PrintStats)
  :LinearSOESolver(SOLVER_TAGS_MumpsParallelSolver),
   theMumpsSOE(0), rank(0), np(0), theComm(MPI_COMM_WORLD)  // Ladruno ADR43 P3a
{
  memset(&id, 0, sizeof(id));

  icntl7 = ICNTL7;
  icntl14 = ICNTL14;
  icntl35 = ICNTL35;   // Ladruno ADR-75 P2 (BLR; 0 = off = stock behavior)
  cntl7 = CNTL7;
  printStats = PrintStats;   // Ladruno ADR-75 P2b
  tAnalysis = -1.0;          // WP-190
  tSolveSum = 0.0;
  nSolveSum = 0;
  init = false;
  needsSetSize = false;
  precision = 0;       // Ladruno WP-194: double unless -precision single
  refineSteps = 0;
  theMixed = 0;
}

MumpsParallelSolver::MumpsParallelSolver(int mpi_comm, int ICNTL7, int ICNTL14,
					 int ICNTL35, double CNTL7,
					 int PrintStats)
  :LinearSOESolver(SOLVER_TAGS_MumpsParallelSolver),
   theMumpsSOE(0), rank(0), np(0), theComm(MPI_COMM_WORLD)  // Ladruno ADR43 P3a
{
  memset(&id, 0, sizeof(id));

  icntl14 = ICNTL14;
  icntl7 = ICNTL7;
  icntl35 = ICNTL35;   // Ladruno ADR-75 P2
  cntl7 = CNTL7;
  printStats = PrintStats;   // Ladruno ADR-75 P2b
  init = false;
  needsSetSize = false;
  precision = 0;       // Ladruno WP-194: double unless -precision single
  refineSteps = 0;
  theMixed = 0;
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
#ifdef _LADRUNO_SMUMPS
  if (theMixed != 0)
    theMixed->terminate();   // Ladruno WP-194: the float instance too
#endif

  theComm = comm;
  needsSetSize = true;

  // A factorization from the old comm is meaningless on the new one; without
  // this, a reconfigure-after-solve with no intervening zeroA() would issue a
  // solve-only JOB=3 on a re-analyzed-but-unfactored instance.
  if (theMumpsSOE != 0)
    theMumpsSOE->factored = false;

  return 0;
}


// WP-190: solves that ran after the last stats print (the final factor's
// lifetime) have no later block to ride in; flush them once at teardown.
void
MumpsParallelSolver::reportSolveTail(void)
{
  if (printStats == 0 || init == false || nSolveSum == 0)
    return;
  // Not collective on purpose: teardown order across ranks is not guaranteed,
  // so this reports rank 0's own time (no max over ranks).
  if (rank == 0)
    opserr << "MUMPS solve tail: t_solve s (rank 0, n) = " << tSolveSum << ", "
	   << nSolveSum << "\n";
  tSolveSum = 0.0;
  nSolveSum = 0;
}

MumpsParallelSolver::~MumpsParallelSolver()
{
  // MPI may already be finalized when a static solver is destroyed.
  int mpiFinal = 0;
  MPI_Finalized(&mpiFinal);
  if (!mpiFinal)
    reportSolveTail();
  id.job=-2; 
  if (init == true)
    dmumps_c(&id); /* Terminate instance */
#ifdef _LADRUNO_SMUMPS
  delete theMixed;   // Ladruno WP-194: terminates the SMUMPS instance
#endif
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
    double tA0 = (printStats != 0) ? MPI_Wtime() : 0.0;   // WP-190
    dmumps_c(&id);
    if (printStats != 0)
      tAnalysis = MPI_Wtime() - tA0;
    
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
    // WP-190: with stats on, JOB=5 is issued as its two halves (JOB=2 then
    // JOB=3, which is what MUMPS defines JOB=5 to be) so each phase can be
    // timed. With stats off the single JOB=5 call is unchanged.
    double tFactor = 0.0;
    if (printStats != 0) {
      id.job = 2;
      double t0 = MPI_Wtime();
      dmumps_c(&id);
      tFactor = MPI_Wtime() - t0;
    } else {
      id.job = 5;
      dmumps_c(&id);
    }
    theMumpsSOE->factored = true;

    // Ladruno ADR-75 P2b: `-stats`. BLR's entire justification is FACTOR
    // MEMORY, which is otherwise completely invisible from OpenSees -- with no
    // way to tell whether compression did anything on a given model. MUMPS
    // convention: a NEGATIVE INFOG entry means -value * 10^6.
    if (printStats != 0) {
      // WP-190: collective max over ranks; outside the rank-0 guard.
      double tmax[3] = { tAnalysis, tFactor, tSolveSum };
      MPI_Allreduce(MPI_IN_PLACE, tmax, 3, MPI_DOUBLE, MPI_MAX, theComm);
      int nSolves = nSolveSum;
      bool hadAnalysis = (tmax[0] >= 0.0);
      tAnalysis = -1.0;
      tSolveSum = 0.0;
      nSolveSum = 0;
      if (rank == 0) {
      double fEntries = (double)id.infog[8];              // INFOG(9)
      if (fEntries < 0) fEntries = -fEntries * 1.0e6;
      opserr << "MUMPS stats: n=" << id.n
	     << " BLR(ICNTL35)=" << icntl35;
      if (icntl35 != 0)
	opserr << " eps(CNTL7)=" << cntl7;
      opserr << "\n";
      opserr << "  factor entries INFOG(9)   = " << fEntries << "\n";
      opserr << "  factor MB/proc INFOG(21)  = " << id.infog[20] << "  (max over procs)\n";
      opserr << "  factor MB tot  INFOG(22)  = " << id.infog[21] << "  (sum over procs)\n";
      opserr << "  elim flops     RINFOG(3)  = " << id.rinfog[2] << "\n";
      // MUMPS reports the BLR-compressed vs full-rank factor/flop comparison in
      // RINFOG(14)/RINFOG(15) when ICNTL(35) is active; printed raw and clearly
      // labelled so the numbers are checkable against the MUMPS guide rather
      // than silently reinterpreted here.
      if (icntl35 != 0) {
	opserr << "  BLR RINFOG(14)            = " << id.rinfog[13] << "\n";
	opserr << "  BLR RINFOG(15)            = " << id.rinfog[14] << "\n";
      }
      // WP-190 wall-clock phase timings (max over ranks). Analysis only when
      // it ran since the previous print; solves = JOB=3 calls since the
      // previous print. A reused (skipped) factorization prints nothing.
      if (hadAnalysis)
	opserr << "  t_analysis s (max)  = " << tmax[0] << "\n";
      opserr << "  t_factor s (max)    = " << tmax[1] << "\n";
      opserr << "  t_solve s (max, n)  = " << tmax[2] << ", " << nSolves << "\n";
      }
    }

    // WP-190: second half of the split JOB=5 (stats on only).
    if (printStats != 0 && id.infog[0] >= 0) {
      id.job = 3;
      double t0 = MPI_Wtime();
      dmumps_c(&id);
      tSolveSum += MPI_Wtime() - t0;
      nSolveSum++;
    }

  } else {

    // Call the MUMPS package to solve the system
    id.job = 3;
    double t0 = (printStats != 0) ? MPI_Wtime() : 0.0;   // WP-190
    dmumps_c(&id);
    if (printStats != 0) {
      tSolveSum += MPI_Wtime() - t0;
      nSolveSum++;
    }
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
	// Ladruno WP-194: -precision single takes its own path; the DMUMPS path
	// below is untouched (and never initialized) when it is on.
	if (precision != 0)
		return solveMixed();

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
  ID icntlData(4);   // WP-190: +printStats (the timing reductions are collective)

  icntlData(0) = icntl7;
  icntlData(1) = icntl14;
  icntlData(2) = icntl35;
  icntlData(3) = printStats;
  theChannel.sendID(0, cTag, icntlData);

  Vector cntlData(1);
  cntlData(0) = cntl7;
  theChannel.sendVector(0, cTag, cntlData);

  // Ladruno WP-194: precision + refine ride in their own ID, after the
  // controls, so the SP subordinate ranks factor in the same precision as P0
  // (a DMUMPS/SMUMPS instance pair across ranks would deadlock or mis-solve).
  ID precData(2);
  precData(0) = precision;
  precData(1) = refineSteps;
  theChannel.sendID(0, cTag, precData);

  return 0;
}

int
MumpsParallelSolver::recvSelf(int ctag,
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  // Ladruno ADR-75 P2: must mirror sendSelf exactly (3 ints + 1 double).
  ID icntlData(4);   // WP-190

  theChannel.recvID(0, ctag, icntlData);

  icntl7 = icntlData(0);
  icntl14 = icntlData(1);
  icntl35 = icntlData(2);
  printStats = icntlData(3);

  Vector cntlData(1);
  theChannel.recvVector(0, ctag, cntlData);
  cntl7 = cntlData(0);

  // Ladruno WP-194: mirror of sendSelf
  ID precData(2);
  theChannel.recvID(0, ctag, precData);
  if (this->setPrecision(precData(0), precData(1)) < 0)
    return -1;

  return 0;
}

int 
MumpsParallelSolver::setLinearSOE(MumpsParallelSOE &theSOE)
{
  theMumpsSOE = &theSOE;
  return 0;
}

// ---------------------------------------------------------------------------
// Ladruno WP-194: mixed precision (LadrunoMumpsMixed.h for the design).

bool
MumpsParallelSolver::singlePrecisionAvailable(void)
{
#ifdef _LADRUNO_SMUMPS
  return true;
#else
  return false;
#endif
}

int
MumpsParallelSolver::setPrecision(int single, int refine)
{
  if (refine < 0)
    refine = 0;
  if (single == 0) {
#ifdef _LADRUNO_SMUMPS
    delete theMixed;
    theMixed = 0;
#endif
    precision = 0;
    refineSteps = 0;
    return 0;
  }
#ifdef _LADRUNO_SMUMPS
  delete theMixed;
  theMixed = new LadrunoMumpsMixed(refine);
  precision = 1;
  refineSteps = refine;
  needsSetSize = true;
  if (theMumpsSOE != 0)
    theMumpsSOE->factored = false;
  return 0;
#else
  opserr << "WARNING MumpsParallelSolver: -precision single needs single-precision MUMPS "
	 << "(libsmumps + smumps_c.h), which this build does not have\n";
  return -1;
#endif
}

int
MumpsParallelSolver::solveMixed(void)
{
#ifdef _LADRUNO_SMUMPS
  if (theMixed == 0 || theMumpsSOE == 0)
    return -1;
  if (needsSetSize) {
    theMixed->requestAnalysis();
    needsSetSize = false;
  }
  LadrunoMumpsControls c;
  c.icntl7 = icntl7;
  c.icntl14 = icntl14;
  c.icntl35 = icntl35;
  c.cntl7 = cntl7;
  c.printStats = printStats;
  bool factored = theMumpsSOE->factored;
  int res = theMixed->solve(theComm, theMumpsSOE->matType, theMumpsSOE->size,
			    theMumpsSOE->nnz, theMumpsSOE->rowA, theMumpsSOE->colA,
			    theMumpsSOE->A, theMumpsSOE->B, theMumpsSOE->X,
			    factored, c);
  theMumpsSOE->factored = factored;
  if (res < 0) {
    // the next solve analyses and factors from scratch
    theMixed->requestAnalysis();
    theMumpsSOE->factored = false;
  }
  return res;
#else
  return -1;
#endif
}













