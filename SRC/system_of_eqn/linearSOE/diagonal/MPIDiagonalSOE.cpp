/* 
 * @author: George petropoulos <gnp>
 *
 * @Description: Sets up a parallel distributed diagonal SOE
 *               all processes are equally used for computation.
 *               Model needs to be appropriatelly distributed.
 *               Responsibility of builder or author of main to account for this
 *               SOE laready assumes appropriate distribution of graph
 *               
 * @Date: 12/05
 *
 * Copyright: ALL RIGHTS RESERVED BY AUTHOR
 *
 */

#include <MPIDiagonalSOE.h>
#include <MPIDiagonalSolver.h>
#include <algorithm>   // Ladruno (ADR-74 setSize fix): std::sort
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

// Ladruno (ADR-74 setSize investigation): dc.s.* sub-brackets so the profiler can
// attribute the ~N^1.9 first-step setSize cost to its internal phases.
#include <profiler/ProfilerMacros.h>



MPIDiagonalSOE::MPIDiagonalSOE(MPIDiagonalSolver &the_Solver, bool lumped)
  :LinearSOE(the_Solver, LinSOE_TAGS_MPIDiagonalSOE),
   size(0), A(0), B(0), X(0), sharedA(0), sharedB(0), maxSharedA(0), maxSharedB(0), isAfactored(false), updateA(true),
   vectX(0), vectB(0), dataShared(0),
   actualNeighbors(0),maxNeighbors(0),myNeighbors(0),myNeighborsSizes(0), 
   myDOFsArray(0), myDOFsSharedArray(0),maxDOFsSharedArray(0),posLocKey(0),
   processID(0), numProcesses(0),
   numChannels(0), theChannels(0), localCol(0),
   myDOFs(0,32), myDOFsShared(0,16), numShared(0), theModel(0), lumpDiagonal(lumped)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  the_Solver.setLinearSOE(*this);
  // cached neighbors
  std::map<int, int*> sharedDOFsMap;
  this->myActualNeighborsSharedDOFs = sharedDOFsMap;
  std::map<int, double*> sharedBsMap;
  this->myActualNeighborsSharedBs = sharedBsMap;
  std::map<int, double*> sharedBsToSendMap;
  this->myActualNeighborsBsToSend = sharedBsToSendMap;

}


MPIDiagonalSOE::~MPIDiagonalSOE()
{
  if (A != 0) delete [] A;
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;
  if (sharedA !=0) delete [] sharedA;
  if (sharedB !=0) delete [] sharedB;
  if (maxSharedA !=0) delete [] maxSharedA;
  if (maxSharedB !=0) delete [] maxSharedB;
  if (vectX !=0) delete vectX;
  if (vectB !=0) delete vectB;
  if (myNeighbors !=0) delete [] myNeighbors;
  if (myNeighborsSizes !=0) delete [] myNeighborsSizes;
  if (dataShared !=0) delete [] dataShared;
  if (myDOFsSharedArray !=0) delete [] myDOFsSharedArray;
  if (myDOFsArray !=0) delete [] myDOFsArray;
  if (maxDOFsSharedArray !=0) delete [] maxDOFsSharedArray;
  if (posLocKey !=0) delete [] posLocKey;


  /// add clean the map contents prior to delete maps
  /// add delete the maps
  if (theChannels != 0)
    delete [] theChannels;
}


int
MPIDiagonalSOE::getNumEqn(void) const
{
  return size;
}


// Ladruno (consistent/Olovsson parallel PCG) ===============================
// Sum vector v's shared-DOF entries across neighbour ranks, in place, so every
// sharing rank ends with the global sum at its shared DOFs. Mirrors the RHS (B)
// exchange in MPIDiagonalSolver::solve() (the "else" / notSet==false branch),
// but on an arbitrary vector. REQUIRES the neighbour maps already built by the
// first solve() (guaranteed: refineAccel runs after the step's diagonal solve);
// a no-op if they are not yet set up. v.Size() must equal the SOE size.
int
MPIDiagonalSOE::assembleSharedSum(Vector &v)
{
  if (v.Size() != size) return -1;
  // The send buffers (myActualNeighborsBsToSend) AND the per-neighbour posloc sizes
  // (myNeighborsSizes) are populated only by the FIRST solve()'s setup pass; a no-op
  // before that (refineAccel always runs post-solve, so this just guards np=1 / a
  // pre-solve call). myActualNeighborsSharedDOFs is allocated earlier in setSize, so
  // it is NOT the right readiness signal.
  if (myActualNeighborsBsToSend.empty()) {
    // Tripwire (adversarial-review hardening): empty buffers with actual neighbours means
    // this was called BEFORE the first solve() set up the exchange — the cross-rank sum
    // would be SILENTLY dropped. Harmless for the shipped consistent integrators (their
    // refineAccel always runs after a solve), but warn once so a future pre-solve caller
    // is not bitten by a silent wrong answer. No neighbours (np=1 / interior) is a genuine
    // no-op and stays silent.
    static bool warnedNotReady = false;
    if (actualNeighbors > 0 && !warnedNotReady) {
      warnedNotReady = true;
      opserr << "WARNING MPIDiagonalSOE::assembleSharedSum - called before the first solve() "
                "set up the neighbour exchange; shared-DOF entries NOT reduced this call.\n";
    }
    return 0;
  }

  double *vv = &v(0);
  int n = (2*actualNeighbors > 0) ? 2*actualNeighbors : 1;
  MPI_Request *theRequests = new MPI_Request[n];
  MPI_Status  *theStatuses = new MPI_Status[n];
  int ct = 0;

  for (int i = 0; i < maxNeighbors; i++)
    if ((i != processID) && (myNeighbors[i] == 1))
      MPI_Irecv((myActualNeighborsSharedBs.find(i)->second), maxShared, MPI_DOUBLE,
                i, MPI_ANY_TAG, MPI_COMM_WORLD, &theRequests[ct++]);

  for (int i = 0; i < maxNeighbors; i++) {
    if ((i != processID) && (myNeighbors[i] == 1)) {
      int upto = myNeighborsSizes[i];
      double *tmpptr = (myActualNeighborsBsToSend.find(i)->second);
      int *posloc = (myActualNeighborsSharedDOFs.find(i)->second);
      for (int jj = 0; jj < upto; jj++)
        tmpptr[jj] = vv[posloc[jj]];
      MPI_Isend(tmpptr, upto, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &theRequests[ct++]);
    }
  }

  MPI_Waitall(ct, theRequests, theStatuses);

  for (int i = 0; i < maxNeighbors; i++) {
    if ((i != processID) && (myNeighbors[i] == 1)) {
      int *posloc = (myActualNeighborsSharedDOFs.find(i)->second);
      double *dat = (myActualNeighborsSharedBs.find(i)->second);
      for (int k = 0; k < myNeighborsSizes[i]; k++)
        vv[posloc[k]] += dat[k];
    }
  }

  delete [] theRequests;
  delete [] theStatuses;
  return 0;
}

// Sum a scalar across all ranks (global inner-product reduction for the PCG).
double
MPIDiagonalSOE::globalReduceSum(double localVal) const
{
  double g = 0.0;
  MPI_Allreduce(&localVal, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return g;
}
// Ladruno end ==============================================================



//
int 
MPIDiagonalSOE::setSize(Graph &theGraph)
{
  /// debug  ////////////////////////////////////////////////////////////////////////////////////////////
  double ts=0.0;
  double te=0.0;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  ts = MPI_Wtime();
  ////////////////////////////////////////////////////////////////////////////////////////////
  
  // Initialize some MPI stuff and some variables
  MPI_Comm_size(MPI_COMM_WORLD,&numProcesses);
  MPI_Status status;
  MPI_Request Srequests[2];
  MPI_Request Rrequests[2];
  int p0Size = 0;
  int result = 0;
  maxNeighbors = numProcesses;
  myNeighbors = new int[maxNeighbors];
  myNeighborsSizes = new int[maxNeighbors];
  for (int i=0; i<maxNeighbors; i++)
    myNeighbors[i] =-1;
  actualNeighbors = 0;
  
  // get the size of the local system
  size = theGraph.getNumVertex();

  if (size == 0) {
    opserr << "WARNING MPIDiagonalSOE::setSize - ZERO size\n";
    exit(-1);
  } 

  int* sharedDOFs = new int[size];
  for (int i=0; i<size; i++)
    sharedDOFs[i] =0;
  
  myDOFsArray = new int[size];

  // get the actual dof;s of the system
  int count = 0;
  { OPS_PROFILE_SCOPE("dc.s.fill");   // Ladruno (ADR-74): graph-vertex walk
  Vertex *theVertex;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    int vertexTag = theVertex->getTag();
    myDOFsArray[count++] = vertexTag;
  }
  }
  static ID otherSize(1);
//  delete &theGraph;


  myDOFs.setData(myDOFsArray, size);

  // sort the dof ID
  { OPS_PROFILE_SCOPE("dc.s.sort");
  // Ladruno (ADR-74 setSize fix): std::sort replaces the hand-rolled leftmost-pivot
  // quicksort (q_sort below, kept for provenance). The DOF graph hands us its
  // std::map-backed vertices in ASCENDING tag order, so this array arrives already
  // sorted — the old q_sort's deterministic O(n^2/2) worst case (measured ~N^1.9:
  // 10.1 / 101.4 / 501.1 s at the 0.25/1/2 M-hex np8 rungs). std::sort produces the
  // identical ascending array (tags unique) at O(n log n) worst case, and also
  // removes q_sort's O(n) recursion depth on sorted input.
  std::sort(myDOFsArray, myDOFsArray + size);
  }

  int* sendSize = new int[1];
  sendSize[0] = size;
  int* recvSize = new int[1];
  recvSize[0] = 0;
  int* tmpsendSize = new int[1];
  *tmpsendSize = size;

  // Broadcast to all other processes size and DOF ID

  OPS_PROFILE_SCOPE("dc.s.rest");   // Ladruno (ADR-74): everything after the sort, refined by the nested scopes below

  int* allSizes = new int[numProcesses];
  for (int i=0; i<numProcesses; i++)
    allSizes[i] = 0;
  int* max = new int[1];             // Ladruno (ADR-74): hoisted above the dc.s.bcast scope (deleted at function end)
  max[0] =0;
  { OPS_PROFILE_SCOPE("dc.s.bcast"); // Ladruno (ADR-74): np-round bcast + two-pointer intersections
  MPI_Allgather(tmpsendSize,1,MPI_INT,allSizes,1,MPI_INT,MPI_COMM_WORLD);
  for (int i=0; i<numProcesses; i++)
    if ( max[0] < allSizes[i] )
      max[0] = allSizes[i];
  
  if (max[0] < size ) {
    opserr << " MPIDiagonalSOE::setSize() L170 : SEVERE ERROR : ABORTING TASK \n";
    exit(-1);
  }


  int* piShared= new int[1];
  int* piData = new int[max[0]];
  ID*  piDOFs = new ID();
  piShared[0] =0;
  
  for (int i=0; i< numProcesses; i++) {
    for (int j=0;j<max[0]; j++)
      piData[j] = 0;
    for (int j=0; j<size; j++)
      piData[j] = myDOFsArray[j];
    MPI_Bcast(piData,allSizes[i],MPI_INT,i,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (i!= processID) {
      piDOFs->setData(piData,allSizes[i]);
      intersections(myDOFs, *piDOFs, size, allSizes[i], piShared, sharedDOFs);
      if ( (*piShared)!=0 & (actualNeighbors < maxNeighbors) ) { 
	myNeighbors[i] = 1;
      }
    }
  }
  
  delete [] piData;
  delete piDOFs;
  delete [] piShared;
  delete[] allSizes;
  }                                  // Ladruno (ADR-74): end dc.s.bcast

  // Receive from all other processes their size and their ID's
  // perform intersection operation and figure out
  // with which processes does myProcessID have to communicate
  // to avoid blocking we call Irecv & wait, to make up for the send's

  
  //figure out total # of shared dofs
  // and allocate appropriate memory
  { OPS_PROFILE_SCOPE("dc.s.shared"); // Ladruno (ADR-74): shared list + allgather + allocations + zero-init
  numShared =0;
  for (int i=0; i<size; i++) 
    if (sharedDOFs[i] == 1)
      numShared++;


  myDOFsSharedArray = new int[numShared];
  
  //copy of all the shared dofs of the current process
  int pos=0;
  for (int i=0; i<size; i++) {
    if (sharedDOFs[i] == 1) {
      myDOFsSharedArray[pos++] = myDOFsArray[i];
    }
  }
  
  myDOFsShared.setData(myDOFsSharedArray, numShared);
  

  //  opserr << " pid " << processID << " has shared dofs " << myDOFsShared << endln;
  // broadcast shared size to all of the processes
  // that myProcessID shared dofs with
  *tmpsendSize = numShared;

  MPI_Allgather(tmpsendSize,1,MPI_INT,myNeighborsSizes,1,MPI_INT,MPI_COMM_WORLD);
  
  this->maxShared =0;
  for(int i=0; i<numProcesses; i++)
    if (maxShared < myNeighborsSizes[i])
      maxShared = myNeighborsSizes[i];

  if (maxShared < numShared ) {
    opserr << " MPIDiagonalSOE::setSize() L247 : SEVERE ERROR : ABORTING TASK \n";
    exit(-1);
  }
  
 // fprintf(stderr, "MAX-SHARED  %d \n", maxShared);
  
  //  cached neighbors 
  actualNeighbors =0;
  //  Now that i know all my neighbors and our shared size maximum
  //   i shall populate my maps of shared dofs so i can cache them
  //   i allocate space for dof id's and dof's
  for (int i=0; i<maxNeighbors; i++) {
    if (myNeighbors[i] == 1) {
      //opserr << " pid : " << processID << " neighbors w/ pid : " << i << endln;
      myActualNeighborsSharedDOFs[i] = new int[maxShared];
      myActualNeighborsSharedBs[i] = new double[maxShared];
      actualNeighbors++;
    }
  }

  // allocations done;
  // cached neighbors ended

  maxDOFsSharedArray = new int[maxShared];
  posLocKey = new int[size];

  // clear up memory and allocate new one
  if (A != 0) delete [] A; A = 0;
  if (B != 0) delete [] B; B = 0;
  if (X != 0) delete [] X; X = 0;
  if (sharedA != 0) delete [] sharedA; sharedA = 0;
  if (sharedB != 0) delete [] sharedB; sharedB = 0;
  if (maxSharedA != 0) delete [] maxSharedA; maxSharedA = 0;
  if (maxSharedB != 0) delete [] maxSharedB; maxSharedB = 0;
  if (dataShared !=0) delete [] dataShared; dataShared =0;
  if (vectX !=0) delete vectX; vectX =0;
  if (vectB !=0) delete vectB; vectB =0;
  
  A = new double[size];
  B = new double[size];
  X = new double[size];
  sharedA = new double[numShared];
  sharedB = new double[numShared];
  maxSharedA = new double[maxShared];
  maxSharedB = new double[maxShared];

  dataShared = new double[1];
  if ( X!= 0)
    vectX = new Vector(X,size);
  if (B != 0)
    vectB = new Vector(B,size);

  if (A == 0 || B == 0 || X == 0 || sharedA == 0 || sharedB == 0 ||maxSharedA == 0 || maxSharedB == 0 || dataShared == 0 || vectX == 0 || vectB == 0 || posLocKey ==0) {
    opserr << "ERROR MPIDiagonalSOE::setSize() - ";
    opserr << " ran out of memory for size: " << size << endln;
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (sharedA != 0) delete [] sharedA;
    if (sharedB != 0) delete [] sharedB;
    if (maxSharedA != 0) delete [] maxSharedA;
    if (maxSharedB != 0) delete [] maxSharedB;
    if (dataShared !=0) delete [] dataShared;
    if (posLocKey !=0) delete [] posLocKey;
    size = 0;
    return -1;
  }

  // initiallize the vectors
  for (int l=0; l<size; l++) {
    A[l] = 0;
    B[l] = 0;
    X[l] = 0;
    posLocKey[l]= -1;
  }
  for (int l=0; l<numShared; l++) {
    sharedA[l] = 0;
    sharedB[l] = 0;
  }
  for (int l=0; l<maxShared; l++) {
    maxSharedA[l] = 0;
    maxSharedB[l] = 0;
  }
  }                                   // Ladruno (ADR-74): end dc.s.shared

  //
  // renumber DOFs 0 through size
  //
  if (theModel == 0) {
    opserr << "WARNING MPIDiagonalSOE::setSize - no AnalysisModel\n";
    exit(-1);
  } 
  else {
    int local_dof =0;
    int nonlocal_dof = size-numShared;

    { OPS_PROFILE_SCOPE("dc.s.localize"); // Ladruno (ADR-74): global->local renumber, binary search per dof
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
      const ID &theID = dofPtr->getID();
      for (int i=0; i<theID.Size(); i++) {
	int dof = theID(i);
	if (dof >= 0) {
	  int newDOF = myDOFs.getLocationOrdered(dof);
	  dofPtr->setID(i, newDOF);
	}
      }
    }
    }
    // iterate through the FE_Element getting them to set their IDs
    { OPS_PROFILE_SCOPE("dc.s.eleid");   // Ladruno (ADR-74): FE_Element ID re-cache
    FE_EleIter &theEle = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
      elePtr->setID();
    }
  }
  delete [] sharedDOFs;
  // invoke setSize() on the Solver
  // this is a hack still for unknown reasons here
  MPIDiagonalSolver *the_Solver = (MPIDiagonalSolver *)this->getSolver();
  the_Solver->setLinearSOE(*this);
  int solverOK = the_Solver->setSize();
  

  if (solverOK < 0) {
    opserr << "WARNING MPIDiagonalSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  delete [] sendSize;
  delete [] recvSize;
  delete [] tmpsendSize;
  delete [] max;
 
 
  te = MPI_Wtime();
  //fprintf(stderr, "SetSize TIME  %g \n", te-ts);

  return result;
}


const Vector&
MPIDiagonalSOE::getpartofA(Vector& At, const ID& ids)
{
  if (A == 0) {
    opserr << "FATAL MPIDiagonalSOE::getA - A == 0";
    exit(-1);
  } 
  else if (isAfactored) {
    for (int i=0; i<ids.Size(); i++)
      At(i) = 1.0/A[ids(i)];
  }
  else {
    for (int i=0; i<ids.Size(); i++)
      At(i) = A[ids(i)];
  }
  return At;
}

int 
MPIDiagonalSOE::addA(const Matrix &m, const ID &id, double fact)
{
 
  //
  // NEED TO COMMENT OUT THE IF
  // if want to update A
  if (isAfactored)
    return 0;
  // check for a quick return 
  if (fact == 0.0)  return 0;
  int dof;
  int loc;
#ifdef _G3DEBUG
  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "MPIDiagonalSOE::addA()	- Matrix and ID not of similar sizes\n";
    return -1;
  }
#endif

  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] += m(i,i);

        // Lump diagonal logic
        if (lumpDiagonal) {
                for (int j = 0; j < i; j++) {
                  A[pos] += m(j, i) ;
                }
                for (int j = i + 1; j < id.Size(); j++) {
                  A[pos] += m(j, i) ;
                }
        }

	dof = myDOFs[pos];
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	  sharedA[loc] = A[pos];
	  posLocKey[pos] = loc;
	}
      }
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] -= m(i,i);

        // Lump diagonal logic
        if (lumpDiagonal) {
                for (int j = 0; j < i; j++) {
                  A[pos] -= m(j, i) ;
                }
                for (int j = i + 1; j < id.Size(); j++) {
                  A[pos] -= m(j, i) ;
                }
        }

	dof = myDOFs[pos];
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	  sharedA[loc] = A[pos];
	  posLocKey[pos] = loc;
	}
      }
    }
  } else {
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] += m(i,i) * fact;

        // Lump diagonal logic
        if (lumpDiagonal) {
                for (int j = 0; j < i; j++) {
                  A[pos] += m(j, i) * fact;
                }
                for (int j = i + 1; j < id.Size(); j++) {
                  A[pos] += m(j, i) * fact;
                }
        }

	dof = myDOFs[pos];
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	  sharedA[loc] = A[pos];
	  posLocKey[pos] = loc;
	}
      }
    }
  }	

  return 0;
}
 
    
int 
MPIDiagonalSOE::addB(const Vector &v, const ID &id, double fact)
{
  if (isAfactored) {
    // check for a quick return 
    if (fact == 0.0)  return 0;
    int dof;
    int loc;
#ifdef _G3DEBUG
    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
      opserr << "MPIDiagonalSOE::addB() -";
      opserr << " Vector and ID not of similar sizes\n";
      return -1;
    }    
#endif
  
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0)
	  B[pos] += v(i);
      }
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0)
	  B[pos] -= v(i);
      }
    } else {
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0)
	  B[pos] += v(i) * fact;
      }
    }
  }  else {
    // check for a quick return 
    if (fact == 0.0)  return 0;
    int dof;
    int loc;
#ifdef _G3DEBUG
    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
      opserr << "MPIDiagonalSOE::addB() -";
      opserr << " Vector and ID not of similar sizes\n";
      return -1;
    }    
#endif
  
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0) {
	  B[pos] += v(i);
	  loc = posLocKey[pos];
	  if ( (loc >=0) && (loc <numShared) )
	    sharedB[loc] = B[pos];
	}
      }
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0) {
	  B[pos] -= v(i);
	  loc = posLocKey[pos];
	  if ( (loc >=0) && (loc <numShared) )
	    sharedB[loc] = B[pos];
	}
      }
    } else {
      for (int i=0; i<id.Size(); i++) {
	int pos = id(i);
	if (pos <size && pos >= 0) {
	  B[pos] += v(i) * fact;
	  loc = posLocKey[pos];
	  if ( (loc >=0) && (loc <numShared) )
	    sharedB[loc] = B[pos];
	}
      }
    }
  }

  return 0;
}


    
int
MPIDiagonalSOE::setB(const Vector &v, double fact)
{
  //  std::cout << " setB " << processID << "\n";
  

  // check for a quick return 
  if (fact == 0.0)  return 0;
  int dof;
  int loc;
  if (v.Size() != size) {
    opserr << "WARNING MPIDiagonalSOE::setB() -";
    opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
    return -1;
  }
  
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<size; i++) {
      B[i] = v(i);
      dof = myDOFs[i];
      /*
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	sharedB[loc] = B[i];
	}
      */
      loc = posLocKey[i];
      if ( (loc >=0) && (loc <numShared) )
	sharedB[loc] = B[i];
      //std::cout << " setB pid : " << processID << " @ pos : " << i << " we got loc : " << loc << "\n";
      //std::cout << " OOOOPSSS \n: ";
    }
  } 
  else if (fact == -1.0) {
    for (int i=0; i<size; i++) {
      B[i] = -v(i);
      dof = myDOFs[i];
      /*
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	sharedB[loc] = B[i];
	}
      */
      loc = posLocKey[i];
      if ( (loc >=0) && (loc < numShared) )
	sharedB[loc] = B[i];
      //std::cout << " setB pid : " << processID << " @ pos : " << i << " we got loc : " << loc << "\n";
      //std::cout << " OOOOPSSS \n: ";
    }
  } 
  else {
    for (int i=0; i<size; i++) {
      B[i] = v(i) * fact;
      dof = myDOFs[i];
      /*
	loc = myDOFsShared.getLocationOrdered(dof);
	if ( (loc >=0 ) && (loc < numShared )) {
	sharedB[loc] = B[i];
	}
      */
      loc = posLocKey[i];
      if ( (loc >=0) && (loc < numShared) )
	sharedB[loc] = B[i];
      //std::cout << " setB pid : " << processID << " @ pos : " << i << " we got loc : " << loc << "\n";
      //std::cout << " OOOOPSSS \n: ";
    }
  }	
  return 0;
}

void 
MPIDiagonalSOE::zeroA(void)
{
  // NEED TO COMMENT OUT THE IF
  // if want to update A
  if (!isAfactored) {
    double *Aptr = A;
    for (int i=0; i<size; i++)
      *Aptr++ = 0;
    Aptr = sharedA;
    for (int i=0; i<numShared; i++)
      *Aptr++ = 0;
    isAfactored = false;
  }
}

void 
MPIDiagonalSOE::zeroB(void)
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
  Bptr = sharedB;
  for (int i=0; i<numShared; i++)
    *Bptr++ = 0;
}


void 
MPIDiagonalSOE::setX(int loc, double value)
{
  if (loc < size && loc >=0)
    X[loc] = value;
}

void 
MPIDiagonalSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
MPIDiagonalSOE::getX(void)
{
  if (vectX == 0) {
    opserr << "FATAL MPIDiagonalSOE::getX - vectX == 0";
    exit(-1);
  }    
  //opserr << "MPIDiagonalSOE::getX(void) : " << vectX->Size() << endln;
  return *vectX;
}

const Vector &
MPIDiagonalSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL MPIDiagonalSOE::getB - vectB == 0";
    exit(-1);
  }        
  return *vectB;
}

double 
MPIDiagonalSOE::normRHS(void)
{
  double norm =0.0;
  for (int i=0; i<size; i++) {
    double Yi = B[i];
    norm += Yi*Yi;
  }
  return sqrt(norm);
  
}    


int
MPIDiagonalSOE::setDiagonalSolver(MPIDiagonalSolver &newSolver)
{
  newSolver.setLinearSOE(*this);
  
  if (size != 0) {
    int solverOK = newSolver.setSize();
    if (solverOK < 0) {
      opserr << "WARNING:MPIDiagonalSOE::setSolver :";
      opserr << "the new solver could not setSeize() - staying with old\n";
      return -1;
    }
  }
  
  return this->setSolver(newSolver);
}


int 
MPIDiagonalSOE::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}


int 
MPIDiagonalSOE::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  return 0;
}

int
MPIDiagonalSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];


  localCol = new ID *[nChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  return 0;
}

int 
MPIDiagonalSOE::setAnalysisModel(AnalysisModel &theAnalysisModel)
{
  //opserr << " MPIDiagonalSOE::setAnalysisModel has been called \n";
  if (&theAnalysisModel !=0) {
    theModel = &theAnalysisModel;
  } else {
    opserr << "MPIDiagonalSOE: FATAL ERROR: theAnalysisModel is 0 \n";
    exit(-1);
  }
  return 0;
}


void 
MPIDiagonalSOE::quickSort(ID& array, int array_size)
{
  q_sort(array, 0, array_size - 1);
}


void 
MPIDiagonalSOE::q_sort(ID& array, int left, int right)
{
  int pivot, l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = array[left];
  while (left < right)
    {
      while ((array[right] >= pivot) && (left < right))
	right--;
      if (left != right)
	{
	  array[left] = array[right];
	  left++;
	}
      while ((array[left] <= pivot) && (left < right))
	left++;
      if (left != right)
	{
	  array[right] = array[left];
	  right--;
	}
    }
  array[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort(array, left, pivot-1);
  if (right > pivot)
    q_sort(array, pivot+1, right);
}



void
MPIDiagonalSOE::intersections(ID& arrayA, ID& arrayB, int sizeA, int sizeB, int* shared, int* sharedDOFs)
{

  /// Assumes a, b arrays are SORTED!!!!!!!!!!!!!!
  /// SO make sure you have sorted your dofs before calling this !!!
  /// shared array is equal to thesize of smallest of those two input arrays

  shared[0]=0; // start by zero shared dofs

  int i =0;
  int j =0;
  //  opserr << " INside interesctions with pID : " << processID << " sizeA : " << sizeA << " sizeB " << sizeB << "\n";
  while ( (i<sizeA) && (j<sizeB) ) {
    if ( arrayA[i] == arrayB[j] ) {
      sharedDOFs[i] = 1;
      shared[0]++;
      i++;
      j++;
    }
    else if (arrayA[i] > arrayB[j] ) {
      j++;
    }
    else {
      i++;
    }
  }
  
}


void MPIDiagonalSOE::dontUpdateA()
{
  updateA = false;
}








