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
// $Date: 2009-05-11 20:56:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.h, revA"

#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

UmfpackGenLinSOE::UmfpackGenLinSOE(UmfpackGenLinSolver &the_Solver)
    :LinearSOE(the_Solver, LinSOE_TAGS_UmfpackGenLinSOE), X(), B(), Ap(), Ai(), Ax(),
     factored(false),  // Ladruno (ADR-40 rank 8/10)
     missWarned(0)     // Ladruno (ADR-75 P1g)
{
    the_Solver.setLinearSOE(*this);
}


UmfpackGenLinSOE::UmfpackGenLinSOE()
    :LinearSOE(LinSOE_TAGS_UmfpackGenLinSOE), X(), B(), Ap(), Ai(), Ax(),
     factored(false),  // Ladruno (ADR-40 rank 8/10)
     missWarned(0)     // Ladruno (ADR-75 P1g)
{
}


UmfpackGenLinSOE::~UmfpackGenLinSOE()
{
}


int
UmfpackGenLinSOE::getNumEqn(void) const
{
    return X.Size();
}

int
UmfpackGenLinSOE::setSize(Graph &theGraph)
{
    // Ladruno (ADR-40 rank 8/10): the sparsity structure is about to be rebuilt
    // (new Symbolic) -> any persisted Numeric is stale. Reset here at the top so
    // even the early-error returns below leave the safe (unfactored) state; the
    // solver's setSize() frees the persisted Numeric before rebuilding Symbolic.
    factored = false;

    int size = theGraph.getNumVertex();
    if (size < 0) {
	opserr<<"size of soe < 0\n";
	return -1;
    }

    // fist itearte through the vertices of the graph to get nnz
    Vertex *theVertex;
    int nnz = 0;
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
	const ID &theAdjacency = theVertex->getAdjacency();
	nnz += theAdjacency.Size() +1; // the +1 is for the diag entry
    }

    // resize A, B, X
    Ap.reserve(size+1);
    Ai.reserve(nnz);
    Ax.resize(nnz,0.0);
    B.resize(size);
    B.Zero();
    X.resize(size);
    X.Zero();

    // fill in Ai and Ap
    Ap.push_back(0);
    for (int a=0; a<size; a++) {

	theVertex = theGraph.getVertexPtr(a);
	if (theVertex == 0) {
	    opserr << "WARNING:UmfpackGenLinSOE::setSize :";
	    opserr << " vertex " << a << " not in graph! - size set to 0\n";
	    size = 0;
	    return -1;
	}

	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	ID col(0,idSize+1);

	// diagonal
	col.insert(theVertex->getTag());

	// now we have to place the entries in the ID into order in Ai
	for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    col.insert(row);
	}

	// copy to Ai
	for (int i=0; i<col.Size(); i++) {
	    Ai.push_back(col(i));
	}

	// set Ap
	Ap.push_back(Ap[a]+col.Size());
    }

    // Ladruno (ADR-75 P1g): ENFORCE the invariant addA's binary search depends on.
    // Each CSC column of Ai must be strictly ascending. It is, by construction —
    // every column is built through ID::insert, a binary-search insertion into a
    // sorted array with dedup. But addA's ops_umfpack_findRow would return a WRONG
    // (or absent) slot on an unsorted column, and the failure mode is a silently
    // wrong tangent, not a crash. So check rather than assume: this is O(nnz) once
    // per setSize, against an assembly loop that runs O(idSize^2 * log collen) per
    // element for the whole analysis — free, and it converts a future silent
    // breakage into a loud one. (Same reasoning, and the same guard, as ADR-75 P1f
    // added to PARDISOGenLinSOE::setSize.)
    for (int a = 0; a < size; a++) {
	for (int k = Ap[a] + 1; k < Ap[a+1]; k++) {
	    if (Ai[k] <= Ai[k-1]) {
		opserr << "WARNING:UmfpackGenLinSOE::setSize : column " << a
		       << " row indices are not strictly ascending at Ai[" << k
		       << "] (" << Ai[k] << " after " << Ai[k-1]
		       << ") - addA's binary search requires ascending CSC\n";
		return -1;
	    }
	}
    }

    // Ladruno (ADR-40 rank 8/10): structure rebuilt (new Symbolic) -> any
    // persisted Numeric is stale; mark unfactored. The solver's setSize() frees
    // the persisted Numeric before rebuilding Symbolic.
    factored = false;

    // invoke setSize() on the Solver
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:UmfpackGenLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }
    return 0;
}

// Ladruno (ADR-75 P1g): binary search for a row index inside one CSC column of Ai.
// Mirrors ops_pardiso_findCol (ADR-75 P1f) — same contract, same shape, so the two
// SOEs stay comparable. Returns the index into Ai/Ax, or -1 if `row` is not in the
// column (a structurally absent entry, which addA skips exactly as the old linear
// scan did by falling out of its loop).
//
// LEGAL because each CSC column of Ai is STRICTLY ASCENDING. That is not an
// assumption: setSize builds every column through ID::insert, which is a
// binary-search insertion into a sorted array with dedup (`if (x == dataMiddle)
// return 1; // already there`) — so the column is a sorted set by construction. It
// is additionally CHECKED at the end of setSize (see there), because "guaranteed by
// construction" silently becomes "wrong answer" the day the construction changes.
static inline int
ops_umfpack_findRow(const int *Ai, int lo, int hi, int row)
{
    while (lo < hi) {
	const int mid = lo + ((hi - lo) >> 1);
	const int r = Ai[mid];
	if (r < row)
	    lo = mid + 1;
	else if (r > row)
	    hi = mid;
	else
	    return mid;
    }
    return -1;
}

int
UmfpackGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return
    if (fact == 0.0) return 0;

    int idSize = id.Size();

    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "UmfpackGenLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }

    // Ladruno (ADR-75 P1g): a FREE miss detector. The `k >= 0` test already exists,
    // so `else missing = 1` costs a never-taken branch. It turns "an element entry
    // has no CSC slot => its stiffness is silently discarded" from undetectable into
    // a once-per-SOE warning; the old linear scan just fell off the end of the
    // column. (Mirrors ADR-75 P1f in PARDISOGenLinSOE.)
    int missing = 0;

    int size = X.Size();
    if (fact == 1.0) { // do not need to multiply
	for (int j=0; j<idSize; j++) {
	    int col = id(j);
	    if (col<0 || col>=size) {
		continue;
	    }
	    for (int i=0; i<idSize; i++) {
		int row = id(i);
		if (row<0 || row>=size) {
		    continue;
		}

		// find place in A (Ladruno ADR-75 P1g: was a linear scan of the
		// whole column for each of idSize^2 entries -> O(idSize^2 * collen))
		const int k = ops_umfpack_findRow(&Ai[0], Ap[col], Ap[col+1], row);
		if (k >= 0) {
		    Ax[k] += m(i,j);
		} else {
		    missing = 1;
		}
	    }
	}
    } else {
	for (int j=0; j<idSize; j++) {
	    int col = id(j);
	    if (col<0 || col>=X.Size()) {
		continue;
	    }
	    for (int i=0; i<idSize; i++) {
		int row = id(i);
		if (row<0 || row>=X.Size()) {
		    continue;
		}

		// find place in A (Ladruno ADR-75 P1g: binary search, see above)
		const int k = ops_umfpack_findRow(&Ai[0], Ap[col], Ap[col+1], row);
		if (k >= 0) {
		    Ax[k] += fact*m(i,j);
		} else {
		    missing = 1;
		}
	    }
	}
    }

    // Ladruno (ADR-75 P1g): reported ONCE per SOE, not per element.
    if (missing == 1 && missWarned == 0) {
	missWarned = 1;
	opserr << "WARNING UmfpackGenLinSOE::addA() - an element entry has no slot in "
	       << "the sparsity pattern; its contribution was DISCARDED. The tangent is "
	       << "wrong. This usually means the DOF connectivity graph setSize() was "
	       << "built from does not cover the element's ID. Reported once.\n";
    }

    return 0;
}


int
UmfpackGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "UmfpackGenLinSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    int size = B.Size();
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] += v(i) * fact;
	}
    }
    
    return 0;
}


int
UmfpackGenLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  {
	B.Zero();
	return 0;
    }

    int size = B.Size();
    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    B[i] = v(i) * fact;
	}
    }
    
    return 0;
}

void
UmfpackGenLinSOE::zeroA(void)
{
    Ax.assign(Ax.size(),0.0);
    // Ladruno (ADR-40 rank 8/10): A is being reassembled (zeroA always precedes
    // the addA loop in IncrementalIntegrator::formTangent) -> the persisted
    // factorization is invalid. Matches BandGenLinSOE::zeroA().
    factored = false;
}

void
UmfpackGenLinSOE::zeroB(void)
{
    B.Zero();
}

void
UmfpackGenLinSOE::setX(int loc, double value)
{
    if (loc<X.Size() && loc>=0) {
	X(loc) = value;
    }
}


void
UmfpackGenLinSOE::setX(const Vector &x)
{
    if (x.Size() == X.Size()) {
	X = x;
    }
}


const Vector &
UmfpackGenLinSOE::getX(void)
{
    return X;
}

const Vector &
UmfpackGenLinSOE::getB(void)
{
    return B;
}

double
UmfpackGenLinSOE::normRHS(void)
{
    return B.Norm();
}

int
UmfpackGenLinSOE::saveSparseA(OPS_Stream& output, int baseIndex)
{
    int size = X.Size();
    if (size == 0) {
        opserr << "WARNING: UmfpackGenLinSOE::saveSparseA() - size is 0\n";
        return -1;
    }
    
    // Assume the header is already written to output stream
    int nnz = Ax.size();
    output << size << " " << size << " " << nnz << "\n";
    
    // Write the sparse matrix entries
    int nnz_written = 0;
    for (int col = 0; col < size; col++) {
        for (int k = Ap[col]; k < Ap[col+1]; k++) {
            int row = Ai[k];
            double value = Ax[k];
            output << (row + baseIndex) << " " << (col + baseIndex) << " " << value << "\n";
            nnz_written++;
        }
    }
    if (nnz_written != nnz) {
        opserr << "WARNING: UmfpackGenLinSOE::saveSparseA() - nnz_written != nnz\n";
        return -1;
    }
    
    return 0;
}

int
UmfpackGenLinSOE::getSparseA(ID& rowIndices, ID& colIndices, Vector& values, int baseIndex)
{
    int size = X.Size();
    if (size == 0) {
        opserr << "WARNING: UmfpackGenLinSOE::getSparseA() - size is 0\n";
        return -1;
    }
    
    int nnz = Ax.size();
    rowIndices.resize(nnz);
    colIndices.resize(nnz);
    values.resize(nnz);
    
    // Fill vectors with non-zero elements
    int idx = 0;
    for (int col = 0; col < size; col++) {
        for (int k = Ap[col]; k < Ap[col+1]; k++) {
            int row = Ai[k];
            double value = Ax[k];
            rowIndices(idx) = row + baseIndex;
            colIndices(idx) = col + baseIndex;
            values(idx) = value;
            idx++;
        }
    }
    
    return 0;
}

int
UmfpackGenLinSOE::getSparseA(std::vector<int>& rowIndices, std::vector<int>& colIndices, std::vector<double>& values, int baseIndex)
{
    int size = X.Size();
    if (size == 0) {
        opserr << "WARNING: UmfpackGenLinSOE::getSparseA() - size is 0\n";
        return -1;
    }
    
    int nnz = Ax.size();
    rowIndices.resize(nnz);
    colIndices.resize(nnz);
    values.resize(nnz);
    
    // Fill vectors with non-zero elements
    int idx = 0;
    for (int col = 0; col < size; col++) {
        for (int k = Ap[col]; k < Ap[col+1]; k++) {
            int row = Ai[k];
            double value = Ax[k];
            rowIndices[idx] = row + baseIndex;
            colIndices[idx] = col + baseIndex;
            values[idx] = value;
            idx++;
        }
    }
    
    return 0;
}

int
UmfpackGenLinSOE::setUmfpackGenLinSolver(UmfpackGenLinSolver &newSolver)
{
    // Ladruno (ADR-40 rank 8/10): a swapped-in solver holds no factorization for
    // this matrix. The solver's Numeric==0 guard already forces a refactor, but
    // reset the flag here too so the factored<->Numeric invariant is locally
    // self-evident (newSolver.setSize() below also frees any Numeric).
    factored = false;
    newSolver.setLinearSOE(*this);
    if (X.Size() != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:UmfpackGenLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    return this->LinearSOE::setSolver(newSolver);
}


int
UmfpackGenLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int
UmfpackGenLinSOE::recvSelf(int cTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
    return 0;
}
