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

// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A


#include <PARDISOGenLinSOE.h>
#include <PARDISOGenLinSolver.h>

#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

PARDISOGenLinSOE::PARDISOGenLinSOE(PARDISOGenLinSolver &the_Solver)
	:LinearSOE(the_Solver, LinSOE_TAGS_PARDISOGenLinSOE),
	size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
	vectX(0), vectB(0),
	Asize(0), Bsize(0),
	factored(false), matType(0)
{
	the_Solver.setLinearSOE(*this);
}


// Ladruno ADR-75 P1d: symmetric variant. matType 1 (SPD) / 2 (symmetric
// indefinite) store only the upper triangle; anything else falls back to the
// full unsymmetric CSR rather than half-storing something PARDISO would then
// read with the wrong mtype.
PARDISOGenLinSOE::PARDISOGenLinSOE(PARDISOGenLinSolver &the_Solver, int _matType)
	:LinearSOE(the_Solver, LinSOE_TAGS_PARDISOGenLinSOE),
	size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
	vectX(0), vectB(0),
	Asize(0), Bsize(0),
	factored(false), matType(_matType)
{
	if (matType < 0 || matType > 2) {
		opserr << "WARNING PARDISOGenLinSOE - unknown matrixType (" << matType
		       << "); unsymmetric storage assumed\n";
		matType = 0;
	}
	the_Solver.setLinearSOE(*this);
}


PARDISOGenLinSOE::~PARDISOGenLinSOE()
{
	if (A != 0) delete[] A;
	if (B != 0) delete[] B;
	if (X != 0) delete[] X;
	if (rowStartA != 0) delete[] rowStartA;
	if (colA != 0) delete[]colA;
	if (vectX != 0) delete vectX;
	if (vectB != 0) delete vectB;
}


int
PARDISOGenLinSOE::getNumEqn(void) const
{
	return size;
}

// Ladruno ADR-75 P1d: the solver reads this to pick `mtype` (see the header).
int
PARDISOGenLinSOE::getMatType(void) const
{
	return matType;
}

int
PARDISOGenLinSOE::setSize(Graph &theGraph)
{

	int result = 0;
	int oldSize = size;
	size = theGraph.getNumVertex();

	// fist itearte through the vertices of the graph to get nnz
	Vertex *theVertex;
	int newNNZ = 0;
	VertexIter &theVertices = theGraph.getVertices();
	while ((theVertex = theVertices()) != 0) {
		const ID &theAdjacency = theVertex->getAdjacency();
		newNNZ += theAdjacency.Size() + 1; // the +1 is for the diag entry
	}

	// Ladruno ADR-75 P1d: half-storage. The connectivity Graph is structurally
	// symmetric (every adjacency appears in both directions), so the strictly
	// off-diagonal count is exactly halved and the diagonal is kept whole. Same
	// arithmetic as MumpsSOE::setSize()'s matType != 0 branch.
	if (matType != 0) {
		newNNZ -= size;
		newNNZ /= 2;
		newNNZ += size;
	}

	nnz = newNNZ;

	if (newNNZ > Asize) { // we have to get more space for A and colA
		if (A != 0)
			delete[] A;
		if (colA != 0)
			delete[] colA;

		A = new double[newNNZ];
		colA = new int[newNNZ];

		if (A == 0 || colA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::PARDISOGenLinSOE :";
			opserr << " ran out of memory for A and colA with nnz = ";
			opserr << newNNZ << " \n";
			size = 0; Asize = 0; nnz = 0;
			result = -1;
		}

		Asize = newNNZ;
	}

	// zero the matrix
	for (int i = 0; i < Asize; i++)
		A[i] = 0;

	factored = false;

	if (size > Bsize) { // we have to get space for the vectors

	// delete the old	
		if (B != 0) delete[] B;
		if (X != 0) delete[] X;
		if (rowStartA != 0) delete[] rowStartA;

		// create the new
		B = new double[size];
		X = new double[size];
		rowStartA = new int[size + 1];

		if (B == 0 || X == 0 || rowStartA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::PARDISOGenLinSOE :";
			opserr << " ran out of memory for vectors (size) (";
			opserr << size << ") \n";
			size = 0; Bsize = 0;
			result = -1;
		}
		else
			Bsize = size;
	}

	// zero the vectors
	for (int j = 0; j < size; j++) {
		B[j] = 0;
		X[j] = 0;
	}

	// create new Vectors objects
	if (size != oldSize) {
		if (vectX != 0)
			delete vectX;

		if (vectB != 0)
			delete vectB;

		vectX = new Vector(X, size);
		vectB = new Vector(B, size);
	}

	// fill in rowStartA and colA
	if (size != 0) {
		rowStartA[0] = 0 + 1;
		int startLoc = 0;
		int lastLoc = 0;
		for (int a = 0; a < size; a++) {

			theVertex = theGraph.getVertexPtr(a);
			if (theVertex == 0) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize :";
				opserr << " vertex " << a << " not in graph! - size set to 0\n";
				size = 0;
				return -1;
			}

			int vertexTag = theVertex->getTag();

			// Ladruno ADR-75 P1d: HARD BOUND on the fill.
			// `nnz` is a computed prediction — for matType != 0 it is
			// (nnz_full - n)/2 + n, which is exact ONLY while the Graph's
			// adjacency is strictly symmetric and self-loop-free. That holds
			// today (Graph::addEdge is two-sided and exit(0)s otherwise), but
			// the failure mode if it ever stops holding is a SILENT HEAP
			// OVERFLOW of colA/A, not a crash — memory corruption that would
			// surface arbitrarily far away. Cheap guard, unbounded payoff.
			// (nnz <= Asize always, so bounding at nnz is the stricter test.)
			if (lastLoc >= nnz) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : CSR fill would "
				          "overrun nnz=" << nnz << " at row " << a
				       << " (matType " << matType << ") - the graph is not the "
				          "symmetric, self-loop-free shape the nnz count assumes\n";
				size = 0;
				return -1;
			}
			colA[lastLoc++] = vertexTag + 1; // place diag in first fortran index start at 1
			const ID &theAdjacency = theVertex->getAdjacency();
			int idSize = theAdjacency.Size();

			// now we have to place the entries in the ID into order in colA
			// (PARDISO requires ASCENDING column indices within each row —
			//  Vertex adjacency is not guaranteed sorted, hence the insertion).
			for (int i = 0; i < idSize; i++) {

				int row = theAdjacency(i);

				// Ladruno ADR-75 P1d: upper triangle only when half-storing.
				// The diagonal already sits at startLoc and is smaller than
				// every entry kept here, so the insertion below can never
				// displace it and the row stays ascending.
				if (matType != 0 && row <= vertexTag)
					continue;

				if (lastLoc >= nnz) {          // see the bound note above
					opserr << "WARNING:PARDISOGenLinSOE::setSize : CSR fill "
					          "would overrun nnz=" << nnz << " in row " << a
					       << "'s adjacency (matType " << matType << ")\n";
					size = 0;
					return -1;
				}

				bool foundPlace = false;
				// find a place in colA for current col
				for (int j = startLoc; j < lastLoc; j++)
					if (colA[j] > row + 1) {
						// move the entries already there one further on
						// and place col in current location
						for (int k = lastLoc; k > j; k--)

							colA[k] = colA[k - 1];
						colA[j] = row + 1;
						foundPlace = true;
						j = lastLoc;
					}
				if (foundPlace == false) // put in at the end
					colA[lastLoc] = row + 1;

				lastLoc++;
			}
			rowStartA[a + 1] = lastLoc + 1;
			startLoc = lastLoc;
		}

		// ---- Ladruno ADR-75 P1d: verify the CSR contract, don't assume it ---
		// MKL PARDISO requires, per row: column indices STRICTLY ASCENDING, and
		// (for the symmetric mtypes especially) the diagonal PRESENT. Violate
		// either and PARDISO does not necessarily complain — it can return a
		// plausible wrong answer, which is the worst failure mode this fork has.
		//
		// The half-storage path relies on an ARGUMENT that the order survives:
		// the diagonal is written first, and only `row > vertexTag` entries are
		// kept, so nothing can sort ahead of it. That argument also silently
		// assumes `theVertex->getTag() == a` (the CSR row index), which the
		// pre-existing code has always assumed but never checked. This loop
		// checks all of it directly instead — O(nnz), negligible next to the
		// O(nnz x rowlen) insertion above, and it validates the unsymmetric
		// path for free.
		if (lastLoc != nnz) {
			opserr << "WARNING:PARDISOGenLinSOE::setSize : filled " << lastLoc
			       << " entries but nnz=" << nnz << " (matType " << matType
			       << ") - the graph is not the shape the nnz count assumes\n";
			size = 0;
			return -1;
		}
		// NB the diagonal is only FIRST in the half-stored case. In the
		// unsymmetric path it is sorted into place like any other column, so
		// the test is "present", not "at rs" — checking position there would
		// reject every correct matrix.
		for (int a = 0; a < size; a++) {
			int rs = rowStartA[a] - 1, re = rowStartA[a + 1] - 1;
			bool sawDiag = false;
			if (rs >= re) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
				       << " is empty (matType " << matType
				       << ") - PARDISO requires a diagonal entry in every row\n";
				size = 0;
				return -1;
			}
			for (int k = rs; k < re; k++) {
				if (k > rs && colA[k] <= colA[k - 1]) {
					opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
					       << " column indices are not strictly ascending at "
					       << k << " (matType " << matType
					       << ") - PARDISO requires ascending CSR\n";
					size = 0;
					return -1;
				}
				if (colA[k] == a + 1)
					sawDiag = true;
			}
			if (sawDiag == false) {
				opserr << "WARNING:PARDISOGenLinSOE::setSize : row " << a
				       << " has no diagonal entry (matType " << matType
				       << ") - PARDISO requires one, even if zero\n";
				size = 0;
				return -1;
			}
		}
	}

	// invoke setSize() on the Solver
	LinearSOESolver *the_Solver = this->getSolver();
	int solverOK = the_Solver->setSize();
	if (solverOK < 0) {
		opserr << "WARNING:PARDISOGenLinSOE::setSize :";
		opserr << " solver failed setSize()\n";
		return solverOK;
	}
	return result;
}

int
PARDISOGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
	// check for a quick return 
	if (fact == 0.0)
		return 0;

	int idSize = id.Size();

	// check that m and id are of similar size
	if (idSize != m.noRows() && idSize != m.noCols()) {
		opserr << "PARDISOGenLinSOE::addA() ";
		opserr << " - Matrix and ID not of similar sizes\n";
		return -1;
	}

	// Ladruno ADR-75 P1d: when half-storing, only the col >= row entries of the
	// element matrix have a home in A. Everything else is identical to the
	// unsymmetric path, so the filter is a single extra comparison rather than
	// a duplicated loop nest (MumpsSOE duplicates; there is no reason to).
	const bool halfStore = (matType != 0);

	if (fact == 1.0) { // do not need to multiply
		for (int i = 0; i < idSize; i++) {
			int row = id(i);
			if (row < size && row >= 0) {
				int startRowLoc = rowStartA[row] - 1;
				int endRowLoc = rowStartA[row + 1] - 1;
				for (int j = 0; j < idSize; j++) {
					int col = id(j);
					if (col < size && col >= 0 && (!halfStore || col >= row)) {
						// find place in A using colA
						for (int k = startRowLoc; k < endRowLoc; k++)
							if (colA[k] == col + 1) {
								A[k] += m(i, j);
								k = endRowLoc;
							}
					}
				}  // for j
			}
		}  // for i
	}
	else {
		for (int i = 0; i < idSize; i++) {
			int row = id(i);
			if (row < size && row >= 0) {
				int startRowLoc = rowStartA[row] - 1;
				int endRowLoc = rowStartA[row + 1] - 1;
				for (int j = 0; j < idSize; j++) {
					int col = id(j);
					if (col < size && col >= 0 && (!halfStore || col >= row)) {
						// find place in A using colA
						for (int k = startRowLoc; k < endRowLoc; k++)
							if (colA[k] == col + 1) {
								A[k] += fact * m(i, j);
								k = endRowLoc;
							}
					}
				}  // for j
			}
		}  // for i
	}
	return 0;
}


int
PARDISOGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
	// check for a quick return 
	if (fact == 0.0)  return 0;

	int idSize = id.Size();
	// check that m and id are of similar size
	if (idSize != v.Size()) {
		opserr << "PARDISOGenLinSOE::addB() ";
		opserr << " - Vector and ID not of similar sizes\n";
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] += v(i);
		}
	}
	else if (fact == -1.0) { // do not need to multiply if fact == -1.0
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] -= v(i);
		}
	}
	else {
		for (int i = 0; i < idSize; i++) {
			int pos = id(i);
			if (pos < size && pos >= 0)
				B[pos] += v(i) * fact;
		}
	}

	return 0;
}


int
PARDISOGenLinSOE::setB(const Vector &v, double fact)
{
	// check for a quick return 
	if (fact == 0.0)  return 0;


	if (v.Size() != size) {
		opserr << "WARNING BandGenLinSOE::setB() -";
		opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i = 0; i < size; i++) {
			B[i] = v(i);
		}
	}
	else if (fact == -1.0) {
		for (int i = 0; i < size; i++) {
			B[i] = -v(i);
		}
	}
	else {
		for (int i = 0; i < size; i++) {
			B[i] = v(i) * fact;
		}
	}
	return 0;
}

void
PARDISOGenLinSOE::zeroA(void)
{
	double *Aptr = A;
	for (int i = 0; i < Asize; i++)
		*Aptr++ = 0;

	factored = false;
}

void
PARDISOGenLinSOE::zeroB(void)
{
	double *Bptr = B;
	for (int i = 0; i < size; i++)
		*Bptr++ = 0;
}

void
PARDISOGenLinSOE::setX(int loc, double value)
{
	if (loc < size && loc >= 0)
		X[loc] = value;
}

void
PARDISOGenLinSOE::setX(const Vector &x)
{
	if (x.Size() == size && vectX != 0)
		*vectX = x;
}

const Vector &
PARDISOGenLinSOE::getX(void)
{
	if (vectX == 0) {
		opserr << "FATAL PARDISOGenLinSOE::getX - vectX == 0";
		exit(-1);
	}
	return *vectX;
}

const Vector &
PARDISOGenLinSOE::getB(void)
{
	if (vectB == 0) {
		opserr << "FATAL PARDISOGenLinSOE::getB - vectB == 0";
		exit(-1);
	}
	return *vectB;
}

double
PARDISOGenLinSOE::normRHS(void)
{
	double norm = 0.0;
	for (int i = 0; i < size; i++) {
		double Yi = B[i];
		norm += Yi * Yi;
	}
	return sqrt(norm);

}


int
PARDISOGenLinSOE::setPARDISOGenLinSolver(PARDISOGenLinSolver &newSolver)
{
	newSolver.setLinearSOE(*this);

	if (size != 0) {
		int solverOK = newSolver.setSize();
		if (solverOK < 0) {
			opserr << "WARNING:PARDISOGenLinSOE::setSolver :";
			opserr << "the new solver could not setSeize() - staying with old\n";
			return -1;
		}
	}

	return this->LinearSOE::setSolver(newSolver);
}


int
PARDISOGenLinSOE::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}

int
PARDISOGenLinSOE::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	return 0;
}
