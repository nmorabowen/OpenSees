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
#include <new>            // Ladruno ADR-75 P1d: std::nothrow (see setSize)

#include <Channel.h>
#include <FEM_ObjectBroker.h>

PARDISOGenLinSOE::PARDISOGenLinSOE(PARDISOGenLinSolver &the_Solver)
	:LinearSOE(the_Solver, LinSOE_TAGS_PARDISOGenLinSOE),
	size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
	vectX(0), vectB(0),
	Asize(0), Bsize(0),
	factored(false), matType(0), asymWarned(0), asymBudget(0), missWarned(0)
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
	factored(false), matType(_matType), asymWarned(0), asymBudget(0), missWarned(0)
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

	// Ladruno ADR-75 P1d: arm the addA symmetry check for the first few TANGENT
	// ASSEMBLIES of this pattern (counted down in zeroA, which runs once per
	// formTangent). Budgeting in passes rather than in addA calls is what makes
	// the cost independent of model size and of how many Newton iterations the
	// run takes — 3 assemblies, not 3 assemblies' worth of a 45-assembly run.
	if (matType != 0 && asymWarned == 0)
		asymBudget = 3;

	// Ladruno ADR-75 P1d (adversarial review): this OOM path was DEAD CODE.
	// Plain `new` THROWS; it never returns 0, so `if (A == 0 ...)` could not
	// fire, and nothing in PythonModule.cpp / tclMain.cpp catches std::bad_alloc
	// — an out-of-memory setSize terminated the process with no OpenSees
	// diagnostic at all. Worse, `Asize = newNNZ` sat OUTSIDE the failure branch
	// and overwrote the `Asize = 0` recovery, so the zeroing loop below would
	// have walked a null pointer had the branch ever been reachable. And the
	// `delete[]` ran BEFORE the allocation, leaving dangling pointers on throw.
	//
	// This matters here specifically: P1c's headline is "UmfPack OOM at 86,490
	// DOF" and this whole lane exists to go past that wall — so the first model
	// that does not fit is exactly the case that must REPORT rather than vanish.
	if (newNNZ > Asize) { // we have to get more space for A and colA
		if (A != 0)
			delete[] A;
		if (colA != 0)
			delete[] colA;
		A = 0; colA = 0;              // never leave these dangling

		A = new (std::nothrow) double[newNNZ];
		colA = new (std::nothrow) int[newNNZ];

		if (A == 0 || colA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::setSize :";
			opserr << " ran out of memory for A and colA with nnz = ";
			opserr << newNNZ << " ("
			       << (newNNZ * (sizeof(double) + sizeof(int))) / (1024.0 * 1024.0)
			       << " MB requested)\n";
			if (A != 0) { delete[] A; A = 0; }
			if (colA != 0) { delete[] colA; colA = 0; }
			size = 0; Asize = 0; nnz = 0;
			return -1;            // do NOT fall through into the zeroing loop
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
		B = 0; X = 0; rowStartA = 0;   // same dangling-on-throw fix as above

		// create the new
		B = new (std::nothrow) double[size];
		X = new (std::nothrow) double[size];
		rowStartA = new (std::nothrow) int[size + 1];

		if (B == 0 || X == 0 || rowStartA == 0) {
			opserr << "WARNING PARDISOGenLinSOE::setSize :";
			opserr << " ran out of memory for vectors (size) (";
			opserr << size << ") \n";
			size = 0; Bsize = 0;
			return -1;             // was: fall through into a null deref
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

// Ladruno ADR-75 P1f: locate `col1` (a ONE-based column index) in the ascending
// run colA[lo, hi). Returns the offset into A, or -1 if the pattern has no such
// entry.
//
// The stock scatter did this with a LINEAR scan, so assembling one element
// matrix cost O(idSize^2 * rowlen): for a 3D brick that is 24*24 lookups each
// rescanning a ~81-entry row. ADR-75b's L3-0 profile measured it at 1699 ms of
// an ~11.9 s step — 1.28x slower than UmfPack's equivalent scatter — and it
// only grew in relative terms as P1a/P1d/P1e made the SOLVE cheaper.
//
// Binary search is legal here because the ascending-CSR invariant is not merely
// assumed, it is ENFORCED: setSize() above rejects the matrix (size = 0,
// return -1) if any row's column indices are not strictly ascending. That check
// is what makes this a safe drop-in rather than a bet.
static inline int
ops_pardiso_findCol(const int *colA, int lo, int hi, int col1)
{
	while (lo < hi) {
		const int mid = lo + ((hi - lo) >> 1);
		const int c = colA[mid];
		if (c < col1)
			lo = mid + 1;
		else if (c > col1)
			hi = mid;
		else
			return mid;
	}
	return -1;
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

	// Ladruno ADR-75 P1d (adversarial review): DETECT the asymmetry we are about
	// to discard. Half-storage on a genuinely unsymmetric tangent produces a
	// converged, plausible, WRONG answer with the upper triangle reflected —
	// the single worst failure mode in this feature, and previously only
	// documented, not caught. The check is nearly free: the lower-triangle
	// entries are visited anyway (to be skipped), and m(j,i) is the mirror
	// already in cache. Reported ONCE per SOE, not per element.
	// BUDGETED (in assembly passes — see setSize/zeroA), and fused into a single
	// lower-triangle pass. This is a diagnostic, not physics: on a correctly
	// symmetric model `asymWarned` never latches, so without a budget the check
	// would tax the assembly loop of a PERFORMANCE feature for the entire
	// analysis.
	//
	// LIMITATION, stated rather than papered over: it samples the first few
	// tangent assemblies. An element whose tangent only turns unsymmetric LATER
	// — contact closing at step 50, a non-associated flow rule after first
	// yield — will NOT be caught. It is a first-pass misconfiguration detector,
	// not a guarantee, and `-matrixType 2` still requires the author to know
	// their tangent is symmetric.
	if (halfStore && asymWarned == 0 && asymBudget > 0 &&
	    idSize == m.noRows() && idSize == m.noCols()) {
		double worstDev = 0.0, scale = 0.0;
		for (int i = 0; i < idSize; i++)
			for (int j = 0; j < i; j++) {
				double u = m(i, j), v = m(j, i);
				double au = u < 0.0 ? -u : u;
				double av = v < 0.0 ? -v : v;
				if (au > scale) scale = au;
				if (av > scale) scale = av;
				double d = u - v;
				if (d < 0.0) d = -d;
				if (d > worstDev) worstDev = d;
			}
		if (scale > 0.0 && worstDev > 1.0e-8 * scale) {
			asymWarned = 1;
			opserr << "WARNING PARDISOGenLinSOE: an assembled element matrix is "
			          "UNSYMMETRIC (max deviation " << worstDev << " vs scale "
			       << scale << ") but `system Pardiso -matrixType " << matType
			       << "` stores only the upper triangle.\n"
			          "     The lower half is being DISCARDED, so this run solves "
			          "the REFLECTED system and may converge to a WRONG answer.\n"
			          "     Use the default `system Pardiso` (-matrixType 0) for "
			          "contact, non-associated flow, follower loads, corotational "
			          "transforms or LadrunoUP. Reported once.\n";
		}
	}

	// Ladruno ADR-75 P1f: a lookup that finds NO home for an element entry means
	// the sparsity pattern is missing a coefficient the assembly needs, i.e. the
	// contribution is silently DISCARDED and the run solves a matrix that is not
	// the assembled tangent. The stock linear scan simply fell off the end of the
	// row in that case, so it could not be observed at all. Should be
	// unreachable — the pattern comes from the same DOF graph as `id` and P1d
	// made a truncated CSR fill a hard error — but this is the same
	// silent-wrong-answer class as the half-storage asymmetry above, and the
	// branch is FREE: the miss test already exists as `k >= 0`.
	int missing = 0;

	if (fact == 1.0) { // do not need to multiply
		for (int i = 0; i < idSize; i++) {
			int row = id(i);
			if (row < size && row >= 0) {
				int startRowLoc = rowStartA[row] - 1;
				int endRowLoc = rowStartA[row + 1] - 1;
				for (int j = 0; j < idSize; j++) {
					int col = id(j);
					if (col < size && col >= 0 && (!halfStore || col >= row)) {
						// find place in A using colA (Ladruno ADR-75 P1f: was a
						// linear scan of the whole row)
						const int k = ops_pardiso_findCol(colA, startRowLoc,
						                                  endRowLoc, col + 1);
						if (k >= 0)
							A[k] += m(i, j);
						else
							missing = 1;
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
						// find place in A using colA (Ladruno ADR-75 P1f)
						const int k = ops_pardiso_findCol(colA, startRowLoc,
						                                  endRowLoc, col + 1);
						if (k >= 0)
							A[k] += fact * m(i, j);
						else
							missing = 1;
					}
				}  // for j
			}
		}  // for i
	}

	// Ladruno ADR-75 P1f: reported ONCE per SOE — a pattern defect is structural,
	// so repeating it per element would bury the run in output.
	if (missing == 1 && missWarned == 0) {
		missWarned = 1;
		opserr << "WARNING PARDISOGenLinSOE::addA() - an element matrix entry has "
		          "NO slot in the CSR pattern, so its stiffness is being "
		          "DISCARDED.\n"
		          "     The assembled system is not the tangent this model "
		          "defines and the answer may be silently wrong. This indicates "
		          "the sparsity pattern\n"
		          "     and the assembly disagree (a setSize/graph bug, not a "
		          "modelling error). Reported once.\n";
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

	// Ladruno ADR-75 P1d: zeroA runs once per tangent assembly, which makes it
	// the only place this SOE can see a pass boundary. Spend one unit of the
	// addA symmetry-check budget here so the check costs a FIXED number of
	// assemblies regardless of model size or Newton-iteration count.
	if (asymBudget > 0)
		asymBudget--;
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
