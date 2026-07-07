/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno: FEAST band-targeted eigensolver SOE (ADR 43, P1). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md.

#include <FeastEigenSOE.h>
#include <FeastEigenSolver.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>

#include <vector>
#include <algorithm>

FeastEigenSOE::FeastEigenSOE(FeastEigenSolver &theSolvr)
:EigenSOE(theSolvr, EigenSOE_TAGS_FeastEigenSOE),
 size(0), nnz(0), csrCapacity(0),
 rowPtr(0), colInd(0), Kvals(0), Mvals(0),
 emin(0.0), emax(0.0), m0(0), nq(8), maxRefine(20), tolExp(12), verbose(false),
 certifyFlag(false)
{
  theSolvr.setEigenSOE(*this);
}

FeastEigenSOE::~FeastEigenSOE()
{
  if (rowPtr != 0) delete [] rowPtr;
  if (colInd != 0) delete [] colInd;
  if (Kvals  != 0) delete [] Kvals;
  if (Mvals  != 0) delete [] Mvals;
}

int
FeastEigenSOE::getNumEqn(void) const
{
  return size;
}

int
FeastEigenSOE::setSize(Graph &theGraph)
{
  int result = 0;
  size = theGraph.getNumVertex();

  // --- pass 1: upper-bound the entry count (adjacency + the diagonal) ------
  // OpenSees Graph vertices do NOT list themselves in their adjacency, so the
  // diagonal is added explicitly below; MKL FEAST requires the diagonal to be
  // present in the pattern for every row.
  int upperNNZ = 0;
  Vertex *vertexPtr;
  VertexIter &theVertices = theGraph.getVertices();
  while ((vertexPtr = theVertices()) != 0)
    upperNNZ += vertexPtr->getAdjacency().Size() + 1;

  // (re)allocate the row-pointer array (n+1) whenever the order changes
  if (rowPtr != 0) delete [] rowPtr;
  rowPtr = new int[size + 1];
  if (rowPtr == 0) {
    opserr << "FeastEigenSOE::setSize() -- ran out of memory for rowPtr, n = " << size << endln;
    size = 0; nnz = 0;
    return -1;
  }

  // (re)allocate the column/value arrays if we need more room
  if (upperNNZ > csrCapacity) {
    if (colInd != 0) delete [] colInd;
    if (Kvals  != 0) delete [] Kvals;
    if (Mvals  != 0) delete [] Mvals;
    colInd = new int[upperNNZ];
    Kvals  = new double[upperNNZ];
    Mvals  = new double[upperNNZ];
    if (colInd == 0 || Kvals == 0 || Mvals == 0) {
      opserr << "FeastEigenSOE::setSize() -- ran out of memory for CSR arrays, nnz = " << upperNNZ << endln;
      csrCapacity = 0; size = 0; nnz = 0;
      return -1;
    }
    csrCapacity = upperNNZ;
  }

  // --- pass 2: build the shared CSR pattern (1-based, ascending columns) ---
  // Each row's column list = {diagonal} U adjacency, de-duplicated & sorted.
  int loc = 0;                 // 0-based running write cursor into colInd
  rowPtr[0] = 1;               // MKL uses 1-based row starts
  std::vector<int> cols;
  for (int a = 0; a < size; a++) {
    vertexPtr = theGraph.getVertexPtr(a);
    if (vertexPtr == 0) {
      opserr << "FeastEigenSOE::setSize() -- vertex " << a << " not in graph\n";
      size = 0; nnz = 0;
      return -1;
    }

    cols.clear();
    cols.push_back(a);                                   // the diagonal
    const ID &theAdjacency = vertexPtr->getAdjacency();
    for (int i = 0; i < theAdjacency.Size(); i++)
      cols.push_back(theAdjacency(i));

    std::sort(cols.begin(), cols.end());
    cols.erase(std::unique(cols.begin(), cols.end()), cols.end());

    for (size_t k = 0; k < cols.size(); k++)
      colInd[loc++] = cols[k] + 1;                       // store 1-based
    rowPtr[a + 1] = loc + 1;                             // 1-based
  }
  nnz = loc;

  // zero the value arrays
  for (int i = 0; i < nnz; i++) {
    Kvals[i] = 0.0;
    Mvals[i] = 0.0;
  }

  // invoke setSize() on the solver (sizes its per-mode workspace)
  EigenSolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "FeastEigenSOE::setSize() -- solver failed in setSize()\n";
    return solverOK;
  }

  return result;
}

int
FeastEigenSOE::findSlot(int row, int col) const
{
  int lo = rowPtr[row] - 1;         // 0-based start of the row
  int hi = rowPtr[row + 1] - 1;     // 0-based end (exclusive)
  int target = col + 1;             // colInd is 1-based
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    int c = colInd[mid];
    if (c < target)      lo = mid + 1;
    else if (c > target) hi = mid;
    else                 return mid;
  }
  return -1;
}

int
FeastEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
  if (fact == 0.0)
    return 0;

  int idSize = id.Size();
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "FeastEigenSOE::addA() -- Matrix and ID not of similar sizes\n";
    return -1;
  }

  for (int i = 0; i < idSize; i++) {
    int row = id(i);
    if (row < 0 || row >= size)
      continue;
    for (int j = 0; j < idSize; j++) {
      int col = id(j);
      if (col < 0 || col >= size)
        continue;
      int slot = findSlot(row, col);
      if (slot >= 0)
        Kvals[slot] += (fact == 1.0) ? m(i, j) : fact * m(i, j);
      else
        opserr << "FeastEigenSOE::addA() -- entry (" << row << "," << col
               << ") not in CSR pattern; dropped\n";
    }
  }

  return 0;
}

int
FeastEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
  if (fact == 0.0)
    return 0;

  int idSize = id.Size();
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "FeastEigenSOE::addM() -- Matrix and ID not of similar sizes\n";
    return -1;
  }

  for (int i = 0; i < idSize; i++) {
    int row = id(i);
    if (row < 0 || row >= size)
      continue;
    for (int j = 0; j < idSize; j++) {
      int col = id(j);
      if (col < 0 || col >= size)
        continue;
      int slot = findSlot(row, col);
      if (slot >= 0)
        Mvals[slot] += (fact == 1.0) ? m(i, j) : fact * m(i, j);
      else
        opserr << "FeastEigenSOE::addM() -- entry (" << row << "," << col
               << ") not in CSR pattern; dropped\n";
    }
  }

  return 0;
}

void
FeastEigenSOE::zeroA(void)
{
  for (int i = 0; i < nnz; i++)
    Kvals[i] = 0.0;
}

void
FeastEigenSOE::zeroM(void)
{
  for (int i = 0; i < nnz; i++)
    Mvals[i] = 0.0;
}

void FeastEigenSOE::setBand(double lo, double hi)   { emin = lo; emax = hi; }
void FeastEigenSOE::setSubspaceSize(int m0_)        { m0 = m0_; }
void FeastEigenSOE::setNumQuad(int nq_)             { nq = nq_; }
void FeastEigenSOE::setMaxRefine(int maxRefine_)    { maxRefine = maxRefine_; }
void FeastEigenSOE::setTol(double eps)             { tolExp = (int)eps; }
void FeastEigenSOE::setVerbose(bool v)             { verbose = v; }
void FeastEigenSOE::setCertify(bool c)             { certifyFlag = c; }

int
FeastEigenSOE::getNumFoundModes(void) const
{
  const FeastEigenSolver *slv =
      static_cast<const FeastEigenSolver *>(theSolver);
  return (slv == 0) ? 0 : slv->getNumFoundModes();
}

double FeastEigenSOE::getEmin(void) const          { return emin; }
double FeastEigenSOE::getEmax(void) const          { return emax; }
int    FeastEigenSOE::getSubspaceSize(void) const  { return m0; }
int    FeastEigenSOE::getNumQuad(void) const       { return nq; }
int    FeastEigenSOE::getMaxRefine(void) const     { return maxRefine; }
int    FeastEigenSOE::getTol(void) const           { return tolExp; }
bool   FeastEigenSOE::getVerbose(void) const       { return verbose; }
bool   FeastEigenSOE::getCertify(void) const       { return certifyFlag; }

int
FeastEigenSOE::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "FeastEigenSOE::sendSelf() -- serial only (ADR 43 P1); "
         << "distributed assembly arrives at P3\n";
  return -1;
}

int
FeastEigenSOE::recvSelf(int commitTag, Channel &theChannel,
                        FEM_ObjectBroker &theBroker)
{
  opserr << "FeastEigenSOE::recvSelf() -- serial only (ADR 43 P1); "
         << "distributed assembly arrives at P3\n";
  return -1;
}
