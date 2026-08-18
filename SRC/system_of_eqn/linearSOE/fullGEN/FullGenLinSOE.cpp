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
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the implementation for FullGenLinSOE


#include <FullGenLinSOE.h>
#include <stdlib.h>

#include <FullGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
using std::nothrow;

FullGenLinSOE::FullGenLinSOE(FullGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), matA(0),
 Asize(0), Bsize(0), 
 factored(false)
{
    theSolvr.setLinearSOE(*this);
}


FullGenLinSOE::FullGenLinSOE(int N, FullGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), matA(0),
 Asize(0), Bsize(0), 
 factored(false)
{
    size = N;

    A = new (nothrow) double[size*size];
	
    if (A == 0) {
	opserr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	opserr << " ran out of memory for A (size,size) (";
	opserr << size <<", " << size << ") \n";
	size = 0; 
    } else {
	// zero the matrix
	Asize = size*size;
	for (int i=0; i<Asize; i++)
	    A[i] = 0;
    
	B = new (nothrow) double[size];
	X = new (nothrow) double[size];
	
	if (B == 0 || X == 0) {
	    opserr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	} else {
	    Bsize = size;
	    // zero the vectors
	    for (int j=0; j<size; j++) {
		B[j] = 0;
		X[j] = 0;
	    }
	}
    }

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);    
    matA  = new Matrix(A, size, size);

    theSolvr.setLinearSOE(*this);
    
    // invoke setSize() on the Solver        
    if (theSolvr.setSize() < 0) {
	opserr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }    
    
}

    
FullGenLinSOE::~FullGenLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;        
    if (matA != 0) delete matA;        
}


int
FullGenLinSOE::getNumEqn(void) const
{
    return size;
}

int 
FullGenLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();

    if (size*size > Asize) { // we have to get another space for A

	if (A != 0) 
	    delete [] A;

	A = new (nothrow) double[size*size];
	
        if (A == 0) {
            opserr << "WARNING FullGenLinSOE::FullGenLinSOE :";
	    opserr << " ran out of memory for A (size,size) (";
	    opserr << size <<", " << size << ") \n";
	    size = 0; Asize = 0;
	    result =  -1;
        } else
	    Asize = size*size;
    }

    // zero the matrix
    for (int i=0; i<Asize; i++)
	A[i] = 0;
	
    factored = false;
    
    if (size > Bsize) { // we have to get space for the vectors
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;

	// create the new
	B = new (nothrow) double[size];
	X = new (nothrow) double[size];
	
        if (B == 0 || X == 0) {
            opserr << "WARNING FullGenLinSOE::FullGenLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	    result =  -1;
        }
	else
	    Bsize = size;
    }

    // zero the vectors
    for (int j=0; j<Bsize; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    // create new Vectors
    // Ladruno: also create when the wrappers are still null. A model whose
    // every DOF is fixed or sp-constrained numbers ZERO equations, so
    // size == oldSize == 0 here and the original condition left vectX/vectB/
    // matA null from the default constructor -- the first getB()/getX() then
    // hit the FATAL exit(-1) path and killed the whole process (silent under
    // the Python module's stream redirection; looked like a native crash).
    if (size != oldSize || vectX == 0) {
	if (vectX != 0)
	    delete vectX;

	if (vectB != 0)
	    delete vectB;

	if (matA != 0)
	    delete matA;
	
	vectX = new Vector(X,Bsize);
	vectB = new Vector(B,Bsize);	
	matA = new Matrix(A,Bsize, Bsize);	
    }

    // invoke setSize() on the Solver    
    LinearSOESolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:FullGenLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    
    
    return result;
}

int 
FullGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "FullGenLinSOE::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i);
		     }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i) * fact;
		     }
		}  // for j
	    } 
	}  // for i
    }    
    return 0;
}



int 
FullGenLinSOE::addColA(const Vector &colData, int col, double fact)
{
  
  if (fact == 0.0)  return 0;
  
  if (colData.Size() != size) {
    opserr << "FullGenLinSOE::addColA() - colData size not equal to n\n";
    return -1;
  }
  
  if (col > size && col < 0) {
    opserr << "FullGenLinSOE::addColA() - col " << col << "outside range 0 to " << size << endln;
    return -1;
  }
  
  
  if (fact == 1.0) { // do not need to multiply 

    double *coliPtr = A + col*size;
    for (int row=0; row<size; row++) {
      *coliPtr += colData(row);
      coliPtr++;
    }

  } else {

    double *coliPtr = A + col*size;
    for (int row=0; row<size; row++) {
      *coliPtr += colData(row) * fact;
      coliPtr++;
    }

  }

  return 0;
}




int 
FullGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "FullGenLinSOE::addB()	- Vector and ID not of similar sizes\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}



int
FullGenLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


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
FullGenLinSOE::zeroA(void)
{
    double *Aptr = A;
    int theSize = size*size;
    for (int i=0; i<theSize; i++)
	*Aptr++ = 0;

    factored = false;
}
	
void 
FullGenLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

int
FullGenLinSOE::formAp(const Vector &p, Vector &Ap)
{
  // Check that p and A are same size
  if (size != p.Size() || size != Ap.Size() || p.Size() != Ap.Size()) {
    opserr << "FullGenLinSOE::formAp -- vectors not of same size\n";
    return -1;
  }

  for (int row = 0; row < size; row++) {
    double sum = 0.0;
    double *APtr = A + row;
    for (int col = 0; col < size; col++) {
      APtr += size;
      sum += (*APtr) * p(col);
    }
    Ap(row) = sum;
  }

  return 0;
}

void 
FullGenLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
FullGenLinSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

// Ladruno: the three accessors below used to `exit(-1)` when their wrapper was
// still null. That is a state a USER can reach from the interpreter -- the
// wrappers are created in setSize(), which the analysis only reaches once
// ConstraintHandler::handle() has succeeded -- so any diagnostic that reads the
// SOE before a successful analyze() KILLED THE PROCESS. From Python that is a
// silent death: exit() unwinds nothing, so there is no traceback and no
// exception, just a dead interpreter (the `printA` report, 2026-08-18).
//
// It is `constraints LadrunoContact` that makes this reachable in practice.
// Under the upstream handlers handle() essentially never fails, but the contact
// handler's whole ADR-78/ADR-85 abort discipline is built on returning -1 from
// handle() for a refused contact -- domainChanged() then fails, setSize() is
// never called, and the very next `printA` (the standard tangent-extraction
// probe, and exactly what a user reaches for to debug the refusal) took the
// interpreter down with it.
//
// A null wrapper now reports itself and returns an EMPTY result. getA() returns
// 0, which is what LinearSOE::getA() itself returns by default and which every
// caller must therefore already handle (LinearSOE::saveSparseA and OPS_printA
// both already branch on it). getX()/getB() return a shared empty Vector: size 0
// is the honest description of an SOE that has never been sized, and OPS_printB
// already has a size == 0 branch. Callers inside the analysis flow are unchanged
// -- they only ever run after a successful domainChanged(), where setSize() has
// allocated the wrappers.
//
// Function-local static (not file-scope) so it is initialized on first use --
// no static-initialization-order dependency on Vector's own machinery.
static const Vector &
ladrunoEmptyVector(void)
{
    static Vector theEmptyVector;   // default ctor => Size() == 0
    return theEmptyVector;
}

const Vector &
FullGenLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "WARNING FullGenLinSOE::getX - the SOE has not been sized yet "
	          "(setSize() has not run, so there is no solution vector); returning "
	          "an empty Vector. Run a successful analyze() first.\n";
	return ladrunoEmptyVector();
    }
    return *vectX;
}

const Vector &
FullGenLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "WARNING FullGenLinSOE::getB - the SOE has not been sized yet "
	          "(setSize() has not run, so there is no right-hand side); returning "
	          "an empty Vector. Run a successful analyze() first.\n";
	return ladrunoEmptyVector();
    }
    return *vectB;
}

const Matrix *
FullGenLinSOE::getA(void)
{
    if (matA == 0) {
	opserr << "WARNING FullGenLinSOE::getA - the SOE has not been sized yet "
	          "(setSize() has not run, so there is no system matrix); returning 0. "
	          "Run a successful analyze() first -- if analyze() reported a failure "
	          "(e.g. a refused contact aborting ConstraintHandler::handle()), fix "
	          "that first: the tangent for this model was never assembled.\n";
	return 0;
    }
    return matA;
}

double 
FullGenLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
FullGenLinSOE::setFullGenSolver(FullGenLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:FullGenLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
FullGenLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int 
FullGenLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}






