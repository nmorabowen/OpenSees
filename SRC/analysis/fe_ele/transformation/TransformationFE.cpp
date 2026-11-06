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
                                                                        
// $Revision: 1.20 $
// $Date: 2008-11-19 23:39:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/transformation/TransformationFE.cpp,v $
                                                                        
// Written: fmk 
// Created: 05/99
//
// Purpose: This file contains the code for implementing the methods
// of the TransformationFE class interface.

#include "TransformationFE.h"
#include <LadrunoThreads.h>   // Ladruno WP-192
#include <TransformationDOF_Group.h>
#include <vector>
#include <algorithm> // Ladruno WP-216: std::max
#include <stdlib.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <TransformationConstraintHandler.h>

#define MAX_NUM_DOF 64

// static variables initialisation
Matrix **TransformationFE::modMatrices; 
Vector **TransformationFE::modVectors;  
Matrix **TransformationFE::theTransformations; 
int TransformationFE::numTransFE(0);           
int TransformationFE::transCounter(0);           
int TransformationFE::sizeTransformations(0);          
double *TransformationFE::dataBuffer = 0;          
double *TransformationFE::localKbuffer = 0;          
int    *TransformationFE::dofData = 0;    ;          
int TransformationFE::sizeBuffer(0);            
// Ladruno WP-216: entries in dofData. It holds one count PER NODE
// (numDOFs.setData(dofData, numGroups)), so it must hold the largest element's
// node count, not MAX_NUM_DOF; file-static to leave the header alone.
static int ladrunoSizeDofData = 0;

//  TransformationFE(Element *, Integrator *theIntegrator);
//	construictor that take the corresponding model element.
TransformationFE::TransformationFE(int tag, Element *ele)
:FE_Element(tag, ele), theDOFs(0), numSPs(0), theSPs(0), modID(0), 
  modTangent(0), modResidual(0), numGroups(0), numTransformedDOF(0)
{
  // set number of original dof at ele
    numOriginalDOF = ele->getNumDOF();

    // create the array of pointers to DOF_Groups
    const ID &nodes = ele->getExternalNodes();
    Domain *theDomain = ele->getDomain();
    int numNodes = nodes.Size();
    theDOFs = new DOF_Group *[numNodes];
    if (theDOFs == 0) {
	opserr << "FATAL TransformationFE::TransformationFE() - out of memory creating ";
	opserr << "array of size : " << numNodes << " for storage of DOF_Group\n";
	exit(-1);
    }

    numGroups = numNodes;

    // now fill the array of DOF_Group pointers
    for (int i=0; i<numNodes; i++) {
	Node *theNode = theDomain->getNode(nodes(i));
	if (theNode == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - no Node with tag: ";
	    opserr << nodes(i) << " in the domain\n";;
	    exit(-1);
	}
	DOF_Group *theDofGroup = theNode->getDOF_GroupPtr();
	if (theDofGroup == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - no DOF_Group : ";
	    opserr << " associated with node: " << nodes(i) << " in the domain\n";;
	    exit(-1);
	}	
	theDOFs[i] = theDofGroup;
    }

    // see if theTransformation array is big enough
    // if not delete the old and create a new one
    if (numNodes > sizeTransformations) {
	if (theTransformations != 0) 
	    delete [] theTransformations;
	
	theTransformations = new Matrix *[numNodes];
	if (theTransformations == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - out of memory ";
	    opserr << "for array of pointers for Transformation matrices of size ";
	    opserr << numNodes;
	    exit(-1);
	}		    
	sizeTransformations = numNodes;
    }	

    // if this is the first element of this type create the arrays for 
    // modified tangent and residual matrices
    if (numTransFE == 0) {

	modMatrices = new Matrix *[MAX_NUM_DOF+1];
	modVectors  = new Vector *[MAX_NUM_DOF+1];
	dataBuffer = new double[MAX_NUM_DOF*MAX_NUM_DOF];
	localKbuffer = new double[MAX_NUM_DOF*MAX_NUM_DOF];
	dofData      = new int[MAX_NUM_DOF];
	ladrunoSizeDofData = MAX_NUM_DOF; // Ladruno WP-216
	sizeBuffer = MAX_NUM_DOF*MAX_NUM_DOF;
	
	if (modMatrices == 0 || modVectors == 0 || dataBuffer == 0 ||
	    localKbuffer == 0 || dofData == 0) {
	    opserr << "TransformationFE::TransformationFE(Element *) ";
	    opserr << " ran out of memory";
	}
	// Ladruno WP-149: <= not <. MAX_NUM_DOF+1 slots, and setID() takes the
	// pooled branch for numTransformedDOF == MAX_NUM_DOF, so slot [MAX_NUM_DOF]
	// was read uninitialized (same off-by-one as FE_Element.cpp).
	for (int i=0; i<=MAX_NUM_DOF; i++) {
	    modMatrices[i] = 0;
	    modVectors[i] = 0;
	}
    }

    // Ladruno WP-216: an element with more than MAX_NUM_DOF nodes (a
    // LadrunoMembrane floor panel) wrote past dofData on every
    // getTangent / getResidual / get*_Force: heap corruption under
    // constraints Transformation. Grow it to the node count.
    if (numNodes > ladrunoSizeDofData) {
	delete [] dofData;
	dofData = new int[numNodes];
	ladrunoSizeDofData = numNodes;
    }

    // increment the number of transformations
    numTransFE++;
}



// ~TransformationFE();    
//	destructor.

TransformationFE::~TransformationFE()
{

    numTransFE--;
    
    if (theDOFs != 0)
	delete [] theDOFs;
    if (theSPs != 0)
	delete [] theSPs;

    int numDOF = 0;    
    if (modID != 0)
	numDOF = modID->Size();
    
    if (modID != 0)
	delete modID;

    if (numDOF > MAX_NUM_DOF) {
	// tangent and residual may have been  created specially
	if (modTangent != 0) delete modTangent;
	if (modResidual != 0) delete modResidual;
    }

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numTransFE == 0) {
	for (int i=0; i<=MAX_NUM_DOF; i++) {   // Ladruno WP-149: <= -- slot [MAX_NUM_DOF] leaked
	    if (modVectors[i] != 0)
		delete modVectors[i];
	    if (modMatrices[i] != 0)
		delete modMatrices[i];
	}
	delete [] modMatrices;
	delete [] modVectors;
	delete [] theTransformations;
	delete [] dataBuffer;
	delete [] localKbuffer;
	delete [] dofData;
	modMatrices = 0;
	modVectors = 0;
	theTransformations = 0;
	dataBuffer = 0;
	localKbuffer = 0;
	dofData = 0;
	ladrunoSizeDofData = 0; // Ladruno WP-216
	sizeTransformations = 0;
	sizeBuffer = 0;
	transCounter = 0;
    }
}    


const ID &
TransformationFE::getDOFtags(void) const 
{
    return this->FE_Element::getDOFtags();
}


const ID &
TransformationFE::getID(void) const
{
    // make sure that it exists
    if (modID == 0) {
	opserr << "FATAL TransformationFE::getID() called before setID()\n";
	exit(-1);
    }
    return *modID;
}


int
TransformationFE::setID(void)
{
    // determine number of DOF
    numTransformedDOF = 0;
    for (int ii=0; ii<numGroups; ii++) {
	DOF_Group *dofPtr = theDOFs[ii];
	numTransformedDOF += dofPtr->getNumDOF();
    }

    // create an ID to hold the array, cannot use existing as 
    // may be different size
    if (modID != 0)
      delete modID;
    modID = 0;

    modID = new ID(numTransformedDOF);
    if (modID == 0 || modID->Size() == 0) {
	opserr << "TransformationFE::setID() ";
	opserr << " ran out of memory for ID of size :";
	opserr << numTransformedDOF << endln;
	exit(-1);
    }

    // fill in the ID
    int current = 0;
    for (int i=0; i<numGroups; i++) {
	DOF_Group *dofPtr = theDOFs[i];
	const ID &theDOFid = dofPtr->getID();

	for (int j=0; j<theDOFid.Size(); j++)  
	    if (current < numTransformedDOF)
		(*modID)(current++) = theDOFid(j);
	    else {
		opserr << "WARNING TransformationFE::setID() - numDOF and";
		opserr << " number of dof at the DOF_Groups\n";
		return -3;
	    }		
    }
    
    // set the pointers to the modified tangent matrix and residual vector
    if (numTransformedDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (modVectors[numTransformedDOF] == 0) {
	    modVectors[numTransformedDOF] = new Vector(numTransformedDOF);
	    modMatrices[numTransformedDOF] = new Matrix(numTransformedDOF,numTransformedDOF);
	    modResidual = modVectors[numTransformedDOF];
	    modTangent = modMatrices[numTransformedDOF];
	    if (modResidual == 0 || modResidual->Size() != numTransformedDOF ||	
		modTangent == 0 || modTangent->noCols() != numTransformedDOF)	{  
		opserr << "TransformationFE::setID() ";
		opserr << " ran out of memory for vector/Matrix of size :";
		opserr << numTransformedDOF << endln;
		exit(-1);
	    }
	} else {
	    modResidual = modVectors[numTransformedDOF];
	    modTangent = modMatrices[numTransformedDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	modResidual = new Vector(numTransformedDOF);
	modTangent = new Matrix(numTransformedDOF, numTransformedDOF);
	if (modResidual == 0 || modResidual->Size() ==0 ||
	    modTangent ==0 || modTangent->noRows() ==0) {
	    
	    opserr << "TransformationFE::setID() ";
	    opserr << " ran out of memory for vector/Matrix of size :";
	    opserr << numTransformedDOF << endln;
	    exit(-1);
	}
    }     

    return 0;
}


// ---------------------------------------------------------------------------
// Ladruno WP-192: the per-thread twin of TransformationFE's class-wide scratch.
// getTangent()/getResidual() build T^t K T / T^t R in buffers that ONE object
// owns for the whole class (ADR-75b 5.4-H1, "one level up" under constraints
// Transformation). A thread inside a WP-192 loop gets its own set, sized like
// the originals (MAX_NUM_DOF), plus its own modTangent/modResidual for this
// FE's transformed size when the FE uses the class-wide pool. FEs larger than
// MAX_NUM_DOF own their modTangent/modResidual and keep them.
// ---------------------------------------------------------------------------
namespace {
  struct LadrunoTFEScratch {
    std::vector<int> dofs;
    std::vector<Matrix *> trans;
    std::vector<double> localK, data;
    std::vector<Matrix *> m;
    std::vector<Vector *> v;
    ~LadrunoTFEScratch() {
      for (size_t i = 0; i < m.size(); i++) delete m[i];
      for (size_t i = 0; i < v.size(); i++) delete v[i];
    }
  };
  thread_local LadrunoTFEScratch ladrunoTFE;
}

void
TransformationFE::ladrunoScratch(int *&dofs, Matrix **&trans, double *&localK,
                                 double *&data, Matrix **modT, Vector **modR) const
{
    LadrunoTFEScratch &t = ladrunoTFE;
    // Ladruno WP-216: one entry per node, like dofData
    const int nDofs = std::max(MAX_NUM_DOF, sizeTransformations);
    if ((int)t.dofs.size() < nDofs)
      t.dofs.resize(nDofs);
    if ((int)t.trans.size() < sizeTransformations)
      t.trans.resize(sizeTransformations, (Matrix *)0);
    if ((int)t.localK.size() < sizeBuffer)
      t.localK.resize(sizeBuffer);
    if ((int)t.data.size() < sizeBuffer)
      t.data.resize(sizeBuffer);
    dofs = t.dofs.data();
    trans = t.trans.data();
    localK = t.localK.data();
    data = t.data.data();

    const int n = numTransformedDOF;
    const bool pooled = (n <= MAX_NUM_DOF);
    if (modT != 0 && pooled) {
      if ((int)t.m.size() <= n)
        t.m.resize(n + 1, (Matrix *)0);
      if (t.m[n] == 0)
        t.m[n] = new Matrix(n, n);
      *modT = t.m[n];
    }
    if (modR != 0 && pooled) {
      if ((int)t.v.size() <= n)
        t.v.resize(n + 1, (Vector *)0);
      if (t.v[n] == 0)
        t.v[n] = new Vector(n);
      *modR = t.v[n];
    }
}


// Ladruno WP-192: a TransformationFE is threadable when its element is (base
// check) and none of its nodes has a TIME-VARYING MP: TransformationDOF_Group::
// getT() rebuilds such a T into the group's own Trans on every call -- a write
// shared by every element on that node. Time-invariant T (rigidDiaphragm,
// equalDOF, rigidLink, SP-only) is read-only and fine.
bool
TransformationFE::ladrunoThreadSafeFE(void)
{
    for (int i = 0; i < numGroups; i++) {
      TransformationDOF_Group *g = dynamic_cast<TransformationDOF_Group *>(theDOFs[i]);
      if (g != 0 && g->ladrunoTimeVaryingT())
        return false;
    }
    return this->FE_Element::ladrunoThreadSafeFE();
}

const Matrix &
TransformationFE::getTangent(Integrator *theNewIntegrator)
{
    const Matrix &theTangent = this->FE_Element::getTangent(theNewIntegrator);

    // Ladruno WP-192: inside a threaded loop every class-wide buffer below
    // (dofData, theTransformations, localKbuffer, dataBuffer, the modTangent
    // pool) is replaced by the calling thread's own copy; outside one, the
    // originals are used and this is the vanilla code path.
    int *dofDataX = dofData;
    Matrix **theTransformationsX = theTransformations;
    double *localKbufferX = localKbuffer;
    double *dataBufferX = dataBuffer;
    Matrix *modTangentX = modTangent;
    if (ladruno_inThreadedRegion())
      ladrunoScratch(dofDataX, theTransformationsX, localKbufferX, dataBufferX,
                     &modTangentX, 0);

    static ID numDOFs(dofData, 1);
    ID ladrunoNumDOFs(dofDataX, numGroups);              // Ladruno WP-192
    ID &numDOFsX = ladruno_inThreadedRegion() ? ladrunoNumDOFs : numDOFs;
    if (!ladruno_inThreadedRegion())
      numDOFs.setData(dofData, numGroups);
    
    // DO THE SP STUFF TO THE TANGENT 
    
    // get the transformation matrix from each dof group & number of local dof
    // for original node.
    int numNode = numGroups;
    for (int a = 0; a<numNode; a++) {
      Matrix *theT = theDOFs[a]->getT();
      theTransformationsX[a] = theT;
      if (theT != 0)
	numDOFsX[a] = theT->noRows(); // T^ 
      else
	numDOFsX[a] = theDOFs[a]->getNumDOF();
    }

    // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
    // where blocks are of size equal to num ele dof at a node

    int startRow = 0;
    int noRowsTransformed = 0;
    int noRowsOriginal = 0;

    static Matrix localKS;                            // Ladruno WP-192: serial path
    Matrix localKT;                                   // Ladruno WP-192: threaded path
    Matrix &localK = ladruno_inThreadedRegion() ? localKT : localKS;

    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {

	int startCol = 0;
	int numDOFi = numDOFsX[i];	
	int noColsOriginal = 0;

	for (int j=0; j<numNode; j++) {

	    const Matrix *Ti = theTransformationsX[i];
	    const Matrix *Tj = theTransformationsX[j];
	    int numDOFj = numDOFsX[j];	
	    localK.setData(localKbufferX, numDOFi, numDOFj);

	    // copy K(i,j) into localK matrix
	    // CHECK SIZE OF BUFFER	    
	    for (int a=0; a<numDOFi; a++)
		for (int b=0; b<numDOFj; b++)
		    localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);

	    // now perform the matrix computation T(i)^T localK T(j)
	    // note: if T == 0 then the Identity is assumed
	    int noColsTransformed = 0;
	    static Matrix localTtKTS;                  // Ladruno WP-192: serial path
	    Matrix localTtKTT;                         // Ladruno WP-192: threaded path
	    Matrix &localTtKT = ladruno_inThreadedRegion() ? localTtKTT : localTtKTS;
	    
	    if (Ti != 0 && Tj != 0) {
		noRowsTransformed = Ti->noCols();
		noColsTransformed = Tj->noCols();
		// CHECK SIZE OF BUFFER
		localTtKT.setData(dataBufferX, noRowsTransformed, noColsTransformed);
		//localTtKT = (*Ti) ^ localK * (*Tj);
		localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
	    } else if (Ti == 0 && Tj != 0) {
		noRowsTransformed = numDOFi;
		noColsTransformed = Tj->noCols();
		// CHECK SIZE OF BUFFER
		localTtKT.setData(dataBufferX, noRowsTransformed, noColsTransformed);
		// localTtKT = localK * (*Tj);	       
		localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
	    } else if (Ti != 0 && Tj == 0) {
		noRowsTransformed = Ti->noCols();
		noColsTransformed = numDOFj;
		localTtKT.setData(dataBufferX, noRowsTransformed, noColsTransformed);
		//localTtKT = (*Ti) ^ localK;
		localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
	    } else {
		noRowsTransformed = numDOFi;
		noColsTransformed = numDOFj;
		localTtKT.setData(dataBufferX, noRowsTransformed, noColsTransformed);
		localTtKT = localK;
	    }
	    // now copy into modTangent the T(i)^t K(i,j) T(j) product
	    for (int c=0; c<noRowsTransformed; c++) 
		for (int d=0; d<noColsTransformed; d++) 
		    (*modTangentX)(startRow+c, startCol+d) = localTtKT(c,d);
	    
	    startCol += noColsTransformed;
	    noColsOriginal += numDOFj;
	}

	noRowsOriginal += numDOFi;
	startRow += noRowsTransformed;
    }

    return *modTangentX;   // Ladruno WP-192
}


const Vector &
TransformationFE::getResidual(Integrator *theNewIntegrator)

{
    const Vector &theResidual = this->FE_Element::getResidual(theNewIntegrator);

    // Ladruno WP-192: the modResidual pool is class-wide; per-thread in a loop.
    Vector *modResidualX = modResidual;
    if (ladruno_inThreadedRegion()) {
      int *d0; Matrix **t0; double *k0; double *b0;
      ladrunoScratch(d0, t0, k0, b0, 0, &modResidualX);
    }
    // DO THE SP STUFF TO THE TANGENT
    
    // perform Tt R  -- as T is block diagonal do T(i)^T R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRowTransformed = 0;
    int startRowOriginal = 0;
    int numNode = numGroups;

    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {
	int noRows = 0;
	int noCols = 0;
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	  noRows = Ti->noCols(); // T^
	  noCols = Ti->noRows();

	  /*
	  Vector orig(noCols);
	  Vector mod(noRows);
	  for (int k=startRowOriginal; k<startRowOriginal+noCols; k++)
	    orig(k-startRowOriginal)= theResidual(k);
	  mod = (*Ti)^orig;
	  for (int k=startRowTransformed; k<startRowTransformed+noRows; k++)
	    (*modResidualX)(k) = mod (k-startRowTransformed);

	  */

	  for (int j=0; j<noRows; j++) {
	    double sum = 0.0;
	    for (int k=0; k<noCols; k++)
	      sum += (*Ti)(k,j) * theResidual(startRowOriginal + k);
	    (*modResidualX)(startRowTransformed +j) = sum;
	  }

	} else {
	  noCols = theDOFs[i]->getNumDOF();
	  noRows = noCols;
	  for (int j=0; j<noRows; j++)
	    (*modResidualX)(startRowTransformed +j) = theResidual(startRowOriginal + j);
	}
	startRowTransformed += noRows;
	startRowOriginal += noCols;
    }

    return *modResidualX;   // Ladruno WP-192
}




// Ladruno (ADR-80 P3) -- the Kratos-style `b -= K * du_D` term.
//
// WHY THIS AND NOT getTangForce(). The stub below has to stay a stub: its only
// input is a GLOBAL vector indexed by myID, and a Transformation-eliminated dof
// has myID == -1, so the prescribed increment it would need is not expressible
// in its argument. The increment lives on the SP_Constraint, reachable only
// through the constrained node's TransformationDOF_Group.
//
// WHAT IT COMPUTES. du in the element's ORIGINAL dof space, zero everywhere
// except the sp'd dofs where it holds (prescribed total - committed nodal
// value); then f = K du with K the element tangent AT THE CURRENT (committed)
// state; then T^t f, so the caller can assemble it with the element's ordinary
// ID -- whose -1 entries drop the eliminated rows for free, leaving exactly the
// free-row coupling force.
//
// Returns 0 (and assembles nothing) when the element touches no sp'd dof, so
// the caller's loop over every FE_Element costs one virtual call per node on
// the overwhelming majority of elements.
const Vector *
TransformationFE::getSPTangentForce(Integrator *theNewIntegrator)
{
    if (modResidual == 0 || numOriginalDOF <= 0)
        return 0;

    // ---- 1. du in ORIGINAL element dof space -----------------------------
    Vector du(numOriginalDOF);
    du.Zero();

    int numSPdof = 0;
    int startOriginal = 0;
    for (int i = 0; i < numGroups; i++) {
        const Matrix *Ti = theDOFs[i]->getT();
        int noCols = (Ti != 0) ? Ti->noRows() : theDOFs[i]->getNumDOF();
        numSPdof += theDOFs[i]->getSPDispIncr(du, startOriginal);
        startOriginal += noCols;
    }

    if (numSPdof == 0)
        return 0;

    // A HOMOGENEOUS sp (plain `fix`) has du == 0 and contributes nothing but a
    // wasted element-tangent formation -- and most boundary elements in any
    // model are of exactly that kind. Reporting "no contribution" for them also
    // lets the caller's fallback recognise a model that has no non-homogeneous
    // sp at all, which is the case the whole route is inert for.
    bool nonzero = false;
    for (int i = 0; i < numOriginalDOF && !nonzero; i++)
        if (du(i) != 0.0)
            nonzero = true;

    if (!nonzero)
        return 0;

    // ---- 2. f = K du, K at the element's current (committed) state --------
    const Matrix &theTangent = this->FE_Element::getTangent(theNewIntegrator);

    Vector f(numOriginalDOF);
    if (f.addMatrixVector(0.0, theTangent, du, 1.0) < 0) {
        opserr << "WARNING TransformationFE::getSPTangentForce() - "
               << "addMatrixVector failed\n";
        return 0;
    }

    // ---- 3. modResidual = T^t f  (block loop mirrors getResidual) ---------
    modResidual->Zero();

    int startRowTransformed = 0;
    int startRowOriginal = 0;
    for (int i = 0; i < numGroups; i++) {
        int noRows = 0;
        int noCols = 0;
        const Matrix *Ti = theDOFs[i]->getT();
        if (Ti != 0) {
            noRows = Ti->noCols(); // T^
            noCols = Ti->noRows();
            for (int j = 0; j < noRows; j++) {
                double sum = 0.0;
                for (int k = 0; k < noCols; k++)
                    sum += (*Ti)(k, j) * f(startRowOriginal + k);
                (*modResidual)(startRowTransformed + j) = sum;
            }
        } else {
            noCols = theDOFs[i]->getNumDOF();
            noRows = noCols;
            for (int j = 0; j < noRows; j++)
                (*modResidual)(startRowTransformed + j) = f(startRowOriginal + j);
        }
        startRowTransformed += noRows;
        startRowOriginal += noCols;
    }

    return modResidual;
}

const Vector &
TransformationFE::getTangForce(const Vector &disp, double fact)
{
    opserr << "TransformationFE::getTangForce() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}

const Vector &
TransformationFE::getK_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addKtToTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}


const Vector &
TransformationFE::getKi_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addKiToTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}

const Vector &
TransformationFE::getM_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addMtoTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}

const Vector &
TransformationFE::getC_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addCtoTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}


void  
TransformationFE::addD_Force(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalD_Force(response, fact);
}   	 

void  
TransformationFE::addM_Force(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalM_Force(response, fact);
}   	 



// CHANGE THE ID SENT
const Vector &
TransformationFE::getLastResponse(void)
{
    Integrator *theLastIntegrator = this->getLastIntegrator();
    if (theLastIntegrator != 0) {
	if (theLastIntegrator->getLastResponse(*modResidual,*modID) < 0) {
	    opserr << "WARNING TransformationFE::getLastResponse(void)";
	    opserr << " - the Integrator had problems with getLastResponse()\n";
	}
    }
    else {
	modResidual->Zero();
	opserr << "WARNING  TransformationFE::getLastResponse()";
	opserr << " No Integrator yet passed\n";
    }
    
    Vector &result = *modResidual;
    return result;
}


int 
TransformationFE::transformResponse(const Vector &modResp, 
				    Vector &unmodResp)
{
    // perform T R  -- as T is block diagonal do T(i) R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRowOriginal = 0;
    int startRowTransformed = 0;
    int numNode = numGroups;
    int noRows = 0;
    int noCols = 0;

    for (int i=0; i<numNode; i++) {
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	    noRows = Ti->noRows();
	    noCols = Ti->noCols();
	    for (int j=0; j<noRows; j++) {
		double sum = 0.0;
		for (int k=0; k<noCols; k++)
		    sum += (*Ti)(j,k) * modResp(startRowTransformed +k) ;
		unmodResp(startRowOriginal + j) = sum;
	    }
	} else {
	    noCols = theDOFs[i]->getNumDOF();
	    noRows = noCols;
	    for (int j=0; j<noCols; j++)
		unmodResp(startRowOriginal + j) = modResp(startRowTransformed +j);
	}
	startRowOriginal += noRows;
	startRowTransformed += noCols;
    }

    return 0;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void  
TransformationFE::addD_ForceSensitivity(int gradNumber, const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalD_ForceSensitivity(gradNumber, response, fact);
}   	 

void  
TransformationFE::addM_ForceSensitivity(int gradNumber, const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalM_ForceSensitivity(gradNumber, response, fact);
}   	 

// AddingSensitivity:END ////////////////////////////////////
