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
                                                                        
// $Revision: 1.19 $
// $Date: 2009-08-25 22:32:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for Element.
// Element is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) Element.h, revA"

#ifndef Element_h
#define Element_h

#include <DomainComponent.h>
#include <ID.h>

class Matrix;
class Vector;
class Info;
class Response;
class ElementalLoad;
class Node;
class Damping;

class Element : public DomainComponent
{
  public:
    Element(int tag, int classTag);    
    virtual ~Element();

    // methods dealing with nodes and number of external dof
    virtual int getNumExternalNodes(void) const =0;
    virtual const ID &getExternalNodes(void)  =0;	
    virtual Node **getNodePtrs(void)  =0;	
    virtual int getNumDOF(void) =0;
    virtual double getCharacteristicLength(void);

    // Ladruno (ADR 20 §9): single source of truth for an element's shape-function
    // weights at a natural coordinate. Given the element's natural coords `xi`
    // (hex: (ξ,η,ζ)∈[-1,1]³; tet: barycentric (L1,L2,L3)), fill `N` with the nodal
    // interpolation weights N_i so that  f(xi) = Σ_i N_i f_i  for any nodal field.
    // Used by LadrunoEmbeddedRebar to embed a rebar node in a non-matching solid
    // host without re-supplying the weights by hand. Default = not implemented
    // (returns -1, leaves N untouched); host elements override and return 0.
    virtual int getInterpolationWeights(const Vector &xi, Vector &N);

    // Ladruno (ADR 23 §3, Phase 2 UR): companion of getInterpolationWeights — the
    // CARTESIAN shape-function gradients at a natural coordinate. Fill `dNdx` (an
    // nNodes x ndm matrix) with dNdx(i,j) = dN_i/dx_j (global coords) so that the
    // displacement gradient at xi is  du_a/dx_j = Sum_i dNdx(i,j) u_i^a. Used by
    // LadrunoEmbeddedNode to tie a constrained node's ROTATIONS to the host's
    // continuum rotation theta = 1/2 curl(u) = skew(grad u) (weights alone cannot —
    // the rotation needs dN/dx, not N). Default = not implemented (returns -1, leaves
    // dNdx untouched); host elements override and return 0. For straight-sided
    // simplex hosts dN/dx is element-constant; for hex/Bezier hosts it varies with xi.
    virtual int getInterpolationGradients(const Vector &xi, Matrix &dNdx);

    // Ladruno (ADR 20 §10.6.1): an element's self-reported explicit critical time
    // step (seconds). A non-negative return FULLY REPLACES the per-element
    // K v = λ M v eigensolve for this element — CriticalTimeStep folds the value
    // into its running minimum and SKIPS the eigensolve (it does NOT take
    // min(self, eigensolve)). Default -1 = "no opinion" (the eigensolve governs).
    // CONTRACT: override (with a value > 0) ONLY when this element's per-element
    // pencil would contribute nothing meaningful on its own — i.e. its real
    // stability bound is invisible to the eigensolve. Returning > 0 while the
    // element ALSO carries a genuine stiff/massive mode would silently drop that
    // mode. The motivating case is LadrunoEmbeddedRebar's bipenalty bound: its host
    // DOFs are massless in the per-element problem, so they slave out the
    // constraint → λ_max=0 (the eigensolve sees no bound), and the self-report is
    // the only signal. An element that sometimes has a real mode should return -1
    // in that regime so the eigensolve governs.
    virtual double getExplicitCriticalTimeStep(void);

    // methods dealing with committed state and update
    virtual int commitState(void);    
    virtual int revertToLastCommit(void) = 0;        
    virtual int revertToStart(void);                
    virtual int update(void);
    virtual bool isSubdomain(void);
    
    // methods to return the current linearized stiffness,
    // damping and mass matrices
    virtual const Matrix &getTangentStiff(void) =0;
    virtual const Matrix &getInitialStiff(void) =0;
    virtual const Matrix &getDamp(void);    
    virtual const Matrix &getMass(void);
    virtual const Matrix &getGeometricTangentStiff();

    // methods for applying loads
    virtual void zeroLoad(void);	
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
    virtual int addLoad(ElementalLoad *theLoad, const Vector &loadFactors);

    virtual int addInertiaLoadToUnbalance(const Vector &accel);
    virtual int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);
    Vector getRayleighDampingFactors() const;


    virtual int setDamping(Domain *theDomain, Damping *theDamping);

    // methods for obtaining resisting force (force includes elemental loads)
    virtual const Vector &getResistingForce(void) =0;
    virtual const Vector &getResistingForceIncInertia(void);        

    // method for obtaining information specific to an element
    virtual Response *setResponse(const char **argv, int argc, 
				  OPS_Stream &theHandler);
    virtual int getResponse(int responseID, Information &eleInformation);
    virtual int getResponseSensitivity(int responseID, int gradIndex,
				       Information &eleInformation);

    virtual int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual int addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    virtual const Vector & getResistingForceSensitivity(int gradIndex);
    virtual const Matrix & getTangentStiffSensitivity(int gradIndex);
    virtual const Matrix & getInitialStiffSensitivity(int gradIndex);
    virtual const Matrix & getCommittedStiffSensitivity(int gradIndex);
    virtual const Matrix & getDampSensitivity(int gradIndex);
    virtual const Matrix & getMassSensitivity(int gradIndex);
    virtual int   commitSensitivity(int gradIndex, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

    virtual int addResistingForceToNodalReaction(int flag);

    virtual int storePreviousK(int numK);
    virtual const Matrix *getPreviousK(int num);

    virtual void onActivate();
    virtual void onDeactivate();

    void activate();
    void deactivate();

    bool isActive();



protected:
	const Vector& getRayleighDampingForces(void);

    double alphaM, betaK, betaK0, betaKc;
    Matrix *Kc; // pointer to hold last committed matrix if needed for rayleigh damping

    Matrix **previousK;
    int numPreviousK;

    int index, nodeIndex;

    static Matrix ** theMatrices; 
    static Vector ** theVectors1; 
    static Vector ** theVectors2; 
    static int numMatrices;

    bool is_this_element_active;

  private:
};


#endif

