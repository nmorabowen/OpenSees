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
                                                                        
// $Revision: 1.15 $
// $Date: 2010-04-23 22:56:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/Brick.h,v $

// Ed "C++" Love
//
// Eight node Brick element 
//

#ifndef BRICK_H
#define BRICK_H


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Damping.h>


class Brick : public Element {

  public :
    
    //null constructor
    Brick();
  
    //full constructor
    Brick(int tag, 
	  int node1,
	  int node2,
	  int node3,
	  int node4,
	  int node5,
	  int node6,
	  int node7,
	  int node8,
	  NDMaterial &theMaterial,
	  double b1 = 0.0, double b2 = 0.0, double b3 = 0.0,
	  Damping *theDamping = 0,
    int massType = 0);
    
    //destructor 
    virtual ~Brick( ) ;

    const char *getClassType(void) const {return "Brick";};

    //set domain
    void setDomain( Domain *theDomain ) ;
    int setDamping(Domain *theDomain, Damping *theDamping);

    // Ladruno WP-192: threaded element loops opt-in. The class / function /
    // shp3d scratch is thread_local; the answer is the materials' (and no
    // Damping object, which is not audited).
    bool ladrunoThreadSafeUpdate(void) const;
    bool ladrunoThreadSafeForm(void) const;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;

    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs(void);

    //return number of dofs
    int getNumDOF( ) ;

    //commit state
    int commitState( ) ;
    
    //revert to last commit 
    int revertToLastCommit( ) ;
    
    //revert to start 
    int revertToStart( ) ;

    // update
    int update(void);

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;
	
    //return stiffness matrix 
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();    
    const Matrix &getMass();    

    void zeroLoad( ) ;
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
      
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);


    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    //plotting 
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

  private : 
    //
    // private attributes
    //

    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[8] ;      //pointers to eight nodes

    //material information
    NDMaterial *materialPointers[8]; //pointers to eight materials

    double b[3];		// Body forces
    double appliedB[3];		// Body forces applied with load
    int applyLoad;

    Vector *load;
    Matrix *Ki;
    int massType;

    //
    // static attributes
    //

    static thread_local Matrix stiff ;   // Ladruno WP-192
    static thread_local Vector resid ;   // Ladruno WP-192
    static thread_local Matrix mass ;   // Ladruno WP-192
    static thread_local Matrix damping ;   // Ladruno WP-192

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static const double sg[2] ;
    static const double wg[8] ;
  
    //local nodal coordinates, three coordinates for each of four nodes
    static thread_local double xl[3][8] ;    // Ladruno WP-192

    //
    // private methods
    //

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute coordinate system
    void computeBasis( ) ;

    //compute B matrix
    const Matrix& computeB( int node, const double shp[4][8] ) ;
  
    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    Damping *theDamping[8];

    // Ladruno WP-195: element tangent cache. The 24x24 stiffness of a
    // small-strain brick is a pure function of the nodal coordinates (xl) and
    // the eight Gauss-point material tangents. When, at a tangent call, the
    // RAW BITS of all eight material tangents and of xl equal those that built
    // the cached stiffness, the cached stiffness is returned instead of the
    // 8 x 8 x 8 B^T D B products -- bit-identical by construction (same
    // inputs, same arithmetic). Elasticity is never assumed: it is observed,
    // on every call, through the material's own getTangent(). Disabled with a
    // Damping object (its stiffness multiplier is not a material tangent) and
    // with LADRUNO_BRICK_TANGENT_CACHE=0 in the environment at construction.
    // Explicitly invalidated by revertToStart, setDomain, setDamping,
    // setParameter/updateParameter and recvSelf. After 16 consecutive
    // misses the comparison pauses (plain path, buffer freed); every 64
    // commits of a paused element it is re-enabled (re-probe). See
    // Ladruno_implementation/195_ladruno_brick_tangent_cache_adr.md.
    double *ladrunoKcache;    // [576 stiff | 24 xl | 36 or 288 tangent bits]
    int ladrunoKcacheState;   // 0 empty, 1 valid (uniform D), 2 valid (per-GP D)
    int ladrunoKcacheCap;     // doubles allocated in ladrunoKcache
    int ladrunoKcacheMisses;  // consecutive misses; the cache pauses at 16
    int ladrunoKcacheIdle;    // commits since the pause; re-probe at 64
    double ladrunoKcacheReprobes; // pauses ended by the periodic re-probe
    bool ladrunoKcacheOn;     // enabled for this element
    void ladrunoInvalidateTangentCache(void);
    double ladrunoKcacheHits;   // tangent calls served from the cache
    double ladrunoKcacheBuilds; // tangent calls that (re)built the cache
    void formResidAndTangentCached(void);

} ; 

#endif

