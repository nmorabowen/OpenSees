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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStrain2D.h,v $
                                                                        
                                                                        
#ifndef ElasticIsotropicPlaneStrain2D_h
#define ElasticIsotropicPlaneStrain2D_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: 
//
// What: "@(#) ElasticIsotropicPlaneStrain2D.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticIsotropicPlaneStrain2D : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropicPlaneStrain2D (int tag, double E, double nu, double rho);
	ElasticIsotropicPlaneStrain2D ();
    ~ElasticIsotropicPlaneStrain2D ();

    const char *getClassType(void) const {return "ElasticIsotropicPlaneStrain2D";};

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);

    const Vector &getStress (void);
    const Vector &getStrain (void);
    double getStressZZ (void);   // Ladruno: sigma_zz = lambda*(eps_xx + eps_yy)

    // Ladruno WP-107 (ADR-75b L3-1): re-entrant ON THE update() PATH, which is
    // the only claim being made.
    //
    // CORRECTED after the red-team review (S2). The first version of this
    // comment said "getTangent/getInitialTangent/getStress write only the
    // instance's own D/sigma". That is FALSE and the contradiction was ten
    // lines below it: `sigma` and `D` are CLASS-WIDE statics (see the private
    // section, and their definitions at the top of the .cpp) and all three of
    // those methods write them. Two instances calling getTangent() on two
    // threads race, full stop.
    //
    // The conclusion survives, for a different and narrower reason: the loop
    // WP-107 threads is Domain::update(), and LadrunoQuad::update() reaches
    // exactly one method on this class -- setTrialStrain(), which writes only
    // the instance's own `epsilon` (.cpp:51-77, four overloads, all
    // unconditionally `return 0`). getTangent / getInitialTangent / getStress
    // are called from getTangentStiff / getResistingForce, which run on the
    // main thread in loops B/C and are NOT threaded by this WP.
    //
    // WHY THE DISTINCTION MATTERS RATHER THAN BEING PEDANTRY: the wrong reason
    // licenses the next mistake. It would authorise allowlisting any element
    // whose update() calls getStress()/getTangent() -- which `LadrunoQuad -eas`
    // (formEAStrue) does, and which is precisely why -eas is NOT allowlisted.
    // The banked lesson of this WP is "a `static` grep is an audit, not a
    // re-entrancy proof"; here the grep had not even been run.
    //
    // So: safe because the update() call graph reaches ONLY setTrialStrain.
    // `D` and `sigma` are class-wide and must never be touched on a worker
    // thread. If loops B/C are ever threaded, this declaration is void.
    //
    // Declared on this LEAF class rather than on ElasticIsotropicMaterial on
    // purpose -- the base has thermal and incremental subclasses that have not
    // been audited, and a base-class `true` would fail OPEN for them.
    virtual bool ladrunoThreadSafeUpdate(void) const { return true; }   // Ladruno WP-107

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

  protected:

  private:
    static Vector sigma;        // Stress vector ... class-wide for returns
    static Matrix D;	        // Elastic constants
    Vector epsilon;	        // Trial strains
    Vector Cepsilon;	        // Committed strains
};


#endif


