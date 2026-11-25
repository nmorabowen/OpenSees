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
                                                                        
// Written: Chris McGann
//          February 2011

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include "DruckerPrager.h"

class DruckerPragerPlaneStrain : public DruckerPrager {

  public : 

  //null constructor
  DruckerPragerPlaneStrain();

  //full constructor
  DruckerPragerPlaneStrain(int tag, double bulk, double shear,
		  double s_y, double r, double r_bar, double Kinfinity, double Kinit, 
		  double d1, double d2, double H, double t, double massDens, double atm);


  //destructor
  ~DruckerPragerPlaneStrain();

  NDMaterial* getCopy();
  const char* getType() const;
  int getOrder() const;

  int setTrialStrain(const Vector &strain_from_element);

  // Unused trialStrain functions
  int setTrialStrain(const Vector &v, const Vector &r);
    
  //send back the strain
  const Vector& getStrain();

  //send back the stress
  const Vector& getStress();

  // Ladruno: out-of-plane sigma_zz from the internal 6-component stress
  double getStressZZ();

  //send back the tangent 
  const Matrix& getTangent();
  const Matrix& getInitialTangent();

  // Ladruno WP-258: re-entrant on DIFFERENT instances (WP-107 / WP-192 contract,
  // Material.h). plastic_integrator() works on per-instance members and stack
  // temporaries, its Matrix::Invert uses the per-thread work area inside a
  // threaded region, and its two diagnostics go through LadrunoOutputGuard.
  // The condensed return buffers below are thread_local (WP-258).
  // Audit: Ladruno_implementation/258_thread_safe_soil_materials.md section 2.1.
  bool ladrunoThreadSafeUpdate(void) const { return true; }
  bool ladrunoThreadSafeForm(void) const { return true; }

  private :

  // static vectors and matrices
  // Ladruno WP-258: thread_local -- getStress/getTangent/getInitialTangent fill
  // and return them, so two instances on two threads would share one buffer.
  static thread_local Vector strain;
  static thread_local Vector stress;
  static thread_local Matrix tangent;


}; 


