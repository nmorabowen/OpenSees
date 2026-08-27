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

// Description: This file contains the class definition for
// LadrunoSANISANDPlaneStrain, the plane-strain wrapper for LadrunoSANISAND.
// Modelled member-for-member on
// SRC/material/nD/UWmaterials/ManzariDafaliasPlaneStrain.{h,cpp}, but
// derived from LadrunoSANISAND (not ManzariDafaliasPlaneStrain) so that
// LadrunoSANISAND's revertToStart/sendSelf/recvSelf/getCopy/Print overrides
// are not lost.

#ifndef LadrunoSANISANDPlaneStrain_h
#define LadrunoSANISANDPlaneStrain_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include <classTags.h> // Ladruno
#include "LadrunoSANISAND.h" // Ladruno

class LadrunoSANISANDPlaneStrain : public LadrunoSANISAND // Ladruno
{
  public :

  //full constructor
  LadrunoSANISANDPlaneStrain(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme = 2,
	int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7,
	double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0); // Ladruno: trailing low-stress-constant args
  //null constructor
  LadrunoSANISANDPlaneStrain();

  //destructor
  ~LadrunoSANISANDPlaneStrain();

  NDMaterial *getCopy(void);
  const char *getType(void) const;
  int getOrder(void) const;

  int setTrialStrain(const Vector &strain_from_element);

  // Unused trialStrain functions
  int setTrialStrain(const Vector &v, const Vector &r);

  //send back the strain
  const Vector& getStrain();

  //send back the stress
  const Vector& getStress();
  const Vector& getStressToRecord();

  //send back the tangent
  const Matrix& getTangent();
  const Matrix& getInitialTangent();

  //send back the state parameters
  const Vector  getState();

  private :

  // static vectors and matrices
  static Vector mSigma_M  ; // mSigma with continuum mechanic sign convention
  static Vector mEpsilon_M; // mEpsilon with continuum mechanic sign convention
  static Vector rSigma;     // Stress for the recorders
  static Matrix mTangent;
  static Matrix mTangent_init;

};

#endif
