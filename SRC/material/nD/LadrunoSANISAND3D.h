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

// Description: This file contains the class definition for LadrunoSANISAND3D,
// the 3D wrapper for LadrunoSANISAND. Modelled member-for-member on
// SRC/material/nD/UWmaterials/ManzariDafalias3D.{h,cpp}, but derived from
// LadrunoSANISAND (not ManzariDafalias3D) so that LadrunoSANISAND's
// revertToStart/sendSelf/recvSelf/getCopy/Print overrides are not lost.
#ifndef LadrunoSANISAND3D_h
#define LadrunoSANISAND3D_h


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <classTags.h> // Ladruno
#include "LadrunoSANISAND.h" // Ladruno

class LadrunoSANISAND3D : public LadrunoSANISAND // Ladruno
{
  public :

  //full constructor
  LadrunoSANISAND3D(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme = 2,
	int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7,
	double Presidual = 0.0, double Pmin = -1.0, int honorTolR = 0); // Ladruno: trailing low-stress-constant args

  //null constructor
  LadrunoSANISAND3D();

  //destructor
  ~LadrunoSANISAND3D();

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
  const Vector getState();


  const Vector& getEStrain();
  const Vector& getPStrain();

  private :

  static Vector mSigma_M  ; // mSigma with continuum mechanic sign convention
  static Vector mEpsilon_M; // mEpsilon with continuum mechanic sign convention

};


#endif
