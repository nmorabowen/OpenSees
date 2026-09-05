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

// Description: This file contains the implementation of the LadrunoSANISAND3D
// class. Ported member-for-member from ManzariDafalias3D.cpp (UW, Pedro
// Arduino, 11.2011); every sign flip below is carried over byte-for-byte --
// the material is compression-positive internally and tension-positive at
// the element face.

#include "LadrunoSANISAND3D.h"

Vector LadrunoSANISAND3D::mEpsilon_M(6);
Vector LadrunoSANISAND3D::mSigma_M(6);

// full constructor
LadrunoSANISAND3D::LadrunoSANISAND3D(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme,
	int tangentType, int JacoType, double TolF, double TolR,
	double Presidual, double Pmin, int honorTolR, int maxSubsteps) // Ladruno
:LadrunoSANISAND(tag, ND_TAG_LadrunoSANISAND3D, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen,
				integrationScheme, tangentType, JacoType, TolF, TolR, Presidual, Pmin, honorTolR,
				maxSubsteps) // Ladruno
{
}


// null constructor
LadrunoSANISAND3D::LadrunoSANISAND3D()
  : LadrunoSANISAND(ND_TAG_LadrunoSANISAND3D) // Ladruno
{

}

// destructor
LadrunoSANISAND3D::~LadrunoSANISAND3D()
{
}

// make a clone of this material
NDMaterial*
LadrunoSANISAND3D::getCopy()
{
    LadrunoSANISAND3D  *clone;
    clone = new LadrunoSANISAND3D();
    *clone = *this;
    return clone;
}

// send back type of material
const char*
LadrunoSANISAND3D::getType() const
{
    return "ThreeDimensional";
}

// send back order of strain
int
LadrunoSANISAND3D::getOrder() const
{
    return 6;
}

// get the strain and integrate plasticity equations
int
LadrunoSANISAND3D::setTrialStrain(const Vector &strain_from_element)
{
	mEpsilon = -1.0 * strain_from_element; // -1.0 is for geotechnical sign convention
	this->integrate();

	// Ladruno (ADR-86b): the ONE line that differs from ManzariDafalias3D here.
	// Vanilla returns a hardcoded 0, so a state determination that did not
	// actually integrate is indistinguishable from one that did -- which is how
	// a 34-minute ModifiedEuler seizure reached the analysis as "converged"
	// (ADR-90 GATE U). ladrunoUpdateStatus() is 0 unless this deck asked for a
	// substep cap AND that cap fired, so an uncapped deck is byte-identical.
	return this->ladrunoUpdateStatus();
}

// unused trial strain functions
int
LadrunoSANISAND3D::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector&
LadrunoSANISAND3D::getStrain()
{
	mEpsilon_M = -1.0 * mEpsilon;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
}

// send back the strain
const Vector&
LadrunoSANISAND3D::getEStrain()
{
	mEpsilon_M = -1.0 * mEpsilonE;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
}

const Vector&
LadrunoSANISAND3D::getPStrain()
{
	mEpsilon_M = -1.0 * (mEpsilon - mEpsilonE);
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
}


// send back the stress
const Vector&
LadrunoSANISAND3D::getStress()
{
	// this->integrate();
	mSigma_M = -1.0 * mSigma;
 	return mSigma_M; // -1.0 is for geotechnical sign convention
}

const Vector&
LadrunoSANISAND3D::getStressToRecord()
{
	return getStress();
}

// send back the tangent
const Matrix&
LadrunoSANISAND3D::getTangent()
{
    if (mTangType == 0)
		return mCe;
	else if (mTangType == 1)
		return mCep;
	else
		return mCep_Consistent;
}

// send back the tangent
const Matrix&
LadrunoSANISAND3D::getInitialTangent()
{
    return mCe;
}

//send back the state parameters
const Vector
LadrunoSANISAND3D::getState()
 {
	 Vector result(26);
	 result.Assemble(mEpsilonE,0,1.0);
	 result.Assemble(mAlpha,6,1.0);
	 result.Assemble(mFabric,12,1.0);
	 result.Assemble(mAlpha_in,18,1.0);
	 result(24) = mVoidRatio;
	 result(25) = mDGamma;

	 return result;
 }
