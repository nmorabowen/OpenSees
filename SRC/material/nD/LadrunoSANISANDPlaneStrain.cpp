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

// Description: This file contains the implementation of the
// LadrunoSANISANDPlaneStrain class. Ported member-for-member from
// ManzariDafaliasPlaneStrain.cpp (UW, Pedro Arduino, 11.2011); every sign
// flip and index mapping below is carried over byte-for-byte -- the
// material is compression-positive internally and tension-positive at the
// element face, and the 3-component plane-strain packing indexes the
// 6-component tensor at {0,1,3} (skipping the out-of-plane normal
// component, index 2).

#include "LadrunoSANISANDPlaneStrain.h"

Vector LadrunoSANISANDPlaneStrain::mEpsilon_M(3);
Vector LadrunoSANISANDPlaneStrain::mSigma_M(3);
Vector LadrunoSANISANDPlaneStrain::rSigma(4);
Matrix LadrunoSANISANDPlaneStrain::mTangent(3,3);
Matrix LadrunoSANISANDPlaneStrain::mTangent_init(3,3);

// full constructor
LadrunoSANISANDPlaneStrain::LadrunoSANISANDPlaneStrain(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme,
	int tangentType, int JacoType, double TolF, double TolR,
	double Presidual, double Pmin, int honorTolR, int maxSubsteps) // Ladruno
:LadrunoSANISAND(tag, ND_TAG_LadrunoSANISANDPlaneStrain, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen,
				integrationScheme, tangentType, JacoType, TolF, TolR, Presidual, Pmin, honorTolR,
				maxSubsteps) // Ladruno
{
}

// null constructor
LadrunoSANISANDPlaneStrain::LadrunoSANISANDPlaneStrain()
  : LadrunoSANISAND(ND_TAG_LadrunoSANISANDPlaneStrain) // Ladruno: vanilla ManzariDafaliasPlaneStrain() calls the base's
                                                        // plain null ctor here, which hardcodes ND_TAG_ManzariDafalias --
                                                        // that classTag bug is deliberately NOT replicated.
// WHY IT MATTERS HERE AND NOT IN VANILLA (do not "align" this with the original):
// nothing ever repairs the tag -- there is no setClassTag call in ManzariDafalias.cpp
// and recvSelf does not restore it -- so a broker- or database-constructed object keeps
// it for life. In vanilla that is nearly invisible, because the usual getCopy() path
// (`*clone = *this`) copies a correct tag over the wrong one. For THIS class it would
// not be: echoLadrunoConstants() and Print() both branch on
// `getClassTag() == ND_TAG_LadrunoSANISAND`, so a restored PlaneStrain carrying the base
// tag would echo ONCE PER GAUSS POINT -- the ~83 MB-of-stderr failure ADR 86 sec.4.4's
// refinement box exists to prevent. See LEDGER_quirks.md; the vanilla fix is one token
// and is deliberately NOT made here (WORKFLOW_GOTCHAS sec.6: ask first).
{
}

// destructor
LadrunoSANISANDPlaneStrain::~LadrunoSANISANDPlaneStrain()
{
}

// make a clone of this material
NDMaterial*
LadrunoSANISANDPlaneStrain::getCopy()
{
    LadrunoSANISANDPlaneStrain  *clone;
    clone = new LadrunoSANISANDPlaneStrain();
    *clone = *this;
    return clone;
}

// send back type of material
const char*
LadrunoSANISANDPlaneStrain::getType() const
{
    return "PlaneStrain";
}

// send back order of strain
int
LadrunoSANISANDPlaneStrain::getOrder() const
{
    return 3;
}

// get the strain and integrate plasticity equations
int
LadrunoSANISANDPlaneStrain::setTrialStrain(const Vector &strain_from_element)
{
	mEpsilon.Zero();
	mEpsilon(0) = -1.0 * strain_from_element(0); // -1.0 is for geotechnical sign convention
	mEpsilon(1) = -1.0 * strain_from_element(1);
	mEpsilon(3) = -1.0 * strain_from_element(2);

    this->integrate();

	// Ladruno (ADR-86b): see the twin note in LadrunoSANISAND3D::setTrialStrain.
	// 0 unless this deck asked for a substep cap AND that cap fired.
	return this->ladrunoUpdateStatus();
}

// unused trial strain functions
int
LadrunoSANISANDPlaneStrain::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector&
LadrunoSANISANDPlaneStrain::getStrain()
{
	mEpsilon_M(0) = -1.0 * mEpsilon(0); // -1.0 is for geotechnical sign convention
	mEpsilon_M(1) = -1.0 * mEpsilon(1);
	mEpsilon_M(2) = -1.0 * mEpsilon(3);

    return mEpsilon_M;
}

// send back the stress
const Vector&
LadrunoSANISANDPlaneStrain::getStress()
{
	mSigma_M(0) = -1.0 * mSigma(0); // -1.0 is for geotechnical sign convention
	mSigma_M(1) = -1.0 * mSigma(1);
	mSigma_M(2) = -1.0 * mSigma(3);

 	return mSigma_M;
}


const Vector&
LadrunoSANISANDPlaneStrain::getStressToRecord()
{
	rSigma(0) = -1.0 * mSigma(0); // -1.0 is for geotechnical sign convention
	rSigma(1) = -1.0 * mSigma(1);
	rSigma(2) = -1.0 * mSigma(2);
	rSigma(3) = -1.0 * mSigma(3);

	return rSigma;
}

// send back the tangent
const Matrix&
LadrunoSANISANDPlaneStrain::getTangent()
{
	Matrix C(6,6);
	if (mTangType == 0)
		C = mCe;
	else if (mTangType == 1)
		C = mCep;
	else
		C = mCep_Consistent;

	mTangent(0,0) = C(0,0);
	mTangent(0,1) = C(0,1);
	mTangent(0,2) = C(0,3);
	mTangent(1,0) = C(1,0);
	mTangent(1,1) = C(1,1);
	mTangent(1,2) = C(1,3);
	mTangent(2,0) = C(3,0);
	mTangent(2,1) = C(3,1);
	mTangent(2,2) = C(3,3);

    return mTangent;
}

// send back the tangent
const Matrix&
LadrunoSANISANDPlaneStrain::getInitialTangent()
{
	mTangent_init(0,0) = mCe(0,0);
	mTangent_init(0,1) = mCe(0,1);
	mTangent_init(0,2) = mCe(0,3);
	mTangent_init(1,0) = mCe(1,0);
	mTangent_init(1,1) = mCe(1,1);
	mTangent_init(1,2) = mCe(1,3);
	mTangent_init(2,0) = mCe(3,0);
	mTangent_init(2,1) = mCe(3,1);
	mTangent_init(2,2) = mCe(3,3);

    return mTangent_init;
}

//send back the state parameters
const Vector
LadrunoSANISANDPlaneStrain::getState()
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
