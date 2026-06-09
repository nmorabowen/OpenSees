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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: Massimo Petracca - ASDEA Software, Italy
// Created: 03/2024

#include <stdlib.h>
#include <math.h>

#include <InitStrainNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <Parameter.h>
#include <OPS_Globals.h>

#include <elementAPI.h>

namespace {

    // a base index to try not to shadow the adpted material's parameters...
    constexpr int param_base_index = 111000;

    // Ladruno: 3D-Voigt {11,22,33,12,23,13} -> dimensional-view component map,
    // keyed by the view order (same convention as LadrunoJ2's vmap). Only the two
    // views the base NDMaterial::getCopy() cannot build (PlaneStrain, AxiSymmetric)
    // are produced natively; everything else (PlaneStress/PlateFiber/BeamFiber/3D)
    // keeps its existing 3D-wrapped path, where the inner stays order 6.
    void initStrainVMap(int ncomp, int *vmap)
    {
        if (ncomp == 3) { vmap[0] = 0; vmap[1] = 1; vmap[2] = 3; }              // PlaneStrain
        else if (ncomp == 4) { vmap[0] = 0; vmap[1] = 1; vmap[2] = 2; vmap[3] = 3; } // AxiSymmetric
        else { for (int a = 0; a < 6; ++a) vmap[a] = a; }                       // ThreeDimensional
    }

}

void
InitStrainNDMaterial::projectEps0(void)
{
    // Ladruno: epsInit (size ncomp) <- epsInit3D (size 6) through the view map.
    int n = epsInit.Size();
    int vmap[6];
    initStrainVMap(n, vmap);
    for (int a = 0; a < n; ++a)
        epsInit(a) = epsInit3D(vmap[a]);
}

void*
OPS_InitStrainNDMaterial(void)
{
    // Pointer to an ND material that will be returned
    NDMaterial* theMaterial = 0;
    NDMaterial* theOtherMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs != 3 && numArgs != 8) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Only 3 or 8 arguments are required, not " << numArgs << endln
            << "nDMaterial InitStrain tag? otherTag? eps0_11? <eps0_22 eps0_33 eps0_12 eps0_23 eps0_13>" << endln;
        return theMaterial;
    }

    int iData[2];
    int numData = 2;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Cannot read $tag and $otherTag" << numArgs << endln;
        return theMaterial;
    }

    theOtherMaterial = OPS_getNDMaterial(iData[1]);
    if (theOtherMaterial == 0) {
        opserr << "Error while creating InitStrain material:" << endln 
            << "Could not find material with tag: " << iData[1] << endln;
        return theMaterial;
    }

    double eps0[6] = { 0,0,0, 0,0,0 };
    numData = numArgs - 2;
    if (OPS_GetDoubleInput(&numData, eps0) != 0) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Could not read initial strain" << endln;
        return theMaterial;
    }

    // Parsing was successful, allocate the material
    if (numArgs == 3) {
        theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, eps0[0]);
    }
    else {
        Vector E0(6);
        for (int i = 0; i < 6; ++i) E0(i) = eps0[i];
        theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, E0);
    }

    if (theMaterial == 0) {
        opserr << "WARNING could not create NDMaterial of type InitStrain\n";
        return theMaterial;
    }

    return theMaterial;
}

InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, const Vector& eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(eps0)
    , epsInit3D(eps0)   // Ladruno: this ctor builds a 3D template -> epsInit == epsInit3D
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // make sure the input strain vector is of correct size
    if (epsInit.Size() != 6) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- input eps0 vector of incorrect size\n";
        exit(-1);
    }
}

InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, double eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(0)
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // initialize epsInit (3D template: isotropic volumetric prestrain on the normals)
    epsInit.resize(6);
    epsInit.Zero();
    for (int i = 0; i < 3; ++i)
        epsInit(i) = eps0;
    epsInit3D = epsInit;   // Ladruno
}

// Ladruno: ADOPTING ctor (dimension-general). Preserves the already-reduced inner
// view; epsInit is view-sized, epsInit3D is the canonical 3D source for eps0_ij.
InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial* innerAdopt,
    const Vector& eps0View, const Vector& eps0Full)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(innerAdopt)
    , epsInit(eps0View)
    , epsInit3D(eps0Full)
{
    // innerAdopt may be null on a failed view copy; callers (getCopy) check.
}

InitStrainNDMaterial::InitStrainNDMaterial()
    : NDMaterial(0, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(6)
    , epsInit3D(6)
{

}

InitStrainNDMaterial::~InitStrainNDMaterial()
{
    if (theMaterial)
        delete theMaterial;
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain)
{
    // Ladruno: size to the inner order (was hardcoded 6) so reduced-order
    // (PlaneStrain/AxiSymmetric) views feed correctly-sized strains. operator=
    // reallocs on size mismatch, so total_strain takes epsInit's order.
    static Vector total_strain;
    total_strain = epsInit;
    total_strain.addVector(1.0, strain, 1.0);
    return theMaterial->setTrialStrain(total_strain);
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrain(strain);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain)
{
    // Ladruno: inner-order sized (was hardcoded 6). getStrain() is the inner total
    // (= eps0 + ele strain); subtract eps0 to recover the element strain, then add
    // the increment and re-impose eps0 via setTrialStrain.
    static Vector strain_from_ele;
    strain_from_ele = theMaterial->getStrain();
    strain_from_ele.addVector(1.0, epsInit, -1.0);
    strain_from_ele.addVector(1.0, strain, 1.0);
    return setTrialStrain(strain_from_ele);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrainIncr(strain);
}

const Vector&
InitStrainNDMaterial::getStress(void)
{
    return theMaterial->getStress();
}

const Matrix&
InitStrainNDMaterial::getTangent(void)
{
    return theMaterial->getTangent();
}

const Matrix&
InitStrainNDMaterial::getInitialTangent(void)
{
    return theMaterial->getInitialTangent();
}

const Vector&
InitStrainNDMaterial::getStrain(void)
{
    return theMaterial->getStrain();
}

int
InitStrainNDMaterial::commitState(void)
{
    return theMaterial->commitState();
}

int
InitStrainNDMaterial::revertToLastCommit(void)
{
    return theMaterial->revertToLastCommit();
}

int
InitStrainNDMaterial::revertToStart(void)
{
    return theMaterial->revertToStart();
}

double
InitStrainNDMaterial::getRho(void)
{
    return theMaterial->getRho();
}

NDMaterial*
InitStrainNDMaterial::getCopy(void)
{
    // Ladruno: clone preserving the CURRENT dimensional view (was: always rebuilt a
    // 3D copy, which corrupted a reduced-order view copy). theMaterial->getCopy() is
    // a same-type copy; epsInit/epsInit3D carry the view-sized and canonical eps0.
    NDMaterial* innerCopy = (theMaterial != nullptr) ? theMaterial->getCopy() : nullptr;
    return new InitStrainNDMaterial(getTag(), innerCopy, epsInit, epsInit3D);
}

NDMaterial *
InitStrainNDMaterial::getCopy(const char *type)
{
    if (strcmp(type, "ThreeDimensional") == 0)
        return getCopy();

    // Ladruno: dimension-general. The base NDMaterial::getCopy() builds PlaneStress/
    // PlateFiber/BeamFiber by wrapping a 3D InitStrain in a condensing view, but it
    // returns 0 for PlaneStrain and AxiSymmetric -- so those (PlaneStrain is the
    // default for LadrunoQuad/CST) used to fail with a null material. Build them
    // natively: inner = the requested view of the stored 3D material; eps0 = the
    // canonical 3D eps0 reduced to that view's component order. Only valid when this
    // instance is the 3D template (epsInit3D sized 6); otherwise fall through.
    if (epsInit3D.Size() == 6) {
        int ncomp = 0;
        if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)
            ncomp = 3;
        else if (strcmp(type, "AxiSymmetric") == 0 || strcmp(type, "AxiSymmetric2D") == 0)
            ncomp = 4;

        if (ncomp != 0) {
            NDMaterial* innerView = theMaterial->getCopy(type);
            if (innerView == nullptr) {
                opserr << "InitStrainNDMaterial::getCopy - the wrapped material cannot "
                       << "provide a " << type << " view\n";
                return nullptr;
            }
            int vmap[6];
            initStrainVMap(ncomp, vmap);
            Vector eps0View(ncomp);
            for (int a = 0; a < ncomp; ++a)
                eps0View(a) = epsInit3D(vmap[a]);
            return new InitStrainNDMaterial(getTag(), innerView, eps0View, epsInit3D);
        }
    }

    // PlaneStress / PlateFiber / BeamFiber / unknown: unchanged legacy path.
    return NDMaterial::getCopy(type);
}

const char*
InitStrainNDMaterial::getType(void) const
{
    return theMaterial->getType();
}

int InitStrainNDMaterial::getOrder(void) const
{
    // Ladruno: delegate to the inner (was hardcoded 6). A native PlaneStrain/
    // AxiSymmetric view reports 3/4 so the element feeds correctly-sized strains;
    // the 3D template (and the InitStrain inside a base PlaneStress/fiber wrapper)
    // still reports 6.
    return (theMaterial != nullptr) ? theMaterial->getOrder() : 6;
}

int
InitStrainNDMaterial::sendSelf(int cTag, Channel& theChannel)
{
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::sendSelf() - theMaterial is null, nothing to send" << endln;
        return -1;
    }

    int dbTag = this->getDbTag();

    // Ladruno: dataID(3) carries the ACTIVE (view) order so recvSelf can size epsInit;
    // the canonical size-6 epsInit3D is what travels on the wire (epsInit is re-derived).
    static ID dataID(4);
    dataID(0) = this->getTag();
    dataID(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        theMaterial->setDbTag(matDbTag);
    }
    dataID(2) = matDbTag;
    dataID(3) = epsInit.Size();
    if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the ID\n";
        return -1;
    }

    if (theChannel.sendVector(dbTag, cTag, epsInit3D) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send epsInit3D\n";
        return -2;
    }

    if (theMaterial->sendSelf(cTag, theChannel) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the Material\n";
        return -3;
    }

    return 0;
}

int
InitStrainNDMaterial::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();

    static ID dataID(4);
    if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the ID\n";
        return -1;
    }
    setTag(dataID(0));

    // as no way to change material, don't have to check classTag of the material
    if (theMaterial == 0) {
        int matClassTag = dataID(1);
        theMaterial = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial == 0) {
            opserr << "InitStrainNDMaterial::recvSelf() - failed to create Material with classTag "
                << matClassTag << endln;
            return -2;
        }
    }
    theMaterial->setDbTag(dataID(2));

    // Ladruno: recover the canonical 3D eps0, then re-project to the active view
    // order (dataID(3)) carried alongside.
    epsInit3D.resize(6);
    if (theChannel.recvVector(dbTag, cTag, epsInit3D) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the epsInit3D vector\n";
        return -3;
    }
    epsInit.resize(dataID(3));
    this->projectEps0();

    if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the Material\n";
        return -4;
    }

    return 0;
}

void
InitStrainNDMaterial::Print(OPS_Stream& s, int flag)
{
    s << "InitStrainNDMaterial tag: " << this->getTag() << endln;
    s << "\tMaterial: " << theMaterial->getTag() << endln;
    s << "\tinitital strain: " << epsInit << endln;
}

int
InitStrainNDMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    // Ladruno: read/write the canonical 3D eps0 (epsInit3D) so the eps0_ij API has
    // stable component indices in every dimensional view; the active (view-sized)
    // epsInit is re-derived in updateParameter via projectEps0().
    if (argc > 0) {
        if (strcmp(argv[0], "eps0") == 0) {
            // eps0 is assumed to impose with a single scalar a volumetric stress = I*eps0
            param.setValue(epsInit3D(0));
            return param.addObject(param_base_index, this);
        }
        else if (strcmp(argv[0], "eps0_11") == 0) {
            param.setValue(epsInit3D(0));
            return param.addObject(param_base_index + 1, this);
        }
        else if (strcmp(argv[0], "eps0_22") == 0) {
            param.setValue(epsInit3D(1));
            return param.addObject(param_base_index + 2, this);
        }
        else if (strcmp(argv[0], "eps0_33") == 0) {
            param.setValue(epsInit3D(2));
            return param.addObject(param_base_index + 3, this);
        }
        else if (strcmp(argv[0], "eps0_12") == 0) {
            param.setValue(epsInit3D(3));
            return param.addObject(param_base_index + 4, this);
        }
        else if (strcmp(argv[0], "eps0_23") == 0) {
            param.setValue(epsInit3D(4));
            return param.addObject(param_base_index + 5, this);
        }
        else if (strcmp(argv[0], "eps0_13") == 0) {
            param.setValue(epsInit3D(5));
            return param.addObject(param_base_index + 6, this);
        }
    }
    return theMaterial->setParameter(argv, argc, param);
}

int InitStrainNDMaterial::updateParameter(int parameterID, Information& info)
{
    // Ladruno: edit the canonical 3D eps0 then re-project onto the active view.
    // Components absent from the current view (e.g. eps0_23 in PlaneStrain) update
    // epsInit3D but project to nothing -- no effect, no out-of-bounds.
    switch (parameterID) {
    case param_base_index:
        epsInit3D(0) = epsInit3D(1) = epsInit3D(2) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 1:
        epsInit3D(0) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 2:
        epsInit3D(1) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 3:
        epsInit3D(2) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 4:
        epsInit3D(3) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 5:
        epsInit3D(4) = info.theDouble;
        this->projectEps0();
        return 0;
    case param_base_index + 6:
        epsInit3D(5) = info.theDouble;
        this->projectEps0();
        return 0;
    default:
        return -1;
    }
}

Response* InitStrainNDMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    return theMaterial->setResponse(argv, argc, output);
}

const Vector&
InitStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
InitStrainNDMaterial::commitSensitivity(const Vector& depsdh,
    int gradIndex, int numGrads)
{
    return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
