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

// Ladruno: ETABS-style shell section stiffness-modifier decorator (ADR 91).
// See LadrunoShellModifierSection.h for the command syntax and the
// mathematical contract (the D' = S*D*S congruence transform).
//
// Written: N. Mora-Bowen (Ladruno), 2026.
// See Ladruno_implementation/91_ladruno_shell_stiffness_modifiers_adr.md

#include <LadrunoShellModifierSection.h>

#include <math.h>
#include <string.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <classTags.h>
#include <elementAPI.h>
#include <Parameter.h>
#include <ID.h>

// static work storage (order 8, fixed)
Vector LadrunoShellModifierSection::sigma(8);
Matrix LadrunoShellModifierSection::D(8, 8);
ID     LadrunoShellModifierSection::code(8);

void* OPS_LadrunoShellModifierSection()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: section LadrunoShellModifier tag? innerSecTag? "
               << "<-f11 v> <-f22 v> <-f12 v> <-m11 v> <-m22 v> <-m12 v> "
               << "<-v13 v> <-v23 v> <-mass v>" << endln;
        return 0;
    }

    int idata[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "WARNING LadrunoShellModifierSection: invalid tag or innerSecTag\n";
        return 0;
    }
    int tag = idata[0];
    int innerTag = idata[1];

    // R4: inner section must already exist
    SectionForceDeformation *theInner = OPS_getSectionForceDeformation(innerTag);
    if (theInner == 0) {
        opserr << "WARNING LadrunoShellModifierSection " << tag
               << ": inner section " << innerTag << " does not exist\n";
        return 0;
    }

    // R1: inner section must be order 8 (a plate/shell section)
    if (theInner->getOrder() != 8) {
        opserr << "WARNING LadrunoShellModifierSection " << tag
               << ": inner section " << innerTag << " has order "
               << theInner->getOrder()
               << " -- LadrunoShellModifier requires an order-8 plate/shell section\n";
        return 0;
    }

    // R2: inner section's response codes must be exactly the plate set,
    // in the standard order -- an order-8 section with different response
    // codes is not a plate section.
    static const int expectedCode[8] = {
        SECTION_RESPONSE_FXX, SECTION_RESPONSE_FYY, SECTION_RESPONSE_FXY,
        SECTION_RESPONSE_MXX, SECTION_RESPONSE_MYY, SECTION_RESPONSE_MXY,
        SECTION_RESPONSE_VXZ, SECTION_RESPONSE_VYZ
    };
    const ID &innerCode = theInner->getType();
    bool typeOk = (innerCode.Size() == 8);
    if (typeOk) {
        for (int i = 0; i < 8; i++) {
            if (innerCode(i) != expectedCode[i]) {
                typeOk = false;
                break;
            }
        }
    }
    if (!typeOk) {
        opserr << "WARNING LadrunoShellModifierSection " << tag
               << ": inner section " << innerTag << " is order-8 but its response "
               << "codes are not exactly {FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ} "
               << "-- it is not a plate/shell section\n";
        return 0;
    }

    double mods[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();

        int idx = -1;
        if (strcmp(opt, "-f11") == 0)        idx = 0;
        else if (strcmp(opt, "-f22") == 0)   idx = 1;
        else if (strcmp(opt, "-f12") == 0)   idx = 2;
        else if (strcmp(opt, "-m11") == 0)   idx = 3;
        else if (strcmp(opt, "-m22") == 0)   idx = 4;
        else if (strcmp(opt, "-m12") == 0)   idx = 5;
        else if (strcmp(opt, "-v13") == 0)   idx = 6;
        else if (strcmp(opt, "-v23") == 0)   idx = 7;
        else if (strcmp(opt, "-mass") == 0)  idx = 8;
        else {
            opserr << "WARNING LadrunoShellModifierSection " << tag
                   << ": unknown option '" << opt << "'\n";
            return 0;
        }

        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING LadrunoShellModifierSection " << tag
                   << ": option '" << opt << "' requires a value\n";
            return 0;
        }

        double val;
        int one = 1;
        if (OPS_GetDoubleInput(&one, &val) < 0) {
            opserr << "WARNING LadrunoShellModifierSection " << tag
                   << ": invalid value for option '" << opt << "'\n";
            return 0;
        }

        // R3: negative modifiers are refused; exactly 0.0 is allowed (warned later)
        if (val < 0.0) {
            opserr << "WARNING LadrunoShellModifierSection " << tag
                   << ": modifier '" << opt << "' = " << val
                   << " is negative -- modifiers must be >= 0.0\n";
            return 0;
        }

        mods[idx] = val;
    }

    return new LadrunoShellModifierSection(tag, theInner, mods);
}

// constructor for blank object that recvSelf needs to be invoked upon
LadrunoShellModifierSection::LadrunoShellModifierSection():
    SectionForceDeformation(0, SEC_TAG_LadrunoShellModifier),
    theSection(0), e(8), otherDbTag(0)
{
    for (int i = 0; i < 9; i++)
        mods[i] = 1.0;
    for (int i = 0; i < 8; i++)
        scale[i] = 1.0;
}

LadrunoShellModifierSection::LadrunoShellModifierSection(int tag,
                                                           SectionForceDeformation *theInnerSection,
                                                           const double m[9]):
    SectionForceDeformation(tag, SEC_TAG_LadrunoShellModifier),
    theSection(0), e(8), otherDbTag(0)
{
    if (theInnerSection == 0) {
        opserr << "LadrunoShellModifierSection::LadrunoShellModifierSection -- null inner section passed\n";
        exit(-1);
    }

    theSection = theInnerSection->getCopy();
    if (theSection == 0) {
        opserr << "LadrunoShellModifierSection::LadrunoShellModifierSection -- failed to copy inner section\n";
        exit(-1);
    }

    for (int i = 0; i < 9; i++)
        mods[i] = m[i];

    this->buildScales();
}

LadrunoShellModifierSection::~LadrunoShellModifierSection()
{
    if (theSection != 0)
        delete theSection;
}

void
LadrunoShellModifierSection::buildScales(void)
{
    static const char* names[9] = {
        "f11", "f22", "f12", "m11", "m22", "m12", "v13", "v23", "mass"
    };

    for (int i = 0; i < 9; i++) {
        if (mods[i] == 0.0) {
            opserr << "WARNING LadrunoShellModifierSection " << this->getTag()
                   << ": modifier '" << names[i] << "' is 0.0 -- the section "
                   << "is SINGULAR in that response mode\n";
        }
    }

    for (int i = 0; i < 8; i++)
        scale[i] = sqrt(mods[i]);
}

int
LadrunoShellModifierSection::setTrialSectionDeformation(const Vector &def)
{
    e = def;

    static Vector eIn(8);
    for (int i = 0; i < 8; i++)
        eIn(i) = scale[i] * def(i);

    return theSection->setTrialSectionDeformation(eIn);
}

const Vector&
LadrunoShellModifierSection::getSectionDeformation(void)
{
    return e;
}

const Vector&
LadrunoShellModifierSection::getStressResultant(void)
{
    const Vector &sIn = theSection->getStressResultant();

    for (int i = 0; i < 8; i++)
        sigma(i) = scale[i] * sIn(i);

    return sigma;
}

const Matrix&
LadrunoShellModifierSection::getSectionTangent(void)
{
    const Matrix &Din = theSection->getSectionTangent();

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            D(i, j) = scale[i] * scale[j] * Din(i, j);

    return D;
}

const Matrix&
LadrunoShellModifierSection::getInitialTangent(void)
{
    const Matrix &Din = theSection->getInitialTangent();

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            D(i, j) = scale[i] * scale[j] * Din(i, j);

    return D;
}

double
LadrunoShellModifierSection::getRho(void)
{
    return mods[8] * theSection->getRho();
}

int
LadrunoShellModifierSection::commitState(void)
{
    return theSection->commitState();
}

int
LadrunoShellModifierSection::revertToLastCommit(void)
{
    return theSection->revertToLastCommit();
}

int
LadrunoShellModifierSection::revertToStart(void)
{
    e.Zero();
    return theSection->revertToStart();
}

SectionForceDeformation*
LadrunoShellModifierSection::getCopy(void)
{
    LadrunoShellModifierSection *theCopy =
        new LadrunoShellModifierSection(this->getTag(), theSection, mods);

    if (theCopy == 0) {
        opserr << "LadrunoShellModifierSection::getCopy -- failed to allocate copy\n";
        exit(-1);
    }

    return theCopy;
}

const ID&
LadrunoShellModifierSection::getType(void)
{
    static bool initialized = false;
    if (!initialized) {
        code(0) = SECTION_RESPONSE_FXX;
        code(1) = SECTION_RESPONSE_FYY;
        code(2) = SECTION_RESPONSE_FXY;
        code(3) = SECTION_RESPONSE_MXX;
        code(4) = SECTION_RESPONSE_MYY;
        code(5) = SECTION_RESPONSE_MXY;
        code(6) = SECTION_RESPONSE_VXZ;
        code(7) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return code;
}

int
LadrunoShellModifierSection::getOrder(void) const
{
    return 8;
}

int
LadrunoShellModifierSection::sendSelf(int cTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();

    static Vector data(10); // tag + 9 modifiers
    data(0) = this->getTag();
    for (int i = 0; i < 9; i++)
        data(1 + i) = mods[i];

    if (theChannel.sendVector(dbTag, cTag, data) < 0) {
        opserr << "LadrunoShellModifierSection::sendSelf() - failed to send data\n";
        return -1;
    }

    // classTag/dbTag of the inner section (ParallelSection pattern)
    ID classTags(2);
    classTags(0) = theSection->getClassTag();
    int secDbTag = theSection->getDbTag();
    if (secDbTag == 0) {
        secDbTag = theChannel.getDbTag();
        if (secDbTag != 0)
            theSection->setDbTag(secDbTag);
    }
    classTags(1) = secDbTag;

    if (theChannel.sendID(dbTag, cTag, classTags) < 0) {
        opserr << "LadrunoShellModifierSection::sendSelf() - failed to send classTags\n";
        return -2;
    }

    if (theSection->sendSelf(cTag, theChannel) < 0) {
        opserr << "LadrunoShellModifierSection::sendSelf() - failed to send inner section\n";
        return -3;
    }

    return 0;
}

int
LadrunoShellModifierSection::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();

    static Vector data(10);
    if (theChannel.recvVector(dbTag, cTag, data) < 0) {
        opserr << "LadrunoShellModifierSection::recvSelf() - failed to receive data\n";
        return -1;
    }

    this->setTag((int)data(0));
    for (int i = 0; i < 9; i++)
        mods[i] = data(1 + i);

    ID classTags(2);
    if (theChannel.recvID(dbTag, cTag, classTags) < 0) {
        opserr << "LadrunoShellModifierSection::recvSelf() - failed to receive classTags\n";
        return -2;
    }

    int matClassTag = classTags(0);
    if (theSection == 0 || theSection->getClassTag() != matClassTag) {
        if (theSection != 0)
            delete theSection;

        SectionForceDeformation *theSectionModel = theBroker.getNewSection(matClassTag);
        if (theSectionModel != 0) {
            theSection = theSectionModel;
            theSection->setDbTag(classTags(1));
        } else {
            opserr << "FATAL LadrunoShellModifierSection::recvSelf() ";
            opserr << "could not get a SectionForceDeformation\n";
            return -3;
        }
    }

    if (theSection->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "LadrunoShellModifierSection::recvSelf() - failed to receive inner section\n";
        return -4;
    }

    this->buildScales();

    return 0;
}

void
LadrunoShellModifierSection::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "LadrunoShellModifierSection, tag: " << this->getTag() << endln;
        s << "  f11 = " << mods[0] << "  f22 = " << mods[1] << "  f12 = " << mods[2] << endln;
        s << "  m11 = " << mods[3] << "  m22 = " << mods[4] << "  m12 = " << mods[5] << endln;
        s << "  v13 = " << mods[6] << "  v23 = " << mods[7] << "  mass = " << mods[8] << endln;
        s << "  Inner section, tag: " << theSection->getTag() << endln;
        if (flag == OPS_PRINT_PRINTMODEL_MATERIAL)
            theSection->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"LadrunoShellModifierSection\", ";
        s << "\"f11\": " << mods[0] << ", \"f22\": " << mods[1] << ", \"f12\": " << mods[2] << ", ";
        s << "\"m11\": " << mods[3] << ", \"m22\": " << mods[4] << ", \"m12\": " << mods[5] << ", ";
        s << "\"v13\": " << mods[6] << ", \"v23\": " << mods[7] << ", \"mass\": " << mods[8] << ", ";
        s << "\"section\": \"" << theSection->getTag() << "\"}";
    }
}

int
LadrunoShellModifierSection::setParameter(const char **argv, int argc, Parameter &param)
{
    if (theSection != 0)
        return theSection->setParameter(argv, argc, param);

    return -1;
}
