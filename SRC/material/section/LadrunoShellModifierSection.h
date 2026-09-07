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
//
//   section LadrunoShellModifier tag? innerSecTag?
//       <-f11 v> <-f22 v> <-f12 v> <-m11 v> <-m22 v> <-m12 v>
//       <-v13 v> <-v23 v> <-mass v>
//
// LadrunoShellModifierSection wraps ANY order-8 plate/shell
// SectionForceDeformation (response codes exactly
// {FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ}, e.g. ElasticMembranePlateSection or
// MembranePlateFiberSection) and applies nine independent ETABS-style
// stiffness modifiers, all defaulting to 1.0:
//
//   m[0]=f11  m[1]=f22  m[2]=f12   membrane modifiers
//   m[3]=m11  m[4]=m22  m[5]=m12   bending modifiers
//   m[6]=v13  m[7]=v23             transverse-shear modifiers
//   m[8]=mass                      scales getRho() only
//
// The transform is the SPD-preserving congruence D' = S*D*S with
// S = diag(sqrt(m[0..7])) applied to the inner section's strain, stress and
// tangent (both consistent and initial). This is exact and energy-consistent
// for an arbitrary (possibly nonlinear) inner section and is NOT a diagonal
// rescale of the tangent -- a diagonal-only scaling can drive the membrane
// block indefinite when a small f11/f22 is combined with an unscaled
// Poisson coupling term. There is deliberately NO weight/aggregate modifier.
//
// Written: N. Mora-Bowen (Ladruno), 2026.
// See Ladruno_implementation/91_ladruno_shell_stiffness_modifiers_adr.md

#ifndef LadrunoShellModifierSection_h
#define LadrunoShellModifierSection_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

class LadrunoShellModifierSection : public SectionForceDeformation
{
  public:
    // constructor for blank object that recvSelf needs to be invoked upon
    LadrunoShellModifierSection();

    // mods[9] = {f11,f22,f12,m11,m22,m12,v13,v23,mass}; theInnerSection is
    // NOT taken over -- a deep copy (theInnerSection->getCopy()) is owned.
    LadrunoShellModifierSection(int tag, SectionForceDeformation *theInnerSection,
                                 const double mods[9]);

    ~LadrunoShellModifierSection();

    const char *getClassType(void) const {return "LadrunoShellModifierSection";};

    int   setTrialSectionDeformation(const Vector &deforms);
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);

    double getRho(void);

    int   commitState(void);
    int   revertToLastCommit(void);
    int   revertToStart(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    // Forwarded to the inner section only -- the nine modifiers themselves
    // are deliberately NOT exposed as settable Parameters (out of scope).
    int setParameter(const char **argv, int argc, Parameter &param);

    // Ladruno NOTE (WP-91): DDM sensitivity is NOT implemented for this
    // wrapper -- getStressResultantSensitivity / getSectionDeformationSensitivity
    // / getSectionTangentSensitivity / getRhoSensitivity / commitSensitivity
    // all fall through to the (non-pure) SectionForceDeformation base-class
    // defaults. Wiring the congruence transform through DDM is future work.

  private:
    // sqrt(mods[0..7]) into scale[], plus one-time opserr warnings for any
    // modifier (including mass) that is exactly 0.0 (allowed, but singular
    // in that response mode).
    void buildScales(void);

    SectionForceDeformation *theSection; // owned deep copy of the inner order-8 section

    double mods[9];   // f11 f22 f12 m11 m22 m12 v13 v23 mass -- all default 1.0
    double scale[8];  // sqrt(mods[0..7]), congruence factors on strain/stress

    Vector e;             // outer (un-scaled) trial section deformation, order 8
    static Vector sigma;  // scaled stress-resultant work vector, order 8
    static Matrix D;      // scaled tangent work matrix, order 8x8
    static ID code;       // section response codes {FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ}

    int otherDbTag; // dbTag reserved for the shipped inner section (ParallelSection pattern)
};

#endif
