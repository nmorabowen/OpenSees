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

// Ladruno WP-254 (program #1131): Model C's reinforcing-steel law as a uniaxial material.
//
//   uniaxialMaterial LadrunoMCFTSteel $tag $fy $fu [-Es $Es] [-esh $esh] [-esu $esu]
//
// The bar law of the Model C shear-critical fibre beam (and the tie law inside
// LadrunoMCFTStrip; both use ladruno_mcft::steelStress):
//
//   backbone env(a), a = strain magnitude:
//     Es a                                   a <= fy/Es
//     fy                                     a <= esh
//     fy + (fu - fy)(1 - (1 - r)^2)          a <= esu,  r = (a - esh)/(esu - esh)
//     1e-3 fy                                a >  esu   (fracture, by TOTAL strain)
//   sigma = clip(Es (eps - ep), -env(max(-eps, fy/Es)), env(max(eps, fy/Es)))
//   committed plastic strain ep = eps - sigma/Es
//
// Elastic unloading and reloading bounded by the backbone of each sign at the
// current total strain: no Bauschinger effect, no buckling. Defaults Es 200000,
// esh 0.008, esu 0.12 (MPa units). Put it in an NDFiber section through
// `nDMaterial LadrunoUniaxialFiber`.

#ifndef LadrunoMCFTSteel_h
#define LadrunoMCFTSteel_h

#include <UniaxialMaterial.h>
#include "../nD/LadrunoMCFTKernel.h"

class LadrunoMCFTSteel : public UniaxialMaterial
{
  public:
    LadrunoMCFTSteel(int tag, const ladruno_mcft::Steel& s);
    LadrunoMCFTSteel();
    ~LadrunoMCFTSteel();

    const char *getClassType(void) const {return "LadrunoMCFTSteel";};

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void) {return epsT;}
    double getStress(void) {return sigT;}
    double getTangent(void) {return EtT;}
    double getInitialTangent(void) {return steel.Es;}

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

  private:
    ladruno_mcft::Steel steel;
    double epC, epsC;        // committed plastic strain, committed strain
    double epsT, sigT, EtT;  // trial
};

#endif
