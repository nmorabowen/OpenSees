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

// Ladruno WP-254 (program #1131): the Model C MCFT strip as an NDFiber beam fibre.
//
//   nDMaterial LadrunoMCFTStrip $tag $fc $fp $ep $e20 -crackSpacing $smx $smy
//       [-ties $rho $fyt $fut] [-tieSteel $Es $esh $esu]
//       [-agg $a] [-betaC $C] [-cvci $c] [-Gxz $G] [-tangentStep $h]
//
// A concrete strip of a shear-critical fibre beam (Sezen & Moehle 2002; Lynn et
// al. 1996 calibration set): a 2D rotating MCFT membrane in the plane of
// bending with ONE internal unknown, the transverse strain eps_y, balanced
// against smeared ties (the law and its root rule: LadrunoMCFTKernel.h).
//
// Views (getCopy):
//   BeamFiber   (order 3): strain (eps_x, gamma_xy, gamma_xz) -> (sx, txy, txz)
//   BeamFiber2d (order 2): strain (eps_x, gamma_xy)           -> (sx, txy)
// Every other view is refused.
//
// gamma_xz. The strip is a membrane in x-y; it has no out-of-plane law. In the
// 3D view txz = Gxz gamma_xz, elastic and uncoupled (default Gxz = Ec/2, the
// membrane's own initial shear modulus; zero Poisson). It gives a 3D section
// its weak-axis shear and torsion stiffness; it carries no damage. Load the
// member in its x-y plane.
//
// Tangent. The condensed strip tangent d(sx, txy)/d(eps_x, gamma_xy) by forward
// differences, each perturbed state re-solving eps_y from the committed memory
// (Model C's own section tangent, done per fibre). Step -tangentStep (default
// 1e-7) on both strains. getInitialTangent is the uncracked one: diag(Ec, Ec/2,
// Gxz).
//
// State. Committed only: eps_y (the root reference), max principal tension and
// compression, tie plastic strain, and a first-commit flag. A trial state is a
// pure function of the committed memory and the trial strain.

#ifndef LadrunoMCFTStrip_h
#define LadrunoMCFTStrip_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include "LadrunoMCFTKernel.h"

class LadrunoMCFTStrip : public NDMaterial
{
  public:
    LadrunoMCFTStrip(int tag, const ladruno_mcft::Strip& strip, double Gxz, double hTan, int order = 3);
    LadrunoMCFTStrip();
    ~LadrunoMCFTStrip();

    const char *getClassType(void) const {return "LadrunoMCFTStrip";};

    int setTrialStrain(const Vector& strain);
    int setTrialStrain(const Vector& strain, const Vector& rate);
    int setTrialStrainIncr(const Vector& strain);
    int setTrialStrainIncr(const Vector& strain, const Vector& rate);
    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);
    const Vector& getStress(void);
    const Vector& getStrain(void);
    double getRho(void) {return 0.0;}

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial* getCopy(void);
    NDMaterial* getCopy(const char* type);
    const char* getType(void) const;
    int getOrder(void) const;

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& matInfo);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    void Print(OPS_Stream& s, int flag = 0);

  private:
    void evaluateTrial(void);

    ladruno_mcft::Strip strip;
    double Gxz;          // elastic out-of-plane shear modulus (order 3)
    double hTan;         // forward-difference step of the tangent
    int order;           // 3 = BeamFiber, 2 = BeamFiber2d

    ladruno_mcft::Hist hC;      // committed memory
    double eyT;                 // trial eps_y
    ladruno_mcft::Out oT;       // trial membrane state
    int foundT;                 // 1: root bracketed near the committed eps_y, 0: full range
    bool tanValid;              // D is current for the trial strain

    Vector eps, epsC, sig;
    Matrix D, Dinit;
};

#endif
