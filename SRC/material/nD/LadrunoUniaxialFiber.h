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

// Ladruno WP-245 (program #1131, rcbf-3): a UniaxialMaterial as an NDFiber
// beam fibre.
//
//   nDMaterial LadrunoUniaxialFiber $tag $uniTag [-rho $rho]
//
// NDFiber sections take only nD materials, so a reinforcing bar (Steel02,
// ReinforcingSteel, LadrunoRebarBuckling, ...) cannot sit in the same section
// as the concrete. This adapter is the bar: the axial strain drives the
// uniaxial law, and the bar carries NO shear (no dowel action, as Model C,
// 06 §1.4):
//
//   BeamFiber   (order 3): strain (eps, gxy, gxz) -> stress (sig(eps), 0, 0),
//                          tangent diag(E_t, 0, 0)
//   BeamFiber2d (order 2): strain (eps, gxy)      -> stress (sig(eps), 0)
//
// Any other view (ThreeDimensional, PlaneStress, PlateFiber, ...) is refused:
// a bar has no meaning there. The fibre's area should be the bar area; the
// concrete fibre it sits in should be net of it (Model C subtracts each bar's
// area from the nearest strip).

#ifndef LadrunoUniaxialFiber_h
#define LadrunoUniaxialFiber_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;

class LadrunoUniaxialFiber : public NDMaterial
{
  public:
    LadrunoUniaxialFiber(int tag, UniaxialMaterial& theUni, int order = 3, double rho = -1.0);
    LadrunoUniaxialFiber();
    ~LadrunoUniaxialFiber();

    const char *getClassType(void) const {return "LadrunoUniaxialFiber";};

    int setTrialStrain(const Vector& strain);
    int setTrialStrain(const Vector& strain, const Vector& rate);
    int setTrialStrainIncr(const Vector& strain);
    int setTrialStrainIncr(const Vector& strain, const Vector& rate);
    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);
    const Vector& getStress(void);
    const Vector& getStrain(void);
    double getRho(void);

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
    UniaxialMaterial* theUni;
    int order;          // 3 = BeamFiber, 2 = BeamFiber2d
    double rho;         // < 0: use the uniaxial material's own getRho()
    Vector eps3, sig3, eps2, sig2;
    Matrix D3, D2;
};

#endif
