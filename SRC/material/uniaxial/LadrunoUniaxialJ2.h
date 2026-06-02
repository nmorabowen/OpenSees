/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Ladruno: uniaxial combined-hardening von Mises (J2) UniaxialMaterial.
//   sigma = E (eps - eps^p)
//   yield: |sigma - X| - sig_y(pbar) <= 0,  X = sum_k X_k  (Chaboche backstress)
//   - Isotropic: Voce saturation + linear  sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso p
//                (shared verbatim with the 3D nDMaterial LadrunoJ2 via
//                 LadrunoHardening.h -- the V7 1e-12 oracle contract)
//   - Kinematic: Chaboche superposed Armstrong-Frederick, X = sum_k X_k,
//                X_k <- (X_k + C_k dp s)/(1 + gamma_k dp)   (gamma_k=0 => linear
//                Prager term; N=1 => single AF; saturates to X_k,sat = C_k/gamma_k)
//   - Return map: implicit backward-Euler. In 1D the flow direction
//                 s = sign(sigma_tr - X_n) is FIXED by the trial state (the
//                 relative stress cannot rotate on a line), so consistency is a
//                 single scalar Newton on dp = accumulated-plastic-strain
//                 increment -- no M(dGamma) tensor-direction solve like the 3D
//                 kernel. Analytic consistent tangent E_alg = E h/(E+h).
//
// This is the uniaxial twin of LadrunoJ2 (nDMaterial). Its primary roles are
// (1) a verification oracle for the 3D Chaboche calibration under uniaxial
// stress, and (2) a fiber material with true multi-backstress AF ratcheting --
// NOT a Steel02/Steel4 replacement. It is the core onto which future steel
// layers (fracture/LCF wrapper, local-buckling flag, weld/HAZ) bolt; v1 exposes
// plasticStrain / plasticWork / backStress responses for that.
//
// classTag MAT_TAG_LadrunoUniaxialJ2 = 33000 (Ladruno private uniaxial band).
// See Ladruno_implementation/13_ladruno_uniaxial_j2_adr.md.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoUniaxialJ2_h
#define LadrunoUniaxialJ2_h

#include <UniaxialMaterial.h>

class LadrunoUniaxialJ2 : public UniaxialMaterial
{
 public:
  static const int MAXBACK = 8;   // max number of Chaboche backstress terms

  LadrunoUniaxialJ2(int tag, double E,
                    double sig0, double Qinf, double bIso, double Hiso,
                    int nBack, const double* C, const double* gam);
  LadrunoUniaxialJ2();
  ~LadrunoUniaxialJ2();

  const char* getClassType(void) const { return "LadrunoUniaxialJ2"; }

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void)          { return Tstrain; }
  double getStress(void)          { return Tstress; }
  double getTangent(void)         { return Ttangent; }
  double getInitialTangent(void)  { return E; }    // elastic — correct for modal

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial* getCopy(void);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int parameterID);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  // material parameters
  double E;                     // Young's modulus
  double sig0;                  // initial yield stress
  double Qinf;                  // Voce saturation increment
  double bIso;                  // Voce saturation rate
  double Hiso;                  // linear isotropic modulus
  int    nBack;                 // number of backstress terms
  double Ckin[MAXBACK];         // AF kinematic moduli C_k
  double gKin[MAXBACK];         // AF recall constants gamma_k

  // committed history
  double CplasticStrain;        // plastic strain eps^p
  double Cebarp;                // accumulated plastic strain pbar
  double CX[MAXBACK];           // backstresses X_k
  double CdGamma;               // last plastic-strain increment dp (IMPL-EX hook)
  double Cwp;                   // accumulated plastic work (fracture/LCF hook)

  // trial history
  double TplasticStrain;
  double Tebarp;
  double TX[MAXBACK];
  double TdGamma;
  double Twp;

  // trial state
  double Tstrain;
  double Tstress;
  double Ttangent;

  // helpers
  double yieldStress(double p) const;   // sig_y(p)   (shared Ladruno law)
  double yieldSlope(double p) const;    // d sig_y/dp (shared Ladruno law)

  int parameterID;
};

#endif
