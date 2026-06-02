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

// Ladruno: combined-hardening von Mises (J2) nDMaterial.
//   - Isotropic: Voce saturation + linear   sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso*p
//   - Kinematic: Chaboche superposed Armstrong-Frederick, alpha = sum_k alpha_k
//                (gamma_k = 0 => linear Prager term; N=1 => single AF)
//   - Return map: implicit backward-Euler, reduced to a SINGLE scalar Newton on
//     dGamma (Kobayashi-Ohno 2002). Analytic consistent tangent (dSNPO 2008
//     eq 7.213 structure + AF dn/deps term).
//
// One self-contained class drives every dimensional view via a `dim` mode: the
// 3D return map runs on the full 6-component tensor and the strain/stress/tangent
// are mapped to the element's reduced ordering. PlaneStress/PlateFiber enforce the
// out-of-plane sigma_22 = 0 by a nested Newton on eps_22 + static condensation of
// the 33 dof (dSNPO 2008 sec 9.2.3 nested-iteration route).
//
// Supersedes J2Plasticity and SimplifiedJ2.
// See Ladruno_implementation/10_ladruno_j2_plasticity.md. classTag 33011.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoJ2_h
#define LadrunoJ2_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LadrunoJ2 : public NDMaterial {
 public:
  static const int MAXBACK = 8;   // max number of Chaboche backstress terms

  // dimensional views (element-facing ordering, engineering shear)
  enum { DIM_3D = 0,         // {00,11,22,01,12,02}      order 6
         DIM_PSTRAIN,        // {00,11,01}  (eps22=0)    order 3
         DIM_AXISYM,         // {00,11,22,01}            order 4
         DIM_PLATEFIBER,     // {00,11,01,12,20} sig22=0 order 5
         DIM_PSTRESS };      // {00,11,01}  sig22=0      order 3

  LadrunoJ2();
  LadrunoJ2(int tag, double K, double G,
            double sig0, double Qinf, double bIso, double Hiso,
            int nBack, const double* C, const double* gam,
            double rho = 0.0, int dimMode = DIM_3D);
  ~LadrunoJ2();

  const char* getClassType(void) const { return "LadrunoJ2"; }

  int setTrialStrain(const Vector& strain);
  int setTrialStrain(const Vector& v, const Vector& r);
  int setTrialStrainIncr(const Vector& v);
  int setTrialStrainIncr(const Vector& v, const Vector& r);

  const Matrix& getTangent(void);
  const Matrix& getInitialTangent(void);
  const Vector& getStress(void);
  const Vector& getStrain(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial* getCopy(void);
  NDMaterial* getCopy(const char* type);
  const char* getType(void) const;
  int getOrder(void) const;

  double getRho(void) { return rho; }

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag = 0);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int paramID);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  // material parameters
  double bulk;     // K
  double shear;    // G
  double sig0;     // initial yield stress
  double Qinf;     // saturation increment (Voce)
  double bIso;     // saturation rate (Voce)
  double Hiso;     // linear isotropic modulus
  double rho;      // mass density
  int    nBack;                 // number of backstress terms
  double Ckin[MAXBACK];         // AF kinematic moduli C_k
  double gKin[MAXBACK];         // AF recall constants gamma_k

  // dimensional view
  int    dim;                   // DIM_*
  int    ncomp;                 // element vector order (3..6)
  int    vmap[6];               // reduced index a -> full 6-comp index
  bool   condense;              // enforce sigma_22 = 0 (PlaneStress / PlateFiber)
  double cEps22;                // committed out-of-plane strain (condensed modes)

  // committed history (symmetric tensors stored as 6 tensor components:
  //   {t00,t11,t22,t01,t12,t02}; shear entries are TRUE tensor components)
  double epsP_n[6];             // plastic strain (deviatoric)
  double ebarP_n;               // accumulated equivalent plastic strain
  double alpha_n[MAXBACK][6];   // backstress terms
  double dGamma_n;              // committed plastic multiplier increment (IMPL-EX hook)

  // trial state
  double strain6[6];            // total strain (tensor components)
  double epsP[6];
  double ebarP;
  double alpha[MAXBACK][6];
  double dGammaTrial;
  double stress6[6];            // stress (tensor comps; shear = true stress)
  double Dtan[6][6];            // algorithmic tangent (engineering 6x6, J2-3D convention)

  // helpers
  void   integrate(void);                   // 3D return map + tangent (dim-agnostic)
  void   buildElasticTangent(double Kt[6][6]) const;
  void   setStateToCommitted(void);
  void   setupDim(void);                    // fill ncomp/vmap/condense + size buffers
  void   condenseTangent(void);             // static condensation of the 33 dof

  int parameterID;

  // element-facing return buffers (sized to ncomp)
  Vector stressOut;
  Vector strainOut;
  Matrix tangentOut;
};

#endif
