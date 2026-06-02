/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: combined-hardening von Mises (J2) nDMaterial.
//   - Isotropic: Voce saturation + linear   sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso*p
//   - Kinematic: Chaboche superposed Armstrong-Frederick, alpha = sum_k alpha_k
//                alpha_k_dot = (2/3) C_k eps_p_dot - gamma_k alpha_k pbar_dot
//                (gamma_k = 0 => linear Prager term; N=1 => single AF)
//   - Return map: implicit backward-Euler, reduced to a SINGLE scalar Newton on
//     dGamma (Kobayashi-Ohno 2002). Direction n = M(dG)/||M(dG)||,
//     M(dG) = s_tr - sum_k alpha_k,n/(1 + sqrt(2/3) gamma_k dG).
//   - Analytic consistent tangent (dSNPO 2008 eq 7.213 structure + AF dn/deps term).
//
// Supersedes J2Plasticity (no kinematic) and SimplifiedJ2 (linear-only + defects).
// v1: ThreeDimensional only. See Ladruno_implementation/10_ladruno_j2_plasticity.md
// classTag: ND_TAG_LadrunoJ2 = 33011.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoJ2_h
#define LadrunoJ2_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LadrunoJ2 : public NDMaterial {
 public:
  static const int MAXBACK = 8;   // max number of Chaboche backstress terms

  LadrunoJ2();
  LadrunoJ2(int tag, double K, double G,
            double sig0, double Qinf, double bIso, double Hiso,
            int nBack, const double* C, const double* gam,
            double rho = 0.0);
  ~LadrunoJ2();

  const char* getClassType(void) const { return "LadrunoJ2"; }

  // strain interface (engineering shear in/out, ordering {11,22,33,12,23,13})
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
  const char* getType(void) const { return "ThreeDimensional"; }
  int getOrder(void) const { return 6; }

  double getRho(void) { return rho; }

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag = 0);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int paramID);

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
  double yieldStress(double pbar) const;   // sig_y(pbar)
  double yieldSlope(double pbar) const;     // d sig_y / d pbar
  void   integrate(void);                   // return map + tangent
  void   buildElasticTangent(double Kt[6][6]) const;
  void   setStateToCommitted(void);

  int parameterID;

  // return buffers
  static Vector stressV;
  static Vector strainV;
  static Matrix tangentM;
};

#endif
