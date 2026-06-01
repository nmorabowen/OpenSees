/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// Description: LogStrainNDMaterial — the logarithmic (Hencky) strain-space
// finite-strain adaptor (seam 3). Lifts an UNCHANGED small-strain 3D NDMaterial
// to finite strain by the spatial multiplicative technique of de Souza Neto,
// Perić & Owen, *Computational Methods for Plasticity* (2008), Box 14.3 (MATISU):
//
//     bᵉᵗʳ = F_Δ bᵉ_n F_Δᵀ  → εᵉᵗʳ = ½ ln bᵉᵗʳ  → inner small-strain law
//        → Kirchhoff τ  → Cauchy σ = J⁻¹ τ                     (Box 14.3 iv)
//        → spatial CONSTITUTIVE modulus c = (1/2J)[D:L:B]      (§14.5)
//
// Implements the FiniteStrainNDMaterial contract: the element calls setTrialF(F)
// with the TOTAL deformation gradient; getStress() returns Cauchy σ and
// getTangent() returns the spatial CONSTITUTIVE modulus c = (1/2J)[D:L:B].
//
// SEAM-3 TANGENT CONTRACT (LOCKED 2026-06-01): getTangent() is the MATERIAL part
// only of the eq.(14.99) spatial tangent — the geometric / initial-stress term
// (−σ_il δ_jk) is NOT included; it is the ELEMENT's K_geo (∫ Gᵀ Σ G dv, built
// from the Cauchy stress and the element's shape-function gradients, which the
// material cannot know). The element forms K = ∫ Bᵀ c B dv + ∫ Gᵀ Σ G dv. This
// keeps the adaptor element-agnostic; see FiniteStrainNDMaterial.h for the full
// contract.
//
// The adaptor owns the committed elastic left Cauchy–Green tensor bᵉ_n, the
// committed F_n, and the committed fed strain ε_feed_n per Gauss point, and wraps
// the inner material's commit / revert / sendSelf / recvSelf.
//
// PLASTIC-INNER STATE PROTOCOL: a stock small-strain plastic NDMaterial treats
// setTrialStrain() input as TOTAL strain and subtracts its OWN committed εᵖ.
// Since the adaptor already carries the elastic history in bᵉ_n, feeding it the
// trial elastic strain directly would double-count εᵖ. The adaptor therefore
// feeds  ε_feed = ε_feed_n + (εᵉᵗʳ − εᵉ_n)  so the inner's subtraction yields
// exactly εᵉᵗʳ (invariant εᵖ_inner = ε_feed_n − εᵉ_n, maintained on commit), then
// recovers εᵉ_{n+1} = Cᵉ:τ (inner elastic compliance) to update bᵉ. This reuses an
// UNCHANGED stateful inner for BOTH elastic and rate-independent plastic laws
// (J2, Drucker–Prager, …) and reproduces de Souza Neto Box 14.3/14.4 exactly;
// hardening persists naturally inside the inner. v1 assumes a LINEAR-elastic inner
// law (so Cᵉ = inv(initial tangent) is constant and the εᵉ recovery is exact).
//
// Reference algorithm + protocol verified against the numpy oracle in
// tests/logstrain_reference.py and tests/test_logstrain_plastic_protocol.py.

#ifndef LogStrainNDMaterial_h
#define LogStrainNDMaterial_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LogStrainNDMaterial : public FiniteStrainNDMaterial
{
 public:
  LogStrainNDMaterial(int tag, NDMaterial &inner);
  LogStrainNDMaterial();
  ~LogStrainNDMaterial();

  const char *getClassType(void) const { return "LogStrainNDMaterial"; }

  // --- the finite-strain seam (FiniteStrainNDMaterial contract) ---
  int setTrialF(const Matrix &F);          // total deformation gradient in
  const Vector &getStress(void);           // Cauchy σ (6-vector)
  const Matrix &getTangent(void);          // spatial constitutive c (6×6, lossy proj.)
  int getSpatialTangentTensor(double c[3][3][3][3]);  // FULL c_ijkl (element channel)
  const Matrix &getInitialTangent(void);
  double getJ(void) const { return Jdet; }

  const Vector &getStrain(void);           // Hencky strain εᵉᵗʳ (6-vector)
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const { return "ThreeDimensional"; }
  int getOrder(void) const { return 6; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);

 private:
  NDMaterial *theMaterial;     // inner small-strain 3D material (deep-copied)

  // committed state (per Gauss point)
  double Fn[9];                // committed total deformation gradient (row-major)
  double Be_n[9];              // committed elastic left Cauchy–Green tensor bᵉ_n
  double epsFeed_n[6];         // committed strain last fed to the inner material
                               // (plastic-inner protocol; engineering-shear Voigt)

  // trial state staged for commit
  double Ftrial9[9];           // trial F
  double BeTrial9[9];          // trial Bᵉᵗʳ (used for the §14.5 spatial tangent)
  double Be_trialUpd[9];       // bᵉ to become bᵉ_n on commit = exp[2 εᵉ_{n+1}]
  double epsFeedTrial[6];      // trial fed strain, staged for commit

  // cached trial outputs
  Vector sigmaCauchy;          // Cauchy stress (6)
  Vector henckyStrain;         // εᵉᵗʳ (6, engineering shear)
  Matrix aTangent;             // spatial tangent (6×6)
  double Jdet;                 // det F

  void setIdentity(double M[9]);
};

#endif
