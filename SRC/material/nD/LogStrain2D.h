/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// Description: LogStrain2D — the 2D (plane) logarithmic (Hencky) finite-strain
// adaptor, the plane sibling of LogStrainNDMaterial. It presents an order-3
// PlaneStrain / PlaneStress face over an UNCHANGED small-strain 3D NDMaterial by
// lifting the element's 2×2 deformation gradient to 3×3 and delegating to an
// internally-composed LogStrainNDMaterial (which owns the full 3D bᵉ state and the
// verified MATISU spatial-multiplicative kernel):
//
//   PLANE STRAIN  — F₃₃ = 1 (no out-of-plane shear); the 3D adaptor is driven
//   verbatim and the in-plane block of σ / c is read back. Because the FULL 3×3
//   elastic left Cauchy–Green tensor is carried, a plastic out-of-plane strain
//   (ε₃₃ ≠ 0 from deviatoric J2 flow) is tracked exactly — plane strain is
//   genuinely the 3D problem with F₃₃ = 1.
//
//   PLANE STRESS  — the thickness stretch λ = F₃₃ is unknown and is solved from
//   σ₃₃(λ) = 0 by a local Newton (finite-difference dσ₃₃/dλ; warm-started from the
//   committed λ). The out-of-plane direction is then statically condensed out of
//   the spatial tangent (δσ₃₃ = 0) to give the consistent 2D modulus. de Souza
//   Neto §14.7. For a linear-elastic inner the condition is linear, so λ recovers
//   the classical thickness contraction ε₃₃ = −ν/(1−ν)(ε₁₁+ε₂₂).
//
// SEAM contract (FiniteStrainND2DMaterial): getStress() = in-plane Cauchy σ
// (3-vector {11,22,12}, engineering shear); getSpatialTangentTensor2D(c) = the
// full in-plane modulus (element forms a = c − σ_il δ_jk); getTangent() = the
// lossy 3×3 projection (output only). Inner MUST be a 3D order-6 material.

#ifndef LogStrain2D_h
#define LogStrain2D_h

#include <FiniteStrainND2DMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LogStrainNDMaterial;

class LogStrain2D : public FiniteStrainND2DMaterial
{
 public:
  enum PlaneType { PlaneStrain = 1, PlaneStress = 2 };

  LogStrain2D(int tag, NDMaterial &inner3D, int planeType);
  LogStrain2D();
  ~LogStrain2D();

  const char *getClassType(void) const { return "LogStrain2D"; }

  // --- the 2D finite-strain seam (FiniteStrainND2DMaterial contract) ---
  int setTrialF(const Matrix &F);              // total 2×2 deformation gradient
  const Vector &getStress(void);               // in-plane Cauchy σ (3-vector)
  const Matrix &getTangent(void);              // in-plane spatial c (3×3, lossy proj.)
  int getSpatialTangentTensor2D(double c[2][2][2][2]);  // FULL in-plane c_ijkl
  const Matrix &getInitialTangent(void);
  double getJ(void) const { return Jarea; }

  const Vector &getStrain(void);               // in-plane Hencky strain (3-vector)
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const {
    return (planeType == PlaneStress) ? "PlaneStress" : "PlaneStrain";
  }
  int getOrder(void) const { return 3; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

 private:
  // internally-composed 3D Hencky adaptor (owns the deep-copied 3D inner + bᵉ_n)
  LogStrainNDMaterial *the3D;
  int planeType;                 // PlaneStrain (1) or PlaneStress (2)

  // committed state
  double lam_n;                  // committed thickness stretch F₃₃ (plane stress)
  // trial state staged for commit
  double lamTrial;

  // cached trial outputs (in-plane, order 3)
  Vector sigma3;                 // in-plane Cauchy σ  {σ11,σ22,σ12}
  Vector hencky3;                // in-plane Hencky strain {ε11,ε22,γ12}
  Matrix tangent3;               // in-plane spatial tangent (3×3, lossy)
  double c2[2][2][2][2];         // FULL in-plane spatial modulus (element channel)
  double Jarea;                  // det(2×2 F)

  // helpers
  int  driveLifted(const Matrix &F2, double F33);   // build F3D, call the3D->setTrialF
  int  buildOutputs(const Matrix &F2, bool condense); // fill sigma3/tangent3/c2; <0 = refuse
  static void voigt3FromTensor(const double c[2][2][2][2], Matrix &D);
};

#endif
