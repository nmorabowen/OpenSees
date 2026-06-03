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

// Ladruno: reinforced-concrete plastic-damage nDMaterial (Phase 1 of the RC
// shell stack, ADR 19_ladruno_rc_shell_adr.md). Thin OpenSees-facing view over
// the header-only LadrunoRCKernel.h (the single source of the constitutive
// update, shared "one core, many views" exactly like LadrunoJ2 + LadrunoJ2Kernel).
//
// One self-contained class drives every dimensional view via a `dim` mode: the
// 3D return runs on the full 6-component tensor; PlaneStress/PlateFiber enforce
// the out-of-plane sigma_33 = 0 by a guarded nested Newton on eps_33 + static
// condensation of the 33 dof. The shell path is LayeredShellFiberSection ->
// getCopy("PlateFiber") (order 5).
//
// Keeps ASDConcrete3D's plastic-damage spine; adds MCFT compression softening
// on the STRENGTH axis (beta = 1/(0.8+170 eps_1)).  betaOn=false => byte-faithful
// to ASDConcrete3D.  ENGINEERING-shear convention (matches ASDConcrete3D), so the
// view does NO tensor<->engineering factor conversions.  classTag 33014.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoRCConcrete_h
#define LadrunoRCConcrete_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include "LadrunoRCKernel.h"

class LadrunoRCConcrete : public NDMaterial {
 public:
  // dimensional views (element-facing ordering, engineering shear)
  enum { DIM_3D = 0,        // {00,11,22,01,12,02}      order 6
         DIM_PSTRESS,       // {00,11,01}  sig33=0      order 3
         DIM_PLATEFIBER };  // {00,11,01,12,02} sig33=0 order 5

  LadrunoRCConcrete();
  LadrunoRCConcrete(int tag, const ladruno_rc_kernel::Params& P, double rho, int dimMode = DIM_3D);
  ~LadrunoRCConcrete();

  const char* getClassType(void) const { return "LadrunoRCConcrete"; }

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

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& matInfo);

 private:
  ladruno_rc_kernel::Params P;       // material constants + backbones (kernel-side)
  ladruno_rc_kernel::RCHist histN;   // committed history
  ladruno_rc_kernel::RCHist histTr;  // trial history (overwritten each integrate)
  double rho;

  // dimensional view
  int    dim;
  int    ncomp;
  int    vmap[6];
  bool   condense;
  double cEps33;                     // committed out-of-plane strain (condensed modes)

  // working state (engineering shear; stress = true stress)
  double strain6[6];
  double stress6[6];
  double Dtan[6][6];
  int    status;

  // helpers
  void setupDim(void);
  void integrate(void);
  void condenseTangent(void);

  // element-facing return buffers
  Vector stressOut;
  Vector strainOut;
  Matrix tangentOut;
};

#endif
