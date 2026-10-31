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

// Ladruno: LadrunoHS -- Hardening Soil (Schanz, Vermeer & Bonnier 1999) and
// Hardening Soil with Small-Strain Stiffness (HSsmall, Benz 2007) nDMaterial
// with selectable Lode-angle cross-section (-lodeShape mc vs -lodeShape mn,
// Panteghini & Lagioia 2014).
//
// One self-contained class drives all dimensional views (ThreeDimensional,
// PlaneStrain, AxiSymmetric, PlaneStress, PlateFiber) via a `dim_` mode.
// Backed by the OpenSees-free C++17 constitutive kernel LadrunoHSKernel.h.
//
// See Ladruno_implementation/210_ladruno_hs_plan.md. classTag 33029.
// Written: N. Mora-Bowen (Ladruno), 2026.

#ifndef LadrunoHS_h
#define LadrunoHS_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <LadrunoHSKernel.h>

class LadrunoHS : public NDMaterial {
 public:
  enum {
    DIM_3D = 0,          // {00,11,22,01,12,02}      order 6
    DIM_PSTRAIN,         // {00,11,01}  (eps22=0)    order 3
    DIM_AXISYM,          // {00,11,22,01}            order 4
    DIM_PLATEFIBER,      // {00,11,01,12,20} sig22=0 order 5
    DIM_PSTRESS          // {00,11,01}  sig22=0      order 3
  };

  LadrunoHS();
  LadrunoHS(int tag,
            const LadrunoHSKernel::HSParams& params,
            double rho = 0.0,
            int dimMode = DIM_3D,
            int stage = 1,
            bool implex = false,
            bool hasSig0 = false,
            const double sig0Comp[6] = nullptr);
  ~LadrunoHS() override;

  const char* getClassType(void) const override { return "LadrunoHS"; }

  int setTrialStrain(const Vector& v) override;
  int setTrialStrain(const Vector& v, const Vector& r) override;
  int setTrialStrainIncr(const Vector& v) override;
  int setTrialStrainIncr(const Vector& v, const Vector& r) override;

  const Matrix& getTangent(void) override;
  const Matrix& getInitialTangent(void) override;
  const Vector& getStress(void) override;
  const Vector& getStrain(void) override;
  double getStressZZ(void) override;

  int commitState(void) override;
  int revertToLastCommit(void) override;
  int revertToStart(void) override;

  NDMaterial* getCopy(void) override;
  NDMaterial* getCopy(const char* type) override;
  const char* getType(void) const override;
  int getOrder(void) const override { return ncomp_; }

  double getRho(void) override { return rho_; }

  int sendSelf(int commitTag, Channel& theChannel) override;
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) override;
  void Print(OPS_Stream& s, int flag = 0) override;

  int setParameter(const char** argv, int argc, Parameter& param) override;
  int updateParameter(int parameterID, Information& info) override;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s) override;
  int getResponse(int responseID, Information& matInfo) override;

 private:
  void setupDim(void);
  void syncStageFromRegistry(void);
  void applyStageTransition(int newStage);
  int  integrate3D(const double epsComp[6]);
  void condenseTangent(double D[6][6]) const;
  void buildInitialTangent6x6(double Kt[6][6]) const;

  LadrunoHSKernel::HSParams params_;
  double rho_;
  int    dim_;
  int    ncomp_;
  int    vmap_[6];
  bool   condense_;
  int    stage_;       // 0 = elastic gravity stage; 1 = elastoplastic HS/HSsmall
  bool   implex_;      // Oliver et al. (2008) IMPL-EX explicit-implicit scheme
  bool   has_sig0_;
  double sig0_comp_[6];

  // Committed state & trial state (compression-positive inside kernel):
  LadrunoHSKernel::HSState st_n_;
  LadrunoHSKernel::HSState st_trial_;
  LadrunoHSKernel::ActiveRegime regime_n_;
  LadrunoHSKernel::ActiveRegime regime_trial_;
  double dsig_p_n_[6];        // committed plastic stress correction for IMPL-EX
  double dsig_p_trial_[6];    // trial plastic stress correction for IMPL-EX
  double sig_impl_trial_[6];  // implicit return-map trial stress for IMPL-EX commit
  double cEps22_comp_;        // committed out-of-plane compressive strain in condensed modes
  double Dtan_[6][6];         // 6x6 trial algorithmic tangent
  double Dtan_n_[6][6];       // 6x6 committed algorithmic tangent
  bool   refused_;            // WP-99 commit-refusal latch when setTrialStrain fails
  bool   registered_in_tag_map_;  // tracks ref-counted registration in ladrunoHSTagStageMap()

  // Element-facing return buffers (sized to ncomp_):
  Vector stressOut_;
  Vector strainOut_;
  Matrix tangentOut_;
  Matrix initTangentOut_;
};

void* OPS_LadrunoHS(void);
void* OPS_LadrunoHSSmall(void);

#endif  // LadrunoHS_h
