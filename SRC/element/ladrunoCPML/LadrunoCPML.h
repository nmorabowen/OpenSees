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

#ifndef LadrunoCPML_h
#define LadrunoCPML_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Node.h>
#include <vector>

#include "LadrunoCPMLKernel.h"

#ifndef ELE_TAG_LadrunoCPML
#define ELE_TAG_LadrunoCPML 33023
#endif

class LadrunoCPML : public Element {
public:
  LadrunoCPML(int tag,
              ladruno_cpml::CPMLTopology topo,
              const ID& nodeTags,
              bool lumped,
              const ladruno_cpml::CPMLMaterialParams& mat,
              const ladruno_cpml::CPMLGradingParams& grading,
              const ladruno_cpml::CPMLGeomConfig& geom,
              int stage = 1,
              double thick = 1.0,
              const ladruno_cpml::CPMLBiotParams& biot = ladruno_cpml::CPMLBiotParams{});
  LadrunoCPML(int tag,
              int ndm,
              const ID& nodeTags,
              int shapeOverride,
              double E,
              double nu,
              double rho,
              double thickness,
              int geomMode,
              const Vector& geomParams,
              double mDeg,
              double R0,
              double kappaMax,
              double alphaMax,
              bool lumped,
              int initialStage,
              const Vector& bodyForces);
  LadrunoCPML(int tag,
              int ndm,
              const ID& nodeTags,
              int shapeOverride,
              double E,
              double nu,
              double rho,
              double thickness,
              int geomMode,
              const Vector& geomParams,
              double mDeg,
              double R0,
              double kappaMax,
              double alphaMax,
              bool lumped,
              int initialStage,
              const Vector& bodyForces,
              const ladruno_cpml::CPMLBiotParams& biot,
              const Vector& fluidBodyForces);
  LadrunoCPML();
  ~LadrunoCPML() override;

  const char* getClassType() const override { return "LadrunoCPML"; }
  static constexpr const char* CLASS_TYPE = "LadrunoCPML";

  int getNumExternalNodes() const override;
  const ID& getExternalNodes() override;
  Node** getNodePtrs() override;
  int getNumDOF() override;
  void setDomain(Domain* theDomain) override;

  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;
  int update() override;

  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc) override;

  const Matrix& getTangentStiff() override;
  const Matrix& getInitialStiff() override;
  const Matrix& getDamp() override;
  const Matrix& getMass() override;

  void zeroLoad() override;
  int addLoad(ElementalLoad* theLoad, double loadFactor) override;
  int addInertiaLoadToUnbalance(const Vector& accel) override;

  const Vector& getResistingForce() override;
  const Vector& getResistingForceIncInertia() override;

  int sendSelf(int commitTag, Channel& theChannel) override;
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) override;
  void Print(OPS_Stream& s, int flag = 0) override;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s) override;
  int getResponse(int responseID, Information& eleInfo) override;

  int setParameter(const char** argv, int argc, Parameter& param) override;
  int updateParameter(int parameterID, Information& info) override;

  bool isValid() const { return valid_; }

protected:
  const Vector& getRayleighDampingForces();

private:
  void ensureBuffers();
  void syncGeomParamsFromGeom();
  void gatherKinematics(bool useCommitted,
                        double* uOut,
                        double* vOut,
                        double* aOut) const;
  double scaleFactor() const { return (ndm_ == 2) ? thick_ : 1.0; }

  ladruno_cpml::CPMLTopology topo_ = ladruno_cpml::CPMLTopology::Q4;
  bool lumped_ = false;
  int stage_ = 1; // 0 = static elasticity only, 1 = active C-PML
  int initialStage_ = 1;
  int shapeOverride_ = 0;
  int geomMode_ = 0;
  Vector geomParams_;
  double thick_ = 1.0;
  double charLen_ = 0.0;
  int ndm_ = 2;
  int nNodes_ = 4;
  int nDOF_ = 8;
  bool valid_ = false;
  double lastDt_ = 0.0;
  bool skipNextUpdate_ = false;

  double b_[3] = {0.0, 0.0, 0.0};
  double appliedB_[3] = {0.0, 0.0, 0.0};
  double fluidBody_[3] = {0.0, 0.0, 0.0};
  double appliedFluidBody_[3] = {0.0, 0.0, 0.0};
  int applyLoad_ = 0;

  ladruno_cpml::CPMLMaterialParams mat_{};
  ladruno_cpml::CPMLGradingParams grading_{};
  ladruno_cpml::CPMLGeomConfig geom_{};
  ladruno_cpml::CPMLBiotParams biot_{};
  ladruno_cpml::CPMLElementKernel kernel_{};

  ID nodeTags_;
  std::vector<Node*> theNodes_;
  std::vector<double> coords_;
  std::vector<double> nodalVols_;
  std::vector<double> fseep_;

  // Stage 0 -> 1 static offset latching: u_0 and f_0 = K_0 * u_0
  std::vector<double> u0_;
  std::vector<double> f0_;
  bool u0Latched_ = false;

  // Trial kinematics buffers
  std::vector<double> uTrial_;
  std::vector<double> vTrial_;
  std::vector<double> aTrial_;

  // OpenSees Matrix/Vector workspaces
  Matrix Kmat_;
  Matrix Cmat_;
  Matrix Mmat_;
  bool MmatValid_ = false;
  int MmatStage_ = -1;
  Matrix* Ki_ = nullptr;
  Vector resid_;
  Vector residInc_;
  Vector basicForce_;
  Vector Qload_;
  Vector raScratch_;
  Vector rayForce_;
};

#endif // LadrunoCPML_h
