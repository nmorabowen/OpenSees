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

#ifndef LadrunoMembrane_h
#define LadrunoMembrane_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <LadrunoMassCache.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class Node;
class Domain;
class Channel;
class FEM_ObjectBroker;
class NDMaterial;
class ElementalLoad;
class Response;
class Information;
class Parameter;
class OPS_Stream;

class LadrunoMembrane : public Element
{
public:
  LadrunoMembrane(int tag,
                  const std::vector<int> &polyNodes,
                  const std::vector<int> &frameTags,
                  const std::vector<int> &supportNodes,
                  const std::vector<std::vector<int>> &holeLoops,
                  const std::string &method,
                  const std::array<double, 3> &spanAxis,
                  bool hasSpanAxis,
                  double orthoRatio,
                  int segments,
                  double meshSize,
                  double thick,
                  NDMaterial *theMat,
                  double E_direct,
                  double nu_direct,
                  double rho,
                  bool lumped,
                  double pressure,
                  double b1,
                  double b2,
                  double b3);

  LadrunoMembrane();
  ~LadrunoMembrane() override;

  const char *getClassType(void) const override { return "LadrunoMembrane"; }

  int getNumExternalNodes(void) const override;
  const ID &getExternalNodes(void) override;
  Node **getNodePtrs(void) override;
  int getNumDOF(void) override;
  void setDomain(Domain *theDomain) override;

  int commitState(void) override;
  int revertToLastCommit(void) override;
  int revertToStart(void) override;
  int update(void) override;

  const Matrix &getTangentStiff(void) override;
  const Matrix &getInitialStiff(void) override;
  const Matrix &getMass(void) override;

  void zeroLoad(void) override;
  int addLoad(ElementalLoad *theLoad, double loadFactor) override;
  int addInertiaLoadToUnbalance(const Vector &accel) override;

  const Vector &getResistingForce(void) override;
  const Vector &getResistingForceIncInertia(void) override;

  int sendSelf(int commitTag, Channel &theChannel) override;
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;
  void Print(OPS_Stream &s, int flag = 0) override;

  Response *setResponse(const char **argv, int argc, OPS_Stream &s) override;
  int getResponse(int responseID, Information &eleInformation) override;

  int setParameter(const char **argv, int argc, Parameter &param) override;
  int updateParameter(int parameterID, Information &info) override;

  bool validateNodesAndBasis(Domain *theDomain);
  bool isGeometryValid(void) const { return m_geometryValid; }

private:
  struct CachedSubLoad {
    int frameIdx = 0;
    double aOverL = 0.0;
    double bOverL = 1.0;
    double waUnit = 0.0;
    double wbUnit = 0.0;
    double length = 0.0;
  };

  struct CachedFrameInfo {
    int frameTag = 0;
    Element *framePtr = nullptr;
    double length = 0.0;
    double aL = 0.0;
    double bL = 1.0;
    double waEqUnit = 0.0;
    double wbEqUnit = 0.0;
    double areaUnit = 0.0;
    double overhangTorqueUnit = 0.0;
    std::array<double, 3> vx = {1.0, 0.0, 0.0};
    std::array<double, 3> vy = {0.0, 1.0, 0.0};
    std::array<double, 3> vz = {0.0, 0.0, 1.0};
    int extNodeI = -1;
    int extNodeJ = -1;
  };

  void rebuildExternalNodeIndex(void);
  bool buildPlaneBasis(void);
  int resolveAndBuildTransfer(void);
  bool hasMass(void) const;
  void computeMeanInPlaneStrain(Vector &eps2d) const;

  std::vector<int> m_polyNodeTags;
  std::vector<int> m_frameTags;
  std::vector<int> m_supportNodeTags;
  std::vector<std::vector<int>> m_holeLoops;

  std::string m_method;
  std::array<double, 3> m_spanAxis;
  bool m_hasSpanAxis;
  double m_orthoRatio;
  int m_segments;
  double m_meshSize;
  double m_thick;
  NDMaterial *m_theMaterial;
  double m_E;
  double m_nu;
  double m_rho;
  bool m_lumped;
  double m_pressure;
  double m_b1;
  double m_b2;
  double m_b3;

  ID m_connectedExternalNodes;
  std::vector<Node *> m_theNodes;
  std::vector<int> m_nodeDofOffsets;
  std::vector<int> m_supportExtIdx;
  std::vector<std::vector<int>> m_holeExtIdx;
  int m_numDOF;
  int m_ndm;

  bool m_geometryValid;
  bool m_framesResolved;
  bool m_hasInPlaneStiff;
  double m_D00_init;
  double m_rhoInit;
  double m_slabArea;
  double m_charLen;
  int m_supportRank;

  std::array<double, 3> m_origin;
  std::array<double, 3> m_e1;
  std::array<double, 3> m_e2;
  std::array<double, 3> m_normal;

  std::vector<CachedFrameInfo> m_frameInfos;
  std::vector<CachedSubLoad> m_unitSubLoads;
  std::vector<double> m_unitSupportAreas;
  std::vector<std::unique_ptr<ElementalLoad>> m_activeBeamLoads;

  Vector m_appliedTotalLoad;
  double m_appliedNormalPressure;
  std::array<double, 3> m_appliedGlobalTraction;
  std::vector<double> m_meshVertices3d;
  std::vector<int> m_meshTriangles;

  Matrix m_Kbase;
  Matrix m_Mbase;
  Matrix m_K;
  Matrix m_M;
  Vector m_P;
  Vector m_Q;
  Matrix *Ki;
  LadrunoMassCache m_massCache;
};

#endif // LadrunoMembrane_h
