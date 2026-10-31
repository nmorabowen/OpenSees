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

#include "LadrunoMembrane.h"
#include "LadrunoMembraneKernel.h"
#include "LadrunoMembraneDiaphragm.h"

#include <Beam2dPartialUniformLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam3dPartialUniformLoad.h>
#include <Beam3dUniformLoad.h>
#include <Channel.h>
#include <Domain.h>
#include <DummyStream.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <LadrunoElementShell.h>
#include <NDMaterial.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <Parameter.h>
#include <SelfWeight.h>
#include <SurfaceLoader.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

LadrunoMembrane::LadrunoMembrane(int tag,
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
                                 double b3,
                                 const std::vector<int> &diaphragmNodes,
                                 const std::array<double, 3> &inPlaneMods)
  : Element(tag, ELE_TAG_LadrunoMembrane),
    m_polyNodeTags(polyNodes),
    m_frameTags(frameTags),
    m_supportNodeTags(supportNodes),
    m_holeLoops(holeLoops),
    m_method(method),
    m_spanAxis(spanAxis),
    m_hasSpanAxis(hasSpanAxis),
    m_orthoRatio(orthoRatio),
    m_segments(segments),
    m_meshSize(meshSize),
    m_thick(thick),
    m_theMaterial(nullptr),
    m_E(E_direct),
    m_nu(nu_direct),
    m_rho(rho),
    m_lumped(lumped),
    m_pressure(pressure),
    m_b1(b1),
    m_b2(b2),
    m_b3(b3),
    m_connectedExternalNodes(1),
    m_numDOF(0),
    m_ndm(3),
    m_geometryValid(false),
    m_framesResolved(false),
    m_hasInPlaneStiff(false),
    m_D00_init(0.0),
    m_rhoInit(0.0),
    m_slabArea(0.0),
    m_charLen(0.0),
    m_supportRank(2),
    m_origin{0.0, 0.0, 0.0},
    m_e1{1.0, 0.0, 0.0},
    m_e2{0.0, 1.0, 0.0},
    m_normal{0.0, 0.0, 1.0},
    m_appliedTotalLoad(3),
    m_appliedNormalPressure(0.0),
    m_appliedGlobalTraction{0.0, 0.0, 0.0},
    m_diaphragmNodeTags(diaphragmNodes),
    m_inPlaneMods(inPlaneMods),
    m_Dmod{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    Ki(0)
{
  if (theMat != nullptr) {
    m_theMaterial = theMat->getCopy("PlaneStress");
    if (m_theMaterial == nullptr) {
      m_theMaterial = theMat->getCopy("PlaneStress2D");
    }
    if (m_theMaterial == nullptr) {
      m_theMaterial = theMat->getCopy();
    }
  }
  this->rebuildExternalNodeIndex();
}

LadrunoMembrane::LadrunoMembrane()
  : Element(0, ELE_TAG_LadrunoMembrane),
    m_method("tributary"),
    m_spanAxis{1.0, 0.0, 0.0},
    m_hasSpanAxis(false),
    m_orthoRatio(1.0),
    m_segments(0),
    m_meshSize(0.0),
    m_thick(0.0),
    m_theMaterial(nullptr),
    m_E(0.0),
    m_nu(0.2),
    m_rho(0.0),
    m_lumped(false),
    m_pressure(0.0),
    m_b1(0.0),
    m_b2(0.0),
    m_b3(0.0),
    m_connectedExternalNodes(1),
    m_numDOF(0),
    m_ndm(3),
    m_geometryValid(false),
    m_framesResolved(false),
    m_hasInPlaneStiff(false),
    m_D00_init(0.0),
    m_rhoInit(0.0),
    m_slabArea(0.0),
    m_charLen(0.0),
    m_supportRank(2),
    m_origin{0.0, 0.0, 0.0},
    m_e1{1.0, 0.0, 0.0},
    m_e2{0.0, 1.0, 0.0},
    m_normal{0.0, 0.0, 1.0},
    m_appliedTotalLoad(3),
    m_appliedNormalPressure(0.0),
    m_appliedGlobalTraction{0.0, 0.0, 0.0},
    m_inPlaneMods{1.0, 1.0, 1.0},
    m_Dmod{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    Ki(0)
{
}

LadrunoMembrane::~LadrunoMembrane()
{
  if (m_theMaterial != nullptr) {
    delete m_theMaterial;
    m_theMaterial = nullptr;
  }
  LadrunoShell::dropKi(Ki);
}

void
LadrunoMembrane::rebuildExternalNodeIndex(void)
{
  std::vector<int> allTags;
  auto addUnique = [&allTags](int tag) -> int {
    for (std::size_t i = 0; i < allTags.size(); ++i) {
      if (allTags[i] == tag) {
        return static_cast<int>(i);
      }
    }
    int idx = static_cast<int>(allTags.size());
    allTags.push_back(tag);
    return idx;
  };

  for (int t : m_polyNodeTags) {
    addUnique(t);
  }

  m_holeExtIdx.clear();
  m_holeExtIdx.resize(m_holeLoops.size());
  for (std::size_t h = 0; h < m_holeLoops.size(); ++h) {
    m_holeExtIdx[h].resize(m_holeLoops[h].size());
    for (std::size_t j = 0; j < m_holeLoops[h].size(); ++j) {
      m_holeExtIdx[h][j] = addUnique(m_holeLoops[h][j]);
    }
  }

  m_supportExtIdx.clear();
  m_supportExtIdx.resize(m_supportNodeTags.size());
  for (std::size_t s = 0; s < m_supportNodeTags.size(); ++s) {
    m_supportExtIdx[s] = addUnique(m_supportNodeTags[s]);
  }

  // WP-204: in-plane-only diaphragm nodes (beam nodes on an edge, wall nodes
  // on a hole edge, ...): retained by the condensation, no gravity share
  m_diaphragmExtIdx.clear();
  m_diaphragmExtIdx.resize(m_diaphragmNodeTags.size());
  for (std::size_t d = 0; d < m_diaphragmNodeTags.size(); ++d) {
    m_diaphragmExtIdx[d] = addUnique(m_diaphragmNodeTags[d]);
  }

  int nExt = static_cast<int>(allTags.size());
  if (nExt > 0) {
    m_connectedExternalNodes.resize(nExt);
    for (int i = 0; i < nExt; ++i) {
      m_connectedExternalNodes(i) = allTags[i];
    }
  } else {
    m_connectedExternalNodes.resize(1);
    m_connectedExternalNodes(0) = 0;
  }
  m_theNodes.assign(nExt, nullptr);
  m_nodeDofOffsets.assign(nExt, 0);
}

int
LadrunoMembrane::getNumExternalNodes(void) const
{
  return m_connectedExternalNodes.Size();
}

const ID &
LadrunoMembrane::getExternalNodes(void)
{
  return m_connectedExternalNodes;
}

Node **
LadrunoMembrane::getNodePtrs(void)
{
  return m_theNodes.empty() ? nullptr : m_theNodes.data();
}

int
LadrunoMembrane::getNumDOF(void)
{
  return m_numDOF;
}

bool
LadrunoMembrane::buildPlaneBasis(void)
{
  int nPoly = static_cast<int>(m_polyNodeTags.size());
  if (nPoly < 3 || static_cast<int>(m_theNodes.size()) < nPoly) {
    opserr << "LadrunoMembrane::buildPlaneBasis - element " << this->getTag()
           << " requires at least 3 polygon nodes\n";
    return false;
  }

  auto getCoord3 = [this](int extIdx) -> std::array<double, 3> {
    const Vector &c = m_theNodes[extIdx]->getCrds();
    double x = (c.Size() >= 1) ? c(0) : 0.0;
    double y = (c.Size() >= 2) ? c(1) : 0.0;
    double z = (c.Size() >= 3) ? c(2) : 0.0;
    return {x, y, z};
  };

  m_origin = getCoord3(0);
  double nx = 0.0, ny = 0.0, nz = 0.0;
  double maxEdge = 0.0;
  for (int i = 0; i < nPoly; ++i) {
    auto pi = getCoord3(i);
    auto pj = getCoord3((i + 1) % nPoly);
    double dx = pj[0] - pi[0];
    double dy = pj[1] - pi[1];
    double dz = pj[2] - pi[2];
    double edgeLen = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (edgeLen <= 1e-12) {
      opserr << "LadrunoMembrane::buildPlaneBasis - zero-length polygon edge on element "
             << this->getTag() << "\n";
      return false;
    }
    if (edgeLen > maxEdge) {
      maxEdge = edgeLen;
    }
    double ux = pi[0] - m_origin[0];
    double uy = pi[1] - m_origin[1];
    double uz = pi[2] - m_origin[2];
    double vx = pj[0] - m_origin[0];
    double vy = pj[1] - m_origin[1];
    double vz = pj[2] - m_origin[2];
    nx += uy * vz - uz * vy;
    ny += uz * vx - ux * vz;
    nz += ux * vy - uy * vx;
  }

  double normN = std::sqrt(nx * nx + ny * ny + nz * nz);
  double area = 0.5 * normN;
  if (maxEdge <= 1e-12 || area / (maxEdge * maxEdge) <= 1e-8) {
    opserr << "LadrunoMembrane::buildPlaneBasis - degenerate/collinear polygon on element "
           << this->getTag() << "\n";
    return false;
  }

  m_slabArea = area;
  m_charLen = std::sqrt(area);
  m_normal = {nx / normN, ny / normN, nz / normN};

  if (m_ndm == 2) {
    if (m_normal[2] >= 0.0) {
      m_e1 = {1.0, 0.0, 0.0};
      m_e2 = {0.0, 1.0, 0.0};
      m_normal = {0.0, 0.0, 1.0};
    } else {
      m_e1 = {1.0, 0.0, 0.0};
      m_e2 = {0.0, -1.0, 0.0};
      m_normal = {0.0, 0.0, -1.0};
    }
  } else {
    std::array<double, 3> ref = {1.0, 0.0, 0.0};
    if (std::abs(m_normal[0]) >= 0.99 && std::abs(m_normal[1]) < 1e-6) {
      ref = {0.0, 1.0, 0.0};
    } else if (std::abs(m_normal[2]) < 0.999999) {
      auto p1 = getCoord3(1);
      ref = {p1[0] - m_origin[0], p1[1] - m_origin[1], p1[2] - m_origin[2]};
    }
    double dotRN = ref[0] * m_normal[0] + ref[1] * m_normal[1] + ref[2] * m_normal[2];
    double e1x = ref[0] - dotRN * m_normal[0];
    double e1y = ref[1] - dotRN * m_normal[1];
    double e1z = ref[2] - dotRN * m_normal[2];
    double normE1 = std::sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
    if (normE1 <= 1e-12) {
      auto p1 = getCoord3(1);
      e1x = p1[0] - m_origin[0];
      e1y = p1[1] - m_origin[1];
      e1z = p1[2] - m_origin[2];
      dotRN = e1x * m_normal[0] + e1y * m_normal[1] + e1z * m_normal[2];
      e1x -= dotRN * m_normal[0];
      e1y -= dotRN * m_normal[1];
      e1z -= dotRN * m_normal[2];
      normE1 = std::sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
    }
    m_e1 = {e1x / normE1, e1y / normE1, e1z / normE1};
    m_e2 = {
      m_normal[1] * m_e1[2] - m_normal[2] * m_e1[1],
      m_normal[2] * m_e1[0] - m_normal[0] * m_e1[2],
      m_normal[0] * m_e1[1] - m_normal[1] * m_e1[0]
    };
  }

  // Check coplanarity of polygon, hole, and support vertices
  int nAll = static_cast<int>(m_theNodes.size());
  for (int i = 0; i < nAll; ++i) {
    auto pi = getCoord3(i);
    double dz = (pi[0] - m_origin[0]) * m_normal[0]
              + (pi[1] - m_origin[1]) * m_normal[1]
              + (pi[2] - m_origin[2]) * m_normal[2];
    if (std::abs(dz) > 1e-2 * maxEdge) {
      opserr << "LadrunoMembrane::buildPlaneBasis - non-planar node "
             << m_connectedExternalNodes(i) << " (|dz|/L = " << std::abs(dz) / maxEdge
             << " > 0.01) on element " << this->getTag() << "\n";
      return false;
    }
  }

  // Validate polygon and holes with PSLG preprocessor
  try {
    std::vector<LadrunoMembraneKernel::Point2D> outer2d(nPoly);
    for (int i = 0; i < nPoly; ++i) {
      auto pi = getCoord3(i);
      double dx = pi[0] - m_origin[0];
      double dy = pi[1] - m_origin[1];
      double dz = pi[2] - m_origin[2];
      outer2d[i] = {
        dx * m_e1[0] + dy * m_e1[1] + dz * m_e1[2],
        dx * m_e2[0] + dy * m_e2[1] + dz * m_e2[2]
      };
    }
    std::vector<std::vector<LadrunoMembraneKernel::Point2D>> holes2d(m_holeExtIdx.size());
    for (std::size_t h = 0; h < m_holeExtIdx.size(); ++h) {
      holes2d[h].resize(m_holeExtIdx[h].size());
      for (std::size_t j = 0; j < m_holeExtIdx[h].size(); ++j) {
        auto ph = getCoord3(m_holeExtIdx[h][j]);
        double dx = ph[0] - m_origin[0];
        double dy = ph[1] - m_origin[1];
        double dz = ph[2] - m_origin[2];
        holes2d[h][j] = {
          dx * m_e1[0] + dy * m_e1[1] + dz * m_e1[2],
          dx * m_e2[0] + dy * m_e2[1] + dz * m_e2[2]
        };
      }
    }
    auto testMesh = LadrunoMembraneKernel::build_pslg_cdt_mesh(outer2d, holes2d, {}, {}, 0.0);
    m_slabArea = testMesh.total_area;
    m_charLen = testMesh.char_length;
  } catch (const std::exception &ex) {
    opserr << "LadrunoMembrane::buildPlaneBasis - invalid PSLG geometry on element "
           << this->getTag() << ": " << ex.what() << "\n";
    return false;
  }

  return true;
}

bool
LadrunoMembrane::validateNodesAndBasis(Domain *theDomain)
{
  if (theDomain == nullptr) {
    return false;
  }
  this->rebuildExternalNodeIndex();
  int nExt = m_connectedExternalNodes.Size();
  m_numDOF = 0;
  m_ndm = 0;

  for (int i = 0; i < nExt; ++i) {
    int nodeTag = m_connectedExternalNodes(i);
    Node *nd = theDomain->getNode(nodeTag);
    if (nd == nullptr) {
      opserr << "LadrunoMembrane::setDomain - node " << nodeTag
             << " does not exist in Domain for element " << this->getTag() << "\n";
      return false;
    }
    m_theNodes[i] = nd;
    int crdSize = nd->getCrds().Size();
    if (i == 0) {
      m_ndm = crdSize;
      if (m_ndm != 2 && m_ndm != 3) {
        opserr << "LadrunoMembrane::setDomain - invalid ndm = " << m_ndm
               << " on node " << nodeTag << " for element " << this->getTag() << "\n";
        return false;
      }
    } else if (crdSize != m_ndm) {
      opserr << "LadrunoMembrane::setDomain - inconsistent nodal dimension on node "
             << nodeTag << " for element " << this->getTag() << "\n";
      return false;
    }
    int ndf = nd->getNumberDOF();
    if (ndf < m_ndm) {
      opserr << "LadrunoMembrane::setDomain - node " << nodeTag << " has ndf = " << ndf
             << " < ndm = " << m_ndm << " for element " << this->getTag() << "\n";
      return false;
    }
    m_nodeDofOffsets[i] = m_numDOF;
    m_numDOF += ndf;
  }

  if (!this->buildPlaneBasis()) {
    return false;
  }
  return true;
}

void
LadrunoMembrane::setDomain(Domain *theDomain)
{
  m_geometryValid = false;
  m_framesResolved = false;
  m_activeBeamLoads.clear();
  LadrunoShell::dropKi(Ki);
  m_massCache.invalidate();

  if (theDomain == nullptr) {
    for (std::size_t i = 0; i < m_theNodes.size(); ++i) {
      m_theNodes[i] = nullptr;
    }
    m_numDOF = 0;
    this->DomainComponent::setDomain(theDomain);
    return;
  }

  if (!this->validateNodesAndBasis(theDomain)) {
    return;
  }

  m_Kbase.resize(m_numDOF, m_numDOF);
  m_Mbase.resize(m_numDOF, m_numDOF);
  m_K.resize(m_numDOF, m_numDOF);
  m_M.resize(m_numDOF, m_numDOF);
  m_P.resize(m_numDOF);
  m_Q.resize(m_numDOF);
  m_Kbase.Zero();
  m_Mbase.Zero();
  m_K.Zero();
  m_M.Zero();
  m_P.Zero();
  m_Q.Zero();

  m_geometryValid = true;
  this->DomainComponent::setDomain(theDomain);
}

int
LadrunoMembrane::resolveAndBuildTransfer(void)
{
  if (m_framesResolved) {
    return 0;
  }
  if (!m_geometryValid) {
    return -1;
  }
  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr) {
    return -1;
  }

  auto project3to2 = [this](double x, double y, double z) -> LadrunoMembraneKernel::Point2D {
    double dx = x - m_origin[0];
    double dy = y - m_origin[1];
    double dz = z - m_origin[2];
    return {
      dx * m_e1[0] + dy * m_e1[1] + dz * m_e1[2],
      dx * m_e2[0] + dy * m_e2[1] + dz * m_e2[2]
    };
  };

  auto getNode3 = [](Node *nd) -> std::array<double, 3> {
    const Vector &c = nd->getCrds();
    return {
      (c.Size() >= 1) ? c(0) : 0.0,
      (c.Size() >= 2) ? c(1) : 0.0,
      (c.Size() >= 3) ? c(2) : 0.0
    };
  };

  int nPoly = static_cast<int>(m_polyNodeTags.size());
  std::vector<LadrunoMembraneKernel::Point2D> outer2d(nPoly);
  for (int i = 0; i < nPoly; ++i) {
    auto p3 = getNode3(m_theNodes[i]);
    outer2d[i] = project3to2(p3[0], p3[1], p3[2]);
  }

  std::vector<std::vector<LadrunoMembraneKernel::Point2D>> holes2d(m_holeExtIdx.size());
  for (std::size_t h = 0; h < m_holeExtIdx.size(); ++h) {
    holes2d[h].resize(m_holeExtIdx[h].size());
    for (std::size_t j = 0; j < m_holeExtIdx[h].size(); ++j) {
      auto p3 = getNode3(m_theNodes[m_holeExtIdx[h][j]]);
      holes2d[h][j] = project3to2(p3[0], p3[1], p3[2]);
    }
  }

  std::vector<LadrunoMembraneKernel::Point2D> support2d(m_supportExtIdx.size());
  for (std::size_t s = 0; s < m_supportExtIdx.size(); ++s) {
    auto p3 = getNode3(m_theNodes[m_supportExtIdx[s]]);
    support2d[s] = project3to2(p3[0], p3[1], p3[2]);
  }

  int nFrames = static_cast<int>(m_frameTags.size());
  std::vector<LadrunoMembraneKernel::Segment2D> frames2d(nFrames);
  m_frameInfos.clear();
  m_frameInfos.resize(nFrames);

  for (int k = 0; k < nFrames; ++k) {
    int fTag = m_frameTags[k];
    Element *fEle = theDomain->getElement(fTag);
    if (fEle == nullptr) {
      opserr << "LadrunoMembrane::resolveAndBuildTransfer - supporting frame element "
             << fTag << " not found in Domain for membrane " << this->getTag() << "\n";
      return -1;
    }
    if (fEle->getNumExternalNodes() < 2) {
      opserr << "LadrunoMembrane::resolveAndBuildTransfer - frame element " << fTag
             << " has fewer than 2 external nodes\n";
      return -1;
    }
    Node **fNodes = fEle->getNodePtrs();
    const ID &fExtIds = fEle->getExternalNodes();
    Node *ndI = (fNodes != nullptr && fNodes[0] != nullptr)
                  ? fNodes[0]
                  : theDomain->getNode(fExtIds(0));
    Node *ndJ = (fNodes != nullptr && fNodes[1] != nullptr)
                  ? fNodes[1]
                  : theDomain->getNode(fExtIds(1));
    if (ndI == nullptr || ndJ == nullptr) {
      opserr << "LadrunoMembrane::resolveAndBuildTransfer - frame element " << fTag
             << " has null endpoint nodes\n";
      return -1;
    }

    auto pI = getNode3(ndI);
    auto pJ = getNode3(ndJ);
    double dx = pJ[0] - pI[0];
    double dy = pJ[1] - pI[1];
    double dz = pJ[2] - pI[2];
    double L_k = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (L_k <= 1e-12) {
      opserr << "LadrunoMembrane::resolveAndBuildTransfer - frame element " << fTag
             << " has zero length\n";
      return -1;
    }

    frames2d[k] = {
      project3to2(pI[0], pI[1], pI[2]),
      project3to2(pJ[0], pJ[1], pJ[2])
    };

    CachedFrameInfo info{};
    info.frameTag = fTag;
    info.framePtr = fEle;
    info.length = L_k;
    info.vx = {dx / L_k, dy / L_k, dz / L_k};

    // Match frame end nodes against connectedExternalNodes if present
    int tagI = ndI->getTag();
    int tagJ = ndJ->getTag();
    for (int idx = 0; idx < m_connectedExternalNodes.Size(); ++idx) {
      if (m_connectedExternalNodes(idx) == tagI) info.extNodeI = idx;
      if (m_connectedExternalNodes(idx) == tagJ) info.extNodeJ = idx;
    }

    if (m_ndm == 2) {
      info.vy = {-info.vx[1], info.vx[0], 0.0};
      info.vz = {0.0, 0.0, 1.0};
    } else {
      bool gotAxes = false;
      DummyStream dummy;
      const char *req[1] = {"localAxes"};
      Response *resp = fEle->setResponse(req, 1, dummy);
      if (resp != nullptr) {
        if (resp->getResponse() == 0) {
          const Vector &la = resp->getInformation().getData();
          if (la.Size() == 9) {
            double vyNorm = std::sqrt(la(3) * la(3) + la(4) * la(4) + la(5) * la(5));
            if (vyNorm > 0.5) {
              info.vx = {la(0), la(1), la(2)};
              info.vy = {la(3), la(4), la(5)};
              info.vz = {la(6), la(7), la(8)};
              gotAxes = true;
            }
          }
        }
        delete resp;
      }
      if (!gotAxes) {
        std::array<double, 3> vxz = {0.0, 0.0, 1.0};
        if (std::abs(info.vx[2]) > 0.99) {
          vxz = {-1.0, 0.0, 0.0};
        }
        double vyx = vxz[1] * info.vx[2] - vxz[2] * info.vx[1];
        double vyy = vxz[2] * info.vx[0] - vxz[0] * info.vx[2];
        double vyz = vxz[0] * info.vx[1] - vxz[1] * info.vx[0];
        double vyLen = std::sqrt(vyx * vyx + vyy * vyy + vyz * vyz);
        if (vyLen > 1e-12) {
          info.vy = {vyx / vyLen, vyy / vyLen, vyz / vyLen};
          info.vz = {
            info.vx[1] * info.vy[2] - info.vx[2] * info.vy[1],
            info.vx[2] * info.vy[0] - info.vx[0] * info.vy[2],
            info.vx[0] * info.vy[1] - info.vx[1] * info.vy[0]
          };
        }
      }
    }
    m_frameInfos[k] = info;
  }

  LadrunoMembraneKernel::Point2D span2d = {1.0, 0.0};
  if (m_hasSpanAxis) {
    double sx = m_spanAxis[0] * m_e1[0] + m_spanAxis[1] * m_e1[1] + m_spanAxis[2] * m_e1[2];
    double sy = m_spanAxis[0] * m_e2[0] + m_spanAxis[1] * m_e2[1] + m_spanAxis[2] * m_e2[2];
    double sNorm = std::hypot(sx, sy);
    if (sNorm <= 1e-12) {
      opserr << "LadrunoMembrane::resolveAndBuildTransfer - spanAxis is orthogonal to membrane plane on element "
             << this->getTag() << "\n";
      return -1;
    }
    span2d = {sx / sNorm, sy / sNorm};
  }

  double E_mat = m_E;
  double nu_mat = m_nu;
  m_D00_init = 0.0;
  if (m_theMaterial != nullptr) {
    const Matrix &D0 = m_theMaterial->getInitialTangent();
    if (D0.noRows() >= 2 && D0.noCols() >= 2 && std::abs(D0(0, 0)) > 1e-15) {
      m_D00_init = D0(0, 0);
      nu_mat = std::max(-0.99, std::min(0.499999, D0(0, 1) / D0(0, 0)));
      E_mat = D0(0, 0) * (1.0 - nu_mat * nu_mat);
    }
  }

  try {
    LadrunoMembraneKernel::DistributionResult distRes;
    bool hasExplicitSupports = (!frames2d.empty() || !support2d.empty());
    if (hasExplicitSupports) {
      distRes = LadrunoMembraneKernel::distribute_surface_load(
          outer2d, 1.0, m_method, frames2d, holes2d, support2d,
          span2d, m_orthoRatio, nu_mat, m_meshSize);
    } else {
      distRes = LadrunoMembraneKernel::distribute_surface_load(
          outer2d, 1.0, "tributary", {}, holes2d, outer2d,
          span2d, m_orthoRatio, nu_mat, m_meshSize);
    }

    m_slabArea = distRes.mesh.total_area;
    m_charLen = distRes.mesh.char_length;
    m_supportRank = distRes.plate_support_rank;

    m_meshVertices3d.resize(3 * distRes.mesh.vertices.size());
    for (std::size_t v = 0; v < distRes.mesh.vertices.size(); ++v) {
      double vx = distRes.mesh.vertices[v].x;
      double vy = distRes.mesh.vertices[v].y;
      m_meshVertices3d[3 * v + 0] = m_origin[0] + vx * m_e1[0] + vy * m_e2[0];
      m_meshVertices3d[3 * v + 1] = m_origin[1] + vx * m_e1[1] + vy * m_e2[1];
      m_meshVertices3d[3 * v + 2] = m_origin[2] + vx * m_e1[2] + vy * m_e2[2];
    }
    m_meshTriangles.resize(3 * distRes.mesh.triangles.size());
    for (std::size_t t = 0; t < distRes.mesh.triangles.size(); ++t) {
      m_meshTriangles[3 * t + 0] = distRes.mesh.triangles[t][0];
      m_meshTriangles[3 * t + 1] = distRes.mesh.triangles[t][1];
      m_meshTriangles[3 * t + 2] = distRes.mesh.triangles[t][2];
    }

    m_meshVertexFlags = distRes.vertex_flags;
    m_meshNodalFieldUnit = distRes.nodal_field;

    auto mapFeatureToTag = [&](int fid) -> double {
      if (hasExplicitSupports) {
        int nF = static_cast<int>(m_frameTags.size());
        int nS = static_cast<int>(m_supportNodeTags.size());
        if (fid >= 0 && fid < nF) {
          return static_cast<double>(m_frameTags[fid]);
        }
        if (fid >= nF && (fid - nF) < nS) {
          return -static_cast<double>(m_supportNodeTags[fid - nF]);
        }
      } else {
        int nP = static_cast<int>(m_polyNodeTags.size());
        if (fid >= 0 && fid < nP) {
          return -static_cast<double>(m_polyNodeTags[fid]);
        }
      }
      return 0.0;
    };

    m_meshTriangleLoadFlowUnit.resize(4 * distRes.mesh.triangles.size(), 0.0);
    for (std::size_t t = 0; t < distRes.mesh.triangles.size(); ++t) {
      int fid = (t < distRes.triangle_support_id.size()) ? distRes.triangle_support_id[t] : -1;
      double fx2 = (t < distRes.triangle_flux_2d.size()) ? distRes.triangle_flux_2d[t].x : 0.0;
      double fy2 = (t < distRes.triangle_flux_2d.size()) ? distRes.triangle_flux_2d[t].y : 0.0;
      m_meshTriangleLoadFlowUnit[4 * t + 0] = mapFeatureToTag(fid);
      m_meshTriangleLoadFlowUnit[4 * t + 1] = fx2 * m_e1[0] + fy2 * m_e2[0];
      m_meshTriangleLoadFlowUnit[4 * t + 2] = fx2 * m_e1[1] + fy2 * m_e2[1];
      m_meshTriangleLoadFlowUnit[4 * t + 3] = fx2 * m_e1[2] + fy2 * m_e2[2];
    }

    m_yieldLines3d.resize(8 * distRes.yield_lines.size(), 0.0);
    for (std::size_t r = 0; r < distRes.yield_lines.size(); ++r) {
      const auto &yl = distRes.yield_lines[r];
      m_yieldLines3d[8 * r + 0] = m_origin[0] + yl.x1 * m_e1[0] + yl.y1 * m_e2[0];
      m_yieldLines3d[8 * r + 1] = m_origin[1] + yl.x1 * m_e1[1] + yl.y1 * m_e2[1];
      m_yieldLines3d[8 * r + 2] = m_origin[2] + yl.x1 * m_e1[2] + yl.y1 * m_e2[2];
      m_yieldLines3d[8 * r + 3] = m_origin[0] + yl.x2 * m_e1[0] + yl.y2 * m_e2[0];
      m_yieldLines3d[8 * r + 4] = m_origin[1] + yl.x2 * m_e1[1] + yl.y2 * m_e2[1];
      m_yieldLines3d[8 * r + 5] = m_origin[2] + yl.x2 * m_e1[2] + yl.y2 * m_e2[2];
      m_yieldLines3d[8 * r + 6] = mapFeatureToTag(yl.fidA);
      m_yieldLines3d[8 * r + 7] = mapFeatureToTag(yl.fidB);
    }

    m_meshSubSegmentsUnit.clear();
    m_unitSubLoads.clear();
    if (hasExplicitSupports) {
      for (std::size_t k = 0; k < distRes.frame_loads.size() && k < m_frameInfos.size(); ++k) {
        const auto &fl = distRes.frame_loads[k];
        m_frameInfos[k].aL = fl.aL;
        m_frameInfos[k].bL = fl.bL;
        m_frameInfos[k].waEqUnit = fl.wa_eq;
        m_frameInfos[k].wbEqUnit = fl.wb_eq;
        m_frameInfos[k].areaUnit = fl.F_total;
        m_frameInfos[k].overhangTorqueUnit = fl.overhang_torque;

        for (const auto &sub : fl.sub_loads) {
          m_meshSubSegmentsUnit.push_back(static_cast<double>(sub.u));
          m_meshSubSegmentsUnit.push_back(static_cast<double>(sub.v));
          m_meshSubSegmentsUnit.push_back(static_cast<double>(m_frameInfos[k].frameTag));
          m_meshSubSegmentsUnit.push_back(sub.aL);
          m_meshSubSegmentsUnit.push_back(sub.bL);
          m_meshSubSegmentsUnit.push_back(0.5 * (sub.wa + sub.wb));
        }

        double L_k = m_frameInfos[k].length;
        if (m_segments > 0 && !fl.sub_loads.empty() && fl.bL > fl.aL + 1e-12) {
          double spanA = fl.aL;
          double spanB = fl.bL;
          double dSpan = (spanB - spanA) / static_cast<double>(m_segments);
          for (int segIdx = 0; segIdx < m_segments; ++segIdx) {
            double t0 = spanA + segIdx * dSpan;
            double t1 = (segIdx + 1 == m_segments) ? spanB : (spanA + (segIdx + 1) * dSpan);
            double L_seg = (t1 - t0) * L_k;
            if (L_seg <= 1e-14) continue;
            double F_bin = 0.0;
            double M_bin = 0.0;
            for (const auto &sub : fl.sub_loads) {
              double sa = std::max(t0, sub.aL);
              double sb = std::min(t1, sub.bL);
              if (sb <= sa + 1e-14 || sub.bL <= sub.aL + 1e-14) continue;
              double r0 = (sa - sub.aL) / (sub.bL - sub.aL);
              double r1 = (sb - sub.aL) / (sub.bL - sub.aL);
              double w0 = sub.wa + r0 * (sub.wb - sub.wa);
              double w1 = sub.wa + r1 * (sub.wb - sub.wa);
              double dL = (sb - sa) * L_k;
              double dF = 0.5 * (w0 + w1) * dL;
              double dM_local = dL * dL * (w0 + 2.0 * w1) / 6.0;
              F_bin += dF;
              M_bin += dF * ((sa - t0) * L_k) + dM_local;
            }
            auto [wa_b, wb_b] = LadrunoMembraneKernel::moment_conserving_trapezoid(F_bin, M_bin, L_seg);
            if (std::abs(wa_b) + std::abs(wb_b) > 1e-15) {
              m_unitSubLoads.push_back({static_cast<int>(k), t0, t1, wa_b, wb_b, L_seg});
            }
          }
        } else {
          for (const auto &sub : fl.sub_loads) {
            double L_seg = (sub.bL - sub.aL) * L_k;
            if (L_seg > 1e-14 && (std::abs(sub.wa) + std::abs(sub.wb) > 1e-15)) {
              m_unitSubLoads.push_back({static_cast<int>(k), sub.aL, sub.bL, sub.wa, sub.wb, L_seg});
            }
          }
        }
      }
      m_unitSupportAreas = distRes.support_node_reactions;
    } else {
      m_unitSupportAreas = distRes.support_node_reactions;
    }

    // In-plane (diaphragm) stiffness and mass if thickness > 0. WP-204: exact
    // Guyan condensation of an LST mesh onto EVERY external node (the ADR-98
    // CST condensation kept the outer polygon only and slaved the rest)
    m_Kbase.Zero();
    m_Mbase.Zero();
    m_hasInPlaneStiff = false;
    double rhoEff = (m_rho > 0.0) ? m_rho : ((m_theMaterial != nullptr) ? m_theMaterial->getRho() : 0.0);
    m_rhoInit = rhoEff;

    bool wantStiff = (m_thick > 0.0 && E_mat > 0.0);
    bool wantMass = (m_thick > 0.0 && rhoEff > 0.0);
    m_dStrainOp.clear();
    m_dMeanStrainOp.clear();
    m_dVerts3d.clear();
    m_dTris.clear();
    m_dTriArea.clear();
    if (wantStiff || wantMass) {
      // WP-204: plane-stress modulus (the material's full initial tangent, or
      // isotropic E/nu) times the ETABS-style f11/f22/f12 modifiers, applied as
      // D_ij * sqrt(f_i f_j) exactly like LadrunoShellModifier's membrane block.
      double Dm[9];
      bool fromMat = false;
      if (wantStiff && m_theMaterial != nullptr) {
        const Matrix &D0 = m_theMaterial->getInitialTangent();
        if (D0.noRows() >= 3 && D0.noCols() >= 3 && std::abs(D0(0, 0)) > 1e-15) {
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
              Dm[3 * i + j] = D0(i, j);
          fromMat = true;
        }
      }
      if (!fromMat) {
        double E_c = wantStiff ? E_mat : 1.0;
        double nu_c = wantStiff ? nu_mat : 0.2;
        double fac = E_c / (1.0 - nu_c * nu_c);
        double iso[9] = {fac, fac * nu_c, 0.0, fac * nu_c, fac, 0.0, 0.0, 0.0, 0.5 * fac * (1.0 - nu_c)};
        std::copy(iso, iso + 9, Dm);
      }
      if (wantStiff) {
        double sf[3] = {std::sqrt(m_inPlaneMods[0]), std::sqrt(m_inPlaneMods[1]),
                        std::sqrt(m_inPlaneMods[2])};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            Dm[3 * i + j] *= sf[i] * sf[j];
      }
      std::copy(Dm, Dm + 9, m_Dmod);

      // The diaphragm mesh: outer loop + holes + every other external node as
      // a mesh vertex (supports and -diaphragmNodes). Frames are not segments
      // here; their end nodes are retained only when they are external nodes.
      int nExt = m_connectedExternalNodes.Size();
      std::vector<LadrunoMembraneKernel::Point2D> extra2d;
      std::vector<LadrunoMembraneKernel::Point2D> ext2d(nExt);
      for (int a = 0; a < nExt; ++a) {
        auto p3 = getNode3(m_theNodes[a]);
        ext2d[a] = project3to2(p3[0], p3[1], p3[2]);
      }
      for (int idx : m_supportExtIdx) extra2d.push_back(ext2d[idx]);
      for (int idx : m_diaphragmExtIdx) extra2d.push_back(ext2d[idx]);
      // (the gravity mesh keeps the kernel default of 1200 vertices; the
      // diaphragm mesh allows 6000 so -meshSize is honoured on large slabs —
      // the sparse LDLt condensation is cheap, the O(n^2) mesher is the cost)
      // WP-208: connection lines, along which the slab edge is straight
      // between consecutive element nodes: every outer / hole loop segment and
      // every supporting frame whose two ends are element nodes (the frame is
      // also a mesher segment, so mesh nodes fall on it). Each line is split
      // at every element node lying on it.
      std::vector<std::pair<int, int>> baseLines;
      for (int i = 0; i < nPoly; ++i)
        baseLines.push_back({i, (i + 1) % nPoly});
      for (const auto &loop : m_holeExtIdx)
        for (std::size_t j = 0; j < loop.size(); ++j)
          baseLines.push_back({loop[j], loop[(j + 1) % loop.size()]});
      std::vector<LadrunoMembraneKernel::Segment2D> linkFrames2d;
      for (const auto &fi : m_frameInfos) {
        if (fi.extNodeI >= 0 && fi.extNodeJ >= 0 && fi.extNodeI != fi.extNodeJ) {
          baseLines.push_back({fi.extNodeI, fi.extNodeJ});
          linkFrames2d.push_back({ext2d[fi.extNodeI], ext2d[fi.extNodeJ]});
        }
      }
      std::vector<LadrunoMembraneDiaphragm::LinkLine> links;
      for (const auto &bl : baseLines) {
        const auto &A = ext2d[bl.first];
        const auto &B = ext2d[bl.second];
        double dx = B.x - A.x, dy = B.y - A.y, len2 = dx * dx + dy * dy;
        if (len2 <= 0.0) continue;
        std::vector<std::pair<double, int>> on = {{0.0, bl.first}, {1.0, bl.second}};
        double tolL = 1e-7 * std::sqrt(len2);
        for (int c = 0; c < nExt; ++c) {
          if (c == bl.first || c == bl.second) continue;
          double s = ((ext2d[c].x - A.x) * dx + (ext2d[c].y - A.y) * dy) / len2;
          if (s <= 1e-9 || s >= 1.0 - 1e-9) continue;
          if (std::hypot(ext2d[c].x - (A.x + s * dx), ext2d[c].y - (A.y + s * dy)) <= tolL)
            on.push_back({s, c});
        }
        std::sort(on.begin(), on.end());
        for (std::size_t k = 0; k + 1 < on.size(); ++k) {
          LadrunoMembraneDiaphragm::LinkLine L;
          L.a = ext2d[on[k].second];
          L.b = ext2d[on[k + 1].second];
          L.slotA = on[k].second;
          L.slotB = on[k + 1].second;
          links.push_back(L);
        }
      }

      auto dmesh = LadrunoMembraneKernel::build_pslg_cdt_mesh(outer2d, holes2d, linkFrames2d, extra2d, m_meshSize, 6000);

      // retained slot a (= external node a) -> mesh vertex
      double tol = 1e-6 * std::max(dmesh.char_length, 1e-12);
      std::vector<int> retained(nExt, -1);
      std::vector<char> taken(dmesh.vertices.size(), 0);
      for (int a = 0; a < nExt; ++a) {
        int best = -1;
        double bestD = tol;
        for (std::size_t v = 0; v < dmesh.vertices.size(); ++v) {
          double d = std::hypot(dmesh.vertices[v].x - ext2d[a].x, dmesh.vertices[v].y - ext2d[a].y);
          if (d <= bestD) { bestD = d; best = static_cast<int>(v); }
        }
        if (best >= 0 && taken[best]) {
          opserr << "WARNING LadrunoMembrane " << this->getTag() << ": node "
                 << m_connectedExternalNodes(a) << " coincides with another node of the"
                 << " panel; it receives no diaphragm stiffness\n";
          best = -1;
        }
        if (best >= 0) taken[best] = 1;
        retained[a] = best;
      }
      auto cond = LadrunoMembraneDiaphragm::condense_lst(
          dmesh.vertices, dmesh.triangles, retained, Dm, m_thick,
          wantMass ? rhoEff * m_thick : 0.0, links);
      for (int slot : cond.isolated) {
        opserr << "WARNING LadrunoMembrane " << this->getTag() << ": node "
               << m_connectedExternalNodes(slot) << " is not on the panel (outside it or"
               << " inside a hole); it receives no diaphragm stiffness\n";
      }
      // a frame node lying on the panel that is not an external node gets no
      // in-plane stiffness: say so (the bridge adds it with -diaphragmNodes)
      for (const auto &fi : m_frameInfos) {
        if (fi.framePtr == nullptr) continue;
        Node **fN = fi.framePtr->getNodePtrs();
        for (int e = 0; e < 2; ++e) {
          int ext = (e == 0) ? fi.extNodeI : fi.extNodeJ;
          if (ext >= 0 || fN == nullptr || fN[e] == nullptr) continue;
          auto p3 = getNode3(fN[e]);
          auto p2 = project3to2(p3[0], p3[1], p3[2]);
          double dz = (p3[0] - m_origin[0]) * m_normal[0] + (p3[1] - m_origin[1]) * m_normal[1] +
                      (p3[2] - m_origin[2]) * m_normal[2];
          if (std::abs(dz) <= 1e-6 * m_charLen &&
              LadrunoMembraneKernel::point_in_polygon(p2, outer2d)) {
            opserr << "WARNING LadrunoMembrane " << this->getTag() << ": node "
                   << fN[e]->getTag() << " of supporting frame " << fi.frameTag
                   << " lies on the panel but is not one of its nodes; it receives no"
                   << " diaphragm stiffness (add it with -diaphragmNodes)\n";
          }
        }
      }

      // keep the recovery operators and the mesh for the force outputs
      m_dStrainOp = cond.strainOp;
      m_dMeanStrainOp = cond.meanStrainOp;
      m_dTriArea = cond.triArea;
      m_dVerts3d.resize(3 * dmesh.vertices.size());
      for (std::size_t v = 0; v < dmesh.vertices.size(); ++v) {
        for (int i = 0; i < 3; ++i)
          m_dVerts3d[3 * v + i] = m_origin[i] + dmesh.vertices[v].x * m_e1[i] + dmesh.vertices[v].y * m_e2[i];
      }
      m_dTris.resize(3 * dmesh.triangles.size());
      for (std::size_t t = 0; t < dmesh.triangles.size(); ++t)
        for (int i = 0; i < 3; ++i)
          m_dTris[3 * t + i] = dmesh.triangles[t][i];

      int nb_dof = 2 * nExt;
      if (wantStiff) {
        m_hasInPlaneStiff = true;
        for (int a = 0; a < nExt; ++a) {
          int off_a = m_nodeDofOffsets[a];
          for (int b = 0; b < nExt; ++b) {
            int off_b = m_nodeDofOffsets[b];
            double k00 = cond.K[(2 * a + 0) * nb_dof + (2 * b + 0)];
            double k01 = cond.K[(2 * a + 0) * nb_dof + (2 * b + 1)];
            double k10 = cond.K[(2 * a + 1) * nb_dof + (2 * b + 0)];
            double k11 = cond.K[(2 * a + 1) * nb_dof + (2 * b + 1)];
            for (int i = 0; i < m_ndm; ++i) {
              double e1i = m_e1[i];
              double e2i = m_e2[i];
              for (int j = 0; j < m_ndm; ++j) {
                double e1j = m_e1[j];
                double e2j = m_e2[j];
                m_Kbase(off_a + i, off_b + j) +=
                    e1i * (k00 * e1j + k01 * e2j) +
                    e2i * (k10 * e1j + k11 * e2j);
              }
            }
          }
        }
      }

      if (wantMass) {
        // the Guyan mass of an LST mesh: its row sums can go negative at a
        // node, so -lumped uses HRZ on the condensed mass (diagonal scaled to
        // the panel mass, always positive)
        double diagSum = 0.0, total = 0.0;
        for (int a = 0; a < nExt; ++a) {
          diagSum += cond.M[(2 * a) * nb_dof + (2 * a)];
          for (int b = 0; b < nExt; ++b)
            total += cond.M[(2 * a) * nb_dof + (2 * b)];
        }
        double hrz = (diagSum > 0.0) ? total / diagSum : 0.0;
        for (int a = 0; a < nExt; ++a) {
          int off_a = m_nodeDofOffsets[a];
          if (m_lumped) {
            double maa = hrz * cond.M[(2 * a) * nb_dof + (2 * a)];
            for (int d = 0; d < m_ndm; ++d)
              m_Mbase(off_a + d, off_a + d) += maa;
            continue;
          }
          for (int b = 0; b < nExt; ++b) {
            int off_b = m_nodeDofOffsets[b];
            double mab = cond.M[(2 * a) * nb_dof + (2 * b)];
            for (int d = 0; d < m_ndm; ++d)
              m_Mbase(off_a + d, off_b + d) += mab;
          }
        }
      }
    }
  } catch (const std::exception &ex) {
    opserr << "LadrunoMembrane::resolveAndBuildTransfer (element " << this->getTag()
           << ") - " << ex.what() << "\n";
    return -1;
  }

  m_framesResolved = true;
  return 0;
}

bool
LadrunoMembrane::hasMass(void) const
{
  double r = (m_rho > 0.0) ? m_rho : ((m_theMaterial != nullptr) ? m_theMaterial->getRho() : 0.0);
  return (r > 0.0 && m_thick > 0.0);
}

// WP-204: (u.e1, u.e2) of every external node, the retained in-plane DOFs
void
LadrunoMembrane::gatherInPlaneDisp(std::vector<double> &ur) const
{
  int nExt = static_cast<int>(m_theNodes.size());
  ur.assign(2 * nExt, 0.0);
  for (int a = 0; a < nExt; ++a) {
    if (m_theNodes[a] == nullptr) continue;
    const Vector &d = m_theNodes[a]->getTrialDisp();
    double ux = (d.Size() >= 1) ? d(0) : 0.0;
    double uy = (d.Size() >= 2) ? d(1) : 0.0;
    double uz = (d.Size() >= 3 && m_ndm == 3) ? d(2) : 0.0;
    ur[2 * a] = ux * m_e1[0] + uy * m_e1[1] + uz * m_e1[2];
    ur[2 * a + 1] = ux * m_e2[0] + uy * m_e2[1] + uz * m_e2[2];
  }
}

// WP-204: the panel-mean in-plane strain of the condensed field (exact area
// average of the LST strain; the stock boundary integral used the outer loop
// only, so holes and interior nodes were ignored)
void
LadrunoMembrane::computeMeanInPlaneStrain(Vector &eps2d) const
{
  eps2d.Zero();
  int n2 = 2 * static_cast<int>(m_theNodes.size());
  if (m_dMeanStrainOp.size() != static_cast<std::size_t>(3 * n2) || n2 == 0) {
    return;
  }
  std::vector<double> ur;
  this->gatherInPlaneDisp(ur);
  for (int r = 0; r < 3; ++r) {
    double s = 0.0;
    for (int j = 0; j < n2; ++j)
      s += m_dMeanStrainOp[r * n2 + j] * ur[j];
    eps2d(r) = s;
  }
}

// ratio of the material's current to initial D00 (1 without a material): the
// nonlinear scaling of the condensed in-plane stiffness (unchanged from ADR-98)
double
LadrunoMembrane::tangentScale(void)
{
  double scale = 1.0;
  if (m_theMaterial != nullptr && std::abs(m_D00_init) > 1e-15) {
    const Matrix &Dt = m_theMaterial->getTangent();
    if (Dt.noRows() >= 2 && Dt.noCols() >= 2) {
      scale = Dt(0, 0) / m_D00_init;
    }
  }
#ifdef LADRUNO_MUTATE_MEMBRANE
  scale = 0.0;
#endif
  return scale;
}

int
LadrunoMembrane::commitState(void)
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "LadrunoMembrane::commitState() - failed in base class\n";
    return retVal;
  }
  if (m_theMaterial != nullptr) {
    retVal += m_theMaterial->commitState();
  }
  return retVal;
}

int
LadrunoMembrane::revertToLastCommit(void)
{
  if (m_theMaterial != nullptr) {
    return m_theMaterial->revertToLastCommit();
  }
  return 0;
}

int
LadrunoMembrane::revertToStart(void)
{
  if (m_theMaterial != nullptr) {
    return m_theMaterial->revertToStart();
  }
  return 0;
}

int
LadrunoMembrane::update(void)
{
  if (!m_geometryValid) {
    return -1;
  }
  if (!m_framesResolved) {
    Domain *theDomain = this->getDomain();
    if (theDomain != nullptr) {
      bool allFramesPresent = true;
      for (int fTag : m_frameTags) {
        if (theDomain->getElement(fTag) == nullptr) {
          allFramesPresent = false;
          break;
        }
      }
      if (allFramesPresent) {
        if (this->resolveAndBuildTransfer() < 0) {
          return -1;
        }
      }
    }
  }
  if (m_theMaterial != nullptr && m_hasInPlaneStiff) {
    Vector eps2d(3);
    this->computeMeanInPlaneStrain(eps2d);
    return m_theMaterial->setTrialStrain(eps2d);
  }
  return 0;
}

const Matrix &
LadrunoMembrane::getTangentStiff(void)
{
  if (!m_framesResolved) {
    this->resolveAndBuildTransfer();
  }
  if (!m_hasInPlaneStiff || m_numDOF == 0) {
    m_K.Zero();
    return m_K;
  }
  double scale = 1.0;
  if (m_theMaterial != nullptr && std::abs(m_D00_init) > 1e-15) {
    const Matrix &Dt = m_theMaterial->getTangent();
    if (Dt.noRows() >= 2 && Dt.noCols() >= 2) {
      scale = Dt(0, 0) / m_D00_init;
    }
  }
#ifdef LADRUNO_MUTATE_MEMBRANE
  scale = 0.0;
#endif
  m_K = m_Kbase;
  if (scale != 1.0) {
    m_K *= scale;
  }
  return m_K;
}

const Matrix &
LadrunoMembrane::getInitialStiff(void)
{
  if (Ki != 0) {
    return *Ki;
  }
  if (!m_framesResolved) {
    this->resolveAndBuildTransfer();
  }
  m_K = m_Kbase;
#ifdef LADRUNO_MUTATE_MEMBRANE
  m_K.Zero();
#endif
  return LadrunoShell::cacheKi(Ki, m_K);
}

const Matrix &
LadrunoMembrane::getMass(void)
{
  if (!m_framesResolved) {
    this->resolveAndBuildTransfer();
  }
  double rhoEff = (m_rho > 0.0) ? m_rho : ((m_theMaterial != nullptr) ? m_theMaterial->getRho() : 0.0);
  if (rhoEff <= 0.0 || m_thick <= 0.0 || m_numDOF == 0) {
    m_M.Zero();
    return m_M;
  }
  const double sig[3] = {rhoEff, m_thick, m_lumped ? 1.0 : 0.0};
  int nNodes = static_cast<int>(m_theNodes.size());
  if (const Matrix *hit = m_massCache.lookup(sig, 3, m_theNodes.data(), nNodes, m_ndm)) {
    return *hit;
  }
  double scale = (m_rhoInit > 0.0) ? (rhoEff / m_rhoInit) : 0.0;
  m_M = m_Mbase;
  if (scale != 1.0) {
    m_M *= scale;
  }
  m_massCache.fill(m_M, sig, 3, m_theNodes.data(), nNodes, m_ndm);
  if (const Matrix *c = m_massCache.cached()) {
    return *c;
  }
  return m_M;
}

void
LadrunoMembrane::zeroLoad(void)
{
  if (m_Q.Size() == m_numDOF) {
    m_Q.Zero();
  }
  m_activeBeamLoads.clear();
  m_appliedTotalLoad.Zero();
  m_appliedNormalPressure = 0.0;
  m_appliedGlobalTraction = {0.0, 0.0, 0.0};
}

int
LadrunoMembrane::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  if (!m_framesResolved) {
    if (this->resolveAndBuildTransfer() < 0) {
      return -1;
    }
  }

  int type = 0;
  const Vector &data = theLoad->getData(type, loadFactor);
  double qx = 0.0, qy = 0.0, qz = 0.0;

  if (type == LOAD_TAG_SelfWeight) {
    double wBase = 1.0;
    double rhoEff = (m_rho > 0.0) ? m_rho : ((m_theMaterial != nullptr) ? m_theMaterial->getRho() : 0.0);
    if (rhoEff > 0.0 && m_thick > 0.0) {
      wBase = rhoEff * m_thick;
    } else if (m_pressure != 0.0) {
      wBase = m_pressure;
    }
    double fx = (data.Size() >= 1) ? data(0) : 0.0;
    double fy = (data.Size() >= 2) ? data(1) : 0.0;
    double fz = (data.Size() >= 3) ? data(2) : 0.0;
    qx = loadFactor * wBase * fx;
    qy = loadFactor * wBase * fy;
    qz = loadFactor * wBase * fz;
  } else if (type == LOAD_TAG_SurfaceLoader) {
    double pBase = (m_pressure != 0.0) ? m_pressure : 1.0;
    qx = loadFactor * (pBase * m_normal[0] + m_b1);
    qy = loadFactor * (pBase * m_normal[1] + m_b2);
    qz = loadFactor * (pBase * m_normal[2] + m_b3);
  } else if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = (data.Size() >= 1) ? data(0) : 0.0;
    double wz = (data.Size() >= 2) ? data(1) : 0.0;
    double wx = (data.Size() >= 3) ? data(2) : 0.0;
    qx = loadFactor * (wx * m_e1[0] + wy * m_e2[0] + wz * m_normal[0]);
    qy = loadFactor * (wx * m_e1[1] + wy * m_e2[1] + wz * m_normal[1]);
    qz = loadFactor * (wx * m_e1[2] + wy * m_e2[2] + wz * m_normal[2]);
  } else if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wy = (data.Size() >= 1) ? data(0) : 0.0;
    double wx = (data.Size() >= 2) ? data(1) : 0.0;
    qx = loadFactor * (wx * m_e1[0] + wy * m_e2[0]);
    qy = loadFactor * (wx * m_e1[1] + wy * m_e2[1]);
    qz = 0.0;
  } else {
    opserr << "LadrunoMembrane::addLoad - load type " << type
           << " unknown for element " << this->getTag() << "\n";
    return -1;
  }

#ifdef LADRUNO_MUTATE_MEMBRANE
  qx = 0.0;
  qy = 0.0;
  qz = 0.0;
#endif

  double qn = qx * m_normal[0] + qy * m_normal[1] + qz * m_normal[2];
  m_appliedNormalPressure += qn;
  m_appliedGlobalTraction[0] += qx;
  m_appliedGlobalTraction[1] += qy;
  m_appliedGlobalTraction[2] += qz;

  // 1. Emit scoped partial trapezoidal beam loads onto supporting frames
  for (const auto &seg : m_unitSubLoads) {
    if (seg.frameIdx < 0 || seg.frameIdx >= static_cast<int>(m_frameInfos.size())) {
      continue;
    }
    const auto &fInfo = m_frameInfos[seg.frameIdx];
    if (fInfo.framePtr == nullptr) {
      continue;
    }

    double wxa_g = seg.waUnit * qx;
    double wya_g = seg.waUnit * qy;
    double wza_g = seg.waUnit * qz;
    double wxb_g = seg.wbUnit * qx;
    double wyb_g = seg.wbUnit * qy;
    double wzb_g = seg.wbUnit * qz;

    double wxa = wxa_g * fInfo.vx[0] + wya_g * fInfo.vx[1] + wza_g * fInfo.vx[2];
    double wya = wxa_g * fInfo.vy[0] + wya_g * fInfo.vy[1] + wza_g * fInfo.vy[2];
    double wza = wxa_g * fInfo.vz[0] + wya_g * fInfo.vz[1] + wza_g * fInfo.vz[2];

    double wxb = wxb_g * fInfo.vx[0] + wyb_g * fInfo.vx[1] + wzb_g * fInfo.vx[2];
    double wyb = wxb_g * fInfo.vy[0] + wyb_g * fInfo.vy[1] + wzb_g * fInfo.vy[2];
    double wzb = wxb_g * fInfo.vz[0] + wyb_g * fInfo.vz[1] + wzb_g * fInfo.vz[2];

    if (m_ndm == 3) {
      auto ptr = std::make_unique<Beam3dPartialUniformLoad>(
          0, wya, wza, wxa, seg.aOverL, seg.bOverL, wyb, wzb, wxb, fInfo.frameTag);
      ElementalLoad *rawPtr = ptr.get();
      m_activeBeamLoads.push_back(std::move(ptr));
      fInfo.framePtr->addLoad(rawPtr, 1.0);
    } else {
      auto ptr = std::make_unique<Beam2dPartialUniformLoad>(
          0, wya, wyb, wxa, wxb, seg.aOverL, seg.bOverL, fInfo.frameTag);
      ElementalLoad *rawPtr = ptr.get();
      m_activeBeamLoads.push_back(std::move(ptr));
      fInfo.framePtr->addLoad(rawPtr, 1.0);
    }

    double areaSeg = 0.5 * (seg.waUnit + seg.wbUnit) * seg.length;
    m_appliedTotalLoad(0) += areaSeg * qx;
    m_appliedTotalLoad(1) += areaSeg * qy;
    m_appliedTotalLoad(2) += areaSeg * qz;
  }

  // 2. Accumulate point support reactions into m_Q
  bool hasExplicitSupports = (!m_frameTags.empty() || !m_supportNodeTags.empty());
  if (hasExplicitSupports) {
    for (std::size_t s = 0; s < m_supportExtIdx.size() && s < m_unitSupportAreas.size(); ++s) {
      double areaS = m_unitSupportAreas[s];
      int extIdx = m_supportExtIdx[s];
      int off = m_nodeDofOffsets[extIdx];
      m_Q(off + 0) += areaS * qx;
      m_Q(off + 1) += areaS * qy;
      if (m_ndm == 3) {
        m_Q(off + 2) += areaS * qz;
      }
      m_appliedTotalLoad(0) += areaS * qx;
      m_appliedTotalLoad(1) += areaS * qy;
      m_appliedTotalLoad(2) += areaS * qz;
    }
  } else {
    int nPoly = static_cast<int>(m_polyNodeTags.size());
    for (int i = 0; i < nPoly && i < static_cast<int>(m_unitSupportAreas.size()); ++i) {
      double areaI = m_unitSupportAreas[i];
      int off = m_nodeDofOffsets[i];
      m_Q(off + 0) += areaI * qx;
      m_Q(off + 1) += areaI * qy;
      if (m_ndm == 3) {
        m_Q(off + 2) += areaI * qz;
      }
      m_appliedTotalLoad(0) += areaI * qx;
      m_appliedTotalLoad(1) += areaI * qy;
      m_appliedTotalLoad(2) += areaI * qz;
    }
  }

  // 3. Transfer any rank < 2 cantilever overhang torque to frame end nodes' rotational DOFs
  if (m_ndm == 3 && m_supportRank < 2) {
    for (std::size_t k = 0; k < m_frameInfos.size(); ++k) {
      double torq = m_frameInfos[k].overhangTorqueUnit * qn;
      if (std::abs(torq) <= 1e-15) {
        continue;
      }
      const auto &vx = m_frameInfos[k].vx;
      std::vector<int> rotRecipients;
      if (m_frameInfos[k].extNodeI >= 0 &&
          m_theNodes[m_frameInfos[k].extNodeI]->getNumberDOF() >= 6) {
        rotRecipients.push_back(m_frameInfos[k].extNodeI);
      }
      if (m_frameInfos[k].extNodeJ >= 0 &&
          m_theNodes[m_frameInfos[k].extNodeJ]->getNumberDOF() >= 6) {
        rotRecipients.push_back(m_frameInfos[k].extNodeJ);
      }
      if (!rotRecipients.empty()) {
        double share = torq / static_cast<double>(rotRecipients.size());
        for (int extIdx : rotRecipients) {
          int off = m_nodeDofOffsets[extIdx];
          m_Q(off + 3) += share * vx[0];
          m_Q(off + 4) += share * vx[1];
          m_Q(off + 5) += share * vx[2];
        }
      }
    }
  }

  return 0;
}

int
LadrunoMembrane::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (!this->hasMass() || m_numDOF == 0) {
    return 0;
  }
  const Matrix &Mq = this->getMass();
  Vector ra(m_numDOF);
  int numNodes = static_cast<int>(m_theNodes.size());
  for (int a = 0; a < numNodes; ++a) {
    const Vector &Raccel = m_theNodes[a]->getRV(accel);
    int ndf_a = m_theNodes[a]->getNumberDOF();
    int off_a = m_nodeDofOffsets[a];
    if (Raccel.Size() != ndf_a) {
      opserr << "LadrunoMembrane::addInertiaLoadToUnbalance - matrix and vector sizes incompatible\n";
      return -1;
    }
    for (int j = 0; j < ndf_a; ++j) {
      ra(off_a + j) = Raccel(j);
    }
  }
  if (m_lumped) {
    for (int i = 0; i < m_numDOF; ++i) {
      m_Q(i) += -Mq(i, i) * ra(i);
    }
  } else {
    m_Q.addMatrixVector(1.0, Mq, ra, -1.0);
  }
  return 0;
}

const Vector &
LadrunoMembrane::getResistingForce(void)
{
  if (!m_framesResolved) {
    this->resolveAndBuildTransfer();
  }
  m_P.Zero();
  if (m_hasInPlaneStiff && m_numDOF > 0) {
    Vector u(m_numDOF);
    int numNodes = static_cast<int>(m_theNodes.size());
    for (int a = 0; a < numNodes; ++a) {
      const Vector &disp = m_theNodes[a]->getTrialDisp();
      int ndf_a = m_theNodes[a]->getNumberDOF();
      int off_a = m_nodeDofOffsets[a];
      for (int j = 0; j < ndf_a && j < disp.Size(); ++j) {
        u(off_a + j) = disp(j);
      }
    }
    double scale = 1.0;
    if (m_theMaterial != nullptr && std::abs(m_D00_init) > 1e-15) {
      const Matrix &Dt = m_theMaterial->getTangent();
      if (Dt.noRows() >= 2 && Dt.noCols() >= 2) {
        scale = Dt(0, 0) / m_D00_init;
      }
    }
#ifdef LADRUNO_MUTATE_MEMBRANE
    scale = 0.0;
#endif
    if (scale != 0.0) {
      m_P.addMatrixVector(0.0, m_Kbase, u, scale);
    }
  }
  m_P.addVector(1.0, m_Q, -1.0);
  return m_P;
}

const Vector &
LadrunoMembrane::getResistingForceIncInertia(void)
{
  this->getResistingForce();
  Vector res(m_numDOF);
  res = m_P;

  if (this->hasMass() && m_numDOF > 0) {
    const Matrix &Mq = this->getMass();
    Vector a(m_numDOF);
    int numNodes = static_cast<int>(m_theNodes.size());
    for (int n = 0; n < numNodes; ++n) {
      const Vector &accel = m_theNodes[n]->getTrialAccel();
      int ndf_n = m_theNodes[n]->getNumberDOF();
      int off_n = m_nodeDofOffsets[n];
      for (int j = 0; j < ndf_n && j < accel.Size(); ++j) {
        a(off_n + j) = accel(j);
      }
    }
    if (m_lumped) {
      for (int i = 0; i < m_numDOF; ++i) {
        res(i) += Mq(i, i) * a(i);
      }
    } else {
      res.addMatrixVector(1.0, Mq, a, 1.0);
    }
  }

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
    res += this->getRayleighDampingForces();
  }
  m_P = res;
  return m_P;
}

int
LadrunoMembrane::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  int methodCode = 0;
  if (m_method == "oneway") methodCode = 1;
  else if (m_method == "fea_poisson") methodCode = 2;
  else if (m_method == "fea_plate") methodCode = 3;

  int totalHoleNodes = 0;
  for (const auto &hl : m_holeLoops) {
    totalHoleNodes += static_cast<int>(hl.size());
  }

  Vector data(28);   // WP-204: +4 (f11, f22, f12, #diaphragm nodes)
  data(0)  = this->getTag();
  data(1)  = static_cast<int>(m_polyNodeTags.size());
  data(2)  = static_cast<int>(m_frameTags.size());
  data(3)  = static_cast<int>(m_supportNodeTags.size());
  data(4)  = static_cast<int>(m_holeLoops.size());
  data(5)  = totalHoleNodes;
  data(6)  = methodCode;
  data(7)  = m_hasSpanAxis ? 1.0 : 0.0;
  data(8)  = m_spanAxis[0];
  data(9)  = m_spanAxis[1];
  data(10) = m_spanAxis[2];
  data(11) = m_orthoRatio;
  data(12) = m_segments;
  data(13) = m_meshSize;
  data(14) = m_thick;
  data(15) = m_E;
  data(16) = m_nu;
  data(17) = m_rho;
  data(18) = m_lumped ? 1.0 : 0.0;
  data(19) = m_pressure;
  data(20) = m_b1;
  data(21) = m_b2;
  data(22) = m_b3;
  data(23) = (m_theMaterial != nullptr) ? 1.0 : 0.0;
  data(24) = m_inPlaneMods[0];
  data(25) = m_inPlaneMods[1];
  data(26) = m_inPlaneMods[2];
  data(27) = static_cast<int>(m_diaphragmNodeTags.size());

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoMembrane::sendSelf - failed to send Vector\n";
    return res;
  }

  int matClassTag = 0;
  int matDbTag = 0;
  if (m_theMaterial != nullptr) {
    matClassTag = m_theMaterial->getClassTag();
    matDbTag = m_theMaterial->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0) {
        m_theMaterial->setDbTag(matDbTag);
      }
    }
  }

  int nPoly = static_cast<int>(m_polyNodeTags.size());
  int nFrames = static_cast<int>(m_frameTags.size());
  int nSupp = static_cast<int>(m_supportNodeTags.size());
  int nHoles = static_cast<int>(m_holeLoops.size());
  int nDiaph = static_cast<int>(m_diaphragmNodeTags.size());
  int idLen = 2 + nPoly + nFrames + nSupp + nHoles + totalHoleNodes + nDiaph;

  ID idData(idLen);
  int idx = 0;
  idData(idx++) = matClassTag;
  idData(idx++) = matDbTag;
  for (int t : m_polyNodeTags) idData(idx++) = t;
  for (int t : m_frameTags) idData(idx++) = t;
  for (int t : m_supportNodeTags) idData(idx++) = t;
  for (const auto &hl : m_holeLoops) idData(idx++) = static_cast<int>(hl.size());
  for (const auto &hl : m_holeLoops) {
    for (int t : hl) idData(idx++) = t;
  }
  for (int t : m_diaphragmNodeTags) idData(idx++) = t;

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoMembrane::sendSelf - failed to send ID\n";
    return res;
  }

  if (m_theMaterial != nullptr) {
    res += m_theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LadrunoMembrane::sendSelf - failed to send material\n";
      return res;
    }
  }
  return res;
}

int
LadrunoMembrane::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  Vector data(28);   // WP-204
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LadrunoMembrane::recvSelf - failed to recv Vector\n";
    return res;
  }

  this->setTag(static_cast<int>(data(0)));
  int nPoly          = static_cast<int>(data(1));
  int nFrames        = static_cast<int>(data(2));
  int nSupp          = static_cast<int>(data(3));
  int nHoles         = static_cast<int>(data(4));
  int totalHoleNodes = static_cast<int>(data(5));
  int methodCode     = static_cast<int>(data(6));
  m_hasSpanAxis      = (data(7) != 0.0);
  m_spanAxis         = {data(8), data(9), data(10)};
  m_orthoRatio       = data(11);
  m_segments         = static_cast<int>(data(12));
  m_meshSize         = data(13);
  m_thick            = data(14);
  m_E                = data(15);
  m_nu               = data(16);
  m_rho              = data(17);
  m_lumped           = (data(18) != 0.0);
  m_pressure         = data(19);
  m_b1               = data(20);
  m_b2               = data(21);
  m_b3               = data(22);
  bool hasMat        = (data(23) != 0.0);
  m_inPlaneMods      = {data(24), data(25), data(26)};
  int nDiaph         = static_cast<int>(data(27));

  if (methodCode == 1) m_method = "oneway";
  else if (methodCode == 2) m_method = "fea_poisson";
  else if (methodCode == 3) m_method = "fea_plate";
  else m_method = "tributary";

  int idLen = 2 + nPoly + nFrames + nSupp + nHoles + totalHoleNodes + nDiaph;
  ID idData(idLen);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LadrunoMembrane::recvSelf - failed to recv ID\n";
    return res;
  }

  int idx = 0;
  int matClassTag = idData(idx++);
  int matDbTag    = idData(idx++);

  m_polyNodeTags.assign(nPoly, 0);
  for (int i = 0; i < nPoly; ++i) m_polyNodeTags[i] = idData(idx++);

  m_frameTags.assign(nFrames, 0);
  for (int i = 0; i < nFrames; ++i) m_frameTags[i] = idData(idx++);

  m_supportNodeTags.assign(nSupp, 0);
  for (int i = 0; i < nSupp; ++i) m_supportNodeTags[i] = idData(idx++);

  std::vector<int> holeSizes(nHoles, 0);
  for (int h = 0; h < nHoles; ++h) holeSizes[h] = idData(idx++);

  m_holeLoops.clear();
  m_holeLoops.resize(nHoles);
  for (int h = 0; h < nHoles; ++h) {
    m_holeLoops[h].assign(holeSizes[h], 0);
    for (int j = 0; j < holeSizes[h]; ++j) {
      m_holeLoops[h][j] = idData(idx++);
    }
  }

  m_diaphragmNodeTags.assign(nDiaph, 0);
  for (int i = 0; i < nDiaph; ++i) m_diaphragmNodeTags[i] = idData(idx++);

  if (hasMat) {
    if (m_theMaterial == nullptr || m_theMaterial->getClassTag() != matClassTag) {
      if (m_theMaterial != nullptr) {
        delete m_theMaterial;
      }
      m_theMaterial = theBroker.getNewNDMaterial(matClassTag);
      if (m_theMaterial == nullptr) {
        opserr << "WARNING LadrunoMembrane::recvSelf - failed to allocate NDMaterial classTag "
               << matClassTag << "\n";
        return -1;
      }
    }
    m_theMaterial->setDbTag(matDbTag);
    res += m_theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING LadrunoMembrane::recvSelf - failed to recv material\n";
      return res;
    }
  } else {
    if (m_theMaterial != nullptr) {
      delete m_theMaterial;
      m_theMaterial = nullptr;
    }
  }

  LadrunoShell::dropKi(Ki);
  m_massCache.invalidate();
  m_framesResolved = false;
  this->rebuildExternalNodeIndex();

  if (this->getDomain() != nullptr) {
    this->setDomain(this->getDomain());
    this->resolveAndBuildTransfer();
  }
  return res;
}

void
LadrunoMembrane::Print(OPS_Stream &s, int flag)
{
  s << "LadrunoMembrane, tag: " << this->getTag() << " method: " << m_method.c_str()
    << " area: " << m_slabArea << endln;
  s << "  Polygon nodes:";
  for (int t : m_polyNodeTags) s << " " << t;
  s << endln;
  if (!m_frameTags.empty()) {
    s << "  Supporting frames:";
    for (int t : m_frameTags) s << " " << t;
    s << endln;
  }
}

Response *
LadrunoMembrane::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType", "LadrunoMembrane");
  output.attr("eleTag", this->getTag());
  int numNodes = m_connectedExternalNodes.Size();
  for (int i = 0; i < numNodes; ++i) {
    char nodeName[32];
    snprintf(nodeName, sizeof(nodeName), "node%d", i + 1);
    output.attr(nodeName, m_connectedExternalNodes(i));
  }

  if (!m_framesResolved && this->getDomain() != nullptr) {
    this->resolveAndBuildTransfer();
  }

  if (argc > 0) {
    if (strcmp(argv[0], "totalLoad") == 0) {
      theResponse = new ElementResponse(this, 1, Vector(3));
    } else if (strcmp(argv[0], "frameLoads") == 0) {
      int nF = static_cast<int>(m_frameTags.size());
      theResponse = new ElementResponse(this, 2, Vector(nF > 0 ? 4 * nF : 1));
    } else if (strcmp(argv[0], "supportReactions") == 0 ||
               strcmp(argv[0], "nodalReactions") == 0) {
      int nS = static_cast<int>(m_supportNodeTags.size());
      theResponse = new ElementResponse(this, 3, Vector(nS > 0 ? nS : 1));
    } else if (strcmp(argv[0], "tributaryAreas") == 0 ||
               strcmp(argv[0], "tributaryWidths") == 0) {
      int nT = static_cast<int>(m_frameTags.size() + m_supportNodeTags.size());
      theResponse = new ElementResponse(this, 4, Vector(nT > 0 ? nT : 1));
    } else if (strcmp(argv[0], "localAxes") == 0) {
      theResponse = new ElementResponse(this, 5, Vector(9));
    } else if (strcmp(argv[0], "area") == 0) {
      theResponse = new ElementResponse(this, 6, Vector(1));
    } else if (strcmp(argv[0], "subSegmentLoads") == 0 ||
               strcmp(argv[0], "beamSubLoads") == 0) {
      int nSub = static_cast<int>(m_unitSubLoads.size());
      theResponse = new ElementResponse(this, 7, Vector(nSub > 0 ? 9 * nSub : 1));
    } else if (strcmp(argv[0], "subSegments") == 0) {
      theResponse = new ElementResponse(this, 8, Vector(1));
    } else if (strcmp(argv[0], "overhangTorques") == 0) {
      int nF = static_cast<int>(m_frameTags.size());
      theResponse = new ElementResponse(this, 9, Vector(nF > 0 ? nF : 1));
    } else if (strcmp(argv[0], "meshStats") == 0) {
      theResponse = new ElementResponse(this, 10, Vector(4));
    } else if (strcmp(argv[0], "meshVertices") == 0) {
      int nV3 = static_cast<int>(m_meshVertices3d.size());
      theResponse = new ElementResponse(this, 11, Vector(nV3 > 0 ? nV3 : 1));
    } else if (strcmp(argv[0], "meshTriangles") == 0) {
      int nT3 = static_cast<int>(m_meshTriangles.size());
      theResponse = new ElementResponse(this, 12, Vector(nT3 > 0 ? nT3 : 1));
    } else if (strcmp(argv[0], "frameTags") == 0) {
      int nF = static_cast<int>(m_frameTags.size());
      theResponse = new ElementResponse(this, 13, Vector(nF > 0 ? nF : 1));
    } else if (strcmp(argv[0], "supportNodeTags") == 0 ||
               strcmp(argv[0], "supportTags") == 0) {
      int nS = static_cast<int>(m_supportNodeTags.size());
      theResponse = new ElementResponse(this, 14, Vector(nS > 0 ? nS : 1));
    } else if (strcmp(argv[0], "meshSubSegments") == 0) {
      int nMSub = static_cast<int>(m_meshSubSegmentsUnit.size());
      theResponse = new ElementResponse(this, 15, Vector(nMSub > 0 ? nMSub : 1));
    } else if (strcmp(argv[0], "meshVertexFlags") == 0) {
      int nVF = static_cast<int>(m_meshVertexFlags.size());
      theResponse = new ElementResponse(this, 16, Vector(nVF > 0 ? nVF : 1));
    } else if (strcmp(argv[0], "meshNodalField") == 0 ||
               strcmp(argv[0], "meshField") == 0) {
      int nNF = static_cast<int>(m_meshNodalFieldUnit.size());
      theResponse = new ElementResponse(this, 17, Vector(nNF > 0 ? nNF : 1));
    } else if (strcmp(argv[0], "meshTriangleLoadFlow") == 0 ||
               strcmp(argv[0], "meshLoadVectors") == 0) {
      int nTF = static_cast<int>(m_meshTriangleLoadFlowUnit.size());
      theResponse = new ElementResponse(this, 18, Vector(nTF > 0 ? nTF : 1));
    } else if (strcmp(argv[0], "yieldLines") == 0) {
      int nYL = static_cast<int>(m_yieldLines3d.size());
      theResponse = new ElementResponse(this, 19, Vector(nYL > 0 ? nYL : 1));
    } else if (strcmp(argv[0], "membraneForces") == 0 || strcmp(argv[0], "membraneForce") == 0) {
      // WP-204: [N11 N22 N12] per diaphragm triangle (centroid, e1/e2 frame)
      int nT = static_cast<int>(m_dTriArea.size());
      for (int t = 0; t < nT; ++t) {
        char nm[48];
        snprintf(nm, sizeof(nm), "N11_%d", t + 1); output.tag("ResponseType", nm);
        snprintf(nm, sizeof(nm), "N22_%d", t + 1); output.tag("ResponseType", nm);
        snprintf(nm, sizeof(nm), "N12_%d", t + 1); output.tag("ResponseType", nm);
      }
      theResponse = new ElementResponse(this, 20, Vector(nT > 0 ? 3 * nT : 1));
    } else if (strcmp(argv[0], "diaphragmVertices") == 0) {
      int n = static_cast<int>(m_dVerts3d.size());
      theResponse = new ElementResponse(this, 21, Vector(n > 0 ? n : 1));
    } else if (strcmp(argv[0], "diaphragmTriangles") == 0) {
      int n = static_cast<int>(m_dTris.size());
      theResponse = new ElementResponse(this, 22, Vector(n > 0 ? n : 1));
    } else if (strcmp(argv[0], "meanMembraneForce") == 0 || strcmp(argv[0], "meanMembraneForces") == 0) {
      output.tag("ResponseType", "N11");
      output.tag("ResponseType", "N22");
      output.tag("ResponseType", "N12");
      theResponse = new ElementResponse(this, 23, Vector(3));
    } else if (strcmp(argv[0], "diaphragmForces") == 0 || strcmp(argv[0], "diaphragmForce") == 0) {
      // WP-204: in-plane (diaphragm) force on each node, global components
      static const char *comp[3] = {"Px", "Py", "Pz"};
      for (int a = 0; a < numNodes; ++a) {
        for (int i = 0; i < m_ndm; ++i) {
          char nm[48];
          snprintf(nm, sizeof(nm), "%s_%d", comp[i], a + 1);
          output.tag("ResponseType", nm);
        }
      }
      theResponse = new ElementResponse(this, 24, Vector(numNodes * m_ndm > 0 ? numNodes * m_ndm : 1));
    } else if ((strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) &&
               m_theMaterial != nullptr && argc > 2) {
      theResponse = m_theMaterial->setResponse(&argv[2], argc - 2, output);
    } else if ((strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0 ||
                strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) &&
               m_theMaterial != nullptr) {
      theResponse = m_theMaterial->setResponse(argv, argc, output);
    }
  }

  return LadrunoShell::finishResponse(this, theResponse, argv, argc, output);
}

int
LadrunoMembrane::getResponse(int responseID, Information &eleInfo)
{
  if (!m_framesResolved && this->getDomain() != nullptr) {
    this->resolveAndBuildTransfer();
  }

  switch (responseID) {
  case 1:
    return eleInfo.setVector(m_appliedTotalLoad);

  case 2: {
    int nF = static_cast<int>(m_frameInfos.size());
    static thread_local Vector vFrameLoads;
    vFrameLoads.resize(nF > 0 ? 4 * nF : 1);
    vFrameLoads.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? m_appliedNormalPressure : 1.0;
    for (int k = 0; k < nF; ++k) {
      vFrameLoads(4 * k + 0) = m_frameInfos[k].aL;
      vFrameLoads(4 * k + 1) = m_frameInfos[k].bL;
      vFrameLoads(4 * k + 2) = m_frameInfos[k].waEqUnit * qScale;
      vFrameLoads(4 * k + 3) = m_frameInfos[k].wbEqUnit * qScale;
    }
    return eleInfo.setVector(vFrameLoads);
  }

  case 3: {
    int nS = static_cast<int>(m_supportNodeTags.size());
    static thread_local Vector vSupp;
    vSupp.resize(nS > 0 ? nS : 1);
    vSupp.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? m_appliedNormalPressure : 1.0;
    for (int s = 0; s < nS && s < static_cast<int>(m_unitSupportAreas.size()); ++s) {
      vSupp(s) = m_unitSupportAreas[s] * qScale;
    }
    return eleInfo.setVector(vSupp);
  }

  case 4: {
    int nF = static_cast<int>(m_frameInfos.size());
    int nS = static_cast<int>(m_supportNodeTags.size());
    int nT = nF + nS;
    static thread_local Vector vTrib;
    vTrib.resize(nT > 0 ? nT : 1);
    vTrib.Zero();
    for (int k = 0; k < nF; ++k) {
      vTrib(k) = m_frameInfos[k].areaUnit;
    }
    for (int s = 0; s < nS && s < static_cast<int>(m_unitSupportAreas.size()); ++s) {
      vTrib(nF + s) = m_unitSupportAreas[s];
    }
    if (nT == 0) {
      vTrib(0) = m_slabArea;
    }
    return eleInfo.setVector(vTrib);
  }

  case 5: {
    static thread_local Vector la(9);
    for (int i = 0; i < 3; ++i) {
      la(i)     = m_e1[i];
      la(i + 3) = m_e2[i];
      la(i + 6) = m_normal[i];
    }
    return eleInfo.setVector(la);
  }

  case 6: {
    static thread_local Vector vArea(1);
    vArea(0) = m_slabArea;
    return eleInfo.setVector(vArea);
  }

  case 7: {
    int nSub = static_cast<int>(m_unitSubLoads.size());
    static thread_local Vector vSubLoads;
    vSubLoads.resize(nSub > 0 ? 9 * nSub : 1);
    vSubLoads.Zero();
    double qNormSq = m_appliedGlobalTraction[0] * m_appliedGlobalTraction[0] +
                     m_appliedGlobalTraction[1] * m_appliedGlobalTraction[1] +
                     m_appliedGlobalTraction[2] * m_appliedGlobalTraction[2];
    double qx = (qNormSq > 1e-30) ? m_appliedGlobalTraction[0] : m_normal[0];
    double qy = (qNormSq > 1e-30) ? m_appliedGlobalTraction[1] : m_normal[1];
    double qz = (qNormSq > 1e-30) ? m_appliedGlobalTraction[2] : m_normal[2];
    for (int s = 0; s < nSub; ++s) {
      const auto &seg = m_unitSubLoads[s];
      int fTag = 0;
      std::array<double, 3> vx = {1.0, 0.0, 0.0};
      std::array<double, 3> vy = {0.0, 1.0, 0.0};
      std::array<double, 3> vz = {0.0, 0.0, 1.0};
      if (seg.frameIdx >= 0 && seg.frameIdx < static_cast<int>(m_frameInfos.size())) {
        const auto &fInfo = m_frameInfos[seg.frameIdx];
        fTag = fInfo.frameTag;
        vx = fInfo.vx;
        vy = fInfo.vy;
        vz = fInfo.vz;
      }
      double wxa_g = seg.waUnit * qx;
      double wya_g = seg.waUnit * qy;
      double wza_g = seg.waUnit * qz;
      double wxb_g = seg.wbUnit * qx;
      double wyb_g = seg.wbUnit * qy;
      double wzb_g = seg.wbUnit * qz;

      double wxa = wxa_g * vx[0] + wya_g * vx[1] + wza_g * vx[2];
      double wya = wxa_g * vy[0] + wya_g * vy[1] + wza_g * vy[2];
      double wza = (m_ndm == 3) ? (wxa_g * vz[0] + wya_g * vz[1] + wza_g * vz[2]) : 0.0;

      double wxb = wxb_g * vx[0] + wyb_g * vx[1] + wzb_g * vx[2];
      double wyb = wxb_g * vy[0] + wyb_g * vy[1] + wzb_g * vy[2];
      double wzb = (m_ndm == 3) ? (wxb_g * vz[0] + wyb_g * vz[1] + wzb_g * vz[2]) : 0.0;

      vSubLoads(9 * s + 0) = static_cast<double>(fTag);
      vSubLoads(9 * s + 1) = seg.aOverL;
      vSubLoads(9 * s + 2) = seg.bOverL;
      vSubLoads(9 * s + 3) = wya;
      vSubLoads(9 * s + 4) = wyb;
      vSubLoads(9 * s + 5) = wza;
      vSubLoads(9 * s + 6) = wzb;
      vSubLoads(9 * s + 7) = wxa;
      vSubLoads(9 * s + 8) = wxb;
    }
    return eleInfo.setVector(vSubLoads);
  }

  case 8: {
    static thread_local Vector vSubCount(1);
    vSubCount(0) = static_cast<double>(m_unitSubLoads.size());
    return eleInfo.setVector(vSubCount);
  }

  case 9: {
    int nF = static_cast<int>(m_frameInfos.size());
    static thread_local Vector vTorq;
    vTorq.resize(nF > 0 ? nF : 1);
    vTorq.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? m_appliedNormalPressure : 1.0;
    for (int k = 0; k < nF; ++k) {
      vTorq(k) = m_frameInfos[k].overhangTorqueUnit * qScale;
    }
    return eleInfo.setVector(vTorq);
  }

  case 10: {
    static thread_local Vector vStats(4);
    vStats(0) = static_cast<double>(m_meshVertices3d.size() / 3);
    vStats(1) = static_cast<double>(m_meshTriangles.size() / 3);
    vStats(2) = static_cast<double>(m_unitSubLoads.size());
    vStats(3) = static_cast<double>(m_supportRank);
    return eleInfo.setVector(vStats);
  }

  case 11: {
    int nV3 = static_cast<int>(m_meshVertices3d.size());
    static thread_local Vector vVerts;
    vVerts.resize(nV3 > 0 ? nV3 : 1);
    vVerts.Zero();
    for (int i = 0; i < nV3; ++i) {
      vVerts(i) = m_meshVertices3d[i];
    }
    return eleInfo.setVector(vVerts);
  }

  case 12: {
    int nT3 = static_cast<int>(m_meshTriangles.size());
    static thread_local Vector vTris;
    vTris.resize(nT3 > 0 ? nT3 : 1);
    vTris.Zero();
    for (int i = 0; i < nT3; ++i) {
      vTris(i) = static_cast<double>(m_meshTriangles[i]);
    }
    return eleInfo.setVector(vTris);
  }

  case 13: {
    int nF = static_cast<int>(m_frameTags.size());
    static thread_local Vector vFTags;
    vFTags.resize(nF > 0 ? nF : 1);
    vFTags.Zero();
    for (int i = 0; i < nF; ++i) {
      vFTags(i) = static_cast<double>(m_frameTags[i]);
    }
    return eleInfo.setVector(vFTags);
  }

  case 14: {
    int nS = static_cast<int>(m_supportNodeTags.size());
    static thread_local Vector vSTags;
    vSTags.resize(nS > 0 ? nS : 1);
    vSTags.Zero();
    for (int i = 0; i < nS; ++i) {
      vSTags(i) = static_cast<double>(m_supportNodeTags[i]);
    }
    return eleInfo.setVector(vSTags);
  }

  case 15: {
    int nMSub = static_cast<int>(m_meshSubSegmentsUnit.size());
    static thread_local Vector vMSub;
    vMSub.resize(nMSub > 0 ? nMSub : 1);
    vMSub.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? m_appliedNormalPressure : 1.0;
    int nSegs = nMSub / 6;
    for (int s = 0; s < nSegs; ++s) {
      vMSub(6 * s + 0) = m_meshSubSegmentsUnit[6 * s + 0];
      vMSub(6 * s + 1) = m_meshSubSegmentsUnit[6 * s + 1];
      vMSub(6 * s + 2) = m_meshSubSegmentsUnit[6 * s + 2];
      vMSub(6 * s + 3) = m_meshSubSegmentsUnit[6 * s + 3];
      vMSub(6 * s + 4) = m_meshSubSegmentsUnit[6 * s + 4];
      vMSub(6 * s + 5) = m_meshSubSegmentsUnit[6 * s + 5] * qScale;
    }
    return eleInfo.setVector(vMSub);
  }

  case 16: {
    int nVF = static_cast<int>(m_meshVertexFlags.size());
    static thread_local Vector vFlags;
    vFlags.resize(nVF > 0 ? nVF : 1);
    vFlags.Zero();
    for (int i = 0; i < nVF; ++i) {
      vFlags(i) = static_cast<double>(m_meshVertexFlags[i]);
    }
    return eleInfo.setVector(vFlags);
  }

  case 17: {
    int nNF = static_cast<int>(m_meshNodalFieldUnit.size());
    static thread_local Vector vField;
    vField.resize(nNF > 0 ? nNF : 1);
    vField.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? m_appliedNormalPressure : 1.0;
    for (int i = 0; i < nNF; ++i) {
      vField(i) = m_meshNodalFieldUnit[i] * qScale;
    }
    return eleInfo.setVector(vField);
  }

  case 18: {
    int nTF = static_cast<int>(m_meshTriangleLoadFlowUnit.size());
    static thread_local Vector vFlow;
    vFlow.resize(nTF > 0 ? nTF : 1);
    vFlow.Zero();
    double qScale = (std::abs(m_appliedNormalPressure) > 1e-15) ? std::abs(m_appliedNormalPressure) : 1.0;
    int nTris = nTF / 4;
    for (int t = 0; t < nTris; ++t) {
      vFlow(4 * t + 0) = m_meshTriangleLoadFlowUnit[4 * t + 0];
      vFlow(4 * t + 1) = m_meshTriangleLoadFlowUnit[4 * t + 1] * qScale;
      vFlow(4 * t + 2) = m_meshTriangleLoadFlowUnit[4 * t + 2] * qScale;
      vFlow(4 * t + 3) = m_meshTriangleLoadFlowUnit[4 * t + 3] * qScale;
    }
    return eleInfo.setVector(vFlow);
  }

  case 19: {
    int nYL = static_cast<int>(m_yieldLines3d.size());
    static thread_local Vector vYL;
    vYL.resize(nYL > 0 ? nYL : 1);
    vYL.Zero();
    for (int i = 0; i < nYL; ++i) {
      vYL(i) = m_yieldLines3d[i];
    }
    return eleInfo.setVector(vYL);
  }

  case 20:
  case 23: {
    // WP-204: N = t * D_mod * eps (times the material tangent ratio)
    int n2 = 2 * static_cast<int>(m_theNodes.size());
    int nT = static_cast<int>(m_dTriArea.size());
    static thread_local Vector vN;
    vN.resize(responseID == 20 ? (nT > 0 ? 3 * nT : 1) : 3);
    vN.Zero();
    if (!m_hasInPlaneStiff || m_dStrainOp.size() != static_cast<std::size_t>(nT) * 3 * n2) {
      return eleInfo.setVector(vN);
    }
    std::vector<double> ur;
    this->gatherInPlaneDisp(ur);
    double ts = m_thick * this->tangentScale();
    auto forces = [&](const double *op, int out) {
      double eps[3];
      for (int r = 0; r < 3; ++r) {
        double s = 0.0;
        for (int j = 0; j < n2; ++j) s += op[r * n2 + j] * ur[j];
        eps[r] = s;
      }
      for (int r = 0; r < 3; ++r)
        vN(out + r) = ts * (m_Dmod[3 * r] * eps[0] + m_Dmod[3 * r + 1] * eps[1] + m_Dmod[3 * r + 2] * eps[2]);
    };
    if (responseID == 20) {
      for (int t = 0; t < nT; ++t)
        forces(&m_dStrainOp[static_cast<std::size_t>(t) * 3 * n2], 3 * t);
    } else if (m_dMeanStrainOp.size() == static_cast<std::size_t>(3 * n2)) {
      forces(m_dMeanStrainOp.data(), 0);
    }
    return eleInfo.setVector(vN);
  }

  case 21: {
    int n = static_cast<int>(m_dVerts3d.size());
    static thread_local Vector v;
    v.resize(n > 0 ? n : 1);
    v.Zero();
    for (int i = 0; i < n; ++i) v(i) = m_dVerts3d[i];
    return eleInfo.setVector(v);
  }

  case 22: {
    int n = static_cast<int>(m_dTris.size());
    static thread_local Vector v;
    v.resize(n > 0 ? n : 1);
    v.Zero();
    for (int i = 0; i < n; ++i) v(i) = static_cast<double>(m_dTris[i]);
    return eleInfo.setVector(v);
  }

  case 24: {
    int nExt = static_cast<int>(m_theNodes.size());
    static thread_local Vector v;
    v.resize(nExt * m_ndm > 0 ? nExt * m_ndm : 1);
    v.Zero();
    if (m_hasInPlaneStiff && m_numDOF > 0) {
      Vector u(m_numDOF);
      for (int a = 0; a < nExt; ++a) {
        const Vector &disp = m_theNodes[a]->getTrialDisp();
        int off = m_nodeDofOffsets[a];
        for (int j = 0; j < disp.Size() && off + j < m_numDOF; ++j) u(off + j) = disp(j);
      }
      Vector f(m_numDOF);
      f.addMatrixVector(0.0, m_Kbase, u, this->tangentScale());
      for (int a = 0; a < nExt; ++a)
        for (int i = 0; i < m_ndm; ++i)
          v(a * m_ndm + i) = f(m_nodeDofOffsets[a] + i);
    }
    return eleInfo.setVector(v);
  }

  default:
    return this->Element::getResponse(responseID, eleInfo);
  }
}

int
LadrunoMembrane::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1) {
    return -1;
  }
  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(m_rho);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0], "pressure") == 0) {
    param.setValue(m_pressure);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0], "thickness") == 0 || strcmp(argv[0], "thick") == 0) {
    param.setValue(m_thick);
    return param.addObject(3, this);
  }
  if (m_theMaterial != nullptr) {
    if (LadrunoShell::isMaterialPointToken(argv[0])) {
      return LadrunoShell::forwardToMaterialPoint(&m_theMaterial, 1, true, argv, argc, param);
    }
    return LadrunoShell::forwardToMaterials(&m_theMaterial, 1, argv, argc, param);
  }
  return -1;
}

int
LadrunoMembrane::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    m_rho = info.theDouble;
    m_massCache.invalidate();
    m_framesResolved = false;
    return 0;
  case 2:
    m_pressure = info.theDouble;
    return 0;
  case 3:
    m_thick = info.theDouble;
    LadrunoShell::dropKi(Ki);
    m_massCache.invalidate();
    m_framesResolved = false;
    return 0;
  default:
    return -1;
  }
}
