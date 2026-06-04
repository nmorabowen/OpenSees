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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// LadrunoIMKBeam2d: 2D concentrated-plasticity beam-column macro element.
// Planar counterpart of LadrunoIMKBeam (3D). Elastic interior in series with
// UNCOUPLED moment-rotation rotational hinges at the ends (one independent
// uniaxial law per end, strong-axis Mz only). No P-M interaction -> immune to
// rigid-diaphragm spurious axial force. Displacement-driven stiffness
// formulation: basic 3x3, exact series F = F_elastic + F_hinge (no n-factor),
// 2x2 internal Newton on the hinge rotations. Column-face offset and geometric
// nonlinearity are delegated to the CrdTransf2d (use geomTransf -jntOffset for
// the rigid offset to the hinge). The per-axis hinge kernel is shared with the
// 3D element via LadrunoIMKHinge.h.
//
// See Ladruno_implementation/14_ladruno_imk_beam.md.

#ifndef LadrunoIMKBeam2d_h
#define LadrunoIMKBeam2d_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class Channel;
class Information;
class Response;
class CrdTransf;
class UniaxialMaterial;
class FEM_ObjectBroker;

class LadrunoIMKBeam2d : public Element
{
 public:
  LadrunoIMKBeam2d();
  LadrunoIMKBeam2d(int tag, int Nd1, int Nd2,
                   double A, double E, double Iz,
                   CrdTransf &theTransf,
                   UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                   double rho = 0.0);
  ~LadrunoIMKBeam2d();

  const char *getClassType(void) const { return "LadrunoIMKBeam2d"; }

  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);
  void setDomain(Domain *theDomain);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  int update(void);
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);

 private:
  // --- geometry / elastic properties of the member ---
  double A, E, Iz;
  double rho;

  // --- hinge materials: [0]=Mz@i, [1]=Mz@j ---
  //     a NULL entry means that end is elastic (no hinge).
  UniaxialMaterial *theMat[2];
  double thetaH[2];        // trial hinge rotations
  double thetaHcommit[2];  // committed hinge rotations

  // --- cached basic-system state, set in update() ---
  Vector q;   // basic forces  [N, Mz_i, Mz_j]
  Matrix kb;  // basic tangent (3x3)
  Vector Q;   // 6-dof element load (inertia unbalance only in v1)

  Node *theNodes[2];
  ID connectedExternalNodes;
  CrdTransf *theCoordTransf;

  static Matrix K;  // 6x6 work matrix (mass / inertia)
  static Vector P;  // 6 work vector

  // per-axis state determination delegated to the shared, dimension-agnostic
  // kernel in LadrunoIMKHinge.h (ladrunoIMKSolveAxis / ladrunoIMKInitBlock).
};

#endif
