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
// LadrunoIMKBeam: 3D concentrated-plasticity beam-column macro element.
// Elastic interior in series with UNCOUPLED moment-rotation rotational hinges
// at the ends (one independent uniaxial law per bending axis, per end). No P-M
// and no biaxial My-Mz interaction -> immune to rigid-diaphragm spurious axial
// force. Displacement-driven stiffness formulation: basic 6x6, exact series
// F = F_elastic + F_hinge (no n-factor), per-axis 2x2 internal Newton on hinge
// rotations. Column-face offset and geometric nonlinearity are delegated to the
// CrdTransf3d (use geomTransf -jntOffset for the rigid offset to the hinge).
//
// See Ladruno_implementation/14_ladruno_imk_beam.md.

#ifndef LadrunoIMKBeam_h
#define LadrunoIMKBeam_h

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
class NDMaterial;
class FEM_ObjectBroker;

class LadrunoIMKBeam : public Element
{
 public:
  LadrunoIMKBeam();
  LadrunoIMKBeam(int tag, int Nd1, int Nd2,
                 double A, double E, double G, double Jx, double Iy, double Iz,
                 CrdTransf &theTransf,
                 UniaxialMaterial *mZi, UniaxialMaterial *mZj,
                 UniaxialMaterial *mYi, UniaxialMaterial *mYj,
                 double rho = 0.0, int cMass = 0,
                 double GAvy = 0.0, double GAvz = 0.0,
                 double P0 = 0.0, double pmAlpha = 0.5,
                 double dI = 0.0, double dJ = 0.0, bool useImplex = false,
                 UniaxialMaterial *mShearY = 0, UniaxialMaterial *mShearZ = 0,
                 NDMaterial *mBiI = 0, NDMaterial *mBiJ = 0,
                 double P0t = 0.0, double Pbal = 0.0, double muBal = 1.0);
  ~LadrunoIMKBeam();

  const char *getClassType(void) const { return "LadrunoIMKBeam"; }

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

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  bool ladrunoThreadSafeUpdate(void) const;
  bool ladrunoThreadSafeForm(void) const;
  bool ladrunoEnergyRayleighOnly(void) const { return true; }

 private:
  // --- geometry / elastic properties of the member ---
  double A, E, G, Jx, Iy, Iz;
  double rho;
  int cMass;
  double GAvy, GAvz;
  double P0, pmAlpha;
  double P0t, Pbal, muBal;
  double dI, dJ;
  bool useImplex;
  double implexErr;

  // --- hinge materials: [0]=Mz@i, [1]=Mz@j, [2]=My@i, [3]=My@j ---
  //     a NULL entry means that end/axis is elastic (no hinge).
  UniaxialMaterial *theMat[4];
  NDMaterial *theBiaxialMat[2];  // optional coupled 2-DOF (Mz, My) hinges at [0]=i, [1]=j
  double thetaH[4];         // trial hinge rotations
  double thetaHcommit[4];   // committed hinge rotations
  double dThetaHprev[4];    // previous committed increment
  double thetaHp[4];        // trial plastic hinge rotations (th - q_c/k0)
  double thetaHpCommit[4];  // committed plastic hinge rotations
  double dThetaHpPrev[4];   // previous committed plastic increment for IMPL-EX
  UniaxialMaterial *theShearMat[2];  // optional series shear hinges: [0]=Vy, [1]=Vz
  double gammaS[2];         // trial shear hinge deformations [gamma_y, gamma_z]
  double gammaScommit[2];   // committed shear hinge deformations
  double dGammaSprev[2];    // previous committed shear increments
  double gammaSp[2];        // trial plastic shear deformations
  double gammaSpCommit[2];  // committed plastic shear deformations
  double dGammaSpPrev[2];   // previous committed plastic shear increments for IMPL-EX
  double dtime_n_commit;    // previous committed time step for IMPL-EX
  bool stepDtLatched;       // latches step dt on first update() of each step
  double dtStep;            // latched step dt
  double dtScaleStep;       // latched dtScale for IMPL-EX within current step

  // --- cached basic-system state, set in update() ---
  Vector q;   // basic forces  [N, Mz_i, Mz_j, My_i, My_j, T]
  Matrix kb;  // basic tangent (6x6)
  Vector Q;   // 12-dof element load (inertia unbalance)
  double q0[5];  // Fixed end forces in basic system
  double p0[5];  // Reactions in basic system
  double vCommit[6];  // Committed basic deformations
  bool q0Dirty;  // True when zeroLoad/addLoad modified q0 since last update()
  bool committedAtCurrentV;  // True immediately after recvSelf at vCommit

  Node *theNodes[2];
  ID connectedExternalNodes;
  CrdTransf *theCoordTransf;

  // per-axis state determination is delegated to the shared, dimension-agnostic
  // kernel in LadrunoIMKHinge.h (ladrunoIMKSolveAxis / ladrunoIMKInitBlock),
  // also used by LadrunoIMKBeam2d.
};

#endif
