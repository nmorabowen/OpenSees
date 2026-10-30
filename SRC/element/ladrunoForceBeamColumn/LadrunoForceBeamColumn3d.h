/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
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

// Standalone Ladruno force-based (flexibility) 3D fiber frame element
// (Ladruno WP-200) combining exact force-based equilibrium state
// determination with Tier-1 per-IP/element/user characteristic length (lch),
// Tier-2 single/dual uncoupled (Mz, My) and coupled biaxial (Mz-My) cohesive
// fracture hinges with 4x4 Schur-complement static condensation, P-M
// interaction, and IMPLEX support.

#ifndef LadrunoForceBeamColumn3d_h
#define LadrunoForceBeamColumn3d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <BeamIntegration.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <vector>
#include <Damping.h>

class Response;
class ElementalLoad;
class UniaxialMaterial;
class NDMaterial;

class LadrunoForceBeamColumn3d : public Element
{
 public:
  LadrunoForceBeamColumn3d();
  LadrunoForceBeamColumn3d(int tag, int nodeI, int nodeJ,
                           int numSections, SectionForceDeformation **sec,
                           BeamIntegration &beamIntegr,
                           CrdTransf &coordTransf, double rho = 0.0,
                           int maxNumIters = 10, double tolerance = 1.0e-12,
                           int maxNumSubdivide = 4, double subdivideFactor = 10.0,
                           Damping *theDamping = 0,
                           int cMass = 0, int lchMode = 0, double userLch = 0.0,
                           UniaxialMaterial *hingeMatZ = 0,
                           UniaxialMaterial *hingeMatY = 0,
                           NDMaterial *hingeMatC = 0,
                           UniaxialMaterial *hingeMatZJ = 0,
                           UniaxialMaterial *hingeMatYJ = 0,
                           NDMaterial *hingeMatCJ = 0,
                           int dualHingeMode = 0,
                            double xiI = 0.0, double xiJ = 1.0,
                            double pmP0 = 0.0, double pmAlpha = 0.5,
                            bool eleImplex = false,
                            double GAvy = 0.0, double GAvz = 0.0,
                            double pmP0t = 0.0, double pmPbal = 0.0,
                            double pmMuBal = 0.0, bool axialRelease = false,
                            int spanLoadMode = 1);

  ~LadrunoForceBeamColumn3d();

  const char *getClassType(void) const { return "LadrunoForceBeamColumn3d"; }

  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);

  int getNumDOF(void);

  void setDomain(Domain *theDomain);
  int setDamping(Domain *theDomain, Damping *theDamping);
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  bool ladrunoThreadSafeUpdate(void) const override;
  bool ladrunoThreadSafeForm(void) const override;

  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getDampingForce(void);
  const Vector &getResistingForceIncInertia(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes = 0, int numModes = 0);

  friend OPS_Stream &operator<<(OPS_Stream &s, LadrunoForceBeamColumn3d &E);
  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);

  double getCharacteristicLength(void);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getKiSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  int getResponseSensitivity(int responseID, int gradNumber,
                             Information &eleInformation);

 protected:
  void setSectionPointers(int numSections, SectionForceDeformation **secPtrs);
  int getInitialFlexibility(Matrix &fe, bool includeHinges = true);
  int getInitialDeformations(Vector &v0);

 private:
  void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
  void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables(void);

  void computeReactions(double *p0);
  void computeSectionForces(Vector &sp, int isec);

  int evalBulkState(const Vector &v_tilde,
                    Vector &Se_out, Matrix &Kbulk_out, bool commitSubdivide = false);
  int solveHingeState(const Vector &v);

  ID connectedExternalNodes;

  BeamIntegration *beamIntegr;
  int numSections;
  SectionForceDeformation **sections;
  CrdTransf *crdTransf;

  double rho;
  int cMass;
  int maxIters;
  double tol;

  int initialFlag;

  Node *theNodes[2];

  Matrix kv;
  Vector Se;

  Matrix kvcommit;
  Vector Secommit;
  Matrix kvbulk;
  Matrix kvbulkCommit;

  Matrix *fs;
  Vector *vs;
  Vector *Ssr;
  Vector *vscommit;
  Matrix *fscommit;
  Vector *Ssrcommit;

  enum { maxNumEleLoads = 100 };
  enum { NDM = 3 };
  enum { NND = 6 };
  enum { NEGD = 12 };
  enum { NEBD = 6 };

  int numEleLoads;
  int sizeEleLoads;
  ElementalLoad **eleLoads;
  double *eleLoadFactors;
  Vector load;
  bool eleLoadDirty;

  Matrix *Ki;

  bool isTorsion;

  Damping *theDamping;

  double current_section_lch;
  int lchMode;
  double userLch;
  double GAvy;
  double GAvz;

  // Tier-2 embedded/series hinges in flexibility space
  int hingeOn;
  int dualHingeMode;
  double hingeXiI;
  double hingeXiJ;
  UniaxialMaterial *theHingeZ;
  UniaxialMaterial *theHingeY;
  NDMaterial *theHingeC;
  UniaxialMaterial *theHingeZJ;
  UniaxialMaterial *theHingeYJ;
  NDMaterial *theHingeCJ;
  double hingeJumpZ;
  double hingeJumpCommitZ;
  double hingeJumpY;
  double hingeJumpCommitY;
  double hingeJumpZJ;
  double hingeJumpCommitZJ;
  double hingeJumpYJ;
  double hingeJumpCommitYJ;
  double hingeMscaleZ;
  double hingeMscaleCommitZ;
  double hingeMscaleY;
  double hingeMscaleCommitY;
  double hingeMscaleZJ;
  double hingeMscaleCommitZJ;
  double hingeMscaleYJ;
  double hingeMscaleCommitYJ;
  double hingeDualKv[6][4];
  double hingeDualKvCommit[6][4];
  double hingeDualKav[4][6];
  double hingeDualKavCommit[4][6];
  double hingeDualKaaInv[4][4];
  double hingeDualKaaInvCommit[4][4];
  double pmP0;
  double pmAlpha;
  double pmP0t;
  double pmPbal;
  double pmMuBal;
  bool axialRelease;
  int spanLoadMode;
  bool useImplex;
  double eleImplexErr;
  double eleImplexErrCommit;
  bool postRestoreSync;
  double vRestore[6];
  bool vRestoreValid;
  Vector vCommitBulk;
  // Bulk deformation the TRIAL state (Se, kvbulk, vs, fs, Ssr) was solved at.
  // evalBulkState iterates from that state, not from the committed one.
  Vector vTrialBulk;
  // vTrialBulk at commitState. Not vCommitBulk: when the last Newton correction was
  // below DBL_EPSILON update() returned early, and the committed state was solved
  // at that earlier deformation, not at the committed one.
  Vector vTrialBulkCommit;
  // Basic deformation and element loads the trial state was last solved for.
  // update() skips a re-solve when neither changed (see solvedStateMatches).
  Vector vSolved;
  bool solvedValid;
  // The last update failed: the trial state is not a solution. update() must re-solve
  // (no early return) and commitState refuses it. Callers that ignore update()'s
  // return code (TransformationConstraintHandler::enforceSPs) would otherwise go on
  // with, and commit, the failed iterate.
  bool trialFailed;
  std::vector<ElementalLoad *> solvedLoads;
  std::vector<double> solvedLoadData;
  // The same record for the committed state, restored by revertToLastCommit.
  Vector vSolvedCommit;
  bool solvedCommitValid;
  std::vector<ElementalLoad *> solvedLoadsCommit;
  std::vector<double> solvedLoadDataCommit;
  bool solvedStateMatches(const Vector &v);
  void recordSolvedState(const Vector &v);

  static thread_local Matrix theMatrix;
  static thread_local Vector theVector;
  static thread_local double workArea[];

  enum { maxNumSections = 20 };

  int maxSubdivisions;
  double subdivideFactor;

  static thread_local Vector vsSubdivide[maxNumSections];
  static thread_local Vector SsrSubdivide[maxNumSections];
  static thread_local Matrix fsSubdivide[maxNumSections];

  int parameterID;
  const Vector &computedqdh(int gradNumber);
  const Matrix &computedfedh(int gradNumber);
  void computeReactionSensitivity(double *dp0dh, int gradNumber);
  void computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber);
};

#endif
