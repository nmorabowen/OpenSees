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

// Standalone Ladruno force-based (flexibility) 2D fiber frame element
// (Ladruno WP-200) combining exact force-based equilibrium state
// determination with Tier-1 per-IP/element/user characteristic length (lch),
// Tier-2 single/dual cohesive fracture hinges with 2x2 Schur-complement
// static condensation, P-M interaction, and IMPLEX support.

#ifndef LadrunoForceBeamColumn2d_h
#define LadrunoForceBeamColumn2d_h

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

class LadrunoForceBeamColumn2d : public Element
{
 public:
  LadrunoForceBeamColumn2d(int tag, int nodeI, int nodeJ,
                           int numSections, SectionForceDeformation **sec,
                           BeamIntegration & beamIntegr,
                           CrdTransf &coordTransf, double massDensPerUnitLength = 0.0,
                           int maxNumIters = 10, double tolerance = 1.0e-12,
                           int cMass = 0,
                           Damping *theDamping = 0,
                           int lchMode = 0, double userLch = 0.0,
                           UniaxialMaterial *hingeMat = 0,
                           UniaxialMaterial *hingeMatJ = 0,
                           int dualHingeMode = 0,
                           double xiI = 0.0, double xiJ = 1.0,
                           double pmP0 = 0.0, double pmAlpha = 0.5,
                           int maxSub = 4, double subFac = 10.0,
                           bool eleImplex = false,
                           double GAvy = 0.0,
                           double pmP0t = 0.0, double pmPbal = 0.0, double pmMuBal = 1.0,
                           bool axialRelease = false, int spanLoadMode = 1);

  LadrunoForceBeamColumn2d();
  ~LadrunoForceBeamColumn2d();

  const char *getClassType(void) const { return "LadrunoForceBeamColumn2d"; }

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

  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getDampingForce(void);
  const Vector &getResistingForceIncInertia(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes = 0, int numModes = 0);

  friend OPS_Stream &operator<<(OPS_Stream &s, LadrunoForceBeamColumn2d &E);
  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  int getResponseSensitivity(int responseID, int gradNumber, Information &eleInformation);

  double getCharacteristicLength(void);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getInitialStiffSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);

  bool ladrunoThreadSafeUpdate(void) const override;
  bool ladrunoThreadSafeForm(void) const override;

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
  void computeReactionSensitivity(double *dp0dh, int gradNumber);
  void computeSectionForces(Vector &sp, int isec);
  void computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber);

  const Vector &computedqdh(int gradNumber);
  const Matrix &computedfedh(int gradNumber);

  int evalBulkState(const Vector &v_tilde,
                    Vector &Se_out, Matrix &Kbulk_out, bool commitSubdivide = false);
  int solveHingeState(const Vector &v);

  enum { maxNumSections = 20 };
  enum { NDM = 2 };
  enum { NND = 3 };
  enum { NEGD = 6 };
  enum { NEBD = 3 };

  int numSections;
  SectionForceDeformation **sections;
  CrdTransf *crdTransf;
  BeamIntegration *beamIntegr;

  ID connectedExternalNodes;
  Node *theNodes[2];

  Vector load;
  double rho;
  int cMass;
  int maxIters;
  double tol;

  int initialFlag;

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
  ElementalLoad *eleLoads[maxNumEleLoads];
  double eleLoadFactors[maxNumEleLoads];
  int numEleLoads;
  bool eleLoadDirty;

  Matrix *Ki;
  int parameterID;
  Damping *theDamping;

  int maxSubdivisions;
  double subdivideFactor;

  double current_section_lch;
  int lchMode;
  double userLch;
  double GAvy;

  // Tier-2 embedded/series hinges in flexibility space
  int hingeOn;
  int dualHingeMode;
  double hingeXiI;
  double hingeXiJ;
  UniaxialMaterial *theHinge;
  UniaxialMaterial *theHingeJ;
  double hingeJump;
  double hingeJumpCommit;
  double hingeJumpJ;
  double hingeJumpJCommit;
  double hingeHinv[2][2];
  double hingeHinvCommit[2][2];
  double hingeKva[3][2];
  double hingeKvaCommit[3][2];
  double hingeKav[2][3];
  double hingeKavCommit[2][3];
  double hingeMscale;
  double hingeMscaleCommit;
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
  double vRestore[3];
  bool vRestoreValid;
  Vector vCommitBulk;
  // Bulk deformation the TRIAL state (Se, kvbulk, vs, fs, Ssr) was solved at.
  // evalBulkState iterates from that state, not from the committed one.
  Vector vTrialBulk;
  // Basic deformation and element loads the trial state was last solved for.
  // update() skips a re-solve when neither changed (see solvedStateMatches).
  Vector vSolved;
  bool solvedValid;
  std::vector<ElementalLoad *> solvedLoads;
  std::vector<double> solvedLoadData;
  bool solvedStateMatches(const Vector &v);
  void recordSolvedState(const Vector &v);

  static thread_local Matrix theMatrix;
  static thread_local Vector theVector;
  static thread_local double workArea[];

  static thread_local Vector vsSubdivide[maxNumSections];
  static thread_local Vector SsrSubdivide[maxNumSections];
  static thread_local Matrix fsSubdivide[maxNumSections];
};

#endif
