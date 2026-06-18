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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for LadrunoDispBeamColumn3d.
// The element displacement field gives rise to constant axial strain,
// linear curvature, and constant twist angle.

#ifndef LadrunoDispBeamColumn3d_h
#define LadrunoDispBeamColumn3d_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Damping.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class BeamIntegration;
class Response;
class UniaxialMaterial;

class LadrunoDispBeamColumn3d : public Element
{
  public:
    LadrunoDispBeamColumn3d(int tag, int nd1, int nd2,
		     int numSections, SectionForceDeformation **s,
		     BeamIntegration &bi, CrdTransf &coordTransf,
             double rho = 0.0, int cMass = 0,
             Damping *theDamping = 0,
             int lchMode = 0, double userLch = 0.0,  // Ladruno (ADR 32): regularization length mode
             int nlGeom = 0,                         // Ladruno (ADR 32): 0=linear basic strain, 1=½(θy²+θz²) bowing (NL)
             UniaxialMaterial *hingeMatZ = 0,        // Ladruno (ADR 33) Tier-2: embedded strong-axis (Mz) cohesive hinge
             UniaxialMaterial *hingeMatY = 0);       // Ladruno (ADR 33 PR-3b): embedded weak-axis (My) cohesive hinge (biaxial)
    LadrunoDispBeamColumn3d();
    ~LadrunoDispBeamColumn3d();

    const char *getClassType(void) const {return "LadrunoDispBeamColumn3d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    int setDamping(Domain *theDomain, Damping *theDamping);

    // Ladruno (ADR 32): per-IP characteristic length for crack-band regularization
    double getCharacteristicLength(void);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getDampingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int            updateParameter(int parameterID, Information &info);
    int            activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getInitialStiffSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    void getBasicStiff(Matrix &kb, int initial = 0);

    // Ladruno (ADR 32): see LadrunoDispBeamColumn2d.h. lchMode 0=ip/1=element/2=user;
    // current_section_lch is transient (set per-IP in update()); nlGeom adds the
    // ½(θy²+θz²) biaxial bowing strain (DispBeamColumnNL3d) when set.
    double current_section_lch;
    int    lchMode;
    double userLch;
    int    nlGeom;

    // Ladruno (ADR 33) Tier-2 PR-3a: embedded strong-axis (Mz) rotation-jump hinge — the
    // literal 2D scalar hinge on the Mz row of the 6-DOF 3D basic system. One scalar jump
    // alpha_z carried by a discrete cohesive law theHingeZ (any UniaxialMaterial). The bulk
    // section sees kappa_z_bulk = (1/L)((6xi-4)v1+(6xi-2)v2) + Gbar_z*alpha_z, Gbar_z=-1/L,
    // so it unloads as the jump grows (no double-count); the cohesive M([[theta]]) carries Gf.
    // alpha_z is converged by an inner Newton in update() and STATICALLY CONDENSED to the
    // 6-DOF basic system BEFORE CorotCrdTransf3d (the PINNED INVARIANT): a guarded rank-1
    // update kb -= hingeKvZ hingeKvZ^T / hingeKaaZ. hingeKvZ is a 6-vector (incl. the
    // cross-tangent rows ks(MY,MZ), ks(T,MZ), ks(P,MZ)) so the condensed off-diagonals are
    // right. ALL gated on hingeOn so the no-hinge path is byte-identical. v1: strong axis
    // only (alpha_y / the coupled 2x2 / torsion are PR-3b+), linear basic strain (XOR -nl),
    // monotonic full-Newton / DisplacementControl. See ADR 33 (built on the 2D ADR 32 PR-2a).
    int              hingeOn;           // 1 if the embedded Mz hinge is active
    UniaxialMaterial *theHingeZ;        // discrete cohesive M_z([[theta_z]]) law (element owns a copy)
    double           hingeJumpZ;        // trial strong-axis rotation jump alpha_z (converged in update())
    double           hingeJumpCommitZ;  // committed alpha_z (inner-Newton warm start + serialization)
    double           hingeKaaZ;         // cached guarded K_alphaalpha_z at the converged state (Z-only rank-1)
    double           hingeKvZ[6];       // cached K_v-alpha_z 6-vector at the converged state
    double           hingeMscaleZ;      // running moment scale (~Mc_z) for a stable inner-Newton tol

    // Ladruno (ADR 33) Tier-2 PR-3b: the WEAK-axis (My) jump alpha_y, making the hinge biaxial.
    // Present (theHingeY != 0) iff -hingeY was given (which requires -hinge). The two bending
    // jumps alpha = [alpha_z, alpha_y] are solved by ONE unified inner Newton (a single
    // setTrialSectionDeformation per IP sets BOTH bulk curvatures kappa_*_bulk = B_*v - alpha_*/L),
    // and condensed by a TRUE coupled 2x2: K_aa carries the off-diagonal (1/L)sum wt*ks(MZ,MY) so
    // the staggered-activation tangent is exact when the section couples bending axes (ks(MZ,MY)!=0,
    // generic for any axially-loaded biaxially-bent section). K_aa^{-1} is an EIGENVALUE-FLOORED
    // inverse (bounds every mode, sign per eigenvalue preserved) — not two rank-1 updates, not a
    // det-floored adjugate. Block-diagonal cohesive (one scalar law per axis); biaxial coupling
    // flows through the bulk ks. See ADR 33 PR-3b.
    UniaxialMaterial *theHingeY;        // discrete cohesive M_y([[theta_y]]) law (element owns a copy)
    double           hingeJumpY;        // trial weak-axis rotation jump alpha_y (converged in update())
    double           hingeJumpCommitY;  // committed alpha_y (inner-Newton warm start + serialization)
    double           hingeKvY[6];       // cached K_v-alpha_y 6-vector at the converged state
    double           hingeMscaleY;      // running moment scale (~Mc_y) for a stable inner-Newton tol
    double           hingeKaaInvZZ;     // cached eigenvalue-floored 2x2 inverse of K_aa (symmetric):
    double           hingeKaaInvZY;     //   [[hingeKaaInvZZ, hingeKaaInvZY],
    double           hingeKaaInvYY;     //    [hingeKaaInvZY, hingeKaaInvYY]] (biaxial rank-2 condensation)

    // Tier-2 helper: inner Newton on the scalar jump alpha_z; sets the 6-DOF section trial
    // deformations with the -alpha_z/L Mz-curvature offset, caches hingeKaaZ/hingeKvZ. (Z-only.)
    int solveHingeJump(const Vector &v, double L);

    // Tier-2 PR-3b helper: unified inner Newton on the biaxial jump [alpha_z, alpha_y]; sets both
    // bulk curvatures in ONE setTrialSectionDeformation/IP, caches hingeKvZ/hingeKvY + the
    // eigenvalue-floored 2x2 inverse hingeKaaInv* for the rank-2 condensation.
    int solveHingeJumpBiaxial(const Vector &v, double L);

	int numSections;
	SectionForceDeformation** theSections; // pointer to the ND material objects
	CrdTransf* crdTransf;        // pointer to coordinate transformation object 
	BeamIntegration* beamInt;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector

    Vector Q;      // Applied nodal loads
    Vector q;      // Basic force
    double q0[5];  // Fixed end forces in basic system (no torsion)
    double p0[5];  // Reactions in basic system (no torsion)

    double rho;    // Mass density per unit length
    int cMass;     // consistent mass flag

    Damping *theDamping;

	int parameterID;

    enum {maxNumSections = 20};

    static double workArea[];
};

#endif

