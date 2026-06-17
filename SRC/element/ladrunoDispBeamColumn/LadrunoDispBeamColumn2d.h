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
// Description: This file contains the class definition for LadrunoDispBeamColumn2d.
// The element displacement field gives rise to constant axial strain and
// linear curvature.

#ifndef LadrunoDispBeamColumn2d_h
#define LadrunoDispBeamColumn2d_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>
#include <Damping.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class Response;
class UniaxialMaterial;

class LadrunoDispBeamColumn2d : public Element
{
  public:
    LadrunoDispBeamColumn2d(int tag, int nd1, int nd2,
		     int numSections, SectionForceDeformation **s,
		     BeamIntegration &bi, CrdTransf &coordTransf,
		     double rho = 0.0, int cMass = 0,
		     Damping *theDamping = 0,
		     int lchMode = 0, double userLch = 0.0,  // Ladruno (ADR 32): regularization length mode
		     int nlGeom = 0,                         // Ladruno (ADR 32): 0=linear basic strain, 1=½θ² bowing (NL)
		     UniaxialMaterial *hingeMat = 0);        // Ladruno (ADR 32) Tier-2: embedded cohesive hinge M([[theta]])
    LadrunoDispBeamColumn2d();
    ~LadrunoDispBeamColumn2d();

    const char *getClassType(void) const {return "LadrunoDispBeamColumn2d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    int setDamping(Domain *theDomain, Damping *theDamping);

    // Ladruno (ADR 32): report the localizing IP's tributary length so crack-band /
    // auto-regularizing materials (ASDConcrete1D -autoRegularization, etc.) regularize
    // over the correct band instead of the whole element (cf. ForceBeamColumn2d).
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
    int getResponseSensitivity(int responseID, int gradNumber,
			       Information &eleInformation);

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

    // Ladruno (ADR 32) Tier-1 regularization-length channel.
    // lchMode: 0 = ip (per-IP tributary wt[i]*L, the safe default),
    //          1 = element (full L, A/B-debug, re-enables multi-IP energy double-counting),
    //          2 = user (a fixed crack-band width = userLch).
    // current_section_lch is transient: set per integration point inside update()
    // immediately before setTrialSectionDeformation(), then read by the active
    // material through ops_TheActiveElement->getCharacteristicLength().
    double current_section_lch;
    int    lchMode;
    double userLch;

    // Ladruno (ADR 32) Stage-1: nlGeom = 0 uses the stock linear basic strain
    // (constant axial + linear curvature); nlGeom = 1 adds the ½θ² bowing term to
    // the axial section strain (DispBeamColumnNL2d formulation), giving the
    // displacement-based P-δ / bowing coupling. Both still rely on the geomTransf
    // for rigid-body rotation; -nl improves the in-element large-deflection response.
    int    nlGeom;

    // Ladruno (ADR 32) Tier-2: embedded strong-discontinuity rotation-jump hinge.
    // hingeOn=1 activates a single scalar rotation jump `alpha` carried by a discrete
    // cohesive law theHinge (any UniaxialMaterial; LadrunoCohesiveHinge is the default).
    // The bulk section at every IP sees the BOUNDED enhanced curvature
    //   kappa_bulk = B*v + Gbar*alpha,   Gbar = -1/L   (constant, orthogonal: integral Gbar dx = -1)
    // so as alpha grows the bulk UNLOADS (no double count); the cohesive M([[theta]])
    // carries all post-peak fracture energy. alpha is converged by an inner Newton in
    // update() and statically condensed to the 3-DOF basic system BEFORE crdTransf
    // (the PINNED INVARIANT): K_basic = K_vv - K_valpha * (1/K_aa) * K_valpha^T, with a
    // GUARDED reciprocal (K_aa is sign-discontinuous at activation and indefinite on the
    // softening branch). All of this is GATED on hingeOn so the no-hinge path is
    // byte-identical to Stage-1. v1: linear basic strain only (mutually exclusive with -nl),
    // single hinge, monotonic full-Newton / DisplacementControl. See ADR 32 (2026-06-17
    // adversarial-review CORRECTIONS) for the derivation and the deferred hardening.
    int              hingeOn;           // 1 if an embedded cohesive hinge is active
    UniaxialMaterial *theHinge;         // the discrete cohesive M([[theta]]) law (element owns a copy)
    double           hingeJump;         // trial rotation jump alpha (converged in update())
    double           hingeJumpCommit;   // committed alpha (inner-Newton warm start + serialization)
    double           hingeKaa;          // cached guarded K_alphaalpha at the converged state
    double           hingeKv[3];        // cached K_v-alpha 3-vector at the converged state
    double           hingeMscale;       // running moment scale (~Mc) for a stable inner-Newton
                                        // convergence tol that does not collapse when the
                                        // fully-broken hinge carries M -> 0 (transient)

    // Tier-2 helper: inner Newton on the scalar jump alpha; sets the section trial
    // deformations with the -alpha/L bulk-curvature offset, caches hingeKaa/hingeKv.
    int solveHingeJump(const Vector &v, double L);

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
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system

    double rho;	   // Mass density per unit length
    int cMass;     // consistent mass flag

    Damping *theDamping;

    enum {maxNumSections = 20};

    static double workArea[];

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif

