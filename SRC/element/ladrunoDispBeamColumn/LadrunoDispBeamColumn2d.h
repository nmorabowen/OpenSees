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
		     UniaxialMaterial *hingeMat = 0,         // Ladruno (ADR 32) Tier-2: embedded cohesive hinge M([[theta]])
		     UniaxialMaterial *hingeMatJ = 0,        // Ladruno: second end-J embedded cohesive hinge
		     int dualHingeMode = 0,                  // 0=legacy constant G=-1/L, 1=Hermite end-rotation modes G_I(xi), G_J(xi)
		     double hingeXiI = 0.0,
		     double hingeXiJ = 1.0,
		     double pmP0 = 0.0,
		     double pmAlpha = 0.5,
		     int spanLoadMode = 1);
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

    // Ladruno (WP-192/WP-196): OpenMP element loop opt-in
    bool ladrunoThreadSafeUpdate(void) const;
    bool ladrunoThreadSafeForm(void) const;
    bool ladrunoEnergyRayleighOnly(void) const { return theDamping == 0; }

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

    // Ladruno (ADR 32) Tier-2: embedded strong-discontinuity rotation-jump hinge(s).
    int              hingeOn;           // 1 if an embedded cohesive hinge is active
    int              dualHingeMode;     // 0 = legacy constant G=-1/L, 1 = Hermite G_I(xi), G_J(xi)
    double           hingeXiI;          // normalized hinge location for end I (default 0.0)
    double           hingeXiJ;          // normalized hinge location for end J (default 1.0)
    UniaxialMaterial *theHinge;         // end-I (or single) discrete cohesive M([[theta]]) law
    UniaxialMaterial *theHingeJ;        // end-J discrete cohesive M([[theta]]) law (optional)
    double           hingeJump;         // trial rotation jump alpha_I (converged in update())
    double           hingeJumpCommit;   // committed alpha_I
    double           hingeJumpJ;        // trial rotation jump alpha_J
    double           hingeJumpJCommit;  // committed alpha_J
    double           hingeKaa;          // cached guarded K_alphaalpha (single-hinge mode)
    double           hingeKaaCommit;    // committed guarded K_alphaalpha
    double           hingeKv[3];        // cached K_v-alpha 3-vector (single-hinge mode)
    double           hingeKvCommit[3];  // committed K_v-alpha 3-vector
    double           hingeHinv[2][2];       // cached 2x2 guarded H^-1 (dual-hinge mode)
    double           hingeHinvCommit[2][2]; // committed 2x2 guarded H^-1
    double           hingeKva[3][2];        // cached 3x2 K_v-alpha (dual-hinge mode)
    double           hingeKvaCommit[3][2];  // committed 3x2 K_v-alpha
    double           hingeMscale;       // running moment scale (~Mc) for inner-Newton tol
    double           hingeMscaleCommit; // committed running moment scale
    double           pmP0;              // P-M interaction squash load P0 (>0 enables scaling for general hinges)
    double           pmAlpha;           // P-M interaction linear/parabolic weight (default 0.5)
    bool             postRestoreSync;   // true right after recvSelf/revertToLastCommit until v or loads change
    bool             vRestoreValid;
    double           vRestore[3];

    // Tier-2 helper: inner Newton on the scalar or dual jump alpha
    int solveHingeJump(const Vector &v, double L);
    int solveDualHingeJump(const Vector &v, double L);
    int solveSpanLoadBubble(const Vector &v, double L);

	int numSections;
	SectionForceDeformation** theSections; // pointer to the ND material objects
	CrdTransf* crdTransf;        // pointer to coordinate transformation object
	BeamIntegration* beamInt;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    Vector Q;      // Applied nodal loads
    Vector q;      // Basic force
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system

    double rho;	   // Mass density per unit length
    int cMass;     // consistent mass flag

    Damping *theDamping;

    enum {maxNumSections = 20};

    // Ladruno (WP-200): Exact span-load internal strain + quartic bubble enhancement
    int    spanLoadMode;                    // 1 = enabled (default), 0 = disabled (-noSpanLoad)
    bool   hasSpanLoad;                     // true when Beam2dUniform/Partial/PointLoad active
    bool   spanLoadDirty;                   // true when zeroLoad/addLoad modified span loads since last update()
    double s0Fixed[maxNumSections][2];      // [i][0]=N_fixed(xi_i), [i][1]=Mz_fixed(xi_i)
    double q0Enh[3];                        // quadrature projection sum_i wt_i * b(xi_i)^T * s0Fixed(xi_i)
    double spanBetaP;                       // trial inelastic bubble axial amplitude
    double spanBetaPCommit;                 // committed inelastic bubble axial amplitude
    double spanBetaZ;                       // trial inelastic bubble curvature amplitude
    double spanBetaZCommit;                 // committed inelastic bubble curvature amplitude
    double spanKbbInv[2][2];                // cached guarded 2x2 bubble flexibility K_bb^-1
    double spanKbbInvCommit[2][2];          // committed guarded 2x2 bubble flexibility K_bb^-1
    double spanKvb[3][2];                   // cached 3x2 coupling K_vb
    double spanKvbCommit[3][2];             // committed 3x2 coupling K_vb

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif

