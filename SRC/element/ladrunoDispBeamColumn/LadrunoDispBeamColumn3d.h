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
class NDMaterial;

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
              UniaxialMaterial *hingeMatZ = 0,        // Ladruno (ADR 33) Tier-2: embedded strong-axis (Mz) cohesive hinge (end I / legacy)
              UniaxialMaterial *hingeMatY = 0,        // Ladruno (ADR 33 PR-3b): embedded weak-axis (My) cohesive hinge (end I)
              NDMaterial *hingeMatC = 0,              // Ladruno (ADR 34): coupled biaxial (Mz-My) cohesive hinge (end I)
              UniaxialMaterial *hingeMatZJ = 0,       // Ladruno (WP-200): end-J strong-axis (Mz) cohesive hinge
              UniaxialMaterial *hingeMatYJ = 0,       // Ladruno (WP-200): end-J weak-axis (My) cohesive hinge
              NDMaterial *hingeMatCJ = 0,             // Ladruno (WP-200): end-J coupled biaxial (Mz-My) cohesive hinge
              int dualMode = 0,                       // Ladruno (WP-200): 1 = Hermite end-rotation mode G_I(xi), G_J(xi)
              double xiI = 0.0,                       // Ladruno (WP-200): normalized location of hinge I in [0, 1]
              double xiJ = 1.0,                       // Ladruno (WP-200): normalized location of hinge J in [0, 1]
              double pmP0 = 0.0,
              double pmAlpha = 0.5,
              int spanLoadMode = 1,
              int pmBalMode = 0,
              double pmP0t = 0.0,
              double pmPbal = 0.0,
              double pmMuBal = 1.0,
              int axialRelease = 0);
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

    // Ladruno (WP-192/WP-196): OpenMP element loop opt-in
    bool ladrunoThreadSafeUpdate(void) const;
    bool ladrunoThreadSafeForm(void) const;
    bool ladrunoEnergyRayleighOnly(void) const { return theDamping == 0; }

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

    int              hingeOn;           // 1 if any embedded hinge is active
    UniaxialMaterial *theHingeZ;        // discrete cohesive M_z([[theta_z]]) law at end I / legacy
    double           hingeJumpZ;        // trial strong-axis rotation jump alpha_z
    double           hingeJumpCommitZ;  // committed alpha_z
    double           hingeKaaZ;         // cached guarded K_alphaalpha_z (Z-only rank-1)
    double           hingeKaaCommitZ;   // committed guarded K_alphaalpha_z
    double           hingeKvZ[6];       // cached K_v-alpha_z 6-vector
    double           hingeKvCommitZ[6]; // committed K_v-alpha_z 6-vector
    double           hingeMscaleZ;      // running moment scale (~Mc_z)
    double           hingeMscaleCommitZ;// committed running moment scale (~Mc_z)

    UniaxialMaterial *theHingeY;        // discrete cohesive M_y([[theta_y]]) law at end I
    double           hingeJumpY;        // trial weak-axis rotation jump alpha_y
    double           hingeJumpCommitY;  // committed alpha_y
    double           hingeKvY[6];       // cached K_v-alpha_y 6-vector
    double           hingeKvCommitY[6]; // committed K_v-alpha_y 6-vector
    double           hingeMscaleY;      // running moment scale (~Mc_y)
    double           hingeMscaleCommitY;// committed running moment scale (~Mc_y)
    double           hingeKaaInvZZ;     // cached eigenvalue-floored 2x2 inverse of K_aa
    double           hingeKaaInvZY;
    double           hingeKaaInvYY;
    double           hingeKaaInvCommitZZ;
    double           hingeKaaInvCommitZY;
    double           hingeKaaInvCommitYY;

    NDMaterial      *theHingeC;         // coupled biaxial cohesive M([[theta_z,theta_y]]) law at end I

    // Ladruno (WP-200): dual end hinges (ends I and J at normalized locations xiI, xiJ)
    int              dualHingeMode;     // 1 if Hermite end-rotation mode is active
    double           hingeXiI;          // normalized location of hinge I in [0, 1] (default 0.0)
    double           hingeXiJ;          // normalized location of hinge J in [0, 1] (default 1.0)
    UniaxialMaterial *theHingeZJ;       // end-J strong-axis (Mz) cohesive law
    UniaxialMaterial *theHingeYJ;       // end-J weak-axis (My) cohesive law
    NDMaterial       *theHingeCJ;       // end-J coupled biaxial (Mz-My) cohesive law
    double           hingeJumpZJ;       // trial end-J strong-axis rotation jump alpha_zJ
    double           hingeJumpCommitZJ; // committed alpha_zJ
    double           hingeJumpYJ;       // trial end-J weak-axis rotation jump alpha_yJ
    double           hingeJumpCommitYJ; // committed alpha_yJ
    double           hingeMscaleZJ;     // running moment scale (~Mc_zJ)
    double           hingeMscaleCommitZJ;
    double           hingeMscaleYJ;     // running moment scale (~Mc_yJ)
    double           hingeMscaleCommitYJ;
    double           hingeDualKv[6][4];         // cached 6x4 K_{v,alpha} for [azI, ayI, azJ, ayJ]
    double           hingeDualKvCommit[6][4];
    double           hingeDualKav[4][6];        // cached 4x6 K_{alpha,v} including P-M derivative
    double           hingeDualKavCommit[4][6];
    double           hingeDualKaaInv[4][4];     // cached 4x4 eigenvalue-floored inverse of K_{alpha,alpha}
    double           hingeDualKaaInvCommit[4][4];
    double           pmP0;                      // P-M interaction squash load P0 (>0 enables scaling for general hinges)
    double           pmAlpha;                   // P-M interaction linear/parabolic weight (default 0.5)
    int              pmBalMode;                 // 1 if 4-param RC balance-point P-M interaction is active
    double           pmP0t;                     // P-M tension capacity (>0)
    double           pmPbal;                    // P-M balance compression load (>0)
    double           pmMuBal;                   // P-M moment ratio Mbal/M0 at P = Pbal
    int              axialRelease;              // 1 = decouple fiber axial strain from basic chord elongation v(0)
    bool             postRestoreSync;           // true right after recvSelf/revertToLastCommit until v or loads change
    bool             vRestoreValid;
    double           vRestore[6];

    // Tier-2 helper: inner Newton on the scalar jump alpha_z (single-hinge Z-only).
    int solveHingeJump(const Vector &v, double L);

    // Tier-2 PR-3b helper: unified inner Newton on the biaxial jump [alpha_z, alpha_y] (single-hinge).
    int solveHingeJumpBiaxial(const Vector &v, double L);

    // WP-200 helper: inner Newton on the 4-DOF dual-hinge jumps [alpha_zI, alpha_yI, alpha_zJ, alpha_yJ].
    int solveDualHingeJump3d(const Vector &v, double L);

    // WP-200 helper: inner Newton on the 2-DOF quartic span-load bubble amplitudes [beta_z, beta_y].
    int solveSpanLoadBubble3d(const Vector &v, double L, const double *xi, const double *wt);

	int numSections;
	SectionForceDeformation** theSections; // pointer to the ND material objects
	CrdTransf* crdTransf;        // pointer to coordinate transformation object 
	BeamIntegration* beamInt;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    Vector Q;      // Applied nodal loads
    Vector q;      // Basic force
    double q0[5];  // Fixed end forces in basic system (no torsion)
    double p0[5];  // Reactions in basic system (no torsion)

    double rho;    // Mass density per unit length
    int cMass;     // consistent mass flag

    Damping *theDamping;

	int parameterID;

    enum {maxNumSections = 20};

    // Ladruno (WP-200): span-load internal strain + 3-DOF (P, Mz, My) quartic bubble enhancement
    int    spanLoadMode;                    // 1 = enabled (default), 0 = legacy nodal-q0-only
    int    hasSpanLoad;                     // 1 if active span load contributed to s0Fixed
    bool   spanLoadDirty;                   // true if addLoad/zeroLoad modified loads since last update()
    double s0Fixed[maxNumSections][3];      // [N_fixed(xi_i), Mz_fixed(xi_i), My_fixed(xi_i)]
    double q0Enh[5];                        // quadrature integral of B^T * s0Fixed
    double spanBetaP;                       // trial quartic bubble axial amplitude beta_p
    double spanBetaPCommit;                 // committed beta_p
    double spanBetaZ;                       // trial quartic bubble amplitude beta_z
    double spanBetaZCommit;                 // committed beta_z
    double spanBetaY;                       // trial quartic bubble amplitude beta_y
    double spanBetaYCommit;                 // committed beta_y
    double spanKbbInv[3][3];                // cached eigenvalue-floored 3x3 inverse of K_bb
    double spanKbbInvCommit[3][3];
    double spanKvb[6][3];                   // cached 6x3 K_{v,beta} coupling matrix
    double spanKvbCommit[6][3];
    double eps0Axial[maxNumSections];       // trial neutral-axis strain for -axialRelease
    double eps0AxialCommit[maxNumSections]; // committed neutral-axis strain for -axialRelease
};

#endif

