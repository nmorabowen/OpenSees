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

#ifndef ElasticTimoshenkoBeam3d_h
#define ElasticTimoshenkoBeam3d_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/13
// Revision: A
//
// Purpose: This file contains the class definition for ElasticTimoshenkoBeam3d.
// ElasticTimoshenkoBeam3d is a frame member that takes shearing deformations
// and rotational inertia effects into account.

#include <Element.h>
#include <Matrix.h>

class CrdTransf;
class SectionForceDeformation;

class ElasticTimoshenkoBeam3d : public Element
{
public:
    // constructor
    ElasticTimoshenkoBeam3d(int tag, int Nd1, int Nd2, double E, double G,
			    double A, double Jx, double Iy, double Iz, double Avy, double Avz,
			    CrdTransf &theTransf, double rho = 0.0, int cMass = 0, int geomNonlinear = 0);
    ElasticTimoshenkoBeam3d(int tag, int Nd1, int Nd2, 
			    SectionForceDeformation &section,
			    CrdTransf &theTransf, double rho = 0.0, int cMass = 0, int geomNonlinear = 0);  
    ElasticTimoshenkoBeam3d();
    
    // destructor
    ~ElasticTimoshenkoBeam3d();
    
    // method to get class type
    const char *getClassType() const {return "ElasticTimoshenkoBeam3d";};
    
    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
    
    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();
    
    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();
    
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag = 0);
    
    // public methods for element recorder
    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
    
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);
    
protected:
    
private:
    // private methods
    void setUp();
    // Ladruno WP-206: basic-system path through the CrdTransf
    void formBasic(void);           // ub -> qb (+ q0, + N*kgb*ub if geomNonlinear)
    void basicToLocal(Vector &qlOut) const;
    void loadsToBasic(void);        // ql0 (local FEF) -> q0b, p0b
    
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;  // contains the tags of the end nodes
    Node *theNodes[2];          // array of nodes
    CrdTransf *theCoordTransf;  // coordinate transformation
    
    // parameters
    double E;     // elastic modulus
    double G;     // shear modulus
    double A;     // cross section (axial) area
    double Jx;    // torsional constant
    double Iy;    // moment of inertia about local y axis
    double Iz;    // moment of inertia about local z axis
    double Avy;   // shear area along local y axis
    double Avz;   // shear area along local z axis
    double rho;   // mass per unit length
    int cMass;    // consistent mass flag
    int nlGeo;    // nonlinear geometry flag
    double phiY;  // ratio of bending to shear stiffness about local y axis
    double phiZ;  // ratio of bending to shear stiffness about local z axis
    double L;     // element length
    
    // state variables
    Vector ul;     // displacements in local system
    Vector ql;     // forces in local system
    Vector ql0;    // fixed end forces due to loads in local system
    Matrix kl;     // stiffness matrix in local system
    Matrix klgeo;  // geometric stiffness matrix in local system
    
    // constant variables
    Matrix Tgl;  // transformation matrix from global to local system
    Matrix Ki;   // initial stiffness matrix in global system
    Matrix M;    // mass matrix in global system
    
    static Matrix theMatrix;  // a class wide Matrix
    static Vector theVector;  // a class wide Vector
    Vector theLoad;

    // Ladruno WP-206: the stock element built its own rotation-only Tgl from
    // the transformation's axes, so a PDelta transf (P-Delta), a Corotational
    // transf and -jntOffset rigid arms were silently ignored. Unless the legacy
    // combination (-geomNonlinear with a Linear transf) is requested, the
    // element now works in the 6-dof basic system and lets the CrdTransf map
    // it to global (offsets, P-Delta, corotational), like ElasticBeam3d.
    bool useBasic;      // true: basic-system path; false: legacy local path
    int geomFlag;       // the user's -geomNonlinear flag (element P-delta)
    Matrix kb;          // 6x6 basic stiffness [N, Mz1, Mz2, My1, My2, T]
    Matrix kgb;         // 6x6 basic P-delta geometric stiffness (per unit N)
    Vector qb;          // basic forces
    Vector q0b;         // basic fixed-end forces from element loads
    Vector p0b;         // 5 reactions (CrdTransf convention) from element loads
};

#endif
