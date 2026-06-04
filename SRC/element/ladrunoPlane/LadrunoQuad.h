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

// LadrunoQuad — unified 4-node bilinear plane (plane-stress/strain) continuum
// element for the Ladruno fork. The 2D sibling of LadrunoBrick: the formulation
// is a parameter, not a class.
//
//   element('LadrunoQuad', tag, n1,n2,n3,n4, matTag,
//           '-formulation', <std|bbar|ssp|eas>,   (default std)
//           '-type', <PlaneStrain|PlaneStress>,    (default PlaneStrain)
//           '-thick', t, '-rho', r, '-body', b1, b2, '-pressure', p)
//
// v1 (ADR 25, Phase 1): small-strain (geometrically linear) kinematics, two
// formulations:
//   std  — full 2x2 Gauss displacement      (reduces to upstream FourNodeQuad)
//   bbar — mean-dilatation B-bar (2D factor 1/2; PlaneStrain only) — cures
//          volumetric locking
// ssp / eas are RESERVED (parser errors with a "not yet implemented" hint) so
// the API is stable; the geometry layer (-geom corot/finite) drops in later via
// the SolidTransformation seam (ADR 25 P4/P5).
//
// See Ladruno_implementation/25_ladruno_plane_elements_adr.md

#ifndef LadrunoQuad_h
#define LadrunoQuad_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class LadrunoQuad : public Element
{
  public:
    enum class Formulation { STD = 0, BBAR = 1, SSP = 2, EAS = 3 };

    LadrunoQuad(int tag, int nd1, int nd2, int nd3, int nd4,
                NDMaterial &m, const char *type, double t,
                Formulation form,
                double rho = 0.0, double b1 = 0.0, double b2 = 0.0,
                double pressure = 0.0);
    LadrunoQuad();
    ~LadrunoQuad();

    const char *getClassType(void) const { return "LadrunoQuad"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // state
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // stiffness / mass / residual
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // parallel / IO
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:

  private:
    NDMaterial **theMaterial;          // 4 material points
    ID connectedExternalNodes;         // tags of the 4 nodes
    Node *theNodes[4];

    static double matrixData[64];
    static Matrix K;                   // 8x8 element matrix (shared scratch)
    static Vector P;                   // 8 resisting force (shared scratch)
    Vector Q;                          // applied nodal loads
    double b[2];                       // body forces
    double appliedB[2];                // body force from load pattern
    int applyLoad;

    Vector pressureLoad;               // consistent pressure load
    double thickness;
    double pressure;
    double rho;
    Formulation formulation;
    int planeType;                     // 1 = PlaneStrain, 2 = PlaneStress

    static double shp[3][4];           // shp[0/1] = dN/dx,dN/dy ; shp[2] = N
    static double shpBar[2][4];        // element-mean gradients (B-bar)
    static double pts[4][2];           // gauss points
    static double wts[4];              // gauss weights

    double shapeFunction(double xi, double eta);  // returns detJ, fills shp
    void computeShapeBar(void);                    // fills shpBar (B-bar)
    void formB(Matrix &B);                         // 3x8 strain-disp from shp/shpBar
    void setPressureLoadAtNodes(void);
    const char *typeString(void) const;

    Matrix *Ki;
};

#endif
