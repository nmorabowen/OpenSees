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

// LadrunoLST — 6-node linear-strain triangle (T6, Ladruno fork, ADR 70 P3).
// Quadratic shape functions ⇒ LINEAR strain, so it carries an inclined shear
// band (no mesh-line snapping, unlike the constant-strain LadrunoCST) and
// converges quadratically — the usable quadratic triangle for the plane family.
//
//   element('LadrunoLST', tag, n1,n2,n3, n4,n5,n6, matTag,
//           '-type', <PlaneStrain|PlaneStress>,    (default PlaneStrain)
//           '-geom', <linear|finite>,              (default linear)
//           '-thick', t, '-rho', r, '-body', b1, b2, '-pressure', p)
//
// Node order matches upstream SixNodeTri: corners n1,n2,n3 (CCW), midsides
// n4=(1-2), n5=(2-3), n6=(3-1). 3-point INTERIOR integration (barycentric
// (2/3,1/6),(1/6,2/3),(1/6,1/6), w=1/6 — the rule SixNodeTri ships).
// std -geom linear reduces to upstream SixNodeTri.
//
// NO bbar / F-bar lane — REFUTED at P3 (ADR 70 §9): constant element-mean
// dilatation (B-bar, and the F-bar scale (J0/J)^{1/2} sampled at the centroid)
// is RANK-DEFICIENT on the T6. The two quadratic CONFORMAL modes (Re/Im z²)
// have identically zero deviatoric strain and a linear dilatation with zero
// element mean, so the averaged operator assigns them ZERO energy: a free
// element shows 5 zero-energy modes (3 RBM + 2 spurious; stacked-B̄ rank 7 of
// the 9 required). Per fork policy no uncontrolled mechanism ships — the T6
// volumetric cure is ADR 70 P4 (nodal F-bar-Patch or P1-projected dilatation).
// '-formulation bbar' is refused at parse with this message.
//
// -geom finite (updated-Lagrangian large-strain, std only) drives a
//   FiniteStrainND2DMaterial (e.g. LogStrain2D) by the TOTAL 2x2 F via
//   setTrialF(F) through the shared LadrunoFiniteStrain2DKernel. PlaneStrain
//   only. corot (large-rotation/small-strain) is out of scope (ADR 70 §1 note).
//
// MASS (deliberate divergence from SixNodeTri): getMass uses HRZ lumping
// (diagonal of the consistent mass, ∫ρN_a² dV, rescaled to the exact total) —
// strictly positive. Upstream's plain N-lumping gives EXACTLY ZERO corner
// masses (negative on distorted elements) on the T6, which is unusable for
// explicit dynamics; same call the fork made for the H20 (ADR-72, HRZ-only).
// The reduce-to-SixNodeTri gate is static, so this does not affect it.
//
// See Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md

#ifndef LadrunoLST_h
#define LadrunoLST_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class LadrunoLST : public Element
{
  public:
    // std only — bbar/F-bar REFUTED on the T6 (rank-deficient, see header note);
    // the enum survives so the serialized formulation slot stays forward-usable.
    enum class Formulation { STD = 0 };
    // geometry axis (ADR 70). Values mirror LadrunoQuad/LadrunoBrick (2 = finite).
    enum class Geom { LINEAR = 0, FINITE = 2 };

    LadrunoLST(int tag, int nd1, int nd2, int nd3,
               int nd4, int nd5, int nd6,
               NDMaterial &m, const char *type, double t,
               Formulation form = Formulation::STD, Geom geom = Geom::LINEAR,
               double rho = 0.0, double b1 = 0.0, double b2 = 0.0,
               double pressure = 0.0,
               double b1bv = 0.0, double b2bv = 0.0);   // Ladruno (W2-E1): bulk-viscosity coeffs
    LadrunoLST();
    ~LadrunoLST();

    const char *getClassType(void) const { return "LadrunoLST"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    double getCharacteristicLength(void);   // crack-band lch = sqrt(2*area)

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

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
    static constexpr int numgp = 3;
    static constexpr int numnodes = 6;

    NDMaterial **theMaterial;          // 3 material points
    ID connectedExternalNodes;         // tags of the 6 nodes
    Node *theNodes[numnodes];

    static double matrixData[(2 * numnodes) * (2 * numnodes)];
    static Matrix K;                   // 12x12 element matrix (shared scratch)
    static Vector P;                   // 12 resisting force (shared scratch)
    Vector Q;                          // applied nodal loads
    double b[2];                       // body forces
    double appliedB[2];                // body force from load pattern
    int applyLoad;

    Vector pressureLoad;               // consistent pressure load
    double thickness;
    double pressure;
    double rho;
    Formulation formulation;
    Geom geom;                         // Ladruno (ADR 70): LINEAR / FINITE geometry axis
    double bulkVisc_b1;                // Ladruno (W2-E1): explicit bulk-viscosity coeffs (0=off)
    double bulkVisc_b2;
    int planeType;                     // 1 = PlaneStrain, 2 = PlaneStress

    static double shp[3][numnodes];    // shp[0/1] = dN/dx,dN/dy ; shp[2] = N
    static double pts[numgp][2];       // integration points (barycentric interior)
    static double wts[numgp];          // integration weights

    double shapeFunction(double xi, double eta);   // returns detJ, fills shp
    void formB(Matrix &B);                         // 3x12 strain-disp from shp

    // --- -geom finite (ADR 70 P3): updated-Lagrangian large strain ---------- //
    bool isFinite(void) const { return geom == Geom::FINITE; }
    int  updateFinite(void);              // per GP: F (+F-bar scale) → setTrialF(F)
    void formFinite(int tang_flag);       // fills static P (resid) + K (tangent) via kernel

    void setPressureLoadAtNodes(void);
    const char *typeString(void) const;

    Matrix *Ki;
};

#endif
