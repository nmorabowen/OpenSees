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

// LadrunoCST — 3-node constant-strain triangle (Ladruno fork). The thin 2D
// sibling of LadrunoQuad: std formulation only (a 1-point triangle is
// rank-sufficient — no hourglass; bbar/ssp/eas have nothing to average). Ships
// as the baseline / triangular-mesh fallback / future E-FEM carrier.
//
//   element('LadrunoCST', tag, n1,n2,n3, matTag,
//           '-type', <PlaneStrain|PlaneStress>, '-geom', <linear|finite>,
//           '-thick', t, '-rho', r, '-body', b1, b2, '-pressure', p)
//
// v1 (ADR 25, Phase 1): small-strain linear; reduces to upstream Tri31.
//
// -geom finite (ADR 70 P2, updated-Lagrangian large-strain) drives a
//   FiniteStrainND2DMaterial (e.g. LogStrain2D) by the TOTAL 2x2 deformation
//   gradient via setTrialF(F) and assembles the spatial internal force + full
//   consistent tangent (a_ijkl = c_ijkl − σ_il δ_jk) through the shared
//   LadrunoFiniteStrain2DKernel. PlaneStrain only. NOTE (the honest baseline,
//   ADR 70 §3.2): element-level F-bar is a NO-OP on the constant-strain T3
//   (J is element-constant ⇒ J0 ≡ J ⇒ F̄ = F — nothing to average), so
//   CST-finite locks volumetrically near incompressibility. The usable
//   finite-strain triangle is LadrunoLST (T6, ADR 70 P3); the T3 cure is
//   nodal F-bar-Patch (ADR 70 P4, not built).
//
// See Ladruno_implementation/25_ladruno_plane_elements_adr.md and
// Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md

#ifndef LadrunoCST_h
#define LadrunoCST_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class LadrunoCST : public Element
{
  public:
    // geometry axis (ADR 70 P2). LINEAR = small strain; FINITE = updated-
    // Lagrangian large strain (setTrialF path). Values mirror LadrunoQuad/Brick.
    enum class Geom { LINEAR = 0, FINITE = 2 };

    LadrunoCST(int tag, int nd1, int nd2, int nd3,
               NDMaterial &m, const char *type, double t,
               Geom geom = Geom::LINEAR,                 // Ladruno (ADR 70)
               double rho = 0.0, double b1 = 0.0, double b2 = 0.0,
               double pressure = 0.0,
               double b1bv = 0.0, double b2bv = 0.0);   // Ladruno (W2-E1): bulk-viscosity coeffs
    LadrunoCST();
    ~LadrunoCST();

    const char *getClassType(void) const { return "LadrunoCST"; }

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
    NDMaterial **theMaterial;
    ID connectedExternalNodes;
    Node *theNodes[3];

    static double matrixData[36];
    static Matrix K;             // 6x6
    static Vector P;             // 6
    Vector Q;
    double b[2];
    double appliedB[2];
    int applyLoad;

    Vector pressureLoad;
    double thickness;
    double pressure;
    double rho;
    Geom geom;                   // Ladruno (ADR 70): LINEAR / FINITE geometry axis
    double bulkVisc_b1;          // Ladruno (W2-E1): explicit bulk-viscosity coeffs (0=off)
    double bulkVisc_b2;
    int planeType;               // 1 = PlaneStrain, 2 = PlaneStress

    static double shp[3][3];
    static double pts[1][2];
    static double wts[1];

    double shapeFunction(double xi, double eta);

    // --- -geom finite (ADR 70 P2): updated-Lagrangian large strain ---------- //
    bool isFinite(void) const { return geom == Geom::FINITE; }
    int  updateFinite(void);              // centroid GP: F → setTrialF(F)
    void formFinite(int tang_flag);       // fills static P (resid) + K (tangent) via kernel

    void setPressureLoadAtNodes(void);
    const char *typeString(void) const;

    Matrix *Ki;

    static constexpr int numgp = 1;
    static constexpr int numnodes = 3;
};

#endif
