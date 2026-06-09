/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
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

// Authors: Nicolas Mora Bowen, Patricio Palacios, Jose Abell, Guppi (Ladruño)
// Created: 03/2026
//
// Description: 6-node quadratic Bézier triangular element for 2D analysis.
//
// Based on:
//   Kadapa, C. "Novel quadratic Bézier triangular and tetrahedral elements
//   using existing mesh generators: Applications to linear nearly
//   incompressible elastostatics and implicit and explicit elastodynamics."
//   International Journal for Numerical Methods in Engineering,
//   2019; 117(5):543-573. doi:10.1002/nme.5967
//
// Features:
//   - Quadratic Bernstein polynomial shape functions (all nonnegative)
//   - Pure displacement and B-bar formulations
//   - Consistent and lumped mass matrices (all-positive lumped mass)
//   - Plane stress and plane strain
//   - Lagrange-to-Bézier control point mapping
//
// v1 limitation: validated for STRAIGHT-SIDED meshes only (mid-edge nodes at
// edge midpoints ⇒ P = X). Curved elements emit a warning in setDomain; the
// Eq. 14 Dirichlet mapping + higher quadrature they need is a follow-up.
//
// Element command:
//   element BezierTri6 $tag $nd1 $nd2 $nd3 $nd4 $nd5 $nd6
//                      $thick $type $matTag
//                      <-bbar> <-cMass> <-pressure $p> <-rho $r> <-bodyForce $b1 $b2>
//
//   Mass: default is the all-positive lumped mass ρVe/6 (Kadapa Eq. 56) for
//   explicit dynamics; pass -cMass for the consistent mass (implicit/eigen).
//
//   $type = "PlaneStress" or "PlaneStrain"
//
// Recorder responses:
//   "forces"              → 12 resisting force components
//   "stresses"            → σ_xx, σ_yy, σ_xy at each GP (9 total)
//   "strains"             → ε_xx, ε_yy, γ_xy at each GP (9 total)
//   "gaussPoint"          → physical x,y of each GP (6 total)
//   "material" $gp <args> → delegate to NDMaterial at GP
//   "stiffness"           → 12×12 tangent stiffness
//   "pressure"            → pressure force vector
//   "charLength"          → element characteristic length lch (crack-band)
//
// Example recorders (Python):
//   ops.recorder('Element', '-file', 'stress.txt', '-ele', 1, 'stresses')
//   ops.recorder('Element', '-file', 'gpxy.txt', '-ele', 1, 'gaussPoint')
//   ops.eleResponse(1, 'material', 1, 'stress')   # GP 1 material stress

#ifndef BezierTri6_h
#define BezierTri6_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

// Class tag ELE_TAG_BezierTri6 is defined in classTags.h (= 33000, ladruno band).

class BezierTri6 : public Element
{
  public:
    // ─── Constructors and Destructor ───────────────────────────

    // Full constructor
    BezierTri6(int tag,
               int nd1, int nd2, int nd3, int nd4, int nd5, int nd6,
               NDMaterial &m, const char *type, double thick,
               double pressure = 0.0, double rho = 0.0,
               double b1 = 0.0, double b2 = 0.0,
               bool useBbar = false, bool cMass = false);

    // Null constructor (for parallel/database reconstruction)
    BezierTri6();

    // Destructor
    ~BezierTri6();

    // ─── Element Interface ────────────────────────────────────

    // Information methods
    const char *getClassType(void) const { return "BezierTri6"; }
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // Element-size characteristic length for crack-band regularization
    // (e.g. ASDConcrete). Overrides Element's min-inter-node-distance default,
    // which on a quadratic element collapses to ~½ the edge length. See .cpp.
    double getCharacteristicLength(void);

    // Lifecycle methods
    void setDomain(Domain *theDomain);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // Stiffness and force methods
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // Serialization
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    // Output
    void Print(OPS_Stream &s, int flag = 0);

    // Response and recording
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    // Parameter interface (stress-control staging / sensitivity)
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    int displaySelf(Renderer &, int mode, float fact,
                    const char **displayModes = 0, int numModes = 0);

  protected:

  private:
    // ─── Constants ────────────────────────────────────────────
    static constexpr int NEN = 6;      // Number of element nodes
    static constexpr int NDOF = 2;     // DOFs per node
    static constexpr int NELD = 12;    // Total element DOFs (6×2)
    static constexpr int NSTRESS = 3;  // Stress/strain components (2D)
    static constexpr int NGAUSS = 3;   // Gauss points for stiffness (degree 2)
    static constexpr int NGAUSS_MASS = 6; // Gauss points for mass (degree 4)

    // ─── Computational Methods ────────────────────────────────

    // Shape functions and derivatives at a point (ξ₁, ξ₂)
    void shapeFunctions(double xi1, double xi2,
                        double N[NEN]) const;

    void shapeDerivatives(double xi1, double xi2,
                          double dN[2][NEN]) const;

    // Compute Jacobian, its determinant, and physical derivatives
    double computeJacobian(const double dN[2][NEN],
                           double J[2][2],
                           double dN_dx[2][NEN]) const;

    // Compute B matrix (strain-displacement) at a point
    void computeBMatrix(const double dN_dx[2][NEN],
                        double B[NSTRESS][NELD]) const;

    // Compute B-bar matrix (modified strain-displacement)
    void computeBBarMatrix(const double dN_dx[2][NEN],
                           const double dN_avg[2][NEN],
                           double Bbar[NSTRESS][NELD]) const;

    // Compute volume-averaged derivatives for B-bar
    void computeVolumeAveragedDerivatives(double dN_avg[2][NEN],
                                           double &volume) const;

    // Map Lagrange node positions to Bézier control points
    void computeControlPoints();

    // Get the element area
    double computeArea() const;

    // ─── Static Gauss Quadrature Data ─────────────────────────
    // 3-point rule (degree 2) for stiffness
    static const double GP3_xi[][2];
    static const double GP3_w[];

    // 6-point rule (degree 4) for mass matrix
    static const double GP6_xi[][2];
    static const double GP6_w[];

    // ─── Static Return Objects (Petracca/Abell pattern) ──────
    // Class-wide scratch matrices to avoid per-call allocation
    static Matrix K_return;   // 12×12 stiffness return matrix
    static Matrix M_return;   // 12×12 mass return matrix
    static Vector P_return;   // 12×1 force return vector

    // ─── Member Data ──────────────────────────────────────────

    // Material at each Gauss point
    NDMaterial **theMaterial;   // Array of NGAUSS material pointers

    // Connectivity
    ID connectedExternalNodes;  // 6 node tags
    Node *theNodes[NEN];        // 6 node pointers

    // Control point coordinates (Bézier, not Lagrange)
    // These may differ from node positions for curved edges
    double controlPts[NEN][2];

    // Material parameters
    double thickness;           // Element thickness
    double rho;                 // Mass density (0 = no mass)
    double pressure;            // Applied surface pressure
    double b[2];                // Body force per unit volume {b1, b2}

    // Formulation flags
    bool useBbar;               // Use B-bar for near-incompressibility
    bool cMass;                 // true = consistent mass; false = lumped (ρVe/6)

    // Applied load vector (from addLoad)
    Vector Q;                   // 12×1 applied load vector

    // Self-weight / body-force state driven by load patterns (addLoad)
    double appliedB[2];         // current body force = loadFactor·data·b
    int applyLoad;              // 1 once a SelfWeight load has been added

    // Internal working storage
    double Ki_data[NELD * NELD]; // Storage for initial stiffness
    Matrix *Ki;                  // Pointer to initial stiffness (lazy)

    // One-time print flag (Abell pattern)
    static int numInstances;
};

#endif // BezierTri6_h
