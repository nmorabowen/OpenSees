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

// LadrunoCSTPair — disjoint 2-triangle F-bar-Patch macro-element (ADR 70 P4a,
// dSNPO §15.1.9). Four nodes n1..n4 (CCW quad), triangles T1=(n1,n2,n3) and
// T2=(n1,n3,n4) split along the n1-n3 diagonal, one centroid GP each. Both
// triangles share the PATCH dilatation
//
//     J̄ = v_patch / V_patch = (J1·V1 + J2·V2) / (V1 + V2)      (eq 15.36)
//     F̄_e = (J̄ / J_e)^{1/2} F_e                                 (eq 15.34, n=2)
//
// — the volumetric cure for the constant-strain T3, on which element-level
// F-bar is a no-op (J ≡ J0; ADR 70 §3.2). The internal force of each triangle
// depends on the OTHER triangle's DOFs through J̄, so the exact consistent
// tangent carries patch-local cross blocks (dSNPO eqs 15.37/15.38) and is
// GENERALLY UNSYMMETRIC (stress-proportional coupling — it vanishes at F = I).
// Because the patch is disjoint and lives inside this single macro-element,
// standard element-local assembly carries the exact Newton operator — the
// P4 design-spike decision (ADR 70 §9).
//
//   element('LadrunoCSTPair', tag, n1,n2,n3,n4, matTag,
//           '-thick', t, '-type', 'PlaneStrain', '-rho', r, '-body', b1, b2)
//
// FINITE-STRAIN ONLY (this element exists solely as the F-bar-Patch cure):
// the material must be a FiniteStrainND2DMaterial (e.g. LogStrain2D), driven
// by the total 2×2 F̄_e via setTrialF. PlaneStrain only (family constraint:
// the finite volume weight omits the thickness stretch). Assembly runs
// through the shared LadrunoFiniteStrain2DKernel in 4-node zero-padded rows;
// the patch coupling reuses addFbarCoupling2D with the volume-weighted patch
// row ḡ = Σ_s (v_s/v_p) g_s in place of the centroid g₀.
//
// Oracle: tests/cstpair_reference.py (FD-exact at finite stress; row form ≡
// dSNPO 15.37/15.38 split; reduce-to-2-CSTs under homogeneous F; 3-RBM rank).
//
// See Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md (§9 P4).

#ifndef LadrunoCSTPair_h
#define LadrunoCSTPair_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class LadrunoCSTPair : public Element
{
  public:
    LadrunoCSTPair(int tag, int nd1, int nd2, int nd3, int nd4,
                   NDMaterial &m, double t,
                   double rho = 0.0, double b1 = 0.0, double b2 = 0.0);
    LadrunoCSTPair();
    ~LadrunoCSTPair();

    const char *getClassType(void) const { return "LadrunoCSTPair"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    double getCharacteristicLength(void);   // mean-triangle sqrt(2*A_mean)

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

  protected:

  private:
    static constexpr int numtri = 2;         // triangles (= materials = GPs)
    static constexpr int numnodes = 4;
    static constexpr int ndf = 2 * numnodes; // 8

    NDMaterial **theMaterial;                // one per triangle
    ID connectedExternalNodes;
    Node *theNodes[numnodes];

    static double matrixData[ndf * ndf];
    static Matrix K;             // 8x8 shared scratch
    static Vector P;             // 8 shared scratch
    Vector Q;
    double b[2];
    double appliedB[2];
    int applyLoad;

    double thickness;
    double rho;

    Matrix *Ki;

    // per-triangle constant reference kinematics in 4-node padded layout:
    // gradRef[t][k*4+a] = dN_a/dX_k (zeros on the off node); area[t] = A_t.
    // Returns min(detJ) (<=0 flags a degenerate/backward triangle).
    double refGradients(double gradRef[numtri][2 * numnodes],
                        double area[numtri]) const;

    int  updatePair(void);          // F_e, J̄, F̄_e -> setTrialF (both triangles)
    void formPair(int tang_flag);   // fills static P / K via the shared kernel
};

#endif
