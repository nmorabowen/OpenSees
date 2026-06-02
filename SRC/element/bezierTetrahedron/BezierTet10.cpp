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
// Created: 05/2026
//
// Description: 10-node quadratic Bézier tetrahedral element for 3D analysis.
//
// Reference:
//   Kadapa, C. "Novel quadratic Bézier triangular and tetrahedral elements
//   using existing mesh generators." Int. J. Numer. Methods Eng.,
//   2019; 117(5):543-573. doi:10.1002/nme.5967 (§5)
//
// ═══════════════════════════════════════════════════════════════════
//  MATHEMATICAL FOUNDATION
// ═══════════════════════════════════════════════════════════════════
//
//  Shape functions: Quadratic Bernstein polynomials on the tetrahedron
//  with barycentric coordinates (L1,L2,L3,L4), L4 = 1-L1-L2-L3.
//
//    N1=L1²  N2=L2²  N3=L3²  N4=L4²              (4 vertices)
//    N5=2L1L2  N6=2L2L3  N7=2L1L3                (mid-edges 1-2, 2-3, 1-3)
//    N8=2L1L4  N9=2L3L4  N10=2L2L4               (mid-edges 1-4, 3-4, 2-4)
//
//  Node order matches OpenSees TenNodeTetrahedron (jaabell/Larenas), so
//  Lagrange T10 / Gmsh meshes feed directly with only the basis swapped.
//
//  Key properties:
//    - Partition of unity: Σ Nₐ = 1
//    - Nonnegativity: Nₐ ≥ 0  ⇒  all-positive lumped mass (explicit dynamics)
//    - Equal integration: ∫ Nₐ dΩ = Vₑ/10 for all a (vertices and mid-edges)
//
//  Voigt order {xx, yy, zz, xy, yz, zx} == Kadapa Eq. 26-27 == TenNodeTet.
//  Isoparametric map uses CONTROL POINTS: x(ξ) = Σ Nₐ(ξ)·Pₐ (P = X if straight).
// ═══════════════════════════════════════════════════════════════════

#include "BezierTet10.h"

#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <OPS_Globals.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <SolidTransformation.h>         // Ladruno — geometry-method layer
#include <SolidTransformationLinear.h>   // Ladruno — default/fallback (identity)

#include <math.h>
#include <string.h>
#include <stdio.h>

// ═══════════════════════════════════════════════════════════════════
//  STATIC DATA
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::numInstances = 0;

Matrix BezierTet10::K_return(NELD, NELD);
Matrix BezierTet10::M_return(NELD, NELD);
Vector BezierTet10::P_return(NELD);

// ─── 4-Point Tetrahedral Gauss Quadrature (degree 2) ──────────
// Symmetric rule; alpha/beta == TenNodeTetrahedron. Each weight = 1/24,
// so the weights sum to 1/6 (the reference tetrahedron volume).
//   a = (5 + 3√5)/20 ≈ 0.585410,  b = (5 - √5)/20 ≈ 0.138197
// Stored as the 3 free barycentric coords (L1,L2,L3); L4 = 1-L1-L2-L3.
const double BezierTet10::GP4_L[][3] = {
    {0.585410196624968, 0.138196601125011, 0.138196601125011},
    {0.138196601125011, 0.585410196624968, 0.138196601125011},
    {0.138196601125011, 0.138196601125011, 0.585410196624968},
    {0.138196601125011, 0.138196601125011, 0.138196601125011}
};
const double BezierTet10::GP4_w[] = {
    1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0
};

// ─── 1D Gauss-Legendre (4-pt) on [0,1] ────────────────────────
// Collapsed (Duffy) into a 4×4×4 rule for the consistent mass: exact to
// degree 7, which covers the degree-6 integrand the (1-a)² Duffy Jacobian
// produces from the degree-4 product NₐNᵦ (ADR O12).
const double BezierTet10::GL4_t[] = {
    0.069431844202974, 0.330009478207572,
    0.669990521792428, 0.930568155797026
};
const double BezierTet10::GL4_w[] = {
    0.173927422568727, 0.326072577431273,
    0.326072577431273, 0.173927422568727
};

// ─── Edge → (vertexA, vertexB) for mid-edge nodes 5..10 (idx 4..9) ──
// Matches TenNodeTetrahedron edge order (jaabell/Larenas N9↔N10 swap).
const int BezierTet10::edgeV[6][2] = {
    {0, 1},   // node 5  : edge (1-2)
    {1, 2},   // node 6  : edge (2-3)
    {0, 2},   // node 7  : edge (1-3)
    {0, 3},   // node 8  : edge (1-4)
    {2, 3},   // node 9  : edge (3-4)
    {1, 3}    // node 10 : edge (2-4)
};


// ═══════════════════════════════════════════════════════════════════
//  CONSTRUCTORS AND DESTRUCTOR
// ═══════════════════════════════════════════════════════════════════

BezierTet10::BezierTet10(int tag,
                         int nd1, int nd2, int nd3, int nd4, int nd5,
                         int nd6, int nd7, int nd8, int nd9, int nd10,
                         NDMaterial &m, double r,
                         double bx, double by, double bz,
                         bool bbar, bool cmass, double press,
                         int geomMethodID)
  : Element(tag, ELE_TAG_BezierTet10),
    theMaterial(0),
    connectedExternalNodes(NEN),
    rho(r), pressure(press),
    useBbar(bbar), cMass(cmass),
    Q(NELD), applyLoad(0), Ki(0), theGeom(0)   // Ladruno — geometry method
{
    // Ladruno: geometry-method layer (linear default / corot). Linear is the
    // identity wrapper so std/bbar run bit-for-bit unchanged.
    theGeom = SolidTransformation::create(geomMethodID);
    if (theGeom == 0)
        theGeom = new SolidTransformationLinear();   // safe fallback (unknown id)

    // One-time attribution print (Abell pattern) — Ladruño banner + authors
    if (numInstances == 0) {
        numInstances++;
        opserr << "\n"
" ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄\n"
"███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███\n"
"███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███\n"
"███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███\n"
"███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███\n"
"███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███\n"
"███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███\n"
"█████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀\n"
"▀                                   ███    ███\n";
        opserr << "BezierTet10 - Quadratic Bezier Tetrahedron\n"
               << "  Authors: Nicolas Mora Bowen, Patricio Palacios, "
                  "Jose Abell, Guppi (Ladruño)\n"
               << "  Ref: Kadapa, IJNME 2019; 117(5):543-573. "
                  "doi:10.1002/nme.5967\n";
    }

    // Store body forces
    b[0] = bx; b[1] = by; b[2] = bz;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;

    // Store connectivity
    int nd[NEN] = {nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, nd9, nd10};
    for (int i = 0; i < NEN; i++)
        connectedExternalNodes(i) = nd[i];

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = 0;
        controlPts[i][0] = controlPts[i][1] = controlPts[i][2] = 0.0;
    }

    for (int i = 0; i < NELD * NELD; i++)
        Ki_data[i] = 0.0;

    // Allocate 3D materials at the NGAUSS Gauss points
    theMaterial = new NDMaterial*[NGAUSS];
    for (int i = 0; i < NGAUSS; i++) {
        theMaterial[i] = m.getCopy("ThreeDimensional");
        if (theMaterial[i] == 0) {
            opserr << "BezierTet10::BezierTet10 - failed to get a 3D copy "
                      "of the material\n";
            exit(-1);
        }
    }
}


BezierTet10::BezierTet10()
  : Element(0, ELE_TAG_BezierTet10),
    theMaterial(0),
    connectedExternalNodes(NEN),
    rho(0.0), pressure(0.0),
    useBbar(false), cMass(false),
    Q(NELD), applyLoad(0), Ki(0), theGeom(0)   // Ladruno — geometry method
{
    // Ladruno: broker/database reconstruction path — recvSelf rebuilds the real
    // method from the serialized id; default to the identity wrapper meanwhile.
    theGeom = new SolidTransformationLinear();

    b[0] = b[1] = b[2] = 0.0;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = 0;
        controlPts[i][0] = controlPts[i][1] = controlPts[i][2] = 0.0;
    }

    for (int i = 0; i < NELD * NELD; i++)
        Ki_data[i] = 0.0;
}


BezierTet10::~BezierTet10()
{
    if (theMaterial != 0) {
        for (int i = 0; i < NGAUSS; i++)
            if (theMaterial[i])
                delete theMaterial[i];
        delete[] theMaterial;
    }
    if (Ki != 0)
        delete Ki;
    if (theGeom != 0)        // Ladruno — geometry-method layer
        delete theGeom;
}


// ═══════════════════════════════════════════════════════════════════
//  INFORMATION METHODS
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::getNumExternalNodes(void) const { return NEN; }
const ID &BezierTet10::getExternalNodes(void) { return connectedExternalNodes; }
Node **BezierTet10::getNodePtrs(void) { return theNodes; }
int BezierTet10::getNumDOF(void) { return NELD; }


// ═══════════════════════════════════════════════════════════════════
//  LIFECYCLE METHODS
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
        for (int i = 0; i < NEN; i++)
            theNodes[i] = 0;
        return;
    }

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
        if (theNodes[i] == 0) {
            opserr << "BezierTet10::setDomain -- node "
                   << connectedExternalNodes(i) << " does not exist\n";
            return;
        }
        if (theNodes[i]->getNumberDOF() != NDOF) {
            opserr << "BezierTet10::setDomain -- node "
                   << connectedExternalNodes(i)
                   << " has " << theNodes[i]->getNumberDOF()
                   << " DOFs, expected " << NDOF << "\n";
            return;
        }
    }

    computeControlPoints();

    // ─── v1 straight-sided guard (ADR 06_bezier_tet10 D9′) ────────
    //  v1 is validated for STRAIGHT-SIDED elements only: all 6 mid-edge
    //  nodes at the edge midpoints ⇒ control points = nodes (P = X).
    //  A curved mesh breaks the interpolatory Dirichlet BCs and under-
    //  integrates BᵀDB. Warn; curved support + Eq. 14 is a follow-up.
    for (int e = 0; e < 6; e++) {
        const Vector &Xc1 = theNodes[edgeV[e][0]]->getCrds();
        const Vector &Xc2 = theNodes[edgeV[e][1]]->getCrds();
        const Vector &Xm  = theNodes[4 + e]->getCrds();
        double L2 = 0.0, d2 = 0.0;
        for (int k = 0; k < 3; k++) {
            double mid = 0.5 * (Xc1(k) + Xc2(k));
            d2 += (Xm(k) - mid) * (Xm(k) - mid);
            L2 += (Xc2(k) - Xc1(k)) * (Xc2(k) - Xc1(k));
        }
        if (L2 > 0.0 && d2 > 1.0e-12 * L2) {   // > ~1e-6 relative deviation
            opserr << "WARNING BezierTet10 " << this->getTag()
                   << " - mid-edge node " << connectedExternalNodes(4 + e)
                   << " is off the edge midpoint (curved element). v1 is "
                      "validated for straight-sided meshes only.\n";
        }
    }

    this->DomainComponent::setDomain(theDomain);
}


int BezierTet10::commitState()
{
    int retVal = 0;
    if ((retVal = this->Element::commitState()) != 0)
        opserr << "BezierTet10::commitState() - failed in base class\n";
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->commitState();
    retVal += theGeom->commitState();   // Ladruno (corot is stateless: no-op)
    return retVal;
}


int BezierTet10::revertToLastCommit()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToLastCommit();
    retVal += theGeom->revertToLastCommit();   // Ladruno (corot stateless)
    return retVal;
}


int BezierTet10::revertToStart()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToStart();
    retVal += theGeom->revertToStart();   // Ladruno (corot stateless)
    if (theNodes[0] != 0)
        computeControlPoints();
    return retVal;
}


int BezierTet10::update()
{
    // Ladruno (seam 0+2): refresh the geometry method and localize the trial
    // displacement into the core frame. For -geom linear this is the identity
    // (uCore == uGlobal), so the strain below is bit-for-bit the direct kernel;
    // for -geom corot it is the de-rotated (Rᵀ) displacement. Strains ε = B·uCore
    // (or B̄·uCore) are then pushed to the materials.
    const Vector &u = this->computeLocalDisp();

    int ret = 0;

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        if (factor < 0.0) {
            opserr << "BezierTet10::update() - degenerate Jacobian at GP " << gp
                   << " (element " << this->getTag() << ")\n";
            return -1;
        }

        Vector strain(NSTRESS);
        for (int i = 0; i < NSTRESS; i++) {
            double sum = 0.0;
            for (int j = 0; j < NELD; j++)
                sum += B[i][j] * u(j);
            strain(i) = sum;
        }

        ret += theMaterial[gp]->setTrialStrain(strain);
    }

    return ret;
}


// ═══════════════════════════════════════════════════════════════════
//  STIFFNESS AND FORCE
// ═══════════════════════════════════════════════════════════════════

const Matrix &BezierTet10::getTangentStiff()
{
    // K = Σᵢ wᵢ |Jᵢ| · Bᵢᵀ Dᵢ Bᵢ   (B → B̄ for B-bar), assembled in the CORE
    // frame in ONE pass that also yields the core internal force, then globalized
    // (seam 3). Ladruno: refresh the geometry method so K and the fCore used for
    // K_geo share one fresh R.
    this->computeLocalDisp();

    // one Gauss pass → core K (into K_return) AND the core internal force fCore.
    static Vector fCore(NELD);
    this->formCore(1, fCore, &K_return);

    // seam 3: globalize the core-frame K and add the corotational geometric
    // stiffness K_geo, which depends on fCore — the SAME force globalizeForce
    // rotates in getResistingForce (both via formCore), so they match by
    // construction. Identity (kGlobal=kCore, K_geo=0) for -geom linear.  // Ladruno
    theGeom->globalizeStiff(K_return, fCore, K_return);

    return K_return;
}


const Matrix &BezierTet10::getInitialStiff()
{
    if (Ki != 0)
        return *Ki;

    K_return.Zero();

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        if (factor < 0.0) {
            opserr << "BezierTet10::getInitialStiff - degenerate Jacobian at GP "
                   << gp << " (element " << this->getTag() << ")\n";
            continue;
        }

        const Matrix &D = theMaterial[gp]->getInitialTangent();

        for (int I = 0; I < NELD; I++) {
            for (int Jc = I; Jc < NELD; Jc++) {
                double sum = 0.0;
                for (int k = 0; k < NSTRESS; k++)
                    for (int l = 0; l < NSTRESS; l++)
                        sum += B[k][I] * D(k, l) * B[l][Jc];
                K_return(I, Jc) += factor * sum;
                if (I != Jc)
                    K_return(Jc, I) += factor * sum;
            }
        }
    }

    // seam 3 (reference config): pin the geometry frame to R = I by refreshing
    // theGeom with cur == ref, THEN globalize with a zero core force. Without the
    // update-to-reference, globalizeStiff would reuse a STALE current-config R if
    // getInitialStiff is queried mid-analysis. R = I is deterministic, so caching
    // Ki below stays valid. Identity for -geom linear.  // Ladruno
    if (theNodes[0] != 0) {
        static Matrix refC(NEN, 3);
        for (int i = 0; i < NEN; i++) {
            const Vector &X = theNodes[i]->getCrds();
            for (int d = 0; d < 3; d++)
                refC(i, d) = X(d);
        }
        theGeom->update(NEN, refC, refC);          // Rmat = I
        static Vector zeroF(NELD);
        zeroF.Zero();
        theGeom->globalizeStiff(K_return, zeroF, K_return);
    }

    Ki = new Matrix(Ki_data, NELD, NELD);
    *Ki = K_return;
    return *Ki;
}


const Matrix &BezierTet10::getMass()
{
    // Lumped (default) or consistent (-cMass). Density: element rho
    // overrides; otherwise the material's own density.
    M_return.Zero();

    double rhoEff = rho;
    if (rhoEff == 0.0) {
        double sum = 0.0;
        for (int i = 0; i < NGAUSS; i++)
            sum += theMaterial[i]->getRho();
        rhoEff = sum / NGAUSS;
    }
    if (rhoEff == 0.0)
        return M_return;

    if (!cMass) {
        // ─── Lumped mass (Kadapa Eq. 57) ──────────────────────
        //   Mᵉ = (ρ Vₑ / 10) · diag[1₁₀, 1₁₀, 1₁₀]
        //  Each of the 10 nodes gets an EQUAL share ρVₑ/10 — exact for the
        //  quadratic Bernstein tet (∫Nₐ = Vₑ/10 ∀a) and ALL POSITIVE, the
        //  reason this element exists for explicit dynamics.
        double Ve = this->computeVolume();
        double m  = rhoEff * Ve / 10.0;
        for (int a = 0; a < NEN; a++) {
            M_return(3*a,     3*a)     = m;
            M_return(3*a + 1, 3*a + 1) = m;
            M_return(3*a + 2, 3*a + 2) = m;
        }
        return M_return;
    }

    // ─── Consistent mass:  M = ∫_Ω ρ Nᵤᵀ Nᵤ dΩ ────────────────
    //  4×4×4 collapsed-Duffy rule (degree 7) integrates the degree-6
    //  integrand exactly for straight edges. All entries ≥ 0 (Nₐ ≥ 0).
    for (int i = 0; i < 4; i++) {
        double a = GL4_t[i];
        for (int j = 0; j < 4; j++) {
            double bb = GL4_t[j];
            for (int k = 0; k < 4; k++) {
                double cc = GL4_t[k];
                double L1 = a;
                double L2 = (1.0 - a) * bb;
                double L3 = (1.0 - a) * (1.0 - bb) * cc;
                double jac = (1.0 - a) * (1.0 - a) * (1.0 - bb);
                double w = GL4_w[i] * GL4_w[j] * GL4_w[k] * jac;

                double N[NEN];
                shapeFunctions(L1, L2, L3, N);

                double dN[3][NEN];
                shapeDerivatives(L1, L2, L3, dN);

                double J[3][3];
                double dN_dx[3][NEN];
                double detJ = computeJacobian(dN, J, dN_dx);

                double factor = w * fabs(detJ) * rhoEff;

                for (int p = 0; p < NEN; p++) {
                    for (int q = 0; q < NEN; q++) {
                        double val = factor * N[p] * N[q];
                        M_return(3*p,     3*q)     += val;
                        M_return(3*p + 1, 3*q + 1) += val;
                        M_return(3*p + 2, 3*q + 2) += val;
                    }
                }
            }
        }
    }

    return M_return;
}


void BezierTet10::zeroLoad()
{
    Q.Zero();
    applyLoad = 0;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;
}


int BezierTet10::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // Supports the standard self-weight (gravity) load so the body force can
    // be ramped by a load pattern / time series (same idiom as FourNodeQuad).
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor * data(0) * b[0];
        appliedB[1] += loadFactor * data(1) * b[1];
        appliedB[2] += loadFactor * data(2) * b[2];
        return 0;
    }

    opserr << "BezierTet10::addLoad - load type " << type
           << " unknown for element " << this->getTag()
           << " (only SelfWeight is supported)\n";
    return -1;
}


int BezierTet10::addInertiaLoadToUnbalance(const Vector &accel)
{
    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;
    if (!hasMass)
        return 0;

    static Vector a(NELD);
    for (int i = 0; i < NEN; i++) {
        const Vector &Raccel = theNodes[i]->getRV(accel);
        if (Raccel.Size() != NDOF) {
            opserr << "BezierTet10::addInertiaLoadToUnbalance - "
                   << "matrix and target sizes mismatch\n";
            return -1;
        }
        a(3*i)     = Raccel(0);
        a(3*i + 1) = Raccel(1);
        a(3*i + 2) = Raccel(2);
    }

    Q.addMatrixVector(1.0, M, a, 1.0);
    return 0;
}


const Vector &BezierTet10::getResistingForce()
{
    // F = globalize(∫_Ω Bᵀσ dΩ)  −  body force  −  pressure  −  Q.
    //
    // Ladruno (corot load-frame contract): the INTERNAL force ∫Bᵀσ is assembled
    // in the core frame and rotated to global by globalizeForce (seam 3). The
    // fixed-direction EXTERNAL loads — body force (gravity), the +z pressure
    // hack, and Q (applied / inertia) — are applied in the GLOBAL frame AFTER
    // globalizeForce so the rotation R never co-rotates them. This keeps the
    // fCore fed to globalizeStiff equal to pure ∫Bᵀσ (so K_geo carries no
    // spurious external-load term) and makes the f/K paths trivially consistent.
    // For -geom linear globalizeForce is the identity and this reproduces the
    // direct kernel.
    this->computeLocalDisp();                     // seam 0+2: refresh frame

    static Vector fInt(NELD);
    this->formCore(0, fInt, 0);                   // core-frame ∫Bᵀσ (no tangent)
    theGeom->globalizeForce(fInt, fInt);          // seam 3: → global
    P_return = fInt;

    // fixed-direction external loads, GLOBAL frame
    double bx = (applyLoad == 1) ? appliedB[0] : b[0];
    double by = (applyLoad == 1) ? appliedB[1] : b[1];
    double bz = (applyLoad == 1) ? appliedB[2] : b[2];
    if (bx != 0.0 || by != 0.0 || bz != 0.0 || pressure != 0.0) {
        for (int gp = 0; gp < NGAUSS; gp++) {
            double w = GP4_w[gp];
            double N[NEN];
            shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);

            double dN[3][NEN];
            shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

            double J[3][3];
            double dN_dx[3][NEN];
            double detJ = computeJacobian(dN, J, dN_dx);
            double factor = w * fabs(detJ);

            for (int a = 0; a < NEN; a++) {
                P_return(3*a)     -= factor * N[a] * bx;
                P_return(3*a + 1) -= factor * N[a] * by;
                P_return(3*a + 2) -= factor * N[a] * bz;
                if (pressure != 0.0)   // +z volume hack (mirrors BezierTri6)
                    P_return(3*a + 2) -= factor * N[a] * pressure;
            }
        }
    }

    P_return.addVector(1.0, Q, -1.0);
    return P_return;
}


// ─── Geometry-method seam helpers (Ladruno) ───────────────────────────
// computeLocalDisp: seam 0+2 — refresh theGeom from the current geometry and
// return the localized (core-frame) trial displacement. Identity for -geom
// linear (uCore == uGlobal). Mirrors LadrunoBrick::computeLocalDisp.
const Vector &BezierTet10::computeLocalDisp(void)
{
    static Matrix refCrds(NEN, 3), curCrds(NEN, 3);
    static Vector uGlobal(NELD), uCore(NELD);

    for (int i = 0; i < NEN; i++) {
        const Vector &X = theNodes[i]->getCrds();
        const Vector &u = theNodes[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) {
            refCrds(i, d)     = X(d);
            curCrds(i, d)     = X(d) + u(d);
            uGlobal(3 * i + d) = u(d);
        }
    }

    theGeom->update(NEN, refCrds, curCrds);
    theGeom->localizeDisp(uGlobal, uCore);
    return uCore;
}

// formBAtGauss: the single guarded strain-displacement assembly shared by every
// Gauss loop (update / formCore / getInitialStiff). Returns w·|detJ|, or −1 if
// the element is degenerate (detJ == 0, where computeJacobian leaves dN_dx unset
// — so the degenerate check MUST live here, not be duplicated/forgotten per loop).
double BezierTet10::formBAtGauss(int gp, const double dN_avg[3][NEN],
                                 double B[NSTRESS][NELD]) const
{
    double dN[3][NEN];
    shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

    double J[3][3];
    double dN_dx[3][NEN];
    double detJ = computeJacobian(dN, J, dN_dx);
    if (fabs(detJ) <= 0.0)
        return -1.0;

    if (useBbar)
        computeBBarMatrix(dN_dx, dN_avg, B);
    else
        computeBMatrix(dN_dx, B);

    return GP4_w[gp] * fabs(detJ);
}

// formCore: ONE Gauss pass → core-frame internal force fInt = ∫Bᵀσ dΩ (always)
// and, when tangFlag, the core tangent K = ∫BᵀDB dΩ. NO body force / pressure /
// Q. getResistingForce calls formCore(0,…) and getTangentStiff formCore(1,…), so
// the fCore globalizeForce rotates equals the fCore globalizeStiff uses for K_geo
// BY CONSTRUCTION — and the tangent path no longer runs a second Bᵀσ loop.
void BezierTet10::formCore(int tangFlag, Vector &fInt, Matrix *K)
{
    fInt.Zero();
    if (tangFlag && K != 0)
        K->Zero();

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        if (factor < 0.0) {
            opserr << "BezierTet10::formCore - degenerate Jacobian at GP " << gp
                   << " (element " << this->getTag() << ")\n";
            continue;
        }

        const Vector &sigma = theMaterial[gp]->getStress();
        for (int I = 0; I < NELD; I++) {
            double s = 0.0;
            for (int k = 0; k < NSTRESS; k++)
                s += B[k][I] * sigma(k);
            fInt(I) += factor * s;
        }

        if (tangFlag && K != 0) {
            const Matrix &D = theMaterial[gp]->getTangent();
            for (int I = 0; I < NELD; I++) {
                for (int Jc = I; Jc < NELD; Jc++) {   // symmetric
                    double sum = 0.0;
                    for (int k = 0; k < NSTRESS; k++)
                        for (int l = 0; l < NSTRESS; l++)
                            sum += B[k][I] * D(k, l) * B[l][Jc];
                    (*K)(I, Jc) += factor * sum;
                    if (I != Jc)
                        (*K)(Jc, I) += factor * sum;
                }
            }
        }
    }
}


const Vector &BezierTet10::getResistingForceIncInertia()
{
    this->getResistingForce();

    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;

    if (hasMass) {
        static Vector a(NELD);
        for (int i = 0; i < NEN; i++) {
            const Vector &accel = theNodes[i]->getTrialAccel();
            a(3*i)     = accel(0);
            a(3*i + 1) = accel(1);
            a(3*i + 2) = accel(2);
        }
        P_return.addMatrixVector(1.0, M, a, 1.0);
    }

    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
        const Vector &v = this->getRayleighDampingForces();
        P_return += v;
    }

    return P_return;
}


// ═══════════════════════════════════════════════════════════════════
//  SHAPE FUNCTIONS AND DERIVATIVES
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::shapeFunctions(double L1, double L2, double L3,
                                 double N[NEN]) const
{
    double L4 = 1.0 - L1 - L2 - L3;

    N[0] = L1 * L1;        // N1  vertex 1
    N[1] = L2 * L2;        // N2  vertex 2
    N[2] = L3 * L3;        // N3  vertex 3
    N[3] = L4 * L4;        // N4  vertex 4
    N[4] = 2.0 * L1 * L2;  // N5  edge (1-2)
    N[5] = 2.0 * L2 * L3;  // N6  edge (2-3)
    N[6] = 2.0 * L1 * L3;  // N7  edge (1-3)
    N[7] = 2.0 * L1 * L4;  // N8  edge (1-4)
    N[8] = 2.0 * L3 * L4;  // N9  edge (3-4)
    N[9] = 2.0 * L2 * L4;  // N10 edge (2-4)
}


void BezierTet10::shapeDerivatives(double L1, double L2, double L3,
                                   double dN[3][NEN]) const
{
    //  ∂Nₐ/∂L1, ∂Nₐ/∂L2, ∂Nₐ/∂L3 with L4 = 1-L1-L2-L3 (∂L4/∂Lk = -1).
    double L4 = 1.0 - L1 - L2 - L3;

    // ∂N/∂L1
    dN[0][0] =  2.0 * L1;          // L1²
    dN[0][1] =  0.0;               // L2²
    dN[0][2] =  0.0;               // L3²
    dN[0][3] = -2.0 * L4;          // L4²
    dN[0][4] =  2.0 * L2;          // 2 L1 L2
    dN[0][5] =  0.0;               // 2 L2 L3
    dN[0][6] =  2.0 * L3;          // 2 L1 L3
    dN[0][7] =  2.0 * (L4 - L1);   // 2 L1 L4
    dN[0][8] = -2.0 * L3;          // 2 L3 L4
    dN[0][9] = -2.0 * L2;          // 2 L2 L4

    // ∂N/∂L2
    dN[1][0] =  0.0;
    dN[1][1] =  2.0 * L2;
    dN[1][2] =  0.0;
    dN[1][3] = -2.0 * L4;
    dN[1][4] =  2.0 * L1;
    dN[1][5] =  2.0 * L3;
    dN[1][6] =  0.0;
    dN[1][7] = -2.0 * L1;
    dN[1][8] = -2.0 * L3;
    dN[1][9] =  2.0 * (L4 - L2);

    // ∂N/∂L3
    dN[2][0] =  0.0;
    dN[2][1] =  0.0;
    dN[2][2] =  2.0 * L3;
    dN[2][3] = -2.0 * L4;
    dN[2][4] =  0.0;
    dN[2][5] =  2.0 * L2;
    dN[2][6] =  2.0 * L1;
    dN[2][7] = -2.0 * L1;
    dN[2][8] =  2.0 * (L4 - L3);
    dN[2][9] = -2.0 * L2;
}


// ═══════════════════════════════════════════════════════════════════
//  JACOBIAN AND B-MATRIX
// ═══════════════════════════════════════════════════════════════════

double BezierTet10::computeJacobian(const double dN[3][NEN],
                                    double J[3][3],
                                    double dN_dx[3][NEN]) const
{
    //  J[i][j] = ∂x_j/∂L_i = Σₐ ∂Nₐ/∂L_i · Pₐ_j   (control points).
    //  dN/dx = J⁻¹ · dN/dL. Returns SIGNED detJ; callers use |detJ| for the
    //  integration measure (orientation-robust — ADR O11).
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            J[i][j] = 0.0;

    for (int a = 0; a < NEN; a++)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                J[i][j] += dN[i][a] * controlPts[a][j];

    double detJ =
        J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]);

    if (fabs(detJ) > 0.0) {
        double inv = 1.0 / detJ;
        double Ji[3][3];
        Ji[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * inv;
        Ji[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]) * inv;
        Ji[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * inv;
        Ji[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]) * inv;
        Ji[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * inv;
        Ji[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]) * inv;
        Ji[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * inv;
        Ji[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]) * inv;
        Ji[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * inv;

        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_dx[i][a] = Ji[i][0]*dN[0][a]
                            + Ji[i][1]*dN[1][a]
                            + Ji[i][2]*dN[2][a];
    }

    return detJ;
}


void BezierTet10::computeBMatrix(const double dN_dx[3][NEN],
                                 double B[NSTRESS][NELD]) const
{
    //  3D strain-displacement (Kadapa Eq. 38). Voigt {xx,yy,zz,xy,yz,zx}.
    //   Bₐ = [ N,x  0    0   ]
    //        [ 0    N,y  0   ]
    //        [ 0    0    N,z ]
    //        [ N,y  N,x  0   ]
    //        [ 0    N,z  N,y ]
    //        [ N,z  0    N,x ]
    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            B[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        double bx = dN_dx[0][a], by = dN_dx[1][a], bz = dN_dx[2][a];
        int cx = 3*a, cy = 3*a + 1, cz = 3*a + 2;

        B[0][cx] = bx;
        B[1][cy] = by;
        B[2][cz] = bz;
        B[3][cx] = by;  B[3][cy] = bx;
        B[4][cy] = bz;  B[4][cz] = by;
        B[5][cx] = bz;  B[5][cz] = bx;
    }
}


void BezierTet10::computeBBarMatrix(const double dN_dx[3][NEN],
                                    const double dN_avg[3][NEN],
                                    double Bbar[NSTRESS][NELD]) const
{
    //  B-bar (Kadapa Eq. 45, full 3D). The three normal-strain rows get the
    //  volumetric (1/3) split with the volume-averaged derivatives B̄; the
    //  three shear rows are unchanged.
    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            Bbar[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        double B1 = dN_dx[0][a], B2 = dN_dx[1][a], B3 = dN_dx[2][a];
        double A1 = dN_avg[0][a], A2 = dN_avg[1][a], A3 = dN_avg[2][a];
        int cx = 3*a, cy = 3*a + 1, cz = 3*a + 2;

        // normal-strain rows (volumetric split)
        Bbar[0][cx] = (A1 + 2.0*B1) / 3.0;
        Bbar[1][cx] = (A1 - B1) / 3.0;
        Bbar[2][cx] = (A1 - B1) / 3.0;
        Bbar[0][cy] = (A2 - B2) / 3.0;
        Bbar[1][cy] = (A2 + 2.0*B2) / 3.0;
        Bbar[2][cy] = (A2 - B2) / 3.0;
        Bbar[0][cz] = (A3 - B3) / 3.0;
        Bbar[1][cz] = (A3 - B3) / 3.0;
        Bbar[2][cz] = (A3 + 2.0*B3) / 3.0;

        // shear rows (unchanged)
        Bbar[3][cx] = B2;  Bbar[3][cy] = B1;
        Bbar[4][cy] = B3;  Bbar[4][cz] = B2;
        Bbar[5][cx] = B3;  Bbar[5][cz] = B1;
    }
}


void BezierTet10::computeVolumeAveragedDerivatives(double dN_avg[3][NEN],
                                                   double &volume) const
{
    //  B̄ₐ = (1/Vₑ) ∫_Ω ∂Nₐ/∂xᵢ dΩ  (4-pt rule, Kadapa Remark 3: full rule).
    for (int i = 0; i < 3; i++)
        for (int a = 0; a < NEN; a++)
            dN_avg[i][a] = 0.0;

    volume = 0.0;

    for (int gp = 0; gp < NGAUSS; gp++) {
        double w = GP4_w[gp];
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

        double J[3][3];
        double dN_dx[3][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        double factor = w * fabs(detJ);
        volume += factor;

        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_avg[i][a] += factor * dN_dx[i][a];
    }

    if (volume > 0.0) {
        double invVol = 1.0 / volume;
        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_avg[i][a] *= invVol;
    }
}


// ═══════════════════════════════════════════════════════════════════
//  LAGRANGE → BÉZIER CONTROL POINT MAPPING
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::computeControlPoints()
{
    //  Kadapa Eq. 11-13, edge-wise:
    //    vertices:  Pᵢ = Xᵢ                                  (i = 1..4)
    //    mid-edge:  P = 2[X_mid - 0.25 Xₐ - 0.25 X_b]        for edge (a,b)
    //  Straight-sided ⇒ P = X for every node.
    double X[NEN][3];
    for (int i = 0; i < NEN; i++) {
        const Vector &crds = theNodes[i]->getCrds();
        X[i][0] = crds(0);
        X[i][1] = crds(1);
        X[i][2] = crds(2);
    }

    // Vertices
    for (int i = 0; i < 4; i++)
        for (int d = 0; d < 3; d++)
            controlPts[i][d] = X[i][d];

    // Mid-edges
    for (int e = 0; e < 6; e++) {
        int a = edgeV[e][0], bnode = edgeV[e][1];
        for (int d = 0; d < 3; d++)
            controlPts[4 + e][d] =
                2.0 * (X[4 + e][d] - 0.25 * X[a][d] - 0.25 * X[bnode][d]);
    }
}


double BezierTet10::computeVolume() const
{
    double vol = 0.0;
    for (int gp = 0; gp < NGAUSS; gp++) {
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);
        double J[3][3];
        double dN_dx[3][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);
        vol += GP4_w[gp] * fabs(detJ);
    }
    return vol;
}


// ═══════════════════════════════════════════════════════════════════
//  CHARACTERISTIC LENGTH  (crack-band regularization, e.g. ASDConcrete)
// ═══════════════════════════════════════════════════════════════════
//
//  Crack-band materials (ASDConcrete3D, etc.) call
//  ops_TheActiveElement->getCharacteristicLength() exactly once — on the
//  first setTrialStrain — to regularize the softening branch by element size.
//
//  Element's base default returns the MINIMUM inter-node distance. On a
//  quadratic element that distance is corner-to-mid-edge ≈ ½ the true edge
//  length, which under-estimates the band width and over-softens the response.
//
//  We instead return an element-size equivalent from the integrated volume:
//  the leg of a right tetrahedron (three mutually perpendicular legs) of
//  equal volume,
//
//      lch = cbrt(6 · V),
//
//  the 3D analogue of BezierTri6's sqrt(2·A). It recovers the true edge length
//  for a right-corner tet and is geometry-true for distorted straight-sided
//  elements; curved Bézier edges shift the factor only in the safe
//  (under-estimating) direction.
double BezierTet10::getCharacteristicLength(void)
{
    double V = this->computeVolume();
    if (V <= 0.0)
        return Element::getCharacteristicLength();  // degenerate: fall back
    return cbrt(6.0 * V);
}


// ═══════════════════════════════════════════════════════════════════
//  SERIALIZATION (sendSelf / recvSelf)
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // tag + 10 nodes + matClassTag + matDbTag + useBbar + cMass + geomID = 16
    // Ladruno: a dedicated slot 15 carries the geometry-method id (the existing
    // 15-slot buffer was full — no packed slot to reuse, so widen to 16).
    static ID iData(16);
    iData(0) = this->getTag();
    for (int i = 0; i < NEN; i++)
        iData(i + 1) = connectedExternalNodes(i);
    iData(11) = theMaterial[0]->getClassTag();
    int matDbTag = theMaterial[0]->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->setDbTag(matDbTag);
    }
    iData(12) = matDbTag;
    iData(13) = useBbar ? 1 : 0;
    iData(14) = cMass ? 1 : 0;
    iData(15) = theGeom->getMethodID();   // Ladruno — geometry method

    res += theChannel.sendID(this->getDbTag(), commitTag, iData);

    // rho, pressure, b[0..2] = 5
    static Vector dData(5);
    dData(0) = rho;
    dData(1) = pressure;
    dData(2) = b[0];
    dData(3) = b[1];
    dData(4) = b[2];

    res += theChannel.sendVector(this->getDbTag(), commitTag, dData);

    for (int i = 0; i < NGAUSS; i++)
        res += theMaterial[i]->sendSelf(commitTag, theChannel);

    return res;
}


int BezierTet10::recvSelf(int commitTag, Channel &theChannel,
                          FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static ID iData(16);
    res += theChannel.recvID(this->getDbTag(), commitTag, iData);

    this->setTag(iData(0));
    for (int i = 0; i < NEN; i++)
        connectedExternalNodes(i) = iData(i + 1);
    int matClassTag = iData(11);
    int matDbTag = iData(12);
    useBbar = (iData(13) == 1);
    cMass = (iData(14) == 1);

    // Ladruno: rebuild the geometry method from the serialized id (Linear
    // fallback for an unknown/zero id) — else a parallel worker silently
    // reverts -geom corot to linear.
    if (theGeom != 0)
        delete theGeom;
    theGeom = SolidTransformation::create(iData(15));
    if (theGeom == 0)
        theGeom = new SolidTransformationLinear();

    static Vector dData(5);
    res += theChannel.recvVector(this->getDbTag(), commitTag, dData);

    rho = dData(0);
    pressure = dData(1);
    b[0] = dData(2);
    b[1] = dData(3);
    b[2] = dData(4);

    if (theMaterial == 0) {
        theMaterial = new NDMaterial*[NGAUSS];
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i] = 0;
    }

    for (int i = 0; i < NGAUSS; i++) {
        if (theMaterial[i] == 0) {
            theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
            if (theMaterial[i] == 0) {
                opserr << "BezierTet10::recvSelf - material creation failed\n";
                return -1;
            }
        }
        theMaterial[i]->setDbTag(matDbTag);
        res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    return res;
}


// ═══════════════════════════════════════════════════════════════════
//  PRINT AND DISPLAY
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        s << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTet10\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "]}";
        return;
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nBezierTet10, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tmass density:  " << rho << endln;
        s << "\tB-bar formulation: " << (useBbar ? "yes" : "no") << endln;
        s << "\tbody forces:  " << b[0] << " " << b[1] << " " << b[2] << endln;

        s << "\tStress (Gauss):" << endln;
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTet10\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "], ";
        s << "\"bbar\": " << (useBbar ? "true" : "false") << ", ";
        s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
    }
}


int BezierTet10::displaySelf(Renderer &theViewer, int displayMode,
                             float fact, const char **modes, int numMode)
{
    // Draw the 6 edges, each as 2 segments through its mid-edge node.
    static Vector v1(3), v2(3);

    for (int e = 0; e < 6; e++) {
        int seq[3] = {edgeV[e][0], 4 + e, edgeV[e][1]};
        for (int seg = 0; seg < 2; seg++) {
            int n1 = seq[seg];
            int n2 = seq[seg + 1];

            const Vector &c1 = theNodes[n1]->getCrds();
            const Vector &c2 = theNodes[n2]->getCrds();
            for (int i = 0; i < 3; i++) {
                v1(i) = c1(i);
                v2(i) = c2(i);
            }

            if (displayMode >= 0) {
                const Vector &d1 = theNodes[n1]->getDisp();
                const Vector &d2 = theNodes[n2]->getDisp();
                for (int i = 0; i < 3; i++) {
                    v1(i) += fact * d1(i);
                    v2(i) += fact * d2(i);
                }
            }

            theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
        }
    }

    return 0;
}


// ═══════════════════════════════════════════════════════════════════
//  RESPONSE AND RECORDING
// ═══════════════════════════════════════════════════════════════════

Response *BezierTet10::setResponse(const char **argv, int argc,
                                   OPS_Stream &output)
{
    Response *theResponse = 0;

    // ─── MPCO_Ladruno geometry self-declaration (contract Part A) ──
    //  topology=tet (simplex), family=bernstein, paramDomain=bary,
    //  rational=0 ⇒ no controlPointWeights. numCtrl=10, numGP=4.
    if (strcmp(argv[0], "basisInfo") == 0) {
        output.tag("ElementBasis");
        output.attr("topology",    "tet");
        output.attr("family",      "bernstein");
        output.attr("paramDomain", "bary");
        output.attr("rational",    0);
        output.attr("numCtrl",     NEN);      // 10 control points
        output.attr("numGP",       NGAUSS);   // 4 Gauss points (result stations)
        output.attr("orderU",      2);        // total polynomial degree (simplex)
        output.endTag();                      // ElementBasis
        return new ElementResponse(this, 101, ID(1));   // non-null sentinel
    }
    if (strcmp(argv[0], "integrationPoints") == 0)
        return new ElementResponse(this, 102, Matrix(NGAUSS, 3));
    if (strcmp(argv[0], "integrationWeights") == 0)
        return new ElementResponse(this, 103, Vector(NGAUSS));

    output.tag("ElementOutput");
    output.attr("eleType", "BezierTet10");
    output.attr("eleTag", this->getTag());
    for (int i = 0; i < NEN; i++) {
        char buf[16];
        sprintf(buf, "node%d", i + 1);
        output.attr(buf, connectedExternalNodes(i));
    }

    // ─── Element resisting forces ─────────────────────────────
    if (strcmp(argv[0], "force") == 0 ||
        strcmp(argv[0], "forces") == 0) {

        for (int i = 1; i <= NEN; i++) {
            char buf[32];
            sprintf(buf, "P1_%d", i); output.tag("ResponseType", buf);
            sprintf(buf, "P2_%d", i); output.tag("ResponseType", buf);
            sprintf(buf, "P3_%d", i); output.tag("ResponseType", buf);
        }
        theResponse = new ElementResponse(this, 1, Vector(NELD));
    }

    // ─── Stiffness matrix ─────────────────────────────────────
    else if (strcmp(argv[0], "stiff") == 0 ||
             strcmp(argv[0], "stiffness") == 0) {
        theResponse = new ElementResponse(this, 2, Matrix(NELD, NELD));
    }

    // ─── Material / integration point response ────────────────
    else if (strcmp(argv[0], "material") == 0 ||
             strcmp(argv[0], "integrPoint") == 0) {

        if (argc < 2) {
            opserr << "BezierTet10::setResponse - need GP number after 'material'\n";
            output.endTag();
            return 0;
        }

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= NGAUSS) {
            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("L1", GP4_L[pointNum-1][0]);
            output.attr("L2", GP4_L[pointNum-1][1]);
            output.attr("L3", GP4_L[pointNum-1][2]);

            double N[NEN];
            shapeFunctions(GP4_L[pointNum-1][0], GP4_L[pointNum-1][1],
                           GP4_L[pointNum-1][2], N);
            double x = 0.0, y = 0.0, z = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
                z += N[a] * controlPts[a][2];
            }
            output.attr("xLoc", x);
            output.attr("yLoc", y);
            output.attr("zLoc", z);

            theResponse = theMaterial[pointNum-1]->setResponse(
                &argv[2], argc-2, output);

            output.endTag();  // GaussPoint
        }
    }

    // ─── Stresses at ALL Gauss points ─────────────────────────
    else if (strcmp(argv[0], "stress") == 0 ||
             strcmp(argv[0], "stresses") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.tag("NdMaterialOutput");
            output.attr("classType", theMaterial[i]->getClassTag());
            output.attr("tag", theMaterial[i]->getTag());
            output.tag("ResponseType", "sigma_xx");
            output.tag("ResponseType", "sigma_yy");
            output.tag("ResponseType", "sigma_zz");
            output.tag("ResponseType", "sigma_xy");
            output.tag("ResponseType", "sigma_yz");
            output.tag("ResponseType", "sigma_zx");
            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }
        theResponse = new ElementResponse(this, 3, Vector(NSTRESS * NGAUSS));
    }

    // ─── Strains at ALL Gauss points ──────────────────────────
    else if (strcmp(argv[0], "strain") == 0 ||
             strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.tag("NdMaterialOutput");
            output.tag("ResponseType", "eps_xx");
            output.tag("ResponseType", "eps_yy");
            output.tag("ResponseType", "eps_zz");
            output.tag("ResponseType", "gamma_xy");
            output.tag("ResponseType", "gamma_yz");
            output.tag("ResponseType", "gamma_zx");
            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }
        theResponse = new ElementResponse(this, 4, Vector(NSTRESS * NGAUSS));
    }

    // ─── Gauss point physical coordinates ─────────────────────
    else if (strcmp(argv[0], "gaussPoint") == 0 ||
             strcmp(argv[0], "gaussCoord") == 0 ||
             strcmp(argv[0], "gpCoord") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("ResponseType", "x");
            output.tag("ResponseType", "y");
            output.tag("ResponseType", "z");
        }
        theResponse = new ElementResponse(this, 5, Vector(3 * NGAUSS));
    }

    // ─── Pressure contribution ────────────────────────────────
    else if (strcmp(argv[0], "pressure") == 0) {
        theResponse = new ElementResponse(this, 6, Vector(NELD));
    }

    // ─── Characteristic length ────────────────────────────────
    // Element-size lch used by crack-band materials (ASDConcrete). Lets the
    // user inspect the value the material regularizes with. See
    // getCharacteristicLength().
    else if (strcmp(argv[0], "charLength") == 0 ||
             strcmp(argv[0], "characteristicLength") == 0) {
        output.tag("ResponseType", "lch");
        theResponse = new ElementResponse(this, 7, Vector(1));
    }

    output.endTag();  // ElementOutput
    return theResponse;
}


int BezierTet10::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {

    case 1:  // resisting forces
        return eleInfo.setVector(this->getResistingForce());

    case 2:  // tangent stiffness
        return eleInfo.setMatrix(this->getTangentStiff());

    case 3: {
        // Stresses at all Gauss points (committed)
        static Vector stressVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &sigma = theMaterial[i]->getStress();
            for (int j = 0; j < NSTRESS; j++)
                stressVec(i * NSTRESS + j) = sigma(j);
        }
        return eleInfo.setVector(stressVec);
    }

    case 4: {
        // Strains at all Gauss points (γ engineering shear)
        static Vector strainVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &eps = theMaterial[i]->getStrain();
            for (int j = 0; j < NSTRESS; j++)
                strainVec(i * NSTRESS + j) = eps(j);
        }
        return eleInfo.setVector(strainVec);
    }

    case 5: {
        // Gauss point physical coordinates  x(ξ) = Σ Nₐ(ξ) Pₐ
        static Vector gpCoords(3 * NGAUSS);
        for (int gp = 0; gp < NGAUSS; gp++) {
            double N[NEN];
            shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);
            double x = 0.0, y = 0.0, z = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
                z += N[a] * controlPts[a][2];
            }
            gpCoords(3*gp)     = x;
            gpCoords(3*gp + 1) = y;
            gpCoords(3*gp + 2) = z;
        }
        return eleInfo.setVector(gpCoords);
    }

    case 6: {
        // Pressure contribution (volume hack, +z), 30-vector
        static Vector pressVec(NELD);
        pressVec.Zero();
        if (pressure != 0.0) {
            for (int gp = 0; gp < NGAUSS; gp++) {
                double w = GP4_w[gp];
                double N[NEN];
                shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);
                double dN[3][NEN];
                shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);
                double J[3][3];
                double dN_dx[3][NEN];
                double detJ = computeJacobian(dN, J, dN_dx);
                double factor = w * fabs(detJ) * pressure;
                for (int a = 0; a < NEN; a++)
                    pressVec(3*a + 2) += factor * N[a];
            }
        }
        return eleInfo.setVector(pressVec);
    }

    case 7: {
        // ─── Characteristic length ─────────────────────────────
        // The element-size lch crack-band materials regularize with.
        static Vector lchVec(1);
        lchVec(0) = this->getCharacteristicLength();
        return eleInfo.setVector(lchVec);
    }

    // ─── MPCO_Ladruno geometry probes (contract Part A) ──────────
    case 101:
        // basisInfo sentinel — metadata emitted via the stream already.
        return 0;

    case 102: {
        // GP_PARAM — barycentric coords (L1,L2,L3) of the 4 GPs; L4=1-ΣL.
        static Matrix L(NGAUSS, 3);
        for (int g = 0; g < NGAUSS; g++) {
            L(g, 0) = GP4_L[g][0];
            L(g, 1) = GP4_L[g][1];
            L(g, 2) = GP4_L[g][2];
        }
        return eleInfo.setMatrix(L);
    }

    case 103: {
        // GP_WEIGHT — weights, same GP order (sum = 1/6 = ref tet volume).
        static Vector w(NGAUSS);
        for (int g = 0; g < NGAUSS; g++)
            w(g) = GP4_w[g];
        return eleInfo.setVector(w);
    }

    default:
        return -1;
    }
}
