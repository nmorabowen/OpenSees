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
// Reference:
//   Kadapa, C. "Novel quadratic Bézier triangular and tetrahedral elements
//   using existing mesh generators: Applications to linear nearly
//   incompressible elastostatics and implicit and explicit elastodynamics."
//   Int. J. Numer. Methods Eng., 2019; 117(5):543-573. doi:10.1002/nme.5967
//
// ═══════════════════════════════════════════════════════════════════
//  MATHEMATICAL FOUNDATION
// ═══════════════════════════════════════════════════════════════════
//
//  Shape functions: Quadratic Bernstein polynomials on triangle
//  with area coordinates (ξ₁, ξ₂), ξ₃ = 1 - ξ₁ - ξ₂
//
//    N₁ = ξ₃²              (corner at ξ₃=1)
//    N₂ = ξ₁²              (corner at ξ₁=1)
//    N₃ = ξ₂²              (corner at ξ₂=1)
//    N₄ = 2ξ₁ξ₃            (mid-edge 1-2)
//    N₅ = 2ξ₁ξ₂            (mid-edge 2-3)
//    N₆ = 2ξ₂ξ₃            (mid-edge 3-1)
//
//  Key properties:
//    - Partition of unity: Σ Nₐ = 1
//    - Nonnegativity: Nₐ ≥ 0  ∀(ξ₁,ξ₂) in the triangle
//    - Equal integration: ∫ Nₐ dΩ = A/6 for all a
//
//  These ensure all-positive lumped mass for explicit dynamics.
//
//  Isoparametric mapping uses CONTROL POINTS, not Lagrange nodes:
//    x(ξ) = Σ Nₐ(ξ) · Pₐ
//
//  For straight-sided elements: Pₐ = Xₐ (control points = nodes).
//  For curved edges: Pₐ are computed from Lagrange nodes via Eq.(13).
//
// ═══════════════════════════════════════════════════════════════════

#include "BezierTri6.h"

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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// ═══════════════════════════════════════════════════════════════════
//  STATIC DATA
// ═══════════════════════════════════════════════════════════════════

// Instance counter for one-time attribution (Abell pattern)
int BezierTri6::numInstances = 0;

// Static return objects — class-wide scratch storage (Petracca pattern)
// Avoids heap allocation on every getTangentStiff()/getMass() call
Matrix BezierTri6::K_return(NELD, NELD);
Matrix BezierTri6::M_return(NELD, NELD);
Vector BezierTri6::P_return(NELD);

// ─── 3-Point Gauss Quadrature (degree 2) ──────────────────────
// Exact for quadratic integrands → sufficient for stiffness (B^T D B
// has a rational integrand, but degree 2 gives good accuracy for
// straight-sided elements where the Jacobian is constant).
//
// Weights sum to 1/2 (area of reference triangle).
const double BezierTri6::GP3_xi[][2] = {
    {1.0/6.0, 1.0/6.0},
    {2.0/3.0, 1.0/6.0},
    {1.0/6.0, 2.0/3.0}
};
const double BezierTri6::GP3_w[] = {
    1.0/6.0, 1.0/6.0, 1.0/6.0
};

// ─── 6-Point Gauss Quadrature (degree 4) ──────────────────────
// Exact for quartic integrands → exact consistent mass matrix
// (product of two quadratic Bernstein polynomials = degree 4).
//
// Dunavant's rule. Raw weights sum to 1.0; we store them divided
// by 2 so they sum to 1/2 (reference triangle area).
const double BezierTri6::GP6_xi[][2] = {
    {0.445948490915965, 0.445948490915965},
    {0.108103018168070, 0.445948490915965},
    {0.445948490915965, 0.108103018168070},
    {0.091576213509771, 0.091576213509771},
    {0.816847572980459, 0.091576213509771},
    {0.091576213509771, 0.816847572980459}
};
const double BezierTri6::GP6_w[] = {
    0.223381589678011 / 2.0,
    0.223381589678011 / 2.0,
    0.223381589678011 / 2.0,
    0.109951743655322 / 2.0,
    0.109951743655322 / 2.0,
    0.109951743655322 / 2.0
};


// ═══════════════════════════════════════════════════════════════════
//  CONSTRUCTORS AND DESTRUCTOR
// ═══════════════════════════════════════════════════════════════════

BezierTri6::BezierTri6(int tag,
                       int nd1, int nd2, int nd3,
                       int nd4, int nd5, int nd6,
                       NDMaterial &m, const char *type, double thick,
                       double press, double r,
                       double bx, double by,
                       bool bbar, bool cmass)
  : Element(tag, ELE_TAG_BezierTri6),
    theMaterial(0),
    connectedExternalNodes(NEN),
    thickness(thick), rho(r), pressure(press),
    useBbar(bbar), cMass(cmass),
    Q(NELD), applyLoad(0), Ki(0)
{
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
        opserr << "BezierTri6 - Quadratic Bezier Triangle\n"
               << "  Authors: Nicolas Mora Bowen, Patricio Palacios, "
                  "Jose Abell, Guppi (Ladruño)\n"
               << "  Ref: Kadapa, IJNME 2019; 117(5):543-573. "
                  "doi:10.1002/nme.5967\n";
    }

    // Store body forces
    b[0] = bx;
    b[1] = by;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;

    // Store connectivity
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    connectedExternalNodes(4) = nd5;
    connectedExternalNodes(5) = nd6;

    // Initialize node pointers
    for (int i = 0; i < NEN; i++)
        theNodes[i] = 0;

    // Initialize control points
    for (int i = 0; i < NEN; i++) {
        controlPts[i][0] = 0.0;
        controlPts[i][1] = 0.0;
    }

    // Zero initial stiffness data
    for (int i = 0; i < NELD * NELD; i++)
        Ki_data[i] = 0.0;

    // Allocate materials at Gauss points
    // Use NGAUSS points for stiffness/stress, NGAUSS_MASS is only
    // used in getMass() where we don't store materials.
    theMaterial = new NDMaterial*[NGAUSS];
    for (int i = 0; i < NGAUSS; i++) {
        theMaterial[i] = m.getCopy(type);
        if (theMaterial[i] == 0) {
            opserr << "BezierTri6::BezierTri6 - failed to get copy of material\n";
            exit(-1);
        }
    }
}


BezierTri6::BezierTri6()
  : Element(0, ELE_TAG_BezierTri6),
    theMaterial(0),
    connectedExternalNodes(NEN),
    thickness(0.0), rho(0.0), pressure(0.0),
    useBbar(false), cMass(false),
    Q(NELD), applyLoad(0), Ki(0)
{
    b[0] = 0.0;
    b[1] = 0.0;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = 0;
        controlPts[i][0] = 0.0;
        controlPts[i][1] = 0.0;
    }

    for (int i = 0; i < NELD * NELD; i++)
        Ki_data[i] = 0.0;
}


BezierTri6::~BezierTri6()
{
    if (theMaterial != 0) {
        for (int i = 0; i < NGAUSS; i++) {
            if (theMaterial[i])
                delete theMaterial[i];
        }
        delete[] theMaterial;
    }

    if (Ki != 0)
        delete Ki;
}


// ═══════════════════════════════════════════════════════════════════
//  INFORMATION METHODS
// ═══════════════════════════════════════════════════════════════════

int BezierTri6::getNumExternalNodes(void) const
{
    return NEN;
}

const ID &BezierTri6::getExternalNodes(void)
{
    return connectedExternalNodes;
}

Node **BezierTri6::getNodePtrs(void)
{
    return theNodes;
}

int BezierTri6::getNumDOF(void)
{
    return NELD;
}


// ═══════════════════════════════════════════════════════════════════
//  LIFECYCLE METHODS
// ═══════════════════════════════════════════════════════════════════

void BezierTri6::setDomain(Domain *theDomain)
{
    // Check domain
    if (theDomain == 0) {
        for (int i = 0; i < NEN; i++)
            theNodes[i] = 0;
        return;
    }

    // Get node pointers
    for (int i = 0; i < NEN; i++) {
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
        if (theNodes[i] == 0) {
            opserr << "BezierTri6::setDomain -- node "
                   << connectedExternalNodes(i) << " does not exist\n";
            return;
        }
        if (theNodes[i]->getNumberDOF() != NDOF) {
            opserr << "BezierTri6::setDomain -- node "
                   << connectedExternalNodes(i)
                   << " has " << theNodes[i]->getNumberDOF()
                   << " DOFs, expected " << NDOF << "\n";
            return;
        }
    }

    // Compute Bézier control points from Lagrange node positions
    computeControlPoints();

    // ─── v1 straight-sided guard (ADR 04_bezier_elements D9) ──────
    //  v1 is validated for STRAIGHT-SIDED elements only: mid-edge nodes at
    //  the edge midpoints ⇒ control points = nodes (P = X), so Dirichlet BCs
    //  stay interpolatory and 3-pt quadrature is exact. A curved mesh silently
    //  breaks both (a prescribed displacement lands on a non-interpolatory
    //  control point; B^T D B is under-integrated). Warn so the user knows.
    //  Curved support + the Eq. 14 control-point↔DOF mapping is a follow-up.
    {
        const int corner[3][2] = {{0,1}, {1,2}, {2,0}};  // edge corners for mids 4,5,6
        for (int e = 0; e < 3; e++) {
            const Vector &Xc1 = theNodes[corner[e][0]]->getCrds();
            const Vector &Xc2 = theNodes[corner[e][1]]->getCrds();
            const Vector &Xm  = theNodes[3 + e]->getCrds();
            double L2 = 0.0, d2 = 0.0;
            for (int k = 0; k < 2; k++) {
                double mid = 0.5 * (Xc1(k) + Xc2(k));
                d2 += (Xm(k) - mid) * (Xm(k) - mid);
                L2 += (Xc2(k) - Xc1(k)) * (Xc2(k) - Xc1(k));
            }
            if (L2 > 0.0 && d2 > 1.0e-12 * L2) {   // > ~1e-6 relative deviation
                opserr << "WARNING BezierTri6 " << this->getTag()
                       << " - mid-edge node " << connectedExternalNodes(3 + e)
                       << " is off the edge midpoint (curved element). v1 is "
                          "validated for straight-sided meshes only; Dirichlet "
                          "BCs and quadrature may be inaccurate.\n";
            }
        }
    }

    // IMPORTANT: call base class setDomain
    this->DomainComponent::setDomain(theDomain);
}


int BezierTri6::commitState()
{
    int retVal = 0;

    // Call Element commitState (handles Rayleigh damping, etc.)
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "BezierTri6::commitState() - failed in base class\n";
    }

    // Commit materials at all Gauss points
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->commitState();

    return retVal;
}


int BezierTri6::revertToLastCommit()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToLastCommit();
    return retVal;
}


int BezierTri6::revertToStart()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToStart();

    // Recompute control points (in case nodes moved via setDomain)
    if (theNodes[0] != 0)
        computeControlPoints();

    return retVal;
}


int BezierTri6::update()
{
    // ═══════════════════════════════════════════════════════════
    //  update() — Set trial strains at each Gauss point
    //
    //  This is the core element computation. The framework has
    //  updated nodal displacements; we extract them, compute
    //  strains via the B matrix (or B-bar), and push them to
    //  the material objects.
    //
    //  ε = B · u   (or ε̄ = B̄ · u for B-bar)
    // ═══════════════════════════════════════════════════════════

    // Extract current nodal displacements
    Vector u(NELD);
    for (int i = 0; i < NEN; i++) {
        const Vector &disp = theNodes[i]->getTrialDisp();
        u(2*i)     = disp(0);
        u(2*i + 1) = disp(1);
    }

    int ret = 0;

    // For B-bar: precompute volume-averaged derivatives
    double dN_avg[2][NEN];
    double volume;
    if (useBbar) {
        computeVolumeAveragedDerivatives(dN_avg, volume);
    }

    // Loop over Gauss points
    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];

        // Compute shape function derivatives
        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        // Compute Jacobian and physical derivatives
        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        if (detJ <= 0.0) {
            opserr << "BezierTri6::update() - non-positive Jacobian\n";
            return -1;
        }

        // Compute B matrix
        double B[NSTRESS][NELD];
        if (useBbar) {
            computeBBarMatrix(dN_dx, dN_avg, B);
        } else {
            computeBMatrix(dN_dx, B);
        }

        // Compute strain: ε = B · u
        //   ε = {ε_xx, ε_yy, γ_xy}  (engineering shear strain)
        Vector strain(NSTRESS);
        for (int i = 0; i < NSTRESS; i++) {
            double sum = 0.0;
            for (int j = 0; j < NELD; j++)
                sum += B[i][j] * u(j);
            strain(i) = sum;
        }

        // Set trial strain on material
        ret += theMaterial[gp]->setTrialStrain(strain);
    }

    return ret;
}


// ═══════════════════════════════════════════════════════════════════
//  STIFFNESS AND FORCE
// ═══════════════════════════════════════════════════════════════════

const Matrix &BezierTri6::getTangentStiff()
{
    // ═══════════════════════════════════════════════════════════
    //  K = ∫_Ω Bᵀ D_tan B dΩ   (tangent stiffness)
    //
    //  Evaluated with NGAUSS-point quadrature:
    //    K = Σᵢ wᵢ |Jᵢ| t · Bᵢᵀ Dᵢ Bᵢ
    //
    //  For B-bar: replace B with B̄
    // ═══════════════════════════════════════════════════════════

    K_return.Zero();

    // B-bar: precompute volume-averaged derivatives
    double dN_avg[2][NEN];
    double volume;
    if (useBbar) {
        computeVolumeAveragedDerivatives(dN_avg, volume);
    }

    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];
        double w = GP3_w[gp];

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        double B[NSTRESS][NELD];
        if (useBbar)
            computeBBarMatrix(dN_dx, dN_avg, B);
        else
            computeBMatrix(dN_dx, B);

        // Get material tangent at this Gauss point
        const Matrix &D = theMaterial[gp]->getTangent();

        // Factor: w × |J| × thickness
        double factor = w * detJ * thickness;

        // K += factor × Bᵀ D B
        // Using explicit loops for clarity and to avoid temporary allocation
        for (int I = 0; I < NELD; I++) {
            for (int J = I; J < NELD; J++) {  // Exploit symmetry
                double sum = 0.0;
                for (int k = 0; k < NSTRESS; k++) {
                    for (int l = 0; l < NSTRESS; l++) {
                        sum += B[k][I] * D(k, l) * B[l][J];
                    }
                }
                K_return(I, J) += factor * sum;
                if (I != J)
                    K_return(J, I) += factor * sum;  // Symmetric
            }
        }
    }

    return K_return;
}


const Matrix &BezierTri6::getInitialStiff()
{
    if (Ki != 0)
        return *Ki;

    // Compute initial stiffness using initial material tangents
    K_return.Zero();

    double dN_avg[2][NEN];
    double volume;
    if (useBbar) {
        computeVolumeAveragedDerivatives(dN_avg, volume);
    }

    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];
        double w = GP3_w[gp];

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        double B[NSTRESS][NELD];
        if (useBbar)
            computeBBarMatrix(dN_dx, dN_avg, B);
        else
            computeBMatrix(dN_dx, B);

        const Matrix &D = theMaterial[gp]->getInitialTangent();
        double factor = w * detJ * thickness;

        for (int I = 0; I < NELD; I++) {
            for (int J = I; J < NELD; J++) {
                double sum = 0.0;
                for (int k = 0; k < NSTRESS; k++)
                    for (int l = 0; l < NSTRESS; l++)
                        sum += B[k][I] * D(k, l) * B[l][J];
                K_return(I, J) += factor * sum;
                if (I != J)
                    K_return(J, I) += factor * sum;
            }
        }
    }

    // Cache the initial stiffness
    Ki = new Matrix(Ki_data, NELD, NELD);
    *Ki = K_return;

    return *Ki;
}


const Matrix &BezierTri6::getMass()
{
    // ═══════════════════════════════════════════════════════════
    //  MASS MATRIX — lumped (default) or consistent (-cMass)
    //
    //  Density: the element -rho overrides; otherwise the material's
    //  own density (NDMaterial::getRho) is used, matching FourNodeQuad.
    // ═══════════════════════════════════════════════════════════

    M_return.Zero();

    // Effective density (element rho overrides material rho)
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
        // ─── Lumped mass (Kadapa Eq. 56) ──────────────────────
        //   Mᵉ = (ρ Vₑ / 6) · diag[1₆, 1₆]
        //
        //  Each of the 6 control points gets an EQUAL share ρVₑ/6.
        //  This is exact for the quadratic Bernstein basis (∫Nₐ = Vₑ/6
        //  for every a) and — crucially — ALL POSITIVE, which is the
        //  whole reason this element exists for explicit dynamics.
        double Ve = this->computeArea() * thickness;   // element volume
        double m  = rhoEff * Ve / 6.0;                 // mass per node
        for (int a = 0; a < NEN; a++) {
            M_return(2*a,   2*a)   = m;
            M_return(2*a+1, 2*a+1) = m;
        }
        return M_return;
    }

    // ─── Consistent mass:  M = ∫_Ω ρ Nᵤᵀ Nᵤ dΩ ────────────────
    //  6-point quadrature (degree 4) integrates the quartic integrand
    //  (product of two quadratic Bernsteins) exactly for straight edges.
    //  All entries ≥ 0 because every Nₐ ≥ 0.
    for (int gp = 0; gp < NGAUSS_MASS; gp++) {
        double xi1 = GP6_xi[gp][0];
        double xi2 = GP6_xi[gp][1];
        double w = GP6_w[gp];

        double N[NEN];
        shapeFunctions(xi1, xi2, N);

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        // Factor: w × |J| × thickness × density
        double factor = w * detJ * thickness * rhoEff;

        // M_{(2a+p)(2b+p)} += factor × Nₐ × Nᵦ
        for (int a = 0; a < NEN; a++) {
            for (int bb = 0; bb < NEN; bb++) {
                double val = factor * N[a] * N[bb];
                M_return(2*a,   2*bb)   += val;  // x-x coupling
                M_return(2*a+1, 2*bb+1) += val;  // y-y coupling
            }
        }
    }

    return M_return;
}


void BezierTri6::zeroLoad()
{
    Q.Zero();
    applyLoad = 0;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
}


int BezierTri6::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // Supports the standard self-weight (gravity) load so body force can
    // be driven by a load pattern / time series instead of being a fixed
    // constructor constant. SelfWeight's data scales the element's base
    // body force b[] (same idiom as FourNodeQuad).
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor * data(0) * b[0];
        appliedB[1] += loadFactor * data(1) * b[1];
        return 0;
    }

    opserr << "BezierTri6::addLoad - load type " << type
           << " unknown for element " << this->getTag()
           << " (only SelfWeight is supported)\n";
    return -1;
}


int BezierTri6::addInertiaLoadToUnbalance(const Vector &accel)
{
    // Mass may come from the element rho OR the material density, so
    // probe the assembled mass matrix rather than testing rho directly.
    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;
    if (!hasMass)
        return 0;

    // Get accelerations at nodes
    static Vector a(NELD);
    for (int i = 0; i < NEN; i++) {
        const Vector &Raccel = theNodes[i]->getRV(accel);
        if (Raccel.Size() != NDOF) {
            opserr << "BezierTri6::addInertiaLoadToUnbalance - "
                   << "matrix and target sizes mismatch\n";
            return -1;
        }
        a(2*i)     = Raccel(0);
        a(2*i + 1) = Raccel(1);
    }

    // Q += M × a
    Q.addMatrixVector(1.0, M, a, 1.0);

    return 0;
}


const Vector &BezierTri6::getResistingForce()
{
    // ═══════════════════════════════════════════════════════════
    //  F^int = ∫_Ω Bᵀ σ dΩ  (internal resisting force)
    //
    //  Also subtracts applied body forces and pressure.
    // ═══════════════════════════════════════════════════════════

    P_return.Zero();

    // B-bar precomputation
    double dN_avg[2][NEN];
    double volume;
    if (useBbar) {
        computeVolumeAveragedDerivatives(dN_avg, volume);
    }

    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];
        double w = GP3_w[gp];

        double N[NEN];
        shapeFunctions(xi1, xi2, N);

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J_mat[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J_mat, dN_dx);

        double B[NSTRESS][NELD];
        if (useBbar)
            computeBBarMatrix(dN_dx, dN_avg, B);
        else
            computeBMatrix(dN_dx, B);

        // Get stress from material
        const Vector &sigma = theMaterial[gp]->getStress();

        double factor = w * detJ * thickness;

        // F^int += factor × Bᵀ σ
        for (int I = 0; I < NELD; I++) {
            double sum = 0.0;
            for (int k = 0; k < NSTRESS; k++)
                sum += B[k][I] * sigma(k);
            P_return(I) += factor * sum;
        }

        // ─── Subtract body forces ─────────────────────────────
        // F^body = -∫ Nᵤᵀ f dΩ  (negative sign: goes to RHS).
        // When a SelfWeight load pattern is active, use the scaled
        // appliedB (rampable); otherwise fall back to the constant b[].
        double bx = (applyLoad == 1) ? appliedB[0] : b[0];
        double by = (applyLoad == 1) ? appliedB[1] : b[1];
        if (bx != 0.0 || by != 0.0) {
            for (int a = 0; a < NEN; a++) {
                P_return(2*a)     -= factor * N[a] * bx;
                P_return(2*a + 1) -= factor * N[a] * by;
            }
        }

        // ─── Subtract pressure load ──────────────────────────
        // Pressure acts as distributed load in the y-direction
        // F^press = -∫ Nᵤᵀ {0, p} dΩ · t
        if (pressure != 0.0) {
            for (int a = 0; a < NEN; a++) {
                P_return(2*a + 1) -= factor * N[a] * pressure;
            }
        }
    }

    // Subtract any applied element loads (from addLoad)
    P_return.addVector(1.0, Q, -1.0);

    return P_return;
}


const Vector &BezierTri6::getResistingForceIncInertia()
{
    // Get static resisting force (fills P_return)
    this->getResistingForce();

    // ─── Add inertia: R += M × a ──────────────────────────────
    // Mass may come from the element rho OR the material density, so
    // probe the assembled mass rather than testing rho directly.
    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;

    if (hasMass) {
        static Vector a(NELD);
        for (int i = 0; i < NEN; i++) {
            const Vector &accel = theNodes[i]->getTrialAccel();
            a(2*i)     = accel(0);
            a(2*i + 1) = accel(1);
        }
        P_return.addMatrixVector(1.0, M, a, 1.0);
    }

    // ─── Add Rayleigh damping if present ──────────────────────
    // Independent of mass: stiffness-proportional damping (betaK) must
    // still be added even when there is no mass.
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
        const Vector &v = this->getRayleighDampingForces();
        P_return += v;
    }

    return P_return;
}


// ═══════════════════════════════════════════════════════════════════
//  SHAPE FUNCTIONS AND DERIVATIVES
// ═══════════════════════════════════════════════════════════════════

void BezierTri6::shapeFunctions(double xi1, double xi2,
                                 double N[NEN]) const
{
    //  Quadratic Bernstein polynomials on a triangle.
    //
    //  Area coordinates: ξ₁, ξ₂, ξ₃ = 1-ξ₁-ξ₂
    //
    //  All Nₐ ≥ 0 — the defining property that makes this element
    //  special for explicit dynamics.

    double xi3 = 1.0 - xi1 - xi2;

    N[0] = xi3 * xi3;          // ξ₃²
    N[1] = xi1 * xi1;          // ξ₁²
    N[2] = xi2 * xi2;          // ξ₂²
    N[3] = 2.0 * xi1 * xi3;   // 2ξ₁ξ₃
    N[4] = 2.0 * xi1 * xi2;   // 2ξ₁ξ₂
    N[5] = 2.0 * xi2 * xi3;   // 2ξ₂ξ₃
}


void BezierTri6::shapeDerivatives(double xi1, double xi2,
                                   double dN[2][NEN]) const
{
    //  ∂Nₐ/∂ξ₁ and ∂Nₐ/∂ξ₂
    //
    //  Using ξ₃ = 1-ξ₁-ξ₂, so ∂ξ₃/∂ξ₁ = ∂ξ₃/∂ξ₂ = -1.

    double xi3 = 1.0 - xi1 - xi2;

    // ∂N/∂ξ₁
    dN[0][0] = -2.0 * xi3;           // ∂(ξ₃²)/∂ξ₁ = -2ξ₃
    dN[0][1] =  2.0 * xi1;           // ∂(ξ₁²)/∂ξ₁ = 2ξ₁
    dN[0][2] =  0.0;                 // ∂(ξ₂²)/∂ξ₁ = 0
    dN[0][3] =  2.0 * (xi3 - xi1);  // ∂(2ξ₁ξ₃)/∂ξ₁ = 2(ξ₃-ξ₁)
    dN[0][4] =  2.0 * xi2;           // ∂(2ξ₁ξ₂)/∂ξ₁ = 2ξ₂
    dN[0][5] = -2.0 * xi2;           // ∂(2ξ₂ξ₃)/∂ξ₁ = -2ξ₂

    // ∂N/∂ξ₂
    dN[1][0] = -2.0 * xi3;           // ∂(ξ₃²)/∂ξ₂ = -2ξ₃
    dN[1][1] =  0.0;                 // ∂(ξ₁²)/∂ξ₂ = 0
    dN[1][2] =  2.0 * xi2;           // ∂(ξ₂²)/∂ξ₂ = 2ξ₂
    dN[1][3] = -2.0 * xi1;           // ∂(2ξ₁ξ₃)/∂ξ₂ = -2ξ₁
    dN[1][4] =  2.0 * xi1;           // ∂(2ξ₁ξ₂)/∂ξ₂ = 2ξ₁
    dN[1][5] =  2.0 * (xi3 - xi2);  // ∂(2ξ₂ξ₃)/∂ξ₂ = 2(ξ₃-ξ₂)
}


// ═══════════════════════════════════════════════════════════════════
//  JACOBIAN AND B-MATRIX
// ═══════════════════════════════════════════════════════════════════

double BezierTri6::computeJacobian(const double dN[2][NEN],
                                    double J[2][2],
                                    double dN_dx[2][NEN]) const
{
    //  Jacobian of isoparametric mapping: x(ξ) = Σ Nₐ(ξ) · Pₐ
    //
    //  J = [∂x/∂ξ₁  ∂y/∂ξ₁]  = Σₐ [∂Nₐ/∂ξ₁·xₐ  ∂Nₐ/∂ξ₁·yₐ]
    //      [∂x/∂ξ₂  ∂y/∂ξ₂]      [∂Nₐ/∂ξ₂·xₐ  ∂Nₐ/∂ξ₂·yₐ]
    //
    //  Uses CONTROL POINTS (not Lagrange nodes!) for the mapping.

    J[0][0] = 0.0; J[0][1] = 0.0;
    J[1][0] = 0.0; J[1][1] = 0.0;

    for (int a = 0; a < NEN; a++) {
        J[0][0] += dN[0][a] * controlPts[a][0];
        J[0][1] += dN[0][a] * controlPts[a][1];
        J[1][0] += dN[1][a] * controlPts[a][0];
        J[1][1] += dN[1][a] * controlPts[a][1];
    }

    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    if (detJ > 0.0) {
        // Compute J⁻¹ and transform derivatives to physical coords
        //   [∂N/∂x] = J⁻¹ [∂N/∂ξ₁]
        //   [∂N/∂y]       [∂N/∂ξ₂]
        double invDetJ = 1.0 / detJ;
        double Jinv[2][2];
        Jinv[0][0] =  J[1][1] * invDetJ;
        Jinv[0][1] = -J[0][1] * invDetJ;
        Jinv[1][0] = -J[1][0] * invDetJ;
        Jinv[1][1] =  J[0][0] * invDetJ;

        for (int a = 0; a < NEN; a++) {
            dN_dx[0][a] = Jinv[0][0] * dN[0][a] + Jinv[0][1] * dN[1][a];
            dN_dx[1][a] = Jinv[1][0] * dN[0][a] + Jinv[1][1] * dN[1][a];
        }
    }

    return detJ;
}


void BezierTri6::computeBMatrix(const double dN_dx[2][NEN],
                                 double B[NSTRESS][NELD]) const
{
    //  Standard strain-displacement matrix for 2D:
    //
    //  B_a = [∂Nₐ/∂x    0     ]
    //        [  0      ∂Nₐ/∂y  ]
    //        [∂Nₐ/∂y  ∂Nₐ/∂x  ]
    //
    //  Voigt ordering: {ε_xx, ε_yy, γ_xy} with engineering shear strain

    // Zero the matrix
    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            B[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        int col_x = 2 * a;
        int col_y = 2 * a + 1;

        B[0][col_x] = dN_dx[0][a];   // ε_xx = ∂u/∂x
        B[1][col_y] = dN_dx[1][a];   // ε_yy = ∂v/∂y
        B[2][col_x] = dN_dx[1][a];   // γ_xy = ∂u/∂y + ∂v/∂x
        B[2][col_y] = dN_dx[0][a];
    }
}


void BezierTri6::computeBBarMatrix(const double dN_dx[2][NEN],
                                    const double dN_avg[2][NEN],
                                    double Bbar[NSTRESS][NELD]) const
{
    //  B-bar formulation (Kadapa Eq. 45, adapted for 2D):
    //
    //  Replace the volumetric (dilatational) part of B with its
    //  volume-weighted average. The deviatoric part stays local.
    //
    //  For node a, column for u_x DOF:
    //    row 0 (ε_xx): (B̄₁ + 2B₁)/3
    //    row 1 (ε_yy): (B̄₁ - B₁)/3
    //    row 2 (γ_xy): B₂             (unchanged)
    //
    //  For node a, column for u_y DOF:
    //    row 0 (ε_xx): (B̄₂ - B₂)/3
    //    row 1 (ε_yy): (B̄₂ + 2B₂)/3
    //    row 2 (γ_xy): B₁             (unchanged)
    //
    //  where B₁ = ∂Nₐ/∂x, B₂ = ∂Nₐ/∂y at the current GP
    //        B̄₁, B̄₂ = volume-averaged derivatives

    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            Bbar[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        double B1 = dN_dx[0][a];       // ∂Nₐ/∂x at this GP
        double B2 = dN_dx[1][a];       // ∂Nₐ/∂y at this GP
        double Bbar1 = dN_avg[0][a];   // volume-averaged ∂Nₐ/∂x
        double Bbar2 = dN_avg[1][a];   // volume-averaged ∂Nₐ/∂y

        int col_x = 2 * a;
        int col_y = 2 * a + 1;

        // u_x DOF column
        Bbar[0][col_x] = (Bbar1 + 2.0 * B1) / 3.0;
        Bbar[1][col_x] = (Bbar1 - B1) / 3.0;
        Bbar[2][col_x] = B2;

        // u_y DOF column
        Bbar[0][col_y] = (Bbar2 - B2) / 3.0;
        Bbar[1][col_y] = (Bbar2 + 2.0 * B2) / 3.0;
        Bbar[2][col_y] = B1;
    }
}


void BezierTri6::computeVolumeAveragedDerivatives(double dN_avg[2][NEN],
                                                    double &volume) const
{
    //  B̄ₐ = (1/V) ∫_Ω ∂Nₐ/∂xᵢ dΩ
    //
    //  Computed with 3-point quadrature.

    for (int i = 0; i < 2; i++)
        for (int a = 0; a < NEN; a++)
            dN_avg[i][a] = 0.0;

    volume = 0.0;

    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];
        double w = GP3_w[gp];

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        double factor = w * detJ * thickness;
        volume += factor;

        for (int a = 0; a < NEN; a++) {
            dN_avg[0][a] += factor * dN_dx[0][a];
            dN_avg[1][a] += factor * dN_dx[1][a];
        }
    }

    // Normalize by volume
    if (volume > 0.0) {
        double invVol = 1.0 / volume;
        for (int a = 0; a < NEN; a++) {
            dN_avg[0][a] *= invVol;
            dN_avg[1][a] *= invVol;
        }
    }
}


// ═══════════════════════════════════════════════════════════════════
//  LAGRANGE → BÉZIER CONTROL POINT MAPPING
// ═══════════════════════════════════════════════════════════════════

void BezierTri6::computeControlPoints()
{
    //  From Kadapa Eq. (13):
    //
    //  Corner control points = corner nodes (interpolatory):
    //    P₁ = X₁,  P₂ = X₂,  P₃ = X₃
    //
    //  Mid-edge control points (generally NOT on the edge):
    //    P₄ = 2[X₄ - 0.25·X₁ - 0.25·X₂]
    //    P₅ = 2[X₅ - 0.25·X₂ - 0.25·X₃]
    //    P₆ = 2[X₆ - 0.25·X₃ - 0.25·X₁]
    //
    //  For straight-sided elements (mid-nodes at midpoints):
    //    P₄ = X₄, P₅ = X₅, P₆ = X₆  (trivially)

    // Get current node coordinates
    double X[NEN][2];
    for (int i = 0; i < NEN; i++) {
        const Vector &crds = theNodes[i]->getCrds();
        X[i][0] = crds(0);
        X[i][1] = crds(1);
    }

    // Corners: control point = node
    for (int d = 0; d < 2; d++) {
        controlPts[0][d] = X[0][d];
        controlPts[1][d] = X[1][d];
        controlPts[2][d] = X[2][d];
    }

    // Mid-edge: Eq. (13d-f)
    for (int d = 0; d < 2; d++) {
        controlPts[3][d] = 2.0 * (X[3][d] - 0.25 * X[0][d] - 0.25 * X[1][d]);
        controlPts[4][d] = 2.0 * (X[4][d] - 0.25 * X[1][d] - 0.25 * X[2][d]);
        controlPts[5][d] = 2.0 * (X[5][d] - 0.25 * X[2][d] - 0.25 * X[0][d]);
    }
}


double BezierTri6::computeArea() const
{
    double area = 0.0;

    for (int gp = 0; gp < NGAUSS; gp++) {
        double xi1 = GP3_xi[gp][0];
        double xi2 = GP3_xi[gp][1];
        double w = GP3_w[gp];

        double dN[2][NEN];
        shapeDerivatives(xi1, xi2, dN);

        double J[2][2];
        double dN_dx[2][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        area += w * detJ;
    }

    return area;
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
//  We instead return an element-size equivalent from the integrated area:
//  the leg of a right-isosceles triangle of equal area,
//
//      lch = sqrt(2 · A),
//
//  which recovers the true edge length for the common right-triangle mesh
//  (a split quad) and is geometry-true for distorted straight-sided elements.
//  Curved Bézier edges shift the factor slightly, but always in the safe
//  (under-estimating) direction. Area only — lch is an in-plane size, so the
//  out-of-plane thickness is deliberately excluded.
double BezierTri6::getCharacteristicLength(void)
{
    double A = this->computeArea();
    if (A <= 0.0)
        return Element::getCharacteristicLength();  // degenerate: fall back
    return sqrt(2.0 * A);
}


// ═══════════════════════════════════════════════════════════════════
//  SERIALIZATION (sendSelf / recvSelf)
// ═══════════════════════════════════════════════════════════════════

int BezierTri6::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // Send integer data
    static ID iData(11);
    iData(0) = this->getTag();
    for (int i = 0; i < NEN; i++)
        iData(i + 1) = connectedExternalNodes(i);
    iData(7) = theMaterial[0]->getClassTag();
    int matDbTag = theMaterial[0]->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->setDbTag(matDbTag);
    }
    iData(8) = matDbTag;
    iData(9) = useBbar ? 1 : 0;
    iData(10) = cMass ? 1 : 0;

    res += theChannel.sendID(this->getDbTag(), commitTag, iData);

    // Send double data
    static Vector dData(6);
    dData(0) = thickness;
    dData(1) = rho;
    dData(2) = pressure;
    dData(3) = b[0];
    dData(4) = b[1];
    dData(5) = 0.0;  // reserved

    res += theChannel.sendVector(this->getDbTag(), commitTag, dData);

    // Send material states
    for (int i = 0; i < NGAUSS; i++)
        res += theMaterial[i]->sendSelf(commitTag, theChannel);

    return res;
}


int BezierTri6::recvSelf(int commitTag, Channel &theChannel,
                          FEM_ObjectBroker &theBroker)
{
    int res = 0;

    // Receive integer data
    static ID iData(11);
    res += theChannel.recvID(this->getDbTag(), commitTag, iData);

    this->setTag(iData(0));
    for (int i = 0; i < NEN; i++)
        connectedExternalNodes(i) = iData(i + 1);
    int matClassTag = iData(7);
    int matDbTag = iData(8);
    useBbar = (iData(9) == 1);
    cMass = (iData(10) == 1);

    // Receive double data
    static Vector dData(6);
    res += theChannel.recvVector(this->getDbTag(), commitTag, dData);

    thickness = dData(0);
    rho = dData(1);
    pressure = dData(2);
    b[0] = dData(3);
    b[1] = dData(4);

    // Allocate materials if needed
    if (theMaterial == 0) {
        theMaterial = new NDMaterial*[NGAUSS];
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i] = 0;
    }

    // Receive material states
    for (int i = 0; i < NGAUSS; i++) {
        if (theMaterial[i] == 0) {
            theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
            if (theMaterial[i] == 0) {
                opserr << "BezierTri6::recvSelf - material creation failed\n";
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

void BezierTri6::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        // JSON-like output
        s << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTri6\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "]}";
        return;
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nBezierTri6, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tthickness:  " << thickness << endln;
        s << "\tmass density:  " << rho << endln;
        s << "\tB-bar formulation: " << (useBbar ? "yes" : "no") << endln;
        s << "\tbody forces:  " << b[0] << " " << b[1] << endln;

        s << "\tStress (Gauss):" << endln;
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTri6\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "], ";
        s << "\"thickness\": " << thickness << ", ";
        s << "\"bbar\": " << (useBbar ? "true" : "false") << ", ";
        s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
    }
}


int BezierTri6::displaySelf(Renderer &theViewer, int displayMode,
                             float fact, const char **modes, int numMode)
{
    // Get end points for each edge and draw
    // For a 6-node triangle, draw the 3 edges with mid-nodes

    static Vector v1(3), v2(3);

    // Draw edges: 0→3→1, 1→4→2, 2→5→0
    int edges[3][3] = {{0,3,1}, {1,4,2}, {2,5,0}};

    for (int e = 0; e < 3; e++) {
        for (int seg = 0; seg < 2; seg++) {
            int n1 = edges[e][seg];
            int n2 = edges[e][seg+1];

            const Vector &c1 = theNodes[n1]->getCrds();
            const Vector &c2 = theNodes[n2]->getCrds();

            for (int i = 0; i < 2; i++) {
                v1(i) = c1(i);
                v2(i) = c2(i);
            }
            v1(2) = 0.0;
            v2(2) = 0.0;

            if (displayMode >= 0) {
                const Vector &d1 = theNodes[n1]->getDisp();
                const Vector &d2 = theNodes[n2]->getDisp();
                for (int i = 0; i < 2; i++) {
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
//
//  This is the critical piece for getting output from the element.
//  The recorder framework works in two phases:
//
//    Phase 1 (setResponse): The recorder asks "can you record X?"
//      → We create a Response object with:
//        - A responseID (integer we pick, passed to getResponse later)
//        - A template Vector/Matrix of the correct size
//        - XML output tags for column headers (essential for file recorders)
//
//    Phase 2 (getResponse): At each committed step, the recorder calls
//      this with the responseID → we fill in the actual data.
//
//  Supported responses:
//    "force" / "forces"       → 12 resisting force components
//    "stiff" / "stiffness"    → 12×12 stiffness matrix
//    "stress" / "stresses"    → σ_xx, σ_yy, σ_xy at each GP (9 values)
//    "strain" / "strains"     → ε_xx, ε_yy, γ_xy at each GP (9 values)
//    "gaussPoint"             → x, y coordinates of each GP (6 values)
//    "material" $gpNum <args> → delegate to material at GP
//    "pressure"               → pressure contribution to resisting force
//
//  For OpenSeesPy usage:
//    ops.eleResponse(tag, 'stresses')
//    ops.eleResponse(tag, 'material', 1, 'stress')
//    ops.recorder('Element', '-file', 'out.txt', '-ele', 1, 'stresses')
//
// ═══════════════════════════════════════════════════════════════════

Response *BezierTri6::setResponse(const char **argv, int argc,
                                   OPS_Stream &output)
{
    Response *theResponse = 0;

    // ═══════════════════════════════════════════════════════════════
    //  the Ladruno recorder geometry self-declaration (element contract Part A)
    //
    //  These probes are answered BEFORE the ElementOutput block so they
    //  emit top-level tags (ElementBasis), not nested under ElementOutput.
    //  They let a self-describing recorder reconstruct the Bernstein
    //  geometric map x(ξ)=ΣNₐ(ξ)Pₐ with NO element-class-specific reader
    //  logic. Without basisInfo the recorder assumes a Lagrange basis and
    //  would mis-place Gauss points on CURVED Bézier elements.
    //
    //  family=bernstein, topology=tri (simplex), paramDomain=bary,
    //  rational=0 (BezierTri6 is non-rational ⇒ no controlPointWeights).
    //  numCtrl=6, numGP=3 (the 3-pt stress/stiffness rule where results
    //  live — NOT the 6-pt mass rule). Control-point and GP order is the
    //  standard FE hierarchical order (vertices then mid-edges).
    // ═══════════════════════════════════════════════════════════════
    if (strcmp(argv[0], "basisInfo") == 0) {
        output.tag("ElementBasis");
        output.attr("topology",    "tri");
        output.attr("family",      "bernstein");
        output.attr("paramDomain", "bary");
        output.attr("rational",    0);
        output.attr("numCtrl",     NEN);      // 6 control points
        output.attr("numGP",       NGAUSS);   // 3 Gauss points (result stations)
        output.attr("orderU",      2);        // total polynomial degree (simplex)
        output.endTag();                      // ElementBasis
        return new ElementResponse(this, 101, ID(1));   // non-null sentinel
    }
    if (strcmp(argv[0], "integrationPoints") == 0)
        return new ElementResponse(this, 102, Matrix(NGAUSS, 2));
    if (strcmp(argv[0], "integrationWeights") == 0)
        return new ElementResponse(this, 103, Vector(NGAUSS));

    output.tag("ElementOutput");
    output.attr("eleType", "BezierTri6");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes(0));
    output.attr("node2", connectedExternalNodes(1));
    output.attr("node3", connectedExternalNodes(2));
    output.attr("node4", connectedExternalNodes(3));
    output.attr("node5", connectedExternalNodes(4));
    output.attr("node6", connectedExternalNodes(5));

    // ─── Element resisting forces ─────────────────────────────
    if (strcmp(argv[0], "force") == 0 ||
        strcmp(argv[0], "forces") == 0) {

        for (int i = 1; i <= NEN; i++) {
            char buf[32];
            sprintf(buf, "P1_%d", i); output.tag("ResponseType", buf);
            sprintf(buf, "P2_%d", i); output.tag("ResponseType", buf);
        }

        theResponse = new ElementResponse(this, 1, Vector(NELD));
    }

    // ─── Stiffness matrix ─────────────────────────────────────
    else if (strcmp(argv[0], "stiff") == 0 ||
             strcmp(argv[0], "stiffness") == 0) {

        theResponse = new ElementResponse(this, 2, Matrix(NELD, NELD));
    }

    // ─── Material/integration point response ──────────────────
    // Usage: material $gpNum <material args>
    //   e.g.: ops.eleResponse(1, 'material', 1, 'stress')
    //         ops.eleResponse(1, 'material', 2, 'strain')
    //         ops.eleResponse(1, 'material', 1, 'tangent')
    else if (strcmp(argv[0], "material") == 0 ||
             strcmp(argv[0], "integrPoint") == 0) {

        if (argc < 2) {
            opserr << "BezierTri6::setResponse - need GP number after 'material'\n";
            output.endTag();
            return 0;
        }

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= NGAUSS) {
            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("eta", GP3_xi[pointNum-1][0]);
            output.attr("neta", GP3_xi[pointNum-1][1]);

            // Compute physical coordinates of this GP
            double N[NEN];
            shapeFunctions(GP3_xi[pointNum-1][0], GP3_xi[pointNum-1][1], N);
            double x = 0.0, y = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
            }
            output.attr("xLoc", x);
            output.attr("yLoc", y);

            // Delegate to material
            theResponse = theMaterial[pointNum-1]->setResponse(
                &argv[2], argc-2, output);

            output.endTag();  // GaussPoint
        }
    }

    // ─── Stresses at ALL Gauss points ─────────────────────────
    // Returns: σ_xx1, σ_yy1, σ_xy1, σ_xx2, σ_yy2, σ_xy2, ...
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
            output.tag("ResponseType", "sigma_xy");

            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }

        theResponse = new ElementResponse(this, 3, Vector(NSTRESS * NGAUSS));
    }

    // ─── Strains at ALL Gauss points ──────────────────────────
    // Returns: ε_xx1, ε_yy1, γ_xy1, ε_xx2, ε_yy2, γ_xy2, ...
    else if (strcmp(argv[0], "strain") == 0 ||
             strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);

            output.tag("NdMaterialOutput");
            output.tag("ResponseType", "eps_xx");
            output.tag("ResponseType", "eps_yy");
            output.tag("ResponseType", "gamma_xy");

            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }

        theResponse = new ElementResponse(this, 4, Vector(NSTRESS * NGAUSS));
    }

    // ─── Gauss point physical coordinates ─────────────────────
    // Returns: x1, y1, x2, y2, x3, y3 (physical coords of 3 GPs)
    // Useful for post-processing: user needs to know WHERE stress is
    else if (strcmp(argv[0], "gaussPoint") == 0 ||
             strcmp(argv[0], "gaussCoord") == 0 ||
             strcmp(argv[0], "gpCoord") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("ResponseType", "x");
            output.tag("ResponseType", "y");
        }

        theResponse = new ElementResponse(this, 5, Vector(2 * NGAUSS));
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


int BezierTri6::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {

    case 1:  // resisting forces
        return eleInfo.setVector(this->getResistingForce());

    case 2:  // tangent stiffness
        return eleInfo.setMatrix(this->getTangentStiff());

    case 3: {
        // ─── Stresses at all Gauss points ─────────────────────
        //
        // The material objects store the committed stress from the
        // last converged step. We just collect them in order:
        //   [σ_xx, σ_yy, σ_xy] at GP1, GP2, GP3
        //
        static Vector stressVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &sigma = theMaterial[i]->getStress();
            for (int j = 0; j < NSTRESS; j++)
                stressVec(i * NSTRESS + j) = sigma(j);
        }
        return eleInfo.setVector(stressVec);
    }

    case 4: {
        // ─── Strains at all Gauss points ──────────────────────
        //
        // Note: γ_xy is engineering shear strain (= 2ε_xy).
        // This matches the Voigt convention used throughout.
        //
        static Vector strainVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &eps = theMaterial[i]->getStrain();
            for (int j = 0; j < NSTRESS; j++)
                strainVec(i * NSTRESS + j) = eps(j);
        }
        return eleInfo.setVector(strainVec);
    }

    case 5: {
        // ─── Gauss point physical coordinates ─────────────────
        //
        // x(ξ) = Σ Nₐ(ξ) · Pₐ  at each GP
        //
        // This is essential for plotting stress contours: the user
        // needs to know the (x,y) location where each σ lives.
        //
        static Vector gpCoords(2 * NGAUSS);
        for (int gp = 0; gp < NGAUSS; gp++) {
            double N[NEN];
            shapeFunctions(GP3_xi[gp][0], GP3_xi[gp][1], N);

            double x = 0.0, y = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
            }
            gpCoords(2*gp)     = x;
            gpCoords(2*gp + 1) = y;
        }
        return eleInfo.setVector(gpCoords);
    }

    case 6: {
        // ─── Pressure contribution ────────────────────────────
        //
        // If a uniform pressure p is applied on the element,
        // the equivalent nodal force is:
        //   F^press = ∫_Ω Nᵀ {0, -p} dΩ · t  (for vertical pressure)
        //
        // For now, pressure acts as a body force in the negative y-direction.
        // A more general implementation would apply it on a specified edge.
        //
        static Vector pressVec(NELD);
        pressVec.Zero();

        if (pressure != 0.0) {
            for (int gp = 0; gp < NGAUSS; gp++) {
                double xi1 = GP3_xi[gp][0];
                double xi2 = GP3_xi[gp][1];
                double w = GP3_w[gp];

                double N[NEN];
                shapeFunctions(xi1, xi2, N);

                double dN[2][NEN];
                shapeDerivatives(xi1, xi2, dN);

                double J[2][2];
                double dN_dx[2][NEN];
                double detJ = computeJacobian(dN, J, dN_dx);

                double factor = w * detJ * thickness * pressure;

                for (int a = 0; a < NEN; a++) {
                    // Pressure acts as body force in y-direction
                    pressVec(2*a + 1) += factor * N[a];
                }
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

    // ─── the Ladruno recorder geometry probes (element contract Part A) ──
    case 101:
        // basisInfo sentinel — metadata was emitted via the stream in
        // setResponse; nothing to fill here.
        return 0;

    case 102: {
        // GP_PARAM — parametric (barycentric) coords of the 3 Gauss
        // points, rows in the element's GP order. paramDomain="bary":
        // store the 2 free area coords (ξ₁, ξ₂); ξ₃ = 1-ξ₁-ξ₂.
        static Matrix xi(NGAUSS, 2);
        for (int g = 0; g < NGAUSS; g++) {
            xi(g, 0) = GP3_xi[g][0];
            xi(g, 1) = GP3_xi[g][1];
        }
        return eleInfo.setMatrix(xi);
    }

    case 103: {
        // GP_WEIGHT — quadrature weights, same GP order as GP_PARAM.
        // (Reference-triangle weights; they sum to 1/2 = ref. area.)
        static Vector w(NGAUSS);
        for (int g = 0; g < NGAUSS; g++)
            w(g) = GP3_w[g];
        return eleInfo.setVector(w);
    }

    default:
        return -1;
    }
}

// ═══════════════════════════════════════════════════════════════════
//  PARAMETER INTERFACE
// ═══════════════════════════════════════════════════════════════════
//
// setParameter is what `parameter` / `addToParameter $ptag element $tag
// <args>` reach. Three routes, mirroring SixNodeTri / LadrunoCST:
//
//   "pressure"            → register THIS element (responseID 2); the
//                           Parameter drives `pressure` through
//                           updateParameter below.
//   "material" $gp <args> → delegate to the NDMaterial at GP $gp (1-based).
//   anything else         → broadcast to every GP material. This is the
//                           stress-control path: commitStressIncrementXX
//                           etc. register the MATERIALS with the Parameter,
//                           which then dispatches updateParameter straight
//                           to them (no element-level hop needed).
//
// Without this override, Element::setParameter returns -1 and
// addToParameter is a SILENT no-op — staged stress control (apeGmsh
// ops.initial_stress / STKO stressControl) does nothing through this
// element while Tri31/SixNodeTri hosts work fine.

int
BezierTri6::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    int res = -1;

    // surface pressure on the element itself
    if (strcmp(argv[0], "pressure") == 0)
        return param.addObject(2, this);

    // a specific Gauss-point material: "material $gp <args>"
    if ((strstr(argv[0], "material") != 0) &&
        (strcmp(argv[0], "materialState") != 0)) {
        if (argc < 3)
            return -1;
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= NGAUSS)
            return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
        return -1;
    }

    // otherwise a forall-material parameter — broadcast to every GP
    for (int i = 0; i < NGAUSS; i++) {
        int matRes = theMaterial[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }

    return res;
}

int
BezierTri6::updateParameter(int parameterID, Information &info)
{
    switch (parameterID) {
    case 2:
        // pressure enters getResistingForce on the fly — no nodal-load
        // recompute needed (unlike SixNodeTri::setPressureLoadAtNodes).
        pressure = info.theDouble;
        return 0;
    default:
        return -1;
    }
}
