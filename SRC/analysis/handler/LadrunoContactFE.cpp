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

// ADR-39 P1a — LadrunoContactFE implementation (empty-connectivity zero adapter).
// See LadrunoContactFE.h for the design rationale.

#include "LadrunoContactFE.h"
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>

LadrunoContactFE::LadrunoContactFE(int tag)
  : FE_Element(tag, /*numDOF_Group=*/0, /*ndof=*/0),
    resid(0), tang(0, 0), mode(EMPTY), theSlave(0), ndm(0), kn(0.0)
{
    // myDOF_Groups and myID are size 0 (empty connectivity): the adapter adds NO
    // edges to the DOF graph, so the numberer permutation is untouched and the
    // result is bitwise-identical to no-contact. (P1a)
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; }
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, int ndm_,
                                   const double p0[3], const double n[3], double kn_)
  : FE_Element(tag, /*numDOF_Group=*/1, /*ndof=*/ndm_),
    resid(ndm_), tang(ndm_, ndm_),
    mode(RIGID_PLANE), theSlave(slaveNode), ndm(ndm_), kn(kn_)
{
    // Connectivity = the slave node's DOF_Group; setID() fills myID with its first
    // ndm equation numbers (the translational DOFs). Same pattern as PenaltySP_FE.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0)
            myDOF_Groups(0) = dg->getTag();
    }
    for (int d = 0; d < 3; d++) { planeP0[d] = p0[d]; planeN[d] = n[d]; }
}

LadrunoContactFE::~LadrunoContactFE()
{
}

double
LadrunoContactFE::rigidPlaneGap(void) const
{
    // g = n . (x_s - p0), x_s = X_s + u_s  (first ndm components)
    const Vector &X = theSlave->getCrds();
    const Vector &u = theSlave->getTrialDisp();
    double g = 0.0;
    for (int d = 0; d < ndm; d++)
        g += planeN[d] * (X(d) + u(d) - planeP0[d]);
    return g;
}

const Vector &
LadrunoContactFE::getResidual(Integrator *)
{
    resid.Zero();
    if (mode == RIGID_PLANE) {
        double g = rigidPlaneGap();
        if (g < 0.0) {                       // penetration -> restoring force +n
            double tn = -kn * g;             // tn = kn*|g| > 0 (Macaulay <-g>_+)
            for (int d = 0; d < ndm; d++)
                resid(d) = tn * planeN[d];   // drives the slave toward g=0 (PenaltySP convention)
        }
    }
    return resid;
}

const Matrix &
LadrunoContactFE::getTangent(Integrator *theIntegrator)
{
    // Route through the integrator so it decides the K/C/M combination. CDL's
    // formEleTangent calls only addMtoTang (no-op here) -> tang stays 0, so the
    // contact stiffness CANNOT pollute the explicit lumped-mass matrix (the bug:
    // a re-formed tangent injecting kn into M makes M_xx ~ kn and a ~ 0, so the
    // node coasts through). Newmark calls addKtToTang(c1) -> c1*K_c.
    tang.Zero();
    if (theIntegrator != 0)
        theIntegrator->formEleTangent(this);   // calls my zeroTangent/addKtToTang/...
    return tang;
}

void
LadrunoContactFE::zeroTangent(void)
{
    tang.Zero();   // formEleTangent zeroes the element tangent before assembling
}

void
LadrunoContactFE::addKtToTang(double fact)
{
    // K_c = kn (n (x) n) over the slave DOFs, only when the pair is penetrating.
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * kn * planeN[i] * planeN[j];
    }
}

void
LadrunoContactFE::addKiToTang(double fact)
{
    // Initial-stiffness algorithms (Newton -initial, ModifiedNewton -initial,
    // HALL_TANGENT) form the LHS from K_initial. For a flat rigid plane the normal
    // is constant, so K_initial == K_current == kn (n (x) n) when penetrating;
    // mirror addKtToTang so the contact stiffness is not silently dropped (the base
    // addKiToTang early-returns on myEle == 0). Residual is unaffected either way.
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * kn * planeN[i] * planeN[j];
    }
}

void
LadrunoContactFE::addCtoTang(double)
{
    // contact carries no damping in P2 -> no-op
}

void
LadrunoContactFE::addMtoTang(double)
{
    // contact pairs contribute no mass (mass lives on the real structural
    // elements at the same nodes); no-op keeps the explicit LHS = lumped M only.
}

// Size-0 returns for the off-hot-path force variants (base would warn on myEle==0).
const Vector &LadrunoContactFE::getTangForce(const Vector &, double) { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getK_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getKi_Force(const Vector &, double)  { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getC_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getM_Force(const Vector &, double)   { resid.Zero(); return resid; }
