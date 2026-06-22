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
#include <Domain.h>                 // Ladruno: ADR-39 P3 (lazy engine re-fetch)
#include <LadrunoContactDomain.h>   // Ladruno: ADR-39 P3 (per-pair friction state)
#include <LadrunoContactKernel.h>   // Ladruno: ADR-39 P2b/P3 (header-only NTS math)

LadrunoContactFE::LadrunoContactFE(int tag)
  : FE_Element(tag, /*numDOF_Group=*/0, /*ndof=*/0),
    resid(0), tang(0, 0), mode(EMPTY), theSlave(0), ndm(0), kn(0.0), nps(0),
    kt(0.0), mu(0.0), theDomain(0), contactTag(0), segIndex(0), consistentTan(false)
{
    // myDOF_Groups and myID are size 0 (empty connectivity): the adapter adds NO
    // edges to the DOF graph, so the numberer permutation is untouched and the
    // result is bitwise-identical to no-contact. (P1a)
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) segNode[i] = 0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, int ndm_,
                                   const double p0[3], const double n[3], double kn_)
  : FE_Element(tag, /*numDOF_Group=*/1, /*ndof=*/ndm_),
    resid(ndm_), tang(ndm_, ndm_),
    mode(RIGID_PLANE), theSlave(slaveNode), ndm(ndm_), kn(kn_), nps(0),
    kt(0.0), mu(0.0), theDomain(0), contactTag(0), segIndex(0), consistentTan(false)
{
    // Connectivity = the slave node's DOF_Group; setID() fills myID with its first
    // ndm equation numbers (the translational DOFs). Same pattern as PenaltySP_FE.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0)
            myDOF_Groups(0) = dg->getTag();
    }
    for (int d = 0; d < 3; d++) { planeP0[d] = p0[d]; planeN[d] = n[d]; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) segNode[i] = 0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes,
                                   int nps_, double kn_, const double odir[3],
                                   double kt_, double mu_, Domain *dom,
                                   int contactTag_, int segIndex_, bool consistentTan_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/3 * (1 + nps_)),
    resid(3 * (1 + nps_)), tang(3 * (1 + nps_), 3 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(3), kn(kn_), nps(nps_),
    kt(kt_), mu(mu_), theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentTan(consistentTan_)
{
    // Connectivity = slave DOF_Group + each segment-node DOF_Group. setID() then
    // fills myID = [slave xyz | seg_1 xyz | ... | seg_nps xyz] (each node ndf==3 ⇒
    // its DOF_Group contributes exactly 3 ⇒ exact ndof match). The B-operator below
    // assumes this layout. The handler guards ndf==3 on every node of the pair.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0) myDOF_Groups(0) = dg->getTag();
    }
    for (int i = 0; i < nps; i++) {
        segNode[i] = segNodes[i];
        if (segNodes[i] != 0) {
            DOF_Group *dg = segNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(1 + i) = dg->getTag();
        }
    }
    for (int i = nps; i < 4; i++) segNode[i] = 0;
    // orientation direction (toward the slave's allowed half-space): the derived
    // normal is flipped to satisfy n·orientDir>0, so it's winding-immune AND stays
    // correct after the slave penetrates (a fixed direction, not a live position).
    for (int d = 0; d < 3; d++) orientDir[d] = odir[d];
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; }
}

LadrunoContactFE::~LadrunoContactFE()
{
}

bool
LadrunoContactFE::segmentActive(double &gap, double n[3], double N[4], double *B,
                                double *gTvec) const
{
    if (mode != SEGMENT || theSlave == 0) return false;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return false;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    if (!LadrunoContactKernel::evalSegment(nps, Xseg, xs, orientDir, gap, n, N))
        return false;
    // gap operator B (1×ndof) over [u_s | u_1..u_nps]: [ nᵀ | −N_i nᵀ ]
    int ndof = 3 * (1 + nps);
    for (int k = 0; k < ndof; k++) B[k] = 0.0;
    for (int d = 0; d < 3; d++) B[d] = n[d];
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) B[3 * (1 + i) + d] = -N[i] * n[d];
    // P3: relative tangential SLIP at the contact point. The closest-point projection
    // makes (x_s − x̄) ∥ n, so POSITIONS carry NO tangential information; the slip is
    // the slave DISPLACEMENT minus the interpolated master DISPLACEMENT at the
    // projection: d = u_s − Σ N_i u_i, tangential part (the engagement origin gT0 is
    // subtracted by the caller). u_s/u_i are read from the same trial config as the
    // projection, at the SAME shape weights N — one consistent contact point.
    if (gTvec != 0) {
        double ubar[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < nps; i++) {
            const Vector &ui = segNode[i]->getTrialDisp();
            for (int d = 0; d < 3; d++) ubar[d] += N[i] * ui(d);
        }
        double drel[3] = { us(0)-ubar[0], us(1)-ubar[1], us(2)-ubar[2] };
        LadrunoContactKernel::tangentPart(drel, n, gTvec);
    }
    return true;
}

void
LadrunoContactFE::addFrictionTang(double fact, const double n[3], const double N[4],
                                  double tn, const double gTeff[3], const double gpT[3],
                                  bool consistent)
{
    // K_fric = Gᵀ K_ss G, G = [I | −N_i I]. Block (a,b) of the ndof×ndof tangent is
    // w_a·w_b·K_ss with the scatter weights w = [1, −N_0, …, −N_{nps−1}] over
    // [slave, seg nodes]. Validated assembly: proto_p35_implicit_tangent.py.
    double Kss[3][3];
    LadrunoContactKernel::frictionTangentBlock(gTeff, gpT, n, tn, kn, kt, mu, consistent, Kss);
    double w[5];
    w[0] = 1.0;
    for (int i = 0; i < nps; i++) w[1 + i] = -N[i];
    int nn = 1 + nps;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < nn; b++) {
            double wab = fact * w[a] * w[b];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tang(3 * a + i, 3 * b + j) += wab * Kss[i][j];
        }
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
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15], gTvec[3];   // ndof <= 3*(1+4) = 15
        if (segmentActive(gap, n, N, B, gTvec)) {
            double tn = LadrunoContactKernel::traction(kn, gap);  // = kn*<-gap>_+ > 0
            int ndof = 3 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;        // r = Bᵀ tn (slave +tn n, master −N_i tn n)

            // --- P3 Coulomb friction (force only; tangent is P3.5) ---
            // mu<=0 SHORT-CIRCUITS before any state touch ⇒ byte-identical to the
            // frictionless P2b path (+ dodges the n̂=tT*/‖tT*‖ 0/0). The engine is
            // re-fetched lazily (wipe deletes it) and null-checked.
            if (mu > 0.0 && theDomain != 0) {
                LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
                if (cd != 0) {
                    int slaveTag = theSlave->getTag();
                    LadrunoContactDomain::FrictionState &st =
                        cd->getOrCreateFrictionState(contactTag, slaveTag, segIndex);
                    // capture the ENGAGEMENT-config tangential origin ONCE at first
                    // activation (else a late-engaging slave's pre-contact tangential
                    // drift becomes a spurious stick traction — design-gate MAJOR-1).
                    if (!st.engaged) {
                        for (int d = 0; d < 3; d++) st.gT0[d] = gTvec[d];
                        st.engaged = true;
                    }
                    double gTeff[3];
                    for (int d = 0; d < 3; d++) gTeff[d] = gTvec[d] - st.gT0[d];
                    double tFric[3], gpTtrial[3];
                    // N for the cone = current penetration force tn (design MINOR-8).
                    LadrunoContactKernel::frictionReturnMap(gTeff, st.gpT, tn, kt, mu,
                                                            tFric, gpTtrial);
                    // trial = PURE fn of committed state (idempotent across the CDL
                    // firstStep double-eval); commit() promotes it.
                    for (int d = 0; d < 3; d++) st.gpTtrial[d] = gpTtrial[d];
                    // mirror the normal block: slave += tFric, master_i += −N_i tFric.
                    // tFric is the already-negated APPLIED force (kernel BLOCKER-1),
                    // so friction OPPOSES the slave motion.
                    for (int d = 0; d < 3; d++) resid(d) += tFric[d];
                    for (int i = 0; i < nps; i++)
                        for (int d = 0; d < 3; d++)
                            resid(3 * (1 + i) + d) += -N[i] * tFric[d];
                }
            }
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
    } else if (mode == SEGMENT) {
        // K_c = kn BᵀB (main NTS term; ∂n/∂u block deferred to P2b-2 — for a FIXED
        // master the slave block kn(n⊗n) is exact and the master DOFs are constrained).
        double gap, n[3], N[4], B[15], gTvec[3];
        if (segmentActive(gap, n, N, B, gTvec)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // P3.5 friction tangent (IMPLICIT only — CDL never calls addKtToTang).
            // Reads COMMITTED gpT (not gpTtrial) so the tangent is the derivative of the
            // residual evaluated at the same state. Default consistentTan=false ⇒
            // symmetric (drop d_TN⊗n, solver-safe); true ⇒ full non-sym consistent.
            if (mu > 0.0 && theDomain != 0) {
                LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
                if (cd != 0) {
                    LadrunoContactDomain::FrictionState &st =
                        cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                    double gTeff[3];
                    for (int d = 0; d < 3; d++)
                        gTeff[d] = st.engaged ? (gTvec[d] - st.gT0[d]) : 0.0;
                    double tn = LadrunoContactKernel::traction(kn, gap);
                    addFrictionTang(fact, n, N, tn, gTeff, st.gpT, consistentTan);
                }
            }
        }
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
    } else if (mode == SEGMENT) {
        // initial-stiffness path: same kn BᵀB (flat segment ⇒ K_initial == K_current)
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // P3.5: the friction INITIAL stiffness is the STICK tangent kt·Gᵀ P_t G
            // (SPD ⇒ a Modified/Initial-Newton contraction; gate Q5). gTeff/gpT do not
            // matter (forced stick ⇒ K_ss = kt·P_t), so no engine slot is needed; the
            // stick tangent is symmetric, so consistent=false here regardless.
            if (mu > 0.0) {
                double tn = LadrunoContactKernel::traction(kn, gap);
                double zero[3] = {0.0, 0.0, 0.0};
                addFrictionTang(fact, n, N, tn, zero, zero, false);
            }
        }
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
