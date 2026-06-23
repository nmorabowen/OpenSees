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

// ADR-39 — LadrunoContactFE: contact FE_Element adapter (NTS penalty + Coulomb friction).
// Three constructor modes, shipped incrementally:
//   EMPTY       (P1a, #345) — empty-connectivity zero adapter; legacy/null path, graph-neutral,
//                             bitwise-identical to no-contact (kept for the byte-identity baseline).
//   RIGID_PLANE (P2a, #346) — slave node vs a rigid analytical plane; conn={slave}, penalty normal.
//   SEGMENT     (P2b/P3/P3.5, #354/#360/#361) — faceted node-to-segment: real connectivity
//                             FE_Element(tag, 1+nps, 3*(1+nps)); getResidual assembles Bᵀ·tn +
//                             Coulomb friction; addKtToTang assembles kn·BᵀB + the consistent
//                             friction tangent. Stateless view — per-pair path state on LadrunoContactDomain.
// See LadrunoContactFE.h for the design rationale.

#include "LadrunoContactFE.h"
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Domain.h>                 // Ladruno: ADR-39 P3 (lazy engine re-fetch)
#include <LadrunoContactDomain.h>   // Ladruno: ADR-39 P3 (per-pair friction state)
#include <LadrunoContactKernel.h>     // Ladruno: ADR-39 P2b (penalty normal-law / traction)
#include <LadrunoContactProjection.h> // Ladruno: ADR-41 A2 (closest-point projection geometry)
#include <LadrunoFrictionKernel.h>    // Ladruno: ADR-41 A1 (Coulomb/Tresca friction return map + tangent)
#include <LadrunoMortarKernel.h>      // Ladruno: ADR-41 C1/C2.1 (clipped-GP mortar D,M,g̃)

LadrunoContactFE::LadrunoContactFE(int tag)
  : FE_Element(tag, /*numDOF_Group=*/0, /*ndof=*/0),
    resid(0), tang(0, 0), mode(EMPTY), theSlave(0), ndm(0), kn(0.0), nps(0),
    kt(0.0), mu(0.0), theDomain(0), contactTag(0), segIndex(0), consistentTan(false)
{
    // myDOF_Groups and myID are size 0 (empty connectivity): the adapter adds NO
    // edges to the DOF graph, so the numberer permutation is untouched and the
    // result is bitwise-identical to no-contact. (P1a)
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0;
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
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0;
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
    for (int i = 0; i < 4; i++) { mortarSlave[i] = 0; mortarMaster[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0;
}

// C2.1/C2.2 — clipped-GP MORTAR contact (one slave facet vs one master facet). theDomain
// (C2.2) reaches the engine for the per-node λ_N + global gap; null ⇒ the C2.1 penalty path.
LadrunoContactFE::LadrunoContactFE(int tag, Node **slaveNodes, int nps_s,
                                   Node **masterNodes, int nps_m, double epsN,
                                   const double odir[3], int contactTag_, int slaveFacetIndex_,
                                   Domain *dom, double mu_, double epsT_, double cohesion_,
                                   double tauMax_, bool consistentTan_)
  : FE_Element(tag, /*numDOF_Group=*/nps_s + nps_m, /*ndof=*/3 * (nps_s + nps_m)),
    resid(3 * (nps_s + nps_m)), tang(3 * (nps_s + nps_m), 3 * (nps_s + nps_m)),
    mode(MORTAR), theSlave(0), ndm(3), kn(epsN), nps(0),
    kt(epsT_), mu(mu_), theDomain(dom), contactTag(contactTag_), segIndex(0),
    consistentTan(consistentTan_), npsS(nps_s), npsM(nps_m), slaveFacetIndex(slaveFacetIndex_),
    mortarCohesion(cohesion_), mortarTauMax(tauMax_)
{
    // Connectivity = each slave-facet DOF_Group then each master-facet DOF_Group.
    // setID() fills myID = [slave_0 xyz | … | master_0 xyz | …]; mortarActive()/the
    // residual + tangent loops below assume this [slave block | master block] layout.
    // Every node is ndf==3 (the handler guards it) ⇒ exact ndof match.
    for (int i = 0; i < nps_s; i++) {
        mortarSlave[i] = slaveNodes[i];
        if (slaveNodes[i] != 0) {
            DOF_Group *dg = slaveNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(i) = dg->getTag();
        }
    }
    for (int i = 0; i < nps_m; i++) {
        mortarMaster[i] = masterNodes[i];
        if (masterNodes[i] != 0) {
            DOF_Group *dg = masterNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(nps_s + i) = dg->getTag();
        }
    }
    for (int i = nps_s; i < 4; i++) mortarSlave[i] = 0;
    for (int i = nps_m; i < 4; i++) mortarMaster[i] = 0;
    for (int i = 0; i < 4; i++) segNode[i] = 0;
    for (int d = 0; d < 3; d++) { orientDir[d] = odir[d]; planeP0[d] = 0.0; planeN[d] = 0.0; }
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
    if (!LadrunoContactProjection::evalSegment(nps, Xseg, xs, orientDir, gap, n, N))
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
        LadrunoFrictionKernel::tangentPart(drel, n, gTvec);
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
    LadrunoFrictionKernel::frictionTangentBlock(gTeff, gpT, n, tn, kn, kt, mu, consistent, Kss);
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

// C2.1: integrate the mortar facet pair at the CURRENT trial config. Reads slave/master
// facet node trial positions (X+u), calls LadrunoMortarKernel::integratePair, and (for a
// non-empty overlap) returns D,M,g̃ + the per-facet master normal n. False ⇒ no overlap.
bool
LadrunoContactFE::mortarActive(double D[4][4], double M[4][4], double g[4], double n[3]) const
{
    double Xs[4][3], Xm[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &X = mortarSlave[i]->getCrds();
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xs[i][d] = X(d) + u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &X = mortarMaster[i]->getCrds();
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xm[i][d] = X(d) + u(d);
    }
    LadrunoMortarKernel::PairResult pr;
    if (LadrunoMortarKernel::integratePair(npsS, Xs, npsM, Xm, orientDir, pr) != 0)
        return false;                                // empty/degenerate overlap
    for (int i = 0; i < 4; i++) {
        g[i] = pr.g[i];
        for (int j = 0; j < 4; j++) { D[i][j] = pr.D[i][j]; M[i][j] = pr.M[i][j]; }
    }
    // per-facet master normal (flat facet ⇒ the per-GP projection normal), oriented
    // toward orientDir — the same n the weighted gap g̃ used inside integratePair.
    return LadrunoMortarKernel::facetNormal(npsM, Xm, orientDir, n);
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
                    LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, tn, kt, mu,
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
    } else if (mode == MORTAR) {
        // C2.2 — augmented Lagrange. The augmented nodal pressure is
        //   p_I = min(0, λ_I + epsN·ḡ_I^facet),  ḡ_I^facet = g̃_I^facet / a_I^facet,
        // where λ_I is the per-GLOBAL-slave-node multiplier (committed-only, on the Domain) and
        // ḡ_I^facet is THIS facet's local normalized gap. Using the LOCAL gap (not a running
        // global one) keeps the per-facet residual sweep DETERMINISTIC — NewtonRaphson forms the
        // residual facet-by-facet, so a shared node reading a running global gap would see a
        // different pressure in each facet (an order-dependent, non-Newton residual). The λ_I
        // term assembles globally (shared ⇒ Σ_facets D_KI^facet λ_I = D_KI^global λ_I), and the
        // penalty term → 0 under augmentation, so the converged state is consistent + epsN-
        // independent (oracle T8). The GLOBAL gap is accumulated here ONLY for the commit-time
        // Uzawa update of λ_I + the ‖ḡ‖ query — never read back into THIS sweep's force.
        // Force along n (self-equilibrating ⇒ Σφ=1): F^s_K = −(D·p)_K n, F^m_L = +(Mᵀ·p)_L n.
        double D[4][4], M[4][4], g[4], n[3];
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (mortarActive(D, M, g, n)) {
            double p[4] = {0, 0, 0, 0};
            for (int I = 0; I < npsS; I++) {
                double aFacet = 0.0;                      // a_I^facet = Σ_J D_IJ = ∫N_I dΓ (this facet)
                for (int J = 0; J < npsS; J++) aFacet += D[I][J];
                if (aFacet <= 1e-300) continue;           // unreferenced slave node
                double lambdaI = 0.0;
                if (cd != 0) {
                    // accumulate the GLOBAL weighted gap (for commit/query) + read λ_I.
                    int nodeTag = mortarSlave[I]->getTag();
                    cd->accumulateMortarGap(contactTag, nodeTag, this->getTag(), g[I], aFacet, kn);
                    lambdaI = cd->getOrCreateMortarNormalState(contactTag, nodeTag).lambdaN;
                }
                double pr = lambdaI + kn * (g[I] / aFacet);   // λ_I + epsN·ḡ_I^facet (kn = epsN)
                p[I] = (pr < 0.0) ? pr : 0.0;             // active iff compression (KKT clamp)
            }
            for (int K = 0; K < npsS; K++) {              // slave block: −(D·p)_K n
                double Dp = 0.0;
                for (int I = 0; I < npsS; I++) Dp += D[K][I] * p[I];
                for (int d = 0; d < 3; d++) resid(3 * K + d) = -Dp * n[d];
            }
            for (int L = 0; L < npsM; L++) {              // master block: +(Mᵀ·p)_L n
                double Mp = 0.0;
                for (int I = 0; I < npsS; I++) Mp += M[I][L] * p[I];
                for (int d = 0; d < 3; d++) resid(3 * (npsS + L) + d) = Mp * n[d];
            }
            // --- C3.1 Coulomb/Tresca friction (force only; tangent is C3.2) ---
            // mu≤0 ∧ c≤0 ∧ τmax≤0 SHORT-CIRCUITS before any slot touch ⇒ byte-identical to the
            // frictionless C2 path (the NTS P3 `mu>0` guard, generalized to the unified cone).
            if ((mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0) && cd != 0)
                addMortarFriction(D, M, n, p, cd);
        } else if (cd != 0) {
            // overlap empty THIS eval (the pair separated / the clip rejected it): zero this
            // facet's contribution so it stops biasing the shared node's accumulated global gap.
            for (int I = 0; I < npsS; I++)
                cd->accumulateMortarGap(contactTag, mortarSlave[I]->getTag(), this->getTag(),
                                        0.0, 0.0, kn);
        }
    }
    return resid;
}

// C3.1 — Coulomb/Tresca friction FORCE on the mortar interface (penalty, λ_T≡0). For each slave
// node in NORMAL contact (p_normal[I] < 0), build the LOCAL weighted tangential gap, run the
// SHIPPED LadrunoFrictionKernel return map with the nodal normal pressure N_I = −p_normal[I], and
// scatter the (already-negated, motion-opposing) tangential traction tFric_I via the D/M operators
// exactly like the normal force: f^s_K += Σ_I D_KI tFric_I, f^m_L += −Σ_I M_IL tFric_I (self-
// equilibrating via Σφ=1 — oracle T2). The slip is per-GLOBAL-slave-node committed state on the
// Domain (gpT/gT0), the C3 analogue of λ_N (mirrors the NTS FrictionState; design [[_adr41_c3_design]]).
void
LadrunoContactFE::addMortarFriction(const double D[4][4], const double M[4][4], const double n[3],
                                    const double p_normal[4], LadrunoContactDomain *cd)
{
    // facet node DISPLACEMENTS (NOT positions). The closest-point projection makes the weighted
    // relative POSITION ∫N_I(x_s − x_m(ξ̄)) purely NORMAL (n·r = g̃), so positions carry NO
    // tangential slip information — the tangential slip is the weighted relative DISPLACEMENT
    // ∫N_I(u_s − u_m(ξ̄)), tangential part (the ADR-39 SEGMENT path uses exactly u_s − ΣN_i u_i).
    double us[4][3], um[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) us[i][d] = u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) um[i][d] = u(d);
    }
    double tFric[4][3] = {{0}};
    for (int I = 0; I < npsS; I++) {
        if (p_normal[I] >= 0.0) continue;             // friction only on in-contact nodes
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;
        double N_I = -p_normal[I];                    // contact pressure magnitude (≥ 0)
        // weighted relative DISPLACEMENT r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K (= ∫N_I(u_s − u_m) dΓ);
        // its tangential part / a_I is the area-normalised mortar tangential SLIP (gT0-referenced
        // below, so a late-engaging node's pre-contact drift is not a spurious stick traction).
        double r[3] = {0, 0, 0};
        for (int J = 0; J < npsS; J++)
            for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
        for (int K = 0; K < npsM; K++)
            for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
        double rn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
        double gbarT[3];
        for (int d = 0; d < 3; d++) gbarT[d] = (r[d] - rn * n[d]) / aFacet;   // tangential, normalised

        LadrunoContactDomain::MortarNormalState &st =
            cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag());
        // engagement origin captured ONCE at first contact (else pre-contact tangential drift
        // becomes a spurious stick traction — the ADR-39 P3 MAJOR-1, reused).
        if (!st.engaged) {
            for (int d = 0; d < 3; d++) st.gT0[d] = gbarT[d];
            st.engaged = true;
        }
        double gTeff[3];
        for (int d = 0; d < 3; d++) gTeff[d] = gbarT[d] - st.gT0[d];
        double tF[3], gpTtrial[3];
        // N for the cone = the nodal normal pressure; epsT rides kt; trial = pure fn of committed
        // gpT ⇒ idempotent across re-evals. Returns the APPLIED (negated) traction opposing motion.
        LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, N_I, kt, mu, tF, gpTtrial,
                                                 mortarCohesion, mortarTauMax);
        for (int d = 0; d < 3; d++) { st.gpTtrial[d] = gpTtrial[d]; tFric[I][d] = tF[d]; }
    }
    // scatter via D (slave) / −M (master), like the normal force but a VECTOR traction
    for (int K = 0; K < npsS; K++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += D[K][I] * tFric[I][d];
            resid(3 * K + d) += s;
        }
    for (int L = 0; L < npsM; L++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += M[I][L] * tFric[I][d];
            resid(3 * (npsS + L) + d) += -s;
        }
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
    } else if (mode == MORTAR) {
        addMortarTang(fact);
    }
}

// C2.1/C2.2 — K_c = epsN · B̃ᵀ diag(act/a) B̃ ⊗ (n⊗n), B̃ = [D, −M] over [slave|master] nodes.
// Material/penalty tangent only — the geometric ∂{D,M,n}/∂u terms are the C1 linearization
// stub, deferred (NTS shipped without ∂n/∂u). SPD on the active set; matches proto_c2_alm.
// C2.2: the active mask uses the SAME per-facet pressure p_I = min(0, λ_I + epsN·ḡ_I^facet)
// the residual used (local gap + the global multiplier λ_I), so the tangent is the derivative
// of the augmented residual at the frozen active set. The λ_I offset only shifts the active
// SET (∂λ/∂u = 0 within a sweep), so the penalty Gram block is unchanged. theDomain==0 ⇒ the
// C2.1 fallback (λ≡0).
void
LadrunoContactFE::addMortarTang(double fact, bool initialStiff)
{
    double D[4][4], M[4][4], g[4], n[3];
    if (!mortarActive(D, M, g, n)) return;
    int nN = npsS + npsM;
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    double W[4] = {0.0, 0.0, 0.0, 0.0};               // act_I / a_I^facet
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;
        double lambdaI = (cd != 0)
            ? cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag()).lambdaN
            : 0.0;
        double pr = lambdaI + kn * (g[I] / aFacet);   // same p_I as the residual (kn = epsN)
        if (pr < 0.0) W[I] = 1.0 / aFacet;            // active iff compression
    }
    for (int A = 0; A < nN; A++) {
        for (int B = 0; B < nN; B++) {
            double Ks = 0.0;                           // K_scalar[A][B]
            for (int I = 0; I < npsS; I++) {
                double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                Ks += bIA * W[I] * bIB;
            }
            Ks *= kn;                                 // epsN
            if (Ks == 0.0) continue;
            for (int dA = 0; dA < 3; dA++)
                for (int dB = 0; dB < 3; dB++)
                    tang(3 * A + dA, 3 * B + dB) += fact * Ks * n[dA] * n[dB];
        }
    }

    // --- C3.2 mortar friction TANGENT (the SYMMETRIC consistent tangent) ---
    // K[A][B] += Σ_I (b_IA b_IB / a_I) K_ss_I, b_IA the B̃=[D,−M] operator row, K_ss_I the 3×3
    // friction block (LadrunoFrictionKernel::frictionTangentBlock) at the SAME per-node slip the
    // residual used. Default SYMMETRIC (drop the non-symmetric Coulomb Csl — solver-safe on any
    // SOE, the ADR-39 P3.5 Q2 rule; the full mortar Csl couples through the normal-gap operator
    // and is deferred with the geometric ∂{D,M,n}/∂u terms). Makes implicit frictional Newton
    // converge (singular without it). Oracle proto_c3_mortar_friction.py T6 FD-checks this assembly.
    if ((mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0) && cd != 0) {
        double us[4][3], um[4][3];
        for (int i = 0; i < npsS; i++) {
            const Vector &u = mortarSlave[i]->getTrialDisp();
            for (int d = 0; d < 3; d++) us[i][d] = u(d);
        }
        for (int i = 0; i < npsM; i++) {
            const Vector &u = mortarMaster[i]->getTrialDisp();
            for (int d = 0; d < 3; d++) um[i][d] = u(d);
        }
        for (int I = 0; I < npsS; I++) {
            double aFacet = 0.0;
            for (int J = 0; J < npsS; J++) aFacet += D[I][J];
            if (aFacet <= 1e-300) continue;
            LadrunoContactDomain::MortarNormalState &st =
                cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag());
            double pr = st.lambdaN + kn * (g[I] / aFacet);
            if (pr >= 0.0) continue;                  // friction only on in-contact nodes
            double N_I = -pr;
            // displacement-based tangential slip (same as the residual), engagement-referenced.
            double r[3] = {0, 0, 0};
            for (int J = 0; J < npsS; J++)
                for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
            for (int K = 0; K < npsM; K++)
                for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
            double rn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
            double gTeff[3];
            for (int d = 0; d < 3; d++) gTeff[d] = (r[d] - rn * n[d]) / aFacet - st.gT0[d];
            // initial-stiffness path ⇒ force the SPD STICK tangent kt·P_t: pass gTeff == gpT so the
            // trial traction is zero (‖tT*‖ ≤ cap ⇒ the kernel returns the stick block). Avoids the
            // rank-deficient slip tangent stalling Modified/Initial-Newton (gate MINOR-1, mirrors SEGMENT).
            const double *gtForKss = initialStiff ? st.gpT : gTeff;
            double Kss[3][3];
            LadrunoFrictionKernel::frictionTangentBlock(gtForKss, st.gpT, n, N_I, kn, kt, mu,
                                                        /*consistent=*/false, Kss,
                                                        mortarCohesion, mortarTauMax);
            // scatter: tang(3A+i,3B+j) += fact·(b_IA b_IB / a_I)·K_ss[i][j]
            for (int A = 0; A < nN; A++) {
                double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                if (bIA == 0.0) continue;
                for (int B = 0; B < nN; B++) {
                    double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                    double w = fact * bIA * bIB / aFacet;
                    if (w == 0.0) continue;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            tang(3 * A + i, 3 * B + j) += w * Kss[i][j];
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
    } else if (mode == MORTAR) {
        addMortarTang(fact, /*initialStiff=*/true);  // penalty K_initial == K_current; friction =
                                                     // the SPD stick tangent (geometric terms deferred)
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
