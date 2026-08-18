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
#include <cmath>                    // Ladruno: ADR-39 B1 (softKt tangent-basis fabs/sqrt)
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <StaticIntegrator.h>       // Ladruno: contact-review P5 (-visc static-integrator gate)
#include <CentralDifferenceLadruno.h> // Ladruno: ADR-39 B1 (explicit dt for the SOFT=1 penalty)
#include <Domain.h>                 // Ladruno: ADR-39 P3 (lazy engine re-fetch)
#include <LadrunoContactDomain.h>   // Ladruno: ADR-39 P3 (per-pair friction state)
#include <LadrunoContactKernel.h>     // Ladruno: ADR-39 P2b (penalty normal-law / traction)
#include <LadrunoContactProjection.h> // Ladruno: ADR-41 A2 (closest-point projection geometry)
#include <LadrunoFrictionKernel.h>    // Ladruno: ADR-41 A1 (Coulomb/Tresca friction return map + tangent)
#include <LadrunoMortarKernel.h>      // Ladruno: ADR-41 C1/C2.1 (clipped-GP mortar D,M,g̃)
#include <LadrunoEdgeKernel.h>        // Ladruno: ADR-57 E2 (segment-segment closest point + edge gap/B)
#include <LadrunoContact2DKernel.h>   // Ladruno: ADR-85 T1b (2D NTS projection / vertex-pair kernel)

LadrunoContactFE::LadrunoContactFE(int tag)
  : FE_Element(tag, /*numDOF_Group=*/0, /*ndof=*/0),
    resid(0), tang(0, 0), mode(EMPTY), theSlave(0), ndm(0), kn(0.0), nps(0),
    kt(0.0), mu(0.0), muc(0.0), theDomain(0), contactTag(0), segIndex(0), consistentTan(false),
    consistentNormal(false), softScale(0.0)
{
    // myDOF_Groups and myID are size 0 (empty connectivity): the adapter adds NO
    // edges to the DOF graph, so the numberer permutation is untouched and the
    // result is bitwise-identical to no-contact. (P1a)
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; edgeNode[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0; isTie = false;
    edgeAlm = false;   // ADR-57 E6: not an edge-edge adapter
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, int ndm_,
                                   const double p0[3], const double n[3], double kn_, double muc_,
                                   double softScale_, Domain *dom, int contactTag_)
  : FE_Element(tag, /*numDOF_Group=*/1, /*ndof=*/ndm_),
    resid(ndm_), tang(ndm_, ndm_),
    mode(RIGID_PLANE), theSlave(slaveNode), ndm(ndm_), kn(kn_), nps(0),
    kt(0.0), mu(0.0), muc(muc_), theDomain(dom), contactTag(contactTag_), segIndex(0),
    consistentTan(false),
    consistentNormal(false), softScale(softScale_ > 0.0 ? softScale_ : 0.0)
{
    // Connectivity = the slave node's DOF_Group; setID() fills myID with its first
    // ndm equation numbers (the translational DOFs). Same pattern as PenaltySP_FE.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0)
            myDOF_Groups(0) = dg->getTag();
    }
    for (int d = 0; d < 3; d++) { planeP0[d] = p0[d]; planeN[d] = n[d]; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; edgeNode[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0; isTie = false;
    edgeAlm = false;   // ADR-57 E6: not an edge-edge adapter
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes,
                                   int nps_, double kn_, const double odir[3],
                                   double kt_, double mu_, Domain *dom,
                                   int contactTag_, int segIndex_, bool consistentTan_, double muc_,
                                   bool consistentNormal_, double softScale_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/3 * (1 + nps_)),
    resid(3 * (1 + nps_)), tang(3 * (1 + nps_), 3 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(3), kn(kn_), nps(nps_),
    kt(kt_), mu(mu_), muc(muc_), theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentTan(consistentTan_), consistentNormal(consistentNormal_),
    softScale(softScale_ > 0.0 ? softScale_ : 0.0)
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
    for (int i = 0; i < 4; i++) edgeNode[i] = 0;   // ADR-57: not an edge-edge adapter
    // orientation direction (toward the slave's allowed half-space): the derived
    // normal is flipped to satisfy n·orientDir>0, so it's winding-immune AND stays
    // correct after the slave penetrates (a fixed direction, not a live position).
    for (int d = 0; d < 3; d++) orientDir[d] = odir[d];
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; }
    for (int i = 0; i < 4; i++) { mortarSlave[i] = 0; mortarMaster[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0; isTie = false;
    edgeAlm = false;   // ADR-57 E6: not an edge-edge adapter
    useSmoothNormal = false;                 // ADR-63 #4a: faceted by default (setSmoothNormals opts in)
    for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) nodalNorm[i][d] = 0.0;
}

// ADR-63 #4a — install this segment's nps FROZEN nodal normals (the engine's per-handle field) and
// flip onto the smoothed-normal path. nn = nps*3 row-major; a null/non-SEGMENT install is a no-op
// (keeps the faceted path).
void
LadrunoContactFE::setSmoothNormals(const double *nn, const int *se)
{
    if (nn == 0 || mode != SEGMENT) return;
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) nodalNorm[i][d] = nn[i*3 + d];
    // ADR-63 P2.1 — install the shared-edge mask for the facet-ownership guard (0 ⇒ leave all-zero
    // ⇒ guard inert). Only the segment's nps flags are meaningful; the rest stay zero.
    if (se != 0)
        for (int k = 0; k < nps; k++) sharedEdge[k] = se[k];
    useSmoothNormal = true;
}

// C2.1/C2.2 — clipped-GP MORTAR contact (one slave facet vs one master facet). theDomain
// (C2.2) reaches the engine for the per-node λ_N + global gap; null ⇒ the C2.1 penalty path.
LadrunoContactFE::LadrunoContactFE(int tag, Node **slaveNodes, int nps_s,
                                   Node **masterNodes, int nps_m, double epsN,
                                   const double odir[3], int contactTag_, int slaveFacetIndex_,
                                   Domain *dom, double mu_, double epsT_, double cohesion_,
                                   double tauMax_, bool consistentTan_, bool isTie_, double muc_,
                                   double softScale_)
  : FE_Element(tag, /*numDOF_Group=*/nps_s + nps_m, /*ndof=*/3 * (nps_s + nps_m)),
    resid(3 * (nps_s + nps_m)), tang(3 * (nps_s + nps_m), 3 * (nps_s + nps_m)),
    mode(MORTAR), theSlave(0), ndm(3), kn(epsN), nps(0),
    kt(epsT_), mu(mu_), muc(muc_), theDomain(dom), contactTag(contactTag_), segIndex(0),
    consistentTan(consistentTan_), consistentNormal(false),
    softScale(softScale_ > 0.0 ? softScale_ : 0.0),   // B2 (P5): SOFT=2 segment-based explicit penalty
    npsS(nps_s), npsM(nps_m), slaveFacetIndex(slaveFacetIndex_),
    mortarCohesion(cohesion_), mortarTauMax(tauMax_), isTie(isTie_)
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
    for (int i = 0; i < 4; i++) { segNode[i] = 0; edgeNode[i] = 0; }
    edgeAlm = false;   // ADR-57 E6: not an edge-edge adapter
    for (int d = 0; d < 3; d++) { orientDir[d] = odir[d]; planeP0[d] = 0.0; planeN[d] = 0.0; }
}

// ADR-57 E2 — EDGE_EDGE: perpendicular edge-edge penalty contact (one slave edge vs one master
// edge). theDomain reaches the engine for the per-pair committed sign (EdgeEdgeState); null ⇒ the
// sign is captured locally each eval from orientDir (deterministic, but not held across re-pairing).
LadrunoContactFE::LadrunoContactFE(int tag, Node *sNodeA, Node *sNodeB, Node *mNodeA, Node *mNodeB,
                                   double epsN, const double odir[3], int contactTag_, Domain *dom,
                                   double mu_, double kt_, double cohesion_, double tauMax_,
                                   bool consistentTan_, double softScale_, bool edgeAlm_)
  : FE_Element(tag, /*numDOF_Group=*/4, /*ndof=*/12),
    resid(12), tang(12, 12),
    mode(EDGE_EDGE), theSlave(0), ndm(3), kn(epsN), nps(0),
    kt(kt_), mu(mu_), muc(0.0), theDomain(dom), contactTag(contactTag_), segIndex(0),
    consistentTan(consistentTan_), consistentNormal(false),
    softScale(softScale_ > 0.0 ? softScale_ : 0.0),   // E5 SOFT (filter ≤0 like the other ctors)
    npsS(0), npsM(0), slaveFacetIndex(0), edgeAlm(edgeAlm_),   // E6 one-scalar ALM (decl-order init)
    mortarCohesion(cohesion_), mortarTauMax(tauMax_),
    isTie(false)
{
    // Connectivity = the 4 edge nodes' DOF_Groups (each ndf==3 ⇒ exact ndof==12). setID() then
    // fills myID = [sa xyz | sb xyz | ma xyz | mb xyz] — the layout edgeGeom()/the residual assume.
    edgeNode[0] = sNodeA; edgeNode[1] = sNodeB; edgeNode[2] = mNodeA; edgeNode[3] = mNodeB;
    for (int i = 0; i < 4; i++)
        if (edgeNode[i] != 0) {
            DOF_Group *dg = edgeNode[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(i) = dg->getTag();
        }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; }
    for (int d = 0; d < 3; d++) { orientDir[d] = odir[d]; planeP0[d] = 0.0; planeN[d] = 0.0; }
}

// ADR-85 T1b -- the 2D NTS SEGMENT (nps == 2) / concave-VERTEX (nps == 1) adapter.
// See the header for the degenerate-segment representation decision (no new Mode
// value, no new class tag -- 2D is a parameterization of SEGMENT).
// ADR-85 T2 -- kt_/mu_/consistentTan_ wire the SCALAR friction return map (WORK ITEM
// 2/3); mu_<=0 (the default) short-circuits before any state touch, exactly the 3D
// SEGMENT ctor's contract, so a frictionless 2D deck is byte-identical to T1b.
LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes, int nps_,
                                   int ndm2, double kn_, double sigma, double LrefIn,
                                   Node *prevFar, Node *nextFar,
                                   double kt_, double mu_, bool consistentTan_,
                                   double muc_, double softScale_,
                                   Domain *dom, int contactTag_, int segIndex_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/ndm2 * (1 + nps_)),
    resid(ndm2 * (1 + nps_)), tang(ndm2 * (1 + nps_), ndm2 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(ndm2), kn(kn_), nps(nps_),
    kt(kt_), mu(mu_), muc(muc_), theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentTan(consistentTan_), consistentNormal(false),
    softScale(softScale_ > 0.0 ? softScale_ : 0.0)
{
    // Connectivity = slave DOF_Group + each master-node DOF_Group. Each node has
    // ndf == ndm == 2 EXACTLY (the handler's ADR-85 equality guards), so setID()
    // fills myID = [slave xy | node_1 xy (| node_2 xy)] -- exactly ndof slots; the
    // 2D B-operators assume this layout (the 3D SEGMENT packing argument, ndm-
    // generalized).
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
    for (int i = 0; i < 4; i++) { mortarSlave[i] = 0; mortarMaster[i] = 0; edgeNode[i] = 0; }
    npsS = npsM = 0; slaveFacetIndex = 0; mortarCohesion = mortarTauMax = 0.0; isTie = false;
    edgeAlm = false;
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; orientDir[d] = 0.0; }
    useSmoothNormal = false;      // -smoothNormal is 3D-only (the handler refuses it in 2D)
    for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) nodalNorm[i][d] = 0.0;
    nts2dSigma   = sigma;
    nts2dLref    = LrefIn;
    // ADR-85 T2 (D5) -- store the TAG, not the pointer (see the header comment on
    // nts2dPrevFarNode/nts2dNextFarNode): the pointer is only used to prove "a
    // neighbour was armed" and to read its tag, both fine to do ONCE here (the
    // handler just fetched it live, this same handle()).
    nts2dHasPrevFar = (prevFar != 0);
    nts2dPrevFarTag = (prevFar != 0) ? prevFar->getTag() : 0;
    nts2dHasNextFar = (nextFar != 0);
    nts2dNextFarTag = (nextFar != 0) ? nextFar->getTag() : 0;
}

LadrunoContactFE::~LadrunoContactFE()
{
}

// ADR-85 T2 (D5) -- live re-resolution of the far-neighbour nodes through theDomain,
// by tag, on every access. Returns 0 both when no neighbour was ever armed AND when
// an armed neighbour no longer resolves (removed by the analysis since this adapter
// was built) -- the caller already treats a null far node as "open terminal end" /
// "not a vertex pair" (segment2DActive's existing null checks), so a removed
// neighbour stands the pair down to that SAME behaviour instead of reading freed
// memory. theDomain==0 (never expected for a live 2D adapter, but defensive) also
// resolves to "no neighbour".
Node *
LadrunoContactFE::nts2dPrevFarNode(void) const
{
    if (!nts2dHasPrevFar || theDomain == 0) return 0;
    return theDomain->getNode(nts2dPrevFarTag);
}

Node *
LadrunoContactFE::nts2dNextFarNode(void) const
{
    if (!nts2dHasNextFar || theDomain == 0) return 0;
    return theDomain->getNode(nts2dNextFarTag);
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
    // ADR-63 #4a: smoothed nodal-normal field (frozen at handle time) when installed; else the
    // shipped faceted normalOriented() (byte-identical). evalSegmentSmooth keeps the SAME facet
    // closest-point projection + in-bounds/penetration gate (active set unchanged ⇒ R5 slide-off
    // untouched); only the normal DIRECTION changes. orientDir is the degenerate-blend fallback.
    bool active = useSmoothNormal
        ? LadrunoContactProjection::evalSegmentSmooth(nps, Xseg, xs, nodalNorm, orientDir, gap, n, N,
                                                      sharedEdge)
        : LadrunoContactProjection::evalSegment(nps, Xseg, xs, orientDir, gap, n, N);
    if (!active)
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

// ADR-85 T1b -- the OPEN-END parametric acceptance window (dimensionless slack on
// xi past a TERMINAL chain side; see the block comment inside segment2DActive).
// Sized ~4 orders above the measured tilt-induced drift of a boundary slave on
// the G-T1b gate decks (1.25e-7 at load; scale pen*tilt/L^2 -- the 6-order figure
// is the drift against tolIn = 1e-9, not against this window) and 3 orders below
// anything a gate measures geometrically (1e-3*L past the block corner). Only
// OPEN ends see it -- interior seams keep the kernel's strict tolIn.
static const double NTS2D_END_SLACK = 1.0e-3;

// ADR-85 T1b -- the 2D NTS narrow phase. Wired VERBATIM to the
// LadrunoContact2DKernel ownership contract (ADR-85 How/1, T1a corrections):
//   1. previous segment SEG2D_IN_BOUNDS  -> the PREVIOUS SEGMENT owns
//   2. else next segment SEG2D_IN_BOUNDS -> the NEXT SEGMENT owns
//   3. else aPrev > 0 && aNext < 0       -> the VERTEX owns (UNSLACKED)
//   4. else                              -> nobody (genuinely off the patch)
// A SEGMENT adapter implements step 1 as a stand-down against its predecessor in
// surface order (nts2dPrevFar -- the convex-corner strip overlap, which in 3D is
// the accepted NTS corner double-count, is CLOSED in 2D by this rule); a VERTEX
// adapter implements steps 1-3 in full. The neighbour reads steer only the
// ACTIVE-SET decision -- forces and tangents involve only the connectivity nodes,
// exactly like the 3D in-bounds test. NO geometry lives here: every predicate is
// a kernel call (the T1a "T1b adds no geometry" contract).
bool
LadrunoContactFE::segment2DActive(double &gap, double n[3], double N[2], double B[6],
                                  double committedSide, double *liveSideOut,
                                  double *gTout) const
{
    using namespace LadrunoContact2DKernel;
    gap = 0.0; n[0] = n[1] = n[2] = 0.0; N[0] = N[1] = 0.0;
    for (int k = 0; k < 6; k++) B[k] = 0.0;
    if (liveSideOut != 0) *liveSideOut = 0.0;
    if (gTout != 0) *gTout = 0.0;
    if (mode != SEGMENT || ndm != 2 || theSlave == 0) return false;

    // ADR-85 T2 (D5) -- re-resolve the far neighbours THROUGH theDomain by tag, not
    // through a raw pointer held since construction (see the header comment on
    // nts2dPrevFarNode/nts2dNextFarNode). A neighbour removed since this adapter was
    // built resolves to 0 here -- the SAME value "never armed" already produces --
    // so every existing null check below stands this pair down exactly as it would
    // for a genuine open end / non-vertex, never dereferencing freed memory.
    Node *prevFarNode = nts2dPrevFarNode();
    Node *nextFarNode = nts2dNextFarNode();

    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    const double xs[2] = { Xs(0) + us(0), Xs(1) + us(1) };

    if (nps == 2) {
        if (segNode[0] == 0 || segNode[1] == 0) return false;
        double X0[2], X1[2];
        {
            const Vector &Xa = segNode[0]->getCrds(); const Vector &ua = segNode[0]->getTrialDisp();
            const Vector &Xb = segNode[1]->getCrds(); const Vector &ub = segNode[1]->getTrialDisp();
            X0[0] = Xa(0) + ua(0); X0[1] = Xa(1) + ua(1);
            X1[0] = Xb(0) + ub(0); X1[1] = Xb(1) + ub(1);
        }
        double xi, g2, n2[2];
        int st = projectSegment2D(X0, X1, xs, nts2dSigma, nts2dLref, xi, g2, n2);
        if (st == SEG2D_DEGENERATE)
            return false;
        // ADR-85 T1b OPEN-END acceptance window (post-gate fix, measured). At an
        // interior (chained) side, out-of-bounds is a REFUSAL, never a clamp --
        // the T1a seam contract, untouched. But a chain's OPEN TERMINAL side has
        // no neighbour and no vertex pair, so the strict tolIn = 1e-9 slack
        // turns TILT-INDUCED projection drift of a boundary slave into a hard
        // force discontinuity: on the G-T1b compression patch the master's end
        // nodes settle ~5e-3 more than its interior node, the end slave's xi
        // drifts to -1.25e-7, both end pairs refuse in the SAME iterate, the top
        // block loses its only contact anchors and Newton explodes (du ~ 1e13,
        // singular block; probe-traced). Accept the projection within
        // NTS2D_END_SLACK parametric of a TERMINAL side only (prevFar == 0 <=>
        // the X0 side is an open end; nextFar == 0 <=> X1 is). The gap is then
        // the kernel's signed distance to the segment's INFINITE LINE and B is
        // still its EXACT first variation (the T1a cross-form identity is
        // algebraic in xi, not restricted to [0,1]), so residual and tangent
        // stay consistent and C0 across the boundary band. Interior seams keep
        // the strict predicate, so the T1a ownership/uniqueness contract is
        // untouched; the honest permanent end treatment (a radial END-CAP
        // vertex, C0 with no window) is deferred by name -- see LEDGER_quirks.
        if (st == SEG2D_OUT_LOW &&
            !(prevFarNode == 0 && xi >= -NTS2D_END_SLACK))
            return false;
        if (st == SEG2D_OUT_HIGH &&
            !(nextFarNode == 0 && xi <= 1.0 + NTS2D_END_SLACK))
            return false;
        // Ordered-ownership step 1: the predecessor owns whenever IT is in-bounds
        // too. Its geometry is read at the SAME trial config; the in-bounds
        // status is sigma-independent, so passing the surface sigma is exact.
        if (prevFarNode != 0) {
            const Vector &Xp = prevFarNode->getCrds();
            const Vector &up = prevFarNode->getTrialDisp();
            const double XP[2] = { Xp(0) + up(0), Xp(1) + up(1) };
            double xiP, gP, nP[2];
            if (projectSegment2D(XP, X0, xs, nts2dSigma, nts2dLref, xiP, gP, nP)
                    == SEG2D_IN_BOUNDS)
                return false;                   // the previous segment owns this slave
        }
        if (g2 >= 0.0) return false;            // owned but separated
        gap = g2; n[0] = n2[0]; n[1] = n2[1];
        shape2D(xi, N);
        bOperatorSegment2D(N[0], N[1], n2, B);
        // ADR-85 T2 -- the ONE scalar tangential relative-DISPLACEMENT (never a
        // 2-vector differenced componentwise -- see the file-level note in
        // LadrunoFrictionKernel.h). t_hat = perp(n) (ADR-85 How/4); drel is the
        // displacement-not-position residual segmentActive's gTvec also uses,
        // just projected to a scalar instead of a tangent-plane 3-vector.
        if (gTout != 0) {
            const Vector &u0 = segNode[0]->getTrialDisp();
            const Vector &u1 = segNode[1]->getTrialDisp();
            const double ubar[2] = { N[0]*u0(0) + N[1]*u1(0), N[0]*u0(1) + N[1]*u1(1) };
            const double drel[2] = { us(0) - ubar[0], us(1) - ubar[1] };
            const double th[2] = { -n2[1], n2[0] };
            *gTout = drel[0]*th[0] + drel[1]*th[1];
        }
        return true;
    }

    // nps == 1: the concave-vertex pair (degenerate-segment representation).
    if (nps != 1 || segNode[0] == 0 || prevFarNode == 0 || nextFarNode == 0)
        return false;
    double XP[2], XV[2], XN[2];
    {
        const Vector &Xp = prevFarNode->getCrds(); const Vector &up = prevFarNode->getTrialDisp();
        const Vector &Xv = segNode[0]->getCrds();  const Vector &uv = segNode[0]->getTrialDisp();
        const Vector &Xn = nextFarNode->getCrds(); const Vector &un = nextFarNode->getTrialDisp();
        XP[0] = Xp(0) + up(0); XP[1] = Xp(1) + up(1);
        XV[0] = Xv(0) + uv(0); XV[1] = Xv(1) + uv(1);
        XN[0] = Xn(0) + un(0); XN[1] = Xn(1) + un(1);
    }
    // Ownership steps 1-2: either adjacent segment in-bounds => a SEGMENT owns.
    {
        double xiT, gT, nT[2];
        if (projectSegment2D(XP, XV, xs, nts2dSigma, nts2dLref, xiT, gT, nT) == SEG2D_IN_BOUNDS)
            return false;
        if (projectSegment2D(XV, XN, xs, nts2dSigma, nts2dLref, xiT, gT, nT) == SEG2D_IN_BOUNDS)
            return false;
    }
    const double tPrev[2] = { XV[0] - XP[0], XV[1] - XP[1] };
    const double tNext[2] = { XN[0] - XV[0], XN[1] - XV[1] };
    WedgeResult w = vertexWedge2D(tPrev, tNext, nts2dSigma, xs, XV, nts2dLref);
    // Step 3, UNSLACKED (this -- not the tolIn-slacked w.inWedge -- is what closes
    // the 5.2-ulp seam band the T1a oracle measured; kernel contract, verbatim).
    if (!(w.aPrev > 0.0 && w.aNext < 0.0))
        return false;
    if (liveSideOut != 0) *liveSideOut = w.sideSign;
    const double side = (committedSide != 0.0) ? committedSide : w.sideSign;
    if (side == 0.0) {
        // AMBIGUOUS (fold-back spike / conditioning gate): DEFER the capture --
        // never guess (kernel handler-flow step 2). Loud once per contact, the
        // WARN_EDGE_SIGN_DEFER precedent.
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VTX2D_DEFER))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": 2D vertex-pair "
                      "side-sign capture DEFERRED -- the corner is (numerically) a fold-back "
                      "spike or the slave sits on the wedge-boundary conditioning gate, so the "
                      "side sign is undecidable; the pair stays inert until the geometry "
                      "disambiguates (ADR-85 How/1).\n";
        return false;
    }
    double g2, n2[2];
    // tolLen = tauPerp*Lref, NOT the kernel comment's suggested tauSeg*Lref
    // (post-gate fix, measured): tauSeg = 1e-8 is the ZERO-LENGTH-SEGMENT gauge,
    // and tauSeg*Lref (1.118e-8 on the G-T1b notch) sits ABOVE the fork's
    // standard 1e-8 seeded penetration -- the corner deck's vertex pair refused
    // its own seed and the slave's y DOF went singular (probe: 1.117e-8 fails,
    // 1.119e-8 passes with the exact vertex answer). The floor's job here is
    // DIRECTION conditioning (r/||r|| of a coincident pair), for which the
    // dimensionless conditioning gauge tauPerp = 1e-12 is the right scale:
    // still RELATIVE to Lref (the ADR-57 P5 unit-trap rule), direction noise at
    // the floor ~eps*|x|/(1e-12*L) ~ 1e-4 rad, and the deferred force there is
    // O(kn*1e-12*L) -- nil. See LEDGER_quirks (ADR-85 T1b gauge collision).
    if (vertexEval2D(XV, xs, side, g2, n2, TAU_PERP_DEFAULT * nts2dLref) != VTX2D_OK)
        return false;
    if (g2 >= 0.0) return false;                // claimed but separated (convex-side graze)
    gap = g2; n[0] = n2[0]; n[1] = n2[1];
    N[0] = 1.0; N[1] = 0.0;                     // |B| weight of the one vertex node (B = [n | -n])
    bOperatorVertex2D(n2, B);
    // ADR-85 T2 -- the vertex pair's scalar tangential slip: the relative
    // DISPLACEMENT between the slave and the ONE vertex node, projected onto
    // t_hat = perp(n) (same construction as the segment branch above; the vertex
    // pair's B = [n | -n] has only one "master" node, so ubar collapses to uv).
    if (gTout != 0) {
        const Vector &uv = segNode[0]->getTrialDisp();
        const double drel[2] = { us(0) - uv(0), us(1) - uv(1) };
        const double th[2] = { -n2[1], n2[0] };
        *gTout = drel[0]*th[0] + drel[1]*th[1];
    }
    return true;
}

// ADR-85 T1b -- the domain-committed vertex side sign (0 when not a 2D vertex
// pair, no engine, or not yet captured; the commit itself is getResidual's).
double
LadrunoContactFE::vtx2DCommittedSide(void) const
{
    if (mode != SEGMENT || ndm != 2 || nps != 1) return 0.0;
    if (theDomain == 0 || theSlave == 0 || segNode[0] == 0) return 0.0;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) return 0.0;
    return cd->getVtx2DSide(contactTag, theSlave->getTag(), segNode[0]->getTag());
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

// ADR-85 T2 -- the 2D sibling of addFrictionTang (WORK ITEM 3). The kernel's
// tangentBlock1D returns two SCALARS (Kss_scalar on the tangential direction,
// dTN_scalar the consistent pressure-coupling coefficient); this method
// reconstructs the 2x2 physical block via the SAME chain rule addFrictionTang's
// 3x3 Kss embeds: d(gTeff)/d(u_rel) = t_hat (this pair's tangential B-operator,
// identical to the projection segment2DActive's gTout uses) and d(tn)/d(u_rel) =
// kn*n (the pair's NORMAL B-operator, i.e. n itself for a unit-weighted node) --
// so Kss_2x2 = Kss_scalar*(t_hat⊗t_hat) + dTN_scalar*(t_hat⊗n) -- dTN_scalar ALREADY
// carries kn (see tangentBlock1D), so kn is NOT reapplied here -- then scattered
// via the SAME G = [I_2 | −N_i I_2] discipline as the 3D block.
void
LadrunoContactFE::addFrictionTang2D(double fact, const double n[3], const double N2[2],
                                    double tn, double gTeff, double gpT, bool consistent)
{
    const double th[2] = { -n[1], n[0] };             // t_hat = perp(n), ADR-85 How/4
    double Kss_s = 0.0, dTN_s = 0.0;
    LadrunoFrictionKernel::tangentBlock1D(gTeff, gpT, tn, kn, kt, mu, consistent, Kss_s, dTN_s);
    double Kss2[2][2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            // ORCHESTRATOR REVIEW FIX (ADR-85 T2): dTN_s ALREADY carries kn --
            // tangentBlock1D returns dTN = -(dcap/dN)*kn*sgn, exactly as the 3D
            // frictionTangentBlock carries it ONCE in `(-dCap_dN*kn*nh[i])*n[j]`.
            // Multiplying by kn again here scaled the consistent cross term as
            // kn^2; with kn ~ 1e6..1e9 that is a catastrophically wrong
            // -consistanttan tangent. (The symmetric DEFAULT was unaffected:
            // dTN_s is identically 0 unless consistent==true.) kn appears ONCE.
            Kss2[i][j] = Kss_s * th[i] * th[j] + dTN_s * th[i] * n[j];
    double w[3];
    w[0] = 1.0;
    for (int i = 0; i < nps; i++) w[1 + i] = -N2[i];
    int nn = 1 + nps;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < nn; b++) {
            double wab = fact * w[a] * w[b];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    tang(2 * a + i, 2 * b + j) += wab * Kss2[i][j];
        }
}

// B3 (P2b-2c) — consistent ∂n/∂u geometric NORMAL tangent. Gathers the slave + master
// CURRENT positions (exactly like segmentActive), calls the oracle-validated pure
// function LadrunoContactProjection::normalGeomTangent (which re-projects deterministically
// ⇒ the SAME ξ̄/n/gap the residual used), and scatters K_geom = kn·gN·∂²gN/∂u² into `tang`.
// Inert unless penetrating + in-bounds; for a flat facet the slave block is identically 0
// (the byte-identity contract). SYMMETRIC ⇒ no -consistanttan / non-sym solver needed.
void
LadrunoContactFE::addNormalGeomTang(double fact)
{
    if (mode != SEGMENT || theSlave == 0) return;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    double Kg[15][15];
    if (!LadrunoContactProjection::normalGeomTangent(nps, Xseg, xs, orientDir, kn, Kg))
        return;                                 // not penetrating / out-of-bounds
    int ndof = 3 * (1 + nps);
    for (int a = 0; a < ndof; a++)
        for (int b = 0; b < ndof; b++)
            tang(a, b) += fact * Kg[a][b];
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

// ADR-57 E2 — EDGE_EDGE geometry at the current trial config. Gathers the 4 edge-node positions
// (X+u), runs the shipped LadrunoEdgeKernel closest-point, and (only for an EE_OK, MARGIN-INTERIOR
// crossing) fills the signed gap gN, the oriented common-perpendicular normal n, the parameters s,t,
// and the B-operator rows. A parallel / zero-length / near-vertex / clamped pair ⇒ returns false
// (no force this eval — routing demoted it). The sign is BODY-FIXED: `committedSign` (±1) is applied
// verbatim; 0 ⇒ orient from orientDir and return the chosen sign in *outSign for the caller to COMMIT.
bool
LadrunoContactFE::edgeGeom(double &gN, double n[3], double &s, double &t, double B[4][3],
                           int committedSign, int *outSign) const
{
    if (mode != EDGE_EDGE) return false;
    double X[4][3];
    for (int i = 0; i < 4; i++) {
        if (edgeNode[i] == 0) return false;
        const Vector &Xi = edgeNode[i]->getCrds();
        const Vector &ui = edgeNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) X[i][d] = Xi(d) + ui(d);
    }
    LadrunoEdgeKernel::ClosestResult cr =
        LadrunoEdgeKernel::closestPtSegSeg(X[0], X[1], X[2], X[3]);
    // only a well-conditioned, MARGIN-INTERIOR crossing is a true edge-edge contact where the
    // cross-product normal is the connector (gN = ±‖w‖). Parallel (EE_PARALLEL), zero-length
    // (EE_DEGENERATE), or near-vertex (!interior) ⇒ inert (routing handles those, ADR §1/§2).
    if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) return false;
    s = cr.s; t = cr.t;
    double d1[3], d2[3];
    for (int d = 0; d < 3; d++) { d1[d] = X[1][d] - X[0][d]; d2[d] = X[3][d] - X[2][d]; }
    // contact-review P5 (H2-style conditioning guard): a first-engagement sign capture with the
    // reference orientDir nearly PERPENDICULAR to n is a numerical coin flip — and the committed
    // body-fixed sign then binds the pair for its whole life. edgeGap returns chosen==0 in that
    // case (and on a degenerate normal): DEFER — treat this eval as no-contact and retry the
    // capture at the next, better-conditioned config. A committed ±1 sign always passes through.
    // The defer is LOUD (gate lens-A #6): on a symmetric rig the reference can stay ⟂ n on EVERY
    // eval (a permanently inert pair would be silent pass-through — worse than the old coin flip),
    // so surface it once per contact and point at the knob that decides it (-outward).
    int chosen = 0;
    gN = LadrunoEdgeKernel::edgeGap(cr.c1, cr.c2, d1, d2, committedSign, orientDir, n, &chosen);
    if (chosen == 0) {
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_EDGE_SIGN_DEFER))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": edge-edge sign "
                      "capture DEFERRED — the orientation reference is (numerically) perpendicular "
                      "to the edge cross-normal, so the body-fixed contact sign is undecidable; the "
                      "pair stays inert until the geometry disambiguates. If this persists, pass an "
                      "-outward with a component along the expected contact normal.\n";
        return false;
    }
    if (outSign != 0) *outSign = chosen;
    LadrunoEdgeKernel::bOperator(n, s, t, B);
    return true;
}

// ADR-57 E3 — tangential slip at the closest point, from DISPLACEMENT (the C3.1 trap). The closest-
// point construction makes the relative POSITION purely normal, so the tangential slip is the
// relative DISPLACEMENT of the two closest points, tangential part: gT = tangentPart((1−s)u_a0 +
// s u_a1 − (1−t)u_b0 − t u_b1, n). Same {s,t,n} edgeGeom returns ⇒ one consistent contact point.
void
LadrunoContactFE::edgeSlip(double s, double t, const double n[3], double gT[3]) const
{
    const Vector &ua0 = edgeNode[0]->getTrialDisp();
    const Vector &ua1 = edgeNode[1]->getTrialDisp();
    const Vector &ub0 = edgeNode[2]->getTrialDisp();
    const Vector &ub1 = edgeNode[3]->getTrialDisp();
    double drel[3];
    for (int d = 0; d < 3; d++)
        drel[d] = ((1.0 - s) * ua0(d) + s * ua1(d)) - ((1.0 - t) * ub0(d) + t * ub1(d));
    LadrunoFrictionKernel::tangentPart(drel, n, gT);
}

// B1 (P4) — the ASSEMBLED translational mass [mx,my,mz] of a node: the engine's pre-computed cache
// (nodal `mass` + Σ_elements diag(M_e) — what the explicit integrator actually inverts) if present,
// else the bare nodal `mass` (Node::getMass diagonal) as a fallback. The handler fills the cache once
// per handle() whenever a SOFT contact exists; without it (e.g. an unexpected call) we still get the
// nodal mass — correct for lumped-mass models, missing only the element-density contribution.
static void
ladrunoNodeMass(Domain *dom, Node *nd, double m[3])
{
    m[0] = m[1] = m[2] = 0.0;
    if (nd == 0) return;
    if (dom != 0) {
        LadrunoContactDomain *cd = dom->getLadrunoContactDomain();
        if (cd != 0 && cd->getNodalMass(nd->getTag(), m)) return;   // assembled (incl. element mass)
    }
    const Matrix &M = nd->getMass();
    int nr = M.noRows();
    for (int d = 0; d < 3 && d < nr; d++) m[d] = M(d, d);
}

// inverse mass PROJECTED on the gap normal n: invMproj = Σ_d n_d²/m_d over the 3 translational DOFs.
// A DOF with no mass (m_d ≤ 0) is treated as INFINITE mass (zero contribution) — the correct gap-mode
// treatment of a FIXED master node (a Dirichlet BC carries no inertia). So a fixed/rigid master
// contributes 0, leaving m_eff = m_slave; a fully massless SLAVE returns 0 ⇒ the caller cannot
// soft-size (handled there).
static double
ladrunoInvMassProj(const double m[3], const double n[3])
{
    double s = 0.0;
    for (int d = 0; d < 3; d++)
        if (m[d] > 0.0) s += n[d] * n[d] / m[d];
    return s;
}

// B1 (P4) — SOFT=1 Courant-stable penalty. The gap-mode generalized mass m_eff = 1/(B M⁻¹ Bᵀ),
// B = [n | −N_i n], reduces (diagonal lumped mass) to m_eff = 1/(invMproj_slave + Σ N_i²·invMproj_i)
// — the harmonic combination of the slave mass and the shape-projected segment mass; for a rigid
// plane (N==0) it collapses to the slave mass. Then the contact's own central-difference stability
// bound ω·dt ≤ 2 (ω = √(k/m_eff)) gives the LARGEST non-throttling penalty 4·m_eff/dt²; SOFT picks
// k_soft = SOFSCL·4·m_eff/dt² (SOFSCL < 1 the safety margin). Validated: proto_b1_soft_penalty.py.
// Only active under the explicit CentralDifferenceLadruno (dynamic_cast — catches its SMS subclasses
// too); any implicit integrator fails the cast ⇒ the configured kn is returned ⇒ byte-identical.
// B1 — gap-mode inverse effective mass B M⁻¹ Bᵀ for unit direction `dir`, B = [dir | −N_i dir] over
// [slave | seg nodes]; N==0 ⇒ RIGID_PLANE (slave only). A DOF with no mass ⇒ ∞ mass ⇒ 0 contribution
// (a FIXED master node drops out, leaving m_eff = m_slave). ≤0 ⇒ massless slave (caller cannot size).
double
LadrunoContactFE::gapModeInvMass(const double dir[3], const double N[4]) const
{
    double ms[3]; ladrunoNodeMass(theDomain, theSlave, ms);
    double invSum = ladrunoInvMassProj(ms, dir);  // 0 ⇒ massless slave
    if (invSum > 0.0 && N != 0)                    // SEGMENT: add the shape-projected segment masses
        for (int i = 0; i < nps; i++) {
            double mi[3]; ladrunoNodeMass(theDomain, segNode[i], mi);
            invSum += N[i] * N[i] * ladrunoInvMassProj(mi, dir);
        }
    return invSum;
}

// contact-review P5 — is the -visc dashpot live this eval? The old "statics-inert" claim
// relied on v ≡ 0 under a StaticIntegrator, but trial velocities are STATE, not rates: a
// static stage AFTER a transient one (or setNodeVel) leaves them non-zero, so the dashpot
// injected an unphysical velocity-proportional force with NO matching tangent
// (StaticIntegrator::formEleTangent never calls addCtoTang). Under a StaticIntegrator the
// dashpot is now OFF (warn once per contact via the engine latch); genuinely-zero
// velocities made the old force exactly 0, so gating it is byte-identical there.
bool
LadrunoContactFE::viscousActive(Integrator *theIntegrator) const
{
    if (muc <= 0.0) return false;
    if (dynamic_cast<StaticIntegrator *>(theIntegrator) == 0) return true;   // transient ⇒ live
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    if (cd == 0 || cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VISC_STATIC))
        opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": -visc (viscous "
                  "contact stabilization) is DISABLED under a static integrator — nodal "
                  "velocities are stale state (not rates) in statics and the damping tangent "
                  "is never assembled there; the dashpot re-arms under a transient integrator.\n";
    return false;
}

double
LadrunoContactFE::softKn(Integrator *theIntegrator, const double n[3], const double N[4]) const
{
    if (softScale <= 0.0 || theIntegrator == 0)
        return kn;
    CentralDifferenceLadruno *cdl = dynamic_cast<CentralDifferenceLadruno *>(theIntegrator);
    if (cdl == 0)
        return kn;                                // implicit / non-CDL ⇒ inert (byte-identical)
    double dt = cdl->getCurrentDeltaT();
    double invSum = gapModeInvMass(n, N);         // 0 ⇒ massless slave (fail below)
    if (dt <= 0.0 || invSum <= 0.0) {            // massless slave or no dt ⇒ cannot soft-size
        // contact-review P5: per-(contact, topic) engine latch, not a process-lifetime static
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (cd == 0 || cd->warnOnce(contactTag, LadrunoContactDomain::WARN_SOFT_NOMASS))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": -soft could not "
                      "size a penalty under the explicit integrator (non-positive dt or a contact "
                      "node with no assembled mass — neither a nodal `mass` nor element-density "
                      "mass); using the configured kn. Give the contact nodes mass for the SOFT=1 "
                      "stable penalty.\n";
        return kn;
    }
    return softScale * 4.0 * (1.0 / invSum) / (dt * dt);   // k_soft = SOFSCL·4·m_eff/dt²
}

// B1 (kt follow-up) — the SOFT=1 Courant-stable TANGENTIAL (Coulomb stick) penalty. softKn sizes only
// the NORMAL kn; under explicit a stiff friction kt still throttles dt_cr via the stick mode
// ω_t = √(kt/m_eff_t) (continuous stick at a high cone DIVERGES when ω_t·dt > 2). softKt sizes
// k_soft_t = softScale·4·m_eff_t/dt² (⇒ ω_t·dt = 2√softScale ≤ 2) from the WORST-CASE (largest inverse
// mass ⇒ smallest m_eff) over the two BASIS tangents t1,t2 to n. EXACT for the isotropic lumped mass
// `system Diagonal` uses (invMproj direction-independent ⇒ m_eff_t == m_eff_n, off-diagonal of the
// tangent-plane mass restriction vanishes); for a genuinely ANISOTROPIC nodal mass the binding direction
// can be a t1/t2 combination, so max(t1,t2) is an approximation (use isotropic lumped mass under explicit)
// — documented, like the B1 row-sum-SOE caveat. The per-node bound is also necessary-not-sufficient under
// multi-node coupling (inherited from softKn); the default SOFSCL=0.10 leaves the absorbing margin. Same
// gate as softKn (soft off / implicit / no dt-or-mass ⇒ the configured kt ⇒ byte-identical; softKn,
// called first in getResidual, already emits the no-mass warning). SEGMENT (NTS) friction only.
double
LadrunoContactFE::softKt(Integrator *theIntegrator, const double n[3], const double N[4]) const
{
    if (softScale <= 0.0 || theIntegrator == 0)
        return kt;
    CentralDifferenceLadruno *cdl = dynamic_cast<CentralDifferenceLadruno *>(theIntegrator);
    if (cdl == 0)
        return kt;                                // implicit / non-CDL ⇒ inert (byte-identical)
    double dt = cdl->getCurrentDeltaT();
    // ADR-85 T1b -- the explicit 2D branch (How/6): the 2D tangent space is ONE-
    // dimensional, t = perp(n) -- EXACT per-mode. The 3D least-aligned-axis
    // construction below would pick the OUT-OF-PLANE axis as t1 for an in-plane
    // n (the wrong space entirely); it stays textually untouched for ndm == 3.
    // ADR-85 T2 (WORK ITEM 4, verified): this branch is now LIVE -- T2 wires 2D
    // friction through getResidual/addKtToTang/addKiToTang (LadrunoContactFE::
    // returnMap1D call sites), and softKt is called only from that friction
    // path, so a `-soft`+`-mu` 2D deck now REACHES this branch under the
    // explicit CentralDifferenceLadruno exactly as coded. The B2 coupled-K_c
    // anisotropic-mass caveat (m_x != m_y makes the per-mode bound necessary-
    // not-sufficient) transfers to 2D verbatim.
    if (ndm == 2) {
        double t1[3] = { -n[1], n[0], 0.0 };      // perp(n), z-padded for the mass helpers
        double invT = gapModeInvMass(t1, N);
        if (dt <= 0.0 || invT <= 0.0)
            return kt;                            // massless / no dt => the configured kt
        return softScale * 4.0 * (1.0 / invT) / (dt * dt);  // k_soft_t = SOFSCL*4*m_eff_t/dt^2
    }
    // two orthonormal tangents to n: t1 ⟂ n built from the coordinate axis least aligned with n, t2=n×t1
    int k = 0; double amin = std::fabs(n[0]);
    if (std::fabs(n[1]) < amin) { amin = std::fabs(n[1]); k = 1; }
    if (std::fabs(n[2]) < amin) { k = 2; }
    double e[3] = {0.0, 0.0, 0.0}; e[k] = 1.0;
    double en = e[0]*n[0] + e[1]*n[1] + e[2]*n[2];
    double t1[3], nrm = 0.0;
    for (int d = 0; d < 3; d++) { t1[d] = e[d] - en * n[d]; nrm += t1[d]*t1[d]; }
    nrm = std::sqrt(nrm);
    if (nrm < 1e-300) return kt;                  // degenerate (n not unit) ⇒ configured kt
    for (int d = 0; d < 3; d++) t1[d] /= nrm;
    double t2[3] = { n[1]*t1[2] - n[2]*t1[1], n[2]*t1[0] - n[0]*t1[2], n[0]*t1[1] - n[1]*t1[0] };
    // worst case = MAX inverse-mass ⇒ MIN m_eff_t ⇒ HIGHEST stick frequency ⇒ the binding stability bound
    double inv1 = gapModeInvMass(t1, N), inv2 = gapModeInvMass(t2, N);
    double invMax = (inv1 > inv2) ? inv1 : inv2;
    if (dt <= 0.0 || invMax <= 0.0)
        return kt;                                // massless / no dt ⇒ the configured kt (softKn warned)
    return softScale * 4.0 * (1.0 / invMax) / (dt * dt);  // k_soft_t = SOFSCL·4·m_eff_t/dt²
}

// ADR-57 E5 — the 4-node edge gap-mode inverse effective mass B M⁻¹ Bᵀ for unit direction `dir`,
// B = [(1−s)dir, s dir, −(1−t)dir, −t dir] over the edge pair [a0,a1,b0,b1]. Closed form
// Σ w_i²·invMproj_i, w=[(1−s),s,(1−t),t] (the signs square away) — the B1 gapModeInvMass generalized
// from the [slave|seg] NTS operator to the 4-node edge operator. A fixed/massless node ⇒ ∞ mass ⇒ 0
// contribution (ladrunoInvMassProj); a fully massless pair ⇒ 0 (the caller cannot soft-size).
double
LadrunoContactFE::edgeGapModeInvMass(const double dir[3], double s, double t) const
{
    double w[4] = { 1.0 - s, s, 1.0 - t, t };     // |weights| (squared ⇒ sign-independent)
    double invSum = 0.0;
    for (int i = 0; i < 4; i++) {
        double mi[3]; ladrunoNodeMass(theDomain, edgeNode[i], mi);
        invSum += w[i] * w[i] * ladrunoInvMassProj(mi, dir);
    }
    return invSum;
}

// ADR-57 E5 — the edge-edge SOFT=1-analogue NORMAL penalty (the B1 softKn for the edge operator).
// Only active under the explicit CentralDifferenceLadruno (dynamic_cast — catches its SMS subclasses);
// any implicit integrator fails the cast ⇒ the configured kn is returned ⇒ byte-identical.
double
LadrunoContactFE::softKnEdge(Integrator *theIntegrator, const double n[3], double s, double t) const
{
    if (softScale <= 0.0 || theIntegrator == 0)
        return kn;
    CentralDifferenceLadruno *cdl = dynamic_cast<CentralDifferenceLadruno *>(theIntegrator);
    if (cdl == 0)
        return kn;                                // implicit / non-CDL ⇒ inert (byte-identical)
    double dt = cdl->getCurrentDeltaT();
    double invSum = edgeGapModeInvMass(n, s, t);  // 0 ⇒ massless pair (fail below)
    if (dt <= 0.0 || invSum <= 0.0) {            // massless / no dt ⇒ cannot soft-size
        // contact-review P5: per-(contact, topic) engine latch, not a process-lifetime static
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (cd == 0 || cd->warnOnce(contactTag, LadrunoContactDomain::WARN_EDGESOFT_NOMASS))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": edge-edge -edgeSoft "
                      "could not size a penalty under the explicit integrator (non-positive dt or an "
                      "edge node with no assembled mass — neither a nodal `mass` nor element-density "
                      "mass); using the configured edgeKn. Give the edge nodes mass for the SOFT "
                      "stable penalty.\n";
        return kn;
    }
    return softScale * 4.0 * (1.0 / invSum) / (dt * dt);   // k_soft = SOFSCL·4·m_eff/dt²
}

// ADR-57 E5 — the edge-edge SOFT TANGENTIAL (Coulomb stick) penalty: the B1-kt n→t rule on the edge
// operator. Sizes k_soft_t = softScale·4·m_eff_t/dt² from the WORST-CASE (largest inverse mass ⇒
// smallest m_eff) over the two basis tangents t1,t2 to n ⇒ ω_t·dt = 2√softScale ≤ 2 (kt never
// throttles dt_cr via the stick mode). Same gate as softKnEdge (soft off / implicit / no dt-or-mass
// ⇒ the configured kt ⇒ byte-identical; softKnEdge, called first in getResidual, emits any no-mass warning).
double
LadrunoContactFE::softKtEdge(Integrator *theIntegrator, const double n[3], double s, double t) const
{
    if (softScale <= 0.0 || theIntegrator == 0)
        return kt;
    CentralDifferenceLadruno *cdl = dynamic_cast<CentralDifferenceLadruno *>(theIntegrator);
    if (cdl == 0)
        return kt;                                // implicit / non-CDL ⇒ inert (byte-identical)
    double dt = cdl->getCurrentDeltaT();
    // two orthonormal tangents to n: t1 ⟂ n from the coordinate axis least aligned with n, t2=n×t1
    int k = 0; double amin = std::fabs(n[0]);
    if (std::fabs(n[1]) < amin) { amin = std::fabs(n[1]); k = 1; }
    if (std::fabs(n[2]) < amin) { k = 2; }
    double e[3] = {0.0, 0.0, 0.0}; e[k] = 1.0;
    double en = e[0]*n[0] + e[1]*n[1] + e[2]*n[2];
    double t1[3], nrm = 0.0;
    for (int d = 0; d < 3; d++) { t1[d] = e[d] - en * n[d]; nrm += t1[d]*t1[d]; }
    nrm = std::sqrt(nrm);
    if (nrm < 1e-300) return kt;                  // degenerate (n not unit) ⇒ configured kt
    for (int d = 0; d < 3; d++) t1[d] /= nrm;
    double t2[3] = { n[1]*t1[2] - n[2]*t1[1], n[2]*t1[0] - n[0]*t1[2], n[0]*t1[1] - n[1]*t1[0] };
    double inv1 = edgeGapModeInvMass(t1, s, t), inv2 = edgeGapModeInvMass(t2, s, t);
    double invMax = (inv1 > inv2) ? inv1 : inv2;  // worst case = MIN m_eff_t = the binding bound
    if (dt <= 0.0 || invMax <= 0.0)
        return kt;                                // massless / no dt ⇒ configured kt (softKnEdge warned)
    return softScale * 4.0 * (1.0 / invMax) / (dt * dt);  // k_soft_t = SOFSCL·4·m_eff_t/dt²
}

const Vector &
LadrunoContactFE::getResidual(Integrator *theIntegrator)
{
    resid.Zero();
    if (mode == RIGID_PLANE) {
        double g = rigidPlaneGap();
        if (g < 0.0) {                       // penetration -> restoring force +n
            // B1: under explicit, kn → the Courant-stable SOFT=1 k_soft (rigid plane ⇒ m_eff = m_s);
            // otherwise the configured kn (byte-identical). N==0 ⇒ slave-only effective mass.
            double knEff = softKn(theIntegrator, planeN, 0);
            double tn = -knEff * g;          // tn = kn*|g| > 0 (Macaulay <-g>_+)
            for (int d = 0; d < ndm; d++)
                resid(d) = tn * planeN[d];   // drives the slave toward g=0 (PenaltySP convention)
            // D2 viscous normal stabilization: f_visc = −muc·(n·v_slave)·n, opposing the approach
            // rate ġ = n·v_slave. Active only while in contact (g<0) and only under a TRANSIENT
            // integrator (contact-review P5: trial velocities are stale STATE under statics — see
            // viscousActive). Force here; the muc·n⊗n damping tangent is in addCtoTang.
            if (theSlave != 0 && viscousActive(theIntegrator)) {
                const Vector &vs = theSlave->getTrialVel();
                double gdot = 0.0;
                for (int d = 0; d < ndm; d++) gdot += planeN[d] * vs(d);
                for (int d = 0; d < ndm; d++) resid(d) += -muc * gdot * planeN[d];
            }
        }
    } else if (mode == SEGMENT && ndm == 2) {
        // ADR-85 T1b/T2 -- the 2D NTS pair (segment nps == 2 / concave vertex
        // nps == 1). mu<=0 is frictionless (byte-identical to the T1b lane,
        // the SAME short-circuit discipline the 3D SEGMENT branch below uses);
        // mu>0 wires the SCALAR return map (ADR-85 T2, below). The shipped 3D
        // SEGMENT branch below is textually untouched -- an ndm == 2 adapter
        // cannot reach it, and every 3D adapter carries ndm == 3 so it skips
        // this arm on one comparison.
        double gap, n[3], N2[2], B[6], liveSide = 0.0, gTnow = 0.0;
        const double committedSide = vtx2DCommittedSide();
        bool act2d = segment2DActive(gap, n, N2, B, committedSide, &liveSide, &gTnow);
        // panel MINOR diagnostic: a committed side sign that DISAGREES with the
        // live wedge classification means the corner has deformed past its
        // reference-config type -- a reference-concave corner gone convex would
        // silently ATTRACT under the committed sign. The committed sign is kept
        // for the pairing epoch (the How/1 contract); this latch makes the
        // configuration loud instead of silent. liveSide is nonzero only when
        // the wedge CLAIMED this eval, so separated pairs never warn.
        if (nps == 1 && committedSide != 0.0 && liveSide != 0.0 &&
            liveSide != committedSide && theDomain != 0) {
            LadrunoContactDomain *cdS = theDomain->getLadrunoContactDomain();
            if (cdS != 0 && cdS->warnOnce(contactTag, LadrunoContactDomain::WARN_VTX2D_SIDE_FLIP))
                opserr << "WARNING LadrunoContactFE - contact " << contactTag
                       << ": 2D vertex-pair COMMITTED side sign disagrees with the LIVE "
                          "wedge classification at vertex node " << segNode[0]->getTag()
                       << " (committed " << committedSide << ", live " << liveSide
                       << ") -- the corner has deformed past its reference-config type "
                          "(a reference-concave corner gone convex would silently attract "
                          "under the committed sign). Inspect the deck; the committed sign "
                          "is kept for this pairing epoch (ADR-85 How/1).\n";
        }
        if (act2d) {
            // FIRST-CAPTURE commit of the vertex side sign (kernel handler-flow
            // step 3: committed at first capture, passed verbatim thereafter --
            // a per-step re-derived sign flips exactly while interpenetrating).
            // For a CONCAVE vertex the captured value is structurally -1
            // (sideSign = -corner inside the wedge), so a capture at a later-
            // rejected implicit trial config commits the same value a clean
            // capture would -- no double-buffer needed (see LadrunoContactDomain).
            if (nps == 1 && committedSide == 0.0 && liveSide != 0.0 && theDomain != 0) {
                LadrunoContactDomain *cdV = theDomain->getLadrunoContactDomain();
                if (cdV != 0)
                    cdV->setVtx2DSide(contactTag, theSlave->getTag(),
                                      segNode[0]->getTag(), liveSide);
            }
            // B1 SOFT=1 is FREE here (ADR-85 How/6): softKn/gapModeInvMass are
            // size-safe for 2-DOF nodes -- n is z-padded (n[2] = 0, so invMproj
            // sums only the in-plane components) and N carries the |B| weights
            // (segment: shape fns; vertex: {1, 0} since B = [n | -n]).
            double N4[4] = { N2[0], N2[1], 0.0, 0.0 };
            double knEff = softKn(theIntegrator, n, N4);
            double tn = LadrunoContactKernel::traction(knEff, gap);   // kn*<-gap>_+ > 0
            int ndof = 2 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;           // r = B^T tn (slave +tn n, master -N_i tn n)

            // NTS force snapshot (the `ladrunoContactForce` query is fed from
            // the SEGMENT branch only -- ADR-85 What). Vertex pairs report under
            // their nSeg+ordinal segIndex, so the per-slave sum stays complete.
            if (theDomain != 0) {
                LadrunoContactDomain *cdF = theDomain->getLadrunoContactDomain();
                if (cdF != 0)
                    cdF->setNtsForce(contactTag, theSlave->getTag(), segIndex, tn);
            }

            // D2 -visc dashpot -- the T1b SEGMENT port (ADR-85 How/6): the
            // shipped stride-3 loops become d < ndm over the 2D B layout
            // (B block of node i starts at ndm*(1+i) = 2*(1+i)); same active
            // set, same static-integrator gate (viscousActive), force-only
            // under the mass-only explicit CDL, force + addCtoTang under
            // implicit -- exactly the 3D contract.
            if (viscousActive(theIntegrator)) {
                double gdot = 0.0;
                const Vector &vs = theSlave->getTrialVel();
                for (int d = 0; d < 2; d++) gdot += B[d] * vs(d);
                for (int i = 0; i < nps; i++) {
                    const Vector &vm = segNode[i]->getTrialVel();
                    for (int d = 0; d < 2; d++) gdot += B[2 * (1 + i) + d] * vm(d);
                }
                for (int k = 0; k < ndof; k++) resid(k) += -muc * gdot * B[k];
            }

            // --- ADR-85 T2 -- 2D NTS Coulomb friction (SCALAR return map) ---
            // mu<=0 SHORT-CIRCUITS before any state touch ⇒ byte-identical to the
            // T1b frictionless 2D pair (the SAME `mu > 0.0` discipline the 3D
            // SEGMENT branch below uses -- no new silence). The engine is
            // re-fetched lazily (wipe deletes it) and null-checked, exactly like
            // the 3D path; the ONLY new state this touches is the EXISTING
            // FrictionState slot 0 (gpT[0]/gpTtrial[0]/gT0[0]/gT0committed[0] --
            // no new members, ADR-85 T2 binding).
            if (mu > 0.0 && theDomain != 0) {
                LadrunoContactDomain *cdFr = theDomain->getLadrunoContactDomain();
                if (cdFr != 0) {
                    int slaveTag = theSlave->getTag();
                    LadrunoContactDomain::FrictionState &st =
                        cdFr->getOrCreateFrictionState(contactTag, slaveTag, segIndex);
                    // capture the ENGAGEMENT-config tangential origin ONCE, as a
                    // SCALAR (gTnow is already the ONE t_hat-projected number
                    // segment2DActive formed from the relative DISPLACEMENT --
                    // never a 2-vector to difference componentwise; the ADR-85
                    // T2 design note's BINDING implementation consequence).
                    if (!st.engaged) { st.gT0[0] = gTnow; st.engaged = true; }
                    double gTeff = gTnow - st.gT0[0];      // scalar - scalar (Sterbenz-exact)
                    // B1 (kt): the SAME softKt this contact's normal block used to
                    // size ktEff (softKt already carries a NAMED 2D branch, How/6 --
                    // now REACHED for the first time, WORK ITEM 4).
                    double ktEff = softKt(theIntegrator, n, N4);
                    double tFricScalar, gpTtrial0;
                    LadrunoFrictionKernel::returnMap1D(gTeff, st.gpT[0], tn, ktEff, mu,
                                                       tFricScalar, gpTtrial0);
                    st.gpTtrial[0] = gpTtrial0;            // SET, not += (idempotent, BLOCKER-2)
                    // scatter tFricScalar*t_hat exactly like the normal block scatters
                    // tn*n (B = [n | -N_i n]); t_hat = perp(n), the SAME direction
                    // gTnow was projected onto.
                    const double th[2] = { -n[1], n[0] };
                    resid(0) += tFricScalar * th[0];
                    resid(1) += tFricScalar * th[1];
                    for (int i = 0; i < nps; i++) {
                        resid(2 * (1 + i) + 0) += -N2[i] * tFricScalar * th[0];
                        resid(2 * (1 + i) + 1) += -N2[i] * tFricScalar * th[1];
                    }
                }
            }
        }
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15], gTvec[3];   // ndof <= 3*(1+4) = 15
        if (segmentActive(gap, n, N, B, gTvec)) {
            // B1: under explicit, kn → the Courant-stable SOFT=1 k_soft (m_eff from the slave +
            // shape-projected segment masses on n); otherwise the configured kn (byte-identical).
            double knEff = softKn(theIntegrator, n, N);
            double tn = LadrunoContactKernel::traction(knEff, gap);  // = kn*<-gap>_+ > 0
            int ndof = 3 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;        // r = Bᵀ tn (slave +tn n, master −N_i tn n)

            // B3: report the normal force into the engine snapshot (the `ladrunoContactForce`
            // query). Pure side-channel — no effect on resid/tang. Overwrites this pair's slot.
            if (theDomain != 0) {
                LadrunoContactDomain *cdF = theDomain->getLadrunoContactDomain();
                if (cdF != 0)
                    cdF->setNtsForce(contactTag, theSlave->getTag(), segIndex, tn);
            }

            // --- D2 viscous normal stabilization (force; tangent in addCtoTang) ---
            // ġ = B·v (normal relative velocity, B=[n|−Nᵢ n]); f_visc = −muc·ġ·B opposes the
            // approach. Active only in contact (segmentActive) and only under a TRANSIENT
            // integrator (contact-review P5 — see viscousActive). Force-only under CDL; force +
            // muc·B Bᵀ tangent under implicit.
            if (viscousActive(theIntegrator)) {
                double gdot = 0.0;
                const Vector &vs = theSlave->getTrialVel();
                for (int d = 0; d < 3; d++) gdot += B[d] * vs(d);
                for (int i = 0; i < nps; i++) {
                    const Vector &vm = segNode[i]->getTrialVel();
                    for (int d = 0; d < 3; d++) gdot += B[3 * (1 + i) + d] * vm(d);
                }
                for (int k = 0; k < ndof; k++) resid(k) += -muc * gdot * B[k];
            }

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
                    // N for the cone = current penetration force tn (design MINOR-8). B1 (kt): under
                    // explicit, ktEff is the Courant-stable SOFT tangential stick penalty (so a stiff
                    // configured kt does not throttle dt_cr via the stick mode); otherwise the
                    // configured kt (byte-identical). The cone μ·tn already rides the SOFT normal tn.
                    double ktEff = softKt(theIntegrator, n, N);
                    LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, tn, ktEff, mu,
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
        // B2 (P5) — SOFT=2 segment-based EXPLICIT penalty. Under the explicit CentralDifferenceLadruno
        // ONLY (the same dynamic_cast gate as B1), replace the per-facet epsN with the Courant-stable
        // k_soft = SOFSCL·4·m_eff/dt² per slave node and assemble a pure single-pass penalty over the
        // clipped overlap (no ALM, no λ accumulation) — catching the corner/edge/T-intersection cases
        // the NTS SOFT=1 lane misses while keeping explicit dt_cr un-throttled. Under any implicit
        // integrator the cast fails ⇒ fall through to the shipped penalty/ALM path with the configured
        // epsN (SOFT=2 is explicit-only ⇒ implicit byte-identical to a plain -mortar penalty).
        // addSoft2Penalty also assembles the SOFT=2 viscous damper (-visc, D2.2) and Courant-stable
        // Coulomb/Tresca friction (-mu/-cohesion/-tauMax) on the same overlap; only -tie is refused.
        if (softScale > 0.0 &&
            dynamic_cast<CentralDifferenceLadruno *>(theIntegrator) != 0) {
            addSoft2Penalty(theIntegrator);
            return resid;
        }
        double D[4][4], M[4][4], g[4], n[3];
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (isTie) {
            // C4 — MESH-TIE: a permanent bond. The active set is frozen ON for every slave node and
            // the FULL 3-vec r_I = ΣD u_s − ΣM u_m is driven to zero (no KKT, no clamp, no friction).
            if (mortarActive(D, M, g, n)) {
                addMortarTieForce(D, M, cd);
            } else if (cd != 0) {
                // overlap empty this eval: zero this facet's r contribution so it stops biasing the
                // shared node's global accumulator (mirrors the normal-gap empty-overlap reset).
                double z[3] = {0.0, 0.0, 0.0};
                for (int I = 0; I < npsS; I++)
                    cd->accumulateMortarTie(contactTag, mortarSlave[I]->getTag(), this->getTag(),
                                            z, 0.0, kn);
            }
            return resid;
        }
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

            // --- D2.2 viscous normal stabilization (force; tangent in addCtoTang) ---
            // Per IN-CONTACT slave node (p[I]<0, the same KKT mask): the weighted normal gap RATE
            // ḡ̇_I = n·(Σ_J D_IJ v_s,J − Σ_K M_IK v_m,K)/a_I (v = getTrialVel), the viscous pressure
            // p_visc_I = μ_c·ḡ̇_I (NO clamp — a dashpot active while in contact), scattered EXACTLY like
            // the normal penalty force (f^s=−(D·p_visc)n, f^m=+(Mᵀ·p_visc)n). It is the C2 normal
            // operator with epsN→μ_c, driven by velocity. Transient integrators only (contact-review
            // P5 — see viscousActive; stale trial velocities are not rates under statics).
            if (viscousActive(theIntegrator)) {
                double vs[4][3], vm[4][3];
                for (int i = 0; i < npsS; i++) {
                    const Vector &v = mortarSlave[i]->getTrialVel();
                    for (int d = 0; d < 3; d++) vs[i][d] = v(d);
                }
                for (int i = 0; i < npsM; i++) {
                    const Vector &v = mortarMaster[i]->getTrialVel();
                    for (int d = 0; d < 3; d++) vm[i][d] = v(d);
                }
                double pv[4] = {0, 0, 0, 0};
                for (int I = 0; I < npsS; I++) {
                    if (p[I] >= 0.0) continue;            // viscous only on in-contact nodes
                    double aFacet = 0.0;
                    for (int J = 0; J < npsS; J++) aFacet += D[I][J];
                    if (aFacet <= 1e-300) continue;
                    double rdot[3] = {0, 0, 0};
                    for (int J = 0; J < npsS; J++)
                        for (int d = 0; d < 3; d++) rdot[d] += D[I][J] * vs[J][d];
                    for (int K = 0; K < npsM; K++)
                        for (int d = 0; d < 3; d++) rdot[d] -= M[I][K] * vm[K][d];
                    double gdot = (rdot[0]*n[0] + rdot[1]*n[1] + rdot[2]*n[2]) / aFacet;
                    pv[I] = muc * gdot;                   // p_visc_I = μ_c·ḡ̇_I
                }
                for (int K = 0; K < npsS; K++) {          // slave: −(D·p_visc)_K n
                    double Dp = 0.0;
                    for (int I = 0; I < npsS; I++) Dp += D[K][I] * pv[I];
                    for (int d = 0; d < 3; d++) resid(3 * K + d) += -Dp * n[d];
                }
                for (int L = 0; L < npsM; L++) {          // master: +(Mᵀ·p_visc)_L n
                    double Mp = 0.0;
                    for (int I = 0; I < npsS; I++) Mp += M[I][L] * pv[I];
                    for (int d = 0; d < 3; d++) resid(3 * (npsS + L) + d) += Mp * n[d];
                }
            }
        } else if (cd != 0) {
            // overlap empty THIS eval (the pair separated / the clip rejected it): zero this
            // facet's contribution so it stops biasing the shared node's accumulated global gap.
            for (int I = 0; I < npsS; I++)
                cd->accumulateMortarGap(contactTag, mortarSlave[I]->getTag(), this->getTag(),
                                        0.0, 0.0, kn);
        }
    } else if (mode == EDGE_EDGE) {
        // ADR-57 E2 — perpendicular edge-edge penalty. Run the closest-point geometry at the trial
        // config, fetch the per-pair body-fixed committed sign (Domain-owned EdgeEdgeState; the §2
        // A-4 fix so gN cannot flip through contact), and assemble the penalty force f = tN·B,
        // tN = εN⟨−gN⟩, B = [(1−s)n | s n | −(1−t)n | −t n]. Active ONLY in penetration (gN<0) and
        // only for a well-conditioned margin-interior crossing (edgeGeom refuses the rest). The
        // master edge gets the equal-and-opposite reaction (Σf = 0). E3 adds the Coulomb/Tresca
        // friction on top (SOFT/ALM are E5/E6).
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        LadrunoContactDomain::EdgeEdgeState *st = 0;
        int committedSign = 0;
        if (cd != 0) {
            st = &cd->getOrCreateEdgeEdgeState(contactTag, edgeNode[0]->getTag(),
                                               edgeNode[1]->getTag(), edgeNode[2]->getTag(),
                                               edgeNode[3]->getTag());
            committedSign = st->signN;
        }
        double gN, n[3], s, t, B[4][3]; int usedSign = 0;
        if (edgeGeom(gN, n, s, t, B, committedSign, &usedSign)) {
            // capture the body-fixed sign ONCE at first engagement, then hold it committed (the §2
            // A-4 fix — applied verbatim every eval, never re-derived from w·n which masks penetration).
            if (st != 0 && st->signN == 0) st->signN = usedSign;
            // E5 — under the explicit CDL with -edgeSoft, kn is REPLACED by the Courant-stable
            // k_soft = SOFSCL·4·m_eff/dt² (the edge gap-mode mass); softScale≤0 / implicit ⇒ kn
            // (byte-identical). Inert in addKtToTang (CDL never calls it; implicit knEff≡kn).
            double knEff = softKnEdge(theIntegrator, n, s, t);
            // E6 — the one-scalar ALM augmented pressure p = min(0, λ_N + εN·gN), traction tN = −p.
            // -edgeAlm OFF ⇒ λ_N≡0 ⇒ p = εN·gN, active iff gN<0 ⇒ the E2/E5 penalty path EXACTLY
            // (byte-identical). ON ⇒ inject the committed per-pair λ_N (the MortarNormalState C2.2
            // Uzawa, point-like) and stash gN + εN for the commit-time Uzawa + the penetration query
            // (committed-only ⇒ written here, mutated into λ_N only in Domain::commit()).
            double lambdaN = (edgeAlm && st != 0) ? st->lambdaN : 0.0;
            double pAug = lambdaN + knEff * gN;          // augmented pressure (≤0 ⇒ active)
            if (edgeAlm && st != 0) { st->gN_committed = gN; st->epsN = knEff; }
            if (pAug < 0.0) {                            // KKT-active (penetration, or λ_N-held contact)
                double tN = -pAug;                       // εN·|gN| (penalty) or |λ_N + εN·gN| (ALM) > 0
                for (int i = 0; i < 4; i++)
                    for (int d = 0; d < 3; d++)
                        resid(3 * i + d) = tN * B[i][d];  // f = tN·B (slave +, master − ⇒ Σf=0)

                // --- E3 Coulomb/Tresca friction (force; tangent is in addKtToTang) ---
                // mu≤0 ∧ c≤0 ∧ τmax≤0 SHORT-CIRCUITS before any slot touch ⇒ byte-identical to the
                // E2 frictionless path (the NTS `mu>0` guard — the kernel's cap≤0 branch returns the
                // RAW elastic, so byte-identity is the GUARD's job, not the kernel's). Needs the engine
                // for the per-pair committed slip (gpT) + engagement origin (gT0); cd==0 ⇒ frictionless.
                if (st != 0 && (mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0)) {
                    double gTvec[3]; edgeSlip(s, t, n, gTvec);
                    // capture the ENGAGEMENT-config tangential origin ONCE at first activation (else a
                    // late-engaging pair's pre-contact tangential drift is a spurious stick traction —
                    // the ADR-39 P3 MAJOR-1, reused).
                    if (!st->engaged) {
                        for (int d = 0; d < 3; d++) st->gT0[d] = gTvec[d];
                        st->engaged = true;
                    }
                    double gTeff[3];
                    for (int d = 0; d < 3; d++) gTeff[d] = gTvec[d] - st->gT0[d];
                    double tFric[3], gpTtrial[3];
                    // N for the cone = the current normal pressure tN; trial slip = PURE fn of committed
                    // gpT ⇒ idempotent across the CDL firstStep double-eval. tFric is already negated
                    // (motion-opposing). The 4th consumer of the ONE shipped return map. E5 — under the
                    // explicit CDL with -edgeSoft the STICK penalty kt is replaced by the Courant-stable
                    // k_soft_t (the B1-kt n→t rule) so a stiff kt never throttles dt_cr; implicit ⇒ kt.
                    double ktEff = softKtEdge(theIntegrator, n, s, t);
                    LadrunoFrictionKernel::frictionReturnMap(gTeff, st->gpT, tN, ktEff, mu, tFric,
                                                             gpTtrial, mortarCohesion, mortarTauMax);
                    for (int d = 0; d < 3; d++) st->gpTtrial[d] = gpTtrial[d];
                    // scatter via the SAME edge weights as the normal force (w = [(1−s),s,−(1−t),−t];
                    // Σw = 0 ⇒ Σf = 0). slave a0/a1 += (1−s)/s·tFric, master b0/b1 += −(1−t)/−t·tFric.
                    double w[4] = { 1.0 - s, s, -(1.0 - t), -t };
                    for (int i = 0; i < 4; i++)
                        for (int d = 0; d < 3; d++)
                            resid(3 * i + d) += w[i] * tFric[d];
                }
            }
        } else if (edgeAlm && st != 0) {
            // E6 — the pair is NOT a current margin-interior crossing (it slid off an endpoint or went
            // parallel — a LARGE-SLIDING case, outside the reference-config MVP scope §7). Drop the
            // STALE committed gap/penalty so the commit Uzawa neither augments λ_N from a dead
            // penetration nor leaves a live εN polluting getMaxEdgePenetration. λ_N is HELD
            // (min(0, λ_N + 0) = λ_N ≤ 0) until the pair re-crosses, then self-corrects (the mortar
            // accumulateMortarGap(0,0) / mortarNormalGCEnd reset pattern, point-like analogue).
            st->gN_committed = 0.0;
            st->epsN = 0.0;
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
        // C3.3 augmented-Lagrange: the OFFSET TRICK gTeff_eff = gTeff + λ_T/epsT injects the
        // committed tangential multiplier so the penalty return map yields the AUGMENTED trial
        // tT* = λ_T + epsT·(gTeff − gpT). λ_T≡0 ⇒ the C3.1/C3.2 penalty friction (epsT = kt).
        double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
        for (int d = 0; d < 3; d++) gTeff[d] = (gbarT[d] - st.gT0[d]) + st.lambdaT[d] * invEpsT;
        double tF[3], gpTtrial[3];
        // N for the cone = the nodal normal pressure; epsT rides kt; trial = pure fn of committed
        // gpT/λ_T ⇒ idempotent across re-evals. Returns the APPLIED (negated) traction opposing motion.
        LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, N_I, kt, mu, tF, gpTtrial,
                                                 mortarCohesion, mortarTauMax);
        // C3.3 Uzawa trial: λ_T ← −tFric (the returned cone-capped traction); committed in commit().
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = gpTtrial[d]; tFric[I][d] = tF[d]; st.lambdaTtrial[d] = -tF[d];
        }
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

// C4 — MESH-TIE force (a permanent bond; the zero-gap limit of contact). For every slave node I of
// the facet pair: build the FULL 3-vec weighted relative DISPLACEMENT r_I = Σ_J D_IJ u_s,J −
// Σ_K M_IK u_m,K (from getTrialDisp — the bond exists from the as-built config, so NO gT0; the C3.1
// displacement-not-position lesson), form t_I = λ_tie,I + epsTie·(r_I/a_I) (NO clamp — equality, all
// 3 components), and scatter via D/−M exactly like the NORMAL force: f^s_K = −Σ_I D_KI t_I,
// f^m_L = +Σ_I M_IL t_I (self-equilibrating, Σφ=1; oracle T2). The LOCAL r_I keeps R(u) deterministic
// per facet (the C2.2 rule); the GLOBAL r is accumulated for the commit Uzawa + ‖r‖ query only.
void
LadrunoContactFE::addMortarTieForce(const double D[4][4], const double M[4][4],
                                    LadrunoContactDomain *cd)
{
    double us[4][3], um[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) us[i][d] = u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) um[i][d] = u(d);
    }
    double t[4][3] = {{0}};
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;             // unreferenced slave node
        // r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K  (full 3-vec weighted relative displacement)
        double r[3] = {0, 0, 0};
        for (int J = 0; J < npsS; J++)
            for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
        for (int K = 0; K < npsM; K++)
            for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
        double lamTie[3] = {0, 0, 0};
        if (cd != 0) {
            int nodeTag = mortarSlave[I]->getTag();
            // accumulate the GLOBAL r (order-independent — the commit Uzawa + ‖r‖ query read it).
            cd->accumulateMortarTie(contactTag, nodeTag, this->getTag(), r, aFacet, kn);
            // the per-GLOBAL-node tie multiplier λ_tie (committed); it assembles globally for free
            // (Σ_facets D_KI λ_tie = D_KI^global λ_tie) — the C2.2 variationally-consistent term.
            const LadrunoContactDomain::MortarNormalState &st =
                cd->getOrCreateMortarNormalState(contactTag, nodeTag);
            for (int d = 0; d < 3; d++) lamTie[d] = st.lambdaTie[d];
        }
        // tie traction t_I = λ_tie,I + epsTie·(r_I^local/a_I) — full 3-vec, NO clamp (kn = epsTie).
        for (int d = 0; d < 3; d++) t[I][d] = lamTie[d] + kn * (r[d] / aFacet);
    }
    // scatter like the NORMAL force: f^s_K = −Σ_I D_KI t_I, f^m_L = +Σ_I M_IL t_I (the buffer is
    // zeroed at the top of getResidual and the tie owns the whole MORTAR residual ⇒ assign with =).
    for (int K = 0; K < npsS; K++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += D[K][I] * t[I][d];
            resid(3 * K + d) = -s;
        }
    for (int L = 0; L < npsM; L++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += M[I][L] * t[I][d];
            resid(3 * (npsS + L) + d) = s;
        }
}

// B2 (P5) — SOFT=2 segment-based EXPLICIT penalty force. Re-integrates the facet pair via the
// shipped clip→Gauss kernel (mortarActive ⇒ D,M,g̃,n over the overlap), then per slave node I sizes
// the Courant-stable penalty from the SEGMENT-BASED gap-mode generalized mass and assembles a pure
// single-pass penalty traction (NO ALM/λ — explicit). The gap operator for the area-normalised node-I
// gap ḡ_I = g̃_I/a_I is B_I : slave node J → (D_IJ/a_I) n, master node K → −(M_IK/a_I) n, so the
// gap-mode generalized mass is m_eff,I = 1/(B_I M⁻¹ B_Iᵀ) = 1/(Σ_J (D_IJ/a_I)² invMproj_sJ +
// Σ_K (M_IK/a_I)² invMproj_mK) from the ASSEMBLED nodal masses projected on n (a fixed/massless node
// ⇒ ∞ mass ⇒ 0 contribution). Then k_soft,I = softScale·4·m_eff,I/dt² ⇒ each contact mode's
// ω_I·dt = 2√softScale ≤ 2 (Courant-stable; the B1 rule generalized from NTS B=[n|−Nᵢ n] to the
// segment B_I=[D,−M]/a_I). p_I = min(0, k_soft,I·ḡ_I) scatters self-equilibratingly along n exactly
// like the C2 normal block: f^s_K = −(D·p)_K n, f^m_L = +(Mᵀ·p)_L n. Validated: proto_b2_soft2_segment.py.
void
LadrunoContactFE::addSoft2Penalty(Integrator *theIntegrator)
{
    double D[4][4], M[4][4], g[4], n[3];
    if (!mortarActive(D, M, g, n))
        return;                                   // no overlap this eval ⇒ zero force (KKT separation)

    double dt = 0.0;
    CentralDifferenceLadruno *cdl = dynamic_cast<CentralDifferenceLadruno *>(theIntegrator);
    if (cdl != 0) dt = cdl->getCurrentDeltaT();   // caller guarantees the cast, but be defensive

    // assembled translational masses of the facet nodes (kept FULL so the friction block below can
    // re-project them on the surface tangents), and their projection on the gap normal n (B1 helpers)
    double msF[4][3] = {{0}}, mmF[4][3] = {{0}};
    double invMs[4] = {0.0, 0.0, 0.0, 0.0}, invMm[4] = {0.0, 0.0, 0.0, 0.0};
    for (int I = 0; I < npsS; I++) {
        ladrunoNodeMass(theDomain, mortarSlave[I], msF[I]);
        invMs[I] = ladrunoInvMassProj(msF[I], n);
    }
    for (int K = 0; K < npsM; K++) {
        ladrunoNodeMass(theDomain, mortarMaster[K], mmF[K]);
        invMm[K] = ladrunoInvMassProj(mmF[K], n);  // fixed/massless master node ⇒ 0 (∞ mass)
    }

    double p[4] = {0.0, 0.0, 0.0, 0.0};
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;                       // a_I = Σ_J D_IJ = ∫ N_I dΓ (this facet)
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;           // unreferenced slave node
        double gbar = g[I] / aFacet;              // area-normalised weighted gap (<0 ⇒ penetration)
        if (gbar >= 0.0) continue;                // separation ⇒ no penalty (KKT clamp)
        // segment-based gap-mode inverse mass for node I (the oracle-validated closed form)
        double invMeff = 0.0;
        for (int J = 0; J < npsS; J++) { double c = D[I][J] / aFacet; invMeff += c * c * invMs[J]; }
        for (int K = 0; K < npsM; K++) { double c = M[I][K] / aFacet; invMeff += c * c * invMm[K]; }
        double knEff;
        if (dt > 0.0 && invMeff > 0.0) {
            knEff = softScale * 4.0 * (1.0 / invMeff) / (dt * dt);   // k_soft = SOFSCL·4·m_eff/dt²
        } else {
            knEff = kn;                           // cannot soft-size ⇒ the configured epsN (+ warn once)
            // contact-review P5: per-(contact, topic) engine latch, not a process-lifetime static
            LadrunoContactDomain *cdW = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
            if (cdW == 0 || cdW->warnOnce(contactTag, LadrunoContactDomain::WARN_SOFT2_NOMASS))
                opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": -mortar -soft "
                          "(SOFT=2) could not size a segment-based penalty under the explicit "
                          "integrator (non-positive dt or a contact node with no assembled mass); "
                          "using the configured epsN. Give the contact nodes mass for the SOFT=2 "
                          "stable penalty.\n";
        }
        p[I] = knEff * gbar;                       // ≤ 0 (compression)
    }

    // scatter the penalty pressure: f^s_K = −(D·p)_K n, f^m_L = +(Mᵀ·p)_L n (self-equilibrating,
    // Σφ=1). resid was zeroed at the top of getResidual and SOFT=2 owns the whole MORTAR residual.
    for (int K = 0; K < npsS; K++) {
        double Dp = 0.0;
        for (int I = 0; I < npsS; I++) Dp += D[K][I] * p[I];
        for (int d = 0; d < 3; d++) resid(3 * K + d) = -Dp * n[d];
    }
    for (int L = 0; L < npsM; L++) {
        double Mp = 0.0;
        for (int I = 0; I < npsS; I++) Mp += M[I][L] * p[I];
        for (int d = 0; d < 3; d++) resid(3 * (npsS + L) + d) = Mp * n[d];
    }

    // --- D2.2 viscous normal stabilization on the SOFT=2 explicit path ---
    // The SAME velocity-proportional normal damper the regular mortar path applies (getResidual,
    // D2.2 block), here on the SOFT=2 active set: per IN-CONTACT slave node (p[I]<0, the soft-penalty
    // KKT mask) the weighted normal gap RATE ḡ̇_I = n·(Σ_J D_IJ v_s,J − Σ_K M_IK v_m,K)/a_I (v =
    // getTrialVel), the viscous pressure p_visc_I = μ_c·ḡ̇_I (NO clamp — a dashpot active while in
    // contact), scattered EXACTLY like the normal penalty force (f^s += −(D·p_visc)n, f^m += +(Mᵀ·p_visc)n).
    // It is the C2 normal operator with epsN→μ_c, driven by velocity — the oracle-validated D2.2 operator
    // (proto_d2_viscous.py T6), unchanged. μ_c=0 ⇒ no viscous term (byte-identical to the plain SOFT=2
    // penalty). Under explicit CDL this is force-only (the tangent is mass-only); on the implicit
    // fall-through the regular mortar path's addCtoTang owns the consistent C tangent.
    if (muc > 0.0) {
        double vs[4][3], vm[4][3];
        for (int i = 0; i < npsS; i++) {
            const Vector &v = mortarSlave[i]->getTrialVel();
            for (int d = 0; d < 3; d++) vs[i][d] = v(d);
        }
        for (int i = 0; i < npsM; i++) {
            const Vector &v = mortarMaster[i]->getTrialVel();
            for (int d = 0; d < 3; d++) vm[i][d] = v(d);
        }
        double pv[4] = {0, 0, 0, 0};
        for (int I = 0; I < npsS; I++) {
            if (p[I] >= 0.0) continue;            // viscous only on in-contact nodes (soft-penalty mask)
            double aFacet = 0.0;
            for (int J = 0; J < npsS; J++) aFacet += D[I][J];
            if (aFacet <= 1e-300) continue;
            double rdot[3] = {0, 0, 0};
            for (int J = 0; J < npsS; J++)
                for (int d = 0; d < 3; d++) rdot[d] += D[I][J] * vs[J][d];
            for (int K = 0; K < npsM; K++)
                for (int d = 0; d < 3; d++) rdot[d] -= M[I][K] * vm[K][d];
            double gdot = (rdot[0]*n[0] + rdot[1]*n[1] + rdot[2]*n[2]) / aFacet;
            pv[I] = muc * gdot;                   // p_visc_I = μ_c·ḡ̇_I
        }
        for (int K = 0; K < npsS; K++) {          // slave: −(D·p_visc)_K n
            double Dp = 0.0;
            for (int I = 0; I < npsS; I++) Dp += D[K][I] * pv[I];
            for (int d = 0; d < 3; d++) resid(3 * K + d) += -Dp * n[d];
        }
        for (int L = 0; L < npsM; L++) {          // master: +(Mᵀ·p_visc)_L n
            double Mp = 0.0;
            for (int I = 0; I < npsS; I++) Mp += M[I][L] * pv[I];
            for (int d = 0; d < 3; d++) resid(3 * (npsS + L) + d) += Mp * n[d];
        }
    }

    // --- frictional SOFT=2: Courant-stable Coulomb/Tresca over the segment overlap ---
    // The SOFT=2 friction analogue of the normal penalty: per IN-CONTACT slave node (soft-penalty mask
    // p_I<0) run the SHIPPED LadrunoFrictionKernel return map with the SOFT normal pressure N_I = |p_I|
    // and a Courant-stable STICK penalty k_soft_t,I sized from the SEGMENT-BASED tangential gap-mode
    // mass m_eff_t,I (the B2 normal m_eff with n→t — the gap operator row B_t,I : slave J → (D_IJ/a_I)t,
    // master K → −(M_IK/a_I)t), worst-cased over the two basis tangents t1,t2 to n so a stiff friction
    // kt never throttles dt_cr via the stick mode ω_t=√(kt/m_eff_t) (the B1-kt rule, generalized to the
    // segment operator; oracle proto_soft2_friction.py). Penalty friction (λ_T≡0, single-pass — the
    // explicit SOFT=2 philosophy): the slip is the DISPLACEMENT-based weighted tangential gap (the C3.1
    // closest-point lesson — positions are purely normal), gT0-referenced. tFric scatters via D/−M like
    // the C3.1 path. μ≤0 ∧ c≤0 ∧ τmax≤0 SHORT-CIRCUITS ⇒ byte-identical to the frictionless SOFT=2.
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    if ((mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0) && cd != 0) {
        // two orthonormal tangents to n (the softKt construction)
        int kx = 0; double amin = std::fabs(n[0]);
        if (std::fabs(n[1]) < amin) { amin = std::fabs(n[1]); kx = 1; }
        if (std::fabs(n[2]) < amin) { kx = 2; }
        double e[3] = {0.0, 0.0, 0.0}; e[kx] = 1.0;
        double en = e[0]*n[0] + e[1]*n[1] + e[2]*n[2];
        double t1[3], nrm = 0.0;
        for (int d = 0; d < 3; d++) { t1[d] = e[d] - en * n[d]; nrm += t1[d]*t1[d]; }
        nrm = std::sqrt(nrm);
        bool haveT = (nrm > 1e-300);
        double t2[3] = {0.0, 0.0, 0.0};
        if (haveT) {
            for (int d = 0; d < 3; d++) t1[d] /= nrm;
            t2[0] = n[1]*t1[2] - n[2]*t1[1]; t2[1] = n[2]*t1[0] - n[0]*t1[2]; t2[2] = n[0]*t1[1] - n[1]*t1[0];
        }
        // facet node DISPLACEMENTS (the slip is displacement-based — the C3.1 lesson)
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
            if (p[I] >= 0.0) continue;            // friction only on in-contact nodes (soft mask)
            double aFacet = 0.0;
            for (int J = 0; J < npsS; J++) aFacet += D[I][J];
            if (aFacet <= 1e-300) continue;
            double N_I = -p[I];                   // soft normal pressure magnitude (≥ 0)
            // per-node Courant-stable stick penalty: worst-case m_eff_t over the two tangents.
            double ktEff = kt;                    // fallback: the configured tangential penalty
            if (haveT && dt > 0.0) {
                double inv1 = 0.0, inv2 = 0.0;
                for (int J = 0; J < npsS; J++) {
                    double c = D[I][J] / aFacet;
                    inv1 += c*c*ladrunoInvMassProj(msF[J], t1); inv2 += c*c*ladrunoInvMassProj(msF[J], t2);
                }
                for (int K = 0; K < npsM; K++) {
                    double c = M[I][K] / aFacet;
                    inv1 += c*c*ladrunoInvMassProj(mmF[K], t1); inv2 += c*c*ladrunoInvMassProj(mmF[K], t2);
                }
                double invMax = (inv1 > inv2) ? inv1 : inv2;   // MAX inv ⇒ MIN m_eff_t ⇒ binding bound
                if (invMax > 0.0) ktEff = softScale * 4.0 * (1.0 / invMax) / (dt * dt);
            }
            // weighted relative DISPLACEMENT r_I, tangential part / a_I (gT0-referenced)
            double r[3] = {0, 0, 0};
            for (int J = 0; J < npsS; J++)
                for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
            for (int K = 0; K < npsM; K++)
                for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
            double rn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
            double gbarT[3];
            for (int d = 0; d < 3; d++) gbarT[d] = (r[d] - rn * n[d]) / aFacet;
            LadrunoContactDomain::MortarNormalState &st =
                cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag());
            if (!st.engaged) {                    // engagement origin captured ONCE (the P3 MAJOR-1)
                for (int d = 0; d < 3; d++) st.gT0[d] = gbarT[d];
                st.engaged = true;
            }
            double gTeff[3];
            for (int d = 0; d < 3; d++) gTeff[d] = gbarT[d] - st.gT0[d];   // penalty friction (λ_T≡0)
            double tF[3], gpTtrial[3];
            LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, N_I, ktEff, mu, tF, gpTtrial,
                                                     mortarCohesion, mortarTauMax);
            for (int d = 0; d < 3; d++) {
                st.gpTtrial[d] = gpTtrial[d]; st.lambdaTtrial[d] = 0.0; tFric[I][d] = tF[d];
            }
        }
        // scatter via D (slave) / −M (master), like the C3.1 friction force (Σφ=1 self-equilibrating)
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
    } else if (mode == SEGMENT && ndm == 2) {
        // ADR-85 T1b -- K_c = kn*B^TB on the SAME active set as the residual
        // (segment2DActive re-evaluates deterministically; the committed vertex
        // side sign is NORMALLY captured by getResidual first, but correctness
        // does not hinge on that order: with no committed sign yet,
        // segment2DActive falls back to the SAME live wedge sign the capture
        // would commit, so an out-of-order tangent classifies identically --
        // the order-safe fallback is what carries it). Uses the configured kn,
        // NOT softKn:
        // SOFT is residual-only and addKtToTang is assembled only under
        // implicit integrators, where softKn == kn (the shipped B1 contract).
        // -geomtan (consistentNormal) is ACCEPTED AS A NO-OP here: the T1a FD
        // gate (proto_t1_nts2d.py family F8) proves the 2D B-operators are the
        // EXACT first variation of the gap, and per the ADR-85 T1b directive no
        // separate geometric term is fabricated for the 2D lane -- kn*B^TB is the
        // assembled tangent with or without -geomtan. (NOTE, disclosed: exact
        // FIRST variation does not make the SECOND variation zero -- a rotating
        // master segment carries a curvature block, omitted here exactly like
        // the shipped default 3D path omits it when -geomtan is off. The no-op
        // costs Newton RATE on rotating masters, never force correctness.)
        double gap, n[3], N2[2], B[6], gTnow = 0.0;
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0, &gTnow)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // ADR-85 T2 -- friction tangent (IMPLICIT only -- CDL never calls
            // addKtToTang). Reads COMMITTED gpT (not gpTtrial), exactly the 3D
            // SEGMENT rule below, so the tangent is the derivative of the
            // residual evaluated at the SAME state. Default consistentTan=false
            // (ct.consistentTan, threaded through the 2D ctor) ⇒ symmetric
            // (drop the dTN cross term, design-gate Q2); true ⇒ the full
            // non-symmetric consistent tangent (needs FullGeneral/UmfPack -- the
            // SAME parse-time WARNING the 3D -consistanttan path already emits,
            // dimension-independent).
            if (mu > 0.0 && theDomain != 0) {
                LadrunoContactDomain *cdT = theDomain->getLadrunoContactDomain();
                if (cdT != 0) {
                    LadrunoContactDomain::FrictionState &st =
                        cdT->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                    double gTeff = st.engaged ? (gTnow - st.gT0[0]) : 0.0;
                    double tn = LadrunoContactKernel::traction(kn, gap);
                    addFrictionTang2D(fact, n, N2, tn, gTeff, st.gpT[0], consistentTan);
                }
            }
        }
    } else if (mode == SEGMENT) {
        // K_c = kn BᵀB (main NTS term — EXACT for a flat/fixed master where n is constant).
        // B3 (P2b-2c): when consistentNormal, ALSO add the consistent ∂n/∂u geometric block
        // kn·gN·∂²gN/∂u² (curvature / large-sliding) so implicit Newton is QUADRATIC on a
        // CURVED interface. SYMMETRIC; for a flat facet the geometric block is 0 ⇒ byte-identical.
        double gap, n[3], N[4], B[15], gTvec[3];
        if (segmentActive(gap, n, N, B, gTvec)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // ADR-63 #4a: the faceted B3 ∂n/∂u block is derived for the per-segment normal, NOT for
            // the smoothed nodal-normal field — applying it under smoothing would be inconsistent with
            // the residual's n_smooth. Suppress it; the smoothed path ships the symmetric frozen-field
            // kn·BᵀB only (ADR-63 D4). The smoothed-normal consistent tangent (∂n_smooth/∂u) is the
            // gated P3. The handler emits a one-time silent-downgrade warning when both flags are set.
            if (consistentNormal && !useSmoothNormal)
                addNormalGeomTang(fact);
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
    } else if (mode == EDGE_EDGE) {
        // ADR-57 E2 — the main (penalty-Gram) tangent K = εN·BᵀB (12×12, symmetric, rank-1, PSD —
        // solver-safe on any system, like the NTS/mortar main blocks). EXACT at the closest point
        // for the frozen-direction linearization; the geometric ∂{n,s,t}/∂u curvature block is E4
        // (gated off). Same active mask as the residual (penetrating + margin-interior crossing).
        // The committed sign was captured by getResidual (which runs first in the iterate); read it.
        // E5 NOTE: this tangent uses `kn` (and friction `kt`) — NOT the SOFT `k_soft`. That is the
        // shipped B1/B2 SOFT contract: SOFT is RESIDUAL-ONLY, and addKtToTang is assembled ONLY under
        // implicit integrators (the explicit CentralDifferenceLadruno is mass-only — it never calls
        // addKtToTang). Under implicit, softKnEdge/softKtEdge return kn/kt (the dynamic_cast to CDL
        // fails), so the residual's knEff/ktEff ≡ kn/kt here ⇒ the tangent IS the residual derivative.
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        LadrunoContactDomain::EdgeEdgeState *st = 0;
        int committedSign = 0;
        if (cd != 0) {
            st = &cd->getOrCreateEdgeEdgeState(contactTag, edgeNode[0]->getTag(),
                                               edgeNode[1]->getTag(), edgeNode[2]->getTag(),
                                               edgeNode[3]->getTag());
            committedSign = st->signN;
        }
        double gN, n[3], s, t, B[4][3]; int usedSign = 0;
        // E6 — the active mask + cone N use the AUGMENTED pressure pAug = λ_N + εN·gN (the C2.2 rule:
        // the tangent is the derivative of the AUGMENTED residual at the frozen active set; the λ_N
        // offset only shifts the active SET, ∂λ/∂u=0 within a sweep ⇒ the penalty-Gram block is still
        // kn·BᵀB). -edgeAlm OFF ⇒ λ_N≡0 ⇒ pAug = kn·gN ⇒ the mask is gN<0 and tN=−kn·gN (byte-identical).
        // Implicit-only here (CDL never calls addKtToTang), so kn ≡ knEff (the E5 SOFT contract).
        double lambdaN = (edgeAlm && st != 0) ? st->lambdaN : 0.0;
        if (edgeGeom(gN, n, s, t, B, committedSign, &usedSign) && lambdaN + kn * gN < 0.0) {
            for (int i = 0; i < 4; i++)
                for (int di = 0; di < 3; di++)
                    for (int j = 0; j < 4; j++)
                        for (int dj = 0; dj < 3; dj++)
                            tang(3 * i + di, 3 * j + dj) += fact * kn * B[i][di] * B[j][dj];

            // --- E3 friction tangent (IMPLICIT only — CDL never calls addKtToTang) ---
            // K_fric[A][B] += w_A w_B·K_ss, w = [(1−s),s,−(1−t),−t], K_ss = frictionTangentBlock at the
            // SAME committed slip the residual used (gpT committed, gTeff at the current config). Default
            // consistentTan=false ⇒ symmetric (drop the Coulomb Csl; solver-safe); true ⇒ the full
            // non-symmetric consistent tangent (needs FullGeneral/UmfPack).
            if (st != 0 && (mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0)) {
                double tN = -(lambdaN + kn * gN);        // the AUGMENTED normal pressure (the cone's N)
                double gTvec[3]; edgeSlip(s, t, n, gTvec);
                double gTeff[3];
                for (int d = 0; d < 3; d++) gTeff[d] = st->engaged ? (gTvec[d] - st->gT0[d]) : 0.0;
                double Kss[3][3];
                LadrunoFrictionKernel::frictionTangentBlock(gTeff, st->gpT, n, tN, kn, kt, mu,
                                                            consistentTan, Kss, mortarCohesion,
                                                            mortarTauMax);
                double w[4] = { 1.0 - s, s, -(1.0 - t), -t };
                for (int A = 0; A < 4; A++)
                    for (int Bn = 0; Bn < 4; Bn++) {
                        double wab = fact * w[A] * w[Bn];
                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++)
                                tang(3 * A + i, 3 * Bn + j) += wab * Kss[i][j];
                    }
            }
        }
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
    if (isTie) {
        // C4 — MESH-TIE tangent: K = epsTie · B̃ᵀ diag(1/a_I) B̃ ⊗ I₃, B̃ = [D, −M]. The active set is
        // frozen ON (W_I = 1/a_I for EVERY slave node — no compression mask) and the per-node block is
        // the FULL identity I₃ (all 3 components tied — normal AND tangential), not the contact n⊗n.
        // SPD, symmetric, no Csl, no active-set switching ⇒ no -consistanttan needed (oracle T2).
        // initialStiff is irrelevant (the penalty Gram is constant; geometric ∂{D,M}/∂u deferred, as C2).
        double W[4] = {0.0, 0.0, 0.0, 0.0};
        for (int I = 0; I < npsS; I++) {
            double aFacet = 0.0;
            for (int J = 0; J < npsS; J++) aFacet += D[I][J];
            if (aFacet > 1e-300) W[I] = 1.0 / aFacet;
        }
        for (int A = 0; A < nN; A++) {
            for (int B = 0; B < nN; B++) {
                double Ks = 0.0;
                for (int I = 0; I < npsS; I++) {
                    double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                    double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                    Ks += bIA * W[I] * bIB;
                }
                Ks *= kn;                              // epsTie
                if (Ks == 0.0) continue;
                for (int d = 0; d < 3; d++)
                    tang(3 * A + d, 3 * B + d) += fact * Ks;   // ⊗ I₃ (diagonal in the 3-space)
            }
        }
        return;
    }
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
            double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
            double gTeff[3];                          // gTeff_eff = gTeff + λ_T/epsT (same as the residual)
            for (int d = 0; d < 3; d++)
                gTeff[d] = (r[d] - rn * n[d]) / aFacet - st.gT0[d] + st.lambdaT[d] * invEpsT;
            // initial-stiffness path ⇒ force the SPD STICK tangent kt·P_t: pass gTeff == gpT so the
            // trial traction is zero (‖tT*‖ ≤ cap ⇒ the kernel returns the stick block). Avoids the
            // rank-deficient slip tangent stalling Modified/Initial-Newton (gate MINOR-1, mirrors SEGMENT).
            const double *gtForKss = initialStiff ? st.gpT : gTeff;
            double Kss[3][3];
            // C3.3: consistentTan ⇒ the full NON-SYMMETRIC Coulomb tangent (the Csl pressure coupling
            // −μ·epsN·t̂⊗n scatters through the normal-gap operator via the same b·b/a — kn=epsN).
            // Needs FullGeneral/UmfPack; default false ⇒ the symmetric tangent (solver-safe). Forced
            // symmetric on the initial-stiffness path (stick, no slip ⇒ no Csl).
            bool useConsistent = consistentTan && !initialStiff;
            LadrunoFrictionKernel::frictionTangentBlock(gtForKss, st.gpT, n, N_I, kn, kt, mu,
                                                        useConsistent, Kss,
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
    } else if (mode == SEGMENT && ndm == 2) {
        // ADR-85 T1b -- initial-stiffness path: K_initial == K_current == kn*B^TB
        // (mirror addKtToTang so Newton -initial / ModifiedNewton -initial do not
        // silently drop the contact stiffness -- the shipped SEGMENT rationale,
        // ndm-generalized).
        double gap, n[3], N2[2], B[6];
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // ADR-85 T2 -- the friction INITIAL stiffness is the SPD STICK
            // tangent kt*(t_hat⊗t_hat) (mirrors the 3D SEGMENT rule below, gate
            // Q5): gTeff==gpT==0.0 forces tangentBlock1D's stick branch
            // regardless of the actual committed state (|kt*(0-0)| = 0 <= any
            // positive cap), so no engine slot is needed; consistent=false is
            // forced too (stick carries no Csl cross term).
            if (mu > 0.0) {
                double tn = LadrunoContactKernel::traction(kn, gap);
                addFrictionTang2D(fact, n, N2, tn, 0.0, 0.0, false);
            }
        }
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
    } else if (mode == EDGE_EDGE) {
        // ADR-57 E2 — K_initial == K_current for the edge-edge penalty (the frozen-direction Gram
        // εN·BᵀB), so mirror addKtToTang; the base addKiToTang early-returns on myEle==0 and would
        // otherwise silently drop the contact stiffness under Newton -initial / ModifiedNewton.
        LadrunoContactDomain::EdgeEdgeState *st = 0;
        int committedSign = 0;
        if (theDomain != 0) {
            LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
            if (cd != 0) {
                st = &cd->getOrCreateEdgeEdgeState(contactTag, edgeNode[0]->getTag(),
                                    edgeNode[1]->getTag(), edgeNode[2]->getTag(),
                                    edgeNode[3]->getTag());
                committedSign = st->signN;
            }
        }
        // E6 — the augmented active mask + cone N (the C2.2 rule, mirroring addKtToTang). -edgeAlm OFF
        // ⇒ λ_N≡0 ⇒ mask is gN<0, tN=−kn·gN (byte-identical).
        double lambdaN = (edgeAlm && st != 0) ? st->lambdaN : 0.0;
        double gN, n[3], s, t, B[4][3]; int usedSign = 0;
        if (edgeGeom(gN, n, s, t, B, committedSign, &usedSign) && lambdaN + kn * gN < 0.0) {
            for (int i = 0; i < 4; i++)
                for (int di = 0; di < 3; di++)
                    for (int j = 0; j < 4; j++)
                        for (int dj = 0; dj < 3; dj++)
                            tang(3 * i + di, 3 * j + dj) += fact * kn * B[i][di] * B[j][dj];
            // E3: the friction INITIAL stiffness is the STICK tangent kt·P_t (SPD ⇒ a Modified/
            // Initial-Newton contraction; the SEGMENT addKiToTang rule). Forced stick ⇒ gTeff==gpT
            // (pass zeros) ⇒ K_ss = kt·P_t, symmetric ⇒ consistent=false regardless; no engine slot.
            if (mu > 0.0 || mortarCohesion > 0.0 || mortarTauMax > 0.0) {
                double tN = -(lambdaN + kn * gN);
                double zero[3] = {0.0, 0.0, 0.0}, Kss[3][3];
                LadrunoFrictionKernel::frictionTangentBlock(zero, zero, n, tN, kn, kt, mu,
                                                            false, Kss, mortarCohesion, mortarTauMax);
                double w[4] = { 1.0 - s, s, -(1.0 - t), -t };
                for (int A = 0; A < 4; A++)
                    for (int Bn = 0; Bn < 4; Bn++) {
                        double wab = fact * w[A] * w[Bn];
                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++)
                                tang(3 * A + i, 3 * Bn + j) += wab * Kss[i][j];
                    }
            }
        }
    }
}

void
LadrunoContactFE::addCtoTang(double fact)
{
    // D2 viscous stabilization — the only contact damping. C_visc = muc·B Bᵀ (the normal-penalty
    // operator with kn→muc), assembled on the SAME active set as the residual viscous force. The
    // transient integrators call formEleTangent → addCtoTang(c2) (Newmark/HHT/CentralDifference),
    // so implicit Newton gets the consistent damping tangent; the fork's mass-only CDL never calls
    // it ⇒ explicit is force-only (the P3 friction explicit/implicit split). muc=0 ⇒ no-op (the
    // shipped no-damping behavior — byte-identical). RIGID_PLANE + SEGMENT (D2.1) + MORTAR (D2.2);
    // a tie (isTie) has no contact-chatter regime ⇒ no viscous (refused at the command surface).
    if (muc <= 0.0)
        return;
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * muc * planeN[i] * planeN[j];
    } else if (mode == SEGMENT && ndm == 2) {
        // ADR-85 T1b -- the T1b dashpot port's damping tangent: C = muc*B^TB on
        // the same active set as the residual dashpot (implicit transient
        // integrators call addCtoTang(c2); the mass-only explicit CDL never
        // does, so explicit stays force-only -- the shipped D2 contract with the
        // stride-3 loops become-ndm'd: ndof = 2*(1+nps)).
        double gap, n[3], N2[2], B[6];
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * muc * B[i] * B[j];
        }
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * muc * B[i] * B[j];
        }
    } else if (mode == MORTAR && !isTie) {
        // D2.2 — C_visc = μ_c · B̃ᵀ diag(W) B̃ ⊗ (n⊗n), B̃=[D,−M], W_I=1/a_I on the contact active
        // set (p_I<0). The C2 normal-penalty Gram block (addMortarTang) with epsN→μ_c, in C not K.
        double D[4][4], M[4][4], g[4], n[3];
        if (!mortarActive(D, M, g, n)) return;
        int nN = npsS + npsM;
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        double W[4] = {0.0, 0.0, 0.0, 0.0};
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
                double Ks = 0.0;
                for (int I = 0; I < npsS; I++) {
                    double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                    double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                    Ks += bIA * W[I] * bIB;
                }
                Ks *= muc;
                if (Ks == 0.0) continue;
                for (int dA = 0; dA < 3; dA++)
                    for (int dB = 0; dB < 3; dB++)
                        tang(3 * A + dA, 3 * B + dB) += fact * Ks * n[dA] * n[dB];
            }
        }
    }
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
