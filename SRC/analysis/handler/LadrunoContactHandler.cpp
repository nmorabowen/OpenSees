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

// ADR-39 P1a — LadrunoContactHandler implementation. See header for design.
// Replicates the PlainHandler DOF_Group/FE loop (mirrors the proven in-fork
// LadrunoProjectionHandler) and appends the contact FE adapter(s).

#include "LadrunoContactHandler.h"
#include "LadrunoContactFE.h"

#include <LadrunoContactDomain.h>    // Ladruno: ADR-39 (adapter count from the engine)
#include <LadrunoContactSurface.h>   // Ladruno: ADR-39 P2a (slave node-set)
#include <LadrunoContactKernel.h>    // Ladruno: ADR-39 P2b-2b (reference normal for -kn auto)
#include <LadrunoContactProjection.h> // Ladruno: ADR-41 A2 (normalOriented projection geometry)
#include <LadrunoContactBucketSort.h>// Ladruno: ADR-39 P2.5 (broad-phase pairing)
#include <LadrunoContact2DKernel.h>  // Ladruno: ADR-85 T1b (2D NTS orientation vote + auto-kn normal)
#include <LadrunoContactReemit.h>   // Ladruno: ADR-60 (finite-sliding re-emit band/metric)
#include <LadrunoContactNormalField.h> // Ladruno: ADR-63 #4a (nodal-normal smoothing status/build)
#include <LadrunoEdgeKernel.h>       // Ladruno: ADR-57 E2 (edge-edge routing + closest point)
#include <Matrix.h>                  // Ladruno: ADR-39 P2b-2b (master getInitialStiff)
#include <Domain.h>
#include <AnalysisModel.h>
#include <Integrator.h>
#include <Node.h>
#include <Vector.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <EQ_Constraint.h>           // Ladruno: contact-review P4 (upstream-parity EQ handling)
#include <EQ_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <FE_Element.h>
#include <Subdomain.h>
#include <ID.h>
#include <classTags.h>
#include <elementAPI.h>
#include <map>
#include <set>
#include <vector>
#include <cmath>

// ADR-78 P1 — the fatal exit lives in LadrunoContactAbort.cpp, NOT here. This TU
// is compiled into the OPS_Analysis OBJECT library, which is built once with no
// parallel define and folded into every target, so an `#ifdef _PARALLEL_*` block
// written HERE compiles to nothing everywhere -- including OpenSeesMP. That was
// measured: a guarded MPI_Abort placed in this file vanished silently and the
// mutation deck hung exactly as it had before. See LadrunoContactAbort.h.
#include "LadrunoContactAbort.h"

// ----------------------------------------------------------------------------
// contact-review P3: pre-flight a contact's surfaces — every referenced node must
// exist and carry 3D coordinates. Previously a typo'd node tag was skipped in
// SILENCE (dead pairs, localized pass-through with no diagnostic), and a 2D node
// (`-ndm 2 -ndf 3` frame) passed the ndf==3 guards while getCrds()(2) read out of
// bounds (unchecked in release ⇒ nondeterministic garbage geometry). Refuse the
// WHOLE contact loudly instead of silently processing a partial surface.
//
// ADR-78 P1 (2026-08-11): the refusal is now an ABORT, not a skip. Warning-and-
// continuing means a declared contact quietly does not happen, and the resulting
// run looks entirely plausible — the P0 mutation battery produced exactly that:
// a model that converged, balanced its reactions, and was wrong by a factor of
// two. Under MPI it is worse still, because a warning on rank 7 of 240 scrolls
// past unread. Callers turn a false return into `return -1` from handle(), which
// StaticAnalysis/DirectIntegrationAnalysis::domainChanged() report as a failure.
// ----------------------------------------------------------------------------
// ADR-85 T0 (2026-08-18): the coordinate rule is now dimension CONSISTENCY, not
// "must be 3D". A surface whose nodes are ALL 2D or ALL 3D passes this
// pre-flight and reports its common dimension through `dimOut`; the LANE gate
// (ladrunoContact2DLaneOk, below) then decides what a 2D surface is allowed to
// do on the lane that referenced it. What stays refused is precisely what the
// pre-3D guard incident was about: a node set that MIXES getCrds() sizes, where
// getCrds()(2) reads out of bounds on the 2D members (unchecked in release =>
// nondeterministic garbage geometry). An all-3D surface takes exactly the
// shipped path and emits exactly the shipped messages.
static bool
ladrunoSurfaceNodesOk(Domain *dom, int ctag, const LadrunoContactSurface *s,
                      int *dimOut = 0)
{
    const ID &tags = s->getNodeTags();
    int dim = 0;                       // ADR-85: the surface's common getCrds() size
    for (int i = 0; i < tags.Size(); i++) {
        Node *n = dom->getNode(tags(i));
        if (n == 0) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " does not exist; ABORTING (a declared contact is never "
                      "silently dropped -- ADR-78 P1)\n";
            return false;
        }
        int nd = n->getCrds().Size();
        if (nd < 2) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " has " << nd << "D coordinates; contact needs 2D or 3D nodes; "
                      "ABORTING (ADR-85 dimension; was ADR-78 P1)\n";
            return false;
        }
        if (dim == 0)
            dim = nd;
        else if (nd != dim) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": contactSurface " << s->getTag() << " MIXES node dimensions -- node "
                   << tags(i) << " has " << nd << "D coordinates but an earlier node of the "
                      "same surface has " << dim << "D. A contact surface must be "
                      "dimension-consistent (this is the historical out-of-bounds "
                      "reproducer: getCrds()(2) on a 2D node is unchecked in release); "
                      "ABORTING (ADR-85 dimension)\n";
            return false;
        }
    }
    if (dimOut != 0)
        *dimOut = dim;               // 0 for an empty surface => the shipped 3D path
    return true;
}

// ----------------------------------------------------------------------------
// ADR-85 T0 -- the 2D LANE gate.
//
// A DECLARATION may now be 2D (`contactSurface -master 2` over 2-coordinate
// nodes), but the 2D NTS kernel + FE SEGMENT path lands in T1b and the 2D mortar
// interval integrator in T3. A 2D surface on one of those lanes is therefore
// refused BY NAME. Silently skipping it is the failure mode this subsystem's
// abort discipline exists to prevent: the run converges, balances its reactions
// and transmits no contact force (ADR-78 P0 -- wrong by a factor of two, and
// entirely plausible-looking).
//
// `dim` is what ladrunoSurfaceNodesOk() wrote. dim != 2 returns true before
// touching anything else, so no 3D deck can reach a single statement below the
// first line of this function.
//
// The ndf check comes FIRST and is EXACT equality (ADR-85 Why + How/8): a
// `-ndm 2 -ndf 3` frame node is the historical reproducer, and FE_Element::setID()
// packs each node's FULL DOF_Group ID sequentially -- an extra rotation DOF
// pushes equations into the contact's myID slots and shifts every later slot
// (silent mis-assembly). Reporting "not yet supported" for such a deck would
// hide a modeling error behind a roadmap message.
// ----------------------------------------------------------------------------
// ADR-85 T1b: the ndf == ndm == 2 EQUALITY gate, split out of
// ladrunoContact2DLaneOk so the now-LIVE 2D NTS lane keeps the equality
// contract (same message, same "ADR-85 ndf" literal) without the
// not-yet-supported refusal that T1b replaces. The T3 lanes still get both,
// through ladrunoContact2DLaneOk below.
static bool
ladrunoContact2DNdfOk(Domain *dom, int ctag, const LadrunoContactSurface *s,
                      int dim, const char *lane)
{
    if (dim != 2)
        return true;                 // 3D (or an empty surface): shipped path, untouched
    const ID &tags = s->getNodeTags();
    for (int i = 0; i < tags.Size(); i++) {
        Node *n = dom->getNode(tags(i));
        if (n != 0 && n->getNumberDOF() != 2) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " has ndf=" << n->getNumberDOF() << " on a 2D contact surface. The "
                   << lane << " lane requires ndf == ndm == 2 EXACTLY (FE_Element::setID() "
                      "packs each node's FULL DOF_Group ID, so an extra DOF pushes equations "
                      "into the contact's myID slots and shifts every later slot -- silent "
                      "mis-assembly). ABORTING (ADR-85 ndf equality contract)\n";
            return false;
        }
    }
    return true;
}

static bool
ladrunoContact2DLaneOk(Domain *dom, int ctag, const LadrunoContactSurface *s,
                       int dim, const char *lane, const char *phase)
{
    if (dim != 2)
        return true;                 // 3D (or an empty surface): shipped path, untouched
    if (!ladrunoContact2DNdfOk(dom, ctag, s, dim, lane))
        return false;
    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
           << ": contactSurface " << s->getTag() << " is 2D and the " << lane
           << " lane is NOT YET SUPPORTED in 2D -- it lands in ADR-85 " << phase
           << ". ABORTING rather than skipping the pair: a silently dropped contact "
              "converges, balances its reactions and transmits nothing (ADR-85 T0; the "
              "ADR-78 P0 failure mode)\n";
    return false;
}

// ----------------------------------------------------------------------------
// ADR-85 T0 -- the PAIR-level dimension gate, run at each shipped
// ladrunoSurfaceNodesOk() call-site pair. This is where an ARITY-FREE surface is
// first seen, and it closes a hole the per-surface gate above cannot:
//
// ladrunoContact2DLaneOk keys off ONE surface, and the `nps == 2` early blocks
// only fire when the MASTER declares a 2-node segment. Neither covers the NTS
// lane's SLAVE, which is a SLAVE_NODES set with NO nodesPerSeg at all. So a 3D
// master (nps 3, 3D nodes) could be paired with an all-2D slave node-set: the
// nps == 2 block would not fire, and the rewritten ladrunoSurfaceNodesOk would
// PASS the consistent all-2D slave where the shipped "< 3 coords" check refused
// it. That 2D slave then flows into the 3D NTS lane, where a `-ndm 2 -ndf 3`
// node clears the `ndf != 3` guard and `getCrds()(2)` reads OUT OF BOUNDS --
// exactly the historical incident this subsystem's guard exists for.
//
// The pair is therefore gated on BOTH dimensions together:
//   * MIXED across the two surfaces -> named FATAL. New territory: the shipped
//     pre-flight refused every 2D surface, so a cross-dimension pair was
//     previously impossible and no 3D deck can reach this branch.
//   * both 2D -> the per-surface lane gate (ndf equality first, then "not yet
//     supported"), same as the nps == 2 blocks.
// An all-3D pair returns on the first line, so the shipped path is untouched and
// the dimension capture that feeds this is write-only (dimOut) -- no side effects.
// ----------------------------------------------------------------------------
// ADR-85 T1b: `lane2DLive` opts a lane's call site into the LIVE 2D path -- a
// well-formed 2D pair then passes (after the ndf-equality gate on both
// surfaces) instead of drawing the T0 not-yet-supported refusal, and
// `pairDimOut` (write-only) reports the pair dimension so the caller can
// branch. Every OTHER gate (existence, per-surface consistency, the
// cross-dimension refusal) is shared verbatim between the live and the T0
// modes. Defaults keep the T3 call sites byte-identical.
static bool
ladrunoContactPairDimOk(Domain *dom, int ctag,
                        const LadrunoContactSurface *ms, int mdim,
                        const LadrunoContactSurface *ss, int sdim,
                        const char *lane, const char *phase,
                        bool lane2DLive = false, int *pairDimOut = 0)
{
    if (pairDimOut != 0)             // ADR-85 T1b (write-only; no 3D effect)
        *pairDimOut = (mdim == 2 || sdim == 2) ? 2 : 3;
    if (mdim == 3 && sdim == 3)
        return true;                 // the shipped 3D path, evaluated no further
    if (mdim != 0 && sdim != 0 && mdim != sdim) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
               << ": this pair MIXES a " << mdim << "D master surface (" << ms->getTag()
               << ") with a " << sdim << "D slave surface (" << ss->getTag()
               << "). Both surfaces of one contact must carry the SAME node dimension: a "
                  "2D node set fed to a 3D lane reads getCrds()(2) out of bounds "
                  "(unchecked in release => nondeterministic garbage geometry), and a "
                  "`-ndm 2 -ndf 3` node clears the ndf == 3 guard on the way in. "
                  "ABORTING (ADR-85 dimension)\n";
        return false;
    }
    // ADR-85 T1b: a LIVE lane keeps the ndf == ndm == 2 equality contract on
    // both (2D) surfaces and then PROCEEDS -- the not-yet-supported refusal is
    // exactly what T1b replaces on the NTS lane. Only reachable with a 2D
    // surface in the pair, so no 3D deck can enter.
    if (lane2DLive)
        return ladrunoContact2DNdfOk(dom, ctag, ms, mdim, lane) &&
               ladrunoContact2DNdfOk(dom, ctag, ss, sdim, lane);
    // dim 0 (an empty surface) and dim 3 pass ladrunoContact2DLaneOk unchanged.
    return ladrunoContact2DLaneOk(dom, ctag, ms, mdim, lane, phase) &&
           ladrunoContact2DLaneOk(dom, ctag, ss, sdim, lane, phase);
}

// ADR-85 T0 -- the whole pre-flight scaffold in one place, so the four call sites
// cannot drift: resolve BOTH surfaces (existence + per-surface dimension
// consistency), then gate the PAIR (cross-dimension refusal, then the per-surface
// lane gates). Returns false with the message ALREADY emitted, so every caller
// reads `return ladrunoContactFatal();`. An all-3D pair costs two dimOut writes
// and one comparison -- no output, no state change.
static bool
ladrunoContactPairPreflight(Domain *dom, int ctag,
                            const LadrunoContactSurface *ms,
                            const LadrunoContactSurface *ss,
                            const char *lane, const char *phase,
                            bool lane2DLive = false, int *pairDimOut = 0)
{
    int mdim = 0, sdim = 0;
    if (!ladrunoSurfaceNodesOk(dom, ctag, ms, &mdim) ||
        !ladrunoSurfaceNodesOk(dom, ctag, ss, &sdim))
        return false;
    return ladrunoContactPairDimOk(dom, ctag, ms, mdim, ss, sdim, lane, phase,
                                   lane2DLive, pairDimOut);
}

// ----------------------------------------------------------------------------
// command factory: constraints LadrunoContact
// ----------------------------------------------------------------------------
void *OPS_LadrunoContactHandler(void)
{
    return new LadrunoContactHandler();
}

// ----------------------------------------------------------------------------
LadrunoContactHandler::LadrunoContactHandler()
  : ConstraintHandler(HANDLER_TAG_LadrunoContactHandler)
{
}

LadrunoContactHandler::~LadrunoContactHandler()
{
}

// P2b-2b: auto-size the penalty kn for ONE (slave, master-segment) pair from the
// OWNING solid element's initial stiffness, sourced GENERICALLY through base-Element
// virtuals (getExternalNodes / getNumDOF / getInitialStiff) so any 3-DOF/node solid
// works with NO element-type coupling and NO vanilla edit. Reduction (validated in
// proto_p2b2b_autokn.py): kn = f_si * mean_over_seg_nodes( nᵀ K_block_node n ), where
// K_block_node is the 3x3 diagonal block at that node's DOFs and n is the reference
// segment normal. This is the element-stiffness form of the LS-DYNA 26.14a penalty
// f·K·A²/V (K_diag ~ E·L ~ E·A²/V), the A²/V geometry absorbed exactly into the
// assembled matrix. f_si = 0.10 (LS-DYNA SLSFAC default). Returns <= 0 on failure
// (no owning solid found / non-3-DOF element / ambiguous reference normal) so the
// caller skips the pair with a warning. The SOFT Courant floor (26.15) is P2b-2c.
static double
ladrunoResolveAutoKn(Domain *theDomain, Node **segNodes, int nps,
                     const double orientDir[3])
{
    const double F_SI = 0.10;
    // reference (undeformed) segment node coords
    double Xref[4][3];
    for (int k = 0; k < nps; k++) {
        const Vector &Xk = segNodes[k]->getCrds();
        for (int d = 0; d < 3; d++) Xref[k][d] = Xk(d);
    }
    // reference outward normal at the segment parametric center, oriented by orientDir
    // (winding-immune; refuses an ambiguous in-plane reference direction).
    double xiC  = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double etaC = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double n[3];
    if (!LadrunoContactProjection::normalOriented(nps, xiC, etaC, Xref, orientDir, n))
        return -1.0;
    int segTags[4];
    for (int k = 0; k < nps; k++) segTags[k] = segNodes[k]->getTag();
    // owning solid = the first Domain element that (a) carries 3 translational DOFs
    // per node and (b) contains ALL the segment's nodes. Brute force over the Domain
    // elements (fine for the gate meshes; the bucket-sort broad phase is P2.5).
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0 || e->getNumDOF() != 3 * nn) continue;   // not a 3-DOF/node solid
        int loc[4]; bool allIn = true;
        for (int k = 0; k < nps; k++) {
            loc[k] = en.getLocation(segTags[k]);
            if (loc[k] < 0) { allIn = false; break; }
        }
        if (!allIn) continue;
        const Matrix &K = e->getInitialStiff();
        if (K.noRows() != 3 * nn || K.noCols() != 3 * nn) return -1.0;
        double acc = 0.0;
        for (int k = 0; k < nps; k++) {
            int c = 3 * loc[k];
            double Kn[3];
            for (int i = 0; i < 3; i++)
                Kn[i] = K(c + i, c + 0) * n[0] + K(c + i, c + 1) * n[1] + K(c + i, c + 2) * n[2];
            acc += n[0] * Kn[0] + n[1] * Kn[1] + n[2] * Kn[2];
        }
        double kn = F_SI * acc / nps;
        return (kn > 0.0) ? kn : -1.0;
    }
    return -1.0;   // no owning solid element found for this segment
}

// ADR-85 T1b -- the 2D sibling of ladrunoResolveAutoKn (ADR-85 Where: auto-kn
// "needs more than the advertised 3*nn -> ndm*nn generalization"). Same LS-DYNA
// element-stiffness reduction kn = f_si * mean_over_seg_nodes(n^T K_block n),
// with (a) the 2-DOF/node owning-solid gate (getNumDOF() == 2*nn -- the fork's
// plane-element family), (b) 2x2 nodal stiffness blocks, and (c) the KERNEL's
// segment normal in place of the tri/quad-only
// LadrunoContactProjection::normalOriented: n = sigma*perp(t)/L at the REFERENCE
// segment (LadrunoContact2DKernel::perp2; sigma is the per-surface How/2 vote,
// already resolved by the caller, so the 3D resolver's ambiguous-orientation
// failure path is structurally absent here -- and n^T K n is sign-independent
// anyway). Returns <= 0 on failure (degenerate segment / no owning 2-DOF solid /
// bad K size); callers convert that to the same FATAL the 3D lane uses.
static double
ladrunoResolveAutoKn2D(Domain *theDomain, Node **segNodes, double sigma)
{
    const double F_SI = 0.10;
    if (segNodes[0] == 0 || segNodes[1] == 0) return -1.0;
    const Vector &XA = segNodes[0]->getCrds();
    const Vector &XB = segNodes[1]->getCrds();
    const double t[2] = { XB(0) - XA(0), XB(1) - XA(1) };
    const double L = std::sqrt(t[0]*t[0] + t[1]*t[1]);
    if (!(L > 0.0)) return -1.0;                       // zero-length reference segment
    double p[2];
    LadrunoContact2DKernel::perp2(t, p);
    const double n[2] = { sigma * p[0] / L, sigma * p[1] / L };
    int segTags[2] = { segNodes[0]->getTag(), segNodes[1]->getTag() };
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0 || e->getNumDOF() != 2 * nn) continue;   // not a 2-DOF/node solid
        int loc[2]; bool allIn = true;
        for (int k = 0; k < 2; k++) {
            loc[k] = en.getLocation(segTags[k]);
            if (loc[k] < 0) { allIn = false; break; }
        }
        if (!allIn) continue;
        const Matrix &K = e->getInitialStiff();
        if (K.noRows() != 2 * nn || K.noCols() != 2 * nn) return -1.0;
        double acc = 0.0;
        for (int k = 0; k < 2; k++) {
            int c = 2 * loc[k];
            double Kn0 = K(c + 0, c + 0) * n[0] + K(c + 0, c + 1) * n[1];
            double Kn1 = K(c + 1, c + 0) * n[0] + K(c + 1, c + 1) * n[1];
            acc += n[0] * Kn0 + n[1] * Kn1;
        }
        double kn = F_SI * acc / 2.0;
        return (kn > 0.0) ? kn : -1.0;
    }
    return -1.0;   // no owning 2-DOF/node solid element found for this segment
}

// ADR-85 T3 -- PURE CODE MOTION (the phase brief's ONE sanctioned shared-line refactor):
// the interface-level orientation vote (How/2), extracted VERBATIM from the T1b NTS
// branch's inline block so the T3 mortar branch can share it without retyping the math.
// EXACT SAME operation order and ALL message strings as the code this replaces -- the NTS
// lane is shipped (54-test suite) and G-T3's dump gate would catch drift, so this is
// literal code motion, not a rewrite.
//
// ONE sigma per interface, from the reference config: refDir = the explicit -outward (2
// components) or the slave-surface reference centroid minus the master-surface reference
// centroid over UNIQUE master nodes (the flat mTags double-counts shared vertices -- the
// smoothSeed lesson; `sTags`' own duplicate-vertex handling is inherited unchanged from
// the NTS call site -- see the T3 implementation report). Every master segment votes
// sigmaFromRef2D(t_ref, refDir); the vote must be UNANIMOUS among nonzero voters.
// Degenerate (no voter, or a centroid datum below the magnitude floor) or SPLIT => the
// message is emitted and this returns false; the caller does `return
// ladrunoContactFatal();`. sigma is FIXED at pairing (handle() re-runs re-derive it).
static bool
ladruno2DOrientationVote(Domain *theDomain, int ctag, bool hasOutward, const double outward[2],
                         const ID &mTags, int nSeg, const ID &sTags, double Lref,
                         double &sigmaOut)
{
    double refDir[2] = { 0.0, 0.0 };
    if (hasOutward) {
        refDir[0] = outward[0]; refDir[1] = outward[1];
    } else {
        double mc2[2] = { 0.0, 0.0 }; int mcn = 0;
        std::set<int> mseen;
        for (int i = 0; i < mTags.Size(); i++) {
            if (mseen.count(mTags(i))) continue;
            mseen.insert(mTags(i));
            const Vector &X = theDomain->getNode(mTags(i))->getCrds();
            mc2[0] += X(0); mc2[1] += X(1); mcn++;
        }
        double sc2[2] = { 0.0, 0.0 }; int scn = 0;
        // REVIEW FIX (T3): dedup the SLAVE list too. At the original NTS call site
        // sTags is a SLAVE_NODES set (duplicates structurally impossible, so this is
        // arithmetically a NO-OP there -- same adds, same order); the T3 mortar call
        // site passes a SLAVE_SEGMENTS flat pair list where every interior chain
        // vertex appears TWICE, and a duplicate-weighted centroid can cross the
        // magnitude floor or, on a curved interface, flip the vote (the smoothSeed
        // lesson, now on the slave side).
        std::set<int> sseen;
        for (int i = 0; i < sTags.Size(); i++) {
            if (sseen.count(sTags(i))) continue;
            sseen.insert(sTags(i));
            const Vector &X = theDomain->getNode(sTags(i))->getCrds();
            sc2[0] += X(0); sc2[1] += X(1); scn++;
        }
        if (mcn > 0 && scn > 0) {
            refDir[0] = sc2[0] / scn - mc2[0] / mcn;
            refDir[1] = sc2[1] / scn - mc2[1] / mcn;
        }
    }
    // ADR-85 T1b panel BLOCKER fix -- the MAGNITUDE gate on the CENTROID datum (see the
    // original NTS site for the full rationale: an unguarded angle-only gate resolves a
    // flush deck's noise-sized separation CONFIDENTLY to the wrong side).
    bool voteMagOk = true;
    if (!hasOutward) {
        const double magFloor = 1.0e-6 * Lref;
        voteMagOk = (refDir[0] * refDir[0] + refDir[1] * refDir[1] > magFloor * magFloor);
    }
    int votePlus = 0, voteMinus = 0;
    for (int seg = 0; seg < nSeg; seg++) {
        const Vector &A = theDomain->getNode(mTags(2 * seg))->getCrds();
        const Vector &Bc = theDomain->getNode(mTags(2 * seg + 1))->getCrds();
        const double t[2] = { Bc(0) - A(0), Bc(1) - A(1) };
        double sg = LadrunoContact2DKernel::sigmaFromRef2D(t, refDir);
        if (sg > 0.0) votePlus++;
        else if (sg < 0.0) voteMinus++;
    }
    if ((votePlus == 0 && voteMinus == 0) || !voteMagOk) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
               << ": the 2D interface orientation vote is DEGENERATE -- "
               << (hasOutward
                       ? "the given -outward direction is zero or (numerically) "
                         "tangent to every master segment"
                       : "the reference-config centroid vote (slave-surface "
                         "centroid minus master-surface centroid) is zero or "
                         "tangent to every master segment (a flush/coincident "
                         "interface)")
               << ". The outward side is genuinely ambiguous and is never "
                  "guessed: pass -outward ox oy with a component along the "
                  "expected contact normal. ABORTING (ADR-85 T1b orientation "
                  "vote)\n";
        return false;
    }
    if (votePlus > 0 && voteMinus > 0) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
               << ": the 2D interface orientation vote SPLIT (" << votePlus
               << " master segments vote +1, " << voteMinus << " vote -1) -- "
                  "the master surface is strongly curved or inconsistently "
                  "wound, so ONE per-surface sigma cannot orient it. Re-wind "
                  "the master segments consistently, split the surface into "
                  "separately declared contacts, or pass -outward ox oy. "
                  "ABORTING (ADR-85 T1b orientation vote)\n";
        return false;
    }
    sigmaOut = (votePlus > 0) ? 1.0 : -1.0;
    return true;
}

// B1 (P4) — build the ASSEMBLED translational nodal-mass cache the SOFT=1 penalty needs. The explicit
// integrator inverts the assembled global diagonal M = Σ_elements diag(M_e) + the nodal `mass`, but
// Node::getMass() holds ONLY the nodal `mass` (element-density mass never reaches the node). So here we
// reconstruct, per node, m[d] = nodal mass(d) + Σ_e diag(M_e) at that node's translational DOFs (the
// SAME diagonal `system Diagonal` inverts) and stash it on the engine for the stateless adapter to
// read. One pass over the elements; called only when a SOFT contact exists (masses are constant
// between domain changes, so rebuild-per-handle is correct). Translation-first DOF order (OpenSees
// nodes order [u | θ]) ⇒ a node's first 3 element DOFs are translational.
static void
ladrunoBuildNodalMass(Domain *theDomain, LadrunoContactDomain *cd)
{
    cd->clearNodalMass();
    // seed every node with its nodal `mass` (the Node::getMass diagonal, translational DOFs)
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0) {
        const Matrix &M = nodPtr->getMass();
        int nr = M.noRows();
        double m[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < 3 && d < nr; d++) m[d] = M(d, d);
        cd->setNodalMass(nodPtr->getTag(), m);
    }
    // add each element's diagonal mass contribution to its nodes
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0) continue;
        int ndof = e->getNumDOF();
        if (ndof <= 0 || ndof % nn != 0) continue;          // non-uniform DOF layout ⇒ cannot map
        int dpn = ndof / nn;                                // DOFs per node (3 solid, 6 beam/shell)
        const Matrix &Me = e->getMass();
        if (Me.noRows() != ndof || Me.noCols() != ndof) continue;   // no/odd mass matrix ⇒ skip
        for (int k = 0; k < nn; k++) {
            double m[3];
            if (!cd->getNodalMass(en(k), m)) { m[0] = m[1] = m[2] = 0.0; }
            for (int d = 0; d < 3 && d < dpn; d++) m[d] += Me(k * dpn + d, k * dpn + d);
            cd->setNodalMass(en(k), m);
        }
    }
}

// ============================================================ ADR-57 E2 edge-edge routing geometry
// Handler-local helpers porting the proto_e1_router oracle (16/16): the facet-normal align_fn regime
// split + the proximity-gated NTS ≻ edge-edge ≻ face-mortar precedence. Built from REFERENCE coords
// at handle() (the §7 scoping — large sliding needs re-emission, the inherited subsystem fence).

// Newell unit normal of a planar convex tri/quad facet (orientation-independent; align_fn = |nS·nM|).
static void ladrunoFacetNormalNewell(int npts, const double X[4][3], double n[3])
{
    n[0] = n[1] = n[2] = 0.0;
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nn > 1e-300) { n[0] /= nn; n[1] /= nn; n[2] /= nn; }
}

// does point x project in-bounds onto the convex facet X (npts) with |gap| ≤ band? gap by ref.
static bool ladrunoPointInFacet(const double x[3], int npts, const double X[4][3],
                                const double nf[3], double band, double &gap)
{
    double c[3] = {0, 0, 0};
    for (int i = 0; i < npts; i++) for (int d = 0; d < 3; d++) c[d] += X[i][d] / npts;
    gap = 0.0;
    for (int d = 0; d < 3; d++) gap += (x[d] - c[d]) * nf[d];
    double proj[3];
    for (int d = 0; d < 3; d++) proj[d] = x[d] - gap * nf[d];
    int sgn = 0;                                   // consistent cross-product sign vs each edge (convex)
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        double e[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        double w[3] = {proj[0]-a[0], proj[1]-a[1], proj[2]-a[2]};
        double cr[3] = {e[1]*w[2]-e[2]*w[1], e[2]*w[0]-e[0]*w[2], e[0]*w[1]-e[1]*w[0]};
        double dp = cr[0]*nf[0] + cr[1]*nf[1] + cr[2]*nf[2];
        int s = (dp > 1e-300) ? 1 : ((dp < -1e-300) ? -1 : 0);
        if (s != 0) { if (sgn == 0) sgn = s; else if (s != sgn) return false; }
    }
    return std::fabs(gap) <= band;
}

// min distance between two facets: vertices-projected-in-other + all boundary edge-edge closest
// points (the proto_e1 min_facet_dist). Reference coords ⇒ the proximity gate.
static double ladrunoMinFacetDist(int npS, const double S[4][3], const double nS[3],
                                  int npM, const double M[4][3], const double nM[3])
{
    double dmin = HUGE_VAL, gap;
    for (int i = 0; i < npS; i++)
        if (ladrunoPointInFacet(S[i], npM, M, nM, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int i = 0; i < npM; i++)
        if (ladrunoPointInFacet(M[i], npS, S, nS, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int ie = 0; ie < npS; ie++)
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            double g = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            if (g < dmin) dmin = g;
        }
    return dmin;
}

// d_band for a master facet: explicit -edgeBand override, else a fraction of the mean edge length.
static double ladrunoEdgeBand(int npM, const double M[4][3], double override_)
{
    if (override_ > 0.0) return override_;
    double el = 0.0;
    for (int k = 0; k < npM; k++) {
        const double *a = M[k]; const double *b = M[(k + 1) % npM];
        double dx = b[0]-a[0], dy = b[1]-a[1], dz = b[2]-a[2];
        el += std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    return 0.25 * (el / npM);
}

// The SINGLE SOURCE OF TRUTH for "edge-edge OWNS this (slave facet, master facet) pair" — the
// proto_e1 routing decision: in proximity, NON-FACING (align_fn < τ_face), NOT owned by an NTS
// slave-vertex-on-face poke, AND ≥1 boundary-edge pair is a well-conditioned margin-interior in-band
// crossing. Called by BOTH the mortar injection loop (which STANDS DOWN when this is true — the ADR §2
// "face-mortar vacates only if edge-edge succeeds", so the mortar clip and the edge-edge force never
// superpose in the oblique band where the clip is non-degenerate) AND the edge-edge injection loop
// (which injects only when this is true). One decision, two consumers ⇒ no double-count (the gate's
// EE-1 fix). Reference coords (§7 scoping).
static bool ladrunoEdgeEdgeOwns(int npS, const double S[4][3], int npM, const double M[4][3],
                                double dBand)
{
    const double TAU_FACE = std::cos(50.0 * 3.14159265358979323846 / 180.0);
    double nS[3], nM[3];
    ladrunoFacetNormalNewell(npS, S, nS);
    ladrunoFacetNormalNewell(npM, M, nM);
    if (ladrunoMinFacetDist(npS, S, nS, npM, M, nM) > dBand) return false;   // separated ⇒ NONE
    double alignFn = std::fabs(nS[0]*nM[0] + nS[1]*nM[1] + nS[2]*nM[2]);
    if (alignFn >= TAU_FACE) return false;                                   // facing ⇒ FACE-MORTAR owns
    double gap;
    for (int k = 0; k < npS; k++)                                           // slave vertex poke ⇒ NTS owns
        if (ladrunoPointInFacet(S[k], npM, M, nM, dBand, gap)) return false;
    for (int ie = 0; ie < npS; ie++)                                        // ≥1 in-band edge crossing
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            if (std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) <= dBand) return true;
        }
    return false;
}

int
LadrunoContactHandler::handle(const ID *nodesLast)
{
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();
    if (theDomain == 0 || theModel == 0 || theIntegrator == 0) {
        opserr << "WARNING LadrunoContactHandler::handle() - setLinks() not called\n";
        return -1;
    }

    // ADR-78 removal lane. handle() re-runs on EVERY domainChanged(), and a
    // collapse analysis (`recorder Collapse` / `remove node`) shrinks the domain
    // underneath a contact declared once at setup. Reconcile before anything
    // reads a surface: prune the removed nodes, retire any interaction left with
    // nothing to act on, both reported by name.
    //
    // This is what makes `recorder Collapse` + contact legal again. P1 had to
    // abort here because a missing node was indistinguishable from a deck typo;
    // it no longer is -- OPS_LadrunoContactSurface() rejects a non-existent tag
    // at DECLARATION time -- so absence now means the analysis removed it, and
    // the ladrunoSurfaceNodesOk() aborts below become genuinely unreachable for
    // that cause while staying as a backstop.
    if (LadrunoContactDomain *cdPrune = theDomain->getLadrunoContactDomain())
        cdPrune->pruneRemovedNodes(theDomain);

    // P1a: MP constraints (equalDOF / rigidLink / rigidDiaphragm) are NOT enforced
    // by the contact handler yet (Plain-style numbering). Warn rather than silently
    // ignore. Delegation to a base handler lands in P1b.
    MP_ConstraintIter &theMPcheck = theDomain->getMPs();
    if (theMPcheck() != 0)
        opserr << "WARNING LadrunoContactHandler::handle() - MP constraints (e.g. equalDOF) "
                  "are present but NOT enforced by the P1a contact handler; use a model with "
                  "SP (fix) constraints only, or 'constraints Transformation' for those models.\n";

    // collect SPs (domain + load-pattern), keyed by node. Like PlainHandler this is a
    // Plain-style handler: it REMOVES a constrained DOF but cannot enforce a NON-homogeneous
    // (imposed-displacement) SP value — that needs the Transformation/Penalty/Lagrange family,
    // which this handler can't be (it must inject the contact FE adapters in handle()). Warn
    // the same way PlainHandler does so a displacement-driven contact model fails LOUD, not
    // silently-zero; point the user to the DisplacementControl integrator (works with this
    // handler — it drives a DOF via the load factor, no imposed SP needed). (B3 follow-up.)
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0) {
        if (theSP->isHomogeneous() == false)
            opserr << "WARNING LadrunoContactHandler::handle() - non-homogeneous SP (imposed "
                      "displacement) at node " << theSP->getNodeTag() << " is NOT enforced by "
                      "the contact handler (Plain-style; homogeneous constraint assumed). Use the "
                      "DisplacementControl integrator for displacement-driven contact.\n";
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));
    }

    // DOF_Groups: free = -2, single-point-constrained = -1.
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    int numDOF = 0, countDOF = 0, count3 = 0;
    while ((nodPtr = theNod()) != 0) {
        DOF_Group *dofPtr = new DOF_Group(numDOF++, nodPtr);
        if (dofPtr == 0) {
            opserr << "WARNING LadrunoContactHandler::handle() - out of memory (DOF_Group)\n";
            return -4;
        }
        const ID &id = dofPtr->getID();
        for (int j = 0; j < id.Size(); j++) {
            dofPtr->setID(j, -2);
            countDOF++;
        }
        int nodeID = nodPtr->getTag();
        std::multimap<int, SP_Constraint *>::iterator first = allSPs.lower_bound(nodeID);
        std::multimap<int, SP_Constraint *>::iterator last = allSPs.upper_bound(nodeID);
        for (std::multimap<int, SP_Constraint *>::iterator it = first; it != last; it++) {
            int dof = it->second->getDOF_Number();
            const ID &cid = dofPtr->getID();
            if (dof >= 0 && dof < cid.Size() && cid(dof) == -2) {
                dofPtr->setID(dof, -1);
                countDOF--;
            } else {
                opserr << "WARNING LadrunoContactHandler::handle() - invalid/duplicate SP at DOF "
                       << dof << " node " << nodeID << endln;
            }
        }
        // contact-review P4 — EQ_Constraint parity with the CURRENT upstream PlainHandler
        // (added upstream after this replica was written; review HIGH-3): enforce a
        // trivial-identity EQ (single-entry constraint == 1.0 ⇒ mark the DOF -4, the
        // numberer's equation-constraint code) and LOUDLY warn-and-ignore any non-trivial
        // one. Without this, every equationConstraint — including the fork's own
        // LadrunoTie (ADR-62) rows — was SILENTLY dropped under `constraints
        // LadrunoContact`: parts of the structure ran unconnected with no diagnostic.
        // Actual EQ ENFORCEMENT stays the projection handler's job (`constraints
        // LadrunoProjection`) — contact + non-trivial ties in one analysis remain
        // unsupported, but now say so instead of silently producing a wrong answer.
        {
            EQ_ConstraintIter &theEQs = theDomain->getEQs();
            EQ_Constraint *eqPtr;
            while ((eqPtr = theEQs()) != 0) {
                if (eqPtr->getNodeConstrained() != nodeID)
                    continue;
                if (eqPtr->isTimeVarying() == true)
                    opserr << "WARNING LadrunoContactHandler::handle() - time-varying "
                              "EQ_Constraint for node " << nodeID << "; non-varying assumed\n";
                const Vector &C = eqPtr->getConstraint();
                if (C.Size() > 1 || C(0) != 1.0) {
                    opserr << "WARNING LadrunoContactHandler::handle() - EQ_Constraint at node "
                           << nodeID << " is not a trivial identity: NOT ENFORCED by the "
                              "contact handler (use `constraints LadrunoProjection` for "
                              "LadrunoTie / equationConstraint models; contact + non-trivial "
                              "ties in one analysis is unsupported)\n";
                } else {
                    int dof = eqPtr->getConstrainedDOFs();
                    const ID &cid = dofPtr->getID();
                    if (dof >= 0 && dof < cid.Size() && cid(dof) == -2) {
                        dofPtr->setID(dof, -4);
                        countDOF--;
                    } else {
                        opserr << "WARNING LadrunoContactHandler::handle() - EQ_Constraint DOF "
                               << dof << " already constrained at node " << nodeID << endln;
                    }
                }
            }
        }
        nodPtr->setDOF_GroupPtr(dofPtr);
        theModel->addDOF_Group(dofPtr);
    }
    theModel->setNumEqn(countDOF);

    // nodesLast (subdomain boundary): mark -3 (as PlainHandler)
    if (nodesLast != 0) {
        for (int i = 0; i < nodesLast->Size(); i++) {
            Node *np = theDomain->getNode((*nodesLast)(i));
            if (np != 0) {
                DOF_Group *dp = np->getDOF_GroupPtr();
                const ID &id = dp->getID();
                for (int j = 0; j < id.Size(); j++)
                    if (id(j) == -2) { dp->setID(j, -3); count3++; }
            }
        }
    }

    // FE_Elements for every Domain element (as PlainHandler). Track numFe so the
    // contact adapter tags continue ABOVE the structural FE tags (addFE_Element
    // silently drops a duplicate tag).
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    int numFe = 0;
    // ADR-60 R6 (D7 serial-only refusal): detect a partitioned host — this Domain is itself a
    // worker Subdomain (dynamic_cast) OR it owns Subdomain elements (the SP coordinator). Under
    // either, the finite-sliding migration trigger would re-sort uncoordinated across ranks, so
    // re-emit reverts to the shipped frozen-config NTS broad phase below (a one-time warning).
    // Sequential ⇒ false ⇒ every re-emit branch is inert ⇒ byte-identical to the shipped path.
    bool hostPartitioned = (dynamic_cast<Subdomain *>(theDomain) != 0);
    while ((elePtr = theEle()) != 0) {
        if (elePtr->isSubdomain() == true) {
            hostPartitioned = true;
            Subdomain *theSub = (Subdomain *)elePtr;
            if (theSub->doesIndependentAnalysis() == false) {
                FE_Element *fePtr = new FE_Element(numFe++, elePtr);
                if (fePtr == 0) return -5;
                theModel->addFE_Element(fePtr);
                theSub->setFE_ElementPtr(fePtr);
            }
        } else {
            FE_Element *fePtr = new FE_Element(numFe++, elePtr);
            if (fePtr == 0) return -5;
            theModel->addFE_Element(fePtr);
        }
    }

    // --- inject the contact FE adapter(s) ---
    // None if no contact engine is attached (pure Plain -> byte-identical to stock).
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();

    // B1 (P4): if ANY NTS contact / rigid plane requests the SOFT=1 penalty, pre-compute the
    // assembled translational nodal-mass cache (nodal `mass` + element-density mass) the adapters
    // read to size k_soft = SOFSCL·4·m_eff/dt². Only when soft is in play ⇒ no cost (and byte-
    // identical) otherwise. Mortar SOFT is out of scope (refused at the command surface).
    if (cd != 0) {
        // ADR-78 P2 review F1: skip RETIRED interactions here too. This flag
        // feeds BOTH the mass-cache build below and the partitioning refusal --
        // a retired soft interaction is inert by the removal lane's own
        // contract ("the run continues WITHOUT it"), so it must neither
        // trigger the cache pass nor MPI_Abort a partitioned job in which no
        // LIVE soft contact exists (pre-fix that abort printed a FATAL with
        // no subject, because the message loops below correctly skip retired).
        bool anySoft = false;
        for (int c = 0; c < cd->getNumContacts() && !anySoft; c++)
            if (!cd->getContact(c).retired &&
                cd->getContact(c).softScale > 0.0) anySoft = true;
        for (int p = 0; p < cd->getNumRigidPlanes() && !anySoft; p++)
            if (!cd->getRigidPlane(p).retired &&
                cd->getRigidPlane(p).softScale > 0.0) anySoft = true;
        // B2 (P5): a SOFT=2 mortar contact also needs the assembled nodal-mass cache (the segment-
        // based m_eff). The cache is over ALL nodes (one element pass), shared with the NTS SOFT=1 lane.
        // ADR-57 E5: a -edgeSoft edge-edge contact needs the SAME cache (the edge gap-mode m_eff).
        for (int c = 0; c < cd->getNumMortarContacts() && !anySoft; c++)
            if (!cd->getMortarContact(c).retired &&
                (cd->getMortarContact(c).softScale > 0.0 ||
                 cd->getMortarContact(c).edgeSoftScale > 0.0)) anySoft = true;
        // ADR-78 P2 (D4): SOFT is REFUSED under partitioning, every lane. k_soft =
        // SOFSCL*4*m_eff/dt^2 reads the ASSEMBLED nodal mass, and this rank's cache
        // holds only ITS elements' mass: a ghosted surface contributes zero, and a
        // partition-BOUNDARY node contributes a PARTIAL mass -- nonzero, so the P1
        // zero-mass guard below cannot see it -- making m_eff too small and the
        // contact silently too soft (the D4 failure class, undetectable locally).
        // np > 1 comes from the per-target MPI TU (LadrunoContactAbort.cpp); an
        // #ifdef probe HERE would be dead code in every build, OpenSeesMP included
        // (the P1 inert-MPI_Abort trap -- see the note at the top of this file).
        //
        // Disclosed over-refusal: an MP parametric sweep whose EVERY rank holds the
        // full model has a complete mass cache and would size SOFT correctly, but
        // that is indistinguishable from manual domain decomposition from inside
        // handle() -- so it refuses too. Named in the message.
        //
        // ADR-85 T2 (D6) ordering fix: this whole -soft pre-scan still runs BEFORE the
        // ADR-85 dimension/pairing gates in the main lane loops below, but the Contact
        // (NTS) mass-check loop now runs `ladrunoContactPairPreflight` on each
        // softScale>0 contact BEFORE inspecting its slave nodes' mass (see below) -- so
        // a deck that is refusable for a dimension/guard reason (mixed dimensions,
        // ndf!=ndm, ...) is refused by THAT named guard first, never mis-attributed to
        // this generic "-soft needs an assembled mass" complaint reached first. (T0
        // left this deliberately alone -- "revisit when T1b rebuilds this loop" -- T1b
        // deferred it again by name as D6, "T2 owns it alongside the loop rebuild for
        // friction"; this is that fix, scoped to the Contact/NTS lane the T0 note named.
        // RigidPlane/Mortar are unchanged here -- neither was named by the T0 note, and
        // this pass did not re-audit their own -soft mass-check loops below.)
        if (anySoft && (hostPartitioned || ladrunoContactNumRanks() > 1)) {
            // Review F4: each offender is its own newline-terminated line (a
            // multi-offender abort used to concatenate the fragments into one
            // run-on string -- the sole diagnostic a torn-down job leaves).
            for (int c = 0; c < cd->getNumContacts(); c++)
                if (!cd->getContact(c).retired && cd->getContact(c).softScale > 0.0)
                    opserr << "FATAL LadrunoContactHandler::handle() - contact "
                           << cd->getContact(c).tag << ": -soft refused under partitioning\n";
            for (int p = 0; p < cd->getNumRigidPlanes(); p++)
                if (!cd->getRigidPlane(p).retired && cd->getRigidPlane(p).softScale > 0.0)
                    opserr << "FATAL LadrunoContactHandler::handle() - contactPlane "
                           << cd->getRigidPlane(p).tag << ": -soft refused under partitioning\n";
            for (int c = 0; c < cd->getNumMortarContacts(); c++) {
                const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
                if (mc.retired) continue;
                if (mc.softScale > 0.0)
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                           << mc.tag << ": -soft refused under partitioning\n";
                if (mc.edgeSoftScale > 0.0)
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                           << mc.tag << ": -edgeSoft refused under partitioning\n";
            }
            opserr << "FATAL LadrunoContactHandler::handle() - SOFT "
                      "is not supported under a partitioned host (np > 1 or "
                      "Subdomain): k_soft = SOFSCL*4*m_eff/dt^2 needs BOTH surfaces' "
                      "ASSEMBLED mass, and this process's cache holds only its own "
                      "elements' share -- a ghosted or partition-boundary node makes "
                      "m_eff silently too small (ADR-78 P2/D4). Use an explicit -kn "
                      "(or -kn auto) instead, or run the SOFT model serial. NOTE: "
                      "this refusal also fires when every rank holds the FULL model "
                      "(an MP parametric sweep) -- completeness of the mass cache is "
                      "not verifiable from one rank; ABORTING (ADR-78 P2)\n";
            return ladrunoContactFatal();
        }
        if (anySoft) {
            ladrunoBuildNodalMass(theDomain, cd);
            // ADR-78 P1, NARROWED after adversarial review. A SOFT contact sizes
            // k_soft = SOFSCL*4*m_eff/dt² from the ASSEMBLED nodal mass built above,
            // and m_eff = 1/(invMproj_slave + Σ Nᵢ²·invMproj_segᵢ).
            //
            // The first version of this guard scanned BOTH surfaces and justified
            // itself with "k_soft collapses to 0 and the contact silently does
            // nothing". BOTH halves of that were wrong, and the adapter says so:
            //
            //  * A massless MASTER is legitimate and supported. gapModeInvMass()
            //    (LadrunoContactFE.cpp:448) treats a DOF with m ≤ 0 as INFINITE mass
            //    -- "the correct gap-mode treatment of a FIXED master node ... a
            //    fixed/rigid master contributes 0, leaving m_eff = m_slave". A rigid
            //    indenter or drum meshed as fixed, element-free master nodes is a
            //    standard pattern; aborting on it rejected a correct model.
            //  * The failure was never silent. softKn() (LadrunoContactFE.cpp:516)
            //    already warns once per (contact, topic) and falls back to the
            //    configured kn, which the kn > 0 guard below guarantees is positive.
            //
            // Only a massless SLAVE is degenerate -- gapModeInvMass returns 0 and the
            // adapter cannot soft-size at all. That is also the partitioned case ADR
            // 0092 INV-3 refuses `soft=` for: a ghosted slave node has no element on
            // the owner rank. So: slave surface only, and the message no longer
            // claims the alternative is silence.
            for (int c = 0; c < cd->getNumContacts(); c++) {
                const LadrunoContactDomain::Contact &ct = cd->getContact(c);
                if (ct.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
                if (ct.softScale <= 0.0) continue;
                const LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
                if (ss == 0) continue;   // surface-undefined is caught, fatally, below
                // ADR-85 T2 (D6): run the SAME dimension/ndf/pairing guard the main
                // NTS loop runs (below) BEFORE this generic mass check, so a deck
                // that is ALSO dimension-broken is named correctly and first. Purely
                // read-only (inspects Domain state, no engine mutation), so calling
                // it here and again in the main loop is redundant but harmless on a
                // valid deck; on an already-fatal deck it just changes which of the
                // two (equivalent, both-correct) messages fires first.
                const LadrunoContactSurface *msPre = cd->getSurface(ct.masterSurfTag);
                if (msPre != 0) {
                    int pairDimPre = 3;
                    if (!ladrunoContactPairPreflight(theDomain, ct.tag, msPre, ss,
                                                     "NTS node-to-segment", "T1b",
                                                     /*lane2DLive=*/true, &pairDimPre))
                        return ladrunoContactFatal();
                }
                const ID &tg = ss->getNodeTags();
                for (int i = 0; i < tg.Size(); i++) {
                    double m[3];
                    bool have = cd->getNodalMass(tg(i), m);
                    if (!have || (m[0] <= 0.0 && m[1] <= 0.0 && m[2] <= 0.0)) {
                        opserr << "FATAL LadrunoContactHandler::handle() - contact "
                               << ct.tag << ": -soft needs an assembled mass at SLAVE node "
                               << tg(i) << ", which carries none (no nodal `mass`, no "
                                  "element on this process). m_eff collapses, so the "
                                  "contact silently falls back to the configured kn and "
                                  "loses its Courant-stable penalty; ABORTING (ADR-78 P1). "
                                  "Under partitioning `-soft` is not supported. A massless "
                                  "MASTER is fine -- it is treated as a fixed/rigid one.\n";
                        return ladrunoContactFatal();
                    }
                }
            }
        }
    }

    // P2b: faceted node-to-segment penalty contact. For each contact definition
    // (MASTER_SEGMENTS surface vs SLAVE_NODES surface) build ONE bound adapter per
    // (slave node, master segment) pair — the gate-sanctioned brute-force pairing
    // (bucket-sort broad phase = P2.5). A degenerate/non-projecting segment yields
    // zero force inside the adapter (so a topology-only contact stays graph-neutral
    // under a connectivity-independent numberer — the P1b regression relies on this).
    if (cd != 0) {
        // P3: rebuild the live friction-state key-set this handle() so the engine can
        // GC slots from a previous (re-meshed/re-paired) analysis. Mark every built
        // frictional pair; sweep at the end.
        cd->frictionGCBegin();
        cd->clearReemit();   // ADR-60: drop last epoch's re-emit anchors (re-registered per opted-in contact below)
        cd->clearNormalFields();  // ADR-63 #4a: drop last handle's nodal-normal fields (rebuilt per opted-in contact below)
        cd->clearVtx2DSide();     // ADR-85 T1b: vertex side signs are per-PAIRING-EPOCH state (sigma is
                                  // fixed at pairing and re-derived on re-pairing, How/2); a 3D deck's
                                  // map is always empty, so this is a no-op there.

        // P3 cross-contact shared-slave warning: a slave node in two contacts' slave
        // sets gets two adapters → double traction (a modeling error the engine can't
        // auto-resolve). Cheap multiset scan, mirrors the ndf!=3 warning.
        {
            std::multimap<int, int> slaveSeen;   // slaveNodeTag -> contactTag
            for (int c = 0; c < cd->getNumContacts(); c++) {
                const LadrunoContactDomain::Contact &ct = cd->getContact(c);
                if (ct.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
                LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
                if (ss == 0 || ss->getKind() != LadrunoContactSurface::SLAVE_NODES)
                    continue;
                const ID &st = ss->getNodeTags();
                for (int i = 0; i < st.Size(); i++) {
                    int tag = st(i);
                    if (slaveSeen.find(tag) != slaveSeen.end())
                        opserr << "WARNING LadrunoContactHandler::handle() - slave node "
                               << tag << " appears in more than one contact's slave set; "
                                  "its contact tractions will be ADDED (check the model)\n";
                    slaveSeen.insert(std::make_pair(tag, ct.tag));
                }
            }
        }

        for (int c = 0; c < cd->getNumContacts(); c++) {
            const LadrunoContactDomain::Contact &ct = cd->getContact(c);
            if (ct.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
            LadrunoContactSurface *ms = cd->getSurface(ct.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": master/slave surface not defined; ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_NODES) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": need a MASTER_SEGMENTS master + SLAVE_NODES slave; "
                          "ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            int nps = ms->getNodesPerSeg();
            // ADR-85 T1b -- the pre-flight now runs ONCE, in LIVE mode, before the
            // arity gate: the same existence / per-surface dimension / cross-
            // dimension / ndf-equality checks as T0 (same functions, same message
            // text), but a well-formed all-2D pair PROCEEDS into the 2D lane below
            // instead of drawing the T0 "lands in T1b" refusal -- that refusal is
            // what this phase replaces, and it is the ONLY T0 message the NTS lane
            // loses. An all-3D pair passes silently (pairDim is a write-only out)
            // and takes the shipped path below, textually untouched. (T0's two call
            // sites -- the nps == 2 early block and the general one after the arity
            // gate, both closing the arity-free-slave hole -- collapse into this
            // single earlier call; on a VALID deck all of them are silent passes,
            // so the collapse is observable only on already-fatal decks, where at
            // most the message ORDER changes.)
            int pairDim = 3;
            if (!ladrunoContactPairPreflight(theDomain, ct.tag, ms, ss,
                                             "NTS node-to-segment", "T1b",
                                             /*lane2DLive=*/true, &pairDim))
                return ladrunoContactFatal();
            if (pairDim == 2) {
                // ==============================================================
                // ADR-85 T1b -- the LIVE 2D NTS lane (frictionless penalty).
                // A 3D deck CANNOT enter: pairDim == 2 requires a surface whose
                // every node carries 2-component coordinates (pre-flight above).
                // Everything in this block is new 2D code ending in `continue`,
                // so the shipped 3D path below never sees a 2D node and stays
                // textually untouched. Wired VERBATIM to the
                // LadrunoContact2DKernel contract (How/1 + the T1a corrections):
                // CONCAVE-only vertex pairs, ORDERED ownership (segment
                // precedence in surface order, then the UNSLACKED vertex claim
                // aPrev > 0 && aNext < 0 -- evaluated inside the adapters), and
                // the bisector side sign committed at first capture. The
                // interface-level orientation vote (How/2) replaces the 3D
                // default's per-pair datum -- never inherited (its ridge-flip +
                // silent-flush-drop hazards are on record).
                // ==============================================================
                if (nps != 2) {
                    // only reachable through a stream restore (the parser pins
                    // arity <-> dimension at declaration); refuse by name --
                    // T0's parser rule, enforced again at handle() time.
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": a 2D NTS contact needs nodesPerSeg == 2 (got " << nps
                           << "); ABORTING (ADR-85 nodesPerSeg)\n";
                    return ladrunoContactFatal();
                }
                // ---- capability fence: what T1b did NOT wire was refused BY
                //      NAME, never silently ignored (the T0 discipline). ADR-85
                //      T2 turns 2D NTS friction (-mu) into the LIVE scalar
                //      return map (LadrunoFrictionKernel::returnMap1D via
                //      LadrunoContactFE::getResidual/addKtToTang/addKiToTang) --
                //      the retired FATAL used to sit here; ct.kt/ct.mu/
                //      ct.consistentTan are threaded into the 2D adapters below.
                //      mu<=0 stays byte-identical to the T1b frictionless lane
                //      (the adapter's own short-circuit, mirroring the 3D
                //      SEGMENT `mu > 0.0` discipline -- no new silence).
                //      -geomtan (ct.consistentNormal) is still an ACCEPTED no-op
                //      (see LadrunoContactFE::addKtToTang's 2D branch), announced
                //      once by the WARN_2D_GEOMTAN_NOOP latch below.
                //      -consistanttan's non-symmetric-solver warning is emitted
                //      at PARSE time (OpenSeesOutputCommands.cpp), dimension-
                //      independent, so 2D needs no extra handle()-time latch --
                //      exactly the 3D NTS lane's own contract (only the MORTAR
                //      lane additionally re-warns per contact at handle() time).
                if (ct.smoothNormal) {
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": -smoothNormal is 3D-only (ADR-63 nodal-normal smoothing is "
                              "built on Newell facet normals) and is fenced out of the 2D "
                              "lane; the 2D corner coverage is the vertex policy, not "
                              "smoothing. Drop -smoothNormal. ABORTING (ADR-85 fence)\n";
                    return ladrunoContactFatal();
                }
                if (ct.enableReemit) {
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": -reemit (finite-sliding re-emission, ADR-60) is 3D-only and "
                              "fenced out of the 2D lane (ADR-85 not-in-scope). Drop -reemit. "
                              "ABORTING (ADR-85 fence)\n";
                    return ladrunoContactFatal();
                }
                // panel fix: -geomtan is ACCEPTED as a documented no-op in 2D,
                // but accepted-and-SILENT hides the deferral from the user --
                // surface it once per contact (the WARN_VTX2D_DEFER idiom).
                if (ct.consistentNormal &&
                    cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_2D_GEOMTAN_NOOP)) {
                    opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": -geomtan is ACCEPTED AS A NO-OP on the 2D NTS lane -- the "
                              "2D B-operator is already the exact FIRST variation of the "
                              "gap, and the geometric curvature term is a named deferral "
                              "(ADR-85 How/5); the assembled tangent is kn*B^T*B with or "
                              "without the flag.\n";
                }
                if (!ct.knAuto && ct.kn <= 0.0) {
                    // the 2D twin of the shipped kn > 0 gate below (same rationale:
                    // dead adapters are a declared contact transmitting nothing).
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": segment contact needs kn > 0 (got " << ct.kn << "); "
                              "ABORTING (ADR-78 P1)\n";
                    return ladrunoContactFatal();
                }
                const ID &mTags = ms->getNodeTags();
                const ID &sTags = ss->getNodeTags();
                int nSeg = mTags.Size() / 2;

                // ---- Lref: the master surface's MIN INITIAL segment length ----
                //      The kernel's caller contract: every relative gauge
                //      (tauSeg zero-length refusal, vertexEval2D's tolLen) hangs
                //      off it, computed ONCE per contact at handle(). Nodes
                //      exist (pre-flight), so the derefs are safe.
                double Lref = 0.0;
                int    LrefBadSeg = -1;
                for (int seg = 0; seg < nSeg; seg++) {
                    const Vector &A = theDomain->getNode(mTags(2 * seg))->getCrds();
                    const Vector &Bc = theDomain->getNode(mTags(2 * seg + 1))->getCrds();
                    double dx = Bc(0) - A(0), dy = Bc(1) - A(1);
                    double L = std::sqrt(dx * dx + dy * dy);
                    if (!(L > 0.0)) { LrefBadSeg = seg; break; }
                    if (Lref <= 0.0 || L < Lref) Lref = L;
                }
                if (nSeg <= 0 || LrefBadSeg >= 0 || !(Lref > 0.0)) {
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": 2D master surface " << ct.masterSurfTag;
                    if (nSeg <= 0)
                        opserr << " has no segments";
                    else
                        opserr << " segment " << LrefBadSeg << " has zero reference length "
                                  "(the kernel's relative zero-length gauge Lref would be "
                                  "disabled)";
                    opserr << "; ABORTING (ADR-85 T1b)\n";
                    return ladrunoContactFatal();
                }

                // ---- CHAIN-INTEGRITY scan (panel MAJOR fix) ------------------
                //      Surface order is the 2D adjacency contract: the ordered-
                //      ownership stand-down and the concave vertex pairs key off
                //      CONSECUTIVE head-to-tail links. That contract was
                //      previously assumed, not validated -- a permuted listing
                //      reproduces the exact T1a double-owner closed form, and a
                //      reversed one silently disarms BOTH the stand-down and the
                //      vertex pairs (probe-confirmed). O(nSeg) rule: a node tag
                //      shared by two segments is legal ONLY as segment k's
                //      SECOND node and segment (k+1)%nSeg's FIRST node (the
                //      wrap covers a closed loop); any other sharing (permuted /
                //      reversed / forked / >2 uses) is a named FATAL.
                {
                    std::map<int, std::vector<std::pair<int, int> > > occ;
                    for (int seg = 0; seg < nSeg; seg++) {
                        occ[mTags(2 * seg)].push_back(std::make_pair(seg, 0));
                        occ[mTags(2 * seg + 1)].push_back(std::make_pair(seg, 1));
                    }
                    for (std::map<int, std::vector<std::pair<int, int> > >::const_iterator
                             it = occ.begin(); it != occ.end(); ++it) {
                        const std::vector<std::pair<int, int> > &v = it->second;
                        if (v.size() <= 1) continue;         // terminal / unshared
                        bool okLink = false;
                        if (v.size() == 2 && v[0].second != v[1].second) {
                            int tailSeg = (v[0].second == 1) ? v[0].first : v[1].first;
                            int headSeg = (v[0].second == 1) ? v[1].first : v[0].first;
                            if (nSeg > 1 && headSeg == (tailSeg + 1) % nSeg)
                                okLink = true;
                        }
                        if (!okLink) {
                            opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                   << ct.tag << ": 2D master surface " << ct.masterSurfTag
                                   << " node " << it->first << " is shared by segments "
                                   << v[0].first << " and " << v[1].first
                                   << " without being a consecutive head-to-tail link "
                                      "(segment k's SECOND node must be segment k+1's "
                                      "FIRST node in declaration order; a closed loop "
                                      "wraps last-to-first). Surface order is the 2D "
                                      "adjacency contract -- the ownership stand-down and "
                                      "the concave vertex pairs key off it, and a "
                                      "permuted/reversed/forked listing would silently "
                                      "re-open the corner double-count. Re-order and "
                                      "re-wind the master segment list head-to-tail. "
                                      "ABORTING (ADR-85 chain integrity)\n";
                            return ladrunoContactFatal();
                        }
                    }
                }

                // ---- the INTERFACE-LEVEL orientation vote (How/2) ------------
                //      ADR-85 T3: extracted to the file-static
                //      ladruno2DOrientationVote() (pure code motion, exact
                //      operation order + message strings preserved) so the T3
                //      mortar branch can share it.
                double sigma = 0.0;
                if (!ladruno2DOrientationVote(theDomain, ct.tag, ct.hasOutward, ct.outward,
                                              mTags, nSeg, sTags, Lref, sigma))
                    return ladrunoContactFatal();

                // ---- CONCAVE-vertex candidates (How/1, the T1a decision) -----
                //      Vertices shared by CONSECUTIVE segments in surface order
                //      (mTags(2i+1) == mTags(2(i+1)), with wrap-around so a
                //      closed loop gets its seam vertex too). Turn sign from
                //      sigma at the REFERENCE config; only CONCAVE (turnSin > 0
                //      == C2D_CONCAVE) corners get a vertex pair -- the T1a
                //      oracle proved a convex vertex pair can NEVER carry force
                //      (its wedge is strictly exterior) while the concave wedge
                //      is the real, entirely-penetrating coverage hole.
                //      SURFACE ORDER IS THE ADJACENCY CONTRACT, and the chain-
                //      integrity scan above ENFORCES it: non-consecutive sharing
                //      is refused by name rather than silently disarming the
                //      stand-down and the vertex pairs -- which, with the
                //      open-end acceptance window, would be an OVERLAP (a
                //      double-count), not a hole.
                struct LadrunoVtx2D {
                    int segPrev, segNext;              // adjacent segment ordinals
                    Node *prevFar, *vtx, *nextFar;     // X_prev0, Xv, X_next1
                };
                std::vector<LadrunoVtx2D> vtx2d;
                for (int seg = 0; nSeg >= 2 && seg < nSeg; seg++) {
                    int nxt = (seg + 1) % nSeg;
                    if (nxt == seg) break;
                    if (mTags(2 * seg + 1) != mTags(2 * nxt)) continue;  // not chained here
                    Node *pF = theDomain->getNode(mTags(2 * seg));
                    Node *vN = theDomain->getNode(mTags(2 * seg + 1));
                    Node *nF = theDomain->getNode(mTags(2 * nxt + 1));
                    const Vector &XP = pF->getCrds();
                    const Vector &XV = vN->getCrds();
                    const Vector &XN = nF->getCrds();
                    const double tPrev[2] = { XV(0) - XP(0), XV(1) - XP(1) };
                    const double tNext[2] = { XN(0) - XV(0), XN(1) - XV(1) };
                    // sign(turn) == the kernel's CornerType seen from the contact
                    // side (normalization dropped -- a sign-only decision).
                    double turn = sigma * (tPrev[0] * tNext[1] - tPrev[1] * tNext[0]);
                    if (turn > 0.0) {                  // C2D_CONCAVE
                        LadrunoVtx2D v;
                        v.segPrev = seg; v.segNext = nxt;
                        v.prevFar = pF; v.vtx = vN; v.nextFar = nF;
                        vtx2d.push_back(v);
                    }
                }

                // ---- broad phase: the STRIDE-3 bucket grid, z = 0 padded -----
                //      LadrunoContactBucketSort stays stride-3 and is reused
                //      UNTOUCHED (ADR-85 Where); the handler feeds it the 2D
                //      reference coords with the third component pinned to 0 --
                //      its degenerate-axis handling collapses NZ to 1. Nodes
                //      exist (pre-flight), so the shipped missing-node backfill
                //      idiom is unnecessary here.
                std::vector<double> segCoords((size_t)nSeg * 2 * 3, 0.0);
                for (int seg = 0; seg < nSeg; seg++) {
                    double *S = &segCoords[(size_t)seg * 2 * 3];
                    for (int k = 0; k < 2; k++) {
                        const Vector &X = theDomain->getNode(mTags(seg * 2 + k))->getCrds();
                        S[k * 3 + 0] = X(0); S[k * 3 + 1] = X(1); S[k * 3 + 2] = 0.0;
                    }
                }
                LadrunoContactBucketSort::Grid grid(nSeg, 2, segCoords.data(),
                                                    ct.cellFrac, 1.0);
                std::vector<int> cand(nSeg);

                // ---- per-slave pair loop -------------------------------------
                //      Ownership note (the panel's checklist question): the
                //      handler builds ONE stateless adapter per (slave, candidate
                //      segment) and per (slave, near concave vertex); the
                //      per-slave ORDERED-OWNERSHIP decision itself is evaluated
                //      at the TRIAL config inside LadrunoContactFE::
                //      segment2DActive on every residual/tangent evaluation --
                //      it cannot live at this loop level because ownership moves
                //      with the trial configuration between two handle() calls.
                //      This loop's job is to ARM it: each segment adapter gets
                //      its surface-order predecessor's far node, each vertex
                //      adapter its two adjacent far nodes.
                for (int si = 0; si < sTags.Size(); si++) {
                    Node *sn = theDomain->getNode(sTags(si));
                    if (sn == 0 || sn->getNumberDOF() != 2) {
                        // the 2D twin of the shipped != 3 backstop (the
                        // pre-flight's ndf-equality gate already refused this;
                        // keep the same defense-in-depth the 3D loop has).
                        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                               << " slave node " << sTags(si) << " ndf="
                               << (sn != 0 ? sn->getNumberDOF() : 0)
                               << " != 2; ABORTING (2D NTS is 2D translational -- ADR-85 "
                                  "ndf equality contract)\n";
                        return ladrunoContactFatal();
                    }
                    const Vector &Xs0 = sn->getCrds();
                    double slavePt[3] = { Xs0(0), Xs0(1), 0.0 };
                    int nCand = grid.candidates(slavePt, cand.data());
                    for (int ci = 0; ci < nCand; ci++) {
                        int seg = cand[ci];
                        Node *segNodes[2];
                        // panel fix: a bool sentinel, NOT `badTag != 0` -- node
                        // tag 0 is legal in OpenSees, and the tag-keyed idiom
                        // would silently skip the refusal for it.
                        bool badNode = false;
                        int badTag = 0, badNdf = 0;
                        for (int k = 0; k < 2; k++) {
                            Node *mn = theDomain->getNode(mTags(seg * 2 + k));
                            if (mn == 0 || mn->getNumberDOF() != 2) {
                                badNode = true;
                                badTag = mTags(seg * 2 + k);
                                badNdf = (mn != 0) ? mn->getNumberDOF() : -1;
                                break;
                            }
                            segNodes[k] = mn;
                        }
                        if (badNode) {
                            opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                   << ct.tag << " master segment " << seg << " node "
                                   << badTag << ": ndf=" << badNdf
                                   << " != 2 (or node absent); ABORTING (2D NTS master "
                                      "segments are 2D translational -- ADR-85 ndf "
                                      "equality contract)\n";
                            return ladrunoContactFatal();
                        }
                        double knUse = ct.kn;
                        if (ct.knAuto) {
                            knUse = ladrunoResolveAutoKn2D(theDomain, segNodes, sigma);
                            if (knUse <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                       << ct.tag << " slave node " << sTags(si) << " segment "
                                       << seg << ": -kn auto could not size a penalty (no "
                                          "owning 2-DOF/node solid element / degenerate "
                                          "segment); ABORTING (ADR-78 P1). Give an explicit "
                                          "-kn, or ensure the master segment's plane element "
                                          "is in the model.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        // the PREDECESSOR segment in surface order (shares this
                        // segment's first node): its far node arms the adapter's
                        // ordered-ownership stand-down (kernel rule step 1 -- the
                        // previous segment owns whenever IT projects the slave
                        // in-bounds). No chained predecessor => no stand-down,
                        // and the X0 side is an OPEN TERMINAL end. The SUCCESSOR
                        // far node mirrors it for the X1 side: non-null <=>
                        // chained (strict seam); null <=> an open end, where the
                        // adapter applies the NTS2D_END_SLACK acceptance window
                        // (the tilt-drift discontinuity fix -- see
                        // LadrunoContactFE::segment2DActive).
                        Node *prevFar = 0, *nextFarSeg = 0;
                        if (nSeg > 1) {
                            int prv = (seg + nSeg - 1) % nSeg;
                            if (prv != seg && mTags(2 * prv + 1) == mTags(2 * seg))
                                prevFar = theDomain->getNode(mTags(2 * prv));
                            int nx2 = (seg + 1) % nSeg;
                            if (nx2 != seg && mTags(2 * seg + 1) == mTags(2 * nx2))
                                nextFarSeg = theDomain->getNode(mTags(2 * nx2 + 1));
                        }
                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sn, segNodes, /*nps=*/2, /*ndm2=*/2, knUse,
                            sigma, Lref, prevFar, nextFarSeg,
                            ct.kt, ct.mu, ct.consistentTan,  // ADR-85 T2: the scalar friction lane
                            ct.muc,        // T1b SEGMENT dashpot port (0 => off)
                            ct.softScale,  // B1 SOFT=1 (free in 2D: size-safe mass helpers)
                            theDomain, ct.tag, seg);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        // ADR-85 T2: mark this segment pair's friction slot LIVE this
                        // handle() (the frictionGCBegin/…/frictionGCEnd sweep at the
                        // end of the NTS loop), mirroring the 3D loop's `ct.mu > 0.0`
                        // mark below. mu<=0 ⇒ no slot is ever created (getResidual's
                        // own short-circuit) ⇒ nothing to mark ⇒ identical to T1b.
                        if (ct.mu > 0.0)
                            cd->frictionGCMark(ct.tag, sTags(si), seg);
                    }
                    // one VERTEX adapter per (slave, concave vertex) whose
                    // adjacent segments are broad-phase candidates of this slave
                    // (superset-preserving: the vertex lies on both segments, so
                    // any slave near the vertex has at least one of them in its
                    // 27-neighbour candidate set). segIndex = nSeg + segPrev
                    // keys the force snapshot / side-sign state uniquely and
                    // rebuild-stably (disjoint from the real segment ordinals).
                    for (size_t vi = 0; vi < vtx2d.size(); vi++) {
                        const LadrunoVtx2D &v = vtx2d[vi];
                        bool nearV = false;
                        for (int ci = 0; ci < nCand && !nearV; ci++)
                            if (cand[ci] == v.segPrev || cand[ci] == v.segNext)
                                nearV = true;
                        if (!nearV) continue;
                        double knV = ct.kn;
                        if (ct.knAuto) {
                            // size the vertex penalty from the PREVIOUS adjacent
                            // segment (deterministic; both segments share Xv).
                            Node *prevSeg[2] = { v.prevFar, v.vtx };
                            knV = ladrunoResolveAutoKn2D(theDomain, prevSeg, sigma);
                            if (knV <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                       << ct.tag << " slave node " << sTags(si)
                                       << " concave vertex node " << v.vtx->getTag()
                                       << ": -kn auto could not size the vertex-pair penalty "
                                          "(no owning 2-DOF/node solid element on the "
                                          "adjacent segment); ABORTING (ADR-78 P1). Give an "
                                          "explicit -kn.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        Node *vNodes[1] = { v.vtx };
                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sn, vNodes, /*nps=*/1, /*ndm2=*/2, knV,
                            sigma, Lref, v.prevFar, v.nextFar,
                            ct.kt, ct.mu, ct.consistentTan,  // ADR-85 T2: the scalar friction lane
                            ct.muc, ct.softScale, theDomain, ct.tag,
                            /*segIndex=*/nSeg + v.segPrev);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        // ADR-85 T2: vertex pairs key their friction slot under the
                        // SAME segIndex = nSeg + segPrev convention setNtsForce/
                        // getOrCreateFrictionState already use (ADR-85 What/T1b) --
                        // mark it live so the GC sweep does not drop it.
                        if (ct.mu > 0.0)
                            cd->frictionGCMark(ct.tag, sTags(si), nSeg + v.segPrev);
                    }
                }
                continue;   // 2D contact done -- the shipped 3D code below never sees it
            }
            if (nps < 3 || nps > 4) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": nodesPerSeg " << nps << " unsupported (need 3 or 4); "
                          "ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            if (!ct.knAuto && ct.kn <= 0.0) {
                // A SEGMENT contact needs a positive penalty (P2b-1 requires -kn).
                // kn == 0 (e.g. `contact ... -outward` with the kn omitted) is inert —
                // it would build dead adapters, i.e. a declared contact that transmits
                // nothing, which ADR-78 P1 makes fatal rather than a warning.
                // `-kn auto` (P2b-2b) carries a 0 placeholder and is resolved per pair
                // below, so it bypasses this guard.
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": segment contact needs kn > 0 (got " << ct.kn << "); "
                          "ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSeg = mTags.Size() / nps;

            // ADR-60 R6 (serial-only): one-time warning when -reemit is requested under a
            // partitioned host. The deformed feed + trigger are disabled below (reemitActive=false
            // ⇒ the shipped frozen-config NTS path), so this only surfaces the silent downgrade.
            if (ct.enableReemit && hostPartitioned &&
                cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_REEMIT_PARTITION)) {
                opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": -reemit is serial-only but the host is partitioned (Subdomain); "
                          "finite-sliding re-emit DISABLED (reverting to the frozen-config NTS "
                          "broad phase) to avoid an uncoordinated cross-rank re-sort.\n";
            }
            // ADR-60 R1 (BLOCKER-MEMBERSHIP): fingerprint this contact's master mTags ORDERING. If
            // it changed since the last handle (element removal / re-mesh reordered the surface),
            // the segment-ordinal friction key would alias onto a different physical segment — drop
            // this contact's friction slots (they re-engage fresh, traction-continuous per D4) and
            // refuse to ARM the trigger this step (named error). Only evaluated when -reemit is on
            // ⇒ OFF leaves theReemitFp untouched ⇒ byte-identical.
            bool membershipChanged = false;
            if (ct.enableReemit) {
                std::vector<int> mt((size_t)mTags.Size());
                for (int i = 0; i < mTags.Size(); i++) mt[i] = mTags(i);
                unsigned long long fp = LadrunoContactReemit::membershipFingerprint(
                    mt.empty() ? 0 : &mt[0], (int)mt.size(), nps);
                if (cd->reemitMembershipChanged(ct.tag, fp)) {
                    membershipChanged = true;
                    opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": master surface node ordering changed since the last handle; the "
                              "segment-keyed re-emit friction history is invalid. Dropping this "
                              "contact's friction state (re-engages fresh) and refusing to arm "
                              "re-emit this step.\n";
                    cd->dropFrictionForContact(ct.tag);
                }
            }
            // The deformed broad-phase feed is active when -reemit is on AND the host is serial
            // (R6). The migration TRIGGER additionally needs a STABLE master membership this
            // handle (R1) — under a membership change we keep the deformed feed (no pass-through)
            // but do not arm an automatic re-handle off the just-changed mesh.
            bool reemitActive = ct.enableReemit && !hostPartitioned;
            bool reemitArmed  = reemitActive && !membershipChanged;

            // P2.5 BROAD PHASE: bucket-sort the master segments (reference coords)
            // once, then per slave query only the 27-neighbour candidate segments
            // instead of all nSeg (brute force). The kept set is a SUPERSET of every
            // near pair, so the contact result is identical to brute force; a huge
            // -cell (=> 1 bucket) reproduces brute force exactly (the equivalence
            // gate). Missing master nodes (malformed model) are filled with the
            // segment's first valid coord — superset-preserving, and the narrow loop
            // below re-fetches + skips them anyway.
            // segCoords is filled with REFERENCE coords first (the shipped P0 path). When re-emit is
            // on it is shifted to the DEFORMED config IN PLACE below — but the R2 search band is taken
            // from these REFERENCE coords (deformation-invariant; it must not drift vs the grid the
            // next re-emit builds — gate D6 + the referenceBand() contract).
            std::vector<double> segCoords((size_t)nSeg * nps * 3, 0.0);
            for (int seg = 0; seg < nSeg; seg++) {
                double *S = &segCoords[(size_t)seg * nps * 3];
                bool haveFirst = false; double first[3] = {0.0, 0.0, 0.0};
                for (int k = 0; k < nps; k++) {
                    Node *mn = theDomain->getNode(mTags(seg * nps + k));
                    if (mn != 0) {
                        const Vector &X = mn->getCrds();
                        S[k*3+0] = X(0); S[k*3+1] = X(1); S[k*3+2] = X(2);
                        if (!haveFirst) { first[0]=X(0); first[1]=X(1); first[2]=X(2); haveFirst = true; }
                    } else {
                        S[k*3+0] = S[k*3+1] = S[k*3+2] = HUGE_VAL;   // mark missing
                    }
                }
                for (int k = 0; k < nps; k++)   // backfill missing with first valid
                    if (S[k*3+0] == HUGE_VAL)
                        for (int d = 0; d < 3; d++) S[k*3+d] = first[d];
            }

            // ADR-60 R2: the deformation-invariant re-emit search band = reference median segment
            // diagonal × cellFrac, from the REFERENCE segCoords ABOVE (before the deformed shift).
            double reemitBand = reemitActive
                ? LadrunoContactReemit::referenceBand(nSeg, nps, segCoords.data(), ct.cellFrac) : 0.0;

            // ADR-60 P1: when re-emit is on, shift segCoords to the DEFORMED (committed) config IN
            // PLACE so the rebuilt candidate set tracks the slid geometry. OFF ⇒ segCoords stays the
            // reference value (the verbatim shipped path) ⇒ byte-identical. A missing-node slot keeps
            // its backfilled reference value (the narrow loop skips that pair anyway).
            if (reemitActive) {
                for (int seg = 0; seg < nSeg; seg++)
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        if (mn == 0) continue;
                        const Vector &u = mn->getDisp();
                        double *S = &segCoords[((size_t)seg * nps + k) * 3];
                        S[0] += u(0); S[1] += u(1); S[2] += u(2);
                    }
            }

            // ADR-60 P1 (BLOCKER-CLIP): the runaway percentile clip targets a STATIC reference
            // outlier; on the deformed feed a legitimately-diverging node would collapse the grid and
            // silently drop real pairs. Disable it (clipPct=0) when re-emit is on — the
            // NX*NY*NZ≤min(nSeg,5000) cap still bounds memory. OFF ⇒ the shipped 1.0 clip (identical).
            // ADR-60 R8: Grid::runawayGuardFired() is intentionally NOT auto-warned. As implemented the
            // guard fires whenever the 1/99-percentile clamp moves a bound — i.e. for ANY mesh with >100
            // segment-centroids (the tails are clipped by design) — so it is not a clean "a node ran away"
            // signal and surfacing it would spam normal large models. The safety Finding-3 cared about (the
            // grid never collapsing during instability) is provided by clipPct=0 here, not by a warning;
            // the accessor stays available for debugging/tests.
            LadrunoContactBucketSort::Grid grid(nSeg, nps, segCoords.data(), ct.cellFrac,
                                                reemitActive ? 0.0 : 1.0);
            std::vector<int> cand(nSeg);

            // ADR-60: register this contact for finite-sliding re-emit (opt-in via -reemit) with the
            // R2 deformation-invariant reference band. Gated on reemitArmed: serial host (R6) AND a
            // stable master membership this handle (R1). OFF ⇒ nothing registered ⇒ byte-identical.
            if (reemitArmed)
                cd->setReemitContact(ct.tag, reemitBand, ct.resortFrac, ct.resortEvery);

            // ADR-63 #4a: build the smooth nodal-normal field for this master surface (opt-in
            // -smoothNormal). Once per handle, from the DEFORMED config (frozen within the step),
            // with the coherent winding cached per the R1 fingerprint. The GLOBAL outward sign is a
            // SURFACE datum captured here ONCE — explicit -outward if given, else the slave-surface
            // REFERENCE centroid minus the master centroid (NEVER a per-slave datum: that was the R3
            // bug). On a refused master (non-manifold / disconnected / non-orientable) a one-time
            // warning fires and smoothFieldOK stays false ⇒ the adapters keep the faceted path.
            // OFF (default) ⇒ nothing built ⇒ byte-identical (BLOCKER-IDENTITY).
            bool smoothFieldOK = false;
            double smoothSeed[3] = {0.0, 0.0, 0.0};   // the GLOBAL outward datum (also the fallback orientDir)
            if (ct.smoothNormal) {
                std::vector<double> segDef((size_t)nSeg * nps * 3, 0.0);
                std::vector<double> segRef((size_t)nSeg * nps * 3, 0.0);   // REFERENCE master (vote basis, F1)
                bool missingNode = false;
                for (int seg = 0; seg < nSeg; seg++)
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        double *S = &segDef[((size_t)seg * nps + k) * 3];
                        double *R = &segRef[((size_t)seg * nps + k) * 3];
                        if (mn != 0) {
                            const Vector &X = mn->getCrds();
                            const Vector &u = mn->getDisp();   // committed/deformed (frozen within step)
                            S[0] = X(0) + u(0); S[1] = X(1) + u(1); S[2] = X(2) + u(2);
                            R[0] = X(0); R[1] = X(1); R[2] = X(2);   // reference (config-independent vote)
                        } else {
                            missingNode = true;   // a malformed master ⇒ refuse the field (review GAP-4):
                                                  // backfilling a missing node makes a garbage facet normal
                                                  // that pollutes the SHARED nodal normals (no per-pair skip
                                                  // here, unlike the bucket grid) → fall back to faceted.
                        }
                    }
                // GLOBAL outward seed (captured once, NOT per-slave): -outward, else slave-ref-centroid
                // − master-ref-centroid over UNIQUE master nodes (the flat mTags double-counts shared
                // nodes, biasing the centroid toward high-valence ridges — review/bring-up bug).
                if (ct.hasOutward) {
                    for (int d = 0; d < 3; d++) smoothSeed[d] = ct.outward[d];
                } else {
                    double mc[3] = {0,0,0}; std::set<int> mseen;
                    for (int i = 0; i < mTags.Size(); i++) {
                        int tg = mTags(i);
                        if (mseen.count(tg)) continue;
                        Node *mn = theDomain->getNode(tg);
                        if (mn != 0) { const Vector &X = mn->getCrds(); for (int d=0;d<3;d++) mc[d]+=X(d); mseen.insert(tg); }
                    }
                    int mcn = (int)mseen.size();
                    double sc[3] = {0,0,0}; int scn = 0;
                    for (int i = 0; i < sTags.Size(); i++) {
                        Node *snd = theDomain->getNode(sTags(i));
                        if (snd != 0) { const Vector &X = snd->getCrds(); for (int d=0;d<3;d++) sc[d]+=X(d); scn++; }
                    }
                    if (mcn > 0 && scn > 0)
                        for (int d = 0; d < 3; d++) smoothSeed[d] = sc[d]/scn - mc[d]/mcn;
                }
                std::vector<int> mt((size_t)mTags.Size());
                for (int i = 0; i < mTags.Size(); i++) mt[i] = mTags(i);
                // ADR-63 auto-sign robustness: on the AUTO path (no -outward) gather the slave
                // REFERENCE coords (config-independent, captured once — D4) for the per-slave majority
                // vote. With -outward we skip them: the sign comes from the explicit aggregate seed.
                std::vector<double> slaveRef;
                if (!ct.hasOutward) {
                    slaveRef.reserve((size_t)sTags.Size() * 3);
                    for (int i = 0; i < sTags.Size(); i++) {
                        Node *snd = theDomain->getNode(sTags(i));
                        if (snd != 0) {
                            const Vector &X = snd->getCrds();
                            slaveRef.push_back(X(0)); slaveRef.push_back(X(1)); slaveRef.push_back(X(2));
                        }
                    }
                }
                int nSlaveVote = (int)slaveRef.size() / 3;
                // contact-review P5: the field is keyed by CONTACT tag (per-contact frozen sign /
                // -outward vote — a second contact sharing the master must not inherit contact A's).
                int st = missingNode ? LadrunoContactNormalField::NON_ORIENTABLE
                       : cd->setNormalField(ct.tag, nps, mt.empty() ? 0 : &mt[0], nSeg,
                                            segDef.data(), smoothSeed,
                                            slaveRef.empty() ? 0 : &slaveRef[0], nSlaveVote,
                                            ct.hasOutward, segRef.data());
                if (st != LadrunoContactNormalField::OK) {
                    if (cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_SMOOTH_REFUSE)) {
                        opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                               << ": -smoothNormal master surface " << ct.masterSurfTag
                               << " could not be coherently oriented (non-manifold / disconnected "
                                  "multi-shell / non-orientable / closed / missing node); falling back to "
                                  "the faceted normal — NOTE the R3 ridge-flip protection is NOT active on "
                                  "this surface, prefer an explicit -outward.\n";
                    }
                } else {
                    smoothFieldOK = true;
                    // ADR-63 review F2/F3/F5: warn once if the FROZEN auto sign is a near coin-flip
                    // (seed ~⟂ the field — a two-sided / saddle / straddle master). Converts a silent
                    // wrong-side field into a loud, actionable signal. Skipped when -outward is explicit.
                    if (!ct.hasOutward) {
                        double conf = cd->getNormalFieldSignConf(ct.tag);
                        if (conf >= 0.0 && conf < 0.1 &&
                            cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_SMOOTH_SIGN)) {
                            opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                                   << ": -smoothNormal auto outward-sign is ILL-CONDITIONED (confidence "
                                   << conf << " ≪ 1 — the slave cloud is nearly tangent to the master, "
                                      "e.g. a two-sided / saddle master). The global sign may be wrong; "
                                      "pass an explicit -outward.\n";
                        }
                    }
                    // ADR-63 (silent-downgrade guard): -smoothNormal + -consistentNormal pre-P3 use the
                    // frozen-field symmetric kn·BᵀB (the faceted B3 ∂n/∂u is suppressed in the adapter);
                    // the smoothed-normal consistent tangent is the gated P3. Surface it once.
                    if (ct.consistentNormal &&
                        cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_SMOOTH_DOWNGRADE)) {
                        opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                               << ": -smoothNormal + -geomtan — the smoothed-normal consistent "
                                  "tangent (∂n_smooth/∂u) is not yet available (ADR-63 P3); using "
                                  "the frozen-field symmetric kn·BᵀB tangent.\n";
                    }
                }
            }

            for (int si = 0; si < sTags.Size(); si++) {
                Node *sn = theDomain->getNode(sTags(si));
                if (sn == 0 || sn->getNumberDOF() != 3) {
                    // ADR-78 P1: dropping a slave node here leaves a hole in an
                    // otherwise-live interface -- the same silently-partial contact
                    // the surface pre-flight exists to prevent, one level down.
                    // (sn == 0 is already impossible past ladrunoSurfaceNodesOk; the
                    // ndf mismatch is the live case.)
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << " slave node " << sTags(si) << " ndf="
                           << (sn != 0 ? sn->getNumberDOF() : 0)
                           << " != 3; ABORTING (P2b is 3D translational -- ADR-78 P1)\n";
                    return ladrunoContactFatal();
                }
                const Vector &Xs0 = sn->getCrds();
                double slavePt[3] = { Xs0(0), Xs0(1), Xs0(2) };
                if (reemitActive) {
                    // ADR-60 P1: query the grid at the slave's DEFORMED (committed) position so a slid
                    // slave finds the segment it is ACTUALLY over (the no-pass-through fix); the SAME
                    // deformed point is the re-emit anchor (drift resets each re-emit ⇒ sane cadence).
                    // OFF ⇒ slavePt stays Xs0 ⇒ byte-identical.
                    const Vector &us = sn->getDisp();
                    slavePt[0] = Xs0(0) + us(0); slavePt[1] = Xs0(1) + us(1); slavePt[2] = Xs0(2) + us(2);
                    // R1: register the anchor only when the trigger is armed (stable membership);
                    // on a membership-change step we keep the deformed query but do not arm a re-sort.
                    if (reemitArmed)
                        cd->addReemitAnchor(ct.tag, sTags(si), slavePt);
                }
                int nCand = grid.candidates(slavePt, cand.data());
                for (int ci = 0; ci < nCand; ci++) {
                    int seg = cand[ci];
                    Node *segNodes[4]; int badTag = 0, badNdf = 0;
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) {
                            badTag = mTags(seg * nps + k);
                            badNdf = (mn != 0) ? mn->getNumberDOF() : -1;
                            break;
                        }
                        segNodes[k] = mn;
                    }
                    if (badTag != 0) {
                        // ADR-78 P1 (added on review): this was `if (!ok) continue;` with NO
                        // message at all -- the last fully SILENT drop in the handler, and on
                        // the master side, which ladrunoSurfaceNodesOk cannot catch because it
                        // checks existence and 3-D coords but NOT ndf. A master node shared
                        // with a beam/shell (pile cap, slab-on-soil) passes pre-flight and then
                        // every pair touching its segments vanished without a word: exactly the
                        // silently-partial interface P1 exists to remove.
                        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                               << " master segment " << seg << " node " << badTag << ": ndf="
                               << badNdf << " != 3 (or node absent); ABORTING (P2b master "
                                  "segments are 3D translational -- ADR-78 P1)\n";
                        return ladrunoContactFatal();
                    }
                    // ADR-60 P1 (BLOCKER-SELF-CONTACT): finite-sliding re-emit can pair a slave with a
                    // segment that CONTAINS it (a folded declared surface). Skip a candidate whose node
                    // set includes the slave node — a self-penetration guard. Re-emit-only ⇒ identical.
                    if (reemitActive) {
                        bool selfSeg = false;
                        for (int k = 0; k < nps; k++)
                            if (mTags(seg * nps + k) == sTags(si)) { selfSeg = true; break; }
                        if (selfSeg) continue;
                    }
                    // orientation direction toward the slave's allowed half-space:
                    // explicit -outward if given, else auto = (slave ref) − (segment
                    // ref centroid). The derived normal is flipped to n·orientDir>0
                    // (winding-immune). Use -outward for just-penetrated starts.
                    double orientDir[3];
                    if (ct.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = ct.outward[d];
                    } else {
                        double cen[3] = {0.0, 0.0, 0.0};
                        for (int k = 0; k < nps; k++) {
                            const Vector &Xk = segNodes[k]->getCrds();
                            for (int d = 0; d < 3; d++) cen[d] += Xk(d);
                        }
                        const Vector &Xs = sn->getCrds();
                        for (int d = 0; d < 3; d++) orientDir[d] = Xs(d) - cen[d] / nps;
                    }
                    // ADR-63 review GAP-2: for a -smoothNormal contact the orientDir is the smoothed
                    // path's degenerate-blend FALLBACK. The per-pair auto direction above is the
                    // R3-buggy datum (it flips on a ridge), so a degenerate query would silently fall
                    // back into R3. Override it with the GLOBAL outward seed (the same surface datum the
                    // field's frozen sign uses) ⇒ the fallback is also global-sign-correct, no R3 hole.
                    if (smoothFieldOK)
                        for (int d = 0; d < 3; d++) orientDir[d] = smoothSeed[d];
                    // P2b-2b: resolve `-kn auto` per (slave, segment) pair from the
                    // owning solid element's stiffness; a fixed `-kn $val` rides ct.kn.
                    double knUse = ct.kn;
                    if (ct.knAuto) {
                        knUse = ladrunoResolveAutoKn(theDomain, segNodes, nps, orientDir);
                        if (knUse <= 0.0) {
                            // ADR-78 P1: was a per-pair skip with a warning, which is the
                            // subtlest of the silent degradations -- the contact still
                            // exists, just with a hole in it, so nothing downstream looks
                            // wrong. It is also the exact failure a partitioned run hits
                            // when the owner rank holds the master NODES but not the
                            // backing solid ELEMENT (apeGmsh ADR 0092 INV-1/INV-4): the
                            // abort is what makes that mis-pick loud instead of silent.
                            opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                   << ct.tag << " slave node " << sTags(si) << " segment "
                                   << seg << ": -kn auto could not size a penalty (no owning "
                                      "solid element / non-3-DOF element / ambiguous normal); "
                                      "ABORTING (ADR-78 P1). Give an explicit -kn, or ensure "
                                      "the master segment's solid element is on this rank.\n";
                            return ladrunoContactFatal();
                        }
                    }
                    // P3: thread kt/mu + the engine ptr + the friction-state key
                    // (contactTag, GLOBAL segment ordinal `seg`). mu<=0 ⇒ the adapter
                    // short-circuits friction (byte-identical to frictionless P2b).
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sn, segNodes, nps, knUse, orientDir,
                                             ct.kt, ct.mu, theDomain, ct.tag, seg,
                                             ct.consistentTan, ct.muc,       // D2 viscous (0 ⇒ off)
                                             ct.consistentNormal,            // B3 ∂n/∂u geom tangent
                                             ct.softScale);                  // B1 SOFT=1 (0 ⇒ off)
                    if (fe == 0) return -5;
                    // ADR-63 #4a: install this segment's FROZEN nodal normals so segmentActive() uses
                    // the smooth field via evalSegmentSmooth (null ⇒ refused/out-of-range ⇒ the adapter
                    // keeps the faceted path). orientDir remains the per-query degenerate-blend fallback.
                    if (smoothFieldOK) {
                        const double *nn = cd->getSegNodalNorm(ct.tag, seg);
                        // ADR-63 P2.1: install the segment's shared-edge mask so evalSegmentSmooth
                        // rejects a projection landing on a SHARED interior edge (facet ownership).
                        const int *se = cd->getSegSharedEdge(ct.tag, seg);
                        if (nn != 0) fe->setSmoothNormals(nn, se);
                    }
                    theModel->addFE_Element(fe);
                    if (ct.mu > 0.0)
                        cd->frictionGCMark(ct.tag, sTags(si), seg);  // live this handle()
                }
            }
        }

        cd->frictionGCEnd();   // prune friction slots no live pair referenced
    }

    // --- ADR-41 C2.1: MORTAR pairing -> ONE clipped-GP penalty adapter per (slave facet,
    //     candidate master facet). Mirrors the NTS loop but FACET-vs-FACET: broad-phase
    //     the master facets, query per slave-facet CENTROID, and let the kernel's overlap
    //     clip reject non-overlapping candidates (status<0 ⇒ the adapter is inert that
    //     eval). Frictionless penalty (λ_N Uzawa = C2.2; friction = C3). No path state. ---
    if (cd != 0) {
        cd->mortarNormalGCBegin();   // C2.2: rebuild the live λ_N node-set this handle()
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            if (mc.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": master/slave surface not defined; ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": need a MASTER_SEGMENTS master + SLAVE_SEGMENTS slave; "
                          "ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            // ADR-85 T0: ONE lane name per iteration, shared by the early block and the
            // pre-flight below (both used to build the same ternary separately).
            const char *laneMort = mc.isTie ? "mortar mesh-tie (-tie)"
                                            : "mortar segment-to-segment";
            // ADR-85 T3 -- the pre-flight now runs in LIVE mode (lane2DLive=true,
            // &pairDim captured): a well-formed 2D pair PROCEEDS into the live 2D
            // mortar lane below instead of drawing the T0 "lands in T3" refusal --
            // exactly the T1b NTS collapse, generalized to the mortar loop. An all-3D
            // pair passes silently (pairDim is a write-only out) and takes the shipped
            // path below, textually untouched.
            int pairDim = 3;
            if (npsM == 2 || npsS == 2) {
                if (!ladrunoContactPairPreflight(theDomain, mc.tag, ms, ss, laneMort, "T3",
                                                 /*lane2DLive=*/true, &pairDim))
                    return ladrunoContactFatal();
            }
            if (pairDim == 2) {
                // ==============================================================
                // ADR-85 T3 -- the LIVE 2D mortar lane (interval clip on the
                // slave-trace measure, ALM, tie, mortar friction). A 3D deck
                // CANNOT enter: pairDim == 2 requires a surface whose every node
                // carries 2-component coordinates (pre-flight above). Everything
                // in this block ends in `continue`, so the shipped 3D mortar
                // path below never sees a 2D node and stays textually untouched.
                // ==============================================================
                if (npsM != 2 || npsS != 2) {
                    // only reachable through a stream restore (the parser pins
                    // arity <-> dimension at declaration); refuse by name.
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": a 2D mortar contact needs nodesPerSeg == 2 on BOTH surfaces "
                              "(got master " << npsM << ", slave " << npsS
                           << "); ABORTING (ADR-85 nodesPerSeg)\n";
                    return ladrunoContactFatal();
                }
                if (mc.softScale > 0.0) {
                    // decision 2b.6 -- SOFT=2 segment-based explicit penalty is a
                    // deferral, not a silent drop: it would ride the T3 mortar
                    // interval exactly as B2 rides the 3D clip, but no 2D demand
                    // has driven it yet (ADR-85 What / not-in-scope list).
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": SOFT=2 segment-based explicit penalty is deferred in 2D "
                              "(ADR-85 deferral row) -- remove -soft or use the 3D lane. "
                              "ABORTING\n";
                    return ladrunoContactFatal();
                }
                if (mc.ngp != 2 &&
                    cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_2D_NGP_NOOP)) {
                    opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": -ngp is ACCEPTED AND IGNORED on the 2D mortar lane -- 2-point "
                              "Gauss is EXACT for the straight-segment slave-trace integrands "
                              "(ADR-85 SS How/3, oracle-proven), so the 2D kernel always uses it.\n";
                }
                const ID &mTags2 = ms->getNodeTags();
                const ID &sTags2 = ss->getNodeTags();
                int nSegM = mTags2.Size() / 2;
                int nSegS = sTags2.Size() / 2;

                // ---- Lref: the min INITIAL segment length over BOTH surfaces --
                //      (mortar integrates BOTH -- the T1b NTS Lref is master-only
                //      because NTS projects onto the master alone; the 2D mortar
                //      interval clip projects the SLAVE endpoints too, so its
                //      relative gauges need the smaller of the two surfaces).
                double Lref = 0.0;
                int LrefBadSurf = -1, LrefBadSeg = -1;
                for (int seg = 0; seg < nSegM; seg++) {
                    const Vector &A = theDomain->getNode(mTags2(2 * seg))->getCrds();
                    const Vector &Bc = theDomain->getNode(mTags2(2 * seg + 1))->getCrds();
                    double dx = Bc(0) - A(0), dy = Bc(1) - A(1);
                    double L = std::sqrt(dx * dx + dy * dy);
                    if (!(L > 0.0)) { LrefBadSurf = 0; LrefBadSeg = seg; break; }
                    if (Lref <= 0.0 || L < Lref) Lref = L;
                }
                for (int seg = 0; LrefBadSeg < 0 && seg < nSegS; seg++) {
                    const Vector &A = theDomain->getNode(sTags2(2 * seg))->getCrds();
                    const Vector &Bc = theDomain->getNode(sTags2(2 * seg + 1))->getCrds();
                    double dx = Bc(0) - A(0), dy = Bc(1) - A(1);
                    double L = std::sqrt(dx * dx + dy * dy);
                    if (!(L > 0.0)) { LrefBadSurf = 1; LrefBadSeg = seg; break; }
                    if (Lref <= 0.0 || L < Lref) Lref = L;
                }
                if (nSegM <= 0 || nSegS <= 0 || LrefBadSeg >= 0 || !(Lref > 0.0)) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": ";
                    if (nSegM <= 0 || nSegS <= 0)
                        opserr << "a 2D mortar surface has no segments";
                    else
                        opserr << (LrefBadSurf == 0 ? "master" : "slave") << " segment "
                               << LrefBadSeg << " has zero reference length (the kernel's "
                                  "relative-gauge caller contract Lref would be disabled)";
                    opserr << "; ABORTING (ADR-85 T3)\n";
                    return ladrunoContactFatal();
                }

                // ---- the SAME interface-level orientation vote as the 2D NTS lane
                //      (ADR-85 T3 -- pure code motion, shared via
                //      ladruno2DOrientationVote; see that helper's block comment).
                double sigma = 0.0;
                if (!ladruno2DOrientationVote(theDomain, mc.tag, mc.hasOutward, mc.outward,
                                              mTags2, nSegM, sTags2, Lref, sigma))
                    return ladrunoContactFatal();

                // ---- ntsCoDeclared union: an NTS contact declared on the SAME
                //      interface (How/3 decision 2b.4). CONTRACT CORRECTION
                //      (coordinator, supersedes the brief's "exact slave surface
                //      tag" rule -- that rule is unimplementable: NTS slave
                //      surfaces are SLAVE_NODES node sets while the mortar lane
                //      requires SLAVE_SEGMENTS, so the two lanes necessarily use
                //      DIFFERENT slave surface tags even on the identical physical
                //      interface, e.g. G-T3(d)'s own gate deck: NTS -slave 1 is
                //      surface 20, mortar -slave-segments 2 1 2 is surface 21).
                //      Binding rule: union the slave node TAGS of every NTS
                //      contact whose MASTER surface tag equals this mortar pair's
                //      master tag; a given mortar SLAVE SEGMENT is "covered" iff
                //      BOTH its node tags are in that union (tested per segment
                //      at injection time below, not once per MortarContact -- a
                //      mortar surface can have segments only partially shadowed
                //      by the NTS slave set).
                std::set<int> ntsSlaveUnion;
                for (int ci = 0; ci < cd->getNumContacts(); ci++) {
                    const LadrunoContactDomain::Contact &ntsCt = cd->getContact(ci);
                    if (ntsCt.retired) continue;
                    if (ntsCt.masterSurfTag != mc.masterSurfTag) continue;
                    LadrunoContactSurface *ntsSlave = cd->getSurface(ntsCt.slaveSurfTag);
                    if (ntsSlave == 0) continue;
                    const ID &ntsSlaveTags = ntsSlave->getNodeTags();
                    for (int ti = 0; ti < ntsSlaveTags.Size(); ti++)
                        ntsSlaveUnion.insert(ntsSlaveTags(ti));
                }

                // penalty: epsN if given, else the kn slot; auto-size per pair if requested
                // (the SAME rule as the 3D lane and the T1b NTS -kn resolution).
                bool epsAuto = mc.epsNAuto || mc.knAuto;
                double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
                if (!epsAuto && epsFixed <= 0.0) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": needs a penalty (-epsN val|auto or kn > 0); "
                              "ABORTING (ADR-78 P1)\n";
                    return ladrunoContactFatal();
                }

                // ---- brute-force double loop (slave segment, master segment) -----
                //      Mirrors the 3D mortar loop's reasoning: EVERY master segment
                //      is a candidate per slave segment (no bucket sort, no pre-cull
                //      -- the 3D loop injects one FE per (sf,seg) pair unconditionally
                //      too and lets the kernel's own overlap clip decide EMPTY/DEGEN/
                //      ALIGN lazily at each getResidual; mirrored faithfully here).
                // REVIEW FIX (T3, efficiency): the 2D auto-penalty resolver's inputs
                // (mNodes, sigma) are INVARIANT across slave segments -- unlike the 3D
                // twin, whose orientDir varies per slave facet -- so it is resolved once
                // per MASTER segment and memoized (-1.0 = not yet resolved).
                std::vector<double> epsAutoPerSeg;
                if (epsAuto) epsAutoPerSeg.assign(nSegM, -1.0);
                for (int sf = 0; sf < nSegS; sf++) {
                    Node *sNodes[2];
                    for (int k = 0; k < 2; k++) {
                        Node *sn = theDomain->getNode(sTags2(sf * 2 + k));
                        if (sn == 0 || sn->getNumberDOF() != 2) {
                            opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave node " << sTags2(sf * 2 + k) << " ndf="
                                   << (sn != 0 ? sn->getNumberDOF() : -1)
                                   << " != 2; ABORTING (2D mortar is 2D translational -- "
                                      "ADR-85 ndf equality contract)\n";
                            return ladrunoContactFatal();
                        }
                        sNodes[k] = sn;
                    }
                    // per-SLAVE-SEGMENT coverage test (the coordinator's binding rule):
                    // this segment's ALIGN refusal is silent iff BOTH its node tags are
                    // in the NTS-slave union -- a mortar surface can have segments only
                    // partially shadowed by the co-declared NTS slave set, so this is
                    // computed per segment, not once per MortarContact.
                    bool ntsCoDeclared = (ntsSlaveUnion.count(sNodes[0]->getTag()) > 0 &&
                                         ntsSlaveUnion.count(sNodes[1]->getTag()) > 0);
                    for (int seg = 0; seg < nSegM; seg++) {
                        Node *mNodes[2];
                        for (int k = 0; k < 2; k++) {
                            Node *mn = theDomain->getNode(mTags2(seg * 2 + k));
                            if (mn == 0 || mn->getNumberDOF() != 2) {
                                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                       << mc.tag << " master segment " << seg << " node "
                                       << mTags2(seg * 2 + k) << ": ndf="
                                       << (mn != 0 ? mn->getNumberDOF() : -1)
                                       << " != 2 (or node absent); ABORTING (2D mortar master "
                                          "segments are 2D translational -- ADR-85 ndf equality "
                                          "contract)\n";
                                return ladrunoContactFatal();
                            }
                            mNodes[k] = mn;
                        }
                        double epsUse = epsFixed;
                        if (epsAuto) {
                            if (epsAutoPerSeg[seg] < 0.0)
                                epsAutoPerSeg[seg] =
                                    ladrunoResolveAutoKn2D(theDomain, mNodes, sigma);
                            epsUse = epsAutoPerSeg[seg];
                            if (epsUse <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                       << mc.tag << " slave segment " << sf << " master segment "
                                       << seg << ": auto penalty could not be sized (no owning "
                                          "2-DOF/node solid element / degenerate segment); "
                                          "ABORTING (ADR-78 P1). Give an explicit -epsN.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        // C3.1 friction: epsT auto => size from the normal penalty; else the
                        // given value (mirrors the 3D mortar lane, MINOR-1 default-and-warn).
                        bool epsTFromEpsN = mc.epsTAuto;   // provenance for the h-scaling below
                        double epsTuse = mc.epsTAuto ? epsUse : mc.epsT;
                        bool wantFric = (mc.mu > 0.0 || mc.cohesion > 0.0 || mc.tauMax > 0.0);
                        if (wantFric && mc.consistentTan &&
                            cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_CSL)) {
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": -consistanttan (non-symmetric Coulomb friction "
                                      "tangent) needs a non-symmetric solver (system FullGeneral/"
                                      "UmfPack/BandGeneral); symmetric solvers will silently "
                                      "corrupt it.\n";
                        }
                        if (wantFric && epsTuse <= 0.0) {
                            epsTuse = epsUse;
                            epsTFromEpsN = true;
                            if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EPST)) {
                                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                       << mc.tag << ": friction requested without -epsT; "
                                          "defaulting the tangential penalty to the normal penalty "
                                          "(-epsT auto). Set -epsT explicitly to control it.\n";
                            }
                        }
                        // ---- -thickness h applied ONCE, HERE, at the 2D injection site
                        //      (ADR-85 SS How/7): epsN, epsT, muc, cohesion, tauMax all
                        //      h-scale; lambda_N/lambda_T/lambda_Tie are NEVER touched
                        //      here -- they inherit h through epsN*h*gbar in the Uzawa
                        //      update (scaling them again is the gated h^2 error).
                        // REVIEW FIX (T3): the AUTO-resolved penalty comes from the owning
                        // element's getInitialStiff(), which ALREADY folds the element's
                        // real thickness into K (the SS How/7 NTS argument: auto absorbs
                        // h automatically) -- h-scaling it again is itself an h^2 error.
                        // Only EXPLICIT per-unit-thickness user inputs h-scale; a
                        // defaulted epsT inherits the provenance of the value it
                        // derives from (epsTFromEpsN above).
                        double h = mc.hThickness;
                        double epsUseH  = epsAuto ? epsUse : epsUse * h;
                        double epsTuseH = epsTFromEpsN ? epsUseH : epsTuse * h;
                        double mucH      = mc.muc      * h;
                        double cohesionH = mc.cohesion  * h;
                        double tauMaxH   = mc.tauMax    * h;

                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sNodes, /*nps_s=*/2, mNodes, /*nps_m=*/2, /*ndm2=*/2,
                            epsUseH, sigma, Lref, mc.tag, sf, theDomain,
                            mc.mu, epsTuseH, cohesionH, tauMaxH, mc.consistentTan, mc.isTie,
                            mucH, ntsCoDeclared);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        // C2.2 twin: this pair's slave nodes have a live lambda_N slot this
                        // handle() (mortarNormalGCMark also covers the friction slot -- they
                        // share the SAME per-(contactTag,slaveNodeTag) MortarNormalState
                        // record, exactly as the 3D mortar loop relies on).
                        cd->mortarNormalGCMark(mc.tag, sNodes[0]->getTag());
                        cd->mortarNormalGCMark(mc.tag, sNodes[1]->getTag());
                    }
                }
                continue;   // 2D mortar contact done -- the shipped 3D code below never sees it
            }
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": nodesPerSeg must be 3 or 4 (slave " << npsS << ", master " << npsM
                       << "); ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            if (mc.hThickness != 1.0) {
                // ADR-85 T3 -- -thickness is a 2D plane-model option; a 3D mortar deck's
                // thickness lives in its elements (SS How/7). Only reachable here (a 3D
                // pair) via a direct API / stream-restore call with hThickness != 1 --
                // the command surface pairs -thickness with -mortar, not with dimension,
                // so this is the dimension-aware backstop.
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": -thickness is a 2D plane-model option; a 3D mortar deck's "
                          "thickness lives in its elements. ABORTING (ADR-85 T3)\n";
                return ladrunoContactFatal();
            }
            // penalty: epsN if given, else the kn slot; auto-size per pair if requested.
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            if (!epsAuto && epsFixed <= 0.0) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": needs a penalty (-epsN val|auto or kn > 0); "
                          "ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            // ADR-78 P1 (missing / bad-dimension node) + ADR-85 T0/T3 pair gate,
            // mirroring the NTS site. BOTH mortar surfaces are faceted (MASTER_SEGMENTS +
            // SLAVE_SEGMENTS), so both carry an arity and the shipped 3/4 check above
            // already implies 3D nodes by the parser's arity<->dimension pairing; this
            // is the defense-in-depth twin that does not depend on that inference
            // holding (e.g. a surface rebuilt by unpackDefinitions, which never sees the
            // parser). Flipped to lane2DLive (ADR-85 T3) for uniformity with the earlier
            // call; reaching here with pairDim == 2 would mean a corrupted stream (arity
            // 3/4 declared over 2D node coordinates). REVIEW-CORRECTED rationale: only
            // the SLAVE-node ndf==3 check below is a named FATAL (the master-node check
            // is a silent per-facet skip) -- the invariant holds because the slave loop
            // runs FIRST, a loop-ordering property, not a per-node gate. If that
            // ordering ever changes, this corrupted-stream case needs its own branch.
            if (!ladrunoContactPairPreflight(theDomain, mc.tag, ms, ss, laneMort, "T3",
                                             /*lane2DLive=*/true, &pairDim))
                return ladrunoContactFatal();
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;

            // C2.1 broad phase = BRUTE FORCE (every master facet a candidate per slave
            // facet); the kernel's overlap clip rejects non-overlaps (status<0 ⇒ inert).
            // The NTS bucket-sort superset contract is proven for POINT slaves (slave
            // NODES) only — a CENTROID query for an EXTENDED slave facet silently drops
            // overlaps when the slave facet is coarser than the master bucket (gate MAJOR,
            // C2.1 review: ~44% area lost on a coarse-slave/fine-master pairing). A
            // slave-aware mortar broad phase (corner-span query / max-diag cell) is a
            // follow-up ("mortar P2.5"); correctness-first brute force here, like NTS P2b-1.
            for (int sf = 0; sf < nSegS; sf++) {
                Node *sNodes[4]; bool ok = true; double scen[3] = {0.0, 0.0, 0.0};
                double Sx[4][3];                                  // ref coords (ADR-57 E2 stand-down)
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() != 3) {
                        // ADR-78 P1 (added on review): was warn-and-skip-the-facet, which is
                        // the mortar twin of the NTS slave-ndf abort. Leaving the two lanes
                        // inconsistent meant a mortar interface could still lose facets while
                        // the NTS one refused to.
                        opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                               << mc.tag << " slave node " << sTags(sf * npsS + k) << " ndf="
                               << (sn != 0 ? sn->getNumberDOF() : -1)
                               << " != 3; ABORTING (mortar is 3D translational -- ADR-78 P1)\n";
                        return ladrunoContactFatal();
                    }
                    sNodes[k] = sn;
                    const Vector &Xk = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = Xk(d); scen[d] += Xk(d); }
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) scen[d] /= npsS;     // slave-facet centroid

                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mNodes[4]; ok = true; double mcen[3] = {0.0, 0.0, 0.0};
                    double Mx[4][3];                              // ref coords (ADR-57 E2 stand-down)
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        mNodes[k] = mn;
                        const Vector &Xk = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = Xk(d); mcen[d] += Xk(d); }
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) mcen[d] /= npsM;
                    // ADR-57 E2 — FACE-MORTAR STANDS DOWN when the edge-edge lane owns this pair (the
                    // gate EE-1 fix): in the oblique band (50°–90°) the mortar clip is NON-degenerate,
                    // so without this the clip pressure and the injected edge-edge force would
                    // SUPERPOSE on the same region. The SAME routing decision the edge-edge loop uses
                    // (one source of truth) ⇒ exactly one lane owns the pair. mc.edgeEdge off ⇒ the
                    // shipped mortar behavior, unchanged.
                    if (mc.edgeEdge &&
                        ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx, ladrunoEdgeBand(npsM, Mx, mc.edgeBand)))
                        continue;
                    // orientation toward the slave's allowed half-space: explicit -outward,
                    // else (slave facet centroid − master facet centroid) at the ref config.
                    double orientDir[3];
                    if (mc.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    } else {
                        for (int d = 0; d < 3; d++) orientDir[d] = scen[d] - mcen[d];
                        // contact-review P5 (the PR-3 gate follow-up): the kernel signs the facet
                        // normal by flipping it toward orientDir, so the hazard is the NORMAL
                        // COMPONENT of the auto orientation vanishing (gate lens-B): a COINCIDENT
                        // conforming pair (scen − mcen ≈ 0) AND a zero-gap STAGGERED non-conforming
                        // pair (scen − mcen tangential, O(h) — the canonical mortar interface) both
                        // leave the flip test to roundoff ⇒ the mortar sign is WINDING LUCK (an
                        // unlucky winding makes the contact attractive, with no diagnostic). Warn
                        // once per contact when |n̂_m·orientDir| ≈ 0 at a PROXIMATE pair (|orientDir|
                        // ≤ 2h; far coplanar pairs the clip rejects anyway stay silent). A -tie pair
                        // is exempt (its full 3-vec bond r→0 is sign-independent — and ties are the
                        // COMMON coincident case). Thresholds RELATIVE to the facet size (unit-safe).
                        if (!mc.isTie) {
                            double d2 = orientDir[0]*orientDir[0] + orientDir[1]*orientDir[1]
                                      + orientDir[2]*orientDir[2];
                            double h2 = 0.0;
                            for (int k = 0; k < npsM; k++) {
                                const double *a = Mx[k]; const double *b = Mx[(k + 1) % npsM];
                                double ex = b[0]-a[0], ey = b[1]-a[1], ez = b[2]-a[2];
                                double e2 = ex*ex + ey*ey + ez*ez;
                                if (e2 > h2) h2 = e2;
                            }
                            double nM[3];
                            ladrunoFacetNormalNewell(npsM, Mx, nM);   // unit (zero if degenerate)
                            double odn = nM[0]*orientDir[0] + nM[1]*orientDir[1]
                                       + nM[2]*orientDir[2];
                            if (odn*odn <= 1e-12 * h2 && d2 <= 4.0 * h2 &&
                                cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_MORTAR_ORIENT)) {
                                opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                       << mc.tag << ": auto orientation is DEGENERATE (the slave-to-"
                                          "master centroid direction has no component along the "
                                          "master facet normal — a coincident or zero-gap staggered "
                                          "interface); the contact sign falls to the facet winding "
                                          "and may be attractive. Pass an explicit -outward.\n";
                            }
                        }
                    }

                    double epsUse = epsFixed;
                    if (epsAuto) {
                        epsUse = ladrunoResolveAutoKn(theDomain, mNodes, npsM, orientDir);
                        if (epsUse <= 0.0) {
                            // ADR-78 P1 (added on review): the SAME resolver and the SAME
                            // failure as the NTS `-kn auto` abort, converted there and left
                            // skipping here. D5.2 named this path fatal, and apeGmsh ADR 0092
                            // INV-1 leans on it: the node-majority owner proxy is only safe
                            // because a mis-pick makes auto-sizing fail LOUDLY. Without this,
                            // that backstop did not exist for mortar mesh-ties -- which is
                            // exactly what ADR-78 P4 plans to run partitioned.
                            opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave facet " << sf << " master facet " << seg
                                   << ": auto penalty could not be sized (no owning solid "
                                      "element / non-3-DOF element / ambiguous normal); "
                                      "ABORTING (ADR-78 P1). Give an explicit -epsN, or ensure "
                                      "the master facet's solid element is on this rank.\n";
                            return ladrunoContactFatal();
                        }
                    }
                    // C3.1 friction: epsT auto ⇒ size from the normal penalty (epsUse); else the
                    // given value. mu/cohesion/tauMax ≤0 ⇒ the adapter short-circuits (frictionless).
                    double epsTuse = mc.epsTAuto ? epsUse : mc.epsT;
                    // gate MINOR-1: friction REQUESTED but no tangential penalty given ⇒ the return
                    // map would stick at ZERO force (silently inert). Default epsT to the normal
                    // penalty and warn ONCE rather than ship dead friction.
                    bool wantFric = (mc.mu > 0.0 || mc.cohesion > 0.0 || mc.tauMax > 0.0);
                    // C3.3: -consistanttan opts into the NON-SYMMETRIC consistent Coulomb friction
                    // tangent (Csl), which REQUIRES a non-symmetric solver (system FullGeneral /
                    // UmfPack / BandGeneral); a symmetric SOE silently drops the lower triangle and
                    // corrupts it. The default (symmetric) tangent is correct on any solver. Warn once.
                    if (wantFric && mc.consistentTan &&
                        cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_CSL)) {
                        opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                               << mc.tag << ": -consistanttan (non-symmetric Coulomb friction "
                                  "tangent) needs a non-symmetric solver (system FullGeneral/"
                                  "UmfPack/BandGeneral); symmetric solvers will silently corrupt it.\n";
                    }
                    if (wantFric && epsTuse <= 0.0) {
                        epsTuse = epsUse;
                        if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EPST)) {
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": friction requested without -epsT; defaulting the "
                                      "tangential penalty to the normal penalty (-epsT auto). Set "
                                      "-epsT explicitly to control it.\n";
                        }
                    }
                    // C4 — a mesh-tie pair (mc.isTie): epsUse is epsTie; friction params are 0
                    // (refused with -tie upstream). The adapter runs the full-3-vec r→0 bond.
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sNodes, npsS, mNodes, npsM, epsUse,
                                             orientDir, mc.tag, sf, theDomain,
                                             mc.mu, epsTuse, mc.cohesion, mc.tauMax, mc.consistentTan,
                                             mc.isTie, mc.muc,    // D2.2 mortar viscous (0 ⇒ off)
                                             mc.softScale);       // B2 SOFT=2 segment-based penalty (0 ⇒ off)
                    if (fe == 0) return -5;
                    theModel->addFE_Element(fe);
                    // C2.2: this pair's slave nodes have a live λ_N slot this handle().
                    for (int k = 0; k < npsS; k++)
                        cd->mortarNormalGCMark(mc.tag, sTags(sf * npsS + k));
                }
            }
        }
        cd->mortarNormalGCEnd();   // prune λ_N slots no live mortar pair referenced
    }

    // --- ADR-57 E2: perpendicular edge-edge routing + adapter injection (the cos_t→0 fallback) ---
    // For each -mortar contact with -edgeedge, run the E1 routing detector over the SAME brute-force
    // (slave facet, master facet) enumeration the mortar lane uses (NOT the bucket sort — that is
    // NTS point-vs-segment only) and, on an EDGE_EDGE route, inject one EDGE_EDGE adapter per kept
    // boundary-edge pair. The mortar loop STANDS DOWN on a pair edge-edge owns (the same
    // ladrunoEdgeEdgeOwns decision), so the mortar clip and the edge-edge force never superpose —
    // even in the oblique band where the clip is non-degenerate (the gate EE-1 fix). Reference-config
    // geometry (§7 scoping; large sliding ⇒ re-emit). -edgeedge absent ⇒ zero adapters ⇒ byte-
    // identical; -edgeedge present but no pair routes ⇒ also zero adapters.
    if (cd != 0) {
        cd->edgeGCBegin();
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            if (mc.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
            if (!mc.edgeEdge) continue;                 // opt-in off ⇒ skip (byte-identical)
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) continue;           // already warned by the mortar loop
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) continue;
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            // ADR-85 T3 -- PERMANENT scope fence (was the ADR-85 T0 "lands in T3"
            // not-yet-supported refusal; T3 has now landed and there is no 2D
            // analogue to land -- naming a future phase here would be a lie). Edge-
            // edge is the 3D vertex-face/face-face/edge-edge trichotomy's THIRD
            // primitive; the 2D taxonomy is segments + concave vertices only (ADR-85
            // "What" taxonomy note) -- there is no 2D feature an edge-edge fallback
            // could ever route to. BEFORE the silent `continue` below (the arity
            // gate says nothing; harmless while nps 2 was impossible, but a STREAM-
            // RESTORED surface -- unpackDefinitions never runs the parser -- could
            // arrive here with nps == 2, and standing down silently is exactly the
            // do-nothing-quietly mode these gates exist to prevent). A 3D pair
            // (nps 3 or 4) never enters this block and reaches the `continue` unchanged.
            if (npsM == 2 || npsS == 2) {
                int pairDimEE = 3;
                if (!ladrunoContactPairPreflight(theDomain, mc.tag, ms, ss,
                                                 "mortar edge-edge", "T3",
                                                 /*lane2DLive=*/true, &pairDimEE))
                    return ladrunoContactFatal();
                if (pairDimEE == 2) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": edge-edge contact has no 2D analogue (ADR-85 scope fence) -- "
                              "the 2D narrow-phase is segments + concave vertices; remove "
                              "-edgeedge from this 2D declaration. ABORTING\n";
                    return ladrunoContactFatal();
                }
            }
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) continue;
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            // ADR-78 P1 (missing / bad-dimension node) + the ADR-85 T0/T3 pair gate.
            // NOTE the control flow: this loop's ARITY gate uses `continue`, but the
            // PRE-FLIGHT result has always been FATAL here, and the dimension gate rides
            // the pre-flight, not the arity -- so shipped control flow is preserved on
            // both counts. Flipped to lane2DLive (ADR-85 T3) for uniformity; reaching
            // here with pairDim == 2 would mean a corrupted stream (arity 3/4 declared
            // over 2D node coordinates). REVIEW-CORRECTED rationale (mirrors the main
            // mortar loop's second flipped site): the named FATAL relies on the SLAVE
            // ndf loop running before the master's silent per-facet skip -- a loop-
            // ordering property, not a per-node gate.
            {
                int pairDimEE2 = 3;
                if (!ladrunoContactPairPreflight(theDomain, mc.tag, ms, ss,
                                                 "mortar edge-edge", "T3",
                                                 /*lane2DLive=*/true, &pairDimEE2))
                    return ladrunoContactFatal();
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;
            // de-dup: a facet edge is shared by two facets, so key each injected adapter by the
            // ORDERED (slave edge nodes, master edge nodes) 4-tuple and inject it once (the §7 step-3
            // set, mirroring the bucket-sort stamp idiom — a NEW handler set, not a reuse).
            std::set<std::vector<int> > injectedEdges;
            for (int sf = 0; sf < nSegS; sf++) {
                Node *sN[4]; bool ok = true; double Scen[3] = {0, 0, 0}, Sx[4][3];
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() != 3) { ok = false; break; }
                    sN[k] = sn;
                    const Vector &X = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = X(d); Scen[d] += X(d); }
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) Scen[d] /= npsS;
                double nS[3]; ladrunoFacetNormalNewell(npsS, Sx, nS);
                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mN[4]; ok = true; double Mcen[3] = {0, 0, 0}, Mx[4][3];
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() != 3) { ok = false; break; }
                        mN[k] = mn;
                        const Vector &X = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = X(d); Mcen[d] += X(d); }
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) Mcen[d] /= npsM;

                    // does edge-edge OWN this pair? (the single-source-of-truth routing — the SAME
                    // decision the mortar stand-down uses, so exactly one lane owns the pair). d_band
                    // from -edgeBand, else the facet edge length.
                    double dBand = ladrunoEdgeBand(npsM, Mx, mc.edgeBand);
                    if (!ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx, dBand)) continue;

                    // orientation toward the slave's allowed half-space (the §2 first-capture sign
                    // reference): explicit -outward, else the SLAVE FACET NORMAL (always non-degenerate,
                    // unlike a centroid difference that can vanish near-coincident centroids — the gate
                    // E2-2 fix; ADR §2 "orient once from {n̂_s, n̂_m}") flipped from master toward slave.
                    double orientDir[3];
                    if (mc.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    } else {
                        double cdv[3] = {Scen[0]-Mcen[0], Scen[1]-Mcen[1], Scen[2]-Mcen[2]};
                        double dc = nS[0]*cdv[0] + nS[1]*cdv[1] + nS[2]*cdv[2];
                        double sgn = (dc >= 0.0) ? 1.0 : -1.0;
                        for (int d = 0; d < 3; d++) orientDir[d] = sgn * nS[d];
                    }
                    // edge penalty: explicit -edgeKn, else auto per master facet (if -edgeKn auto or the
                    // mortar penalty is auto), else the resolved fixed mortar penalty (the ADR default).
                    double edgeKnUse = mc.edgeKn;
                    if (edgeKnUse <= 0.0) {
                        if (mc.edgeKnAuto || epsAuto)
                            edgeKnUse = ladrunoResolveAutoKn(theDomain, mN, npsM, orientDir);
                        else
                            edgeKnUse = epsFixed;
                    }
                    if (edgeKnUse <= 0.0) {
                        // ADR-78 P1 (added on review): third lane of the same auto-penalty
                        // failure. Converting NTS and mortar but not edge-edge would leave
                        // the fallback lane -- the one that catches what the face clip
                        // degenerates on -- quietly dropping pairs.
                        opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                               << " slave facet " << sf << " master facet " << seg
                               << ": -edgeedge penalty could not be sized; ABORTING (ADR-78 P1). "
                                  "Give an explicit -edgeKn.\n";
                        return ladrunoContactFatal();
                    }
                    // E3 friction: edgeKt (tangential penalty). Friction requested but no -edgeKt ⇒
                    // the return map would stick at ZERO force (silently inert); default it to the
                    // edge normal penalty and warn ONCE (mirrors the mortar -epsT auto-default).
                    double edgeKtUse = mc.edgeKt;
                    bool edgeFric = (mc.edgeMu > 0.0 || mc.edgeCohesion > 0.0 || mc.edgeTauMax > 0.0);
                    if (edgeFric && edgeKtUse <= 0.0) {
                        edgeKtUse = edgeKnUse;
                        if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EDGEKT)) {
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                                   << ": edge-edge friction requested without -edgeKt; defaulting the "
                                      "tangential penalty to the edge normal penalty. Set -edgeKt to control it.\n";
                        }
                    }

                    // enumerate the (≤4)×(≤4) boundary-edge pairs; inject the margin-interior, well-
                    // conditioned, in-band crossings (the E0 kernel decides interior + non-parallel).
                    for (int ie = 0; ie < npsS; ie++) {
                        Node *sa = sN[ie], *sb = sN[(ie + 1) % npsS];
                        for (int je = 0; je < npsM; je++) {
                            Node *ma = mN[je], *mb = mN[(je + 1) % npsM];
                            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                                Sx[ie], Sx[(ie + 1) % npsS], Mx[je], Mx[(je + 1) % npsM]);
                            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
                            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
                            double g = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                            if (g > dBand) continue;
                            int sat = sa->getTag(), sbt = sb->getTag();
                            int mat = ma->getTag(), mbt = mb->getTag();
                            std::vector<int> key(4);
                            key[0] = (sat < sbt) ? sat : sbt; key[1] = (sat < sbt) ? sbt : sat;
                            key[2] = (mat < mbt) ? mat : mbt; key[3] = (mat < mbt) ? mbt : mat;
                            if (injectedEdges.find(key) != injectedEdges.end()) continue;
                            injectedEdges.insert(key);
                            LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sa, sb, ma, mb,
                                                       edgeKnUse, orientDir, mc.tag, theDomain,
                                                       mc.edgeMu, edgeKtUse, mc.edgeCohesion,
                                                       mc.edgeTauMax, mc.edgeConsistentTan,
                                                       mc.edgeSoftScale,    // E5 explicit SOFT (0 ⇒ off)
                                                       mc.edgeAlm);         // E6 one-scalar ALM (off ⇒ penalty)
                            if (fe == 0) return -5;
                            theModel->addFE_Element(fe);
                            cd->edgeGCMark(mc.tag, sat, sbt, mat, mbt);  // live this handle()
                        }
                    }
                }
            }
        }
        cd->edgeGCEnd();   // prune edge-edge slots no live pair referenced
    }

    // P2a: rigid analytical-plane contacts -> ONE bound adapter per slave node
    // (connectivity = that node). ndm derived per-node from its coordinates.
    if (cd != 0) {
        for (int p = 0; p < cd->getNumRigidPlanes(); p++) {
            const LadrunoContactDomain::RigidPlane &rp = cd->getRigidPlane(p);
            if (rp.retired) continue;   // ADR-78 removal lane: retired by pruneRemovedNodes()
            LadrunoContactSurface *surf = cd->getSurface(rp.slaveSurfTag);
            if (surf == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - contactPlane " << rp.tag
                       << ": slave surface " << rp.slaveSurfTag
                       << " not defined; ABORTING (ADR-78 P1)\n";
                return ladrunoContactFatal();
            }
            const ID &nodeTags = surf->getNodeTags();
            for (int k = 0; k < nodeTags.Size(); k++) {
                Node *sn = theDomain->getNode(nodeTags(k));
                if (sn == 0) {
                    opserr << "FATAL LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " not in domain; ABORTING (ADR-78 P1)\n";
                    return ladrunoContactFatal();
                }
                int nd = sn->getCrds().Size();
                if (sn->getNumberDOF() < nd) {
                    // The adapter couples the slave's first nd translational DOFs;
                    // a node with fewer DOFs than coordinates (ndf < ndm) would let
                    // base setID() leave trailing myID slots at 0 (contact assembled
                    // onto equation 0) and rigidPlaneGap() read past getTrialDisp().
                    // ADR-78 P1: abort rather than drop the node -- skipping left a
                    // rigid plane that silently held back only some of its slaves.
                    opserr << "FATAL LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " has ndf=" << sn->getNumberDOF()
                           << " < ndm=" << nd << "; ABORTING (ADR-78 P1)\n";
                    return ladrunoContactFatal();
                }
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn,
                                                            rp.muc,         // D2 viscous (0 ⇒ off)
                                                            rp.softScale,   // B1 SOFT=1 (0 ⇒ off)
                                                            theDomain,      // B1 engine mass cache
                                                            rp.tag);        // P5 warn-latch key
                if (fe == 0) return -5;
                theModel->addFE_Element(fe);
            }
        }
    }

    return count3;
}

void
LadrunoContactHandler::clearAll(void)
{
    // reset node DOF_Group pointers (as PlainHandler); the AnalysisModel owns and
    // deletes the DOF_Groups + FE_Elements (incl. our contact adapters).
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return;
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0)
        nodPtr->setDOF_GroupPtr(0);
}

int
LadrunoContactHandler::sendSelf(int, Channel &)
{
    // Deliberately empty, and staying that way (ADR-78 P2). The handler is
    // STATELESS -- every definition lives on the Domain's contact engine, and
    // that freight now rides Domain::sendSelf -> LadrunoContactDomain::sendSelf
    // (shipped in P2), the same channel stream that carries nodes/elements. A
    // handler shipped through the DDM analysis aggregate therefore has nothing
    // to add: its remote twin rebuilds everything from the received Domain.
    return 0;
}

int
LadrunoContactHandler::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}
