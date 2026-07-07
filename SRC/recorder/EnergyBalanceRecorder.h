/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
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


#ifndef EnergyBalanceRecorder_h
#define EnergyBalanceRecorder_h

// Written: jaabell (12/2024)
// Revision: C  (per-region + residual/closure, nodal mass & damping,
//               trapezoidal integration, pointer caching; nmb 05/2026)
//
// Description: records a structural-dynamics energy balance for the whole
// model and, optionally, for one or more MeshRegions. Each step writes
// EBR_NUM_ENERGY_COMPONENTS columns per block:
//
//   KE   kinetic energy            1/2 v^T M v          (instantaneous)
//   IE   internal energy           int F^T v dt         (accumulated)
//   DW   damping work              int v^T C v dt       (accumulated)
//   ULW  unbalanced-load work      int v^T Pext dt      (accumulated)
//   RES  residual                  ULW - (KE + IE + DW) (closure error)
//   ERR  relative residual [%]     100 * RES / Eref     (Eref = running
//                                                        max total energy)
//
// Column layout written each step:
//   [time]  (KE IE DW ULW RES ERR)_model  (KE IE DW ULW RES ERR)_region0  ...
//
// M and C include BOTH element matrices (Element::getMass/getDamp, which
// already carries element Rayleigh aM*M + bK*K) AND nodal contributions
// (Node::getMass, Node::getDamp = aM*nodalMass). Pext is the external
// applied *nodal* load (Node::getUnbalancedLoad after Domain::applyLoad).
//
// Known coverage gaps (RES exposes their effect, so they are visible):
//   * Modal damping (applied inside the integrator's solve, not via any
//     element/node getDamp) is NOT captured by DW.
//   * Element / distributed loads (eleLoad) assemble into the element
//     unbalance, not Node::unbalLoad, so they are NOT captured by ULW.
//   * Under _PARALLEL_PROCESSING each rank records only its own partition
//     (per-partition totals, ULW double-counted on shared nodes). Reduce
//     offline for a global balance. A warning is emitted at construction.
//
// Time integration of the three work integrals is trapezoidal
// (1/2 (rate_{n-1}+rate_n) dt); the first record() only seeds the rate so
// no spurious increment is taken over the initial time jump.
//
// Ladruno ADR-69 (-v2, opt-in): channel-aware layout. Model block becomes
//   KE_ele KE_nod IE DW_ele DW_nod ULW [E_inject] [E_lnvd] [E_modal] [E_hg]
//   RES ERR%
// where the element/nodal split quarantines the MPI-unsafe nodal terms
// (element terms reduce correctly across ranks; nodal terms are multiply-
// counted on shared boundary nodes), E_inject is the absorbing-boundary
// incident injection (derived from the published resisting-sign leak L:
// IE shown = IE_raw - L - E_hg, E_inject = -L), E_lnvd is the FLAC local
// non-viscous damping dissipation, E_modal (P2) is the integrator-published
// modal-damping dissipation (published from IncrementalIntegrator::commit();
// integrators that override commit() without chaining - HHT family, *_TP
// explicit - do not publish and their modal dissipation stays in RES), and
// E_hg (P2) is the recorder-PULLED hourglass-stabilization energy (summed
// element "hourglassEnergy" responses; present iff the model contains an
// element exposing that response at first record). Channel columns appear
// only when a producer declared them (LadrunoEnergyChannels.h) or, for
// E_hg, when the response probe hit. Region blocks keep the legacy 6
// columns in both modes. Without -v2 the output is byte-identical to v1.
// E_inject covers LysmerTriangle and ASDAbsorbingBoundary2D/3D (bottom +
// lateral) in full.
//
// Ladruno ADR-69 P2 (-ownedNodes <regionTag>, opt-in): count NODAL terms
// (nodal-mass KE, nodal alphaM Rayleigh DW, nodal-load ULW) only for nodes
// in the given MeshRegion. Under flat-per-rank MPI (openseesmp) shared
// boundary nodes are constructed on every touching rank and their nodal
// terms would be multiply-counted in a cross-rank sum; passing each rank
// the set of nodes it OWNS (canonical-rank convention - apeGmsh emits
// this) makes the per-rank nodal blocks sum exactly. Applies to BOTH
// layouts' nodal terms. Absent flag = all nodes owned (byte-identical).
// A missing region warns and counts NO nodes (loud, visible misconfig).
//
// Ladruno ADR-69 P2: under a detected MPI launcher (PMI_SIZE/
// OMPI_COMM_WORLD_SIZE/SLURM_NTASKS > 1) file-based outputs are suffixed
// per rank ("energy.txt" -> "energy.part-<rank>.txt") so ranks never race
// on one file; stitch offline (element block sums exactly; nodal block
// per the -ownedNodes contract).
//
// -v2 COLUMN LAYOUT REFLECTS PROCESS HISTORY, not just the current model:
// channel declarations are process-sticky (monotone, no reset on wipe — the
// recorder baselines energy DELTAS so totals never corrupt a later model,
// but the declared set only grows). In a multi-model interpreter session a
// later -v2 recorder shows a benign all-zero column for any channel a
// previous model declared (e.g. E_inject after an earlier Lysmer model).
// Parse -v2 output by the echoed/XML column names, never by position.



#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <TimeSeries.h>
#include <EnergyBalanceKernel.h>   // ebkernel::EnergyAccumulator + per-entity math

#include <vector>
#include <unordered_map>
#include <unordered_set>   // Ladruno ADR-69 P2: -ownedNodes membership

class Domain;
class FE_Datastore;
class Node;
class Element;
class MeshRegion;
class Response;   // Ladruno ADR-69 P2: hourglassEnergy pull cache

#define EBR_NUM_ENERGY_COMPONENTS 6

class EnergyBalanceRecorder: public Recorder
{
public:
    EnergyBalanceRecorder();
    EnergyBalanceRecorder(
        Domain &theDomain,
        OPS_Stream &theOutputHandler,
        bool echoTimeFlag = true,
        const ID &regions = ID(),
        bool v2Layout = false,        // Ladruno ADR-69: channel-aware split layout
        int ownedNodesRegion = -1);   // Ladruno ADR-69 P2: -ownedNodes region (-1 = all)

    ~EnergyBalanceRecorder();

    int record(int commitTag, double timeStamp);
    int flush();

    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);


protected:

private:
    int initialize(void);
    void buildCache(void);

    Domain *theDomain;
    OPS_Stream *theOutputHandler;

    Vector response;

    bool echoTimeFlag;       // include the time column in the output
    bool initializationDone; // response/accumulators sized and zeroed
    bool firstRecord;        // first record() only seeds the rate

    double time_last;

    // ---- whole-model accumulator (work integrals + closure/ERR) ----
    // The KE/IE/DW/ULW/RES/ERR math lives ONLY in ebkernel (one source of
    // truth, ADR D8); the recorder owns the per-scope accumulator state.
    ebkernel::EnergyAccumulator model_acc;

    // ---- Ladruno ADR-69 v2 layout (opt-in via -v2) ----
    // Legacy layout above stays byte-identical; the v2 path uses its own
    // accumulator with the element/nodal split + energy channels. Which
    // channel columns exist is fixed at initialize() from the registry's
    // declared set (producers declare at construction/setDomain — define the
    // integrator before the first analyze so -lnvd is visible by then).
    bool v2Layout;
    bool chInject;               // ABSORB_LEAK declared -> E_inject column
    bool chLnvd;                 // LNVD_WORK declared  -> E_lnvd column
    bool chModal;                // Ladruno ADR-69 P2: MODAL_WORK -> E_modal
    bool chHourglass;            // Ladruno ADR-69 P2: response probe hit -> E_hg
    int  numModelCols;           // model-block width (6 legacy or v2 width)
    ebkernel::EnergyAccumulatorV2 model_acc2;

    // Ladruno ADR-69 P2: hourglass "energy" pull (mechanism C) - Response*
    // objects resolved by probing every element for a "hourglassEnergy"
    // response; rebuilt (and the old ones deleted) with the cache.
    // P2.1 removal safety: each entry keeps the element TAG so record() can
    // verify getElement(tag) still returns the cached pointer before the
    // virtual getResponse() call - `remove element` deletes the Element the
    // cached ElementResponse wraps, and nothing routes domainChanged() to
    // recorders (Domain::hasDomainChanged() is stateful/consumed by the
    // Analysis, so the recorder cannot poll it; re-validate by tag instead,
    // the LadrunoMonitorRecorder pattern).
    struct HgEntry {
        int tag;
        Element *ele;
        Response *res;
    };
    std::vector<HgEntry> hgCache;

    // Ladruno ADR-69 P2: -ownedNodes. Region tag (-1 = inactive) and the
    // resolved node set (by TAG - removal/pointer-reuse safe, P2.1); nodal
    // terms are only counted for members.
    int ownedRegionTag;
    std::unordered_set<int> ownedNodes;

    // ---- per-region accumulators ----
    ID regionTags;
    int numRegions;
    std::vector<ebkernel::EnergyAccumulator> region_acc;

    // ---- caching (rebuilt on domainChanged): region membership + scratch --
    // Ladruno ADR-69 P2.1: keyed by TAG, not pointer - a removed entity's
    // freed pointer could be reused by a new allocation and silently inherit
    // the old region binning; tags cannot alias. cachedNumElements/Nodes are
    // structural staleness sentinels checked every record() (domainChanged()
    // is never delivered to recorders; see HgEntry note above).
    bool cacheValid;
    std::unordered_map<int, std::vector<int> > elemRegions;
    std::unordered_map<int, std::vector<int> > nodeRegions;
    int cachedNumElements;
    int cachedNumNodes;
    int maxNumDOF;
    Vector velScratch;
};

#endif
