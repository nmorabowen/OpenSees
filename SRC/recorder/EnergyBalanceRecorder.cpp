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

// Written: jaabell (12/2024); per-region + closure rewrite nmb (05/2026)
//
// Description: see EnergyBalanceRecorder.h for the energy-balance definition,
// column layout, and the documented coverage gaps (modal damping, element
// loads, parallel per-partition behaviour).

#include <EnergyBalanceRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include "Ladruno_LaunchEnv.h"   // Ladruno WP-163 M4/MP-3
#include "LadrunoEnergyBalanceMPI.h"   // Ladruno: optional -v2 columns agreed across ranks
#include <MeshRegion.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <RelativeMassDamping.h>   // Ladruno (relativeMassDamping)

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <TCP_Stream.h>

#include <elementAPI.h>
#include <classTags.h>   // Ladruno: RECORDER_TAGS_EnergyBalanceRecorder

#include <Response.h>      // Ladruno ADR-69 P2: hourglassEnergy pull
#include <Information.h>   // Ladruno ADR-69 P2: response data access
#include <DummyStream.h>   // Ladruno ADR-69 P2: silent setResponse probe

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string>    // Ladruno ADR-69 P2: per-rank filename suffix
#include <sstream>   // Ladruno ADR-69 P2
#include <cstdlib>   // Ladruno ADR-69 P2: std::getenv
#include <profiler/ProfilerMacros.h>   // Ladruno WP-196: commit-phase split
#include <profiler/Profiler.h>         // Ladruno WP-196: deep() gate
#include <LadrunoThreads.h>            // Ladruno WP-196: the WP-192 parallel-for
#include <OPS_Globals.h>               // Ladruno WP-196: ops_TheActiveElement

void*
OPS_EnergyBalanceRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: recorder EnergyBalance -file <fileName> <-time> "
                  "<-region tag ...> <-v2>\n";
        return 0;
    }

    OPS_Stream *theOutputStream = 0;
    const char* filename = 0;

    const int STANDARD_STREAM = 0;
    const int DATA_STREAM = 1;
    const int XML_STREAM = 2;
    const int BINARY_STREAM = 4;
    const int DATA_STREAM_CSV = 5;
    const int TCP_STREAM = 6;
    const int DATA_STREAM_ADD = 7;

    int eMode = STANDARD_STREAM;

    bool echoTimeFlag = false;
    bool doScientific = false;
    int precision = 6;
    bool closeOnWrite = false;
    bool v2Layout = false;        // Ladruno ADR-69
    int ownedNodesRegion = -1;    // Ladruno ADR-69 P2: -ownedNodes region tag
    bool relMassSplit = false;     // Ladruno (relativeMassDamping)

    const char *inetAddr = 0;
    int inetPort = 0;

    ID regions(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-v2") == 0) {
            // Ladruno ADR-69: channel-aware layout with element/nodal split
            v2Layout = true;
        }
        else if (strcmp(option, "-relativeMassSplit") == 0) {
            // Ladruno (relativeMassDamping): DW_rel / W_ref columns (-v2)
            relMassSplit = true;
        }
        else if (strcmp(option, "-scientific") == 0) {
            doScientific = true;
        }
        else if (strcmp(option, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM;
        }
        else if (strcmp(option, "-fileAdd") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM_ADD;
        }
        else if (strcmp(option, "-closeOnWrite") == 0) {
            closeOnWrite = true;
        }
        else if (strcmp(option, "-csv") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = DATA_STREAM_CSV;
        }
        else if (strcmp(option, "-tcp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                inetAddr = OPS_GetString();
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &inetPort) < 0) {
                    opserr << "WARNING: EnergyBalance -tcp failed to read inetPort\n";
                    return 0;
                }
            }
            eMode = TCP_STREAM;
        }
        else if (strcmp(option, "-xml") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = XML_STREAM;
        }
        else if (strcmp(option, "-binary") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
            eMode = BINARY_STREAM;
        }
        else if (strcmp(option, "-precision") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &precision) < 0) {
                    opserr << "WARNING: EnergyBalance -precision failed to read value\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-region") == 0) {
            // -region <tag> ; repeatable, e.g. ... -region 1 -region 2
            // Each named region adds one (KE IE DW ULW RES ERR) block after
            // the whole-model total. Region energy is summed over the region's
            // elements (KE, IE, DW) and nodes (KE, DW, ULW).
            int tag;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING: EnergyBalance -region requires an integer region tag\n";
                    return 0;
                }
                regions[regions.Size()] = tag;
            }
        }
        else if (strcmp(option, "-ownedNodes") == 0) {
            // Ladruno ADR-69 P2: count nodal terms (KE_nod, DW_nod, ULW)
            // only for nodes in this MeshRegion - the per-rank ownership
            // gate for flat-per-rank MPI (see header). -1 stays "all".
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &ownedNodesRegion) < 0) {
                    opserr << "WARNING: EnergyBalance -ownedNodes requires an "
                              "integer region tag\n";
                    return 0;
                }
            }
        }
    }

    // Ladruno ADR-69 P2: per-rank filename suffix under a detected MPI
    // launcher, so flat-per-rank (openseesmp) processes never race on one
    // file. Same env probe order as LadrunoRecorder (Intel/MS-MPI, OpenMPI,
    // then SLURM last - sbatch exports SLURM_NTASKS into the batch shell
    // itself, so a real launcher's own variables take precedence).
    // "energy.txt" -> "energy.part-<rank>.txt" (suffix before the last
    // extension; appended when there is none).
    std::string partFilename;
    if (filename != 0) {
        // WP-163 M4/MP-3: the shared probe (Ladruno_LaunchEnv.h) — SLURM only
        // inside an srun step; SIZE > 1 without a valid RANK is refused.
        int rank = 0, size = 1;
        std::string source, error;
        const ladruno::launch::Status st =
            ladruno::launch::detectRank(rank, size, source, error);
        if (st == ladruno::launch::Inconsistent) {
            opserr << "WARNING: EnergyBalance: " << error.c_str()
                   << " - cannot tell which per-rank file this process owns\n";
            return 0;
        }
        if (st == ladruno::launch::Launched) {
            std::string stem(filename), ext;
            const size_t dot = stem.find_last_of('.');
            if (dot != std::string::npos && dot > 0) {
                ext = stem.substr(dot);
                stem.erase(dot);
            }
            std::stringstream ss;
            ss << stem << ".part-" << rank << ext;
            partFilename = ss.str();
            filename = partFilename.c_str();
        }
    }

    // data handler
    if (eMode == DATA_STREAM && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_ADD && filename != 0)
        theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_CSV && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
    else if (eMode == XML_STREAM && filename != 0)
        theOutputStream = new XmlFileStream(filename);
    else if (eMode == BINARY_STREAM && filename != 0)
        theOutputStream = new BinaryFileStream(filename);
    else if (eMode == TCP_STREAM && inetAddr != 0)
        theOutputStream = new TCP_Stream(inetPort, inetAddr);
    else
        theOutputStream = new StandardStream();

    theOutputStream->setPrecision(precision);

    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        delete theOutputStream;
        return 0;
    }

#ifdef _PARALLEL_PROCESSING
    opserr << "EnergyBalanceRecorder: WARNING under parallel processing the energy "
              "balance is computed PER PARTITION (local domain only) - totals are not "
              "reduced across ranks and ULW is double-counted on shared nodes. Reduce "
              "the per-rank output offline for a global balance.\n";
#endif

    EnergyBalanceRecorder* recorder = new EnergyBalanceRecorder(
        *domain, *theOutputStream, echoTimeFlag, regions, v2Layout,
        ownedNodesRegion);
    if (relMassSplit) {
        if (!v2Layout)
            opserr << "EnergyBalanceRecorder: WARNING -relativeMassSplit needs -v2; ignored" << endln;
        else
            recorder->setRelativeMassSplit(true);
    }

    return recorder;
}


EnergyBalanceRecorder::EnergyBalanceRecorder()
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(0), theOutputHandler(0),
      response(EBR_NUM_ENERGY_COMPONENTS),
      echoTimeFlag(true), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      v2Layout(false), chInject(false), chLnvd(false),
      chModal(false), chHourglass(false), chMassScale(false),
      numModelCols(EBR_NUM_ENERGY_COMPONENTS), model_acc2(),
      hgCache(), ownedRegionTag(-1), ownedNodes(),
      regionTags(), numRegions(0),
      region_acc(),
      cacheValid(false), cachedNumElements(-1), cachedNumNodes(-1),
      maxNumDOF(0), velScratch(),
      sweepNumSafe(0), sweepAnnounced(0)   // Ladruno WP-196
{

}


EnergyBalanceRecorder::EnergyBalanceRecorder(
    Domain &theDomain_,
    OPS_Stream &theOutputHandler_,
    bool echoTimeFlag_,
    const ID &regions,
    bool v2Layout_,
    int ownedNodesRegion)
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(&theDomain_), theOutputHandler(&theOutputHandler_),
      response((echoTimeFlag_ ? 1 : 0) + EBR_NUM_ENERGY_COMPONENTS * (1 + regions.Size())),
      echoTimeFlag(echoTimeFlag_), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      v2Layout(v2Layout_), chInject(false), chLnvd(false),
      chModal(false), chHourglass(false), chMassScale(false),
      numModelCols(EBR_NUM_ENERGY_COMPONENTS), model_acc2(),
      hgCache(), ownedRegionTag(ownedNodesRegion), ownedNodes(),
      regionTags(regions), numRegions(regions.Size()),
      region_acc((size_t)regions.Size()),
      cacheValid(false), cachedNumElements(-1), cachedNumNodes(-1),
      maxNumDOF(0), velScratch(),
      sweepNumSafe(0), sweepAnnounced(0)   // Ladruno WP-196
{
    response.Zero();
}


EnergyBalanceRecorder::~EnergyBalanceRecorder()
{
    // Ladruno ADR-69 P2: the hourglass pull cache owns its Response objects
    for (size_t i = 0; i < hgCache.size(); ++i)
        delete hgCache[i].res;
    hgCache.clear();

    if (theOutputHandler != 0) {
        if (initializationDone)
            theOutputHandler->endTag(); // closes the "Data" tag opened in initialize()
        delete theOutputHandler;
    }
}


// NOTE (ADR D8): the per-entity energy math (addElementEnergy / addNodeEnergy)
// and the per-scope trapezoidal/closure arithmetic now live in ONE place,
// ebkernel (EnergyBalanceKernel.h), shared with ladruno::EnergyBalanceSource.
// This recorder calls ebkernel::addElementEnergy / addNodeEnergy in the sweep
// below and ebkernel::EnergyAccumulator::step() for the integration + closure.


void
EnergyBalanceRecorder::buildCache(void)
{
    elemRegions.clear();
    nodeRegions.clear();

    // region membership: element/node pointer -> list of region indices
    for (int r = 0; r < numRegions; ++r) {
        MeshRegion *region = theDomain->getRegion(regionTags(r));
        if (region == 0) {
            opserr << "EnergyBalanceRecorder - WARNING region " << regionTags(r)
                   << " not found in the domain; its energy columns will be zero\n";
            continue;
        }
        const ID &elems = region->getElements();
        for (int i = 0; i < elems.Size(); ++i) {
            // keyed by TAG (removal-safe, no pointer reuse aliasing)
            if (theDomain->getElement(elems(i)) != 0)
                elemRegions[elems(i)].push_back(r);
        }
        const ID &rnodes = region->getNodes();
        for (int i = 0; i < rnodes.Size(); ++i) {
            if (theDomain->getNode(rnodes(i)) != 0)
                nodeRegions[rnodes(i)].push_back(r);
        }
    }

    // size the velocity scratch to the largest element DOF count
    // Ladruno ADR-69 P2: in the same sweep, probe every element for a
    // "hourglassEnergy" response (mechanism C recorder-pull); the non-null
    // hits form the E_hg cache. The probe stream is a DummyStream so the
    // setResponse tag machinery never touches the real output.
    for (size_t i = 0; i < hgCache.size(); ++i)
        delete hgCache[i].res;
    hgCache.clear();
    static DummyStream probeStream;
    static const char *hgArgv[1] = { "hourglassEnergy" };
    maxNumDOF = 0;
    Element *ele;
    ElementIter &elements = theDomain->getElements();
    while ((ele = elements()) != 0) {
        // Ladruno WP-163 M1: skip the SP ShadowSubdomains P0's iterator yields
        // (no node pointers; a setResponse probe would go to the remote actor).
        if (ele->isSubdomain())
            continue;
        const int n = ele->getNumDOF();
        if (n > maxNumDOF)
            maxNumDOF = n;
        if (v2Layout) {
            Response *hgRes = ele->setResponse(hgArgv, 1, probeStream);
            if (hgRes != 0) {
                HgEntry e = { ele->getTag(), ele, hgRes };
                hgCache.push_back(e);
            }
        }
    }
    if (velScratch.Size() < maxNumDOF)
        velScratch.resize(maxNumDOF);

    // Ladruno ADR-69 P2.1: structural staleness sentinels - record() compares
    // these against the live domain every call and rebuilds on any change
    // (recorders never receive domainChanged(); see header).
    cachedNumElements = theDomain->getNumElements();
    cachedNumNodes = theDomain->getNumNodes();

    // Ladruno ADR-69 P2: resolve the -ownedNodes region into a node set.
    // A missing region warns and owns NO nodes (loud, visible misconfig -
    // the nodal columns go to zero rather than silently double-counting).
    ownedNodes.clear();
    if (ownedRegionTag >= 0) {
        MeshRegion *owned = theDomain->getRegion(ownedRegionTag);
        if (owned == 0) {
            opserr << "EnergyBalanceRecorder - WARNING -ownedNodes region "
                   << ownedRegionTag << " not found; NO nodal terms will be "
                      "counted on this rank\n";
        }
        else {
            const ID &onodes = owned->getNodes();
            for (int i = 0; i < onodes.Size(); ++i) {
                if (theDomain->getNode(onodes(i)) != 0)
                    ownedNodes.insert(onodes(i));   // by TAG (P2.1)
            }
        }
    }

    cacheValid = true;
}


int
EnergyBalanceRecorder::record(int commitTag, double timeStamp)
{
    if (theDomain == 0)
        return 0;

    if (theOutputHandler == 0) {
        opserr << "EnergyBalanceRecorder::record() - no output handler has been set\n";
        return -1;
    }

    if (initializationDone == false) {
        if (this->initialize() != 0) {
            opserr << "EnergyBalanceRecorder::record() - failed in initialize()\n";
            return -1;
        }
    }

    // Ladruno ADR-69 P2.1: removal/addition safety. Nothing routes
    // domainChanged() to recorders (Domain::record() calls them directly;
    // the analysis propagates domain changes only to its own components),
    // and Domain::hasDomainChanged() is stateful (the Analysis consumes it),
    // so re-validate structurally on every record: entity-count change or
    // any hourglass cache entry whose tag no longer resolves to the SAME
    // Element* (removed, or removed-and-replaced) forces a cache rebuild
    // BEFORE the hgCache getResponse() virtual calls below could touch a
    // freed element.
    if (cacheValid == true) {
        if (theDomain->getNumElements() != cachedNumElements ||
            theDomain->getNumNodes() != cachedNumNodes)
            cacheValid = false;
        else
            for (size_t i = 0; i < hgCache.size(); ++i)
                if (theDomain->getElement(hgCache[i].tag) != hgCache[i].ele) {
                    cacheValid = false;
                    break;
                }
    }

    if (cacheValid == false)
        this->buildCache();

    const double dT = timeStamp - time_last;

    int timeOffset = 0;
    if (echoTimeFlag == true) {
        timeOffset = 1;
        response(0) = timeStamp;
    }

    //
    // single sweep: compute each element/node contribution once (via the shared
    // ebkernel per-entity math), accumulate into the whole-model total and bin
    // into the regions it belongs to.
    //
    double g_ke = 0., g_internal_rate = 0., g_damping_rate = 0., g_unbalanced_rate = 0.;
    // Ladruno ADR-69 (-v2): element/nodal split sums, accumulated ALONGSIDE
    // the legacy variables (which keep their exact accumulation order so the
    // legacy layout stays byte-identical). v2_ms_ke (P3) is the fictitious
    // added-mass KE share, accumulated by the kernel through the optional
    // out-param (null on the legacy path -> zero extra work there).
    double v2_ke_ele = 0., v2_ke_nod = 0., v2_dw_ele = 0., v2_dw_nod = 0.;
    double v2_ms_ke = 0.;
    double rr_rel_rate = 0., rr_ref_rate = 0.;   // Ladruno (relativeMassDamping)
    std::vector<double> step_region_ke((size_t)numRegions, 0.0);
    std::vector<double> step_region_internal_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_damping_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_unbalanced_rate((size_t)numRegions, 0.0);

    {
        // Ladruno WP-196: the per-element terms are formed first (threaded when
        // allowed, see sweepElements) and reduced here in the live element
        // order, so every sum below sees the same operands in the same order
        // as the former one-pass loop. ebkernel::elementEnergyTerms is bitwise
        // addElementEnergy on zeroed locals; t.ms is what it added to *msKE.
        OPS_PROFILE_SCOPE("eb.elements");
        this->sweepElements();
        const bool msActive = v2Layout &&
            Ladruno::MassScalingEnergyRegistry::instance().active();
        for (size_t e = 0; e < sweepEles.size(); ++e) {
            Element *ele = sweepEles[e];
            const ebkernel::ElementEnergyTerms &t = sweepTerms[e];
            const double ke = t.ke, ie = t.ie, dw = t.dw;
            if (msActive)
                v2_ms_ke += t.ms;
            g_ke += ke; g_internal_rate += ie; g_damping_rate += dw;
            if (v2Layout) { v2_ke_ele += ke; v2_dw_ele += dw; }

            if (numRegions > 0) {
                std::unordered_map<int, std::vector<int> >::const_iterator it =
                    elemRegions.find(ele->getTag());
                if (it != elemRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke[(size_t)r] += ke;
                        step_region_internal_rate[(size_t)r] += ie;
                        step_region_damping_rate[(size_t)r] += dw;
                    }
            }
        }
    }

    {
        OPS_PROFILE_SCOPE("eb.nodes");   // Ladruno WP-196
        Node *node;
        NodeIter &nodes = theDomain->getNodes();
        while ((node = nodes()) != 0) {
            // Ladruno ADR-69 P2: -ownedNodes gate - nodal terms are only
            // counted for owned nodes (both layouts; see header contract)
            if (ownedRegionTag >= 0 &&
                ownedNodes.find(node->getTag()) == ownedNodes.end())
                continue;
            double ke = 0., dw = 0., ulw = 0.;
            ebkernel::addNodeEnergy(node, ke, dw, ulw,
                                    v2Layout ? &v2_ms_ke : 0,
                                    chRelMass ? &rr_rel_rate : 0,
                                    chRelMass ? &rr_ref_rate : 0);
            g_ke += ke; g_damping_rate += dw; g_unbalanced_rate += ulw;
            if (v2Layout) { v2_ke_nod += ke; v2_dw_nod += dw; }

            if (numRegions > 0) {
                std::unordered_map<int, std::vector<int> >::const_iterator it =
                    nodeRegions.find(node->getTag());
                if (it != nodeRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke[(size_t)r] += ke;
                        step_region_damping_rate[(size_t)r] += dw;
                        step_region_unbalanced_rate[(size_t)r] += ulw;
                    }
            }
        }
    }

    //
    // integrate the work rates + close the balance through the shared kernel
    // accumulator (trapezoidal, first-record seeding, per-scope running-max ERR
    // are all defined ONCE in ebkernel::EnergyAccumulator::step). Fill the
    // response: whole-model block then one block per region.
    //
    if (v2Layout) {
        // Ladruno ADR-69: channel-aware model block. Channel totals are
        // cumulative registry reads; the accumulator snapshots its baseline
        // on the first record so a fresh recorder starts from zero.
        const Ladruno::EnergyChannelRegistry &reg =
            Ladruno::EnergyChannelRegistry::instance();

        // Ladruno ADR-69 P2: pull the hourglass-stabilization energy total
        // (mechanism C) - the element accumulates at commit cadence, the
        // recorder just samples the running scalar, so no aliasing.
        double hg_total = 0.0;
        for (size_t i = 0; i < hgCache.size(); ++i) {
            Response *hgRes = hgCache[i].res;
            if (hgRes->getResponse() >= 0) {
                const Information &info = hgRes->getInformation();
                if (info.theVector != 0 && info.theVector->Size() > 0)
                    hg_total += (*info.theVector)(0);
            }
        }

        double out[ebkernel::NUM_V2_SLOTS];
        model_acc2.step(dT, firstRecord,
                        v2_ke_ele, v2_ke_nod,
                        g_internal_rate, v2_dw_ele, v2_dw_nod,
                        g_unbalanced_rate,
                        reg.energy(Ladruno::EnergyChannelRegistry::ABSORB_LEAK),
                        reg.energy(Ladruno::EnergyChannelRegistry::LNVD_WORK),
                        reg.energy(Ladruno::EnergyChannelRegistry::MODAL_WORK),
                        hg_total,
                        v2_ms_ke,
                        out);
        int col = timeOffset;
        response(col++) = out[ebkernel::V2_KE_ELE];
        response(col++) = out[ebkernel::V2_KE_NOD];
        response(col++) = out[ebkernel::V2_IE];
        response(col++) = out[ebkernel::V2_DW_ELE];
        response(col++) = out[ebkernel::V2_DW_NOD];
        response(col++) = out[ebkernel::V2_ULW];
        if (chInject)
            response(col++) = out[ebkernel::V2_E_INJECT];
        if (chLnvd)
            response(col++) = out[ebkernel::V2_E_LNVD];
        if (chModal)
            response(col++) = out[ebkernel::V2_E_MODAL];
        if (chHourglass)
            response(col++) = out[ebkernel::V2_E_HG];
        if (chMassScale)
            response(col++) = out[ebkernel::V2_KE_MS];
        if (chRelMass) {
            // Ladruno (relativeMassDamping): advisory split of DW_nod
            if (firstRecord) {
                rr_dwRel += rr_rel_rate * dT;
                rr_wRef  += rr_ref_rate * dT;
            } else {
                rr_dwRel += 0.5 * (rr_prevRel + rr_rel_rate) * dT;
                rr_wRef  += 0.5 * (rr_prevRef + rr_ref_rate) * dT;
            }
            rr_prevRel = rr_rel_rate;
            rr_prevRef = rr_ref_rate;
            response(col++) = rr_dwRel;
            response(col++) = rr_wRef;
        }
        response(col++) = out[ebkernel::V2_RES];
        response(col++) = out[ebkernel::V2_ERR];
    }
    else {
        double out[ebkernel::NUM_ENERGY_COMPONENTS];
        model_acc.step(dT, firstRecord,
                       g_ke, g_internal_rate, g_damping_rate, g_unbalanced_rate, out);
        int col = timeOffset;
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            response(col + c) = out[c];
    }

    for (int r = 0; r < numRegions; ++r) {
        double out[ebkernel::NUM_ENERGY_COMPONENTS];
        region_acc[(size_t)r].step(dT, firstRecord,
                                   step_region_ke[(size_t)r],
                                   step_region_internal_rate[(size_t)r],
                                   step_region_damping_rate[(size_t)r],
                                   step_region_unbalanced_rate[(size_t)r], out);
        // Region blocks keep the legacy 6 columns in both layouts; in v2 the
        // model block may be wider, so offset by its actual width.
        const int col = timeOffset + numModelCols + EBR_NUM_ENERGY_COMPONENTS * r;
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            response(col + c) = out[c];
    }

    firstRecord = false;

    { OPS_PROFILE_SCOPE("eb.write");   // Ladruno WP-196
    theOutputHandler->write(response); }

    time_last = timeStamp;

    return 0;
}


// Ladruno WP-196: the element sweep. Snapshot the live element order, form
// every element's terms into its own slot, and leave the reduction to record().
// Threaded runs follow the WP-192 shape: maximal runs of thread-safe elements
// are formed on the threads, every other element alone on this thread in its
// place, so an element that is not audited sees exactly its serial context.
void
EnergyBalanceRecorder::sweepElements(void)
{
    sweepEles.clear();
    {
        Element *ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0)
            sweepEles.push_back(ele);
    }
    const int nEle = (int)sweepEles.size();
    if ((int)sweepTerms.size() != nEle)
        sweepTerms.resize((size_t)nEle);
    const bool legacy = ebkernel::legacyEnergyCheck();

    const int nThreads = this->sweepThreads(nEle);
    if (nThreads <= 1) {
        for (int i = 0; i < nEle; ++i)
            ebkernel::elementEnergyTerms(sweepEles[i], legacy, sweepTerms[i]);
        return;
    }

    struct Ctx {
        Element **e;
        ebkernel::ElementEnergyTerms *t;
        bool legacy;
    };
    struct Body {
        static int run(int i, void *c) {
            Ctx *x = (Ctx *)c;
            ebkernel::elementEnergyTerms(x->e[i], x->legacy, x->t[i]);
            return 0;
        }
    };
    Element *activeBefore = ops_TheActiveElement;
    int i = 0;
    while (i < nEle) {
        if (!sweepSafe[(size_t)i]) {
            ebkernel::elementEnergyTerms(sweepEles[i], legacy, sweepTerms[i]);
            i++;
            continue;
        }
        int j = i;
        while (j < nEle && sweepSafe[(size_t)j])
            j++;
        Ctx ctx = { &sweepEles[(size_t)i], &sweepTerms[(size_t)i], legacy };
        ladruno_parallelFor(j - i, nThreads, 16, &Body::run, &ctx, 0);
        i = j;
    }
    // the serial sweep does not touch ops_TheActiveElement; neither may this
    ops_TheActiveElement = activeBefore;
}


// Ladruno WP-196: how many threads the element sweep may use (1 = serial).
// The WP-192 gates (requested count, OpenMP, the MPI hybrid acknowledgement,
// the LADRUNO_THREADS_LOOPS mask bit 16, a plain Domain, no deep profiler),
// plus: no consistent selective-mass-scaling registry (its per-element
// lookups are not audited), and at least two thread-safe elements.
int
EnergyBalanceRecorder::sweepThreads(int nEle)
{
    const char *why = 0;
    const int n = ladruno_loopThreads(&why);
    if (n <= 1)
        return 1;
    if (!ladruno_loopEnabled(LADRUNO_LOOP_ENERGY))
        return 1;
    if (theDomain->ladrunoThreadedUpdateAllowed() == false)
        return 1;
    if (ops_profiler::theProfiler().deep())
        return 1;
    if (Ladruno::MassScalingEnergyRegistry::instance().active())
        return 1;

    // per-element opt-in, re-derived whenever the element list changed
    if (safeFor != sweepEles) {
        safeFor = sweepEles;
        sweepSafe.assign((size_t)nEle, 0);
        sweepNumSafe = 0;
        for (int i = 0; i < nEle; ++i) {
            Element *e = sweepEles[(size_t)i];
            if (e->isSubdomain() || e->getNodePtrs() == 0)
                continue;
            // ladrunoPrepareThreadedForm allocates the Rayleigh pool slot here,
            // serially, as the WP-192 form loops do
            if (e->ladrunoThreadSafeForm() && e->ladrunoPrepareThreadedForm()) {
                sweepSafe[(size_t)i] = 1;
                sweepNumSafe++;
            }
        }
    }
    if (sweepNumSafe < 2)
        return 1;
    if (sweepAnnounced != n) {
        opserr << "EnergyBalanceRecorder: element sweep threaded on " << n
               << " threads: " << sweepNumSafe << " of " << nEle
               << " elements formed in parallel, summed in serial order (WP-196)\n";
        sweepAnnounced = n;
    }
    return n;
}


int
EnergyBalanceRecorder::setDomain(Domain &theDom)
{
    theDomain = &theDom;
    time_last = theDomain->getCurrentTime();
    firstRecord = true;
    cacheValid = false;
    return 0;
}


int
EnergyBalanceRecorder::domainChanged(void)
{
    // region membership / element DOF sizes may have changed: rebuild lazily.
    cacheValid = false;
    return 0;
}


int
EnergyBalanceRecorder::sendSelf(int commitTag, Channel &theChannel)
{
    // Recorders are not serialized across processes in OpenSees - each rank
    // constructs its own recorder locally via OPS_EnergyBalanceRecorder.
    return 0;
}


int
EnergyBalanceRecorder::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
    return 0;
}


int
EnergyBalanceRecorder::initialize(void)
{
    if (theDomain == 0) {
        opserr << "EnergyBalanceRecorder::initialize() - the domain has not been set\n";
        return -1;
    }

    // Ladruno ADR-69 (-v2): fix the channel columns from the registry's
    // declared set. Producers declare at construction/setDomain (or, for
    // MODAL_WORK, at the first base commit - which precedes this initialize
    // since the first record() runs inside that same commitDomain), so
    // anything built before analyze is visible here. E_hg is decided from
    // the response-probe cache instead of the registry - the column
    // reflects the CURRENT model, so it must be probed before the widths
    // are fixed (buildCache here; record() skips its lazy rebuild).
    if (v2Layout) {
        if (cacheValid == false)
            this->buildCache();
        const Ladruno::EnergyChannelRegistry &reg =
            Ladruno::EnergyChannelRegistry::instance();
        chInject = reg.declared(Ladruno::EnergyChannelRegistry::ABSORB_LEAK);
        chLnvd   = reg.declared(Ladruno::EnergyChannelRegistry::LNVD_WORK);
        chModal  = reg.declared(Ladruno::EnergyChannelRegistry::MODAL_WORK);
        chHourglass = !hgCache.empty();
        // Ladruno ADR-69 P3: KE_ms gated on the SMS registry's CURRENT blocks
        // (pull-style, like E_hg — the registry clears on integrator teardown,
        // so this reflects the live model, not process history). The SMS
        // integrator's domainChanged has published by the first record; an
        // integrator defined after the first record misses its column (same
        // documented limitation as the registry channels).
        chMassScale =
            Ladruno::MassScalingEnergyRegistry::instance().anyActive();
        chRelMass = forceRelMass || (OPS_numRelativeMassDampings() > 0);   // Ladruno (relativeMassDamping)
        // Ladruno: the flags above are RANK-LOCAL (a rank that owns no DRM /
        // absorbing element has no ABSORB_LEAK producer, and so on). Agree them
        // across ranks so every interpreter-per-rank file and echo has the same
        // layout; a rank with no local contribution writes 0 there. The kernel
        // folds every channel into RES whether or not its column is written, so
        // RES/ERR% are unchanged. Identity outside openseesmp-style MPI.
        {
            int fl[6] = { chInject ? 1 : 0, chLnvd ? 1 : 0, chModal ? 1 : 0,
                          chHourglass ? 1 : 0, chMassScale ? 1 : 0,
                          chRelMass ? 1 : 0 };
            ladrunoEnergyBalanceAgreeMax(fl, 6);
            chInject = fl[0] != 0; chLnvd = fl[1] != 0; chModal = fl[2] != 0;
            chHourglass = fl[3] != 0; chMassScale = fl[4] != 0;
            chRelMass = fl[5] != 0;
        }
        numModelCols = 6 + (chInject ? 1 : 0) + (chLnvd ? 1 : 0)
                     + (chModal ? 1 : 0) + (chHourglass ? 1 : 0)
                     + (chMassScale ? 1 : 0) + (chRelMass ? 2 : 0) + 2;

        // Coverage (ADR-69 P1.5/P1.6): ASDAbsorbingBoundary2D/3D publish
        // both their BOTTOM compliant-base injection (addBaseActions) and
        // their LATERAL free-field transfer (addRffToSoil) to E_inject,
        // same channel as LysmerTriangle - no unaccounted ASD injection
        // path remains, so no coverage warning is needed here.
    }
    else {
        numModelCols = EBR_NUM_ENERGY_COMPONENTS;
    }

    // size and zero the response and all per-scope accumulators
    const int timeOffset = (echoTimeFlag == true) ? 1 : 0;
    response.resize(timeOffset + numModelCols
                    + EBR_NUM_ENERGY_COMPONENTS * numRegions);
    response.Zero();

    model_acc.reset();
    model_acc2.reset();
    rr_dwRel = rr_wRef = rr_prevRel = rr_prevRef = 0.0;   // Ladruno (relativeMassDamping)
    region_acc.assign((size_t)numRegions, ebkernel::EnergyAccumulator());

    // self-describing column metadata (consumed by XML/binary streams; text
    // streams ignore tags - the layout is echoed once below and documented
    // in the header).
    static const char *comp[EBR_NUM_ENERGY_COMPONENTS] =
        { "KE", "IE", "DW", "ULW", "RES", "ERR%" };

    if (echoTimeFlag == true) {
        theOutputHandler->tag("TimeOutput");
        theOutputHandler->tag("ResponseType", "time");
        theOutputHandler->endTag();
    }

    theOutputHandler->tag("EnergyOutput");
    theOutputHandler->attr("region", "model");
    if (v2Layout) {
        theOutputHandler->tag("ResponseType", "KE_ele");
        theOutputHandler->tag("ResponseType", "KE_nod");
        theOutputHandler->tag("ResponseType", "IE");
        theOutputHandler->tag("ResponseType", "DW_ele");
        theOutputHandler->tag("ResponseType", "DW_nod");
        theOutputHandler->tag("ResponseType", "ULW");
        if (chInject) theOutputHandler->tag("ResponseType", "E_inject");
        if (chLnvd)   theOutputHandler->tag("ResponseType", "E_lnvd");
        if (chModal)  theOutputHandler->tag("ResponseType", "E_modal");
        if (chHourglass) theOutputHandler->tag("ResponseType", "E_hg");
        if (chMassScale) theOutputHandler->tag("ResponseType", "KE_ms");
        if (chRelMass) { theOutputHandler->tag("ResponseType", "DW_rel"); theOutputHandler->tag("ResponseType", "W_ref"); }
        theOutputHandler->tag("ResponseType", "RES");
        theOutputHandler->tag("ResponseType", "ERR%");
    } else {
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            theOutputHandler->tag("ResponseType", comp[c]);
    }
    theOutputHandler->endTag();

    for (int r = 0; r < numRegions; ++r) {
        theOutputHandler->tag("EnergyOutput");
        theOutputHandler->attr("region", regionTags(r));
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            theOutputHandler->tag("ResponseType", comp[c]);
        theOutputHandler->endTag();
    }

    theOutputHandler->tag("Data");

    // echo the column layout once for users of the plain-text sidecar
    opserr << "EnergyBalanceRecorder: columns =";
    if (echoTimeFlag == true) opserr << " time";
    if (v2Layout) {
        opserr << " [model: KE_ele KE_nod IE DW_ele DW_nod ULW";
        if (chInject) opserr << " E_inject";
        if (chLnvd)   opserr << " E_lnvd";
        if (chModal)  opserr << " E_modal";
        if (chHourglass) opserr << " E_hg";
        if (chMassScale) opserr << " KE_ms";
        if (chRelMass)    opserr << " DW_rel W_ref";
        opserr << " RES ERR%]";
    } else {
        opserr << " [model: KE IE DW ULW RES ERR%]";
    }
    for (int r = 0; r < numRegions; ++r)
        opserr << " [region " << regionTags(r) << ": KE IE DW ULW RES ERR%]";
    opserr << endln;

    initializationDone = true;

    return 0;
}


int EnergyBalanceRecorder::flush(void) {
    if (theOutputHandler != 0)
        return theOutputHandler->flush();
    return 0;
}
