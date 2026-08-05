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

// Ladruno: live analysis-monitor recorder. See LadrunoMonitorRecorder.h and
// Ladruno_implementation/08_analysis_monitor.md.

#include "LadrunoMonitorRecorder.h"
#include "LadrunoMonitorSink.h"

#include <Domain.h>
#include <Node.h>
#include <MeshRegion.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>
#include <LoadPattern.h>              // Ladruno (ADR-73 P4): -overlay mode
#include <LadrunoPorousOverlay.h>

#include <vector>

#include <string.h>
#include <stdio.h>

static int respToDataFlag(const char *resp)
{
    if (resp == 0 || strcmp(resp, "disp") == 0) return 0;
    if (strcmp(resp, "vel") == 0) return 1;
    if (strcmp(resp, "accel") == 0) return 2;
    if (strcmp(resp, "reaction") == 0) return 7;
    // reaction including inertia (mass*accel); + Rayleigh damping. Names + flag
    // meaning mirror NodeRecorder (dataFlag 8 -> calculateNodalReactions(1),
    // dataFlag 9 -> (2)) so a dynamic run streams true reactions, not the
    // static part only.
    if (strcmp(resp, "reactionIncInertia") == 0 ||
        strcmp(resp, "reactionInclInertia") == 0 ||
        strcmp(resp, "reactionIncludingInertia") == 0) return 8;
    if (strcmp(resp, "rayleighForces") == 0 ||
        strcmp(resp, "rayleighDampingForces") == 0) return 9;
    return -1; // unsupported in v1
}

LadrunoMonitorRecorder::LadrunoMonitorRecorder()
    : Recorder(RECORDER_TAGS_LadrunoMonitorRecorder),
      theDofs(0), theNodalTags(0), dataFlag(0), overlayMode(false), overlayTag(-1),
      respName("disp"), filename(), theDomain(0),
      everyK(1), minIntervalSec(0.0), recordCallCount(0),
      haveEmitted(false), sink(0), initialized(false), initFailed(false)
{
}

LadrunoMonitorRecorder::LadrunoMonitorRecorder(const ID &dofs,
                                               const ID &nodalTags,
                                               const char *resp,
                                               Domain &dom,
                                               const char *fname,
                                               int every,
                                               double hz)
    : Recorder(RECORDER_TAGS_LadrunoMonitorRecorder),
      theDofs(dofs), theNodalTags(nodalTags),
      dataFlag(respToDataFlag(resp)), overlayMode(false), overlayTag(-1),
      respName(resp ? resp : "disp"),
      filename(fname ? fname : ""), theDomain(&dom),
      everyK(every < 1 ? 1 : every),
      minIntervalSec(hz > 0.0 ? 1.0 / hz : 0.0),
      recordCallCount(0), haveEmitted(false),
      sink(0), initialized(false), initFailed(false)
{
}

// Overlay-pressure mode (ADR-73 P4). dataFlag is unused in this mode (p is a
// nodal scalar); overlayMode selects the branch in initialize()/record().
LadrunoMonitorRecorder::LadrunoMonitorRecorder(int ovTag,
                                               const ID &nodalSubsetTags,
                                               Domain &dom,
                                               const char *fname,
                                               int every,
                                               double hz)
    : Recorder(RECORDER_TAGS_LadrunoMonitorRecorder),
      theDofs(0), theNodalTags(nodalSubsetTags),
      dataFlag(0), overlayMode(true), overlayTag(ovTag), respName("p"),
      filename(fname ? fname : ""), theDomain(&dom),
      everyK(every < 1 ? 1 : every),
      minIntervalSec(hz > 0.0 ? 1.0 / hz : 0.0),
      recordCallCount(0), haveEmitted(false),
      sink(0), initialized(false), initFailed(false)
{
}

LadrunoMonitorRecorder::~LadrunoMonitorRecorder()
{
    if (sink != 0) {
        sink->close();
        delete sink;
        sink = 0;
    }
}

int LadrunoMonitorRecorder::initialize(void)
{
    initialized = true; // attempt once, regardless of outcome

    if (theDomain == 0) {
        opserr << "LadrunoMonitorRecorder: no domain set\n";
        initFailed = true;
        return -1;
    }

    // -------- overlay pressure mode (ADR-73 P4) --------
    // Columns overlay<tag>.p.node<n>, one scalar per region node (subset or all).
    // The subset was validated against the region nodes at parse time; here we
    // just resolve the overlay and freeze the SWMR columns.
    if (overlayMode) {
        if (filename.empty()) {
            opserr << "LadrunoMonitorRecorder: no -sink file given\n";
            initFailed = true;
            return -1;
        }
        LoadPattern *lp = theDomain->getLoadPattern(overlayTag);
        if (lp == 0 || lp->getClassTag() != PATTERN_TAG_LadrunoPorousOverlay) {
            opserr << "LadrunoMonitorRecorder: -overlay " << overlayTag
                   << " is not a LadrunoPorousOverlay in the domain\n";
            initFailed = true;
            return -1;
        }
        LadrunoPorousOverlay *ov = (LadrunoPorousOverlay *)lp;
        if (!ov->snapshotReady()) {
            opserr << "LadrunoMonitorRecorder: overlay " << overlayTag
                   << " has no completed geometry snapshot; cannot stream p "
                      "(build the model/pattern before the recorder)\n";
            initFailed = true;
            return -1;
        }

        std::vector<std::string> columns;
        theNodeTags.clear();
        char buf[128];
        if (theNodalTags.Size() > 0) {
            for (int i = 0; i < theNodalTags.Size(); ++i) {
                int tag = theNodalTags(i);
                theNodeTags.push_back(tag);
                snprintf(buf, sizeof(buf), "overlay%d.p.node%d", overlayTag, tag);
                columns.push_back(buf);
            }
        } else {
            const std::vector<int> &rn = ov->getRegionNodeTags();
            for (size_t i = 0; i < rn.size(); ++i) {
                theNodeTags.push_back(rn[i]);
                snprintf(buf, sizeof(buf), "overlay%d.p.node%d", overlayTag, rn[i]);
                columns.push_back(buf);
            }
        }
        if (columns.empty()) {
            opserr << "LadrunoMonitorRecorder: overlay " << overlayTag
                   << " has no region nodes to stream\n";
            initFailed = true;
            return -1;
        }

        frameBuf.assign(columns.size(), 0.0);
        sink = new ladruno::MonitorSink(filename, columns);
        if (sink->open() != 0) {
            opserr << "LadrunoMonitorRecorder: could not open sink "
                   << filename.c_str() << "\n";
            delete sink;
            sink = 0;
            initFailed = true;
            return -1;
        }
        opserr << "LadrunoMonitorRecorder: streaming " << (int)columns.size()
               << " overlay-p channels to " << filename.c_str() << " (SWMR)\n";
        return 0;
    }

    if (dataFlag < 0) {
        opserr << "LadrunoMonitorRecorder: unsupported -resp '" << respName.c_str()
               << "' (supports disp|vel|accel|reaction|reactionIncInertia|rayleighForces)\n";
        initFailed = true;
        return -1;
    }
    if (filename.empty()) {
        opserr << "LadrunoMonitorRecorder: no -sink file given\n";
        initFailed = true;
        return -1;
    }

    // Resolve valid nodes and build self-describing channel labels. Channels
    // are ordered node-major, dof-minor -- the exact order record() emits.
    std::vector<std::string> columns;
    theNodeTags.clear();
    char buf[128];
    for (int i = 0; i < theNodalTags.Size(); ++i) {
        int tag = theNodalTags(i);
        Node *nd = theDomain->getNode(tag);
        if (nd == 0) {
            opserr << "LadrunoMonitorRecorder: node " << tag << " not found, skipped\n";
            continue;
        }
        theNodeTags.push_back(tag);
        for (int j = 0; j < theDofs.Size(); ++j) {
            snprintf(buf, sizeof(buf), "node%d.%s.dof%d",
                     tag, respName.c_str(), theDofs(j) + 1);
            columns.push_back(buf);
        }
    }

    if (columns.empty()) {
        opserr << "LadrunoMonitorRecorder: no valid node/dof channels\n";
        initFailed = true;
        return -1;
    }

    frameBuf.assign(columns.size(), 0.0);

    sink = new ladruno::MonitorSink(filename, columns);
    if (sink->open() != 0) {
        opserr << "LadrunoMonitorRecorder: could not open sink " << filename.c_str() << "\n";
        delete sink;
        sink = 0;
        initFailed = true;
        return -1;
    }

    opserr << "LadrunoMonitorRecorder: streaming " << (int)columns.size()
           << " channels to " << filename.c_str() << " (SWMR)\n";
    return 0;
}

int LadrunoMonitorRecorder::record(int commitTag, double timeStamp)
{
    if (!initialized)
        this->initialize();   // lazy: open the file on the first commit
    if (initFailed || sink == 0)
        return 0;             // inert: never error-spam, never stall the solve

    ++recordCallCount;

    // -every K decimation: emit on calls 1, 1+K, 1+2K, ...
    if (everyK > 1 && ((recordCallCount - 1) % everyK) != 0)
        return 0;

    // -hz wall-clock throttle (always allow the very first emitted frame).
    if (minIntervalSec > 0.0 && haveEmitted) {
        std::chrono::duration<double> dt =
            std::chrono::steady_clock::now() - lastEmit;
        if (dt.count() < minIntervalSec)
            return 0;
    }

    // Overlay pressure mode (ADR-73 P4): gather committed p at each column node
    // via the overlay pattern. Re-resolved by tag each emit (same removal-robust
    // idiom as the node path); an absent overlay streams 0.0 rather than dangling.
    if (overlayMode) {
        LoadPattern *lp = theDomain->getLoadPattern(overlayTag);
        LadrunoPorousOverlay *ov =
            (lp != 0 && lp->getClassTag() == PATTERN_TAG_LadrunoPorousOverlay)
                ? (LadrunoPorousOverlay *)lp : 0;
        for (size_t n = 0; n < theNodeTags.size(); ++n)
            frameBuf[n] = (ov != 0) ? ov->pAtNodeTag(theNodeTags[n]) : 0.0;
        sink->append(commitTag, timeStamp, frameBuf);
        lastEmit = std::chrono::steady_clock::now();
        haveEmitted = true;
        return 0;
    }

    // Reaction modes: ask the domain to (re)compute nodal reactions. flag 0 =
    // static part only; 1 = + inertia (mass*accel); 2 = + inertia + Rayleigh
    // damping. Mirrors NodeRecorder so a dynamic run streams true reactions.
    if (dataFlag == 7)
        theDomain->calculateNodalReactions(0);
    else if (dataFlag == 8)
        theDomain->calculateNodalReactions(1);
    else if (dataFlag == 9)
        theDomain->calculateNodalReactions(2);

    // Gather selected scalars in channel order (node-major, dof-minor). Nodes
    // are re-resolved by tag on every emit, so a node removed mid-run
    // (progressive collapse / `remove node`) can never leave a dangling Node*:
    // its columns stay (SWMR fixes the channel set at open) and stream 0.0 until
    // the node reappears. This replaces the old cached-pointer path whose
    // domainChanged() no-op would dereference freed memory after a removal.
    size_t k = 0;
    const int numDOF = theDofs.Size();
    for (size_t n = 0; n < theNodeTags.size(); ++n) {
        Node *nd = theDomain->getNode(theNodeTags[n]);
        const Vector *r = 0;
        if (nd != 0) {
            switch (dataFlag) {
            case 0: r = &nd->getTrialDisp();  break;
            case 1: r = &nd->getTrialVel();   break;
            case 2: r = &nd->getTrialAccel(); break;
            case 7: case 8: case 9: r = &nd->getReaction(); break;
            default: break;
            }
        }
        for (int j = 0; j < numDOF; ++j) {
            int dof = theDofs(j);
            frameBuf[k++] = (r != 0 && r->Size() > dof) ? (*r)(dof) : 0.0;
        }
    }

    sink->append(commitTag, timeStamp, frameBuf);

    lastEmit = std::chrono::steady_clock::now();
    haveEmitted = true;
    return 0;
}

int LadrunoMonitorRecorder::flush(void)
{
    return 0; // SWMR sink flushes per frame
}

int LadrunoMonitorRecorder::restart(void)
{
    return 0;
}

int LadrunoMonitorRecorder::domainChanged(void)
{
    // No-op by design: record() re-resolves every monitored node by tag on each
    // emit, so a domain change (node/element add or remove) needs no cached-state
    // refresh here. The channel set is fixed at first record() (SWMR cannot add
    // or drop channels), so a removed node's columns remain and stream 0.0.
    return 0;
}

int LadrunoMonitorRecorder::setDomain(Domain &theDom)
{
    theDomain = &theDom;
    return 0;
}

// v1 is sequential-only; the parallel path is a follow-up (08_analysis_monitor.md).
int LadrunoMonitorRecorder::sendSelf(int, Channel &)
{
    return 0;
}

int LadrunoMonitorRecorder::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

// ---------------------------------------------------------------------------
// Interpreter parser: invoked via  recorder Monitor -node ... -sink file ...
// Registered in OpenSeesOutputCommands.cpp's recordersMap as "Monitor".
// ---------------------------------------------------------------------------
void *OPS_LadrunoMonitorRecorder()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0)
        return 0;

    std::vector<int> nodes;
    std::vector<int> dofs;   // 0-based
    const char *resp = "disp";
    const char *filename = 0;
    int everyK = 1;
    double hz = 0.0;

    // -overlay (ADR-73 P4): pore-pressure mode (sawOverlay once seen).
    // Exclusive of -node/-region/-dof/-resp; -nodes is repurposed as the
    // region-node subset. The `saw*` flags detect illegal mixing (gate d).
    int overlayTag = -1;
    bool sawOverlay = false;
    bool sawNodeSingular = false;   // -node
    bool sawRegion = false;
    bool sawDof = false;
    bool sawResp = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();

        if (strcmp(opt, "-node") == 0 || strcmp(opt, "-nodes") == 0) {
            // -node is a normal-mode flag (illegal under -overlay); -nodes is the
            // overlay region-node subset AND the normal-mode node list.
            if (strcmp(opt, "-node") == 0) sawNodeSingular = true;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, v;
                int oldn = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&num, &v) < 0) {
                    // un-consume the flag only if the failed read consumed it
                    // (openseespy does; the classic Tcl elementAPI does not)
                    if (OPS_GetNumRemainingInputArgs() < oldn) OPS_ResetCurrentInputArg(-1);
                    break;
                }
                nodes.push_back(v);
            }
        }
        else if (strcmp(opt, "-region") == 0) {
            sawRegion = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, tag;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING recorder Monitor: bad -region tag\n";
                    return 0;
                }
                MeshRegion *reg = theDomain->getRegion(tag);
                if (reg == 0) {
                    opserr << "WARNING recorder Monitor: region " << tag << " not found\n";
                    return 0;
                }
                const ID &rn = reg->getNodes();
                for (int i = 0; i < rn.Size(); ++i)
                    nodes.push_back(rn(i));
            }
        }
        else if (strcmp(opt, "-dof") == 0 || strcmp(opt, "-dofs") == 0) {
            sawDof = true;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, d;
                int oldn = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&num, &d) < 0) {
                    // un-consume the flag only if the failed read consumed it
                    // (openseespy does; the classic Tcl elementAPI does not)
                    if (OPS_GetNumRemainingInputArgs() < oldn) OPS_ResetCurrentInputArg(-1);
                    break;
                }
                dofs.push_back(d - 1);
            }
        }
        else if (strcmp(opt, "-resp") == 0 || strcmp(opt, "-response") == 0) {
            sawResp = true;
            if (OPS_GetNumRemainingInputArgs() > 0)
                resp = OPS_GetString();
        }
        else if (strcmp(opt, "-overlay") == 0) {
            sawOverlay = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, tag;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING recorder Monitor: -overlay requires an "
                              "int overlay pattern tag\n";
                    return 0;
                }
                overlayTag = tag;
            }
            else {
                opserr << "WARNING recorder Monitor: -overlay requires an int "
                          "overlay pattern tag\n";
                return 0;
            }
        }
        else if (strcmp(opt, "-sink") == 0 || strcmp(opt, "-file") == 0 ||
                 strcmp(opt, "-filename") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                filename = OPS_GetString();
        }
        else if (strcmp(opt, "-every") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &everyK) < 0) everyK = 1;
            }
        }
        else if (strcmp(opt, "-hz") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetDoubleInput(&num, &hz) < 0) hz = 0.0;
            }
        }
        // unknown options are ignored (forward-compatible)
    }

    if (filename == 0) {
        opserr << "WARNING recorder Monitor: missing -sink <file>\n";
        return 0;
    }

    // -------- overlay pressure mode (ADR-73 P4) --------
    if (sawOverlay) {
        // exclusive of -node/-region/-dof/-resp (SWMR columns are frozen at open,
        // keep the modes disjoint). -nodes is allowed as the region-node subset.
        if (sawNodeSingular || sawRegion || sawDof || sawResp) {
            opserr << "WARNING recorder Monitor: -overlay is exclusive of "
                      "-node/-region/-dof/-resp; use -nodes for the region-node "
                      "subset\n";
            return 0;
        }
        LoadPattern *lp = theDomain->getLoadPattern(overlayTag);
        if (lp == 0 || lp->getClassTag() != PATTERN_TAG_LadrunoPorousOverlay) {
            opserr << "WARNING recorder Monitor: -overlay " << overlayTag
                   << " is not a LadrunoPorousOverlay in the domain\n";
            return 0;
        }
        LadrunoPorousOverlay *ov = (LadrunoPorousOverlay *)lp;
        // -nodes subset (if given) must be region nodes (fatal otherwise);
        // empty = all region nodes.
        if (!nodes.empty()) {
            if (!ov->snapshotReady()) {
                opserr << "WARNING recorder Monitor: overlay " << overlayTag
                       << " has no completed snapshot; cannot validate the "
                          "-nodes subset (issue the recorder after the "
                          "model/pattern is built)\n";
                return 0;
            }
            const std::vector<int> &rn = ov->getRegionNodeTags();
            for (size_t i = 0; i < nodes.size(); ++i) {
                bool found = false;
                for (size_t j = 0; j < rn.size(); ++j)
                    if (rn[j] == nodes[i]) { found = true; break; }
                if (!found) {
                    opserr << "WARNING recorder Monitor: -overlay -nodes tag "
                           << nodes[i] << " is not a region node of overlay "
                           << overlayTag << "\n";
                    return 0;
                }
            }
        }
        ID subsetID((int)nodes.size());
        for (int i = 0; i < (int)nodes.size(); ++i) subsetID(i) = nodes[i];
        return new LadrunoMonitorRecorder(overlayTag, subsetID, *theDomain,
                                          filename, everyK, hz);
    }

    if (nodes.empty()) {
        opserr << "WARNING recorder Monitor: no -node/-region nodes given\n";
        return 0;
    }
    if (dofs.empty()) {
        opserr << "WARNING recorder Monitor: no -dof given\n";
        return 0;
    }

    ID nodeID((int)nodes.size());
    for (int i = 0; i < (int)nodes.size(); ++i) nodeID(i) = nodes[i];
    ID dofID((int)dofs.size());
    for (int i = 0; i < (int)dofs.size(); ++i) dofID(i) = dofs[i];

    return new LadrunoMonitorRecorder(dofID, nodeID, resp, *theDomain,
                                      filename, everyK, hz);
}
