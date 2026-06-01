/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

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

#include <string.h>
#include <stdio.h>

static int respToDataFlag(const char *resp)
{
    if (resp == 0 || strcmp(resp, "disp") == 0) return 0;
    if (strcmp(resp, "vel") == 0) return 1;
    if (strcmp(resp, "accel") == 0) return 2;
    if (strcmp(resp, "reaction") == 0) return 7;
    return -1; // unsupported in v1
}

LadrunoMonitorRecorder::LadrunoMonitorRecorder()
    : Recorder(RECORDER_TAGS_LadrunoMonitorRecorder),
      theDofs(0), theNodalTags(0), dataFlag(0), respName("disp"),
      filename(), theDomain(0),
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
      dataFlag(respToDataFlag(resp)),
      respName(resp ? resp : "disp"),
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
    if (dataFlag < 0) {
        opserr << "LadrunoMonitorRecorder: unsupported -resp '" << respName.c_str()
               << "' (v1 supports disp|vel|accel|reaction)\n";
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
    theNodes.clear();
    char buf[128];
    for (int i = 0; i < theNodalTags.Size(); ++i) {
        int tag = theNodalTags(i);
        Node *nd = theDomain->getNode(tag);
        if (nd == 0) {
            opserr << "LadrunoMonitorRecorder: node " << tag << " not found, skipped\n";
            continue;
        }
        theNodes.push_back(nd);
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

    if (dataFlag == 7)
        theDomain->calculateNodalReactions(0);

    // Gather selected scalars in channel order (node-major, dof-minor).
    size_t k = 0;
    const int numDOF = theDofs.Size();
    for (size_t n = 0; n < theNodes.size(); ++n) {
        Node *nd = theNodes[n];
        const Vector *r = 0;
        switch (dataFlag) {
        case 0: r = &nd->getTrialDisp();  break;
        case 1: r = &nd->getTrialVel();   break;
        case 2: r = &nd->getTrialAccel(); break;
        case 7: r = &nd->getReaction();   break;
        default: break;
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

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();

        if (strcmp(opt, "-node") == 0 || strcmp(opt, "-nodes") == 0) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, v;
                if (OPS_GetIntInput(&num, &v) < 0) { OPS_ResetCurrentInputArg(-1); break; }
                nodes.push_back(v);
            }
        }
        else if (strcmp(opt, "-region") == 0) {
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
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1, d;
                if (OPS_GetIntInput(&num, &d) < 0) { OPS_ResetCurrentInputArg(-1); break; }
                dofs.push_back(d - 1);
            }
        }
        else if (strcmp(opt, "-resp") == 0 || strcmp(opt, "-response") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0)
                resp = OPS_GetString();
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
