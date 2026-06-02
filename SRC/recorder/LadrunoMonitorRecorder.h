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

// Ladruno: live analysis-monitor recorder.
//
// Written for the Ladruno fork (N. Mora-Bowen). See
// Ladruno_implementation/08_analysis_monitor.md.
//
// A LadrunoMonitorRecorder is a normal Recorder added to the domain, so its
// record() fires at every commit. Instead of writing an at-rest results file,
// it streams the selected nodal scalars (disp/vel/accel/reaction) to a
// SWMR-HDF5 sink that a separate viewer process tails live. Two gates bound
// the data rate at explicit step counts so streaming can never perturb the
// solve (resolved decision #2):
//   * -every K : emit only every Kth committed step,
//   * -hz   H  : emit at most H frames/second (steady_clock wall throttle).
//
// v1 scope: nodal scalar time-histories only, sequential runs. Element /
// region-response channels and the parallel path are follow-ups.

#ifndef LadrunoMonitorRecorder_h
#define LadrunoMonitorRecorder_h

#include <Recorder.h>
#include <ID.h>

#include <chrono>
#include <string>
#include <vector>

class Domain;
class Node;

namespace ladruno { class MonitorSink; }

class LadrunoMonitorRecorder : public Recorder {
public:
    LadrunoMonitorRecorder();
    LadrunoMonitorRecorder(const ID &theDofs,
                           const ID &theNodalTags,
                           const char *respName,
                           Domain &theDomain,
                           const char *filename,
                           int everyK,
                           double hz);
    ~LadrunoMonitorRecorder();

    int record(int commitTag, double timeStamp);
    int flush();

    int restart(void);
    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:
    int initialize(void);

    ID theDofs;            // 0-based dof indices
    ID theNodalTags;
    std::vector<Node *> theNodes;   // resolved valid nodes (parallel to labels)
    int dataFlag;          // 0 disp, 1 vel, 2 accel, 7 reaction
    std::string respName;
    std::string filename;
    Domain *theDomain;

    // rate gates
    int everyK;                 // >=1
    double minIntervalSec;      // from -hz; <=0 => no wall throttle
    long long recordCallCount;  // counts every record() call (for -every)
    std::chrono::steady_clock::time_point lastEmit;
    bool haveEmitted;

    ladruno::MonitorSink *sink; // owned
    bool initialized;
    bool initFailed;
    std::vector<double> frameBuf;
};

#endif
