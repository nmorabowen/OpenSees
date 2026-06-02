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
#include <MeshRegion.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <TCP_Stream.h>

#include <elementAPI.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

void*
OPS_EnergyBalanceRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: recorder EnergyBalance -file <fileName> <-time> "
                  "<-region tag ...>\n";
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

    const char *inetAddr = 0;
    int inetPort = 0;

    ID regions(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
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
        *domain, *theOutputStream, echoTimeFlag, regions);

    return recorder;
}


EnergyBalanceRecorder::EnergyBalanceRecorder()
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(0), theOutputHandler(0),
      response(EBR_NUM_ENERGY_COMPONENTS),
      echoTimeFlag(true), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      regionTags(), numRegions(0),
      region_acc(),
      cacheValid(false), maxNumDOF(0), velScratch()
{

}


EnergyBalanceRecorder::EnergyBalanceRecorder(
    Domain &theDomain_,
    OPS_Stream &theOutputHandler_,
    bool echoTimeFlag_,
    const ID &regions)
    : Recorder(RECORDER_TAGS_EnergyBalanceRecorder),
      theDomain(&theDomain_), theOutputHandler(&theOutputHandler_),
      response((echoTimeFlag_ ? 1 : 0) + EBR_NUM_ENERGY_COMPONENTS * (1 + regions.Size())),
      echoTimeFlag(echoTimeFlag_), initializationDone(false), firstRecord(true),
      time_last(0.),
      model_acc(),
      regionTags(regions), numRegions(regions.Size()),
      region_acc((size_t)regions.Size()),
      cacheValid(false), maxNumDOF(0), velScratch()
{
    response.Zero();
}


EnergyBalanceRecorder::~EnergyBalanceRecorder()
{
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
            Element *ele = theDomain->getElement(elems(i));
            if (ele != 0)
                elemRegions[ele].push_back(r);
        }
        const ID &rnodes = region->getNodes();
        for (int i = 0; i < rnodes.Size(); ++i) {
            Node *node = theDomain->getNode(rnodes(i));
            if (node != 0)
                nodeRegions[node].push_back(r);
        }
    }

    // size the velocity scratch to the largest element DOF count
    maxNumDOF = 0;
    Element *ele;
    ElementIter &elements = theDomain->getElements();
    while ((ele = elements()) != 0) {
        const int n = ele->getNumDOF();
        if (n > maxNumDOF)
            maxNumDOF = n;
    }
    if (velScratch.Size() < maxNumDOF)
        velScratch.resize(maxNumDOF);

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
    std::vector<double> step_region_ke((size_t)numRegions, 0.0);
    std::vector<double> step_region_internal_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_damping_rate((size_t)numRegions, 0.0);
    std::vector<double> step_region_unbalanced_rate((size_t)numRegions, 0.0);

    {
        Element *ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0) {
            double ke = 0., ie = 0., dw = 0.;
            ebkernel::addElementEnergy(ele, velScratch, ke, ie, dw);
            g_ke += ke; g_internal_rate += ie; g_damping_rate += dw;

            if (numRegions > 0) {
                std::unordered_map<Element*, std::vector<int> >::const_iterator it =
                    elemRegions.find(ele);
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
        Node *node;
        NodeIter &nodes = theDomain->getNodes();
        while ((node = nodes()) != 0) {
            double ke = 0., dw = 0., ulw = 0.;
            ebkernel::addNodeEnergy(node, ke, dw, ulw);
            g_ke += ke; g_damping_rate += dw; g_unbalanced_rate += ulw;

            if (numRegions > 0) {
                std::unordered_map<Node*, std::vector<int> >::const_iterator it =
                    nodeRegions.find(node);
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
    {
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
        const int col = timeOffset + EBR_NUM_ENERGY_COMPONENTS * (1 + r);
        for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
            response(col + c) = out[c];
    }

    firstRecord = false;

    theOutputHandler->write(response);

    time_last = timeStamp;

    return 0;
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

    // size and zero the response and all per-scope accumulators
    const int timeOffset = (echoTimeFlag == true) ? 1 : 0;
    response.resize(timeOffset + EBR_NUM_ENERGY_COMPONENTS * (1 + numRegions));
    response.Zero();

    model_acc.reset();
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
    for (int c = 0; c < EBR_NUM_ENERGY_COMPONENTS; ++c)
        theOutputHandler->tag("ResponseType", comp[c]);
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
    opserr << " [model: KE IE DW ULW RES ERR%]";
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
