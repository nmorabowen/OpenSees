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
      internal_energy(0.), damping_work(0.), unbalanced_load_work(0.),
      prev_internal_rate(0.), prev_damping_rate(0.), prev_unbalanced_rate(0.),
      eref_global(0.),
      regionTags(), numRegions(0),
      region_internal_energy(), region_damping_work(), region_unbalanced_load_work(),
      prev_region_internal_rate(), prev_region_damping_rate(), prev_region_unbalanced_rate(),
      region_eref(),
      step_region_ke(), step_region_internal_rate(),
      step_region_damping_rate(), step_region_unbalanced_rate(),
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
      internal_energy(0.), damping_work(0.), unbalanced_load_work(0.),
      prev_internal_rate(0.), prev_damping_rate(0.), prev_unbalanced_rate(0.),
      eref_global(0.),
      regionTags(regions), numRegions(regions.Size()),
      region_internal_energy(regions.Size()), region_damping_work(regions.Size()),
      region_unbalanced_load_work(regions.Size()),
      prev_region_internal_rate(regions.Size()), prev_region_damping_rate(regions.Size()),
      prev_region_unbalanced_rate(regions.Size()), region_eref(regions.Size()),
      step_region_ke(regions.Size()), step_region_internal_rate(regions.Size()),
      step_region_damping_rate(regions.Size()), step_region_unbalanced_rate(regions.Size()),
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


namespace {

// Accumulate the energy contributions of a single element. vel is a caller-
// supplied scratch vector (>= numDOF) so the hot loop performs no allocation.
//   kinetic       += 1/2 v^T M v   (instantaneous kinetic energy)
//   internal_rate += F^T v         (internal power; caller integrates over dt)
//   damping_rate  += v^T C v       (damping power;  caller integrates over dt)
static void addElementEnergy(Element *ele, Vector &vel,
                             double &kinetic,
                             double &internal_rate,
                             double &damping_rate)
{
    const Matrix &M = ele->getMass();
    const Matrix &C = ele->getDamp();
    const Vector &F = ele->getResistingForce();
    const int numExternalNodes = ele->getNumExternalNodes();
    const int numDOF = ele->getNumDOF();
    Node **elenodes = ele->getNodePtrs();

    int cnt = 0;
    for (int i = 0; i < numExternalNodes && cnt < numDOF; ++i) {
        const Vector &node_vel = elenodes[i]->getVel();
        for (int j = 0; j < node_vel.Size() && cnt < numDOF; ++j)
            vel(cnt++) = node_vel(j);
    }
    for (; cnt < numDOF; ++cnt)
        vel(cnt) = 0.0;

    // Some elements return an empty (0x0) mass or damping matrix when they
    // carry none; only accumulate when the matrix is the expected size.
    const bool haveMass = (M.noRows() == numDOF && M.noCols() == numDOF);
    const bool haveDamp = (C.noRows() == numDOF && C.noCols() == numDOF);
    if (haveMass)
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);
    if (haveDamp)
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                damping_rate += C(i, j) * vel(i) * vel(j);

    if (F.Size() == numDOF) {
        double dot = 0.0;
        for (int i = 0; i < numDOF; ++i)
            dot += F(i) * vel(i);
        internal_rate += dot;
    }
}

// Accumulate the nodal contributions: lumped-mass kinetic energy, nodal
// Rayleigh damping (alphaM * nodalMass), and the work rate of the external
// applied nodal load. M is fully consumed before getDamp() is called because
// Node::getMass()/getDamp() can share the same scratch matrix when massless.
static void addNodeEnergy(Node *node,
                          double &kinetic,
                          double &damping_rate,
                          double &unbalanced_rate)
{
    const Vector &vel = node->getVel();
    const int n = vel.Size();

    const Matrix &M = node->getMass();
    if (M.noRows() == n && M.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);

    const Matrix &C = node->getDamp();
    if (C.noRows() == n && C.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                damping_rate += C(i, j) * vel(i) * vel(j);

    const Vector &Pext = node->getUnbalancedLoad();
    if (Pext.Size() == n)
        unbalanced_rate += vel ^ Pext;
}

} // namespace


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
    // single sweep: compute each element/node contribution once, accumulate
    // into the whole-model total and bin into the regions it belongs to.
    //
    double g_ke = 0., g_internal_rate = 0., g_damping_rate = 0., g_unbalanced_rate = 0.;
    step_region_ke.Zero();
    step_region_internal_rate.Zero();
    step_region_damping_rate.Zero();
    step_region_unbalanced_rate.Zero();

    {
        Element *ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0) {
            double ke = 0., ie = 0., dw = 0.;
            addElementEnergy(ele, velScratch, ke, ie, dw);
            g_ke += ke; g_internal_rate += ie; g_damping_rate += dw;

            if (numRegions > 0) {
                std::unordered_map<Element*, std::vector<int> >::const_iterator it =
                    elemRegions.find(ele);
                if (it != elemRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke(r) += ke;
                        step_region_internal_rate(r) += ie;
                        step_region_damping_rate(r) += dw;
                    }
            }
        }
    }

    {
        Node *node;
        NodeIter &nodes = theDomain->getNodes();
        while ((node = nodes()) != 0) {
            double ke = 0., dw = 0., ulw = 0.;
            addNodeEnergy(node, ke, dw, ulw);
            g_ke += ke; g_damping_rate += dw; g_unbalanced_rate += ulw;

            if (numRegions > 0) {
                std::unordered_map<Node*, std::vector<int> >::const_iterator it =
                    nodeRegions.find(node);
                if (it != nodeRegions.end())
                    for (size_t k = 0; k < it->second.size(); ++k) {
                        const int r = it->second[k];
                        step_region_ke(r) += ke;
                        step_region_damping_rate(r) += dw;
                        step_region_unbalanced_rate(r) += ulw;
                    }
            }
        }
    }

    //
    // integrate the work rates (trapezoidal). The first record() only seeds
    // the previous rate so no increment is taken over the initial time jump.
    //
    if (firstRecord) {
        // first recorded step: integrate with the rectangle rule (there is no
        // previous-step rate yet). Skipping it would drop one increment of work
        // and leave a constant residual offset.
        internal_energy      += g_internal_rate   * dT;
        damping_work         += g_damping_rate    * dT;
        unbalanced_load_work += g_unbalanced_rate * dT;
        for (int r = 0; r < numRegions; ++r) {
            region_internal_energy(r)      += step_region_internal_rate(r)   * dT;
            region_damping_work(r)         += step_region_damping_rate(r)    * dT;
            region_unbalanced_load_work(r) += step_region_unbalanced_rate(r) * dT;
        }
    } else {
        // trapezoidal rule
        internal_energy      += 0.5 * (prev_internal_rate   + g_internal_rate)   * dT;
        damping_work         += 0.5 * (prev_damping_rate    + g_damping_rate)    * dT;
        unbalanced_load_work += 0.5 * (prev_unbalanced_rate + g_unbalanced_rate) * dT;
        for (int r = 0; r < numRegions; ++r) {
            region_internal_energy(r)      += 0.5 * (prev_region_internal_rate(r)   + step_region_internal_rate(r))   * dT;
            region_damping_work(r)         += 0.5 * (prev_region_damping_rate(r)    + step_region_damping_rate(r))    * dT;
            region_unbalanced_load_work(r) += 0.5 * (prev_region_unbalanced_rate(r) + step_region_unbalanced_rate(r)) * dT;
        }
    }
    firstRecord = false;

    prev_internal_rate   = g_internal_rate;
    prev_damping_rate    = g_damping_rate;
    prev_unbalanced_rate = g_unbalanced_rate;
    for (int r = 0; r < numRegions; ++r) {
        prev_region_internal_rate(r)   = step_region_internal_rate(r);
        prev_region_damping_rate(r)    = step_region_damping_rate(r);
        prev_region_unbalanced_rate(r) = step_region_unbalanced_rate(r);
    }

    //
    // fill the response: whole-model block then one block per region.
    // RES = ULW - (KE + IE + DW); ERR% = 100 * RES / running-max-total-energy.
    //
    {
        const double sum = g_ke + internal_energy + damping_work;
        const double res = unbalanced_load_work - sum;
        // energy scale = running max of summed component magnitudes (does not
        // collapse to ~0 when KE and incremental IE cancel, e.g. free vibration)
        const double mag = fabs(g_ke) + fabs(internal_energy) + fabs(damping_work)
                         + fabs(unbalanced_load_work);
        if (mag > eref_global) eref_global = mag;
        const double err = (eref_global > 1.0e-16) ? 100.0 * res / eref_global : 0.0;

        int col = timeOffset;
        response(col + 0) = g_ke;
        response(col + 1) = internal_energy;
        response(col + 2) = damping_work;
        response(col + 3) = unbalanced_load_work;
        response(col + 4) = res;
        response(col + 5) = err;
    }

    for (int r = 0; r < numRegions; ++r) {
        const double ke  = step_region_ke(r);
        const double ie  = region_internal_energy(r);
        const double dw  = region_damping_work(r);
        const double ulw = region_unbalanced_load_work(r);
        const double sum = ke + ie + dw;
        const double res = ulw - sum;
        const double mag = fabs(ke) + fabs(ie) + fabs(dw) + fabs(ulw);
        if (mag > region_eref(r)) region_eref(r) = mag;
        const double err = (region_eref(r) > 1.0e-16) ? 100.0 * res / region_eref(r) : 0.0;

        const int col = timeOffset + EBR_NUM_ENERGY_COMPONENTS * (1 + r);
        response(col + 0) = ke;
        response(col + 1) = ie;
        response(col + 2) = dw;
        response(col + 3) = ulw;
        response(col + 4) = res;
        response(col + 5) = err;
    }

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

    // size and zero the response and all per-region accumulators
    const int timeOffset = (echoTimeFlag == true) ? 1 : 0;
    response.resize(timeOffset + EBR_NUM_ENERGY_COMPONENTS * (1 + numRegions));
    response.Zero();

    region_internal_energy.resize(numRegions);
    region_damping_work.resize(numRegions);
    region_unbalanced_load_work.resize(numRegions);
    prev_region_internal_rate.resize(numRegions);
    prev_region_damping_rate.resize(numRegions);
    prev_region_unbalanced_rate.resize(numRegions);
    region_eref.resize(numRegions);
    step_region_ke.resize(numRegions);
    step_region_internal_rate.resize(numRegions);
    step_region_damping_rate.resize(numRegions);
    step_region_unbalanced_rate.resize(numRegions);

    region_internal_energy.Zero();
    region_damping_work.Zero();
    region_unbalanced_load_work.Zero();
    prev_region_internal_rate.Zero();
    prev_region_damping_rate.Zero();
    prev_region_unbalanced_rate.Zero();
    region_eref.Zero();

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
