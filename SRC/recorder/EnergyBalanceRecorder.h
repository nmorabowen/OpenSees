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



#include <Recorder.h>
#include <ID.h>
#include <Vector.h>
#include <TimeSeries.h>

#include <vector>
#include <unordered_map>

class Domain;
class FE_Datastore;
class Node;
class Element;
class MeshRegion;

#define EBR_NUM_ENERGY_COMPONENTS 6

class EnergyBalanceRecorder: public Recorder
{
public:
    EnergyBalanceRecorder();
    EnergyBalanceRecorder(
        Domain &theDomain,
        OPS_Stream &theOutputHandler,
        bool echoTimeFlag = true,
        const ID &regions = ID());

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

    // ---- whole-model accumulators (work integrals) and closure helpers ----
    double internal_energy;
    double damping_work;
    double unbalanced_load_work;
    double prev_internal_rate;     // previous-step rates, for trapezoidal rule
    double prev_damping_rate;
    double prev_unbalanced_rate;
    double eref_global;            // running max total energy, for ERR%

    // ---- per-region accumulators (mirrors of the whole-model ones) ----
    ID regionTags;
    int numRegions;
    Vector region_internal_energy;
    Vector region_damping_work;
    Vector region_unbalanced_load_work;
    Vector prev_region_internal_rate;
    Vector prev_region_damping_rate;
    Vector prev_region_unbalanced_rate;
    Vector region_eref;
    // per-step scratch (sized numRegions), members to avoid per-step alloc
    Vector step_region_ke;
    Vector step_region_internal_rate;
    Vector step_region_damping_rate;
    Vector step_region_unbalanced_rate;

    // ---- caching (rebuilt on domainChanged): region membership + scratch --
    bool cacheValid;
    std::unordered_map<Element*, std::vector<int> > elemRegions;
    std::unordered_map<Node*, std::vector<int> > nodeRegions;
    int maxNumDOF;
    Vector velScratch;
};

#endif
