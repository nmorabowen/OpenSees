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

// ADR-30 P1 — LadrunoProjectionHandler (HANDLER_TAG 33001, first fork handler).
//
// Plain-style assembly that does NOT eliminate MP slave DOFs: every node DOF keeps
// its own equation and its own diagonal mass (homogeneous SPs excluded, as Plain).
// The multi-point constraints are instead enforced in explicit dynamics by the
// LadrunoConstraintProjector (M-orthogonal acceleration projection), pushed to a
// projection-aware integrator (CentralDifferenceLadruno) via LadrunoProjectionConsumer.
//
// At handle() it builds the connected-component constraint groups (union-find over
// DOF vertices), runs the topological diagnostics (chain / double-slave / zero-free)
// and the lumped-mass guard, and stores the groups; at doneNumberingDOF() it freezes
// equation indices into the projector and pushes it to the integrator.
//
// v1 scope: equalDOF / equalDOF-style identity MPs. rigidLink/rigidDiaphragm transport
// (general C) = P2. Command: constraints LadrunoProjection <-verbose> <-projectICs>
// <-icTol $tol>. See Ladruno_implementation/30_..._adr.md and _adr30_p1_design.md.

#ifndef LadrunoProjectionHandler_h
#define LadrunoProjectionHandler_h

#include <ConstraintHandler.h>
#include <vector>
#include <set>
#include <utility>

class LadrunoConstraintProjector;
class FE_Element;
class DOF_Group;
class Node;
class SP_Constraint;

class LadrunoProjectionHandler : public ConstraintHandler
{
  public:
    LadrunoProjectionHandler();
    LadrunoProjectionHandler(bool verbose, bool projectICs, double icTol);
    ~LadrunoProjectionHandler();

    int handle(const ID *nodesNumberedLast = 0);
    int doneNumberingDOF(void);
    int applyLoad(void);          // P4b: enforce prescribed-motion (non-homog SP) displacement
    void clearAll(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // P3 tie-force query: the constraint force f_tie = M(a_raw - a_proj) at (nodeTag, dof)
    // from the last projection step. 0 if no projector, unknown node/dof, or untied DOF.
    double getTieForce(int nodeTag, int dof) const;

  protected:

  private:
    // --- per-slave record (one constrained DOF and the retained masters it ties to)
    struct SlaveRec {
        int vtx;                                 // vertex id of the constrained DOF
        std::vector<int> masterVtx;              // retained vertices it ties to
        std::vector<double> coeff;               // C row entries (= 1 for equalDOF)
        double delta;                            // Uc0(i) - sum_j C(i,j) Ur0(j)
    };
    // --- a connected-component constraint group (filled at handle())
    struct RawGroup {
        std::vector<int> masterVtx;              // retained DOF vertices (cols of L)
        std::vector<SlaveRec> slaves;            // constrained rows of L
    };

    // --- P4b prescribed-motion record (one non-homogeneous SP / imposedMotion DOF) ---
    // The DOF stays SP-excluded (eqn = -1, like a fix); applyLoad() imposes its prescribed
    // displacement on the node each step (mirrors TransformationDOF_Group::enforceSPs).
    // ImposedMotionSP sets vel/accel itself in domain applyLoad (before handler applyLoad);
    // a plain SP supplies only the displacement (vel/accel = 0 — same as Transformation).
    struct PrescribedDOF {
        int nodeTag;                             // resolved fresh via getNode() in applyLoad()
        int dof;                                 // local dof number
        SP_Constraint *sp;                       // for getValue()+getInitialValue()
    };
    std::vector<PrescribedDOF> prescribedDOFs;
    std::set<std::pair<int,int> > prescribedKey; // (nodeTag,dof) of every prescribed DOF

    // --- P4c derived-prescribed slave (Tier 2: prescribed MASTER) ---------------------
    // A slave DOF tied ONLY to prescribed master(s) (no free retained master) is itself
    // fully determined by the prescribed motion: u_c = sum_k C_k u_{p_k} + delta, and
    // likewise v_c / a_c (no delta). Rather than the (acceleration-only, drifting) known-RHS
    // projection, such a slave is KINEMATICALLY imposed: excluded from the equation set
    // (eqn = -1, like its masters) and its u/v/a set from the masters each step in applyLoad
    // (after the masters' own disp is imposed) -> exact tracking, no drift, like Transformation.
    // A slave tied to BOTH a free and a prescribed master (mixed) is refused (handle()).
    struct DerivedSlave {
        int nodeTag, dof;                        // the constrained DOF (eqn-excluded)
        std::vector<int> masterNode, masterDof;  // its prescribed master DOFs
        std::vector<double> coeff;               // C row entries on those masters
        double delta;                            // Uc0 - sum_k C_k Ur0_k (disp offset only)
    };
    std::vector<DerivedSlave> derivedSlaves;
    std::set<std::pair<int,int> > derivedKey;    // (nodeTag,dof) of every derived-prescribed slave

    // vertex bookkeeping (a vertex = one (node,dof) pair)
    std::vector<int> vtxNode;                    // node tag per vertex
    std::vector<int> vtxDof;                     // local dof per vertex
    std::vector<RawGroup> theGroups;

    LadrunoConstraintProjector *theProjector;    // owned
    bool verbose;
    bool projectICs;
    double icTol;

    int buildGroups(void);                       // handle()-time grouping + diagnostics
    int classifyDerivedSlaves(int &countDOF);    // P4c: exclude prescribed-driven slaves, refuse mixed
    int consistentMassGuard(void);               // refuse cMass on a tied DOF
};

#endif
