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

// Ladruno: LadrunoDynamicRelaxation — fork-authored TransientIntegrator
// (classTag 33005). Quasi-static solver by EXPLICIT, MATRIX-FREE dynamic
// relaxation: a fictitious-mass leap-frog march with kinetic (Cundall) damping
// driven to STATIC equilibrium. Runs under the stock DirectIntegrationAnalysis
// driver with `algorithm Linear` + `system Diagonal` (exactly as
// CentralDifferenceLadruno does) — NO analysis-core surgery.
//
//   M* a_n = f_ext - f_int(u_n) ;  v_{n+1/2} = v_{n-1/2} + dt a_n ;
//   u_{n+1} = u_n + dt v_{n+1/2}
//   kinetic damping: zero ALL velocities at each kinetic-energy peak.
//
// M* is the integrator's OWN fictitious lumped mass (default Gershgorin row-sum
// of the element stiffness, scale-free / critically-fast), NOT Element::getMass()
// (which is consistent-by-default + density-scaled + zero on zero-density models
// — the ADR-20 §2.5 BLOCKER). It is assembled onto the diagonal SOE by an
// overridden formTangent(), so theSOE->solve() degenerates to a M*^{-1} apply and
// K is NEVER assembled/factorized in the stepping loop.
//
// Termination is script-owned (the stock transient driver has no convergence
// early-exit): a Ladruno_scripts relax loop runs analyze(chunk,dt) and watches
// the velocity/displacement decay (Ladruno_scripts/relax_to_static.py).
//
// See Ladruno_implementation/21_ladruno_dynamic_relaxation_adr.md.

#ifndef LadrunoDynamicRelaxation_h
#define LadrunoDynamicRelaxation_h

#include <TransientIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;

class LadrunoDynamicRelaxation : public TransientIntegrator
{
  public:
    // massMode: 0 = gershgorin (default), 1 = lumped (real mass row-sum * scale),
    //           2 = unity
    // All args defaulted, so this doubles as the broker default ctor
    // (new LadrunoDynamicRelaxation() -> all defaults, then recvSelf).
    LadrunoDynamicRelaxation(int massMode = 0, double dtPseudo = 1.0,
                             double massScale = 1.0, int recomputeEvery = 0,
                             bool interp = false, double divergenceFactor = 0.0,
                             bool verbose = false,
                             int dampMode = 0, double zetaTarget = 1.0,
                             bool autoRefresh = true);
    ~LadrunoDynamicRelaxation();

    // TransientIntegrator contract
    int newStep(double deltaT);
    int update(const Vector &deltaU);
    int commit(void);
    int domainChanged(void);

    // matrix-free LHS: assemble the cached fictitious M* onto the diagonal SOE
    int formTangent(int statusFlag = CURRENT_TANGENT);
    int formEleTangent(FE_Element *theEle);   // unused (formTangent overridden)
    int formNodTangent(DOF_Group *theDof);    // unused (formTangent overridden)

    const Vector &getVel(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    int buildFictitiousMass(void);   // populate *Mstar in domainChanged()

    // --- options ---
    int    massMode;                 // 0 gershgorin | 1 lumped | 2 unity
    double dtPseudo;                 // fictitious pseudo-time step (default 1)
    double massScale;               // scale for -mass lumped
    int    recomputeEvery;           // refresh M* every N steps (0 = never)
    bool   interp;                   // parabolic KE-peak interpolation (reserved)
    double divergenceFactor;         // KE runaway breaker (0 = off)
    bool   verbose;

    // --- fictitious mass + leap-frog state (CDL skeleton) ---
    Vector *Mstar;                   // integrator-owned diagonal fictitious mass
    Vector *Ut, *Vhalf, *Aprev, *Vfull, *Azero;
    bool   firstStep;
    int    updateCount;              // one solve/step guard
    int    size;
    int    stepCount;                // -recompute cadence

    // --- kinetic damping / termination bookkeeping ---
    double prevKE;                   // last kinetic energy (peak detection)
    double kineticEnergy;            // current 1/2 v^T M* v
    double residualNorm;             // current ||f_ext - f_int|| = ||M* . a||_inf
    double deltaT;

    // --- v2: damping mode + auto-refresh ---
    int    dampMode;                 // 0 = kinetic/Cundall (default) | 1 = viscous-critical
    double zetaTarget;               // viscous damping ratio (default 1.0 = critical)
    double cVisc;                    // realized mass-proportional coeff (C* = cVisc*M*)
    bool   autoRefresh;              // rebuild M* at KE peaks (gershgorin) — knob-free snap-through
};

#endif
