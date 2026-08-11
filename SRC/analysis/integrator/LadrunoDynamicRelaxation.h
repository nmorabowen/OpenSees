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
// STABILITY MARGIN (-massSafety f, note 83 §3). The raw Gershgorin mass
// m_i = (dt²/4)Σ_j|K_ij| bounds ω_max·dt ≤ 2 with EQUALITY, i.e. it marches
// exactly ON the central-difference stability boundary, where the amplification
// matrix is defective and round-off is amplified instead of damped — measured, a
// zero-push hold on an EXACT equilibrium walked to an 87 kN residual on a 300 kN
// problem with no error reported. `-massSafety f` divides the prefactor by f²
// (m_i = (dt²/(4f²))Σ_j|K_ij|) so ω_max·dt ≤ 2f, at the cost of f² relaxation per
// step. Default 0.5 (see the .cpp for the measured justification).
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

// Gershgorin-mass safety factor f: omega_max*dt <= 2f. The DEFAULT is deliberately
// below 1 — f = 1 is the bare stability boundary and is a measured silent-wrong-
// answer generator (note 83 §3). Named so the tests and the ledger quote ONE number.
#define LADRUNO_DR_DEFAULT_MASS_SAFETY 0.5

// Step cadence of the stability-margin probe. One extra element-tangent pass every
// N steps — the same work a `-recompute` refresh does, but WITHOUT replacing M*, so
// it measures rather than hides drift. At 500 it is well under 1% of a DR step
// budget; `-marginEvery 0` switches it off and makes `stabilityMargin` report -2
// ("not measured") rather than a stale number.
#define LADRUNO_DR_DEFAULT_MARGIN_EVERY 500

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
                             bool autoRefresh = true,
                             double massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY,
                             int marginEvery = LADRUNO_DR_DEFAULT_MARGIN_EVERY);
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

    // --- Layer-1.5 quasi-staticness / micro-burst signals (ADR-31) ---
    // The robust-solve driver reads these to decide when DR has settled to a
    // static solution. `residualNorm` is the FORCE-based unbalance
    // ||f_ext - f_int|| -- the mass-free settling criterion (R-DR-ENERGY: the
    // EnergyBalance recorder's physical-mass KE is ~0 / wrong on DR's pseudo-mass
    // models, so the settling gate must be force-based, not a KE ratio).
    // `kineticEnergy` is the pseudo-mass KE for the micro-burst classifier.
    // Read-only on existing in-memory state (no sendSelf change, serial v1).
    double getResidualNorm(void) const;     // ||f_ext - f_int|| (settling gate)
    double getKineticEnergy(void) const;    // 1/2 v^T M* v (micro-burst signal)

    // Gershgorin stability margin of the mass we are MARCHING with, measured
    // against the CURRENT tangent: max_i (dt^2/4)*sum_j|K_ij| / M*_i, which is
    // (omega_max*dt/2)^2. <= 1 is stable; it equals massSafety^2 whenever the
    // tangent has not changed since M* was built, and exceeds 1 once the model
    // has stiffened past the mass — the silent-divergence detector.
    //
    // It is the WORST value since the last domainChanged, not the latest: a mass
    // rebuild resets the instantaneous margin to massSafety^2 by construction, so a
    // "latest" reading would let the refresh that papers over an excursion also
    // erase the evidence of it. This is the number robust_drive logs to audit a
    // whole rung-5 relaxation after the fact.
    //
    // Sampled TWICE over, deliberately: free at every M* rebuild (dense, and timed
    // at KE peaks where the model moves fastest) AND by an independent probe every
    // -marginEvery steps. Neither alone is enough — the rebuild sample is blind
    // under `-noAutoRefresh` (measured: 1.07 vs a flat 1.00 on the same diverging
    // run) and the fixed cadence is too coarse to catch a transient excursion.
    //
    // NEGATIVE means NOT MEASURED, and must never be read as "safe": -1 = no
    // gershgorin mass (lumped/unity carry no such bound), -2 = probe disabled
    // (-marginEvery 0).
    double getStabilityMargin(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    int buildFictitiousMass(void);   // populate *Mstar in domainChanged()
    // probe the LIVE tangent against the mass in use, without touching M*. Runs on
    // its own step cadence so the diagnostic never depends on the refresh policy.
    int measureStabilityMargin(void);
    // record + announce one margin sample (live-tangent mass vs marched mass)
    void applyMargin(const Vector &live, const Vector &march);

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
    Vector *Mprobe;                  // scratch row-sums for the margin probe
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

    // --- v3: explicit stability margin (note 83 §3) ---
    double massSafety;               // f in omega_max*dt <= 2f (gershgorin only)
    double stabMargin;               // last measured (omega_max*dt/2)^2; <0 = n/a
    bool   marginWarned;             // one-time boundary-crossing diagnostic
    int    marginEvery;              // probe cadence in steps (0 = off)
    int    marginCount;              // steps since the last probe
    bool   haveBuiltOnce;            // M* exists (own flag: stabMargin is a VALUE)
};

#endif
