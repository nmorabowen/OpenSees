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

#ifndef ExplicitBatheLNVD_h
#define ExplicitBatheLNVD_h

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 XII 2024
// Revision: B  (local damping wired + activated, validation, resync; nmb 05/2026)
//
// Description: ExplicitBatheLNVD performs a transient analysis using the Noh-Bathe
// explicit scheme (see ExplicitBathe) with FLAC-style Local Non-Viscous Damping
// for dynamic relaxation / quasi-static solving. At every solve the assembled
// residual r = f_ext - f_int is modified to
//
//     r_i  <-  r_i - alpha * |r_i| * sign(v_i)
//
// i.e. an artificial force proportional to the local unbalanced-force magnitude
// and opposing the local velocity. It removes kinetic energy (damps oscillation
// toward equilibrium) and vanishes as v -> 0, so it does NOT bias the converged
// static solution. The damping is applied symmetrically to BOTH Noh-Bathe
// sub-steps (via the formUnbalance() override, which the solution algorithm
// calls for sub-step 1 and update() calls for sub-step 2).
//
//   alpha (FLAC local damping coefficient): 0 <= alpha < 1, classic default 0.8.
//   Use LNVD to drive a model to static equilibrium; watch the reported
//   unbalanced-force norm -> 0 to judge convergence. For true dynamics use the
//   undamped ExplicitBathe instead.
//
// Reference for base algorithm:
//   Gunwoo Noh, Klaus-Jürgen Bathe, "An explicit time integration scheme for the
//   analysis of wave propagations", Computers & Structures 129 (2013) 178-193,
//   https://doi.org/10.1016/j.compstruc.2013.06.007.

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>   // CTSLumping, CTSResult, computeCriticalTimeStep()

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitBatheLNVD : public TransientIntegrator
{
public:
    // constructors
    ExplicitBatheLNVD();
    ExplicitBatheLNVD(double p, double alpha_flac, int compute_critical_timestep_ = 0,
                      bool verbose = false, bool cflAbort = false,
                      double divergenceFactor = 0.0,
                      bool cflUseTangent = false, int cflRecomputeEvery = 0,
                      CTSLumping lumping = CTSLumping::RowSum);

    // destructor
    ~ExplicitBatheLNVD();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    // assemble the residual and add the FLAC local non-viscous damping
    int formUnbalance(void);

    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    const Vector &getVel(void);

    // conservative (central-difference) critical time step; <=0 if not computed
    double getCriticalTimeStep(void) const override;
    // infinity-norm of the most recent unbalanced force (dynamic-relaxation
    // convergence indicator); <0 until the first solve
    double getUnbalanceNorm(void) const;

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    // add the FLAC local non-viscous damping to the assembled SOE residual
    void addLocalDamping(void);

    double deltaT;

    // Explicit Bathe parameters
    double p;
    double q0;
    double q1;
    double q2;

    // FLAC local non-viscous damping coefficient (0 <= alpha < 1)
    double alpha_flac;

    // State variables
    Vector *U_t;                    // Displacement at time t
    Vector *V_t;                    // Velocity at time t
    Vector *A_t;                    // Acceleration at time t
    Vector *U_tpdt;                 // Displacement at time t + p*Dt
    Vector *V_tpdt;                 // Velocity at time t + p*Dt
    Vector *V_fake;                 // Trial velocity for setResponse
    Vector *A_tpdt;                 // Acceleration at time t + p*Dt
    Vector *U_tdt;                  // Displacement at time t + Dt
    Vector *V_tdt;                  // Velocity at time t + Dt
    Vector *A_tdt;                  // Acceleration at time t + Dt

    int updateCount;

    double a0, a1, a2, a3, a4, a5, a6, a7;

    int compute_critical_timestep;
    double damped_minimum_critical_timestep;
    double undamped_minimum_critical_timestep;
    int damped_critical_element_tag;
    int undamped_critical_element_tag;

    // Diagnostics / run control
    bool verbose;
    bool cflAbort;
    double divergenceFactor;
    double prevKE;
    bool firstStep;
    double lastUnbalanceNorm;       // |r|_inf at the most recent solve

    // critical-time-step refresh policy (see ExplicitBathe)
    bool cflUseTangent;
    int cflRecomputeEvery;
    int cflStepCount;
    bool cflFirstComputation;
    CTSLumping lumping;       // element-mass lumping for the dt_cr pencil (-lump)
};

#endif
