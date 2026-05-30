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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024
// Revision: A
//
// Description: This file contains the class definition for ExplicitBathe.
// ExplicitBathe is an algorithmic class for performing a transient analysis
// using the explicit Bathe time integration scheme. 
//
// This algorithm is similar to a central difference scheme but perturbed to 
// introduce numerical damping. It is a second-order accurate explicit scheme. 
// Unlike the CentralDifference class, this one only assembles the mass matrix 
// on the right hand side, making it much easier to use with diagonal mass matrices. 
//
// The time-step required for stability is approximately twice that of
// standard explicit central difference schemes, making it more efficient for
// many applications. The optional critical-time-step computation below reports
// the CONSERVATIVE central-difference limit dt_cd = 2/omega_max (a guaranteed-
// safe lower bound for this scheme) together with the Noh-Bathe limit
// EB_NB_STABILITY_FACTOR * dt_cd (the published ~2x advantage at the optimal p).
//
// Recommended configuration (explicit recipe):
//   - lumped (diagonal) element mass, e.g. -lMass / "-mass" on elements
//   - system Diagonal      (trivial M^{-1}; consistent mass + sparse solver is
//                           wrong/pointless here)
//   - algorithm Linear     (the scheme needs exactly 2 solves/step)
//   - integrator ExplicitBathe 0.54
// Pair with `recorder EnergyBalance ... -file energy.txt` to watch for spurious
// energy growth (the instability signature of an over-large dt).
//
// Critical-time-step (-cfl / -tangent / -recompute / criticalTimeStep()) caveats:
//   - It is a per-element estimate using row-sum-lumped element mass. This is
//     exact/conservative for translational DOFs (bars, solids) but row-sum
//     lumping can be non-conservative for ROTATIONAL DOFs (beams/shells); treat
//     dt_cr as a guide there, not a guarantee.
//   - It ignores constraints (equalDOF/rigid diaphragms/MP) and pure nodal mass,
//     so a constrained or nodal-mass-only model may report a non-binding or
//     <=0 (criticalTimeStep() returns <=0) value.
//   - -recompute N / -tangent only UPDATE the reported dt_cr; they do NOT change
//     the analysis dt. For stiffening/contact, either pair with -cflAbort or
//     re-query criticalTimeStep() from the driver and sub-divide (see the
//     adaptive_substep example). It is O(N_elements) DGGEV per refresh -- use a
//     large N.
//
// Key features:
// - Second-order accurate
// - Conditionally stable (dt <= 2/omega_max central-difference reference)
// - Built-in numerical damping (controlled by parameter p, p=0.54 is a good choice)
// - Only mass matrix on RHS (no tangent matrix assembly required)
// - Optional automatic critical time step calculation
//
// Reference:
// Gunwoo Noh, Klaus-Jürgen Bathe,
// "An explicit time integration scheme for the analysis of wave propagations",
// Computers & Structures, Volume 129, 2013, Pages 178-193,
// ISSN 0045-7949,
// https://doi.org/10.1016/j.compstruc.2013.06.007.

#ifndef ExplicitBathe_h
#define ExplicitBathe_h

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>   // CTSLumping, CTSResult, computeCriticalTimeStep()

// Published Noh-Bathe stability advantage over central difference (~2x at the
// optimal sub-step parameter). dt_cr_NB = EB_NB_STABILITY_FACTOR * (2/omega_max).
// Reference: Noh & Bathe (2013), Computers & Structures 129, 178-193.
#define EB_NB_STABILITY_FACTOR 2.0

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitBathe : public TransientIntegrator
{
public:
    // Constructors
    ExplicitBathe();

    ExplicitBathe(double p, int compute_critical_timestep_ = 0,
                  bool verbose = false, bool cflAbort = false,
                  double divergenceFactor = 0.0,
                  bool cflUseTangent = false, int cflRecomputeEvery = 0,
                  CTSLumping lumping = CTSLumping::RowSum);

    // Destructor
    ~ExplicitBathe();

    // Conservative (central-difference) critical time step, available after a
    // step has run with compute_critical_timestep enabled; <=0 if not computed.
    double getCriticalTimeStep(void) const override;
    
    // Methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Method to obtain current velocity
    const Vector &getVel(void);
    
    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    // Time step
    double deltaT;

    // Explicit Bathe parameters
    double p;           // Sub-step parameter (0 < p < 1), controls damping
    double q0;          // Integration coefficient q0
    double q1;          // Integration coefficient q1
    double q2;          // Integration coefficient q2

    // State vectors at different time levels
    // The scheme has two sub-steps per time step:
    // Step 1: t -> t + p*dt (using acceleration A_tpdt)
    // Step 2: t + p*dt -> t + dt (using accelerations A_tpdt and A_tdt)
    
    Vector *U_t;        // Displacement at time t
    Vector *V_t;        // Velocity at time t 
    Vector *A_t;        // Acceleration at time t 
    
    Vector *U_tpdt;     // Displacement at time t + p*dt
    Vector *V_tpdt;     // Velocity at time t + p*dt
    Vector *V_fake;     // Temporary velocity for setting response
    Vector *A_tpdt;     // Acceleration at time t + p*dt
    
    Vector *U_tdt;      // Displacement at time t + dt
    Vector *V_tdt;      // Velocity at time t + dt
    Vector *A_tdt;      // Acceleration at time t + dt
    Vector *R_tdt;      // Forces at time t + dt (unused, kept for compatibility)

    // Update counter
    int updateCount;    // Counts updates per step (should be exactly 2)

    // Integration coefficients (computed from p)
    double a0;          // = p * dt
    double a1;          // = (p * dt)^2 / 2
    double a2;          // = a0 / 2
    double a3;          // = (1 - p) * dt
    double a4;          // = ((1 - p) * dt)^2 / 2
    double a5;          // = q0 * a3
    double a6;          // = (0.5 + q1) * a3
    double a7;          // = q2 * a3

    // Critical time step computation
    int compute_critical_timestep;          // Control flag for critical dt computation
    double damped_minimum_critical_timestep;    // Minimum critical dt (damped)
    double undamped_minimum_critical_timestep;  // Minimum critical dt (undamped)
    int damped_critical_element_tag;            // Element with minimum damped dt
    int undamped_critical_element_tag;          // Element with minimum undamped dt

    // Diagnostics / run control
    bool verbose;             // gate the per-step opserr reporting (default off)
    bool cflAbort;            // abort the step if dt exceeds the Noh-Bathe limit
    double divergenceFactor;  // >0: abort if kinetic energy grows by this factor
                              //     in one step (spurious-energy circuit breaker)
    double prevKE;            // previous-step kinetic-energy proxy (0.5*v.v)
    bool firstStep;           // first newStep() seeds prevKE / checks cold start

    // critical-time-step refresh policy
    bool cflUseTangent;       // use getTangentStiff() instead of getInitialStiff()
    int cflRecomputeEvery;    // recompute dt_cr every N steps (0 = once); tracks
                              // stiffness changes in nonlinear runs
    int cflStepCount;         // step counter for periodic recompute
    bool cflFirstComputation; // gate the detailed report to the first computation
    CTSLumping lumping;       // element-mass lumping for the dt_cr pencil (-lump)
};

#endif
