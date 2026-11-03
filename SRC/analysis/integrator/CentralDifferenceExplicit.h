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

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell
//
// References:
//   Belytschko, T., Liu, W.K., Moran, B., "Nonlinear Finite Elements for
//     Continua and Structures", Wiley, 2000 (central difference, stable step).
//   Hughes, T.J.R., "The Finite Element Method", Prentice-Hall, 1987.
//   Goudreau, G.L., Hallquist, J.O., "Recent developments in large scale finite
//     element Lagrangian hydrocode technology", Computer Methods in Applied
//     Mechanics and Engineering, 30, 1982.
//
// Description: explicit leap-frog central difference. Mass-only left hand side,
// viscous damping lagged at the known half-step velocity:
//
//     a_n       = M^-1 ( P_n - C v_{n-1/2} - F_int(u_n) )
//     v_{n+1/2} = v_{n-1/2} + dt a_n
//     u_{n+1}   = u_n + dt v_{n+1/2}
//
// The scheme is Newmark(beta=0, gamma=1/2): second order, no algorithmic
// dissipation, period shortening of order (omega dt)^2/24.
//
// Velocities. getVel() returns the half-step velocity v_{n-1/2}, which is the
// velocity seen by the residual and by modal damping. The nodes receive the
// full-step velocity v_n = v_{n-1/2} + dt/2 a_n after each solve, so recorders
// read a velocity synchronized with the displacement.
//
// First step. a_0 = M^-1 ( P_0 - C v_0 - F_int(u_0) ) and v_{-1/2} = v_0 - dt/2 a_0
// are evaluated on the first newStep(), where dt is known and the system of
// equations can be formed. The result is u_1 = u_0 + dt v_0 + dt^2/2 a_0.
//
// Stability. The limit is dt_cr = (2/omega_max) (sqrt(1+xi^2) - xi), with
// xi = 0 giving 2/omega_max. getCriticalTimeStep() reports the estimate of
// CriticalTimeStep.h. A damped limit well below the undamped limit is the
// signature of stiffness-proportional Rayleigh damping, for which dt_cr scales
// with 2/(betaK omega_max^2); the class prints a note in that case. Refusal is
// opt-in through -cflAbort.
//
// Required recipe: lumped or diagonal mass, `system Diagonal`, `algorithm Linear`
// (exactly one solve per step), dt below dt_cr. For a coupled, consistent-mass
// or implicitly damped central difference `integrator NewmarkExplicit 0.5`
// remains the appropriate class.

#ifndef CentralDifferenceExplicit_h
#define CentralDifferenceExplicit_h

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>

class DOF_Group;
class FE_Element;
class Vector;
class LinearSOE;

class CentralDifferenceExplicit : public TransientIntegrator
{
public:
    CentralDifferenceExplicit();
    CentralDifferenceExplicit(int compute_critical_timestep,
                              bool verbose, bool cflAbort,
                              double divergenceFactor,
                              bool cflUseTangent, int cflRecomputeEvery,
                              CTSLumping lumping);
    ~CentralDifferenceExplicit();

    // Damped limit when damping is active, 2/omega_max otherwise; <=0 when no
    // element produced an estimate.
    double getCriticalTimeStep(void) const override;

    // Mass only on the left hand side; no K, no C, no dt.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    // The assembled tangent is pure mass, so after the first assembly the
    // factored system is reused until domainChanged() or revertToLastStep()
    // invalidates it. An updateParameter on a density changes M without a domain
    // change: such a model needs -noMassCache.
    int formTangent(int statFlag);

    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Called after Domain::revertToLastCommit() on a failed step. Re-seeds the
    // private leap-frog state from the committed node state and re-arms the
    // first-step starter, so a retry centers v_{-1/2} for the dt it uses. The
    // restart is machine-exact without damping and second order with damping,
    // because the starter evaluates C v at the full-step velocity.
    int revertToLastStep(void);

    const Vector &getVel(void);

    // -noMassCache restores the mass assembly at every step.
    void setMassCacheEnabled(bool on) { useMassCache = on; if (!on) massTangentValid = false; }

    // -commitSolveState skips the second element state determination in update().
    // Element state is then committed at the strain u_{n+1} and the half-step
    // velocity v_{n+1/2} used by the solve. Bit-identical for rate-independent
    // materials; materials with rate-dependent memory commit a slightly different
    // state, hence opt-in.
    void setCommitSolveState(bool on) { commitSolveState = on; }

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

private:
    double deltaT;

    Vector *Ut;        // displacement u_n, advanced to u_{n+1} during newStep
    Vector *Vhalf;     // half-step velocity v_{n-1/2}
    Vector *Aprev;     // acceleration a_n carried between steps
    Vector *Vfull;     // full-step output velocity
    Vector *Azero;     // zero vector: the residual carries no inertia term

    int updateCount;   // exactly one solve per step
    bool firstStep;    // first-step starter pending

    int    compute_critical_timestep;     // 0 = off, 1 = (re)compute, 2 = computed
    double damped_minimum_critical_timestep;
    double undamped_minimum_critical_timestep;
    int    damped_critical_element_tag;
    int    undamped_critical_element_tag;

    bool   verbose;
    bool   cflAbort;
    double divergenceFactor;   // >0: abort on kinetic energy proxy growth
    double prevKE;             // running maximum of the kinetic energy proxy. A
                               // previous-step baseline trips falsely at the
                               // kinetic energy troughs of free vibration.
    double committedPrevKE;    // prevKE at newStep() entry, for revertToLastStep

    bool   cflUseTangent;      // getTangentStiff() instead of getInitialStiff()
    int    cflRecomputeEvery;  // refresh dt_cr every N steps (0 = once)
    int    cflStepCount;
    int    committedCflStepCount;
    bool   cflFirstComputation;
    CTSLumping lumping;
    bool   betaKWarned;

    void seedFromCommitted(AnalysisModel *theModel);

    bool   useMassCache;
    bool   massTangentValid;
    bool   massCacheNoted;

    bool   commitSolveState;
    bool   commitSolveStateNoted;
};

#endif
