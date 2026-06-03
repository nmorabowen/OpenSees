/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
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

// Ladruno: LadrunoIndirectControl — indirect / CMOD displacement-control static
// integrator (classTag 33006). Generalizes stock DisplacementControl (a single
// nodal DOF) to a weighted MULTI-DOF control quantity
//     zeta = c . U = sum_k coef_k * U(node_k, dof_k)
// e.g. a crack-mouth-opening displacement (CMOD = u_A - u_B). Such a relative
// quantity stays MONOTONE through snap-back -- where the load factor AND the
// individual nodal DOFs both reverse -- so it follows equilibrium branches that
// geometric arc-length and single-DOF DisplacementControl cannot. This is the
// de Borst "indirect displacement control" answer to snap-back; ADR-20 section 8
// follow-up #3. See Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.
//
// Constraint (per step):  c . dU = Delta-zeta   (the control increment `incr`)
//   predictor (newStep): solve K dUhat = phat ;  dLambda = Delta-zeta/(c.dUhat) ;
//                        dU = dLambda dUhat
//   corrector (update):  solve K dUhat = phat ;  dLambda = -(c.dUbar)/(c.dUhat) ;
//                        dU = dUbar + dLambda dUhat
// where phat is the reference load vector and c the participation vector (coef_k
// at the controlled global equation numbers). Mirrors DisplacementControl's
// dUhat/dUbar bookkeeping verbatim, with the scalar control component
// dUahat = dUhat(theDofID) replaced by the dot product c.dUhat. Optional
// Ramm-style increment adaptation (-iter numIter dmin dmax), as DisplacementControl
// carries. No DDM sensitivity machinery (out of scope, as in LadrunoArcLength).

#ifndef LadrunoIndirectControl_h
#define LadrunoIndirectControl_h

#include <StaticIntegrator.h>
#include <Vector.h>
#include <ID.h>

class LinearSOE;
class AnalysisModel;
class Domain;
class FE_Element;

class LadrunoIndirectControl : public StaticIntegrator
{
  public:
    // nodes/dofs/coefs define the control quantity c.U (dofs are 1-based on input,
    // stored 0-based). incr = control increment Delta-zeta. numIncr>0 enables the
    // Ramm increment adaptation; min/max clamp the (signed-magnitude) increment.
    LadrunoIndirectControl(const ID &nodes, const ID &dofs, const Vector &coefs,
                           double incr, int numIncr = 1,
                           double min = 0.0, double max = 0.0);
    LadrunoIndirectControl();   // for the FEM_ObjectBroker (must recvSelf)
    ~LadrunoIndirectControl();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    double controlDot(const Vector &v) const;   // c.v over the controlled eqns

    // --- control definition (input) ---
    ID     ctrlNode;            // node tags
    ID     ctrlDof;             // 0-based dof indices
    Vector ctrlCoef;           // participation coefficients
    double theIncrement;        // control increment Delta-zeta (adapted if numIncr>0)
    int    specNumIncrStep;     // desired corrector iterations/step (Ramm target)
    int    numIncrLastStep;     // corrector iterations of the last step (Jlast)
    double minIncrement;        // increment-magnitude clamp
    double maxIncrement;

    // --- resolved at domainChanged ---
    ID     ctrlEqn;             // global equation number per control entry (-1 if free)

    // --- arc-length-style per-step state (mirrors DisplacementControl) ---
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;               // reference load vector
    double  deltaLambdaStep, currentLambda;
};

#endif
