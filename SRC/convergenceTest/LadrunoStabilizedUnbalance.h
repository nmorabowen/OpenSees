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

// Ladruno: LadrunoStabilizedUnbalance — a true-equilibrium NormUnbalance
// convergence test for the viscous-stabilized LadrunoArcLength -stabilize mode
// (ADR-20 §8 follow-up #4). classTag CONVERGENCE_TEST_LadrunoStabilizedUnbalance.
//
// In -stabilize mode the integrator regularizes the residual with an artificial
// viscous force f_v, so the SOE B that a stock CTestNormUnbalance norms is the
// POLLUTED `lambda p - f_int - f_v`: Newton is "converged" while the TRUE static
// unbalance `||lambda p - f_int||` is still nonzero by `||f_v||`. This test norms
// the TRUE unbalance instead -- a STRICTER criterion that drives the artificial
// force below tolerance, so the converged state is a genuine equilibrium point
// (not merely a regularized one).
//
// Mechanism (zero analysis-core surgery): setEquiSolnAlgo() captures BOTH the
// LinearSOE and the active IncrementalIntegrator. If that integrator is a
// LadrunoArcLength currently stabilizing, test() uses its getStabilizedTrueResidual()
// (the `||lambda p - f_int||` recorded in formUnbalance before f_v is subtracted);
// otherwise it degrades GRACEFULLY to the stock SOE `||B||` norm, so the test is a
// drop-in replacement for NormUnbalance with any integrator. Structure mirrors
// CTestNormUnbalance verbatim. See
// Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.

#ifndef LadrunoStabilizedUnbalance_h
#define LadrunoStabilizedUnbalance_h

#include <ConvergenceTest.h>
#include <Vector.h>

class EquiSolnAlgo;
class LinearSOE;
class IncrementalIntegrator;

class LadrunoStabilizedUnbalance : public ConvergenceTest
{
  public:
    LadrunoStabilizedUnbalance();
    LadrunoStabilizedUnbalance(double tol, int maxNumIter, int printFlag,
                               int normType = 2, double maxTol = OPS_MAXTOL);
    ~LadrunoStabilizedUnbalance();

    ConvergenceTest *getCopy(int iterations);

    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);

    int test(void);
    int start(void);

    int getNumTests(void);
    int getMaxNumTests(void);
    double getRatioNumToMax(void);
    const Vector &getNorms(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:
    // returns the residual norm to test against: the LadrunoArcLength true
    // unbalance when stabilizing, else the SOE ||B||. Sets usedTrue accordingly.
    double currentResidualNorm(bool &usedTrue) const;

  private:
    LinearSOE             *theSOE;
    IncrementalIntegrator *theIntegrator;   // for the stabilized true residual
    double tol;
    double maxTol;
    int    maxNumIter;
    int    currentIter;
    int    printFlag;
    int    nType;
    Vector norms;
};

#endif
