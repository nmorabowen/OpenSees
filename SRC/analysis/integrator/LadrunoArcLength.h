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

// Ladruno: LadrunoArcLength — fork-authored static integrator (classTag 33004).
//
// A strict superset of stock ArcLength that adds, default-OFF, the Ramm
// desired-iteration radius adaptation ("Layer A") that stock LoadControl /
// DisplacementControl / MinUnbalDispNorm carry but ArcLength never received.
// With no -adapt flag the class is bit-identical to ArcLength.
//
// Constraint (same as ArcLength):
//   i=1   dU^T dU + alpha^2 dLambda^2 = arcLength^2
//   i>1   the ArcLength quadratic in dLambda (positive-theta root).
//
// Layer A (when -adapt Jd ellMin ellMax [-p exp]): at the start of each step,
//   arcLength <- clamp( arcLength * (Jd / Jlast)^p , ellMin, ellMax )
// where Jlast = corrector iterations of the last converged step (Jlast<-max(1,.)
// to guard a predictor-only convergence). See
// Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.
//
// v1 deliberately omits the DDM sensitivity machinery of stock ArcLength
// (sensitivity is out of scope here); the viscous -stabilize mode and the
// Layer-B reduceStep/revert mutators are deferred to follow-up PRs.

#ifndef LadrunoArcLength_h
#define LadrunoArcLength_h

#include <StaticIntegrator.h>
#include <Vector.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

class LadrunoArcLength : public StaticIntegrator
{
  public:
    // new args are all defaulted so the broker `new LadrunoArcLength(1.0)` survives.
    LadrunoArcLength(double arcLength, double alpha = 1.0,
                     bool adapt = false, int Jd = 5,
                     double ellMin = 0.0, double ellMax = 0.0,
                     double pExp = 1.0,
                     bool stabilize = false, double fTarget = 2.0e-4,
                     bool adaptStab = false, int massMode = 0,
                     double massScale = 1.0, double cVisc = 0.0,
                     double dtPseudo = 1.0);
    ~LadrunoArcLength();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    // -stabilize seam: regularize K and inject the artificial viscous force.
    int formTangent(int statusFlag = CURRENT_TANGENT);   // K (+ cOverDt*M* if stabilize)
    int formUnbalance(void);                             // B (- f_v if stabilize)
    int commit(void);                                    // watchdog + adaptStab + calibrate

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    int buildArtificialMass(void);     // integrator-owned diagonal M* (lifted from DR)
    // --- stock ArcLength state (verbatim) ---
    double arcLength2;                 // arc radius, SQUARED
    double alpha2;                     // load-factor scaling, squared
    double a, b, c, b24ac;             // quadratic coefficients
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;                      // reference load vector
    double deltaLambdaStep, currentLambda;
    int    signLastDeltaLambdaStep;

    // --- Layer A: Ramm desired-iteration radius adaptation ---
    bool   adapt;                      // master gate (false => bit-identical)
    int    Jd;                         // desired corrector iterations / step
    double ellMin, ellMax;             // radius clamp
    double pExp;                       // exponent (1 = LoadControl, 0.5 = Crisfield)
    int    numIncrLastStep;            // corrector-iteration counter (Jlast)

    // --- -stabilize: Abaqus-STABILIZE-style viscous regularization ---
    // MUTUALLY EXCLUSIVE with the arc-length quadratic: -stabilize => the class
    // runs viscous-regularized incremental LOAD CONTROL (no constraint quadratic).
    bool   stabilize;                  // master gate (false => identity preserved)
    bool   adaptStab;                  // rescale cVisc each commit to hold f
    double fTarget;                    // dissipated-energy fraction (default 2e-4)
    double cVisc;                      // viscous factor (supplied via -cVisc, or calibrated)
    double dtPseudo;                   // fictitious dt (folded into cOverDt)
    double cOverDt;                    // = cVisc/dtPseudo, the one knob applied to the math
    bool   cCalibrated;                // one-shot calibration latch
    double Estrain0;                   // calibration baseline strain energy
    double dissipVisc;                 // cumulative viscous work (watchdog)
    double residualTrueNorm;           // ||lambda p - f_int|| WITHOUT f_v (watchdog)
    double dLambda0;                   // load increment for the stabilized predictor
    int    massMode;                   // 0 gershgorin (no dt^2/4) | 1 lumped | 2 unity
    double massScale;                  // for -mass lumped
    int    nEqn;                       // cached numEqn for the diagonal poke
    Vector *Mstar;                     // integrator-owned artificial diagonal mass
};

#endif
