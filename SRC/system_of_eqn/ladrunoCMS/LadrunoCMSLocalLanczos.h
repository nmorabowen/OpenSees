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

#ifndef LadrunoCMSLocalLanczos_h
#define LadrunoCMSLocalLanczos_h

#include "LadrunoCMSOptions.h"

#include <string>
#include <vector>

namespace ladruno_cms {

struct LocalLanczosControls {
    double tolerance = 1.0e-10;
    int maximumOperatorApplications = 500;
    int maximumRestarts = 20;
    // Zero selects min(n,max(8*k+2*b,80)); positive values are a test/debug cap.
    int maximumBasisDimension = 0;
    int level = 0;
    int subdomain = 0;
};

struct LocalEigenResult {
    std::vector<double> eigenvalues;
    // Column-major matrix, dimension x numberOfModes.
    std::vector<double> eigenvectors;
    std::vector<double> normalizedResiduals;
    int operatorApplications = 0;
    int restarts = 0;
    int breakdowns = 0;
    int maximumBasisUsed = 0;
    int maximumRetainedAtRestart = 0;
    double maximumMOrthogonalityLoss = 0.0;
    // ADR-1000 P4 section 4. The section-27/29 profiles put ~65% of hierarchy
    // time inside this routine, so it gets its own attribution rather than
    // another round of optimizing by inspection. Seconds.
    double rayleighRitzSeconds = 0.0;
    double orthonormalizeSeconds = 0.0;
    double operatorSeconds = 0.0;
    double residualSeconds = 0.0;
    int rayleighRitzCalls = 0;
    int blockSolves = 0;
};

// Ritz values must be ascending. Returns -1 for invalid inputs.
int clusterSafeRetentionCount(
    const std::vector<double> &ascendingRitzValues,
    int nominalRetained,
    double relativeGapTolerance);

int solveLocalLanczos(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    int numberOfModes,
    const LocalLanczosControls &controls,
    LocalEigenResult &result,
    std::string &message);

} // namespace ladruno_cms

#endif
