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

#ifndef LadrunoCMSSubspaceRefiner_h
#define LadrunoCMSSubspaceRefiner_h

#include "LadrunoCMSOptions.h"

#include <string>
#include <vector>

namespace ladruno_cms {

struct SubspaceRefinementControls {
    int globalDimension = 0;
    int numberOfModes = 0;
    int numberOfIterationVectors = 0;
    int maximumIterations = 160;
    double tolerance = 1.0e-8;
    double massRankTolerance = 0.0;
};

struct SubspaceIterationDiagnostic {
    int iteration = 0;
    double maximumResidual = 0.0;
    double massOrthogonalityLoss = 0.0;
};

struct SubspaceRefinementResult {
    std::vector<double> eigenvalues;
    // Column-major globalDimension x numberOfModes, complete on every rank.
    std::vector<double> eigenvectors;
    std::vector<double> normalizedResiduals;
    std::vector<double> residualNorms;
    std::vector<double> stiffnessActionNorms;
    std::vector<double> massActionNorms;
    std::vector<int> maximumResidualEquation;
    std::vector<SubspaceIterationDiagnostic> history;
    int iterations = 0;
    int numberOfIterationVectors = 0;
    bool converged = false;
    double maximumMassOrthogonalityLoss = 0.0;
    double factorizationSeconds = 0.0;
    double solveSeconds = 0.0;
    double inverseRefinementSeconds = 0.0;
    double massActionSeconds = 0.0;
    double rayleighRitzSeconds = 0.0;
};

// Collective on MPI_COMM_WORLD.  Each CSR contains uniquely owned entries of
// the original pencil, while initialVectors is complete on every rank.
int refineDistributedSubspace(
    const SymmetricCSR &localStiffness,
    const SymmetricCSR &localMass,
    const std::vector<double> &initialVectors,
    const SubspaceRefinementControls &controls,
    SubspaceRefinementResult &result,
    std::string &message);

} // namespace ladruno_cms

#endif
