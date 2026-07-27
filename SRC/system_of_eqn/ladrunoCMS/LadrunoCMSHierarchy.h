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

#ifndef LadrunoCMSHierarchy_h
#define LadrunoCMSHierarchy_h

#include "LadrunoCMSOptions.h"

#include <string>
#include <vector>

namespace ladruno_cms {

struct FineSubdomainInput {
    int fine = -1;
    int coarse = -1;
    // Local row -> original global equation.
    std::vector<int> equations;
    SymmetricCSR stiffness;
    SymmetricCSR mass;
    // Negative means a complete fixed-interface basis.
    int modesToKeep = -1;
};

struct TwoLevelHierarchyInput {
    SymmetricCSR globalStiffness;
    SymmetricCSR globalMass;
    std::vector<FineSubdomainInput> fineSubdomains;
    int numberOfCoarseSubdomains = 0;
    // Negative means a complete level-1 fixed-interface basis.
    int modesLevel1 = -1;
    int numberOfModes = 0;
    int denseMax = 2000;
    double localTolerance = 1.0e-10;
    int maximumOperatorApplications = 500;
    int maximumRestarts = 20;
    double massRtol = 1.0e-12;
    double massAtol = 1.0e-14;
    double assemblyRtol = 1.0e-12;
    double assemblyAtol = 1.0e-14;
    // Test/diagnostic aid only. Production reconstruction acts directly on the
    // requested modes and does not materialize the global n x r transformation.
    bool storeTotalTransformation = false;
    std::size_t maximumTransformationEntries = 1000000u;
    // BENCHMARK DIAGNOSTIC ONLY -- omit the level-1 Craig-Bampton reduction (T1)
    // and let S1 assemble the un-reduced group pencils, i.e. the "level-2 only"
    // ablation of ADR-1000 P4 section 2b. It exists to attribute cost/reduction
    // between the two levels and is NEVER an accepted solver configuration or a
    // fallback: no command option can set it, and LadrunoCMSEigenSolver refuses
    // an input that carries it. A run with this flag reports
    // diagnostics.appliedT1 == false and diagnostics.ablatedLevel1 == true so it
    // can never be mistaken for the mandatory T2 -> S2 -> T1 -> S1 chain.
    bool diagnosticAblateLevel1 = false;
};

struct HierarchyDiagnostics {
    bool appliedT2 = false;
    bool appliedS2 = false;
    bool appliedT1 = false;
    bool appliedS1 = false;
    // True only for the benchmark-only level-1 ablation above. Any consumer that
    // treats a result as the mandatory chain must reject this.
    bool ablatedLevel1 = false;
    int originalDimension = 0;
    int afterLevel2BeforeCompatibility = 0;
    int afterLevel2Compatibility = 0;
    int afterLevel1BeforeCompatibility = 0;
    int finalRawDimension = 0;
    int finalMasslessDimension = 0;
    int finalDynamicDimension = 0;
    int nestedEnrichmentDimension = 0;
    int originalMasslessDimension = 0;
    double maximumDuplicateJump = 0.0;
    std::vector<int> maximumResidualEquation;
    // 0: fine interior, 1: level-2 interface, 2: level-1 interface.
    std::vector<int> maximumResidualRole;
    std::vector<int> maximumResidualMassless;
    std::vector<double> residualNorm;
    std::vector<double> stiffnessActionNorm;
    std::vector<double> massActionNorm;
    std::vector<int> level2CompatibilityCounts;
    std::vector<int> level1CompatibilityCounts;
    // ADR-1000 P4 section 1 -- per-phase wall clock of the DISTRIBUTED
    // hierarchy, so `hierarchySeconds` can be attributed instead of guessed at.
    // Seconds, measured on this rank; zero on paths that never run the phase
    // (the serial two-level entry point leaves them all zero).
    double partitionSeconds = 0.0;        // owner counts + interior/boundary split
    double fineModesSeconds = 0.0;        // T2, the local Craig-Bampton reduction
    double compatibilitySeconds = 0.0;    // S2, the group merge
    double level1Seconds = 0.0;           // T1, the coarse reduction
    double globalSolveSeconds = 0.0;      // S1 + the leader/global solve
    double backSubstitutionSeconds = 0.0; // reconstruction to global coordinates
    double publicationSeconds = 0.0;      // gathering/publishing the eigenvectors
    // Peak resident set of THIS rank, bytes; 0 when the platform query fails.
    std::size_t peakResidentBytes = 0;
};

struct TwoLevelHierarchyResult {
    std::vector<double> eigenvalues;
    // Column-major originalDimension x numberOfModes.
    std::vector<double> eigenvectors;
    std::vector<double> normalizedResiduals;
    // Optional, column-major originalDimension x finalRawDimension.
    std::vector<double> totalTransformation;
    SymmetricCSR finalStiffness;
    SymmetricCSR finalMass;
    HierarchyDiagnostics diagnostics;
};

struct DistributedHierarchyInput {
    int globalDimension = 0;
    int fine = -1;
    int coarse = -1;
    int numberOfCoarseSubdomains = 0;
    std::vector<int> localEquations;
    SymmetricCSR localStiffness;
    SymmetricCSR localMass;
    const std::vector<AssemblyRecord> *ownedStiffnessContributions = nullptr;
    const std::vector<AssemblyRecord> *ownedMassContributions = nullptr;
    int modesLevel2 = 0;
    int modesLevel1 = 0;
    int numberOfModes = 0;
    int denseMax = 2000;
    double tolerance = 1.0e-8;
    int maximumOperatorApplications = 500;
    // Restart budget of the local fixed-interface Lanczos (T2). See
    // Options::maxRestarts -- was hard-coded to 20 here.
    int maximumRestarts = 20;
    double massRtol = 1.0e-12;
    double massAtol = 1.0e-14;
};

struct DistributedHierarchyResult {
    std::vector<double> eigenvalues;
    // Column-major globalDimension x numberOfModes, complete on every rank.
    std::vector<double> eigenvectors;
    std::vector<double> normalizedResiduals;
    HierarchyDiagnostics diagnostics;
};

int solveTwoLevelHierarchy(
    const TwoLevelHierarchyInput &input,
    TwoLevelHierarchyResult &result,
    std::string &message);

int solveDistributedHierarchy(
    const DistributedHierarchyInput &input,
    DistributedHierarchyResult &result,
    std::string &message);

// Forms an explicitly nested Rayleigh--Ritz space from an accepted result and
// a newly computed CMS result.  The input vectors are complete on every rank;
// K and M actions remain distributed through uniquely owned contributions.
int solveNestedRitzUnion(
    int globalDimension,
    int numberOfModes,
    const std::vector<AssemblyRecord> &ownedStiffnessContributions,
    const std::vector<AssemblyRecord> &ownedMassContributions,
    const DistributedHierarchyResult &previous,
    DistributedHierarchyResult &candidate,
    std::string &message);

} // namespace ladruno_cms

#endif
