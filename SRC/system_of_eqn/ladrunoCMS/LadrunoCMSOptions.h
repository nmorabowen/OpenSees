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

#ifndef LadrunoCMSOptions_h
#define LadrunoCMSOptions_h

#include <cstddef>
#include <string>
#include <vector>

namespace ladruno_cms {

enum class HierarchyMode { Auto, Logical };
enum class AssemblyVerification { Off, Signature, Full };
enum class RefinementMode { Subspace, None };
enum class ContributionKind { Stiffness, Mass };
enum class EquationRole { Unseen, Interior, Level2Interface, Level1Interface };

struct AssemblyMetrics {
    double maxAbs = 0.0;
    double sum = 0.0;
    double sumAbs = 0.0;
    double sumSquares = 0.0;
    double weightedProjection = 0.0;
    std::size_t numberOfValues = 0;
};

struct AssemblyRecord {
    std::size_t ordinal = 0;
    ContributionKind kind = ContributionKind::Stiffness;
    int rows = 0;
    int columns = 0;
    std::vector<int> originalIDs;
    std::vector<int> activeIDs;
    std::vector<double> upperValues;
    AssemblyMetrics metrics;
};

struct SymmetricCSR {
    int dimension = 0;
    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<double> upperValues;

    int validate(std::string &message) const;
    std::vector<double> toDense() const;
    std::vector<double> multiply(const std::vector<double> &x) const;
};

struct EffectiveOwnership {
    std::vector<std::vector<int>> ownersByEquation;
    std::vector<int> interiorOwner;
    std::vector<EquationRole> role;
};

struct CoordinateMassClassification {
    std::vector<int> dynamic;
    std::vector<int> massless;
    double threshold = 0.0;
};

struct Options {
    HierarchyMode hierarchy = HierarchyMode::Logical;
    int level1 = 0;
    int level2 = 0;
    int modesLevel2 = 0;
    int modesLevel1 = 0;
    double tolerance = 1.0e-8;
    int maxEnrich = 4;
    int maxIterations = 500;
    RefinementMode refinement = RefinementMode::Subspace;
    // Zero selects Bathe's q=max(p+8,2p) rule.
    int iterationVectors = 0;
    // Building 1A needs roughly O(10^2) inverse actions with the automatic
    // q=max(p+8,2p) block; the acceptance gate remains the original residual.
    int maxRefineIterations = 160;
    int denseMax = 2000;
    double assemblyRtol = 1.0e-12;
    double assemblyAtol = 1.0e-14;
    double massRtol = 1.0e-12;
    double massAtol = 1.0e-14;
    AssemblyVerification verifyAssembly = AssemblyVerification::Signature;
    std::size_t verifyFullMaxBytes = 268435456u;
    bool verbose = false;

    int resolvedIterationVectors(int numModes, std::string &message) const;
    int validate(int worldSize, int numModes, std::string &message) const;
};

struct DenseMemoryEstimate {
    std::size_t matricesAndVectors = 0;
    std::size_t staticTile = 0;
    std::size_t packedGather = 0;
    std::size_t knownLowerBound = 0;
};

int parseCommandOptions(
    const std::vector<std::string> &arguments,
    Options &options,
    int &numberOfModes,
    std::string &message);

DenseMemoryEstimate estimateDenseRootMemory(
    const int *localDimensions,
    int numberOfGroups,
    int dynamicDimension,
    int masslessDimension,
    int numModes,
    int staticBlockSize = 32);

AssemblyMetrics computeAssemblyMetrics(
    const double *rowMajorValues,
    int rows,
    int columns,
    double factor);

int makeAssemblyRecord(
    ContributionKind kind,
    std::size_t ordinal,
    const double *rowMajorValues,
    int rows,
    int columns,
    const std::vector<int> &originalIDs,
    int globalDimension,
    double factor,
    double relativeTolerance,
    double absoluteTolerance,
    AssemblyRecord &result,
    std::string &message);

bool compareAssemblyRecords(
    const AssemblyRecord &reference,
    const AssemblyRecord &candidate,
    double relativeTolerance,
    double absoluteTolerance,
    std::string &message);

int assignContributionOwner(
    const std::vector<int> &activeIDs,
    const std::vector<int> &fineLabels,
    std::string &message);

int classifyEffectiveOwnership(
    int numberOfEquations,
    const std::vector<AssemblyRecord> &contributions,
    const std::vector<int> &contributionOwners,
    const std::vector<int> &coarseOfFine,
    EffectiveOwnership &result,
    std::string &message);

int buildSymmetricCSR(
    int dimension,
    const std::vector<AssemblyRecord> &contributions,
    ContributionKind kind,
    SymmetricCSR &result,
    std::string &message);

int buildOwnedLocalPencil(
    const std::vector<AssemblyRecord> &stiffnessContributions,
    const std::vector<AssemblyRecord> &massContributions,
    std::vector<int> &globalEquations,
    SymmetricCSR &stiffness,
    SymmetricCSR &mass,
    std::string &message);

int classifyCoordinateMass(
    const SymmetricCSR &mass,
    double relativeTolerance,
    double absoluteTolerance,
    CoordinateMassClassification &result,
    std::string &message);

} // namespace ladruno_cms

#endif
