#include "LadrunoCMSHierarchy.h"
#include "LadrunoCMSMumps.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#ifdef _WIN32
extern "C" int DSYGVX(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_TEST_DSYGVX DSYGVX
#else
extern "C" int dsygvx_(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_TEST_DSYGVX dsygvx_
#endif

namespace {

int failures = 0;

void require(bool condition, const std::string &message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}
// Evaluates the call FIRST, then builds the diagnostic. As two arguments of
// require() their evaluation order is unspecified in C++, and MSVC was building
// the message string BEFORE the call filled `message` -- every failure printed
// an empty diagnostic. See LEDGER_quirks (2026-07-26).
#define REQUIRE_CALL(status, text)                                                 do {                                                                               const bool ladrunoCallOk = (status);                                           require(ladrunoCallOk, text);                                              } while (0)


bool close(double left, double right, double tolerance = 1.0e-8)
{
    return std::fabs(left - right) <= tolerance *
        std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
}

ladruno_cms::SymmetricCSR csrFromDense(const std::vector<double> &dense, int order)
{
    ladruno_cms::SymmetricCSR result;
    result.dimension = order;
    result.rowOffsets.push_back(0);
    for (int row = 0; row < order; ++row) {
        for (int column = row; column < order; ++column) {
            const double value = dense[static_cast<std::size_t>(row) * order + column];
            if (value != 0.0) {
                result.columnIndices.push_back(column);
                result.upperValues.push_back(value);
            }
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    return result;
}

struct Fixture {
    ladruno_cms::TwoLevelHierarchyInput input;
    std::vector<double> stiffness;
    std::vector<double> mass;
};

Fixture makeFixture(bool completeBasis)
{
    Fixture fixture;
    constexpr int order = 9;
    fixture.stiffness.assign(order * order, 0.0);
    fixture.mass.assign(order * order, 0.0);
    const int equations[4][3] = {{0, 1, 2}, {2, 3, 4}, {4, 5, 6}, {6, 7, 8}};
    for (int fine = 0; fine < 4; ++fine) {
        const double shift = 0.15 * fine;
        const std::vector<double> localStiffness = {
            2.2 + shift, -1.0, 0.0,
            -1.0, 2.6 + shift, -0.8,
            0.0, -0.8, 1.9 + shift};
        const std::vector<double> localMass = {
            0.8 + 0.05 * fine, 0.04, 0.0,
            0.04, 1.0 + 0.05 * fine, 0.03,
            0.0, 0.03, 1.2 + 0.05 * fine};
        ladruno_cms::FineSubdomainInput local;
        local.fine = fine;
        local.coarse = fine < 2 ? 0 : 1;
        local.equations.assign(equations[fine], equations[fine] + 3);
        local.stiffness = csrFromDense(localStiffness, 3);
        local.mass = csrFromDense(localMass, 3);
        local.modesToKeep = completeBasis ? -1 : 1;
        fixture.input.fineSubdomains.push_back(local);
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                const std::size_t global =
                    static_cast<std::size_t>(equations[fine][row]) * order +
                    equations[fine][column];
                fixture.stiffness[global] += localStiffness[static_cast<std::size_t>(row) * 3 + column];
                fixture.mass[global] += localMass[static_cast<std::size_t>(row) * 3 + column];
            }
        }
    }
    fixture.input.globalStiffness = csrFromDense(fixture.stiffness, order);
    fixture.input.globalMass = csrFromDense(fixture.mass, order);
    fixture.input.numberOfCoarseSubdomains = 2;
    fixture.input.modesLevel1 = completeBasis ? -1 : 2;
    fixture.input.numberOfModes = completeBasis ? 5 : 3;
    fixture.input.denseMax = 20;
    fixture.input.localTolerance = 1.0e-10;
    fixture.input.maximumOperatorApplications = 300;
    fixture.input.maximumRestarts = 20;
    fixture.input.storeTotalTransformation = completeBasis;
    return fixture;
}

std::vector<double> directEigenvalues(
    const ladruno_cms::SymmetricCSR &stiffnessCSR,
    const ladruno_cms::SymmetricCSR &massCSR,
    int requested)
{
    const int order = stiffnessCSR.dimension;
    const std::vector<double> stiffnessRows = stiffnessCSR.toDense();
    const std::vector<double> massRows = massCSR.toDense();
    std::vector<double> stiffness(static_cast<std::size_t>(order) * order);
    std::vector<double> mass(static_cast<std::size_t>(order) * order);
    for (int row = 0; row < order; ++row) {
        for (int column = 0; column < order; ++column) {
            stiffness[static_cast<std::size_t>(column) * order + row] =
                stiffnessRows[static_cast<std::size_t>(row) * order + column];
            mass[static_cast<std::size_t>(column) * order + row] =
                massRows[static_cast<std::size_t>(row) * order + column];
        }
    }
    int itype = 1;
    char jobz = 'V';
    char range = 'I';
    char uplo = 'U';
    int n = order;
    int leading = order;
    double lower = 0.0;
    double upper = 0.0;
    int first = 1;
    int last = requested;
    double absoluteTolerance = 0.0;
    int found = 0;
    std::vector<double> eigenvalues(static_cast<std::size_t>(order));
    std::vector<double> eigenvectors(static_cast<std::size_t>(order) * requested);
    int leadingVectors = order;
    int workspaceSize = 8 * order;
    std::vector<double> workspace(static_cast<std::size_t>(workspaceSize));
    std::vector<int> integerWorkspace(static_cast<std::size_t>(5 * order));
    std::vector<int> failed(static_cast<std::size_t>(order));
    int info = 0;
    LADRUNO_TEST_DSYGVX(
        &itype, &jobz, &range, &uplo, &n, stiffness.data(), &leading,
        mass.data(), &leading, &lower, &upper, &first, &last,
        &absoluteTolerance, &found, eigenvalues.data(), eigenvectors.data(),
        &leadingVectors, workspace.data(), &workspaceSize,
        integerWorkspace.data(), failed.data(), &info);
    require(info == 0 && found == requested, "direct LAPACK reference failed");
    eigenvalues.resize(static_cast<std::size_t>(std::max(0, found)));
    return eigenvalues;
}

void checkCompleteTwoLevelFlow()
{
    Fixture fixture = makeFixture(true);
    const std::vector<double> reference = directEigenvalues(
        fixture.input.globalStiffness, fixture.input.globalMass, 5);
    ladruno_cms::TwoLevelHierarchyResult result;
    std::string message;
    REQUIRE_CALL(ladruno_cms::solveTwoLevelHierarchy(fixture.input, result, message) == 0,
            "complete two-level hierarchy failed: " + message);
    require(result.diagnostics.appliedT2 && result.diagnostics.appliedS2 &&
            result.diagnostics.appliedT1 && result.diagnostics.appliedS1,
            "one of T2/S2/T1/S1 was not applied");
    require(result.diagnostics.afterLevel2BeforeCompatibility == 12,
            "unexpected pre-S2 dimension");
    require(result.diagnostics.afterLevel2Compatibility == 10,
            "S2 did not merge both within-group interfaces");
    require(result.diagnostics.afterLevel1BeforeCompatibility == 10,
            "complete T1 changed the assembled coarse dimension");
    require(result.diagnostics.finalRawDimension == 9,
            "S1 did not merge the cross-group interface");
    require(*std::max_element(result.diagnostics.level2CompatibilityCounts.begin(),
                              result.diagnostics.level2CompatibilityCounts.end()) == 2,
            "S2 compatibility multiplicity was not observed");
    require(*std::max_element(result.diagnostics.level1CompatibilityCounts.begin(),
                              result.diagnostics.level1CompatibilityCounts.end()) == 2,
            "S1 compatibility multiplicity was not observed");
    require(result.totalTransformation.size() == 81u,
            "bounded diagnostic transformation was not reconstructed");
    require(result.diagnostics.maximumDuplicateJump < 2.0e-10,
            "shared equations are discontinuous after back-substitution");
    for (std::size_t mode = 0; mode < reference.size(); ++mode) {
        require(close(result.eigenvalues[mode], reference[mode], 2.0e-8),
                "complete CMS eigenvalue differs from direct LAPACK");
        require(result.normalizedResiduals[mode] < 2.0e-8,
                "complete CMS mode has an excessive global residual");
    }
}

void checkTruncatedTwoLevelFlow()
{
    Fixture fixture = makeFixture(false);
    const std::vector<double> reference = directEigenvalues(
        fixture.input.globalStiffness, fixture.input.globalMass, 3);
    ladruno_cms::TwoLevelHierarchyResult result;
    std::string message;
    REQUIRE_CALL(ladruno_cms::solveTwoLevelHierarchy(fixture.input, result, message) == 0,
            "truncated two-level hierarchy failed: " + message);
    require(result.diagnostics.finalRawDimension == 5,
            "truncated hierarchy has an unexpected final dimension");
    require(result.totalTransformation.empty(),
            "production path materialized the global transformation");
    require(result.diagnostics.maximumDuplicateJump < 2.0e-10,
            "truncated back-substitution violates compatibility");
    for (std::size_t mode = 0; mode < reference.size(); ++mode) {
        require(result.eigenvalues[mode] + 1.0e-9 >= reference[mode],
                "truncated Ritz eigenvalue violates the upper-bound property");
        require(std::isfinite(result.normalizedResiduals[mode]),
                "truncated hierarchy returned a non-finite residual");
    }
}

void checkCoordinateMassCondensationAcrossHierarchy()
{
    Fixture fixture = makeFixture(true);
    const std::vector<double> localMass = {
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.03,
        0.0, 0.03, 1.2};
    fixture.input.fineSubdomains[0].mass = csrFromDense(localMass, 3);
    fixture.mass[0] = 0.0;
    fixture.mass[1] = 0.0;
    fixture.mass[9] = 0.0;
    fixture.input.globalMass = csrFromDense(fixture.mass, 9);

    ladruno_cms::StaticCondensation directCondensation;
    std::string message;
    REQUIRE_CALL(ladruno_cms::condenseCoordinateMass(
                fixture.input.globalStiffness, fixture.input.globalMass,
                fixture.input.massRtol, fixture.input.massAtol,
                directCondensation, message) == 0,
            "direct coordinate-mass condensation failed: " + message);
    const std::vector<double> reference = directEigenvalues(
        directCondensation.reducedStiffness,
        directCondensation.reducedMass, 5);

    ladruno_cms::TwoLevelHierarchyResult result;
    message.clear();
    REQUIRE_CALL(ladruno_cms::solveTwoLevelHierarchy(fixture.input, result, message) == 0,
            "hierarchical coordinate-mass condensation failed: " + message);
    require(result.diagnostics.finalRawDimension == 8,
            "hierarchy did not eliminate the finite-inertia coordinate at T2");
    for (std::size_t mode = 0; mode < reference.size(); ++mode) {
        require(close(result.eigenvalues[mode], reference[mode], 3.0e-8),
                "hierarchical finite eigenvalue differs after mass condensation");
        require(result.normalizedResiduals[mode] < 3.0e-8,
                "reconstructed finite mode has an excessive residual");
    }
}

void checkExplicitMemoryGuards()
{
    Fixture denseFixture = makeFixture(true);
    denseFixture.input.numberOfModes = 3;
    denseFixture.input.denseMax = 4;
    denseFixture.input.storeTotalTransformation = false;
    ladruno_cms::TwoLevelHierarchyResult result;
    std::string message;
    require(ladruno_cms::solveTwoLevelHierarchy(
                denseFixture.input, result, message) < 0,
            "denseMax guard accepted an oversized final pencil");
    require(message.find("denseMax") != std::string::npos,
            "denseMax rejection did not name its controlling option");

    Fixture transformationFixture = makeFixture(true);
    transformationFixture.input.maximumTransformationEntries = 80u;
    message.clear();
    require(ladruno_cms::solveTwoLevelHierarchy(
                transformationFixture.input, result, message) < 0,
            "diagnostic transformation entry guard was bypassed");
    require(message.find("entry limit") != std::string::npos,
            "diagnostic transformation guard was not actionable");
}

void checkDistributedFourRankFlow(int rank, int size)
{
    if (size != 4)
        return;
    Fixture fixture = makeFixture(true);
    const std::vector<double> reference = directEigenvalues(
        fixture.input.globalStiffness, fixture.input.globalMass, 5);
    const ladruno_cms::FineSubdomainInput &local =
        fixture.input.fineSubdomains[static_cast<std::size_t>(rank)];
    const std::vector<double> stiffnessDense = local.stiffness.toDense();
    const std::vector<double> massDense = local.mass.toDense();
    ladruno_cms::AssemblyRecord stiffnessRecord;
    ladruno_cms::AssemblyRecord massRecord;
    std::string message;
    REQUIRE_CALL(ladruno_cms::makeAssemblyRecord(
                ladruno_cms::ContributionKind::Stiffness, 0,
                stiffnessDense.data(), 3, 3, local.equations, 9, 1.0,
                1.0e-12, 1.0e-14, stiffnessRecord, message) == 0,
            "distributed stiffness record construction failed: " + message);
    REQUIRE_CALL(ladruno_cms::makeAssemblyRecord(
                ladruno_cms::ContributionKind::Mass, 0,
                massDense.data(), 3, 3, local.equations, 9, 1.0,
                1.0e-12, 1.0e-14, massRecord, message) == 0,
            "distributed mass record construction failed: " + message);

    ladruno_cms::DistributedHierarchyInput input;
    input.globalDimension = 9;
    input.fine = rank;
    input.coarse = rank / 2;
    input.numberOfCoarseSubdomains = 2;
    input.localEquations = local.equations;
    input.localStiffness = local.stiffness;
    input.localMass = local.mass;
    std::vector<ladruno_cms::AssemblyRecord> ownedStiffness{stiffnessRecord};
    std::vector<ladruno_cms::AssemblyRecord> ownedMass{massRecord};
    input.ownedStiffnessContributions = &ownedStiffness;
    input.ownedMassContributions = &ownedMass;
    input.modesLevel2 = 3;
    input.modesLevel1 = 10;
    input.numberOfModes = 5;
    input.denseMax = 20;
    input.tolerance = 1.0e-8;
    input.maximumOperatorApplications = 300;
    ladruno_cms::DistributedHierarchyResult result;
    message.clear();
    REQUIRE_CALL(ladruno_cms::solveDistributedHierarchy(input, result, message) == 0,
            "distributed T2/S2/T1/S1 flow failed: " + message);
    require(result.diagnostics.appliedT2 && result.diagnostics.appliedS2 &&
            result.diagnostics.appliedT1 && result.diagnostics.appliedS1,
            "distributed flow skipped a mandatory transformation");
    require(result.diagnostics.afterLevel2BeforeCompatibility == 12 &&
            result.diagnostics.afterLevel2Compatibility == 10 &&
            result.diagnostics.afterLevel1BeforeCompatibility == 10 &&
            result.diagnostics.finalRawDimension == 9,
            "distributed hierarchy dimensions differ from the serial oracle");
    require(result.eigenvectors.size() == 45u,
            "distributed reconstruction did not publish complete vectors");
    require(result.diagnostics.maximumDuplicateJump < 2.0e-10,
            "distributed reconstruction has incompatible interface copies");
    for (std::size_t mode = 0; mode < reference.size(); ++mode) {
        require(close(result.eigenvalues[mode], reference[mode], 3.0e-8),
                "distributed eigenvalue differs from direct LAPACK");
        require(result.normalizedResiduals[mode] < 3.0e-8,
                "distributed complete-basis residual is excessive");
    }

    // Regression: recomputed local CMS bases are not inherently nested.  The
    // production enrichment must therefore preserve the preceding Ritz space
    // explicitly before checking Rayleigh--Ritz monotonicity.
    ladruno_cms::DistributedHierarchyResult enriched = result;
    for (int mode = 0; mode < input.numberOfModes; ++mode) {
        const int perturbationRow = (mode + 5) % input.globalDimension;
        enriched.eigenvectors[
            static_cast<std::size_t>(mode) * input.globalDimension +
            perturbationRow] += 2.0e-2;
    }
    message.clear();
    REQUIRE_CALL(ladruno_cms::solveNestedRitzUnion(
                input.globalDimension, input.numberOfModes,
                ownedStiffness, ownedMass, result, enriched, message) == 0,
            "nested distributed enrichment failed: " + message);
    require(enriched.diagnostics.nestedEnrichmentDimension > input.numberOfModes,
            "nested enrichment discarded every independent new direction");
    for (int mode = 0; mode < input.numberOfModes; ++mode) {
        const double allowance = 1.0e-10 *
            std::max(1.0, std::fabs(result.eigenvalues[static_cast<std::size_t>(mode)]));
        require(enriched.eigenvalues[static_cast<std::size_t>(mode)] <=
                    result.eigenvalues[static_cast<std::size_t>(mode)] + allowance,
                "nested enrichment increased a Ritz eigenvalue");
        require(enriched.normalizedResiduals[static_cast<std::size_t>(mode)] < 5.0e-8,
                "nested enrichment damaged an exact distributed mode");
    }
    std::vector<double> rootVectors = result.eigenvectors;
    MPI_Bcast(rootVectors.data(), static_cast<int>(rootVectors.size()),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (std::size_t entry = 0; entry < rootVectors.size(); ++entry)
        require(close(result.eigenvectors[entry], rootVectors[entry], 1.0e-12),
                "published eigenvectors differ between ranks");
}


// ---------------------------------------------------------------------------
// P3d -- rank/partition-permutation invariance.
//
// Which rank happens to own which fine subdomain is an accident of the
// partitioner. Permuting that mapping must leave the eigenvalues, the residuals
// AND the subspace (MAC ~ 1) untouched; if it does not, the hierarchy is
// smuggling rank identity into the mathematics. Two permutations are exercised:
// a reversal (coarse groups stay contiguous rank ranges) and an interleave
// (coarse group 0 becomes ranks {0,2}), which is the one that would catch an
// implementation assuming a coarse group is a contiguous block of ranks.
// ---------------------------------------------------------------------------

// `subdomain` selects which fine subdomain's DATA this rank carries. The `fine`
// LABEL is always the rank: solveDistributedHierarchy validates `input.fine ==
// rank` (LadrunoCMSHierarchy.cpp:1293), so the label is not a free index and a
// permuted assignment can only be expressed by moving the data, which is what a
// different partitioner would hand you anyway.
bool runDistributedHierarchy(
    int rank,
    int subdomain,
    int coarse,
    ladruno_cms::DistributedHierarchyResult &result,
    std::string &message)
{
    Fixture fixture = makeFixture(true);
    const ladruno_cms::FineSubdomainInput &local =
        fixture.input.fineSubdomains[static_cast<std::size_t>(subdomain)];
    const std::vector<double> stiffnessDense = local.stiffness.toDense();
    const std::vector<double> massDense = local.mass.toDense();
    ladruno_cms::AssemblyRecord stiffnessRecord;
    ladruno_cms::AssemblyRecord massRecord;
    if (ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Stiffness,
            static_cast<std::size_t>(subdomain), stiffnessDense.data(), 3, 3,
            local.equations, 9, 1.0, 1.0e-12, 1.0e-14, stiffnessRecord,
            message) != 0)
        return false;
    if (ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Mass,
            static_cast<std::size_t>(subdomain), massDense.data(), 3, 3,
            local.equations, 9, 1.0, 1.0e-12, 1.0e-14, massRecord,
            message) != 0)
        return false;

    ladruno_cms::DistributedHierarchyInput input;
    input.globalDimension = 9;
    input.fine = rank;
    input.coarse = coarse;
    input.numberOfCoarseSubdomains = 2;
    input.localEquations = local.equations;
    input.localStiffness = local.stiffness;
    input.localMass = local.mass;
    std::vector<ladruno_cms::AssemblyRecord> ownedStiffness{stiffnessRecord};
    std::vector<ladruno_cms::AssemblyRecord> ownedMass{massRecord};
    input.ownedStiffnessContributions = &ownedStiffness;
    input.ownedMassContributions = &ownedMass;
    input.modesLevel2 = 3;
    input.modesLevel1 = 10;
    input.numberOfModes = 5;
    input.denseMax = 20;
    input.tolerance = 1.0e-8;
    input.maximumOperatorApplications = 300;
    message.clear();
    return ladruno_cms::solveDistributedHierarchy(input, result, message) == 0;
}

double modalAssuranceCriterion(
    const std::vector<double> &left,
    const std::vector<double> &right,
    int dimension,
    int mode)
{
    const std::size_t offset = static_cast<std::size_t>(mode) * dimension;
    double cross = 0.0;
    double leftNorm = 0.0;
    double rightNorm = 0.0;
    for (int row = 0; row < dimension; ++row) {
        const double a = left[offset + static_cast<std::size_t>(row)];
        const double b = right[offset + static_cast<std::size_t>(row)];
        cross += a * b;
        leftNorm += a * a;
        rightNorm += b * b;
    }
    if (leftNorm <= 0.0 || rightNorm <= 0.0)
        return 0.0;
    return (cross * cross) / (leftNorm * rightNorm);
}

void checkRankPermutationInvariance(int rank, int size)
{
    if (size != 4) {
        if (rank == 0)
            std::cout << "rank-permutation invariance SKIPPED: the fixture has "
                         "4 fine subdomains and needs exactly 4 ranks, got "
                      << size << '\n';
        return;
    }
    // The coarse grouping stays structural (ranks 0,1 -> group 0; 2,3 -> group
    // 1); what the permutation changes is WHICH subdomain's data lands in which
    // group. With a complete local basis both runs are exact, so the spectra and
    // the subspace must agree regardless.
    ladruno_cms::DistributedHierarchyResult baseline;
    std::string message;
    // NOTE: run first, THEN build the failure string. Passing the call and
    // "..." + message as two arguments of require() leaves their evaluation
    // order unspecified, and the message was being captured BEFORE the call ran
    // -- which reported every failure with an empty diagnostic.
    const bool baselineOk =
        runDistributedHierarchy(rank, rank, rank / 2, baseline, message);
    require(baselineOk, "baseline distributed hierarchy failed: " + message);

    const int permutations[2][4] = {
        {3, 2, 1, 0},   // reversal
        {0, 2, 1, 3}};  // interleave: group 0 now holds subdomains {0,2}
    const char *labels[2] = {"reversal", "interleave"};

    for (int which = 0; which < 2; ++which) {
        const int subdomain = permutations[which][rank];
        const std::string label(labels[which]);
        ladruno_cms::DistributedHierarchyResult permuted;
        message.clear();
        const bool permutedOk = runDistributedHierarchy(
            rank, subdomain, rank / 2, permuted, message);
        require(permutedOk,
                "permuted (" + label + ") distributed hierarchy failed: " +
                    message);
        if (!permutedOk ||
            permuted.eigenvalues.size() != baseline.eigenvalues.size())
            continue;

        // Teeth check: the permutation must actually have moved data, or the
        // comparisons below are comparing a run with itself.
        const Fixture layout = makeFixture(true);
        const int moved =
            layout.input.fineSubdomains[static_cast<std::size_t>(subdomain)]
                        .equations !=
                    layout.input.fineSubdomains[static_cast<std::size_t>(rank)]
                        .equations
                ? 1
                : 0;
        int movedAnywhere = 0;
        MPI_Allreduce(&moved, &movedAnywhere, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
        require(movedAnywhere == 1,
                "permutation (" + label + ") left every rank holding its "
                "original subdomain -- the invariance check is vacuous");

        for (std::size_t mode = 0; mode < baseline.eigenvalues.size(); ++mode) {
            const double reference = baseline.eigenvalues[mode];
            const double candidate = permuted.eigenvalues[mode];
            const double relative = std::fabs(candidate - reference) /
                std::max(1.0, std::fabs(reference));
            require(relative <= 1.0e-9,
                    "permutation (" + label + ") moved eigenvalue " +
                        std::to_string(mode) + " by a relative " +
                        std::to_string(relative));
            require(permuted.normalizedResiduals[mode] < 3.0e-8,
                    "permutation (" + label + ") damaged the residual of mode " +
                        std::to_string(mode));
            const double mac = modalAssuranceCriterion(
                baseline.eigenvectors, permuted.eigenvectors, 9,
                static_cast<int>(mode));
            require(mac >= 1.0 - 1.0e-9,
                    "permutation (" + label + ") rotated the subspace: MAC of "
                    "mode " + std::to_string(mode) + " is " +
                        std::to_string(mac));
        }
        require(permuted.diagnostics.finalRawDimension ==
                    baseline.diagnostics.finalRawDimension,
                "permutation (" + label + ") changed the reduced dimension");
    }
}

// ---------------------------------------------------------------------------
// P3d -- the three interface topologies.
//
// The shared fixture is one shape: a chain whose subdomains share equations both
// INSIDE a coarse group (a level-2 interface) and ACROSS groups (a level-1
// interface). That exercises S2 and S1 together and can hide a compatibility
// path that only works when the other one is also active. The two degenerate
// shapes are built here as a SEPARATE fixture so the shared makeFixture -- and
// the seven tests pinned to its 9/12/10/10 dimensions -- stay untouched.
//
//   combined    subdomains {0,1,2} {2,3,4} {4,5,6} {6,7,8}      order 9
//               shares 2 (in group 0), 4 (across groups), 6 (in group 1)
//   level2Only  {0,1,2} {2,3,4} | {5,6,7} {7,8,9}               order 10
//               shares only within each coarse group; the groups are decoupled
//   level1Only  {0,1,2} {3,4,5} | {2,6,7} {5,8,9}               order 10
//               no sharing inside a group; groups couple through 2 and 5
// ---------------------------------------------------------------------------

struct TopologyFixture {
    int order = 0;
    std::vector<std::vector<int>> equations;
    std::vector<ladruno_cms::SymmetricCSR> stiffness;
    std::vector<ladruno_cms::SymmetricCSR> mass;
    ladruno_cms::SymmetricCSR globalStiffness;
    ladruno_cms::SymmetricCSR globalMass;
};

TopologyFixture makeTopologyFixture(int topology)
{
    static const int tables[3][4][3] = {
        {{0, 1, 2}, {2, 3, 4}, {4, 5, 6}, {6, 7, 8}},   // combined
        {{0, 1, 2}, {2, 3, 4}, {5, 6, 7}, {7, 8, 9}},   // level-2 only
        {{0, 1, 2}, {3, 4, 5}, {2, 6, 7}, {5, 8, 9}}};  // level-1 only
    TopologyFixture fixture;
    fixture.order = topology == 0 ? 9 : 10;
    std::vector<double> globalStiffness(
        static_cast<std::size_t>(fixture.order) * fixture.order, 0.0);
    std::vector<double> globalMass(globalStiffness.size(), 0.0);
    for (int fine = 0; fine < 4; ++fine) {
        const double shift = 0.15 * fine;
        const std::vector<double> localStiffness = {
            2.2 + shift, -1.0, 0.0,
            -1.0, 2.6 + shift, -0.8,
            0.0, -0.8, 1.9 + shift};
        const std::vector<double> localMass = {
            0.8 + 0.05 * fine, 0.04, 0.0,
            0.04, 1.0 + 0.05 * fine, 0.03,
            0.0, 0.03, 1.2 + 0.05 * fine};
        std::vector<int> ids(tables[topology][fine], tables[topology][fine] + 3);
        fixture.equations.push_back(ids);
        fixture.stiffness.push_back(csrFromDense(localStiffness, 3));
        fixture.mass.push_back(csrFromDense(localMass, 3));
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                const std::size_t global =
                    static_cast<std::size_t>(ids[static_cast<std::size_t>(row)]) *
                        fixture.order +
                    ids[static_cast<std::size_t>(column)];
                globalStiffness[global] +=
                    localStiffness[static_cast<std::size_t>(row) * 3 + column];
                globalMass[global] +=
                    localMass[static_cast<std::size_t>(row) * 3 + column];
            }
        }
    }
    fixture.globalStiffness = csrFromDense(globalStiffness, fixture.order);
    fixture.globalMass = csrFromDense(globalMass, fixture.order);
    return fixture;
}

void checkInterfaceTopologies(int rank, int size)
{
    if (size != 4) {
        if (rank == 0)
            std::cout << "interface-topology sweep SKIPPED: needs exactly 4 "
                         "ranks, got " << size << '\n';
        return;
    }
    const char *labels[3] = {"combined", "level2Only", "level1Only"};
    for (int topology = 0; topology < 3; ++topology) {
        const std::string label(labels[topology]);
        const TopologyFixture fixture = makeTopologyFixture(topology);
        const std::vector<double> reference =
            directEigenvalues(fixture.globalStiffness, fixture.globalMass, 5);
        const std::vector<double> stiffnessDense =
            fixture.stiffness[static_cast<std::size_t>(rank)].toDense();
        const std::vector<double> massDense =
            fixture.mass[static_cast<std::size_t>(rank)].toDense();

        std::string message;
        ladruno_cms::AssemblyRecord stiffnessRecord;
        ladruno_cms::AssemblyRecord massRecord;
        const bool stiffnessOk = ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Stiffness,
            static_cast<std::size_t>(rank), stiffnessDense.data(), 3, 3,
            fixture.equations[static_cast<std::size_t>(rank)], fixture.order,
            1.0, 1.0e-12, 1.0e-14, stiffnessRecord, message) == 0;
        require(stiffnessOk,
                "topology " + label + ": stiffness record failed: " + message);
        message.clear();
        const bool massOk = ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Mass,
            static_cast<std::size_t>(rank), massDense.data(), 3, 3,
            fixture.equations[static_cast<std::size_t>(rank)], fixture.order,
            1.0, 1.0e-12, 1.0e-14, massRecord, message) == 0;
        require(massOk, "topology " + label + ": mass record failed: " + message);
        if (!stiffnessOk || !massOk)
            continue;

        ladruno_cms::DistributedHierarchyInput input;
        input.globalDimension = fixture.order;
        input.fine = rank;
        input.coarse = rank / 2;
        input.numberOfCoarseSubdomains = 2;
        input.localEquations = fixture.equations[static_cast<std::size_t>(rank)];
        input.localStiffness = fixture.stiffness[static_cast<std::size_t>(rank)];
        input.localMass = fixture.mass[static_cast<std::size_t>(rank)];
        std::vector<ladruno_cms::AssemblyRecord> ownedStiffness{stiffnessRecord};
        std::vector<ladruno_cms::AssemblyRecord> ownedMass{massRecord};
        input.ownedStiffnessContributions = &ownedStiffness;
        input.ownedMassContributions = &ownedMass;
        input.modesLevel2 = 3;
        input.modesLevel1 = 10;
        input.numberOfModes = 5;
        input.denseMax = 20;
        input.tolerance = 1.0e-8;
        input.maximumOperatorApplications = 300;

        ladruno_cms::DistributedHierarchyResult result;
        message.clear();
        const bool solved =
            ladruno_cms::solveDistributedHierarchy(input, result, message) == 0;
        require(solved,
                "topology " + label + ": distributed hierarchy failed: " +
                    message);
        if (!solved)
            continue;

        // The T2 -> S2 -> T1 -> S1 chain is mandatory per the ADR, including on
        // a topology where one of the interface sets is empty.
        require(result.diagnostics.appliedT2 && result.diagnostics.appliedS2 &&
                    result.diagnostics.appliedT1 && result.diagnostics.appliedS1,
                "topology " + label + ": the mandatory transformation chain was "
                "not applied end to end");
        require(result.eigenvectors.size() ==
                    static_cast<std::size_t>(fixture.order) * 5,
                "topology " + label + ": incomplete published eigenvectors");
        require(result.diagnostics.maximumDuplicateJump < 2.0e-10,
                "topology " + label + ": incompatible interface copies");
        for (std::size_t mode = 0; mode < reference.size(); ++mode) {
            require(close(result.eigenvalues[mode], reference[mode], 3.0e-8),
                    "topology " + label + ": eigenvalue " +
                        std::to_string(mode) + " differs from direct LAPACK");
            require(result.normalizedResiduals[mode] < 3.0e-8,
                    "topology " + label + ": residual of mode " +
                        std::to_string(mode) + " is excessive");
        }
    }
}

// ---------------------------------------------------------------------------
// P3d -- the level-1 ablation, allowed ONLY as a labelled diagnostic.
//
// ADR-1000 P4 section 2b wants a "level-2 only" ablation (omit T1) to attribute
// cost and reduction between the two levels, and P3d requires that it can never
// be an accepted configuration or a fallback. Three things are checked:
//   1. it is a REAL ablation -- the final space is strictly larger than the full
//      chain's on a fixture where T1 actually truncates (else it is a no-op
//      dressed up as a benchmark);
//   2. it is LABELLED -- appliedT1 == false and ablatedLevel1 == true, so no
//      consumer can mistake it for the mandatory chain;
//   3. it is UNREACHABLE from the production path -- the distributed hierarchy
//      the solver actually calls always reports the full chain.
// The fourth guard (the option parser refuses an ablation flag) lives in
// ladruno_cms_core_check, next to the rest of the parser tests.
// ---------------------------------------------------------------------------

void checkLevelAblationIsDiagnosticOnly()
{
    // A truncating fixture: modesLevel1 = 2, so T1 genuinely removes coordinates
    // and the ablation must show up as a larger final space.
    Fixture full = makeFixture(false);
    ladruno_cms::TwoLevelHierarchyResult fullResult;
    std::string message;
    REQUIRE_CALL(
        ladruno_cms::solveTwoLevelHierarchy(full.input, fullResult, message) == 0,
        "full-chain hierarchy failed: " + message);

    Fixture ablated = makeFixture(false);
    ablated.input.diagnosticAblateLevel1 = true;
    ladruno_cms::TwoLevelHierarchyResult ablatedResult;
    message.clear();
    REQUIRE_CALL(
        ladruno_cms::solveTwoLevelHierarchy(
            ablated.input, ablatedResult, message) == 0,
        "level-1 ablation failed to run as a diagnostic: " + message);

    std::cout << "level-1 ablation diagnostic: finalRawDimension full="
              << fullResult.diagnostics.finalRawDimension << " ablated="
              << ablatedResult.diagnostics.finalRawDimension
              << " (T1 removes "
              << ablatedResult.diagnostics.finalRawDimension -
                     fullResult.diagnostics.finalRawDimension
              << " coordinates)\n";

    // 1. a real ablation, not a no-op
    require(ablatedResult.diagnostics.finalRawDimension >
                fullResult.diagnostics.finalRawDimension,
            "the level-1 ablation did not enlarge the final space (" +
                std::to_string(ablatedResult.diagnostics.finalRawDimension) +
                " vs " +
                std::to_string(fullResult.diagnostics.finalRawDimension) +
                ") -- it is a no-op, so it cannot attribute anything");

    // 2. labelled, and only the level-1 stage is skipped
    require(!ablatedResult.diagnostics.appliedT1,
            "the ablated run still reports appliedT1");
    require(ablatedResult.diagnostics.ablatedLevel1,
            "the ablated run is not labelled as ablated");
    require(ablatedResult.diagnostics.appliedT2 &&
                ablatedResult.diagnostics.appliedS2 &&
                ablatedResult.diagnostics.appliedS1,
            "the ablation skipped more than the level-1 reduction");
    require(fullResult.diagnostics.appliedT1 &&
                !fullResult.diagnostics.ablatedLevel1,
            "the full chain mislabelled itself as ablated");

    // The diagnostic must still be a usable measurement: a larger retained space
    // cannot be less accurate than the truncated one.
    const std::size_t modes = std::min(ablatedResult.eigenvalues.size(),
                                       fullResult.eigenvalues.size());
    for (std::size_t mode = 0; mode < modes; ++mode) {
        const double allowance = 1.0e-9 *
            std::max(1.0, std::fabs(fullResult.eigenvalues[mode]));
        require(ablatedResult.eigenvalues[mode] <=
                    fullResult.eigenvalues[mode] + allowance,
                "the ablated (larger) space produced a HIGHER Ritz value at mode " +
                    std::to_string(mode) + " -- the reduction is inconsistent");
    }

    // 3. default-off: nothing gets ablated by accident
    const ladruno_cms::TwoLevelHierarchyInput defaults;
    require(!defaults.diagnosticAblateLevel1,
            "level-1 ablation is ON by default -- it must be opt-in");
}

void checkProductionPathIsNeverAblated(int rank, int size)
{
    if (size != 4)
        return;
    ladruno_cms::DistributedHierarchyResult result;
    std::string message;
    const bool ok = runDistributedHierarchy(rank, rank, rank / 2, result, message);
    require(ok, "distributed hierarchy failed: " + message);
    if (!ok)
        return;
    // The distributed path is the one LadrunoCMSEigenSolver calls. It has no
    // ablation input at all, and the solver refuses a result that claims one.
    require(!result.diagnostics.ablatedLevel1,
            "the production distributed hierarchy reported an ablated level 1");
    require(result.diagnostics.appliedT1 && result.diagnostics.appliedS1,
            "the production distributed hierarchy skipped the mandatory chain");
}

void checkAssemblyMismatchIsRejected()
{
    Fixture fixture = makeFixture(true);
    fixture.input.fineSubdomains[0].stiffness.upperValues[0] += 0.1;
    ladruno_cms::TwoLevelHierarchyResult result;
    std::string message;
    require(ladruno_cms::solveTwoLevelHierarchy(fixture.input, result, message) < 0,
            "inconsistent local/global assembly was accepted");
    require(message.find("do not reproduce") != std::string::npos,
            "assembly mismatch did not produce an actionable diagnostic");
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    checkCompleteTwoLevelFlow();
    checkTruncatedTwoLevelFlow();
    checkCoordinateMassCondensationAcrossHierarchy();
    checkExplicitMemoryGuards();
    checkAssemblyMismatchIsRejected();
    checkDistributedFourRankFlow(rank, size);
    checkRankPermutationInvariance(rank, size);
    checkInterfaceTopologies(rank, size);
    if (rank == 0)
        checkLevelAblationIsDiagnosticOnly();
    checkProductionPathIsNeverAblated(rank, size);
    MPI_Finalize();
    if (failures != 0)
        return 1;
    std::cout << "two-level hierarchy checks passed\n";
    return 0;
}
