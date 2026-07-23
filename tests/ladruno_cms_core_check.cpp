#include "LadrunoCMSOptions.h"

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

int failures = 0;

void require(bool condition, const std::string &message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool close(double left, double right, double tolerance = 1.0e-13)
{
    return std::fabs(left - right) <= tolerance;
}

ladruno_cms::AssemblyRecord record(
    ladruno_cms::ContributionKind kind,
    std::vector<int> ids,
    std::vector<double> upper)
{
    ladruno_cms::AssemblyRecord value;
    value.kind = kind;
    value.activeIDs = std::move(ids);
    value.upperValues = std::move(upper);
    return value;
}

void testOptions()
{
    ladruno_cms::Options options;
    options.level1 = 2;
    options.level2 = 2;
    options.modesLevel2 = 3;
    options.modesLevel1 = 2;
    std::string message;
    require(options.validate(4, 1, message) == 0, "valid logical hierarchy rejected");
    require(options.validate(2, 1, message) < 0, "world-size mismatch accepted");
    options.domainMode = ladruno_cms::DomainMode::Physical;
    options.verifyAssembly = ladruno_cms::AssemblyVerification::Off;
    require(options.validate(4, 1, message) == 0, "valid physical mode rejected");
    options.verifyAssembly = ladruno_cms::AssemblyVerification::Signature;
    require(options.validate(4, 1, message) < 0,
            "physical mode accepted replicated assembly verification");
}

void testAssemblySignatures()
{
    const double values[] = {1.0, -2.0, -2.0, 5.0};
    const double nearValues[] = {1.0, -2.0, -2.0, 5.0 + 1.0e-12};
    const double farValues[] = {1.0, -2.0, -2.0, 5.0 + 1.0e-4};
    ladruno_cms::AssemblyRecord reference;
    reference.ordinal = 7;
    reference.rows = 2;
    reference.columns = 2;
    reference.originalIDs = {-1, 3};
    reference.activeIDs = {3};
    reference.metrics = ladruno_cms::computeAssemblyMetrics(values, 2, 2, 0.5);

    ladruno_cms::AssemblyRecord candidate = reference;
    candidate.metrics = ladruno_cms::computeAssemblyMetrics(nearValues, 2, 2, 0.5);
    std::string message;
    require(
        ladruno_cms::compareAssemblyRecords(reference, candidate, 1.0e-10, 1.0e-13, message),
        "numerically equivalent assembly signature rejected");
    candidate.metrics = ladruno_cms::computeAssemblyMetrics(farValues, 2, 2, 0.5);
    require(
        !ladruno_cms::compareAssemblyRecords(reference, candidate, 1.0e-10, 1.0e-13, message),
        "different assembly values accepted");
    candidate = reference;
    candidate.originalIDs = {-1, 4};
    require(
        !ladruno_cms::compareAssemblyRecords(reference, candidate, 1.0e-10, 1.0e-13, message),
        "different assembly IDs accepted");

    const double filteredValues[] = {
        10.0, 7.0, 8.0,
        9.0, 2.0, -1.0,
        6.0, -1.0, 3.0};
    ladruno_cms::AssemblyRecord filtered;
    require(
        ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Mass,
            4,
            filteredValues,
            3,
            3,
            {-1, 2, 2},
            3,
            0.5,
            1.0e-12,
            1.0e-14,
            filtered,
            message) == 0,
        "assembly record construction failed");
    require(
        filtered.originalIDs == std::vector<int>({-1, 2, 2}) &&
            filtered.activeIDs == std::vector<int>({2, 2}),
        "negative-ID filtering or repeated IDs are incorrect");
    require(
        filtered.upperValues.size() == 3u && close(filtered.upperValues[0], 1.0) &&
            close(filtered.upperValues[1], -0.5) && close(filtered.upperValues[2], 1.5),
        "factor or packed active contribution is incorrect");

    const double nonsymmetric[] = {1.0, 2.0, 3.0, 4.0};
    require(
        ladruno_cms::makeAssemblyRecord(
            ladruno_cms::ContributionKind::Stiffness,
            0,
            nonsymmetric,
            2,
            2,
            {0, 1},
            2,
            1.0,
            1.0e-12,
            1.0e-14,
            filtered,
            message) < 0,
        "non-symmetric active contribution accepted");
}

void testOwnership()
{
    std::string message;
    require(
        ladruno_cms::assignContributionOwner({0, 1, 2}, {2, 1, 1}, message) == 1,
        "majority owner is incorrect");
    require(
        ladruno_cms::assignContributionOwner({0, 1}, {1, 0}, message) == 0,
        "owner tie did not select the lowest label");
    require(
        ladruno_cms::assignContributionOwner({3}, {0, 1}, message) < 0,
        "out-of-range active equation accepted");

    const std::vector<ladruno_cms::AssemblyRecord> contributions = {
        record(ladruno_cms::ContributionKind::Stiffness, {0, 1}, {}),
        record(ladruno_cms::ContributionKind::Stiffness, {1, 2}, {}),
        record(ladruno_cms::ContributionKind::Stiffness, {2}, {})};
    ladruno_cms::EffectiveOwnership ownership;
    require(
        ladruno_cms::classifyEffectiveOwnership(
            3, contributions, {0, 1, 2}, {0, 0, 1}, ownership, message) == 0,
        "effective ownership classification failed");
    require(
        ownership.role[0] == ladruno_cms::EquationRole::Interior &&
            ownership.interiorOwner[0] == 0,
        "single-owner coordinate is not interior");
    require(
        ownership.role[1] == ladruno_cms::EquationRole::Level2Interface,
        "same-coarse shared coordinate is not an S2 interface");
    require(
        ownership.role[2] == ladruno_cms::EquationRole::Level1Interface,
        "cross-coarse shared coordinate is not an S1 interface");
}

void testCSR()
{
    const std::vector<ladruno_cms::AssemblyRecord> contributions = {
        record(ladruno_cms::ContributionKind::Stiffness, {0, 1}, {2.0, -1.0, 2.0}),
        record(ladruno_cms::ContributionKind::Stiffness, {1, 2}, {3.0, -2.0, 4.0})};
    ladruno_cms::SymmetricCSR csr;
    std::string message;
    require(
        ladruno_cms::buildSymmetricCSR(
            3, contributions, ladruno_cms::ContributionKind::Stiffness, csr, message) == 0,
        "CSR assembly failed");
    const std::vector<double> dense = csr.toDense();
    const std::vector<double> expected = {
        2.0, -1.0, 0.0,
        -1.0, 5.0, -2.0,
        0.0, -2.0, 4.0};
    require(dense.size() == expected.size(), "dense CSR size is incorrect");
    for (std::size_t index = 0; index < dense.size() && index < expected.size(); ++index)
        require(close(dense[index], expected[index]), "CSR does not match dense assembly");
    const std::vector<double> product = csr.multiply({1.0, 2.0, 3.0});
    require(
        product.size() == 3u && close(product[0], 0.0) &&
            close(product[1], 3.0) && close(product[2], 8.0),
        "symmetric CSR multiplication is incorrect");

    ladruno_cms::SymmetricCSR repeated;
    require(
        ladruno_cms::buildSymmetricCSR(
            1,
            {record(ladruno_cms::ContributionKind::Stiffness, {0, 0}, {1.0, 2.0, 3.0})},
            ladruno_cms::ContributionKind::Stiffness,
            repeated,
            message) == 0,
        "repeated-ID CSR assembly failed");
    require(
        repeated.toDense().size() == 1u && close(repeated.toDense()[0], 8.0),
        "repeated IDs lost an off-diagonal contribution");
}

void testCoordinateMass()
{
    std::string message;
    ladruno_cms::CoordinateMassClassification classification;
    ladruno_cms::SymmetricCSR mass;
    mass.dimension = 3;
    mass.rowOffsets = {0, 1, 2, 3};
    mass.columnIndices = {0, 1, 2};
    mass.upperValues = {2.0, 0.0, 3.0};
    require(
        ladruno_cms::classifyCoordinateMass(
            mass, 1.0e-12, 1.0e-14, classification, message) == 0,
        "coordinate-aligned semidefinite mass rejected");
    require(
        classification.dynamic == std::vector<int>({0, 2}) &&
            classification.massless == std::vector<int>({1}),
        "dynamic/massless coordinate sets are incorrect");

    mass.rowOffsets = {0, 2, 3, 4};
    mass.columnIndices = {0, 1, 1, 2};
    mass.upperValues = {2.0, 1.0e-2, 0.0, 3.0};
    require(
        ladruno_cms::classifyCoordinateMass(
            mass, 1.0e-12, 1.0e-14, classification, message) < 0,
        "non-coordinate-aligned mass nullspace accepted");

    mass.rowOffsets = {0, 1, 2, 3};
    mass.columnIndices = {0, 1, 2};
    mass.upperValues = {2.0, -1.0, 3.0};
    require(
        ladruno_cms::classifyCoordinateMass(
            mass, 1.0e-12, 1.0e-14, classification, message) < 0,
        "negative mass diagonal accepted");

    mass.upperValues = {0.0, 0.0, 0.0};
    require(
        ladruno_cms::classifyCoordinateMass(
            mass, 1.0e-12, 1.0e-14, classification, message) < 0,
        "mass without dynamic coordinates accepted");
}

void testCommandParserAndLocalPencil()
{
    ladruno_cms::Options options;
    int modes = 0;
    std::string message;
    const std::vector<std::string> valid = {
        "-hierarchy", "logical", "-level1", "2", "-level2", "2",
        "-modesL2", "8", "-modesL1", "12", "-tol", "1e-8",
        "-refinement", "subspace", "-iterationVectors", "auto",
        "-maxRefineIter", "12",
        "-partition", "metis", "-localEigen", "lanczos",
        "-globalSolver", "dense", "-denseMax", "50", "-verbose", "5"};
    require(ladruno_cms::parseCommandOptions(valid, options, modes, message) == 0,
            "valid command options were rejected: " + message);
    require(modes == 5 && options.level1 == 2 && options.level2 == 2 &&
            options.modesLevel2 == 8 && options.modesLevel1 == 12 && options.verbose,
            "valid command options were parsed incorrectly");
    require(options.resolvedIterationVectors(modes, message) == 13 &&
            options.maxRefineIterations == 12,
            "subspace-refinement options were parsed incorrectly");
    require(options.domainMode == ladruno_cms::DomainMode::ReplicatedReference,
            "default command mode no longer preserves the P2 reference");
    require(ladruno_cms::parseCommandOptions(
                {"-domainMode", "physical", "-hierarchy", "logical",
                 "-level1", "2", "-level2", "2", "-modesL2", "4",
                 "-modesL1", "4", "2"},
                options, modes, message) == 0 &&
            options.domainMode == ladruno_cms::DomainMode::Physical &&
            options.verifyAssembly == ladruno_cms::AssemblyVerification::Off,
            "physical command mode was not parsed with local assembly semantics");
    require(ladruno_cms::parseCommandOptions(
                {"-domainMode", "physical", "-hierarchy", "logical",
                 "-level1", "2", "-level2", "2", "-modesL2", "4",
                 "-modesL1", "4", "-verifyAssembly", "signature", "2"},
                options, modes, message) < 0,
            "physical mode accepted replicated signature verification");
    require(ladruno_cms::parseCommandOptions(
                {"-hierarchy", "logical", "-level1", "2", "-level1", "2",
                 "-level2", "2", "-modesL2", "4", "-modesL1", "4", "2"},
                options, modes, message) < 0,
            "duplicate command option was accepted");
    require(ladruno_cms::parseCommandOptions(
                {"-hierarchy", "auto", "-level1", "1", "-modesL2", "4",
                 "-modesL1", "4", "2"},
                options, modes, message) < 0,
            "auto hierarchy accepted explicit topology values");
    require(ladruno_cms::parseCommandOptions(
                {"-hierarchy", "auto", "-modesL2", "4", "-modesL1", "4",
                 "-globalSolver", "unknown", "2"},
                options, modes, message) < 0,
            "unsupported final backend was accepted");

    const std::vector<ladruno_cms::AssemblyRecord> stiffness = {
        record(ladruno_cms::ContributionKind::Stiffness, {5, 2}, {2.0, -1.0, 3.0})};
    const std::vector<ladruno_cms::AssemblyRecord> mass = {
        record(ladruno_cms::ContributionKind::Mass, {5, 2}, {1.0, 0.0, 2.0})};
    std::vector<int> equations;
    ladruno_cms::SymmetricCSR localStiffness;
    ladruno_cms::SymmetricCSR localMass;
    require(ladruno_cms::buildOwnedLocalPencil(
                stiffness, mass, equations, localStiffness, localMass, message) == 0,
            "owned local pencil construction failed: " + message);
    require(equations == std::vector<int>({2, 5}),
            "owned local equations are not stable and sorted");
    const std::vector<double> dense = localStiffness.toDense();
    require(dense.size() == 4u && close(dense[0], 3.0) && close(dense[1], -1.0) &&
            close(dense[2], -1.0) && close(dense[3], 2.0),
            "global-to-local contribution remapping changed the pencil");
}

} // namespace

int main()
{
    testOptions();
    testAssemblySignatures();
    testOwnership();
    testCSR();
    testCoordinateMass();
    testCommandParserAndLocalPencil();
    if (failures != 0)
        return 1;
    std::cout << "Ladruno CMS core checks passed\n";
    return 0;
}
