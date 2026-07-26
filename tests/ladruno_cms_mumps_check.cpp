#include "LadrunoCMSMumps.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void require(bool condition, const char *message)
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

void require(bool condition, const std::string &message)
{
    require(condition, message.c_str());
}

bool close(double left, double right, double tolerance = 1.0e-11)
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
            if (value == 0.0)
                continue;
            result.columnIndices.push_back(column);
            result.upperValues.push_back(value);
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    return result;
}

std::vector<double> multiplyDense(
    const std::vector<double> &matrix,
    int order,
    const std::vector<double> &vectors,
    int numberOfColumns)
{
    std::vector<double> result(
        static_cast<std::size_t>(order) * numberOfColumns, 0.0);
    for (int column = 0; column < numberOfColumns; ++column) {
        for (int row = 0; row < order; ++row) {
            for (int inner = 0; inner < order; ++inner)
                result[static_cast<std::size_t>(column) * order + row] +=
                    matrix[static_cast<std::size_t>(row) * order + inner] *
                    vectors[static_cast<std::size_t>(column) * order + inner];
        }
    }
    return result;
}

void testMumpsMultiRHS()
{
    const std::vector<double> dense = {4.0, 1.0, 1.0, 3.0};
    const ladruno_cms::SymmetricCSR matrix = csrFromDense(dense, 2);
    ladruno_cms::MumpsSPD solver;
    std::string message;
    require(solver.factorize(matrix, message) == 0, "SPD MUMPS factorization failed");
    require(solver.isFactorized() && solver.dimension() == 2, "MUMPS state is inconsistent");
    const std::vector<double> rightHandSides = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> solution = rightHandSides;
    require(solver.solve(solution, 2, message) == 0, "MUMPS multi-RHS solve failed");
    const std::vector<double> recovered = multiplyDense(dense, 2, solution, 2);
    for (std::size_t index = 0; index < recovered.size(); ++index)
        require(close(recovered[index], rightHandSides[index]), "MUMPS solve residual is too large");

    ladruno_cms::MumpsSPD indefiniteSolver;
    require(
        indefiniteSolver.factorize(csrFromDense({1.0, 0.0, 0.0, -1.0}, 2), message) < 0,
        "MUMPS SPD guard accepted a negative pivot");
}

// Ranks 0 and 1 carry every entry, so this fixture is valid for any size >= 2;
// it used to be gated to size == 4, which silently skipped the whole
// distributed path at np=2 (ADR-1000 section 20.5).
void testDistributedMumps(int rank, int size)
{
    if (size != 2 && size != 4)
        return;
    ladruno_cms::SymmetricCSR local;
    local.dimension = 2;
    local.rowOffsets = {0, 0, 0};
    if (rank == 0) {
        local.rowOffsets = {0, 2, 2};
        local.columnIndices = {0, 1};
        local.upperValues = {4.0, 1.0};
    } else if (rank == 1) {
        local.rowOffsets = {0, 0, 1};
        local.columnIndices = {1};
        local.upperValues = {3.0};
    }
    ladruno_cms::DistributedMumpsSPD solver;
    std::string message;
    require(solver.factorize(local, message) == 0,
            "distributed SPD MUMPS factorization failed");
    std::vector<double> rightHandSides = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> solution = rightHandSides;
    require(solver.solve(solution, 2, message) == 0,
            "distributed MUMPS multi-RHS solve failed");
    const std::vector<double> recovered = multiplyDense(
        {4.0, 1.0, 1.0, 3.0}, 2, solution, 2);
    for (std::size_t index = 0; index < recovered.size(); ++index)
        require(close(recovered[index], rightHandSides[index]),
                "distributed MUMPS solve residual is too large");
}

// ---------------------------------------------------------------------------
// Part-0 fixture (ADR-1000 P3 execution plan, section 20.5 of the ADR).
//
// The fixture above is a 2x2 pencil: it proves the distributed wiring, but
// MUMPS never leaves its trivial path on it, so it cannot see the np=2 analysis
// failure Part 0 is about ("orderMinPriority: no valid number of stages in
// multisector", 15.6 GiB spent in ANALYSIS, np=4 and np=6 unaffected). Part 0
// needs a distributed pencil big enough that MUMPS makes a real ordering
// decision, run at BOTH np=2 (must pass) and np=4 (must not regress).
//
// Two shapes, because the reported failure came from a rank-local CMS pencil
// that is DENSE-as-CSR, not sparse (see the P4 plan's memory-debt section):
//   * Laplace -- a 3-D 7-point Laplacian, genuinely sparse;
//   * Dense   -- a diagonally dominant full block, the shape reduceCraigBampton
//                actually hands MUMPS today.
// Both are SPD by strict diagonal dominance. Only the PATTERN drives the
// analysis phase, so the values are kept trivially well conditioned.
//
// Sizes are overridable (LADRUNO_CMS_CHECK_SIDE, LADRUNO_CMS_CHECK_DENSE_ORDER)
// so the same binary can be cranked up to hunt an ordering failure without
// making the default run expensive. Defaults are CI-cheap.
// ---------------------------------------------------------------------------

int environmentSize(const char *name, int fallback)
{
    const char *raw = std::getenv(name);
    if (raw == nullptr)
        return fallback;
    const long parsed = std::strtol(raw, nullptr, 10);
    if (parsed < 2 || parsed > 1000000)
        return fallback;
    return static_cast<int>(parsed);
}

struct ModelMatrix {
    bool dense = false;
    int side = 0;   // grid side of the Laplacian; unused when dense
    int order = 0;

    // Upper-triangular entries of `row` (columns >= row), ascending.
    void upperRow(int row, std::vector<int> &columns,
                  std::vector<double> &values) const
    {
        columns.clear();
        values.clear();
        columns.push_back(row);
        if (dense) {
            values.push_back(static_cast<double>(order) + 1.0);
            for (int column = row + 1; column < order; ++column) {
                columns.push_back(column);
                values.push_back(
                    1.0 / (1.0 + static_cast<double>(column - row)));
            }
            return;
        }
        values.push_back(6.0);
        const int x = row % side;
        const int y = (row / side) % side;
        const int z = row / (side * side);
        if (x + 1 < side) {
            columns.push_back(row + 1);
            values.push_back(-1.0);
        }
        if (y + 1 < side) {
            columns.push_back(row + side);
            values.push_back(-1.0);
        }
        if (z + 1 < side) {
            columns.push_back(row + side * side);
            values.push_back(-1.0);
        }
    }
};

// This rank's contiguous row block, in the distributed-assembled form
// factorizeDistributed expects: GLOBAL dimension, only owned rows populated.
ladruno_cms::SymmetricCSR buildLocalBlock(
    const ModelMatrix &model, int rank, int size)
{
    ladruno_cms::SymmetricCSR result;
    result.dimension = model.order;
    result.rowOffsets.assign(1, 0);
    const long long order = model.order;
    const int begin = static_cast<int>((order * rank) / size);
    const int end = static_cast<int>((order * (rank + 1)) / size);
    std::vector<int> columns;
    std::vector<double> values;
    for (int row = 0; row < model.order; ++row) {
        if (row >= begin && row < end) {
            model.upperRow(row, columns, values);
            result.columnIndices.insert(
                result.columnIndices.end(), columns.begin(), columns.end());
            result.upperValues.insert(
                result.upperValues.end(), values.begin(), values.end());
        }
        result.rowOffsets.push_back(
            static_cast<int>(result.columnIndices.size()));
    }
    return result;
}

std::vector<double> multiplySymmetric(
    const ModelMatrix &model,
    const std::vector<double> &vectors,
    int numberOfColumns)
{
    std::vector<double> result(vectors.size(), 0.0);
    std::vector<int> columns;
    std::vector<double> values;
    for (int row = 0; row < model.order; ++row) {
        model.upperRow(row, columns, values);
        for (std::size_t entry = 0; entry < columns.size(); ++entry) {
            const int column = columns[entry];
            const double value = values[entry];
            for (int rhs = 0; rhs < numberOfColumns; ++rhs) {
                const std::size_t base =
                    static_cast<std::size_t>(rhs) * model.order;
                result[base + static_cast<std::size_t>(row)] +=
                    value * vectors[base + static_cast<std::size_t>(column)];
                if (column != row)
                    result[base + static_cast<std::size_t>(column)] +=
                        value * vectors[base + static_cast<std::size_t>(row)];
            }
        }
    }
    return result;
}

void testDistributedMumpsAtScale(int rank, int size)
{
    if (size != 2 && size != 4)
        return;

    ModelMatrix laplace;
    laplace.side = environmentSize("LADRUNO_CMS_CHECK_SIDE", 16);
    laplace.order = laplace.side * laplace.side * laplace.side;
    ModelMatrix denseBlock;
    denseBlock.dense = true;
    denseBlock.order = environmentSize("LADRUNO_CMS_CHECK_DENSE_ORDER", 1000);

    const ModelMatrix models[2] = {laplace, denseBlock};
    const char *labels[2] = {"laplace", "dense"};
    const int numberOfColumns = 2;

    for (int index = 0; index < 2; ++index) {
        const ModelMatrix &model = models[index];
        const std::string label(labels[index]);
        if (rank == 0)
            std::cout << "distributed MUMPS at scale: shape=" << label
                      << " order=" << model.order << " ranks=" << size << '\n';

        std::vector<double> expected(
            static_cast<std::size_t>(model.order) * numberOfColumns, 0.0);
        for (int row = 0; row < model.order; ++row) {
            expected[static_cast<std::size_t>(row)] =
                1.0 + 0.25 * static_cast<double>(row % 7);
            expected[static_cast<std::size_t>(model.order + row)] =
                (row % 2 == 0) ? 1.0 : -1.0;
        }
        std::vector<double> solution =
            multiplySymmetric(model, expected, numberOfColumns);

        const ladruno_cms::SymmetricCSR local =
            buildLocalBlock(model, rank, size);
        ladruno_cms::DistributedMumpsSPD solver;
        std::string message;
        const int factorized = solver.factorize(local, message);
        require(factorized == 0,
                "distributed SPD factorization failed at scale (" + label +
                    ", np=" + std::to_string(size) + "): " + message);
        if (factorized != 0)
            continue;
        REQUIRE_CALL(solver.solve(solution, numberOfColumns, message) == 0,
                "distributed multi-RHS solve failed at scale (" + label +
                    ", np=" + std::to_string(size) + "): " + message);

        double worst = 0.0;
        for (std::size_t entry = 0; entry < solution.size(); ++entry)
            worst = std::max(
                worst, std::fabs(solution[entry] - expected[entry]));
        require(worst <= 1.0e-9,
                "distributed solve at scale (" + label + ", np=" +
                    std::to_string(size) + ") recovered the wrong vector, "
                    "worst componentwise error " + std::to_string(worst));
    }
}

// The rank-local Craig-Bampton pencils go through MumpsSPD (MPI_COMM_SELF), not
// the distributed class, so the same ordering trap has to be probed there too:
// same dense shape, one rank, no collectives.
void testSerialMumpsAtScale(int rank)
{
    if (rank != 0)
        return;
    ModelMatrix model;
    model.dense = true;
    model.order = environmentSize("LADRUNO_CMS_CHECK_DENSE_ORDER", 1000);
    std::cout << "serial MUMPS at scale: shape=dense order=" << model.order
              << '\n';

    ladruno_cms::SymmetricCSR matrix;
    matrix.dimension = model.order;
    matrix.rowOffsets.assign(1, 0);
    std::vector<int> columns;
    std::vector<double> values;
    for (int row = 0; row < model.order; ++row) {
        model.upperRow(row, columns, values);
        matrix.columnIndices.insert(
            matrix.columnIndices.end(), columns.begin(), columns.end());
        matrix.upperValues.insert(
            matrix.upperValues.end(), values.begin(), values.end());
        matrix.rowOffsets.push_back(
            static_cast<int>(matrix.columnIndices.size()));
    }

    std::vector<double> expected(static_cast<std::size_t>(model.order), 0.0);
    for (int row = 0; row < model.order; ++row)
        expected[static_cast<std::size_t>(row)] =
            1.0 + 0.25 * static_cast<double>(row % 7);
    std::vector<double> solution = multiplySymmetric(model, expected, 1);

    ladruno_cms::MumpsSPD solver;
    std::string message;
    const int factorized = solver.factorize(matrix, message);
    require(factorized == 0,
            "serial SPD factorization failed at scale (dense, order " +
                std::to_string(model.order) + "): " + message);
    if (factorized != 0)
        return;
    REQUIRE_CALL(solver.solve(solution, 1, message) == 0,
            "serial solve failed at scale: " + message);
    double worst = 0.0;
    for (std::size_t entry = 0; entry < solution.size(); ++entry)
        worst = std::max(worst, std::fabs(solution[entry] - expected[entry]));
    require(worst <= 1.0e-9,
            "serial solve at scale recovered the wrong vector, worst "
            "componentwise error " + std::to_string(worst));
}

void testStaticCondensation()
{
    const std::vector<double> stiffness = {
        6.0, 1.0, 1.0, 0.5,
        1.0, 5.0, 2.0, -1.0,
        1.0, 2.0, 4.0, 1.0,
        0.5, -1.0, 1.0, 3.0};
    const std::vector<double> mass = {
        2.0, 0.0, 0.0, 0.0,
        0.0, 3.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0};
    ladruno_cms::StaticCondensation condensation;
    std::string message;
    require(
        ladruno_cms::condenseCoordinateMass(
            csrFromDense(stiffness, 4),
            csrFromDense(mass, 4),
            1.0e-12,
            1.0e-14,
            condensation,
            message) == 0,
        "coordinate-mass condensation failed");
    require(
        condensation.dynamic == std::vector<int>({0, 1}) &&
            condensation.massless == std::vector<int>({2, 3}),
        "condensation D/Z sets are incorrect");
    require(
        condensation.persistentCrossNonzeros() == 4u,
        "K_ZD was not retained in sparse form");

    const std::vector<double> reducedVectors = {1.0, 0.0, 0.0, 1.0};
    std::vector<double> fullVectors;
    require(
        ladruno_cms::reconstructStaticCoordinates(
            condensation, reducedVectors, 2, fullVectors, message) == 0,
        "static reconstruction failed");
    const std::vector<double> stiffnessTimesG =
        multiplyDense(stiffness, 4, fullVectors, 2);
    for (int column = 0; column < 2; ++column) {
        require(
            close(stiffnessTimesG[static_cast<std::size_t>(column) * 4 + 2], 0.0) &&
                close(stiffnessTimesG[static_cast<std::size_t>(column) * 4 + 3], 0.0),
            "reconstructed vector violates massless-coordinate equilibrium");
    }

    const std::vector<double> reducedStiffness = condensation.reducedStiffness.toDense();
    for (int column = 0; column < 2; ++column) {
        for (int row = 0; row < 2; ++row) {
            double fromFullCongruence = 0.0;
            for (int physical = 0; physical < 4; ++physical)
                fromFullCongruence +=
                    fullVectors[static_cast<std::size_t>(row) * 4 + physical] *
                    stiffnessTimesG[static_cast<std::size_t>(column) * 4 + physical];
            require(
                close(reducedStiffness[static_cast<std::size_t>(row) * 2 + column],
                      fromFullCongruence),
                "condensed stiffness does not equal G^T K G");
        }
    }
    const std::vector<double> reducedMass = condensation.reducedMass.toDense();
    require(
        reducedMass.size() == 4u && close(reducedMass[0], 2.0) &&
            close(reducedMass[1], 0.0) && close(reducedMass[2], 0.0) &&
            close(reducedMass[3], 3.0),
        "condensed mass is not M_DD");

    std::vector<double> singularStiffness = stiffness;
    singularStiffness[2u * 4u + 2u] = 0.0;
    singularStiffness[2u * 4u + 3u] = 0.0;
    singularStiffness[3u * 4u + 2u] = 0.0;
    singularStiffness[3u * 4u + 3u] = 0.0;
    require(
        ladruno_cms::condenseCoordinateMass(
            csrFromDense(singularStiffness, 4),
            csrFromDense(mass, 4),
            1.0e-12,
            1.0e-14,
            condensation,
            message) < 0,
        "singular K_ZZ passed the condensation guard");

    std::vector<double> indefiniteMass = mass;
    indefiniteMass[0] = 1.0;
    indefiniteMass[1] = 2.0;
    indefiniteMass[4] = 2.0;
    indefiniteMass[5] = 1.0;
    require(
        ladruno_cms::condenseCoordinateMass(
            csrFromDense(stiffness, 4),
            csrFromDense(indefiniteMass, 4),
            1.0e-12,
            1.0e-14,
            condensation,
            message) < 0,
        "indefinite M_DD passed the condensation guard");
}

void testBlockedSparseCondensation()
{
    constexpr int dynamic = 40;
    constexpr int massless = 2;
    constexpr int order = dynamic + massless;
    std::vector<double> stiffness(static_cast<std::size_t>(order) * order, 0.0);
    std::vector<double> mass(static_cast<std::size_t>(order) * order, 0.0);
    for (int equation = 0; equation < dynamic; ++equation) {
        stiffness[static_cast<std::size_t>(equation) * order + equation] = 10.0;
        mass[static_cast<std::size_t>(equation) * order + equation] =
            1.0 + 0.01 * equation;
        const int z = dynamic + equation % massless;
        const double coupling = 0.05 + 0.001 * equation;
        stiffness[static_cast<std::size_t>(equation) * order + z] = coupling;
        stiffness[static_cast<std::size_t>(z) * order + equation] = coupling;
    }
    stiffness[static_cast<std::size_t>(dynamic) * order + dynamic] = 4.0;
    stiffness[static_cast<std::size_t>(dynamic + 1) * order + dynamic + 1] = 5.0;
    stiffness[static_cast<std::size_t>(dynamic) * order + dynamic + 1] = 0.2;
    stiffness[static_cast<std::size_t>(dynamic + 1) * order + dynamic] = 0.2;

    ladruno_cms::StaticCondensation condensation;
    std::string message;
    require(
        ladruno_cms::condenseCoordinateMass(
            csrFromDense(stiffness, order),
            csrFromDense(mass, order),
            1.0e-12,
            1.0e-14,
            condensation,
            message) == 0,
        "blocked sparse condensation failed");
    require(
        condensation.persistentCrossNonzeros() == static_cast<std::size_t>(dynamic),
        "blocked condensation densified sparse K_ZD");

    std::vector<double> identity(static_cast<std::size_t>(dynamic) * dynamic, 0.0);
    for (int column = 0; column < dynamic; ++column)
        identity[static_cast<std::size_t>(column) * dynamic + column] = 1.0;
    std::vector<double> transformation;
    require(
        ladruno_cms::reconstructStaticCoordinates(
            condensation, identity, dynamic, transformation, message) == 0,
        "blocked static reconstruction failed");
    const std::vector<double> stiffnessTimesTransformation =
        multiplyDense(stiffness, order, transformation, dynamic);
    for (int column = 0; column < dynamic; ++column) {
        for (int z = dynamic; z < order; ++z)
            require(
                close(stiffnessTimesTransformation[
                    static_cast<std::size_t>(column) * order + z], 0.0),
                "blocked reconstruction violates a Z equilibrium row");
    }
    const std::vector<double> reducedStiffness = condensation.reducedStiffness.toDense();
    for (int column = 0; column < dynamic; ++column) {
        for (int row = 0; row < dynamic; ++row) {
            double congruence = 0.0;
            for (int physical = 0; physical < order; ++physical)
                congruence +=
                    transformation[static_cast<std::size_t>(row) * order + physical] *
                    stiffnessTimesTransformation[
                        static_cast<std::size_t>(column) * order + physical];
            require(
                close(reducedStiffness[static_cast<std::size_t>(row) * dynamic + column],
                      congruence),
                "blocked Khat does not equal G^T K G");
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    {
        testMumpsMultiRHS();
        testDistributedMumps(rank, size);
        testDistributedMumpsAtScale(rank, size);
        testSerialMumpsAtScale(rank);
        testStaticCondensation();
        testBlockedSparseCondensation();
    }
    MPI_Finalize();
    if (failures != 0)
        return 1;
    std::cout << "Ladruno CMS MUMPS and condensation checks passed\n";
    return 0;
}
