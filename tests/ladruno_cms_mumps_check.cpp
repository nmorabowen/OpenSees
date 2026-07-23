#include "LadrunoCMSMumps.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
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

void testDistributedMumps(int rank, int size)
{
    if (size != 4)
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
        testStaticCondensation();
        testBlockedSparseCondensation();
    }
    MPI_Finalize();
    if (failures != 0)
        return 1;
    std::cout << "Ladruno CMS MUMPS and condensation checks passed\n";
    return 0;
}
