#include "LadrunoCMSSubspaceRefiner.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
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

bool close(double left, double right, double tolerance = 1.0e-9)
{
    return std::fabs(left - right) <= tolerance *
        std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
}

ladruno_cms::SymmetricCSR rootOwnedDiagonal(
    const std::vector<double> &diagonal,
    int rank)
{
    ladruno_cms::SymmetricCSR result;
    result.dimension = static_cast<int>(diagonal.size());
    result.rowOffsets.assign(diagonal.size() + 1u, 0);
    if (rank != 0)
        return result;
    for (std::size_t row = 0; row < diagonal.size(); ++row) {
        if (diagonal[row] != 0.0) {
            result.columnIndices.push_back(static_cast<int>(row));
            result.upperValues.push_back(diagonal[row]);
        }
        result.rowOffsets[row + 1u] =
            static_cast<int>(result.columnIndices.size());
    }
    return result;
}

void checkOriginalPencilConvergence(int rank)
{
    constexpr int order = 12;
    constexpr int requested = 3;
    constexpr int block = 6;
    std::vector<double> stiffnessDiagonal(order);
    std::vector<double> massDiagonal(order, 1.0);
    for (int row = 0; row < order; ++row)
        stiffnessDiagonal[static_cast<std::size_t>(row)] = row + 1.0;
    std::vector<double> initial(static_cast<std::size_t>(order) * block, 0.0);
    for (int column = 0; column < block; ++column) {
        initial[static_cast<std::size_t>(column) * order + column] = 1.0;
        initial[static_cast<std::size_t>(column) * order + column + block] = 0.2;
    }
    ladruno_cms::SubspaceRefinementControls controls;
    controls.globalDimension = order;
    controls.numberOfModes = requested;
    controls.numberOfIterationVectors = block;
    controls.maximumIterations = 24;
    controls.tolerance = 1.0e-10;
    controls.massRankTolerance = 1.0e-12;
    ladruno_cms::SubspaceRefinementResult result;
    std::string message;
    const int code = ladruno_cms::refineDistributedSubspace(
        rootOwnedDiagonal(stiffnessDiagonal, rank),
        rootOwnedDiagonal(massDiagonal, rank),
        initial, controls, result, message);
    require(code == 0, "original-pencil refinement failed: " + message);
    require(result.converged, "original-pencil refinement did not converge");
    require(result.history.size() > 1u,
            "test fixture did not exercise an inverse subspace iteration");
    if (result.history.size() > 1u)
        require(result.history.back().maximumResidual <
                    result.history.front().maximumResidual,
                "original-pencil residual did not decrease");
    for (int mode = 0; mode < requested &&
         mode < static_cast<int>(result.eigenvalues.size()); ++mode) {
        require(close(result.eigenvalues[static_cast<std::size_t>(mode)], mode + 1.0),
                "refined eigenvalue is inaccurate");
        require(result.normalizedResiduals[static_cast<std::size_t>(mode)] <=
                    controls.tolerance,
                "accepted mode exceeds the original-pencil tolerance");
    }
}

void checkCoordinateSemidefiniteMass(int rank)
{
    const std::vector<double> stiffness = {2.0, 3.0, 4.0, 5.0};
    const std::vector<double> mass = {1.0, 1.0, 0.0, 0.0};
    const std::vector<double> initial = {
        1.0, 0.0, 0.4, 0.0,
        0.0, 1.0, 0.0, -0.3};
    ladruno_cms::SubspaceRefinementControls controls;
    controls.globalDimension = 4;
    controls.numberOfModes = 2;
    controls.numberOfIterationVectors = 2;
    controls.maximumIterations = 4;
    controls.tolerance = 1.0e-12;
    controls.massRankTolerance = 1.0e-12;
    ladruno_cms::SubspaceRefinementResult result;
    std::string message;
    const int code = ladruno_cms::refineDistributedSubspace(
        rootOwnedDiagonal(stiffness, rank), rootOwnedDiagonal(mass, rank),
        initial, controls, result, message);
    require(code == 0, "coordinate-semidefinite mass was rejected: " + message);
    require(result.eigenvalues.size() == 2u && close(result.eigenvalues[0], 2.0) &&
                close(result.eigenvalues[1], 3.0),
            "coordinate-semidefinite eigenvalues are incorrect");
}

void checkMassRankLossIsRejected(int rank)
{
    const std::vector<double> initial = {0.0, 1.0};
    ladruno_cms::SubspaceRefinementControls controls;
    controls.globalDimension = 2;
    controls.numberOfModes = 1;
    controls.numberOfIterationVectors = 1;
    controls.maximumIterations = 2;
    controls.tolerance = 1.0e-10;
    ladruno_cms::SubspaceRefinementResult result;
    std::string message;
    require(ladruno_cms::refineDistributedSubspace(
                rootOwnedDiagonal({2.0, 3.0}, rank),
                rootOwnedDiagonal({1.0, 0.0}, rank),
                initial, controls, result, message) < 0,
            "an M-null starting block was accepted");
}

void checkSingularStiffnessIsRejected(int rank)
{
    const std::vector<double> initial = {1.0, 1.0};
    ladruno_cms::SubspaceRefinementControls controls;
    controls.globalDimension = 2;
    controls.numberOfModes = 1;
    controls.numberOfIterationVectors = 1;
    controls.maximumIterations = 2;
    controls.tolerance = 1.0e-12;
    ladruno_cms::SubspaceRefinementResult result;
    std::string message;
    require(ladruno_cms::refineDistributedSubspace(
                rootOwnedDiagonal({0.0, 2.0}, rank),
                rootOwnedDiagonal({1.0, 1.0}, rank),
                initial, controls, result, message) < 0,
            "a singular original stiffness was accepted by the SPD inverse action");
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    checkOriginalPencilConvergence(rank);
    checkCoordinateSemidefiniteMass(rank);
    checkMassRankLossIsRejected(rank);
    checkSingularStiffnessIsRejected(rank);
    int globalFailures = 0;
    MPI_Allreduce(&failures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Finalize();
    if (globalFailures != 0)
        return 1;
    if (rank == 0)
        std::cout << "Ladruno CMS original-pencil subspace checks passed\n";
    return 0;
}
