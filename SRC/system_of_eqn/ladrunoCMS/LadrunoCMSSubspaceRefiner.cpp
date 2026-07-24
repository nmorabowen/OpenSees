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

#include "LadrunoCMSSubspaceRefiner.h"

#include "LadrunoCMSMumps.h"

#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#ifdef _WIN32
extern "C" int DSYGVX(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_REFINER_DSYGVX DSYGVX
#else
extern "C" int dsygvx_(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_REFINER_DSYGVX dsygvx_
#endif

namespace ladruno_cms {
namespace {

bool checkedProduct(int left, int right, std::size_t &value)
{
    if (left < 0 || right < 0)
        return false;
    const std::size_t first = static_cast<std::size_t>(left);
    const std::size_t second = static_cast<std::size_t>(right);
    if (first != 0u && second > std::numeric_limits<std::size_t>::max() / first)
        return false;
    value = first * second;
    return true;
}

double dotColumn(
    const std::vector<double> &left,
    const std::vector<double> &right,
    int dimension,
    int leftColumn,
    int rightColumn)
{
    double value = 0.0;
    const std::size_t leftOffset =
        static_cast<std::size_t>(leftColumn) * dimension;
    const std::size_t rightOffset =
        static_cast<std::size_t>(rightColumn) * dimension;
    for (int row = 0; row < dimension; ++row)
        value += left[leftOffset + static_cast<std::size_t>(row)] *
            right[rightOffset + static_cast<std::size_t>(row)];
    return value;
}

int globalAction(
    const SymmetricCSR &localMatrix,
    const std::vector<double> &vectors,
    int numberOfColumns,
    std::vector<double> &actions,
    std::string &message)
{
    const int dimension = localMatrix.dimension;
    std::size_t entries = 0;
    if (!checkedProduct(dimension, numberOfColumns, entries) ||
        entries > static_cast<std::size_t>(INT_MAX) ||
        vectors.size() != entries) {
        message = "distributed subspace action dimensions exceed the MPI range";
        return -1;
    }
    std::vector<double> local(entries, 0.0);
    for (int column = 0; column < numberOfColumns; ++column) {
        const auto begin = vectors.begin() +
            static_cast<std::ptrdiff_t>(column) * dimension;
        std::vector<double> vector(begin, begin + dimension);
        const std::vector<double> product = localMatrix.multiply(vector);
        if (product.size() != static_cast<std::size_t>(dimension)) {
            message = "distributed subspace sparse action failed";
            return -2;
        }
        std::copy(product.begin(), product.end(),
                  local.begin() + static_cast<std::ptrdiff_t>(column) * dimension);
    }
    actions.assign(entries, 0.0);
    MPI_Allreduce(local.data(), actions.data(), static_cast<int>(entries),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return 0;
}

int massOrthonormalize(
    const SymmetricCSR &localMass,
    const std::vector<double> &input,
    int numberOfColumns,
    double rankTolerance,
    std::vector<double> &basis,
    std::vector<double> &massBasis,
    double &orthogonalityLoss,
    double &massActionSeconds,
    std::string &message)
{
    const int dimension = localMass.dimension;
    const double actionStart = MPI_Wtime();
    std::vector<double> inputMass;
    if (globalAction(localMass, input, numberOfColumns, inputMass, message) < 0)
        return -1;
    massActionSeconds += MPI_Wtime() - actionStart;
    basis.assign(input.size(), 0.0);
    massBasis.assign(input.size(), 0.0);
    for (int column = 0; column < numberOfColumns; ++column) {
        const std::size_t offset = static_cast<std::size_t>(column) * dimension;
        std::copy(input.begin() + static_cast<std::ptrdiff_t>(offset),
                  input.begin() + static_cast<std::ptrdiff_t>(offset + dimension),
                  basis.begin() + static_cast<std::ptrdiff_t>(offset));
        std::copy(inputMass.begin() + static_cast<std::ptrdiff_t>(offset),
                  inputMass.begin() + static_cast<std::ptrdiff_t>(offset + dimension),
                  massBasis.begin() + static_cast<std::ptrdiff_t>(offset));
        const double originalMassNorm =
            dotColumn(basis, massBasis, dimension, column, column);
        if (!std::isfinite(originalMassNorm) || originalMassNorm <= 0.0) {
            message = "subspace iteration block lost rank in the M product";
            return -2;
        }
        for (int pass = 0; pass < 2; ++pass) {
            for (int previous = 0; previous < column; ++previous) {
                const double projection =
                    dotColumn(basis, massBasis, dimension, previous, column);
                const std::size_t previousOffset =
                    static_cast<std::size_t>(previous) * dimension;
                for (int row = 0; row < dimension; ++row) {
                    const std::size_t current = offset + static_cast<std::size_t>(row);
                    const std::size_t old = previousOffset + static_cast<std::size_t>(row);
                    basis[current] -= projection * basis[old];
                    massBasis[current] -= projection * massBasis[old];
                }
            }
        }
        const double massNorm =
            dotColumn(basis, massBasis, dimension, column, column);
        const double threshold = rankTolerance *
            std::max(1.0, std::fabs(originalMassNorm));
        if (!std::isfinite(massNorm) || massNorm <= threshold) {
            message = "subspace iteration block lost rank in the M product";
            return -3;
        }
        const double inverseNorm = 1.0 / std::sqrt(massNorm);
        for (int row = 0; row < dimension; ++row) {
            const std::size_t current = offset + static_cast<std::size_t>(row);
            basis[current] *= inverseNorm;
            massBasis[current] *= inverseNorm;
        }
    }

    orthogonalityLoss = 0.0;
    for (int column = 0; column < numberOfColumns; ++column)
        for (int row = 0; row < numberOfColumns; ++row) {
            const double expected = row == column ? 1.0 : 0.0;
            orthogonalityLoss = std::max(
                orthogonalityLoss,
                std::fabs(dotColumn(
                    basis, massBasis, dimension, row, column) - expected));
        }
    return 0;
}

int solveProjected(
    std::vector<double> projectedStiffness,
    int order,
    std::vector<double> &eigenvalues,
    std::vector<double> &coordinates,
    std::string &message)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int localFailure = 0;
    if (rank == 0) {
        std::vector<double> projectedMass(
            static_cast<std::size_t>(order) * order, 0.0);
        for (int diagonal = 0; diagonal < order; ++diagonal)
            projectedMass[static_cast<std::size_t>(diagonal) * order + diagonal] = 1.0;
        int itype = 1;
        char jobz = 'V';
        char range = 'I';
        char uplo = 'U';
        int leading = order;
        double lowerValue = 0.0;
        double upperValue = 0.0;
        int first = 1;
        int last = order;
        double absoluteTolerance = 0.0;
        int found = 0;
        eigenvalues.assign(static_cast<std::size_t>(order), 0.0);
        coordinates.assign(static_cast<std::size_t>(order) * order, 0.0);
        int workspaceSize = std::max(1, 8 * order);
        std::vector<double> workspace(static_cast<std::size_t>(workspaceSize));
        std::vector<int> integerWorkspace(static_cast<std::size_t>(5 * order));
        std::vector<int> failed(static_cast<std::size_t>(order));
        int info = 0;
        LADRUNO_REFINER_DSYGVX(
            &itype, &jobz, &range, &uplo, &order,
            projectedStiffness.data(), &leading,
            projectedMass.data(), &leading,
            &lowerValue, &upperValue, &first, &last, &absoluteTolerance,
            &found, eigenvalues.data(), coordinates.data(), &leading,
            workspace.data(), &workspaceSize, integerWorkspace.data(),
            failed.data(), &info);
        if (info != 0 || found != order) {
            message = "subspace Rayleigh-Ritz dsygvx failed with INFO=" +
                std::to_string(info) + " found=" + std::to_string(found);
            localFailure = 1;
        }
    }
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (localFailure) {
        if (rank != 0)
            message = "subspace Rayleigh-Ritz failed on rank 0";
        return -1;
    }
    if (rank != 0) {
        eigenvalues.resize(static_cast<std::size_t>(order));
        coordinates.resize(static_cast<std::size_t>(order) * order);
    }
    MPI_Bcast(eigenvalues.data(), order, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coordinates.data(), order * order, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return 0;
}

int rayleighRitz(
    const SymmetricCSR &localStiffness,
    const SymmetricCSR &localMass,
    const std::vector<double> &input,
    int numberOfColumns,
    double rankTolerance,
    std::vector<double> &eigenvalues,
    std::vector<double> &vectors,
    std::vector<double> &stiffnessActions,
    std::vector<double> &massActions,
    double &orthogonalityLoss,
    double &massActionSeconds,
    std::string &message)
{
    const int dimension = localStiffness.dimension;
    std::vector<double> basis;
    std::vector<double> basisMass;
    if (massOrthonormalize(
            localMass, input, numberOfColumns, rankTolerance,
            basis, basisMass, orthogonalityLoss,
            massActionSeconds, message) < 0)
        return -1;
    std::vector<double> basisStiffness;
    if (globalAction(
            localStiffness, basis, numberOfColumns,
            basisStiffness, message) < 0)
        return -2;
    std::vector<double> projected(
        static_cast<std::size_t>(numberOfColumns) * numberOfColumns, 0.0);
    for (int column = 0; column < numberOfColumns; ++column)
        for (int row = 0; row <= column; ++row) {
            const double forward = dotColumn(
                basis, basisStiffness, dimension, row, column);
            const double reverse = dotColumn(
                basis, basisStiffness, dimension, column, row);
            const double value = 0.5 * (forward + reverse);
            projected[static_cast<std::size_t>(column) * numberOfColumns + row] = value;
            projected[static_cast<std::size_t>(row) * numberOfColumns + column] = value;
        }
    std::vector<double> coordinates;
    if (solveProjected(
            projected, numberOfColumns,
            eigenvalues, coordinates, message) < 0)
        return -3;

    std::size_t entries = 0;
    if (!checkedProduct(dimension, numberOfColumns, entries)) {
        message = "subspace Ritz vectors exceed the supported size";
        return -4;
    }
    vectors.assign(entries, 0.0);
    stiffnessActions.assign(entries, 0.0);
    massActions.assign(entries, 0.0);
    for (int mode = 0; mode < numberOfColumns; ++mode)
        for (int source = 0; source < numberOfColumns; ++source) {
            const double coefficient = coordinates[
                static_cast<std::size_t>(mode) * numberOfColumns + source];
            for (int row = 0; row < dimension; ++row) {
                const std::size_t target =
                    static_cast<std::size_t>(mode) * dimension + row;
                const std::size_t origin =
                    static_cast<std::size_t>(source) * dimension + row;
                vectors[target] += coefficient * basis[origin];
                stiffnessActions[target] += coefficient * basisStiffness[origin];
                massActions[target] += coefficient * basisMass[origin];
            }
        }
    return 0;
}

// MUMPS disables its built-in iterative refinement when NRHS>1.  Subspace
// iteration necessarily solves a block, so perform the recommended two
// correction steps explicitly while reusing the same distributed factor.
int solveInverseBlock(
    DistributedMumpsSPD &factor,
    const SymmetricCSR &localStiffness,
    const std::vector<double> &rightHandSides,
    int numberOfColumns,
    std::vector<double> &solution,
    double &solveSeconds,
    double &refinementSeconds,
    std::string &message)
{
    solution = rightHandSides;
    double start = MPI_Wtime();
    if (factor.solve(solution, numberOfColumns, message) < 0)
        return -1;
    solveSeconds += MPI_Wtime() - start;
    const double refinementStart = MPI_Wtime();
    for (int step = 0; step < 2; ++step) {
        std::vector<double> stiffnessSolution;
        if (globalAction(
                localStiffness, solution, numberOfColumns,
                stiffnessSolution, message) < 0)
            return -2;
        std::vector<double> correction(rightHandSides.size(), 0.0);
        for (std::size_t entry = 0; entry < correction.size(); ++entry)
            correction[entry] = rightHandSides[entry] - stiffnessSolution[entry];
        start = MPI_Wtime();
        if (factor.solve(correction, numberOfColumns, message) < 0)
            return -3;
        solveSeconds += MPI_Wtime() - start;
        for (std::size_t entry = 0; entry < solution.size(); ++entry)
            solution[entry] += correction[entry];
    }
    refinementSeconds += MPI_Wtime() - refinementStart;
    return 0;
}

void evaluateResiduals(
    const std::vector<double> &eigenvalues,
    const std::vector<double> &stiffnessActions,
    const std::vector<double> &massActions,
    int dimension,
    int numberOfModes,
    SubspaceRefinementResult &result)
{
    result.normalizedResiduals.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    result.residualNorms.assign(static_cast<std::size_t>(numberOfModes), 0.0);
    result.stiffnessActionNorms.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    result.massActionNorms.assign(static_cast<std::size_t>(numberOfModes), 0.0);
    result.maximumResidualEquation.assign(
        static_cast<std::size_t>(numberOfModes), -1);
    for (int mode = 0; mode < numberOfModes; ++mode) {
        const double eigenvalue = eigenvalues[static_cast<std::size_t>(mode)];
        double residualSquare = 0.0;
        double stiffnessSquare = 0.0;
        double massSquare = 0.0;
        double maximum = -1.0;
        int maximumEquation = -1;
        const std::size_t offset = static_cast<std::size_t>(mode) * dimension;
        for (int row = 0; row < dimension; ++row) {
            const std::size_t index = offset + static_cast<std::size_t>(row);
            const double stiffness = stiffnessActions[index];
            const double mass = massActions[index];
            const double residual = stiffness - eigenvalue * mass;
            residualSquare += residual * residual;
            stiffnessSquare += stiffness * stiffness;
            massSquare += mass * mass;
            if (std::fabs(residual) > maximum) {
                maximum = std::fabs(residual);
                maximumEquation = row;
            }
        }
        const double residualNorm = std::sqrt(std::max(0.0, residualSquare));
        const double stiffnessNorm = std::sqrt(std::max(0.0, stiffnessSquare));
        const double massNorm = std::sqrt(std::max(0.0, massSquare));
        result.residualNorms[static_cast<std::size_t>(mode)] = residualNorm;
        result.stiffnessActionNorms[static_cast<std::size_t>(mode)] = stiffnessNorm;
        result.massActionNorms[static_cast<std::size_t>(mode)] = massNorm;
        result.maximumResidualEquation[static_cast<std::size_t>(mode)] =
            maximumEquation;
        // Combined backward error ||Kx - lambda*Mx|| / (||Kx|| + |lambda|*||Mx||)
        // — NOT ||Kx||-only. Matches the ADR (section 3.6) and the numpy oracle;
        // at convergence (||Kx|| ~ |lambda|*||Mx||) it is ~2x tighter than the
        // ||Kx||-only form would report.
        result.normalizedResiduals[static_cast<std::size_t>(mode)] =
            residualNorm / std::max(
                std::numeric_limits<double>::min(),
                stiffnessNorm + std::fabs(eigenvalue) * massNorm);
    }
}

void retainRequestedModes(
    const std::vector<double> &allEigenvalues,
    const std::vector<double> &allVectors,
    int dimension,
    int numberOfModes,
    SubspaceRefinementResult &result)
{
    result.eigenvalues.assign(
        allEigenvalues.begin(), allEigenvalues.begin() + numberOfModes);
    result.eigenvectors.assign(
        allVectors.begin(),
        allVectors.begin() + static_cast<std::ptrdiff_t>(dimension) * numberOfModes);
}

} // namespace

int refineDistributedSubspace(
    const SymmetricCSR &localStiffness,
    const SymmetricCSR &localMass,
    const std::vector<double> &initialVectors,
    const SubspaceRefinementControls &controls,
    SubspaceRefinementResult &result,
    std::string &message)
{
    result = SubspaceRefinementResult{};
    message.clear();
    std::size_t expected = 0;
    int localFailure =
        controls.globalDimension < 1 ||
        controls.numberOfModes < 1 ||
        controls.numberOfIterationVectors < controls.numberOfModes ||
        controls.maximumIterations < 0 ||
        !(controls.tolerance > 0.0) ||
        localStiffness.dimension != controls.globalDimension ||
        localMass.dimension != controls.globalDimension ||
        localStiffness.validate(message) < 0 ||
        localMass.validate(message) < 0 ||
        !checkedProduct(
            controls.globalDimension,
            controls.numberOfIterationVectors,
            expected) ||
        initialVectors.size() != expected;
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (globalFailure) {
        if (!localFailure)
            message = "subspace-refinement inputs are invalid on another rank";
        else if (message.empty())
            message = "subspace-refinement inputs are invalid";
        return -1;
    }
    const double rankTolerance = controls.massRankTolerance > 0.0
        ? controls.massRankTolerance
        : 1024.0 * std::numeric_limits<double>::epsilon();
    result.numberOfIterationVectors = controls.numberOfIterationVectors;

    std::vector<double> allEigenvalues;
    std::vector<double> allVectors;
    std::vector<double> stiffnessActions;
    std::vector<double> massActions;
    double orthogonalityLoss = 0.0;
    double ritzStart = MPI_Wtime();
    if (rayleighRitz(
            localStiffness, localMass, initialVectors,
            controls.numberOfIterationVectors, rankTolerance,
            allEigenvalues, allVectors, stiffnessActions, massActions,
            orthogonalityLoss, result.massActionSeconds, message) < 0)
        return -2;
    result.rayleighRitzSeconds += MPI_Wtime() - ritzStart;
    evaluateResiduals(
        allEigenvalues, stiffnessActions, massActions,
        controls.globalDimension, controls.numberOfModes, result);
    double maximumResidual = *std::max_element(
        result.normalizedResiduals.begin(), result.normalizedResiduals.end());
    result.maximumMassOrthogonalityLoss = orthogonalityLoss;
    result.history.push_back({0, maximumResidual, orthogonalityLoss});
    result.converged = maximumResidual <= controls.tolerance;
    if (result.converged) {
        retainRequestedModes(
            allEigenvalues, allVectors, controls.globalDimension,
            controls.numberOfModes, result);
        return 0;
    }

    DistributedMumpsSPD factor;
    const double factorStart = MPI_Wtime();
    if (factor.factorize(localStiffness, message) < 0) {
        message = "original stiffness factorization failed: " + message;
        return -3;
    }
    result.factorizationSeconds = MPI_Wtime() - factorStart;

    for (int iteration = 1; iteration <= controls.maximumIterations; ++iteration) {
        std::vector<double> inverseBlock;
        if (solveInverseBlock(
                factor, localStiffness, massActions,
                controls.numberOfIterationVectors, inverseBlock,
                result.solveSeconds, result.inverseRefinementSeconds,
                message) < 0) {
            message = "subspace inverse action failed: " + message;
            return -4;
        }

        ritzStart = MPI_Wtime();
        if (rayleighRitz(
                localStiffness, localMass, inverseBlock,
                controls.numberOfIterationVectors, rankTolerance,
                allEigenvalues, allVectors, stiffnessActions, massActions,
                orthogonalityLoss, result.massActionSeconds, message) < 0)
            return -5;
        result.rayleighRitzSeconds += MPI_Wtime() - ritzStart;
        evaluateResiduals(
            allEigenvalues, stiffnessActions, massActions,
            controls.globalDimension, controls.numberOfModes, result);
        maximumResidual = *std::max_element(
            result.normalizedResiduals.begin(), result.normalizedResiduals.end());
        result.maximumMassOrthogonalityLoss = std::max(
            result.maximumMassOrthogonalityLoss, orthogonalityLoss);
        result.history.push_back({iteration, maximumResidual, orthogonalityLoss});
        result.iterations = iteration;
        if (maximumResidual <= controls.tolerance) {
            result.converged = true;
            retainRequestedModes(
                allEigenvalues, allVectors, controls.globalDimension,
                controls.numberOfModes, result);
            return 0;
        }
    }

    retainRequestedModes(
        allEigenvalues, allVectors, controls.globalDimension,
        controls.numberOfModes, result);
    std::ostringstream diagnostic;
    diagnostic << std::scientific << std::setprecision(8)
               << "subspace refinement did not converge: maximum residual="
               << maximumResidual << " tolerance=" << controls.tolerance
               << " iterations=" << controls.maximumIterations;
    message = diagnostic.str();
    return -6;
}

} // namespace ladruno_cms
