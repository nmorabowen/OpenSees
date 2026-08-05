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

#include "LadrunoCMSLocalLanczos.h"

#include <mpi.h>

#include "LadrunoCMSMumps.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

#ifdef _WIN32
extern "C" void DSYEV(
    char *, char *, int *, double *, int *, double *, double *, int *, int *);
#define LADRUNO_DSYEV DSYEV
#else
extern "C" void dsyev_(
    char *, char *, int *, double *, int *, double *, double *, int *, int *);
#define LADRUNO_DSYEV dsyev_
#endif

namespace ladruno_cms {
namespace {

double dot(const double *left, const double *right, int size)
{
    double value = 0.0;
    for (int index = 0; index < size; ++index)
        value += left[index] * right[index];
    return value;
}

double norm2(const std::vector<double> &values)
{
    return std::sqrt(std::max(0.0, dot(values.data(), values.data(),
                                      static_cast<int>(values.size()))));
}

std::vector<double> column(
    const std::vector<double> &matrix,
    int rows,
    int columnIndex)
{
    const std::size_t begin = static_cast<std::size_t>(columnIndex) * rows;
    return std::vector<double>(matrix.begin() + static_cast<std::ptrdiff_t>(begin),
                               matrix.begin() + static_cast<std::ptrdiff_t>(begin + rows));
}

void appendColumn(std::vector<double> &matrix, const std::vector<double> &values)
{
    matrix.insert(matrix.end(), values.begin(), values.end());
}

bool finiteVector(const std::vector<double> &values)
{
    return std::all_of(
        values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

int mOrthonormalize(
    std::vector<double> &candidate,
    const std::vector<double> &basis,
    int dimension,
    const SymmetricCSR &mass,
    double &orthogonalityLoss)
{
    const int basisColumns = static_cast<int>(basis.size() / static_cast<std::size_t>(dimension));
    std::vector<double> massCandidate = mass.multiply(candidate);
    if (massCandidate.size() != candidate.size())
        return -1;
    const double initialSquaredNorm = dot(candidate.data(), massCandidate.data(), dimension);
    if (!std::isfinite(initialSquaredNorm) || initialSquaredNorm <= 0.0)
        return 0;
    const double initialNorm = std::sqrt(initialSquaredNorm);

    for (int pass = 0; pass < 2; ++pass) {
        massCandidate = mass.multiply(candidate);
        if (massCandidate.size() != candidate.size())
            return -1;
        for (int basisColumn = 0; basisColumn < basisColumns; ++basisColumn) {
            const double *basisVector =
                basis.data() + static_cast<std::size_t>(basisColumn) * dimension;
            const double projection = dot(basisVector, massCandidate.data(), dimension);
            for (int row = 0; row < dimension; ++row)
                candidate[static_cast<std::size_t>(row)] -= projection * basisVector[row];
        }
    }

    massCandidate = mass.multiply(candidate);
    if (massCandidate.size() != candidate.size())
        return -1;
    const double squaredNorm = dot(candidate.data(), massCandidate.data(), dimension);
    const double breakdownThreshold = 100.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, initialNorm);
    if (!std::isfinite(squaredNorm) || squaredNorm <= breakdownThreshold * breakdownThreshold)
        return 0;
    const double inverseNorm = 1.0 / std::sqrt(squaredNorm);
    for (double &value : candidate)
        value *= inverseNorm;

    massCandidate = mass.multiply(candidate);
    double localLoss = std::fabs(dot(candidate.data(), massCandidate.data(), dimension) - 1.0);
    for (int basisColumn = 0; basisColumn < basisColumns; ++basisColumn) {
        const double *basisVector =
            basis.data() + static_cast<std::size_t>(basisColumn) * dimension;
        localLoss = std::max(
            localLoss, std::fabs(dot(basisVector, massCandidate.data(), dimension)));
    }
    orthogonalityLoss = std::max(orthogonalityLoss, localLoss);
    return 1;
}

std::vector<double> deterministicVector(
    int dimension,
    int level,
    int subdomain,
    int restart,
    int attempt)
{
    std::vector<double> values(static_cast<std::size_t>(dimension), 0.0);
    const double seed = 1.0 + 17.0 * level + 31.0 * subdomain +
        43.0 * restart + 59.0 * attempt;
    for (int row = 0; row < dimension; ++row) {
        const double coordinate = static_cast<double>(row + 1);
        values[static_cast<std::size_t>(row)] =
            std::sin(0.6180339887498949 * seed * coordinate) +
            0.5 * std::cos(0.4142135623730950 * (seed + 3.0) * coordinate);
    }
    return values;
}

int applyOperator(
    const SymmetricCSR &mass,
    MumpsSPD &stiffnessFactor,
    const std::vector<double> &input,
    std::vector<double> &output,
    std::string &message)
{
    output = mass.multiply(input);
    if (output.size() != input.size()) {
        message = "mass multiplication failed in the local Lanczos operator";
        return -1;
    }
    if (stiffnessFactor.solve(output, 1, message) < 0) {
        message = "local Lanczos inverse action failed: " + message;
        return -2;
    }
    return 0;
}

int symmetricEigen(
    std::vector<double> &matrix,
    int order,
    std::vector<double> &eigenvalues,
    std::string &message)
{
    eigenvalues.assign(static_cast<std::size_t>(order), 0.0);
    char jobz = 'V';
    char uplo = 'U';
    int leadingDimension = order;
    int info = 0;
    int workspaceSize = -1;
    double workspaceQuery = 0.0;
    LADRUNO_DSYEV(
        &jobz, &uplo, &order, matrix.data(), &leadingDimension,
        eigenvalues.data(), &workspaceQuery, &workspaceSize, &info);
    if (info != 0 || !std::isfinite(workspaceQuery)) {
        message = "LAPACK workspace query failed in Rayleigh-Ritz";
        return -1;
    }
    if (workspaceQuery > static_cast<double>(std::numeric_limits<int>::max())) {
        message = "Rayleigh-Ritz workspace exceeds LAPACK integer range";
        return -2;
    }
    workspaceSize = std::max(1, static_cast<int>(std::ceil(workspaceQuery)));
    std::vector<double> workspace(static_cast<std::size_t>(workspaceSize));
    LADRUNO_DSYEV(
        &jobz, &uplo, &order, matrix.data(), &leadingDimension,
        eigenvalues.data(), workspace.data(), &workspaceSize, &info);
    if (info != 0 || !finiteVector(matrix) || !finiteVector(eigenvalues)) {
        message = "LAPACK eigensolve failed in Rayleigh-Ritz with INFO=" +
            std::to_string(info);
        return -3;
    }
    return 0;
}

// ADR-1000 P4 section 4: `massOperatorBasis` is M*operatorBasis, maintained
// incrementally by the caller. This routine used to recompute it here -- one
// sparse mat-vec plus two vector allocations per basis column, on EVERY call,
// for a basis that had only grown by one block since the previous call. The
// values are identical; they are just computed once per column instead of once
// per column per iteration.
int rayleighRitz(
    const std::vector<double> &basis,
    const std::vector<double> &massOperatorBasis,
    int dimension,
    std::vector<double> &ritzValues,
    std::vector<double> &ritzCoordinates,
    std::string &message)
{
    const int columns = static_cast<int>(basis.size() / static_cast<std::size_t>(dimension));
    if (massOperatorBasis.size() !=
        static_cast<std::size_t>(dimension) * static_cast<std::size_t>(columns)) {
        message = "cached M*operator basis is out of step with the Krylov basis";
        return -1;
    }
    ritzCoordinates.assign(static_cast<std::size_t>(columns) * columns, 0.0);
    for (int columnIndex = 0; columnIndex < columns; ++columnIndex) {
        const double *massOperatorColumn =
            massOperatorBasis.data() + static_cast<std::size_t>(columnIndex) * dimension;
        for (int rowIndex = 0; rowIndex < columns; ++rowIndex) {
            const double *basisColumn =
                basis.data() + static_cast<std::size_t>(rowIndex) * dimension;
            ritzCoordinates[static_cast<std::size_t>(columnIndex) * columns + rowIndex] =
                dot(basisColumn, massOperatorColumn, dimension);
        }
    }
    for (int columnIndex = 0; columnIndex < columns; ++columnIndex) {
        for (int rowIndex = 0; rowIndex < columnIndex; ++rowIndex) {
            const double average = 0.5 * (
                ritzCoordinates[static_cast<std::size_t>(columnIndex) * columns + rowIndex] +
                ritzCoordinates[static_cast<std::size_t>(rowIndex) * columns + columnIndex]);
            ritzCoordinates[static_cast<std::size_t>(columnIndex) * columns + rowIndex] = average;
            ritzCoordinates[static_cast<std::size_t>(rowIndex) * columns + columnIndex] = average;
        }
    }
    return symmetricEigen(ritzCoordinates, columns, ritzValues, message);
}

std::vector<double> linearCombination(
    const std::vector<double> &basis,
    int dimension,
    const double *coordinates,
    int numberOfBasisColumns)
{
    std::vector<double> result(static_cast<std::size_t>(dimension), 0.0);
    for (int basisColumn = 0; basisColumn < numberOfBasisColumns; ++basisColumn) {
        const double coefficient = coordinates[basisColumn];
        const double *basisVector =
            basis.data() + static_cast<std::size_t>(basisColumn) * dimension;
        for (int row = 0; row < dimension; ++row)
            result[static_cast<std::size_t>(row)] += coefficient * basisVector[row];
    }
    return result;
}

double normalizedResidual(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    const std::vector<double> &vector,
    double eigenvalue)
{
    const std::vector<double> stiffnessVector = stiffness.multiply(vector);
    const std::vector<double> massVector = mass.multiply(vector);
    std::vector<double> residual(vector.size(), 0.0);
    for (std::size_t index = 0; index < vector.size(); ++index)
        residual[index] = stiffnessVector[index] - eigenvalue * massVector[index];
    const double denominator = norm2(stiffnessVector) +
        std::fabs(eigenvalue) * norm2(massVector);
    return norm2(residual) / std::max(denominator, std::numeric_limits<double>::min());
}

} // namespace

int clusterSafeRetentionCount(
    const std::vector<double> &ascendingRitzValues,
    int nominalRetained,
    double relativeGapTolerance)
{
    const int numberOfValues = static_cast<int>(ascendingRitzValues.size());
    if (nominalRetained < 1 || nominalRetained > numberOfValues ||
        relativeGapTolerance < 0.0 || !std::isfinite(relativeGapTolerance) ||
        !finiteVector(ascendingRitzValues) ||
        !std::is_sorted(ascendingRitzValues.begin(), ascendingRitzValues.end()))
        return -1;
    int retained = nominalRetained;
    while (retained < numberOfValues) {
        const double retainedValue =
            ascendingRitzValues[static_cast<std::size_t>(numberOfValues - retained)];
        const double nextValue =
            ascendingRitzValues[static_cast<std::size_t>(numberOfValues - retained - 1)];
        const double gap = std::fabs(retainedValue - nextValue) /
            std::max(
                std::numeric_limits<double>::min(),
                std::max(std::fabs(retainedValue), std::fabs(nextValue)));
        if (gap > relativeGapTolerance)
            break;
        ++retained;
    }
    return retained;
}

int solveLocalLanczos(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    int numberOfModes,
    const LocalLanczosControls &controls,
    LocalEigenResult &result,
    std::string &message)
{
    message.clear();
    result = LocalEigenResult{};
    if (stiffness.dimension != mass.dimension || stiffness.validate(message) < 0 ||
        mass.validate(message) < 0 || stiffness.dimension < 1) {
        if (message.empty())
            message = "local Lanczos requires non-empty stiffness and mass of equal size";
        return -1;
    }
    const int dimension = stiffness.dimension;
    if (numberOfModes < 1 || numberOfModes > dimension ||
        !(controls.tolerance > 0.0) || controls.maximumOperatorApplications < 1 ||
        controls.maximumRestarts < 0 || controls.maximumBasisDimension < 0) {
        message = "invalid local Lanczos dimensions or controls";
        return -2;
    }

    {
        MumpsSPD massGuard;
        if (massGuard.factorize(mass, message) < 0) {
            message = "local Lanczos mass is not SPD: " + message;
            return -3;
        }
    }
    MumpsSPD stiffnessFactor;
    if (stiffnessFactor.factorize(stiffness, message) < 0) {
        message = "local Lanczos stiffness is not SPD: " + message;
        return -4;
    }

    const int blockSize = std::min(dimension, std::max(1, std::min(8, numberOfModes)));
    const long long proposedBasis = std::max(
        8LL * numberOfModes + 2LL * blockSize, 80LL);
    const int defaultBasis = static_cast<int>(
        std::min<long long>(dimension, proposedBasis));
    const int maximumBasis = controls.maximumBasisDimension > 0
        ? std::min(dimension, controls.maximumBasisDimension)
        : defaultBasis;
    if (maximumBasis < numberOfModes + 1 && dimension > numberOfModes) {
        message = "local Lanczos basis cap leaves no expansion direction";
        return -5;
    }

    std::vector<double> basis;
    std::vector<double> operatorBasis;
    // M * operatorBasis, kept in step column by column (see rayleighRitz).
    std::vector<double> massOperatorBasis;
    int activeFirstColumn = 0;
    int activeColumns = 0;
    for (int initial = 0; initial < blockSize &&
         result.operatorApplications < controls.maximumOperatorApplications &&
         static_cast<int>(basis.size() / static_cast<std::size_t>(dimension)) < maximumBasis;
         ++initial) {
        std::vector<double> candidate = deterministicVector(
            dimension, controls.level, controls.subdomain, 0, initial);
        const double orthoStarted = MPI_Wtime();
        const int accepted = mOrthonormalize(
            candidate, basis, dimension, mass, result.maximumMOrthogonalityLoss);
        result.orthonormalizeSeconds += MPI_Wtime() - orthoStarted;
        if (accepted < 0) {
            message = "local Lanczos failed during initial M-orthogonalization";
            return -6;
        }
        if (accepted == 0) {
            ++result.breakdowns;
            continue;
        }
        appendColumn(basis, candidate);
        std::vector<double> operatorColumn;
        const double operatorStarted = MPI_Wtime();
        if (applyOperator(mass, stiffnessFactor, candidate, operatorColumn, message) < 0)
            return -7;
        result.operatorSeconds += MPI_Wtime() - operatorStarted;
        appendColumn(operatorBasis, operatorColumn);
        appendColumn(massOperatorBasis, mass.multiply(operatorColumn));
        ++result.operatorApplications;
        ++activeColumns;
    }
    if (activeColumns == 0) {
        message = "local Lanczos could not construct an initial M-orthonormal block";
        return -8;
    }

    std::vector<double> ritzValues;
    std::vector<double> ritzCoordinates;
    while (true) {
        const int basisColumns = static_cast<int>(
            basis.size() / static_cast<std::size_t>(dimension));
        result.maximumBasisUsed = std::max(result.maximumBasisUsed, basisColumns);
        const double ritzStarted = MPI_Wtime();
        if (rayleighRitz(
                basis, massOperatorBasis, dimension,
                ritzValues, ritzCoordinates, message) < 0)
            return -9;
        result.rayleighRitzSeconds += MPI_Wtime() - ritzStarted;
        ++result.rayleighRitzCalls;

        const int availableModes = std::min(numberOfModes, basisColumns);
        result.eigenvalues.assign(static_cast<std::size_t>(availableModes), 0.0);
        result.eigenvectors.assign(
            static_cast<std::size_t>(dimension) * availableModes, 0.0);
        result.normalizedResiduals.assign(static_cast<std::size_t>(availableModes),
                                          std::numeric_limits<double>::infinity());
        int converged = 0;
        for (int mode = 0; mode < availableModes; ++mode) {
            const int ritzColumn = basisColumns - 1 - mode;
            const double inverseEigenvalue = ritzValues[static_cast<std::size_t>(ritzColumn)];
            if (!(inverseEigenvalue > 0.0) || !std::isfinite(inverseEigenvalue))
                continue;
            std::vector<double> vector = linearCombination(
                basis,
                dimension,
                ritzCoordinates.data() + static_cast<std::size_t>(ritzColumn) * basisColumns,
                basisColumns);
            const double eigenvalue = 1.0 / inverseEigenvalue;
            const double residualStarted = MPI_Wtime();
            const double residual = normalizedResidual(stiffness, mass, vector, eigenvalue);
            result.residualSeconds += MPI_Wtime() - residualStarted;
            result.eigenvalues[static_cast<std::size_t>(mode)] = eigenvalue;
            result.normalizedResiduals[static_cast<std::size_t>(mode)] = residual;
            std::copy(
                vector.begin(), vector.end(),
                result.eigenvectors.begin() + static_cast<std::ptrdiff_t>(
                    static_cast<std::size_t>(mode) * dimension));
            if (residual <= controls.tolerance)
                ++converged;
        }
        if (availableModes == numberOfModes && converged == numberOfModes)
            return 0;
        if (result.operatorApplications >= controls.maximumOperatorApplications) {
            const double maximumResidual = result.normalizedResiduals.empty()
                ? std::numeric_limits<double>::infinity()
                : *std::max_element(result.normalizedResiduals.begin(),
                                    result.normalizedResiduals.end());
            const double minimumResidual = result.normalizedResiduals.empty()
                ? std::numeric_limits<double>::infinity()
                : *std::min_element(result.normalizedResiduals.begin(),
                                    result.normalizedResiduals.end());
            std::ostringstream detail;
            detail.setf(std::ios::scientific);
            detail.precision(6);
            detail << "local Lanczos exhausted maximum operator applications: level="
                   << controls.level << ", subdomain=" << controls.subdomain
                   << ", applications=" << result.operatorApplications
                   << ", restarts=" << result.restarts
                   << ", basis=" << basisColumns
                   << ", converged=" << converged << '/' << numberOfModes
                   << ", minResidual=" << minimumResidual
                   << ", maxResidual=" << maximumResidual;
            message = detail.str();
            return -10;
        }

        if (basisColumns >= maximumBasis) {
            if (basisColumns >= dimension) {
                message = "local Lanczos full-space Ritz residual did not meet tolerance";
                return -11;
            }
            if (result.restarts >= controls.maximumRestarts) {
                message = "local Lanczos exhausted maximum restarts";
                return -12;
            }
            const int nominalRetained = std::min(basisColumns, numberOfModes + blockSize);
            const double clusterTolerance = std::sqrt(controls.tolerance);
            const int retained = clusterSafeRetentionCount(
                ritzValues, nominalRetained, clusterTolerance);
            if (retained < 0) {
                message = "local Lanczos received invalid Ritz values at restart";
                return -13;
            }
            if (retained >= maximumBasis) {
                message = "local Lanczos cannot retain an uncut cluster within the basis cap";
                return -14;
            }
            result.maximumRetainedAtRestart = std::max(
                result.maximumRetainedAtRestart, retained);
            std::vector<double> retainedBasis;
            std::vector<double> retainedOperatorBasis;
            retainedBasis.reserve(static_cast<std::size_t>(dimension) * retained);
            retainedOperatorBasis.reserve(static_cast<std::size_t>(dimension) * retained);
            for (int kept = 0; kept < retained; ++kept) {
                const int ritzColumn = basisColumns - 1 - kept;
                const double *coordinates = ritzCoordinates.data() +
                    static_cast<std::size_t>(ritzColumn) * basisColumns;
                appendColumn(
                    retainedBasis,
                    linearCombination(basis, dimension, coordinates, basisColumns));
                appendColumn(
                    retainedOperatorBasis,
                    linearCombination(operatorBasis, dimension, coordinates, basisColumns));
            }
            basis.swap(retainedBasis);
            operatorBasis.swap(retainedOperatorBasis);
            massOperatorBasis.clear();
            for (int kept = 0; kept < retained; ++kept)
                appendColumn(
                    massOperatorBasis,
                    mass.multiply(column(operatorBasis, dimension, kept)));
            activeFirstColumn = 0;
            activeColumns = std::min(blockSize, retained);
            ++result.restarts;
        }

        const int currentColumns = static_cast<int>(
            basis.size() / static_cast<std::size_t>(dimension));
        const int room = maximumBasis - currentColumns;
        const int requestedExpansion = std::min(blockSize, room);
        std::vector<double> candidates;
        for (int columnIndex = 0; columnIndex < activeColumns; ++columnIndex) {
            const int source = activeFirstColumn + columnIndex;
            if (source >= currentColumns)
                break;
            appendColumn(candidates, column(operatorBasis, dimension, source));
        }

        const int expansionStart = currentColumns;
        int acceptedColumns = 0;
        int candidateIndex = 0;
        int deterministicAttempt = 0;
        // M*x for each accepted column, solved as one block below.
        std::vector<double> pendingOperator;
        pendingOperator.reserve(
            static_cast<std::size_t>(dimension) * requestedExpansion);
        // acceptedColumns counts applications not yet added to the running total,
        // so the budget must include them or a block could overshoot the cap.
        while (acceptedColumns < requestedExpansion &&
               result.operatorApplications + acceptedColumns <
                   controls.maximumOperatorApplications) {
            std::vector<double> candidate;
            const int candidateColumns = static_cast<int>(
                candidates.size() / static_cast<std::size_t>(dimension));
            if (candidateIndex < candidateColumns) {
                candidate = column(candidates, dimension, candidateIndex++);
            } else {
                candidate = deterministicVector(
                    dimension,
                    controls.level,
                    controls.subdomain,
                    result.restarts,
                    1000 + deterministicAttempt++);
                ++result.breakdowns;
            }
            const double blockOrthoStarted = MPI_Wtime();
            const int accepted = mOrthonormalize(
                candidate, basis, dimension, mass, result.maximumMOrthogonalityLoss);
            result.orthonormalizeSeconds += MPI_Wtime() - blockOrthoStarted;
            if (accepted < 0) {
                message = "local Lanczos failed during block M-orthogonalization";
                return -15;
            }
            if (accepted == 0) {
                ++result.breakdowns;
                if (deterministicAttempt > 4 * dimension) {
                    message = "local Lanczos could not recover from block breakdown";
                    return -16;
                }
                continue;
            }
            appendColumn(basis, candidate);
            // ADR-1000 P4 section 4: the inverse action is DEFERRED. Each
            // accepted column only needs M*x here; the K^-1 solves for the whole
            // block are issued as ONE multi-RHS MUMPS call below. Measured on the
            // sheet profile, a single-RHS solve cost 0.77 ms against ~31 us for a
            // sparse mat-vec on the same matrix -- i.e. the cost was per-CALL
            // overhead, not arithmetic, so batching removes most of it. Nothing
            // downstream in this loop reads operatorBasis, so deferring is safe:
            // the next block's candidates come from the columns this expansion is
            // about to produce, and they are consumed on the NEXT iteration.
            const double blockOperatorStarted = MPI_Wtime();
            std::vector<double> massCandidate = mass.multiply(candidate);
            result.operatorSeconds += MPI_Wtime() - blockOperatorStarted;
            if (massCandidate.size() != candidate.size()) {
                message = "mass multiplication failed in the local Lanczos operator";
                return -17;
            }
            pendingOperator.insert(
                pendingOperator.end(), massCandidate.begin(), massCandidate.end());
            ++acceptedColumns;
        }
        if (acceptedColumns == 0) {
            message = "local Lanczos added no vector during expansion";
            return -18;
        }
        {
            const double blockSolveStarted = MPI_Wtime();
            if (stiffnessFactor.solve(pendingOperator, acceptedColumns, message) < 0) {
                message = "local Lanczos inverse action failed: " + message;
                return -17;
            }
            result.operatorSeconds += MPI_Wtime() - blockSolveStarted;
            operatorBasis.insert(
                operatorBasis.end(), pendingOperator.begin(), pendingOperator.end());
            for (int added = 0; added < acceptedColumns; ++added)
                appendColumn(
                    massOperatorBasis,
                    mass.multiply(column(pendingOperator, dimension, added)));
            result.operatorApplications += acceptedColumns;
            ++result.blockSolves;
            pendingOperator.clear();
        }
        activeFirstColumn = expansionStart;
        activeColumns = acceptedColumns;
    }
}

} // namespace ladruno_cms
