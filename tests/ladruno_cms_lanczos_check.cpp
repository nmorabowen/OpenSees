#include "LadrunoCMSLocalLanczos.h"
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
            if (value == 0.0)
                continue;
            result.columnIndices.push_back(column);
            result.upperValues.push_back(value);
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    return result;
}

std::vector<double> multiply(
    const std::vector<double> &left,
    const std::vector<double> &right,
    int rows,
    int inner,
    int columns,
    bool transposeRight = false)
{
    std::vector<double> result(static_cast<std::size_t>(rows) * columns, 0.0);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            for (int index = 0; index < inner; ++index) {
                const double rightValue = transposeRight
                    ? right[static_cast<std::size_t>(column) * inner + index]
                    : right[static_cast<std::size_t>(index) * columns + column];
                result[static_cast<std::size_t>(row) * columns + column] +=
                    left[static_cast<std::size_t>(row) * inner + index] * rightValue;
            }
        }
    }
    return result;
}

struct Pencil {
    int order = 0;
    std::vector<double> stiffness;
    std::vector<double> mass;
};

Pencil transformedDiagonalPencil(const std::vector<double> &eigenvalues)
{
    const int order = static_cast<int>(eigenvalues.size());
    std::vector<double> lower(static_cast<std::size_t>(order) * order, 0.0);
    for (int row = 0; row < order; ++row) {
        lower[static_cast<std::size_t>(row) * order + row] = 1.0 + 0.02 * row;
        if (row > 0)
            lower[static_cast<std::size_t>(row) * order + row - 1] = 0.08;
        if (row > 2)
            lower[static_cast<std::size_t>(row) * order + row - 3] = -0.03;
    }
    std::vector<double> diagonal(static_cast<std::size_t>(order) * order, 0.0);
    for (int index = 0; index < order; ++index)
        diagonal[static_cast<std::size_t>(index) * order + index] = eigenvalues[index];
    Pencil result;
    result.order = order;
    result.mass = multiply(lower, lower, order, order, order, true);
    const std::vector<double> lowerDiagonal =
        multiply(lower, diagonal, order, order, order);
    result.stiffness = multiply(lowerDiagonal, lower, order, order, order, true);
    return result;
}

std::vector<double> lapackEigenvalues(const Pencil &pencil, int requested)
{
    const int order = pencil.order;
    std::vector<double> stiffnessColumnMajor(static_cast<std::size_t>(order) * order);
    std::vector<double> massColumnMajor(static_cast<std::size_t>(order) * order);
    for (int row = 0; row < order; ++row) {
        for (int column = 0; column < order; ++column) {
            stiffnessColumnMajor[static_cast<std::size_t>(column) * order + row] =
                pencil.stiffness[static_cast<std::size_t>(row) * order + column];
            massColumnMajor[static_cast<std::size_t>(column) * order + row] =
                pencil.mass[static_cast<std::size_t>(row) * order + column];
        }
    }
    int itype = 1;
    char jobz = 'V';
    char range = 'I';
    char uplo = 'U';
    int n = order;
    int leading = order;
    double lowerValue = 0.0;
    double upperValue = 0.0;
    int first = 1;
    int last = requested;
    double absoluteTolerance = 0.0;
    int found = 0;
    std::vector<double> values(static_cast<std::size_t>(order), 0.0);
    std::vector<double> vectors(static_cast<std::size_t>(order) * requested, 0.0);
    int leadingVectors = order;
    int workspaceSize = std::max(1, 8 * order);
    std::vector<double> workspace(static_cast<std::size_t>(workspaceSize));
    std::vector<int> integerWorkspace(static_cast<std::size_t>(5 * order));
    std::vector<int> failed(static_cast<std::size_t>(order));
    int info = 0;
    LADRUNO_TEST_DSYGVX(
        &itype, &jobz, &range, &uplo, &n,
        stiffnessColumnMajor.data(), &leading,
        massColumnMajor.data(), &leading,
        &lowerValue, &upperValue, &first, &last, &absoluteTolerance,
        &found, values.data(), vectors.data(), &leadingVectors,
        workspace.data(), &workspaceSize, integerWorkspace.data(), failed.data(), &info);
    require(info == 0 && found == requested, "LAPACK reference eigensolve failed");
    values.resize(static_cast<std::size_t>(std::max(0, found)));
    return values;
}

void checkMOrthogonality(
    const ladruno_cms::LocalEigenResult &result,
    const ladruno_cms::SymmetricCSR &mass,
    int order,
    int modes)
{
    for (int column = 0; column < modes; ++column) {
        const std::vector<double> vector(
            result.eigenvectors.begin() + static_cast<std::ptrdiff_t>(column * order),
            result.eigenvectors.begin() + static_cast<std::ptrdiff_t>((column + 1) * order));
        const std::vector<double> massVector = mass.multiply(vector);
        for (int row = 0; row < modes; ++row) {
            double product = 0.0;
            for (int index = 0; index < order; ++index)
                product += result.eigenvectors[static_cast<std::size_t>(row) * order + index] *
                    massVector[static_cast<std::size_t>(index)];
            require(
                close(product, row == column ? 1.0 : 0.0, 2.0e-7),
                "Lanczos modes are not M-orthonormal");
        }
    }
}

void testDefaultBasisAgainstLapack()
{
    const Pencil pencil = transformedDiagonalPencil(
        {1.0, 1.0 + 1.0e-9, 2.0, 3.0, 4.0, 5.0,
         6.0, 7.0, 8.0, 9.0, 10.0, 11.0});
    const ladruno_cms::SymmetricCSR stiffness = csrFromDense(pencil.stiffness, pencil.order);
    const ladruno_cms::SymmetricCSR mass = csrFromDense(pencil.mass, pencil.order);
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-10;
    controls.maximumOperatorApplications = 300;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    require(
        ladruno_cms::solveLocalLanczos(
            stiffness, mass, 4, controls, result, message) == 0,
        "default local Lanczos failed: " + message);
    const std::vector<double> reference = lapackEigenvalues(pencil, 4);
    require(result.eigenvalues.size() == reference.size(), "Lanczos mode count differs from LAPACK");
    for (std::size_t mode = 0; mode < reference.size() && mode < result.eigenvalues.size(); ++mode) {
        require(close(result.eigenvalues[mode], reference[mode], 2.0e-8),
                "Lanczos eigenvalue differs from LAPACK");
        require(result.normalizedResiduals[mode] <= controls.tolerance,
                "Lanczos accepted a residual above tolerance");
    }
    checkMOrthogonality(result, mass, pencil.order, 4);
}

void testForcedThickRestart()
{
    std::vector<double> spectrum(30);
    for (int index = 0; index < 30; ++index)
        spectrum[static_cast<std::size_t>(index)] = 1.0 + 0.15 * index;
    // k+b retains four Ritz vectors; place a cluster across the 4/5 boundary.
    spectrum[4] = spectrum[3] * (1.0 + 1.0e-10);
    const Pencil pencil = transformedDiagonalPencil(spectrum);
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-8;
    controls.maximumOperatorApplications = 400;
    controls.maximumBasisDimension = 8;
    controls.maximumRestarts = 20;
    controls.level = 2;
    controls.subdomain = 7;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    const ladruno_cms::SymmetricCSR stiffness = csrFromDense(pencil.stiffness, pencil.order);
    const ladruno_cms::SymmetricCSR mass = csrFromDense(pencil.mass, pencil.order);
    require(
        ladruno_cms::solveLocalLanczos(
            stiffness, mass, 2, controls, result, message) == 0,
        "forced-restart local Lanczos failed: " + message);
    require(result.restarts > 0, "forced-restart fixture did not execute a restart");
    require(result.maximumBasisUsed <= controls.maximumBasisDimension,
            "Lanczos exceeded its basis cap");
    const std::vector<double> reference = lapackEigenvalues(pencil, 2);
    for (int mode = 0; mode < 2; ++mode)
        require(close(result.eigenvalues[static_cast<std::size_t>(mode)],
                      reference[static_cast<std::size_t>(mode)], 5.0e-7),
                "restarted Lanczos eigenvalue differs from LAPACK");
    checkMOrthogonality(result, mass, pencil.order, 2);
}

void testClusterSafeRetentionRule()
{
    // Ascending inverse Ritz values. The nominal top four end at 0.5; the
    // adjacent 0.5*(1-1e-10) must be retained with a 1e-5 relative-gap rule.
    const std::vector<double> ritz = {
        0.1, 0.2, 0.5 * (1.0 - 1.0e-10), 0.5, 0.7, 0.8, 0.9};
    require(
        ladruno_cms::clusterSafeRetentionCount(ritz, 4, 1.0e-5) == 5,
        "thick-restart retention cut a Ritz cluster");
    require(
        ladruno_cms::clusterSafeRetentionCount(ritz, 4, 1.0e-12) == 4,
        "thick-restart retention merged a resolved Ritz gap");
}

void testBreakdownFailsCleanlyAtImpossibleTolerance()
{
    std::vector<double> spectrum(10, 2.0);
    spectrum[0] = 1.0;
    const Pencil pencil = transformedDiagonalPencil(spectrum);
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-30;
    controls.maximumOperatorApplications = 30;
    controls.maximumBasisDimension = 8;
    controls.maximumRestarts = 2;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    require(
        ladruno_cms::solveLocalLanczos(
            csrFromDense(pencil.stiffness, pencil.order),
            csrFromDense(pencil.mass, pencil.order),
            1,
            controls,
            result,
            message) < 0,
        "impossible Lanczos tolerance was silently accepted");
    require(result.breakdowns > 0,
            "invariant-subspace fixture did not exercise breakdown recovery");
    require(message.find("exhausted") != std::string::npos ||
                message.find("full-space") != std::string::npos ||
                message.find("cluster") != std::string::npos,
            "breakdown path did not return an explicit terminal diagnostic");
}

void testCoordinateCondensationFlow()
{
    const Pencil pencil{
        4,
        {6.0, 1.0, 1.0, 0.5,
         1.0, 5.0, 2.0, -1.0,
         1.0, 2.0, 4.0, 1.0,
         0.5, -1.0, 1.0, 3.0},
        {2.0, 0.0, 0.0, 0.0,
         0.0, 3.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0}};
    ladruno_cms::StaticCondensation condensation;
    std::string message;
    require(
        ladruno_cms::condenseCoordinateMass(
            csrFromDense(pencil.stiffness, 4), csrFromDense(pencil.mass, 4),
            1.0e-12, 1.0e-14, condensation, message) == 0,
        "pre-Lanczos coordinate condensation failed");
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-10;
    ladruno_cms::LocalEigenResult reduced;
    require(
        ladruno_cms::solveLocalLanczos(
            condensation.reducedStiffness, condensation.reducedMass,
            2, controls, reduced, message) == 0,
        "Lanczos on condensed pencil failed: " + message);
    std::vector<double> fullVectors;
    require(
        ladruno_cms::reconstructStaticCoordinates(
            condensation, reduced.eigenvectors, 2, fullVectors, message) == 0,
        "post-Lanczos static reconstruction failed");
    const ladruno_cms::SymmetricCSR fullStiffness = csrFromDense(pencil.stiffness, 4);
    const ladruno_cms::SymmetricCSR fullMass = csrFromDense(pencil.mass, 4);
    for (int mode = 0; mode < 2; ++mode) {
        const std::vector<double> vector(
            fullVectors.begin() + static_cast<std::ptrdiff_t>(mode * 4),
            fullVectors.begin() + static_cast<std::ptrdiff_t>((mode + 1) * 4));
        const std::vector<double> stiffnessVector = fullStiffness.multiply(vector);
        const std::vector<double> massVector = fullMass.multiply(vector);
        double residualNorm = 0.0;
        double denominator = 0.0;
        for (int row = 0; row < 4; ++row) {
            const double residual = stiffnessVector[static_cast<std::size_t>(row)] -
                reduced.eigenvalues[static_cast<std::size_t>(mode)] *
                    massVector[static_cast<std::size_t>(row)];
            residualNorm += residual * residual;
            denominator += stiffnessVector[static_cast<std::size_t>(row)] *
                stiffnessVector[static_cast<std::size_t>(row)];
        }
        require(std::sqrt(residualNorm / denominator) < 5.0e-10,
                "reconstructed condensed eigenvector has a large full residual");
    }
}

void testWideSpectrumAgainstLapack()
{
    std::vector<double> spectrum = {
        1.0e-4, 1.00001e-4, 5.0e-4, 0.1, 0.5, 1.0, 2.0, 4.0,
        8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0};
    const Pencil pencil = transformedDiagonalPencil(spectrum);
    const ladruno_cms::SymmetricCSR stiffness = csrFromDense(pencil.stiffness, pencil.order);
    const ladruno_cms::SymmetricCSR mass = csrFromDense(pencil.mass, pencil.order);
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-9;
    controls.maximumOperatorApplications = 300;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    require(
        ladruno_cms::solveLocalLanczos(
            stiffness, mass, 3, controls, result, message) == 0,
        "wide-spectrum local Lanczos failed: " + message);
    const std::vector<double> reference = lapackEigenvalues(pencil, 3);
    for (int mode = 0; mode < 3; ++mode) {
        require(close(result.eigenvalues[static_cast<std::size_t>(mode)],
                      reference[static_cast<std::size_t>(mode)], 2.0e-7),
                "wide-spectrum Lanczos eigenvalue differs from LAPACK");
        require(result.normalizedResiduals[static_cast<std::size_t>(mode)] <=
                    controls.tolerance,
                "wide-spectrum Lanczos residual exceeds tolerance");
    }
}

void testInvalidMassRejected()
{
    const ladruno_cms::SymmetricCSR stiffness = csrFromDense({2.0, 0.0, 0.0, 3.0}, 2);
    const ladruno_cms::SymmetricCSR mass = csrFromDense({1.0, 2.0, 2.0, 1.0}, 2);
    ladruno_cms::LocalLanczosControls controls;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    require(
        ladruno_cms::solveLocalLanczos(
            stiffness, mass, 1, controls, result, message) < 0,
        "local Lanczos accepted an indefinite mass");
}

void testOperatorBudgetIsStrict()
{
    const Pencil pencil = transformedDiagonalPencil({1.0, 2.0, 3.0, 4.0, 5.0});
    ladruno_cms::LocalLanczosControls controls;
    controls.tolerance = 1.0e-12;
    controls.maximumOperatorApplications = 1;
    ladruno_cms::LocalEigenResult result;
    std::string message;
    require(
        ladruno_cms::solveLocalLanczos(
            csrFromDense(pencil.stiffness, pencil.order),
            csrFromDense(pencil.mass, pencil.order),
            3,
            controls,
            result,
            message) < 0,
        "local Lanczos ignored an insufficient operator budget");
    require(result.operatorApplications <= controls.maximumOperatorApplications,
            "local Lanczos exceeded maximumOperatorApplications");
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    {
        testDefaultBasisAgainstLapack();
        testForcedThickRestart();
        testClusterSafeRetentionRule();
        testBreakdownFailsCleanlyAtImpossibleTolerance();
        testCoordinateCondensationFlow();
        testWideSpectrumAgainstLapack();
        testInvalidMassRejected();
        testOperatorBudgetIsStrict();
    }
    MPI_Finalize();
    if (failures != 0)
        return 1;
    std::cout << "Ladruno CMS local Lanczos checks passed\n";
    return 0;
}
