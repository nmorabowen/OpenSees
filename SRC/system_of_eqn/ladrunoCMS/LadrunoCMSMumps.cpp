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

#include "LadrunoCMSMumps.h"

#include <dmumps_c.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>
#include <type_traits>
#include <utility>

namespace ladruno_cms {
namespace {

static_assert(std::is_integral<MUMPS_INT>::value, "MUMPS_INT must be integral");
static_assert(
    std::numeric_limits<MUMPS_INT>::max() >= std::numeric_limits<int>::max(),
    "the linked MUMPS integer width cannot represent OpenSees equation indices");

bool checkedProduct(std::size_t left, std::size_t right, std::size_t &result)
{
    if (left != 0u && right > std::numeric_limits<std::size_t>::max() / left)
        return false;
    result = left * right;
    return true;
}

int csrFromEntries(
    int dimension,
    const std::map<std::pair<int, int>, double> &entries,
    SymmetricCSR &result,
    std::string &message)
{
    result = SymmetricCSR{};
    if (dimension < 0 ||
        entries.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        message = "reduced CSR dimensions exceed the integer index range";
        return -1;
    }
    result.dimension = dimension;
    result.rowOffsets.assign(static_cast<std::size_t>(dimension) + 1u, 0);
    for (const auto &entry : entries) {
        const int row = entry.first.first;
        const int column = entry.first.second;
        if (row < 0 || row > column || column >= dimension || !std::isfinite(entry.second)) {
            message = "invalid entry while constructing a reduced CSR matrix";
            return -2;
        }
        ++result.rowOffsets[static_cast<std::size_t>(row) + 1u];
    }
    for (int row = 0; row < dimension; ++row)
        result.rowOffsets[static_cast<std::size_t>(row) + 1u] +=
            result.rowOffsets[static_cast<std::size_t>(row)];
    result.columnIndices.reserve(entries.size());
    result.upperValues.reserve(entries.size());
    for (const auto &entry : entries) {
        result.columnIndices.push_back(entry.first.second);
        result.upperValues.push_back(entry.second);
    }
    return result.validate(message);
}

int splitPrincipalBlocks(
    const SymmetricCSR &matrix,
    const std::vector<int> &dynamicMap,
    const std::vector<int> &masslessMap,
    std::map<std::pair<int, int>, double> &dynamicEntries,
    std::map<std::pair<int, int>, double> &masslessEntries,
    std::map<std::pair<int, int>, double> &masslessDynamicEntries,
    std::string &message)
{
    for (int row = 0; row < matrix.dimension; ++row) {
        for (int position = matrix.rowOffsets[static_cast<std::size_t>(row)];
             position < matrix.rowOffsets[static_cast<std::size_t>(row) + 1u]; ++position) {
            const int column = matrix.columnIndices[static_cast<std::size_t>(position)];
            const double value = matrix.upperValues[static_cast<std::size_t>(position)];
            const int dynamicRow = dynamicMap[static_cast<std::size_t>(row)];
            const int dynamicColumn = dynamicMap[static_cast<std::size_t>(column)];
            const int masslessRow = masslessMap[static_cast<std::size_t>(row)];
            const int masslessColumn = masslessMap[static_cast<std::size_t>(column)];
            if (dynamicRow >= 0 && dynamicColumn >= 0) {
                dynamicEntries[std::make_pair(dynamicRow, dynamicColumn)] += value;
            } else if (masslessRow >= 0 && masslessColumn >= 0) {
                masslessEntries[std::make_pair(masslessRow, masslessColumn)] += value;
            } else if (masslessRow >= 0 && dynamicColumn >= 0) {
                masslessDynamicEntries[
                    std::make_pair(masslessRow, dynamicColumn)] += value;
            } else if (dynamicRow >= 0 && masslessColumn >= 0) {
                masslessDynamicEntries[
                    std::make_pair(masslessColumn, dynamicRow)] += value;
            } else {
                message = "coordinate split did not classify a matrix entry";
                return -2;
            }
        }
    }
    return 0;
}

} // namespace

struct MumpsSPD::Impl {
    DMUMPS_STRUC_C control;
    std::vector<MUMPS_INT> rows;
    std::vector<MUMPS_INT> columns;
    std::vector<double> values;
    int order = 0;
    bool initialized = false;
    bool factorized = false;

    Impl()
    {
        std::memset(&control, 0, sizeof(control));
    }

    void terminate()
    {
        if (!initialized)
            return;
        int mpiInitialized = 0;
        int mpiFinalized = 0;
        MPI_Initialized(&mpiInitialized);
        if (mpiInitialized)
            MPI_Finalized(&mpiFinalized);
        if (mpiInitialized && !mpiFinalized) {
            control.job = -2;
            dmumps_c(&control);
        }
        initialized = false;
        factorized = false;
    }
};

struct DistributedMumpsSPD::Impl {
    DMUMPS_STRUC_C control;
    std::vector<MUMPS_INT> rows;
    std::vector<MUMPS_INT> columns;
    std::vector<double> values;
    int order = 0;
    bool initialized = false;
    bool factorized = false;

    Impl()
    {
        std::memset(&control, 0, sizeof(control));
    }

    void terminate()
    {
        if (!initialized)
            return;
        int mpiInitialized = 0;
        int mpiFinalized = 0;
        MPI_Initialized(&mpiInitialized);
        if (mpiInitialized)
            MPI_Finalized(&mpiFinalized);
        if (mpiInitialized && !mpiFinalized) {
            control.job = -2;
            dmumps_c(&control);
        }
        initialized = false;
        factorized = false;
    }
};

struct StaticCondensation::Impl {
    std::vector<int> crossRowOffsets;
    std::vector<int> crossColumns;
    std::vector<double> crossValues;
    std::unique_ptr<MumpsSPD> masslessStiffnessFactor;
};

StaticCondensation::StaticCondensation()
    : impl_(new Impl())
{
}

StaticCondensation::~StaticCondensation() = default;
StaticCondensation::StaticCondensation(StaticCondensation &&) noexcept = default;
StaticCondensation &StaticCondensation::operator=(StaticCondensation &&) noexcept = default;

std::size_t StaticCondensation::persistentCrossNonzeros() const
{
    return impl_ ? impl_->crossValues.size() : 0u;
}

MumpsSPD::MumpsSPD()
    : impl_(new Impl())
{
}

MumpsSPD::~MumpsSPD()
{
    impl_->terminate();
}

DistributedMumpsSPD::DistributedMumpsSPD()
    : impl_(new Impl())
{
}

DistributedMumpsSPD::~DistributedMumpsSPD()
{
    impl_->terminate();
}

int DistributedMumpsSPD::factorize(
    const SymmetricCSR &localMatrix,
    std::string &message)
{
    message.clear();
    impl_->terminate();
    impl_->rows.clear();
    impl_->columns.clear();
    impl_->values.clear();
    impl_->order = 0;

    int localFailure = localMatrix.validate(message) < 0 || localMatrix.dimension < 1;
    int minimumDimension = localMatrix.dimension;
    int maximumDimension = localMatrix.dimension;
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &minimumDimension, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maximumDimension, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (globalFailure || minimumDimension != maximumDimension) {
        if (!localFailure)
            message = minimumDimension != maximumDimension
                ? "distributed MUMPS matrix dimensions differ between ranks"
                : "distributed MUMPS received an invalid local matrix on another rank";
        return -1;
    }
    if (localMatrix.columnIndices.size() >
        static_cast<std::size_t>(std::numeric_limits<MUMPS_INT>::max())) {
        localFailure = 1;
        message = "distributed MUMPS local nonzero count exceeds its integer width";
    }
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (globalFailure) {
        if (!localFailure)
            message = "distributed MUMPS integer range failed on another rank";
        return -2;
    }

    impl_->rows.reserve(localMatrix.columnIndices.size());
    impl_->columns.reserve(localMatrix.columnIndices.size());
    impl_->values.reserve(localMatrix.upperValues.size());
    for (int row = 0; row < localMatrix.dimension; ++row) {
        for (int position = localMatrix.rowOffsets[static_cast<std::size_t>(row)];
             position < localMatrix.rowOffsets[static_cast<std::size_t>(row + 1)];
             ++position) {
            const int column = localMatrix.columnIndices[static_cast<std::size_t>(position)];
            impl_->rows.push_back(static_cast<MUMPS_INT>(column) + 1);
            impl_->columns.push_back(static_cast<MUMPS_INT>(row) + 1);
            impl_->values.push_back(localMatrix.upperValues[static_cast<std::size_t>(position)]);
        }
    }

    std::memset(&impl_->control, 0, sizeof(impl_->control));
    impl_->control.par = 1;
    impl_->control.sym = 1;
    impl_->control.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    impl_->control.job = -1;
    dmumps_c(&impl_->control);
    if (impl_->control.infog[0] < 0) {
        message = "distributed MUMPS initialization failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]);
        return -3;
    }
    impl_->initialized = true;
    impl_->control.icntl[0] = -1;
    impl_->control.icntl[1] = -1;
    impl_->control.icntl[2] = -1;
    impl_->control.icntl[3] = 0;
    // ICNTL(7)=0 (AMD), NOT 7 (auto). Under auto, MUMPS selects PORD for a
    // DENSE-as-CSR pencil -- exactly what reduceCraigBampton feeds here -- and
    // PORD dies in ANALYSIS with "Error in function orderMinPriority: no valid
    // number of stages in multisector". Reproduced deterministically by the
    // order-12000 dense fixture in ladruno_cms_mumps_check at BOTH np=2 and
    // np=4; a sparse order-64000 pencil at the same rank counts is unaffected,
    // so the trigger is the dense pattern, not the process count. AMD has no
    // multisector concept and is always built into MUMPS. ICNTL(28) (sequential
    // vs parallel analysis) was measured to be irrelevant: forcing
    // ICNTL(28)=1 while leaving ICNTL(7)=7 still fails identically.
    // ADR-1000 Part 0.
    impl_->control.icntl[6] = 0;
    impl_->control.icntl[12] = 1;
    impl_->control.icntl[13] = 35;
    impl_->control.icntl[17] = 3; // distributed assembled matrix
    impl_->control.icntl[23] = 1;
    impl_->control.n = static_cast<MUMPS_INT>(localMatrix.dimension);
    impl_->control.nz_loc = static_cast<MUMPS_INT>(impl_->values.size());
    impl_->control.irn_loc = impl_->rows.data();
    impl_->control.jcn_loc = impl_->columns.data();
    impl_->control.a_loc = impl_->values.data();
    impl_->control.job = 4;
    dmumps_c(&impl_->control);
    int attempts = 0;
    while ((impl_->control.infog[0] == -8 || impl_->control.infog[0] == -9) &&
           attempts < 5) {
        const MUMPS_INT doubled =
            impl_->control.icntl[13] * static_cast<MUMPS_INT>(2);
        impl_->control.icntl[13] = std::min(
            static_cast<MUMPS_INT>(4000),
            std::max(static_cast<MUMPS_INT>(100), doubled));
        impl_->control.job = 4;
        dmumps_c(&impl_->control);
        ++attempts;
    }
    if (impl_->control.infog[0] < 0) {
        message = "distributed MUMPS SPD factorization failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]);
        impl_->terminate();
        return -4;
    }
    if (impl_->control.infog[11] != 0 || impl_->control.infog[27] != 0) {
        message = "distributed MUMPS rejected the SPD guard: negative pivots=" +
            std::to_string(impl_->control.infog[11]) + " null pivots=" +
            std::to_string(impl_->control.infog[27]);
        impl_->terminate();
        return -5;
    }
    impl_->order = localMatrix.dimension;
    impl_->factorized = true;
    return 0;
}

int DistributedMumpsSPD::solve(
    std::vector<double> &rightHandSides,
    int numberOfColumns,
    std::string &message)
{
    message.clear();
    if (!impl_->factorized || numberOfColumns < 1) {
        message = "distributed MUMPS solve requires a factorized matrix and RHS";
        return -1;
    }
    std::size_t expected = 0;
    if (!checkedProduct(static_cast<std::size_t>(impl_->order),
                        static_cast<std::size_t>(numberOfColumns), expected) ||
        rightHandSides.size() != expected || expected >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        message = "distributed MUMPS RHS dimensions exceed the supported range";
        return -2;
    }
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    impl_->control.nrhs = static_cast<MUMPS_INT>(numberOfColumns);
    impl_->control.lrhs = static_cast<MUMPS_INT>(impl_->order);
    impl_->control.rhs = rank == 0 ? rightHandSides.data() : nullptr;
    impl_->control.job = 3;
    dmumps_c(&impl_->control);
    if (impl_->control.infog[0] < 0) {
        message = "distributed MUMPS solve failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]);
        return -3;
    }
    MPI_Bcast(rightHandSides.data(), static_cast<int>(expected), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    return 0;
}

int DistributedMumpsSPD::dimension() const
{
    return impl_->order;
}

bool DistributedMumpsSPD::isFactorized() const
{
    return impl_->factorized;
}

int MumpsSPD::factorize(const SymmetricCSR &matrix, std::string &message)
{
    message.clear();
    impl_->terminate();
    impl_->rows.clear();
    impl_->columns.clear();
    impl_->values.clear();
    impl_->order = 0;

    if (matrix.validate(message) < 0 || matrix.dimension < 1) {
        if (message.empty())
            message = "MUMPS requires a non-empty symmetric CSR matrix";
        return -1;
    }
    int mpiInitialized = 0;
    int mpiFinalized = 0;
    MPI_Initialized(&mpiInitialized);
    if (mpiInitialized)
        MPI_Finalized(&mpiFinalized);
    if (!mpiInitialized || mpiFinalized) {
        message = "MUMPS factorization requires an active MPI runtime";
        return -2;
    }
    if (matrix.columnIndices.size() >
        static_cast<std::size_t>(std::numeric_limits<MUMPS_INT>::max())) {
        message = "MUMPS nonzero count exceeds its integer width";
        return -3;
    }

    impl_->rows.reserve(matrix.columnIndices.size());
    impl_->columns.reserve(matrix.columnIndices.size());
    impl_->values.reserve(matrix.upperValues.size());
    for (int row = 0; row < matrix.dimension; ++row) {
        for (int position = matrix.rowOffsets[static_cast<std::size_t>(row)];
             position < matrix.rowOffsets[static_cast<std::size_t>(row) + 1u]; ++position) {
            const int column = matrix.columnIndices[static_cast<std::size_t>(position)];
            // CSR is upper triangular; MUMPS receives the equivalent lower triangle.
            impl_->rows.push_back(static_cast<MUMPS_INT>(column) + 1);
            impl_->columns.push_back(static_cast<MUMPS_INT>(row) + 1);
            impl_->values.push_back(matrix.upperValues[static_cast<std::size_t>(position)]);
        }
    }

    std::memset(&impl_->control, 0, sizeof(impl_->control));
    impl_->control.par = 1;
    impl_->control.sym = 1;
    impl_->control.comm_fortran = MPI_Comm_c2f(MPI_COMM_SELF);
    impl_->control.job = -1;
    dmumps_c(&impl_->control);
    if (impl_->control.infog[0] < 0) {
        message = "MUMPS initialization failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]);
        return -4;
    }
    impl_->initialized = true;

    impl_->control.icntl[0] = -1;
    impl_->control.icntl[1] = -1;
    impl_->control.icntl[2] = -1;
    impl_->control.icntl[3] = 0;
    impl_->control.icntl[4] = 0;
    // ICNTL(7)=0 (AMD) for the same reason as the distributed factorization
    // above -- and this site matters more: the rank-local Craig-Bampton pencils
    // are what reach MumpsSPD, and they are the DENSE-as-CSR ones. Measured on
    // the order-12000 dense fixture: with ICNTL(7)=7 this path dies in PORD
    // ("orderMinPriority ... multisector") even on a single rank, while the
    // distributed path next to it (already switched to AMD) completes.
    // ADR-1000 Part 0.
    impl_->control.icntl[6] = 0;
    impl_->control.icntl[12] = 1; // exact inertia/null-pivot behavior on the root
    impl_->control.icntl[13] = 35;
    impl_->control.icntl[17] = 0;
    impl_->control.icntl[23] = 1; // detect null pivot rows
    impl_->control.n = static_cast<MUMPS_INT>(matrix.dimension);
    impl_->control.nz = static_cast<MUMPS_INT>(impl_->values.size());
    impl_->control.irn = impl_->rows.data();
    impl_->control.jcn = impl_->columns.data();
    impl_->control.a = impl_->values.data();

    impl_->control.job = 4;
    dmumps_c(&impl_->control);
    int attempts = 0;
    while ((impl_->control.infog[0] == -8 || impl_->control.infog[0] == -9) &&
           attempts < 5) {
        const MUMPS_INT doubled = impl_->control.icntl[13] * static_cast<MUMPS_INT>(2);
        impl_->control.icntl[13] = std::min(
            static_cast<MUMPS_INT>(4000),
            std::max(static_cast<MUMPS_INT>(100), doubled));
        impl_->control.job = 4;
        dmumps_c(&impl_->control);
        ++attempts;
    }
    if (impl_->control.infog[0] < 0) {
        message = "MUMPS analysis/factorization failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]) + " INFOG(2)=" +
            std::to_string(impl_->control.infog[1]);
        impl_->terminate();
        return -5;
    }
    if (impl_->control.infog[11] != 0 || impl_->control.infog[27] != 0) {
        message = "MUMPS rejected the SPD guard: negative pivots=" +
            std::to_string(impl_->control.infog[11]) + " null pivots=" +
            std::to_string(impl_->control.infog[27]);
        impl_->terminate();
        return -6;
    }
    impl_->order = matrix.dimension;
    impl_->factorized = true;
    return 0;
}

int MumpsSPD::solve(
    std::vector<double> &rightHandSides,
    int numberOfColumns,
    std::string &message)
{
    message.clear();
    if (!impl_->factorized || numberOfColumns < 1) {
        message = "MUMPS solve requires a factorization and at least one RHS";
        return -1;
    }
    std::size_t expected = 0;
    if (!checkedProduct(
            static_cast<std::size_t>(impl_->order),
            static_cast<std::size_t>(numberOfColumns),
            expected) ||
        rightHandSides.size() != expected) {
        message = "MUMPS RHS dimensions are inconsistent";
        return -2;
    }
    if (!std::all_of(
            rightHandSides.begin(), rightHandSides.end(),
            [](double value) { return std::isfinite(value); })) {
        message = "MUMPS RHS contains a non-finite value";
        return -3;
    }
    impl_->control.rhs = rightHandSides.data();
    impl_->control.nrhs = static_cast<MUMPS_INT>(numberOfColumns);
    impl_->control.lrhs = static_cast<MUMPS_INT>(impl_->order);
    impl_->control.job = 3;
    dmumps_c(&impl_->control);
    if (impl_->control.infog[0] < 0) {
        message = "MUMPS solve failed with INFOG(1)=" +
            std::to_string(impl_->control.infog[0]);
        return -4;
    }
    if (!std::all_of(
            rightHandSides.begin(), rightHandSides.end(),
            [](double value) { return std::isfinite(value); })) {
        message = "MUMPS solve produced a non-finite value";
        return -5;
    }
    return 0;
}

int MumpsSPD::dimension() const
{
    return impl_->order;
}

bool MumpsSPD::isFactorized() const
{
    return impl_->factorized;
}

int condenseCoordinateMass(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    double relativeTolerance,
    double absoluteTolerance,
    StaticCondensation &result,
    std::string &message)
{
    message.clear();
    result = StaticCondensation{};
    if (stiffness.dimension != mass.dimension || stiffness.validate(message) < 0 ||
        mass.validate(message) < 0) {
        if (message.empty())
            message = "stiffness and mass dimensions differ";
        return -1;
    }
    CoordinateMassClassification classification;
    if (classifyCoordinateMass(
            mass, relativeTolerance, absoluteTolerance, classification, message) < 0)
        return -2;

    result.originalDimension = stiffness.dimension;
    result.dynamic = classification.dynamic;
    result.massless = classification.massless;
    std::vector<int> dynamicMap(static_cast<std::size_t>(stiffness.dimension), -1);
    std::vector<int> masslessMap(static_cast<std::size_t>(stiffness.dimension), -1);
    for (std::size_t index = 0; index < result.dynamic.size(); ++index)
        dynamicMap[static_cast<std::size_t>(result.dynamic[index])] = static_cast<int>(index);
    for (std::size_t index = 0; index < result.massless.size(); ++index)
        masslessMap[static_cast<std::size_t>(result.massless[index])] = static_cast<int>(index);

    std::map<std::pair<int, int>, double> stiffnessDD;
    std::map<std::pair<int, int>, double> stiffnessZZ;
    std::map<std::pair<int, int>, double> stiffnessZD;
    if (splitPrincipalBlocks(
            stiffness,
            dynamicMap,
            masslessMap,
            stiffnessDD,
            stiffnessZZ,
            stiffnessZD,
            message) < 0)
        return -3;

    std::map<std::pair<int, int>, double> massDD;
    std::map<std::pair<int, int>, double> ignoredMassZZ;
    std::map<std::pair<int, int>, double> ignoredMassZD;
    if (splitPrincipalBlocks(
            mass,
            dynamicMap,
            masslessMap,
            massDD,
            ignoredMassZZ,
            ignoredMassZD,
            message) < 0)
        return -4;
    if (csrFromEntries(
            static_cast<int>(result.dynamic.size()),
            massDD,
            result.reducedMass,
            message) < 0)
        return -5;

    {
        MumpsSPD massGuard;
        if (massGuard.factorize(result.reducedMass, message) < 0) {
            message = "dynamic mass is not SPD: " + message;
            return -6;
        }
    }

    if (result.massless.empty()) {
        if (csrFromEntries(
                static_cast<int>(result.dynamic.size()),
                stiffnessDD,
                result.reducedStiffness,
                message) < 0)
            return -7;
        return 0;
    }

    SymmetricCSR stiffnessMassless;
    if (csrFromEntries(
            static_cast<int>(result.massless.size()),
            stiffnessZZ,
            stiffnessMassless,
            message) < 0)
        return -8;
    result.impl_->masslessStiffnessFactor.reset(new MumpsSPD());
    if (result.impl_->masslessStiffnessFactor->factorize(stiffnessMassless, message) < 0) {
        message = "massless stiffness block K_ZZ is not SPD: " + message;
        return -9;
    }

    const int numberOfDynamic = static_cast<int>(result.dynamic.size());
    const int numberOfMassless = static_cast<int>(result.massless.size());
    result.impl_->crossRowOffsets.assign(result.massless.size() + 1u, 0);
    for (const auto &entry : stiffnessZD) {
        if (entry.first.first < 0 || entry.first.first >= numberOfMassless ||
            entry.first.second < 0 || entry.first.second >= numberOfDynamic ||
            !std::isfinite(entry.second)) {
            message = "invalid sparse K_ZD entry";
            return -10;
        }
        ++result.impl_->crossRowOffsets[static_cast<std::size_t>(entry.first.first) + 1u];
    }
    for (int row = 0; row < numberOfMassless; ++row)
        result.impl_->crossRowOffsets[static_cast<std::size_t>(row) + 1u] +=
            result.impl_->crossRowOffsets[static_cast<std::size_t>(row)];
    result.impl_->crossColumns.reserve(stiffnessZD.size());
    result.impl_->crossValues.reserve(stiffnessZD.size());
    for (const auto &entry : stiffnessZD) {
        result.impl_->crossColumns.push_back(entry.first.second);
        result.impl_->crossValues.push_back(entry.second);
    }

    std::size_t correctionSize = 0;
    if (!checkedProduct(result.dynamic.size(), result.dynamic.size(), correctionSize)) {
        message = "static Schur correction overflows storage";
        return -11;
    }
    std::vector<double> correction(correctionSize, 0.0);
    constexpr int staticBlockSize = 32;
    for (int firstColumn = 0; firstColumn < numberOfDynamic;
         firstColumn += staticBlockSize) {
        const int blockColumns = std::min(staticBlockSize, numberOfDynamic - firstColumn);
        std::size_t tileSize = 0;
        if (!checkedProduct(
                result.massless.size(), static_cast<std::size_t>(blockColumns), tileSize)) {
            message = "static-response tile overflows storage";
            return -12;
        }
        std::vector<double> response(tileSize, 0.0);
        for (int z = 0; z < numberOfMassless; ++z) {
            for (int position = result.impl_->crossRowOffsets[static_cast<std::size_t>(z)];
                 position < result.impl_->crossRowOffsets[static_cast<std::size_t>(z) + 1u];
                 ++position) {
                const int dynamicColumn =
                    result.impl_->crossColumns[static_cast<std::size_t>(position)];
                if (dynamicColumn >= firstColumn &&
                    dynamicColumn < firstColumn + blockColumns) {
                    response[static_cast<std::size_t>(dynamicColumn - firstColumn) *
                                 numberOfMassless + z] =
                        result.impl_->crossValues[static_cast<std::size_t>(position)];
                }
            }
        }
        if (result.impl_->masslessStiffnessFactor->solve(
                response, blockColumns, message) < 0) {
            message = "K_ZZ blocked multi-RHS solve failed: " + message;
            return -13;
        }
        for (double &value : response)
            value = -value;
        for (int localColumn = 0; localColumn < blockColumns; ++localColumn) {
            const int column = firstColumn + localColumn;
            for (int z = 0; z < numberOfMassless; ++z) {
                const double responseValue =
                    response[static_cast<std::size_t>(localColumn) * numberOfMassless + z];
                for (int position = result.impl_->crossRowOffsets[static_cast<std::size_t>(z)];
                     position < result.impl_->crossRowOffsets[static_cast<std::size_t>(z) + 1u];
                     ++position) {
                    const int row =
                        result.impl_->crossColumns[static_cast<std::size_t>(position)];
                    correction[static_cast<std::size_t>(column) * numberOfDynamic + row] +=
                        result.impl_->crossValues[static_cast<std::size_t>(position)] *
                        responseValue;
                }
            }
        }
    }

    double maximumDifference = 0.0;
    double maximumMagnitude = 0.0;
    int maximumRow = -1;
    int maximumColumn = -1;
    for (int column = 0; column < numberOfDynamic; ++column) {
        for (int row = 0; row <= column; ++row) {
            const double correctionRowColumn =
                correction[static_cast<std::size_t>(column) * numberOfDynamic + row];
            const double correctionColumnRow =
                correction[static_cast<std::size_t>(row) * numberOfDynamic + column];
            const double difference =
                std::fabs(correctionRowColumn - correctionColumnRow);
            maximumMagnitude = std::max(
                maximumMagnitude,
                std::max(std::fabs(correctionRowColumn),
                         std::fabs(correctionColumnRow)));
            if (difference > maximumDifference) {
                maximumDifference = difference;
                maximumRow = row;
                maximumColumn = column;
            }
        }
    }
    const double roundoffRelativeTolerance =
        256.0 * std::numeric_limits<double>::epsilon() *
        std::max(1, std::max(numberOfDynamic, numberOfMassless));
    const double symmetryRelativeTolerance =
        std::max(relativeTolerance, roundoffRelativeTolerance);
    const double symmetryScale = std::max(1.0, maximumMagnitude);
    const double symmetryThreshold =
        absoluteTolerance + symmetryRelativeTolerance * symmetryScale;
    if (maximumDifference > symmetryThreshold) {
        std::ostringstream detail;
        detail.setf(std::ios::scientific);
        detail.precision(6);
        detail << "static Schur correction lost symmetry at (" << maximumRow
               << ',' << maximumColumn << "): absDiff=" << maximumDifference
               << ", normScale=" << symmetryScale
               << ", threshold=" << symmetryThreshold;
        message = detail.str();
        return -14;
    }
    for (int column = 0; column < numberOfDynamic; ++column) {
        for (int row = 0; row <= column; ++row) {
            const double correctionRowColumn =
                correction[static_cast<std::size_t>(column) * numberOfDynamic + row];
            const double correctionColumnRow =
                correction[static_cast<std::size_t>(row) * numberOfDynamic + column];
            stiffnessDD[std::make_pair(row, column)] +=
                0.5 * (correctionRowColumn + correctionColumnRow);
        }
    }
    if (csrFromEntries(
            numberOfDynamic, stiffnessDD, result.reducedStiffness, message) < 0)
        return -15;
    return 0;
}

int reconstructStaticCoordinates(
    const StaticCondensation &condensation,
    const std::vector<double> &reducedVectors,
    int numberOfColumns,
    std::vector<double> &fullVectors,
    std::string &message)
{
    message.clear();
    fullVectors.clear();
    if (condensation.originalDimension < 0 || numberOfColumns < 1 ||
        condensation.dynamic.size() + condensation.massless.size() !=
            static_cast<std::size_t>(condensation.originalDimension)) {
        message = "invalid static-reconstruction dimensions";
        return -1;
    }
    std::size_t reducedSize = 0;
    std::size_t fullSize = 0;
    if (!checkedProduct(condensation.dynamic.size(), static_cast<std::size_t>(numberOfColumns),
                        reducedSize) ||
        !checkedProduct(static_cast<std::size_t>(condensation.originalDimension),
                        static_cast<std::size_t>(numberOfColumns), fullSize) ||
        reducedVectors.size() != reducedSize) {
        message = "reduced vectors do not match static-reconstruction dimensions";
        return -2;
    }
    if (!condensation.impl_ ||
        (!condensation.massless.empty() &&
         !condensation.impl_->masslessStiffnessFactor)) {
        message = "static reconstruction has no K_ZZ factorization";
        return -3;
    }
    fullVectors.assign(fullSize, 0.0);
    for (int column = 0; column < numberOfColumns; ++column) {
        for (std::size_t d = 0; d < condensation.dynamic.size(); ++d) {
            const double value = reducedVectors[
                static_cast<std::size_t>(column) * condensation.dynamic.size() + d];
            fullVectors[static_cast<std::size_t>(column) * condensation.originalDimension +
                        condensation.dynamic[d]] = value;
        }
    }
    constexpr int staticBlockSize = 32;
    for (int firstColumn = 0; firstColumn < numberOfColumns; firstColumn += staticBlockSize) {
        const int blockColumns = std::min(staticBlockSize, numberOfColumns - firstColumn);
        std::size_t tileSize = 0;
        if (!checkedProduct(
                condensation.massless.size(), static_cast<std::size_t>(blockColumns), tileSize)) {
            message = "static-reconstruction tile overflows storage";
            return -4;
        }
        if (tileSize == 0u)
            continue;
        std::vector<double> response(tileSize, 0.0);
        for (std::size_t z = 0; z < condensation.massless.size(); ++z) {
            for (int position = condensation.impl_->crossRowOffsets[z];
                 position < condensation.impl_->crossRowOffsets[z + 1u]; ++position) {
                const int dynamicColumn =
                    condensation.impl_->crossColumns[static_cast<std::size_t>(position)];
                const double stiffnessValue =
                    condensation.impl_->crossValues[static_cast<std::size_t>(position)];
                for (int localColumn = 0; localColumn < blockColumns; ++localColumn)
                    response[static_cast<std::size_t>(localColumn) * condensation.massless.size() +
                             z] -= stiffnessValue *
                        reducedVectors[
                            static_cast<std::size_t>(firstColumn + localColumn) *
                                condensation.dynamic.size() +
                            dynamicColumn];
            }
        }
        if (condensation.impl_->masslessStiffnessFactor->solve(
                response, blockColumns, message) < 0) {
            message = "static-reconstruction K_ZZ solve failed: " + message;
            return -5;
        }
        for (int localColumn = 0; localColumn < blockColumns; ++localColumn) {
            for (std::size_t z = 0; z < condensation.massless.size(); ++z)
                fullVectors[
                    static_cast<std::size_t>(firstColumn + localColumn) *
                        condensation.originalDimension +
                    condensation.massless[z]] =
                    response[static_cast<std::size_t>(localColumn) *
                                 condensation.massless.size() + z];
        }
    }
    if (!std::all_of(
            fullVectors.begin(), fullVectors.end(),
            [](double value) { return std::isfinite(value); })) {
        message = "static reconstruction produced a non-finite value";
        return -6;
    }
    return 0;
}

} // namespace ladruno_cms
