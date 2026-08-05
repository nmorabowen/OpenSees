// ADR-1000 P3b -- the distributed assembly oracle.
//
// The gate the P3 execution plan calls "the single most important P3 gate":
// prove that "physically distributed" means the ALGEBRA is a genuine split, not
// a replicated matrix that happens to give the right answer. Everything above it
// in the hierarchy (T2/S2/T1/S1, refinement) consumes the rank-local pencils
// built here, so if the split double-counts or drops a contribution, every
// downstream gate is measuring the wrong operator.
//
// Method. One structured 2-D grid of bilinear cells, one DOF per node, bottom
// row constrained (global equation -1, so the constrained-DOF path is exercised
// too). The same element contributions are assembled twice through the SAME
// production calls -- makeAssemblyRecord + buildSymmetricCSR:
//   * reference: every element, on one rank, no MPI;
//   * distributed: each rank only its owned column strip.
// Then the two are compared as OPERATORS (K.x and M.x, with the distributed
// action formed exactly as LadrunoCMSSubspaceRefiner::globalAction does it --
// SymmetricCSR::multiply per rank + MPI_Allreduce) and ENTRYWISE (sum of the
// rank-local matrices vs the reference matrix).
//
// The matvec is the measuring instrument, not the subject: the identical kernel
// runs on both sides, so any difference is attributable to the assembly split.
//
// Deliberately size-generic: runs its distributed legs at ANY size >= 2. A check
// that guards its collective leg on `size == N` silently tests nothing at other
// sizes -- that is exactly how Part 0's np=2 acceptance gate came to be vacuous
// (LEDGER_quirks, 2026-07-26).

#include "LadrunoCMSOptions.h"

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

// Evaluates the call FIRST, then builds the diagnostic. As two arguments of
// require() their evaluation order is unspecified in C++, and MSVC was building
// the message string BEFORE the call filled `message` -- every failure printed
// an empty diagnostic. See LEDGER_quirks (2026-07-26).
#define REQUIRE_CALL(status, text)                                                 do {                                                                               const bool ladrunoCallOk = (status);                                           require(ladrunoCallOk, text);                                              } while (0)

// ---------------------------------------------------------------------------
// Model: nx by ny bilinear cells, one DOF per node, bottom row constrained.
// ---------------------------------------------------------------------------

struct GridModel {
    int nx = 6;
    int ny = 5;

    int nodesX() const { return nx + 1; }
    int nodesY() const { return ny + 1; }
    int numberOfElements() const { return nx * ny; }

    // Global equation of node (i,j); -1 where constrained (j == 0).
    int equation(int i, int j) const
    {
        if (j == 0)
            return -1;
        return (j - 1) * nodesX() + i;
    }

    int globalDimension() const { return ny * nodesX(); }

    // Counter-clockwise node equations of cell (ex,ey).
    std::vector<int> elementIDs(int ex, int ey) const
    {
        return {equation(ex, ey), equation(ex + 1, ey),
                equation(ex + 1, ey + 1), equation(ex, ey + 1)};
    }

    // Contiguous column strips, so every internal strip boundary is a real
    // interface shared by two ranks.
    int owner(int ex, int size) const
    {
        return static_cast<int>((static_cast<long long>(ex) * size) / nx);
    }
};

// Bilinear-quad Laplacian stiffness and consistent mass, row major.
const std::vector<double> &elementStiffness()
{
    static const std::vector<double> values = {
         4.0 / 6.0, -1.0 / 6.0, -2.0 / 6.0, -1.0 / 6.0,
        -1.0 / 6.0,  4.0 / 6.0, -1.0 / 6.0, -2.0 / 6.0,
        -2.0 / 6.0, -1.0 / 6.0,  4.0 / 6.0, -1.0 / 6.0,
        -1.0 / 6.0, -2.0 / 6.0, -1.0 / 6.0,  4.0 / 6.0};
    return values;
}

const std::vector<double> &elementMass()
{
    static const std::vector<double> values = {
        4.0 / 36.0, 2.0 / 36.0, 1.0 / 36.0, 2.0 / 36.0,
        2.0 / 36.0, 4.0 / 36.0, 2.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 2.0 / 36.0, 4.0 / 36.0, 2.0 / 36.0,
        2.0 / 36.0, 1.0 / 36.0, 2.0 / 36.0, 4.0 / 36.0};
    return values;
}

// Assemble the elements selected by `keep` into one rank-local pencil, through
// the production calls. `keep(elementIndex) -> bool`.
template <typename Selector>
bool assemblePencil(
    const GridModel &model,
    const Selector &keep,
    ladruno_cms::SymmetricCSR &stiffness,
    ladruno_cms::SymmetricCSR &mass,
    std::vector<std::size_t> &ordinals,
    std::string &message)
{
    std::vector<ladruno_cms::AssemblyRecord> stiffnessRecords;
    std::vector<ladruno_cms::AssemblyRecord> massRecords;
    ordinals.clear();
    for (int ey = 0; ey < model.ny; ++ey) {
        for (int ex = 0; ex < model.nx; ++ex) {
            const int element = ey * model.nx + ex;
            if (!keep(element))
                continue;
            const std::size_t ordinal = static_cast<std::size_t>(element);
            const std::vector<int> ids = model.elementIDs(ex, ey);
            ladruno_cms::AssemblyRecord stiffnessRecord;
            ladruno_cms::AssemblyRecord massRecord;
            if (ladruno_cms::makeAssemblyRecord(
                    ladruno_cms::ContributionKind::Stiffness, ordinal,
                    elementStiffness().data(), 4, 4, ids,
                    model.globalDimension(), 1.0, 1.0e-12, 1.0e-14,
                    stiffnessRecord, message) != 0)
                return false;
            if (ladruno_cms::makeAssemblyRecord(
                    ladruno_cms::ContributionKind::Mass, ordinal,
                    elementMass().data(), 4, 4, ids,
                    model.globalDimension(), 1.0, 1.0e-12, 1.0e-14,
                    massRecord, message) != 0)
                return false;
            stiffnessRecords.push_back(stiffnessRecord);
            massRecords.push_back(massRecord);
            ordinals.push_back(ordinal);
        }
    }
    if (ladruno_cms::buildSymmetricCSR(
            model.globalDimension(), stiffnessRecords,
            ladruno_cms::ContributionKind::Stiffness, stiffness, message) != 0)
        return false;
    return ladruno_cms::buildSymmetricCSR(
               model.globalDimension(), massRecords,
               ladruno_cms::ContributionKind::Mass, mass, message) == 0;
}

// The production distributed action: SymmetricCSR::multiply per rank, summed
// across MPI_COMM_WORLD (LadrunoCMSSubspaceRefiner's globalAction, inlined
// because that helper is file-static).
std::vector<double> distributedAction(
    const ladruno_cms::SymmetricCSR &localMatrix,
    const std::vector<double> &vector)
{
    const std::vector<double> local = localMatrix.multiply(vector);
    std::vector<double> global(local.size(), 0.0);
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

double worstRelativeDifference(
    const std::vector<double> &left, const std::vector<double> &right)
{
    double scale = 0.0;
    for (double value : right)
        scale = std::max(scale, std::fabs(value));
    if (scale == 0.0)
        scale = 1.0;
    double worst = 0.0;
    for (std::size_t entry = 0; entry < left.size(); ++entry)
        worst = std::max(worst, std::fabs(left[entry] - right[entry]) / scale);
    return worst;
}

// Deterministic pseudo-random probes: a fixed-seed LCG, so the vectors are
// identical on every rank and every platform without <random>'s freedom.
std::vector<double> probeVector(const GridModel &model, int which)
{
    const int dimension = model.globalDimension();
    std::vector<double> vector(static_cast<std::size_t>(dimension), 0.0);
    if (which == 0) {                       // e1
        vector[0] = 1.0;
    } else if (which == 1) {                // ones
        std::fill(vector.begin(), vector.end(), 1.0);
    } else if (which == 2) {                // ramp
        for (int row = 0; row < dimension; ++row)
            vector[static_cast<std::size_t>(row)] =
                static_cast<double>(row + 1);
    } else {                                // fixed-seed pseudo-random
        unsigned long long state =
            1469598103934665603ull + static_cast<unsigned long long>(which);
        for (int row = 0; row < dimension; ++row) {
            state = state * 6364136223846793005ull + 1442695040888963407ull;
            const double unit =
                static_cast<double>((state >> 11) & 0xFFFFFFFFull) /
                static_cast<double>(0xFFFFFFFFull);
            vector[static_cast<std::size_t>(row)] = 2.0 * unit - 1.0;
        }
    }
    return vector;
}

// ---------------------------------------------------------------------------

void checkDistributedAssembly(int rank, int size)
{
    if (size < 2)
        return;
    const GridModel model;
    const int dimension = model.globalDimension();
    std::string message;

    ladruno_cms::SymmetricCSR referenceStiffness;
    ladruno_cms::SymmetricCSR referenceMass;
    std::vector<std::size_t> referenceOrdinals;
    REQUIRE_CALL(assemblePencil(
                model, [](int) { return true; }, referenceStiffness,
                referenceMass, referenceOrdinals, message),
            "reference assembly failed: " + message);
    require(referenceOrdinals.size() ==
                static_cast<std::size_t>(model.numberOfElements()),
            "reference assembly did not cover every element");

    ladruno_cms::SymmetricCSR localStiffness;
    ladruno_cms::SymmetricCSR localMass;
    std::vector<std::size_t> ownedOrdinals;
    const auto owned = [&model, rank, size](int element) {
        return model.owner(element % model.nx, size) == rank;
    };
    REQUIRE_CALL(assemblePencil(model, owned, localStiffness, localMass,
                           ownedOrdinals, message),
            "distributed assembly failed: " + message);

    // Exact element coverage: every element owned by exactly one rank.
    long long ownedHere = static_cast<long long>(ownedOrdinals.size());
    long long ownedEverywhere = 0;
    MPI_Allreduce(&ownedHere, &ownedEverywhere, 1, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    require(ownedEverywhere == model.numberOfElements(),
            "owned elements summed to " + std::to_string(ownedEverywhere) +
                " instead of " + std::to_string(model.numberOfElements()) +
                " -- the partition drops or duplicates elements");
    require(ownedHere > 0, "a rank owns no elements; the partition is degenerate");

    // Contribution locality: the rank-local pencil must NOT span every equation.
    // (A replicated matrix would touch all of [0,n) on every rank.)
    std::vector<char> touched(static_cast<std::size_t>(dimension), 0);
    for (int row = 0; row < dimension; ++row) {
        for (int entry = localStiffness.rowOffsets[static_cast<std::size_t>(row)];
             entry < localStiffness.rowOffsets[static_cast<std::size_t>(row + 1)];
             ++entry) {
            touched[static_cast<std::size_t>(row)] = 1;
            touched[static_cast<std::size_t>(
                localStiffness.columnIndices[static_cast<std::size_t>(entry)])] = 1;
        }
    }
    const long long touchedHere =
        std::count(touched.begin(), touched.end(), static_cast<char>(1));
    require(touchedHere < dimension,
            "the rank-local pencil touches every one of the " +
                std::to_string(dimension) +
                " equations -- this looks replicated, not split");

    // Every equation must be touched by at least one rank, and the interface
    // (equations touched by >= 2 ranks) must be non-empty, or the partition is
    // not exercising shared DOFs at all.
    std::vector<int> localTouch(static_cast<std::size_t>(dimension), 0);
    for (int row = 0; row < dimension; ++row)
        localTouch[static_cast<std::size_t>(row)] =
            touched[static_cast<std::size_t>(row)] ? 1 : 0;
    std::vector<int> touchCount(static_cast<std::size_t>(dimension), 0);
    MPI_Allreduce(localTouch.data(), touchCount.data(), dimension, MPI_INT,
                  MPI_SUM, MPI_COMM_WORLD);
    const int unreached =
        static_cast<int>(std::count(touchCount.begin(), touchCount.end(), 0));
    require(unreached == 0,
            std::to_string(unreached) +
                " equations are owned by no rank -- contributions were dropped");
    const int shared = static_cast<int>(std::count_if(
        touchCount.begin(), touchCount.end(), [](int count) { return count > 1; }));
    require(shared > 0,
            "no equation is shared between ranks -- the partition has no "
            "interface, so this fixture cannot detect a double count");

    // ---- the Kx / Mx oracle -------------------------------------------------
    for (int which = 0; which < 6; ++which) {
        const std::vector<double> probe = probeVector(model, which);
        const std::vector<double> referenceStiffnessAction =
            referenceStiffness.multiply(probe);
        const std::vector<double> referenceMassAction =
            referenceMass.multiply(probe);
        const std::vector<double> distributedStiffnessAction =
            distributedAction(localStiffness, probe);
        const std::vector<double> distributedMassAction =
            distributedAction(localMass, probe);
        const double stiffnessDifference = worstRelativeDifference(
            distributedStiffnessAction, referenceStiffnessAction);
        const double massDifference =
            worstRelativeDifference(distributedMassAction, referenceMassAction);
        require(stiffnessDifference <= 1.0e-12,
                "distributed K.x differs from the serial K.x on probe " +
                    std::to_string(which) + " by " +
                    std::to_string(stiffnessDifference));
        require(massDifference <= 1.0e-12,
                "distributed M.x differs from the serial M.x on probe " +
                    std::to_string(which) + " by " +
                    std::to_string(massDifference));
    }

    // ---- entrywise: sum of the rank-local matrices == the reference ---------
    const std::vector<double> localDense = localStiffness.toDense();
    const std::vector<double> referenceDense = referenceStiffness.toDense();
    std::vector<double> summedDense(localDense.size(), 0.0);
    MPI_Allreduce(localDense.data(), summedDense.data(),
                  static_cast<int>(localDense.size()), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    require(worstRelativeDifference(summedDense, referenceDense) <= 1.0e-12,
            "the rank-local stiffness matrices do not sum to the serial matrix");

    // ---- negative control 1: replicated input (the double-count trap) -------
    // Every rank assembles EVERY element. The oracle must reject this; if it
    // does not, the checks above prove nothing.
    ladruno_cms::SymmetricCSR replicatedStiffness;
    ladruno_cms::SymmetricCSR replicatedMass;
    std::vector<std::size_t> replicatedOrdinals;
    REQUIRE_CALL(assemblePencil(model, [](int) { return true; },
                           replicatedStiffness, replicatedMass,
                           replicatedOrdinals, message),
            "replicated assembly fixture failed: " + message);
    const std::vector<double> probe = probeVector(model, 3);
    const double replicatedDifference = worstRelativeDifference(
        distributedAction(replicatedStiffness, probe),
        referenceStiffness.multiply(probe));
    require(replicatedDifference > 1.0e-9,
            "NEGATIVE CONTROL FAILED: a replicated (double-counted) pencil was "
            "not distinguished from the correct split -- the oracle is blind");

    // ---- negative control 2: one element owned twice ------------------------
    // Rank 1 additionally claims element 0, which rank 0 already owns.
    const auto corrupted = [&owned, rank](int element) {
        return owned(element) || (rank == 1 && element == 0);
    };
    ladruno_cms::SymmetricCSR corruptedStiffness;
    ladruno_cms::SymmetricCSR corruptedMass;
    std::vector<std::size_t> corruptedOrdinals;
    REQUIRE_CALL(assemblePencil(model, corrupted, corruptedStiffness,
                           corruptedMass, corruptedOrdinals, message),
            "corrupted-owner fixture failed: " + message);
    const double corruptedDifference = worstRelativeDifference(
        distributedAction(corruptedStiffness, probe),
        referenceStiffness.multiply(probe));
    require(corruptedDifference > 1.0e-9,
            "NEGATIVE CONTROL FAILED: one element counted twice was not "
            "distinguished from the correct split");
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2 && rank == 0)
        std::cout << "distributed assembly oracle SKIPPED: needs >= 2 ranks, "
                     "got " << size << '\n';
    checkDistributedAssembly(rank, size);
    MPI_Finalize();
    if (failures != 0)
        return 1;
    std::cout << "Ladruno CMS distributed assembly oracle passed\n";
    return 0;
}
