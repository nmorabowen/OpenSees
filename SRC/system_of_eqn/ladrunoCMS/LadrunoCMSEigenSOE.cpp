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

#include "LadrunoCMSEigenSOE.h"

#include "LadrunoCMSEigenSolver.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <ID.h>
#include <Matrix.h>
#include <MetisWrapper.h>
#include <OPS_Globals.h>
#include <Vertex.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <climits>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

namespace {

int collectiveFailure(int localFailure, MPI_Comm communicator)
{
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, communicator);
    return globalFailure;
}

void appendStructure(
    const std::vector<ladruno_cms::AssemblyRecord> &records,
    std::vector<long long> &encoded)
{
    encoded.push_back(static_cast<long long>(records.size()));
    for (const auto &record : records) {
        encoded.push_back(static_cast<long long>(record.ordinal));
        encoded.push_back(static_cast<long long>(record.kind));
        encoded.push_back(record.rows);
        encoded.push_back(record.columns);
        encoded.push_back(static_cast<long long>(record.originalIDs.size()));
        for (int equation : record.originalIDs)
            encoded.push_back(equation);
        encoded.push_back(static_cast<long long>(record.activeIDs.size()));
        for (int equation : record.activeIDs)
            encoded.push_back(equation);
        encoded.push_back(static_cast<long long>(record.upperValues.size()));
        encoded.push_back(static_cast<long long>(record.metrics.numberOfValues));
    }
}

void appendMetrics(
    const std::vector<ladruno_cms::AssemblyRecord> &records,
    std::vector<double> &metrics)
{
    for (const auto &record : records) {
        metrics.push_back(record.metrics.maxAbs);
        metrics.push_back(record.metrics.sum);
        metrics.push_back(record.metrics.sumAbs);
        metrics.push_back(record.metrics.sumSquares);
        metrics.push_back(record.metrics.weightedProjection);
    }
}

void appendValues(
    const std::vector<ladruno_cms::AssemblyRecord> &records,
    std::vector<double> &values)
{
    for (const auto &record : records)
        values.insert(values.end(), record.upperValues.begin(), record.upperValues.end());
}

bool closeCollectiveValue(
    double local,
    double reference,
    double relativeTolerance,
    double absoluteTolerance)
{
    return std::fabs(local - reference) <= absoluteTolerance +
        relativeTolerance * std::max(1.0, std::max(std::fabs(local), std::fabs(reference)));
}

} // namespace

LadrunoCMSEigenSOE::LadrunoCMSEigenSOE(
    LadrunoCMSEigenSolver &solver,
    const ladruno_cms::Options &options)
    : EigenSOE(solver, EigenSOE_TAGS_LadrunoCMS), options_(options)
{
    solver.setEigenSOE(*this);
}

int LadrunoCMSEigenSOE::getNumEqn() const
{
    return size_;
}

int LadrunoCMSEigenSOE::setSize(Graph &graph)
{
    size_ = graph.getNumVertex();
    stiffness_.clear();
    mass_.clear();
    stiffnessSignatures_.clear();
    massSignatures_.clear();
    nextStiffnessOrdinal_ = 0;
    nextMassOrdinal_ = 0;
    stiffnessVerificationBytes_ = 0;
    massVerificationBytes_ = 0;
    std::string message;
    if (buildTopology(graph, message) < 0) {
        opserr << "LadrunoCMSEigenSOE::setSize - " << message.c_str() << endln;
        return -1;
    }
    return getSolver()->setSize();
}

int LadrunoCMSEigenSOE::addA(const Matrix &matrix, const ID &id, double factor)
{
    return addContribution(ladruno_cms::ContributionKind::Stiffness, matrix, id, factor);
}

int LadrunoCMSEigenSOE::addM(const Matrix &matrix, const ID &id, double factor)
{
    return addContribution(ladruno_cms::ContributionKind::Mass, matrix, id, factor);
}

void LadrunoCMSEigenSOE::zeroA()
{
    stiffness_.clear();
    stiffnessSignatures_.clear();
    nextStiffnessOrdinal_ = 0;
    stiffnessVerificationBytes_ = 0;
}

void LadrunoCMSEigenSOE::zeroM()
{
    mass_.clear();
    massSignatures_.clear();
    nextMassOrdinal_ = 0;
    massVerificationBytes_ = 0;
}

int LadrunoCMSEigenSOE::sendSelf(int, Channel &)
{
    opserr << "LadrunoCMSEigenSOE::sendSelf - unsupported in ADR 1000 v1\n";
    return -1;
}

int LadrunoCMSEigenSOE::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    opserr << "LadrunoCMSEigenSOE::recvSelf - unsupported in ADR 1000 v1\n";
    return -1;
}

const ladruno_cms::Options &LadrunoCMSEigenSOE::getOptions() const
{
    return options_;
}

const std::vector<ladruno_cms::AssemblyRecord> &
LadrunoCMSEigenSOE::getStiffnessContributions() const
{
    return stiffness_;
}

const std::vector<ladruno_cms::AssemblyRecord> &
LadrunoCMSEigenSOE::getMassContributions() const
{
    return mass_;
}

const std::vector<int> &LadrunoCMSEigenSOE::getFineLabels() const
{
    return fineLabels_;
}

const std::vector<int> &LadrunoCMSEigenSOE::getCoarseOfFine() const
{
    return coarseOfFine_;
}

int LadrunoCMSEigenSOE::getWorldRank() const
{
    return worldRank_;
}

int LadrunoCMSEigenSOE::getWorldSize() const
{
    return worldSize_;
}

int LadrunoCMSEigenSOE::getNumberOfCoarseSubdomains() const
{
    return numberOfCoarseSubdomains_;
}

int LadrunoCMSEigenSOE::getFineSubdomainsPerCoarse() const
{
    return fineSubdomainsPerCoarse_;
}

int LadrunoCMSEigenSOE::buildTopology(Graph &graph, std::string &message)
{
    message.clear();
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized)
        MPI_Finalized(&finalized);
    if (!initialized || finalized) {
        message = "an active MPI runtime is required";
        return -1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank_);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);
    if (size_ < 1 || worldSize_ < 1) {
        message = "the equation graph and MPI world must be nonempty";
        return -2;
    }

    MPI_Comm groupComm = MPI_COMM_NULL;
    int group = -1;
    int groupRank = -1;
    int groupSize = 0;
    if (options_.hierarchy == ladruno_cms::HierarchyMode::Logical) {
        if (options_.level1 < 1 || options_.level2 < 1 ||
            static_cast<long long>(options_.level1) * options_.level2 != worldSize_) {
            message = "logical hierarchy requires level1*level2 == MPI world size";
            return -3;
        }
        numberOfCoarseSubdomains_ = options_.level1;
        fineSubdomainsPerCoarse_ = options_.level2;
        group = worldRank_ / fineSubdomainsPerCoarse_;
        MPI_Comm_split(MPI_COMM_WORLD, group, worldRank_, &groupComm);
    } else {
        MPI_Comm_split_type(
            MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, worldRank_, MPI_INFO_NULL, &groupComm);
    }
    MPI_Comm_rank(groupComm, &groupRank);
    MPI_Comm_size(groupComm, &groupSize);

    if (options_.hierarchy == ladruno_cms::HierarchyMode::Auto) {
        int minimumGroupSize = 0;
        int maximumGroupSize = 0;
        MPI_Allreduce(&groupSize, &minimumGroupSize, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&groupSize, &maximumGroupSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        int leader = groupRank == 0 ? 1 : 0;
        MPI_Allreduce(&leader, &numberOfCoarseSubdomains_, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        fineSubdomainsPerCoarse_ = groupSize;
        int localFailure = minimumGroupSize != maximumGroupSize ? 1 : 0;
        if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
            MPI_Comm_free(&groupComm);
            message = "auto hierarchy requires the same rank count on every shared-memory node";
            return -4;
        }
        MPI_Comm leaderComm = MPI_COMM_NULL;
        MPI_Comm_split(
            MPI_COMM_WORLD, groupRank == 0 ? 0 : MPI_UNDEFINED,
            worldRank_, &leaderComm);
        if (groupRank == 0)
            MPI_Comm_rank(leaderComm, &group);
        MPI_Bcast(&group, 1, MPI_INT, 0, groupComm);
        if (leaderComm != MPI_COMM_NULL)
            MPI_Comm_free(&leaderComm);
    }

    coarseOfFine_.assign(static_cast<std::size_t>(worldSize_), -1);
    MPI_Allgather(&group, 1, MPI_INT, coarseOfFine_.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> groupWorldRanks(static_cast<std::size_t>(groupSize), -1);
    MPI_Allgather(&worldRank_, 1, MPI_INT,
                  groupWorldRanks.data(), 1, MPI_INT, groupComm);

    std::vector<int> coarseLabels(static_cast<std::size_t>(size_), -1);
    int localFailure = 0;
    if (worldRank_ == 0) {
        if (numberOfCoarseSubdomains_ == 1) {
            std::fill(coarseLabels.begin(), coarseLabels.end(), 0);
        } else if (size_ < numberOfCoarseSubdomains_) {
            localFailure = 1;
        } else {
            Metis partitioner(numberOfCoarseSubdomains_);
            if (partitioner.partition(graph, numberOfCoarseSubdomains_) < 0) {
                localFailure = 1;
            } else {
                for (int equation = 0; equation < size_; ++equation) {
                    Vertex *vertex = graph.getVertexPtr(equation);
                    if (vertex == nullptr || vertex->getColor() < 1 ||
                        vertex->getColor() > numberOfCoarseSubdomains_) {
                        localFailure = 1;
                        break;
                    }
                    coarseLabels[static_cast<std::size_t>(equation)] =
                        vertex->getColor() - 1;
                }
            }
        }
    }
    if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
        MPI_Comm_free(&groupComm);
        message = "coarse METIS partition failed or produced invalid labels";
        return -5;
    }
    MPI_Bcast(coarseLabels.data(), size_, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> proposedFineLabels(static_cast<std::size_t>(size_), -1);
    localFailure = 0;
    if (groupRank == 0) {
        std::vector<int> equations;
        equations.reserve(static_cast<std::size_t>(size_));
        std::vector<int> localOfGlobal(static_cast<std::size_t>(size_), -1);
        for (int equation = 0; equation < size_; ++equation) {
            if (coarseLabels[static_cast<std::size_t>(equation)] == group) {
                localOfGlobal[static_cast<std::size_t>(equation)] =
                    static_cast<int>(equations.size());
                equations.push_back(equation);
            }
        }
        if (equations.size() < static_cast<std::size_t>(groupSize)) {
            localFailure = 1;
        } else if (groupSize == 1) {
            for (int equation : equations)
                proposedFineLabels[static_cast<std::size_t>(equation)] = groupWorldRanks[0];
        } else {
            Graph induced(static_cast<int>(equations.size()));
            for (std::size_t local = 0; local < equations.size(); ++local) {
                if (!induced.addVertex(new Vertex(
                        static_cast<int>(local), equations[local]), false)) {
                    localFailure = 1;
                    break;
                }
            }
            if (!localFailure) {
                for (std::size_t local = 0; local < equations.size(); ++local) {
                    Vertex *source = graph.getVertexPtr(equations[local]);
                    if (source == nullptr) {
                        localFailure = 1;
                        break;
                    }
                    const ID &adjacency = source->getAdjacency();
                    for (int edge = 0; edge < adjacency.Size(); ++edge) {
                        const int adjacent = adjacency(edge);
                        if (adjacent < 0 || adjacent >= size_)
                            continue;
                        const int other = localOfGlobal[static_cast<std::size_t>(adjacent)];
                        if (other > static_cast<int>(local))
                            induced.addEdge(static_cast<int>(local), other);
                    }
                }
            }
            Metis partitioner(groupSize);
            if (!localFailure && partitioner.partition(induced, groupSize) < 0)
                localFailure = 1;
            if (!localFailure) {
                std::vector<int> partCounts(static_cast<std::size_t>(groupSize), 0);
                for (std::size_t local = 0; local < equations.size(); ++local) {
                    Vertex *vertex = induced.getVertexPtr(static_cast<int>(local));
                    const int localPart = vertex == nullptr ? -1 : vertex->getColor() - 1;
                    if (localPart < 0 || localPart >= groupSize) {
                        localFailure = 1;
                        break;
                    }
                    ++partCounts[static_cast<std::size_t>(localPart)];
                    proposedFineLabels[static_cast<std::size_t>(equations[local])] =
                        groupWorldRanks[static_cast<std::size_t>(localPart)];
                }
                if (std::find(partCounts.begin(), partCounts.end(), 0) != partCounts.end())
                    localFailure = 1;
            }
        }
    }
    if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
        MPI_Comm_free(&groupComm);
        message = "fine METIS partition failed or produced an empty subdomain";
        return -6;
    }
    fineLabels_.assign(static_cast<std::size_t>(size_), -1);
    MPI_Allreduce(proposedFineLabels.data(), fineLabels_.data(), size_, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    if (std::find(fineLabels_.begin(), fineLabels_.end(), -1) != fineLabels_.end()) {
        MPI_Comm_free(&groupComm);
        message = "the hierarchical partition left an equation without a fine owner";
        return -7;
    }
    MPI_Comm_free(&groupComm);
    return 0;
}

int LadrunoCMSEigenSOE::verifyReplicatedAssembly(std::string &message) const
{
    message.clear();
    if (options_.verifyAssembly == ladruno_cms::AssemblyVerification::Off)
        return 0;
    std::vector<long long> localStructure;
    appendStructure(stiffnessSignatures_, localStructure);
    appendStructure(massSignatures_, localStructure);
    long long referenceStructureSize = worldRank_ == 0
        ? static_cast<long long>(localStructure.size())
        : 0;
    MPI_Bcast(&referenceStructureSize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    int localFailure = referenceStructureSize < 0 || referenceStructureSize > INT_MAX;
    std::vector<long long> referenceStructure;
    if (!localFailure)
        referenceStructure.resize(static_cast<std::size_t>(referenceStructureSize));
    if (worldRank_ == 0 && !localFailure)
        referenceStructure = localStructure;
    if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
        message = "assembly structure exceeds MPI count range";
        return -1;
    }
    MPI_Bcast(referenceStructure.data(), static_cast<int>(referenceStructureSize),
              MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    localFailure = localStructure != referenceStructure;

    std::vector<double> localMetrics;
    appendMetrics(stiffnessSignatures_, localMetrics);
    appendMetrics(massSignatures_, localMetrics);
    long long referenceMetricSize = worldRank_ == 0
        ? static_cast<long long>(localMetrics.size())
        : 0;
    MPI_Bcast(&referenceMetricSize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    if (referenceMetricSize < 0 || referenceMetricSize > INT_MAX)
        localFailure = 1;
    std::vector<double> referenceMetrics;
    if (referenceMetricSize >= 0 && referenceMetricSize <= INT_MAX)
        referenceMetrics.resize(static_cast<std::size_t>(referenceMetricSize));
    if (worldRank_ == 0 && referenceMetrics.size() == localMetrics.size())
        referenceMetrics = localMetrics;
    if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
        message = "assembly metrics exceed MPI count range or structure differs";
        return -2;
    }
    MPI_Bcast(referenceMetrics.data(), static_cast<int>(referenceMetricSize),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (localMetrics.size() != referenceMetrics.size()) {
        localFailure = 1;
    } else {
        for (std::size_t index = 0; index < localMetrics.size(); ++index) {
            const std::size_t recordIndex = index / 5u;
            const bool stiffnessRecord = recordIndex < stiffnessSignatures_.size();
            const std::size_t localRecord = stiffnessRecord
                ? recordIndex
                : recordIndex - stiffnessSignatures_.size();
            const auto &record = stiffnessRecord
                ? stiffnessSignatures_[localRecord]
                : massSignatures_[localRecord];
            const double scaledAtol = options_.assemblyAtol * static_cast<double>(
                std::max<std::size_t>(1u, record.metrics.numberOfValues));
            if (!closeCollectiveValue(
                    localMetrics[index], referenceMetrics[index],
                    options_.assemblyRtol, scaledAtol)) {
                localFailure = 1;
                break;
            }
        }
    }

    if (options_.verifyAssembly == ladruno_cms::AssemblyVerification::Full) {
        std::vector<double> localValues;
        appendValues(stiffnessSignatures_, localValues);
        appendValues(massSignatures_, localValues);
        long long referenceValueSize = worldRank_ == 0
            ? static_cast<long long>(localValues.size())
            : 0;
        MPI_Bcast(&referenceValueSize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (referenceValueSize < 0 || referenceValueSize > INT_MAX)
            localFailure = 1;
        std::vector<double> referenceValues;
        if (referenceValueSize >= 0 && referenceValueSize <= INT_MAX)
            referenceValues.resize(static_cast<std::size_t>(referenceValueSize));
        if (worldRank_ == 0 && referenceValues.size() == localValues.size())
            referenceValues = localValues;
        if (collectiveFailure(localFailure, MPI_COMM_WORLD)) {
            message = "full assembly values exceed MPI count range or signatures differ";
            return -3;
        }
        MPI_Bcast(referenceValues.data(), static_cast<int>(referenceValueSize),
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (localValues.size() != referenceValues.size()) {
            localFailure = 1;
        } else {
            for (std::size_t index = 0; index < localValues.size(); ++index) {
                if (!closeCollectiveValue(
                        localValues[index], referenceValues[index],
                        options_.assemblyRtol, options_.assemblyAtol)) {
                    localFailure = 1;
                    break;
                }
            }
        }
    }
    int failingRank = localFailure ? worldRank_ : worldSize_;
    int firstFailingRank = worldSize_;
    MPI_Allreduce(&failingRank, &firstFailingRank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (firstFailingRank != worldSize_) {
        message = "replicated assembly differs first on MPI rank " +
            std::to_string(firstFailingRank);
        return -4;
    }
    return 0;
}

int LadrunoCMSEigenSOE::addContribution(
    ladruno_cms::ContributionKind kind,
    const Matrix &matrix,
    const ID &id,
    double factor)
{
    if (factor == 0.0)
        return 0;
    const int order = id.Size();
    if (order < 0 || matrix.noRows() != order || matrix.noCols() != order) {
        opserr << "LadrunoCMSEigenSOE::addContribution - Matrix and ID sizes differ\n";
        return -1;
    }
    const std::size_t unsignedOrder = static_cast<std::size_t>(order);
    if (unsignedOrder != 0u &&
        unsignedOrder > std::numeric_limits<std::size_t>::max() / unsignedOrder) {
        opserr << "LadrunoCMSEigenSOE::addContribution - Matrix storage size overflows\n";
        return -2;
    }
    std::vector<int> originalIDs;
    originalIDs.reserve(unsignedOrder);
    for (int local = 0; local < order; ++local) {
        originalIDs.push_back(id(local));
    }

    std::vector<double> rowMajor(
        unsignedOrder * unsignedOrder);
    for (int row = 0; row < order; ++row) {
        for (int column = 0; column < order; ++column)
            rowMajor[static_cast<std::size_t>(row) * order + column] = matrix(row, column);
    }
    ladruno_cms::AssemblyRecord contribution;
    std::string message;
    std::size_t &nextOrdinal = kind == ladruno_cms::ContributionKind::Stiffness
        ? nextStiffnessOrdinal_
        : nextMassOrdinal_;
    const int status = ladruno_cms::makeAssemblyRecord(
        kind,
        nextOrdinal,
        rowMajor.data(),
        matrix.noRows(),
        matrix.noCols(),
        originalIDs,
        size_,
        factor,
        options_.assemblyRtol,
        options_.assemblyAtol,
        contribution,
        message);
    if (status < 0) {
        opserr << "LadrunoCMSEigenSOE::addContribution - " << message.c_str() << '\n';
        return status;
    }
    ++nextOrdinal;
    if (contribution.activeIDs.empty())
        return 0;

    const int owner = ladruno_cms::assignContributionOwner(
        contribution.activeIDs, fineLabels_, message);
    if (owner < 0) {
        opserr << "LadrunoCMSEigenSOE::addContribution - " << message.c_str() << '\n';
        return -3;
    }

    if (options_.verifyAssembly != ladruno_cms::AssemblyVerification::Off) {
        ladruno_cms::AssemblyRecord signature = contribution;
        std::size_t &verificationBytes =
            kind == ladruno_cms::ContributionKind::Stiffness
                ? stiffnessVerificationBytes_
                : massVerificationBytes_;
        if (options_.verifyAssembly == ladruno_cms::AssemblyVerification::Full) {
            const std::size_t valueBytes = signature.upperValues.size() * sizeof(double);
            const std::size_t idBytes =
                (signature.originalIDs.size() + signature.activeIDs.size()) * sizeof(int);
            const std::size_t maximum = std::numeric_limits<std::size_t>::max();
            const bool recordOverflow = valueBytes > maximum - idBytes;
            const std::size_t recordBytes = recordOverflow ? maximum : valueBytes + idBytes;
            const bool currentOverflow =
                stiffnessVerificationBytes_ > maximum - massVerificationBytes_;
            const std::size_t currentBytes = currentOverflow
                ? maximum
                : stiffnessVerificationBytes_ + massVerificationBytes_;
            if (recordOverflow || currentOverflow ||
                currentBytes > maximum - recordBytes ||
                currentBytes + recordBytes > options_.verifyFullMaxBytes) {
                opserr << "LadrunoCMSEigenSOE::addContribution - full assembly "
                       << "verification exceeds verifyFullMaxBytes\n";
                return -4;
            }
            verificationBytes += valueBytes + idBytes;
        } else {
            signature.upperValues.clear();
            signature.upperValues.shrink_to_fit();
        }
        if (kind == ladruno_cms::ContributionKind::Stiffness)
            stiffnessSignatures_.push_back(std::move(signature));
        else
            massSignatures_.push_back(std::move(signature));
    }

    if (owner == worldRank_) {
        if (kind == ladruno_cms::ContributionKind::Stiffness)
            stiffness_.push_back(std::move(contribution));
        else
            mass_.push_back(std::move(contribution));
    }
    return 0;
}
