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

#include "LadrunoCMSOptions.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace ladruno_cms {
namespace {

bool checkedAdd(std::size_t left, std::size_t right, std::size_t &value)
{
    if (right > std::numeric_limits<std::size_t>::max() - left)
        return false;
    value = left + right;
    return true;
}

bool checkedMultiply(std::size_t left, std::size_t right, std::size_t &value)
{
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
        return false;
    value = left * right;
    return true;
}

bool parseInteger(const std::string &text, int &value)
{
    errno = 0;
    char *end = nullptr;
    const long parsed = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' ||
        parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max())
        return false;
    value = static_cast<int>(parsed);
    return true;
}

bool parseSize(const std::string &text, std::size_t &value)
{
    errno = 0;
    char *end = nullptr;
    const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' || text.front() == '-' ||
        parsed > std::numeric_limits<std::size_t>::max())
        return false;
    value = static_cast<std::size_t>(parsed);
    return true;
}

bool parseReal(const std::string &text, double &value)
{
    errno = 0;
    char *end = nullptr;
    value = std::strtod(text.c_str(), &end);
    return errno == 0 && end != text.c_str() && *end == '\0' && std::isfinite(value);
}

} // namespace

int parseCommandOptions(
    const std::vector<std::string> &arguments,
    Options &options,
    int &numberOfModes,
    std::string &message)
{
    options = Options{};
    numberOfModes = 0;
    message.clear();
    if (arguments.empty() || !parseInteger(arguments.back(), numberOfModes) ||
        numberOfModes < 1) {
        message = "the last argument must be a positive numModes integer";
        return -1;
    }
    std::set<std::string> seen;
    bool hierarchySpecified = false;
    bool verifyAssemblySpecified = false;
    bool modesLevel2Specified = false;
    bool modesLevel1Specified = false;
    bool level1Specified = false;
    bool level2Specified = false;
    std::size_t position = 0;
    const std::size_t optionEnd = arguments.size() - 1u;
    const auto takeValue = [&](const std::string &flag, std::string &value) -> bool {
        if (position >= optionEnd) {
            message = flag + " requires a value";
            return false;
        }
        value = arguments[position++];
        return true;
    };
    while (position < optionEnd) {
        const std::string flag = arguments[position++];
        if (flag.empty() || flag.front() != '-') {
            message = "unexpected positional argument before numModes: " + flag;
            return -2;
        }
        if (!seen.insert(flag).second) {
            message = "duplicate option: " + flag;
            return -3;
        }
        if (flag == "-verbose") {
            options.verbose = true;
            continue;
        }
        std::string value;
        if (!takeValue(flag, value))
            return -4;
        if (flag == "-domainMode") {
            if (value == "replicatedReference")
                options.domainMode = DomainMode::ReplicatedReference;
            else if (value == "physical")
                options.domainMode = DomainMode::Physical;
            else {
                message = "-domainMode accepts only physical or replicatedReference";
                return -5;
            }
        } else if (flag == "-hierarchy") {
            hierarchySpecified = true;
            if (value == "logical")
                options.hierarchy = HierarchyMode::Logical;
            else if (value == "auto")
                options.hierarchy = HierarchyMode::Auto;
            else {
                message = "-hierarchy accepts only auto or logical";
                return -5;
            }
        } else if (flag == "-level1") {
            level1Specified = parseInteger(value, options.level1);
            if (!level1Specified) {
                message = "-level1 requires an integer";
                return -6;
            }
        } else if (flag == "-level2") {
            level2Specified = parseInteger(value, options.level2);
            if (!level2Specified) {
                message = "-level2 requires an integer";
                return -7;
            }
        } else if (flag == "-modesL2") {
            modesLevel2Specified = parseInteger(value, options.modesLevel2);
            if (!modesLevel2Specified) {
                message = "-modesL2 requires an integer";
                return -8;
            }
        } else if (flag == "-modesL1") {
            modesLevel1Specified = parseInteger(value, options.modesLevel1);
            if (!modesLevel1Specified) {
                message = "-modesL1 requires an integer";
                return -9;
            }
        } else if (flag == "-tol") {
            if (!parseReal(value, options.tolerance)) {
                message = "-tol requires a finite number";
                return -10;
            }
        } else if (flag == "-maxEnrich") {
            if (!parseInteger(value, options.maxEnrich)) {
                message = "-maxEnrich requires an integer";
                return -11;
            }
        } else if (flag == "-maxIter") {
            if (!parseInteger(value, options.maxIterations)) {
                message = "-maxIter requires an integer";
                return -12;
            }
        } else if (flag == "-refinement") {
            if (value == "subspace")
                options.refinement = RefinementMode::Subspace;
            else if (value == "none")
                options.refinement = RefinementMode::None;
            else {
                message = "-refinement accepts only subspace or none";
                return -13;
            }
        } else if (flag == "-iterationVectors") {
            if (value == "auto") {
                options.iterationVectors = 0;
            } else if (!parseInteger(value, options.iterationVectors)) {
                message = "-iterationVectors requires auto or an integer";
                return -14;
            }
        } else if (flag == "-maxRefineIter") {
            if (!parseInteger(value, options.maxRefineIterations)) {
                message = "-maxRefineIter requires an integer";
                return -15;
            }
        } else if (flag == "-denseMax") {
            if (!parseInteger(value, options.denseMax)) {
                message = "-denseMax requires an integer";
                return -16;
            }
        } else if (flag == "-assemblyRtol") {
            if (!parseReal(value, options.assemblyRtol)) {
                message = "-assemblyRtol requires a finite number";
                return -14;
            }
        } else if (flag == "-assemblyAtol") {
            if (!parseReal(value, options.assemblyAtol)) {
                message = "-assemblyAtol requires a finite number";
                return -15;
            }
        } else if (flag == "-massRtol") {
            if (!parseReal(value, options.massRtol)) {
                message = "-massRtol requires a finite number";
                return -16;
            }
        } else if (flag == "-massAtol") {
            if (!parseReal(value, options.massAtol)) {
                message = "-massAtol requires a finite number";
                return -17;
            }
        } else if (flag == "-verifyFullMaxBytes") {
            if (!parseSize(value, options.verifyFullMaxBytes)) {
                message = "-verifyFullMaxBytes requires a nonnegative integer";
                return -18;
            }
        } else if (flag == "-verifyAssembly") {
            verifyAssemblySpecified = true;
            if (value == "off")
                options.verifyAssembly = AssemblyVerification::Off;
            else if (value == "signature")
                options.verifyAssembly = AssemblyVerification::Signature;
            else if (value == "full")
                options.verifyAssembly = AssemblyVerification::Full;
            else {
                message = "-verifyAssembly accepts only off, signature or full";
                return -19;
            }
        } else if (flag == "-partition") {
            if (value != "metis") {
                message = "v1 supports only -partition metis";
                return -20;
            }
        } else if (flag == "-localEigen") {
            if (value != "lanczos") {
                message = "v1 supports only -localEigen lanczos";
                return -21;
            }
        } else if (flag == "-globalSolver") {
            if (value != "dense") {
                message = "v1 supports only -globalSolver dense";
                return -22;
            }
        } else {
            message = "unknown option: " + flag;
            return -23;
        }
    }
    if (!hierarchySpecified || !modesLevel2Specified || !modesLevel1Specified) {
        message = "-hierarchy, -modesL2 and -modesL1 are required";
        return -24;
    }
    if (options.hierarchy == HierarchyMode::Logical &&
        (!level1Specified || !level2Specified)) {
        message = "logical hierarchy requires -level1 and -level2";
        return -25;
    }
    if (options.hierarchy == HierarchyMode::Auto &&
        (level1Specified || level2Specified)) {
        message = "auto hierarchy derives level1/level2 and rejects explicit values";
        return -26;
    }
    if (options.domainMode == DomainMode::Physical) {
        if (!verifyAssemblySpecified)
            options.verifyAssembly = AssemblyVerification::Off;
        else if (options.verifyAssembly != AssemblyVerification::Off) {
            message = "physical domain mode requires -verifyAssembly off";
            return -27;
        }
    }
    return 0;
}

int Options::validate(int worldSize, int numModes, std::string &message) const
{
    message.clear();
    if (worldSize < 1) {
        message = "MPI world size must be positive";
        return -1;
    }
    if (numModes < 1) {
        message = "numModes must be positive";
        return -2;
    }
    if (hierarchy == HierarchyMode::Logical) {
        if (level1 < 1 || level2 < 1) {
            message = "logical hierarchy requires positive level1 and level2";
            return -3;
        }
        const long long product = static_cast<long long>(level1) * level2;
        if (product != worldSize) {
            message = "logical hierarchy requires level1*level2 == MPI world size";
            return -4;
        }
    }
    if (modesLevel2 < 1 || modesLevel1 < 1) {
        message = "modesL2 and modesL1 must be positive";
        return -5;
    }
    if (!(tolerance > 0.0) || maxEnrich < 0 || maxIterations < 1 ||
        maxRefineIterations < 0) {
        message = "invalid convergence controls";
        return -6;
    }
    const int initialVectors = resolvedIterationVectors(numModes, message);
    if (initialVectors < 0)
        return -7;
    if (denseMax < initialVectors) {
        message = "denseMax must be at least the number of CMS starting vectors";
        return -7;
    }
    if (assemblyRtol < 0.0 || assemblyAtol < 0.0 ||
        massRtol < 0.0 || massAtol < 0.0) {
        message = "assembly and mass tolerances must be nonnegative";
        return -8;
    }
    if (verifyAssembly == AssemblyVerification::Full && verifyFullMaxBytes == 0) {
        message = "full assembly verification requires a positive byte limit";
        return -9;
    }
    if (domainMode == DomainMode::Physical &&
        verifyAssembly != AssemblyVerification::Off) {
        message = "physical domain mode cannot verify replicated assembly";
        return -10;
    }
    return 0;
}

int Options::resolvedIterationVectors(int numModes, std::string &message) const
{
    message.clear();
    if (numModes < 1) {
        message = "numModes must be positive";
        return -1;
    }
    if (refinement == RefinementMode::None)
        return numModes;
    if (iterationVectors < 0) {
        message = "iterationVectors must be auto or a positive integer";
        return -2;
    }
    if (iterationVectors > 0) {
        if (iterationVectors < numModes) {
            message = "iterationVectors must be at least numModes";
            return -3;
        }
        return iterationVectors;
    }
    if (numModes > (std::numeric_limits<int>::max() - 8) ||
        numModes > std::numeric_limits<int>::max() / 2) {
        message = "automatic iteration-vector count exceeds the integer range";
        return -4;
    }
    return std::max(numModes + 8, 2 * numModes);
}

DenseMemoryEstimate estimateDenseRootMemory(
    const int *localDimensions,
    int numberOfGroups,
    int dynamicDimension,
    int masslessDimension,
    int numModes,
    int staticBlockSize)
{
    DenseMemoryEstimate result;
    if (numberOfGroups < 0 || dynamicDimension < 0 || masslessDimension < 0 ||
        numModes < 0 || staticBlockSize < 1 ||
        (numberOfGroups > 0 && localDimensions == nullptr)) {
        result.knownLowerBound = std::numeric_limits<std::size_t>::max();
        return result;
    }

    const std::size_t rD = static_cast<std::size_t>(dynamicDimension);
    const std::size_t rZ = static_cast<std::size_t>(masslessDimension);
    const std::size_t modes = static_cast<std::size_t>(numModes);
    const std::size_t block = static_cast<std::size_t>(
        dynamicDimension < staticBlockSize ? dynamicDimension : staticBlockSize);

    std::size_t square = 0;
    std::size_t twoMatrices = 0;
    std::size_t vectors = 0;
    std::size_t tile = 0;
    if (!checkedMultiply(rD, rD, square) ||
        !checkedMultiply(square, 2u * sizeof(double), twoMatrices) ||
        !checkedMultiply(rD, modes, vectors) ||
        !checkedMultiply(vectors, sizeof(double), vectors) ||
        !checkedMultiply(rZ, block, tile) ||
        !checkedMultiply(tile, sizeof(double), tile)) {
        result.knownLowerBound = std::numeric_limits<std::size_t>::max();
        return result;
    }
    if (!checkedAdd(twoMatrices, vectors, result.matricesAndVectors)) {
        result.knownLowerBound = std::numeric_limits<std::size_t>::max();
        return result;
    }
    result.staticTile = tile;

    std::size_t packedEntries = 0;
    for (int group = 0; group < numberOfGroups; ++group) {
        if (localDimensions[group] < 0) {
            result.knownLowerBound = std::numeric_limits<std::size_t>::max();
            return result;
        }
        const std::size_t dimension = static_cast<std::size_t>(localDimensions[group]);
        std::size_t product = 0;
        if (!checkedMultiply(dimension, dimension + 1u, product) ||
            !checkedAdd(packedEntries, product / 2u, packedEntries)) {
            result.knownLowerBound = std::numeric_limits<std::size_t>::max();
            return result;
        }
    }
    const std::size_t bytesPerEntry = 2u * sizeof(int) + 2u * sizeof(double);
    if (!checkedMultiply(packedEntries, bytesPerEntry, result.packedGather) ||
        !checkedAdd(result.matricesAndVectors, result.staticTile, result.knownLowerBound) ||
        !checkedAdd(result.knownLowerBound, result.packedGather, result.knownLowerBound)) {
        result.knownLowerBound = std::numeric_limits<std::size_t>::max();
    }
    return result;
}

AssemblyMetrics computeAssemblyMetrics(
    const double *rowMajorValues,
    int rows,
    int columns,
    double factor)
{
    AssemblyMetrics result;
    if (rows < 0 || columns < 0 || !std::isfinite(factor) ||
        ((rows > 0 && columns > 0) && rowMajorValues == nullptr)) {
        result.maxAbs = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    std::size_t count = 0;
    if (!checkedMultiply(
            static_cast<std::size_t>(rows),
            static_cast<std::size_t>(columns),
            count)) {
        result.maxAbs = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    result.numberOfValues = count;
    for (std::size_t index = 0; index < count; ++index) {
        const double value = factor * rowMajorValues[index];
        if (!std::isfinite(value)) {
            result.maxAbs = std::numeric_limits<double>::quiet_NaN();
            return result;
        }
        const double magnitude = std::fabs(value);
        result.maxAbs = std::max(result.maxAbs, magnitude);
        result.sum += value;
        result.sumAbs += magnitude;
        result.sumSquares += value * value;
        const double weight =
            (std::fmod(static_cast<double>(index) * 131.0 + 17.0, 257.0) + 1.0) /
            257.0;
        result.weightedProjection += weight * value;
    }
    return result;
}

int makeAssemblyRecord(
    ContributionKind kind,
    std::size_t ordinal,
    const double *rowMajorValues,
    int rows,
    int columns,
    const std::vector<int> &originalIDs,
    int globalDimension,
    double factor,
    double relativeTolerance,
    double absoluteTolerance,
    AssemblyRecord &result,
    std::string &message)
{
    message.clear();
    result = AssemblyRecord{};
    if (rows < 0 || columns < 0 || rows != columns ||
        originalIDs.size() != static_cast<std::size_t>(rows)) {
        message = "assembly matrix must be square and match its ID list";
        return -1;
    }
    if (globalDimension < 0 || relativeTolerance < 0.0 || absoluteTolerance < 0.0 ||
        !std::isfinite(factor)) {
        message = "invalid assembly dimension, factor, or tolerance";
        return -2;
    }

    result.kind = kind;
    result.ordinal = ordinal;
    result.rows = rows;
    result.columns = columns;
    result.originalIDs = originalIDs;
    result.metrics = computeAssemblyMetrics(rowMajorValues, rows, columns, factor);
    if (!std::isfinite(result.metrics.maxAbs)) {
        message = "assembly contains a non-finite value or invalid storage size";
        return -3;
    }

    std::vector<int> activePositions;
    activePositions.reserve(originalIDs.size());
    result.activeIDs.reserve(originalIDs.size());
    for (int local = 0; local < rows; ++local) {
        const int equation = originalIDs[static_cast<std::size_t>(local)];
        if (equation < 0)
            continue;
        if (equation >= globalDimension) {
            message = "assembly equation is outside the global range";
            return -4;
        }
        activePositions.push_back(local);
        result.activeIDs.push_back(equation);
    }

    std::size_t packedSize = 0;
    if (!checkedMultiply(activePositions.size(), activePositions.size() + 1u, packedSize)) {
        message = "active assembly order overflows packed storage";
        return -5;
    }
    result.upperValues.reserve(packedSize / 2u);
    for (std::size_t column = 0; column < activePositions.size(); ++column) {
        const int localColumn = activePositions[column];
        for (std::size_t row = 0; row <= column; ++row) {
            const int localRow = activePositions[row];
            const double upper = factor * rowMajorValues[
                static_cast<std::size_t>(localRow) * columns + localColumn];
            const double lower = factor * rowMajorValues[
                static_cast<std::size_t>(localColumn) * columns + localRow];
            const double scale = std::max(1.0, std::max(std::fabs(upper), std::fabs(lower)));
            const double tolerance = absoluteTolerance + relativeTolerance * scale;
            if (!std::isfinite(upper) || !std::isfinite(lower) ||
                std::fabs(upper - lower) > tolerance) {
                message = "active assembly contribution is non-symmetric or non-finite";
                return -6;
            }
            result.upperValues.push_back(0.5 * (upper + lower));
        }
    }
    return 0;
}

bool compareAssemblyRecords(
    const AssemblyRecord &reference,
    const AssemblyRecord &candidate,
    double relativeTolerance,
    double absoluteTolerance,
    std::string &message)
{
    message.clear();
    if (relativeTolerance < 0.0 || absoluteTolerance < 0.0) {
        message = "assembly tolerances must be nonnegative";
        return false;
    }
    if (reference.ordinal != candidate.ordinal ||
        reference.kind != candidate.kind ||
        reference.rows != candidate.rows ||
        reference.columns != candidate.columns ||
        reference.originalIDs != candidate.originalIDs ||
        reference.metrics.numberOfValues != candidate.metrics.numberOfValues) {
        message = "assembly structure differs";
        return false;
    }

    const double refMetrics[] = {
        reference.metrics.maxAbs,
        reference.metrics.sum,
        reference.metrics.sumAbs,
        reference.metrics.sumSquares,
        reference.metrics.weightedProjection};
    const double candidateMetrics[] = {
        candidate.metrics.maxAbs,
        candidate.metrics.sum,
        candidate.metrics.sumAbs,
        candidate.metrics.sumSquares,
        candidate.metrics.weightedProjection};
    const double countScale = static_cast<double>(
        std::max<std::size_t>(1u, reference.metrics.numberOfValues));
    for (int metric = 0; metric < 5; ++metric) {
        if (!std::isfinite(refMetrics[metric]) || !std::isfinite(candidateMetrics[metric])) {
            message = "assembly metric is non-finite";
            return false;
        }
        const double scale = std::max(
            1.0, std::max(std::fabs(refMetrics[metric]), std::fabs(candidateMetrics[metric])));
        const double tolerance = absoluteTolerance * countScale + relativeTolerance * scale;
        if (std::fabs(refMetrics[metric] - candidateMetrics[metric]) > tolerance) {
            message = "assembly numerical metrics differ";
            return false;
        }
    }
    return true;
}

int assignContributionOwner(
    const std::vector<int> &activeIDs,
    const std::vector<int> &fineLabels,
    std::string &message)
{
    message.clear();
    if (activeIDs.empty()) {
        message = "a contribution needs at least one active equation";
        return -1;
    }
    std::map<int, int> counts;
    for (int equation : activeIDs) {
        if (equation < 0 || static_cast<std::size_t>(equation) >= fineLabels.size()) {
            message = "active equation is outside the fine partition map";
            return -1;
        }
        const int owner = fineLabels[static_cast<std::size_t>(equation)];
        if (owner < 0) {
            message = "fine partition labels must be nonnegative";
            return -1;
        }
        ++counts[owner];
    }
    int selected = -1;
    int selectedCount = -1;
    for (const auto &entry : counts) {
        if (entry.second > selectedCount) {
            selected = entry.first;
            selectedCount = entry.second;
        }
    }
    return selected;
}

int classifyEffectiveOwnership(
    int numberOfEquations,
    const std::vector<AssemblyRecord> &contributions,
    const std::vector<int> &contributionOwners,
    const std::vector<int> &coarseOfFine,
    EffectiveOwnership &result,
    std::string &message)
{
    message.clear();
    result = EffectiveOwnership{};
    if (numberOfEquations < 0 || contributions.size() != contributionOwners.size() ||
        coarseOfFine.empty()) {
        message = "invalid effective-ownership dimensions";
        return -1;
    }
    std::vector<std::set<int>> ownerSets(static_cast<std::size_t>(numberOfEquations));
    for (std::size_t index = 0; index < contributions.size(); ++index) {
        const int owner = contributionOwners[index];
        if (owner < 0 || static_cast<std::size_t>(owner) >= coarseOfFine.size()) {
            message = "contribution owner is outside coarseOfFine";
            return -2;
        }
        std::set<int> uniqueEquations;
        for (int equation : contributions[index].activeIDs) {
            if (equation < 0 || equation >= numberOfEquations) {
                message = "contribution equation is outside the global range";
                return -3;
            }
            uniqueEquations.insert(equation);
        }
        for (int equation : uniqueEquations)
            ownerSets[static_cast<std::size_t>(equation)].insert(owner);
    }

    result.ownersByEquation.resize(static_cast<std::size_t>(numberOfEquations));
    result.interiorOwner.assign(static_cast<std::size_t>(numberOfEquations), -1);
    result.role.assign(static_cast<std::size_t>(numberOfEquations), EquationRole::Unseen);
    for (int equation = 0; equation < numberOfEquations; ++equation) {
        const std::set<int> &owners = ownerSets[static_cast<std::size_t>(equation)];
        std::vector<int> &stored = result.ownersByEquation[static_cast<std::size_t>(equation)];
        stored.assign(owners.begin(), owners.end());
        if (owners.empty())
            continue;
        if (owners.size() == 1u) {
            result.role[static_cast<std::size_t>(equation)] = EquationRole::Interior;
            result.interiorOwner[static_cast<std::size_t>(equation)] = *owners.begin();
            continue;
        }
        std::set<int> coarseOwners;
        for (int owner : owners) {
            const int coarse = coarseOfFine[static_cast<std::size_t>(owner)];
            if (coarse < 0) {
                message = "coarse labels must be nonnegative";
                return -4;
            }
            coarseOwners.insert(coarse);
        }
        result.role[static_cast<std::size_t>(equation)] = coarseOwners.size() > 1u
            ? EquationRole::Level1Interface
            : EquationRole::Level2Interface;
    }
    return 0;
}

int SymmetricCSR::validate(std::string &message) const
{
    message.clear();
    const std::size_t expectedRows = dimension < 0
        ? 0u
        : static_cast<std::size_t>(dimension) + 1u;
    if (dimension < 0 || rowOffsets.size() != expectedRows ||
        rowOffsets.empty() || rowOffsets.front() != 0) {
        message = "invalid CSR dimensions or row offsets";
        return -1;
    }
    if (columnIndices.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        columnIndices.size() != upperValues.size() ||
        rowOffsets.back() != static_cast<int>(columnIndices.size())) {
        message = "CSR index and value counts differ";
        return -2;
    }
    for (int row = 0; row < dimension; ++row) {
        const int begin = rowOffsets[static_cast<std::size_t>(row)];
        const int end = rowOffsets[static_cast<std::size_t>(row + 1)];
        if (begin < 0 || end < begin || end > static_cast<int>(columnIndices.size())) {
            message = "CSR row offsets are not monotone";
            return -3;
        }
        int previous = row - 1;
        for (int position = begin; position < end; ++position) {
            const int column = columnIndices[static_cast<std::size_t>(position)];
            if (column < row || column >= dimension || column <= previous ||
                !std::isfinite(upperValues[static_cast<std::size_t>(position)])) {
                message = "CSR upper triangle is unsorted or invalid";
                return -4;
            }
            previous = column;
        }
    }
    return 0;
}

std::vector<double> SymmetricCSR::toDense() const
{
    std::string message;
    if (validate(message) < 0)
        return {};
    std::vector<double> dense(
        static_cast<std::size_t>(dimension) * static_cast<std::size_t>(dimension), 0.0);
    for (int row = 0; row < dimension; ++row) {
        for (int position = rowOffsets[static_cast<std::size_t>(row)];
             position < rowOffsets[static_cast<std::size_t>(row + 1)]; ++position) {
            const int column = columnIndices[static_cast<std::size_t>(position)];
            const double value = upperValues[static_cast<std::size_t>(position)];
            dense[static_cast<std::size_t>(row) * dimension + column] += value;
            if (column != row)
                dense[static_cast<std::size_t>(column) * dimension + row] += value;
        }
    }
    return dense;
}

std::vector<double> SymmetricCSR::multiply(const std::vector<double> &x) const
{
    std::string message;
    if (x.size() != static_cast<std::size_t>(dimension) || validate(message) < 0)
        return {};
    std::vector<double> result(static_cast<std::size_t>(dimension), 0.0);
    for (int row = 0; row < dimension; ++row) {
        for (int position = rowOffsets[static_cast<std::size_t>(row)];
             position < rowOffsets[static_cast<std::size_t>(row + 1)]; ++position) {
            const int column = columnIndices[static_cast<std::size_t>(position)];
            const double value = upperValues[static_cast<std::size_t>(position)];
            result[static_cast<std::size_t>(row)] += value * x[static_cast<std::size_t>(column)];
            if (column != row)
                result[static_cast<std::size_t>(column)] += value * x[static_cast<std::size_t>(row)];
        }
    }
    return result;
}

int buildSymmetricCSR(
    int dimension,
    const std::vector<AssemblyRecord> &contributions,
    ContributionKind kind,
    SymmetricCSR &result,
    std::string &message)
{
    message.clear();
    result = SymmetricCSR{};
    if (dimension < 0) {
        message = "CSR dimension must be nonnegative";
        return -1;
    }
    std::map<std::pair<int, int>, double> entries;
    for (const AssemblyRecord &contribution : contributions) {
        if (contribution.kind != kind)
            continue;
        const std::size_t order = contribution.activeIDs.size();
        std::size_t packedSize = 0;
        if (!checkedMultiply(order, order + 1u, packedSize) ||
            contribution.upperValues.size() != packedSize / 2u) {
            message = "packed contribution size is invalid";
            return -2;
        }
        std::size_t packed = 0;
        for (std::size_t column = 0; column < order; ++column) {
            const int globalColumn = contribution.activeIDs[column];
            if (globalColumn < 0 || globalColumn >= dimension) {
                message = "contribution column is outside CSR dimension";
                return -3;
            }
            for (std::size_t row = 0; row <= column; ++row, ++packed) {
                const int globalRow = contribution.activeIDs[row];
                const double value = contribution.upperValues[packed];
                if (globalRow < 0 || globalRow >= dimension || !std::isfinite(value)) {
                    message = "contribution row or value is invalid";
                    return -4;
                }
                const int first = std::min(globalRow, globalColumn);
                const int second = std::max(globalRow, globalColumn);
                entries[std::make_pair(first, second)] +=
                    row != column && globalRow == globalColumn ? 2.0 * value : value;
            }
        }
    }

    result.dimension = dimension;
    if (entries.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        message = "CSR has more entries than its integer index type can represent";
        return -5;
    }
    result.rowOffsets.assign(static_cast<std::size_t>(dimension) + 1u, 0);
    for (const auto &entry : entries)
        ++result.rowOffsets[static_cast<std::size_t>(entry.first.first + 1)];
    for (int row = 0; row < dimension; ++row)
        result.rowOffsets[static_cast<std::size_t>(row + 1)] +=
            result.rowOffsets[static_cast<std::size_t>(row)];
    result.columnIndices.reserve(entries.size());
    result.upperValues.reserve(entries.size());
    for (const auto &entry : entries) {
        result.columnIndices.push_back(entry.first.second);
        result.upperValues.push_back(entry.second);
    }
    return result.validate(message);
}

int buildOwnedLocalPencil(
    const std::vector<AssemblyRecord> &stiffnessContributions,
    const std::vector<AssemblyRecord> &massContributions,
    std::vector<int> &globalEquations,
    SymmetricCSR &stiffness,
    SymmetricCSR &mass,
    std::string &message)
{
    message.clear();
    globalEquations.clear();
    stiffness = SymmetricCSR{};
    mass = SymmetricCSR{};
    std::set<int> equations;
    for (const auto &record : stiffnessContributions)
        equations.insert(record.activeIDs.begin(), record.activeIDs.end());
    for (const auto &record : massContributions)
        equations.insert(record.activeIDs.begin(), record.activeIDs.end());
    if (equations.empty() || *equations.begin() < 0) {
        message = "an owned subdomain has no valid active equations";
        return -1;
    }
    globalEquations.assign(equations.begin(), equations.end());
    std::map<int, int> localOfGlobal;
    for (std::size_t local = 0; local < globalEquations.size(); ++local)
        localOfGlobal[globalEquations[local]] = static_cast<int>(local);

    const auto remap = [&localOfGlobal](
                           const std::vector<AssemblyRecord> &source,
                           std::vector<AssemblyRecord> &target) -> bool {
        target = source;
        for (AssemblyRecord &record : target) {
            for (int &equation : record.activeIDs) {
                const auto found = localOfGlobal.find(equation);
                if (found == localOfGlobal.end())
                    return false;
                equation = found->second;
            }
        }
        return true;
    };
    std::vector<AssemblyRecord> localStiffness;
    std::vector<AssemblyRecord> localMass;
    if (!remap(stiffnessContributions, localStiffness) ||
        !remap(massContributions, localMass)) {
        message = "owned contribution remapping failed";
        return -2;
    }
    const int dimension = static_cast<int>(globalEquations.size());
    if (buildSymmetricCSR(
            dimension, localStiffness, ContributionKind::Stiffness,
            stiffness, message) < 0) {
        message = "owned stiffness assembly failed: " + message;
        return -3;
    }
    if (buildSymmetricCSR(
            dimension, localMass, ContributionKind::Mass,
            mass, message) < 0) {
        message = "owned mass assembly failed: " + message;
        return -4;
    }
    return 0;
}

int classifyCoordinateMass(
    const SymmetricCSR &mass,
    double relativeTolerance,
    double absoluteTolerance,
    CoordinateMassClassification &result,
    std::string &message)
{
    message.clear();
    result = CoordinateMassClassification{};
    if (relativeTolerance < 0.0 || absoluteTolerance < 0.0 || mass.validate(message) < 0) {
        if (message.empty())
            message = "mass tolerances must be nonnegative";
        return -1;
    }
    double scale = 1.0;
    for (double value : mass.upperValues)
        scale = std::max(scale, std::fabs(value));
    result.threshold = absoluteTolerance + relativeTolerance * scale;

    std::vector<double> diagonal(static_cast<std::size_t>(mass.dimension), 0.0);
    for (int row = 0; row < mass.dimension; ++row) {
        for (int position = mass.rowOffsets[static_cast<std::size_t>(row)];
             position < mass.rowOffsets[static_cast<std::size_t>(row + 1)]; ++position) {
            if (mass.columnIndices[static_cast<std::size_t>(position)] == row)
                diagonal[static_cast<std::size_t>(row)] +=
                    mass.upperValues[static_cast<std::size_t>(position)];
        }
    }
    std::vector<bool> isMassless(static_cast<std::size_t>(mass.dimension), false);
    for (int equation = 0; equation < mass.dimension; ++equation) {
        const double value = diagonal[static_cast<std::size_t>(equation)];
        if (value < -result.threshold) {
            message = "mass has a negative diagonal coordinate";
            return -2;
        }
        if (std::fabs(value) <= result.threshold) {
            isMassless[static_cast<std::size_t>(equation)] = true;
            result.massless.push_back(equation);
        } else {
            result.dynamic.push_back(equation);
        }
    }
    if (result.dynamic.empty()) {
        message = "mass has no dynamic coordinates";
        return -3;
    }
    for (int row = 0; row < mass.dimension; ++row) {
        for (int position = mass.rowOffsets[static_cast<std::size_t>(row)];
             position < mass.rowOffsets[static_cast<std::size_t>(row + 1)]; ++position) {
            const int column = mass.columnIndices[static_cast<std::size_t>(position)];
            const double value = mass.upperValues[static_cast<std::size_t>(position)];
            if ((isMassless[static_cast<std::size_t>(row)] ||
                 isMassless[static_cast<std::size_t>(column)]) &&
                std::fabs(value) > result.threshold) {
                message = "mass nullspace is not coordinate aligned";
                return -4;
            }
        }
    }
    return 0;
}

} // namespace ladruno_cms
